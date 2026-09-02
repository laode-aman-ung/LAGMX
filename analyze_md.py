#!/usr/bin/env python3
"""LAGMX post-production analysis.

LAGMX takes a protein-ligand system from PDB to a finished production
trajectory. What it has never done is tell you anything about that trajectory:
every project ended with an .xtc and a fresh round of hand-typed `gmx rms`,
`gmx sasa`, `gmx covar` invocations, each with its own interactive group
prompts and its own chance of fitting on the wrong group.

This script closes that gap. It runs the standard analysis battery over every
complex_*/ directory that has a finished production run, writes CSV and PNG
per analysis, and reduces each complex to one row of summary numbers so that
several systems can be compared side by side.

    cd run_matrix && python3 ../analyze_md.py

Like LAGMX.py, every path is resolved against the directory you run it from,
not against the location of this file.

Analyses
    rmsd      complex stability, protein backbone and ligand separately
    rmsf      per-residue flexibility
    rg        radius of gyration, protein compactness
    sasa      solvent accessible surface, protein and complex
    hbond     protein-ligand hydrogen bonds over time
    contacts  per-residue protein-ligand contact occupancy
    pca       essential dynamics, PC1/PC2 projection
    fel       free energy landscape over PC1/PC2
    mmpbsa    binding free energy via gmx_MMPBSA (MM/GBSA and/or MM/PBSA)

Configuration comes from gmx_config.txt; see ANALYSIS_DEFAULTS below for the
keys this script adds. All of them are optional.
"""

import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Keys this script adds to gmx_config.txt. Everything has a working default so
# that an existing config file keeps working untouched.
ANALYSIS_DEFAULTS = {
    "analysis": "rmsd,rmsf,rg,sasa,hbond,contacts,pca,fel",
    "analysis_skip_ns": "0",          # drop this much from the start of the run
    "analysis_contact_cutoff": "0.4",  # nm, protein-ligand contact definition
    "analysis_gmx": "gmx",            # which gmx to analyse with
    "mmpbsa_python": "",              # env holding gmx_MMPBSA; empty = disabled
    "mmpbsa_method": "gb",            # gb, pb, or both
    "mmpbsa_frames": "100",           # frames sampled from the production run
    "mmpbsa_igb": "5",
    "mmpbsa_salt": "0.150",
}

ALL_ANALYSES = ["rmsd", "rmsf", "rg", "sasa", "hbond", "contacts", "pca", "fel", "mmpbsa"]

GMX = "gmx"


# --------------------------------------------------------------------------
# process helpers
# --------------------------------------------------------------------------

def gmx_run(args, stdin=None, cwd=None, quiet=True):
    """Call a gmx module with an argv list rather than a shell string.

    Group selection goes in on stdin as plain group numbers, which is what the
    interactive prompts read. Returns (ok, combined_output) instead of raising:
    one analysis failing must not take the other eight down with it.
    """
    cmd = [GMX] + [str(a) for a in args]
    env = dict(os.environ, GMX_MAXBACKUP="-1")
    try:
        proc = subprocess.run(
            cmd, input=stdin, cwd=cwd, capture_output=True, text=True,
            timeout=7200, env=env
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, f"{exc}"
    out = (proc.stdout or "") + (proc.stderr or "")
    if not quiet:
        print(out)
    return proc.returncode == 0, out


def say(msg, indent=0):
    print(f"{' ' * indent}{msg}", flush=True)


# --------------------------------------------------------------------------
# file format readers
# --------------------------------------------------------------------------

def read_xvg(path):
    """Return (data, legends) from a Grace .xvg file.

    Legends are pulled from the @ s0 legend lines so that plots and CSV
    headers carry the names GROMACS chose, instead of "column 1".
    """
    rows, legends, xlabel, ylabel = [], [], "", ""
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith("@"):
                m = re.match(r'@\s+s\d+\s+legend\s+"(.*)"', line)
                if m:
                    legends.append(m.group(1))
                m = re.match(r'@\s+xaxis\s+label\s+"(.*)"', line)
                if m:
                    xlabel = m.group(1)
                m = re.match(r'@\s+yaxis\s+label\s+"(.*)"', line)
                if m:
                    ylabel = m.group(1)
                continue
            if line.startswith(("#", "&")):
                continue
            parts = line.split()
            if not parts:
                continue
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    if not rows:
        return np.empty((0, 0)), {"legends": legends, "xlabel": xlabel, "ylabel": ylabel}
    width = min(len(r) for r in rows)
    data = np.array([r[:width] for r in rows])
    return data, {"legends": legends, "xlabel": xlabel, "ylabel": ylabel}


def read_xpm(path):
    """Parse a GROMACS .xpm matrix into (values, x_axis, y_axis).

    gmx sham writes the free energy landscape as an XPM colour map. The kJ/mol
    value of each colour lives in the C comment after the colour definition --

        "A  c #000000 " /* "0" */,

    -- not inside the quoted colour string itself, so the file has to be
    decoded rather than read as an image. The header line gives the geometry
    and, crucially, how many characters encode one pixel.
    """
    header_re = re.compile(r'^"\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*"')
    nx = ny = ncolours = nchar = None
    colours, rows, xaxis, yaxis = {}, [], [], []

    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("/* x-axis:"):
                xaxis += [float(v) for v in re.findall(r"[-\d.eE+]+", line[10:])]
                continue
            if line.startswith("/* y-axis:"):
                yaxis += [float(v) for v in re.findall(r"[-\d.eE+]+", line[10:])]
                continue
            if not line.startswith('"'):
                continue

            if nx is None:
                m = header_re.match(line)
                if m:
                    nx, ny, ncolours, nchar = (int(g) for g in m.groups())
                continue

            if len(colours) < ncolours:
                m = re.match(r'^"(.{%d})\s+c\s+\S+\s*"\s*/\*\s*"([^"]*)"' % nchar, line)
                if m:
                    try:
                        colours[m.group(1)] = float(m.group(2))
                    except ValueError:
                        colours[m.group(1)] = np.nan
                continue

            body = line[1:line.rfind('"')]
            if len(body) >= nx * nchar:
                rows.append(body)

    if not rows or not colours:
        return None, None, None
    values = np.array([
        [colours.get(row[i:i + nchar], np.nan) for i in range(0, nx * nchar, nchar)]
        for row in rows
    ])
    # XPM rows run top to bottom, the y axis runs bottom to top.
    values = values[::-1]
    return values, np.array(xaxis), np.array(yaxis)


# --------------------------------------------------------------------------
# index groups
# --------------------------------------------------------------------------

def build_index(tpr, out_ndx, workdir, merged_group="Protein_LIG"):
    """Write an index file and return {group name: group number}.

    The default make_ndx groups already carry everything needed: "Other" is
    exactly the set of non-protein, non-water, non-ion residues, which for a
    LAGMX system is exactly the ligands, however many there are and whatever
    they are called. A merged Protein+ligand group is appended for centring
    and for MM/PBSA.
    """
    ok, out = gmx_run(["make_ndx", "-f", tpr, "-o", out_ndx], stdin="q\n", cwd=workdir)
    if not ok:
        return None
    groups = {name: int(num) for num, name in re.findall(r"^\s*(\d+)\s+(\S+)\s*:", out, re.M)}
    if "Other" not in groups:
        say("no 'Other' group: system has no ligand, ligand analyses will be skipped", 6)
    else:
        new_index = max(groups.values()) + 1
        ok, out2 = gmx_run(
            ["make_ndx", "-f", tpr, "-n", out_ndx, "-o", out_ndx],
            stdin=f'"Protein" | "Other"\nname {new_index} {merged_group}\nq\n', cwd=workdir,
        )
        if ok:
            groups[merged_group] = new_index
    return groups


# --------------------------------------------------------------------------
# trajectory preparation
# --------------------------------------------------------------------------

def prepare_trajectory(cdir, adir, tpr, xtc, groups, merged_group, skip_ps):
    """Undo periodic boundary artefacts before anything is measured.

    Three passes, in this order and no other: make molecules whole, stop them
    jumping between images, then centre the complex in a compact box. Skipping
    this is the classic way to get an RMSD trace with a cliff in it that looks
    like unbinding and is really the ligand crossing the box edge.
    """
    system = groups.get("System", 0)
    centre = groups.get(merged_group, groups.get("Protein", 1))

    whole = os.path.join(adir, "_whole.xtc")
    nojump = os.path.join(adir, "_nojump.xtc")
    final = os.path.join(adir, "md_center.xtc")

    steps = [
        (["trjconv", "-s", tpr, "-f", xtc, "-o", whole, "-pbc", "whole"], f"{system}\n"),
        (["trjconv", "-s", tpr, "-f", whole, "-o", nojump, "-pbc", "nojump"], f"{system}\n"),
        (["trjconv", "-s", tpr, "-f", nojump, "-o", final,
          "-pbc", "mol", "-ur", "compact", "-center"], f"{centre}\n{system}\n"),
    ]
    if skip_ps > 0:
        steps[-1][0].extend(["-b", str(skip_ps)])

    for args, stdin in steps:
        ok, out = gmx_run(args + ["-n", os.path.join(adir, "analysis.ndx")],
                          stdin=stdin, cwd=cdir)
        if not ok:
            say(f"trjconv failed: {out.strip().splitlines()[-1] if out.strip() else '?'}", 6)
            return None

    for tmp in (whole, nojump):
        if os.path.exists(tmp):
            os.remove(tmp)

    # A single reference frame, used as the -s for analyses and for viewing.
    gmx_run(["trjconv", "-s", tpr, "-f", final, "-o", os.path.join(adir, "start.pdb"),
             "-n", os.path.join(adir, "analysis.ndx"), "-dump", "0"],
            stdin=f"{system}\n", cwd=cdir)
    return final


# --------------------------------------------------------------------------
# plotting and export
# --------------------------------------------------------------------------

def export(data, meta, adir, name, title, xlabel=None, ylabel=None, columns=None):
    """Write one analysis to CSV and PNG."""
    if data.size == 0:
        return
    csv_path = os.path.join(adir, f"{name}.csv")
    header = columns or ([meta.get("xlabel") or "x"] +
                         (meta.get("legends") or
                          [f"y{i}" for i in range(1, data.shape[1])]))
    header = header[:data.shape[1]]
    while len(header) < data.shape[1]:
        header.append(f"y{len(header)}")
    np.savetxt(csv_path, data, delimiter=",", header=",".join(header), comments="", fmt="%.6g")

    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    for col in range(1, data.shape[1]):
        ax.plot(data[:, 0], data[:, col], lw=1.0, label=header[col])
    ax.set_xlabel(xlabel or meta.get("xlabel") or "")
    ax.set_ylabel(ylabel or meta.get("ylabel") or "")
    ax.set_title(title)
    if data.shape[1] > 2:
        ax.legend(fontsize=8)
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(adir, f"{name}.png"), dpi=150)
    plt.close(fig)


def stat(data, col=1):
    """Mean and standard deviation of one column, ignoring empty input."""
    if data.size == 0 or data.shape[1] <= col:
        return None, None
    values = data[:, col]
    return float(np.mean(values)), float(np.std(values))


# --------------------------------------------------------------------------
# individual analyses
# --------------------------------------------------------------------------

def analyse_rmsd(ctx):
    """Backbone RMSD for stability, ligand RMSD for whether it stayed put.

    The ligand curve is deliberately fitted on the protein backbone, not on
    the ligand itself: fitting a ligand to its own reference measures internal
    conformational change and hides the thing you actually want to see, which
    is the ligand drifting out of the pocket.
    """
    results = {}
    backbone = ctx["groups"].get("Backbone")
    if backbone is None:
        return results

    ok, _ = gmx_run(["rms", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                     "-o", ctx["a"]("rmsd_protein.xvg"), "-tu", "ns"],
                    stdin=f"{backbone}\n{backbone}\n", cwd=ctx["cdir"])
    if ok:
        data, meta = read_xvg(ctx["a"]("rmsd_protein.xvg"))
        export(data, meta, ctx["adir"], "rmsd_protein", "RMSD protein backbone",
               "waktu (ns)", "RMSD (nm)", ["waktu_ns", "rmsd_nm"])
        mean, sd = stat(data)
        results["rmsd_protein_nm"], results["rmsd_protein_sd"] = mean, sd

    ligand = ctx["groups"].get("Other")
    if ligand is not None:
        ok, _ = gmx_run(["rms", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                         "-o", ctx["a"]("rmsd_ligand.xvg"), "-tu", "ns"],
                        stdin=f"{backbone}\n{ligand}\n", cwd=ctx["cdir"])
        if ok:
            data, meta = read_xvg(ctx["a"]("rmsd_ligand.xvg"))
            export(data, meta, ctx["adir"], "rmsd_ligand", "RMSD ligand (fit on protein)",
                   "waktu (ns)", "RMSD (nm)", ["waktu_ns", "rmsd_nm"])
            mean, sd = stat(data)
            results["rmsd_ligand_nm"], results["rmsd_ligand_sd"] = mean, sd
    return results


def analyse_rmsf(ctx):
    """Per-residue fluctuation, C-alpha only, averaged over each residue."""
    ca = ctx["groups"].get("C-alpha")
    if ca is None:
        return {}
    ok, _ = gmx_run(["rmsf", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                     "-o", ctx["a"]("rmsf.xvg"), "-oq", ctx["a"]("rmsf_bfactor.pdb"), "-res"],
                    stdin=f"{ca}\n", cwd=ctx["cdir"])
    if not ok:
        return {}
    data, meta = read_xvg(ctx["a"]("rmsf.xvg"))
    export(data, meta, ctx["adir"], "rmsf", "RMSF per residu",
           "residu", "RMSF (nm)", ["residu", "rmsf_nm"])
    mean, _ = stat(data)
    out = {"rmsf_mean_nm": mean}
    if data.size:
        top = data[np.argsort(-data[:, 1])][:10]
        out["rmsf_residu_teratas"] = " ".join(f"{int(r[0])}:{r[1]:.2f}" for r in top)
    return out


def analyse_rg(ctx):
    """Radius of gyration. Falls back to the legacy tool: gmx gyrate was
    reimplemented in recent releases and older builds only ship one name."""
    protein = ctx["groups"].get("Protein")
    if protein is None:
        return {}
    for tool in ("gyrate", "gyrate-legacy"):
        ok, _ = gmx_run([tool, "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                         "-o", ctx["a"]("rg.xvg")], stdin=f"{protein}\n", cwd=ctx["cdir"])
        if ok:
            break
    else:
        return {}
    data, meta = read_xvg(ctx["a"]("rg.xvg"))
    if data.size:
        data = data[:, :2]
    export(data, meta, ctx["adir"], "rg", "Radius of gyration",
           "waktu (ps)", "Rg (nm)", ["waktu_ps", "rg_nm"])
    mean, sd = stat(data)
    return {"rg_nm": mean, "rg_sd": sd}


def analyse_sasa(ctx):
    """Solvent accessible surface of the protein and of the whole complex.

    The difference between the two, against the ligand's own free surface, is
    the buried area -- a cheap proxy for how deep the ligand sits.
    """
    out = {}
    jobs = [("sasa_protein", 'group "Protein"')]
    if "Other" in ctx["groups"]:
        jobs.append(("sasa_ligand", 'group "Other"'))
        jobs.append(("sasa_complex", 'group "Protein" or group "Other"'))
    for name, sel in jobs:
        ok, _ = gmx_run(["sasa", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                         "-surface", sel, "-o", ctx["a"](f"{name}.xvg")], cwd=ctx["cdir"])
        if not ok:
            continue
        data, meta = read_xvg(ctx["a"](f"{name}.xvg"))
        if data.size:
            data = data[:, :2]
        export(data, meta, ctx["adir"], name, name.replace("_", " ").upper(),
               "waktu (ps)", "SASA (nm^2)", ["waktu_ps", "sasa_nm2"])
        mean, sd = stat(data)
        out[f"{name}_nm2"], out[f"{name}_sd"] = mean, sd

    if all(k in out for k in ("sasa_protein_nm2", "sasa_ligand_nm2", "sasa_complex_nm2")):
        out["sasa_terkubur_nm2"] = (out["sasa_protein_nm2"] + out["sasa_ligand_nm2"]
                                    - out["sasa_complex_nm2"])
    return out


def analyse_hbond(ctx):
    """Protein-ligand hydrogen bonds per frame.

    gmx hbond was rewritten with a selection interface in recent GROMACS, while
    older builds prompt for two index groups. Try the modern call first, fall
    back to the legacy one, so this works on both.
    """
    if "Other" not in ctx["groups"]:
        return {}
    target = ctx["a"]("hbond_num.xvg")

    ok, _ = gmx_run(["hbond", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                     "-ref", 'group "Protein"', "-sel", 'group "Other"',
                     "-num", target], cwd=ctx["cdir"])
    if not ok or not os.path.exists(target):
        protein, ligand = ctx["groups"]["Protein"], ctx["groups"]["Other"]
        for tool in ("hbond", "hbond-legacy"):
            ok, _ = gmx_run([tool, "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                             "-num", target], stdin=f"{protein}\n{ligand}\n", cwd=ctx["cdir"])
            if ok and os.path.exists(target):
                break
    if not os.path.exists(target):
        return {}

    data, meta = read_xvg(target)
    if data.size:
        data = data[:, :2]
    export(data, meta, ctx["adir"], "hbond", "Ikatan hidrogen protein-ligan",
           "waktu (ps)", "jumlah H-bond", ["waktu_ps", "n_hbond"])
    mean, sd = stat(data)
    out = {"hbond_rerata": mean, "hbond_sd": sd}
    if data.size:
        out["hbond_maks"] = float(np.max(data[:, 1]))
        out["hbond_frac_ada"] = float(np.mean(data[:, 1] > 0))
    return out


def _atom_residue_map(pdb_path):
    """Map 1-based atom serial to (chain, resid, resname) from a dumped frame."""
    mapping, serial = {}, 0
    with open(pdb_path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                serial += 1
                mapping[serial] = (line[21].strip() or "-", line[22:26].strip(),
                                   line[17:20].strip())
    return mapping


def analyse_contacts(ctx):
    """Which residues actually touch the ligand, and for what fraction of the run.

    A single contact map from the final frame says nothing about persistence.
    What matters for picking key residues is occupancy: the share of frames in
    which a residue sits within the cutoff of the ligand.
    """
    if "Other" not in ctx["groups"]:
        return {}
    cutoff = ctx["contact_cutoff"]
    sel = (f'group "Protein" and same residue as within {cutoff} of group "Other"')
    size_xvg, idx_dat = ctx["a"]("contacts_size.xvg"), ctx["a"]("contacts_index.dat")

    ok, _ = gmx_run(["select", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                     "-select", sel, "-os", size_xvg, "-oi", idx_dat], cwd=ctx["cdir"])
    if not ok or not os.path.exists(idx_dat):
        return {}

    data, meta = read_xvg(size_xvg)
    if data.size:
        data = data[:, :2]
    export(data, meta, ctx["adir"], "contacts_count",
           f"Atom protein dalam {cutoff} nm dari ligan",
           "waktu (ps)", "jumlah atom", ["waktu_ps", "n_atom"])

    amap = _atom_residue_map(ctx["a"]("start.pdb")) if os.path.exists(ctx["a"]("start.pdb")) else {}
    counts, frames = {}, 0
    with open(idx_dat, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            frames += 1
            seen = set()
            for token in parts[2:]:
                try:
                    key = amap.get(int(token))
                except ValueError:
                    continue
                if key and key not in seen:
                    seen.add(key)
                    counts[key] = counts.get(key, 0) + 1

    if not frames or not counts:
        return {}
    rows = sorted(counts.items(), key=lambda kv: -kv[1])
    with open(ctx["a"]("contacts_residue.csv"), "w") as fh:
        fh.write("rantai,nomor,residu,frame_kontak,okupansi\n")
        for (chain, resid, resname), n in rows:
            fh.write(f"{chain},{resid},{resname},{n},{n / frames:.4f}\n")

    top = rows[:15][::-1]
    fig, ax = plt.subplots(figsize=(7.5, max(3.0, 0.28 * len(top))))
    ax.barh([f"{r[0][2]}{r[0][1]}:{r[0][0]}" for r in top],
            [r[1] / frames for r in top], color="#4c78a8")
    ax.set_xlabel("okupansi kontak")
    ax.set_title(f"Residu kontak ligan (cutoff {cutoff} nm)")
    ax.grid(axis="x", alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig(ctx["a"]("contacts_residue.png"), dpi=150)
    plt.close(fig)

    mean, _ = stat(data)
    return {
        "kontak_atom_rerata": mean,
        "kontak_residu_kunci": " ".join(
            f"{c[2]}{c[1]}:{c[0]}({n / frames:.2f})" for c, n in rows[:8]),
    }


def analyse_pca(ctx):
    """Essential dynamics: covariance of C-alpha motion, projected on PC1/PC2."""
    ca = ctx["groups"].get("C-alpha")
    if ca is None:
        return {}
    eigvec = ctx["a"]("eigenvec.trr")
    ok, _ = gmx_run(["covar", "-s", ctx["tpr"], "-f", ctx["xtc"], "-n", ctx["ndx"],
                     "-o", ctx["a"]("eigenval.xvg"), "-v", eigvec,
                     "-av", ctx["a"]("average.pdb"), "-l", ctx["a"]("covar.log")],
                    stdin=f"{ca}\n{ca}\n", cwd=ctx["cdir"])
    if not ok:
        return {}

    vals, meta = read_xvg(ctx["a"]("eigenval.xvg"))
    out = {}
    if vals.size:
        total = float(np.sum(vals[:, 1]))
        if total > 0:
            out["pca_pc1_persen"] = round(100 * vals[0, 1] / total, 2)
            out["pca_pc2_persen"] = round(100 * vals[1, 1] / total, 2) if len(vals) > 1 else None
            out["pca_pc1_pc2_persen"] = round(
                100 * float(np.sum(vals[:2, 1])) / total, 2)
        export(vals[:20], meta, ctx["adir"], "pca_eigenvalue",
               "Eigenvalue PCA (20 mode pertama)", "mode", "eigenvalue (nm^2)",
               ["mode", "eigenvalue_nm2"])

    proj = ctx["a"]("pca_proj_1_2.xvg")
    ok, _ = gmx_run(["anaeig", "-v", eigvec, "-s", ctx["tpr"], "-f", ctx["xtc"],
                     "-n", ctx["ndx"], "-first", "1", "-last", "2", "-2d", proj],
                    stdin=f"{ca}\n{ca}\n", cwd=ctx["cdir"])
    if ok and os.path.exists(proj):
        data, _ = read_xvg(proj)
        if data.size >= 2:
            fig, ax = plt.subplots(figsize=(5.4, 5.0))
            sc = ax.scatter(data[:, 0], data[:, 1], c=np.arange(len(data)),
                            cmap="viridis", s=6)
            ax.set_xlabel("PC1 (nm)")
            ax.set_ylabel("PC2 (nm)")
            ax.set_title("Proyeksi PCA")
            fig.colorbar(sc, ax=ax, label="urutan frame")
            fig.tight_layout()
            fig.savefig(ctx["a"]("pca_proj.png"), dpi=150)
            plt.close(fig)
            np.savetxt(ctx["a"]("pca_proj_1_2.csv"), data[:, :2], delimiter=",",
                       header="pc1_nm,pc2_nm", comments="", fmt="%.6g")
    return out


def analyse_fel(ctx):
    """Free energy landscape over the PC1/PC2 projection.

    Runs only after PCA, because it is that projection that gets binned. The
    minimum of the surface is reported so that the lowest-energy conformer can
    be pulled out of the trajectory afterwards.
    """
    proj = ctx["a"]("pca_proj_1_2.xvg")
    if not os.path.exists(proj):
        return {}
    xpm = ctx["a"]("fel_gibbs.xpm")
    ok, _ = gmx_run(["sham", "-f", proj, "-ls", xpm, "-notime",
                     "-lsh", ctx["a"]("fel_enthalpy.xpm"),
                     "-lss", ctx["a"]("fel_entropy.xpm")], cwd=ctx["cdir"])
    if not ok or not os.path.exists(xpm):
        return {}

    values, xaxis, yaxis = read_xpm(xpm)
    if values is None:
        return {}
    np.savetxt(ctx["a"]("fel_gibbs.csv"), values, delimiter=",", fmt="%.6g")

    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    extent = None
    if xaxis is not None and yaxis is not None and len(xaxis) and len(yaxis):
        extent = [xaxis.min(), xaxis.max(), yaxis.min(), yaxis.max()]
    im = ax.imshow(values, origin="lower", aspect="auto", extent=extent, cmap="jet")
    ax.set_xlabel("PC1 (nm)")
    ax.set_ylabel("PC2 (nm)")
    ax.set_title("Free energy landscape")
    fig.colorbar(im, ax=ax, label="G (kJ/mol)")
    fig.tight_layout()
    fig.savefig(ctx["a"]("fel_gibbs.png"), dpi=150)
    plt.close(fig)

    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {}
    return {"fel_min_kJmol": float(np.min(finite)),
            "fel_max_kJmol": float(np.max(finite))}


def _resolve_mmpbsa(setting):
    """Accept an env prefix, a python interpreter, or the executable itself."""
    if not setting:
        return shutil.which("gmx_MMPBSA")
    if os.path.isdir(setting):
        for candidate in (os.path.join(setting, "bin", "gmx_MMPBSA"),
                          os.path.join(setting, "gmx_MMPBSA")):
            if os.access(candidate, os.X_OK):
                return candidate
        return None
    if os.access(setting, os.X_OK):
        if os.path.basename(setting).startswith("python"):
            sibling = os.path.join(os.path.dirname(setting), "gmx_MMPBSA")
            return sibling if os.access(sibling, os.X_OK) else None
        return setting
    return None


def _frame_count(ctx):
    """Number of frames in the prepared trajectory, via the RMSD trace if it
    exists (cheap) and gmx check otherwise."""
    rmsd = ctx["a"]("rmsd_protein.xvg")
    if os.path.exists(rmsd):
        data, _ = read_xvg(rmsd)
        if data.size:
            return len(data)
    ok, out = gmx_run(["check", "-f", ctx["xtc"]], cwd=ctx["cdir"])
    m = re.search(r"Step\s+(\d+)", out or "")
    return int(m.group(1)) if m else 0


def analyse_mmpbsa(ctx):
    """Binding free energy with gmx_MMPBSA.

    Kept opt-in behind mmpbsa_python because it needs its own conda
    environment and, unlike everything else here, it costs real time: an
    end-state calculation over a few hundred frames is minutes to hours, not
    seconds.
    """
    exe = ctx["mmpbsa_exe"]
    if not exe:
        return {}
    if "Other" not in ctx["groups"] or "Protein" not in ctx["groups"]:
        return {}
    topol = os.path.join(ctx["cdir"], "topol.top")
    if not os.path.exists(topol):
        say("topol.top not found, MM/PBSA skipped", 6)
        return {}

    total = _frame_count(ctx)
    if total < 2:
        return {}
    wanted = max(1, min(ctx["mmpbsa_frames"], total))
    interval = max(1, total // wanted)

    method = ctx["mmpbsa_method"].lower()
    blocks = [
        "&general",
        f'sys_name="{os.path.basename(ctx["cdir"])}",',
        f"startframe=1, endframe={total}, interval={interval},",
        "verbose=2,",
        "/",
    ]
    if method in ("gb", "both", "gbsa"):
        blocks += ["&gb", f"igb={ctx['mmpbsa_igb']}, saltcon={ctx['mmpbsa_salt']},", "/"]
    if method in ("pb", "both", "pbsa"):
        blocks += ["&pb", f"istrng={ctx['mmpbsa_salt']}, inp=2, radiopt=0,", "/"]
    with open(ctx["a"]("mmpbsa.in"), "w") as fh:
        fh.write("\n".join(blocks) + "\n")

    env = dict(os.environ)
    env["PATH"] = os.path.dirname(os.path.abspath(GMX)) + os.pathsep + env.get("PATH", "")
    cmd = [exe, "-O", "-i", ctx["a"]("mmpbsa.in"),
           "-cs", ctx["tpr"], "-ci", ctx["ndx"],
           "-cg", str(ctx["groups"]["Protein"]), str(ctx["groups"]["Other"]),
           "-ct", ctx["xtc"], "-cp", "topol.top",
           "-o", ctx["a"]("mmpbsa_results.dat"),
           "-eo", ctx["a"]("mmpbsa_frames.csv"), "-nogui"]
    say(f"gmx_MMPBSA: {total} frame, interval {interval}, metode {method}", 6)
    try:
        proc = subprocess.run(cmd, cwd=ctx["cdir"], capture_output=True,
                              text=True, env=env, timeout=86400)
    except (OSError, subprocess.TimeoutExpired) as exc:
        say(f"gmx_MMPBSA failed: {exc}", 6)
        return {}
    if proc.returncode != 0:
        tail = (proc.stdout + proc.stderr).strip().splitlines()[-3:]
        say("gmx_MMPBSA failed: " + " | ".join(tail), 6)
        return {}

    out, dat = {}, ctx["a"]("mmpbsa_results.dat")
    if os.path.exists(dat):
        section = None
        for line in open(dat, errors="replace"):
            if "GENERALIZED BORN" in line:
                section = "gb"
            elif "POISSON BOLTZMANN" in line:
                section = "pb"
            m = re.match(r"^\s*(?:Delta\s+)?TOTAL\s+([-\d.]+)\s+([-\d.]+)", line)
            if m and section:
                out[f"mmpbsa_{section}_dG_kcal"] = float(m.group(1))
                out[f"mmpbsa_{section}_sd"] = float(m.group(2))
    return out


ANALYSIS_FUNCS = {
    "rmsd": analyse_rmsd,
    "rmsf": analyse_rmsf,
    "rg": analyse_rg,
    "sasa": analyse_sasa,
    "hbond": analyse_hbond,
    "contacts": analyse_contacts,
    "pca": analyse_pca,
    "fel": analyse_fel,
    "mmpbsa": analyse_mmpbsa,
}


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def read_config(directory):
    """Same parser LAGMX.py uses, plus the analysis defaults."""
    values = dict(ANALYSIS_DEFAULTS)
    path = os.path.join(directory, "gmx_config.txt")
    if not os.path.exists(path):
        sys.exit(
            f"LAGMX analyze: no 'gmx_config.txt' in {directory}\n\n"
            "analyze_md.py reads its configuration and complex_*/ directories "
            "from the\ndirectory you run it from, exactly like LAGMX.py. Run it "
            "from the same\nplace you ran the simulation:\n\n"
            "    cd run_matrix && python3 ../analyze_md.py\n"
        )
    for line in open(path, errors="replace"):
        if line.strip().startswith("#"):
            continue
        try:
            name, value = line.strip().split(": ", 1)
        except ValueError:
            continue
        values[name.strip()] = value.strip()
    return values


def analyse_complex(cdir, cfg, requested):
    """Run the requested analyses on one complex directory."""
    name = os.path.basename(cdir)
    tpr, xtc = os.path.join(cdir, "md.tpr"), os.path.join(cdir, "md.xtc")
    if not (os.path.exists(tpr) and os.path.exists(xtc)):
        say(f"{name}: no finished production run (md.tpr/md.xtc missing), skipped")
        return None

    adir = os.path.join(cdir, "analysis")
    os.makedirs(adir, exist_ok=True)
    say(f"{name}")

    merged = cfg.get("merged_group", "Protein_LIG")
    groups = build_index(tpr, os.path.join(adir, "analysis.ndx"), cdir, merged)
    if not groups:
        say("make_ndx failed, skipped", 6)
        return None

    skip_ps = float(cfg["analysis_skip_ns"]) * 1000.0
    prepared = prepare_trajectory(cdir, adir, tpr, xtc, groups, merged, skip_ps)
    if not prepared:
        return None

    ctx = {
        "cdir": cdir, "adir": adir, "tpr": tpr, "xtc": prepared,
        "ndx": os.path.join(adir, "analysis.ndx"), "groups": groups,
        "a": lambda f: os.path.join(adir, f),
        "contact_cutoff": float(cfg["analysis_contact_cutoff"]),
        "mmpbsa_exe": _resolve_mmpbsa(cfg.get("mmpbsa_python", "").strip()),
        "mmpbsa_method": cfg.get("mmpbsa_method", "gb"),
        "mmpbsa_frames": int(cfg.get("mmpbsa_frames", "100")),
        "mmpbsa_igb": cfg.get("mmpbsa_igb", "5"),
        "mmpbsa_salt": cfg.get("mmpbsa_salt", "0.150"),
    }

    summary = {"complex": name}
    for key in requested:
        func = ANALYSIS_FUNCS.get(key)
        if func is None:
            continue
        say(f"{key} ...", 6)
        try:
            summary.update(func(ctx) or {})
        except Exception as exc:                                # noqa: BLE001
            say(f"{key} failed: {exc}", 8)

    with open(os.path.join(adir, "summary.csv"), "w") as fh:
        fh.write("besaran,nilai\n")
        for key, value in summary.items():
            fh.write(f"{key},{value}\n")
    return summary


def main():
    global GMX
    directory = os.getcwd()
    cfg = read_config(directory)

    GMX = cfg.get("analysis_gmx", "gmx").strip() or "gmx"
    resolved = shutil.which(GMX) or (GMX if os.access(GMX, os.X_OK) else None)
    if not resolved:
        sys.exit(f"LAGMX analyze: '{GMX}' not found on PATH. Set analysis_gmx in gmx_config.txt.")
    GMX = resolved

    requested = [a.strip().lower() for a in cfg["analysis"].split(",") if a.strip()]
    if "all" in requested:
        requested = list(ALL_ANALYSES)
    unknown = [a for a in requested if a not in ANALYSIS_FUNCS]
    if unknown:
        sys.exit(f"LAGMX analyze: unknown analysis {unknown}; choose from {ALL_ANALYSES}")

    ok, version = gmx_run(["--version"])
    banner = next((l.strip() for l in version.splitlines() if "GROMACS version" in l), "?")
    say(f"gmx      : {GMX}")
    say(f"           {banner}")
    say(f"analisis : {', '.join(requested)}")
    say(f"buang    : {cfg['analysis_skip_ns']} ns pertama")

    complex_dirs = sorted(d for d in glob.glob(os.path.join(directory, "complex*"))
                          if os.path.isdir(d))
    if not complex_dirs:
        sys.exit("LAGMX analyze: no complex*/ directories here.")
    say(f"kompleks : {', '.join(os.path.basename(d) for d in complex_dirs)}\n")

    rows = [r for r in (analyse_complex(d, cfg, requested) for d in complex_dirs) if r]
    if not rows:
        say("\nTidak ada kompleks yang bisa dianalisis.")
        return 1

    columns, seen = [], set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                columns.append(key)
    out_csv = os.path.join(directory, "analysis_summary.csv")
    with open(out_csv, "w") as fh:
        fh.write(",".join(columns) + "\n")
        for row in rows:
            fh.write(",".join(str(row.get(c, "")) for c in columns) + "\n")

    # Printed transposed: metrics down the side, complexes across. With four
    # or five systems a metric-per-row table is the one you can actually read.
    headline = [
        ("rmsd_protein_nm", "RMSD protein (nm)"),
        ("rmsd_ligand_nm", "RMSD ligan (nm)"),
        ("rmsf_mean_nm", "RMSF rerata (nm)"),
        ("rg_nm", "Rg (nm)"),
        ("sasa_terkubur_nm2", "SASA terkubur (nm2)"),
        ("hbond_rerata", "H-bond rerata"),
        ("kontak_atom_rerata", "kontak atom rerata"),
        ("pca_pc1_pc2_persen", "PC1+PC2 (%)"),
        ("fel_min_kJmol", "FEL min (kJ/mol)"),
        ("mmpbsa_gb_dG_kcal", "MM/GBSA dG (kcal/mol)"),
        ("mmpbsa_pb_dG_kcal", "MM/PBSA dG (kcal/mol)"),
    ]
    width = max(14, *(len(r["complex"]) + 2 for r in rows))
    say("\n=================== RINGKASAN ===================")
    say(f"{'besaran':<24}" + "".join(f"{r['complex']:>{width}}" for r in rows))
    say("-" * (24 + width * len(rows)))
    for key, label in headline:
        if key not in columns:
            continue
        cells = []
        for row in rows:
            v = row.get(key, "")
            cells.append(f"{v:>{width}.3f}" if isinstance(v, float) else f"{str(v):>{width}}")
        say(f"{label:<24}" + "".join(cells))

    for row in rows:
        key_res = row.get("kontak_residu_kunci")
        if key_res:
            say(f"\nresidu kunci {row['complex']}: {key_res}")
    say(f"\n-> {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
