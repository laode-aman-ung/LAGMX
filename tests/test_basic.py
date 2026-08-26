"""Lightweight tests that don't need GROMACS/AmberTools installed.

These check repository structure, config/mdp sanity, and the pure-Python
logic in fix_ligand_valence.py. They do NOT run the actual MD pipeline --
that needs GROMACS/AmberTools and is verified manually against run_matrix/
(see CONTRIBUTING.md).
"""
import glob
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PY_FILES = [
    "pyAutoGMX.py",
    "fix_structure.py",
    "run_matrix/pyAutoGMX.py",
    "run_matrix/fix_ligand_valence.py",
    "run_matrix/fix_structure.py",
]

# complex_6broken is intentionally malformed and expected to fail the real
# pipeline; it's excluded from checks that assume a working scenario.
KNOWN_BROKEN_SCENARIOS = {"complex_6broken"}


def test_python_files_exist_and_compile():
    for rel in PY_FILES:
        path = os.path.join(ROOT, rel)
        assert os.path.exists(path), f"missing {rel}"
        subprocess.run([sys.executable, "-m", "py_compile", path], check=True)


def test_top_level_and_run_matrix_scripts_match():
    # pyAutoGMX.py exists in two places (repo root and run_matrix/) and is
    # meant to be kept identical; this test exists specifically because the
    # two copies drifted apart once during development (see README).
    with open(os.path.join(ROOT, "pyAutoGMX.py")) as f:
        root_version = f.read()
    with open(os.path.join(ROOT, "run_matrix", "pyAutoGMX.py")) as f:
        matrix_version = f.read()
    assert root_version == matrix_version, (
        "pyAutoGMX.py and run_matrix/pyAutoGMX.py have diverged -- "
        "keep them identical or this project ends up publishing "
        "inconsistent claims about what the tool can do"
    )


def test_run_matrix_scenarios_have_inputs():
    matrix = os.path.join(ROOT, "run_matrix")
    complex_dirs = sorted(
        d for d in os.listdir(matrix)
        if os.path.isdir(os.path.join(matrix, d)) and d.startswith("complex_")
    )
    assert len(complex_dirs) >= 6, "expected at least 6 scenarios in run_matrix/"
    for d in complex_dirs:
        p = os.path.join(matrix, d)
        rec = glob.glob(os.path.join(p, "rec*.pdb"))
        lig = glob.glob(os.path.join(p, "lig*.pdb")) + glob.glob(os.path.join(p, "lig*.mol2"))
        assert rec, f"{d}: no receptor file (rec*.pdb)"
        assert lig, f"{d}: no ligand file (lig*.pdb / lig*.mol2)"


def test_gmx_config_has_required_keys():
    required = ["box_type", "distance", "solvent", "charge", "atom_type"]
    for cfg in ["gmx_config.txt", "run_matrix/gmx_config.txt"]:
        path = os.path.join(ROOT, cfg)
        with open(path) as f:
            text = f.read()
        for key in required:
            assert f"{key}:" in text, f"{cfg}: missing required key '{key}'"


def test_mdp_templates_have_required_params():
    for mdp_dir in ["", "run_matrix"]:
        for name in ["em.mdp", "nvt.mdp", "npt.mdp", "md.mdp", "ions.mdp"]:
            path = os.path.join(ROOT, mdp_dir, name)
            assert os.path.exists(path), f"missing {path}"
            with open(path) as f:
                text = f.read()
            assert "integrator" in text, f"{path}: missing 'integrator ='"


def test_fix_ligand_valence_detects_missing_hydrogen(tmp_path):
    sys.path.insert(0, os.path.join(ROOT, "run_matrix"))
    import fix_ligand_valence as flv

    # Minimal synthetic MOL2: atom 1 is C.3 (expects 4 neighbours) but only
    # has 3 -- the same class of defect fix_ligand_valence.py was built to
    # repair (see README: complex_4mol2 needed exactly this fix).
    mol2 = """@<TRIPOS>MOLECULE
test
4 3 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1  C1        0.000    0.000    0.000 C.3     1  UNL1    0.0000
      2  C2        1.500    0.000    0.000 C.3     1  UNL1    0.0000
      3  N1       -0.750    1.299    0.000 N.3     1  UNL1    0.0000
      4  C3       -0.750   -1.299    0.000 C.3     1  UNL1    0.0000
@<TRIPOS>BOND
     1     1     2    1
     2     1     3    1
     3     1     4    1
"""
    p = tmp_path / "synthetic.mol2"
    p.write_text(mol2)

    _, _, _, _, atom_lines, bond_lines = flv.parse_mol2(str(p))
    result = flv.find_missing_h_atom(atom_lines, bond_lines)

    assert result is not None, "expected atom 1 (C.3, only 3 neighbours) to be flagged"
    atom_id, deficit = result
    assert atom_id == 1
    assert deficit == 1
