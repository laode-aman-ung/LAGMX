#!/usr/bin/env python3
"""Tambal atom MOL2 yang kehilangan hidrogen (valensi kurang dari yang
diharapkan tipe SYBYL-nya), lalu verifikasi lewat antechamber acdoctor.

TERPISAH dari pyAutoGMX.py -- ini alat perbaikan DATA ligand, dijalankan
manual sebelum pyAutoGMX, bukan dipanggil oleh pipeline itu sendiri.
pyAutoGMX.py tidak diubah sama sekali; begitu file .mol2 lig_* di folder
kompleks sudah benar, gentop_gmx() akan memprosesnya seperti biasa.

Kenapa bukan `obabel -h`: dicoba dulu, tapi obabel gagal mengekulisasi
cincin aromatik pada file ini dan malah menambah masalah baru (atom N
lain jadi salah valensi). Perbaikan di sini hanya menambah HIDROGEN YANG
HILANG pada atom yang kekurangan tetangga dibanding tipe SYBYL-nya --
tidak mengubah/menafsirkan ulang bagian molekul yang lain.

Pemakaian:
    python3 fix_ligand_valence.py <input.mol2> [output.mol2]

Kalau output.mol2 tidak diberikan, input ditimpa di tempat (backup asli
otomatis disimpan sebagai <input>.orig kalau belum ada).
"""
import sys
import os
import shutil
import subprocess
import re
import math

# Jumlah tetangga (heavy atom + H) yang diharapkan per tipe SYBYL, dipakai
# untuk mendeteksi atom yang valensinya kurang. Hanya tipe yang benar-benar
# ditemui kasus "Weird atomic valence" oleh antechamber yang perlu akurat;
# tipe lain di sini sekadar jaga-jaga untuk kasus serupa di ligand lain.
EXPECTED_DEGREE = {
    'C.3': 4, 'C.2': 3, 'C.1': 2, 'C.ar': 3, 'C.cat': 3,
    'N.3': 3, 'N.2': 3, 'N.1': 1, 'N.pl3': 3, 'N.4': 4,
    'O.3': 2, 'O.2': 1,
    'S.3': 2, 'S.2': 1,
}

BOND_ORDER_H = 1  # bond H baru selalu single


def parse_mol2(path):
    with open(path) as f:
        lines = f.readlines()
    idx = {}
    for i, line in enumerate(lines):
        if line.startswith('@<TRIPOS>'):
            idx[line.strip()] = i
    mol_i = idx['@<TRIPOS>MOLECULE']
    atom_i = idx['@<TRIPOS>ATOM']
    bond_i = idx['@<TRIPOS>BOND']
    n_atoms, n_bonds = (int(x) for x in lines[mol_i + 2].split()[:2])
    atom_lines = lines[atom_i + 1: atom_i + 1 + n_atoms]
    bond_lines = lines[bond_i + 1: bond_i + 1 + n_bonds]
    return lines, mol_i, atom_i, bond_i, atom_lines, bond_lines


def find_missing_h_atom(atom_lines, bond_lines):
    """Cari SATU atom dengan derajat (jumlah tetangga) kurang dari yang
    diharapkan tipe SYBYL-nya. Kembalikan (atom_id, kekurangan) atau None."""
    degree = {}
    for bl in bond_lines:
        parts = bl.split()
        a, b = int(parts[1]), int(parts[2])
        degree[a] = degree.get(a, 0) + 1
        degree[b] = degree.get(b, 0) + 1

    for al in atom_lines:
        parts = al.split()
        aid = int(parts[0])
        atype = parts[5]
        expected = EXPECTED_DEGREE.get(atype)
        if expected is None:
            continue
        cur = degree.get(aid, 0)
        if cur < expected:
            return aid, expected - cur
    return None


def neighbor_positions(atom_lines, bond_lines, atom_id):
    coords = {}
    for al in atom_lines:
        p = al.split()
        coords[int(p[0])] = (float(p[2]), float(p[3]), float(p[4]))
    neighbors = []
    for bl in bond_lines:
        p = bl.split()
        a, b = int(p[1]), int(p[2])
        if a == atom_id:
            neighbors.append(b)
        elif b == atom_id:
            neighbors.append(a)
    return coords[atom_id], [coords[n] for n in neighbors]


def add_missing_h(path_in, path_out):
    lines, mol_i, atom_i, bond_i, atom_lines, bond_lines = parse_mol2(path_in)
    result = find_missing_h_atom(atom_lines, bond_lines)
    if result is None:
        return False, None
    atom_id, deficit = result
    if deficit != 1:
        raise RuntimeError(
            f"atom {atom_id} kekurangan {deficit} tetangga -- perbaikan "
            "otomatis di script ini cuma menangani kekurangan 1 H, cek manual"
        )

    center, neigh = neighbor_positions(atom_lines, bond_lines, atom_id)
    # Arah H baru: sp3/sp2 completion -- negatif jumlah vektor unit ke
    # tetangga yang sudah ada, dinormalisasi ke panjang ikatan C-H standar.
    vx = vy = vz = 0.0
    for nx, ny, nz in neigh:
        dx, dy, dz = nx - center[0], ny - center[1], nz - center[2]
        norm = math.sqrt(dx * dx + dy * dy + dz * dz)
        vx += dx / norm
        vy += dy / norm
        vz += dz / norm
    dirlen = math.sqrt(vx * vx + vy * vy + vz * vz)
    if dirlen < 1e-6:
        raise RuntimeError(f"geometri tetangga atom {atom_id} degenerate, tidak bisa dihitung arah H")
    bond_len = 1.09  # A, panjang ikatan C-H khas
    hx = center[0] - (vx / dirlen) * bond_len
    hy = center[1] - (vy / dirlen) * bond_len
    hz = center[2] - (vz / dirlen) * bond_len

    n_atoms = len(atom_lines)
    n_bonds = len(bond_lines)
    new_atom_id = n_atoms + 1
    subst_name = atom_lines[atom_id - 1].split()[7] if len(atom_lines[atom_id - 1].split()) > 7 else 'UNK1'
    new_atom_line = (
        f"{new_atom_id:>7} H_fix{new_atom_id:<4}{hx:>11.4f}{hy:>10.4f}{hz:>10.4f} "
        f"H       1 {subst_name}       0.0000\n"
    )
    new_bond_line = f"{n_bonds + 1:>6}{atom_id:>6}{new_atom_id:>6}    1\n"

    new_lines = list(lines)
    new_lines[mol_i + 2] = f" {n_atoms + 1} {n_bonds + 1} 0 0 0\n"
    insert_atom_at = atom_i + 1 + n_atoms
    new_lines.insert(insert_atom_at, new_atom_line)
    # offset bond section start by 1 karena satu baris atom baru disisipkan
    insert_bond_at = bond_i + 1 + 1 + n_bonds
    new_lines.insert(insert_bond_at, new_bond_line)

    with open(path_out, 'w') as f:
        f.writelines(new_lines)
    return True, atom_id


def verify_with_antechamber(mol2_path, workdir):
    # cwd dipaksa ke workdir supaya file sisa antechamber tidak berserakan di
    # tempat lain; karena itu -i/-o HARUS relatif terhadap workdir (basename
    # saja), bukan mol2_path apa adanya (yang bisa masih menyertakan prefix
    # direktori dari path asli, menyebabkan path ganda).
    mol2_basename = os.path.basename(mol2_path)
    base = os.path.splitext(mol2_basename)[0]
    out_basename = f"_verify_{base}.mol2"
    out = os.path.join(workdir, out_basename)
    cmd = [
        "antechamber", "-i", mol2_basename, "-fi", "mol2",
        "-o", out_basename, "-fo", "mol2", "-c", "gas", "-nc", "0", "-s", "2", "-at", "gaff2",
    ]
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    ok = proc.returncode == 0
    if not ok:
        m = re.search(r"Fatal Error.*?\n.*", proc.stdout + proc.stderr)
        detail = m.group(0) if m else (proc.stdout + proc.stderr)[-400:]
    else:
        detail = None
    for f in os.listdir(workdir):
        if f.startswith(("ANTECHAMBER", "ATOMTYPE", "sqm", "NEWPDB", "PREP")) or f.endswith((".AC", ".AC0", ".INF", ".INT")):
            try:
                os.remove(os.path.join(workdir, f))
            except OSError:
                pass
    if os.path.exists(out):
        os.remove(out)
    return ok, detail


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    src = sys.argv[1]
    dst = sys.argv[2] if len(sys.argv) > 2 else src
    workdir = os.path.dirname(os.path.abspath(src)) or "."

    backup = src + ".orig"
    if src == dst and not os.path.exists(backup):
        shutil.copy(src, backup)
        print(f"Backup asli disimpan di {backup}")

    tmp = src + ".fixtmp"
    shutil.copy(src, tmp)

    max_rounds = 5
    fixed_atoms = []
    for _ in range(max_rounds):
        ok, detail = verify_with_antechamber(tmp, workdir)
        if ok:
            break
        changed, atom_id = add_missing_h(tmp, tmp)
        if not changed:
            print(f"GAGAL: antechamber masih menolak tapi tidak ada atom kekurangan-H yang terdeteksi.\nDetail: {detail}")
            os.remove(tmp)
            sys.exit(2)
        fixed_atoms.append(atom_id)
    else:
        print(f"GAGAL setelah {max_rounds} percobaan. Detail terakhir: {detail}")
        os.remove(tmp)
        sys.exit(2)

    shutil.move(tmp, dst)
    print(f"OK: {src} -> {dst}; H ditambahkan pada atom asal ID {fixed_atoms}; verifikasi antechamber lolos.")


if __name__ == "__main__":
    main()
