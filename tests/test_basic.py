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
    "run_matrix/fix_ligand_valence.py",
]

# complex_6broken is intentionally malformed and expected to fail the real
# pipeline; it's excluded from checks that assume a working scenario.
KNOWN_BROKEN_SCENARIOS = {"complex_6broken"}


def test_python_files_exist_and_compile():
    for rel in PY_FILES:
        path = os.path.join(ROOT, rel)
        assert os.path.exists(path), f"missing {rel}"
        subprocess.run([sys.executable, "-m", "py_compile", path], check=True)


def test_only_one_copy_of_pyautogmx_exists():
    # pyAutoGMX.py used to be duplicated into run_matrix/ and the two
    # copies drifted apart during development (one had checkpoint-resume,
    # the other had SEQRES curation, neither had both). There must be
    # exactly one canonical copy from here on -- run_matrix/run_pyautogmx.sh
    # invokes it via a relative path (../pyAutoGMX.py) instead of a second
    # copy living next to the test scenarios.
    assert os.path.exists(os.path.join(ROOT, "pyAutoGMX.py"))
    assert not os.path.exists(os.path.join(ROOT, "run_matrix", "pyAutoGMX.py")), (
        "run_matrix/pyAutoGMX.py should not exist -- there must be exactly "
        "one copy of the script (repo root); see run_pyautogmx.sh, which "
        "invokes ../pyAutoGMX.py"
    )
    assert not os.path.exists(os.path.join(ROOT, "run_matrix", "fix_structure.py")), (
        "run_matrix/fix_structure.py should not exist -- fix_structure.py "
        "is resolved relative to pyAutoGMX.py's own location, not cwd, so "
        "it only needs to exist once, at the repo root"
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
    path = os.path.join(ROOT, "run_matrix", "gmx_config.txt")
    with open(path) as f:
        text = f.read()
    for key in required:
        assert f"{key}:" in text, f"run_matrix/gmx_config.txt: missing required key '{key}'"


def test_mdp_templates_have_required_params():
    for name in ["em.mdp", "nvt.mdp", "npt.mdp", "md.mdp", "ions.mdp"]:
        path = os.path.join(ROOT, "run_matrix", name)
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
