# pyAutoGMX

Automated preparation and simulation of single- and multi-ligand protein
complexes in GROMACS.

**Status: prototype (TRL 4).** Validated in a controlled matrix of test
scenarios on one workstation. Not yet validated by users outside this
project. See [Scope, evidence and limits](#scope-evidence-and-limits) below.

## What it does

Preparing one protein-ligand system for GROMACS normally takes nine manual
stages, several of them interactive (curating the receptor structure,
parameterizing each ligand, merging topologies, solvating, adding ions,
minimizing, choosing index groups, ...). Every additional ligand multiplies
the ways this can go wrong: identical ligand copies need distinct
moleculetype names or GROMACS silently corrupts the topology, restraint
includes have to live inside the right moleculetype block, atom types have
to be deduplicated by name rather than string-matched, and each ligand's net
charge has to be detected and passed through correctly.

`pyAutoGMX.py` runs the whole thing as one script, in three stages:

1. **Preparation** -- optional structure curation (missing internal loops,
   validated against a reference sequence), ligand parameterization (GAFF2 +
   AM1-BCC or a faster charge method), topology merging, solvation, ion
   addition, energy minimization, index-group construction.
2. **Equilibration** -- short NVT then NPT for every prepared complex,
   before any production run starts.
3. **Production** -- queued one complex at a time. Target length is set
   once in `gmx_config.txt` (`production_ns`); an interrupted or
   deliberately extended run resumes from checkpoint instead of restarting.

A single invocation processes every `complex_*` directory it finds, so it
works the same way for one ligand or several, identical or distinct,
neutral or charged.

## Requirements

- GROMACS (built with CUDA if you want GPU-accelerated `mdrun`)
- AmberTools (`antechamber`, `parmchk2`, `tleap`) for ligand parameterization
- Python 3 with `parmed`
- Optional: a separate Python environment with `pdbfixer`/OpenMM, only
  needed if `fix_structure: yes` is set (structure curation runs through a
  configurable interpreter so pdbfixer does not need to share an
  environment with AmberTools/parmed)

## Repository layout

There is exactly **one** copy of the tool: [`pyAutoGMX.py`](pyAutoGMX.py)
(and its companion [`fix_structure.py`](fix_structure.py)) at the repo
root. Everything else it needs -- `gmx_config.txt`, the `.mdp` templates,
and the `complex_*/` directories it processes -- is read relative to
whatever directory you run it *from*, not relative to where the script
lives. [`run_matrix/`](run_matrix/) is the one maintained, runnable example
in this repo: it has its own `gmx_config.txt` and `.mdp` files (tuned for
the test matrix) plus seven `complex_*/` scenarios. Use it as the template
for a new project: copy `run_matrix/gmx_config.txt` and `run_matrix/*.mdp`
into your own working directory, add your own `complex_*/` folders next to
them, and point a launcher script at `pyAutoGMX.py` the same way
[`run_matrix/run_pyautogmx.sh`](run_matrix/run_pyautogmx.sh) does
(`python3 ../pyAutoGMX.py`, run from inside your working directory).

## Usage

```bash
cd run_matrix
./run_pyautogmx.sh
```

Each `complex_*/` directory must contain:
- one or more receptor files: `rec*.pdb`
- one or more ligand files: `lig*.pdb` or `lig*.mol2`
- optionally `ref*.pdb` (a reference structure with SEQRES records, used
  only if the receptor itself has none and structure curation is enabled)

All simulation options are read from `gmx_config.txt` in the working
directory you run `pyAutoGMX.py` from -- nothing is hardcoded in the
script. See the comments in
[`run_matrix/gmx_config.txt`](run_matrix/gmx_config.txt) for every option,
including:

- `charge` -- `bcc` (AM1-BCC, accurate, slow) or `gas` (Gasteiger, fast, for
  quick iteration)
- `production_ns` -- target production length in nanoseconds; leave empty
  to fall back to the `nsteps` already in `md.mdp`
- `mdrun_options` -- passed straight through to `gmx mdrun` (GPU/thread
  selection)
- `fix_structure` / `fixer_python` / `seqres_reference` -- optional
  structure curation

## Test matrix

[`run_matrix/`](run_matrix/) is a small, version-controlled set of test
scenarios exercising the parts of the pipeline that are easy to get wrong:

| Directory | Exercises |
|---|---|
| `complex_1single` | baseline: one receptor chain, one ligand |
| `complex_2triple` | receptor split across three files/chains, three distinct ligands |
| `complex_3onefile` | one receptor file with three chains, three *identical* ligand copies (moleculetype deduplication) |
| `complex_4mol2` | ligand supplied as MOL2 instead of PDB |
| `complex_5charged` | a charged ligand (ADP3-) alongside a neutral one |
| `complex_6broken` | deliberately malformed ligand -- should fail cleanly, not silently produce a wrong topology |
| `complex_7multi` | five ligands in one system: two identical, three distinct, one charged |

All six valid scenarios have been run end to end (preparation through 5 ns
production) on a single RTX 3070, with zero atom-name mismatches and zero
LINCS warnings across all of them.

`run_matrix/fix_ligand_valence.py` is a standalone repair tool (not called
by `pyAutoGMX.py` itself) for MOL2 ligand files that are missing an explicit
hydrogen -- it detects atoms with fewer neighbors than their SYBYL atom type
requires, completes the missing bond by tetrahedral/trigonal geometry, and
verifies the result through `antechamber`.

## Scope, evidence and limits

**Already demonstrated:** identical ligand copies handled without collision;
up to five mixed ligands (identical, distinct, charged) in one system on the
same code path; net charge detected and verified against the built
topology, including a charged ligand; single-file multi-chain receptors and
MOL2 ligand input running end to end; a wrong reference sequence rejected
(below 90% identity) before it can corrupt a receptor; production length
resumed/extended from checkpoint via a single config value.

**Still being validated:** 100 ns production per system (5 ns reached so
far); RMSD, RMSF and binding-site occupancy analysis.

**Stated limits:**
- 5 ns demonstrates that the pipeline runs correctly, not scientific
  convergence.
- `grompp` still runs with `-maxwarn` during equilibration and production.
- Index-group selection assumes ligands are the only non-protein,
  non-solvent species in the system.
- MOL2 input is not guaranteed complete -- one test ligand needed a missing
  hydrogen restored (see `fix_ligand_valence.py`) before `antechamber`
  would accept it.

The contribution here is orchestration, not new physics: established tools
(GROMACS, AmberTools/antechamber, GAFF2, AM1-BCC) connected into a sequence
that is reproducible for any number of ligands and any number of complexes.

## Citation

If you use pyAutoGMX, please also cite the tools it wraps:

- Abraham, M.J. et al. GROMACS. *SoftwareX* 1-2, 19-25 (2015).
- Wang, J. et al. Development and testing of a general Amber force field.
  *J. Comput. Chem.* 25, 1157-1174 (2004).
- Jakalian, A. et al. Fast, efficient generation of high-quality atomic
  charges. AM1-BCC method.

## License

MIT -- see [LICENSE](LICENSE).

## Authors

La Ode Aman<sup>1</sup>, Arfan<sup>2</sup>, Aiyi Asnawi<sup>3</sup>,
Netty Ino Ischak<sup>4</sup>, Dizky Ramadani Putri Papeo<sup>1</sup>,
Hamsidar Hasan<sup>1</sup>, A. Mu'thi Andy Suryadi<sup>1</sup>

1. Department of Pharmacy, Faculty of Sports and Health, Universitas Negeri Gorontalo, Indonesia
2. Faculty of Pharmacy, Universitas Halu Oleo, Kendari, Southeast Sulawesi, Indonesia
3. Faculty of Pharmacy, Universitas Bhakti Kencana, Bandung, West Java, Indonesia
4. Department of Chemistry, Faculty of Mathematics and Natural Sciences, Universitas Negeri Gorontalo, Indonesia

Corresponding author: La Ode Aman (laode_aman@ung.ac.id)
