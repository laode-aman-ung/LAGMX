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
- Python 3.10+ with `parmed`
- Optional: a separate Python environment with `pdbfixer`/OpenMM, only
  needed if `fix_structure: yes` is set (structure curation runs through a
  configurable interpreter so pdbfixer does not need to share an
  environment with AmberTools/parmed)

## Installation

All of the above are available from conda-forge and bioconda:

```bash
conda env create -f environment.yml
conda activate pyautogmx
```

Use `mamba` if you have it; solving this environment with plain `conda` can
take several minutes. See [`environment.yml`](environment.yml) for how to
create the optional `pdbfixer`/OpenMM environment.

## Quickstart

To check that the pipeline works on your machine before committing to a real
run:

```bash
./quickstart.sh
```

This runs the smallest scenario from the test matrix -- one receptor chain,
one ligand -- through the full pipeline in a few minutes, using Gasteiger
charges and 10 ps of equilibration and production. It is a smoke test, not a
scientific result: it answers "does this work here?", nothing more. Output
lands in `quickstart_run/`, which is gitignored and safe to delete. If it
fails, the log it points you at is what to attach to an issue.

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

This uses whatever `gmx`, `antechamber` and `python3` are on your PATH and
lets GROMACS choose the hardware it finds, so it runs on a CPU-only machine
as well as a GPU one. To reproduce the stricter setup used for the runs
reported in the paper -- refusing to start unless the `gmx` binary is a CUDA
build and an NVIDIA GPU is visible, so that `mdrun` fails loudly instead of
falling back to the CPU -- pass `--require-gpu`:

```bash
./run_pyautogmx.sh --require-gpu
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

## Preparing your own system

`quickstart.sh` proves the pipeline runs on your machine. To simulate your
own complexes, create a working directory, copy the configuration and `.mdp`
templates into it, and put one directory per complex next to them:

```bash
mkdir -p ~/my-project && cd ~/my-project
cp /path/to/pyAutoGMX/run_matrix/gmx_config.txt .
cp /path/to/pyAutoGMX/run_matrix/*.mdp .
mkdir complex_egfr
# ... put your receptor and ligand files in complex_egfr/ ...
python3 /path/to/pyAutoGMX/pyAutoGMX.py
```

Every `complex_*/` directory in the working directory is processed in turn,
so several systems are simply several directories.

### What each complex directory must contain

File names are matched by prefix, in sorted order:

| Pattern | Meaning |
| --- | --- |
| `rec*.pdb` | Receptor. One file per chain, or one file containing several chains. |
| `lig*.pdb` or `lig*.mol2` | Ligands. One file per ligand copy, including identical copies. |
| `ref*.pdb` | Optional. A reference structure carrying SEQRES records, used only when the receptor has none and `fix_structure: yes` is set. |

Anything that does not match `rec*` or `lig*` is ignored, so a ligand named
`ATP.pdb` will not be picked up — rename it `lig_ATP.pdb`.

### What the inputs have to satisfy

These are requirements of the pipeline, not style preferences:

- **Each ligand copy is its own file.** Three copies of the same compound
  means three `lig*` files. pyAutoGMX gives each one a distinct
  moleculetype name; it does not split a single file containing three
  copies.
- **Ligands must be chemically complete**, hydrogens included. `antechamber`
  rejects a ligand whose valences do not add up, and this is the most common
  reason a run stops during preparation. `run_matrix/fix_ligand_valence.py`
  repairs the specific case of MOL2 atoms missing hydrogens relative to
  their SYBYL atom type:

  ```bash
  python3 /path/to/pyAutoGMX/run_matrix/fix_ligand_valence.py lig_mine.mol2
  ```

  Note that `obabel -h` is **not** a safe substitute here: on our test
  ligands it failed to equalise aromatic rings and introduced new valence
  errors elsewhere in the molecule.
- **Formal charge must be readable from the file.** With `net_charge: auto`
  (the default), the charge of a PDB ligand is read from the formal-charge
  column (columns 79-80, e.g. `N1+`, `O1-`), and the charge of a MOL2 ligand
  is the sum of its partial charges. If that column is blank on a charged
  ligand, `antechamber` will treat the molecule as neutral and produce wrong
  charges without reporting an error. Set `net_charge` explicitly if you are
  unsure.
- **The receptor should contain protein only.** Index-group selection
  assumes ligands are the only non-protein, non-solvent species present, so
  crystallographic waters, ions and cofactors you do not intend to simulate
  should be removed from `rec*.pdb` first.

### Then check the settings

Open `gmx_config.txt` and set at least `charge` (`gas` for a fast first
attempt, `bcc` for production accuracy), `production_ns`, and
`mdrun_options` if you want to pin threads or a GPU. Every option is
documented in the file itself.

## Running on a GPU

`mdrun_options` is empty by default, so `gmx mdrun` runs with GROMACS' own
`-nb auto`: it uses a GPU when it detects a compatible one, and the CPU
otherwise. On a machine with a working CUDA build and an NVIDIA card,
nothing needs to be configured — the GPU is used automatically.

### "GPU support" in the binary does not mean the GPU is used

A GROMACS binary can report GPU support and still run entirely on the CPU if
the runtime driver is not visible to it. This is what that looks like in
`md.log`:

```
GPU support         : OpenCL
Running on 1 node with total 20 cores, 40 processing units
    (GPU detection failed: No valid OpenCL driver found)
```

The run completes normally and reports no error — it is simply slow. A
100 ns production run can finish on the CPU without anyone noticing.

### Check that the GPU was actually used

After a run:

```bash
grep -E "GPU detection|GPU selected|using .* GPU" complex_*/md.log
```

`1 GPU selected` means the GPU was used. `GPU detection failed` means it was
not, whatever the binary claims to support.

### Make a missing GPU an error instead

If you expect a GPU and would rather the run stop than quietly fall back,
request it explicitly in `gmx_config.txt`:

```
mdrun_options: -ntmpi 1 -ntomp 16 -nb gpu
```

With `-nb gpu`, `mdrun` fails instead of using the CPU. `-ntmpi 1` is
required whenever `-ntomp` is set. The launcher has a matching preflight
that refuses to start unless the binary is a CUDA build and an NVIDIA GPU is
visible:

```bash
cd run_matrix
./run_pyautogmx.sh --require-gpu
```

### Getting a CUDA build

The `gromacs` package installed by `environment.yml` is the generic build,
which reports OpenCL rather than CUDA. For CUDA, ask for that build
explicitly:

```bash
conda install -c conda-forge "gromacs=*=nompi_cuda*"
```

See the comments in [`environment.yml`](environment.yml) for why the CUDA
build is not selected by default even on an NVIDIA machine.

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
