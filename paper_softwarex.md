# pyAutoGMX: Automated preparation and simulation of single- and multi-ligand protein complexes in GROMACS

**La Ode Aman**^1\*, **Arfan**^2 (ORCID: 0000-0003-3004-7101), **Aiyi Asnawi**^3 (ORCID: 0000-0002-8179-0520), **Netty Ino Ischak**^4 (ORCID: 0000-0002-7693-8842), **Dizky Ramadani Putri Papeo**^1 (ORCID: 0009-0001-3842-4189), **Hamsidar Hasan**^1 (ORCID: 0000-0003-0148-5233), **A. Mu'thi Andy Suryadi**^1 (ORCID: 0000-0002-7367-5661)

^1 Department of Pharmacy, Faculty of Sports and Health, Universitas Negeri Gorontalo, Jl. Jenderal Sudirman No. 6, Kota Gorontalo, Gorontalo 96128, Indonesia
^2 Faculty of Pharmacy, Universitas Halu Oleo, Kendari, Southeast Sulawesi, Indonesia
^3 Faculty of Pharmacy, Universitas Bhakti Kencana, Bandung, West Java, Indonesia
^4 Department of Chemistry, Faculty of Mathematics and Natural Sciences, Universitas Negeri Gorontalo, Gorontalo, Indonesia

\* Corresponding author: La Ode Aman (laode_aman@ung.ac.id)

Author emails: laode_aman@ung.ac.id (La Ode Aman); arfan09@uho.ac.id (Arfan);
aiyi.asnawi@bku.ac.id (Aiyi Asnawi); nettyischak@ung.ac.id (Netty Ino Ischak);
dizky@ung.ac.id (Dizky Ramadani Putri Papeo); hamsidar.hasan@ung.ac.id
(Hamsidar Hasan); a.muthi@ung.ac.id (A. Mu'thi Andy Suryadi)

<!--
NOTE FOR SUBMISSION: SoftwareX requires this manuscript in their Elsevier
Word/LaTeX template, not plain Markdown. This file is the content source --
copy each section into the template as-is. Fields marked TODO below need
your input before submission.
-->

## Code metadata

| Nr. | Code metadata description | Please fill in this column |
|---|---|---|
| C1 | Current code version | v1.0.0 <!-- TODO: create this tag, see note below --> |
| C2 | Permanent link to code/repository used for this code version | <https://github.com/laode-aman-ung/pyAutoGMX> |
| C3 | Permanent link to reproducible capsule | TODO (optional; e.g. Code Ocean, if prepared) |
| C4 | Legal code license | MIT |
| C5 | Code versioning system used | git |
| C6 | Software code languages, tools and services used | Python 3; GROMACS; AmberTools (`antechamber`, `parmchk2`, `tleap`); `parmed`; optionally `pdbfixer`/OpenMM |
| C7 | Compilation requirements, operating environments and dependencies | Linux (developed and tested on Ubuntu); GROMACS (CUDA build required for GPU-accelerated `mdrun`, optional otherwise); AmberTools; Python 3 with `parmed` |
| C8 | If available, link to developer documentation/manual | `README.md` and `CONTRIBUTING.md` in the repository (C2) |
| C9 | Support email for questions | laode_aman@ung.ac.id |

<!--
TODO before submission: tag the commit this paper describes as a release,
e.g.:
    git tag -a v1.0.0 -m "Initial public release" cf7e93571f65795b6fdcd5cd5863c7ae930d26c5
    git push origin v1.0.0
so C1/C2 point at something immutable rather than a moving branch.
-->

## Abstract

`pyAutoGMX` is a single-script pipeline that automates the preparation,
equilibration, and production stages of a GROMACS molecular dynamics (MD)
simulation for one or more protein-ligand complexes at once. It is built
specifically to handle multiple ligands per complex correctly -- identical
copies, distinct compounds, and charged species -- a case that introduces
correctness problems (moleculetype name collisions, misplaced restraint
scoping, duplicated atom types, incorrectly assumed-neutral charge) that
existing preparation tools either don't address or leave to the user to
discover. A configuration-driven production stage resumes from checkpoint
automatically, whether a run was interrupted or deliberately extended. The
tool is validated against a version-controlled matrix of seven scenarios,
six of which run end to end with zero atom-name mismatches and zero LINCS
constraint warnings.

## 1. Motivation and significance

Preparing a single protein-ligand system for a GROMACS MD simulation
normally requires roughly nine manual stages: curating the receptor
structure, parameterizing the ligand, merging topologies, solvating the
system, adding ions, running energy minimization, constructing index
groups, and two stages of equilibration before production can even start.
Several of these stages require interactive input -- selecting a group in
`gmx genion`, choosing a histidine protonation state in `pdb2gmx` -- which
makes the process slow to repeat and easy to get subtly wrong.

Adding more than one ligand to the system does not make this proportionally
harder; it makes it harder in ways that are easy to miss until a run
silently produces incorrect results. Four problems recur specifically in
multi-ligand systems:

1. **Moleculetype collisions.** If two ligand copies in the same complex
   share a residue name (a common case: the same compound bound at two
   sites, or multiple copies of a homo-oligomer's ligand), GROMACS needs
   distinct `moleculetype` names, restraint files, and `[ molecules ]`
   entries for each copy -- reusing the same name silently corrupts the
   topology rather than raising an error.
2. **Restraint scope.** Each ligand's position-restraint `#include` has to
   live inside that specific ligand's `moleculetype` block. Placed at the
   top of `topol.top` instead -- a natural place to put it by analogy with
   a single-ligand system -- it ends up applying only to whichever molecule
   happens to be defined last.
3. **Duplicated atom types.** Merging `[ atomtypes ]` sections from
   multiple ligand topologies naively (e.g. by string-matching lines)
   produces duplicate or conflicting entries once ligand count grows past
   one or two.
4. **Net charge.** Antechamber assumes a neutral molecule unless told
   otherwise, and gives no error message if that assumption is wrong -- so
   a charged ligand (e.g. a nucleotide such as ADP, charge -3) silently
   gets the wrong partial charges unless the net charge is detected from
   the input and passed through explicitly, per ligand.

Existing tools each address part of the overall problem but not the
combination of "runs end to end" and "handles the multi-ligand case
correctly by default":

- **CHARMM-GUI** [@jo2008charmmgui] is a widely used web service that
  generates a ready-to-run input file set for a protein-ligand system,
  including multi-ligand systems. It does not run the simulation itself:
  preparation happens on CHARMM-GUI's servers, and the user still launches
  and manages equilibration and production on their own machine.
- **acpype** [@sousadasilva2012acpype] automates ligand topology generation
  (a wrapper around `antechamber` [@wang2006antechamber]) but stops there --
  it does not merge topologies, solvate, equilibrate, or run production.
- **BioExcel Building Blocks (BioBB)** [@andrio2019biobb] provides a large
  library of composable Python wrappers around individual GROMACS/AmberTools
  commands. It gives a user the parts to build a pipeline with, not a
  pipeline that runs correctly by default; assembling those parts into a
  working, multi-ligand-safe sequence -- including the four problems above --
  is left to the user.

`pyAutoGMX` targets that gap: one script, driven by a single configuration
file, that a user with a directory of receptor and ligand files can invoke
once and get a fully prepared, equilibrated, and simulated system for
however many ligands (and however many complexes) they have -- with the
multi-ligand-specific failure modes handled by the tool rather than
discovered by the user after the fact.

## 2. Software description

### 2.1. Software architecture

`pyAutoGMX.py` processes every `complex_*/` directory it finds in its
working directory, in three stages, run for *all* complexes before moving
to the next stage (so equilibration for every system completes before the
first production run starts, and a single system's late failure doesn't
waste GPU time better spent on the rest of the batch):

1. **Preparation.** Optional structure curation via `pdbfixer` (run through
   a separately configurable Python interpreter, so `pdbfixer`/OpenMM never
   needs to share an environment with AmberTools/`parmed`), including a
   sequence-identity check: a reference SEQRES sequence is only attached
   (and missing residues modelled) if it matches at least 90% of the
   observed residues, preventing a mismatched reference from silently
   introducing the wrong residues. Ligand parameterization runs through
   `antechamber` with a configurable charge method (AM1-BCC
   [@jakalian2002am1bcc] for production accuracy, or Gasteiger charges for
   fast iteration) and GAFF2 atom types [@wang2004gaff], followed by
   `parmchk2` and `tleap`. Topologies are merged with moleculetype names
   deduplicated automatically, restraint includes placed inside the correct
   moleculetype block, and atom types deduplicated by name. The system is
   then solvated, ionized, and energy-minimized, and index groups are
   constructed generically (via GROMACS's own "Other" group, which
   captures every ligand present without the tool needing to know their
   names or count).
2. **Equilibration.** Short NVT then NPT runs, executed for every prepared
   complex before production begins for any of them.
3. **Production.** Queued one complex at a time. Target length is a single
   configuration value (`production_ns`). If a complex already has a prior
   production checkpoint, `pyAutoGMX` detects it, skips preparation and
   equilibration entirely, and extends the existing run to the new target
   length via `gmx convert-tpr -until` (an *absolute* end time, deliberately
   used instead of `mdrun -nsteps`/`convert-tpr -nsteps`, both of which are
   relative to whatever step the checkpoint is already at and would silently
   overshoot the intended target). More generally, the tool detects the
   furthest pipeline stage each complex directory has actually completed
   (none / energy-minimized / NVT-complete / NPT-complete / production) and
   resumes from there, so an interrupted run only repeats the specific
   stage it was interrupted in, not the whole pipeline.

All simulation options -- box geometry, solvent model, charge method, atom
type, force field, water model, `mdrun` resource flags, production length,
and structure-curation settings -- are read from a `gmx_config.txt` file in
the working directory; nothing is hardcoded in the script.

### 2.2. Software functionalities

- End-to-end preparation, equilibration, and production for an arbitrary
  number of complexes and an arbitrary number of ligands per complex, from
  one invocation.
- Automatic handling of identical ligand copies, distinct ligands, and
  charged ligands within the same system.
- Structure curation with reference-sequence validation (rejects a
  mismatched reference rather than modelling the wrong residues).
- Configurable ligand charge method (AM1-BCC or Gasteiger) as a
  speed/accuracy trade-off.
- Configuration-driven, checkpoint-resumable production length: extending
  a completed run, or resuming an interrupted one, requires changing one
  number, not re-deriving `nsteps` or re-running preparation.
- A standalone repair utility (`fix_ligand_valence.py`, deliberately kept
  separate from the main pipeline) for MOL2 ligand input missing an
  explicit hydrogen: it detects atoms with fewer neighbors than their
  SYBYL atom type requires, completes the missing bond by
  tetrahedral/trigonal geometry, and verifies the result by re-running it
  through `antechamber`.
- A version-controlled test matrix (Section 3) that doubles as a
  regression suite for the multi-ligand-specific behavior described in
  Section 1.

### 2.3. Sample usage

```bash
cd my_project/          # contains gmx_config.txt, *.mdp, and complex_*/ dirs
python3 /path/to/pyAutoGMX.py
```

A minimal `gmx_config.txt` excerpt controlling ligand charge method and
production length:

```
charge: bcc            # AM1-BCC (accurate); use "gas" for fast iteration
production_ns: 5        # resumes/extends automatically on rerun
mdrun_options: -ntmpi 1 -ntomp 16 -nb gpu
```

Each `complex_*/` directory needs only its input structures:

```
complex_mysystem/
├── rec_receptor.pdb
├── lig_ligandA.pdb
└── lig_ligandB.mol2
```

## 3. Illustrative example

The repository's `run_matrix/` directory is both the tool's documentation
by example and its test suite: seven scenarios, each isolating one
behavior that is easy to get wrong (Table 1).

**Table 1.** Test matrix scenarios.

| Scenario | Exercises |
|---|---|
| `complex_1single` | Baseline: one receptor chain, one ligand |
| `complex_2triple` | Receptor split across three files/chains, three distinct ligands |
| `complex_3onefile` | One receptor file with three chains, three *identical* ligand copies (moleculetype deduplication) |
| `complex_4mol2` | Ligand supplied as MOL2 instead of PDB |
| `complex_5charged` | A charged ligand (ADP\textsuperscript{3-}) alongside a neutral one |
| `complex_6broken` | Deliberately malformed ligand -- expected to fail cleanly, not silently produce a wrong topology |
| `complex_7multi` | Five ligands in one system: two identical, three distinct, one charged |

The most demanding of these, `complex_7multi`, combines every failure mode
described in Section 1 in a single system: two copies of the same ligand
(triggering automatic moleculetype renaming), three additional distinct
ligands, and one charged species, all merged into one topology and index
group. It is invoked exactly the same way as `complex_1single` -- the tool
does not need to be told how many ligands a complex has, or which ones are
identical.

**Table 2.** Verification results, all six valid scenarios, 5 ns
production on a single consumer GPU (RTX 3070).

| Metric | 1single | 2triple | 3-chain, 1 file | MOL2 input | Charged | 5-ligand mixed |
|---|---|---|---|---|---|---|
| Ligands | 1 | 3 | 3 (identical) | 3 (identical) | 2 | 5 (2 identical + 3 distinct, 1 charged) |
| Total atoms | 56,201 | 75,819 | 78,555 | 78,499 | 64,453 | 152,944 |
| Temperature, NPT (K) | 299.92 ± 1.38 | 300.03 ± 1.21 | 299.94 ± 1.25 | 299.98 ± 1.15 | 299.91 ± 1.24 | 300.05 ± 0.78 |
| Density, NPT (kg/m^3) | 1002.7 ± 1.5 | 1018.0 ± 0.6 | 1020.2 ± 1.0 | 1016.1 ± 0.7 | 1000.1 ± 1.3 | 986.9 ± 0.6 |
| Production performance | 212 ns/day | 160 ns/day | 156 ns/day | 154 ns/day | 187 ns/day | 83 ns/day |
| Atom-name mismatches | 0 | 0 | 0 | 0 | 0 | 0 |
| LINCS warnings | 0 | 0 | 0 | 0 | 0 | 0 |

The seventh scenario, `complex_6broken`, is a deliberately malformed
ligand and is *expected* to fail; it exists to confirm that a bad input is
rejected during preparation rather than silently propagating into a wrong
topology. `complex_4mol2`'s ligand input additionally needed one missing
hydrogen restored (via `fix_ligand_valence.py`, Section 2.2) before
`antechamber` would accept it -- documented in the repository as a stated
limitation of MOL2 input handling, not silently worked around.

## 4. Impact

`pyAutoGMX` is aimed at researchers running structure-based MD who need to
prepare more than one protein-ligand system -- or more than one ligand per
system -- without re-deriving the multi-ligand-specific fixes in Section 1
by hand each time. Concretely, it turns what would otherwise be nine manual
stages per system, several requiring interactive choices, into one
config-driven invocation that processes an arbitrary batch of complexes
unattended. This is most useful in the early stages of structure-based drug
design or hit-to-lead work, where the same receptor is screened against
many candidate ligands, or where a single binding site must be evaluated
with several ligands present simultaneously (e.g. cofactor plus candidate,
or multiple binding sites on a multimeric receptor).

The tool's contribution is orchestration -- established, individually
validated components (GROMACS [@abraham2015gromacs], AmberTools,
GAFF2/AM1-BCC) connected into a sequence that is reproducible for any
number of ligands and any number of complexes -- rather than a new
simulation method. Its value is in removing a specific, recurring source of
manual error, evidenced by the fact that the multi-ligand correctness
problems in Section 1 were each discovered as real failures during this
project's own development, not anticipated in the abstract.

`pyAutoGMX` is at an early stage: validated on a controlled matrix of seven
scenarios on one workstation, not yet used or validated by anyone outside
this project, and 5 ns of production demonstrates that the pipeline
executes correctly rather than that any given system's dynamics have
converged. These limitations, and others (`grompp` is still run with
`-maxwarn`; index-group selection assumes ligands are the only
non-protein, non-solvent species present), are stated explicitly in the
repository and are the subject of ongoing work.

## 5. Conclusions

`pyAutoGMX` automates the GROMACS preparation-through-production pipeline
for protein-ligand complexes, with correctness for multi-ligand systems --
identical copies, distinct compounds, charged species -- as a first-class
design goal rather than an afterthought. It is validated against a
version-controlled, seven-scenario test matrix, six of which run end to
end with zero topology or constraint-stability warnings across systems
ranging from one ligand to five. The software, its test matrix, and its
stated limitations are publicly available under the MIT license at
<https://github.com/laode-aman-ung/pyAutoGMX>.

## Acknowledgements

This work was presented in preliminary form at the Gorontalo International
Multidisciplinary Health and Wellness Conference (GIMHWELT).

## Declaration of competing interest

The authors declare no competing financial interests or personal
relationships that could have appeared to influence the work reported in
this paper.

## References
