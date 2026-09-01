---
title: 'pyAutoGMX: Automated preparation and simulation of single- and multi-ligand protein complexes in GROMACS'
tags:
  - Python
  - molecular dynamics
  - GROMACS
  - drug discovery
  - computational chemistry
  - structural bioinformatics
authors:
  - name: La Ode Aman
    orcid: 0000-0003-4478-6423
    affiliation: 1
    corresponding: true
  - name: Arfan
    orcid: 0000-0003-3004-7101
    affiliation: 2
  - name: Aiyi Asnawi
    orcid: 0000-0002-8179-0520
    affiliation: 3
  - name: Purnawan Pontana Putra
    orcid: 0000-0001-9466-4569
    affiliation: 4
  - name: Netty Ino Ischak
    orcid: 0000-0002-7693-8842
    affiliation: 5
  - name: Dizky Ramadani Putri Papeo
    orcid: 0009-0001-3842-4189
    affiliation: 1
  - name: Hamsidar Hasan
    orcid: 0000-0003-0148-5233
    affiliation: 1
  - name: A. Mu'thi Andy Suryadi
    orcid: 0000-0002-7367-5661
    affiliation: 1
affiliations:
  - name: Department of Pharmacy, Faculty of Sports and Health, Universitas Negeri Gorontalo, Jl. Jenderal Sudirman No. 6, Kota Gorontalo, Gorontalo 96128, Indonesia
    index: 1
  - name: Faculty of Pharmacy, Universitas Halu Oleo, Kendari, Southeast Sulawesi, Indonesia
    index: 2
  - name: Faculty of Pharmacy, Universitas Bhakti Kencana, Bandung, West Java, Indonesia
    index: 3
  - name: Department of Pharmaceutical Chemistry, Faculty of Pharmacy, Universitas Andalas, Padang 25163, Indonesia
    index: 4
  - name: Department of Chemistry, Faculty of Mathematics and Natural Sciences, Universitas Negeri Gorontalo, Gorontalo, Indonesia
    index: 5
date: 27 August 2026
bibliography: paper.bib
---

## Summary

`pyAutoGMX` is a single Python script that automates the preparation,
equilibration, and production stages of a molecular dynamics (MD) simulation
of a protein-ligand complex in GROMACS [@abraham2015gromacs]. Preparing one
such system by hand normally requires nine manual stages -- receptor
structure curation, ligand parameterization, topology merging, solvation,
ion addition, energy minimization, index-group selection, and two stages of
equilibration -- several of which require interactive input. `pyAutoGMX`
runs all of this, plus a configurable production run, from one invocation,
for every complex directory it is pointed at.

The tool is built specifically to handle **multiple ligands per complex
correctly**, whether those ligands are identical copies, distinct
compounds, or charged species. Correctness problems specific to this
setting -- moleculetype name collisions between identical ligand copies,
position-restraint includes attached at the wrong scope in the merged
topology, duplicated `[ atomtypes ]` entries, and per-ligand net charge that
is silently assumed to be neutral -- are handled automatically rather than
left to the user to discover through failed or silently incorrect runs.
Ligand parameterization uses GAFF2 with a configurable charge method
(AM1-BCC via `antechamber` [@wang2006antechamber] for production accuracy,
or Gasteiger charges for fast iteration), and force-field parameters follow
the GAFF methodology [@wang2004gaff; @jakalian2002am1bcc].

Production length is a single configuration value; an interrupted or
deliberately extended run resumes from the GROMACS checkpoint rather than
restarting the whole pipeline, and the tool detects and skips whichever
pipeline stages (preparation, energy minimization, NVT, NPT) a given
complex has already completed.

## Statement of need

Protein-ligand MD is routine in structure-based drug discovery, but the
setup cost per system is high and scales badly with the number of ligands.
A researcher screening a series of candidate compounds, or studying a
receptor with several co-bound cofactors, repeats nine largely manual
stages for every complex, several of which prompt for interactive input.
The multi-ligand failure modes are quiet rather than loud: GROMACS accepts
a topology in which two identical ligand copies share a moleculetype name,
and `antechamber` accepts a charged ligand as neutral. Both produce a
simulation that runs to completion and is wrong.

`pyAutoGMX` is aimed at computational chemists and pharmacy researchers who
need many such systems prepared reproducibly, and who would rather not
maintain a private pipeline of shell glue to get them. It reduces the setup
of an arbitrary number of complexes to editing one configuration file and
issuing one command.

## State of the field

Several existing tools address parts of this workflow, but none combine
automated, config-driven, end-to-end execution with first-class support for
an arbitrary number of ligands per complex:

- **CHARMM-GUI** [@jo2008charmmgui] is a widely used web service that
  generates ready-to-run input files for a protein-ligand system, including
  multi-ligand systems. It does not itself run the MD pipeline: preparation
  happens on CHARMM-GUI's servers, and equilibration/production still have
  to be launched and managed by the user on their own machine.
- **acpype** [@sousadasilva2012acpype] automates ligand topology generation
  (a wrapper around `antechamber`), but does not merge topologies, solvate,
  equilibrate, or run production -- it addresses one stage of the pipeline,
  not the pipeline itself.
- **BioExcel Building Blocks (BioBB)** [@andrio2019biobb] provides a large
  library of composable Python wrappers around individual GROMACS/AmberTools
  commands (`gmx solvate`, `gmx genion`, `gmx grompp`, and others). It gives
  users the building blocks to construct a pipeline, rather than a pipeline
  that runs by default; assembling and validating that pipeline -- including
  the multi-ligand-specific correctness issues above -- is left to the user.

`pyAutoGMX` targets the gap between these: a single, config-driven script
that a user with a set of receptor and ligand files can point at a
directory and get a fully prepared, equilibrated, and simulated system out
the other end, with the multi-ligand-specific failure modes handled by the
tool rather than discovered by the user.

## Software design

`pyAutoGMX` is deliberately a single ~900-line script rather than a package,
depending only on GROMACS, AmberTools, ParmEd, and the Python standard
library. The target user typically works on a shared HPC node or a lab
workstation where installing a package tree is friction rather than
convenience, and a file that can be copied next to the data and run fits
that setting better. The cost is that the tool is harder to unit-test and to
reuse programmatically; this is a prototype trade-off we expect to revisit.

Configuration lives in one plain-text `gmx_config.txt` rather than in
command-line flags. Every parameter that affects the physics -- force field,
water model, box type and margin, charge method, atom type, production
length, `mdrun` resource options -- is recorded in a single file that sits
beside the results, so the settings that produced a trajectory remain
recoverable months later without reconstructing a shell history. Flags would
have made the common case shorter to type and the reproducible case harder.

Stage resumption is derived from output files rather than from a state
database. `detect_stage()` decides how far a complex has progressed by
testing for the *final* artefact of each stage (`em.gro`, `nvt.gro`,
`npt.gro`, `md.cpt`), so a run interrupted midway through equilibration
repeats only that stage and leaves the expensive topology construction and
minimization intact. This keeps the tool stateless and robust to being
killed, at the price of assuming the directory is not edited by hand
between runs.

Two correctness decisions are worth naming. `build_ligand_ids()` assigns
every ligand file a unique moleculetype name up front, appending a numeric
suffix on collision, so identical ligand copies cannot silently share a
`[ moleculetype ]` block. `detect_net_charge()` reads the formal-charge
column of a PDB, or sums the partial charges of a MOL2 file, and passes the
result to `antechamber -nc` instead of letting it default to neutral. Both
run on the input files before any GROMACS command executes, so a
mis-specified system fails during preparation rather than yielding a
trajectory that is wrong but plausible.

## Validation

`pyAutoGMX` ships with a version-controlled test matrix (`run_matrix/`) of
seven scenarios exercising the parts of the pipeline most likely to break:
a single receptor chain with one ligand; a receptor split across three
files/chains with three distinct ligands; a single receptor file with three
chains and three *identical* ligand copies (moleculetype deduplication); a
ligand supplied as MOL2 instead of PDB; a charged ligand (ADP3-)
alongside a neutral one; a deliberately malformed ligand, which is expected
to fail cleanly rather than produce a silently wrong topology; and a
five-ligand system mixing identical, distinct, and charged species.

All six valid scenarios have been run end to end -- preparation through 5 ns
of production -- on a single consumer GPU (RTX 3070), with zero atom-name
mismatches and zero LINCS constraint warnings across all of them. Energy
minimization converged in every case; NPT temperature and density are
reported per scenario in the project README.

## Research impact statement

`pyAutoGMX` has not yet been used outside the authors' group, and we make no
claim of realized external impact. Its case rests on near-term significance
and on being ready for others to pick up.

The validation matrix described above is the primary evidence. All seven
scenarios, including their receptor and ligand inputs and their `.mdp`
parameter files, are committed to the repository, so any reader can
reproduce the reported runs from a clean checkout rather than take the
numbers on trust. The deliberately malformed scenario is part of that
matrix: it documents the boundary at which the tool refuses to proceed, so
users can distinguish a genuine input problem from a tool failure.

For community readiness, the software is MIT-licensed, its test suite runs
in continuous integration across Python 3.10-3.12, and it ships contribution
guidelines. The workflow it automates -- multi-ligand protein-ligand MD in
GROMACS -- is common enough in structure-based drug discovery that we expect
the audience for a working end-to-end pipeline to be broader than our own
group. Realized impact remains to be demonstrated, and we regard growing it
as the main task ahead of this software rather than a settled result.

## Limitations

`pyAutoGMX` is at an early stage (prototype, not yet used outside this
project). Five nanoseconds of production demonstrates that the pipeline
executes correctly, not that any given system's dynamics have converged;
`grompp` is still run with `-maxwarn` during equilibration and production;
index-group selection assumes ligands are the only non-protein,
non-solvent species present; and MOL2 ligand input is not guaranteed to be
chemically complete -- one test ligand required a missing hydrogen to be
restored by a separate repair script before `antechamber` would accept it.
These limitations are stated explicitly in the project README and are the
subject of ongoing work.

## AI usage disclosure

Generative AI tools were used in the preparation of both the software and
this manuscript, and are disclosed here in full.

**Software.** AI coding assistants were used for parts of `pyAutoGMX.py`,
including debugging, drafting individual functions, and clarifying GROMACS
and AmberTools command usage. The pipeline architecture and the majority of
the implementation were written by the authors. All AI-assisted code was
exercised by the authors against the validation matrix described above; the
results reported here were produced by running that code, not by inspecting
it.

**Manuscript.** AI assistance was used for parts of the prose throughout the
manuscript and the project README, for drafting and for editing. In the
revision that produced this version, an AI coding assistant (Claude,
Anthropic) drafted the State of the field, Software design, Research impact
statement, and AI usage disclosure sections, rewrote the Statement of need,
and authored the continuous integration workflow that compiles the
manuscript.

**Verification.** All AI-assisted output was reviewed by the authors before
inclusion. The design decisions reported in Software design were checked
against the source code; the validation results were checked against the
simulation outputs; and each reference was checked against its DOI. The
authors take full responsibility for the content of this paper and of the
software.

## Acknowledgements

This research received no external funding; it was self-funded by the authors.

## References
