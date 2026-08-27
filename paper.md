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
    affiliation: 1
    corresponding: true
  - name: Arfan
    orcid: 0000-0003-3004-7101
    affiliation: 2
  - name: Aiyi Asnawi
    orcid: 0000-0002-8179-0520
    affiliation: 3
  - name: Netty Ino Ischak
    orcid: 0000-0002-7693-8842
    affiliation: 4
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
  - name: Department of Chemistry, Faculty of Mathematics and Natural Sciences, Universitas Negeri Gorontalo, Gorontalo, Indonesia
    index: 4
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

## Validation

`pyAutoGMX` ships with a version-controlled test matrix (`run_matrix/`) of
seven scenarios exercising the parts of the pipeline most likely to break:
a single receptor chain with one ligand; a receptor split across three
files/chains with three distinct ligands; a single receptor file with three
chains and three *identical* ligand copies (moleculetype deduplication); a
ligand supplied as MOL2 instead of PDB; a charged ligand (ADP\textsuperscript{3-})
alongside a neutral one; a deliberately malformed ligand, which is expected
to fail cleanly rather than produce a silently wrong topology; and a
five-ligand system mixing identical, distinct, and charged species.

All six valid scenarios have been run end to end -- preparation through 5 ns
of production -- on a single consumer GPU (RTX 3070), with zero atom-name
mismatches and zero LINCS constraint warnings across all of them. Energy
minimization converged in every case; NPT temperature and density are
reported per scenario in the project README.

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

## Acknowledgements

This research received no external funding; it was self-funded by the authors.

## References
