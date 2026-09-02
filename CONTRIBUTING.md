# Contributing to LAGMX

Thanks for considering a contribution. This project is a research tool at
an early (prototype) stage, so the process below is intentionally light.

## Reporting a bug

Open a GitHub issue and include:

- your GROMACS and AmberTools versions (`gmx --version`, `antechamber -v`)
- the `gmx_config.txt` you used
- the failing command and the relevant excerpt of its output (the actual
  error, not just "it failed")
- if the problem is specific to one complex, the receptor/ligand file(s)
  involved (or a minimal reproduction) -- most bugs in this project turn
  out to depend on the exact structure of the input, not just the config

## Proposing a change

1. Fork the repository and branch from `main`.
2. Keep the change scoped to one thing. Don't bundle an unrelated
   reformat, rename, or "while I was in there" cleanup with a functional
   change -- it makes the diff hard to review and hard to revert if wrong.
3. Run the test suite before opening a pull request:

   ```bash
   pip install pytest
   pytest tests/ -v
   ```

   CI runs the same tests on every push/PR. They check repository
   structure and Python syntax/logic that doesn't require GROMACS or
   AmberTools to be installed -- they do **not** run the actual MD
   pipeline, since GitHub-hosted runners don't have GROMACS/AmberTools
   available.
4. If your change touches pipeline behavior (anything in `LAGMX.py`,
   the `.mdp` templates, or `gmx_config.txt` parsing), also run at least
   one scenario from `run_matrix/` end to end locally, and say what you ran
   and what happened in the pull request description. This is currently
   the only way changes to the actual simulation behavior get verified --
   please don't skip it for anything beyond a trivial/cosmetic change.
5. If you're fixing something that a new `run_matrix/complex_*` scenario
   would catch in the future, consider adding one (see below) in the same
   PR.

## Adding a test scenario to `run_matrix/`

`LAGMX.py` discovers any `complex_*` directory automatically -- no
registration step needed. To add one:

1. Create `run_matrix/complex_<name>/`.
2. Add one or more receptor files: `rec*.pdb`.
3. Add one or more ligand files: `lig*.pdb` or `lig*.mol2`.
4. That's it. `tests/test_basic.py` will pick it up automatically (it
   asserts every `complex_*` directory has at least one receptor and one
   ligand file); running it through the actual pipeline still requires
   GROMACS/AmberTools locally, as described above.

Prefer a scenario that isolates one specific behavior (like the existing
ones do: identical ligand copies, MOL2 input, a charged ligand, a
deliberately broken ligand) over one that exercises many things at once --
it makes failures much easier to diagnose.

## Code style

- Identifiers and function/variable names are in English; comments
  explaining *why* something is done a certain way (not *what* the code
  does) are written in Indonesian throughout the existing codebase --
  that's intentional and fine to continue, no need to translate.
- No comment for anything the code already makes obvious.
- Prefer extending an existing function's generic handling (e.g. "works for
  N ligands") over adding a special case for a specific input.

## Code of conduct

Be respectful and assume good faith. Disagree with the code, not the
person. Reports of abusive or harassing behavior can be sent directly to
the maintainer via the contact info on their GitHub profile.
