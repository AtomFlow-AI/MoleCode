# Dependencies

The skill wraps the repository's `molecode` package; the CLI adds the repo root
to `sys.path` automatically, so no install is required *if* the runtime already
has the scientific dependencies.

## Required Python

Python 3.9 or newer.

## Runtime dependencies

| Package | Needed for |
| --- | --- |
| `rdkit` | every SMILES/PSMILES conversion and `validate` |
| `networkx` | the `compare` command (markush graph isomorphism) |

Pure graph text editing needs nothing; conversion to/from SMILES needs RDKit.

## If dependencies are missing

`python scripts/molecode_convert.py doctor` reports what is present. Prefer the
least-invasive option, and **ask the user before installing into their active
environment**:

1. Use a Python environment that already has `rdkit` + `networkx`.
2. Conda/mamba:

   ```bash
   conda install -c conda-forge rdkit networkx
   ```

3. Pip, when wheels exist for the platform:

   ```bash
   python -m pip install rdkit networkx
   ```

4. Install the repository itself (pulls both via `pyproject.toml`):

   ```bash
   python -m pip install -e .      # from the repository root
   ```

5. Isolated local virtualenv (does not touch the user's env):

   ```bash
   python -m venv .molecode-venv
   .molecode-venv/bin/python -m pip install rdkit networkx
   .molecode-venv/bin/python scripts/molecode_convert.py doctor
   ```

Do not silently install packages into the user's active environment.
