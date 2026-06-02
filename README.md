# SEARCH

SEARCH is a CLI-first, config-driven bioinformatics project scaffold for local
CDD/CD-Search preparation, PSSM-based analysis, and mutation-site analysis.

This initial version focuses on project setup, environment initialization,
local package installation, and preparation of the local CDD RPS-BLAST database.
The full pipeline stages will be added in later feature commits.

## Overview

SEARCH is designed as a modular command-line project for running
bioinformatics workflows through configuration files.

The project is organized around three principles:

- **CLI-first execution**: workflows are intended to be executed from the command line.
- **Config-driven workflow design**: pipeline parameters and paths are stored in YAML files.
- **Reproducible local setup**: external databases and environments are prepared through documented setup steps.

At this stage, the repository provides the initial project scaffold:

- Python package metadata through `pyproject.toml`
- Conda environment setup through `environment.yml`
- Initial SEARCH CLI scaffold
- Shared utility modules
- Example configuration layout
- CDD database setup script

## Repository structure

```text
SEARCH/
├── configs/
│   ├── examples/
│   │   └── stages/
│   │       └── example_stage.example.yaml
│   ├── inputs/
│   │   └── .gitkeep
│   └── stages/
│       └── .gitkeep
├── data/
│   ├── external/
│   │   └── .gitkeep
│   └── processed/
│       └── .gitkeep
├── scripts/
│   └── setup_cdd_database.sh
├── src/
│   └── search/
│       ├── __init__.py
│       ├── cli.py
│       ├── exceptions.py
│       └── utils/
│           ├── __init__.py
│           ├── config_utils.py
│           └── fs_utils.py
├── environment.yml
├── pyproject.toml
├── README.md
└── .gitignore
```

The formal pipeline configuration files under `configs/stages/` and the
pipeline implementations under `src/search/pipelines/` and
`src/search/processing/` will be added in later commits.

## 1. Create the Conda environment

The Conda environment is defined in:

```text
environment.yml
```

Create the environment:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate search
```

The environment name is expected to be:

```text
search
```

This is important because the CDD setup script activates the `search`
environment internally before running `makeprofiledb`.

## 2. Install the SEARCH package locally

The Python package metadata and CLI entry point are defined in:

```text
pyproject.toml
```

Install the project in editable mode:

```bash
pip install -e .
```

This makes the `search` command available in the active Conda environment.

Check that the CLI is available:

```bash
search --help
```

For the initial scaffold, you can list registered pipelines:

```bash
search list
```

At this initial stage, no real pipeline stages may be registered yet. The
full pipeline registry will be added together with the actual pipeline
implementations in later feature commits.

## 3. Validate an example config

An example stage configuration is provided at:

```text
configs/examples/stages/example_stage.example.yaml
```

You can validate that the YAML config can be loaded:

```bash
search validate-config --config configs/examples/stages/example_stage.example.yaml
```

This step verifies that:

- the SEARCH CLI is installed correctly
- YAML loading works
- the basic project scaffold is functional

## 4. Prepare the CDD archive

The local CD-Search workflow depends on the NCBI Conserved Domain Database
(CDD). The CDD archive is not tracked by Git because it is a large external
resource.

Create the expected archive directory:

```bash
mkdir -p data/external/cdd/archive
```

Download the CDD archive:

```bash
wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz \
  -O data/external/cdd/archive/cdd.tar.gz
```

The archive must be placed exactly at:

```text
data/external/cdd/archive/cdd.tar.gz
```

The setup script expects this exact path.

## 5. Build the local CDD RPS-BLAST database

After the Conda environment is created and the CDD archive is downloaded,
run:

```bash
bash scripts/setup_cdd_database.sh .
```

The script will:

1. Check that `data/external/cdd/archive/cdd.tar.gz` exists
2. Extract CDD `.smp` profile files into:

   ```text
   data/external/cdd/profiles/
   ```

3. Generate a profile list at:

   ```text
   data/external/cdd/rpsblast/Cdd.pn
   ```

4. Build a local RPS-BLAST database with `makeprofiledb`
5. Save the database files under:

   ```text
   data/external/cdd/rpsblast/
   ```

The final RPS-BLAST database basename is:

```text
data/external/cdd/rpsblast/Cdd
```

Future CD-Search pipeline configuration files should point to this basename,
not to a single database file.

## 6. Data and Git tracking policy

Large external resources and generated outputs should not be committed to Git.

Examples of files and directories that should remain untracked:

```text
data/external/cdd/archive/cdd.tar.gz
data/external/cdd/profiles/
data/external/cdd/rpsblast/
results/
results_tmp/
*.xml
*.tsv
```

The repository should track only lightweight source code, configuration
examples, scripts, and small placeholder files such as `.gitkeep`.

## Development status

This first version is an initialization scaffold.

Included:

- project metadata
- Conda environment definition
- local package installation support
- initial CLI scaffold
- shared utility modules
- example config layout
- CDD setup script

Planned next steps:

- add formal pipeline YAML files under `configs/stages/`
- add pipeline orchestration modules under `src/search/pipelines/`
- add processing logic under `src/search/processing/`
- register executable pipeline stages in the SEARCH CLI

## Initial commit scope

The intended first commit should contain only the project scaffold and setup
files. Formal pipeline stages and analysis logic should be committed
separately in later feature commits.
