# SEARCH

This is the official codebase for **SEARCH**.

> Project type: **cli-tool**
> Python package: **search**
> CLI command: **search**

## Overview

This project follows a CLI-first, config-driven, and reproducible project structure.

It is designed for command-line tools, data processing pipelines, and reproducible workflow execution.

## Installation

Create the environment:

```bash
conda env create -f environment.yml
conda activate search
```

Install the package in editable mode:

```bash
pip install -e .[dev]
```

Or use:

```bash
make install
```

## Quick Start

Show help:

```bash
search --help
```

Run a pipeline stage:

```bash
search run <stage_name> --config configs/stages/<stage_name>.yaml
```

Example:

```bash
search run example.stage --config configs/examples/stages/example_stage.example.yaml
```

Or use:

```bash
make run-example
```

## Project Structure

```text
configs/      Configuration files and stage examples
data/         Data lifecycle directories
experiments/  Reproducible run outputs
results/      Aggregated outputs, reports, figures, and tables
src/          Source code package
tests/        Unit and integration tests
docs/         Documentation
notebooks/    Exploratory notebooks
```

## Development

Run tests:

```bash
make test
```

Run lint checks:

```bash
make lint
```

Format code:

```bash
make format
```

Clean cache files:

```bash
make clean
```

## Documentation

See the `docs/` directory for usage, configuration, pipeline, output, and development notes.

## Contributing

Contributions are welcome. Please open an issue or submit a pull request.

## Acknowledgements

We sincerely thank the authors of the open-source projects used in this work.
