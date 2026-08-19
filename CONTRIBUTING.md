# Contributing to LeukoQuant

Thank you for your interest in contributing. This document covers how to report bugs, request features, and submit code changes.

---

## Reporting bugs

Use the [Bug Report template](https://github.com/stylianosc/leukoquant/issues/new?template=bug_report.md).

Please include:
- The exact command you ran
- The full error message or log output
- Your OS and Python version (`python3 --version`)
- Whether you are running locally or on an SGE cluster

---

## Requesting features

Use the [Feature Request template](https://github.com/stylianosc/leukoquant/issues/new?template=feature_request.md).

---

## Development setup

```bash
git clone https://github.com/stylianosc/leukoquant.git
cd leukoquant

# Create a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install in editable mode with dev dependencies
pip install -e ".[dev]"
```

---

## Running tests

```bash
pytest leukoquant/tests/ -v
```

Unit tests run without containers or sample data. Integration and smoke tests require additional setup - see the [User Guide](https://stylianosc.github.io/leukoquant/user_guide/) for details.

---

## Code style

This project uses [Black](https://black.readthedocs.io/) for formatting and [Flake8](https://flake8.pycqa.org/) for linting.

```bash
black leukoquant/
flake8 leukoquant/
```

Line length is 100 characters (configured in `pyproject.toml`).

---

## Pull requests

1. Fork the repository and create a feature branch from `main`.
2. Make your changes, add tests where appropriate.
3. Ensure `black` and `flake8` pass with no errors.
4. Open a pull request with a clear description of the change and why it is needed.

For substantial changes, please open an issue first to discuss the approach before investing time in an implementation.

---

## Licence

By contributing you agree that your contributions will be licensed under the [Apache 2.0 License](LICENSE).
