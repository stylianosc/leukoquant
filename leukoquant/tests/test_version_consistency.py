"""Every place that states the package version must state the same one.

Version 2.3.0 shipped with the README citation block still saying 2.2.0 until this check existed."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _find(path, pattern):
    f = ROOT / path
    if not f.exists():          # an installed copy has no repository files around it
        return None
    m = re.search(pattern, f.read_text(), re.M)
    assert m, f"{path}: no version found"
    return m.group(1)


def test_all_version_statements_agree():
    versions = {
        "pyproject.toml": _find("pyproject.toml", r'^version\s*=\s*"([^"]+)"'),
        "setup.py": _find("setup.py", r'version\s*=\s*"([^"]+)"'),
        "CITATION.cff": _find("CITATION.cff", r"^version:\s*(\S+)"),
        "CITATION.bib": _find("CITATION.bib", r"version\s*=\s*\{([^}]+)\}"),
        "README.md": _find("README.md", r"version\s*=\s*\{([^}]+)\}"),
    }
    versions = {k: v for k, v in versions.items() if v is not None}
    assert len(set(versions.values())) <= 1, versions
