"""Shared utilities for subject resolution across all processors."""

import glob
import os
from pathlib import Path
from typing import List, Union


def read_subjects(subject_input: str) -> List[str]:
    """Return a list of subject IDs from a bare ID or a text file (one per line)."""
    p = Path(subject_input)
    if p.exists() and p.is_file():
        return [line.strip() for line in p.read_text().splitlines() if line.strip()]
    return [subject_input]


def resolve_subject_pattern(
    pattern: str, subject: str, all: bool = False
) -> Union[str, List[str]]:
    """Replace {subject} in a glob pattern and return match(es).

    When all=False (default): returns the first sorted match as a string;
    raises FileNotFoundError if nothing matches - existing behaviour unchanged.

    When all=True: returns a sorted list of all matching absolute paths.
    Returns an empty list (does not raise) when nothing matches, so callers
    can decide whether to warn/skip or raise an error.

    If no {subject} placeholder is present but the pattern contains glob
    characters, the glob is still resolved (e.g. bare paths like I*.nii.gz).
    If no glob characters and no {subject}, the path is returned as-is (literal).
    """
    if "{subject}" in pattern:
        resolved = pattern.format(subject=subject)
    elif any(c in pattern for c in ("*", "?", "[")):
        resolved = pattern
    else:
        # Literal path - no substitution or glob needed.
        p = os.path.abspath(pattern)
        if all:
            return [p] if os.path.exists(p) else []
        return p

    candidates = sorted(glob.glob(resolved))
    if all:
        return [os.path.abspath(c) for c in candidates]
    if not candidates:
        raise FileNotFoundError(
            f"No files found for subject '{subject}' with pattern: {pattern}"
        )
    return os.path.abspath(candidates[0])
