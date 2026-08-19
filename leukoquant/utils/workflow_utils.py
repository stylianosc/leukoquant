"""Workflow path-translation helpers shared by Snakemake workflows.

These helpers mirror the logic used in the z-score workflow but live in
a dedicated utils module so the Snakemake files remain declarative.
"""
from pathlib import Path
from typing import Dict

from leukoquant.utils.z_score_utils import translate_path as _translate_path


def translate_path(path: str, bind_map: Dict[str, str]) -> str:
    return _translate_path(path, bind_map)


def translate_pattern(pattern: str, bind_map: Dict[str, str]) -> str:
    """Translate the base part of a base:glob pattern to container path.

    If pattern does not contain a ':', translate the whole string.
    """
    if not pattern:
        return pattern
    if ":" in pattern:
        base, glob_part = pattern.split(":", 1)
        base_t = _translate_path(base, bind_map)
        return f"{base_t}:{glob_part}"
    return _translate_path(pattern, bind_map)


def translate_mask_spec(mask_spec: str, bind_map: Dict[str, str]) -> str:
    """Translate each mask mapping's base path in a mask_spec string.

    mask_spec format: name=base:glob[:space],comma-separated.
    """
    if not mask_spec:
        return ""
    parts = []
    for item in [seg.strip() for seg in mask_spec.split(",") if seg.strip()]:
        if "=" not in item:
            parts.append(item)
            continue
        name, rest = item.split("=", 1)
        # preserve optional ':space' suffix
        if ":" in rest:
            pattern_part, suffix = rest.rsplit(":", 1)
            translated = translate_pattern(pattern_part, bind_map)
            parts.append(f"{name}={translated}:{suffix}")
        else:
            translated = translate_pattern(rest, bind_map)
            parts.append(f"{name}={translated}")
    return ",".join(parts)

