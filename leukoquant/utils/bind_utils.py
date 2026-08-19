"""Utilities for building Singularity/Apptainer bind-mount lists.

All processors that pass file inputs into a container accumulate bind entries
of the form ``host_path:container_path``.  These helpers:

  1. Build *directory-level* binds (one entry per unique host directory rather
     than one per file), keeping the total entry count proportional to unique
     directories rather than unique files.
  2. Consolidate the resulting list to a single common-ancestor bind per
     filesystem root, reducing even large per-subject lists (800+ subjects ×
     multiple modalities) to a handful of entries that are well within the
     kernel's ARG_MAX limit.

Typical usage
-------------
::

    from leukoquant.utils.bind_utils import (
        dir_level_bind_files,
        dir_level_bind_file_lists,
        consolidate_bind_entries,
    )

    t1_sing, t1_binds  = dir_level_bind_files(t1_files,  "input_t1")
    dwi_sing, dwi_binds = dir_level_bind_file_lists(dwi_lists, "input_dwi")

    # Consolidate ALL per-subject input binds together.
    input_binds, remap = consolidate_bind_entries(t1_binds + dwi_binds)

    # Remap container paths to match the consolidated mount points.
    t1_sing  = [remap(p) for p in t1_sing]
    dwi_sing = [[remap(p) for p in subj] for subj in dwi_sing]

    bind_parts = (
        input_binds
        + [f"{leukoquant_parent}:/leukoquant"]
        + [f"{out_path}:/output"]
    )
    singularity_bind = "--bind " + ",".join(bind_parts)
"""

import os
from collections import defaultdict
from pathlib import Path
from typing import Callable, Dict, List, Tuple


def dir_level_bind_files(
    files: List[str], prefix: str
) -> Tuple[List[str], List[str]]:
    """Return (singularity_paths, bind_entries) binding parent directories.

    Each unique host directory is bound once under a ``/{prefix}_N`` container
    mount point.  All files within that directory share the same container
    directory path.  This keeps the total bind count proportional to unique
    host directories rather than to file count - critical for large cohorts
    where per-file binds produce thousands of entries.

    Empty strings in ``files`` are passed through as empty strings with no
    corresponding bind entry.

    Args:
        files:  Flat list of host file paths (absolute or resolvable).
                Empty strings are treated as absent and skipped.
        prefix: Prefix for container directory names, e.g. ``"input_t1"``.

    Returns:
        singularity_paths: Container paths parallel to ``files``.
        bind_entries:      ``"host_dir:container_dir"`` strings, one per unique dir.
    """
    sing_paths: List[str] = []
    bind_entries: List[str] = []
    dir_map: Dict[str, str] = {}

    for f in files:
        if not f:
            sing_paths.append("")
            continue
        p = Path(f).absolute()
        host_dir = str(p.parent)
        if host_dir not in dir_map:
            container_dir = f"/{prefix}_{len(dir_map)}"
            dir_map[host_dir] = container_dir
            bind_entries.append(f"{host_dir}:{container_dir}")
        sing_paths.append(f"{dir_map[host_dir]}/{p.name}")

    return sing_paths, bind_entries


def dir_level_bind_file_lists(
    file_lists: List[List[str]], prefix: str
) -> Tuple[List[List[str]], List[str]]:
    """Like :func:`dir_level_bind_files` but for nested per-subject file lists.

    The outer list is indexed by subject; each inner list holds one or more
    file paths for that subject's session.  A single shared directory map
    ensures that the same host directory is bound only once even if files from
    multiple subjects reside in the same location.

    Args:
        file_lists: Outer list indexed by subject; inner list contains file paths.
        prefix:     Prefix for container directory names.

    Returns:
        sing_lists:   Nested list mirroring ``file_lists`` with container paths.
        bind_entries: ``"host_dir:container_dir"`` strings, one per unique dir.
    """
    sing_lists: List[List[str]] = []
    bind_entries: List[str] = []
    dir_map: Dict[str, str] = {}

    for files in file_lists:
        sing_list: List[str] = []
        for f in files:
            if not f:
                sing_list.append("")
                continue
            p = Path(f).absolute()
            host_dir = str(p.parent)
            if host_dir not in dir_map:
                container_dir = f"/{prefix}_{len(dir_map)}"
                dir_map[host_dir] = container_dir
                bind_entries.append(f"{host_dir}:{container_dir}")
            sing_list.append(f"{dir_map[host_dir]}/{p.name}")
        sing_lists.append(sing_list)

    return sing_lists, bind_entries


def consolidate_bind_entries(
    bind_entries: List[str],
) -> Tuple[List[str], Callable[[str], str]]:
    """Reduce per-directory bind entries to one common-ancestor bind per filesystem root.

    All host directories that share a common prefix are consolidated into a
    single ``ancestor_dir:container_root`` bind.  Container paths are remapped
    by the returned ``remap`` callable so existing singularity path lists stay
    consistent after consolidation.

    Multiple filesystem roots (e.g. different drives on Windows, or truly
    unrelated top-level paths) are handled by grouping entries and assigning
    ``/input_N`` per group (``/input`` for the single-group case).

    Passing an empty list is safe and returns ``([], identity)``.

    Args:
        bind_entries: ``"host_dir:container_dir"`` strings.  The host side must
                      be an absolute directory path (i.e. produced by
                      :func:`dir_level_bind_files` or equivalent).

    Returns:
        consolidated_entries: Reduced bind list (one per filesystem root group).
        remap:                 ``remap(old_container_path) -> new_container_path``.
                               Paths not covered by the consolidation are returned
                               unchanged, so it is always safe to apply to any path.
    """
    if not bind_entries:
        return bind_entries, lambda p: p

    host_dirs = [e.split(":", 1)[0] for e in bind_entries]

    # Group by filesystem root: drive letter on Windows, "/" on Linux/macOS.
    root_groups: defaultdict[str, List[str]] = defaultdict(list)
    for d in host_dirs:
        drive = os.path.splitdrive(d)[0] or "/"
        root_groups[drive].append(d)

    old_to_new: Dict[str, str] = {}
    new_bind_entries: List[str] = []

    for group_idx, (_, dirs_in_group) in enumerate(sorted(root_groups.items())):
        ancestor = os.path.commonpath(dirs_in_group)
        container_root = "/input" if len(root_groups) == 1 else f"/input_{group_idx}"
        new_bind_entries.append(f"{ancestor}:{container_root}")

        for entry in bind_entries:
            host_dir, old_container_dir = entry.split(":", 1)
            if host_dir not in dirs_in_group:
                continue
            rel = os.path.relpath(host_dir, ancestor)
            new_container_dir = (
                container_root if rel == "." else f"{container_root}/{rel}"
            )
            old_to_new[old_container_dir] = new_container_dir

    def remap(path: str) -> str:
        """Translate an old container path to its new consolidated equivalent."""
        if not path:
            return path
        for old_dir, new_dir in old_to_new.items():
            if path.startswith(old_dir + "/"):
                return new_dir + path[len(old_dir):]
            if path == old_dir:
                return new_dir
        return path

    return new_bind_entries, remap
