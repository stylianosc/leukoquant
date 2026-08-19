"""
container_utils.py

Checks whether a Singularity .sif container exists on disk.
If not, downloads it from Hugging Face Hub automatically.

Usage (at the top of any workflow .smk file or processor):

    from leukoquant.utils.container_utils import (
        ensure_container,
        MINICONDA_SIF_FILENAME,
        FREESURFER_SIF_FILENAME,
        MINICONDA_SIF_HF_PATH,
        FREESURFER_SIF_HF_PATH,
    )

    ensure_container(
        sif_path="/path/to/containers/miniconda_unified_container.sif",
        filename=MINICONDA_SIF_HF_PATH,
    )
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

# Hugging Face dataset repo hosting the pre-built .sif files
HF_REPO = "stylianosc/leukoquant"

# Canonical paths inside the HF repo (subfolder/filename)
MINICONDA_SIF_HF_PATH  = "containers/miniconda_unified_container.sif"
FREESURFER_SIF_HF_PATH = "containers/freesurfer_unified_container.sif"

# Canonical local filenames (basename only)
MINICONDA_SIF_FILENAME  = "miniconda_unified_container.sif"
FREESURFER_SIF_FILENAME = "freesurfer_unified_container.sif"

# Dedicated Hugging Face dataset repo hosting the GIF anatomical atlas
# database (T1/FLAIR/labels, MiDeFace-de-identified). Private for now; will
# be flipped public once data-sharing/ethics governance for redistribution
# is confirmed.
GIF_DB_HF_REPO = "stylianosc/gif-database"

# Tarball filenames inside GIF_DB_HF_REPO. Each extracts to a directory of
# the same basename (minus .tar.gz), matching the GIF database layout
# leukoquant expects (db.xml + labels.xml + GroupMask.nii.gz + image folder
# + labels/).
GIF_DB_T1_FILENAME    = "db_mideface.tar.gz"
GIF_DB_FLAIR_FILENAME = "db_FLAIR_mideface.tar.gz"

# How long a caller will wait for another process's in-progress download
# before giving up. Generous, since large containers can legitimately take
# several minutes; also bounds how long a genuinely stale lock (left behind
# by a crashed process) blocks everyone else.
_LOCK_WAIT_TIMEOUT_S = 900
_LOCK_POLL_INTERVAL_S = 2


def ensure_container(
    sif_path: str,
    hf_repo: str = HF_REPO,
    filename: str | None = None,
) -> None:
    """
    Ensure the Singularity .sif file exists at `sif_path`.

    If the file is missing, download it from `hf_repo` on Hugging Face.

    Safe under concurrent invocation: when hundreds of SGE tasks launch at
    once (e.g. a fresh multi-dataset run), they can all call this function
    for the same shared container within the same few seconds. The
    unlocked check-then-download-then-move sequence this replaced let
    concurrent callers race: a `singularity exec` running in one process
    could observe `sif_path` mid-write by another (missing or truncated),
    which breaks that container's internal bind-mounts for that one
    invocation - e.g. a `ModuleNotFoundError` for a package that's only
    importable inside the container.

    Two things make this safe:
      1. An exclusive lock file (`sif_path + ".lock"`), created via
         O_CREAT|O_EXCL, serialises the download - only one process
         downloads, everyone else waits and then finds the file already
         present. O_EXCL-based lockfiles are reliable on NFS in a way
         `fcntl.flock()` historically is not.
      2. The file is downloaded to a temp path inside the SAME directory
         as `sif_path`, then placed with `os.replace()`, which is an
         atomic rename within one filesystem - no reader ever observes a
         partially-written file at `sif_path`.

    Parameters
    ----------
    sif_path:
        Absolute or relative path where the .sif file should exist.
    hf_repo:
        Hugging Face dataset repo ID, e.g. "stylianosc/leukoquant".
    filename:
        Path inside the HF repo including subfolder,
        e.g. "containers/miniconda_unified_container.sif".
        Defaults to "containers/<basename of sif_path>".
    """
    if os.path.isfile(sif_path):
        return  # already present, nothing to do

    dest_dir = os.path.dirname(sif_path) or "."
    os.makedirs(dest_dir, exist_ok=True)

    if filename is None:
        filename = "containers/" + os.path.basename(sif_path)

    lock_path = sif_path + ".lock"
    _acquire_lock(lock_path, lambda: os.path.isfile(sif_path))
    try:
        # Re-check now that we hold the lock: another process may have
        # already finished the download while we were waiting for it.
        if os.path.isfile(sif_path):
            return
        _download_container(sif_path, dest_dir, hf_repo, filename)
    finally:
        try:
            os.remove(lock_path)
        except OSError:
            pass  # already removed, or never fully created - not fatal


def _acquire_lock(lock_path: str, is_ready) -> None:
    """Block until this process holds the exclusive download lock.

    Polls with O_CREAT|O_EXCL rather than fcntl.flock(), since flock()'s
    behaviour on NFS-mounted files is not consistently reliable across
    clients, whereas O_EXCL file creation is.

    `is_ready` is a zero-arg callable checked on every poll: if another
    process already finished downloading the target artifact (a .sif file
    for `ensure_container()`, a GIF db.xml for `ensure_gif_db()`) while we
    were waiting, there's nothing left to wait for.
    """
    deadline = time.monotonic() + _LOCK_WAIT_TIMEOUT_S
    while True:
        try:
            fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, str(os.getpid()).encode())
            os.close(fd)
            return
        except FileExistsError:
            pass

        if is_ready():
            return

        if time.monotonic() >= deadline:
            raise TimeoutError(
                f"[container_utils] Timed out after {_LOCK_WAIT_TIMEOUT_S}s waiting for "
                f"download lock: {lock_path}\n"
                f"[container_utils] If no other process is actually downloading "
                f"right now (check for a crashed job), delete this lock file and retry."
            )
        time.sleep(_LOCK_POLL_INTERVAL_S)


def _download_container(sif_path: str, dest_dir: str, hf_repo: str, filename: str) -> None:
    print(
        f"[container_utils] Container not found: {sif_path}\n"
        f"[container_utils] Downloading '{filename}' from "
        f"Hugging Face repo '{hf_repo}' ...",
        flush=True,
    )

    try:
        from huggingface_hub import hf_hub_download
        import logging as _logging
        # Suppress huggingface_hub and httpx verbose logs
        _logging.getLogger("huggingface_hub").setLevel(_logging.ERROR)
        _logging.getLogger("huggingface_hub.utils._http").setLevel(_logging.ERROR)
        _logging.getLogger("httpx").setLevel(_logging.ERROR)
    except ImportError:
        print(
            "[container_utils] ERROR: 'huggingface_hub' is not installed.\n"
            "Install it with:  pip install huggingface_hub",
            file=sys.stderr,
        )
        sys.exit(1)

    import warnings
    # Suppress the unauthenticated HF token warning
    os.environ.setdefault("HF_HUB_DISABLE_PROGRESS_BARS", "0")
    warnings.filterwarnings("ignore", message=".*unauthenticated.*")

    sif_name = os.path.basename(sif_path)
    print(f"[container_utils] Downloading {sif_name} - this may take a while for large containers...", flush=True)

    # Download into a per-call temp directory INSIDE dest_dir (not directly
    # to sif_path, and not relying on hf_hub_download's own cache location,
    # which may be on a different filesystem to dest_dir) so the final
    # placement below is guaranteed to be a same-filesystem, atomic rename.
    tmp_dir = tempfile.mkdtemp(dir=dest_dir, prefix=".hf_download_")
    try:
        downloaded = hf_hub_download(
            repo_id=hf_repo,
            filename=filename,
            repo_type="dataset",
            local_dir=tmp_dir,
            token=os.environ.get("HF_TOKEN", None),
        )

        if os.path.islink(downloaded):
            # hf_hub_download symlinked into its shared cache (e.g.
            # ~/.cache/huggingface), which may be a different filesystem
            # than dest_dir. Copy the real bytes in before the atomic
            # rename, since os.replace() requires same-filesystem paths.
            real_copy = os.path.join(tmp_dir, "container.sif")
            shutil.copyfile(os.path.realpath(downloaded), real_copy)
            downloaded = real_copy

        os.replace(downloaded, sif_path)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    print(
        f"[container_utils] Download complete: {sif_path}",
        flush=True,
    )


def ensure_gif_db(db_dir: str, filename: str, hf_repo: str = GIF_DB_HF_REPO) -> None:
    """
    Ensure a GIF atlas database directory exists at `db_dir` (i.e. that
    `db_dir/db.xml` is present).

    If missing, downloads a `.tar.gz` named `filename` from `hf_repo` on
    Hugging Face and extracts it in place. The tarball's top-level entry
    must be a single directory sharing `db_dir`'s basename (e.g.
    `db_mideface.tar.gz` extracts to a top-level `db_mideface/`) - this is
    how the GIF atlas tarballs are packaged.

    Mirrors `ensure_container()`'s concurrency-safety approach (O_CREAT|O_EXCL
    lock + atomic same-filesystem placement via `os.replace()`), adapted for
    a directory-of-files payload instead of a single file: the tarball is
    downloaded and extracted into temp directories inside `db_dir`'s parent
    (same filesystem as the final destination), then the extracted directory
    is moved into place with a single atomic rename - no reader ever
    observes a partially-extracted `db_dir`.

    Parameters
    ----------
    db_dir:
        Absolute or relative path where the GIF database directory should
        exist, e.g. ".../GIF/db_mideface".
    filename:
        Tarball filename inside the HF repo, e.g. "db_mideface.tar.gz".
    hf_repo:
        Hugging Face dataset repo ID. Defaults to the dedicated GIF database repo.
    """
    db_path = Path(db_dir)
    db_xml = db_path / "db.xml"
    if db_xml.is_file():
        return  # already present, nothing to do

    parent_dir = db_path.parent
    parent_dir.mkdir(parents=True, exist_ok=True)

    lock_path = str(db_path) + ".lock"
    _acquire_lock(lock_path, lambda: db_xml.is_file())
    try:
        # Re-check now that we hold the lock: another process may have
        # already finished the download while we were waiting for it.
        if db_xml.is_file():
            return
        _download_gif_db(db_path, parent_dir, hf_repo, filename)
    finally:
        try:
            os.remove(lock_path)
        except OSError:
            pass  # already removed, or never fully created - not fatal


def _download_gif_db(db_path: Path, parent_dir: Path, hf_repo: str, filename: str) -> None:
    print(
        f"[container_utils] GIF database not found: {db_path}\n"
        f"[container_utils] Downloading '{filename}' from "
        f"Hugging Face repo '{hf_repo}' ...",
        flush=True,
    )

    try:
        from huggingface_hub import hf_hub_download
        import logging as _logging
        _logging.getLogger("huggingface_hub").setLevel(_logging.ERROR)
        _logging.getLogger("huggingface_hub.utils._http").setLevel(_logging.ERROR)
        _logging.getLogger("httpx").setLevel(_logging.ERROR)
    except ImportError:
        print(
            "[container_utils] ERROR: 'huggingface_hub' is not installed.\n"
            "Install it with:  pip install huggingface_hub",
            file=sys.stderr,
        )
        sys.exit(1)

    import tarfile
    import warnings
    os.environ.setdefault("HF_HUB_DISABLE_PROGRESS_BARS", "0")
    warnings.filterwarnings("ignore", message=".*unauthenticated.*")

    print(f"[container_utils] Downloading {filename} - this may take a while for large archives...", flush=True)

    tmp_download_dir = tempfile.mkdtemp(dir=parent_dir, prefix=".hf_download_")
    try:
        downloaded = hf_hub_download(
            repo_id=hf_repo,
            filename=filename,
            repo_type="dataset",
            local_dir=tmp_download_dir,
            token=os.environ.get("HF_TOKEN", None),
        )

        if os.path.islink(downloaded):
            # hf_hub_download symlinked into its shared cache, which may be
            # on a different filesystem than parent_dir. Copy the real bytes
            # in before extraction, since the tarfile module needs a real
            # (or at least readable) file, and the later os.replace() below
            # requires the extracted directory to share parent_dir's filesystem.
            real_copy = os.path.join(tmp_download_dir, os.path.basename(filename))
            shutil.copyfile(os.path.realpath(downloaded), real_copy)
            downloaded = real_copy

        tmp_extract_dir = tempfile.mkdtemp(dir=parent_dir, prefix=".hf_extract_")
        try:
            with tarfile.open(downloaded, "r:gz") as tf:
                tf.extractall(tmp_extract_dir)

            extracted_root = Path(tmp_extract_dir) / db_path.name
            if not (extracted_root / "db.xml").is_file():
                raise RuntimeError(
                    f"[container_utils] Downloaded archive '{filename}' did not extract to the "
                    f"expected top-level directory '{db_path.name}/' containing db.xml. "
                    f"Found instead: {sorted(p.name for p in Path(tmp_extract_dir).iterdir())}"
                )

            os.replace(extracted_root, db_path)
        finally:
            shutil.rmtree(tmp_extract_dir, ignore_errors=True)
    finally:
        shutil.rmtree(tmp_download_dir, ignore_errors=True)

    print(
        f"[container_utils] Download complete: {db_path}",
        flush=True,
    )
