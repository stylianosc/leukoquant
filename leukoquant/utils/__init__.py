"""Utilities for leukoquant.

This module exposes utility entry points without importing submodules at
package import time to avoid import-order/circular-import warnings when
executing submodules with ``python -m``.

The actual implementations are imported lazily when the public wrappers are
called (keeps import side-effects local and deterministic).
"""

from typing import Any

__all__ = ["convert_atlas_to_freesurfer"]


def convert_atlas_to_freesurfer(*args: Any, **kwargs: Any) -> Any:
    """Lazy wrapper for :func:`leukoquant.utils.atlas_converter.convert_atlas_to_freesurfer`.

    Import the heavy submodule only when the function is actually called,
    preventing it from appearing in ``sys.modules`` during package import.
    """
    from .atlas_converter import convert_atlas_to_freesurfer as _convert

    return _convert(*args, **kwargs)
