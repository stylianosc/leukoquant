"""
Core functionality for leukoquant.

This module contains the main processing functions for:
- Image preprocessing and registration
- Tractography analysis
- Lesion segmentation
- Metrics extraction
- Quality control
- GIF/BaMoS processing
"""

# Import only the modules that currently exist
# In leukoquant/core/__init__.py
from .gif_processor import apply_gif, GIFProcessor
from .bamos_processor import apply_bamos, BaMoSProcessor

__all__ = [
    'apply_gif', 'GIFProcessor', 'apply_bamos', 'BaMoSProcessor'
]
