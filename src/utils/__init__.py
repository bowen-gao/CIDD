"""
CIDD Utilities Module

This module contains utility functions for data processing and analysis.
"""

try:
    from .unimap_ret import ret_fragments
    __all__ = ["ret_fragments"]
except ImportError:
    __all__ = []