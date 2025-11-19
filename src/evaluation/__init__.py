"""
CIDD Evaluation Module

This module contains functions for evaluating generated molecules,
including docking scores, drug-likeness metrics, and other properties.
"""

try:
    from .docking_vina import VinaDockingTask
    from .docking_qvina import QVinaDockingTask
    from .scoring_func import *
    __all__ = ["VinaDockingTask", "QVinaDockingTask"]
except ImportError:
    # Some modules may not be available if optional dependencies are not installed
    __all__ = []