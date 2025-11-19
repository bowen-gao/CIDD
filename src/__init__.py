"""
CIDD: Collaborative Intelligence for Structure-Based Drug Design Empowered by LLMs

This package provides tools for AI-guided molecular generation and optimization
using large language models combined with structure-based drug design methods.
"""

__version__ = "1.0.0"
__author__ = "Bowen Gao"
__email__ = "bowen-gao@example.com"  # Update with actual email
__license__ = "MIT"

from . import generation
from . import evaluation
from . import utils

__all__ = ["generation", "evaluation", "utils"]