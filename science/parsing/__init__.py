"""Contains function for parsing data and structure (pdb, gro) files."""

__all__ = ['loadxvg', 'loadCol', 'loadVal', 'Structure']
__author__ = 'Anton Jansen'

from .parsing import loadxvg
from .parsing import loadCol
from .parsing import loadVal
from .structure import Structure
