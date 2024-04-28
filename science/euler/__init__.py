"""Contains functions for number theory, project Euler, etc."""

__all__ = [
    "numDigits",
    "sumDigits",
    "countDigits",
    "nthDigit",
    "firstNdigits",
    "lastNdigits",
    "reverseNumber",
    "isSquare",
    "isCoprime",
    "isPerfect",
    "isJuf",
    "isPalindrome",
    "isPermutation",
    "isAutomorphic",
    "genLucky",
    "partition",
    "fibonacci",
    "genPanDigital",
    "genPrimTriples",
    "genAllTriples",
    "radix",
    "dec2roman",
    "roman2dec",
    "gcd",
    "leastCommonMultiple",
]

__author__ = "Anton Jansen"

from .euler import *
from .fastmath import *
from .primes import *
from .combinatorics import *
