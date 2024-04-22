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
    "partition",
    "fibonacci",
    "genPanDigital",
    "genPrimTriples",
    "genAllTriples",
    "num2vec",
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
