from math import factorial as _factorial
from math import comb as _comb
from itertools import permutations as _permutations
from itertools import combinations as _combinations

from .euler import num2vec as _num2vec


# Helper function.
def _list2num(array: list) -> int:
    """Convert a list to a number.

    Args:
        array (list): input.

    Returns:
        int: number.
    """

    number = 0
    for digit in array:
        number = number * 10 + digit

    return number


def numperms(array: list):
    """Number of permutations for array. This is len(array)!.

    Args:
        array (list): input.

    Returns:
        int: number of permutations.
    """

    return _factorial(len(array))


def genperms(array: list, asint: bool = False):
    """Generate all permutations of array.
    Note: to go in reverse, simply use array[::-1].

    Args:
        array (list): input.
        asint (bool): return the permutation as an integer.

    Yields:
        tuple or int: permutation.
    """

    for permutation in _permutations(array):
        yield _list2num(permutation) if asint else permutation


def numcombs(array: list, r: int) -> int:
    """Number of combinations for array. This is binom(len(array), r).

    Args:
        array (list): input.
        r (int): number of elements to choose.

    Returns:
        int: number of combinations.
    """

    return _comb(len(array), r)


def gencombs(array: list, r: int, asint: bool = False):
    """Generate combinations of array.

    Args:
        array (list): input.
        r (int): number of elements to choose.
        asint (bool): return the combination as an integer.

    Yields:
        tuple or int: combination.
    """

    for combination in _combinations(array, r):
        yield _list2num(combination) if asint else combination


def genallcombs(array: list, asint: bool = False):
    """Generate all combinations of array.
    Note: also generates the identity combination () (0 if asint=True).

    Args:
        array (list): input.
        asint (bool): return the combination as an integer.

    Yields:
        tuple or int: combination.
    """

    for r in range(len(array) + 1):
        for combination in _combinations(array, r):
            yield _list2num(combination) if asint else combination


def numbersplit(number: int):
    """Split a number into all possible combinations of numbers.
    E.g. 1234 --> 1234, 12 34, 1 234, 123 4, 1 23 4, 1 2 3 4.

    Args:
        number (int): number.

    Yields:
        list: split.
    """

    n = _num2vec(number)

    for i in range(1 << (len(n) - 1)):
        result = []
        temp = [n[0]]

        for j in range(len(n) - 1):

            if (i >> j) & 1:
                result.append(_list2num(temp))
                temp = [n[j + 1]]

            else:
                temp.append(n[j + 1])

        result.append(_list2num(temp))

        yield result