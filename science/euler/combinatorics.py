from typing import Generator
from math import factorial as _factorial, comb as _comb
from itertools import permutations as _permutations, combinations as _combinations
from sympy.utilities.iterables import multiset_permutations as _multiset_permutations
from .utility import num2list as _num2list, list2num as _list2num


def numperms(array: list, unique: bool = False) -> int:
    """Number of permutations for array. This is n! for non-unique permutations
    and n! / (n1! * n2! * ... * nk!) for unique permutations.

    Args:
        array (list): input.
        unique (bool): unique permutations? Defaults to False.

    Returns:
        int: number of permutations.
    """

    n = _factorial(len(array))

    if unique:
        for count in [array.count(element) for element in set(array)]:
            n //= _factorial(count)

    return n


def genperms(
    array: list, unique: bool = False, asint: bool = False
) -> Generator[tuple[int] | int, None, None]:
    """Generate all permutations of array.
    Note: to go in reverse, simply use array[::-1].

    Args:
        array (list): input.
        unique (bool): only unique permutations.
        asint (bool): return the permutation as an integer.

    Yields:
        tuple or int: permutation.
    """

    if unique:
        for permutation in _multiset_permutations(array):
            yield _list2num(permutation) if asint else tuple(permutation)
    else:
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


def gencombs(
    array: list, r: int, asint: bool = False
) -> Generator[tuple | int, None, None]:
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


def genallcombs(
    array: list, asint: bool = False
) -> Generator[tuple[int] | int, None, None]:
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


def numbersplit(number: int) -> Generator[tuple[int], None, None]:
    """Split a number into all possible combinations of numbers.
    E.g. 1234 --> 1234, 12 34, 1 234, 123 4, 1 23 4, 1 2 3 4.

    Args:
        number (int): number.

    Yields:
        tuple: split.
    """

    n = _num2list(number)

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

        yield tuple(result)


def countPartitions(money: int, coins: list[int] = []) -> int:
    """Returns the number of ways money can be divided (optionally: by coins provided in the coins array).

    Args:
        money (int): starting number.
        coins (list, optional): coins provided in a coin vector. Defaults to [1, 2, ..., money + 1].

    Returns:
        int: number of ways.
    """

    if not coins:
        coins = range(1, money + 1)

    ways = [0] * (money + 1)
    ways[0] = 1

    for ii in range(0, len(coins)):
        for jj in range(coins[ii], money + 1):
            ways[jj] += ways[jj - coins[ii]]

    return ways[money]


def genPartitions(money: int) -> Generator[list[int], None, None]:

    # Base case of recursion: zero is the sum of the empty list.
    if money == 0:
        yield []
        return

    # Modify partitions of n-1 to form partitions of n.
    for part in genPartitions(money - 1):
        yield [1] + part
        if part and (len(part) < 2 or part[1] > part[0]):
            yield [part[0] + 1] + part[1:]
