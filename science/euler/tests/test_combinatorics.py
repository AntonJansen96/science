from science.euler.combinatorics import *


def test_numperms():
    assert numperms([1]) == 1
    assert numperms([1, 2, 3]) == 6
    assert numperms([1, 2, 3, 4]) == 24


def test_genperms():
    # fmt: off
    assert list(genperms([1, 2, 3])) == [(1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)]
    assert list(genperms([1, 2, 3], asint=True)) == [123, 132, 213, 231, 312, 321]
    assert list(genperms([0, 0, 1])) == [(0, 0, 1), (0, 1, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0), (1, 0, 0)]
    assert list(genperms([0, 0, 1], unique=True)) == [(0, 0, 1), (0, 1, 0), (1, 0, 0)]
    # fmt: on


def test_numcombs():
    assert numcombs([1, 2, 3], 2) == 3
    assert numcombs([1, 2, 3, 4], 2) == 6


def test_gencombs():
    # fmt: off
    assert list(gencombs([1, 2, 3, 4], r=2)) == [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    assert list(gencombs([1, 2, 3, 4], r=3, asint=True)) == [123, 124, 134, 234]
    # fmt: on


def test_genallcombs():
    # fmt: off
    assert list(genallcombs([1, 2, 3])) == [(), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]
    # fmt: on


def test_numbersplit():
    # fmt: off
    assert list(numbersplit(1234)) == [(1234,), (1, 234), (12, 34), (1, 2, 34), (123, 4), (1, 23, 4), (12, 3, 4), (1, 2, 3, 4)]
    # fmt: on


def test_countPartitions():
    assert countPartitions(5) == 7
    assert countPartitions(6) == 11


def test_genPartitions():
    # fmt: off
    assert list(genPartitions(5)) == [[1, 1, 1, 1, 1], [1, 1, 1, 2], [1, 2, 2], [1, 1, 3], [2, 3], [1, 4], [5]]
    # fmt: on
