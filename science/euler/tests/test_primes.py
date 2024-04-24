#!/usr/bin/env python3

from science.euler.primes import Primes


def test_init():
    primes = Primes()
    assert primes._nPrimes == 11
    assert primes._max == 1000
    assert primes._sieveArray == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]


def test_expand():
    primes = Primes()
    primes._expand()
    assert primes._max == 4000
    assert len(primes._sieveArray) == 19
    assert primes._sieveArray[-1] == 61


def test_sieve():
    primes = Primes()
    assert primes.sieve(40) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
