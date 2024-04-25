#!/usr/bin/env python3

from science.euler.primes import Primes

primes = Primes()


def test_init():

    assert primes._nPrimes == 11
    assert primes._max == 1000
    assert primes._sieveArray == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]


def test_expand():
    primes._expand()
    assert primes._max == 4000
    assert len(primes._sieveArray) == 19
    assert primes._sieveArray[-1] == 61


def test_sieve():
    assert Primes.sieve(40) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]


def test_genprime():
    prime = Primes.genprime()
    assert next(prime) == 2
    assert next(prime) == 3
    assert next(prime) == 5
    assert next(prime) == 7
    assert next(prime) == 11


def test_iscircularprime():
    assert Primes.iscircularprime(197) == True
    assert Primes.iscircularprime(227) == False
