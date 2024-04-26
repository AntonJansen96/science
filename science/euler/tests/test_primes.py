#!/usr/bin/env python3

import pytest
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
    assert primes._sieveArray[-1] == 67


def test_nextprime():
    prime = primes.nextprime()
    assert next(prime) == 2
    assert next(prime) == 3
    assert next(prime) == 5
    prime = primes.nextprime(3)
    assert next(prime) == 3
    prime = primes.nextprime(96)
    assert next(prime) == 97


def test_prevprime():
    prime = primes.prevprime(2)
    with pytest.raises(StopIteration):
        next(prime)
    prime = primes.prevprime(3)
    assert next(prime) == 2
    with pytest.raises(StopIteration):
        next(prime)
    prime = primes.prevprime(5)
    assert next(prime) == 3
    assert next(prime) == 2
    with pytest.raises(StopIteration):
        next(prime)


def test_sieve():
    assert Primes.sieve(40) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]


def test_iscircularprime():
    assert Primes.iscircularprime(197) == True
    assert Primes.iscircularprime(227) == False


def test_primerange():
    assert Primes.primerange(200, 250) == [211, 223, 227, 229, 233, 239, 241]


def test_gentwinprimes():
    twinprime = primes.gentwinprimes()
    assert next(twinprime) == (3, 5)
    assert next(twinprime) == (5, 7)
    assert next(twinprime) == (11, 13)
    assert next(twinprime) == (17, 19)


def test_istwinprime():
    assert primes.istwinprime(3) == True
    assert primes.istwinprime(5) == True
    assert primes.istwinprime(97) == False
