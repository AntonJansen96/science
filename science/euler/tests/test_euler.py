#!/usr/bin/env python3

from science.euler.euler import isAutomorphic


def test_isAutomorphic():
    assert isAutomorphic(76) == True
    assert isAutomorphic(20) == False
