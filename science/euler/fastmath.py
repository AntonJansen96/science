from ._fastpow import fastpow


def intlog10(num: int) -> int:
    """Returns the integer part of the base 10 logarithm of a number.

    Args:
        num (int): number.

    Returns:
        int: integer part of the logarithm.
    """

    n = 0

    while num >= 10000000000000000:
        n += 16
        num //= 10000000000000000

    if num >= 100000000:
        n += 8
        num //= 100000000

    if num >= 10000:
        n += 4
        num //= 10000

    if num >= 100:
        n += 2
        num //= 100

    if num >= 10:
        n += 1
        num //= 10

    return n
