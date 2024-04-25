def mulmod(a: int, b: int, m: int) -> int:
    """Returns the product of two numbers modulo m, i.e. (a * b) % m = a * b % m.

    Args:
        a (int): number.
        b (int): number.
        m (int): modulus.

    Returns:
        int: product modulo m.
    """

    # This is already highly optimized in Python under the hood.
    return a * b % m


def powmod(a: int, b: int, m: int) -> int:
    """Returns the power of a number modulo m, i.e. a^b % m.

    Args:
        a (int): number.
        b (int): exponent.
        m (int): modulus.

    Returns:
        int: power modulo m.
    """

    # This is already highly optimized in Python under the hood.
    return pow(a, b, m)


def modinverse(a: int, m: int) -> int:
    """Returns the modular inverse of a number.

    Args:
        a (int): number.
        m (int): modulus.

    Returns:
        int: modular inverse.
    """

    # This is already highly optimized in Python under the hood.
    return pow(a, -1, m)


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


def fastpow(n: int, m: int) -> int:
    """Exponentiate num [0, 10] with power.
    Note: result must be smaller than int64, if it is not, return 0.

    Args:
        n (int): base.
        m (int): exponent.

    Returns:
        int: result.
    """

    power2 = [
        1,
        2,
        4,
        8,
        16,
        32,
        64,
        128,
        256,
        512,
        1024,
        2048,
        4096,
        8192,
        16384,
        32768,
        65536,
        131072,
        262144,
        524288,
        1048576,
        2097152,
        4194304,
        8388608,
        16777216,
        33554432,
        67108864,
        134217728,
        268435456,
        536870912,
        1073741824,
        2147483648,
        4294967296,
        8589934592,
        17179869184,
        34359738368,
        68719476736,
        137438953472,
        274877906944,
        549755813888,
        1099511627776,
        2199023255552,
        4398046511104,
        8796093022208,
        17592186044416,
        35184372088832,
        70368744177664,
        140737488355328,
        281474976710656,
        562949953421312,
        1125899906842624,
        2251799813685248,
        4503599627370496,
        9007199254740992,
        18014398509481984,
        36028797018963968,
        72057594037927936,
        144115188075855872,
        288230376151711744,
        576460752303423488,
        1152921504606846976,
        2305843009213693952,
        4611686018427387904,
    ]

    power3 = [
        1,
        3,
        9,
        27,
        81,
        243,
        729,
        2187,
        6561,
        19683,
        59049,
        177147,
        531441,
        1594323,
        4782969,
        14348907,
        43046721,
        129140163,
        387420489,
        1162261467,
        3486784401,
        10460353203,
        31381059609,
        94143178827,
        282429536481,
        847288609443,
        2541865828329,
        7625597484987,
        22876792454961,
        68630377364883,
        205891132094649,
        617673396283947,
        1853020188851841,
        5559060566555523,
        16677181699666569,
        50031545098999707,
        150094635296999121,
        450283905890997363,
        1350851717672992089,
        4052555153018976267,
    ]

    power4 = [
        1,
        4,
        16,
        64,
        256,
        1024,
        4096,
        16384,
        65536,
        262144,
        1048576,
        4194304,
        16777216,
        67108864,
        268435456,
        1073741824,
        4294967296,
        17179869184,
        68719476736,
        274877906944,
        1099511627776,
        4398046511104,
        17592186044416,
        70368744177664,
        281474976710656,
        1125899906842624,
        4503599627370496,
        18014398509481984,
        72057594037927936,
        288230376151711744,
        1152921504606846976,
        4611686018427387904,
    ]

    power5 = [
        1,
        5,
        25,
        125,
        625,
        3125,
        15625,
        78125,
        390625,
        1953125,
        9765625,
        48828125,
        244140625,
        1220703125,
        6103515625,
        30517578125,
        152587890625,
        762939453125,
        3814697265625,
        19073486328125,
        95367431640625,
        476837158203125,
        2384185791015625,
        11920928955078125,
        59604644775390625,
        298023223876953125,
        1490116119384765625,
        7450580596923828125,
    ]

    power6 = [
        1,
        6,
        36,
        216,
        1296,
        7776,
        46656,
        279936,
        1679616,
        10077696,
        60466176,
        362797056,
        2176782336,
        13060694016,
        78364164096,
        470184984576,
        2821109907456,
        16926659444736,
        101559956668416,
        609359740010496,
        3656158440062976,
        21936950640377856,
        131621703842267136,
        789730223053602816,
        4738381338321616896,
    ]

    power7 = [
        1,
        7,
        49,
        343,
        2401,
        16807,
        117649,
        823543,
        5764801,
        40353607,
        282475249,
        1977326743,
        13841287201,
        96889010407,
        678223072849,
        4747561509943,
        33232930569601,
        232630513987207,
        1628413597910449,
        11398895185373143,
        79792266297612001,
        558545864083284007,
        3909821048582988049,
    ]

    power8 = [
        1,
        8,
        64,
        512,
        4096,
        32768,
        262144,
        2097152,
        16777216,
        134217728,
        1073741824,
        8589934592,
        68719476736,
        549755813888,
        4398046511104,
        35184372088832,
        281474976710656,
        2251799813685248,
        18014398509481984,
        144115188075855872,
        1152921504606846976,
    ]

    power9 = [
        1,
        9,
        81,
        729,
        6561,
        59049,
        531441,
        4782969,
        43046721,
        387420489,
        3486784401,
        31381059609,
        282429536481,
        2541865828329,
        22876792454961,
        205891132094649,
        1853020188851841,
        16677181699666569,
        150094635296999121,
        1350851717672992089,
    ]

    power10 = [
        1,
        10,
        100,
        1000,
        10000,
        100000,
        1000000,
        10000000,
        100000000,
        1000000000,
        10000000000,
        100000000000,
        1000000000000,
        10000000000000,
        100000000000000,
        1000000000000000,
        10000000000000000,
        100000000000000000,
        1000000000000000000,
    ]

    def overflow():
        print(f"fastpow: {m}^{n} larger than int64, using regular pow.")
        return pow(n, m, 1)

    if n == 0:
        return 0

    if n == 1:
        return 1

    if n == 2:
        return power2[m] if m < 63 else overflow()

    if n == 3:
        return power3[m] if m < 40 else overflow()

    if n == 4:
        return power4[m] if m < 32 else overflow()

    if n == 5:
        return power5[m] if m < 28 else overflow()

    if n == 6:
        return power6[m] if m < 25 else overflow()

    if n == 7:
        return power7[m] if m < 23 else overflow()

    if n == 8:
        return power8[m] if m < 21 else overflow()

    if n == 9:
        return power9[m] if m < 20 else overflow()

    if n == 10:
        return power10[m] if m < 19 else overflow()

    print(f"fastpow: {n} not in range [0, 10], using regular pow.")
    return pow(n, m, 1)
