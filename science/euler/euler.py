# NUMBERS AND DIGITS ###########################################################


def numDigits(num: int) -> int:
    """Returns the number of digits of a number.

    Args:
        num (int): number.

    Returns:
        int: number of digits.
    """

    from math import log10

    return int(log10(num)) + 1


def sumDigits(num: int) -> int:
    """Returns the digit sum of a number.

    Args:
        num (int): number.

    Returns:
        int: sum of digits.
    """

    count = 0

    while num:

        count += num % 10
        num //= 10

    return count


def countDigits(num: int, digit: int) -> int:
    """Returns the number of of occurances of digit in num.

    Args:
        num (int): number.

    Returns:
        int: digit.
    """

    count = 0

    while num:

        if num % 10 == digit:
            count += 1

        num //= 10

    return count


def nthDigit(num: int, digit: int) -> int:
    """Returns the nth digit of a number (starting at the least significant digit).

    Args:
        num (int): number.
        digit (int): digit.

    Returns:
        int: count.
    """

    for _ in range(digit - 1):
        num //= 10

    return num % 10


def firstNdigits(num: int, N: int) -> int:
    """Returns the first N digits of a number.

    Args:
        num (int): number.

    Returns:
        int: count.
    """

    while num >= 10**N:
        num //= 10

    return num


# def firstNdigits(base: int, exp: int, N: int) -> int:
#     """Returns the first N digits of base^exp.

#     Args:
#         base (int): number.
#         exp (int): exponent.

#     Returns:
#         int: count.
#     """

#     from math import log10, floor

#     num = exp * log10(base)
#     num -= floor(num)
#     num = 10**num

#     return int(floor(10 ** (N - 1) * num))


def lastNdigits(num: int, N: int) -> int:
    """Returns the last N digits of a number.

    Args:
        num (int): number.
        N (int): number of digits.

    Returns:
        int: last N digits.
    """

    return num % 10**N


def reverseNumber(num: int) -> int:
    """Reverse a number. Example: 1234 -> 4321.

    Args:
        num (int): number.

    Returns:
        int: reversed number.
    """

    invnum = 0

    while num:

        invnum = invnum * 10 + (num % 10)
        num //= 10

    return invnum


# CHECK FOR CERTAIN PROPERTIES #################################################


def isSquare(num: int) -> bool:
    """Check if a number is a square.

    Args:
        num (int): number.

    Returns:
        bool: True if number is a square, False otherwise.
    """

    from math import isqrt

    h = num & 0xF  # Last hexadecimal "digit".

    if h > 9:
        return False

    if h not in [2, 3, 5, 6, 7, 8]:
        introot = isqrt(num)
        return introot * introot == num

    return False


def isCoprime(a: int, b: int) -> bool:
    """Check if a and b are coprime.

    Args:
        a (int): number.
        b (int): number.

    Returns:
        bool: True if a and b are coprime, False otherwise.
    """

    return gcd(a, b) == 1


def isPerfect(num: int) -> bool:
    """Check if a number is perfect (assumes num fits in int64).

    Args:
        num (int): number.

    Returns:
        bool: True if number is perfect, False otherwise.
    """

    if num in [
        6,
        28,
        496,
        8128,
        33550336,
        8589869056,
        137438691328,
        2305843008139952128,
    ]:
        return True

    return False


def isJuf(num: int) -> bool:
    """Check if a number is a Juf number.

    Args:
        num (int): number.

    Returns:
        bool: True if number is a Juf number, False otherwise.
    """

    return (num % 7 == 0) or ("7" in str(num)) or (num > 10 and isPalindrome(num))


def isPalindrome(num: int) -> bool:
    """Check if a number is a palindrome.

    Args:
        num (int): number.

    Returns:
        bool: True if number is a palindrome, False otherwise.
    """

    return num == reverseNumber(num)


def isPermutation(a: int, b: int) -> bool:
    """Check if a and b are permutations of each other.

    Args:
        a (int): number.
        b (int): number.

    Returns:
        bool: True if a and b are permutations of each other, False otherwise.
    """

    def __fingerprint(num):
        result = 0
        while num > 0:
            digit = num % 10
            num //= 10
            result += 1 << (5 * digit)
        return result

    return __fingerprint(a) == __fingerprint(b)


# MISCELLANEOUS FUNCTIONS ######################################################


def gcd(a: int, b: int) -> int:
    """Euclidian algorithm. Returns greatest common divisor of a and b.

    Args:
        a (int): number.
        b (int): number.

    Returns:
       int: greatest common denominator.
    """

    from math import gcd

    return gcd(a, b)
