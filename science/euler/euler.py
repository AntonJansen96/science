from typing import List, Generator, Tuple
from math import isqrt as _isqrt
from .fastmath import intlog10 as _intlog10, mulmod as _mulmod
from .combinatorics import genperms as _genperms

# CONTINUED FRACTIONS ##########################################################


class ContinuedFraction:
    """Class to convert between continued fractions and fractions/square roots."""

    def __init__(self) -> None:
        """Initializes ContinuedFraction object."""

    @staticmethod
    def fromfrac(p: int, q: int) -> list:
        """Converts a fraction p/q to a continued fraction representation.

        Args:
            p (int): numerator.
            q (int): denominator.

        Returns:
            list: continued fraction representation.
        """

        cf = []
        while q != 0:
            cf.append(p // q)
            p, q = q, p % q
        return cf

    @staticmethod
    def tofrac(cf: list, dec: bool = True) -> float | tuple:
        """Converts a continued fraction back to a fraction p/q.

        Args:
            cf (list): continued fraction representation.

        Returns:
            float | tuple: (p, q) or p/q.
        """

        p, q = 1, 0
        for coeff in reversed(cf):
            p, q = q, p
            p += coeff * q
        return p / q if dec else (p, q)

    @staticmethod
    def fromsqrt(N: int, lim: int = 0) -> list:
        """Converts a square root N to a continued fraction representation.

        Args:
            N (int): sqrt(N).
            lim (int, optional): Limit for the number of coefficients in the
            continued fraction. If set to zero, will compute until full period
            is found. Defaults to 0.

        Returns:
            list: continued fraction representation.
        """

        from math import isqrt

        if isqrt(N) ** 2 == N:
            return [isqrt(N)]

        cf = []
        m, d, a = 0, 1, isqrt(N)
        cf.append(a)

        if lim:
            for _ in range(1, lim):

                m = d * a - m
                d = (N - m * m) // d
                a = (isqrt(N) + m) // d
                cf.append(a)

            return cf

        else:
            states = set()
            while (m, d, a) not in states:
                states.add((m, d, a))

                m = d * a - m
                d = (N - m * m) // d
                a = (isqrt(N) + m) // d
                cf.append(a)

            return cf[:-1]

    @staticmethod
    def tosqrt(cf: list, iter: int = 2, dec: bool = True) -> float | tuple:
        """Converts a continued fraction back to a square root sqrt(N).
        Note that unlike the fraction case, the continued fraction representation
        of a square root is periodic. Therefore, the number N resulting from this
        function will be an approximation.

        Args:
            cf (list): continued fraction representation of the root.
            iter (int, optional): Number of times to cycle through the periodic
            part of the continued fraction of a square root when computing back
            the decimal. Defaults to 2.

        Returns:
            float | tuple: decimal or fractional approximation.
        """

        assert iter >= 0

        # We don't make cf #iter times as large but rather reloop the current cf.
        def helper(cf: list, iter: int):
            for _ in range(iter):
                for idx in range(len(cf) - 1, 0, -1):
                    yield cf[idx]
            yield cf[0]

        p, q = 1, 0
        for coeff in helper(cf, iter=iter):
            p, q = q, p
            p += coeff * q

        return (p * p) / (q * q) if dec else (p, q)


# NUMBERS AND DIGITS ###########################################################


def numDigits(num: int) -> int:
    """Returns the number of digits of a number.

    Args:
        num (int): number.

    Returns:
        int: number of digits.
    """

    return _intlog10(num) + 1


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

    h = num & 0xF  # Last hexadecimal "digit".

    if h > 9:
        return False

    if h not in [2, 3, 5, 6, 7, 8]:
        introot = _isqrt(num)
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

    # This (from LeetCode) is almost twice as fast as the code below.
    stringnum = str(num)
    return stringnum[::-1] == stringnum

    return num == reverseNumber(num)


def isPermutation(a: int, b: int) -> bool:
    """Check if a and b are permutations of each other.

    Args:
        a (int): number.
        b (int): number.

    Returns:
        bool: True if a and b are permutations of each other, False otherwise.
    """

    def fingerprint(num):
        result = 0
        while num > 0:
            digit = num % 10
            num //= 10
            result += 1 << (5 * digit)
        return result

    return fingerprint(a) == fingerprint(b)


def isAutomorphic(num: int) -> bool:
    """Check if a number is automorphic.
    Automorphic means that the square of the number ends with the number itself.

    Args:
        num (int): number.

    Returns:
        bool: True if number is automorphic, False otherwise.
    """

    m = 1
    temp = num
    while temp > 0:
        m *= 10
        temp //= 10

    return num == _mulmod(num, num, m)


def genLucky(lim: int) -> List[int]:
    """Returns a list of the lucky numbers up to and including lim.
    Note: because of the nature of lucky numbers, a separate isLucky function
    does not make sense. See also https://en.wikipedia.org/wiki/Lucky_number.

    Args:
        lim (int): upper limit.

    Returns:
        List[int]: The lucky numbers up to lim.
    """

    # This could probably be optimized but good for now.

    lucky = list(range(lim + 1))
    for v in lucky:
        # Get a sice of the lucky list starting from the index v & ~1
        # (which is always even) up to the end of the list,
        # stepping by v if v is not zero, otherwise by 2.
        del lucky[v & ~1 :: v or 2]

    return lucky


# GENERATORS ###################################################################


def fibonacci() -> Generator[int, None, None]:
    """Generator for Fibonacci numbers.

    Yields:
        int: Fibonacci number.
    """

    a, b = 0, 1
    while True:
        yield a
        a, b = b, a + b


def genPanDigital(a: int, b: int) -> Generator[int, None, None]:
    """Generates all numbers that are a to b pandigital.

    Args:
        a (int): start.
        b (int): end.

    Yields:
        num: pandigital number.
    """

    for permutation in _genperms([digit for digit in range(a, b + 1)], asint=True):
        if numDigits(permutation) == b - a + 1:
            yield permutation


def genPrimTriples(perimLim: int) -> Tuple[List[int], List[int], List[int]]:
    """Generates primitive Pythagorean triples with perimeter (a + b + c) < perimLim.

    Args:
        perimLim (int): perimeter limit.

    Returns:
        three lists A, B, C.
    """

    A = []
    B = []
    C = []

    for m in range(1, _isqrt(perimLim) + 1):
        for n in range(m + 1, _isqrt(perimLim) + 1):

            if not (m & 1 and n & 1) and isCoprime(m, n):
                a = n * n - m * m
                b = 2 * m * n
                c = n * n + m * m

                if a + b + c <= perimLim:
                    A.append(a)
                    B.append(b)
                    C.append(c)

    return A, B, C


def genAllTriples(perimLim: int) -> Tuple[List[int], List[int], List[int]]:
    """Generates all Pythagorean triples with perimeter (a + b + c) < perimLim.

    Args:
        perimLim (int): perimeter limit.

    Returns:
        three lists A, B, C.
    """

    A = []
    B = []
    C = []

    for m in range(1, _isqrt(perimLim) + 1):
        for n in range(m + 1, _isqrt(perimLim) + 1):

            if not (m & 1 and n & 1) and isCoprime(m, n):
                a = n * n - m * m
                b = 2 * m * n
                c = n * n + m * m

                k = 1
                while k * (a + b + c) <= perimLim:
                    print("hello")
                    A.append(k * a)
                    B.append(k * b)
                    C.append(k * c)
                    k += 1

    return A, B, C


# CONVERSIONS ##################################################################


def radix(value: int, radix: int) -> int:
    """Converts a number to a different radix.

    Args:
        value (int): number.
        radix (int): radix.

    Returns:
        int: number in different radix.
    """

    result = 0
    power = 1

    while value:
        result += power * (value % radix)
        value //= radix
        power *= 10

    return result


def dec2roman(dec: int) -> str:
    """Converts a decimal number to a Roman numeral.

    Args:
        dec (int): number.

    Returns:
        str: Roman numeral.
    """

    M = ["", "M", "MM", "MMM", "MMMM"]
    H = ["", "C", "CC", "CCC", "CD", "D", "DC", "DCC", "DCCC", "CM"]
    T = ["", "X", "XX", "XXX", "XL", "L", "LX", "LXX", "LXXX", "XC"]
    O = ["", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"]

    a = int(dec / 1000)
    b = int(dec % 1000 / 100)
    c = int(dec % 100 / 10)
    d = int(dec % 10)

    return M[a] + H[b] + T[c] + O[d]


def roman2dec(roman: str) -> int:
    """Converts a Roman numeral to a decimal number.

    Args:
        roman (str): Roman numeral.

    Returns:
        int: number.
    """

    chars = {"M": 1000, "D": 500, "C": 100, "L": 50, "X": 10, "V": 5, "I": 1}

    decimal = 0

    for ii in range(len(roman) - 1):

        if chars[roman[ii]] >= chars[roman[ii + 1]]:
            decimal += chars[roman[ii]]

        else:
            decimal -= chars[roman[ii]]

    decimal += chars[roman[-1]]

    return decimal


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


def leastCommonMultiple(a: int, b: int) -> int:
    """Returns the least common multiple of a and b.

    Args:
        a (int): number.
        b (int): number.

    Returns:
        int: least common multiple.
    """

    return a * b // gcd(a, b)
