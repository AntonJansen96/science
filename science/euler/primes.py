from math import isqrt
from .fastmath import powmod, mulmod


class Primes:
    """Class for doing stuff with prime numbers."""

    def __init__(self, setMax: int = 1000) -> None:
        """Initialize the class."""

        self.sieveArray = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
        self.nPrimes = 11  # Number of primes in sieveArray.
        self.max = 1000  # Maximum number.

        while self.max < setMax:
            self.__expand()

    def __expand(self) -> None:
        """Expand the sieveArray."""

        self.max *= 4  # Use 4 because root(4) = 2 is factor to increase primes with.
        next = self.sieveArray[-1]

        while next < isqrt(self.max) + 1:

            for prime in self.sieveArray:

                if next % prime == 0:
                    break

                if prime > isqrt(next):
                    self.sieveArray.append(next)
                    break

            next += 2

        self.nPrimes = len(self.sieveArray)

    def sieve(self, limit: int) -> list:
        """Generate all prime numbers up to a limit using the sieve of Eratosthenes.

        Args:
            limit (int): limit.

        Returns:
            list: prime numbers.
        """

        sieve = [2]
        next = 3

        while next <= limit:

            for prime in sieve:

                if next % prime == 0:
                    break

                if prime > isqrt(next):
                    sieve.append(next)
                    break

            next += 2

        return sieve

    def isPrime(self, num: int) -> bool:
        """Check if a number is prime.

        Args:
            num (int): number.

        Returns:
            bool: True if prime, False otherwise.
        """

        # Optimization: miller Rabin is faster than implementation below (on M1).
        return self.__millerRabin(num)

        if (num & 1) == 0:  # if num is even.
            return num == 2  # only even prime is 2.

        idx = 1
        while self.sieveArray[idx] <= isqrt(num):
            # If divisible by a prime, smaller than
            # the root of num, num is not prime.

            if num % self.sieveArray[idx] == 0:
                return False

            idx += 1

            if idx == self.nPrimes:
                self.__expand()

        return num > 1  # No number smaller than 2 is prime.

    def primePi(self, num: int) -> int:
        """Count the number of primes less than or equal to num.

        Args:
            num (int): number.

        Returns:
            int: number of primes.
        """

        # Will contain all primes up to sqrt(num).
        primes = []

        # If primes contains enough elements then run a fast binary search:
        if primes and primes[-1] > num:
            # Find smallest number larger than num:
            nextPrime = next((x for x in primes if x > num), None)
            return primes.index(nextPrime)

        v = isqrt(num)
        # About sqrt(num) * 12 bytes, for num = 10^12 => 12 MByte plus primes[].
        higher = [0] * (v + 2)
        lower = [0] * (v + 2)
        used = [False] * (v + 2)

        result = num - 1  # Assume all numbers are primes.
        # The remaining lines subtract composites until result contains #primes.

        for p in range(2, v + 1):  # set up lower and upper bound.
            lower[p] = p - 1
            higher[p] = num // p - 1

        for p in range(2, v + 1):
            if lower[p] == lower[p - 1]:  # Composite?
                continue

            if (
                not primes or p > primes[-1]
            ):  # Store prime numbers (if not already existing).
                primes.append(p)

            temp = lower[p - 1]  # Remove more composites.
            result -= higher[p] - temp

            pSquare = p * p
            end = min(v, num // pSquare)

            j = 1 + (p & 1)  # Alternate between 1 and 2.

            for i in range(p + j, end + 2, j):  # Adjust upper bound.
                if used[i]:
                    continue

                d = i * p
                if d <= v:
                    higher[i] -= higher[d] - temp
                else:
                    higher[i] -= lower[num // d] - temp

            for i in range(v, pSquare - 1, -1):  # Adjust lower bound.
                lower[i] -= lower[i // p] - temp

            for i in range(pSquare, end + 1, p * j):  # Cross off multiples.
                used[i] = True

        return result

    def __millerRabin(self, num: int) -> bool:
        """Check if a number is prime using the Miller-Rabin primality test.
        Guaranteed to be correct for numbers less than 2^64.

        Args:
            num (int): number.

        Returns:
            bool: True if prime, False otherwise.
        """

        # Only guaranteed to work for num < 2^64.
        assert num < 18446744073709551616

        bitmask_primes_2_to_31 = (  # Trivial cases.
            (1 << 2)
            | (1 << 3)
            | (1 << 5)
            | (1 << 7)
            | (1 << 11)
            | (1 << 13)
            | (1 << 17)
            | (1 << 19)
            | (1 << 23)
            | (1 << 29)
        )  # = 0x208A2Ac

        if num < 31:
            return bitmask_primes_2_to_31 & (1 << num) != 0

        if (  # divisible by a small prime.
            num % 2 == 0
            or num % 3 == 0
            or num % 5 == 0
            or num % 7 == 0
            or num % 11 == 0
            or num % 13 == 0
            or num % 17 == 0
        ):
            return False

        # We filtered all composite numbers < 17 * 19, all
        # others below 17 * 19 must be prime.
        if num < 17 * 19:
            return True

        # Test num against those numbers ("whitnesses").
        # Good bases can be found at http://miller-rabin.appspot.com.
        STOP = 0
        TestAgainst1 = [377687, STOP]
        TestAgainst2 = [31, 73, STOP]
        TestAgainst3 = [2, 7, 61, STOP]
        # First three sequences are good up to 2^32.
        TestAgainst4 = [2, 13, 1662803, STOP]
        TestAgainst7 = [2, 325, 9375, 28178, 450775, 9780504, 1795265022, STOP]

        # Good up to 2^64.
        testAgainst = TestAgainst7
        if num < 5329:  # Use less tests if feasible...
            testAgainst = TestAgainst1
        elif num < 9080191:
            testAgainst = TestAgainst2
        elif num < 4759123141:
            testAgainst = TestAgainst3
        elif num < 1122004669633:
            testAgainst = TestAgainst4

        d = num - 1  # Find num - 1 = d * 2^j.
        shift = 0
        d >>= 1
        while d & 1 == 0:
            shift += 1
            d >>= 1

        testAgainstIndex = 0
        while testAgainst[testAgainstIndex] != STOP:  # Test num against all bases.
            x = powmod(testAgainst[testAgainstIndex], d, num)
            testAgainstIndex += 1

            if x == 1 or x == num - 1:  # Is test^d % num == 1 || -1 ?
                continue

            maybePrime = False  # Now either prime or strong pseudo-prime.
            for r in range(shift):  # Check test^(d*2^r) for 0 <= r < shift:
                x = mulmod(x, x, num)  # x^2 % num (initial x was test^d)

                if x == 1:  # If x % num == 1 => not prime.
                    return False

                if (
                    x == num - 1
                ):  # If x % num == -1 => not prime or an even stronger pseudo-prime.
                    maybePrime = True
                    break  # Next iteration.

            if not maybePrime:  # Not prime.
                return False

        return True  # Prime.
