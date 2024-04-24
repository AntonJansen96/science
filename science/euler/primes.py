from typing import List, Set, Union
from math import isqrt as _isqrt
from .fastmath import powmod as _powmod, mulmod as _mulmod


class Primes:
    """Class for doing stuff with prime numbers and factoring."""

    def __init__(self, setMax: int = 1000, debug: bool = False) -> None:
        """Initialize the Primes class.

        Args:
            setMax (int, optional): upper limit of numbers to check. Defaults to 1000.
            Note: will expand the sieveArray automatically if setMax is exceeded.
        """

        self._debug = debug
        self._sieveArray = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
        self._nPrimes = 11  # Number of primes in sieveArray.
        self._max = 1000  # Maximum number.

        # Handle logging.
        if self._debug:
            import logging # Lazy import.

            logging.basicConfig(level=logging.INFO)
            self._logger = logging.getLogger(__name__)

        # Call expand if request max prime size is larger than current max.
        while self._max < setMax:
            self._log("_expand() called from __init__.")
            self._expand()

    def _log(self, message: str) -> None:
        """Print a logging message using Python's logging module."""

        if self._debug:
            self._logger.info(message)

    def _expand(self) -> None:
        """Expand the sieveArray."""

        self._max *= 4  # Use 4 because root(4) = 2 is factor to increase primes with.
        self._log(f"_max is now {self._max}.")

        next = self._sieveArray[-1]

        while next < _isqrt(self._max) + 1:

            for prime in self._sieveArray:

                if next % prime == 0:
                    break

                if prime > _isqrt(next):
                    self._sieveArray.append(next)
                    break

            next += 2

        self._log(
            "Length of _sieveArray is now {}, largest prime is {}.".format(
                len(self._sieveArray), self._sieveArray[-1]
            )
        )
        self._nPrimes = len(self._sieveArray)

    @staticmethod
    def sieve(limit: int) -> List[int]:
        """Generate all prime numbers up to a limit using the sieve of Eratosthenes.
        Note: does not use self._sieveArray but starts from the beginning.

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

                if prime > _isqrt(next):
                    sieve.append(next)
                    break

            next += 2

        return sieve

    @staticmethod
    def isprime(num: int) -> bool:
        """Check if a number is prime.

        Args:
            num (int): number.

        Returns:
            bool: True if prime, False otherwise.
        """

        # Optimization: miller Rabin is faster than implementation below (on M1).
        # If we want to use the code below, isprime can no-longer be a @staticmethod.
        return Primes._millerrabin(num)

        if (num & 1) == 0:  # if num is even.
            return num == 2  # only even prime is 2.

        idx = 1
        while self._sieveArray[idx] <= _isqrt(num):
            # If divisible by a prime, smaller than
            # the root of num, num is not prime.

            if num % self._sieveArray[idx] == 0:
                return False

            idx += 1

            if idx == self._nPrimes:
                self._expand()

        return num > 1  # No number smaller than 2 is prime.

    @staticmethod
    def primepi(num: int) -> int:
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

        v = _isqrt(num)
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

    def primefactors(self, num: int, multiplicity: bool = True) -> List[int]:
        """Find the prime factors of a number.

        Args:
            num (int): number.
            multiplicity (bool, optional): include multiplicity. Defaults to True.

        Returns:
            list: prime factors.
        """

        # Speedup.
        if self.isprime(num):
            return [num]

        factors = []
        idx = 1  # Start at second prime in primes (3) because we already checked 2.

        add = True  # Handle multiplicity.
        while (num >> 1) << 1 == num:  # while num is divisible by 2
            if add:
                factors.append(2)
                if not multiplicity:
                    add = False
            num //= 2  # num = num // 2

        # Only have to check up until root of num.
        while self._sieveArray[idx] <= _isqrt(num):
            prime = self._sieveArray[idx]

            add = True
            while num % prime == 0:  # While num is divisible by prime ...
                if add:
                    factors.append(prime)  # ... add that prime to factors,
                    if not multiplicity:
                        add = False
                num //= prime  # and divide.

            idx += 1

            # If we have reached the end of the current primes list,
            # call expand() (i.e. double the list).
            if idx == self._nPrimes:
                self._expand()

        # If at the end we're left with a number larger than two,
        # it must be a prime so add to factors.
        if num > 2:
            factors.append(num)

        return factors

    def factors(self, num: int, proper: bool = False) -> Union[Set[int], List[int]]:
        """Find all factors of a number.
        Note: factors are an unsorted set while proper factors are a sorted list.

        Args:
            num (int): number.
            proper (bool): exclude num itself.

        Returns:
            list: factors.
        """

        assert num > 0, "num <= 0"

        # Do not first check if a number is prime, because
        # we're already doing that in self.primefactors().
        primefactors = self.primefactors(num)

        # If prime, return immediately.
        if len(primefactors) == 1:
            return [1] if proper else {1, num}

        # Iterative approach with bit manipulation technique.
        factors = {1}
        for prime in primefactors:
            factors.update({factor * prime for factor in factors})

        # We need to sort for proper factors because
        # the number itself is not always at the end.
        return sorted(factors)[:-1] if proper else factors

    # def factors(self, num: int, proper: bool = False) -> List[int]:
    #     """Find all factors of a number. Note: list is not sorted.

    #     Args:
    #         num (int): number.
    #         proper (bool): exclude num itself.

    #     Returns:
    #         list: factors.
    #     """

    #     assert num > 0, "num <= 0"

    #     factors = []

    #     # This is depth-first search is slow, maybe want to find a faster method.
    #     def findfactors(num, primefactors, idx, factor):

    #         if idx == len(primefactors):

    #             factors.append(factor)
    #             return

    #         while True:

    #             findfactors(num, primefactors, idx + 1, factor)
    #             factor *= primefactors[idx]

    #             if num % factor != 0:
    #                 break

    #     findfactors(num, self.primefactors(num, multiplicity=False), 0, 1)

    #     return factors[:-1] if proper else factors

    def largestfactor(self, num: int) -> int:
        """Find the largest nontrivial factor of a number.
        Returns 1 if a number is prime.

        Args:
            num (int): number.

        Returns:
            int: largest factor.
        """

        return max(self.factors(num, proper=True))

    def totient(self, num: int) -> int:
        """Calculate Euler's totient function for a number.

        Args:
            num (int): number.

        Returns:
            int: totient.
        """

        totient = num

        if (num >> 1) << 1 == num:  # num % 2 == 0
            totient -= totient // 2

            while (num >> 1) << 1 == num:  # while num % 2 = 0, num /= 2
                num >>= 1

        idx = 1
        while self._sieveArray[idx] <= _isqrt(num):
            prime = self._sieveArray[idx]

            if num % prime == 0:
                totient -= totient // prime

                while num % prime == 0:
                    num //= prime

            idx += 1

            if idx == self._nPrimes:  # If we have reached the end of the
                self._expand()  # current primes list, call expand().

        if num > 2:
            totient -= totient // num

        return totient

    def isamicable(self, num: int) -> bool:
        """Check if a number is amicable.

        Args:
            num (int): number.

        Returns:
            bool: True if amicable, False otherwise.
        """

        assert num > 1, "num <= 1"

        x = sum(self.factors(num, proper=True))
        y = sum(self.factors(x, proper=True))

        return x != y and num == y

    def isabundant(self, num: int) -> bool:
        """Check if a number is abundant.

        Args:
            num (int): number.

        Returns:
            bool: True if abundant, False otherwise.
        """

        # Abundant if sum of factors is greater than 2 * num.
        return sum(self.factors(num)) > 2 * num

    @staticmethod
    def _millerrabin(num: int) -> bool:
        """Check if a number is prime using the Miller-Rabin primality test.
        Guaranteed to be correct for numbers less than 2^64.

        Args:
            num (int): number.

        Returns:
            bool: True if prime, False otherwise.
        """

        # Only guaranteed to work for num < 2^64.
        assert num < 18446744073709551616, "num > 2^64"

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
            x = _powmod(testAgainst[testAgainstIndex], d, num)
            testAgainstIndex += 1

            if x == 1 or x == num - 1:  # Is test^d % num == 1 || -1 ?
                continue

            maybePrime = False  # Now either prime or strong pseudo-prime.
            for r in range(shift):  # Check test^(d*2^r) for 0 <= r < shift:
                x = _mulmod(x, x, num)  # x^2 % num (initial x was test^d)

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
