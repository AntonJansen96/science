from math import isqrt


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
