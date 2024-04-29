from typing import Generator


def cycleCount(generator: Generator[int, None, None], maxCycle: int) -> int:
    """Detect cycle length in a sequence by shifting and overlapping the lists.

    Args:
        generator (Generator[int, None, None]): sequence generator.
        maxCycle (int): maximum cycle length to search for.

    Returns:
        int: cycle length. Returns 0 if no cycle was found.
    """

    lst = []
    for _ in range(2 * maxCycle):
        lst.append(next(generator))
        length = len(lst)
        for i in range(1, length // 2 + 1):
            if lst[:i] * (length // i) == lst[: i * (length // i)]:
                return i

    return 0


def pisanoPeriod(n: int) -> int:
    """Returns the Pisano period for the Fibonacci sequence modulo n.

    Args:
        n (int): modulo n.

    Returns:
        int: Pisano period.
    """

    a = 0
    b = 1
    c = a + b
    for i in range(0, n * n):
        c = (a + b) % n
        a, b = b, c
        if a == 0 and b == 1:
            return i + 1
