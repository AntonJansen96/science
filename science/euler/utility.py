def num2list(num: int) -> list:
    """Converts a number to a vector of digits.

    Args:
        num (int): number.

    Returns:
        list: vector of digits.
    """

    numbers = []

    while num:
        numbers.append(num % 10)
        num //= 10

    return numbers[::-1]


def list2num(array: list) -> int:
    """Convert a list to a number.

    Args:
        array (list): input.

    Returns:
        int: number.
    """

    number = 0
    for digit in array:
        number = number * 10 + digit

    return number
