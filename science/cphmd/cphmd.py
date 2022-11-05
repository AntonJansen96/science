def protonation(xList, cutoff=0.8) -> float:
    """Returns the average protonation, i.e. the fraction of frames in which
    the lambda-coordinate is < 1 - cutoff.

    Args:
        xList (list): coordinate list.
        cutoff (float, optional): Defaults to 0.8.

    Returns:
        float: the average protonation.
    """

    lambda_proto   = 0
    lambda_deproto = 0

    for x in xList:

        if x > cutoff:
            lambda_deproto += 1

        if x < 1 - cutoff:
            lambda_proto += 1

    if lambda_proto + lambda_deproto == 0:
        fraction = 0
    else:
        fraction = float(lambda_proto) / (lambda_proto + lambda_deproto)

    return fraction


def deprotonation(xList, cutoff=0.8) -> float:
    """Returns the average deprotonation, i.e. the fraction of frames in which
    the lambda-coordinate is > cutoff.

    Args:
        xList (list): coordinate list.
        cutoff (float, optional): Defaults to 0.8.

    Returns:
        float: the average deprotonation.
    """

    lambda_proto   = 0
    lambda_deproto = 0

    for x in xList:

        if x > cutoff:
            lambda_deproto += 1

        if x < 1 - cutoff:
            lambda_proto += 1

    if lambda_proto + lambda_deproto == 0:
        fraction = 0
    else:
        fraction = float(lambda_deproto) / (lambda_proto + lambda_deproto)

    return fraction


def movingDeprotonation(xList, cutoff=0.8):
    """Returns a list containing the moving average deprotonation.

    Args:
        xList (list): coordinate list.
        cutoff (float, optional): Defaults to 0.8.

    Returns:
        list: list containing the moving average deprotonation.
    """
    Av = len(xList) * [0]
    lambda_proto = 1
    lambda_deproto = 0

    for idx in range(0, len(xList)):
        if xList[idx] > cutoff:
            lambda_deproto += 1
        elif xList[idx] < 1 - cutoff:
            lambda_proto += 1

        Av[idx] = float(lambda_deproto) / (lambda_proto + lambda_deproto)

    return Av
