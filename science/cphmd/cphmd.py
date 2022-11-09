import MDAnalysis


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


def getLambdaFileIndices(structure, resid):
    """Returns an array containing the lambda-file indices for the specified resid.
    Only takes into account ASPT, GLUT, HSPT.

    Args:
        structure (string): pdb file name.
        resid (int): residue id.

    Returns:
       list: List of lambda indices.
    """
    u                  = MDAnalysis.Universe(structure)
    numChains          = len(u.segments) - 1
    segmentAatoms      = u.segments[0].atoms
    titratableAtoms    = segmentAatoms.select_atoms('resname ASPT GLUT HSPT')
    titratableResnames = list(titratableAtoms.residues.resnames)
    titratableResids   = list(titratableAtoms.residues.resids)
    targetidx          = titratableResids.index(resid)

    numASPTGLUT        = len(segmentAatoms.select_atoms('resname ASPT GLUT').residues)
    numHSPT            = len(segmentAatoms.select_atoms('resname HSPT').residues)
    factor             = numASPTGLUT + 3 * numHSPT

    count = 1
    for idx in range(0, len(titratableResnames)):

        if idx == targetidx:
            array = []
            for ii in range(0, numChains):
                array.append(count + ii * factor)
            return array

        if titratableResnames[idx] in ['ASPT', 'GLUT']:
            count += 1

        elif titratableResnames[idx] == 'HSPT':
            count += 3
