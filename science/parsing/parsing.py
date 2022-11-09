import pandas


def loadxvg(fname, col=[0, 1], dt=1, b=0):
    """Loads an .xvg file into a list of lists.
    May also be used to load float columns from files in general.

    Args:
        fname (string): file name.
        col (list, optional): Columns to load. Defaults to [0, 1].
        dt (int, optional): Step size. Defaults to 1.
        b (int, optional): Starting point. Defaults to 0.

    Returns:
        list of lists : contains the columns that were loaded.
    """

    count = -1
    data = [[] for _ in range(len(col))]
    for stringLine in open(fname).read().splitlines():
        if stringLine[0] in ['@', '#', '&']:
            continue
        # THIS IS FOR THE dt PART.
        count += 1
        if count % dt != 0:
            continue

        listLine = stringLine.split()
        # AND THIS IS FOR THE b PART.
        if b != 0 and float(listLine[col[0]]) < b:
            continue

        for idx in col:
            data[idx].append(float(listLine[col[idx]]))
    return data


def loadCol(fname, col=1, header=None):
    """Loads a column in a file into a list. Automatically determines the type
    for the values in the list. Defaults to strings if multiple types detected
    in column.

    Args:
        fname (string): file name.
        col (int): the column number (starts at 1). Defaults to 1.
        header (int): row numbers to consider as header. header=None means no header, header=0 means 1 header line, header=1 means 2 header lines etc. Defaults to None.

    Returns:
        list: The column loaded into a list.
    """

    df = pandas.read_table(fname, header=header, delim_whitespace=True, na_filter=False)

    return list(df.iloc[:, col - 1])


def loadVal(fname, row, col, sep=None):
    """Retrieves a value from a file and returns it as a string.
    By default uses whitespace as a delimiter.

    Args:
        fname (string): file name.
        row (int): row number (starts at 1).
        col (int): column number (starts at 1).
        sep (str, optional): separator (white space by default). Defaults to None.

    Returns:
        string: the value.
    """

    try:
        for x, y in enumerate(open(fname)):
            if x == row - 1:
                return y.split(sep)[col - 1]

    except IndexError:
        return None
