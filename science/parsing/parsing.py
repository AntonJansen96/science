import sys
import os
import pandas

class Sanitize:
    """Class for sanitizing user input of various types."""

    def __init__(self, var, name: str = '', v: bool = True, exit: bool = True) -> None:
        """Initialize Sanitize object.

        Args:
            var (Any): variable under consideration.
            name (str, optional): Additional description of variable (for user message).
            v (bool, optional): verbose. Defaults to True.
            exit (bool, optional): exit python upon encountering an error. Defaults to False.
        """

        self.var = var
        self.__name = name
        self.__verbose = v
        self.__exit = exit
        self.__good = True

    def __error(self, msg: str) -> None:
        """Handles the error message program logic. Sets exitcode to 1 when executed.

        Args:
            msg (str): Information on the error.
        """

        pre = 'SanitizeError: '  # Pre-string for error messages.

        if self.__verbose:

            if self.__name:
                if type(self.var) == str:
                    print(f'{pre}[\'{self.var}\'] specified for [{self.__name}] {msg}.')
                else:
                    print(f'{pre}[{self.var}] specified for [{self.__name}] {msg}.')

            else:
                if type(self.var) == str:
                    print(f'{pre}[\'{self.var}\'] {msg}.')
                else:
                    print(f'{pre}[{self.var}] {msg}.')

        self.__good = False

    def __endbehavior(self):
        # if exit=True and we've had an error, exit python.
        if self.__exit and not self.__good:
            sys.exit(1)
        # If exit=False we return var, regardless of whether we've had an error.
        return self.var

    def num(self, Type=None, Range: list = [], signed: bool = False):
        """Sanitize numerical types (int, float, bool).

        Args:
            Type (Any, optional): (list of) acceptable type(s). Defaults to None.
            Range (list, optional): acceptable interval. Defaults to [].
            signed (bool, optional): Signed or unsigned. Defaults to False.

        Returns:
            int: exitcode.
        """

        # If Type was user-specified in any way, turn it into a list.
        if Type is not None and type(Type) != list:
            Type = [Type]

        # First of all, if you use this function var should be an int, float, bool.
        if type(self.var) not in [int, float, bool]:
            self.__error(f'not a number (should be in {[int, float, bool]})')

        else:
            # Second, if Type was user-specified in any way, check for match.
            if Type is not None and type(self.var) not in Type:
                self.__error(f'should be of type {Type} (found {type(self.var)})')

            # Check range.
            if len(Range) and (self.var < Range[0] or self.var > Range[1]):
                self.__error(f'outside of acceptable interval {Range}')

            # Check signed versus unsigned.
            if signed and self.var < 0:
                self.__error('cannot be negative')

        return self.__endbehavior()

    def string(self, Range: list = [], upper: bool = False, lower: bool = False, ws: bool = True):
        """sanitize strings (str).

        Args:
            Range (list, optional): acceptable string length. Defaults to [].
            upper (bool, optional): all uppercase. Defaults to False.
            lower (bool, optional): all lowercase. Defaults to False.
            ws (bool, optional): accept whitespace. Defaults to True.

        Returns:
            int: exitcode.
        """
        # Confirm whether var is a string.
        if type(self.var) != str:
            self.__error(f'should be of type {str}')

        # Check range.
        if len(Range) and (len(self.var) < Range[0] or len(self.var) > Range[1]):
            self.__error(f'outside of acceptable length {Range}')

        if upper and not self.var.isupper():
            self.__error('should be all uppercase')

        if lower and not self.var.islower():
            self.__error('should be all lowercase')

        if (not ws) and (self.var.count(' ') > 0):
            self.__error('cannot contain whitespace')

        return self.__endbehavior()

    def path(self, ext: str = '', out: bool = False, abs: bool = False):
        """Sanitize file paths (str).

        Args:
            ext (str, optional): acceptable extension(s). Defaults to ''.
            out (bool, optional): path meant for creation/output? Defaults to False.
            abs (bool, optional): absolute path required. Defaults to False.

        Returns:
            int: exitcode.
        """
        # Only strings can represent file paths.
        if type(self.var) != str:
            self.__error(f'should be of type {str}')
            return int(not self.__good)

        # File paths cannot be empty strings.
        if self.var == '':
            self.__error('cannot be empty')

        # File paths cannot be directories.
        elif os.path.isdir(self.var):
            self.__error('is a directory')

        # (input) file path should correspond to an exisiting file.
        elif not os.path.exists(self.var) and not out:
            self.__error('corresponding file does not exist')

        # (output) file DIRECTORY should already exist:
        if out and not os.path.isdir(os.path.split(os.path.abspath(self.var))[0]):
            self.__error('directory for output does not exist')

        # Check extension.
        if ext:
            tail = os.path.split(self.var)[1]

            if type(ext) != list:
                ext = [ext]

            if tail.count(' ') != 0:
                self.__error('cannot contain whitespace')

            if (tail.count('.') == 0) or (tail.count('.') > 1) or (tail[-1] == '.'):
                self.__error(f'ambiguous extension (should be {ext})')

            elif tail.index('.') == 0:
                self.__error('cannot be just an extension')

            elif tail[tail.index('.'):] not in ext:
                self.__error(f'should have extension {ext}')

        # Check absolute file path.
        if abs and not os.path.isabs(self.var):
            self.__error('should be an absolute file path')

        return self.__endbehavior()


def loadxvg(fname: str, col: list = [0, 1], dt: int = 1, b: int = 0):
    """Loads an .xvg file into a list of lists.
    May also be used to load float columns from files in general.

    Args:
        fname (str): file name.
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


def loadCol(fname: str, col: int = 1, header=None):
    """Loads a column in a file into a list. Automatically determines the type
    for the values in the list. Defaults to strings if multiple types detected
    in column.

    Args:
        fname (str): file name.
        col (int): the column number (starts at 1). Defaults to 1.
        header (int): row numbers to consider as header. header=None means no header, header=0 means 1 header line, header=1 means 2 header lines etc. Defaults to None.

    Returns:
        list: The column loaded into a list.
    """
    try:
        df = pandas.read_table(fname, header=header, delim_whitespace=True, na_filter=False)

    # Bug fix for when file / column is empty.
    except pandas.errors.EmptyDataError:
        print(f'loadCol: warning, \'{fname}\' is empty.')
        return []

    return list(df.iloc[:, col - 1])


def loadVal(fname: str, row: int, col: int, sep=None):
    """Retrieves a value from a file and returns it as a string.
    By default uses whitespace as a delimiter.

    Args:
        fname (str): file name.
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
