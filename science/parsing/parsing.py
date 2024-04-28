import os
import sys
from datetime import datetime as _datetime
from typing import List, Union


class Sanitize:
    """Class for sanitizing user input of various types."""

    def __init__(self, var, name: str = "", v: bool = True, exit: bool = True) -> None:
        """Initialize Sanitize object.

        Args:
            var (Any): variable under consideration.
            name (str, optional): additional description of variable (for user message).
            v (bool, optional): verbose. Defaults to True.
            exit (bool, optional): exit python upon encountering an error. Defaults to False.
        """

        self.var = var
        self._name = name
        self._verbose = v
        self._exit = exit
        self._good = True

    def _error(self, msg: str) -> None:
        """Handles the error message program logic. Sets exitcode to 1 when executed.

        Args:
            msg (str): Information on the error.
        """

        pre = "SanitizeError: "  # Pre-string for error messages.

        if self._verbose:

            if self._name:
                if type(self.var) == str:
                    print(f"{pre}['{self.var}'] specified for [{self._name}] {msg}.")
                else:
                    print(f"{pre}[{self.var}] specified for [{self._name}] {msg}.")

            else:
                if type(self.var) == str:
                    print(f"{pre}['{self.var}'] {msg}.")
                else:
                    print(f"{pre}[{self.var}] {msg}.")

        self._good = False

    def _endbehavior(self):
        # if exit=True and we've had an error, exit python.
        if self._exit and not self._good:
            sys.exit(1)
        # if exit=False and we've had an error, return None.
        elif not self._good:
            return None
        # If exit=False and we haven't had an error we're good.
        else:
            return self.var

    def num(self, Type=None, Range: List[Union[int, float]] = [], signed: bool = False):
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
            self._error(f"not a number (should be in {[int, float, bool]})")

        else:
            # Second, if Type was user-specified in any way, check for match.
            if Type is not None and type(self.var) not in Type:
                self._error(f"should be of type {Type} (found {type(self.var)})")

            # Check range.
            if len(Range) and (self.var < Range[0] or self.var > Range[1]):
                self._error(f"outside of acceptable interval {Range}")

            # Check signed versus unsigned.
            if signed and self.var < 0:
                self._error("cannot be negative")

        return self._endbehavior()

    def string(
        self,
        Range: list = [int],
        upper: bool = False,
        lower: bool = False,
        ws: bool = True,
    ):
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
            self._error(f"should be of type {str}")

        # Check range.
        if len(Range) and (len(self.var) < Range[0] or len(self.var) > Range[1]):
            self._error(f"outside of acceptable length {Range}")

        if upper and not self.var.isupper():
            self._error("should be all uppercase")

        if lower and not self.var.islower():
            self._error("should be all lowercase")

        if (not ws) and (self.var.count(" ") > 0):
            self._error("cannot contain whitespace")

        return self._endbehavior()

    def path(self, ext: str = "", out: bool = False, abs: bool = False):
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
            self._error(f"should be of type {str}")
            return int(not self._good)

        # File paths cannot be empty strings.
        if self.var == "":
            self._error("cannot be empty")

        # File paths cannot be directories.
        elif os.path.isdir(self.var):
            self._error("is a directory")

        # (input) file path should correspond to an exisiting file.
        elif not os.path.exists(self.var) and not out:
            self._error("corresponding file does not exist")

        # (output) file DIRECTORY should already exist:
        if out and not os.path.isdir(os.path.split(os.path.abspath(self.var))[0]):
            self._error("directory for output does not exist")

        # Check extension.
        if ext:
            tail = os.path.split(self.var)[1]

            if type(ext) != list:
                ext = [ext]

            # Add "." to extensions if not present.
            for idx in range(len(ext)):
                if ext[idx][0] != ".":
                    ext[idx] = "." + ext[idx]

            # Comment this because file names can contain whitespace in linux.
            # if tail.count(" ") != 0:
            #     self._error("cannot contain whitespace")

            if "." not in tail or (tail.index(".") == 0 and tail.count(".") == 1):
                self._error(f"does not have a file extension (should be {ext})")

            elif "." + tail.split(".")[-1] not in ext:
                self._error(f"should have extension {ext}")

        # Check absolute file path.
        if abs and not os.path.isabs(self.var):
            self._error("should be an absolute file path")

        return self._endbehavior()


class User:
    """Provides logging and formats user updates, warnings, and errors."""

    def __init__(
        self,
        baseName: str = "",
        verbosity: int = 2,
        logFileName: str = "",
        maxLineLength: int = 90,
    ):
        """Initialize User object.

        Args:
            baseName (str): base message string. Defaults to ''.
            verbosity (int): verbosity value. 0 = supress all, 1 = only warnings and errors, 2 = regular, 3 = verbose. Defaults to 2.
            logFileName (str): name of log file (not specified means no logging).
            maxLineLength (int): maximum length of sentence before a newline (does not consider length of preMessage etc). Defaults to 90.
        """

        self._baseName = baseName
        self._verbosity = verbosity
        self._logFileName = logFileName
        self._maxLineLength = maxLineLength

    def _log(self, message: str):
        if self._logFileName:
            with open(self._logFileName, "a+") as logfile:
                time = _datetime.now().strftime("%Y/%m/%d|%H:%M:%S")
                print(time + " " + message, file=logfile)

    def _output(self, message: str):
        """Handles output to terminal and (optionally) logging.

        Args:
            message (str): message to be printed/logged.
        """

        self._log(message)
        print(message)

    def doInput(self) -> str:
        """Handles input from terminal and (optionally) logging."""

        string = input(self._baseName)
        self._log(string)
        return string

    def _base(self, message: str, preMessage: str = ""):
        """Base method for handeling user messages.

        Args:
            message (str): message.
            preMessage (str, optional): pre-message. Defaults to ''.
        """

        # If we input a dictionary of list as a message, print and return.
        if isinstance(message, list) or isinstance(message, dict):
            self._output(message)
            return

        firstLineWritten = False
        charCount = 0
        currentLine = ""

        for word in message.split(" "):

            charCount += len(word) + 1

            if charCount < self._maxLineLength:
                currentLine += word + " "
            else:
                self._output(f"{self._baseName}{preMessage}{currentLine}")
                firstLineWritten = True
                charCount = len(word) + 1
                currentLine = word + " "

            # This we add so that not every subsequent line has preMessage,
            # but it is aligned nonetheless.
            if firstLineWritten and preMessage != "":
                preMessage = " ".ljust(len(preMessage))

        self._output(f"{self._baseName}{preMessage}{currentLine.lstrip()}")  # Flush.

    def verbose(self, message: str):
        """Print message when verbosity is high.

        Args:
            message (str): message.
        """

        if self._verbosity > 2:
            self._base(message)

    def update(self, message: str):
        """Print default message.

        Args:
            message (str): message.
        """

        if self._verbosity > 1:
            self._base(message)

    def warning(self, message: str):
        """Print warning message.

        Args:
            message (str): message.
        """

        if self._verbosity > 0:
            self._output(self._baseName)
            self._base(message, preMessage="WARNING - ")
            self._output(self._baseName)

    def error(self, message: str):
        """Print error message and sys.exit() the program.

        Args:
            message (str): message.
        """

        if self._verbosity > 0:
            self._output(self._baseName)
            self._base(message, preMessage="ERROR - ")
            self._output(self._baseName)

        sys.exit(1)

    def inputOptionHandler(self, message: str, options: List[str]) -> int:
        """Handles user input when options are required.

        Args:
            message (str): message.
            options (list): list of strings describing the options. Starts counting from 0.

        Returns:
            int: number corresponding to the option selected by the user.
        """

        valids = []
        msgstring = f"{self._baseName}{message}"

        # Loop through the options list and create string for display
        for idx in range(0, len(options)):
            msgstring += f"\n{self._baseName}{idx}. {options[idx]}"
            valids.append(str(idx))

        while True:
            self._output(msgstring)
            val = input(f"{self._baseName}Type a number: ")

            if val in valids:
                self._output(f"{self._baseName}selected {val}")
                self._output("")
                return int(val)

            self._output(
                f"{self._baseName}{val} is not a valid option, please try again:\n"
            )


def loadxvg(fname: str, col: List[int] = [0, 1], dt: int = 1, b: int = 0):
    """Loads an .xvg file into a list of lists.
    May also be used to load float columns from files in general.

    Args:
        fname (str): file name.
        col (list, optional): columns to load. Defaults to [0, 1].
        dt (int, optional): step size. Defaults to 1.
        b (int, optional): starting point. Defaults to 0.

    Returns:
        list of lists : contains the columns that were loaded.
    """

    count = -1
    data = [[] for _ in range(len(col))]
    for stringLine in open(fname).read().splitlines():
        if stringLine[0] in ["@", "#", "&"]:
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
        header (int): row numbers to consider as header. header=None means no header,
        header=0 means 1 header line, header=1 means 2 header lines etc. Defaults to None.

    Returns:
        list: The column loaded into a list.
    """

    # Lazy/deferred import.
    import pandas

    try:
        df = pandas.read_table(
            fname, header=header, delim_whitespace=True, na_filter=False
        )

    # Bug fix for when file / column is empty.
    except pandas.errors.EmptyDataError:
        print(f"loadCol: warning, '{fname}' is empty.")
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


def pickleDump(var, file: str, protocol: int = 5) -> None:
    """Dumps a variable to a file name file using protocol 5 (most efficient).

    Args:
        var (any): variable name.
        file (str): name of file to dump to.
    """

    # Lazy/deferred import.
    from pickle import dump

    dump(var, open(file, "wb"), protocol)


def pickleLoad(file: str):
    """Loads a pickled variable from a file.

    Args:
        file (str): name of file to load from.

    Returns:
        any: variable that was loaded.
    """

    # Lazy/deferred import.
    from pickle import load

    return load(open(file, "rb"))
