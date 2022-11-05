import subprocess
import os
import time


class Stopwatch:
    """Stopwatch class. Helpful for profiling.
    """
    def __init__(self, description=''):
        """Initialize Stopwatch object. Count starts upon object initialization.

        Args:
            description (str, optional): Stopwatch description. Defaults to ''.
        """

        self.description = description
        self.starttime = time.time()

    def time(self):
        stoptime = time.time() - self.starttime
        if self.description == '':
            print('{:.3f}'.format(stoptime))
        else:
            print('{} took {:.3f} s'.format(self.description, stoptime))


def gromacs(command, stdin=[], terminal=True, basepath='/usr/local/gromacs_constantph', logFile='gromacs.log'):
    """Python function for handeling calls to GROMACS.

    Args:
        command (string): GROMACS command, e.g. 'make_ndx -f protein.pdb'.
        stdin (list, optional): list of input arguments. Defaults to [].
        terminal (bool, optional): print output to terminal. Defaults to True.
        basepath (str, optional): base path for GROMACS version to be used. Defaults to 'usr/local/gromacs_constantph'.
        logFile (str, optional): log file to use when terminal=False. Defaults to 'gromacs.log'.

    Returns:
        int: return code (0 if things were successful).
    """

    # If we don't pass envvars to subprocess (which happens by default) this works.
    path_to_gmx = os.path.normpath(basepath + '/' + 'bin/gmx')
    command = "{} {}".format(path_to_gmx, command)

    # Use the << EOF method to handle user input. Prepare string.
    if stdin:
        xstr = ' << EOF\n'
        for val in stdin:
            xstr += '{}\n'.format(val)
        command += xstr + 'EOF'

    if terminal:
        process = subprocess.run(command, shell=True, env={})
    else:
        with open(logFile, 'a+') as file:
            process = subprocess.run(command, shell=True, stdout=file, stderr=file, env={})

    if process.returncode != 0:
        if terminal:
            print("Failed to run \"{}\" (exitcode {}).".format(command, process.returncode))

        print("Failed to run \"{}\" (exitcode {}). Check your logfile ({})".format(command, process.returncode, logFile))

    return process.returncode


def inputOptionHandler(message, options):
    """Handles input options when prompting a user.

    Args:
        message (string): message for the user.
        options (list): list of options.

    Returns:
        int: number corresponding to selected option.
    """

    valids = []
    msgstring = "{}:".format(message)

    # Loop through the options list and create string for display.
    for idx in range(0, len(options)):
        msgstring += "\n{}. {}".format(idx, options[idx])
        valids.append(str(idx))

    while True:
        print(msgstring)
        val = input("Type a number: ")

        if val in valids:
            print()
            return int(val)

        print("{} is not a valid option, please try again:\n".format(val))


def exists(path):
    """Returns True if the path exists. Works for directories as well as files.

    Args:
        path (string): path to a directory or file.

    Returns:
        bool: does path exist?
    """

    return bool(os.path.exists(path))
