import time
import os
import subprocess
import MDAnalysis
import scipy.stats as stats
import copy


class Stopwatch:
    """Stopwatch class. Helpful for profiling."""

    def __init__(self, description=''):
        """Initialize Stopwatch object. Count starts upon object initialization.

        Args:
            description (str, optional): Stopwatch description. Defaults to ''.
        """

        self.description = description
        self.starttime = time.time()

    def time(self):
        """Print elapsed time (in seconds).
        """
        stoptime = time.time() - self.starttime
        if self.description == '':
            print('{:.3f}'.format(stoptime))
        else:
            print('{} took {:.3f} s'.format(self.description, stoptime))


def gromacs(command: str, stdin: list = [], terminal: bool = True, basepath: str = '/usr/local/gromacs_constantph', logFile: str = 'gromacs.log'):
    """Python function for handeling calls to GROMACS.

    Args:
        command (str): GROMACS command, e.g. 'make_ndx -f protein.pdb'.
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
        else:
            print("Failed to run \"{}\" (exitcode {}). Check your logfile ({})".format(command, process.returncode, logFile))

    return process.returncode


def createIndexFile(inputFile: str, outputFile: str, groups: list) -> str:
    """Creates an index (.ndx) file for GROMACS for a list of specified groups.
    These groups are strings that use the MDAnalysis selection syntax.

    Args:
        inputFile (str): input structure file (pdb / gro) name.
        outputFile (str): output file (ndx) name.
        groups (list): list of strings. These strings use the MDAnalysis selection syntax.

    Returns:
        str: outputFile.
    """

    u = MDAnalysis.Universe(inputFile)

    with open(outputFile, 'w') as file:

        print('Creating index file \'{}\' containing {} index groups:'.format(outputFile, len(groups)))

        groupidx = 0

        for group in groups:

            lineCount = 0

            # Write header. We replace any white space with underscores.
            header = group.replace(' ', '_')
            file.write('[ {} ]\n'.format(header))

            indices = list(u.select_atoms(group).atoms.indices)

            for idx in indices:
                # It's idx + 1 because we want atom numbers, not indices.
                file.write('{:<7d}'.format(idx + 1))

                # This is to prevent having too many atom numbers on one line.
                lineCount += 1
                if lineCount == 13:
                    file.write('\n')
                    lineCount = 0

            file.write('\n\n')

            print('group {} \'{}\' contains {} atoms'.format(groupidx, header, len(indices)))
            groupidx += 1

    return outputFile


def inputOptionHandler(message: str, options: list):
    """Handles input options when prompting a user.

    Args:
        message (str): message for the user.
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


def triplet2letter(triplet: str):
    """Converts a resname triplet (e.g. GLU) to a single letter (e.g. E).
    Note: also does NA->Na+ and CL->Cl-.

    Args:
        triplet (str): triplet.

    Returns:
        char: single letter.
    """

    dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
            'ASPT': 'D', 'GLUT': 'E', 'HSPT': 'H', 'NA': 'Na+', 'CL': 'Cl-'}

    return dict[triplet]


def ttestPass(sample1: list, sample2: list, alpha: float = 0.05):
    """Returns True if the means of sample1 and sample2 differ SIGNIFICANTLY.
    That is, with a confidence interval of 1 - alpha %. Uses Welch's t-test.

    Args:
        sample1 (list): Some sample.
        sample2 (list): Sample to compare to.
        alpha (float, optional): Significance of the test. Defaults to 0.05.

    Returns:
        bool: Whether or not the two sample means differ signifcantly.
    """

    pvalue = stats.ttest_ind(sample1, sample2, equal_var=False)[1]

    return bool(pvalue < alpha)


def makeSuperDict(keyLists: list):
    """Initializes a (recursively) nested dictionary of dictionaries based on a list
    of lists of keys. Helpful for organizing large amounts of data.

    Example:
        keyLists = [['A', 'B'], [1, 2], 0.33] -> {'A': {1: 0.33, 2: 0.33}, 'B': {1: 0.33, 2: 0.33}}.

    Args:
        keyLists (list): the list of lists of nested keys. The last element corresponds
        to the value of the innermost key-value pair.

    Returns:
        dict: desired nested dictionary or 'superDict' structure.
    """

    # Assertions
    assert isinstance(keyLists[0], list)
    assert len(keyLists) > 1

    # Calculate size of the superDict and provide user update.
    size = 1
    for idx in range(0, len(keyLists) - 1):
        size *= len(keyLists[idx])
    print('Created superDict of size {} (recursion depth {}).'.format(size, len(keyLists) - 1))

    # EXPLANATION THROUGH THE FOLLOWING EXAMPLE:

    # A = []
    # B = {}
    # for key in metrics:
    #     B[key] = copy.deepcopy(A)
    # C = {}
    # for key in chains:
    #     C[key] = copy.deepcopy(B)
    # D = {}
    # for key in reps:
    #     D[key] = copy.deepcopy(C)
    # E = {}
    # for key in sims:
    #     E[key] = copy.deepcopy(D)

    # List  = [sims, reps, chains, metrics, []]
    # array = [{}, {}, {}, {}, List[4]]

    # for key in List[3]:
    #     array[3][key] = copy.deepcopy(array[4])

    # for key in List[2]:
    #     array[2][key] = copy.deepcopy(array[3])

    # for key in List[1]:
    #     array[1][key] = copy.deepcopy(array[2])

    # for key in List[0]:
    #     array[0][key] = copy.deepcopy(array[1])

    array = [copy.deepcopy({}) for _ in range(0, len(keyLists) - 1)]
    array += [copy.deepcopy(keyLists[-1])]

    for idx in range(1, len(keyLists))[::-1]:
        for key in keyLists[idx - 1]:
            array[idx - 1][key] = copy.deepcopy(array[idx])

    return array[0]


def genRestraints(pdb: str, fname: str, atomSelection: str):
    """Write a position restraint file for specified atomSelection.
    Note: funct = 1, fcx = fcy = fcz = 1000 (hardcoded).

    Args:
        pdb (str): structure file.
        fname (str): output file (.itp) name.
        atomSelection (str): MDAnalysis style selection string.
    """

    u = MDAnalysis.Universe(pdb)
    sel = u.select_atoms(atomSelection).atoms.indices

    with open(fname, 'w+') as file:
        # Write header
        file.write('; Position restraints for {}\n'.format(pdb))
        file.write('; Selection is \'{}\'\n\n'.format(atomSelection))
        file.write('[ position_restraints ]\n')
        file.write(';  i funct       fcx        fcy        fcz\n')

        # Write position restraints
        for idx in list(sel):
            file.write('{}  1  1000  1000  1000\n'.format(idx + 1))

    # User update and reminder
    print('Generated position restraints file \'{}\''.format(fname))
    print('Don\'t forget to add the following to (the bottom of) your topology:\n')
    print('#ifdef POSRES_NAME\n#include \"{}\"'.format(fname))
    print('#endif\n')
    print('And the following to your .mdp file:\n')
    print('define = -DPOSRES_NAME\n')


def backup(name: str, verbose: bool = True):
    """Create a GROMACS-style (e.g. '#MD.log.1#') backup of file name.

    Args:
        name (str): (base) file name to backup.
        verbose (bool, optional): provide a user update. Defaults to True.
    """
    count = 1
    while os.path.isfile(name):
        if os.path.isfile(f'#{name}.{count}#'):
            count += 1
        else:
            os.system(f'mv {name} \#{name}.{count}\#')  # noqa
            if verbose:
                print(f'Backed up {name} to #{name}.{count}#')
            return
