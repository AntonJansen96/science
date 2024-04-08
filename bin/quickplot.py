#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argparse
import argcomplete
import matplotlib.pyplot as plt

# import science.cphmd
from science.parsing import Sanitize
# import science

print(Sanitize("hello world").path(ext="txt"))

def parsecmdline() -> argparse.Namespace:
    """Perform command line parsing using python argparse module.

    Returns:
        argparse.Namespace: parsed command line arguments.
    """

    # Main tool description.
    desc1 = "quickplot make quick python matplotlib plots from the command line."

    # Epilogue.
    desc2 = ""

    parser = argparse.ArgumentParser(prog="quickplot", description=desc1, epilog=desc2)

    parser.add_argument(
        "-f",
        required=True,
        dest="file",
        action="store",
        help="Specify data file.",
    )

    # Required for auto-completing using argcomplete.
    argcomplete.autocomplete(parser)

    # Do the actual parsing.
    CLI = parser.parse_args()

    # Object containing all the parsed key-value pairs.
    return CLI


def parsexvgfile():
    pass


if __name__ == "__main__":
    CLI = parsecmdline()
    print(CLI)

    print(science.cphmd.deprotonation([1.0, 1.0, 0.0]))

# cmdline = sys.argv[1:]

# if not cmdline:
#     print("Usage: quickplot.py file.xvg:1,2,3")
#     sys.exit(1)

# for arg in cmdline:
#     denomidx = arg.index(":")
#     fname = arg[:denomidx]
#     cols = [int(val) for val in arg[denomidx + 1:].split(',')]

#     data = loadxvg(fname=fname, col=cols)
#     # first column is always the x-axis, rest of columns are y-axis.
#     for idx in range(1, len(cols)):
#         plt.plot(data[0], data[idx], label=f"{fname} col {cols[idx]}")

# plt.legend()
# plt.show()
# # plt.savefig("quick.png")
# # os.system("code quick.png")
