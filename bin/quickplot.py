#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argparse
import argcomplete

# Set default parameters.
default_cropmarin = 10
default_dpi = 150
default_fontsize = 12
default_linewidth = 0.5
default_output = "out.png"


def parsecmdline() -> argparse.Namespace:
    """Perform command line parsing using the Python argparse module.

    Returns:
        argparse.Namespace: parsed command line arguments.
    """

    description = "quickplot make quick python matplotlib plots from the command line."
    parser = argparse.ArgumentParser(prog="quickplot", description=description)

    # The input file name is a positional argument.
    parser.add_argument(
        "file",
        nargs="*",
        default="nothingspecified",
        help="Input data file.",
    )

    parser.add_argument(
        "-nl",
        required=False,
        dest="nolegend",
        action="store_const",
        const=1,
        help="Do not plot the legend.",
    )

    parser.add_argument(
        "-s",
        required=False,
        dest="save",
        nargs="?",
        const=default_output,
        help=f"Plot output file name (defaults to {default_output} if no file specified).",
    )

    parser.add_argument(
        "-dpi",
        required=False,
        dest="dpi",
        action="store",
        help=f"Set dpi (defaults to {default_dpi}).",
        type=int,
        default=default_dpi,
    )

    parser.add_argument(
        "-fo",
        required=False,
        dest="fontsize",
        action="store",
        help=f"Set font size (defaults to {default_fontsize}).",
        type=float,
        default=default_fontsize,
    )

    parser.add_argument(
        "-xl",
        required=False,
        dest="xlabel",
        action="store",
        help="Set x-axis label. This will override values found in the input data (eg .xvg).",
    )

    parser.add_argument(
        "-yl",
        required=False,
        dest="ylabel",
        action="store",
        help="Set y-axis label. This will override values found in the input data (eg .xvg).",
    )

    parser.add_argument(
        "-axis",
        required=False,
        dest="axis",
        action="store",
        help="Set axis limits. Format: xmin,xmax,ymin,ymax.",
    )

    # Required for auto-completing using argcomplete.
    argcomplete.autocomplete(parser)

    # Do the actual parsing.
    CLI = parser.parse_args()

    # If no input file is specified, raise an error.
    if CLI.file == "nothingspecified":
        parser.error("No input file specified.")

    # Make the output file a .png if no extension is provided.
    if CLI.save is not None and "." not in CLI.save:
        CLI.save = CLI.save + ".png"

    return CLI


class QuickPlot:
    """QuickPlot class to handle the plotting of data using matplotlib."""

    def __init__(self, CLI: argparse.Namespace) -> None:
        """Initializes the QuickPlot class and sanitizes the input arguments.

        Args:
            CLI (argparse.Namespace): the parsed command line arguments.
        """

        from science.parsing import Sanitize  # Lazy import.

        # Parse the input string and perform sanitization.
        self.inputData = {}
        for path in CLI.file:
            if path.count(":") == 0 or path[-1] == ":":
                self.inputData[Sanitize(path).path()] = [0, 1]
            elif path.count(":") == 1:
                head = Sanitize(path.split(":")[0]).path()
                tail = path.split(":")[1]
                self.inputData[head] = [int(num) for num in tail.split(",")]
            else:
                print("Invalid input format")

        # updateStr = "quickplot: plotting "
        # for key in self.inputData:
        #     updateStr += f"{key} col {self.inputData[key]}, "
        # print(updateStr[:-2])

        # Set nolegend to True if the flag is set.
        self.nolegend: bool = CLI.nolegend != None

        # Set save to False if the flag is not set, otherwise sanitize the path.
        self.save = (
            False
            if CLI.save is None
            else Sanitize(CLI.save, "output").path(
                out=True,
                ext=[".pdf", ".png", ".jpeg", ".svg", ".eps", ".tiff", ".bmp", ".gif"],
            )
        )

        # Set dpi and check if it is within the range [10, 1000].
        self.dpi: int = Sanitize(CLI.dpi, "dpi").num(Range=[10, 1000])

        # Initialize axis to False or array of floats and check type and length.
        # Also check if the values are in the correct order.
        if CLI.axis is not None:
            try:
                self.axis = [float(val) for val in CLI.axis.split(",")]
                if (
                    len(self.axis) != 4
                    or (self.axis[0] > self.axis[1])
                    or (self.axis[2] > self.axis[3])
                ):
                    print("invalid axis values. ignoring...")
                    self.axis = False
            except:
                print("invalid axis values. ignoring...")
                self.axis = False
        else:
            self.axis = False

        # Set xlabel and ylabel to False if the flag is not set.
        self.xlabel = False if CLI.xlabel is None else CLI.xlabel
        self.ylabel = False if CLI.ylabel is None else CLI.ylabel

        # Set fontsize and check if it is within the range [1, 100].
        self.fontsize = Sanitize(CLI.fontsize, "fontsize").num(signed=True)

    def run(self) -> None:
        """Run the QuickPlot class to plot the data using matplotlib."""

        import matplotlib.pyplot as plt  # Lazy import.

        # Set DPI
        plt.figure(dpi=self.dpi)

        # Set font size.
        plt.rcParams.update({"font.size": self.fontsize})

        for key in self.inputData:
            data, xlabel, ylabel, legendList = self.loadxvgComplete(
                fname=key, col=self.inputData[key]
            )

            # Plot the data.
            t = data[0]
            for idx in range(1, len(data)):
                plt.plot(
                    t, data[idx], label=legendList[idx - 1], linewidth=default_linewidth
                )

        # Set labels.
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        # Set axis limits.
        if self.axis:
            plt.axis(self.axis)

        # Set legend.
        if self.xlabel:
            plt.xlabel(self.xlabel)

        # Set ylabel.
        if self.ylabel:
            plt.ylabel(self.ylabel)

        # Disable legend if nolegend is set.
        if not self.nolegend:
            plt.legend()

        # Save the plot if the save flag is set.
        if self.save:
            import os

            plt.savefig(self.save)
            os.system(
                f"convert {self.save} -trim -bordercolor white -border {default_cropmarin}x{default_cropmarin} {self.save}"
            )
            os.system(f"code {self.save}")
        else:  # Show the plot if the save flag is not set.
            plt.show()

    def loadxvgComplete(self, fname: str, col: list):
        """Load data from an xvg file and return the data, xlabel, ylabel, and legendList.

        Args:
            fname (str): file name.
            col (list): list of columns to load.

        Returns:
            tuple: data, xlabel, ylabel, legendList.
        """

        data = self.loadxvg(fname=fname, col=col)

        # Load the x-axis and y-axis labels and the legend list.
        legendList = []
        for line in open(fname):
            if line[0] == "@":
                if "xaxis" in line:
                    xlabel = line[19:-2]
                elif "yaxis" in line:
                    ylabel = line[19:-2]
                elif 'legend "' in line:
                    legendList.append(line[13:-2])
            elif line[0] not in ["#", "@"]:
                break

        return data, xlabel, ylabel, legendList


    def loadxvg(self, fname: str, col: list = [0, 1], dt: int = 1, b: int = 0):
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

if __name__ == "__main__":
    # Run the QuickPlot class with the parsed command line arguments.
    QuickPlot(parsecmdline()).run()
