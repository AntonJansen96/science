import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis
import numpy as np
from scipy.special import erf, erfinv
from science.parsing import loadCol


class GenPotentials:
    """Constructs the bias and pH potentials used in CpHMD in Python."""

    def __init__(self, dwpE: float, pKa: float, pH: float) -> None:
        assert dwpE in [0, 7.5], "currently only dwpE 7.5 or 0.0 supported!"

        self.h = dwpE       # User level input parameters.
        self.pKa = pKa
        self.pH = pH

        self.s = 0.30
        self.w = 1000  # Checked in constant_ph.cpp.
        e1 = erfinv(1 - 2.0 / self.w)
        e10 = erfinv(1 - 20.0 / self.w)
        sig0 = 0.02  # Checked in constant_ph.cpp.
        self.r = (e1 - e10) / (2 * sig0)
        self.m = 2.0 * sig0 * (2.0 * e1 - e10) / (e1 - e10)

        # self.k = 0.5 * self.h   # Initial values for Umin parameters.
        # self.a = 0.05
        # self.b = -0.1

        # Iteratively find Umin parameters k, a, b.
        # Anton: this is not implemented (yet).

        self.k = 4.7431   # Final value for dwpE = 7.5 kJ/mol.
        self.a = 0.0435   # Final value for dwpE = 7.5 kJ/mol.
        self.b = 0.0027   # Final value for dwpE = 7.5 kJ/mol.

    def __Uwall(self, lamda: float) -> float:
        A = 1 - erf(self.r * (lamda + self.m))
        B = 1 + erf(self.r * (lamda - 1 - self.m))
        return 0.5 * self.w * (A + B)

    def __Umin(self, lamda: float) -> float:
        A = -(lamda - 1 - self.b)**2 / (2 * self.a**2)
        B = -(lamda + self.b)**2 / (2 * self.a**2)
        return -self.k * (np.exp(A) + np.exp(B))

    def __Ubarrier(self, lamda: float) -> float:
        d = 0.5 * self.h
        return d * np.exp(-(lamda - 0.5)**2 / (2 * self.s**2))

    def U_bias(self, lamda: float) -> float:
        """Bias potential.

        Args:
            lamda (float): lambda coordinate value.

        Returns:
            float: potential (kJ/mol).
        """
        return self.__Uwall(lamda) + self.__Umin(lamda) + self.__Ubarrier(lamda)

    def U_pH(self, lamda: float) -> float:
        """Henderson-Hasselbalch based pH potential.

        Args:
            lamda (float): lambda coordinate value.

        Returns:
            float: potential (kJ/mol).
        """
        ph_prefactor = 0.001 * 8.3145 * 300 * np.log(10) * (self.pKa - self.pH)    
        k1 = 2.5 * self.r
        x0 = 2 * self.a

        if self.pH > self.pKa:
            ph_fraction = 1 / (1 + np.exp(-2 * k1 * (lamda - 1 + x0 * self.s)))
        else:
            ph_fraction = 1 / (1 + np.exp(-2 * k1 * (lamda - x0)))

        return ph_prefactor * ph_fraction


class InverseBoltzmann:
    """For performing Inverse-Boltzmann related tasks such as adding a (bias)
    potential to a replicas histogram or for adjusting dwpE."""

    def __init__(self, baseName: str, coordsArray: list, Nbins: int = 35, Nrange=(-0.10, 1.10)) -> None:

        matplotlib.rcParams.update({'font.size': 24})

        self.baseName: str = baseName
        self.Nbins: int = Nbins
        self.Nrange: tuple = Nrange

        self.binsArray: list = []
        self.histArray: list = []
        self.enerArray: list = []

        for coords in coordsArray:
            hist, bins = np.histogram(coords, bins=Nbins, range=self.Nrange, density=True)
            bins = [(bins[i] + bins[i + 1]) / 2 for i in range(0, len(bins) - 1)]
            self.histArray.append(hist)
            self.binsArray.append(bins)

            U = []
            for p in hist:
                U.append(8.3145 * 300 * -np.log(p))
            self.enerArray.append([E / 1000. for E in U])  # J to kJ.

        self.added: list = len(self.binsArray[0]) * [0]

    def plot(self, name: str) -> None:

        plt.figure(dpi=200)
        for idx in range(len(self.binsArray)):
            plt.plot(self.binsArray[idx], self.histArray[idx])
        plt.xlim(self.Nrange)
        plt.xlabel(r'$\lambda$-coordinate')
        plt.ylabel('Probability density')
        plt.tight_layout()
        plt.savefig(f"{self.baseName}_{name}_hist.png")

        plt.figure(dpi=200)
        for idx in range(len(self.binsArray)):
            plt.plot(self.binsArray[idx], self.enerArray[idx])
        plt.hlines(y=0, xmin=-0.1, xmax=1.1, linestyles='--', color='black', linewidth=0.5)
        plt.hlines(y=min(self.enerArray[idx]), xmin=-0.1, xmax=1.1, linestyles='--', linewidth=0.5)
        plt.xlim(self.Nrange)
        plt.ylim(-5, 20)
        plt.xlabel(r'$\lambda$-coordinate')
        plt.ylabel('Energy (kJ/mol)')
        plt.tight_layout()
        plt.savefig(f"{self.baseName}_{name}_ener.png")

    def addTestPotential(self, dwpE=5):
        for val in range(len(self.binsArray[0])):
            if self.binsArray[0][val] > 0.2 and self.binsArray[0][val] < 0.8:
                self.added[val] -= dwpE
                for idx in range(len(self.binsArray)):
                    self.enerArray[idx][val] -= dwpE
                    self.histArray[idx][val] = np.exp(-1000 * self.enerArray[idx][val] / (8.3145 * 300))

        plt.figure(dpi=200)
        plt.plot(self.binsArray[0], self.added)
        plt.hlines(y=0, xmin=-0.1, xmax=1.1, linestyles='--', color='black', linewidth=0.5)
        plt.xlim(self.Nrange)
        plt.ylim(-10, 20)
        plt.xlabel(r'$\lambda$-coordinate')
        plt.ylabel('Energy (kJ/mol)')
        plt.tight_layout()
        plt.savefig(f"{self.baseName}_added.png")

    def addBiasPotential(self, dwpE=7.5):
        bias = BiasPotential(dwpE)

        for val in range(len(self.binsArray[0])):
            E = bias.potential(self.binsArray[0][val])
            self.added[val] += E
            for idx in range(len(self.binsArray)):
                self.enerArray[idx][val] += E
                self.histArray[idx][val] = np.exp(-1000 * self.enerArray[idx][val] / (8.3145 * 300))

        plt.figure(dpi=200)
        plt.plot(self.binsArray[0], self.added)
        plt.hlines(y=0, xmin=-0.1, xmax=1.1, linestyles='--', color='black', linewidth=0.5)
        plt.xlim(self.Nrange)
        plt.ylim(-5, 20)
        plt.xlabel(r'$\lambda$-coordinate')
        plt.ylabel('Energy (kJ/mol)')
        plt.tight_layout()
        plt.savefig(f"{self.baseName}_added.png")

    def addpHPotential(self, pH: float, pKa: float):
        for val in range(len(self.binsArray[0])):
            E = 8.3145 * 300 * np.log(10)
            r = 13.51   # From bias potential 7.5 kJ/mol.
            a = 0.0435  # From bias potential 7.5 kJ/mol.

            k1 = 2.5 * r    # Where r comes from the bias potential parameters...
            x0 = 2 * a      # Where a comes from the bias potential parameters...
            if pKa > pH:
                E *= 1 / (1 + np.exp(-2 * k1 * (self.binsArray[0][val] - 1 + x0)))
            else:
                E *= 1 / (1 + np.exp(-2 * k1 * (self.binsArray[0][val] - x0)))
            E /= 1000.0

            self.added[val] += E
            for idx in range(len(self.binsArray)):
                self.enerArray[idx][val] += E
                self.histArray[idx][val] = np.exp(-1000 * self.enerArray[idx][val] / (8.3145 * 300))

        plt.figure(dpi=200)
        plt.plot(self.binsArray[0], self.added)
        plt.hlines(y=0, xmin=-0.1, xmax=1.1, linestyles='--', color='black', linewidth=0.5)
        plt.xlim(self.Nrange)
        # plt.ylim(-5, 20)
        plt.xlabel(r'$\lambda$-coordinate')
        plt.ylabel('Energy (kJ/mol)')
        plt.tight_layout()
        plt.savefig(f"{self.baseName}_added.png")


def protonation(xList: list, cutoff: float = 0.8) -> float:
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


def deprotonation(xList: list, cutoff: float = 0.8) -> float:
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


def movingDeprotonation(xList: list, cutoff: float = 0.8):
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


def getLambdaFileIndices(universe: MDAnalysis.Universe, resid: int, numChains: int = 5):
    """Returns an array containing the lambda-file indices for the specified resid.
    Takes into account ASPT, GLUT, HSPT, ARGT, LYST, TYRT.

    Args:
        universe (MDAnalysis.Universe): MDAnalysis universe variable containing the structure.
        resid (int): titratable residue resid.

    Returns:
       list: corresponding cphmd-coord numbers.
    """

    segmentAatoms      = universe.select_atoms("chainID A")
    titratableAtoms    = segmentAatoms.select_atoms('resname ASPT GLUT HSPT ARGT LYST TYRT')
    titratableResnames = list(titratableAtoms.residues.resnames)
    titratableResids   = list(titratableAtoms.residues.resids)
    targetidx          = titratableResids.index(resid)

    numASPTGLUT        = len(segmentAatoms.select_atoms('resname ASPT GLUT ARGT LYST TYRT').residues)
    numHSPT            = len(segmentAatoms.select_atoms('resname HSPT').residues)
    factor             = numASPTGLUT + 3 * numHSPT

    count = 1
    for idx in range(0, len(titratableResnames)):

        if idx == targetidx:
            array = []
            for ii in range(0, numChains):
                array.append(count + ii * factor)
            return array

        if titratableResnames[idx] in ['ASPT', 'GLUT', 'ARGT', 'LYST', 'TYRT']:
            count += 1

        elif titratableResnames[idx] == 'HSPT':
            count += 3


def theoreticalProtonation(pH: float, pKa: float) -> float:
    """Returns theoretical protonation fraction as calculated by the Henderson-Hasselbach equation,
    i.e. protonation = 1 / (1 + 10^(pH - pKa)).

    Args:
        pH (float): (solvent) pH.
        pKa (float): macroscopic pKa.

    Returns:
        float: protonation fraction.
    """

    return 1 / (1 + 10**(pH - pKa))


def theoreticalMicropKa(pH: float, protonation: float) -> float:
    """Return the theoretical microscopic pKa as calculated by the Henderson-Hasselbalch equation,
    i.e. pKa = pH - log_10(1 / f_p - 1)

    Args:
        pH (float): (solvent) pH.
        protonation (float): protonation fraction.

    Returns:
        float: theoretical microscopic pKa.
    """

    return pH - np.log10(1 / protonation - 1)


def extractCharges(proto: str, depro: str) -> None:
    """Extract the titratable atoms and charges by comparing two .itp files.
    Assumes all headers ([atoms] etc.) are removed before running.

    Args:
        proto (str): .itp file name of protonated structure.
        depro (str): .itp file name of deprotonated structure.
    """

    atoms_p   = loadCol(proto, 5)
    charge_p  = loadCol(proto, 7)
    atoms_dp  = loadCol(depro, 5)
    charge_dp = loadCol(depro, 7)

    protoDict = {}
    for idx in range(0, len(atoms_p)):
        protoDict[atoms_p[idx]] = charge_p[idx]

    deproDict = {}
    for idx in range(0, len(atoms_dp)):
        deproDict[atoms_dp[idx]] = charge_dp[idx]

    atoms = []
    A = []
    B = []
    for key in protoDict:
        try:
            if protoDict[key] != deproDict[key]:
                atoms.append(key)
                A.append(protoDict[key])
                B.append(deproDict[key])
        except KeyError:
            atoms.append(key)
            A.append(protoDict[key])
            B.append(0.000)

    print("titratable atoms", atoms)
    print("qqA", A, sum(A))
    print("qqB_1", B, sum(B))


def plotdVdl(dvdl: list, name: list = ['curve'], Nrange: tuple = (-0.1, 1.1)) -> None:
    """Makes a (combined) plot of the specified dV/dl coefficients.

    Args:
        dvdl (list): list of dV/dl coefficients in lambdagrouptypes.dat/.mdp format.
        name (list): list of names belonging to the sets of dV/dl coefficients.
        Nrange (tuple): desired lambda coordinate range.
    """

    matplotlib.rcParams.update({'font.size': 18})
    plt.figure(dpi=200)

    coords = np.arange(Nrange[0], Nrange[1], 0.01)

    for ii in range(0, len(dvdl)):

        coeffs = [float(val) for val in dvdl[ii].split()][::-1]
        energies = []

        for coord in coords:
            val = 0
            for idx in range(0, len(coeffs)):
                val += coeffs[idx] * coord**idx
            energies.append(val)

        plt.plot(coords, energies, label=name[ii])

    plt.legend()
    plt.xlabel(r"$\lambda$-coordinate")
    plt.ylabel("Energy (kJ/mol)")
    plt.tight_layout()
    plt.savefig('dvdlplot.png')
