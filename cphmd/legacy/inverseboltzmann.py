import matplotlib.pyplot as plt
import numpy as np
import math

from phbuilder.structure import Structure

def inverseBoltzmann(fileName, resid):

    def barrier(l):
        k = 4.7431      # Coded this for nothing: this is only at 7.5!
        a = 0.0435
        b = 0.0027
        d = 3.75
        s = 0.30
        w = 1000.0
        r = 13.5
        m = 0.2019
    
        A = np.exp(- (l - 1 - b)**2 / (2 * a**2) )
        B = np.exp(- (l + b)**2     / (2 * a**2) )
        C = np.exp(- (l - 0.5)**2   / (2 * s**2) )
        D = 1 - math.erf(r * (l + m))
        E = 1 + math.erf(r * (l - 1 - m))

        return -k * (A + B) + d * C + 0.5 * w * (D + E)

    numChains = 5
    resPerChain = 43

    # GET THE DATA

    pdb = Structure(fileName, verbosity=3)

    idx = 1
    data = []
    for residue in pdb.d_residues:
        if residue.d_chain == 'A':
            if residue.d_resid == resid:
                for jj in range(0, numChains):
                    print("loading lambda_{}.dat...".format(idx + resPerChain * jj))
                    data += loadCol('lambda_{}.dat'.format(idx + resPerChain * jj), 2)[50000:]
                break

            if residue.d_resname in ['ASPT', 'GLUT']:
                idx += 1

            elif residue.d_resname == 'HSPT':
                idx += 3

    # GET THE HISTOGRAM AND CORRESPONDING LAMBDA LISTS
    
    bins1 = [x / 1000. for x in range(-100, 1100, 3)] # 400 bins
    hist, bins2 = np.histogram(data, bins=bins1, density=True)
    bins2 = bins2[1:] # remove first element to make arrays equal size

    # DO THE ACTUAL BOLTZMANN INVERSION

    R = 8.3145  # gas constant (not kb because we have kJ/mol as a unit)
    T = 300     # K

    energyList = []
    for p in hist:
        energyList.append(R * T * -np.log(p))

    energyList = [E / 1000. for E in energyList] # From J to kJ.

    # PLOT AND EYEBALL

    plt.plot(bins2, energyList)
    plt.grid()
    plt.ylabel("Energy (kJ/mol)")
    plt.xlabel(r"$\lambda$-coordinate")
    plt.show()
