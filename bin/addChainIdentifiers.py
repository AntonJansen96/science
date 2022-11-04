#!/bin/python3

import sys, shelve, os, string

inputFile = sys.argv[1]
outputFile = sys.argv[2]

# Set/update variable to universe
def add(varName, value):
    shelve.open('universe')[varName] = value

# Check whether universe contains a certain varName
def has(varName):
    return varName in shelve.open('universe')

# Retrieve variable from universe
def get(varName):
    if has(varName):
        return shelve.open('universe')[varName]

    data = eval(input("couldn't retrieve var \"{0}\" from  Enter manually: ".format(varName)))
    print("add {0} = {1} {2}".format(varName, data, type(data)))
    add(varName, data)
    return data

# Stores the information for a residue.
class Residue:
    def __init__(self, atoms, resname, chain, resid, x, y, z):
        self.d_atoms   = atoms      # list      holds atom types
        self.d_resname = resname    # string    holds residue name
        self.d_chain   = chain      # string    holds chain name (A, B, etc.)
        self.d_resid   = resid      # int       holds residue number
        self.d_x       = x          # list      holds x-coordinates
        self.d_y       = y          # list      holds y-coordinates
        self.d_z       = z          # list      holds z-coordinates

# Stores the information for the periodic box.
class Crystal:
    def __init__(self, a, b, c, alpha, beta, gamma, space, Z):
        self.d_a       = a          # Angstroms
        self.d_b       = b          # Angstroms
        self.d_c       = c          # Angstroms
        self.d_alpha   = alpha      # degrees
        self.d_beta    = beta       # degrees
        self.d_gamma   = gamma      # degrees
        self.d_space   = space      # Space group
        self.d_Z       = Z          # Z value

# Loads structure (pdb/gro) into d_residues.
def load(name):
    extension = os.path.splitext(name)[1]

    if (extension == ".pdb"):
        read_pdb(name)

# Writes d_residues to a structure (pdb/gro) file.
def write(name):
    extension = os.path.splitext(name)[1]

    if (extension == ".pdb"):
        write_pdb(name)

# Load a .pdb file into d_residues.
def read_pdb(name):

    with open(name) as file:

        atomLines = []

        for line in file.readlines():

            # Get title.
            if (line[0:6] == "TITLE "):
                add('d_title', line[7:80].strip())

            # Get periodic box information (if any).
            elif (line[0:6] == "CRYST1"):
                add('d_box', Crystal(float(line[6:15]), float(line[15:24]), float(line[24:33]), float(line[33:40]), float(line[40:47]), float(line[47:54]), line[55:66], int(line[66:70])))

            # Get the ATOM lines.
            elif (line[0:6] == "ATOM  "):
                atomLines.append(line)

            # Only import the first MODEL...
            elif (line[0:6] == "ENDMDL"):
                break

    # Loop through the atomLines and create a list of Residue objects.

    d_residues = []
    atoms      = []
    x          = []
    y          = []
    z          = []
    lastLine   = False

    for idx in range(0, len(atomLines)):

        atoms.append(atomLines[idx][12:16].strip())
        x.append(float(atomLines[idx][30:38]))
        y.append(float(atomLines[idx][38:46]))
        z.append(float(atomLines[idx][46:54]))

        try:
            currentResID = int(atomLines[idx][22:27]) # fix
            nextResID    = int(atomLines[idx + 1][22:27])
        except IndexError:
            lastLine = True

        if (currentResID != nextResID or lastLine):
            
            currentResName = atomLines[idx][17:21].strip()
            currentChain   = atomLines[idx][21:22]
            
            # Create the Residue object.
            d_residues.append(Residue(atoms, currentResName, currentChain, currentResID, x, y, z))

            # Reset.
            atoms = []
            x     = []
            y     = []
            z     = []

    # Add the list of Residues to universe.
    add('d_residues', d_residues)

def write_pdb(name):
    with open(name, 'w') as file:
        if has('d_title'):
            file.write("TITLE     {0}\n".format(get('d_title')))

        if has('d_box'):
            cryst = get('d_box')
            file.write("CRYST1{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f} {:11s}{:>4d}\n".format(cryst.d_a, cryst.d_b, cryst.d_c, cryst.d_alpha, cryst.d_beta, cryst.d_gamma, cryst.d_space, cryst.d_Z))

        file.write("MODEL {:8d}\n".format(1))

        atomNumber = 1
        for residue in get('d_residues'):
            for idx in range(0, len(residue.d_atoms)):
                
                atom = residue.d_atoms[idx]
                if len(atom) == 3:
                    atom = ' ' + atom

                file.write("{:6s}{:5d} {:^4s}{:1s}{:4s}{:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', atomNumber % 100000, atom, '', residue.d_resname, residue.d_chain, residue.d_resid % 10000, '', residue.d_x[idx], residue.d_y[idx], residue.d_z[idx]))
                atomNumber += 1

        file.write("TER\nENDMDL\n")

# MAIN CODE

letters = string.ascii_uppercase + string.ascii_lowercase

load(inputFile)
residues = get('d_residues')

jj = 0

try:
    for idx in range(0, len(residues)):
        if residues[idx].d_resname not in ['SOL', 'BUF', 'POPC', 'NA', 'CL', 'SOD', 'TIP3', 'CLA']:
            residues[idx].d_chain = letters[jj]

        if residues[idx].d_resid + 1 != residues[idx + 1].d_resid and jj != 51:
            jj += 1
except:
    IndexError

add('d_residues', residues)
write(outputFile)

os.remove('universe')
