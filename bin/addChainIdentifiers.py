#!/usr/bin/env python3

import sys
import string

from science.parsing import Sanitize
from science.parsing import Structure

assert len(sys.argv) == 3, 'Specify source file and target file'

sourceFile = Sanitize(sys.argv[1]).path(ext=['.gro', '.pdb'])
outputFile = Sanitize(sys.argv[2]).path(ext=['.pdb'], out=True)

letters = string.ascii_uppercase + string.ascii_lowercase

pdb = Structure(sourceFile)

ii = 0
try:
    for idx in range(0, len(pdb.d_residues)):
        if pdb.d_residues[idx].d_resname not in ['SOL', 'BUF', 'POPC', 'NA', 'CL', 'SOD', 'TIP3', 'CLA']:
            pdb.d_residues[idx].d_chain = letters[ii]

        if pdb.d_residues[idx].d_resid + 1 != pdb.d_residues[idx + 1].d_resid and ii != 51:
            ii += 1

except IndexError:
    pass

pdb.write(outputFile)
