#!/usr/bin/env python3

# For processing camera location+angle (viewpoints) in VMD.

import sys

assert len(sys.argv) == 3, 'Specify source file and target file'

sourceFile = sys.argv[1]
targetFile = sys.argv[2]

# Get viewpoint from source file

with open(sourceFile) as file:
    for line in file.read().splitlines():
        if 'set viewpoints(' in line:
            camera = line
            break

def writeViewBlock(file):
    file.write('\n# VIEWPOINT BLOCK ##############################################################\n')
    file.write('set viewplist {}\n')
    file.write(f'{camera}\n')
    file.write('lappend viewplist [molinfo top]\n')
    file.write('foreach v $viewplist {\n')
    file.write('  molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)\n}\n')
    file.write('unset viewplist\n')

# Write viewpoint block to target file:
with open(targetFile, 'a+') as file:
    writeViewBlock(file)

# Rewrite source file to only consist of the view block:
with open(sourceFile, 'w+') as file:
    writeViewBlock(file)
