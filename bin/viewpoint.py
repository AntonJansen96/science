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

# Process camera line further

for idx in range(0, len(camera)):
    if camera[idx] == '{':
        camera = camera[idx:]
        break

def writeViewBlock(file):
    file.write('\n# VIEWPOINT\n')
    file.write(f'molinfo [molinfo top] set {{center_matrix rotate_matrix scale_matrix global_matrix}} {camera}\n')

# Write viewpoint block to target file:
with open(targetFile, 'a+') as file:
    writeViewBlock(file)

# Rewrite source file to only consist of the view block:
# with open(sourceFile, 'w+') as file:
#     writeViewBlock(file)
