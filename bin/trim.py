#!/usr/bin/env python3

import os
import sys

from science.parsing import Sanitize

# Ensure the user provides a file.
assert len(sys.argv) == 2, "Usage: trim.py <file>"

# Get the file path and sanitize it.
file = Sanitize(sys.argv[1]).path(ext=[".png", ".jpg", ".jpeg"])

# Trim whitespace from an image.
os.system(f"convert {file} -trim {file}")
