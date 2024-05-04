#!/usr/bin/env python3

import os
import sys
import subprocess
from science.parsing import Sanitize

"""Generate a python file with a shebang line and make it executable"""

assert len(sys.argv) == 2, "Usage: pyfile.py <filename>"

filename = Sanitize(sys.argv[1]).path(ext=".py", out=True)

open(filename, 'w+').write("#!/usr/bin/env python3\n\n")
os.chmod(filename, 0o755)
subprocess.run(["code", filename])
