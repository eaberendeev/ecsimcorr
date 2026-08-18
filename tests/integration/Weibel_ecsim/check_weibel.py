#!/usr/bin/env python3
"""Draws the report figures and checks the numbers of a finished Weibel run
(plain ecsim scheme, this test's gen_config.py).

Usage:  python3 check_weibel.py [workdir] [--outdir DIR] [--quiet]

  workdir is optional: by default the run directory is derived from this
  test's gen_config.py (DirName_Dx_.._np_.._Dt_..) and looked up in the
  current directory and in the repo root.

Exit code: 0 = all checks passed, 1 = some check failed, 2 = bad input.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SHARED = os.path.join(os.path.dirname(HERE), "weibel_check.py")

sys.exit(subprocess.call(
    [sys.executable, SHARED] + sys.argv[1:]
    + ["--config", os.path.join(HERE, "gen_config.py")]))