#!/usr/bin/env python3
"""Analyze the full ECSIM-sized symplectic Weibel experiment."""
import runpy
from pathlib import Path
import sys

shared = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(shared))
runpy.run_path(str(shared / "weibel_check.py"), run_name="__main__")
