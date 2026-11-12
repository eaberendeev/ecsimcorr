#!/usr/bin/env python3
"""Checks the fast magnetosonic frequency of the run (wave_cases/fms_check.py)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "wave_cases"))
import fms_check  # noqa: E402

if __name__ == "__main__":
    sys.exit(fms_check.main("fms_symplectic"))
