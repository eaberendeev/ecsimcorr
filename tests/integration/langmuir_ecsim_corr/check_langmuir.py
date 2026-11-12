#!/usr/bin/env python3
"""Compares the Langmuir wave of the run with kinetic theory (wave_cases/langmuir.py)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "wave_cases"))
import langmuir  # noqa: E402

if __name__ == "__main__":
    sys.exit(langmuir.check(sys.argv[1] if len(sys.argv) > 1 else "."))
