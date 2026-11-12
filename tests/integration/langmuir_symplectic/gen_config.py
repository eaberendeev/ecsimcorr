#!/usr/bin/env python3
"""Langmuir wave (Bohm-Gross frequency, Landau damping) with the symplectic scheme.

The case and its parameters are in tests/integration/wave_cases/langmuir.py.
"""

import sys
from pathlib import Path

# The integration runner executes a copy of this file from the project root.
_root = next(path for path in Path(__file__).resolve().parents if (path / "srcBeren").is_dir())
sys.path.insert(0, str(_root / "tests" / "integration" / "wave_cases"))
import langmuir  # noqa: E402

# 1M markers x 490 steps of the explicit scheme
NumProcs = 8
DirName = "Res_Langmuir_symplectic"

system_config, particles_config = langmuir.configuration("symplectic")


def generate_config():
    langmuir.write_configuration("symplectic")
    return DirName


if __name__ == "__main__":
    generate_config()
