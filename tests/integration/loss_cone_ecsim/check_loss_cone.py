#!/usr/bin/env python3
"""Checks the escaped fraction of the run (mirror_cases/loss_cone_check.py)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "mirror_cases"))
import loss_cone_check  # noqa: E402

if __name__ == "__main__":
    sys.exit(loss_cone_check.main("loss_cone_ecsim"))
