import os
import sys

if __package__ in (None, ""):  # `python3 path/to/bplot`
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from bplot.cli import main
else:
    from .cli import main

main()
