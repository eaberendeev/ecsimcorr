"""bplot — plots of a beren3d run.

  bplot info                         what the run has written
  bplot maps Ez Bz n                 2D maps of plane slices, one PNG per output step
  bplot profiles n Ephi --profile r  1D profiles of all steps and a space-time diagram
  bplot series energyField* '*Particles'
  bplot spectra Electrons

The run directory is the current one; give another as the first argument:
  bplot maps ../Res_Test Ez
"""
import argparse
import os
import sys

from . import quantities as quantity_names

EXAMPLES = """examples:
  bplot maps Ez Bz n                       fields and every species density, all steps
  bplot maps Ez=0.03 Bz=0.17:0.23 --avg    fixed colour limits, averaged files
  bplot maps Er Ephi n@Ions --profile r    each map with its azimuthal mean underneath
  bplot maps n --plane x --time 100:400 --every 5 --movie
"""


def _common_time(p):
    p.add_argument("--time", metavar="A:B", help="time range A:B, A:, :B; one value = nearest step; 'last'")
    p.add_argument("--every", type=int, default=1, metavar="K", help="take every K-th of the selected steps")


def _common_out(p):
    p.add_argument("--out", help="output file (directory for maps); default RUN/Plots/...")
    p.add_argument("--dpi", type=int, default=120)


def _common_slices(p):
    p.add_argument("quantities", nargs="+", metavar="Q", help="quantities, see the list above")
    p.add_argument("--plane", help="x, y or z (plane nearest the middle), z=6.0 (nearest), z3 (file index); default: z")
    p.add_argument("--avg", action="store_true", help="use the time-averaged slices (average_steps)")
    p.add_argument("-j", "--jobs", type=int, default=os.cpu_count(), help="parallel processes")
    _common_time(p)
    _common_out(p)


def parser():
    top = argparse.ArgumentParser(prog="bplot", description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    sub = top.add_subparsers(dest="command", required=True)

    sub.add_parser("info", help="list planes, steps, table columns, probes and spectra")

    names = quantity_names.__doc__ + "\n" + EXAMPLES
    maps = sub.add_parser("maps", help="2D maps of plane slices", description=names,
                          formatter_class=argparse.RawTextHelpFormatter)
    _common_slices(maps)
    maps.add_argument("--profile", metavar="KIND", help="add a 1D profile under each map: r, x, y, z, xc, yc, zc")
    maps.add_argument("--only-profile", action="store_true", help="draw the profiles without the maps")
    maps.add_argument("--scale", choices=("run", "frame"), default="run",
                      help="colour limits shared by all frames (run) or per frame")
    maps.add_argument("--cols", type=int, help="panels per row")
    maps.add_argument("--movie", action="store_true", help="also join the frames into movie.mp4")
    maps.add_argument("--fps", type=float, default=8)

    profiles = sub.add_parser("profiles", help="1D profiles over time and a space-time diagram", description=names,
                              formatter_class=argparse.RawTextHelpFormatter)
    _common_slices(profiles)
    profiles.add_argument("--profile", default="r", metavar="KIND", help="r (default), x, y, z, xc, yc, zc")
    profiles.add_argument("--data", action="store_true", help="also save the profiles to .npz")
    profiles.add_argument("--show", action="store_true", help="open a window")

    series = sub.add_parser("series", help="time series from energy.txt, boundary.txt, continuity.txt, probes")
    series.add_argument("columns", nargs="+", metavar="COLUMN", help="column names or patterns: 'energyField*' "
                        "'*Particles'; probes as Ex@p0")
    _common_time(series)
    series.add_argument("--avg", type=int, default=1, metavar="N", help="mean over blocks of N rows")
    series.add_argument("--smooth", type=int, default=1, metavar="N", help="running mean over N rows")
    series.add_argument("--rel", action="store_true", help="subtract the first value")
    series.add_argument("--log", action="store_true", help="log scale of |value|")
    series.add_argument("--split", action="store_true", help="one panel per column")
    series.add_argument("--x", choices=("time", "step"), default="time")
    series.add_argument("--show", action="store_true", help="open a window")
    _common_out(series)

    spectra = sub.add_parser("spectra", help="energy spectra coloured by time")
    spectra.add_argument("species", nargs="*", help="species (default: all)")
    _common_time(spectra)
    spectra.add_argument("--rebin", type=int, default=1, metavar="N", help="merge N neighbouring energy bins")
    spectra.add_argument("--logx", action="store_true")
    spectra.add_argument("--show", action="store_true", help="open a window")
    _common_out(spectra)
    return top


def _split_run(argv):
    """Take the run directory off the command line: `bplot RUN maps ...` or `bplot maps RUN ...`."""
    for i in (0, 1):
        if len(argv) > i and os.path.isdir(argv[i]):
            return argv[i], argv[:i] + argv[i + 1:]
    return ".", argv


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    path, argv = _split_run(argv)
    args = parser().parse_args(argv)
    if not os.path.isfile(os.path.join(path, "system_config.json")):
        raise SystemExit(f"{os.path.abspath(path)} is not a run directory (no system_config.json)")
    import matplotlib
    if not getattr(args, "show", False):
        matplotlib.use("Agg")
    from .run import Run
    run = Run(path)
    if args.command == "info":
        from .info import info
        info(run)
        return
    if args.command in ("maps", "profiles"):
        from . import maps
        try:
            selected = quantity_names.parse(args.quantities, run.species)
        except ValueError as error:
            raise SystemExit(str(error))
        plane = pick_plane(run, selected, args.plane, args.avg)
        result = (maps.maps if args.command == "maps" else maps.profiles)(run, selected, plane, args)
    else:
        from . import curves
        result = (curves.series(run, args.columns, args) if args.command == "series"
                  else curves.spectra(run, args.species, args))
    if getattr(args, "show", False) and result is not None:
        import matplotlib.pyplot as plt
        plt.show()


def pick_plane(run, selected, text, averaged):
    q = selected[0]
    planes = [key[2] for key in run.slices if key[:2] == (q.owner, q.base) and key[3] == averaged]
    middle = [0.5 * n * d for n, d in zip(run.cells, run.step)]
    planes.sort(key=lambda p: ("zyx".index(p.axis), abs(run.plane_coordinate(p) - middle["xyz".index(p.axis)])))
    if not planes:
        raise SystemExit(f"no {'averaged ' if averaged else ''}slices of '{q.label}'; see `bplot info`")
    if not text:
        return planes[0]
    axis, rest = text[0].lower(), text[1:]
    same = [p for p in planes if p.axis == axis]
    if not same:
        raise SystemExit(f"no {axis} plane for '{q.label}'; planes: {' '.join(p.tag for p in planes)}")
    if not rest:
        return same[0]
    if rest.startswith("="):
        return min(same, key=lambda p: abs(run.plane_coordinate(p) - float(rest[1:])))
    found = [p for p in same if p.index == int(rest)]
    if not found:
        raise SystemExit(f"no plane {text}; planes: {' '.join(p.tag for p in planes)}")
    return found[0]
