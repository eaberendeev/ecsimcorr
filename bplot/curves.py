"""Time series from the text tables and energy spectra."""
import fnmatch
import os

import numpy as np

from . import style
from .run import Frame, select


def match_columns(columns, patterns):
    names = []
    for pattern in patterns:
        found = [c for c in columns if fnmatch.fnmatchcase(c, pattern) and c not in names]
        if not found:
            close = [c for c in columns if pattern.lower().strip("*") in c.lower()][:8]
            hint = f"; similar: {' '.join(close)}" if close else "; see `bplot info`"
            raise SystemExit(f"no column matches '{pattern}'{hint}")
        names += found
    return names


def block_mean(x, y, n):
    if n <= 1:
        return x, y
    m = len(y) // n * n
    return x[:m].reshape(-1, n).mean(1), y[:m].reshape(-1, n).mean(1)


def running_mean(y, n):
    if n <= 1:
        return y
    kernel = np.ones(n) / n
    pad = np.pad(y, (n // 2, n - 1 - n // 2), mode="edge")
    return np.convolve(pad, kernel, mode="valid")


def series(run, patterns, args):
    import matplotlib.pyplot as plt
    columns = run.tables()
    names = match_columns(columns, patterns)
    panels = [[n] for n in names] if args.split else [names]
    fig, axes = plt.subplots(len(panels), 1, figsize=(10, 2.6 + 2.2 * len(panels)), sharex=True,
                             layout="constrained", squeeze=False)
    for ax, group in zip(axes[:, 0], panels):
        for name in group:
            steps, y = columns[name]
            frames = select([Frame(i, int(s), int(s)) for i, s in enumerate(steps)], args.time, args.every, run)
            index = np.array([f.seq for f in frames], dtype=int)
            x = steps[index] if args.x == "step" else run.time(steps[index])
            y = y[index]
            if args.rel:
                y = y - y[0]
            x, y = block_mean(x, y, args.avg)
            y = running_mean(y, args.smooth)
            ax.plot(x, np.abs(y) if args.log else y, lw=1.2, label=name)
        if args.log:
            ax.set_yscale("log")
        else:
            style._sci(ax.yaxis)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, ncol=max(1, len(group) // 8))
    axes[-1, 0].set_xlabel("step" if args.x == "step" else "t")
    notes = [f"block mean {args.avg}" if args.avg > 1 else "", f"running mean {args.smooth}" if args.smooth > 1 else "",
             "minus initial value" if args.rel else "", "|value|" if args.log else ""]
    fig.suptitle(", ".join(n for n in notes if n) or os.path.basename(run.path))
    tag = "_".join(p.replace("*", "") for p in patterns)[:60] or "series"
    return _save(fig, run, args.out, f"series_{tag}.png", args.dpi)


def spectra(run, species, args):
    import matplotlib.pyplot as plt
    species = species or run.species
    fig, axes = plt.subplots(1, len(species), figsize=(6 * len(species), 4.5), layout="constrained", squeeze=False)
    for ax, name in zip(axes[0], species):
        records = run.spectra(name)
        if not records:
            raise SystemExit(f"no energy spectra for '{name}'")
        frames = select([Frame(i, r[0], r[0]) for i, r in enumerate(records)], args.time, args.every, run)
        chosen = [records[f.seq] for f in frames]
        times = np.array([run.time(r[0]) for r in chosen])
        style.draw_curves(fig, ax, [rebin(r[1], r[2], args.rebin) for r in chosen], times, "energy", name)
        ax.set_yscale("log")
        if args.logx:
            ax.set_xscale("log")
        ax.set_ylabel("count")
    return _save(fig, run, args.out, "spectra.png", args.dpi)


def rebin(energy, count, n):
    """Merge n neighbouring bins: mean energy, total count."""
    m = len(count) // n * n
    if n <= 1 or m == 0:
        return energy, count
    return energy[:m].reshape(-1, n).mean(1), count[:m].reshape(-1, n).sum(1)


def _save(fig, run, out, name, dpi):
    path = out or os.path.join(run.path, "Plots", name)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.savefig(path, dpi=dpi)
    print(path)
    return fig
