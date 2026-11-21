"""2D maps of plane slices, 1D profiles taken from them, and space-time diagrams."""
import math
import multiprocessing
import os
import shutil
import subprocess

import numpy as np

from . import style
from .quantities import values
from .run import AXES, IN_PLANE, select


# ---- frame selection ------------------------------------------------------
def aligned_frames(run, quantities, plane, averaged, time_range, every):
    """[(first, last, {quantity: Frame})] for steps present for every quantity."""
    per_quantity = []
    for q in quantities:
        frames = run.frames(q.owner, q.base, plane, averaged)
        if not frames:
            kind = "averaged " if averaged else ""
            raise SystemExit(f"no {kind}{q.base} files for '{q.label}' on plane {plane.tag}; see `bplot info`")
        per_quantity.append({(f.first, f.last): f for f in frames})
    common = sorted(set.intersection(*(set(f) for f in per_quantity)))
    if not common:
        raise SystemExit("the selected quantities are written at different steps; plot them separately")
    reference = select([per_quantity[0][k] for k in common], time_range, every, run)
    return [(f.first, f.last, {q: frames[(f.first, f.last)] for q, frames in zip(quantities, per_quantity)})
            for f in reference]


# ---- profiles -------------------------------------------------------------
def profile(data, kind, run, plane):
    """(coordinate, values, axis label) of a 1D profile of a plane slice."""
    coords = run.plane_axes(plane)
    axes = [AXES[a] for a in IN_PLANE[plane.axis]]
    if kind == "r":
        if plane.axis != "z":
            raise SystemExit("profile r needs a z plane")
        x, y = (c - c0 for c, c0 in zip(coords, run.center))
        r = np.hypot(*np.meshgrid(x, y, indexing="ij"))
        dr = min(run.step[0], run.step[1])
        rmax = run.config.get("CylinderDomain", {}).get("radius", min(abs(x[0]), x[-1], abs(y[0]), y[-1]))
        bins = np.rint(r / dr).astype(int).ravel()
        count = np.bincount(bins)
        total = np.bincount(bins, data.ravel())
        keep = np.flatnonzero((count > 0) & (np.arange(count.size) * dr <= rmax + 1e-9))
        return keep * dr, total[keep] / count[keep], "r"
    along, cut = kind[0], kind[1:] == "c"
    if along not in axes or kind[1:] not in ("", "c"):
        raise SystemExit(f"profile '{kind}' does not fit plane {plane.tag} (in-plane axes {''.join(axes)})")
    dim = axes.index(along)
    other = 1 - dim
    if cut:
        centre = (run.center + [0.5 * run.cells[2] * run.step[2]])[IN_PLANE[plane.axis][other]]
        index = int(np.argmin(abs(coords[other] - centre)))
        line = np.take(data, index, axis=other)
        return coords[dim], line, f"{along}   ({axes[other]} = {coords[other][index]:g})"
    return coords[dim], data.mean(axis=other), f"{along}   (mean over {axes[other]})"


# ---- workers --------------------------------------------------------------
_job = {}


def _init(job):
    _job.update(job)


def _measure(item):
    """Robust colour and profile ranges of every quantity on one frame."""
    _, _, frames = item
    run, plane, averaged, kind = _job["run"], _job["plane"], _job["averaged"], _job["profile"]
    cache, ranges = {}, {}
    for q, frame in frames.items():
        data = values(q, run, plane, frame, averaged, cache)
        line = profile(data, kind, run, plane)[1] if kind else data
        ranges[q] = (style.percentiles(data), (float(np.nanmin(line)), float(np.nanmax(line))))
    return ranges


def _render(task):
    """Draw one frame. Each worker lays its figure out once and then only swaps the data."""
    (first, last, frames), limits = task
    run, plane, averaged, kind = _job["run"], _job["plane"], _job["averaged"], _job["profile"]
    if "figure" not in _job:
        _job["figure"] = _layout(list(frames))
    fig, title, panels = _job["figure"]
    cache = {}
    for q, frame in frames.items():
        image, line = panels[q]
        data = values(q, run, plane, frame, averaged, cache)
        if image:
            clim = q.limits or (limits[q][0] if limits else style.robust_limits(data, q.signed))
            image.set_data(data.T)
            image.set(clim=clim, cmap=style.cmap(clim))
        if line:
            line.set_ydata(profile(data, kind, run, plane)[1])
            if limits:
                line.axes.set_ylim(*limits[q][1])
            else:
                line.axes.relim()
                line.axes.autoscale_view()
    title.set_text(f"{run.frame_label(next(iter(frames.values())))}     plane {plane.axis} = "
                   f"{run.plane_coordinate(plane):g}")
    path = os.path.join(_job["out"], f"{first:07d}.png")
    fig.savefig(path, dpi=_job["dpi"])
    fig.set_layout_engine("none")  # keep the first frame's layout
    return path


def _layout(quantities):
    """Figure with an empty map and/or profile per quantity: (figure, title, {q: (image, line)})."""
    import matplotlib.pyplot as plt
    run, plane, kind, only_profile = _job["run"], _job["plane"], _job["profile"], _job["only_profile"]
    n = len(quantities)
    cols = _job["cols"] or (n if n <= 3 else math.ceil(math.sqrt(n)))
    rows = math.ceil(n / cols)
    a1, a2 = run.plane_axes(plane)
    aspect = min(max((a2[-1] - a2[0]) / max(a1[-1] - a1[0], 1e-30), 0.3), 3.0)
    map_h = 0 if only_profile else 3.6 * aspect
    prof_h = 2.0 if kind else 0
    fig = plt.figure(figsize=(4.4 * cols, rows * (map_h + prof_h + 0.9) + 0.6), layout="constrained")
    grid = fig.add_gridspec(rows, cols)
    blank = np.zeros((len(a1), len(a2)))
    panels = {}
    for i, q in enumerate(quantities):
        heights = [h for h in (map_h, prof_h) if h]
        sub = grid[i // cols, i % cols].subgridspec(len(heights), 1, height_ratios=heights)
        image = line = None
        if not only_profile:
            image = style.draw_map(fig, fig.add_subplot(sub[0]), blank, a1, a2, IN_PLANE[plane.axis], q, -1, 1)
        if kind:
            x, values_, label = profile(blank, kind, run, plane)
            line = style.draw_profile(fig.add_subplot(sub[-1]), x, values_, label, q.label if only_profile else None)
        panels[q] = (image, line)
    return fig, fig.suptitle("", fontsize=13), panels


def _pool(jobs, job, count):
    jobs = max(1, min(jobs, math.ceil(count / 4)))  # a worker pays ~1 s to start matplotlib
    if jobs == 1:
        _init(job)
        return None
    return multiprocessing.get_context("fork").Pool(jobs, _init, (job,))


def _map(pool, fn, items):
    return list(pool.imap(fn, items)) if pool else [fn(i) for i in items]


def _shared_limits(pool, items, quantities):
    measured = _map(pool, _measure, items)
    limits = {}
    for q in quantities:
        lo = min(m[q][0][0] for m in measured)
        hi = max(m[q][0][1] for m in measured)
        plo = min(m[q][1][0] for m in measured)
        phi = max(m[q][1][1] for m in measured)
        pad = 0.05 * (phi - plo) or abs(phi) * 0.05 or 1.0
        limits[q] = (style.limits(lo, hi, q.signed), (plo - pad, phi + pad))
    return limits


# ---- commands -------------------------------------------------------------
def maps(run, quantities, plane, args):
    items = aligned_frames(run, quantities, plane, args.avg, args.time, args.every)
    out = args.out or os.path.join(run.path, "Plots", f"maps_{plane.tag}{'_avg' if args.avg else ''}")
    os.makedirs(out, exist_ok=True)
    job = dict(run=run, plane=plane, averaged=args.avg, profile=args.profile, only_profile=args.only_profile,
               cols=args.cols, dpi=args.dpi, out=out)
    pool = _pool(args.jobs, job, len(items))
    try:
        limits = _shared_limits(pool, items, quantities) if args.scale == "run" and len(items) > 1 else None
        tasks = [(item, limits) for item in items]
        paths = []
        for path in (pool.imap(_render, tasks) if pool else map(_render, tasks)):
            paths.append(path)
            print(f"\r{len(paths)}/{len(items)} frames", end="", flush=True)
        print()
    finally:
        if pool:
            pool.close()
    print(f"frames: {out}")
    if args.movie:
        print(f"movie:  {movie(paths, os.path.join(out, 'movie.mp4'), args.fps)}")


def profiles(run, quantities, plane, args):
    """Profiles of all selected frames: curves coloured by time and a space-time map."""
    import matplotlib.pyplot as plt
    items = aligned_frames(run, quantities, plane, args.avg, args.time, args.every)
    job = dict(run=run, plane=plane, averaged=args.avg, profile=args.profile)
    pool = _pool(args.jobs, job, len(items))
    try:
        lines = _map(pool, _profiles_of, items)
    finally:
        if pool:
            pool.close()
    times = np.array([run.time(0.5 * (first + last)) for first, last, _ in items])
    fig, axes = plt.subplots(len(quantities), 2, figsize=(12, 3.4 * len(quantities)), layout="constrained",
                             squeeze=False)
    for row, q in enumerate(quantities):
        x, label = lines[0][q][0], lines[0][q][2]
        stack = np.array([line[q][1] for line in lines])
        style.draw_curves(fig, axes[row, 0], [(x, line) for line in stack], times, label, q.label)
        style.draw_spacetime(fig, axes[row, 1], x, times, stack, label, q)
    span = f"{run.frame_label(items[0][2][quantities[0]])}  …  {run.frame_label(items[-1][2][quantities[0]])}"
    fig.suptitle(f"profile {args.profile}, plane {plane.axis} = {run.plane_coordinate(plane):g}\n{span}")
    out = args.out or os.path.join(run.path, "Plots",
                                   f"profiles_{args.profile}_{plane.tag}{'_avg' if args.avg else ''}.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=args.dpi)
    if args.data:
        np.savez(os.path.splitext(out)[0] + ".npz", time=times,
                 **{f"{q.label.replace(' ', '_')}": np.array([line[q][1] for line in lines]) for q in quantities},
                 **{f"{q.label.replace(' ', '_')}_coord": lines[0][q][0] for q in quantities})
    print(out)
    return fig


def _profiles_of(item):
    _, _, frames = item
    cache = {}
    return {q: profile(values(q, _job["run"], _job["plane"], f, _job["averaged"], cache), _job["profile"],
                       _job["run"], _job["plane"])
            for q, f in frames.items()}


def movie(paths, target, fps):
    if not shutil.which("ffmpeg"):
        return "ffmpeg not found, frames only"
    listing = target + ".txt"
    with open(listing, "w", encoding="utf-8") as f:
        for path in paths + paths[-1:]:
            f.write(f"file '{os.path.abspath(path)}'\nduration {1 / fps}\n")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", listing,
                    "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p", target], check=True)
    os.remove(listing)
    return target
