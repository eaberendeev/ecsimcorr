"""Colormaps carried over from PlotScripts and the drawing primitives."""
import matplotlib
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import numpy as np

from .run import AXES

# Diverging blue-white-red map for signed fields (white band around zero).
FIELD = colors.LinearSegmentedColormap("fields", {
    "red": ((0.0, 0, 0), (0.35, 0.586, 0.586), (0.4, 0, 0), (0.45, 0, 0), (0.495, 1, 1), (0.505, 1, 1),
            (0.55, 1, 1), (0.65, 1, 1), (1.0, 0.586, 0.586)),
    "green": ((0.0, 0, 0), (0.35, 0, 0), (0.4, 0, 0), (0.45, 0.586, 0.586), (0.495, 1, 1), (0.505, 1, 1),
              (0.55, 0.586, 0.586), (0.65, 0, 0), (1.0, 0, 0)),
    "blue": ((0.0, 0.586, 0.586), (0.35, 1, 1), (0.4, 1, 1), (0.45, 1, 1), (0.495, 1, 1), (0.505, 1, 1),
             (0.55, 0, 0), (0.65, 0.586, 0.586), (1.0, 0, 0)),
}, N=1024)
# Sequential white-blue-orange-red map for densities and other positive values.
DENSITY = colors.LinearSegmentedColormap("dens", {
    "red": ((0.0, 1, 1), (0.15, 0, 0), (0.35, 0, 0), (0.55, 1, 1), (0.75, 1, 1), (1.0, 0.586, 0.586)),
    "green": ((0.0, 1, 1), (0.15, 0.586, 0.586), (0.35, 0, 0), (0.55, 0.586, 0.586), (0.75, 0, 0), (1.0, 0, 0)),
    "blue": ((0.0, 1, 1), (0.15, 1, 1), (0.35, 1, 1), (0.55, 0, 0), (0.75, 0.586, 0.586), (1.0, 0, 0)),
}, N=1024)


def percentiles(data):
    """0.1 % and 99.9 % percentiles: the colour range ignores single hot cells."""
    finite = data[np.isfinite(data)]
    return tuple(float(v) for v in np.percentile(finite, (0.1, 99.9))) if finite.size else (0.0, 0.0)


def limits(lo, hi, signed):
    """Colour limits, symmetric about zero when signed data changes sign or is zero.
    Symmetric limits select the diverging map (see cmap)."""
    if signed and lo <= 0 <= hi:
        m = max(-lo, hi) or 1.0
        return -m, m
    if hi <= lo:
        hi = lo + (abs(lo) or 1.0)
    return lo, hi


def robust_limits(data, signed):
    return limits(*percentiles(data), signed)


def cmap(limits):
    """Diverging map with white at zero for limits symmetric about zero, sequential otherwise."""
    return FIELD if limits[0] == -limits[1] else DENSITY


def _sci(axis):
    fmt = ticker.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-2, 3))
    axis.set_major_formatter(fmt)


def draw_map(fig, ax, data, a1, a2, axes, q, vmin, vmax):
    im = ax.imshow(data.T, origin="lower", extent=(a1[0], a1[-1], a2[0], a2[-1]),
                   cmap=cmap((vmin, vmax)), vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_xlabel(AXES[axes[0]])
    ax.set_ylabel(AXES[axes[1]])
    ax.set_title(q.label)
    bar = fig.colorbar(im, ax=ax, shrink=0.9, pad=0.02)
    bar.locator = ticker.MaxNLocator(5)
    _sci(bar.ax.yaxis)
    bar.update_ticks()
    return im


def draw_profile(ax, x, line, xlabel, title):
    (curve,) = ax.plot(x, line, lw=1.4, color="C0")
    ax.set_xlabel(xlabel)
    if title:
        ax.set_title(title)
    ax.grid(alpha=0.3)
    _sci(ax.yaxis)
    return curve


def draw_curves(fig, ax, curves, times, xlabel, title):
    """Curves [(x, y)] coloured by their time."""
    norm = colors.Normalize(times[0], times[-1] if times[-1] > times[0] else times[0] + 1)
    cmap = matplotlib.colormaps["viridis"]
    for (x, y), t in zip(curves, times):
        ax.plot(x, y, lw=1.0, color=cmap(norm(t)))
    fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), ax=ax, label="t", pad=0.02)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    _sci(ax.yaxis)


def draw_spacetime(fig, ax, x, times, stack, xlabel, q):
    vmin, vmax = q.limits or robust_limits(stack, q.signed)
    edges_t = _edges(times)
    im = ax.pcolormesh(_edges(x), edges_t, stack, cmap=cmap((vmin, vmax)), vmin=vmin, vmax=vmax, shading="flat")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("t")
    ax.set_title(q.label)
    bar = fig.colorbar(im, ax=ax, pad=0.02)
    _sci(bar.ax.yaxis)


def _edges(centres):
    c = np.asarray(centres, dtype=float)
    if c.size == 1:
        return np.array([c[0] - 0.5, c[0] + 0.5])
    mid = 0.5 * (c[1:] + c[:-1])
    return np.concatenate(([2 * c[0] - mid[0]], mid, [2 * c[-1] - mid[-1]]))
