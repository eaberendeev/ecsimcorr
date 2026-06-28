import numpy as np
import matplotlib.colors as col
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

_FIELD_CDICT = {
    "red":   ((0.0, 0.0, 0.0), (0.35, 0.586, 0.586), (0.4, 0, 0),
              (0.45, 0, 0), (0.495, 1, 1), (0.505, 1, 1),
              (0.55, 1, 1), (0.65, 1, 1), (1.0, 0.586, 0.586)),
    "green": ((0.0, 0.0, 0.0), (0.35, 0, 0), (0.4, 0, 0),
              (0.45, 0.586, 0.586), (0.495, 1, 1), (0.505, 1, 1),
              (0.55, 0.586, 0.586), (0.65, 0, 0), (1.0, 0, 0)),
    "blue":  ((0.0, 0.586, 0.586), (0.35, 1, 1), (0.4, 1, 1),
              (0.45, 1, 1), (0.495, 1, 1), (0.505, 1, 1),
              (0.55, 0, 0), (0.65, 0.586, 0.586), (1.0, 0, 0)),
}
_DENS_CDICT = {
    "red":   ((0.0, 1, 1), (0.15, 0, 0), (0.35, 0, 0),
              (0.55, 1, 1), (0.75, 1, 1), (1.0, 0.586, 0.586)),
    "green": ((0.0, 1, 1), (0.15, 0.586, 0.586), (0.35, 0, 0),
              (0.55, 0.586, 0.586), (0.75, 0, 0), (1.0, 0, 0)),
    "blue":  ((0.0, 1, 1), (0.15, 1, 1), (0.35, 1, 1),
              (0.55, 0, 0), (0.75, 0.586, 0.586), (1.0, 0, 0)),
}
CMAP_FIELD = col.LinearSegmentedColormap("fields", _FIELD_CDICT, N=1024)
CMAP_DENS = col.LinearSegmentedColormap("dens", _DENS_CDICT, N=1024)

DIVERGING_QUANTITIES = {
    "Ex", "Ey", "Ez", "Er", "Ephi",
    "Bx", "By", "Bz", "Br", "Bphi",
    "Jx", "Jy", "Jz", "Jr", "Jphi",
    "Prp", "Prz", "Pzp",
}


def auto_cmap(quantity):
    if quantity in DIVERGING_QUANTITIES:
        return CMAP_FIELD
    return CMAP_DENS


def auto_title(quantity, species):
    if species:
        return f"{species} {quantity}"
    return quantity


def render_2d(ax, fig, data, extent, cmap, vmin, vmax, title, plane_labels):
    plotdata = np.swapaxes(data, 0, 1)
    kwargs = {"cmap": cmap, "origin": "lower", "extent": extent, "aspect": "auto"}
    if vmin is not None and vmax is not None:
        kwargs["vmin"] = vmin
        kwargs["vmax"] = vmax
    im = ax.imshow(plotdata, **kwargs)
    ax.set_xlabel(f"${plane_labels[0]}$")
    ax.set_ylabel(f"${plane_labels[1]}$")
    ax.set_title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="6%", pad=0.45)
    cbar = fig.colorbar(im, cax=cax, orientation="horizontal")
    cbar.locator = ticker.MaxNLocator(nbins=5)
    cbar.formatter.set_powerlimits((0, 0))
    cbar.ax.xaxis.set_ticks_position("bottom")
    cbar.update_ticks()
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))


def render_1d(ax, x_axis, profile, xlabel, title):
    ax.plot(x_axis, profile, color="C0", lw=1.5)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
    ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
    ax.grid(True, alpha=0.3)
