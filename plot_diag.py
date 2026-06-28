#!/usr/bin/env python3
import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

from diagplot import (
    DataCache,
    compute_1d_profile,
    auto_cmap,
    auto_title,
    render_2d,
    render_1d,
    CMAP_FIELD,
    CMAP_DENS,
)

# ============================================================================
#  CONFIGURATION — edit this section
# ============================================================================
CONFIG = {
    "workdir": "./",
    "plane": "Z",
    "slice_num": "003",
    "t_auto": True,
    "t0": 0,
    "tmax": 9999,
    "tstep": 1,
    "figsize": None,
    "dpi": 150,
    "outdir": "./Anime/2D",

    "nrows": 3,
    "ncols": 3,

    # Each panel: row, col, quantity, species (optional), type, vmin/vmax, ...
    #
    # type: "2d"    — 2D color plot
    #       "1d"    — 1D profile only
    #       "2d+1d" — 2D color + 1D profile side-by-side
    #
    # 1d_mode (for "1d" and "2d+1d"):
    #   "phi_avg" — azimuthal average (Z-plane only), default for Z
    #   "slice_0" — central slice, profile along axis 0, default for X/Y
    #   "slice_1" — central slice, profile along axis 1
    #   "avg_0"   — average over axis 1, profile along axis 0
    #   "avg_1"   — average over axis 0, profile along axis 1
    #
    # quantity (fields, no species needed):
    #   Ex, Ey, Ez, Er, Ephi, Bx, By, Bz, Br, Bphi
    #
    # quantity (per-species, set "species"):
    #   Density, Jx, Jy, Jz, Jr, Jphi, Prr, Ppp, Pzz, Prp, Prz, Pzp
    #
    # cmap: "field" (diverging) or "dens" (sequential); auto if omitted
    #       field cmap auto-symmetrizes vmin/vmax around 0
    # abs: true to plot |data| (Density always abs)
    # vmin, vmax: color limits; omit for auto-scale
    # title: override auto-generated title

    "panels": [
        {"row": 0, "col": 0, "quantity": "Er",   "type": "2d+1d", "vmin": -0.03, "vmax": 0.03},
        {"row": 0, "col": 1, "quantity": "Ephi", "type": "2d+1d", "vmin": -0.03, "vmax": 0.03},
        {"row": 0, "col": 2, "quantity": "Ez",   "type": "2d",    "vmin": -0.03, "vmax": 0.03},

        {"row": 1, "col": 0, "quantity": "Density", "species": "Electrons", "type": "2d+1d"},
        {"row": 1, "col": 1, "quantity": "Jphi",    "species": "Electrons", "type": "2d+1d"},
        {"row": 1, "col": 2, "quantity": "Bz",      "type": "2d+1d", "vmin": 0.17, "vmax": 0.23},

        {"row": 2, "col": 0, "quantity": "Density", "species": "Ions", "type": "2d+1d"},
        {"row": 2, "col": 1, "quantity": "Jphi",    "species": "Ions", "type": "2d+1d"},
        {"row": 2, "col": 2, "quantity": "Prr",     "species": "Ions", "type": "2d+1d"},
    ],
}

# ============================================================================

VALID_PLANES = {"X", "Y", "Z"}
VALID_TYPES = {"2d", "1d", "2d+1d"}
VALID_MODES = {"phi_avg", "slice_0", "slice_1", "avg_0", "avg_1"}
VALID_CMAPS = {"field", "dens"}
FIELD_QUANTITIES = {"Ex", "Ey", "Ez", "Er", "Ephi", "Bx", "By", "Bz", "Br", "Bphi"}
SPECIES_QUANTITIES = {"Density", "Jx", "Jy", "Jz", "Jr", "Jphi",
                      "Prr", "Ppp", "Pzz", "Prp", "Prz", "Pzp"}
ALL_QUANTITIES = FIELD_QUANTITIES | SPECIES_QUANTITIES
CYLINDRICAL_QUANTITIES = {"Er", "Ephi", "Br", "Bphi", "Jr", "Jphi"}

W_2D = 3
W_1D = 2

MODE_TITLES = {
    "phi_avg": lambda t: r"$\langle$" + t + r"$\rangle_\varphi$",
    "slice_0": lambda t: t + " (mid)",
    "slice_1": lambda t: t + " (mid)",
    "avg_0":   lambda t: r"$\langle$" + t + r"$\rangle$",
    "avg_1":   lambda t: r"$\langle$" + t + r"$\rangle$",
}


def get_plane_info(plane, sys_cfg):
    nx = int(sys_cfg["NumCellsX"])
    ny = int(sys_cfg["NumCellsY"])
    nz = int(sys_cfg["NumCellsZ"])
    dx = float(sys_cfg["Dx"])
    dy = float(sys_cfg["Dy"])
    dz = float(sys_cfg["Dz"])
    if plane == "Z":
        return nx, ny, dx, dy, ("x", "y"), [0, nx * dx, 0, ny * dy]
    if plane == "Y":
        return nx, nz, dx, dz, ("x", "z"), [0, nx * dx, 0, nz * dz]
    if plane == "X":
        return ny, nz, dy, dz, ("y", "z"), [0, ny * dy, 0, nz * dz]
    raise ValueError(f"Unknown plane: {plane}")


def default_1d_mode(plane):
    return "phi_avg" if plane == "Z" else "slice_0"


def validate_config(cfg, species_list):
    errors = []
    warnings = []

    plane = cfg.get("plane", "")
    if plane not in VALID_PLANES:
        errors.append(f'plane="{plane}" invalid, must be one of {sorted(VALID_PLANES)}')

    nrows = cfg.get("nrows", 0)
    ncols = cfg.get("ncols", 0)
    panels = cfg.get("panels", [])

    if not panels:
        errors.append("panels list is empty")

    for i, p in enumerate(panels):
        tag = f"panels[{i}]"

        for key in ("row", "col", "quantity"):
            if key not in p:
                errors.append(f'{tag}: missing required key "{key}"')

        row = p.get("row", 0)
        col = p.get("col", 0)
        if row >= nrows:
            errors.append(f'{tag}: row={row} >= nrows={nrows}')
        if col >= ncols:
            errors.append(f'{tag}: col={col} >= ncols={ncols}')

        q = p.get("quantity", "")
        if q and q not in ALL_QUANTITIES:
            errors.append(
                f'{tag}: quantity="{q}" unknown. '
                f'Valid: {sorted(ALL_QUANTITIES)}'
            )

        species = p.get("species")
        if q in SPECIES_QUANTITIES and not species:
            errors.append(f'{tag}: quantity="{q}" requires "species"')
        if species and species not in species_list:
            errors.append(
                f'{tag}: species="{species}" not found. '
                f'Available: {species_list}'
            )
        if q in FIELD_QUANTITIES and species:
            warnings.append(f'{tag}: quantity="{q}" is a field, "species" ignored')

        ptype = p.get("type", "2d")
        if ptype not in VALID_TYPES:
            errors.append(
                f'{tag}: type="{ptype}" invalid, '
                f'must be one of {sorted(VALID_TYPES)}'
            )

        mode = p.get("1d_mode")
        if mode and mode not in VALID_MODES:
            errors.append(
                f'{tag}: 1d_mode="{mode}" invalid, '
                f'must be one of {sorted(VALID_MODES)}'
            )
        eff_mode = mode or default_1d_mode(plane)
        if eff_mode == "phi_avg" and plane != "Z" and ptype in ("1d", "2d+1d"):
            errors.append(
                f'{tag}: 1d_mode="phi_avg" only works on plane="Z", '
                f'use slice_0/slice_1/avg_0/avg_1'
            )

        if q in CYLINDRICAL_QUANTITIES and plane != "Z":
            warnings.append(
                f'{tag}: quantity="{q}" is cylindrical, '
                f'meaningful only on plane="Z"'
            )

        cmap = p.get("cmap")
        if cmap and cmap not in VALID_CMAPS:
            errors.append(
                f'{tag}: cmap="{cmap}" invalid, '
                f'must be one of {sorted(VALID_CMAPS)}'
            )

    for w in warnings:
        print(f"  WARNING: {w}")

    if errors:
        print("\nConfig errors:")
        for e in errors:
            print(f"  ERROR: {e}")
        raise SystemExit(1)


def main():
    cfg = CONFIG
    workdir = cfg["workdir"]

    sys_path = os.path.join(workdir, "system_config.json")
    part_path = os.path.join(workdir, "particles_config.json")
    for path in (sys_path, part_path):
        if not os.path.isfile(path):
            print(f"ERROR: config not found: {path}")
            raise SystemExit(1)

    with open(sys_path, "r", encoding="utf-8") as f:
        sys_cfg = json.load(f)
    with open(part_path, "r", encoding="utf-8") as f:
        part_cfg = json.load(f)

    species_list = [p["Name"] for p in part_cfg.get("particles", [])]

    validate_config(cfg, species_list)

    grid_nx, grid_ny, dx1, dx2, plane_labels, extent = get_plane_info(
        cfg["plane"], sys_cfg,
    )

    tdelay = (int(sys_cfg["diagnostics"]["TimeStepDelayDiag2D"])
              * float(sys_cfg["Dt"]))

    os.makedirs(cfg["outdir"], exist_ok=True)

    cache = DataCache(workdir, cfg["plane"], cfg["slice_num"], species_list)

    nrows = cfg["nrows"]
    ncols = cfg["ncols"]
    panels = cfg["panels"]
    def_mode = default_1d_mode(cfg["plane"])

    col_ratios = []
    col_offsets = []
    col_wide = []
    gs_offset = 0
    for c in range(ncols):
        col_offsets.append(gs_offset)
        col_panels = [p for p in panels if p["col"] == c]
        is_wide = any(p.get("type", "2d") == "2d+1d" for p in col_panels)
        col_wide.append(is_wide)
        if is_wide:
            col_ratios.extend([W_2D, W_1D])
            gs_offset += 2
        else:
            col_ratios.append(W_2D)
            gs_offset += 1

    if cfg.get("figsize"):
        figsize = cfg["figsize"]
    else:
        figsize = (sum(col_ratios) * 1.8 + 2, nrows * 5.5 + 1.5)

    if cfg.get("t_auto", False):
        pattern = f"FieldE_Plane{cfg['plane']}_{cfg['slice_num']}_"
        diag_dir = os.path.join(workdir, "Fields", "Diag2D")
        if not os.path.isdir(diag_dir):
            print(f"ERROR: directory not found: {diag_dir}")
            raise SystemExit(1)
        timesteps = sorted(
            int(f.replace(pattern, ""))
            for f in os.listdir(diag_dir)
            if f.startswith(pattern)
        )
        if not timesteps:
            print(f'ERROR: no files matching "{pattern}*" in {diag_dir}')
            avail = {f.split("_Plane")[1].split("_")[0]
                     for f in os.listdir(diag_dir)
                     if f.startswith("FieldE_Plane")}
            if avail:
                slices = {f.split("_")[1] for f in os.listdir(diag_dir)
                          if f.startswith(f"FieldE_Plane{cfg['plane']}_")}
                print(f'  Available planes: {sorted(avail)}')
                if slices:
                    print(f'  Available slices for plane {cfg["plane"]}: {sorted(slices)}')
            raise SystemExit(1)
        print(f"Auto-detected {len(timesteps)} timesteps: "
              f"{timesteps[0]}..{timesteps[-1]}")
    else:
        timesteps = list(range(cfg["t0"], cfg["tmax"], cfg["tstep"]))

    print(f"Figure size: {figsize[0]:.1f} x {figsize[1]:.1f} inches")

    for t in timesteps:
        ts = str(t).zfill(4)
        cache.clear()
        print(f"Timestep {ts}")

        fig = plt.figure(figsize=figsize)
        total_gs_cols = len(col_ratios)
        gs = gridspec.GridSpec(
            nrows + 1, total_gs_cols,
            height_ratios=[0.06] + [1] * nrows,
            width_ratios=col_ratios,
            hspace=0.4, wspace=0.35,
        )

        rmap = None
        rendered = 0

        for panel in panels:
            row = panel["row"]
            c = panel["col"]
            quantity = panel["quantity"]
            species = panel.get("species", None)
            ptype = panel.get("type", "2d")
            vmin = panel.get("vmin", None)
            vmax = panel.get("vmax", None)
            use_abs = panel.get("abs", False)
            title = panel.get("title", None) or auto_title(quantity, species)
            cmap_name = panel.get("cmap", None)
            mode_1d = panel.get("1d_mode", def_mode)

            try:
                data = cache.get(quantity, species, ts)
            except FileNotFoundError:
                print(f"  WARNING: file not found for {quantity}"
                      f"{'/' + species if species else ''}, skipping")
                continue
            rendered += 1

            if use_abs or quantity == "Density":
                data = np.abs(data)

            if cmap_name == "field":
                cmap = CMAP_FIELD
            elif cmap_name == "dens":
                cmap = CMAP_DENS
            else:
                cmap = auto_cmap(quantity)

            if cmap is CMAP_FIELD:
                if vmin is None and vmax is None:
                    abs_max = float(np.nanmax(np.abs(data)))
                    vmin, vmax = -abs_max, abs_max
                elif vmin is None:
                    vmin = -abs(vmax)
                elif vmax is None:
                    vmax = abs(vmin)

            gs_col = col_offsets[c]

            if ptype == "2d":
                ax = fig.add_subplot(gs[row + 1, gs_col])
                render_2d(ax, fig, data, extent, cmap, vmin, vmax,
                          title, plane_labels)

            elif ptype == "1d":
                ax = fig.add_subplot(gs[row + 1, gs_col])
                x_axis, profile, xlabel, rmap = compute_1d_profile(
                    data, mode_1d, dx1, dx2, plane_labels, rmap,
                )
                title_1d = MODE_TITLES[mode_1d](title)
                render_1d(ax, x_axis, profile, xlabel, title_1d)

            elif ptype == "2d+1d":
                ax2d = fig.add_subplot(gs[row + 1, gs_col])
                ax1d = fig.add_subplot(gs[row + 1, gs_col + 1])
                render_2d(ax2d, fig, data, extent, cmap, vmin, vmax,
                          title, plane_labels)
                x_axis, profile, xlabel, rmap = compute_1d_profile(
                    data, mode_1d, dx1, dx2, plane_labels, rmap,
                )
                title_1d = MODE_TITLES[mode_1d](title)
                render_1d(ax1d, x_axis, profile, xlabel, title_1d)

        if rendered == 0:
            plt.close(fig)
            continue

        phys_time = round(t * tdelay, 2)
        fig.text(0.5, 0.98, rf"$t \cdot \omega_p = {phys_time}$",
                 fontsize=16, ha="center", va="top",
                 bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

        outname = os.path.join(cfg["outdir"],
                               f"Diag_{cfg['plane']}_{ts}.png")
        fig.savefig(outname, format="png", dpi=cfg["dpi"],
                    bbox_inches="tight")
        plt.close(fig)
        print(f"  -> {outname}")


if __name__ == "__main__":
    main()
