import numpy as np


def yee_grid_to_common(Ex, Ey):
    Ex_interp = np.empty_like(Ex)
    Ey_interp = np.empty_like(Ey)
    Ex_interp[0, :] = Ex[0, :]
    Ex_interp[1:, :] = 0.5 * (Ex[:-1, :] + Ex[1:, :])
    Ey_interp[:, 0] = Ey[:, 0]
    Ey_interp[:, 1:] = 0.5 * (Ey[:, :-1] + Ey[:, 1:])
    return Ex_interp, Ey_interp


def cartesian_to_cylindrical(Fx, Fy, Nx, Ny):
    x0 = (Nx - 1) / 2.0
    y0 = (Ny - 1) / 2.0
    ix = np.arange(Nx)[:, None] - x0
    iy = np.arange(Ny)[None, :] - y0
    r = np.sqrt(ix**2 + iy**2)
    r_safe = np.where(r > 1e-10, r, 1.0)
    cos_phi = np.where(r > 1e-10, ix / r_safe, 0.0)
    sin_phi = np.where(r > 1e-10, iy / r_safe, 0.0)
    Fr = Fx * cos_phi + Fy * sin_phi
    Fphi = -Fx * sin_phi + Fy * cos_phi
    return Fr, Fphi


def build_rmap(Nx, Ny):
    x0 = (Nx - 1) / 2.0
    y0 = (Ny - 1) / 2.0
    rmax = int(min(x0, y0))
    rmap = [[] for _ in range(rmax)]
    for x in range(Nx):
        for y in range(Ny):
            r2 = (x - x0)**2 + (y - y0)**2
            ri = int(np.sqrt(r2))
            if ri < rmax:
                rmap[ri].append((x, y))
    return [np.array(r).T if r else np.array([[], []], dtype=int) for r in rmap]


def phi_averaged(data, rmap):
    result = np.zeros(len(rmap))
    for i, r in enumerate(rmap):
        if r.shape[1] > 0:
            result[i] = data[r[0], r[1]].mean()
    return result


def compute_1d_profile(data, mode, dx1, dx2, plane_labels=("", ""), rmap=None):
    nx, ny = data.shape

    if mode == "phi_avg":
        if rmap is None:
            rmap = build_rmap(nx, ny)
        profile = phi_averaged(data, rmap)
        dr = min(dx1, dx2)
        x_axis = np.arange(len(profile)) * dr
        xlabel = "$r$"
        return x_axis, profile, xlabel, rmap

    if mode == "slice_0":
        profile = data[:, ny // 2]
        x_axis = np.arange(nx) * dx1
        xlabel = f"${plane_labels[0]}$"
        return x_axis, profile, xlabel, rmap

    if mode == "slice_1":
        profile = data[nx // 2, :]
        x_axis = np.arange(ny) * dx2
        xlabel = f"${plane_labels[1]}$"
        return x_axis, profile, xlabel, rmap

    if mode == "avg_0":
        profile = np.mean(data, axis=1)
        x_axis = np.arange(nx) * dx1
        xlabel = f"${plane_labels[0]}$"
        return x_axis, profile, xlabel, rmap

    if mode == "avg_1":
        profile = np.mean(data, axis=0)
        x_axis = np.arange(ny) * dx2
        xlabel = f"${plane_labels[1]}$"
        return x_axis, profile, xlabel, rmap

    raise ValueError(f"Unknown 1d_mode: {mode}")
