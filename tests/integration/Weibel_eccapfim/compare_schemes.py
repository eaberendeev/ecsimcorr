#!/usr/bin/env python3
"""Compare matched ECCAPFIM/ECSIM-Corr Weibel calculations.

Each ROOT/{scheme}_dt{dt} contains configs, energy.txt, z-resolved field
slices, Performance/globalFunctions.time and runtime.json. The field diagnostic
must contain every unique physical Z node; averaging isolates kz=0.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from check_weibel import MODES, FIT_WINDOW, fit_growth, gamma_theory, read_plane


def analyze_case(directory):
    cfg = json.loads((directory / "system_config.json").read_text())
    pcfg = json.loads((directory / "particles_config.json").read_text())
    runtime = json.loads((directory / "runtime.json").read_text())
    rows = np.atleast_1d(np.genfromtxt(directory / "energy.txt", names=True))
    timing = np.atleast_1d(np.genfromtxt(directory / "Performance/globalFunctions.time", names=True))
    if any(not np.all(np.isfinite(rows[n])) for n in rows.dtype.names):
        raise ValueError(f"Nonfinite diagnostics in {directory}")
    if not np.all(np.isfinite(timing["Total"])) or np.any(timing["Total"] <= 0):
        raise ValueError(f"Invalid timing in {directory}")
    if not np.isclose(rows["Time"][-1], cfg["LastTimestep"] * cfg["Dt"]):
        raise ValueError(f"Incomplete calculation in {directory}")
    nx, ny, nz = (cfg[n] for n in ("NumCellsX", "NumCellsY", "NumCellsZ"))
    field_cfg = next(o for o in cfg["diagnostics"]["outputs"] if o["type"] == "fields_2d")
    snapshots = {}
    for path in (directory / "Fields/Diag2D").glob("FieldB_PlaneZ_*_*"):
        plane, seq = path.name.split("_")[-2:]
        if not seq.isascii() or not seq.isdecimal():
            continue
        snapshots.setdefault(int(seq), {})[int(plane)] = path
    times = np.array(sorted(snapshots)) * field_cfg["interval"] * cfg["Dt"]
    modes = {mode: [] for mode in MODES}
    for snapshot in sorted(snapshots):
        planes = snapshots[snapshot]
        # Beren3D stores g..g+Nz-1 as unique nodes; g is inferred from the data shape.
        if len(planes) != nz or sorted(planes) != list(range(1, nz + 1)):
            raise ValueError(f"Need all physical z planes, got {sorted(planes)} in {directory}")
        average = np.mean([read_plane(path, nx, ny) for path in planes.values()], axis=0)
        ft = np.fft.fft2(average) / (nx * ny)
        for mode in MODES:
            modes[mode].append(float(np.linalg.norm(ft[:2, mode[0], mode[1]])))
    fits = {}
    for mode, amplitude in modes.items():
        fit = fit_growth(times, np.array(amplitude))
        k = 2 * np.pi * np.hypot(mode[0] / (nx * cfg["Dx"]), mode[1] / (ny * cfg["Dy"]))
        theory = float(gamma_theory(k))
        fit.update(theory=theory, relative_error=fit["gamma"] / theory - 1)
        fits[str(mode)] = fit
    total = rows["Electrons"] + rows["energyFieldE"] + rows["energyFieldB"]
    anisotropy = 2 * rows["ElectronsZ"] / rows["ElectronsXY"]
    linear_timing = timing["Total"][(timing["Time"] * cfg["Dt"] >= FIT_WINDOW[0]) &
                                    (timing["Time"] * cfg["Dt"] <= FIT_WINDOW[1])]
    warm_timing = timing["Total"][timing["Time"] > 5]
    count = rows["ElectronsParticles"]
    if not np.all(count == count[0]):
        raise ValueError(f"Particle count changed in {directory}")
    metrics = {
        "scheme": cfg["Scheme"], "dt": cfg["Dt"], "grid": [nx, ny, nz],
        "particles": int(count[0]), "ppc": pcfg["particles"][0]["NumPartPerCell"],
        "threads": runtime["threads"], "final_time": float(rows["Time"][-1]),
        "initial_kinetic": float(rows["Electrons"][0]),
        "initial_kinetic_xy": float(rows["ElectronsXY"][0]),
        "initial_kinetic_z": float(rows["ElectronsZ"][0]),
        "wall_seconds": runtime["wall_seconds"],
        "mean_step_seconds": float(np.mean(warm_timing)),
        "median_step_seconds": float(np.median(warm_timing)),
        "linear_mean_step_seconds": float(np.mean(linear_timing)),
        "seconds_per_physical_time": float(np.mean(linear_timing) / cfg["Dt"]),
        "relative_energy_drift": float(np.max(np.abs(total - total[0])) / total[0]),
        "maximum_step_energy_residual": float(np.max(rows["energyConserve"])),
        "magnetic_peak": float(np.max(rows["energyFieldB"])),
        "magnetic_final": float(rows["energyFieldB"][-1]),
        "electric_final": float(rows["energyFieldE"][-1]),
        "anisotropy_initial": float(anisotropy[0]), "anisotropy_final": float(anisotropy[-1]),
        "fits": fits,
    }
    if "nonlinearIterations" in rows.dtype.names:
        metrics["mean_nonlinear_iterations"] = float(np.mean(rows["nonlinearIterations"][1:]))
    return metrics, rows, times, modes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--partial", action="store_true", help="Read only completed cases")
    args = parser.parse_args()
    data = {}
    for scheme in ("ecsim_corr", "eccapfim"):
        for dt in (1.2, 0.6):
            name = f"{scheme}_dt{dt:g}"
            directory = args.root / name
            if args.partial and not (directory / "runtime.json").exists():
                continue
            data[name] = analyze_case(directory)
    if not data:
        raise ValueError("No completed calculations")
    metrics = {name: values[0] for name, values in data.items()}
    particle_inputs = [json.loads((args.root / name / "particles_config.json").read_text()) for name in data]
    if any(value != particle_inputs[0] for value in particle_inputs[1:]):
        raise ValueError("Particle configurations differ between methods or time steps")
    shared_keys = ("NumCellsX", "NumCellsY", "NumCellsZ", "Dx", "Dy", "Dz", "Boundary_conditions")
    configs = [json.loads((args.root / name / "system_config.json").read_text()) for name in data]
    if any(any(c[key] != configs[0][key] for key in shared_keys) for c in configs):
        raise ValueError("Grid or boundaries differ between comparison cases")
    initial = [m["initial_kinetic"] for m in metrics.values()]
    if max(initial)-min(initial) > 1e-11 * max(initial):
        raise ValueError("Different initial kinetic energies: samples are not matched")
    directory = args.root / "Analysis"
    directory.mkdir(exist_ok=True)
    (directory / "comparison.json").write_text(json.dumps(metrics, indent=2, allow_nan=False) + "\n")
    convergence = {}
    for scheme in ("ecsim_corr", "eccapfim"):
        if not all(f"{scheme}_dt{dt:g}" in data for dt in (1.2, 0.6)):
            continue
        coarse, fine = (data[f"{scheme}_dt{dt:g}"] for dt in (1.2, 0.6))
        mask = (fine[2] >= FIT_WINDOW[0]) & (fine[2] <= FIT_WINDOW[1])
        temporal = {}
        for mode in MODES:
            reference = np.array(fine[3][mode])[mask]
            approximation = np.interp(fine[2][mask], coarse[2], coarse[3][mode])
            temporal[str(mode)] = {
                "relative_amplitude_difference": float(np.linalg.norm(approximation-reference) / np.linalg.norm(reference)),
                "gamma_difference": coarse[0]["fits"][str(mode)]["gamma"] - fine[0]["fits"][str(mode)]["gamma"],
            }
        convergence[scheme] = temporal
    (directory / "dt_comparison.json").write_text(json.dumps(convergence, indent=2, allow_nan=False) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    colors = {"ecsim_corr": "tab:orange", "eccapfim": "tab:blue"}
    for name, (m, rows, times, modes) in data.items():
        style = "-" if m["dt"] == 1.2 else "--"
        kw = {"color": colors[m["scheme"]], "ls": style, "label": f'{m["scheme"]}, dt={m["dt"]}'}
        axes[0, 0].semilogy(rows["Time"], rows["energyFieldB"] / m["initial_kinetic"], **kw)
        axes[0, 1].semilogy(times, modes[(1, 1)], **kw)
        total = rows["Electrons"] + rows["energyFieldE"] + rows["energyFieldB"]
        axes[1, 0].plot(rows["Time"], (total-total[0])/total[0], **kw)
        axes[1, 1].plot(rows["Time"], 2*rows["ElectronsZ"]/rows["ElectronsXY"], **kw)
    axes[0, 0].set_ylabel("Magnetic energy / initial kinetic energy")
    axes[0, 1].set_ylabel("Mode (1,1) amplitude, kz=0")
    axes[1, 0].set_ylabel("Relative total energy drift")
    axes[1, 1].set_ylabel("T_parallel / T_perp")
    for ax in axes.flat:
        ax.set_xlabel("omega_p t"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.savefig(directory / "comparison.png", dpi=160)
    plt.close(fig)
    for name, m in metrics.items():
        f = m["fits"]["(1, 1)"]
        print(f'{name}: {m["linear_mean_step_seconds"]:.4f} s/step, '
              f'gamma11={f["gamma"]:.5f} +/- {f["stderr"]:.5f}, '
              f'error={f["relative_error"]:.1%}, energy drift={m["relative_energy_drift"]:.2e}')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
