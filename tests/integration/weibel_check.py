#!/usr/bin/env python3
"""Automated Weibel-instability test for the new ecsim_corr runs (Sec. 3.4
of the catcode paper), new-code format only.

Draws the report figures and checks the numbers of a finished run:

  1. energies         : log W_B, W_E vs t, theory growth-rate line + PIC fit
  2. field maps       : Bx and Ez at t ~ 252 and at the last snapshot
  3. mode analysis    : per-mode growth rates vs theory (printed table)
  4. polarization     : B perp k, |Fx|=|Fy|, phase locked, Bz small (printed)
  5. energy balance   : W_B saturation, W_E behavior, conservation (printed)

The scheme is read from system_config.json ("Scheme"). Two tolerance sets:
  - ecsim_corr  : W_E must stay at thermal noise, fast growth, high saturation
  - no-corr     : W_E must GROW, slower growth, lower saturation (the known
                  signatures of the missing charge-correction step)

Usage:  python3 weibel_check.py [workdir] [--config gen_config.py]
                [--outdir DIR] [--quiet]

  workdir is optional: by default it is derived from the test's gen_config.py
  (DirName_Dx_.._np_.._Dt_.., the same naming as generate_config()), which is
  taken from --config, or from the adjacent Weibel_ecsim_corr / Weibel_ecsim
  test dirs, or from gen_config.py in the current directory. The run is looked
  up in the current directory and in the repo root (parent of tests/).

Output:
  - text report  -> <outdir>/test_weibel_report.txt (and stdout)
  - figures      -> <outdir>/test_energies.png, <outdir>/test_maps.png
  - exit code    : 0 = all checks passed, 1 = some check failed, 2 = bad input
"""
import os
import sys
import json
import argparse
import importlib.util
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import erfcx
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# expected physics (theory, eq. (25) of the paper; T_perp = 0.1 keV, T_par = 1 keV)
# ---------------------------------------------------------------------------
T_PERP_KEV = 0.1
T_PAR_KEV = 1.0
R = T_PAR_KEV / T_PERP_KEV
BETA_PERP = 511.0 / (2.0 * T_PERP_KEV)
EXP_GAMMA_MAX = 0.0236          # theory Gamma_max at k = 1.16
TOL_THEORY = 0.002

GAMMA_FIT_WIN = (72.0, 270.0)   # linear-stage window for the mode fit
GAMMA_WB_FIT_WIN = (60.0, 260.0)
PHASE_FIT_WIN = (126.0, 270.0)  # window for the phase-drift check
MAP_TIME = 252.0                # time of the "linear stage" maps

# tolerances (name, lower, upper) in the relevant units; per scheme:
#   gamma_rel  : Gamma(1,1) / theory value (0.0232)
#   gamma_wb   : Gamma from ln W_B slope / Gamma_max (0.0236)
#   we_factor  : max(W_E late) / mean(W_E early)
#   wb_max     : W_B saturation level
#   wb_end_min : W_B at the last snapshot (lower bound only)
TOL = {
    "ecsim_corr": dict(
        gamma_rel=(0.70, 1.15),
        gamma_wb=(0.50, 1.10),
        we_factor=(0.0, 4.0),
        wb_max=(0.15, 1.00),
        wb_end_min=0.05,
    ),
    "nocorr": dict(
        gamma_rel=(0.40, 1.15),
        gamma_wb=(0.25, 1.00),
        we_factor=(2.0, 500.0),
        wb_max=(0.03, 0.25),
        wb_end_min=0.02,
    ),
}
TOL_POL_RATIO = (0.5, 2.0)      # |Fx(1,1)|/|Fy(1,1)|
TOL_POL_PHASE = 0.4             # |arg(Fx/Fy) - pi| in rad
TOL_BZ_RATIO = 0.5              # Bz_rms / Bxy_rms
TOL_APERIODIC = 0.6             # phase drift of mode (1,1) over linear stage
TOL_ENERGY_CONS = 1e-10         # |dE| during the run
TOL_DN_RMS = 0.05               # density fluctuation rms

# ---------------------------------------------------------------------------
# theory: dispersion of the aperiodic Weibel instability (paper, eq. 25)
# ---------------------------------------------------------------------------
def gamma_theory(k):
    if k <= 0.0 or k >= np.sqrt(R - 1.0):
        return np.nan
    def F(g, kk):
        alpha = g * np.sqrt(BETA_PERP) / kk
        return g ** 2 + kk ** 2 - R + 1.0 + R * np.sqrt(np.pi) * alpha * erfcx(alpha)
    try:
        return brentq(F, 1e-9, 10.0, args=(k,))
    except ValueError:
        return np.nan


# ---------------------------------------------------------------------------
# data access (new-code format only)
# ---------------------------------------------------------------------------
FIELD_SNAPS = 12                # fields dumped every 12 steps
DT = 1.5


def read_field_plane(wd, base, snap):
    """FieldE/FieldB_PlaneZ_003_XXXX: header (2 floats) + 3 arrays of 33x33."""
    path = os.path.join(wd, "Fields", "Diag2D", f"{base}_PlaneZ_003_{snap:04d}")
    with open(path, "rb") as fh:
        nx, ny = np.fromfile(fh, dtype=np.float32, count=2)
        raw = np.fromfile(fh, dtype=np.float32)
    n = int(nx) * int(ny)
    return [raw[i * n:(i + 1) * n].reshape(int(nx), int(ny))[:30, :30]
            for i in range(3)]


def mode_amplitudes(wd, modes):
    snaps = sorted(int(f.replace("FieldB_PlaneZ_003_", ""))
                   for f in os.listdir(os.path.join(wd, "Fields", "Diag2D"))
                   if f.startswith("FieldB_PlaneZ_003_"))
    if not snaps:
        sys.exit("error: no FieldB_PlaneZ_003_* files")
    t = np.array(snaps) * FIELD_SNAPS * DT
    amp = {m: np.zeros(len(snaps)) for m in modes}
    for i, s in enumerate(snaps):
        bx, by, _ = read_field_plane(wd, "FieldB", s)
        fxb, fyb = np.fft.fft2(bx), np.fft.fft2(by)
        for (n, m) in modes:
            amp[(n, m)][i] = np.sqrt(abs(fxb[n, m]) ** 2 + abs(fyb[n, m]) ** 2)
    return t, amp, snaps


def fit_gamma(t, a, win):
    msk = (t >= win[0]) & (t <= win[1])
    if msk.sum() < 5:
        return np.nan, np.nan, t[msk]
    p, cov = np.polyfit(t[msk], np.log(a[msk]), 1, cov=True)
    return p[0], np.sqrt(cov[0, 0]), t[msk]


class Tee:
    """writes to a file and to the original stdout."""
    def __init__(self, path):
        self.file = open(path, "w")
        self.orig = sys.stdout
    def write(self, s):
        self.file.write(s)
        self.orig.write(s)
    def flush(self):
        self.file.flush()
        self.orig.flush()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("workdir", nargs="?", default=None,
                    help="run directory (default: from gen_config.py naming)")
    ap.add_argument("--config", default=None,
                    help="gen_config.py used for the default workdir name")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    wd = args.workdir
    if wd is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        cfg = args.config
        if cfg is None:
            for cand in ([os.path.join(script_dir, ts, "gen_config.py")
                          for ts in ("Weibel_ecsim_corr", "Weibel_ecsim")]
                         + ["gen_config.py"]):
                if os.path.isfile(cand):
                    cfg = cand
                    break
        if cfg is None:
            print("error: no gen_config.py found; pass --config "
                  "(or the workdir explicitly)", file=sys.stderr)
            return 2
        spec = importlib.util.spec_from_file_location("gen_config", cfg)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        name = (f"{mod.DirName}_Dx_{mod.Dx}_np_{mod.NumPartPerCell}"
                f"_Dt_{mod.Dt}")
        for base in (os.getcwd(),
                     os.path.abspath(os.path.join(script_dir, "..", ".."))):
            base = os.path.abspath(base)
            wd = os.path.join(base, name)
            if os.path.isdir(wd):
                break
        if not os.path.isdir(wd):
            print(f"error: no run '{name}'\n"
                  f"  (looked in {os.getcwd()} and "
                  f"{os.path.abspath(os.path.join(script_dir, '..', '..'))})\n"
                  f"  run the integration test first, or pass the workdir "
                  f"explicitly", file=sys.stderr)
            return 2
        print(f"[weibel_check] workdir (default name from {cfg}): {wd}")
    wd = os.path.abspath(wd)

    with open(os.path.join(wd, "system_config.json")) as f:
        sys_cfg = json.load(f)
    scheme = sys_cfg.get("Scheme", "")
    tol = TOL["ecsim_corr"] if "corr" in scheme.lower() else TOL["nocorr"]

    outdir = args.outdir or os.path.join(wd, "Analysis")
    os.makedirs(outdir, exist_ok=True)
    sys.stdout = Tee(os.path.join(outdir, "test_weibel_report.txt"))

    fails = passed = 0

    def report(name, value, ok, fmt="%.4g"):
        nonlocal fails, passed
        fails += int(not ok)
        passed += int(ok)
        tag = "PASS" if ok else "FAIL"
        print(f"  [{tag}] {name} = {fmt % value}")
        return ok

    d = np.loadtxt(os.path.join(wd, "energy.txt"), skiprows=1)
    t, w_e, w_b, cons = d[:, 0], d[:, 9], d[:, 10], d[:, 12]

    print("=" * 78)
    print("Weibel test report  (run: %s, scheme: %s)" % (wd, scheme))
    print("=" * 78)
    g_max = gamma_theory(1.16)
    g11_th = gamma_theory(2.0 * np.pi * np.sqrt(2.0) / 9.0)
    print(f"theory: Gamma_max = {g_max:.4f} (k=1.16), "
          f"Gamma(0.99) = {g11_th:.4f} (mode (1,1))")
    report("theory Gamma_max reproduced",
           g_max, abs(g_max - EXP_GAMMA_MAX) < TOL_THEORY)

    # ---------------- modes ----------------
    modes = [(1, 1), (1, 0), (0, 1), (2, 0), (0, 2), (2, 1), (1, 2), (0, 3)]
    t_s, amp, snaps = mode_amplitudes(wd, modes)
    k_unit = 2.0 * np.pi / 9.0

    g11, g11e, win = fit_gamma(t_s, amp[(1, 1)], GAMMA_FIT_WIN)
    ok = (not np.isnan(g11) and
          tol["gamma_rel"][0] * g11_th < g11 < tol["gamma_rel"][1] * g11_th)
    report(f"Gamma(mode (1,1)) = {g11:.4f} +/- {g11e:.4f} in t={win[0]:.0f}..{win[-1]:.0f}"
           f"  (theory {g11_th:.4f})", g11, ok)

    print("  mode analysis (Gamma vs theory):")
    for m in modes:
        if m == (1, 1):
            continue
        gm, ge, w = fit_gamma(t_s, amp[m], GAMMA_FIT_WIN)
        th = gamma_theory(np.sqrt(m[0] ** 2 + m[1] ** 2) * k_unit)
        if not np.isnan(gm) and np.isfinite(th):
            print(f"    k=({m[0]},{m[1]}) k={np.sqrt(m[0]**2+m[1]**2)*k_unit:.3f}: "
                  f"Gamma={gm:.4f} +/- {ge:.4f}  theory={th:.4f}")

    # ---------------- energies ----------------
    msk = (t >= GAMMA_WB_FIT_WIN[0]) & (t <= GAMMA_WB_FIT_WIN[1])
    p_b, cov_b = np.polyfit(t[msk], np.log(w_b[msk]), 1, cov=True)
    gw = p_b[0] / 2.0
    ok = tol["gamma_wb"][0] * g_max < gw < tol["gamma_wb"][1] * g_max
    report(f"Gamma from ln W_B slope (t={GAMMA_WB_FIT_WIN[0]:.0f}..{GAMMA_WB_FIT_WIN[1]:.0f}) = {gw:.4f}"
           f"  (theory max {g_max:.4f})", gw, ok)

    wb_max, i_max = w_b.max(), np.argmax(w_b)
    ok = tol["wb_max"][0] < wb_max < tol["wb_max"][1]
    report(f"W_B saturation = {wb_max:.3f} at t={t[i_max]:.0f}"
           f" (expected {tol['wb_max']})", wb_max, ok)
    ok = w_b[-1] > tol["wb_end_min"]
    report(f"W_B at the end = {w_b[-1]:.3f} (>= {tol['wb_end_min']})", w_b[-1], ok)

    we_early = w_e[(t > 50) & (t < 150)].mean()
    we_late = w_e[t > 300].max()
    fact = we_late / we_early
    ok = tol["we_factor"][0] < fact < tol["we_factor"][1]
    report(f"W_E late/early = {we_late:.2e} / {we_early:.2e} (factor {fact:.2f},"
           f" expected {tol['we_factor']})", fact, ok)

    decay = np.abs(cons).max()
    ok = decay < TOL_ENERGY_CONS
    report(f"energy conservation max|dE| = {decay:.1e} (< {TOL_ENERGY_CONS:.0e})", decay, ok)

    # ---------------- polarization & aperiodicity ----------------
    i_pk = np.argmin(np.abs(t_s - t[i_max]))
    bx, by, bz = read_field_plane(wd, "FieldB", snaps[i_pk])
    print(f"  polarization check at W_B peak (t={t_s[i_pk]:.0f}, snap {snaps[i_pk]}):")
    Fx, Fy = np.fft.fft2(bx)[1, 1], np.fft.fft2(by)[1, 1]
    r_amp, ph = abs(Fx / Fy), abs(np.angle(Fx / Fy))
    ok = (TOL_POL_RATIO[0] < r_amp < TOL_POL_RATIO[1] and
          abs(ph - np.pi) < TOL_POL_PHASE)
    report(f"  polarization mode (1,1): |Fx/Fy|={r_amp:.3f}, arg(Fx/Fy)={ph:.3f} rad"
           " (expect 1, pi)", r_amp, ok)
    bz_r = bz.std() / np.sqrt(bx.std() ** 2 + by.std() ** 2)
    ok = bz_r < TOL_BZ_RATIO
    report(f"  Bz_rms / Bxy_rms = {bz_r:.3f} (< {TOL_BZ_RATIO})", bz_r, ok)

    i0 = np.argmin(np.abs(t_s - PHASE_FIT_WIN[0]))
    i1 = np.argmin(np.abs(t_s - PHASE_FIT_WIN[1]))
    phas = np.array([np.angle(np.fft.fft2(read_field_plane(wd, "FieldB", s)[0])[1, 1])
                     for s in snaps[i0:i1 + 1]])
    drift = max(phas) - min(phas)
    ok = drift < TOL_APERIODIC
    report(f"  phase(1,1) drift over linear stage = {drift:.3f} rad"
           f" (< {TOL_APERIODIC}, aperiodic mode)", drift, ok)

    # ---------------- density fluctuation ----------------
    dens_snaps = sorted(int(f.replace("Density_PlaneZ_003_", ""))
                        for f in os.listdir(os.path.join(wd, "Particles",
                                                         "Electrons", "Diag2D"))
                        if f.startswith("Density_PlaneZ_003_"))
    if dens_snaps:
        s = dens_snaps[np.argmin(np.abs(np.array(dens_snaps) * 18.0 - MAP_TIME))]
        with open(os.path.join(wd, "Particles", "Electrons", "Diag2D",
                               f"Density_PlaneZ_003_{s:04d}"), "rb") as fh:
            nx, ny = np.fromfile(fh, dtype=np.float32, count=2)
            raw = np.fromfile(fh, dtype=np.float32)
        n = int(nx) * int(ny)
        rho = raw[:n].reshape(int(nx), int(ny))[:30, :30]
        dn = abs(-(rho + 1.0)).std()           # rho/n0 ~ -1; fluctuation rms
        ok = dn < TOL_DN_RMS
        report(f"density fluctuation rms = {dn * 100:.2f}% (< {TOL_DN_RMS * 100:.0f}%)",
               dn, ok)

    # ---------------- figures ----------------
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.semilogy(t, w_b, lw=2, color="tab:red", label=r"$W_B$")
    ax.semilogy(t, w_e, lw=2, color="tab:blue", ls="--", label=r"$W_E$")
    ax.semilogy(t[msk], np.exp(np.polyval(p_b, t[msk])), "k:", lw=2,
                label=r"fit, $\Gamma=%.4f$" % gw)
    tm = t[msk][len(t[msk]) // 2]
    a_ref = w_b[np.argmin(np.abs(t - tm))]
    tt = np.linspace(t[msk][0], 400.0, 60)
    ax.semilogy(tt, a_ref * np.exp(2.0 * g_max * (tt - tm)), "k--", lw=1.5,
                label=r"theory, $\Gamma_{max}=%.4f$" % g_max)
    ax.set_xlabel(r"$t\,\omega_p$"); ax.set_ylabel("Energy")
    ax.set_ylim(1e-5, 1.5); ax.legend(fontsize=10); ax.grid(alpha=0.3)
    ax.set_title(f"Weibel field energies, scheme = {scheme}")
    fig.savefig(os.path.join(outdir, "test_energies.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    def near(tt):
        i = np.argmin(np.abs(t_s - tt))
        return snaps[i], t_s[i]

    def plot_map(ax, data, title, cmap, vmin, vmax):
        im = ax.imshow(data.T, origin="lower", extent=[0, 9, 0, 9], cmap=cmap,
                       vmin=vmin, vmax=vmax)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(r"$x\,c/\omega_p$"); ax.set_ylabel(r"$y\,c/\omega_p$")
        return im

    s_lin, t_lin = near(MAP_TIME)
    s_end, t_end = near(t_s[-1])
    bxl, byl, bzl = read_field_plane(wd, "FieldB", s_lin)
    bxr, byr, bzr = read_field_plane(wd, "FieldB", s_end)
    exl, eyl, ezl = read_field_plane(wd, "FieldE", s_lin)
    exr, eyr, ezr = read_field_plane(wd, "FieldE", s_end)
    b_amp = max(np.abs(bxl).max(), np.abs(bxr).max(), np.abs(byl).max())
    e_amp = max(np.abs(ezl).max(), np.abs(ezr).max())

    fig = plt.figure(figsize=(11.0, 9.0))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.25)
    panels = [
        (fig.add_subplot(gs[0, 0]), bxl, r"$B_x$, $t\,\omega_p=%.0f$" % t_lin, b_amp),
        (fig.add_subplot(gs[0, 1]), bxr, r"$B_x$, $t\,\omega_p=%.0f$" % t_end, b_amp),
        (fig.add_subplot(gs[1, 0]), ezl, r"$E_z$, $t\,\omega_p=%.0f$" % t_lin, e_amp),
        (fig.add_subplot(gs[1, 1]), ezr, r"$E_z$, $t\,\omega_p=%.0f$" % t_end, e_amp),
    ]
    for ax, data, title, lim in panels:
        im = plot_map(ax, data, title, "RdBu_r", -lim, lim)
        fig.colorbar(im, ax=ax, shrink=0.85)
    fig.suptitle("Weibel: Bx and Ez maps (z=0.75), linear stage and end")
    fig.savefig(os.path.join(outdir, "test_maps.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # ---------------- summary ----------------
    print("=" * 78)
    status = "ALL CHECKS PASSED" if fails == 0 else f"{fails} CHECK(S) FAILED"
    print(f"Summary: {passed} passed, {fails} failed  ->  {status}")
    print("Figures: %s" % os.path.abspath(outdir))
    print("Report:  %s" % os.path.abspath(os.path.join(outdir, "test_weibel_report.txt")))
    sys.stdout.file.close()
    sys.stdout = sys.stdout.orig
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())