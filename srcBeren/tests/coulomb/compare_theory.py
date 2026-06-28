#!/usr/bin/env python3
"""Compare the Coulomb collision tests against analytic relaxation theory.

Theory:
  * single-species temperature isotropization  -> NRL/Trubnikov rate nu_T
  * two-species temperature equilibration       -> Glinskiy et al. Eq.17 (finite M)

Both rates are evaluated in CGS (Gaussian) units and converted to the code time
unit (1/omega_pe).  The simulation uses velocities in units of c, masses in m_e,
density in n0, so a code temperature T_code maps to T[erg] = T_code * m_e * c^2.

The Takizuka-Abe variance implemented in the code,
    <delta^2> = 2*pi (q1 q2)^2 n lnL dt / (m_red^2 u^3)   [Gaussian CGS],
is the standard TA77 form, so these analytic curves should match the points
within statistical noise.

Run coulomb.exe first (it writes <test>.csv and <test>_params.json), then:
    python3 compare_theory.py
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')   # headless: write PNGs without a display
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Physical constants (Gaussian CGS), matching srcBeren/utils/sgs.h
# ---------------------------------------------------------------------------
ME = 9.10938356e-28      # g
QE = 4.80320427e-10      # esu
C = 2.99792458e10        # cm/s
C2 = C * C
PI = np.pi


def omega_pe(n0):
    return np.sqrt(4.0 * PI * n0 * QE * QE / ME)


def rk4(deriv, y0, t_max, nsteps=8000):
    dt = t_max / nsteps
    y = np.array(y0, dtype=float)
    ts = np.empty(nsteps + 1)
    ys = np.empty((nsteps + 1, len(y0)))
    ts[0] = 0.0
    ys[0] = y
    t = 0.0
    for i in range(nsteps):
        k1 = deriv(t, y)
        k2 = deriv(t + dt / 2, y + dt / 2 * k1)
        k3 = deriv(t + dt / 2, y + dt / 2 * k2)
        k4 = deriv(t + dt, y + dt * k3)
        y = y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
        ts[i + 1] = t
        ys[i + 1] = y
    return ts, ys


def fit_exp_rate(t, y):
    """Fit y = y0 * exp(-k t) over the range where y stays above the noise floor.
    Returns the decay rate k (>0).  The lower cut (20% of the initial value)
    keeps the fit out of the statistical noise plateau."""
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)
    mask = y > 0.2 * y[0]
    if mask.sum() < 3:
        mask = y > 0
    slope, _ = np.polyfit(t[mask], np.log(y[mask]), 1)
    return -slope


def traj_rel_err(t_sim, sim, t_th, th):
    """RMS relative deviation of the simulated curve from the analytic ODE curve,
    sampled at the simulation times."""
    th_i = np.interp(t_sim, t_th, th)
    return np.sqrt(np.mean(((sim - th_i) / th_i) ** 2))


def cumtrapz(y, x):
    """Cumulative trapezoidal integral, c[0]=0 (no scipy dependency)."""
    d = np.diff(x)
    avg = 0.5 * (y[1:] + y[:-1])
    return np.concatenate(([0.0], np.cumsum(d * avg)))


# ---------------------------------------------------------------------------
# Isotropization theory (NRL temperature-isotropization rate)
#   dT_par/dt  = -2 nu_T (T_par - T_perp)
#   dT_perp/dt = +  nu_T (T_par - T_perp)
#   nu_T = 2 sqrt(pi) e^4 n lnL / (sqrt(m) T_par^{3/2})
#          * A^{-2} [ -3 + (A+3) f(A) ],   A = T_perp/T_par - 1
#   f(A) = atan(sqrt(A))/sqrt(A)   (A>0),  atanh(sqrt(-A))/sqrt(-A)  (A<0),  1 (A=0)
# ---------------------------------------------------------------------------
def make_iso_theory(p):
    n0 = p['n0']
    wpe = omega_pe(n0)
    lnL = p['coulomb_log']
    e = QE * abs(p['charge'])
    m = p['mass'] * ME
    n = p['n_density'] * n0

    def nu_T(Tpar_code, Tperp_code):
        Tpar = Tpar_code * ME * C2  # erg
        A = Tperp_code / Tpar_code - 1.0
        if abs(A) < 1e-8:
            br = 4.0 / 15.0
        elif A > 0:
            sa = np.sqrt(A)
            br = (-3.0 + (A + 3.0) * np.arctan(sa) / sa) / (A * A)
        else:
            sa = np.sqrt(-A)
            br = (-3.0 + (A + 3.0) * np.arctanh(sa) / sa) / (A * A)
        nu_phys = 2.0 * np.sqrt(PI) * e**4 * n * lnL / (np.sqrt(m) * Tpar**1.5) * br
        return nu_phys / wpe

    def deriv(t, y):
        Tpar, Tperp = y
        D = Tpar - Tperp
        nu = nu_T(Tpar, Tperp)
        return np.array([-2.0 * nu * D, +nu * D])

    return deriv, nu_T


# ---------------------------------------------------------------------------
# Equilibration theory, finite mass ratio (Glinskiy, Timofeev & Prikhodko,
# "Drift-kinetic PIC model ...", Eqs. 17-19).  Same TA77 collisions, same code
# units (t in 1/omega_pe); the only equilibration theory used in the tests.
#   t = tau_col * Int_{Te(0)/Tinf}^{Te(t)/Tinf}  [2/M + y(1-1/M)]^{3/2}/(1-y) dy
#   tau_col = (3 (2 pi)^{3/2} / (4 lnL)) (n_e c^3 / omega_pe^3) M (Tinf/m_e c^2)^{3/2}
#   Tinf = (Te(0)+Ti(0))/2,   M = m_i/m_e
# ---------------------------------------------------------------------------
def make_eq_theory_glinskij(p):
    n0 = p['n0']
    wpe = omega_pe(n0)
    lnL = p['coulomb_log']
    if p['mass1'] <= p['mass2']:
        e_idx, i_idx = 1, 2
    else:
        e_idx, i_idx = 2, 1
    M = p[f'mass{i_idx}'] / p[f'mass{e_idx}']
    Te0 = p[f'T{e_idx}_0']                 # code units (T / m_e c^2)
    Ti0 = p[f'T{i_idx}_0']
    ne = p[f'n{e_idx}_density'] * n0       # actual electron density [cm^-3]
    Tinf = (Te0 + Ti0) / 2.0               # code units (= T_inf / m_e c^2)

    # n_e c^3 / omega_pe(n_e)^3 = c^3 m_e^2 / (16 pi^2 n_e e^4) makes tau_col density-
    # independent of the (boosted) test density; the result below is the physical
    # tau_col [s], converted to code time (1/omega_pe(n0)) with omega_pe at n0.
    tau_col_phys = (3.0 * (2.0 * PI) ** 1.5 / (4.0 * lnL)) * \
        (C**3 * ME**2 / (16.0 * PI**2 * ne * QE**4)) * M * Tinf**1.5
    tau_col = tau_col_phys * wpe           # code units (wpe = omega_pe(n0))
    k_lin = 1.0 / (tau_col * (1.0 + 1.0 / M) ** 1.5)  # near-eq |Te-Ti| decay rate
    info = {'e_idx': e_idx, 'i_idx': i_idx, 'tau_col': tau_col, 'Tinf': Tinf, 'M': M, 'k_lin': k_lin}

    def curve(t_max_code):
        y0 = Te0 / Tinf
        eps = 1e-5
        y_end = 1.0 + eps if y0 > 1.0 else 1.0 - eps
        y = np.linspace(y0, y_end, 20000)
        g = (2.0 / M + y * (1.0 - 1.0 / M)) ** 1.5 / (1.0 - y)
        t = tau_col * cumtrapz(g, y)        # monotonically increasing, >= 0
        Te = y * Tinf
        Ti = (Te0 + Ti0) - Te
        mask = t <= t_max_code
        return t[mask], Te[mask], Ti[mask]

    return curve, info


# ---------------------------------------------------------------------------
def load(test):
    csv = test + '.csv'
    pjson = test + '_params.json'
    if not (os.path.exists(csv) and os.path.exists(pjson)):
        print(f"  Skipping {test}: missing {csv} or {pjson}")
        return None, None
    d = np.genfromtxt(csv, delimiter=',', names=True)
    with open(pjson) as f:
        p = json.load(f)
    return d, p


def plot_isotropization(test, ax_T, ax_aniso, ax_E):
    d, p = load(test)
    if d is None:
        return
    t = d['time']
    deriv, nu_T = make_iso_theory(p)
    ts, ys = rk4(deriv, [p['T_par0'], p['T_perp0']], t[-1])

    ax_T.plot(t, d['T_par'], 'rs', ms=3, label=r'sim $T_\parallel$')
    ax_T.plot(t, d['T_perp'], 'bo', ms=3, label=r'sim $T_\perp$')
    ax_T.plot(ts, ys[:, 0], 'r-', lw=1.5, label=r'theory $T_\parallel$')
    ax_T.plot(ts, ys[:, 1], 'b-', lw=1.5, label=r'theory $T_\perp$')
    ax_T.set_xlabel(r'$t\,\omega_{pe}$')
    ax_T.set_ylabel(r'$T$ [$m_e c^2$]')
    ax_T.set_title('Isotropization: temperatures')
    ax_T.legend(fontsize=8)
    ax_T.grid(True, alpha=0.3)

    D_sim = d['T_par'] - d['T_perp']
    D_th = ys[:, 0] - ys[:, 1]
    k_fit = fit_exp_rate(t, D_sim)
    k_th = 3.0 * nu_T(p['T_par0'], p['T_perp0'])
    err = traj_rel_err(t, d['T_par'], ts, ys[:, 0])
    ax_aniso.semilogy(t, D_sim, 'ks', ms=3, label='sim')
    ax_aniso.semilogy(ts, D_th, 'g-', lw=1.5, label='theory (ODE)')
    ax_aniso.set_xlabel(r'$t\,\omega_{pe}$')
    ax_aniso.set_ylabel(r'$T_\parallel - T_\perp$ [$m_e c^2$]')
    ax_aniso.set_title('Anisotropy decay')
    ax_aniso.legend(fontsize=8)
    ax_aniso.grid(True, alpha=0.3, which='both')
    ax_aniso.text(0.03, 0.05,
                  f'k_fit  = {k_fit:.3e}\n'
                  f'3 nu_T = {k_th:.3e}\n'
                  f'RMS err($T_\\parallel$) = {err * 100:.1f}%',
                  transform=ax_aniso.transAxes, fontsize=9, va='bottom',
                  bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

    ax_E.plot(t, d['Etotal_rel_err'], 'k-')
    ax_E.set_xlabel(r'$t\,\omega_{pe}$')
    ax_E.set_ylabel(r'$\Delta E / E_0$')
    ax_E.set_title('Energy conservation')
    ax_E.grid(True, alpha=0.3)
    ax_E.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    print(f"[isotropization] fitted rate k = {k_fit:.4e}, theory 3*nu_T(0) = {k_th:.4e}, "
          f"ratio = {k_fit / k_th:.3f}, RMS(T_par) = {err * 100:.1f}%")


def plot_equilibration(test, ax_T, ax_diff, ax_E):
    """Simulation vs Glinskiy Eq.17 (finite mass ratio M = m_i/m_e)."""
    d, p = load(test)
    if d is None:
        return
    t = d['time']
    curve_g, info_g = make_eq_theory_glinskij(p)
    tg, Te_g, Ti_g = curve_g(t[-1])
    e_idx = info_g['e_idx']

    def to_cols(Te, Ti):
        return (Te, Ti) if e_idx == 1 else (Ti, Te)

    g1, g2 = to_cols(Te_g, Ti_g)
    lab1, lab2 = (r'$T_e$', r'$T_i$') if e_idx == 1 else (r'$T_i$', r'$T_e$')
    sim_e = d['T1'] if e_idx == 1 else d['T2']

    ax_T.plot(t, d['T1'], 'rs', ms=3, label='sim ' + lab1)
    ax_T.plot(t, d['T2'], 'bo', ms=3, label='sim ' + lab2)
    ax_T.plot(tg, g1, 'r-', lw=1.6, label='Glinskij ' + lab1)
    ax_T.plot(tg, g2, 'b-', lw=1.6, label='Glinskij ' + lab2)
    ax_T.set_xlabel(r'$t\,\omega_{pe}$')
    ax_T.set_ylabel(r'$T$ [$m_e c^2$]')
    ax_T.set_title('Equilibration: sim vs Glinskij Eq.17 (M=%.0f)' % info_g['M'])
    ax_T.legend(fontsize=8)
    ax_T.grid(True, alpha=0.3)

    D_sim = np.abs(d['T1'] - d['T2'])
    err_g = traj_rel_err(t, sim_e, tg, Te_g)
    ax_diff.semilogy(t, D_sim, 'ks', ms=3, label='sim')
    ax_diff.semilogy(tg, np.abs(Te_g - Ti_g), 'm-', lw=1.8, label='Glinskij Eq.17')
    ax_diff.set_xlabel(r'$t\,\omega_{pe}$')
    ax_diff.set_ylabel(r'$|T_e - T_i|$ [$m_e c^2$]')
    ax_diff.set_title(r'$|T_e-T_i|$ decay  ($\tau_{col}=%.2e$)' % info_g['tau_col'])
    ax_diff.legend(fontsize=8)
    ax_diff.grid(True, alpha=0.3, which='both')
    ax_diff.text(0.03, 0.05,
                 f'RMS Glinskij = {err_g * 100:.1f}%',
                 transform=ax_diff.transAxes, fontsize=9, va='bottom',
                 bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

    ax_E.plot(t, d['Etotal_rel_err'], 'k-')
    ax_E.set_xlabel(r'$t\,\omega_{pe}$')
    ax_E.set_ylabel(r'$\Delta E / E_0$')
    ax_E.set_title('Energy conservation')
    ax_E.grid(True, alpha=0.3)
    ax_E.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    print(f"[equilibration]  M={info_g['M']:.0f}  tau_col(Glinskij)={info_g['tau_col']:.3e}  "
          f"RMS(Te) Glinskij={err_g * 100:.1f}%")


def plot_dt_convergence(tests, ax_T, ax_diff, ax_rms):
    """Timestep convergence at fixed density, N and t_max.  The TA77 finite-
    <delta^2> bias scales ~ dt, so RMS(Te) vs Glinskij theory should fall
    ~linearly with dt and the curves converge to the analytic solution as
    dt -> 0."""
    colors = ['C3', 'C0', 'C2', 'C4']
    results = []
    tg = Te_g = Ti_g = None
    e_idx = 1
    info_g = None
    for test, col in zip(tests, colors):
        d, p = load(test)
        if d is None:
            continue
        if tg is None:
            curve_g, info_g = make_eq_theory_glinskij(p)
            tg, Te_g, Ti_g = curve_g(d['time'][-1])
            e_idx = info_g['e_idx']
        t = d['time']
        sim_e = d['T1'] if e_idx == 1 else d['T2']
        sim_i = d['T2'] if e_idx == 1 else d['T1']
        dt = p['dt']
        err_g = traj_rel_err(t, sim_e, tg, Te_g)
        results.append((dt, err_g))
        ax_T.plot(t, sim_e, 'o', color=col, ms=3, label=f'dt={dt:g} $T_e$')
        ax_T.plot(t, sim_i, 's', color=col, ms=3, mfc='none')
        ax_diff.semilogy(t, np.abs(d['T1'] - d['T2']), 'o', color=col, ms=3, label=f'dt={dt:g}')

    if tg is None:
        ax_T.axis('off')
        ax_T.text(0.5, 0.5, 'missing dt-convergence data', ha='center', va='center', transform=ax_T.transAxes)
        return None

    ax_T.plot(tg, Te_g, 'k-', lw=1.4, label=r'Glinskij $T_e$')
    ax_T.plot(tg, Ti_g, 'g:', lw=1.8, label=r'Glinskij $T_i$')
    ax_T.set_xlabel(r'$t\,\omega_{pe}$')
    ax_T.set_ylabel(r'$T$ [$m_e c^2$]')
    ax_T.set_title('dt convergence: temperatures (M=%.0f)' % info_g['M'])
    ax_T.legend(fontsize=7, ncol=2)
    ax_T.grid(True, alpha=0.3)

    ax_diff.semilogy(tg, np.abs(Te_g - Ti_g), 'k-', lw=1.6, label='Glinskij')
    ax_diff.set_xlabel(r'$t\,\omega_{pe}$')
    ax_diff.set_ylabel(r'$|T_e - T_i|$ [$m_e c^2$]')
    ax_diff.set_title(r'$|T_e-T_i|$ decay vs $dt$')
    ax_diff.legend(fontsize=7)
    ax_diff.grid(True, alpha=0.3, which='both')

    # ---- RMS(Te) vs dt ----
    results.sort(key=lambda r: r[0])
    dts = np.array([r[0] for r in results])
    eg = np.array([r[1] * 100 for r in results])
    ax_rms.plot(dts, eg, 's-', color='g', label='vs Glinskij')
    ref = eg[-1] * dts / dts[-1]
    ax_rms.plot(dts, ref, 'r--', alpha=0.7, label=r'$\propto dt$ (TA77 bias)')
    ax_rms.set_xlabel(r'$dt$ [$1/\omega_{pe}$]')
    ax_rms.set_ylabel('RMS(Te) vs theory [%]')
    ax_rms.set_title('Convergence of RMS error vs dt')
    ax_rms.legend(fontsize=8)
    ax_rms.grid(True, alpha=0.3)
    ax_rms.set_ylim(bottom=0)
    ax_rms.set_xlim(left=0)

    print("[dt-convergence] dt      RMS(Te) vs Glinskij")
    for dt, e_g in results:
        print(f"[dt-convergence] {dt:5.2f}    {e_g * 100:6.2f}%")
    return results


def make_master_figure():
    """Single figure with every comparison."""
    fig = plt.figure(figsize=(16, 21))
    gs = fig.add_gridspec(5, 3, hspace=0.55, wspace=0.35, top=0.955, bottom=0.035, left=0.07, right=0.97)
    fig.suptitle('Coulomb collisions (Takizuka-Abe): all comparisons', fontsize=16, fontweight='bold')

    plot_isotropization('iso_electrons',
                        fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2]))
    plot_isotropization('iso_electrons_r5',
                        fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2]))
    plot_equilibration('eq_two_species',
                       fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1]), fig.add_subplot(gs[2, 2]))
    plot_equilibration('eq_two_species_real',
                       fig.add_subplot(gs[3, 0]), fig.add_subplot(gs[3, 1]), fig.add_subplot(gs[3, 2]))

    ax_txt = fig.add_subplot(gs[4, :])
    ax_txt.axis('off')
    ax_txt.text(0.0, 1.0,
                'Rows:\n'
                '1  isotropization  T_par = 2 T_perp   (TA77 Fig.10)\n'
                '2  isotropization  T_par = 5 T_perp\n'
                '3  equilibration   m_i/m_e = 100\n'
                '4  equilibration   m_i/m_e = 1836  (real)\n\n'
                'Theory:\n'
                '  isotropization -> NRL / TA77 Eq.20\n'
                '  equilibration  -> Glinskij Eq.17  (finite M)\n\n'
                r'T_par = m<v_x^2>,  T_perp = (m/2)(<v_y^2>+<v_z^2>)' '\n'
                'units: T in m_e c^2, t in 1/omega_pe\n'
                'energy conserved to ~1e-15 in every run',
                transform=ax_txt.transAxes, fontsize=11, va='top', ha='left', family='monospace',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

    fig.savefig('coulomb_all_comparisons.png', dpi=140, bbox_inches='tight')
    print("Saved: coulomb_all_comparisons.png")


if __name__ == '__main__':
    plt.rcParams.update({'figure.dpi': 120, 'axes.titlesize': 11, 'axes.labelsize': 10,
                         'legend.fontsize': 8, 'xtick.labelsize': 9, 'ytick.labelsize': 9})

    # ---- main figure: all four tests, simulation vs analytic theory ----
    fig = plt.figure(figsize=(16, 17))
    gs = fig.add_gridspec(4, 3, hspace=0.5, wspace=0.35, top=0.95, bottom=0.04, left=0.07, right=0.97)
    fig.suptitle('Coulomb collisions: simulation vs analytic theory', fontsize=15, fontweight='bold')

    plot_isotropization('iso_electrons',
                        fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2]))
    plot_isotropization('iso_electrons_r5',
                        fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2]))
    plot_equilibration('eq_two_species',
                       fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1]), fig.add_subplot(gs[2, 2]))
    plot_equilibration('eq_two_species_real',
                       fig.add_subplot(gs[3, 0]), fig.add_subplot(gs[3, 1]), fig.add_subplot(gs[3, 2]))

    plt.savefig('coulomb_theory_comparison.png', dpi=150, bbox_inches='tight')
    print("Saved: coulomb_theory_comparison.png")

    # ---- timestep convergence study (dt = 3 vs dt = 1) ----
    figd = plt.figure(figsize=(16, 4.7))
    gsd = figd.add_gridspec(1, 3, wspace=0.32, top=0.84, bottom=0.17, left=0.06, right=0.98)
    figd.suptitle('Timestep convergence at fixed density, N and t_max (M=100, e-i equilibration)',
                  fontsize=14, fontweight='bold')
    plot_dt_convergence(['eq_two_species', 'eq_two_species_dt1'],
                        figd.add_subplot(gsd[0, 0]), figd.add_subplot(gsd[0, 1]), figd.add_subplot(gsd[0, 2]))
    figd.savefig('coulomb_dt_convergence.png', dpi=150, bbox_inches='tight')
    print("Saved: coulomb_dt_convergence.png")

    # ---- single combined figure with every comparison ----
    make_master_figure()
