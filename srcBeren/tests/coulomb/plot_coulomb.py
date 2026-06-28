#!/usr/bin/env python3
"""Raw diagnostics for the Coulomb collision tests.
Run coulomb first to generate the CSV files, then run this script.
For comparison against analytic theory use compare_theory.py instead."""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')   # headless: write PNGs without a display
import matplotlib.pyplot as plt


def plot_isotropization(filename, title, ax_temp, ax_energy):
    if not os.path.exists(filename):
        print(f"  Skipping {filename} (not found)")
        return
    d = np.genfromtxt(filename, delimiter=',', names=True)
    t = d['time']
    ax_temp.plot(t, d['T_par'], 'r-', label=r'$T_\parallel$', linewidth=1.5)
    ax_temp.plot(t, d['T_perp'], 'b-', label=r'$T_\perp$', linewidth=1.5)
    ax_temp.plot(t, d['T_total'], 'k--', label=r'$T_{tot}$', linewidth=1.0, alpha=0.7)
    ax_temp.set_xlabel(r'$t\,\omega_{pe}$')
    ax_temp.set_ylabel(r'$T$ [$m_e c^2$]')
    ax_temp.set_title(title + ': temperature isotropization')
    ax_temp.legend()
    ax_temp.grid(True, alpha=0.3)

    ax_energy.plot(t, d['Etotal_rel_err'], 'k-', linewidth=1)
    ax_energy.set_xlabel(r'$t\,\omega_{pe}$')
    ax_energy.set_ylabel(r'$\Delta E / E_0$')
    ax_energy.set_title(title + ': energy conservation')
    ax_energy.grid(True, alpha=0.3)
    ax_energy.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))


def plot_equilibration(filename, title, ax_temp, ax_energy):
    if not os.path.exists(filename):
        print(f"  Skipping {filename} (not found)")
        return
    d = np.genfromtxt(filename, delimiter=',', names=True)
    t = d['time']
    ax_temp.plot(t, d['T1'], 'r-', label=r'$T_1$ (species 1)', linewidth=1.5)
    ax_temp.plot(t, d['T2'], 'b-', label=r'$T_2$ (species 2)', linewidth=1.5)
    ax_temp.set_xlabel(r'$t\,\omega_{pe}$')
    ax_temp.set_ylabel(r'$T$ [$m_e c^2$]')
    ax_temp.set_title(title + ': temperature equilibration')
    ax_temp.legend()
    ax_temp.grid(True, alpha=0.3)

    ax_energy.plot(t, d['Etotal_rel_err'], 'k-', linewidth=1)
    ax_energy.set_xlabel(r'$t\,\omega_{pe}$')
    ax_energy.set_ylabel(r'$\Delta E / E_0$')
    ax_energy.set_title(title + ': energy conservation')
    ax_energy.grid(True, alpha=0.3)
    ax_energy.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))


if __name__ == '__main__':
    plt.rcParams.update({'figure.dpi': 120, 'axes.titlesize': 12, 'axes.labelsize': 10,
                         'legend.fontsize': 9, 'xtick.labelsize': 9, 'ytick.labelsize': 9})

    fig = plt.figure(figsize=(16, 8))
    gs = fig.add_gridspec(2, 2, hspace=0.4, wspace=0.3, top=0.9, bottom=0.1, left=0.08, right=0.96)
    fig.suptitle('Coulomb Collision Tests (Takizuka-Abe) - raw diagnostics', fontsize=15, fontweight='bold')

    plot_isotropization('iso_electrons.csv', 'e-e', fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]))
    plot_equilibration('eq_two_species.csv', 'e-i', fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]))

    plt.savefig('coulomb_tests.png', dpi=150, bbox_inches='tight')
    print("Saved: coulomb_tests.png")
