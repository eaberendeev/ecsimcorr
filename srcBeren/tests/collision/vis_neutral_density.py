#!/usr/bin/env python3
"""Visualise neutral density decay for different collision channels.

The script reads:
- simulation output (time in 1/omega_p and neutral counts),
- a JSON test config with densities/energies/masses,
- an optional choice for neutral velocity distribution (Maxwellian or mono beam).

It then reconstructs the theoretical neutral density decay n(t)/n0 = exp(-nu/omega_p * t)
based on the selected collision option and compares it to the simulation.

Example run (monoenergetic neutrals, charged Maxwellian):
    python vis_neutral_density.py test_charge_exchange_m_1836_dt_300.000000.txt \\
        --config test_config.json --neutral-distribution mono --samples 6000 --output plot.png
"""
from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Optional

import matplotlib.pyplot as plt
import numpy as np

# Physical constants (CGS units where applicable) matching the C++ codebase.
ME = 9.10938356e-28        # Electron mass, g
QE = 4.80320427e-10        # Elementary charge, statC
C_LIGHT = 2.99792458e10    # Speed of light, cm/s
MC2_KEV = 511.0            # Electron rest energy, keV
EV_TO_ERG = 1.602176634e-12
AMU_H = 1.00783

# Mass ratios (in electron masses)
M_E_RATIO = 1.0
M_I_RATIO = 1836.0

# Cross-section model constants (mirroring cross_section.h)
A = 1.53e-25
B = 0.6
C = 0.56

MC2_EV = 511000.0
EV_TO_MC2 = (1.0 / MC2_EV)
HYDROGEN_IONIZATION_EV = 13.5984346

P_THRESHOLD = HYDROGEN_IONIZATION_EV * EV_TO_MC2

# P_THRESHOLD = 2.6614481409001956e-05

A1 = 3.2345
A2 = 235.88
A3 = 0.038371
A4 = 3.8068e-6
A5 = 1.1832e-10
A6 = 2.3713

B1 = 12.899
B2 = 61.897
B3 = 9.2731e3
B4 = 4.9749e-4
B5 = 3.9890e-2
B6 = -1.5900
B7 = 3.1834
B8 = -3.7154


@dataclass
class SimulationMeta:
    species: str
    temperature_eV: float
    base_density: float
    reaction: str
    mass_ratio: float


@dataclass
class TestSettings:
    name: str
    n0: float
    neutral_relative_density: float
    num_particles: int
    num_steps: int
    neutrals_energy_kev: float
    particles_energy_kev: float
    charged_mass: float
    neutral_mass: float
    dt: float
    collision_scheme: str
    reaction: str

    @property
    def charged_density(self) -> float:
        """Density of charged species in cm^-3 (neutral_relative_density is given in units of n0)."""
        return self.n0 

    @property
    def charged_species(self) -> str:
        if self.reaction == 'electron_ionization':
            return 'electron'
        return 'proton'


def sigma_e_cm2(energy_mc2: np.ndarray) -> np.ndarray:
    """Electron-impact ionisation cross section in cm^2."""
    energy = np.asarray(energy_mc2, dtype=float)
    result = np.zeros_like(energy)
    mask = energy > P_THRESHOLD
    if np.any(mask):
        ratio = energy[mask] / P_THRESHOLD
        s1 = B * np.exp(-C * (ratio - 1.0))
        s2 = A * np.log(ratio) / (energy[mask] * P_THRESHOLD)
        result[mask] = s2 * (1.0 - s1)
    return result


def sigma_p_cm2(energy_mc2: np.ndarray) -> np.ndarray:
    """Proton-impact ionisation cross section in cm^2."""
    energy = np.asarray(energy_mc2, dtype=float)
    result = np.zeros_like(energy)
    mask = energy > P_THRESHOLD
    if np.any(mask):
        e_kev_per_amu = energy[mask] * MC2_KEV / AMU_H
        left = np.exp(-B2 / e_kev_per_amu) * np.log1p(B3 * e_kev_per_amu) / e_kev_per_amu
        right = B4 * np.exp(-B5 * e_kev_per_amu) / (
            np.power(e_kev_per_amu, B6) + B7 * np.power(e_kev_per_amu, B8)
        )
        result[mask] = 1e-16 * B1 * (left + right)
    return result


def sigma_cx_cm2(energy_mc2: np.ndarray) -> np.ndarray:
    """Charge-exchange cross section in cm^2."""
    energy = np.asarray(energy_mc2, dtype=float)
    e_kev_per_amu = energy * MC2_KEV / AMU_H
    numerator = A1 * np.log(A2 / e_kev_per_amu + A6)
    denominator = 1.0 + A3 * e_kev_per_amu + A4 * np.power(e_kev_per_amu, 3.5) + A5 * np.power(e_kev_per_amu, 5.4)
    return 1e-16 * numerator / denominator


def maxwellian_sigma_v(temperature_eV: float, mass_ratio: float, sigma_func, nodes: int = 80) -> float:
    """Maxwellian average of sigma*v using Gauss-Laguerre quadrature."""
    if temperature_eV <= 0.0:
        return 0.0
    theta = temperature_eV * EV_TO_ERG / (ME * C_LIGHT ** 2)
    if theta <= 0.0:
        return 0.0
    x, w = np.polynomial.laguerre.laggauss(nodes)
    energies = theta * x
    sigma_vals = sigma_func(energies)
    integral = np.sum(w * sigma_vals * x)
    prefactor = C_LIGHT * (2.0 / math.sqrt(math.pi)) * math.sqrt(2.0 * theta / mass_ratio)
    return prefactor * integral


def extract_dt_from_filename(path: Path) -> Optional[float]:
    match = re.search(r"_dt_([0-9]+(?:\.[0-9]+)?)", path.stem)
    if not match:
        return None
    return float(match.group(1))


def detect_reaction(options: dict[str, bool]) -> str:
    active = [name for name, enabled in options.items() if enabled]
    if len(active) != 1:
        raise ValueError(f"Expected exactly one collision option to be true, got: {options}")
    reaction = active[0]
    if reaction not in {'electron_ionization', 'proton_charge_exchange', 'proton_ionization'}:
        raise ValueError(f"Unsupported reaction type: {reaction}")
    return reaction


def load_test_settings(config_path: Path, data_path: Path, preferred_name: Optional[str]) -> TestSettings:
    with config_path.open() as f:
        config = json.load(f)
    tests: Iterable[dict[str, Any]] = config.get('tests', [])
    tests = list(tests)
    if not tests:
        raise ValueError(f"No tests found in config {config_path}")

    dt_from_file = extract_dt_from_filename(data_path)
    stem = data_path.stem

    def score(entry: dict[str, Any]) -> int:
        current_score = 0
        prefix = str(entry.get('name_prefix', ''))
        if preferred_name:
            if preferred_name in prefix or prefix in preferred_name or preferred_name in stem:
                current_score += 2
        if prefix and prefix in stem:
            current_score += 1
        if dt_from_file is not None and 'dt' in entry:
            try:
                if math.isclose(float(entry['dt']), dt_from_file, rel_tol=1e-6, abs_tol=1e-3):
                    current_score += 1
            except Exception:
                pass
        return current_score

    scored = sorted(((score(t), t) for t in tests), key=lambda x: x[0], reverse=True)
    best_score, best_entry = scored[0]
    if best_score == 0:
        best_entry = tests[0]

    reaction = detect_reaction(best_entry['collision_options'])
    settings = TestSettings(
        name=str(best_entry.get('name_prefix', stem)),
        n0=float(best_entry['n0']),
        neutral_relative_density=float(best_entry.get('neutral_relative_density', 1.0)),
        num_particles=int(best_entry.get('num_particles', 0)),
        num_steps=int(best_entry.get('num_steps', 0)),
        neutrals_energy_kev=float(best_entry['neutrals_energy_kev']),
        particles_energy_kev=float(best_entry['particles_energy_kev']),
        charged_mass=float(best_entry['charged_mass']),
        neutral_mass=float(best_entry.get('neutral_mass', M_I_RATIO + 1.0)),
        dt=float(best_entry.get('dt', dt_from_file or 0.0)),
        collision_scheme=str(best_entry.get('collision_scheme', '')),
        reaction=reaction,
    )
    return settings


def load_time_and_neutral(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load time (omega_p^-1) and neutral counts from a text log."""
    data = np.loadtxt(path, comments='#', usecols=(1, 4))
    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError('Expected two columns (time, neutral_count) when reading simulation file.')
    time = data[:, 0]
    neutral_count = data[:, 1]
    return time, neutral_count


def compute_plasma_frequency(density_cm3: float) -> float:
    return math.sqrt(4.0 * math.pi * density_cm3 * QE ** 2 / ME)


def sample_maxwellian_velocity(temperature_eV: float, mass_ratio: float, size: int, rng: np.random.Generator) -> np.ndarray:
    """Sample 3D velocities from a Maxwellian with given temperature and mass ratio in me units."""
    kT = temperature_eV * EV_TO_ERG
    mass = mass_ratio * ME
    if kT <= 0.0 or mass <= 0.0:
        return np.zeros((size, 3))
    sigma = math.sqrt(kT / mass)
    return rng.normal(loc=0.0, scale=sigma, size=(size, 3))


def mono_beam_velocity(energy_eV: float, mass_ratio: float, direction: np.ndarray | None = None) -> np.ndarray:
    """Return a fixed velocity vector for a monoenergetic beam."""
    mass = mass_ratio * ME
    if direction is None:
        direction = np.array([0.0, 0.0, 1.0])
    norm = np.linalg.norm(direction)
    unit_dir = direction / norm if norm > 0 else np.array([0.0, 0.0, 1.0])
    speed = math.sqrt(max(0.0, 2.0 * energy_eV * EV_TO_ERG / mass))
    return unit_dir * speed


def sigma_v_monte_carlo(
    charged_temp_eV: float,
    neutral_temp_eV: float,
    neutral_energy_eV: float,
    charged_mass_ratio: float,
    neutral_mass_ratio: float,
    sigma_func,
    neutral_distribution: str,
    samples: int = 4000,
    seed: int = 13,
) -> float:
    """Monte Carlo average of sigma*v for Maxwellian charged species and optional neutral beam."""
    rng = np.random.default_rng(seed)
    charged_vel = sample_maxwellian_velocity(charged_temp_eV, charged_mass_ratio, samples, rng)

    if neutral_distribution == 'maxwell':
        neutral_vel = sample_maxwellian_velocity(neutral_temp_eV, neutral_mass_ratio, samples, rng)
    elif neutral_distribution == 'mono':
        beam = mono_beam_velocity(neutral_energy_eV, neutral_mass_ratio)
        neutral_vel = np.tile(beam, (samples, 1))
    else:
        raise ValueError("neutral_distribution must be 'maxwell' or 'mono'")

    rel_vel = charged_vel - neutral_vel
    speeds = np.linalg.norm(rel_vel, axis=1)
    energy_mc2 = 0.5 * charged_mass_ratio * np.square(speeds / C_LIGHT)
    sigma_vals = sigma_func(energy_mc2)
    return float(np.mean(sigma_vals * speeds))


def annotate_plot(ax, meta: SimulationMeta, settings: TestSettings, sigma_v: float, nu: float, gamma: float, neutral_distribution: str, samples: int) -> None:
    reaction_tex = meta.reaction.replace('_', r'\_')
    tau_text = r"\infty" if gamma <= 0 else f"{1.0 / gamma:.3e}"
    lines = [
        rf"$\mathrm{{Reaction~type:}}\ \mathrm{{{reaction_tex}}}$",
        rf"$\mathrm{{Charged~species:}}\ \mathrm{{{meta.species}}},\ m/m_{{\mathrm{{e}}}} = {meta.mass_ratio:.0f}$",
        rf"$T_{{\mathrm{{charged}}}} = {meta.temperature_eV/1e3:.3g}\,\mathrm{{keV}}$",
        rf"$T_{{\mathrm{{neutral}}}} = {settings.neutrals_energy_kev:.3g}\,\mathrm{{keV}}\ \mathrm{{({neutral_distribution})}}$",
        rf"$n_0 = {settings.n0:.3e}\,\mathrm{{cm^{{-3}}}}",
        rf"$n_{{neutrals}} = {meta.base_density*settings.neutral_relative_density:.3e}\,\mathrm{{cm^{{-3}}}}$",
        rf"$\langle\sigma v\rangle = {sigma_v:.3e}\,\mathrm{{cm^3\,s^{{-1}}}}$",
        rf"$\nu = {nu:.3e}\,\mathrm{{s^{{-1}}}}$",
        rf"$\nu/\omega_p = {gamma:.3e}$",
        rf"$\tau\omega_p = {tau_text}$",
        rf"$\mathrm{{samples}} = {samples}$",
    ]
    text = '\n'.join(lines)
    ax.text(
        0.52,
        0.98,
        text,
        transform=ax.transAxes,
        va='top',
        ha='left',
        fontsize=9,
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'),
    )


def select_sigma_function(reaction: str):
    if reaction == 'electron_ionization':
        return sigma_e_cm2
    if reaction == 'proton_ionization':
        return sigma_p_cm2
    if reaction == 'proton_charge_exchange':
        return sigma_cx_cm2
    raise ValueError(f'Unknown reaction type: {reaction}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='Path to neutral density data file.')
    parser.add_argument('--config', type=Path, default=Path('test_config.json'), help='Path to test_config.json.')
    parser.add_argument('--test-name', help='Optional test name prefix to pick a specific entry from the config.')
    parser.add_argument('--neutral-distribution', choices=['maxwell', 'mono'], default='maxwell', help='Neutral velocity distribution; charged particles are always Maxwellian.')
    parser.add_argument('--samples', type=int, default=4000, help='Monte Carlo samples for <sigma v>.')
    parser.add_argument('--seed', type=int, default=13, help='Random seed for sampling velocities.')
    parser.add_argument('--output', type=Path, help='Save plot to this path instead of showing interactively.')
    args = parser.parse_args()

    settings = load_test_settings(args.config, args.input, args.test_name)
    time_wp, neutral_density = load_time_and_neutral(args.input)
    if len(time_wp) == 0:
        raise ValueError('No data points found in input file.')

    sigma_func = select_sigma_function(settings.reaction)
    charged_temp_eV = settings.particles_energy_kev * 1e3
    neutral_energy_eV = settings.neutrals_energy_kev * 1e3
    sigma_v = sigma_v_monte_carlo(
        charged_temp_eV=charged_temp_eV,
        neutral_temp_eV=neutral_energy_eV,
        neutral_energy_eV=neutral_energy_eV,
        charged_mass_ratio=settings.charged_mass,
        neutral_mass_ratio=settings.neutral_mass,
        sigma_func=sigma_func,
        neutral_distribution=args.neutral_distribution,
        samples=args.samples,
        seed=args.seed,
    )

    plasma_frequency = compute_plasma_frequency(settings.charged_density)
    print(sigma_v)
    collision_frequency = settings.charged_density * sigma_v
    gamma = collision_frequency / plasma_frequency if plasma_frequency > 0 else 0.0

    time_shifted = time_wp - time_wp[0]
    neutral_norm = neutral_density / neutral_density[0]
    theory_norm = np.exp(-gamma * time_shifted)

    meta = SimulationMeta(
        species=settings.charged_species,
        temperature_eV=charged_temp_eV,
        base_density=settings.charged_density,
        reaction=settings.reaction,
        mass_ratio=settings.charged_mass,
    )

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(time_wp, neutral_norm, label='Model', linewidth=1.5, color='green')
    ax.plot(time_wp, theory_norm, label='Theory', linestyle='--', linewidth=1.5, color='black')
    ax.set_xlabel(r't, 1/$\omega_p$')
    ax.set_ylabel(r'$n_n$/$n_n$$_0$')
    ax.grid(True, which='both', alpha=0.3)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(time_wp[0], time_wp[-1])
    ax.legend()
    annotate_plot(ax, meta, settings, sigma_v, collision_frequency, gamma, args.neutral_distribution, args.samples)

    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=200)
    else:
        plt.show()


if __name__ == '__main__':
    main()
