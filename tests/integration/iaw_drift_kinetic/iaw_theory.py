"""Kinetic dispersion of the ion-acoustic wave along a uniform B.

Both species respond along the field only (drift-kinetic, k perpendicular
to nothing: k || B, T_perp = 0), so the electrostatic susceptibility is the
unmagnetized one-dimensional Maxwellian response
    chi_s = (omega_ps^2 / v_ts^2) (1 + zeta Z(zeta)),  zeta = omega / (sqrt(2) k v_ts),
with Z(zeta) = i sqrt(pi) w(zeta) (Faddeeva). On the grid the electric field
is gathered with S1 and the charge deposited with the Esirkepov shape, which
together with the finite-difference Gauss law multiply the field response by
sinc(k dz / 2)^4 (as in xpic drift_kinetic_ex5_iaw):
    k^2 / sinc^4 + sum_s chi_s(omega) = 0.
Units: c = 1, n = 1, |q| = 1, temperatures in keV with m_e c^2 = 511 keV.
"""

import math

import numpy as np
from scipy.optimize import fsolve
from scipy.special import wofz

MEC2_KEV = 511.0


def plasma_z(zeta):
    return 1j * math.sqrt(math.pi) * wofz(zeta)


def dispersion(omega, k, species, dz):
    sinc = math.sin(0.5 * k * dz) / (0.5 * k * dz) if dz > 0 else 1.0
    total = k * k / sinc**4
    for charge, mass, density, temperature_kev in species:
        vt = math.sqrt(temperature_kev / MEC2_KEV / mass)
        zeta = omega / (math.sqrt(2.0) * k * vt)
        total += density * charge**2 / mass / vt**2 * (1.0 + zeta * plasma_z(zeta))
    return total


def ion_acoustic_root(k, species, dz):
    """Least-damped forward root omega = omega_r - i gamma."""
    # Fluid guess: omega = k c_s, c_s^2 = (T_e + 3 T_i) / m_i.
    te = next(t for q, m, n, t in species if q < 0) / MEC2_KEV
    ti = next(t for q, m, n, t in species if q > 0) / MEC2_KEV
    mi = next(m for q, m, n, t in species if q > 0)
    guess = k * math.sqrt((te + 3 * ti) / mi)

    def residual(x):
        value = dispersion(x[0] + 1j * x[1], k, species, dz)
        return [value.real, value.imag]

    re, im = fsolve(residual, [guess, -0.1 * guess], xtol=1e-14)
    return re, -im


if __name__ == "__main__":
    species = [(-1.0, 1.0, 1.0, 5.0), (1.0, 25.0, 1.0, 0.1)]
    omega, gamma = ion_acoustic_root(2 * math.pi / 32, species, 1.0)
    print(f"omega = {omega:.16e}  gamma = {gamma:.16e}")
    print("xpic:   3.9304083806753808e-03          4.7291325876081436e-04")
