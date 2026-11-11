"""Perpendicular fast magnetosonic wave: hot magnetized kinetic ions and
drift-kinetic (isothermal-pressure) electrons, from xpic
hybrid_drift_kinetic_ex2_fms.py. Units c = 1, n = 1, |q| = 1, m_e = 1,
temperatures in units of m_e c^2.
"""

import math

import numpy as np


class Parameters:
    def __init__(self, B0, n, mi, ti, te, k, dx, dt):
        self.B0, self.n, self.mi, self.Ti, self.Te = B0, n, mi, ti, te
        self.k, self.dx, self.dt = k, dx, dt
        self.vA2 = B0**2 / (n * mi)
        self.cs2 = 2 * (te + ti) / mi
        self.Omega_i = B0 / mi
        self.beta = 2 * n * (te + ti) / B0**2


def omega_fluid(p, k):
    return k * math.sqrt((p.vA2 + p.cs2) / (1 + p.vA2))


def bessel_i_scaled(order, x, terms=60):
    """exp(-x) I_n(x) by the power series, enough for x = (k rho_i)^2 << 1."""
    order = abs(order)
    term = (0.5 * x)**order / math.factorial(order)
    total = 0.0
    for m in range(terms):
        total += term
        term *= (0.5 * x)**2 / ((m + 1) * (m + 1 + order))
    return math.exp(-x) * total


def dispersion(p, k, omega, harmonics=20):
    lam = k**2 * p.Ti / (p.mi * p.Omega_i**2)  # (k rho_i)^2
    wp2 = p.n / p.mi
    exx, D, eyy = 1.0, 1 / (omega * p.B0), 1.0
    for m in range(-harmonics, harmonics + 1):
        I = bessel_i_scaled(m, lam)
        dI = bessel_i_scaled(1, lam) if m == 0 else 0.5 * (
            bessel_i_scaled(m - 1, lam) + bessel_i_scaled(m + 1, lam))
        r = wp2 / omega / (omega - m * p.Omega_i)
        exx -= r * m * m * I / lam
        D += r * m * (dI - I)
        eyy -= r * (m * m * I / lam + 2 * lam * (I - dI))
    eyy -= 2 * p.n * p.Te * k**2 / (p.B0 * omega)**2
    return exx * (eyy - (k / omega)**2) - D**2


def omega_kinetic(p, k):
    lo, hi = 0.5 * omega_fluid(p, k), min(1.5 * omega_fluid(p, k), 0.99 * p.Omega_i)
    f_lo = dispersion(p, k, lo)
    if not np.isfinite(f_lo) or f_lo * dispersion(p, k, hi) >= 0:
        raise ValueError("The kinetic FMS root is not bracketed")
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        f_mid = dispersion(p, k, mid)
        if np.sign(f_mid) == np.sign(f_lo):
            lo, f_lo = mid, f_mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def omega_discrete(p, omega):
    """Yee curls replace k by 2/dx sin(k dx/2); the implicit midpoint rule
    maps omega to 2/dt atan(omega dt/2)."""
    k = 2 / p.dx * math.sin(0.5 * p.k * p.dx)
    omega *= k / p.k
    return 2 / p.dt * math.atan(0.5 * omega * p.dt)
