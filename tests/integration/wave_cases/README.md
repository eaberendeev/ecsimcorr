# Wave tests of the kinetic schemes

Two linear waves with a known kinetic dispersion relation, run by every
kinetic scheme on the same input. The cases live here; each scheme has a thin
`tests/integration/<case>_<scheme>/` directory that selects the scheme.

| Case | Input | What is compared |
|---|---|---|
| `langmuir.py` | 64 x 1 x 1 cells, `dx = 0.1 c/omega_pe`, `k lambda_De = 0.4`, 8192 markers per cell, immobile ions | frequency (2 %) and Landau damping rate (20 %) against the kinetic root |
| `fms_kinetic.py` + `fms_check.py` | 32 x 1 x 1 cells, `dx = 10 c/omega_pe`, `B0 = 0.2`, `m_i/m_e = 100`, `Te = Ti = 1 keV`, 64 markers per cell, `Omega_e dt = 0.4` (`dt = 0.5` for the explicit scheme) | frequency of the fast magnetosonic wave (2 %) against the hot-ion kinetic root (`fms_theory.py`) |

Exactly conserving schemes must also keep the step energy balance
(`1e-6` and `1e-8` of the kinetic energy).

```
python3 tests/integration/run_integration.py --test fms_ecsim_corr
python3 tests/integration/run_integration.py --test all
```

## Results

Errors against the kinetic theory, 2026-11-12, Release build.

| Scheme | Langmuir: omega | Langmuir: Landau rate | Magnetosonic: omega | Energy balance (Langmuir / magnetosonic) |
|---|---:|---:|---:|---:|
| `ecsim` | 0.58 % | 6.4 % | 0.37 % | 7e-15 / 2e-11 |
| `ecsim_corr`, `CorrectorPasses: 1` in the magnetosonic case | 0.57 % | 7.3 % | 0.35 % | 8e-15 / 1e-13 |
| `symplectic` | 0.64 % | 5.3 % | 1.16 % | drift 4e-7 / 4 % |

The explicit symplectic scheme does not conserve energy exactly: on the
magnetosonic grid (`dx` = 226 Debye lengths) the electrons heat by 4 % over
61 thousand steps; the frequency is still right.

## ECSIM-Corr needs a corrector pass at `Omega_e dt > 0.3`

The single-pass `ecsim_corr` fails the magnetosonic case: the frequency is
6.7 % off, the ions lose 18 % of their energy and the electrons gain 15 %, at
an exact total balance. `ecsim` on the same input is clean.

What the runs showed (uniform plasma and the wave, 64 to 1024 markers per cell):

- It is an instability, not noise. With 1024 markers per cell the wave runs
  cleanly for half of the run; then `E_x`, `E_y` at `k dx` = 1.5-2 grow
  exponentially (rate 2e-4) at `omega = 0.189`, just below the numerical
  electron cyclotron frequency `2 atan(Omega_e dt / 2) / dt = 0.197`.
- The velocity rescaling is not the source. With `lambda = 1` the fields grow
  alike, the electrons heat out of nothing (+13 %) and the field solver fails
  in the end. With `lambda` the balance closes, but the ions pay.
- `Omega_e dt` controls it, not `omega_pe dt`:

  | `Omega_e dt` | steps per gyration | ions | electrons | `E` field |
  |---:|---:|---:|---:|---|
  | 0.20 | 31 | +0.5 % | 0 % | steady |
  | 0.25 | 25 | +0.4 % | +0.3 % | steady |
  | 0.30 | 21 | +0.4 % | +0.3 % | steady |
  | 0.35 | 18 | +0.5 % | +0.1 % | steady |
  | 0.40 | 16 | -5.4 % | +4.9 % | grows 100 times |

  (`dt = 2`, `B0` varied, 1024 markers per cell, 30 thousand time units.) At
  the same `Omega_e dt = 0.4` with `dt = 1`, `B0 = 0.4` the ions lose 79 %.
  The bound is from this one setup (`dx = 10 c/omega_pe`,
  `omega_pe / Omega_e` = 5-10), not a general criterion.
- It is not the one-cell periodic axes (32 x 2 x 2 and 32 x 4 x 4 cells behave
  alike), not the temperature (0.25 keV: the same), and no particle crosses a
  cell (`v_te dt / dx` = 0.009).

The cause is that the particles are pushed by the predicted field `E'` while
the field advances with `(E^n + E^{n+1}) / 2` of their Esirkepov current.
Why the mismatch pumps the cyclotron mode has not been derived; a cold linear
model of the predictor-corrector pair is only damped.

### `CorrectorPasses`

Each pass moves `E'` by one quasi-Newton step towards `(E^n + E^{n+1}) / 2`,

    (I - M + L) dE = (I - M) ((E^n + E^{n+1}) / 2 - E'),

with the ECSIM response `L` in place of the derivative of the Esirkepov
current, returns the particles to `x^{n+1/2}`, `v^n` and repeats the push. The
mismatch falls about 15 times per pass. The final current is still the
Esirkepov current of the final paths, so charge continuity holds to rounding
(residual 5e-15 with 0, 1 and 2 passes) and the energy balance stays exact.

| Magnetosonic case | no pass | 1 pass | 2 passes |
|---|---|---|---|
| 64 markers per cell, `dt = 2`: omega error | 6.7 %, fails | 0.35 % | 0.50 % |
| same: ions / electrons | -17.7 % / +14.6 % | -3.0 % / +2.8 % | -0.6 % / +0.3 % |
| 1024 markers per cell, `dt = 2`: omega error | 0.74 % | 0.10 % | |
| same: ions / electrons | -5.4 % / +4.9 % | +0.2 % / 0 % | |
| 1024 markers, `dt = 5` (`Omega_e dt = 1`): omega error | 14.5 %, fails | 0.15 % | 0.12 % |
| 1024 markers, `dt = 10` (`Omega_e dt = 2`): omega error | | 0.08 % | 0.10 % |
| 1024 markers, `dt = 1`, `B0 = 0.4`: ions | -79 % | +0.2 % | |

### What the passes change elsewhere

Checked on 2026-11-13 with 0, 1 and 2 passes; `ecsim` is the reference where
it is listed. Magnetosonic case with 1024 markers per cell unless noted.

| Time step, `B0 = 0.2` | `Omega_e dt` | omega error: 0 / 1 / 2 passes | ions: 0 / 1 passes |
|---:|---:|---|---|
| 0.5 | 0.1 | 0.16 % / 0.07 % / | +0.2 % / +0.2 % |
| 1 | 0.2 | 0.16 % / 0.08 % / | +0.1 % / +0.2 % |
| 2 | 0.4 | 0.74 % / 0.10 % / 0.04 % | -5.4 % / +0.2 % |
| 5 | 1 | 14.5 % / 0.15 % / 0.12 % | +1.1 % (electrons -24 %) / +0.2 % |
| 10 | 2 | run aborts at 40 % / 0.08 % / 0.10 % | -17 % / +0.2 % |

| Field, `dt = 2` | `Omega_e dt` | omega error: 0 / 1 passes / `ecsim` | ions: 0 / 1 passes / `ecsim` |
|---:|---:|---|---|
| 0.05 | 0.1 | 3.2 % / 0.34 % / | +1.0 % / 0 % / |
| 0.1 | 0.2 | 0.31 % / 0.05 % / | 0 % / 0 % / |
| 0.4 | 0.8 | 10 % / 0.11 % / 0.16 % | -18 % / +0.8 % / +0.9 % |
| 0.8 | 1.6 | 0.18 % / 0.17 % / 0.18 % | +7.5 % / +9.4 % / +9.1 % |

The ion heating at `B0 = 0.8` is the damping of the wave: `ecsim` gives the same.

- Where the single-pass scheme is stable (`Omega_e dt` <= 0.2) a pass changes
  nothing but the frequency error, which halves.
- Langmuir wave (no magnetic field): the frequency and damping errors with 1
  and 2 passes are those of the single-pass scheme to five digits.
- Uniform plasma, 64 markers per cell, 4000 steps: with passes the noise goes
  to the `ecsim` level. Without a field the magnetic noise energy is 1.4e-3
  / 9.5e-5 / 1.3e-4 for 0 / 1 / 2 passes (`ecsim` 1.9e-4); at `B0 = 0.2` the
  ions lose 3.7 % / 1.2 % / 0.6 % (`ecsim` 0.6 %).

Weibel instability (`Weibel_ecsim_corr`, no external field): all 11 checks
pass with 0, 1 and 2 passes. The linear stage is the same; the nonlinear one
is not:

| | 0 passes, `dt = 1.5` | 0 passes, `dt = 0.75` | 1 pass, `dt = 1.5` | 2 passes, `dt = 1.5` | 1 pass, `dt = 0.75` |
|---|---:|---:|---:|---:|---:|
| growth rate of mode (1,1), theory 0.0232 | 0.0197 | 0.0205 | 0.0199 | 0.0199 | 0.0206 |
| peak magnetic energy | 0.393 | 0.404 | 0.498 | 0.498 | 0.495 |
| time of the peak | 358 | 358 | 627 | 621 | 622 |
| magnetic energy at t = 750 | 0.240 | 0.289 | 0.460 | 0.461 | 0.458 |
| `Bz_rms / Bxy_rms` at the peak | 0.125 | 0.125 | 0.384 | 0.375 | 0.384 |

With a pass the nonlinear stage does not depend on the time step or on the
number of passes (within 1 %); without it the field decays after the first
maximum, less so at the smaller step. The symplectic scheme (`dt = 0.15`,
`Weibel_symplectic/README.md`) has its maximum 0.438 at t = 620 and 0.402 at
t = 750. Which nonlinear level is right is not established: the corrected
ECSIM-Corr is 14 % above the symplectic one, and no resolution study was made.

Cost: a copy of the particles of every charged species, one more push with
deposition and two more field solves per pass. On the 64-marker case a step
takes 18 % longer with one pass (single thread).

Limits seen: at `dt = 25` and `dt = 50` (`c dt / dx` = 2.5 and 5) the
unpreconditioned field solver stops converging for `ecsim` and `ecsim_corr`
alike, so larger steps were not tested.
