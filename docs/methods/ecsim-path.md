# ECSIM-Path

`"Scheme": "ecsim_path"` is ECSIM-Corr without the global energy correction.
Instead of rescaling the velocities after the step, each particle is pushed by
the field it actually crosses and deposits the path it actually takes, so the
work of the field is done particle by particle, where the particle is. The code
is in `srcBeren/simulation/ecsim_path/`; the class derives from
`SimulationEcsimCorr` and replaces four of its steps.

## The step

1. **Predictor current.** The ECSIM response `L` stays the point one at
   `x^{n+1/2}`, but the free current `J0` is the Esirkepov current of the
   straight path the particle is predicted to take, `x^n -> x^n + dt alpha v^n`
   (`pathGather::predict_current`). A full path response was tried: the same
   energy error at 17 s per step and 12 GB.
2. **Path push.** `v^n -> v^{n+1}` with the mean predicted field on the straight
   step path `x^n -> x^{n+1} = x^n + dt v_mid` (two fixed-point passes,
   `push_on_path_and_deposit`). The particle deposits that path (minus the first
   half-path, already deposited by the drift). Its energy then changes by the
   work of the field on its own current.
3. **One corrector pass**, always: the step is repeated with the corrected field,
   as `CorrectorPasses: 1` of ECSIM-Corr. Without it the energy error is larger
   and the fields are noisy at `dt = 10` (strong mirror: the field energy grows
   three times); `CorrectorPasses: 0` is rejected.
4. **Magnetization.** The straight path is the Crank-Nicolson chord of the Boris
   rotation, so it carries only `mu / (1 + theta^2)` of a gyrating particle's
   magnetic moment, `theta = |q/m| B dt / 2`. The rest goes to the field
   equations as a magnetization current: each particle adds
   `dM = -w m u^2 / (2B) theta^2 / (1 + theta^2) b` at `x^{n+1/2}`
   (`u` the gyration speed), `J_M = curl dM`. After the field advance its
   gyration energy changes by `-dM . (B^{n+1} - B^n)`; the field loses the same
   amount. At `Omega dt << 1` and at `B = 0`, `dM` vanishes. Diagnostics:
   `magnetizationWork`, and `magnetizationRemainder` for the energy a particle
   could not take (its perpendicular energy would go negative); the remainder
   is given to the species by a common velocity scale.

`EnergyCorrection` is rejected: there is nothing to correct.

## What holds and what does not

- Charge continuity is exact (residual 1e-15), far moves included.
- The step energy balance is not exact: the predictor freezes the path weights.
  The error falls with the step as about `dt^4` (5 keV plasma, 0.1 cell per
  unit time: 1e-5 at `dt = 2`, 1e-8 at `dt = 0.5`, `tests/unit/ecsim_path.cpp`).
- The magnetization exchange is exact: at `Omega_e dt = 2.5` the balance error
  is below 1e-2 of the magnetization work.

## Limits

- **Thermal path per step**: `v_te dt / dx <= 0.5`. Above it the frozen path
  weights make the step unstable (cold-plasma tests; a warning is printed at the
  start of the run). `omega_pe dt` alone does not limit it: a 1.25 keV plasma
  runs at `omega_pe dt = 10`.
- **Cells much larger than the Debye length**: on the magnetosonic case
  (`dx` = 226 Debye lengths) `omega_pe dt = 2`, `Omega_e dt = 0.4` heats the
  electrons through a slow electrostatic mode; `dt` halved or a second corrector
  pass removes it ([wave tests](../../tests/integration/wave_cases/README.md)).
- Boris itself rotates by `2 atan(Omega dt / 2)` instead of `Omega dt`, so at
  `Omega_e dt >> 1` the gyrophase and the orbit radius (`rho sqrt(1 + theta^2)`)
  are wrong, while `mu`, the grad-B drift and the mirror force stay right.

## Mirror runs

Weak mirror, 256 markers per cell, `t = 3000`, the reference is a run with a small
step:

| | electron energy | central density | central `Bz` |
|---|---:|---:|---:|
| reference | 76.84 | 0.627 | 0.4146 |
| `ecsim_path`, `dt = 5` | 76.34 | 0.618 | 0.4153 |
| `ecsim_path`, `dt = 10` | 79.4 | 0.610 | 0.4146 |
| `ecsim_corr`, 1 pass, `dt = 5` | 75.65 | 0.564 | 0.4144 |
| `ecsim_path` without magnetization, `dt = 5` | | | 0.4201 |

Without the magnetization the plasma loses its diamagnetism. In the strong
mirror (`t = 1500`, reference electron energy 92.18) `ecsim_path` gives 93.07
at `dt = 5` and 93.72 at `dt = 10`; the energy correction of ECSIM-Corr cools
the electrons to 86.3 and 67.0.

## Cost

1.15 to 1.5 times a step of ECSIM-Corr with one corrector pass: the path
weights are found by a trial Esirkepov pass per particle, twice in the push and
once in the predictor.
