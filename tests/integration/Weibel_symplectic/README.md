# Symplectic Weibel validations

## The existing ECSIM-sized experiment

`gen_config.py` imports the physical configuration from `../Weibel_ecsim/gen_config.py`:
30 x 30 x 5 cells, dx=dy=dz=0.3, 500 electrons/cell (2,250,000 particles),
T_x=T_y=0.1 keV, T_z=1 keV, periodic boundaries, initial E=B=0, end time 750.
There is no quiet start or imposed magnetic seed in this scenario. The same
random initialization code and four-thread initialization are used.

Only the method, time step and output cadence change. The symplectic method
uses dt=0.15 (Maxwell CFL number 0.866), requiring 5,000 steps. The fields and
density are sampled every 14.4 normalized time units, exactly as in the ECSIM
fixture, by changing the slice interval from 12 to 96 steps. Recovery is added
every 500 steps (75 normalized time units); each backup overwrites the previous
one. No test directories are removed. `IntegrationHeavy=True` excludes this
case from the default integration sweep.

Build and run in a new directory, using the supervisor instead of the generic
integration wrapper (which assumes the ECSIM energy-work diagnostic):

```bash
cmake -S srcBeren -B /tmp/beren-symplectic-release -DBUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release -DUSE_TIMERS=ON
cmake --build /tmp/beren-symplectic-release --target beren3d -j 4
python3 tests/integration/Weibel_symplectic/supervise.py /path/to/new/run \
  --binary /tmp/beren-symplectic-release/beren3d
```

The supervisor caps OpenMP at four threads, checks diagnostics every 30 seconds,
writes `status.json` atomically, and imposes a 12-hour wall-time limit. It stops
the process group on malformed/non-finite diagnostics, relative Hamiltonian
change above 5%, Gauss-defect drift above 1e-8, or div B above 1e-8. These are
execution guards, not physical validation criteria. A three-step full-size
pilot on 2026-10-28 used 14.5 seconds and about 1.4 GB peak RSS; the complete
unoptimized sequential deposit is expected to take approximately 5–7 hours.

At completion `check_weibel.py` invokes the shared Weibel analysis, using actual
configured snapshot times. It checks the same mode growth, magnetic saturation,
polarization, phase drift and density fluctuation criteria as the corrected
ECSIM benchmark. The energy criterion for this method is instead
`max|H-H0|/H0 < 0.005`, since symplectic splitting does not exactly conserve H.
Gauss drift and div B must also stay below 1e-8. These criteria are fixed before
the production run; completing the run does not imply passing the criteria.
The existing ECSIM energy criterion is unchanged.

The reference contains electrons only and starts from E=0, without a Poisson
projection. Its initial `div E-rho` includes the omitted positive background
and particle-sampling fluctuations. To preserve exactly the requested initial
conditions, this test does not insert ions or project E. It monitors the change
of that initial defect (`gaussDrift`), not its absolute value (`gaussError`).
This matters when interpreting a comparison with the corrected ECSIM method.

Live monitoring and final analysis:

```bash
cat /path/to/run/status.json
tail /path/to/run/energy.txt
python3 tests/integration/Weibel_symplectic/check_weibel.py /path/to/run --quiet
```

Final plots and a pass/fail report are written to `Analysis/`; `run.log` and
`analysis.log` retain process output. `status.json` distinguishes a running
calculation, a process failure and a failed physical validation.

## Smaller independent 1D3V mode validation

`run.py` is a separate quiet-start test, not a substitute for the full fixture.
The C++ `symplectic_weibel` executable subclasses `SimulationSymplectic` only to
set initial data. All time advances use the production implementation.
It initializes a periodic plasma with one B_z(kx) mode, T_y/T_x=10,
sigma_y=0.1c, k=1.16 omega_pe/c and almost immobile neutralizing ions
(m_i/m_e=1e12). Equal-weight normal-quantile tensor sampling is deterministic:
32 v_x nodes, 8 v_y nodes and 2 passive v_z signs per cell. Nodes are rescaled
to the requested second moments. v_z does not enter the seeded E_y/B_z mode.
The finite velocity quadrature is an approximation to a Maxwellian.

The kinetic root is obtained from

```
gamma^2 + k^2 + 1 - A + A*sqrt(pi)*eta*erfcx(eta) = 0
eta = gamma/(sqrt(2)*k*sigma_x), A=T_y/T_x
```

with c=omega_pe=1. This is the transverse dispersion relation in
[Kuldinow & Hara (2025), Eq. (2.2)](https://doi.org/10.1017/S0022377825000303),
specialized to the y polarization. The theory function is independently checked
against the cold-parallel fluid polynomial and the instability cutoff.

```bash
cmake --build /tmp/beren-symplectic-release --target symplectic_weibel -j 4
python3 tests/integration/Weibel_symplectic/run.py /tmp/new-linear-run \
  --binary /tmp/beren-symplectic-release/tests/symplectic_weibel/symplectic_weibel
python3 tests/integration/Weibel_symplectic/run.py /tmp/new-linear-run --status
```

A fixed window t=60..110 is used for ln|B_k| fitting. Acceptance requires growth
within 10% of kinetic theory, fit R²>0.98, amplification above 30, relative H
error below 0.005, absolute Gauss error below 1e-10 and div B below 1e-12.
A control with A=1 requires amplification below three. These are numerical
assertions, not just CSV/plot generation. `weibel.csv`, `report.json` and
`weibel.png` retain the measurements. Reported velocity ratios are **raw second
moments**, not local thermal pressures with bulk flow removed.

Completed on 2026-10-28, theory gamma=0.05328678:

| Cells | dt | End time | Measured gamma | Error vs theory | Wall time |
|---|---|---|---|---|---|
| 32 | 0.1 | 160 | 0.05250598 | 1.465% | 52.6 s |
| 32 | 0.05 | 160 | 0.05250801 | 1.461% | 101.1 s |
| 64 | 0.05 | 110 | 0.05245994 | 1.552% | 138.5 s |
| 32, A=1 | 0.1 | 120 | no amplification above the seed | stable control | 39.2 s |

All four passed their criteria. The first run amplified the magnetic mode by
1231; its maximum relative energy error was 1.20e-5. Maximum absolute Gauss
error across the four runs was 1.98e-14; div B was zero for this polarization.
Time-step halving changed fitted gamma by 0.0039%; grid doubling changed it by
0.092%. The fits include residual electromagnetic transients and finite velocity
quadrature error; these measurements do not establish asymptotic spatial order
or the nonlinear saturation of the full 2D benchmark.
