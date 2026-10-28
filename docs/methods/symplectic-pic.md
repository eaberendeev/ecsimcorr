# Symplectic Vlasov–Maxwell PIC

`Scheme: "symplectic"` implements the symmetric second-order Hamiltonian
splitting of Xiao, Qin, Liu, He, Zhang and Sun, *Physics of Plasmas* **22**,
112504 (2015), [DOI](https://doi.org/10.1063/1.4935904),
[author manuscript](https://arxiv.org/abs/1510.06972).
The local library contains the published paper as `sympl_PIC.Id_10093.pdf`.
The implementation uses the **lowest-order cubical Whitney forms** described
in the appendix; it does not implement the paper's higher-order spatial shapes.

## Discretization and time levels

The evolved state is `(E_n, B_n, x_n, v_n)`. Positions and velocities are
synchronous, unlike the staggered particle positions in ECSIM. All coefficients
use the existing normalization `c=1`; velocities are nonrelativistic.
The Yee locations are exactly `Domain::get_node_position`: charge at nodes,
`E_d` on edges shifted by half a cell in `d`, `B_d` on faces shifted by half a
cell in the transverse directions. Only `num_cells` unique periodic nodes
are evolved. The endpoint and ghosts are copies, excluded from field energies.

Let `lambda_i(u)` be the piecewise linear nodal hat and `chi_i(u)` the
piecewise constant indicator of cell `[i,i+1)`. The reconstruction is

- charge: `lambda_i(x) lambda_j(y) lambda_k(z)`;
- `E_x`: `chi_i(x) lambda_j(y) lambda_k(z)`;
- `B_x`: `lambda_i(x) chi_j(y) chi_k(z)`;
- the other components follow cyclic permutation.

Here coordinates in each factor are scaled by that axis's cell size through
`Domain::to_cell_coordinates`. These reconstructions commute with discrete
curl. Ordinary CIC interpolation centered at every Yee component does not
provide this property. Faraday uses forward differences; Ampère uses their
transpose, equivalently backward differences on the periodic mesh.

`ParticlesArray::mpw()` is `density / NumPartPerCell`. Thus nodal charge is
`sum(q*mpw*shape)`, and integrated edge current is `sum(q*mpw*integral(dx*shape))`.
There is no extra division by cell volume. Energies, as elsewhere in Beren,
are sums of normalized energy density contributions: multiply **both** field
and particle sums by `cell_volume` to obtain the volume integral.

## Exact flows

The Hamiltonian is split into electric energy, magnetic energy and the three
Cartesian particle kinetic energies, Eqs. (26)–(29). Each stage is an exact
solution of its subsystem:

- `H_E`: `B -= dt curl E`, `v += dt (q/m) gather_E(E,x)`; E and x are fixed.
- `H_B`: `E += dt curl^T B`; B, x and v are fixed.
- `H_a`: move only `x_a` by `dt*v_a`, update transverse velocities by
  `(q/m) integral(dx_a e_a cross B(x))`, and subtract the integrated edge
  current from E. B and `v_a` are fixed.

A full step follows Eq. (55):

```
Hx(dt/2), Hy(dt/2), Hz(dt/2), HB(dt/2), HE(dt),
HB(dt/2), Hz(dt/2), Hy(dt/2), Hx(dt/2)
```

Each axis trajectory is split at every grid plane. On each piece the required
transverse B and edge-current shapes are constant along the direction of
motion, so the line integrals are exact. Coordinates are wrapped only after
the signed path has been integrated. This includes multiple cells and multiple
periods, in either direction. Final particle buckets are updated by extracting
only particles whose destination cell differs and inserting them after the scan.
The existing `update_cells` requires neighboring-cell moves.

`currentOnGrid` accumulates `integral J dt` in both halves, is divided by the
**full** dt once after the step and is then the step-averaged current. Density
is deposited with the same nodal shape at initialization and after each step.
Neutrals drift and retain number density but contribute no charge/current.
Shared deposits use one scratch field per OpenMP thread. Each H_axis stage
sums these currents in thread-index order, then updates E and the species'
integrated current. This is valid because B is fixed and E does not enter the
particle motion within H_axis. Density uses separate thread-local buffers.
Buffers are cleared before every stage/species and are not checkpointed.
Changing the thread count changes summation order and may change round-off;
particle initialization and RNG streams are unchanged by this optimization.
No RNG is used in the time advance.

## Supported model and constraints

This implementation supports a periodic rectangular domain on all three axes,
without collisions, injection, absorbing layers or imposed time-dependent
fields. Unsupported configurations fail during initialization. Uniform initial
E and B are set via

```json
"Symplectic": {
  "initial_electric": [0.01, 0, 0],
  "initial_magnetic": [0, 0, 0.1]
}
```

B is the **total** evolved magnetic field, including uniform magnetic flux.
`ExternalFieldE`/`ExternalFieldB` must be absent or empty. The initial vectors
are ignored on recovery. Recovery stores synchronous E, B, x and v with the
existing output factory and file format. Only restart from this scheme with
the same grid/species/normalization is valid; ECSIM checkpoints have different
time levels. Diagnostic reference energy and Gauss defect are reset on restart.

A conservative Maxwell bound is enforced:
`Dt * sqrt(1/Dx^2 + 1/Dy^2 + 1/Dz^2) < 1`.
This bound follows from the vacuum Maxwell part of the implemented splitting;
it is not a quoted general stability theorem from the paper. In the published
paper the Bernstein-wave example uses `Dt = Dl/(2c)`, Eq. (62).

The exactness of each split flow for arbitrary dt does not imply stability of
the composed step. For a homogeneous cold electron plasma, with immobile ions
and `omega_pe=1`, the `(E,v)` map is

```
[ E_new ]   [ 1-dt^2/2   dt*(1-dt^2/4) ] [ E ]
[ v_new ] = [   -dt         1-dt^2/2   ] [ v ]
```

Its determinant is one and its trace is `2-dt^2`. It is linearly stable for
`0 < omega_pe*dt < 2`; the endpoint 2 is marginal and generally has a growing
Jordan mode. The numerical oscillation frequency in the stable range is
`omega_num = 2/dt * asin(omega_pe*dt/2)`. Thus `omega_pe*dt=1.5` is stable in
this test but gives about 13.1% frequency error. A step greater than the plasma
period `2*pi/omega_pe` is outside this stability interval. The unit test checks
bounded energy at 1.5 and 1.9 and exponential numerical growth at 2.1, using a
coarse grid that satisfies Maxwell CFL in all three cases. These are tests of
the homogeneous cold limit, not sufficient conditions for a general PIC run.

This does not guarantee stability at an unresolved plasma or gyrofrequency.
Unlike the implicit ECSIM schemes, the explicit scheme requires temporal
resolution of those scales as well. Symplecticity does not imply exact energy
conservation: bounded oscillatory energy error is expected for resolved motion.

The method preserves the **initial** Gauss defect `div E - rho`. It does not
solve Poisson's equation or repair inconsistent initial data. A physical
initial state must satisfy Gauss's law, including any specified immobile
neutralizing background; the code does not insert such a background. In
particular, uniform E requires pointwise charge neutrality. The example uses
co-located equal-weight opposite species. Initial B must be divergence-free.

## Diagnostics and reproducible checks

Existing `diagnostics.outputs` controls files; absent outputs write no diagnostic
files. `energy_balance` includes per-species energies, `energyFieldE`,
`energyFieldB`, `energyTotal`, signed `energyError` relative to initialization,
`gaussError` (maximum magnitude of `div E-rho`),
`gaussDrift` (maximum change of `div E-rho`), and `divB` (maximum magnitude).
`energyError` is a total change since initialization, not a per-step ECSIM
work residual. Initial Gauss error may be nonzero even when `gaussDrift=0`.
The console summary uses this signed total change for its `dE` column.
Generic probes/pressure outputs use their existing interpolation/moment rules,
not the Whitney gather used by the advance. Field dumps and slices store the
actual Yee coefficients.

Build and run independent numerical assertions:

```bash
cmake -S srcBeren -B /tmp/beren-symplectic -DBUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Debug -DUSE_TIMERS=ON
cmake --build /tmp/beren-symplectic -j 4
cmake --build /tmp/beren-symplectic --target run_symplectic_invariants
cmake --build /tmp/beren-symplectic --target run_tests
```

The new unit target checks local continuity for positive/negative, zero,
multiple-cell and multiple-period paths on an anisotropic grid; exact uniform-B
Lorentz integrals; commuting curl reconstruction; `div curl=0`; reversal of the
full particle-field map; preservation of the full noncanonical Poisson tensor
by a finite-difference Jacobian; second-order error against analytical cold
plasma oscillations and a semidiscrete vacuum Maxwell mode; the cold-plasma
energy envelope over 20,000 steps; factory selection, neutral drift, complete
step continuity/Gauss preservation, particle rebinning, invalid configurations,
actual energy output, repeated finalize and continuation through disk recovery.

A small runnable electron–positron experiment (creates a **new** directory):

```bash
python3 examples/symplectic_config.py /tmp/beren-symplectic-example --steps 200
cd /tmp/beren-symplectic-example
OMP_NUM_THREADS=2 /tmp/beren-symplectic/beren3d
```

This experiment exercises the production executable and outputs, but does not
replace the quantitative assertions. Landau damping and Bernstein dispersion
from the paper remain larger physical validations; they are not claimed here.
The old history (`7efb974`, `4de0c69`) was inspected but not restored: its last
version injects a hardcoded test particle, prints debug arrays and exits after
two steps. The new implementation directly follows the exact flows above.

### Verification on 2026-10-28

GCC Debug with timers: all 17 tests labelled `unit` passed. The new invariant
executable also passed separately with `OMP_NUM_THREADS=1` (CTest uses 2).
Its Poisson-map maximum defect was `1.32443e-10`; halving dt gave error ratios
`4.00291` for cold plasma and `4.0062` for the vacuum mode. The cold-plasma
relative energy envelope was `0.000625` over 20,000 steps at dt=0.05.
The 200-step production example wrote 201 energy rows; maximum absolute Gauss
error and drift were `1.60e-15`, maximum div B `2.16e-16`, and maximum energy
change divided by initial electric energy `1.22e-3`. All table values were finite.
Heavy validation targets, full integrations, Landau damping and Bernstein-wave
benchmarks were not run.

## Weibel experiments

See [the two Weibel validations](../../tests/integration/Weibel_symplectic/README.md).
The quiet-start 1D3V mode test passed a kinetic growth-rate comparison, half-dt
repeat, doubled-grid repeat and an equal-second-moment stable control. The
original ECSIM-sized 30x30x5, 500-particles/cell fixture is configured separately
with a durable supervisor; completion and pass/fail are recorded in the run's
`status.json`.


## Particle-kernel optimization

The original tensor-product gather, charge and kinetic-flow kernels remain as
`gather_reference`, `deposit_charge_reference` and `kinetic_flow_reference`.
The production kernels wrap each needed node once, omit identically zero
weights, reconstruct only the transverse B components needed by an axis flow,
and cache the fixed transverse stencil during that flow. The time composition,
line-integral segmentation, particle shapes and energy definition are unchanged.

`step_timing` measures prepare+advance separately from diagnostics, excluding
initialization and two warm-up steps. It emits every sample to `step_timing.csv`
and coarse stage profiling to `step_profile.json`. Run from a directory with
the experiment's two JSON files, for example:

```bash
cmake --build /tmp/beren-symplectic-opt --target step_timing -j 4
OMP_NUM_THREADS=4 /tmp/beren-symplectic-opt/tests/step_timing/step_timing 2 10
```

On 2026-10-28, Release with timers, the ECSIM-sized Weibel population
(2,250,000 particles) gave these seconds per prepare+advance step:

| Implementation | Threads | Mean seconds per step, 10 measured steps |
|---|---|---|
| Original symplectic | 4 (serial deposits) | 4.721 |
| Optimized symplectic | 4 | 0.394; repeat under build load 0.509 |
| Optimized symplectic | 1 | 1.242 |
| ECSIM | 4 | 0.848; repeat under build load 0.972 |

The original profile spent about 70% in kinetic flows and 22% in the electric
particle kick; diagnostics were below 1%. The optimized method is faster **per
step** here, while it still needs eight times as many steps as ECSIM at
Dt=1.2 versus 0.15. These timings do not imply it is faster per unit physical
time or that the two methods have the same error at their respective Dt.

Checks include optimized/reference comparisons on anisotropic grids, exact
plane starts and multi-period trajectories; the unchanged Poisson-map,
continuity, reversal, convergence, energy-envelope and recovery tests; and a
25,000-particle mixed charged/neutral test comparing serial with 4->2->4
threads, including neutral multi-period rebinning. All 17 Debug unit tests
passed. The repeated 1D3V Weibel test measured gamma=0.052505981684798,
matching the pre-optimization 0.052505981684793 to round-off. Full-size
12-step energy/constraint diagnostics also agreed to round-off. The full
nonlinear ECSIM-sized optimized Weibel validation is supervised separately.

The original and optimized full-size runs were also compared through step 732
(t=109.8): maximum differences in the printed electric/magnetic energies were
1.01e-14 and in signed Hamiltonian error 1.24e-14. This checks optimization
agreement over many cell crossings; physical acceptance of the full nonlinear
run remains a separate pending validation.
