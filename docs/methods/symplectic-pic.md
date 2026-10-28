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
periods, in either direction. Final particle buckets are rebuilt explicitly:
the existing `update_cells` requires neighboring-cell moves.

`currentOnGrid` accumulates `integral J dt` in both halves, is divided by the
**full** dt once after the step and is then the step-averaged current. Density
is deposited with the same nodal shape at initialization and after each step.
Neutrals drift and retain number density but contribute no charge/current.
Shared deposits are sequential in this first implementation; thread count
does not change their summation order. No RNG is used in the time advance.

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
