# Desktop performance workloads

These are short workloads for locating bottlenecks, not physical validation runs.
All defaults use **50 × 50 × 30 cells**, **dx = dy = dz = 0.5**, and **exactly
5 steps**. ECSIM and Corr use dt = 1.2; symplectic uses dt = 0.25 (Maxwell
CFL = 0.866 < 1). Compare cost per step; these runs cover different physical
durations. All populations are initialized once; there is no ongoing injection.

## Cases

PPC means macro-particles **per species per occupied cell**, averaged over random
positions. Electrons and hydrogen ions have masses 1 and 1836, temperatures
0.1 and 0.001 keV, and equal macro weights. n0 = 10¹⁸ cm⁻³ throughout, including
the controls for collision cases. No imposed fields, damping, checkpoints,
field slices, spectra or pressure diagnostics are enabled. `VerboseStep` is
enabled for the existing wall-time and solver log messages.

| `--case` | Initial population | Approximate total particles | Main target |
|---|---|---:|---|
| `dense` | Periodic, uniform e/i, 200 PPC each | 30 million | Particle gather, push, deposition, response assembly |
| `sparse` | Same periodic box, 2 PPC each | 300 thousand | Grid and matrix work relative to particle work |
| `mixed` | Periodic, 2 PPC background + 198 PPC extra in central 1/8 volume | 4.01 million | Uneven occupancy, cache behavior, OpenMP load imbalance |
| `cylinder` | Open cylinder, e/i, 200 PPC each | 15.08 million | Open-end losses, particle removal, cylindrical geometry |
| `cylinder-sparse` | Same cylinder, 2 PPC each | 151 thousand | Sparse occupancy and fixed grid costs |
| `cylinder-coulomb` | Same as `cylinder`, with e–e, i–i, e–i collisions | 15.08 million | Coulomb pairing, scattering, random numbers |
| `cylinder-neutral` | Same as `cylinder`, plus 200 PPC neutral hydrogen | 22.62 million | Control for neutral transport cost |
| `cylinder-ionization` | Same as `cylinder-neutral`, electron impact ionization | 22.62 million initially | Neutral lookup, cross sections, reactions, particle creation |

Only `dense`, `sparse`, and `mixed` support all three schemes. Cylinder cases
support `ecsim` and `ecsim_corr`; the current symplectic implementation requires
a collisionless periodic box.

The default cylindrical wall has radius **12**, center (12.5, 12.5), and open
ends z = 0, 15. Plasma starts at radius **10**, across the full length
**z = 0…15**, reaching both open ends. Its volume occupies about 50% of the
Cartesian box. Particles can leave through the ends from the first step;
the measured workload includes boundary handling and a declining population.
There is a radial vacuum gap of 2 to the cylindrical wall. Cylinder
collisions use ln Λ = 15; the ionization case uses `physical_only` electron
impact on H, with proton ionization and charge exchange disabled. Compare it
against `cylinder-neutral`, and Coulomb against `cylinder`.

`mixed` changes computational occupancy and physical density together: the
background is 0.01 n0 and the core is n0, with equal macro weights throughout.
It is a neutral initial density contrast, not an equilibrium. `sparse` retains
physical density n0 with fewer, heavier macro-particles. Counts are expected
global counts, not a guarantee of equal population in every cell.

## Build once, then run

From the repository root:

```bash
cmake -S srcBeren -B _build_desktop -DCMAKE_BUILD_TYPE=Release \
  -DUSE_TIMERS=ON -DMEMORY_TIMERS_LEVEL=0 -DBUILD_TESTS=OFF
cmake --build _build_desktop -j 6

python3 configs/desktop/run.py --binary _build_desktop/beren3d \
  --case dense --scheme ecsim_corr --threads 6 \
  --output test_results/desktop/dense_corr_t6_r1

python3 configs/desktop/run.py --binary _build_desktop/beren3d \
  --case cylinder-ionization --scheme ecsim --threads 6 \
  --output test_results/desktop/ionization_ecsim_t6_r1
```

Each output directory must be new. It contains the input JSONs, `phys.par`,
`run.json` (binary hash, platform and selected environment), `stdout.log`,
`beren3d.log`, `profile_it_5.json`, and `summary.json`. The runner uses the
production `beren3d` loop, so **collisions are included**. The existing
`step_timing` executable bypasses `collision_step` and is unsuitable here.

Step 1 warms up allocations, solver state and the reference assembly check;
only steps **2–5** enter the profile and reported mean/median. The runner sets
`VALIDATION_PERIODICITY` to the run length so reference checking stays on the
first step, including with `--steps` overrides. Initialization and first-step
cost remain in process elapsed time, separately from measured step time.
Step times come from the production `Step wall time` log: preparation,
advancement, collisions and diagnostics are included; profile serialization
and the wall-time log message itself are outside it. Other verbose solver
messages remain inside the measured step. `master_thread_timers` are inclusive
timers and **must not be added together** (nested scopes overlap); inspect the
Chrome trace for worker threads and parallel imbalance. Profiling itself adds
overhead. The runner rejects missing profiles and incomplete measured steps.

Start with sparse, then mixed, then dense. Dense e/i alone carry at least
1.44 GB of position/velocity data; cell-vector capacity, field matrices,
reference assembly and per-thread buffers add substantial memory. Use the
same thread count and binary when comparing cases. Repeat each case at least
three times in separate directories, then compare medians; use 1/2/4/6/etc.
threads for scaling. Five steps are a quick bottleneck probe; use `--steps 20`
or more to reduce noise after identifying an interesting case. Keep the
machine otherwise idle and record changes to solver environment options.

## Generate inputs without running

```bash
python3 configs/desktop/gen_config.py --case mixed --scheme symplectic \
  --output test_results/desktop/mixed_symplectic_inputs
```

`--ppc N` overrides dense/core PPC (sparse defaults to 2; mixed needs N ≥ 3).
`--cells NX NY NZ` scales geometry with the box; for a quick startup check use
`--cells 8 8 6 --ppc 8`. Such runs check configuration wiring, not desktop
performance. `build.py --config configs/desktop/gen_config.py --timers 1`
also accepts the default dense Corr configuration; use the runner above for
the fixed profiling window and structured result summary.
