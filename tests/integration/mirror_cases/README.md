# Mirror tests of the kinetic schemes

| Case | Input | What is compared |
|---|---|---|
| `loss_cone.py` + `loss_cone_check.py` | two coils, mirror ratio 5.67, open walls, 16 x 16 x 64 cells, 5 keV test electrons (16k markers, density 1e-6) at the centre, `Omega_e dt` = 2 at the centre and 11 in the throats | escaped fraction against `1 - sqrt(1 - B0 / Bmax)` (10 %); at most 1 % may leave in the second half of the run |

Each scheme has a thin `tests/integration/loss_cone_<scheme>/` directory:

```
python3 tests/integration/run_integration.py --test loss_cone_ecsim_path
```

## Results

2026-11-24, Release build; theory 0.0933 (averaged over the cloud).

| Scheme | `dt = 0.5` | `dt = 2` | `dt = 5` (`Omega_e dt` = 28 in the throats) |
|---|---:|---:|---:|
| `ecsim` | 0.0919 | 0.0917 | 0.0914 |
| `ecsim_corr`, 1 pass | 0.0919 | 0.0919 | 0.0917 |
| `ecsim_path` | 0.0919 | 0.0919 | 0.0917 |

The leak over the second half is 0.35 % in every run, independent of the step:
these are slow electrons of the cone still on their way out, not scattering.
All schemes push with Boris, which keeps the magnetic moment and the mirror
force at any `Omega_e dt` (it gets the gyrophase and the orbit radius wrong),
so the case checks the coupling of the push to the field solve, the far
moves and the open walls rather than distinguishing the schemes. The
self-consistent magnetization (diamagnetism) is not tested here: the plasma is
a test population.
