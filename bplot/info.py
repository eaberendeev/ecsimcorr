"""`bplot info`: what a run has written."""
import os
import textwrap

from .run import _read_table


def info(run):
    c = run.config
    cells = " × ".join(map(str, run.cells))
    print(f"run      {run.path}")
    print(f"grid     {cells} cells, d = {' '.join(f'{d:g}' for d in run.step)}, dt = {run.dt:g}, "
          f"scheme {c.get('Scheme', '?')}")
    if "CylinderDomain" in c:
        print(f"cylinder r = {c['CylinderDomain']['radius']:g} around ({', '.join(f'{v:g}' for v in run.center)})")
    print(f"species  {' '.join(run.species) or '-'}")

    print("\n2D slices: step range, interval; 'avg' = averaged windows (use --avg)")
    for plane in run.planes():
        print(f"  plane {plane.tag:<5} {plane.axis} = {run.plane_coordinate(plane):g}")
        owners = sorted({k[0] for k in run.slices if k[2] == plane}, key=lambda o: (o is not None, o or ""))
        for owner in owners:
            parts = []
            for base in sorted({k[1] for k in run.slices if k[2] == plane and k[0] == owner}):
                parts.append(f"{base} {_span(run, owner, base, plane)}")
            print(f"    {owner or 'fields':<10} " + "\n               ".join(parts))

    for name in ("energy.txt", "boundary.txt", "continuity.txt"):
        path = os.path.join(run.path, name)
        if os.path.isfile(path):
            header, data = _read_table(path)
            last = f", last step {int(data[-1, 0])}" if data.size and header[0] == "Step" else ""
            print(f"\n{name}: {len(data)} rows{last}")
            print(textwrap.fill(" ".join(h for h in header if h not in ("Step", "Time")), 100,
                                initial_indent="  ", subsequent_indent="  "))
    points = run.probe_points()
    if points:
        print("\nprobes: columns as Ex@p0 ...")
        for i, p in enumerate(points):
            print(f"  p{i}  ({', '.join(f'{v:g}' for v in p)})")
    for species in run.species:
        records = run.spectra(species)
        if records:
            print(f"\nspectra {species}: {len(records)} steps {records[0][0]} … {records[-1][0]}")


def _span(run, owner, base, plane):
    text = []
    for averaged in (False, True):
        frames = run.frames(owner, base, plane, averaged)
        if not frames:
            continue
        if averaged:
            text.append(f"avg {len(frames)} × {frames[0].last - frames[0].first + 1} steps")
        else:
            step = frames[1].first - frames[0].first if len(frames) > 1 else 0
            text.append(f"{len(frames)} steps {frames[0].first}…{frames[-1].first} by {step}")
        if len(text) == 2:
            break
    return "(" + ", ".join(text) + ")"
