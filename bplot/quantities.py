"""Names of plottable 2D quantities and their values on the physical nodes of a plane.

NAME[@OWNER][=LIMITS]
  vector families  E B En Bn Eext Binit J  + component x y z r phi (none = magnitude)
  scalars          rho (total charge density), n (species density, |rho_s|),
                   Prr Ppp Pzz Prp Prz Pzp (species second moments)
  OWNER            a species name, '*' (all species, default) or 'total' (J only)
  LIMITS           a:b, or m for -m:m
"""
import re
from dataclasses import dataclass

import numpy as np

from .run import GHOST, IN_PLANE, read_slice

VECTORS = {  # family -> (file basename, Yee placement, belongs to a species)
    "Eext": ("FieldE_external", "edge", False), "Binit": ("FieldBInit", "face", False),
    "En": ("FieldEn", "edge", False), "Bn": ("FieldBn", "face", False),
    "E": ("FieldE", "edge", False), "B": ("FieldB", "face", False), "J": ("Current", "edge", True),
}
PRESSURES = ("Prr", "Ppp", "Pzz", "Prp", "Prz", "Pzp")
TOKEN = re.compile(r"^(?P<name>[A-Za-z]+)(?:@(?P<owner>[\w*]+))?(?:=(?P<lim>[-+.\deE:]+))?$")
VECTOR = re.compile(r"^(?P<family>Eext|Binit|En|Bn|E|B|J)(?P<comp>x|y|z|r|phi)?$")


@dataclass(frozen=True)
class Quantity:
    label: str
    owner: str | None  # species, or None for grid fields
    base: str          # file basename
    component: str     # x y z r phi abs, or '' for a scalar
    placement: str     # edge, face or node
    signed: bool       # diverging colormap, symmetric limits
    limits: tuple | None = None

    @property
    def needs_center(self):
        return self.component in ("r", "phi")


def parse(tokens, species):
    """Expand command-line tokens into quantities; raise ValueError with a hint."""
    result = []
    for token in tokens:
        m = TOKEN.match(token)
        if not m:
            raise ValueError(f"cannot read quantity '{token}'")
        name, owner, limits = m["name"], m["owner"], _limits(m["lim"])
        vector = VECTOR.match(name)
        if vector:
            base, placement, per_species = VECTORS[vector["family"]]
            component = vector["comp"] or "abs"
            if vector["family"] == "J" and owner == "total":
                per_species, base, owner = False, "TotalCurrent", None
            entry = (name, base, component, placement, component != "abs")
        elif name == "rho":
            per_species, entry = False, (name, "ChargeDensity", "", "node", True)
        elif name == "n":
            per_species, entry = True, (name, "Density", "abs", "node", False)
        elif name in PRESSURES:
            per_species, entry = True, (name, name, "", "node", name[1] != name[2])
        else:
            raise ValueError(f"unknown quantity '{name}'; see `bplot maps --help`")
        label, base, component, placement, signed = entry
        if not per_species:
            if owner not in (None, "total"):
                raise ValueError(f"'{name}' is a grid field and takes no species")
            result.append(Quantity(label, None, base, component, placement, signed, limits))
            continue
        owners = species if owner in (None, "*") else [owner]
        for s in owners:
            if s not in species:
                raise ValueError(f"unknown species '{s}'; run has {', '.join(species) or 'none'}")
            result.append(Quantity(f"{label} {s}", s, base, component, placement, signed, limits))
    return result


def _limits(text):
    if not text:
        return None
    if ":" in text:
        a, b = text.split(":", 1)
        return float(a), float(b)
    m = abs(float(text))
    return -m, m


def values(q, run, plane, frame, averaged, cache):
    """Values of q on the physical nodes of plane, shape (n1, n2)."""
    path = run.slice_path(q.owner, q.base, plane, averaged, frame.seq)
    if path not in cache:
        cache[path] = read_slice(path)
    raw = cache[path]
    axes = IN_PLANE[plane.axis]
    sizes = [run.cells[a] + 1 for a in axes]
    physical = (slice(GHOST, GHOST + sizes[0]), slice(GHOST, GHOST + sizes[1]))

    def component(c):
        return _to_nodes(raw[c], _shifted(q.placement, c), axes)[physical]

    if not q.component:
        return component(0)
    if q.component in "xyz":
        return component("xyz".index(q.component))
    if q.component == "abs" and q.placement == "node":
        return np.abs(component(0))
    if q.component == "abs":
        return np.sqrt(sum(component(c) ** 2 for c in range(3)))
    if plane.axis != "z":
        raise ValueError(f"{q.label}: r/phi components need a z plane")
    x, y = (coords - c for coords, c in zip(run.plane_axes(plane), run.center))
    x, y = np.meshgrid(x, y, indexing="ij")
    r = np.hypot(x, y)
    r[r == 0] = np.inf
    fx, fy = component(0), component(1)
    return (fx * x + fy * y) / r if q.component == "r" else (fy * x - fx * y) / r


def _shifted(placement, c):
    """Global axes along which component c sits half a cell above the node."""
    return {"edge": {c}, "face": {0, 1, 2} - {c}}.get(placement, set())


def _to_nodes(a, shifted, axes):
    for dim, axis in enumerate(axes):
        if axis in shifted:
            lower = np.roll(a, 1, axis=dim)
            a = 0.5 * (a + lower)  # node i = mean of the half nodes i - 1/2 and i + 1/2
    return a
