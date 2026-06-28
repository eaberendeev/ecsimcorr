#!/usr/bin/env python3
"""Convert old-format particle backup files (with initCoord+initVelocity) to new format.

Old Particle layout:  coord(3d) + velocity(3d) + initCoord(3d) + initVelocity(3d) = 96 bytes
New Particle layout:  coord(3d) + velocity(3d)                                    = 48 bytes

Usage:
    python convert_backup.py <Recovery_dir>
    python convert_backup.py data/Recovery

Converts all .backup files under <dir>/Particles/ in-place (old files saved as .backup.old).
"""

import struct
import sys
import shutil
from pathlib import Path

OLD_PARTICLE_SIZE = 12 * 8  # 96 bytes
NEW_PARTICLE_SIZE = 6 * 8   # 48 bytes


def convert_file(path: Path):
    data = path.read_bytes()
    if len(data) < 4:
        print(f"  SKIP {path} (too small)")
        return

    n_particles = struct.unpack_from("i", data, 0)[0]
    expected = 4 + n_particles * OLD_PARTICLE_SIZE

    if len(data) != expected:
        if len(data) == 4 + n_particles * NEW_PARTICLE_SIZE:
            print(f"  SKIP {path} (already new format, {n_particles} particles)")
            return
        print(f"  ERROR {path}: expected {expected} bytes, got {len(data)}")
        return

    backup_path = path.with_suffix(".backup.old")
    shutil.copy2(path, backup_path)
    print(f"  saved original -> {backup_path}")

    with open(path, "wb") as f:
        f.write(struct.pack("i", n_particles))
        offset = 4
        for _ in range(n_particles):
            f.write(data[offset : offset + NEW_PARTICLE_SIZE])
            offset += OLD_PARTICLE_SIZE

    new_size = 4 + n_particles * NEW_PARTICLE_SIZE
    print(f"  converted {path}: {n_particles} particles, {len(data)} -> {new_size} bytes")


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <Recovery_dir>")
        sys.exit(1)

    recovery_dir = Path(sys.argv[1])
    particles_dir = recovery_dir / "Particles"

    if not particles_dir.is_dir():
        print(f"Error: {particles_dir} not found")
        sys.exit(1)

    for backup_file in sorted(particles_dir.rglob("*.backup")):
        convert_file(backup_file)

    print("Done.")


if __name__ == "__main__":
    main()
