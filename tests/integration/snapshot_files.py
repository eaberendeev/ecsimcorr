"""Snapshot filename selection, independent of plotting dependencies."""
import os


def instantaneous_snapshots(directory, prefix):
    suffixes = (name[len(prefix):] for name in os.listdir(directory) if name.startswith(prefix))
    return sorted(int(suffix) for suffix in suffixes
                  if suffix and all("0" <= char <= "9" for char in suffix))
