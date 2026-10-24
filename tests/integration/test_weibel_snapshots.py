import tempfile
import unittest
from pathlib import Path

from snapshot_files import instantaneous_snapshots


class SnapshotSelectionTest(unittest.TestCase):
    def test_average_and_unrelated_files_are_excluded(self):
        for prefix in ("FieldB_PlaneZ_003_", "Density_PlaneZ_003_"):
            with self.subTest(prefix=prefix), tempfile.TemporaryDirectory() as directory:
                for suffix in ("0012", "0000", "avg0000", "avg0001", "", "2.backup", "-1", "１２"):
                    (Path(directory) / (prefix + suffix)).touch()
                (Path(directory) / "unrelated_0001").touch()
                self.assertEqual(instantaneous_snapshots(directory, prefix), [0, 12])


if __name__ == "__main__":
    unittest.main()
