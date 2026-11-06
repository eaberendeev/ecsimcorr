"""Check that the dependency checker accepts valid code and reports deliberate errors."""

from pathlib import Path
import tempfile
import unittest

from check_architecture import check


class ArchitectureChecks(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.write('core/vector3.h')
        self.write('grid/grid.h', '#include "core/vector3.h"\n')
        self.write('simulation/ecsim/step.h')
        self.write('simulation/boundary_conditions.h')

    def write(self, name, text=''):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def test_allowed_geometry_and_external_includes(self):
        self.write('core/axis.h')
        self.write('core/vector3.h', '#include "axis.h"\n#include <vector>\n')
        self.write('3rd_party/nlohmann/json.hpp')
        self.write('config/data.h', '#include "nlohmann/json.hpp"\n')
        self.assertEqual(check(self.root), [])

    def test_rejects_lower_layer_depending_on_scheme(self):
        for module in ('core', 'grid', 'runtime', 'utils', 'particles', 'fields', 'coupling'):
            with self.subTest(module=module):
                name = f'{module}/invalid.h'
                self.write(name, '#include "simulation/ecsim/step.h"\n')
                self.assertTrue(any('forbidden dependency' in error for error in check(self.root)))
                (self.root / name).unlink()

    def test_boundary_exception_is_limited_to_exact_include(self):
        self.write('diagnostics/DensityMeasurements.h', '#include "simulation/boundary_conditions.h"\n')
        self.assertEqual(check(self.root), [])
        self.write('diagnostics/DensityMeasurements.h', '#include "simulation/ecsim/step.h"\n')
        self.assertTrue(any('forbidden dependency' in error for error in check(self.root)))

    def test_rejects_noncanonical_and_missing_project_includes(self):
        self.write('fields/invalid.h', '#include "../grid/grid.h"\n#include "missing.h"\n')
        errors = check(self.root)
        self.assertTrue(any('use explicit project path' in error for error in errors))
        self.assertTrue(any('unresolved project include' in error for error in errors))

    def test_recognizes_legacy_angle_project_includes(self):
        self.write('utils/containers.h')
        self.write('coupling/invalid.h', '#include <containers.h>\n#include <grid/missing.h>\n')
        errors = check(self.root)
        self.assertTrue(any('use explicit project path: utils/containers.h' in error for error in errors))
        self.assertTrue(any('unresolved project include: grid/missing.h' in error for error in errors))

    def test_sibling_angle_include_needs_root_path(self):
        self.write('particles/particle_array.h')
        self.write('particles/consumer.h', '#include <particle_array.h>\n')
        self.assertTrue(any('use explicit project path' in error for error in check(self.root)))
        self.write('particles/consumer.h', '#include "particle_array.h"\n')
        self.assertEqual(check(self.root), [])

    def test_rejects_directory_scoped_build_configuration(self):
        self.write('CMakeLists.txt', 'include_directories(include)\nset(LIBRARIES fields)\n')
        self.assertEqual(len(check(self.root)), 2)

    def test_shared_flags_are_allowed_only_in_root_cmake(self):
        self.write('CMakeLists.txt', 'add_compile_options(-Wall)\nadd_compile_definitions(USE_TIMERS)\nadd_link_options(--coverage)\n')
        self.assertEqual(check(self.root), [])
        self.write('fields/CMakeLists.txt', 'add_compile_options(-Wall)\n')
        self.assertEqual(len(check(self.root)), 1)

    def test_vendor_sources_are_excluded(self):
        self.write('3rd_party/vendor/file.h', '#include "nonexistent.h"\n')
        self.write('3rd_party/vendor/CMakeLists.txt', 'include_directories(include)\n')
        self.assertEqual(check(self.root), [])


if __name__ == '__main__':
    unittest.main()
