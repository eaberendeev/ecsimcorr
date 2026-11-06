#!/usr/bin/env python3
"""Check explicit project includes and the established module boundaries.

This checks source-level includes, including inactive preprocessor branches.
It does not infer C++ symbol dependencies or claim that the whole project is a DAG.
"""

import argparse
from pathlib import Path
import re
import sys

INCLUDE = re.compile(r'^\s*#\s*include\s*(["<])([^">]+)[">]', re.MULTILINE)
SOURCE_SUFFIXES = {'.h', '.hpp', '.cpp', '.cc'}
LEGACY_INCLUDE_DIRS = ('include', 'diagnostics', 'fields', 'particles', 'utils',
                       'simulation', 'simulation/ecsim', 'simulation/ecsim_corr')
LOWER = {'core', 'grid', 'include', 'runtime', 'utils'}
ALLOWED = {
    'core': {'core'},
    'grid': {'core', 'grid'},
    'include': {'core', 'include'},
    'runtime': {'core', 'runtime'},
    'utils': {'core', 'include', 'runtime', 'utils'},
    'fields': LOWER | {'fields'},
    'particles': LOWER | {'particles'},
    'collisions': LOWER | {'particles', 'collisions'},
    'coupling': LOWER | {'particles', 'coupling'},
    'io': LOWER | {'particles', 'io'},
    'config': LOWER | {'config', 'fields', 'particles', 'collisions'},
    'diagnostics': LOWER | {'config', 'particles', 'coupling', 'io', 'diagnostics'},
}
# Existing configuration construction and the shared boundary API. These exact
# edges are retained; they do not permit dependencies on concrete simulation steps.
EXCEPTIONS = {
    ('particles/particle_array.cpp', 'config/particle_config.h'),
    ('config/boundary_config.cpp', 'simulation/boundary_conditions.h'),
    ('diagnostics/DiagnosticMeasurements.h', 'simulation/boundary_conditions.h'),
    ('diagnostics/DensityMeasurements.h', 'simulation/boundary_conditions.h'),
}
SHARED_BUILD_COMMANDS = {'add_compile_options', 'add_compile_definitions', 'add_link_options'}
SCOPED_COMMANDS = re.compile(
    r'^\s*(include_directories|link_libraries|add_compile_options|add_compile_definitions|add_link_options)\s*\(',
    re.MULTILINE | re.IGNORECASE,
)
GLOBAL_LIBRARIES = re.compile(r'\$\{LIBRARIES\}|\bset\s*\(\s*LIBRARIES\b', re.IGNORECASE)


def check(source_dir):
    root = Path(source_dir).resolve()
    errors = []
    for path in sorted(root.rglob('*')):
        relative = path.relative_to(root)
        if '3rd_party' in relative.parts or not path.is_file():
            continue
        content = None
        if path.name == 'CMakeLists.txt' or path.suffix == '.cmake':
            content = path.read_text()
            for pattern in (SCOPED_COMMANDS, GLOBAL_LIBRARIES):
                for match in pattern.finditer(content):
                    # Common flags are set once before all project subdirectories.
                    if (relative.as_posix() == 'CMakeLists.txt' and pattern is SCOPED_COMMANDS
                            and match.group(1).lower() in SHARED_BUILD_COMMANDS):
                        continue
                    line = content.count('\n', 0, match.start()) + 1
                    errors.append(f'{relative}:{line}: use target-scoped build requirements: {match.group().strip()}')
        if path.suffix not in SOURCE_SUFFIXES:
            continue
        if content is None:
            content = path.read_text()
        for match in INCLUDE.finditer(content):
            delimiter, name = match.groups()
            line = content.count('\n', 0, match.start()) + 1
            location = f'{relative}:{line}'
            # Quoted includes first search beside the source; angle includes do not.
            candidates = ([path.parent / name] if delimiter == '"' else []) + [root / name]
            target = next((p.resolve() for p in candidates if p.is_file()), None)
            if target is None:
                # Recognize legacy angle includes such as <containers.h>, too.
                target = next(((root / directory / name).resolve()
                               for directory in LEGACY_INCLUDE_DIRS
                               if (root / directory / name).is_file()), None)
            if target is None:
                project_path = (root / Path(name).parts[0]).is_dir()
                if (delimiter == '"' or project_path) and not any(
                    (root / vendor / name).is_file()
                    for vendor in ('3rd_party', '3rd_party/eigen', '3rd_party/amgcl')
                ):
                    errors.append(f'{location}: unresolved project include: {name}')
                continue
            try:
                target_relative = target.relative_to(root)
            except ValueError:
                errors.append(f'{location}: include escapes source directory: {name}')
                continue
            if (delimiter != '"' or target.parent != path.parent) and name != target_relative.as_posix():
                errors.append(f'{location}: use explicit project path: {target_relative.as_posix()}')
            module = relative.parts[0]
            destination = target_relative.parts[0]
            edge = (relative.as_posix(), target_relative.as_posix())
            if module in ALLOWED and destination not in ALLOWED[module] and edge not in EXCEPTIONS:
                errors.append(f'{location}: forbidden dependency {module} -> {target_relative.as_posix()}')
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'srcBeren')
    args = parser.parse_args()
    if not (args.source_dir / 'core').is_dir():
        parser.error('--source-dir must contain the Beren3D source tree')
    errors = check(args.source_dir)
    for error in errors:
        print(error, file=sys.stderr)
    if errors:
        return 1
    print('Project includes and module boundaries: OK')
    return 0


if __name__ == '__main__':
    sys.exit(main())
