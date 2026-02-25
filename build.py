#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import sys
from gen_config import generate_config

default_workdir = "./_build/bin"

def run(cmd, cwd=None):
    print(f"+ {cmd}")
    subprocess.check_call(cmd, shell=True, cwd=cwd)


def main():
    parser = argparse.ArgumentParser(description="Build beren3d with CMake.")
    parser.add_argument(
        "--type",
        default="Release",
        choices=["Debug", "Release"],
        help="CMAKE_BUILD_TYPE",
    )
    parser.add_argument(
        "--rebuild", action="store_true", help="Remove _build before configure"
    )
    parser.add_argument("--rerun", action="store_true", help="Remove workdir")
    parser.add_argument(
        "--eigen",
        default=os.environ.get("EIGEN_PATH", ""),
        help="Path to Eigen (overrides EIGEN_PATH env)",
    )
    parser.add_argument(
        "--amgcl",
        default=os.environ.get("AMGCL_PATH", ""),
        help="Path to AMGCL (overrides AMGCL_PATH env)",
    )
    parser.add_argument(
        "--jobs", type=int, default=os.cpu_count() or 8, help="Parallel build jobs"
    )
    args = parser.parse_args()

    root = os.path.abspath(os.path.dirname(__file__))
    build_dir = os.path.join(root, "_build")
    src_dir = os.path.join(root, "srcBeren")

    if args.rebuild and os.path.isdir(build_dir):
        print(f"Removing build dir: {build_dir}")
        shutil.rmtree(build_dir)

    os.makedirs(build_dir, exist_ok=True)

    cmake_config = [
        f"-DPATH_TO_EIGEN={args.eigen}",
        f"-DPATH_TO_AMGCL={args.amgcl}",
        f"-DCMAKE_BUILD_TYPE={args.type}",
        src_dir,
    ]

    run(f"cmake {' '.join(cmake_config)}", cwd=build_dir)
    run(f"cmake --build . -j{args.jobs}", cwd=build_dir)
    run(f"cmake --install .", cwd=build_dir)

    bin_path = os.path.join(build_dir, "bin", "beren3d")
    if not os.path.exists(bin_path):
        print("Error: binary not found after build:", bin_path, file=sys.stderr)
        sys.exit(1)

    print("\nBuild finished.")
    print("Binary:", bin_path)
    print(
        "Tip: For Debug run directly from _build/bin, for runs use gen_config + run_local.sh."
    )

    workdir = generate_config()

    if args.rerun:
        print(f"Removing work directory: {workdir}")
        shutil.rmtree(workdir, ignore_errors=True)

    if args.type == "Debug":
        print("Running in Debug mode.")
        workdir = default_workdir
    else:
        print("Running in Release mode.")
        if os.path.exists(workdir):
            if args.rerun:
                shutil.rmtree(workdir)
            else:
                print(
                    f"Ошибка: рабочая директория {workdir} уже существует. Используйте --rerun для удаления.",
                    file=sys.stderr,
                )
                sys.exit(1)
        os.makedirs(workdir)
        shutil.copy(build_dir + "/bin/" + "beren3d", workdir)
        shutil.copytree("srcBeren", workdir + "/srcBeren")
        shutil.copytree("PlotScripts", workdir + "/PlotScripts")
        shutil.copy("run.sh", workdir)
        shutil.copy("build.py", workdir)
        shutil.copy("gen_config.py", workdir)

    for fname in ["system_config.json", "particles_config.json", "phys.par"]:
        src = fname
        if os.path.isfile(src):
            shutil.copy(src, workdir)
        else:
            print(f"Предупреждение: файл {fname} не найден")
            sys.exit(1)

    f = open("workdir.tmp", "w")
    f.write(workdir)
    f.close()


if __name__ == "__main__":
    main()
