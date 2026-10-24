#!/usr/bin/env python3
"""Summarize GCC JSON coverage for project sources, excluding tests/vendor code.

Usage: python3 tests/coverage_summary.py --build-dir /tmp/beren-coverage
Run the instrumented tests before collecting coverage. This reports executable
line coverage of that build, not branch coverage or physical validation.
"""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import tempfile


def summarize(build, repo):
    gcov = shutil.which("gcov")
    if not gcov:
        raise RuntimeError("gcov is not installed")
    graphs = sorted(build.rglob("*.gcno"))
    if not graphs:
        raise RuntimeError("No .gcno files: configure with ENABLE_COVERAGE=ON and build the tests")
    lines = {}
    with tempfile.TemporaryDirectory(prefix="beren-gcov-") as temporary:
        for index, graph in enumerate(graphs):
            directory = Path(temporary) / str(index)
            directory.mkdir()
            result = subprocess.run([gcov, "--json-format", str(graph)], cwd=str(directory),
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode:
                raise RuntimeError("gcov failed for {}: {}".format(graph, result.stderr.decode(errors="replace")))
            for report in directory.glob("*.gcov.json.gz"):
                with gzip.open(str(report), "rt") as stream:
                    data = json.load(stream)
                for item in data["files"]:
                    path = Path(item["file"])
                    if not path.is_absolute():
                        path = Path(data["current_working_directory"]) / path
                    try:
                        relative = str(path.resolve().relative_to(repo / "srcBeren"))
                    except ValueError:
                        continue
                    if relative.startswith(("3rd_party/", "tests/")):
                        continue
                    counts = lines.setdefault(relative, {})
                    for line in item["lines"]:
                        number = line["line_number"]
                        counts[number] = counts.get(number, False) or line["count"] > 0
    if not lines:
        raise RuntimeError("No project coverage found")

    def metric(covered, total):
        return {"covered": covered, "total": total, "percent": round(100 * covered / total, 1) if total else 0}

    groups, files = {}, {}
    for path, counts in lines.items():
        covered, total = sum(counts.values()), len(counts)
        files[path] = metric(covered, total)
        group = groups.setdefault(path.split("/")[0], [0, 0])
        group[0] += covered
        group[1] += total
    return {"groups": {name: metric(*counts) for name, counts in sorted(groups.items())},
            "files": files,
            "total": metric(sum(sum(counts.values()) for counts in lines.values()),
                            sum(len(counts) for counts in lines.values()))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    build = args.build_dir.resolve()
    try:
        report = summarize(build, Path(__file__).resolve().parents[1])
    except (RuntimeError, OSError) as error:
        parser.error(str(error))
    output = args.output or build / "coverage-summary.json"
    output.write_text(json.dumps(report, indent=2))
    for name, counts in report["groups"].items():
        print("{}: {}/{} ({}%)".format(name, counts["covered"], counts["total"], counts["percent"]))
    print("Total: {}% — {}".format(report["total"]["percent"], output))


if __name__ == "__main__":
    main()
