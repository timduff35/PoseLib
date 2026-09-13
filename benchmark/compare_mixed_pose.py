#!/usr/bin/env python3
"""Extract unmodified reference solvers and build a paired comparison benchmark.

Example:
  python3 benchmark/compare_mixed_pose.py --reference <base-commit> --build build
  build/benchmark/mixed_pose_benchmark 100000 1 generic 2 > generic.csv

Additional CMake arguments can follow --, e.g. -- -DCMAKE_CXX_COMPILER=clang++.
"""
import argparse
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--reference", required=True, help="Git revision containing the original E3Q3 solvers")
parser.add_argument("--build", default="build")
parser.add_argument("cmake_args", nargs=argparse.REMAINDER)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
reference = subprocess.check_output(["git", "rev-parse", "--verify", args.reference + "^{commit}"], cwd=root, text=True).strip()
build = Path(args.build).resolve()
source = build / "mixed_pose_reference"
source.mkdir(parents=True, exist_ok=True)
for name in ("p1p2ll", "p2p1ll"):
    original = subprocess.check_output(["git", "show", f"{reference}:PoseLib/solvers/{name}.cc"], cwd=root, text=True)
    if "re3q3_rotation" not in original:
        raise RuntimeError(f"{reference} does not contain the original E3Q3 {name} solver")
    # Only rename the exported symbol; preserve the complete reference algorithm.
    original = original.replace(f"int {name}(", f"int reference_{name}(")
    (source / f"{name}.cc").write_text(original)
(source / "REVISION").write_text(reference + "\n")
extra = args.cmake_args[1:] if args.cmake_args[:1] == ["--"] else args.cmake_args
subprocess.run(["cmake", "-S", str(root), "-B", str(build), "-DWITH_BENCHMARK=ON",
                "-DCMAKE_BUILD_TYPE=Release", f"-DMIXED_POSE_REFERENCE_DIR={source}", *extra], check=True)
subprocess.run(["cmake", "--build", str(build), "--target", "mixed_pose_benchmark", "-j", "6"], check=True)
print(f"Reference: {reference}\nExecutable: {build / 'benchmark/mixed_pose_benchmark'}")
