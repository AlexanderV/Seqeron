#!/usr/bin/env python3
"""Single-line mutation runner for the ESA reference kernels (REVIEW-S1a, docs/refactoring/suffix-tree-upgrade/adr/REVIEW-S1a.md).

Applies each mutant of mutants.json to a file of src/SuffixTree/Algorithms/SuffixTree.Esa, builds the Esa test project and runs the stages
default -> fuzz (ESA_FUZZ_SCALE=0.1) -> vectors -> exhaustive until one fails; prints one JSON line per mutant (status KILLED / SURVIVED /
BUILD-FAIL / BAD-PATTERN) and restores the production file afterwards (also on error). Run it on a quiet machine, nothing else editing the tree.

  python3 run_mutants.py [M01 M02 ...]        SEQERON_BENCH_DIR=<regenerated bench dir with data/*.txt> for the Vectors stage
  python3 run_mutants.py --kernel <dir>        mutate another source directory (for example the future optimized kernels); the old/new patterns
                                               of mutants.json must then be adapted to that source.
Mutant fields: id, desc, file (relative to the source directory), old (exactly one occurrence), new, class (killed | equivalent | killed-by-added-test).
"""
import json, os, re, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
SRC = REPO / "src/SuffixTree/Algorithms/SuffixTree.Esa"
PROJ = "tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj"
ENV = dict(os.environ, DOTNET_CLI_TELEMETRY_OPTOUT="1")
STAGES = [
    ("default", ["dotnet", "test", PROJ, "-c", "Release", "--nologo", "--no-build"], {}),
    ("fuzz", ["dotnet", "test", PROJ, "-c", "Release", "--nologo", "--no-build", "--filter", "Category=Fuzz"], {"ESA_FUZZ_SCALE": "0.1"}),
    ("vectors", ["dotnet", "test", PROJ, "-c", "Release", "--nologo", "--no-build", "--filter", "Category=Vectors"], {}),
    ("exhaustive", ["dotnet", "test", PROJ, "-c", "Release", "--nologo", "--no-build", "--filter", "Category=Exhaustive"], {}),
]


def run(cmd, extra=None, timeout=1200):
    env = dict(ENV, **(extra or {}))
    try:
        p = subprocess.run(cmd, cwd=REPO, env=env, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        return 124, "TIMEOUT"


def main(argv):
    src = SRC
    if "--kernel" in argv:
        i = argv.index("--kernel"); src = Path(argv[i + 1]); del argv[i:i + 2]
    mutants = json.loads((HERE / "mutants.json").read_text(encoding="utf-8"))
    for m in mutants:
        if argv and m["id"] not in argv:
            continue
        path = src / m["file"]
        orig = path.read_text(encoding="utf-8")
        res = {"id": m["id"], "desc": m["desc"], "file": m["file"]}
        if orig.count(m["old"]) != 1:
            res["status"] = f"BAD-PATTERN occurrences={orig.count(m['old'])}"
            print(json.dumps(res), flush=True)
            continue
        t0 = time.time()
        try:
            path.write_text(orig.replace(m["old"], m["new"]), encoding="utf-8", newline="")
            rc, out = run(["dotnet", "build", PROJ, "-c", "Release", "--nologo", "-v", "q"])
            if rc != 0:
                res["status"] = "BUILD-FAIL"; res["detail"] = [l for l in out.splitlines() if "error" in l][:3]
            else:
                res["status"] = "SURVIVED"; res["stages"] = {}
                for name, cmd, extra in STAGES:
                    rc, out = run(cmd, extra)
                    tail = [l for l in out.splitlines() if l.startswith(("Passed!", "Failed!"))]
                    res["stages"][name] = tail[-1] if tail else out[-200:]
                    if rc != 0:
                        failed = re.findall(r"^\s*Failed (\S+)", out, re.M)
                        res.update(status="KILLED", killed_by=name, failed_tests=failed[:6], failed_count=len(failed))
                        break
        finally:
            path.write_text(orig, encoding="utf-8", newline="")
        res["seconds"] = round(time.time() - t0)
        print(json.dumps(res), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
