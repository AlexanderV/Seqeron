#!/usr/bin/env python3
"""The orchestrator's own mutants of the S2.2 CI scripts (P1-P18: killed = the Python unit tests of analysis/ci fail) and of the workflow (W1-W16: killed = check_deep_workflow.py
exits non-zero). Hand-made, applied to a scratch copy. usage: python3 -I run_s22_orch_mutants.py ; the table is results_s22_orch.tsv (together with the C# mutants of the guard G1-G13,
which were run by hand against BuildExpectationTests and BufAndArenaTests: see the table).
The first pass of this list left 6 survivors (P2, P15, P16, P17, W13, W14); they are closed by tests / checker rules (commit of the gap closing) and the list is 34 of 34 killed."""
import os, shutil, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, *[".."] * 5))   # the repository root (this file is docs/refactoring/suffix-tree-upgrade/analysis/mutation/...)
MUT = os.path.join(tempfile.gettempdir(), "s22_orch_mut")  # a scratch copy of the CI scripts, the workflow and skips.tsv; the repository is never touched
CI = "docs/refactoring/suffix-tree-upgrade/analysis/ci"
WF = ".github/workflows/suffix-tree-deep.yml"
SK = "tests/SuffixTree/SuffixTree.Esa.Tests/skips.tsv"

def fresh():
    if os.path.exists(MUT): shutil.rmtree(MUT)
    os.makedirs(os.path.join(MUT, CI)); os.makedirs(os.path.join(MUT, ".github/workflows")); os.makedirs(os.path.join(MUT, os.path.dirname(SK)))
    for f in os.listdir(os.path.join(SRC, CI)):
        p = os.path.join(SRC, CI, f)
        if os.path.isfile(p): shutil.copy(p, os.path.join(MUT, CI, f))
    shutil.copy(os.path.join(SRC, WF), os.path.join(MUT, WF)); shutil.copy(os.path.join(SRC, SK), os.path.join(MUT, SK))

PY = [  # python mutants: killed when the unit tests of analysis/ci fail
 ("p1", "skip_accounting.py", "if r.outcome in trx_lib.EXECUTED:\n            counts[\"executed\"] += 1", "if r.outcome == \"Passed\":\n            counts[\"executed\"] += 1"),
 ("p2", "skip_accounting.py", "if r.key in explicit or r.class_name in explicit_class:", "if r.key in explicit or r.class_name in explicit_class or r.message.startswith(\"OneTimeSetUp\"):"),
 ("p3", "skip_accounting.py", "if mode == \"checked\" and row.in_builds == \"default\":", "if mode == \"checked\" and row.in_builds == \"both\":"),
 ("p4", "skip_accounting.py", "            if not partial:\n                failures.append(f\"LISTED BUT ABSENT", "            if True:\n                failures.append(f\"LISTED BUT ABSENT"),
 ("p5", "skip_accounting.py", "elif row.in_builds == \"both\" or (row.in_builds == \"default\" and mode == \"default\"):", "elif row.in_builds == \"both\":"),
 ("p6", "skip_accounting.py", "        if row.kind in (\"explicit\", \"explicit-class\"):\n            if ran:", "        if row.kind in (\"explicit\", \"explicit-class\"):\n            if False:"),
 ("p7", "skip_accounting.py", "        print(f\"{len(report.failures)} discrepancy(ies) between the run and {os.path.basename(args.list_path)}\")\n        return 1", "        print(f\"{len(report.failures)} discrepancy(ies) between the run and {os.path.basename(args.list_path)}\")\n        return 0"),
 ("p8", "skip_accounting.py", "        print(f\"CANNOT ACCOUNT: {e}\", file=sys.stderr)\n        return 2", "        print(f\"CANNOT ACCOUNT: {e}\", file=sys.stderr)\n        return 1"),
 ("p9", "skip_accounting.py", "        if row.test in seen:", "        if False:"),
 ("p10", "skip_accounting.py", "if row.kind != \"ignore\" and row.in_builds != \"both\":", "if False:"),
 ("p11", "trx_lib.py", "    if not results:\n        raise TrxError", "    if False:\n        raise TrxError"),
 ("p12", "trx_lib.py", "EXECUTED = frozenset({\"Passed\", \"Failed\",", "EXECUTED = frozenset({\"Passed\", \"Failed\", \"NotExecuted\","),
 ("p13", "trx_lib.py", "classes[unit.get(\"id\")] = method.get(\"className\") or \"\"", "classes[unit.get(\"id\")] = method.get(\"name\") or \"\""),
 ("p14", "trx_summary.py", "unexpected = [r for r in failed if r.key not in expected]", "unexpected = [r for r in failed if r.key in expected]"),
 ("p15", "trx_summary.py", "failed = [r for r in results if r.outcome in (\"Failed\", \"Error\", \"Timeout\", \"Aborted\")]", "failed = [r for r in results if r.outcome in (\"Failed\",)]"),
 ("p16", "skip_accounting.py", "in_builds = \"default\" if \"checked build\" in r.message or \"EsaBounds\" in r.message else \"both\"", "in_builds = \"both\""),
 ("p17", "skip_accounting.py", "            if r.message.startswith(\"OneTimeSetUp:\"):", "            if False:"),
 ("p18", "skip_accounting.py", "    if args.mode not in MODES:\n        print(f\"--mode must be", "    if False:\n        print(f\"--mode must be"),
]
WFM = [  # workflow mutants: killed when check_deep_workflow.py exits non-zero
 ("w1", "on:\n  pull_request:", "on:\n  pull_request_target:"),
 ("w2", "  push:\n    branches: [master]", "  push:\n    branches: [master, 'claude/**']"),
 ("w3", "permissions:\n  contents: read", "permissions:\n  contents: write"),
 ("w4", "run: dotnet build \"$ESA_TESTS\" -c Release --no-restore -p:EsaBounds=true", "run: dotnet build \"$ESA_TESTS\" -c Release --no-restore"),
 ("w5", "ESA_EXPECT_BOUNDS: 'true'", "ESA_EXPECT_BOUNDS: 'false'"),
 ("w6", "timeout-minutes: 75 ", "timeout-minutes: 7 "),
 ("w7", "      github.event_name == 'pull_request' ||\n      github.event_name == 'push' ||\n", "      github.event_name == 'pull_request' ||\n"),
 ("w8", "  fuzz:\n    name: Fuzz\n", "  fuzz:\n    name: Fuzz\n    continue-on-error: true\n"),
 ("w9", "uses: actions/checkout@v4", "uses: actions/checkout@main"),
 ("w10", "    - cron: '43 2 * * 6'   # weekly, Saturday (UTC): large", "    - cron: '43 2 * * 0'   # weekly, Saturday (UTC): large"),
 ("w11", "ESA_FUZZ_SEED: ${{ steps.seed.outputs.seed }}", "ESA_FUZZ_SEED: ${{ secrets.FUZZ_SEED }}"),
 ("w12", "--filter \"Category=Exhaustive\"", "--filter \"Category=Exhaustive|Category=Fuzz\""),
 ("w13", "      - name: Generate the benchmark inputs (Python 3 standard library only)\n", "      - name: Generate the benchmark inputs (skipped)\n        if: ${{ false }}\n"),
 ("w14", "      SEED_INPUT: ${{ inputs.seed }}", "      SEED_INPUT: ''"),
 ("w15", "--filter \"FullyQualifiedName~FullDifferentialFuzzTests\" --logger \"trx;LogFileName=harness.trx\"", "--filter \"FullyQualifiedName~FuzzTests\" --logger \"trx;LogFileName=harness.trx\""),
 ("w16", "          if [[ \"$SEED_INPUT\" =~ ^[0-9]{1,18}$ ]]; then SEED=\"$SEED_INPUT\"; fi\n", ""),
]

def run(cmd, cwd):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=300)

def main():
    fresh(); base_unit = run(["python3", "-m", "unittest", "discover", "-s", os.path.join(MUT, CI), "-p", "test_*.py"], MUT)
    base_wf = run(["python3", os.path.join(MUT, CI, "check_deep_workflow.py"), os.path.join(MUT, WF)], MUT)
    print("baseline unit rc", base_unit.returncode, "workflow check rc", base_wf.returncode, flush=True)
    if base_unit.returncode or base_wf.returncode: sys.exit(3)
    out = []
    for mid, f, old, new in PY:
        fresh(); p = os.path.join(MUT, CI, f); t = open(p, encoding="utf-8").read()
        if t.count(old) != 1: out.append((mid, "NOAPPLY", t.count(old))); print(mid, "NOAPPLY", t.count(old), flush=True); continue
        open(p, "w", encoding="utf-8").write(t.replace(old, new))
        r = run(["python3", "-m", "unittest", "discover", "-s", os.path.join(MUT, CI), "-p", "test_*.py"], MUT)
        out.append((mid, "KILLED" if r.returncode else "SURVIVED", r.returncode)); print(mid, out[-1][1], flush=True)
    for mid, old, new in WFM:
        fresh(); p = os.path.join(MUT, WF); t = open(p, encoding="utf-8").read()
        if t.count(old) < 1: out.append((mid, "NOAPPLY", 0)); print(mid, "NOAPPLY", flush=True); continue
        open(p, "w", encoding="utf-8").write(t.replace(old, new, 1))
        r = run(["python3", os.path.join(MUT, CI, "check_deep_workflow.py"), p], MUT)
        out.append((mid, "KILLED" if r.returncode else "SURVIVED", r.returncode)); print(mid, out[-1][1], flush=True)
    with open(os.path.join(HERE, "results_s22_orch_py.tsv"), "w") as fh:
        for row in out: fh.write("\t".join(str(c) for c in row) + "\n")
    print("killed", sum(1 for r in out if r[1] == "KILLED"), "survived", sum(1 for r in out if r[1] == "SURVIVED"), "noapply", sum(1 for r in out if r[1] == "NOAPPLY"))
main()
