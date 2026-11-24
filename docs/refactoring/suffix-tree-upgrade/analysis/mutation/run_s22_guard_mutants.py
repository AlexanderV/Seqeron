#!/usr/bin/env python3
"""Hand-made mutants of the guard of the checked build of work package S2.2 (prefix s22_): BuildExpectation.cs (the rule) and the live check in BufAndArenaTests.cs. Each mutant is one
textual change, built in the Release configuration (analyzers off) and run against the guard's test classes in the DEFAULT CI configuration (no ESA_EXPECT_BOUNDS):

    dotnet test <project> -c Release --no-build --filter "FullyQualifiedName~BuildExpectationTests|FullyQualifiedName~BufAndArenaTests.BuildMode"

It is KILLED when a test fails (a build error is not a result: the mutant is reported `build-error` and has to be rewritten). Then, for information, the same mutant is run with
ESA_EXPECT_BOUNDS=true in the default build (the dishonest configuration): the evidence says whether the live guard (BufAndArenaTests.BuildMode_...) still fails there, as it must, or
PASSES a dishonest build (the mutant has made the CI guard blind; the tests of the rule are what kill it in the default configuration). Run it in a separate git worktree, never in the working checkout; the source files are restored after every mutant.

usage: run_s22_guard_mutants.py <worktree> <result.tsv>      (environment: MUTANT_BUILD_LIMIT s, default 400)
"""
import os
import re
import signal
import subprocess
import sys

ROOT = os.path.abspath(sys.argv[1])
OUT = sys.argv[2]
BUILD_LIMIT = int(os.environ.get("MUTANT_BUILD_LIMIT", "400"))
PROJ = ROOT + "/tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj"
RULE = "tests/SuffixTree/SuffixTree.Esa.Tests/Safety/BuildExpectation.cs"
LIVE = "tests/SuffixTree/SuffixTree.Esa.Tests/Safety/BufAndArenaTests.cs"
FILTER = "FullyQualifiedName~BuildExpectationTests|FullyQualifiedName~BufAndArenaTests.BuildMode"
FLAGS = ["-p:RunAnalyzers=false", "-p:EnforceCodeStyleInBuild=false", "-p:TreatWarningsAsErrors=false"]
env = dict(os.environ, PATH=os.environ["PATH"] + ":/root/.dotnet", DOTNET_ROOT="/root/.dotnet", DOTNET_CLI_TELEMETRY_OPTOUT="1", DOTNET_NOLOGO="1")
env.pop("ESA_EXPECT_BOUNDS", None)

MUTANTS = [
    ("G01", RULE, 'string.Equals(trimmed, "true", StringComparison.OrdinalIgnoreCase)', 'string.Equals(trimmed, "ture", StringComparison.OrdinalIgnoreCase)', "compares the value to a wrong string (`ture`)"),
    ("G02", RULE, 'string.Equals(trimmed, "true", StringComparison.OrdinalIgnoreCase)', 'string.Equals(trimmed, "true", StringComparison.Ordinal)', "`true` is case sensitive (`True` and `TRUE` are garbage)"),
    ("G03", RULE, '        if (value is null) return null;\n        bool expectChecked;', '        if (value is null || true) return null;\n        bool expectChecked;', "ignores the variable (every value is read as not set)"),
    ("G04", RULE, 'if (expectChecked == boundsChecked) return null;', 'if (expectChecked == true) return null;', "ignores BoundsChecked (only the expectation counts)"),
    ("G05", RULE, 'if (expectChecked == boundsChecked) return null;', 'if (expectChecked != boundsChecked) return null;', "inverts the condition"),
    ("G06", RULE, 'internal const string VariableName = "ESA_EXPECT_BOUNDS";', 'internal const string VariableName = "ESA_EXPECT_BOUND";', "reads another variable"),
    ("G07", RULE, 'Verdict(Environment.GetEnvironmentVariable(VariableName), EsaBuild.BoundsChecked)', 'Verdict(Environment.GetEnvironmentVariable(VariableName), !EsaBuild.BoundsChecked)', "the process verdict inverts the build mode"),
    ("G08", LIVE, 'if (verdict is not null) Assert.Fail(verdict);', 'if (verdict is null) Assert.Fail(verdict);', "the live guard fails when the rule holds (inverted)"),
    ("G09", LIVE, '        if (verdict is not null) Assert.Fail(verdict);\n', '', "the live guard drops the failure (the verdict is read and ignored)"),
    ("G10", RULE, '        else return $"{VariableName} is \'{value}\': it must be unset, \'true\' or \'false\'";', '        else return null;', "garbage is read as not set"),
    ("G11", RULE, '        if (value is null) return null;\n        bool expectChecked;', '        if (string.IsNullOrEmpty(value)) return null;\n        bool expectChecked;', "the empty string is read as not set (a workflow expression that expanded to nothing)"),
    ("G12", RULE, 'string trimmed = value.Trim();', 'string trimmed = value;', "blanks around the value are not ignored"),
    ("G13", RULE, 'else if (string.Equals(trimmed, "false", StringComparison.OrdinalIgnoreCase)) expectChecked = false;', 'else if (string.Equals(trimmed, "false", StringComparison.OrdinalIgnoreCase)) return null;', "`false` never fails (the expectation is only one way)"),
    ("G14", RULE, 'Build with -p:EsaBounds=true (dotnet build ... -p:EsaBounds=true, then dotnet test --no-build), or unset the variable.', 'Build it differently, or unset the variable.', "the failure of a default build does not say what to do"),
    ("G15", RULE, '$"{VariableName} is \'{value}\': it must be unset', '$"the variable is \'{value}\': it must be unset', "the garbage message does not name the variable"),
    ("G16", RULE, 'string.Equals(trimmed, "false", StringComparison.OrdinalIgnoreCase)', 'string.Equals(trimmed, "false", StringComparison.Ordinal)', "`false` is case sensitive"),
    ("G17", RULE, 'Verdict(Environment.GetEnvironmentVariable(VariableName), EsaBuild.BoundsChecked)', 'Verdict(Environment.GetEnvironmentVariable(VariableName), true)', "the process verdict ignores the build of the assembly"),
    ("G18", RULE, '        return expectChecked\n            ? $"{VariableName}=true but SuffixTree.Esa was compiled WITHOUT ESA_BOUNDS: this is a default build, in which the 19 tests that need the checked accessor are skipped (and a skipped test is green). "',
     '        return expectChecked\n            ? $"{VariableName}=true but SuffixTree.Esa was compiled WITHOUT ESA_BOUNDS: this is a default build. "', "the failure of a default build does not say why (skipped tests are green)"),
]


def sh(args, limit, extra_env=None):
    e = dict(env)
    if extra_env:
        e.update(extra_env)
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=e, start_new_session=True)
    try:
        txt, _ = p.communicate(timeout=limit)
        return txt, False
    except subprocess.TimeoutExpired:
        os.killpg(p.pid, signal.SIGKILL)
        txt, _ = p.communicate()
        return txt, True


def build():
    txt, hung = sh(["dotnet", "build", PROJ, "-c", "Release", "--nologo", "-v", "q"] + FLAGS, BUILD_LIMIT)
    if hung:
        return "BUILD-TIMEOUT"
    if re.search(r"error (CS|IDE|CA|S)\d+", txt) or "Build FAILED" in txt:
        m = re.search(r"error ((CS|IDE|CA|S)\d+)[^\n]*", txt)
        return "BUILD-ERROR " + (m.group(0)[:120] if m else "")
    return ""


def test(extra_env=None):
    txt, hung = sh(["dotnet", "test", PROJ, "-c", "Release", "--no-build", "--nologo", "--filter", FILTER], 300, extra_env)
    m = re.search(r"(Passed|Failed)!\s+- Failed:\s+(\d+), Passed:\s+(\d+), Skipped:\s+(\d+), Total:\s+(\d+)", txt)
    failed = re.findall(r"^\s+Failed (\S+)", txt, re.M)
    return (int(m.group(2)) if m else -1, int(m.group(3)) if m else -1, failed, hung)


def main():
    originals = {}
    rows = []
    try:
        err = build()
        if err:
            sys.exit("the unmutated tree does not build: " + err)
        failed, passed, names, hung = test()
        base_red = test({"ESA_EXPECT_BOUNDS": "true"})
        rows.append(("C01", "control-ok" if failed == 0 and passed > 0 else "CONTROL-FAILED", "the unmutated guard, default configuration", f"{passed} passed, {failed} failed"))
        rows.append(("C02", "control-ok" if base_red[0] == 1 else "CONTROL-FAILED", "the unmutated guard, ESA_EXPECT_BOUNDS=true in the default build", f"{base_red[1]} passed, {base_red[0]} failed: {base_red[2]}"))
        for mid, rel, old, new, what in MUTANTS:
            path = os.path.join(ROOT, rel)
            if path not in originals:
                originals[path] = open(path, encoding="utf-8").read()
            text = originals[path]
            if text.count(old) != 1:
                rows.append((mid, "NOT-APPLIED", what, f"the old text occurs {text.count(old)} times"))
                continue
            open(path, "w", encoding="utf-8").write(text.replace(old, new))
            try:
                err = build()
                if err:
                    rows.append((mid, "build-error", what, err))
                    continue
                failed, passed, names, hung = test()
                verdict = "killed" if (failed > 0 or hung) else "SURVIVED"
                evidence = f"{failed} failed of {failed + passed}: " + ", ".join(n.replace('SuffixTree.Esa.Tests.Safety.', '') for n in names[:3]) + (" ..." if len(names) > 3 else "")
                red = test({"ESA_EXPECT_BOUNDS": "true"})
                live = any(n.endswith("BuildMode_IsReported_AndMatchesTheExpectedMode") for n in red[2])
                evidence += f" | default build + ESA_EXPECT_BOUNDS=true: the live guard {'still fails (red)' if live else 'PASSES a dishonest build'}"
                rows.append((mid, verdict, what, evidence))
            finally:
                open(path, "w", encoding="utf-8").write(originals[path])
    finally:
        for path, text in originals.items():
            open(path, "w", encoding="utf-8").write(text)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("id\tresult\tmutant\tevidence\n")
        for r in rows:
            f.write("\t".join(c.replace("\t", " ").replace("\n", " ") for c in r) + "\n")
    mutants = [r for r in rows if not r[0].startswith("C")]
    print(f"{len(mutants)} mutants: {sum(1 for r in mutants if r[1] == 'killed')} killed, {sum(1 for r in mutants if r[1] == 'SURVIVED')} survived, "
          f"{sum(1 for r in mutants if r[1] in ('build-error', 'NOT-APPLIED'))} not run")
    for r in rows:
        if r[1] not in ("killed", "control-ok"):
            print("  ", r)


if __name__ == "__main__":
    main()
