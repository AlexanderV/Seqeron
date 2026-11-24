#!/usr/bin/env python3
"""Structural check of .github/workflows/suffix-tree-deep.yml (work package S2.2).

actionlint is not available on the machine that wrote the workflow, and nothing can start the workflow before the owner merges it, so this script reads
the YAML (PyYAML, duplicate keys refused) and checks what the policy of PROGRESS-S2.md "Test strategy" promises:

  * the triggers: pull_request and push (master only, never another branch) with exactly the three paths, a nightly and a weekly cron, a manual
    dispatch with the `set` input; no other event;
  * the `if:` of every job, EVALUATED for every scenario (an event, the cron, the input) with a small interpreter of the subset of the GitHub expression
    language that the file uses: which jobs run on a pull request, on a push to master, at night, weekly, and for every value of `set`;
  * permissions (contents: read, nowhere more), the pinned actions (checkout, setup-dotnet, upload-artifact at the majors of dotnet.yml, no other),
    a timeout on every job (within the band that the measured times allow), no secret, no network tool, the concurrency group, the DOTNET_* variables,
    the checked build's flags and variable, continue-on-error only on the `large` job, the filters of the test steps, the seed rule of the harness.

usage: check_deep_workflow.py [path of the workflow]       exit 0 = all checks hold, 1 = at least one failed (each is printed), 2 = cannot run
Only the Python standard library and PyYAML.
"""
import re
import sys

try:
    import yaml
except ImportError:  # pragma: no cover
    print("PyYAML is not installed (python3 -c 'import yaml'); this check needs it", file=sys.stderr)
    sys.exit(2)

DEFAULT_PATH = ".github/workflows/suffix-tree-deep.yml"

PATHS = ["src/SuffixTree/**", "tests/SuffixTree/**", ".github/workflows/suffix-tree-deep.yml"]
NIGHTLY = "17 1 * * *"
WEEKLY = "43 2 * * 6"
NIGHTLY_JOBS = ["checked", "exhaustive", "fuzz", "harness", "vectors"]
ALL_JOBS = NIGHTLY_JOBS + ["large"]
SET_OPTIONS = ["nightly", "all", "checked", "exhaustive", "fuzz", "harness", "vectors", "large"]
ALLOWED_ACTIONS = {"actions/checkout@v4", "actions/setup-dotnet@v4", "actions/upload-artifact@v4"}
# timeout band per job, in minutes: not below 2x the measured time on 4 vCPU plus the compile (a hosted runner may have 2 vCPU), not above 4 hours
TIMEOUT_BAND = {"checked": (45, 240), "exhaustive": (50, 240), "fuzz": (45, 240), "harness": (40, 240), "vectors": (15, 240), "large": (15, 240)}
FILTERS = {"exhaustive": "Category=Exhaustive", "fuzz": "Category=Fuzz", "vectors": "Category=Vectors", "large": "Category=Large",
           "harness": "FullyQualifiedName~FullDifferentialFuzzTests"}
FORBIDDEN_RUN = re.compile(r"(?m)(\b(curl|wget|pip|pip3|npm|npx|apt|apt-get|docker|nuget)\b|git\s+(clone|push|fetch|pull)|^\s*gh\s)")
FORBIDDEN_RAW = re.compile(r"secrets\.|GITHUB_TOKEN|\$\{\{\s*vars\.")

failures = []


def check(ok, message):
    if not ok:
        failures.append(message)
    return ok


class UniqueKeyLoader(yaml.SafeLoader):
    """SafeLoader that refuses a duplicate key (PyYAML keeps the last one silently)."""


def _construct_mapping(loader, node, deep=False):
    seen = set()
    for key_node, _ in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in seen:
            raise yaml.constructor.ConstructorError(None, None, f"duplicate key {key!r}", key_node.start_mark)
        seen.add(key)
    return yaml.SafeLoader.construct_mapping(loader, node, deep)


UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping)


# ---------------------------------------------------------------------------------------------------------------------------------------------
# the interpreter of the expression subset: == != && || ! ( ) , string literals, dotted names, contains(), fromJSON(), format()
# ---------------------------------------------------------------------------------------------------------------------------------------------

TOKEN = re.compile(r"\s*(?:(==|!=|&&|\|\||!|\(|\)|,)|'((?:[^']|'')*)'|([A-Za-z_][A-Za-z0-9_.\-]*)|(\d+))")


def tokenize(text):
    pos = 0
    tokens = []
    text = text.strip()
    while pos < len(text):
        m = TOKEN.match(text, pos)
        if not m or m.end() == pos:
            raise ValueError(f"cannot tokenize {text[pos:pos + 20]!r}")
        if m.group(1):
            tokens.append(("op", m.group(1)))
        elif m.group(2) is not None:
            tokens.append(("str", m.group(2).replace("''", "'")))
        elif m.group(3):
            tokens.append(("name", m.group(3)))
        else:
            tokens.append(("num", int(m.group(4))))
        pos = m.end()
    return tokens


class Evaluator:
    def __init__(self, tokens, context):
        self.t = tokens
        self.i = 0
        self.ctx = context

    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else (None, None)

    def take(self, kind=None, value=None):
        tok = self.peek()
        if tok[0] is None or (kind and tok[0] != kind) or (value is not None and tok[1] != value):
            raise ValueError(f"unexpected {tok} at token {self.i}")
        self.i += 1
        return tok

    def parse(self):
        value = self.or_()
        if self.i != len(self.t):
            raise ValueError(f"trailing tokens {self.t[self.i:]}")
        return value

    def or_(self):
        left = self.and_()
        while self.peek() == ("op", "||"):
            self.take()
            right = self.and_()
            left = left if truthy(left) else right
        return left

    def and_(self):
        left = self.eq()
        while self.peek() == ("op", "&&"):
            self.take()
            right = self.eq()
            left = right if truthy(left) else left
        return left

    def eq(self):
        left = self.unary()
        while self.peek() in (("op", "=="), ("op", "!=")):
            op = self.take()[1]
            right = self.unary()
            same = str(left if left is not None else "").lower() == str(right if right is not None else "").lower()
            left = same if op == "==" else not same
        return left

    def unary(self):
        if self.peek() == ("op", "!"):
            self.take()
            return not truthy(self.unary())
        return self.primary()

    def primary(self):
        kind, value = self.peek()
        if kind == "op" and value == "(":
            self.take()
            inner = self.or_()
            self.take("op", ")")
            return inner
        if kind == "str":
            self.take()
            return value
        if kind == "num":
            self.take()
            return value
        if kind == "name":
            self.take()
            if self.peek() == ("op", "("):
                self.take()
                args = []
                if self.peek() != ("op", ")"):
                    args.append(self.or_())
                    while self.peek() == ("op", ","):
                        self.take()
                        args.append(self.or_())
                self.take("op", ")")
                return self.call(value, args)
            if value in ("true", "false"):
                return value == "true"
            return self.ctx.get(value)  # an unknown or unset name is null, as in GitHub
        raise ValueError(f"unexpected token {kind} {value}")

    @staticmethod
    def call(name, args):
        import json
        lowered = name.lower()
        if lowered == "fromjson":
            return json.loads(args[0])
        if lowered == "contains":
            haystack, needle = args
            if isinstance(haystack, list):
                return any(str(x).lower() == str(needle if needle is not None else "").lower() for x in haystack)
            return str(needle if needle is not None else "").lower() in str(haystack if haystack is not None else "").lower()
        if lowered == "format":
            out = args[0]
            for k, a in enumerate(args[1:]):
                out = out.replace("{" + str(k) + "}", str(a))
            return out
        raise ValueError(f"function {name} is not in the subset")


def truthy(value):
    return bool(value) and value != "0"


def evaluate(expression, context):
    text = expression.strip()
    m = re.fullmatch(r"\$\{\{(.*)\}\}", text, re.S)
    if m:
        text = m.group(1)
    return Evaluator(tokenize(text), context).parse()


def scenarios():
    """(label, context, the jobs that must run)."""
    def ctx(event, schedule="", set_=None):
        c = {"github.event_name": event, "github.event.schedule": schedule}
        if set_ is not None:
            c["inputs.set"] = set_
        return c

    yield "pull_request", ctx("pull_request"), ["checked"]
    yield "push (master)", ctx("push"), ["checked"]
    yield "schedule nightly", ctx("schedule", NIGHTLY), NIGHTLY_JOBS
    yield "schedule weekly", ctx("schedule", WEEKLY), ["large"]
    yield "schedule of another cron", ctx("schedule", "0 0 * * *"), []
    for event in ("issues", "release", "workflow_run", "repository_dispatch", "pull_request_target"):
        yield f"event {event}", ctx(event), []
    yield "dispatch nightly", ctx("workflow_dispatch", "", "nightly"), NIGHTLY_JOBS
    yield "dispatch all", ctx("workflow_dispatch", "", "all"), ALL_JOBS
    for job in ALL_JOBS:
        yield f"dispatch {job}", ctx("workflow_dispatch", "", job), [job]
    yield "dispatch with an empty set", ctx("workflow_dispatch", "", ""), []
    yield "dispatch with an unknown set", ctx("workflow_dispatch", "", "everything"), []


# ---------------------------------------------------------------------------------------------------------------------------------------------
# the checks
# ---------------------------------------------------------------------------------------------------------------------------------------------

def steps_of(job):
    return job.get("steps") or []


def run_text(job):
    return "\n".join(str(s.get("run", "")) for s in steps_of(job))


def main(path):
    with open(path, encoding="utf-8") as f:
        raw = f.read()
    try:
        doc = yaml.load(raw, Loader=UniqueKeyLoader)
    except yaml.YAMLError as e:
        print(f"FAIL: the YAML does not parse: {e}")
        return 1
    if not check(isinstance(doc, dict), "the document is not a mapping"):
        return report()

    on = doc.get("on", doc.get(True))  # PyYAML (YAML 1.1) reads the key `on` as the boolean True
    check(isinstance(on, dict), "no `on:` mapping")
    on = on if isinstance(on, dict) else {}
    check(sorted(on) == ["pull_request", "push", "schedule", "workflow_dispatch"], f"the events are {sorted(map(str, on))}: exactly pull_request, push, schedule, workflow_dispatch")
    pr = on.get("pull_request") or {}
    push = on.get("push") or {}
    check(pr.get("paths") == PATHS, f"pull_request.paths is {pr.get('paths')!r}, expected {PATHS!r}")
    check(push.get("paths") == PATHS, f"push.paths is {push.get('paths')!r}, expected {PATHS!r}")
    check(push.get("branches") == ["master"], f"push.branches is {push.get('branches')!r}: the push trigger must name master and nothing else (the owner's CI minutes)")
    check(set(push) <= {"branches", "paths"}, f"push has the keys {sorted(push)}: only branches and paths (no tags, no branches-ignore, no paths-ignore)")
    check(set(pr) <= {"paths"}, f"pull_request has the keys {sorted(pr)}: only paths")
    crons = [s.get("cron") for s in (on.get("schedule") or [])]
    check(crons == [NIGHTLY, WEEKLY], f"the crons are {crons!r}, expected {[NIGHTLY, WEEKLY]!r}")
    dispatch = on.get("workflow_dispatch") or {}
    set_input = (dispatch.get("inputs") or {}).get("set") or {}
    check(set_input.get("options") == SET_OPTIONS, f"workflow_dispatch input `set` options are {set_input.get('options')!r}, expected {SET_OPTIONS!r}")
    check(set_input.get("default") == "nightly", "the default of the input `set` is not nightly")
    check("seed" in (dispatch.get("inputs") or {}), "no `seed` input for the harness")

    check(doc.get("permissions") == {"contents": "read"}, f"permissions is {doc.get('permissions')!r}, expected {{'contents': 'read'}}")
    check(isinstance(doc.get("concurrency"), dict) and "group" in doc["concurrency"], "no concurrency group")
    env = doc.get("env") or {}
    for var in ("CI", "DOTNET_NOLOGO", "DOTNET_CLI_TELEMETRY_OPTOUT"):
        check(var in env, f"env.{var} is missing (as in dotnet.yml)")
    check(not any(k.startswith("DOTNET_GC") for k in env), "a DOTNET_GC* variable is set")

    jobs = doc.get("jobs") or {}
    check(sorted(jobs) == sorted(ALL_JOBS), f"the jobs are {sorted(jobs)}, expected {sorted(ALL_JOBS)}")

    for name, job in jobs.items():
        check(job.get("runs-on") == "ubuntu-latest", f"{name}: runs-on is {job.get('runs-on')!r}")
        low, high = TIMEOUT_BAND.get(name, (1, 240))
        timeout = job.get("timeout-minutes")
        check(isinstance(timeout, int) and low <= timeout <= high, f"{name}: timeout-minutes is {timeout!r}, expected an integer in [{low}, {high}]")
        check("permissions" not in job or job["permissions"] == {"contents": "read"}, f"{name}: permissions override {job.get('permissions')!r}")
        check(("continue-on-error" in job) == (name == "large") and (name != "large" or job.get("continue-on-error") is True),
              f"{name}: continue-on-error is {job.get('continue-on-error')!r}; only `large` has it, and it is true")
        for step in steps_of(job):
            if "uses" in step:
                check(step["uses"] in ALLOWED_ACTIONS, f"{name}: the action {step['uses']!r} is not one of {sorted(ALLOWED_ACTIONS)}")
        uploads = [s for s in steps_of(job) if str(s.get("uses", "")).startswith("actions/upload-artifact")]
        check(len(uploads) == 1, f"{name}: {len(uploads)} upload-artifact steps (expected one, for the .trx)")
        for u in uploads:
            check("trx" in str((u.get("with") or {}).get("path", "")), f"{name}: the upload does not carry the .trx files")
            check("cancelled" in str(u.get("if", "")), f"{name}: the upload does not run after a failure (if: !cancelled())")
            check((u.get("with") or {}).get("name") == f"trx-{name}", f"{name}: the artifact is named {(u.get('with') or {}).get('name')!r}, expected trx-{name} (v4 artifact names must be unique in a run)")
        summaries = [s for s in steps_of(job) if "trx_summary.py" in str(s.get("run", ""))]
        check(len(summaries) == 1 and "GITHUB_STEP_SUMMARY" in str(summaries[0].get("run", "")), f"{name}: the job summary step is missing")
        check(all("cancelled" in str(s.get("if", "")) for s in summaries), f"{name}: the summary step does not run after a failure")
        text = run_text(job)
        check("dotnet restore" in text and "dotnet build" in text and "dotnet test" in text, f"{name}: restore, build and test are not all there")
        check("--no-build" in text, f"{name}: the test step rebuilds (a second build can drop the build flags)")
        check("--results-directory" in text and "trx;LogFileName=" in text, f"{name}: the .trx is not written to a known directory")
        if name in FILTERS:
            test_steps = [st for st in steps_of(job) if "dotnet test" in str(st.get("run", "")) and "--no-build" in str(st.get("run", "")) and "--logger" in str(st.get("run", ""))
                          and "trx_summary" not in str(st.get("run", ""))]
            check(len(test_steps) == 1 and f'--filter "{FILTERS[name]}"' in test_steps[0]["run"], f"{name}: the test step does not run with --filter \"{FILTERS[name]}\"")
        if name != "checked":
            check("EsaBounds" not in text, f"{name}: builds with EsaBounds (only the checked job does)")

    checked = jobs.get("checked") or {}
    check((checked.get("env") or {}).get("ESA_EXPECT_BOUNDS") == "true", "checked: env ESA_EXPECT_BOUNDS must be the string 'true'")
    build_steps = [s for s in steps_of(checked) if "dotnet build" in str(s.get("run", ""))]
    check(len(build_steps) == 1 and "-p:EsaBounds=true" in build_steps[0]["run"], "checked: the build step lacks -p:EsaBounds=true")
    accounting = [s for s in steps_of(checked) if "skip_accounting.py" in str(s.get("run", ""))]
    check(len(accounting) == 1 and "--mode checked" in accounting[0]["run"] and "skips.tsv" in accounting[0]["run"] and "cancelled" in str(accounting[0].get("if", "")),
          "checked: the skip accounting step is missing, not in checked mode, not on skips.tsv or not run after a failure")
    selftest = [s for s in steps_of(checked) if "unittest discover" in str(s.get("run", ""))]
    check(len(selftest) == 1 and "test_*.py" in selftest[0]["run"], "checked: no self-test step of the CI scripts (python3 -m unittest discover ... test_*.py)")
    for name in jobs:
        if name != "checked":
            check("ESA_EXPECT_BOUNDS" not in str(jobs[name].get("env", "")), f"{name}: sets ESA_EXPECT_BOUNDS")

    harness = jobs.get("harness") or {}
    htext = "\n".join(str(s.get("run", "")) for s in steps_of(harness)) + "\n" + "\n".join(str(s.get("env", "")) for s in steps_of(harness))
    check("date -u +%Y%m%d" in htext, "harness: the seed is not derived from the UTC date")
    seed_steps = [s for s in steps_of(harness) if s.get("id") == "seed"]
    check(len(seed_steps) == 1 and "GITHUB_STEP_SUMMARY" in str(seed_steps[0].get("run", "")) and "$SEED" in seed_steps[0]["run"], "harness: no step `seed` that writes the seed to the job summary")
    test_steps = [s for s in steps_of(harness) if "dotnet test" in str(s.get("run", "")) and "--no-build" in str(s.get("run", ""))]
    henv = (test_steps[0].get("env") if test_steps else None) or {}
    check(henv.get("ESA_FUZZ_BUDGET_S") == "600", f"harness: budget is {henv.get('ESA_FUZZ_BUDGET_S')!r}, expected '600'")
    check("steps.seed.outputs.seed" in str(henv.get("ESA_FUZZ_SEED", "")), "harness: ESA_FUZZ_SEED does not come from the seed step")
    check(henv.get("ESA_FUZZ_HANG_S") == "900", "harness: ESA_FUZZ_HANG_S is not 900")
    check(re.search(r"\^\[0-9\]\{1,18\}\$", htext) is not None, "harness: the `seed` input is not validated as digits (shell injection)")
    check("inputs.seed" not in "\n".join(str(s.get("run", "")) for s in steps_of(harness)), "harness: inputs.seed is expanded inside a run script (use an env variable)")
    check("inputs.seed" in str((harness.get("env") or {}).get("SEED_INPUT", "")), "harness: the `seed` input does not reach the seed step (SEED_INPUT must be ${{ inputs.seed }})")
    # no step may be switched off or made conditional in a way the file's own policy does not name: the only `if` of a step is `!cancelled()` (summary, upload, accounting)
    for job_name, job in jobs.items():
        for step in steps_of(job):
            if "if" in step:
                check(re.sub(r"\s+", "", str(step["if"])) == "${{!cancelled()}}", f"{job_name}: step {step.get('name')!r} has `if: {step['if']}` (the only step condition allowed is ${{{{ !cancelled() }}}})")

    vectors = jobs.get("vectors") or {}
    vtext = run_text(vectors)
    gen = vtext.find("gen_inputs.py")
    check(gen >= 0 and gen < vtext.find("dotnet test"), "vectors: gen_inputs.py does not run before the tests")
    check("SEQERON_BENCH_DIR=$RUNNER_TEMP/bench" in vtext and "GITHUB_ENV" in vtext, "vectors: SEQERON_BENCH_DIR is not set (through GITHUB_ENV: the runner context is not available in a job-level env)")
    check(vtext.find("SEQERON_BENCH_DIR=") < gen, "vectors: SEQERON_BENCH_DIR is set after the generator")
    check(not any("runner." in str(j.get("env", "")) for j in jobs.values()) and "runner." not in str(doc.get("env", "")),
          "a job-level or workflow-level env uses the runner context (not available there)")

    # the `if:` of every job, evaluated
    for label, context, expected in scenarios():
        running = []
        for name, job in jobs.items():
            condition = job.get("if")
            if condition is None:
                running.append(name)
                continue
            try:
                value = evaluate(str(condition), context)
            except ValueError as e:
                check(False, f"{name}: cannot evaluate the if: {e}")
                continue
            if truthy(value):
                running.append(name)
        check(sorted(running) == sorted(expected), f"scenario {label}: jobs {sorted(running)} would run, expected {sorted(expected)}")

    # the literals of the crons in the `if:` are the ones of the schedule
    all_ifs = " ".join(str(j.get("if", "")) for j in jobs.values())
    for cron in (NIGHTLY, WEEKLY):
        check(cron in all_ifs, f"no job `if:` mentions the cron {cron!r}")
    check(set(re.findall(r"github\.event\.schedule == '([^']*)'", all_ifs)) == {NIGHTLY, WEEKLY}, "the cron literals of the `if:` differ from the schedule")

    # the text rules
    check(FORBIDDEN_RAW.search(raw) is None, "the workflow mentions a secret, the token or a repository variable")
    all_run = "\n".join(run_text(j) for j in jobs.values())
    forbidden = FORBIDDEN_RUN.search(all_run)
    check(forbidden is None, f"a run script contains {forbidden.group(0)!r}: no network step except dotnet restore and setup-dotnet" if forbidden else "")
    check("S2.8" in raw, "the comment on the continue-on-error of `large` does not name S2.8")
    check("dotnet.yml" in raw and "SuffixTree.EmProto.Tests" in raw and "SuffixTree.Persistent.Tests" in raw, "the comment that EmProto and Persistent are in the default CI is missing")
    check("default branch" in raw.lower(), "the comment that a schedule runs only from the default branch is missing")
    return report()


def report():
    if failures:
        for f in failures:
            print("FAIL:", f)
        print(f"{len(failures)} check(s) failed")
        return 1
    print("all checks hold")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_PATH))
