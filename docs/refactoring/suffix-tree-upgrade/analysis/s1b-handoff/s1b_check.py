#!/usr/bin/env python3
"""Result checker of the hand-off package s1b-handoff (work package 1b.9f). Reads the form that the scripts filled in (form.json) and prints PASS / FAIL / WARN / N/A / NOT RUN per gate.
The gates are fixed HERE, not read from the form: a form that carries a different gate number is told so, and judged by the numbers below.
  R1  RssAnon of the build as the build samples it itself (every attempt, every phase peak is at most the whole peak)   <= 320 MB (decision of 2026-11-18, A1 relaxed: every n, default parameters, four threads)
  R2  RssAnon of the build as an independent sampler reads /proc/<pid>/status                                          <= 320 MB; a figure taken from `ps` (macOS, no /proc) is RSS, not RssAnon: not judged
  B1  the build finished (exit code 0, the last attempt too) and the container has the size the build reported
  I1  identity of the container: its SHA-256 equals the reference of reference_hashes.tsv for the same text / size / PSI, and equals the container of every other run of the same text in the form
  I2  the input text has the SHA-256 of the reference (when there is one)
  V1  e8-verify --container: exit code 0 and all eleven checks are `pass` (a skipped check is not a pass)
  S1  the Structure level of the reader: pass, or not applicable above the reader's limit (2^31 - 3 symbols), or not run (optional)
  Q1  the query set: every answer of the search over the suffix array and of the engine equals the oracle's
  D1  peak scratch within 13 bytes per symbol (a warning, not a failure)
  P1  the process runs on the architecture it was built for (no emulation)
  X1  the CRC-32C path: the check value, the table path equal to the hardware path (a CPU without the instruction is a warning)
  T1  the test sets: no failed test; a count that differs from the reference head's is a warning
usage: s1b_check.py <form.json> [--reference reference_hashes.tsv] [--require-tests] [--run LABEL]...      s1b_check.py --gates
  --run   judge only the named runs of the form (the others still serve as references of the identity gate); the environment and the test sets are judged always
exit code: 0 PASS, 1 FAIL, 2 the file is not a readable form, 3 INCOMPLETE (nothing failed, a required gate has no data)"""
import json
import os
import sys

RSSANON_GATE_MB = 320
READER_LIMIT_N = 2_147_483_645            # the largest text the reader (Structure level, the engine) takes: SesaVerifier and EsaSuffixTree, ADR-004 1b.7f
QUERY_MIN_TOTAL = 500
SCRATCH_BYTES_PER_SYMBOL = 13.0
FORM_VERSION = 1
VERIFY_CHECKS = ['files', 'meta', 'text', 'input', 'sa_permutation', 'sa_order', 'lcp_values', 'lcp_width', 'stats', 'child', 'sha256']
# the reference head's totals (the orchestrator's runs of 2026-11-22 at head 522f8c5f, plus the 20 tests of 1b.9f in the two Esa builds); a different count is a warning, a failed test is a failure
EXPECTED_TESTS = {
    'esa-default': {'passed': 4846, 'skipped': 19},
    'esa-exhaustive': {'passed': 82, 'skipped': 0},
    'esa-fuzz': {'passed': 39, 'skipped': 0},
    'esa-checked': {'passed': 4847, 'skipped': 18},
    'emproto': {'passed': 535, 'skipped': 0},
    'persistent': {'passed': 672, 'skipped': 0},
}


class Report:
    def __init__(self):
        self.counts = {'PASS': 0, 'FAIL': 0, 'WARN': 0, 'N/A': 0, 'NOT RUN': 0}
        self.incomplete = False

    def line(self, tag, gate, text):
        self.counts[tag] += 1
        print(f'  {tag:<7} {gate} {text}')

    def need(self, gate, text):
        """a required gate that has no data: INCOMPLETE"""
        self.counts['NOT RUN'] += 1
        self.incomplete = True
        print(f'  NOT RUN {gate} {text}')


def num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def load_reference(path):
    table = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            head = fh.readline().rstrip('\n').split('\t')
            for raw in fh:
                row = dict(zip(head, raw.rstrip('\n').split('\t')))
                if row.get('kind'):
                    table[f"{row['kind']}:{row['n']}:{row['psi']}"] = row
    return table


def check_environment(env, rep):
    print('environment')
    if not isinstance(env, dict):
        rep.line('FAIL', 'P1', "missing section 'environment'")
        return
    emulation = str(env.get('emulation', 'unknown'))
    arch, image, proc, osarch = env.get('arch'), env.get('image_arch'), env.get('process_arch'), env.get('os_arch')
    same = {'aarch64': 'arm64', 'arm64': 'arm64', 'x86_64': 'amd64', 'amd64': 'amd64'}
    if emulation != 'no':
        rep.line('FAIL', 'P1', f'emulation: {emulation}. Timings under emulation mean nothing and the arm64 path is not exercised; build and run the image for the Docker daemon\'s architecture')
    elif proc and osarch and proc != osarch:
        rep.line('FAIL', 'P1', f'the .NET process is {proc} on an {osarch} OS')
    elif image not in (None, 'unknown') and arch and same.get(str(arch).lower()) != same.get(str(image).lower()):
        rep.line('FAIL', 'P1', f'the image is for {image} but the kernel reports {arch}')
    else:
        rep.line('PASS', 'P1', f'native: kernel {arch}, image {image}, .NET process {proc}')
    if env.get('crc32c_check_value') != 'ok':
        rep.line('FAIL', 'X1', f"CRC-32C of 123456789 is not E3069283 on this machine ({env.get('crc32c_check_value')})")
    elif env.get('crc32c_table_equals_hardware') is not True:
        rep.line('FAIL', 'X1', 'the CRC-32C of the table path and of the chosen path differ')
    elif env.get('crc32c_hardware') is True:
        rep.line('PASS', 'X1', f"CRC-32C on the CPU instruction ({env.get('crc32c_instruction')}), equal to the table path")
    else:
        rep.line('WARN', 'X1', 'no CRC-32C instruction was taken: the managed table path runs (same values, slower); the arm64 instruction path is not exercised by this run')
    print(f"  info    {env.get('cores')} cores, {env.get('mem_total_bytes')} bytes of memory, cgroup {env.get('cgroup_version')}, {env.get('fs_type')}, .NET {env.get('dotnet')}, RssAnon source: {env.get('rssanon_source')}")


def check_run(label, run, rep, table, others):
    n = run.get('n') if num(run.get('n')) else None
    print(f"run {label} ({run.get('kind')}, n = {n:,}, resume {'on' if run.get('resume') else 'off'})" if n else f'run {label}')
    for section in ('input', 'build', 'rss', 'container', 'verify', 'structure', 'queries'):
        if not isinstance(run.get(section), dict):
            rep.line('FAIL', section[:1].upper(), f"missing section '{section}' in run {label}")
    if n is None:
        rep.line('FAIL', 'B1', f"missing 'n' in run {label}")
        return
    build, rss, cont, ver = run.get('build') or {}, run.get('rss') or {}, run.get('container') or {}, run.get('verify') or {}
    # B1
    if build:
        attempts = build.get('attempts') or []
        last = attempts[-1].get('exit_code') if attempts and isinstance(attempts[-1], dict) else None
        if build.get('exit_code') != 0 or (attempts and last != 0):
            rep.line('FAIL', 'B1', f"the build did not finish: exit code {build.get('exit_code')}, last attempt {last}")
        elif cont and cont.get('bytes') != build.get('container_bytes'):
            rep.line('FAIL', 'B1', f"the container has {cont.get('bytes')} bytes, the build reported {build.get('container_bytes')}")
        else:
            rep.line('PASS', 'B1', f"finished in {len(attempts) or 1} attempt(s), {build.get('total_s_sum')} s inside the build, container {build.get('container_bytes'):,} bytes ({build.get('container_bytes') / n:.2f} per symbol)")
    # R1 / R2
    if rss:
        if rss.get('gate_mb') not in (None, RSSANON_GATE_MB):
            rep.line('WARN', 'R0', f"the form carries the gate {rss.get('gate_mb')} MB; this checker judges by {RSSANON_GATE_MB} MB")
        inproc = [rss.get('in_process_peak_mb')] + [a.get('rssanon_in_process_mb') for a in (build.get('attempts') or []) if isinstance(a, dict)]
        if not num(rss.get('in_process_peak_mb')):
            rep.line('FAIL', 'R1', "missing field 'rss.in_process_peak_mb'")
        else:
            worst = max(v for v in inproc if num(v))
            if worst > RSSANON_GATE_MB:
                rep.line('FAIL', 'R1', f'RssAnon peak {worst:.1f} MB (the build\'s own sampler, 10 ms) is above the gate of {RSSANON_GATE_MB} MB')
            else:
                rep.line('PASS', 'R1', f'RssAnon peak {worst:.1f} MB (the build\'s own sampler, 10 ms) <= {RSSANON_GATE_MB} MB; margin {RSSANON_GATE_MB - worst:.1f} MB')
        if rss.get('external_source') != 'proc':
            rep.line('N/A', 'R2', f"not judged: the independent figure comes from '{rss.get('external_source')}' (RSS with the mapped file pages, not RssAnon); the gate is a Linux RssAnon gate, run it in the Linux container")
        elif not num(rss.get('external_peak_mb')):
            rep.line('FAIL', 'R2', "missing field 'rss.external_peak_mb'")
        else:
            ext = [rss.get('external_peak_mb')] + [a.get('rssanon_external_mb') for a in (build.get('attempts') or []) if isinstance(a, dict)]
            worst = max(v for v in ext if num(v))
            if worst > RSSANON_GATE_MB:
                rep.line('FAIL', 'R2', f'RssAnon peak {worst:.1f} MB (independent sampler of /proc/<pid>/status) is above the gate of {RSSANON_GATE_MB} MB')
            else:
                rep.line('PASS', 'R2', f'RssAnon peak {worst:.1f} MB (independent sampler of /proc/<pid>/status) <= {RSSANON_GATE_MB} MB')
    # I1 / I2
    if cont:
        sha = cont.get('sha256')
        key = f"{run.get('kind')}:{n}:{'psi' if run.get('psi', True) else 'nopsi'}"
        refs = []
        if key in table:
            refs.append(('the reference table', table[key]['container_sha256']))
        if cont.get('reference_sha256'):
            refs.append(('the reference in the form', cont['reference_sha256']))
        refs += [(f'the other run {other}', o['container']['sha256']) for other, o in others if isinstance(o.get('container'), dict) and o['container'].get('sha256')]
        bad = [(src, h) for src, h in refs if h != sha]
        if not sha:
            rep.line('FAIL', 'I1', "missing field 'container.sha256'")
        elif bad:
            rep.line('FAIL', 'I1', f"container SHA-256 {sha[:16]}... differs from {bad[0][0]} ({bad[0][1][:16]}...)")
        elif refs:
            rep.line('PASS', 'I1', f"container SHA-256 {sha[:16]}... equals {' and '.join(src for src, _ in refs)}")
        else:
            rep.line('N/A', 'I1', f"no reference for {key} (a new size): the hash {sha[:16]}... is recorded; a second run, or a resume after a kill, must give the same one")
    inp = run.get('input') or {}
    if inp and run.get('kind') and f"{run.get('kind')}:{n}:psi" in table:
        want = table[f"{run.get('kind')}:{n}:psi"].get('input_sha256')
        if want and inp.get('sha256') != want:
            rep.line('FAIL', 'I2', f"the input text has SHA-256 {str(inp.get('sha256'))[:16]}..., the reference text {want[:16]}...: another text was indexed (the generator differs)")
        elif want:
            rep.line('PASS', 'I2', 'the input text equals the reference text')
    # V1
    if ver:
        if ver.get('status') == 'not_run':
            rep.need('V1', 'e8-verify --container was not run')
        else:
            checks = ver.get('checks') or {}
            notpass = [c for c in VERIFY_CHECKS if checks.get(c) != 'pass']
            if ver.get('exit_code') != 0 or ver.get('all') != 'pass' or notpass:
                rep.line('FAIL', 'V1', f"e8-verify --container: exit code {ver.get('exit_code')}, all = {ver.get('all')}, not a pass: {', '.join(notpass) or '-'}")
            else:
                rep.line('PASS', 'V1', f"e8-verify --container: 11 of 11 checks pass in {ver.get('seconds')} s")
    # S1
    st = run.get('structure')
    if isinstance(st, dict):
        status = st.get('status')
        if status == 'pass':
            if n > READER_LIMIT_N:
                rep.line('WARN', 'S1', 'the Structure level passed above the reader\'s documented limit: the limit in this checker is out of date')
            else:
                rep.line('PASS', 'S1', f"Structure level of the reader: pass in {st.get('seconds')} s, RssAnon peak {st.get('rssanon_peak_mb')} MB")
        elif status == 'fail':
            rep.line('FAIL', 'S1', f"Structure level of the reader: {st.get('message') or 'refused the container'}")
        elif status == 'not_applicable':
            if n > READER_LIMIT_N:
                rep.line('N/A', 'S1', f'the reader takes at most {READER_LIMIT_N:,} symbols (a 32-bit surface, stage S2): not applicable at n = {n:,}; e8-verify --container is the check at this size')
            else:
                rep.line('FAIL', 'S1', f'the reader refused the Structure level at n = {n:,}, which is within its limit of {READER_LIMIT_N:,}: {st.get("message")}')
        elif status == 'not_run':
            rep.line('NOT RUN', 'S1', 'Structure level not run (an optional step)')
        else:
            rep.line('FAIL', 'S1', f"unknown status '{status}'")
    # Q1
    q = run.get('queries')
    if isinstance(q, dict):
        if q.get('status') == 'not_run':
            rep.need('Q1', 'the query set was not run')
        else:
            problems = []
            if not num(q.get('total')) or q['total'] < QUERY_MIN_TOTAL:
                problems.append(f"only {q.get('total')} queries (at least {QUERY_MIN_TOTAL})")
            if q.get('sa_mismatches') != 0:
                problems.append(f"{q.get('sa_mismatches')} answers of the suffix-array search differ from the oracle")
            engine = q.get('engine_status')
            if engine == 'ok' and q.get('engine_mismatches') != 0:
                problems.append(f"{q.get('engine_mismatches')} answers of the engine differ from the oracle")
            if engine == 'not_applicable' and n <= READER_LIMIT_N:
                problems.append(f'the engine refused the container at n = {n:,}, within its limit')
            if engine not in ('ok', 'not_applicable', 'skipped'):
                problems.append(f"engine status '{engine}'")
            if q.get('oracle_present_without_hit') != 0:
                problems.append(f"{q.get('oracle_present_without_hit')} queries cut from the text were not found by the oracle")
            if problems:
                rep.line('FAIL', 'Q1', '; '.join(problems))
            elif engine == 'skipped':
                rep.line('WARN', 'Q1', f"{q.get('total')} queries equal the oracle through the suffix-array search; the engine was switched off")
            else:
                via = 'the suffix-array search and the engine' if engine == 'ok' else 'the suffix-array search (the engine cannot open a container above its 32-bit surface)'
                rep.line('PASS', 'Q1', f"{q.get('total')} queries (count, exists, locate) answered by {via} equal the oracle")
    # D1
    disk = run.get('disk')
    if isinstance(disk, dict) and isinstance(build, dict) and num(build.get('scratch_peak_bytes')):
        per = build['scratch_peak_bytes'] / n
        if per > SCRATCH_BYTES_PER_SYMBOL:
            rep.line('WARN', 'D1', f'peak scratch {per:.2f} bytes per symbol is above the {SCRATCH_BYTES_PER_SYMBOL:.0f} that the pre-flight plans for')
        else:
            rep.line('PASS', 'D1', f"peak scratch {per:.2f} bytes per symbol (<= {SCRATCH_BYTES_PER_SYMBOL:.0f}); sampled peak of everything on the volume {disk.get('total_peak_bytes_sampled'):,} bytes")


def check_tests(tests, rep, require):
    print('tests')
    for name, expected in EXPECTED_TESTS.items():
        t = (tests or {}).get(name)
        if not isinstance(t, dict) or t.get('status') in (None, 'not_run'):
            if require:
                rep.need('T1', f'{name}: not run')
            continue
        failed, passed, skipped = t.get('failed'), t.get('passed'), t.get('skipped')
        if t.get('status') == 'fail' or (num(failed) and failed > 0) or not num(passed):
            rep.line('FAIL', 'T1', f"{name}: {failed} failed, {passed} passed, {skipped} skipped")
        elif passed != expected['passed'] or skipped != expected['skipped']:
            rep.line('WARN', 'T1', f"{name}: 0 failed, {passed} passed, {skipped} skipped (the reference head has {expected['passed']} + {expected['skipped']} skipped: a platform may skip or add tests)")
        else:
            rep.line('PASS', 'T1', f"{name}: 0 failed, {passed} passed, {skipped} skipped in {t.get('seconds')} s")
    for name in (tests or {}):
        if name not in EXPECTED_TESTS:
            rep.line('WARN', 'T1', f"{name}: a set the checker does not know")


def main(argv):
    if '--gates' in argv:
        print(f'RSSANON_GATE_MB = {RSSANON_GATE_MB}\nREADER_LIMIT_N = {READER_LIMIT_N}\nQUERY_MIN_TOTAL = {QUERY_MIN_TOTAL}\nSCRATCH_BYTES_PER_SYMBOL = {SCRATCH_BYTES_PER_SYMBOL}')
        return 0
    args = [a for a in argv if not a.startswith('--')]
    ref_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reference_hashes.tsv')
    if '--reference' in argv:
        ref_path = argv[argv.index('--reference') + 1]
        args = [a for a in args if a != ref_path]
    only = [argv[i + 1] for i, a in enumerate(argv[:-1]) if a == '--run']
    args = [a for a in args if a not in only]
    if len(args) != 1:
        print(__doc__)
        return 2
    try:
        with open(args[0]) as fh:
            form = json.load(fh)
    except (OSError, ValueError) as e:
        print(f'cannot read the form {args[0]}: {e}')
        return 2
    if not isinstance(form, dict) or form.get('form') != 's1b-handoff':
        print('this is not a form of the s1b-handoff package')
        return 2
    if form.get('form_version') != FORM_VERSION:
        print(f"form version {form.get('form_version')}: this checker reads version {FORM_VERSION}")
        return 2
    table = load_reference(ref_path)
    rep = Report()
    print(f"S1B result check of {args[0]} (generated {form.get('generated_utc')}, package commit {form.get('package_commit')})")
    check_environment(form.get('environment'), rep)
    runs = form.get('runs') if isinstance(form.get('runs'), dict) else {}
    if not runs:
        rep.need('B1', 'the form has no run')
    keys = {}
    for label, run in runs.items():
        if isinstance(run, dict):
            keys[label] = f"{run.get('kind')}:{run.get('n')}:{'psi' if run.get('psi', True) else 'nopsi'}"
    if only and [x for x in only if x not in runs]:
        print('the form has no run ' + ', '.join(x for x in only if x not in runs))
        return 2
    for label, run in runs.items():
        if only and label not in only:
            continue
        if not isinstance(run, dict):
            rep.line('FAIL', 'B1', f'run {label} is not an object')
            continue
        others = [(other, o) for other, o in runs.items() if other != label and isinstance(o, dict) and keys.get(other) == keys.get(label)]
        check_run(label, run, rep, table, others)
    check_tests(form.get('tests'), rep, '--require-tests' in argv)
    c = rep.counts
    print(f"SUMMARY: {c['PASS']} PASS, {c['FAIL']} FAIL, {c['WARN']} WARN, {c['N/A']} N/A, {c['NOT RUN']} NOT RUN")
    if c['FAIL']:
        print('RESULT: FAIL')
        return 1
    if rep.incomplete:
        print('RESULT: INCOMPLETE (nothing failed; a required gate has no data)')
        return 3
    print('RESULT: PASS')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
