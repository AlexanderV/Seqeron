#!/usr/bin/env python3
"""Fills the form of the hand-off package (work package 1b.9f): reads what the scripts left in the data directory (the rows of the program, the verifier's rows, the sampler's summaries, the hashes, the
test results) and writes form.json (what the checker reads) and form.md (the same numbers for a human). It computes nothing that matters: every gate is judged by s1b_check.py. The files are small
(tens of kilobytes); the owner sends form.json back by hand; nothing in it names a person, a host, a user or a folder (paths are not copied; the kernel release and the .NET version are).
usage: s1b_form.py <data directory> [--out form.json] [--md form.md] [--commit HASH]"""
import glob
import json
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FORM_VERSION = 1
VERIFY_CHECKS = ['files', 'meta', 'text', 'input', 'sa_permutation', 'sa_order', 'lcp_values', 'lcp_width', 'stats', 'child', 'sha256']
PHASES = ['scan', 'create', 'pack', 'ph1', 'ph23', 'ph45', 'bwt', 'psi', 'crc', 'seal']


def read_text(path, default=''):
    try:
        with open(path, errors='replace') as fh:
            return fh.read()
    except OSError:
        return default


def read_json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def read_kv(path):
    out = {}
    for line in read_text(path).splitlines():
        k, _, v = line.partition('\t')
        if k:
            out[k] = v
    return out


def rows(path):
    """metric -> (value, unit, notes) of the last row of each metric; rows of the Bench (R op metric value unit notes) and of the verifier (R E8 case metric value unit notes)"""
    out = {}
    for line in read_text(path).splitlines():
        f = line.split('\t')
        if len(f) < 5 or f[0] != 'R':
            continue
        if f[1] == 'E8' and len(f) >= 7:
            out[f[3]] = (f[4], f[5], f[6])
        else:
            out[f[2]] = (f[3], f[4] if len(f) > 4 else '', f[5] if len(f) > 5 else '')
    return out


def num(v, cast=float):
    try:
        return cast(v)
    except (TypeError, ValueError):
        return None


def value(r, metric, cast=float):
    return num(r[metric][0], cast) if metric in r else None


def build_section(runroot):
    attempts, phases = [], {p: 0.0 for p in PHASES}
    agg = {'counted_read_bytes': 0, 'counted_write_bytes': 0, 'scratch_peak_bytes': 0, 'arena_peak_bytes': 0}
    container_bytes = lcp_kind = overflow = crc = None
    verdict = None
    ext_peaks, in_peaks, sources, frees, scratch_s, total_s, free_start = [], [], set(), [], 0, 0, None
    for out in sorted(glob.glob(os.path.join(runroot, 'build.*.out')), key=lambda p: int(re.findall(r'build\.(\d+)\.out', p)[0])):
        k = re.findall(r'build\.(\d+)\.out', out)[0]
        r = rows(out)
        rc = num(read_text(os.path.join(runroot, f'build.{k}.rc')).strip(), int)
        samp = read_json(os.path.join(runroot, f'build.{k}.sampler.json')) or {}
        att = {'exit_code': rc if rc is not None else -1, 'total_s': value(r, 'total_s'), 'wall_s': value(r, 'wall_s') or num(read_text(os.path.join(runroot, f'build.{k}.wall_s')).strip()),
               'skipped': (r.get('skipped') or ('', '', ''))[0], 'state_found': (r.get('state_found') or ('', '', ''))[0], 'rssanon_in_process_mb': value(r, 'rssanon_peak_mb'),
               'rssanon_external_mb': samp.get('rssanon_peak_mb'), 'samples': samp.get('samples'), 'rssanon_phase_mb': {p: value(r, f'rssanon_{p}_mb') for p in ('ph0', 'ph1', 'ph23', 'ph45', 'ph6', 'seal') if f'rssanon_{p}_mb' in r}}
        attempts.append(att)
        for p in PHASES:
            phases[p] += value(r, f'{p}_s') or 0.0
        for key, metric in (('counted_read_bytes', 'counted_read'), ('counted_write_bytes', 'counted_write')):
            agg[key] += value(r, metric, int) or 0
        for key, metric in (('scratch_peak_bytes', 'scratch_peak'), ('arena_peak_bytes', 'arena_peak')):
            agg[key] = max(agg[key], value(r, metric, int) or 0)
        if 'file_size' in r:
            container_bytes = value(r, 'file_size', int)
            m = re.search(r'LCP (\w+)', r['file_size'][2])
            lcp_kind = m.group(1) if m else lcp_kind
        overflow = value(r, 'overflow_pairs', int) if 'overflow_pairs' in r else overflow
        crc = (r.get('crc_hardware') or (crc, '', ''))[0]
        if 'rssanon_peak_mb' in r:
            m = re.search(r'verdict (\w+)', r['rssanon_peak_mb'][2])
            verdict = m.group(1) if m else verdict
        if att['rssanon_in_process_mb'] is not None:
            in_peaks.append(att['rssanon_in_process_mb'])
        if samp.get('rssanon_peak_mb') is not None:
            ext_peaks.append(samp['rssanon_peak_mb'])
        if samp.get('source'):
            sources.add(samp['source'])
        scratch_s = max(scratch_s, samp.get('scratch_peak_bytes') or 0)
        total_s = max(total_s, samp.get('total_peak_bytes') or 0)
        if samp.get('free_min_bytes') is not None:
            frees.append(samp['free_min_bytes'])
        if free_start is None:
            free_start = samp.get('free_start_bytes')
    finished = [a for a in attempts if a['exit_code'] == 0]
    build = {'exit_code': attempts[-1]['exit_code'] if attempts else None, 'attempts': attempts, 'phases_s': phases,
             'total_s_sum': round(sum(a['total_s'] or 0 for a in attempts), 3), 'wall_s_sum': round(sum(a['wall_s'] or 0 for a in attempts), 3), 'container_bytes': container_bytes,
             'lcp_kind': lcp_kind, 'overflow_pairs': overflow, 'crc_hardware': crc, 'finished_attempts': len(finished), **agg}
    rss = {'gate_mb': 320, 'in_process_peak_mb': max(in_peaks) if in_peaks else None, 'external_peak_mb': max(ext_peaks) if ext_peaks else None,
           'external_source': ('ps-rss' if 'ps-rss' in sources else 'proc') if sources else None, 'verdict': verdict}
    disk = {'scratch_peak_bytes_sampled': scratch_s, 'total_peak_bytes_sampled': total_s, 'free_min_bytes': min(frees) if frees else None, 'free_start_bytes': free_start}
    return build, rss, disk


def verify_section(runroot):
    path = os.path.join(runroot, 'verify.out')
    if not os.path.exists(path):
        return {'status': 'not_run'}
    r = rows(path)
    rc = num(read_text(os.path.join(runroot, 'verify.rc')).strip(), int)
    checks = {c: (r.get(f'verify_{c}') or ('missing', '', ''))[0] for c in VERIFY_CHECKS}
    facts = {k[len('verify_fact_'):]: v[0] for k, v in r.items() if k.startswith('verify_fact_')}
    return {'status': 'pass' if rc == 0 and (r.get('verify_all') or ('',))[0] == 'pass' else 'fail', 'exit_code': rc if rc is not None else -1, 'seconds': value(r, 'verify_wall_s'),
            'all': (r.get('verify_all') or ('missing',))[0], 'checks': checks, 'facts': facts}


def structure_section(runroot):
    path = os.path.join(runroot, 'structure.out')
    if not os.path.exists(path):
        return {'status': 'not_run'}
    r = rows(path)
    result = (r.get('open_result') or ('missing', '', ''))
    status = {'ok': 'pass', 'not_applicable': 'not_applicable', 'corrupt': 'fail', 'error': 'error'}.get(result[0], 'fail')
    return {'status': status, 'seconds': value(r, 'open_seconds'), 'structure_s': value(r, 'open_structure_s'), 'sections_s': value(r, 'open_sections_s'),
            'rssanon_peak_mb': value(r, 'open_rssanon_peak_mb'), 'message': result[2]}


def queries_section(runroot):
    path = os.path.join(runroot, 'queries.out')
    if not os.path.exists(path):
        return {'status': 'not_run'}
    r = rows(path)
    return {'status': 'pass' if (r.get('q_result') or ('',))[0] == 'PASS' else 'fail', 'total': value(r, 'q_total', int), 'present': value(r, 'q_present', int), 'absent': value(r, 'q_absent', int),
            'sa_mismatches': value(r, 'q_sa_mismatches', int), 'engine_status': (r.get('q_engine_status') or ('missing',))[0], 'engine_mismatches': value(r, 'q_engine_mismatches', int),
            'oracle_present_without_hit': value(r, 'q_oracle_present_without_hit', int), 'oracle_hits_total': value(r, 'q_oracle_hits_total', int), 'set_sha256': (r.get('q_set_sha256') or ('',))[0],
            'oracle_scan_s': value(r, 'q_oracle_scan_s'), 'sa_s': value(r, 'q_sa_s'), 'engine_s': value(r, 'q_engine_s'), 'locate_skipped_over_cap': value(r, 'q_locate_skipped_over_cap', int)}


def run_section(data, label, ref):
    runroot = os.path.join(data, 'runs', label)
    meta = read_json(os.path.join(runroot, 'run.json')) or {}
    sums = {}
    for line in read_text(os.path.join(runroot, 'input.sums')).splitlines():
        k, _, v = line.partition('=')
        sums[k] = v
    sha_line = read_text(os.path.join(runroot, 'container.sha256')).split()
    build, rss, disk = build_section(runroot)
    run = {'kind': meta.get('kind'), 'n': meta.get('n'), 'threads': meta.get('threads'), 'psi': meta.get('psi', True), 'resume': meta.get('resume'), 'fsync': meta.get('fsync'),
           'input': {'sha256': sums.get('sha256'), 'generator': sums.get('generator') or meta.get('generator'), 'bytes': num(sums.get('n'), int)},
           'build': build, 'rss': rss, 'disk': disk,
           'container': {'bytes': num(sha_line[1], int) if len(sha_line) > 1 else None, 'sha256': sha_line[0] if sha_line else None,
                         'reference_key': f"{meta.get('kind')}:{meta.get('n')}:{'psi' if meta.get('psi', True) else 'nopsi'}", 'reference_sha256': None},
           'verify': verify_section(runroot), 'structure': structure_section(runroot), 'queries': queries_section(runroot)}
    row = ref.get(run['container']['reference_key'])
    if row:
        run['container']['reference_sha256'] = row['container_sha256']
    if not sha_line:
        del run['container']
    return run


def load_reference(path):
    table = {}
    try:
        with open(path) as fh:
            head = fh.readline().rstrip('\n').split('\t')
            for raw in fh:
                r = dict(zip(head, raw.rstrip('\n').split('\t')))
                if r.get('kind'):
                    table[f"{r['kind']}:{r['n']}:{r['psi']}"] = r
    except OSError:
        pass
    return table


def environment_section(data):
    kv = read_kv(os.path.join(data, 'env', 'env.tsv'))
    h = rows(os.path.join(data, 'env', 'handoff.out'))
    yes = lambda m: (h.get(m) or ('',))[0] == 'yes'  # noqa: E731
    env = {'arch': kv.get('arch'), 'image_arch': kv.get('image_arch'), 'process_arch': (h.get('env_process_arch') or (None,))[0], 'os_arch': (h.get('env_os_arch') or (None,))[0],
           'emulation': kv.get('emulation', 'unknown'), 'kernel': kv.get('kernel'), 'cores': num(kv.get('cores'), int), 'cpu_limit': kv.get('cpu_limit'), 'mem_total_bytes': num(kv.get('mem_total_bytes'), int),
           'cgroup_version': kv.get('cgroup_version'), 'cgroup_memory_limit': kv.get('cgroup_memory_limit'), 'swap_total_bytes': num(kv.get('swap_total_bytes'), int), 'fs_type': kv.get('fs_type'),
           'disk_free_bytes': num(kv.get('disk_free_bytes'), int), 'disk_size_bytes': num(kv.get('disk_size_bytes'), int), 'dotnet': kv.get('dotnet'), 'rssanon_source': kv.get('rssanon_source', 'proc'),
           'preflight_result': kv.get('preflight_result'), 'in_container': kv.get('in_container'), 'crc32c_hardware': yes('env_crc32c_hardware'), 'crc32c_instruction': (h.get('env_crc32c_instruction') or (None,))[0],
           'crc32c_check_value': (h.get('env_crc32c_check_value') or (None,))[0], 'crc32c_table_equals_hardware': yes('env_crc32c_table_equals_hardware'), 'vector128': yes('env_vector128_hardware'),
           'vector256': yes('env_vector256_hardware'), 'os': (h.get('env_os') or (None,))[0]}
    return env


def md(form):
    L = ['# Hand-off form (package s1b-handoff, work package 1b.9f)', '', f"Generated {form['generated_utc']} by package commit `{form['package_commit']}`. The numbers below are what the scripts measured; send `form.json` back (by hand: nothing is sent anywhere).", '']
    e = form['environment']
    L += ['## Environment', '', '| item | value |', '|---|---|']
    for k in ('arch', 'image_arch', 'process_arch', 'emulation', 'kernel', 'cores', 'cpu_limit', 'mem_total_bytes', 'cgroup_version', 'cgroup_memory_limit', 'fs_type', 'disk_free_bytes', 'dotnet', 'crc32c_hardware', 'crc32c_instruction', 'crc32c_table_equals_hardware', 'rssanon_source', 'preflight_result'):
        L.append(f'| {k} | {e.get(k)} |')
    for label, run in form['runs'].items():
        b, r, d, v, s, q = run.get('build') or {}, run.get('rss') or {}, run.get('disk') or {}, run.get('verify') or {}, run.get('structure') or {}, run.get('queries') or {}
        n = run.get('n') or 0
        L += ['', f"## Run {label}: {run.get('kind')}, n = {n:,}, resume {'on' if run.get('resume') else 'off'}, {run.get('threads')} threads", '', '| item | value |', '|---|---|']
        L.append(f"| attempts (exit codes) | {', '.join(str(a['exit_code']) for a in b.get('attempts', []))} |")
        L.append(f"| seconds inside the build (sum) / wall (sum) | {b.get('total_s_sum')} / {b.get('wall_s_sum')} |")
        L.append('| phase seconds | ' + ', '.join(f'{k} {v_:.1f}' for k, v_ in (b.get('phases_s') or {}).items() if v_) + ' |')
        L.append(f"| container bytes | {b.get('container_bytes')} ({(b.get('container_bytes') or 0) / n:.2f} per symbol) |" if n else '| container bytes | |')
        L.append(f"| peak scratch bytes (build counter) / sampled on disk | {b.get('scratch_peak_bytes')} / {d.get('scratch_peak_bytes_sampled')} |")
        L.append(f"| peak of everything on the volume (sampled) | {d.get('total_peak_bytes_sampled')} |")
        L.append(f"| RssAnon peak MB: build's own sampler / independent sampler ({r.get('external_source')}) | {r.get('in_process_peak_mb')} / {r.get('external_peak_mb')} (gate 320) |")
        L.append(f"| container SHA-256 | `{(run.get('container') or {}).get('sha256')}` |")
        L.append(f"| e8-verify --container | {v.get('status')}, {v.get('seconds')} s |")
        L.append(f"| Structure level | {s.get('status')}, {s.get('seconds')} s |")
        L.append(f"| queries | {q.get('status')}: {q.get('total')} queries, search mismatches {q.get('sa_mismatches')}, engine {q.get('engine_status')} / mismatches {q.get('engine_mismatches')} |")
    L += ['', '## Tests', '', '| set | status | passed | failed | skipped | seconds |', '|---|---|---|---|---|---|']
    for name, t in form['tests'].items():
        L.append(f"| {name} | {t.get('status')} | {t.get('passed')} | {t.get('failed')} | {t.get('skipped')} | {t.get('seconds')} |")
    return '\n'.join(L) + '\n'


def main(argv):
    if not argv or argv[0].startswith('--'):
        print(__doc__)
        return 2
    data = argv[0]
    opt = dict(zip(argv[1::2], argv[2::2]))
    ref = load_reference(os.path.join(HERE, 'reference_hashes.tsv'))
    runs = {}
    for d in sorted(glob.glob(os.path.join(data, 'runs', '*'))):
        if os.path.isdir(d) and os.path.exists(os.path.join(d, 'run.json')):
            runs[os.path.basename(d)] = run_section(data, os.path.basename(d), ref)
    tests = {}
    for t in sorted(glob.glob(os.path.join(data, 'tests', '*.json'))):
        j = read_json(t)
        if isinstance(j, dict):
            tests[os.path.basename(t)[:-5]] = j
    form = {'form': 's1b-handoff', 'form_version': FORM_VERSION, 'generated_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            'package_commit': opt.get('--commit') or os.environ.get('S1B_COMMIT') or 'unknown', 'environment': environment_section(data), 'runs': runs, 'tests': tests, 'notes': []}
    out = opt.get('--out', os.path.join(data, 'form.json'))
    with open(out, 'w') as fh:
        json.dump(form, fh, indent=1, sort_keys=False)
    with open(opt.get('--md', os.path.join(data, 'form.md')), 'w') as fh:
        fh.write(md(form))
    print(f'form: {out} ({len(runs)} run(s), {len(tests)} test set(s))')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
