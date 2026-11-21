#!/usr/bin/env python3
"""Tables of rows_1b9d_cache.tsv (run_1b9d_cache.sh): one line per cell (input, limit): the runs that finished, the median [min-max] of the time of the whole build and of its phases, the slowdown against the unlimited
cell of the same input and configuration, the effective cache (time-weighted resident fraction of the TEXT section and the rank file over PH3), the bytes read from the device (the process's read_bytes and the device's
sectors), the dirty-page amplification (pages dirtied by the process = proc_write_bytes over the bytes that the build wrote: explicit + mapped + the zeros of the reservation) and the device's, RssAnon, memory.failcnt,
the abandon / thrash flags of E7 and the identity of the container. usage: summarize_1b9d_cache.py rows_1b9d_cache.tsv"""
import collections, statistics, sys

rows = collections.defaultdict(lambda: collections.defaultdict(dict))  # case -> run -> metric -> value
notes = {}
for line in open(sys.argv[1]):
    p = line.rstrip('\n').split('\t')
    if len(p) < 6 or p[0] != 'R' or p[1] != '1b9d':
        continue
    rows[p[2]][p[3]][p[4]] = p[5]
    notes.setdefault((p[2], p[4]), p[6] if len(p) > 6 else '')

def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None

def med(vals, fmt='%.3g'):
    v = [x for x in vals if x is not None]
    if not v:
        return '-'
    if len(v) == 1:
        return fmt % v[0]
    return (fmt + ' [' + fmt + '-' + fmt + ']') % (statistics.median(v), min(v), max(v))

def metric(case, name, scale=1.0):
    return [f(r.get(name)) * scale if f(r.get(name)) is not None else None for r in rows[case].values() if r.get('run_done') == 'yes']

def logical(r):
    try:
        return float(r['counted_write']) + float(r['counted_mapped_write']) + float(r['counted_reserved_write'])
    except (KeyError, ValueError):
        return None

cases = sorted(rows)
unl = {}
for c in cases:
    inp = c.split('.')[0]
    if rows[c] and next(iter(rows[c].values())).get('settings', '').startswith('limit=max'):
        t = [x for x in metric(c, 'total_s') if x is not None]
        if t:
            unl[inp] = statistics.median(t)

print('| cell | runs finished / started | total s | slowdown | PH1 s | PH2+PH3 s | PH4+PH5 s | PH6 s | CRC + seal s | effective cache (PH3) | device reads MB | dirty amplification (pages / device) | RssAnon peak MB | failcnt | killed | container identical |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for c in cases:
    runs = rows[c]
    done = [r for r in runs.values() if r.get('run_done') == 'yes' and r.get('exit_code') == '0']
    inp = c.split('.')[0]
    total = [f(r.get('total_s')) for r in done]
    slow = [t / unl[inp] if t is not None and inp in unl else None for t in total]
    ph6 = [(f(r.get('bwt_s')) or 0) + (f(r.get('psi_s')) or 0) for r in done]
    tail = [(f(r.get('crc_s')) or 0) + (f(r.get('seal_s')) or 0) for r in done]
    ph23 = [(f(r.get('ph23_s')) or 0) for r in done]
    ph45 = [(f(r.get('ph45_s')) or 0) for r in done]
    dev_read = [(f(r.get('dev_sectors_read')) or 0) * 512 / 1e6 for r in done]
    amp = []
    amp_dev = []
    for r in done:
        lg = logical(r)
        if lg:
            amp.append(float(r['proc_write_bytes']) / lg)
            amp_dev.append(float(r['dev_sectors_written_synced']) * 512 / lg)
    killed = ','.join(sorted({r.get('killed', '?') for r in runs.values() if r.get('killed', 'no') != 'no'})) or 'no'
    ident = ' '.join(sorted({r.get('cmp_container', '-') for r in done})) if done else '-'
    print('| %s | %d / %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s / %s | %s | %s | %s | %s |' % (
        c, len(done), len(runs), med(total, '%.1f'), med(slow, '%.2f'), med([f(r.get('ph1_s')) for r in done], '%.2f'), med(ph23, '%.1f'), med(ph45, '%.1f'), med(ph6, '%.1f'), med(tail, '%.2f'),
        med([f(r.get('eff_resident_ph3')) for r in done], '%.3f'), med(dev_read, '%.0f'), med(amp, '%.2f'), med(amp_dev, '%.2f'),
        med([f(r.get('ext_rssanon_peak')) for r in done], '%.0f'), med([f(r.get('cg_failcnt_final')) for r in runs.values()], '%.3g'), killed, ident))
