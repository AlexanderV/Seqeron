#!/usr/bin/env python3
"""fresh-process A/B of the whole production build: alternates the Bench of two builds (sesa build, 4 threads, psi, no fsync); usage: ab_fresh.py <text> <pairs> <label> [binA binB]"""
import subprocess, sys, os, statistics, re, shutil
text, pairs, label = sys.argv[1], int(sys.argv[2]), sys.argv[3]
bins = {'old': '/tmp/stbench/b9h/bin_old/SuffixTree.Esa.Bench/SuffixTree.Esa.Bench.dll', 'new': '/tmp/stbench/b9h/bin_new/SuffixTree.Esa.Bench/SuffixTree.Esa.Bench.dll'}
env = dict(os.environ, PATH=os.environ['PATH'] + ':/root/.dotnet', DOTNET_ROOT='/root/.dotnet', DOTNET_CLI_TELEMETRY_OPTOUT='1', DOTNET_NOLOGO='1', SEQERON_BENCH_DIR='/tmp/stbench')
env.pop('DOTNET_GCHeapHardLimit', None)
work = '/tmp/stbench/b9h/w_ab'
res = {'old': [], 'new': [], 'old_rss': [], 'new_rss': []}
def run(v):
    shutil.rmtree(work, ignore_errors=True); os.makedirs(work)
    out = subprocess.run(['dotnet', bins[v], 'sesa', 'build', '--input', text, '--work', work, '--threads', '4', '--psi', 'true', '--fsync', 'false'], env=env, capture_output=True, text=True).stdout
    t = float(re.search(r'total_s\t([\d.]+)', out).group(1)); r = float(re.search(r'rssanon_peak_mb\t([\d.]+)', out).group(1))
    return t, r
for i in range(pairs + 1):
    order = ['old', 'new'] if i % 2 == 0 else ['new', 'old']
    for v in order:
        t, r = run(v)
        if i > 0:
            res[v].append(t); res[v + '_rss'].append(r)
        print('%s pair %d %s total_s %.3f rssanon_mb %.1f' % (label, i, v, t, r), flush=True)
ratios = sorted(n / o for o, n in zip(res['old'], res['new']))
print('%s: old median %.3f s [%.3f-%.3f], new median %.3f s [%.3f-%.3f], new/old per pair median %.4f [%.4f-%.4f] over %d pairs; saving median %.3f s; rssanon peak old %.1f new %.1f MB (max over the pairs)' % (
    label, statistics.median(res['old']), min(res['old']), max(res['old']), statistics.median(res['new']), min(res['new']), max(res['new']), statistics.median(ratios), ratios[0], ratios[-1], len(ratios),
    statistics.median([o - n for o, n in zip(res['old'], res['new'])]), max(res['old_rss']), max(res['new_rss'])))
shutil.rmtree(work, ignore_errors=True)
