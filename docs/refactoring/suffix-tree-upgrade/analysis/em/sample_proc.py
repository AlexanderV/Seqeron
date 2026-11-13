#!/usr/bin/env python3
# Independent /proc sampler of the E5 protocol (the orchestrator's idea: "RssAnon sampling as in E0 plus an independent /proc sampler, 20 ms"): it does not share a line of code with the program under test
# (which samples RssAnon at 10 ms from /proc/self/status and Dirty at 20 ms from /proc/meminfo inside the process).
# usage: python3 -I sample_proc.py <pidfile> <rows file> <case> <run>
# It takes the baseline of /proc/meminfo Dirty first (median of 5 samples, 20 ms apart), waits for <pidfile> to appear, then samples RssAnon of /proc/<pid>/status and Dirty / Writeback of /proc/meminfo every 20 ms
# until the process is gone, and appends rows  R <TAB> E5 <TAB> case <TAB> run <TAB> metric <TAB> value <TAB> unit <TAB> notes.
import os, sys, time

def meminfo():
    d = w = 0
    with open('/proc/meminfo', 'rb') as f:
        for line in f:
            if line.startswith(b'Dirty:'):
                d = int(line.split()[1])
            elif line.startswith(b'Writeback:'):
                w = int(line.split()[1])
    return d, w

def rssanon(pid):
    try:
        with open('/proc/%d/status' % pid, 'rb') as f:
            for line in f:
                if line.startswith(b'RssAnon:'):
                    return int(line.split()[1])
    except OSError:
        return None
    return 0

pidfile, rows, case, run = sys.argv[1:5]
base = sorted(meminfo()[0] for _ in range(5) if not time.sleep(0.02))[2]
while not os.path.exists(pidfile):
    time.sleep(0.005)
time.sleep(0.01)
pid = int(open(pidfile).read().strip())
period = 0.02
t = time.monotonic()
max_rss = max_dirty = max_wb = 0
n = 0
sum_delta = 0
while True:
    r = rssanon(pid)
    if r is None:
        break
    d, w = meminfo()
    max_rss = max(max_rss, r)
    max_dirty = max(max_dirty, d)
    max_wb = max(max_wb, w)
    sum_delta += d - base
    n += 1
    t += period
    delay = t - time.monotonic()
    if delay > 0:
        time.sleep(delay)
    else:
        t = time.monotonic()
with open(rows, 'a') as out:
    def row(metric, value, unit, notes=''):
        out.write('\t'.join(['R', 'E5', case, run, metric, str(value), unit, notes]) + '\n')
    row('ext_rssanon_peak', '%.3f' % (max_rss * 1024 / 1e6), 'MB', 'independent sampler: max RssAnon of /proc/<pid>/status over %d samples at 20 ms' % n)
    row('ext_dirty_peak', '%.3f' % ((max_dirty - base) * 1024 / 1e6), 'MB', 'independent sampler: max /proc/meminfo Dirty minus the baseline (%.2f MB, median of 5 samples before the start); gate G3 <= 64 MB; absolute max %.2f MB' % (base * 1024 / 1e6, max_dirty * 1024 / 1e6))
    row('ext_dirty_mean', '%.3f' % (sum_delta / max(n, 1) * 1024 / 1e6), 'MB', 'independent sampler: mean of Dirty minus the baseline')
    row('ext_writeback_peak', '%.3f' % (max_wb * 1024 / 1e6), 'MB', 'independent sampler: max /proc/meminfo Writeback')
