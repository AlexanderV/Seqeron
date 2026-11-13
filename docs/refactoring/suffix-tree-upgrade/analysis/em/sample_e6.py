#!/usr/bin/env python3
# Independent /proc sampler of the E6 protocol: it does not share a line of code with the program under test (which samples RssAnon at 10 ms from /proc/self/status and Dirty at 20 ms from /proc/meminfo inside the process).
# usage: python3 -I sample_e6.py <pidfile> <rows file> <case> <run> <scratch dir>
# It takes the baseline of /proc/meminfo Dirty first (median of 5 samples, 20 ms apart), waits for <pidfile> to appear, then samples every 20 ms: RssAnon of /proc/<pid>/status, Dirty and Writeback of /proc/meminfo, the byte counters of
# /proc/<pid>/io (the last values seen before the process exits are its totals), and every 100 ms the size of the scratch directory (sum of st_size = apparent, sum of st_blocks * 512 = allocated; both the peak).
# It appends rows  R <TAB> E6 <TAB> case <TAB> run <TAB> metric <TAB> value <TAB> unit <TAB> notes.
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

def procio(pid):
    try:
        out = {}
        with open('/proc/%d/io' % pid, 'rb') as f:
            for line in f:
                k, v = line.split(b':')
                out[k.decode()] = int(v)
        return out
    except OSError:
        return None

def scratch_size(path):
    apparent = allocated = files = 0
    stack = [path]
    while stack:
        p = stack.pop()
        try:
            with os.scandir(p) as it:
                for e in it:
                    try:
                        if e.is_dir(follow_symlinks=False):
                            stack.append(e.path)
                        else:
                            st = e.stat(follow_symlinks=False)
                            apparent += st.st_size
                            allocated += st.st_blocks * 512
                            files += 1
                    except OSError:
                        pass
        except OSError:
            pass
    return apparent, allocated, files

pidfile, rows, case, run, scratch = sys.argv[1:6]
base = sorted(meminfo()[0] for _ in range(5) if not time.sleep(0.02))[2]
while not os.path.exists(pidfile):
    time.sleep(0.005)
time.sleep(0.01)
pid = int(open(pidfile).read().strip())
period = 0.02
t = time.monotonic()
next_scratch = t
max_rss = max_dirty = max_wb = 0
max_app = max_alloc = max_files = 0
n = 0
sum_delta = 0
last_io = {}
while True:
    r = rssanon(pid)
    if r is None:
        break
    d, w = meminfo()
    io = procio(pid)
    if io:
        last_io = io
    max_rss = max(max_rss, r)
    max_dirty = max(max_dirty, d)
    max_wb = max(max_wb, w)
    sum_delta += d - base
    n += 1
    now = time.monotonic()
    if now >= next_scratch:
        a, al, fl = scratch_size(scratch)
        max_app = max(max_app, a)
        max_alloc = max(max_alloc, al)
        max_files = max(max_files, fl)
        next_scratch = now + 0.1
    t += period
    delay = t - time.monotonic()
    if delay > 0:
        time.sleep(delay)
    else:
        t = time.monotonic()
with open(rows, 'a') as out:
    def row(metric, value, unit, notes=''):
        out.write('\t'.join(['R', 'E6', case, run, metric, str(value), unit, notes]) + '\n')
    row('ext_rssanon_peak', '%.3f' % (max_rss * 1024 / 1e6), 'MB', 'independent sampler: max RssAnon of /proc/<pid>/status over %d samples at 20 ms' % n)
    row('ext_dirty_peak', '%.3f' % ((max_dirty - base) * 1024 / 1e6), 'MB', 'independent sampler: max /proc/meminfo Dirty minus the baseline (%.2f MB, median of 5 samples before the start); absolute max %.2f MB' % (base * 1024 / 1e6, max_dirty * 1024 / 1e6))
    row('ext_dirty_mean', '%.3f' % (sum_delta / max(n, 1) * 1024 / 1e6), 'MB', 'independent sampler: mean of Dirty minus the baseline')
    row('ext_writeback_peak', '%.3f' % (max_wb * 1024 / 1e6), 'MB', 'independent sampler: max /proc/meminfo Writeback')
    row('ext_scratch_apparent_peak', max_app, 'bytes', 'independent sampler: peak sum of st_size of the files under the scratch directory (every 100 ms; the peak lasts seconds)')
    row('ext_scratch_allocated_peak', max_alloc, 'bytes', 'independent sampler: peak sum of st_blocks * 512 under the scratch directory (sparse files count only the written blocks)')
    row('ext_scratch_files_peak', max_files, 'n', 'independent sampler: most files under the scratch directory at one sample')
    for k in ('rchar', 'wchar', 'read_bytes', 'write_bytes'):
        row('ext_proc_' + k, last_io.get(k, 0), 'bytes', 'independent sampler: /proc/<pid>/io %s, last value before the process exited (whole process incl. runtime start)' % k)
