#!/usr/bin/env python3
# E7 cache-test sampler (part C, "the true 0.5x test"): independent of the program under test.
# usage: python3 -I sample_e7_cache.py <pidfile> <rows file> <case> <run> <scratch dir> <packed text file> <cgroup dir or -> <series file>
# Waits for <pidfile>, then until the process exits:
#   every 20 ms  RssAnon of /proc/<pid>/status (peak, and the mean over the window below), /proc/<pid>/io (last values)
#   every 200 ms the page residency of the packed text file and of the DCS rank file <scratch>/ranks.bin (mincore on a fresh read-only mapping: nothing is touched, nothing is charged), the cgroup's
#                memory.usage_in_bytes, memory.stat (cache, rss = anonymous, mapped_file, dirty, writeback), memory.failcnt, and the phase markers: ranks.bin exists (PH1 end .. PH3 end), <scratch>/ph3 exists (PH3 running)
# Rows: eff_resident_ph3 = time-weighted mean over PH3 (ranks.bin and ph3/ both exist) of resident bytes of (text + rank file) / their size; eff_resident_ranks = the same over the whole life of the rank file
# (PH1 end .. PH3 end); eff_budget_ph3 = (cgroup limit - mean cgroup anonymous memory over PH3) / (text + rank file) (the E3 style budget; limit max: not defined); the series (time, flags, residency, cgroup counters) goes to <series file>.
import ctypes, mmap, os, sys, time

libc = ctypes.CDLL(None, use_errno=True)
libc.mmap.restype = ctypes.c_void_p
libc.mmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_long]
libc.munmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
libc.mincore.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p]
PAGE = os.sysconf('SC_PAGE_SIZE')

def residency(path):
    """(resident bytes, size) of a file by mincore, or (0, 0) when it does not exist (yet)."""
    try:
        fd = os.open(path, os.O_RDONLY)
    except OSError:
        return 0, 0
    try:
        size = os.fstat(fd).st_size
        if size == 0:
            return 0, 0
        addr = libc.mmap(None, size, 1, 1, fd, 0)  # PROT_READ, MAP_SHARED
        if addr in (None, ctypes.c_void_p(-1).value):
            return 0, size
        try:
            pages = (size + PAGE - 1) // PAGE
            vec = (ctypes.c_ubyte * pages)()
            if libc.mincore(addr, size, vec) != 0:
                return 0, size
            res = sum(1 for b in vec if b & 1)
            return min(size, res * PAGE), size
        finally:
            libc.munmap(addr, size)
    finally:
        os.close(fd)

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

def readint(path):
    try:
        with open(path) as f:
            return int(f.read().strip())
    except (OSError, ValueError):
        return None

def cgstat(cg):
    out = {}
    try:
        with open(cg + '/memory.stat') as f:
            for line in f:
                k, v = line.split()
                out[k] = int(v)
    except OSError:
        pass
    return out

pidfile, rows, case, run, scratch, text, cg, series = sys.argv[1:9]
cg = None if cg == '-' else cg
while not os.path.exists(pidfile):
    time.sleep(0.005)
time.sleep(0.01)
pid = int(open(pidfile).read().strip())
ranks = os.path.join(scratch, 'ranks.bin')
ph3 = os.path.join(scratch, 'ph3')
t0 = time.monotonic()
next_slow = t0
period = 0.02
t = t0
max_rss = 0
last_io = {}
samples = []   # (time, ranks exist, ph3 exist, resident bytes, size bytes, cgroup anon, cgroup cache, usage)
rss_win = []   # RssAnon at 20 ms inside PH3
n = 0
f = open(series, 'w')
f.write('t_s\tranks\tph3\tresident_text\tsize_text\tresident_rank\tsize_rank\tcg_usage\tcg_rss\tcg_cache\tcg_mapped\tcg_dirty\tcg_writeback\tcg_failcnt\trssanon_kb\n')
while True:
    r = rssanon(pid)
    if r is None:
        break
    io = procio(pid)
    if io:
        last_io = io
    max_rss = max(max_rss, r)
    now = time.monotonic()
    if os.path.exists(ph3) and os.path.exists(ranks):
        rss_win.append(r)
    if now >= next_slow:
        ranks_ex = os.path.exists(ranks)
        ph3_ex = os.path.exists(ph3)
        rt, st = residency(text)
        rr, sr = residency(ranks) if ranks_ex else (0, 0)
        usage = readint(cg + '/memory.usage_in_bytes') if cg else None
        fail = readint(cg + '/memory.failcnt') if cg else None
        s = cgstat(cg) if cg else {}
        samples.append((now - t0, ranks_ex, ph3_ex, rt, st, rr, sr, s.get('rss'), s.get('cache'), usage))
        f.write('%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n' % (now - t0, ranks_ex, ph3_ex, rt, st, rr, sr, usage, s.get('rss'), s.get('cache'), s.get('mapped_file'), s.get('dirty'), s.get('writeback'), fail, r))
        next_slow = now + 0.2
    n += 1
    t += period
    delay = t - time.monotonic()
    if delay > 0:
        time.sleep(delay)
    else:
        t = time.monotonic()
f.close()

def window(pred):
    """time-weighted mean of the resident fraction of (text + rank file) over the samples where pred holds, plus the mean cgroup anonymous memory there"""
    num = den = 0.0
    anon = []
    cnt = 0
    for i in range(len(samples) - 1):
        s = samples[i]
        if not pred(s):
            continue
        dt = samples[i + 1][0] - s[0]
        size = s[4] + s[6]
        if size == 0:
            continue
        num += dt * (s[3] + s[5]) / size
        den += dt
        cnt += 1
        if s[7] is not None:
            anon.append(s[7])
    return (num / den if den > 0 else float('nan')), den, (sum(anon) / len(anon) if anon else None), cnt

e3, dur3, anon3, c3 = window(lambda s: s[1] and s[2])
er, durr, anonr, cr = window(lambda s: s[1])
size_rank = max((s[6] for s in samples), default=0)
size_text = max((s[4] for s in samples), default=0)
limit = readint(cg + '/memory.limit_in_bytes') if cg else None
with open(rows, 'a') as out:
    def row(metric, value, unit, notes=''):
        out.write('\t'.join(['R', 'E7', case, run, metric, str(value), unit, notes]) + '\n')
    row('ext_rssanon_peak', '%.3f' % (max_rss * 1024 / 1e6), 'MB', 'independent sampler: max RssAnon of /proc/<pid>/status at 20 ms')
    if rss_win:
        row('ext_rssanon_mean_ph3', '%.3f' % (sum(rss_win) / len(rss_win) * 1024 / 1e6), 'MB', 'mean RssAnon over PH3 (rank file and ph3 directory exist), %d samples' % len(rss_win))
    row('eff_resident_ph3', '%.4f' % e3, '', 'time-weighted resident fraction of (packed text + rank file) over PH3 (mincore every 200 ms, %d samples, %.1f s)' % (c3, dur3))
    row('eff_resident_ranks', '%.4f' % er, '', 'the same over the whole life of the rank file (PH1 end .. PH3 end, %d samples, %.1f s)' % (cr, durr))
    row('text_bytes', size_text, 'bytes', 'packed text file')
    row('rank_bytes', size_rank, 'bytes', 'DCS rank file')
    if limit is not None:
        row('cg_limit', limit, 'bytes', 'memory.limit_in_bytes of the sub-group')
        row('cg_failcnt', readint(cg + '/memory.failcnt'), 'n', 'memory.failcnt at the end')
        row('cg_max_usage', readint(cg + '/memory.max_usage_in_bytes'), 'bytes', 'memory.max_usage_in_bytes at the end')
        if anon3 is not None and size_text + size_rank > 0 and limit < 1 << 60:
            row('eff_budget_ph3', '%.4f' % ((limit - anon3) / (size_text + size_rank)), '', 'E3-style budget ratio (limit - mean cgroup anonymous memory over PH3 (%.1f MB)) / (text + rank file)' % (anon3 / 1e6))
        if anon3 is not None:
            row('cg_rss_mean_ph3', '%.3f' % (anon3 / 1e6), 'MB', 'mean memory.stat rss (anonymous) of the sub-group over PH3')
    for k in ('rchar', 'wchar', 'read_bytes', 'write_bytes'):
        row('ext_proc_' + k, last_io.get(k, 0), 'bytes', 'independent sampler: /proc/<pid>/io %s, last value before the process exited (whole process incl. runtime start)' % k)
