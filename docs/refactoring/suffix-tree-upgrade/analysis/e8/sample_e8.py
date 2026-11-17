#!/usr/bin/env python3
# E8 sampler: independent of the program under test (it shares no line of code with it). UNBUILT AND UNTESTED IN A CONTAINER, ON macOS AND ON arm64 (v1 cgroup path tested on the machine of the author, v2 only against a synthetic tree).
# usage: python3 -I sample_e8.py <pidfile> <rows file> <case> <run> <scratch dir> <out prefix> <input text> <data dir> <series file>
# Waits for <pidfile>, then until the process exits:
#   every 20 ms   RssAnon of /proc/<pid>/status, /proc/<pid>/io (the last values are the totals)
#   every 200 ms  page residency (mincore on a fresh read-only mapping: nothing is touched) of the packed text <out prefix>.pk and of the DCS rank file <scratch>/ranks.bin; the memory counters of the cgroup
#                 (cgroup v1 AND v2): usage, anonymous, file, mapped, dirty, writeback, the failure counter (v1 failcnt, v2 events max / oom_kill); the phase markers: ph1, ranks.bin, ph3, ph4, ph5 under the scratch directory
#   every 1 s     the size of the scratch directory, of the outputs and of the input (apparent: st_size, allocated: st_blocks x 512), the free space of the volume of the data directory, /proc/pressure/memory
# Rows: R<TAB>E8<TAB>case<TAB>run<TAB>metric<TAB>value<TAB>unit<TAB>notes   (appended to <rows file>); the 200 ms series goes to <series file>.
# E8_CGROUP_ROOT (default /sys/fs/cgroup) and E8_PROC (default /proc) move the roots the kernel's files are read from (tests).
import ctypes, os, sys, time

CG_ROOT = os.environ.get('E8_CGROUP_ROOT', '/sys/fs/cgroup')
PROC = os.environ.get('E8_PROC', '/proc')
PAGE = os.sysconf('SC_PAGE_SIZE')


def _libc():
    libc = ctypes.CDLL(None, use_errno=True)
    libc.mmap.restype = ctypes.c_void_p
    libc.mmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_long]
    libc.munmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    libc.mincore.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p]
    return libc


def residency(libc, path):
    """(resident bytes, size) of a file by mincore, (0, 0) when it does not exist (yet)."""
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
            return min(size, sum(1 for b in vec if b & 1) * PAGE), size
        finally:
            libc.munmap(addr, size)
    finally:
        os.close(fd)


def read_text(path):
    try:
        with open(path) as f:
            return f.read()
    except OSError:
        return None


def read_int(path):
    t = read_text(path)
    if t is None:
        return None
    t = t.strip()
    if t == 'max':
        return None
    try:
        return int(t)
    except ValueError:
        return None


def cg_find():
    """(version, directory) of the memory cgroup of this container: 2 / 1 / (0, None). A container with its own cgroup namespace sees its group at the mount root."""
    own = read_text(os.path.join(PROC, 'self/cgroup')) or ''
    if os.path.exists(os.path.join(CG_ROOT, 'cgroup.controllers')):
        sub = ''
        for line in own.splitlines():
            parts = line.split(':', 2)
            if len(parts) == 3 and parts[0] == '0':
                sub = parts[2]
        for d in (CG_ROOT + sub, CG_ROOT):
            if os.path.exists(os.path.join(d, 'memory.max')) or os.path.exists(os.path.join(d, 'memory.current')):
                return 2, d
        return 0, None
    base = os.path.join(CG_ROOT, 'memory')
    top = base if os.path.isdir(base) else CG_ROOT
    sub = ''
    for line in own.splitlines():
        parts = line.split(':', 2)
        if len(parts) == 3 and 'memory' in parts[1].split(','):
            sub = parts[2]
    for d in (top + sub, top):
        if os.path.exists(os.path.join(d, 'memory.limit_in_bytes')):
            return 1, d
    return 0, None


def cg_limit(version, d):
    if version == 2:
        return read_int(os.path.join(d, 'memory.max'))
    if version == 1:
        v = read_int(os.path.join(d, 'memory.limit_in_bytes'))
        return None if v is None or v >= 1 << 60 else v
    return None


def parse_kv(text):
    out = {}
    for line in (text or '').splitlines():
        p = line.split()
        if len(p) == 2:
            try:
                out[p[0]] = int(p[1])
            except ValueError:
                pass
    return out


def cg_snapshot(version, d):
    """unified counters: usage, anon, file, mapped, dirty, writeback, fail (v1 failcnt / v2 events max), oom_kill"""
    s = {}
    if version == 2:
        st = parse_kv(read_text(os.path.join(d, 'memory.stat')))
        ev = parse_kv(read_text(os.path.join(d, 'memory.events')))
        s = {'usage': read_int(os.path.join(d, 'memory.current')), 'anon': st.get('anon'), 'file': st.get('file'), 'mapped': st.get('file_mapped'),
             'dirty': st.get('file_dirty'), 'writeback': st.get('file_writeback'), 'fail': ev.get('max'), 'oom_kill': ev.get('oom_kill')}
    elif version == 1:
        st = parse_kv(read_text(os.path.join(d, 'memory.stat')))
        s = {'usage': read_int(os.path.join(d, 'memory.usage_in_bytes')), 'anon': st.get('rss'), 'file': st.get('cache'), 'mapped': st.get('mapped_file'),
             'dirty': st.get('dirty'), 'writeback': st.get('writeback'), 'fail': read_int(os.path.join(d, 'memory.failcnt')), 'oom_kill': None}
        s['oom_kill'] = parse_kv(read_text(os.path.join(d, 'memory.oom_control'))).get('oom_kill')
    return s


def cg_peak(version, d):
    if version == 2:
        return read_int(os.path.join(d, 'memory.peak'))
    if version == 1:
        return read_int(os.path.join(d, 'memory.max_usage_in_bytes'))
    return None


def rssanon_kb(pid):
    try:
        with open('%s/%d/status' % (PROC, pid), 'rb') as f:
            for line in f:
                if line.startswith(b'RssAnon:'):
                    return int(line.split()[1])
    except OSError:
        return None
    return 0


def procio(pid):
    try:
        out = {}
        with open('%s/%d/io' % (PROC, pid), 'rb') as f:
            for line in f:
                k, v = line.split(b':')
                out[k.decode()] = int(v)
        return out
    except OSError:
        return None


def tree_size(path):
    apparent = allocated = 0
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
                    except OSError:
                        pass
        except OSError:
            pass
    return apparent, allocated


def files_size(prefix):
    apparent = allocated = 0
    d, base = os.path.split(prefix)
    try:
        for name in os.listdir(d or '.'):
            if name.startswith(base + '.'):
                try:
                    st = os.stat(os.path.join(d, name))
                    apparent += st.st_size
                    allocated += st.st_blocks * 512
                except OSError:
                    pass
    except OSError:
        pass
    return apparent, allocated


def volume_free(path):
    try:
        v = os.statvfs(path)
        return v.f_bavail * v.f_frsize, (v.f_blocks - v.f_bfree) * v.f_frsize
    except OSError:
        return None, None


def psi_total(kind):
    t = read_text(os.path.join(PROC, 'pressure/memory'))
    if not t:
        return None
    for line in t.splitlines():
        if line.startswith(kind + ' '):
            for part in line.split():
                if part.startswith('total='):
                    return int(part[6:])
    return None


def main():
    pidfile, rows, case, run, scratch, prefix, inp, data, series = sys.argv[1:10]
    libc = _libc()
    version, cgdir = cg_find()
    limit0 = cg_limit(version, cgdir) if cgdir else None
    text_path = prefix + '.pk'
    ranks = os.path.join(scratch, 'ranks.bin')
    markers = {'ph1': os.path.join(scratch, 'ph1'), 'ranks': ranks, 'ph3': os.path.join(scratch, 'ph3'), 'ph4': os.path.join(scratch, 'ph4'), 'ph5': os.path.join(scratch, 'ph5')}
    while not os.path.exists(pidfile):
        time.sleep(0.005)
    time.sleep(0.01)
    pid = int(open(pidfile).read().strip())
    in_size = os.path.getsize(inp) if os.path.exists(inp) else 0
    free0, used0 = volume_free(data)
    psi_some0, psi_full0 = psi_total('some'), psi_total('full')
    t0 = time.monotonic()
    next_slow = next_disk = t0
    period = 0.02
    t = t0
    max_rss = 0
    rss_sum = rss_n = 0
    last_io = {}
    samples = []        # (time, ranks exist, ph3 exist, resident text, size text, resident rank, size rank, cg anon)
    first_seen, last_seen = {}, {}
    max_cg = {'usage': 0, 'anon': 0, 'file': 0, 'mapped': 0, 'dirty': 0, 'writeback': 0}
    cg_last = {}
    app_peak = alloc_peak = total_app_peak = total_alloc_peak = 0
    free_min = free0
    used_peak = 0
    f = open(series, 'w')
    f.write('t_s\tranks\tph3\tresident_text\tsize_text\tresident_rank\tsize_rank\tcg_usage\tcg_anon\tcg_file\tcg_mapped\tcg_dirty\tcg_writeback\tcg_fail\trssanon_kb\tscratch_apparent\n')
    scratch_app = 0
    rss_ph3 = []
    while True:
        r = rssanon_kb(pid)
        if r is None:
            break
        io = procio(pid)
        if io:
            last_io = io
        max_rss = max(max_rss, r)
        rss_sum += r
        rss_n += 1
        now = time.monotonic()
        if now >= next_slow:
            flags = {k: os.path.exists(v) for k, v in markers.items()}
            for k, v in flags.items():
                if v:
                    first_seen.setdefault(k, now - t0)
                    last_seen[k] = now - t0
            rt, st = residency(libc, text_path)
            rr, sr = residency(libc, ranks) if flags['ranks'] else (0, 0)
            snap = cg_snapshot(version, cgdir) if cgdir else {}
            for k in max_cg:
                if snap.get(k) is not None:
                    max_cg[k] = max(max_cg[k], snap[k])
            cg_last = snap or cg_last
            if flags['ranks'] and flags['ph3']:
                rss_ph3.append(r)
            samples.append((now - t0, flags['ranks'], flags['ph3'], rt, st, rr, sr, snap.get('anon')))
            f.write('%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n' % (now - t0, flags['ranks'], flags['ph3'], rt, st, rr, sr, snap.get('usage'), snap.get('anon'), snap.get('file'),
                                                                                         snap.get('mapped'), snap.get('dirty'), snap.get('writeback'), snap.get('fail'), r, scratch_app))
            next_slow = now + 0.2
        if now >= next_disk:
            a, al = tree_size(scratch)
            oa, oal = files_size(prefix)
            scratch_app = a
            app_peak, alloc_peak = max(app_peak, a), max(alloc_peak, al)
            total_app_peak = max(total_app_peak, a + oa + in_size)
            total_alloc_peak = max(total_alloc_peak, al + oal + in_size)
            fr, us = volume_free(data)
            if fr is not None:
                free_min = fr if free_min is None else min(free_min, fr)
                if us is not None and used0 is not None:
                    used_peak = max(used_peak, us - used0)
            next_disk = now + 1.0
        t += period
        delay = t - time.monotonic()
        if delay > 0:
            time.sleep(delay)
        else:
            t = time.monotonic()
    f.close()

    def window(pred):
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
    final = cg_snapshot(version, cgdir) if cgdir else {}
    with open(rows, 'a') as out:
        def row(metric, value, unit, notes=''):
            out.write('\t'.join(['R', 'E8', case, run, metric, str(value), unit, notes]) + '\n')
        row('ext_rssanon_peak', '%.3f' % (max_rss * 1024 / 1e6), 'MB', 'independent sampler: max RssAnon of /proc/<pid>/status over %d samples at 20 ms (the dotnet process of e6)' % rss_n)
        row('ext_rssanon_mean', '%.3f' % (rss_sum / max(rss_n, 1) * 1024 / 1e6), 'MB', 'mean RssAnon over the run')
        if rss_ph3:
            row('ext_rssanon_mean_ph3', '%.3f' % (sum(rss_ph3) / len(rss_ph3) * 1024 / 1e6), 'MB', 'mean RssAnon while the rank file and the ph3 directory exist, %d samples' % len(rss_ph3))
        for k in ('rchar', 'wchar', 'read_bytes', 'write_bytes'):
            row('ext_proc_' + k, last_io.get(k, 0), 'bytes', 'independent sampler: /proc/<pid>/io %s, last value before the process exited (whole process incl. runtime start)' % k)
        row('eff_resident_ph3', '%.4f' % e3, '', 'time-weighted resident fraction of (packed text + rank file) over PH3 (rank file and ph3 directory exist): mincore every 200 ms, %d samples, %.1f s' % (c3, dur3))
        row('eff_resident_ranks', '%.4f' % er, '', 'the same over the whole life of the rank file (end of PH1 .. end of PH3), %d samples, %.1f s' % (cr, durr))
        row('text_bytes', size_text, 'bytes', 'packed text file <prefix>.pk')
        row('rank_bytes', size_rank, 'bytes', 'DCS rank file (largest size seen)')
        for k in markers:
            if k in first_seen:
                row('marker_%s_first_s' % k, '%.2f' % first_seen[k], 's', 'first sample at which %s exists under the scratch directory (200 ms resolution), seconds from the start of the sampler' % k)
                row('marker_%s_last_s' % k, '%.2f' % last_seen[k], 's', 'last sample at which it exists')
        row('ext_scratch_apparent_peak', app_peak, 'bytes', 'peak sum of st_size under the scratch directory (every 1 s; the peak of PH4 lasts a minute or more at 1 Gbp)')
        row('ext_scratch_allocated_peak', alloc_peak, 'bytes', 'peak sum of st_blocks x 512 under the scratch directory')
        row('ext_disk_total_apparent_peak', total_app_peak, 'bytes', 'peak of scratch + outputs (<prefix>.*) + input text, apparent sizes (gate G3: at most 70 GB at 3 Gbp)')
        row('ext_disk_total_allocated_peak', total_alloc_peak, 'bytes', 'the same with allocated blocks')
        if free_min is not None:
            row('ext_disk_free_min', free_min, 'bytes', 'smallest free space of the volume of the data directory (statvfs every 1 s); at the start %s' % free0)
            row('ext_disk_used_peak_delta', used_peak, 'bytes', 'largest increase of the used space of the volume over its value at the start')
        row('cg_version', version, '', 'cgroup hierarchy read (0 = none found)')
        if limit0 is not None:
            row('cg_limit', limit0, 'bytes', 'memory limit of the container at the start of the run')
        else:
            row('cg_limit', 'max', '', 'no memory limit')
        if cgdir:
            peak = cg_peak(version, cgdir)
            row('cg_usage_max_sampled', max_cg['usage'], 'bytes', 'largest memory.current / usage_in_bytes seen at 200 ms')
            if peak is not None:
                row('cg_peak', peak, 'bytes', 'kernel peak of the group (v1 max_usage_in_bytes, v2 memory.peak)')
            row('cg_anon_max', max_cg['anon'], 'bytes', 'largest anonymous memory of the group (v1 rss, v2 anon)')
            row('cg_file_max', max_cg['file'], 'bytes', 'largest page cache of the group (v1 cache, v2 file)')
            row('cg_mapped_max', max_cg['mapped'], 'bytes', 'largest mapped file memory of the group')
            row('cg_dirty_max', max_cg['dirty'], 'bytes', 'largest dirty memory of the group')
            row('cg_writeback_max', max_cg['writeback'], 'bytes', 'largest writeback memory of the group')
            row('cg_failcnt', final.get('fail'), 'n', 'times the group hit its limit (v1 memory.failcnt, v2 memory.events max) at the end of the run')
            if final.get('oom_kill') is not None:
                row('cg_oom_kill', final.get('oom_kill'), 'n', 'OOM kills of the group at the end (v2 memory.events, v1 oom_control when reported)')
            if anon3 is not None:
                row('cg_anon_mean_ph3', '%.3f' % (anon3 / 1e6), 'MB', 'mean anonymous memory of the group over PH3')
                if limit0 is not None and size_text + size_rank > 0:
                    row('eff_budget_ph3', '%.4f' % ((limit0 - anon3) / (size_text + size_rank)), '', 'E3-style budget ratio: (limit - mean anonymous memory over PH3) / (text + rank file)')
        s1, f1 = psi_total('some'), psi_total('full')
        if s1 is not None and psi_some0 is not None:
            row('psi_memory_some_us', s1 - psi_some0, 'us', 'microseconds some task stalled on memory during the run (/proc/pressure/memory, VM-wide)')
            if f1 is not None and psi_full0 is not None:
                row('psi_memory_full_us', f1 - psi_full0, 'us', 'microseconds all tasks stalled on memory during the run (VM-wide)')


if __name__ == '__main__':
    main()
