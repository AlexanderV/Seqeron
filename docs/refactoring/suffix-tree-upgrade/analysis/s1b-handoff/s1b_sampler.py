#!/usr/bin/env python3
"""Independent sampler of the hand-off package (work package 1b.9f): while a process runs it reads, every --interval-ms, RssAnon / RssFile / VmHWM / Threads from /proc/<pid>/status (the figure that the
320 MB gate of the production build is about: anonymous resident pages, the mapped file pages of the container are not in it), and every --disk-interval-ms the bytes that the directories of the run hold:
the scratch (a directory whose name ends in .tmp.d), the container (a file ending in .tmp or .sesa) and the rest (the input text, the logs), as allocated blocks, and the free space of the volume. It ends when the
process ends or the --stop-file appears, and writes a series (TSV) and a summary (JSON). It sends nothing anywhere.
Where /proc has no RssAnon line (macOS, or S1B_PROC points at nothing) it asks `ps -o rss=` instead and says so (`source` = ps-rss): that is the resident size WITH the mapped file pages, not RssAnon, and the
checker does not judge the gate by it.
usage: s1b_sampler.py --pid PID --out series.tsv --summary summary.json [--interval-ms 100] [--disk-interval-ms 1000] [--dir DIR]... [--stop-file FILE]"""
import json
import os
import subprocess
import sys
import time

PROC = os.environ.get('S1B_PROC', '/proc')


def read_status(pid):
    """(rssanon_kb, rssfile_kb, vmhwm_kb, threads) from /proc/<pid>/status, or None when the process or the RssAnon line is not there"""
    try:
        with open(f'{PROC}/{pid}/status') as fh:
            text = fh.read()
    except OSError:
        return None
    got = {}
    for line in text.splitlines():
        key, _, rest = line.partition(':')
        if key in ('RssAnon', 'RssFile', 'VmHWM', 'Threads'):
            try:
                got[key] = int(rest.split()[0])
            except (ValueError, IndexError):
                pass
    if 'RssAnon' not in got:
        return None
    return got['RssAnon'], got.get('RssFile', 0), got.get('VmHWM', 0), got.get('Threads', 0)


def read_ps(pid):
    try:
        out = subprocess.run(['ps', '-o', 'rss=', '-p', str(pid)], capture_output=True, text=True, timeout=5).stdout.split()
        return (int(out[0]), 0, 0, 0) if out else None
    except (OSError, ValueError, subprocess.SubprocessError):
        return None


def alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    # a zombie that nobody has waited for is gone for our purposes
    try:
        with open(f'/proc/{pid}/stat') as fh:
            return fh.read().rsplit(')', 1)[1].split()[0] != 'Z'
    except (OSError, IndexError):
        return True


def disk_usage(dirs):
    """allocated bytes of scratch, container and the rest under dirs (st_blocks * 512; a sparse file counts what it holds)"""
    scratch = container = other = 0
    for top in dirs:
        stack = [(top, False)]
        while stack:
            path, in_scratch = stack.pop()
            try:
                with os.scandir(path) as it:
                    for e in it:
                        try:
                            if e.is_dir(follow_symlinks=False):
                                stack.append((e.path, in_scratch or e.name.endswith('.tmp.d')))
                            elif e.is_file(follow_symlinks=False):
                                b = e.stat(follow_symlinks=False).st_blocks * 512
                                if in_scratch:
                                    scratch += b
                                elif e.name.endswith(('.tmp', '.sesa')):
                                    container += b
                                else:
                                    other += b
                        except OSError:
                            pass
            except OSError:
                pass
    return scratch, container, other


def free_bytes(path):
    try:
        st = os.statvfs(path)
        return st.f_bavail * st.f_frsize
    except OSError:
        return None


def main(argv):
    opts = {'dir': []}
    i = 0
    while i < len(argv):
        k = argv[i]
        if not k.startswith('--') or i + 1 >= len(argv):
            print(__doc__)
            return 2
        if k == '--dir':
            opts['dir'].append(argv[i + 1])
        else:
            opts[k[2:]] = argv[i + 1]
        i += 2
    if 'pid' not in opts or 'out' not in opts or 'summary' not in opts:
        print(__doc__)
        return 2
    pid = int(opts['pid'])
    interval = int(opts.get('interval-ms', 100)) / 1000.0
    disk_interval = int(opts.get('disk-interval-ms', 1000)) / 1000.0
    stop_file = opts.get('stop-file')
    dirs = opts['dir']
    t0 = time.monotonic()
    source = 'proc'
    first = read_status(pid)
    if first is None and alive(pid):
        first = read_ps(pid)
        source = 'ps-rss'
    s = {'source': source, 'interval_ms': int(interval * 1000), 'samples': 0, 'rssanon_peak_kb': 0, 'rssfile_peak_kb': 0, 'vmhwm_kb': 0, 'threads_max': 0,
         'scratch_peak_bytes': 0, 'container_peak_bytes': 0, 'total_peak_bytes': 0, 'free_min_bytes': None, 'free_start_bytes': free_bytes(dirs[0]) if dirs else None}
    next_disk = 0.0
    last_alive = time.monotonic()
    with open(opts['out'], 'w') as series:
        series.write('t_s\trssanon_kb\trssfile_kb\tvmhwm_kb\tthreads\tscratch_bytes\tcontainer_bytes\tother_bytes\tfree_bytes\n')
        cur = first
        while True:
            now = time.monotonic() - t0
            if cur is not None:
                s['samples'] += 1
                s['rssanon_peak_kb'] = max(s['rssanon_peak_kb'], cur[0])
                s['rssfile_peak_kb'] = max(s['rssfile_peak_kb'], cur[1])
                s['vmhwm_kb'] = max(s['vmhwm_kb'], cur[2])
                s['threads_max'] = max(s['threads_max'], cur[3])
                last_alive = time.monotonic()
            cols = ['', '', '', '']
            if dirs and now >= next_disk:
                next_disk = now + disk_interval
                sc, co, ot = disk_usage(dirs)
                fr = free_bytes(dirs[0])
                s['scratch_peak_bytes'] = max(s['scratch_peak_bytes'], sc)
                s['container_peak_bytes'] = max(s['container_peak_bytes'], co)
                s['total_peak_bytes'] = max(s['total_peak_bytes'], sc + co + ot)
                if fr is not None:
                    s['free_min_bytes'] = fr if s['free_min_bytes'] is None else min(s['free_min_bytes'], fr)
                cols = [str(sc), str(co), str(ot), '' if fr is None else str(fr)]
            if cur is not None:
                series.write(f'{now:.3f}\t{cur[0]}\t{cur[1]}\t{cur[2]}\t{cur[3]}\t' + '\t'.join(cols) + '\n')
            if stop_file and os.path.exists(stop_file):
                break
            if not alive(pid):
                break
            time.sleep(interval)
            cur = read_status(pid) if source == 'proc' else read_ps(pid)
            if cur is None and not alive(pid):
                break
    s['wall_s'] = round(time.monotonic() - t0, 3)
    s['rssanon_peak_mb'] = s['rssanon_peak_kb'] * 1024 / 1e6      # MB = 10^6 bytes, as the builder's own RssAnonMb
    s['rssfile_peak_mb'] = s['rssfile_peak_kb'] * 1024 / 1e6
    with open(opts['summary'], 'w') as fh:
        json.dump(s, fh, indent=1)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
