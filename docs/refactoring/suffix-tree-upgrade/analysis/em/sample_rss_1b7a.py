#!/usr/bin/env python3
"""Independent RssAnon sampler for gate G3 of work package 1b.7a: polls /proc every 20 ms for a process whose command line names testhost.dll (the test host of dotnet test), keeps the maximum of its RssAnon
   (anonymous resident pages, kB) and appends one row '#ext <label> <variant> <max RssAnon kB> <samples>' to the rows file when the process has gone (or after 180 s). No code of the measured program.
   usage: sample_rss_1b7a.py <rows file> <label> <variant>"""
import os, sys, time

rows, label, variant = sys.argv[1:4]
deadline = time.time() + 180
seen = False
best = 0
samples = 0
while time.time() < deadline:
    found = False
    for pid in os.listdir('/proc'):
        if not pid.isdigit():
            continue
        try:
            with open('/proc/%s/cmdline' % pid, 'rb') as f:
                cmd = f.read()
            if b'testhost.dll' not in cmd:
                continue
            with open('/proc/%s/status' % pid) as f:
                for line in f:
                    if line.startswith('RssAnon:'):
                        best = max(best, int(line.split()[1]))
                        samples += 1
                        found = True
                        break
        except (OSError, ValueError):
            continue
    if found:
        seen = True
    elif seen:
        break
    time.sleep(0.02)
with open(rows, 'a') as out:
    out.write('#ext\t%s\t%s\t%d\t%d\n' % (label, variant, best, samples))
