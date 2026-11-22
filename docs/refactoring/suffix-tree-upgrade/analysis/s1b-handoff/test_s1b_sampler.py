#!/usr/bin/env python3
"""Tests of the independent sampler s1b_sampler.py (work package 1b.9f): it reads RssAnon from /proc/<pid>/status of a process that holds a known amount of anonymous memory and sees it at least as
large as the amount; it separates the scratch directory, the container and the other files of the work directory by name; it ends when the process ends and when the stop file appears; it falls back
to `ps` (RSS, not RssAnon, and says so) where /proc has no RssAnon line. usage: python3 -I test_s1b_sampler.py"""
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLER = os.path.join(HERE, 's1b_sampler.py')

HOLD = r'''
import sys, time, os
n = int(sys.argv[1]) * 1024 * 1024
buf = bytearray(n)
for i in range(0, n, 4096):
    buf[i] = 1                      # touch every page: anonymous resident memory
work = sys.argv[2]
os.makedirs(os.path.join(work, 'out.sesa.tmp.d', 'ph1'), exist_ok=True)
with open(os.path.join(work, 'out.sesa.tmp.d', 'ph1', 'a.bin'), 'wb') as f:
    f.write(b'x' * (3 * 1024 * 1024))
with open(os.path.join(work, 'out.sesa.tmp'), 'wb') as f:
    f.write(b'y' * (5 * 1024 * 1024))
with open(os.path.join(work, 'text.txt'), 'wb') as f:
    f.write(b'z' * (1 * 1024 * 1024))
time.sleep(float(sys.argv[3]))
'''


def sample(hold_mb, seconds, extra=(), env=None):
    with tempfile.TemporaryDirectory() as work:
        child = subprocess.Popen([sys.executable, '-I', '-c', HOLD, str(hold_mb), work, str(seconds)])
        out = os.path.join(work, 'series.tsv')
        summ = os.path.join(work, 'summary.json')
        p = subprocess.run([sys.executable, '-I', SAMPLER, '--pid', str(child.pid), '--out', out, '--summary', summ, '--interval-ms', '20', '--disk-interval-ms', '100', '--dir', work, *extra],
                           capture_output=True, text=True, env=dict(os.environ, **(env or {})), timeout=60)
        child.wait()
        with open(summ) as fh, open(out) as sh:
            return p.returncode, json.load(fh), sh.read()


class SamplerTests(unittest.TestCase):
    def test_it_sees_the_anonymous_memory_of_the_process(self):
        code, s, series = sample(150, 1.5)
        self.assertEqual(code, 0)
        self.assertEqual(s['source'], 'proc')
        self.assertGreaterEqual(s['rssanon_peak_mb'], 150 * 1.048576 * 0.98, s)
        self.assertLess(s['rssanon_peak_mb'], 150 * 1.048576 + 60, s)
        self.assertGreater(s['samples'], 10)
        self.assertEqual(s['interval_ms'], 20)
        self.assertTrue(series.startswith('t_s\t'))

    def test_it_separates_scratch_container_and_the_rest_by_name(self):
        code, s, _ = sample(8, 1.5)
        self.assertEqual(code, 0)
        self.assertGreaterEqual(s['scratch_peak_bytes'], 3 * 1024 * 1024)
        self.assertLess(s['scratch_peak_bytes'], 4 * 1024 * 1024 + 65536)
        self.assertGreaterEqual(s['container_peak_bytes'], 5 * 1024 * 1024)
        self.assertLess(s['container_peak_bytes'], 6 * 1024 * 1024)
        self.assertGreaterEqual(s['total_peak_bytes'], 9 * 1024 * 1024)
        self.assertGreater(s['free_min_bytes'], 0)

    def test_it_ends_with_the_process_and_says_how_long_that_took(self):
        t = time.time()
        code, s, _ = sample(8, 0.8)
        self.assertEqual(code, 0)
        self.assertLess(time.time() - t, 15)
        self.assertGreater(s['wall_s'], 0.3)

    def test_the_stop_file_ends_it_while_the_process_still_runs(self):
        with tempfile.TemporaryDirectory() as work:
            child = subprocess.Popen([sys.executable, '-I', '-c', 'import time; time.sleep(30)'])
            stop = os.path.join(work, 'stop')
            sampler = subprocess.Popen([sys.executable, '-I', SAMPLER, '--pid', str(child.pid), '--out', os.path.join(work, 's.tsv'), '--summary', os.path.join(work, 's.json'), '--interval-ms', '20', '--stop-file', stop])
            time.sleep(0.5)
            open(stop, 'w').close()
            self.assertEqual(sampler.wait(timeout=10), 0)
            child.kill()
            child.wait()
            self.assertTrue(os.path.exists(os.path.join(work, 's.json')))

    def test_without_a_proc_rssanon_line_it_falls_back_to_ps_and_says_it_is_not_rssanon(self):
        code, s, _ = sample(40, 1.0, env={'S1B_PROC': '/nonexistent-proc'})
        self.assertEqual(code, 0)
        self.assertEqual(s['source'], 'ps-rss')
        self.assertGreaterEqual(s['rssanon_peak_mb'], 40)       # the RSS figure of ps (it includes mapped file pages), kept under the same name but flagged by `source`

    def test_a_pid_that_does_not_exist_gives_an_empty_summary_not_a_crash(self):
        with tempfile.TemporaryDirectory() as work:
            summ = os.path.join(work, 's.json')
            p = subprocess.run([sys.executable, '-I', SAMPLER, '--pid', '999999999', '--out', os.path.join(work, 's.tsv'), '--summary', summ, '--interval-ms', '20'], capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(summ) as fh:
                s = json.load(fh)
            self.assertEqual(s['samples'], 0)


if __name__ == '__main__':
    unittest.main(verbosity=1)
