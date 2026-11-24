#!/usr/bin/env python3
"""Tests of trx_summary.py (work package S2.2): the job summary of the deep workflow. `python3 -m unittest discover -s <this directory> -p "test_*.py"`."""
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from test_skip_accounting import make_trx  # noqa: E402
import trx_summary  # noqa: E402
import trx_lib  # noqa: E402

LARGE = "SuffixTree.Esa.Tests.Large.LargeRedTests"


class Summary(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)

    def trx(self, results):
        path = os.path.join(self.dir.name, "r.trx")
        with open(path, "w", encoding="utf-8") as f:
            f.write(make_trx(results))
        return path

    def expected_file(self, *names):
        path = os.path.join(self.dir.name, "red.txt")
        with open(path, "w", encoding="utf-8") as f:
            f.write("# c\n" + "".join(f"{n}\tS2.7\n" for n in names))
        return path

    def run_cli(self, *args):
        p = subprocess.run([sys.executable, os.path.join(HERE, "trx_summary.py"), *args], capture_output=True, text=True)
        return p.returncode, p.stdout

    def test_counts_and_the_failed_tests(self):
        path = self.trx([("A", "p1", "Passed", ""), ("A", "p2", "Passed", ""), ("A", "f1", "Failed", "boom\nsecond line"), ("A", "n1", "NotExecuted", "x")])
        code, out = self.run_cli("--title", "T", path)
        self.assertEqual(code, 0)
        self.assertIn("### T", out)
        self.assertIn("2 passed, 1 failed, 1 not executed", out)
        self.assertIn("`A.f1`: boom second line", out)
        self.assertIn("10 min 00 s", out)

    def test_an_error_a_timeout_and_an_aborted_test_are_failures_too(self):
        path = self.trx([("A", "p1", "Passed", ""), ("A", "e1", "Error", "crashed"), ("A", "t1", "Timeout", "too slow"), ("A", "a1", "Aborted", "host died")])
        code, out = self.run_cli(path)
        self.assertEqual(code, 0)
        self.assertIn("1 passed, 3 failed, 0 not executed", out)
        for name, message in (("e1", "crashed"), ("t1", "too slow"), ("a1", "host died")):
            self.assertIn(f"`A.{name}`: {message}", out)

    def test_the_red_tests_are_split_into_expected_unexpected_and_green(self):
        path = self.trx([(LARGE, "Red1", "Failed", "refused"), (LARGE, "Red2", "Passed", ""), (LARGE, "NewFailure", "Failed", "a real finding"), ("Other", "Fine", "Passed", "")])
        listed = self.expected_file(f"{LARGE}.Red1", f"{LARGE}.Red2", f"{LARGE}.Gone")
        code, out = self.run_cli("--expected-red", listed, path)
        self.assertEqual(code, 0)
        self.assertIn("Red as expected (1 of 3 listed)", out)
        self.assertIn(f"`{LARGE}.Red1` (green at S2.7)", out)
        self.assertIn("UNEXPECTED failures (1)", out)
        self.assertIn(f"`{LARGE}.NewFailure`: a real finding", out)
        self.assertIn(f"Listed as red but green now (1)", out)
        self.assertIn(f"`{LARGE}.Red2`", out.split("Listed as red but green now")[1])
        self.assertIn(f"Listed but not in the run (1)", out)

    def test_no_unexpected_failure_says_none(self):
        path = self.trx([(LARGE, "Red1", "Failed", "refused")])
        code, out = self.run_cli("--expected-red", self.expected_file(f"{LARGE}.Red1"), path)
        self.assertIn("UNEXPECTED failures (0):** none", out)

    def test_the_grep_lines_of_the_output(self):
        path = os.path.join(self.dir.name, "g.trx")
        text = make_trx([("A", "t1", "Passed", "noise\nLARGE-REPORT test=t1 seconds=1\nmore noise")])
        with open(path, "w", encoding="utf-8") as f:
            f.write(text)
        code, out = self.run_cli("--grep", "LARGE-REPORT", path)
        self.assertIn("LARGE-REPORT test=t1 seconds=1", out)
        self.assertNotIn("noise", out)

    def test_a_missing_trx_is_said_and_is_not_an_error_of_the_summary_step(self):
        code, out = self.run_cli("--title", "T", os.path.join(self.dir.name, "missing.trx"))
        self.assertEqual(code, 0)
        self.assertIn("no usable test result file", out)

    def test_the_committed_expected_list_names_the_six_red_tests(self):
        # nine at the end of S2.2; S2.3 turned Width4_TheViewsOfAnOpenedContainerExist_NoNullAtAnySize green and removed its row; S2.4 turned Width4_TheStatisticsOfTheHeaderAreNotNarrowedToInt and
        # Width4_TheDescentOverTheContainerAnswersContains_AtTheEndsOfTheText green and removed their rows (S2.7 removes four, S2.8 the last two)
        expected = trx_summary.read_expected(os.path.join(HERE, "large_expected_red.txt"))
        self.assertEqual(len(expected), 6)
        self.assertEqual(sorted(set(expected.values())), ["S2.7", "S2.8"])
        self.assertTrue(all(k.startswith(LARGE + ".") for k in expected))


if __name__ == "__main__":
    unittest.main()
