#!/usr/bin/env python3
"""Tests of skip_accounting.py and trx_lib.py (work package S2.2): `python3 -m unittest discover -s docs/refactoring/suffix-tree-upgrade/analysis/ci -p "test_*.py"`.

The trx files are built in memory in the format of `dotnet test --logger trx` (the NUnit adapter): a result per test with its outcome and, for a test that did not
run, the message (an `Assert.Ignore` reason, or the reason of an [Explicit] test). Python standard library only; no network, no dotnet.
"""
import os
import subprocess
import sys
import tempfile
import unittest
import uuid
from xml.sax.saxutils import escape, quoteattr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import skip_accounting as sa  # noqa: E402
import trx_lib  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, *[".."] * 5))
REAL_LIST = os.path.join(REPO, "tests", "SuffixTree", "SuffixTree.Esa.Tests", "skips.tsv")
NS = "http://microsoft.com/schemas/VisualStudio/TeamTest/2010"
HEADER = "kind\ttest\tin_builds\treason\thow_to_unskip"


def make_trx(results, start="2026-11-24T05:00:00.0000000+00:00", finish="2026-11-24T05:10:00.0000000+00:00"):
    """results: (class_name, test_name, outcome, message) -> the text of a trx."""
    defs, res = [], []
    ids = {}
    for class_name, test_name, outcome, message in results:
        unit = ids.setdefault((class_name, test_name), str(uuid.uuid4()))
        defs.append(f'<UnitTest name={quoteattr(test_name)} id="{unit}"><TestMethod className={quoteattr(class_name)} name={quoteattr(test_name)} /></UnitTest>')
        out = ""
        if message:
            out = f"<Output><StdOut>{escape(message)}</StdOut><ErrorInfo><Message>{escape(message)}</Message></ErrorInfo></Output>"
        res.append(f'<UnitTestResult testId="{unit}" testName={quoteattr(test_name)} outcome="{outcome}" duration="00:00:00.0010000">{out}</UnitTestResult>')
    return (f'<?xml version="1.0" encoding="utf-8"?><TestRun xmlns="{NS}"><Times start="{start}" finish="{finish}" />'
            f'<Results>{"".join(res)}</Results><TestDefinitions>{"".join(defs)}</TestDefinitions></TestRun>')


C = "SuffixTree.Esa.Tests."


def run_rows():
    """A small list and the run that matches it in the default build."""
    rows = [
        ("ignore", C + "External.BoundsTests.OnlyChecked1", "default", "bounds: only the checked build", "the checked build"),
        ("ignore", C + "External.BoundsTests.OnlyChecked2(7)", "default", "bounds: only the checked build", "the checked build"),
        ("ignore", C + "Contract.ConformanceTests(ukkonen-adapter).NoNavigation", "both", "backend: no navigation", "never"),
        ("ignore", C + "External.FullDiskTests.NeedsMount", "host", "host: could not mount a tmpfs", "root"),
        ("explicit", C + "Review.ReviewTests.ReportOnly", "both", "report only", "--filter Category=ReviewIntrinsics"),
        ("explicit-class", C + "Contract.ConformanceTests(esa-file-w5)", "both", "POS_WIDTH 5 container", "--filter Category=Width5"),
    ]
    text = "\n".join(["# a list", HEADER] + ["\t".join(r) for r in rows]) + "\n"
    results = [
        (C + "External.BoundsTests", "OnlyChecked1", "NotExecuted", "only the checked build"),
        (C + "External.BoundsTests", "OnlyChecked2(7)", "NotExecuted", "only the checked build"),
        (C + "External.BoundsTests", "AnOrdinaryTest", "Passed", ""),
        (C + "Contract.ConformanceTests(ukkonen-adapter)", "NoNavigation", "NotExecuted", "no navigation"),
        (C + "Contract.ConformanceTests(esa-file)", "NoNavigation", "Passed", ""),
        (C + "External.FullDiskTests", "NeedsMount", "Passed", ""),
        (C + "Review.ReviewTests", "ReportOnly", "NotExecuted", "report only"),
        (C + "Review.ReviewTests", "AnotherTest", "Passed", ""),
        (C + "Contract.ConformanceTests(esa-file-w5)", "A", "NotExecuted", "OneTimeSetUp: POS_WIDTH 5 container"),
        (C + "Contract.ConformanceTests(esa-file-w5)", "B", "NotExecuted", "OneTimeSetUp: POS_WIDTH 5 container"),
    ]
    return text, results


class Base(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)

    def write(self, name, text):
        path = os.path.join(self.dir.name, name)
        with open(path, "w", encoding="utf-8") as f:
            f.write(text)
        return path

    def account(self, list_text, results, mode="default", partial=False):
        rows = sa.load_list(self.write("skips.tsv", list_text))
        return sa.account(trx_lib.read_trx(self.write("run.trx", make_trx(results))), rows, mode, partial)

    def failures_of(self, report, tag):
        return [f for f in report.failures if f.startswith(tag)]


class TheRunThatMatchesTheList(Base):
    def test_a_default_run_is_accounted_without_a_failure(self):
        text, results = run_rows()
        report = self.account(text, results)
        self.assertEqual(report.failures, [])
        self.assertEqual(report.counts["skipped"], 3)  # two bounds tests and the backend test: the host test ran
        self.assertEqual(report.counts["explicit"], 3)  # the explicit method and the two tests of the explicit fixture
        self.assertEqual(report.counts["executed"], 4)

    def test_a_checked_run_expects_the_bounds_tests_to_run_and_the_backend_test_to_skip(self):
        text, results = run_rows()
        results = [(c, n, "Passed", "") if c.endswith("BoundsTests") and n.startswith("OnlyChecked") else (c, n, o, m) for c, n, o, m in results]
        report = self.account(text, results, mode="checked")
        self.assertEqual(report.failures, [])
        self.assertEqual(report.counts["skipped"], 1)

    def test_an_explicit_test_is_not_a_skip(self):
        text, results = run_rows()
        report = self.account(text, results)
        self.assertEqual(report.counts["skipped"] + report.counts["explicit"], 6)
        self.assertNotIn("ReportOnly", " ".join(report.failures))

    def test_the_host_row_may_be_skipped_or_run_in_either_mode(self):
        text, results = run_rows()
        skipped = [(c, n, "NotExecuted", "could not mount") if n == "NeedsMount" else (c, n, o, m) for c, n, o, m in results]
        self.assertEqual(self.account(text, skipped).failures, [])
        self.assertEqual(self.account(text, results).failures, [])
        checked = [(c, n, "Passed", "") if n.startswith("OnlyChecked") else (c, n, o, m) for c, n, o, m in skipped]
        self.assertEqual(self.account(text, checked, mode="checked").failures, [])

    def test_a_failed_test_is_not_a_skip_and_the_failure_is_left_to_the_test_step(self):
        text, results = run_rows()
        results = [(c, n, "Failed", "boom") if n == "AnOrdinaryTest" else (c, n, o, m) for c, n, o, m in results]
        self.assertEqual(self.account(text, results).failures, [])


class AnUnlistedSkip(Base):
    def test_an_assert_ignore_x_in_any_test_is_a_failure(self):
        text, results = run_rows()
        results = [(c, n, "NotExecuted", "x") if n == "AnOrdinaryTest" else (c, n, o, m) for c, n, o, m in results]
        report = self.account(text, results)
        failures = self.failures_of(report, "UNLISTED SKIP")
        self.assertEqual(len(failures), 1, report.failures)
        self.assertIn("External.BoundsTests.AnOrdinaryTest", failures[0])
        self.assertIn("x", failures[0])

    def test_an_unlisted_skip_in_the_checked_mode_is_a_failure_too(self):
        text, results = run_rows()
        results = [(c, n, "Passed", "") if n.startswith("OnlyChecked") else (c, n, o, m) for c, n, o, m in results]
        results = [(c, n, "NotExecuted", "x") if n == "AnotherTest" else (c, n, o, m) for c, n, o, m in results]
        self.assertEqual(len(self.failures_of(self.account(text, results, mode="checked"), "UNLISTED SKIP")), 1)

    def test_an_inconclusive_test_is_a_skip_that_must_be_listed(self):
        text, results = run_rows()
        results = [(c, n, "Inconclusive", "the counter was noisy") if n == "AnOrdinaryTest" else (c, n, o, m) for c, n, o, m in results]
        self.assertEqual(len(self.failures_of(self.account(text, results), "UNLISTED SKIP")), 1)

    def test_a_new_explicit_method_is_an_unlisted_skip_until_it_is_listed(self):
        text, results = run_rows()
        results.append((C + "Review.ReviewTests", "AFreshExplicit", "NotExecuted", "Heavy; run with --filter Category=Heavy"))
        self.assertEqual(len(self.failures_of(self.account(text, results), "UNLISTED SKIP")), 1)

    def test_a_skip_in_a_class_that_only_looks_like_a_listed_one_is_unlisted(self):
        text, results = run_rows()
        # the listed fixture is `esa-file-w5`; `esa-file-w5-trusted` shares its prefix and is another fixture
        results.append((C + "Contract.ConformanceTests(esa-file-w5-trusted)", "A", "NotExecuted", "x"))
        results.append((C + "Contract.ConformanceTests(esa-file)", "NoNavigation2", "NotExecuted", "x"))
        self.assertEqual(len(self.failures_of(self.account(text, results), "UNLISTED SKIP")), 2)

    def test_a_name_that_only_starts_like_a_listed_one_is_unlisted(self):
        text, results = run_rows()
        # a listed ignore row `OnlyChecked1` / `OnlyChecked2(7)`; a listed explicit fixture `...FuzzTests` (written without arguments)
        text += "\t".join(["explicit-class", C + "Fuzz.SomeFuzzTests", "both", "r", "n"]) + "\n"
        results.append((C + "Fuzz.SomeFuzzTests", "T", "NotExecuted", "OneTimeSetUp: heavy"))
        results.append((C + "External.BoundsTests", "OnlyChecked10", "NotExecuted", "x"))
        results.append((C + "External.BoundsTests", "OnlyChecked2(70)", "NotExecuted", "x"))
        results.append((C + "Fuzz.SomeFuzzTestsExtra", "T", "NotExecuted", "x"))
        failures = self.failures_of(self.account(text, results), "UNLISTED SKIP")
        self.assertEqual(len(failures), 3, failures)
        self.assertTrue(any("OnlyChecked10" in f for f in failures))
        self.assertTrue(any("OnlyChecked2(70)" in f for f in failures))
        self.assertTrue(any("SomeFuzzTestsExtra" in f for f in failures))

    def test_the_same_name_in_another_class_is_unlisted(self):
        text, results = run_rows()
        results.append((C + "Other.BoundsTests", "OnlyChecked1", "NotExecuted", "x"))
        self.assertEqual(len(self.failures_of(self.account(text, results), "UNLISTED SKIP")), 1)

    def test_a_fixture_that_is_skipped_in_its_setup_and_is_not_listed_fails_and_is_proposed_once_as_a_class_row(self):
        # an [Explicit] fixture (or a OneTimeSetUp ignore) that nobody listed: every one of its results is an unlisted skip; the proposal is ONE explicit-class row
        text, results = run_rows()
        results.append((C + "Fuzz.NewHeavyTests", "T1", "NotExecuted", "OneTimeSetUp: heavy"))
        results.append((C + "Fuzz.NewHeavyTests", "T2", "NotExecuted", "OneTimeSetUp: heavy"))
        report = self.account(text, results)
        self.assertEqual(len(self.failures_of(report, "UNLISTED SKIP")), 2, report.failures)
        self.assertEqual(report.proposals, [f"explicit-class\t{C}Fuzz.NewHeavyTests\tboth\theavy\tTODO"])

    def test_the_proposal_of_a_test_that_needs_the_checked_build_is_a_default_row_and_any_other_skip_a_both_row(self):
        text, results = run_rows()
        results.append((C + "External.NewBounds", "ByTheBuild", "NotExecuted", "only the checked build keeps the list"))
        results.append((C + "External.NewBounds", "ByTheFlag", "NotExecuted", "needs -p:EsaBounds=true"))
        results.append((C + "External.NewBounds", "ByTheBackend", "NotExecuted", "the backend has no such surface"))
        report = self.account(text, results)
        rows = {p.split("\t")[1].rsplit(".", 1)[1]: p.split("\t")[2] for p in report.proposals}
        self.assertEqual(rows, {"ByTheBuild": "default", "ByTheFlag": "default", "ByTheBackend": "both"})

    def test_the_failure_names_the_reason_and_proposes_the_row(self):
        text, results = run_rows()
        results = [(c, n, "NotExecuted", "because of reasons") if n == "AnOrdinaryTest" else (c, n, o, m) for c, n, o, m in results]
        report = self.account(text, results)
        self.assertTrue(any("because of reasons" in f for f in report.failures))
        self.assertEqual(report.proposals, [f"ignore\t{C}External.BoundsTests.AnOrdinaryTest\tboth\tbecause of reasons\tTODO"])


class AStaleRow(Base):
    def test_a_both_row_whose_test_ran_is_stale(self):
        text, results = run_rows()
        results = [(c, n, "Passed", "") if n == "NoNavigation" else (c, n, o, m) for c, n, o, m in results]
        failures = self.failures_of(self.account(text, results), "STALE")
        self.assertEqual(len(failures), 1)
        self.assertIn("ukkonen-adapter", failures[0])

    def test_a_default_row_whose_test_ran_is_stale_in_a_default_run(self):
        text, results = run_rows()
        results = [(c, n, "Passed", "") if n == "OnlyChecked1" else (c, n, o, m) for c, n, o, m in results]
        self.assertEqual(len(self.failures_of(self.account(text, results), "STALE")), 1)

    def test_a_default_row_that_is_still_skipped_in_a_checked_run_is_the_failure_the_job_exists_for(self):
        text, results = run_rows()  # the checked run, but the bounds tests were skipped: the build was not checked
        failures = self.failures_of(self.account(text, results, mode="checked"), "CHECKED-BUILD TEST SKIPPED")
        self.assertEqual(len(failures), 2)

    def test_a_both_row_that_ran_in_a_checked_run_is_stale(self):
        text, results = run_rows()
        results = [(c, n, "Passed", "") if n.startswith("OnlyChecked") or n == "NoNavigation" else (c, n, o, m) for c, n, o, m in results]
        self.assertEqual(len(self.failures_of(self.account(text, results, mode="checked"), "STALE")), 1)

    def test_an_explicit_row_whose_test_ran_is_stale(self):
        text, results = run_rows()
        results = [(c, n, "Passed", "") if n == "ReportOnly" else (c, n, o, m) for c, n, o, m in results]
        self.assertEqual(len(self.failures_of(self.account(text, results), "STALE")), 1)

    def test_an_explicit_class_with_one_test_that_ran_is_stale(self):
        text, results = run_rows()
        results = [(c, n, "Passed", "") if c.endswith("(esa-file-w5)") and n == "B" else (c, n, o, m) for c, n, o, m in results]
        failures = self.failures_of(self.account(text, results), "STALE")
        self.assertEqual(len(failures), 1)
        self.assertIn("esa-file-w5", failures[0])

    def test_a_row_with_two_results_of_one_name_needs_both_skipped(self):
        text, results = run_rows()
        results.append((C + "Contract.ConformanceTests(ukkonen-adapter)", "NoNavigation", "Passed", ""))
        self.assertEqual(len(self.failures_of(self.account(text, results), "STALE")), 1)

    def test_a_row_without_a_result_is_absent_unless_the_run_is_partial(self):
        text, results = run_rows()
        results = [r for r in results if r[1] != "NoNavigation" and r[1] != "NeedsMount" and r[1] != "ReportOnly" and not r[0].endswith("(esa-file-w5)")]
        report = self.account(text, results)
        self.assertEqual(len(self.failures_of(report, "LISTED BUT ABSENT")), 4, report.failures)
        self.assertEqual(self.account(text, results, partial=True).failures, [])

    def test_a_partial_run_still_fails_on_an_unlisted_skip_and_on_a_stale_row(self):
        text, results = run_rows()
        results = [r for r in results if r[1] != "ReportOnly"]
        results = [(c, n, "NotExecuted", "x") if n == "AnOrdinaryTest" else (c, n, o, m) for c, n, o, m in results]
        results = [(c, n, "Passed", "") if n == "NoNavigation" else (c, n, o, m) for c, n, o, m in results]
        report = self.account(text, results, partial=True)
        self.assertEqual(len(self.failures_of(report, "UNLISTED SKIP")), 1)
        self.assertEqual(len(self.failures_of(report, "STALE")), 1)


class TheList(Base):
    def load(self, text):
        return sa.load_list(self.write("l.tsv", text))

    def test_a_duplicate_test_name_is_refused(self):
        row = ["ignore", C + "A.B.C", "both", "r", "never"]
        with self.assertRaises(sa.ListError) as e:
            self.load("\n".join([HEADER, "\t".join(row), "\t".join(row)]) + "\n")
        self.assertIn("duplicate", str(e.exception))

    def test_the_same_name_with_two_kinds_is_a_duplicate(self):
        with self.assertRaises(sa.ListError):
            self.load("\n".join([HEADER, "\t".join(["ignore", "A.B", "both", "r", "n"]), "\t".join(["explicit", "A.B", "both", "r", "n"])]) + "\n")

    def test_a_wrong_header_a_wrong_column_count_and_unknown_words_are_refused(self):
        good = ["ignore", "A.B", "both", "r", "n"]
        for text in (
            "kind\ttest\tbuilds\treason\thow\n" + "\t".join(good) + "\n",
            HEADER + "\n" + "\t".join(good[:4]) + "\n",
            HEADER + "\n" + "\t".join(good + ["x"]) + "\n",
            HEADER + "\n" + "\t".join(["skip", "A.B", "both", "r", "n"]) + "\n",
            HEADER + "\n" + "\t".join(["ignore", "A.B", "sometimes", "r", "n"]) + "\n",
            HEADER + "\n" + "\t".join(["explicit", "A.B", "default", "r", "n"]) + "\n",
            HEADER + "\n" + "\t".join(["ignore", "", "both", "r", "n"]) + "\n",
            HEADER + "\n" + "\t".join(["ignore", "A.B", "both", "", "n"]) + "\n",
            HEADER + "\n" + "\t".join(["ignore", "A.B", "both", "r", ""]) + "\n",
            "# only comments\n",
            "",
        ):
            with self.subTest(text=text[:60]):
                with self.assertRaises(sa.ListError):
                    self.load(text)

    def test_comments_and_blank_lines_are_ignored_and_the_rows_are_typed(self):
        rows = self.load("# c\n\n" + HEADER + "\n# between\n" + "\t".join(["ignore", "A.B", "default", "r", "n"]) + "\n\n" + "\t".join(["explicit-class", "A.C(x)", "both", "r", "n"]) + "\n")
        self.assertEqual([(r.kind, r.test, r.in_builds) for r in rows], [("ignore", "A.B", "default"), ("explicit-class", "A.C(x)", "both")])

    def test_the_committed_list_has_the_skips_of_the_default_run(self):
        rows = sa.load_list(REAL_LIST)
        ignore = [r for r in rows if r.kind == "ignore"]
        self.assertEqual(len([r for r in ignore if r.in_builds in ("default", "both")]), 37, "the 37 skipped tests of the default run")
        self.assertEqual(len([r for r in ignore if r.in_builds == "both"]), 18, "the 18 that the checked build skips as well")
        self.assertEqual(len([r for r in ignore if r.in_builds == "default"]), 19, "the 19 that need the checked accessor")
        self.assertEqual(len([r for r in ignore if r.in_builds == "host"]), 9, "the tests that mount a tmpfs")
        self.assertEqual(len([r for r in rows if r.kind == "explicit"]), 4)
        self.assertEqual(len([r for r in rows if r.kind == "explicit-class"]), 44, "40 at the end of S2.2, plus the fixture GeneratedViewLargeTests of S2.3, plus WideEngineSparseLargeTests and WideEngineLargeTests and CorruptionCorpusS24Tests of S2.4")


class TheTrx(Base):
    def test_an_empty_run_a_broken_file_and_a_foreign_xml_are_refused(self):
        empty = make_trx([])
        for name, text in (("empty.trx", empty), ("broken.trx", "<TestRun"), ("foreign.trx", "<html/>"), ("none.trx", None)):
            with self.subTest(name=name):
                path = self.write(name, text) if text is not None else os.path.join(self.dir.name, "missing.trx")
                with self.assertRaises(trx_lib.TrxError):
                    trx_lib.read_trx(path)

    def test_the_names_with_arguments_quotes_and_markup_are_read_back_as_they_are(self):
        results = [(C + "Contract.IndexConformanceTests(esa-file-w5)", 'SaAndLcp("a<b>&c","d,e")', "NotExecuted", 'reason <with> "markup" & more')]
        parsed = trx_lib.read_trx(self.write("a.trx", make_trx(results)))
        self.assertEqual(parsed[0].class_name, C + "Contract.IndexConformanceTests(esa-file-w5)")
        self.assertEqual(parsed[0].test_name, 'SaAndLcp("a<b>&c","d,e")')
        self.assertEqual(parsed[0].message, 'reason <with> "markup" & more')
        self.assertEqual(parsed[0].key, C + 'Contract.IndexConformanceTests(esa-file-w5).SaAndLcp("a<b>&c","d,e")')
        self.assertEqual(parsed[0].outcome, "NotExecuted")

    def test_a_result_without_a_definition_is_refused(self):
        text = make_trx([("A", "B", "Passed", "")]).replace("<TestDefinitions>", "<TestDefinitions><!--").replace("</TestDefinitions>", "--></TestDefinitions>")
        with self.assertRaises(trx_lib.TrxError):
            trx_lib.read_trx(self.write("nodef.trx", text))

    def test_the_wall_clock_and_the_durations(self):
        path = self.write("a.trx", make_trx([("A", "B", "Passed", "")]))
        self.assertAlmostEqual(trx_lib.wall_seconds(path), 600.0)
        self.assertAlmostEqual(trx_lib.read_trx(path)[0].seconds, 0.001)

    def test_only_the_outcomes_that_ran_count_as_executed(self):
        self.assertEqual(sorted(trx_lib.EXECUTED), sorted(["Passed", "Failed", "Error", "Timeout", "Aborted", "PassedButRunAborted"]))


class TheCommandLine(Base):
    def run_cli(self, *args):
        p = subprocess.run([sys.executable, os.path.join(HERE, "skip_accounting.py"), *args], capture_output=True, text=True)
        return p.returncode, p.stdout + p.stderr

    def test_exit_codes_and_the_output(self):
        text, results = run_rows()
        lst = self.write("skips.tsv", text)
        trx = self.write("run.trx", make_trx(results))
        code, out = self.run_cli("--mode", "default", "--list", lst, "--trx", trx)
        self.assertEqual(code, 0, out)
        self.assertIn("holds", out)
        bad = self.write("bad.trx", make_trx([(c, n, "NotExecuted", "x") if n == "AnOrdinaryTest" else (c, n, o, m) for c, n, o, m in results]))
        code, out = self.run_cli("--mode", "default", "--list", lst, "--trx", bad)
        self.assertEqual(code, 1, out)
        self.assertIn("UNLISTED SKIP", out)
        self.assertIn("AnOrdinaryTest", out)
        code, out = self.run_cli("--mode", "default", "--list", lst, "--trx", bad, "--propose")
        self.assertEqual(code, 1, out)
        self.assertIn("ignore\t" + C + "External.BoundsTests.AnOrdinaryTest\tboth\tx\tTODO", out)
        code, out = self.run_cli("--mode", "checked", "--list", lst, "--trx", trx)
        self.assertEqual(code, 1, out)
        self.assertIn("CHECKED-BUILD TEST SKIPPED", out)
        code, out = self.run_cli("--mode", "default", "--list", lst, "--trx", os.path.join(self.dir.name, "missing.trx"))
        self.assertEqual(code, 2, out)
        code, out = self.run_cli("--mode", "default", "--list", self.write("empty.tsv", "# nothing\n"), "--trx", trx)
        self.assertEqual(code, 2, out)
        code, out = self.run_cli("--mode", "nightly", "--list", lst, "--trx", trx)
        self.assertEqual(code, 2, out)
        partial = self.write("partial.trx", make_trx([r for r in results if r[1] != "ReportOnly"]))
        self.assertEqual(self.run_cli("--mode", "default", "--list", lst, "--trx", partial)[0], 1)
        self.assertEqual(self.run_cli("--mode", "default", "--list", lst, "--trx", partial, "--partial")[0], 0)


if __name__ == "__main__":
    unittest.main()
