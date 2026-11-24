"""Reads a Visual Studio .trx file (the format of `dotnet test --logger trx`) for the CI scripts of work package S2.2. Only the Python standard library.

A result is (class_name, test_name, outcome, message, stdout, seconds):
  * class_name  the class of the test as the trx names it: with the arguments of a parametrised fixture, `Ns.Class(esa-file-w5)`;
  * test_name   the name of the test as the trx names it: with the arguments of a test case, `Method(1,2)`, or the name that SetName gave;
  * outcome     the trx outcome: Passed, Failed, NotExecuted (an `Assert.Ignore`, an `Explicit` test that was not selected, an `Assume`), Inconclusive, ...
  * key         `class_name.test_name`, what the skip list names.
The NUnit adapter writes an ignored test and an Explicit test that was not selected the same way (NotExecuted, with a message): the trx cannot tell them
apart, which is why the skip list names both kinds (skip_accounting.py).
"""
import re
import xml.etree.ElementTree as ET
from typing import NamedTuple

NS = "{http://microsoft.com/schemas/VisualStudio/TeamTest/2010}"
# the outcomes of a test that RAN (it may have failed, but it was executed); everything else (NotExecuted, Inconclusive, Pending, ...) did not pass
EXECUTED = frozenset({"Passed", "Failed", "Error", "Timeout", "Aborted", "PassedButRunAborted"})


class Result(NamedTuple):
    class_name: str
    test_name: str
    outcome: str
    message: str
    stdout: str
    seconds: float

    @property
    def key(self) -> str:
        return f"{self.class_name}.{self.test_name}"


class TrxError(Exception):
    """The file is not a usable trx (unreadable, no results, a result without a definition)."""


def _seconds(duration: str) -> float:
    m = re.fullmatch(r"(?:(\d+)\.)?(\d+):(\d+):(\d+(?:\.\d+)?)", duration or "")
    if not m:
        return 0.0
    days, h, mi, s = m.groups()
    return int(days or 0) * 86400 + int(h) * 3600 + int(mi) * 60 + float(s)


def read_trx(path: str) -> list:
    """All the results of the file, in file order. Raises TrxError for a file that cannot be a test run (no results at all counts: an empty run must not pass)."""
    try:
        root = ET.parse(path).getroot()
    except (ET.ParseError, OSError) as e:
        raise TrxError(f"cannot read {path}: {e}") from e
    if root.tag != NS + "TestRun":
        raise TrxError(f"{path}: the root element is {root.tag!r}, not a TestRun")
    classes = {}
    for unit in root.iter(NS + "UnitTest"):
        method = unit.find(NS + "TestMethod")
        if method is not None:
            classes[unit.get("id")] = method.get("className") or ""
    results = []
    for node in root.iter(NS + "UnitTestResult"):
        test_id = node.get("testId")
        if test_id not in classes:
            raise TrxError(f"{path}: the result of {node.get('testName')!r} has no test definition ({test_id})")
        message = node.find(f"{NS}Output/{NS}ErrorInfo/{NS}Message")
        stdout = node.find(f"{NS}Output/{NS}StdOut")
        results.append(Result(classes[test_id], node.get("testName") or "", node.get("outcome") or "", (message.text or "") if message is not None else "",
                              (stdout.text or "") if stdout is not None else "", _seconds(node.get("duration") or "")))
    if not results:
        raise TrxError(f"{path}: the run has no result at all")
    return results


def wall_seconds(path: str) -> float:
    """The span of the run (Times start to finish) in seconds; 0 when the file does not say."""
    from datetime import datetime
    try:
        times = ET.parse(path).getroot().find(NS + "Times")
        if times is None:
            return 0.0
        # the trx writes seven fractional digits; fromisoformat of Python before 3.11 takes three or six
        trim = lambda t: re.sub(r"(\.\d{6})\d+", r"\1", t)
        start = datetime.fromisoformat(trim(times.get("start")))
        finish = datetime.fromisoformat(trim(times.get("finish")))
        return max(0.0, (finish - start).total_seconds())
    except (ET.ParseError, OSError, ValueError, TypeError):
        return 0.0
