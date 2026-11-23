#!/usr/bin/env python3
"""Self-test of the prose tool-count checker of gen-catalog.py (audit R14-5, B09): a front page that quotes a stale total must be reported.

  python3 -m unittest scripts/skills/test_gen_catalog_counts.py
"""
import importlib.util
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SPEC = importlib.util.spec_from_file_location("gen_catalog", os.path.join(HERE, "gen-catalog.py"))
GEN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GEN)


class ProseTotalMentions(unittest.TestCase):
    def numbers(self, text):
        return [n for n, _ in GEN.prose_total_mentions(text)]

    def test_every_quoted_form_of_the_total_is_found(self):
        text = "\n".join([
            "[![MCP](https://img.shields.io/badge/MCP-483%20tools-6f42c1)](#3)",
            "chain **483 strict-schema tools** on demand",
            "With **504 tools across 11 servers**",
            "without loading all 504 MCP tool schemas.",
            "Registering all 474 tools would load every schema",
            "(find the right tool among 486 without loading schemas)",
            "the ledger of the 642 MCP tools",
            "~483 tools",
        ])
        self.assertEqual(self.numbers(text), [483, 483, 504, 504, 474, 486, 642, 483])

    def test_unrelated_numbers_are_not_counts(self):
        text = "11 servers, 247 algorithm docs, 170,370 source lines (~483x less), 21 skills, 427 reports, 25 339 tests"
        self.assertEqual(self.numbers(text), [])

    def test_stale_total_is_reported_and_the_right_one_is_not(self):
        import tempfile
        with tempfile.TemporaryDirectory() as root:
            for rel, text in (("README.md", "With **483 tools** the model drowns"), (".claude/skills/README.md", "With **642 tools** it does not")):
                path = os.path.join(root, *rel.split("/"))
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "w", encoding="utf-8") as f:
                    f.write(text)
            saved = GEN.REPO_ROOT
            GEN.REPO_ROOT = root
            try:
                problems = GEN.check_prose_total_counts(642)
            finally:
                GEN.REPO_ROOT = saved
        self.assertEqual(problems, ["README.md:1 says 483 tools, the tool docs hold 642"])


if __name__ == "__main__":
    unittest.main()
