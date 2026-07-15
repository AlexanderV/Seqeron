from __future__ import annotations

import sys
import subprocess
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import wiki_graph_extract  # noqa: E402
import wiki_graph_lint  # noqa: E402
import wiki_lint  # noqa: E402
import wiki_markdown  # noqa: E402
import wiki_search  # noqa: E402
import wiki_search_benchmark  # noqa: E402
import wiki_stats  # noqa: E402


class MarkdownLinkTests(unittest.TestCase):
    def test_extracts_aliases_and_escaped_table_aliases(self):
        body = "[[plain]] [[shown|Alias]] [[table\\|Table alias]] [[#local]]"
        self.assertEqual(
            ["plain", "shown", "table"],
            wiki_markdown.extract_wikilinks(body),
        )

    def test_ignores_inline_code(self):
        body = "`[[one]]` ``code ` [[two]]`` [[real]]"
        self.assertEqual(["real"], wiki_markdown.extract_wikilinks(body))

    def test_ignores_four_character_fence_containing_triple_fence(self):
        body = "````md\n[[fake]]\n```\n````\n[[real]]\n"
        self.assertEqual(["real"], wiki_markdown.extract_wikilinks(body))

    def test_accepts_longer_closing_fence(self):
        body = "```md\n[[fake]]\n`````\n[[real]]\n"
        self.assertEqual(["real"], wiki_markdown.extract_wikilinks(body))

    def test_ignores_tilde_fence(self):
        body = "~~~~\n[[fake]]\n~~~~\n[[real]]\n"
        self.assertEqual(["real"], wiki_markdown.extract_wikilinks(body))

    def test_unclosed_fence_consumes_rest_of_document(self):
        body = "[[before]]\n```\n[[inside]]\n"
        self.assertEqual(["before"], wiki_markdown.extract_wikilinks(body))

    def test_indented_code_and_tab_indented_code_are_ignored(self):
        body = "    [[space-code]]\n\n\t[[tab-code]]\n\n[[real]]"
        self.assertEqual(["real"], wiki_markdown.extract_wikilinks(body))

    def test_indentation_does_not_interrupt_a_paragraph(self):
        body = "paragraph\n    [[continuation]]\n"
        self.assertEqual(["continuation"], wiki_markdown.extract_wikilinks(body))

    def test_fences_inside_blockquotes_and_lists_are_ignored(self):
        body = (
            "> ~~~~\n> [[quoted]]\n> ~~~~~\n\n"
            "10. item\n    ```\n    [[listed]]\n    ````\n\n"
            "[[real]]"
        )
        self.assertEqual(["real"], wiki_markdown.extract_wikilinks(body))

    def test_fences_in_mixed_list_and_blockquote_containers_are_ignored(self):
        body = (
            "- > ```\n  > [[list-quote]]\n  > ```\n\n"
            "> - ~~~\n>   [[quote-list]]\n>   ~~~\n\n"
            "[[real]]"
        )
        self.assertEqual(["real"], wiki_markdown.extract_wikilinks(body))

    def test_excess_list_padding_creates_indented_code(self):
        body = "-     [[indented-code]]\n\n- [[visible-item]]"
        self.assertEqual(["visible-item"], wiki_markdown.extract_wikilinks(body))

    def test_unclosed_container_fence_ends_with_its_container(self):
        body = "> ```\n> [[quoted]]\n\n[[outside]]"
        self.assertEqual(["outside"], wiki_markdown.extract_wikilinks(body))

    def test_code_spans_do_not_cross_paragraph_boundaries(self):
        body = "`open\n\n[[visible]]\n\nclose` [[also-visible]]"
        self.assertEqual(
            ["visible", "also-visible"],
            wiki_markdown.extract_wikilinks(body),
        )

    def test_escaped_backticks_do_not_create_code_spans(self):
        body = r"\`[[visible]]\` [[real]]"
        self.assertEqual(["visible", "real"], wiki_markdown.extract_wikilinks(body))

    def test_pre_and_comment_blocks_are_ignored(self):
        body = (
            "<pre>\n[[pre-code]]\n</pre>\n\n"
            "<!--\n[[comment]]\n-->\n\n"
            "[[real]]"
        )
        self.assertEqual(["real"], wiki_markdown.extract_wikilinks(body))

    def test_inline_html_comments_are_ignored(self):
        body = "before <!-- [[comment]] --> [[real]]"
        self.assertEqual(["real"], wiki_markdown.extract_wikilinks(body))

    def test_cross_page_anchor_resolves_to_page_and_multiline_link_is_rejected(self):
        body = "[[page#section|Section]] [[broken\nlink]] [[real]]"
        self.assertEqual(["page", "real"], wiki_markdown.extract_wikilinks(body))

    def test_root_and_blockquote_fence_matrix(self):
        for prefix in ("", "> ", "> > "):
            for character in ("`", "~"):
                for opening_length in range(3, 7):
                    for extra_closing in range(3):
                        with self.subTest(prefix=prefix, character=character,
                                          opening=opening_length, extra=extra_closing):
                            body = (
                                f"{prefix}{character * opening_length}\n"
                                f"{prefix}[[fake]]\n"
                                f"{prefix}{character * (opening_length + extra_closing)}\n\n"
                                "[[real]]"
                            )
                            self.assertEqual(["real"], wiki_markdown.extract_wikilinks(body))

    def test_all_tools_share_the_same_extractor(self):
        for module in (
            wiki_lint,
            wiki_search,
            wiki_graph_extract,
            wiki_graph_lint,
            wiki_stats,
        ):
            self.assertIs(wiki_markdown.extract_wikilinks, module.extract_wikilinks)


class ProvenanceDedupTests(unittest.TestCase):
    @staticmethod
    def page(slug: str, page_type: str, source: str = "docs/shared.md") -> dict:
        return {
            "slug": slug,
            "meta": {"type": page_type, "sources": [source]},
        }

    def test_distinct_concepts_with_same_source_are_preserved(self):
        scored = [
            (10.0, self.page("concept-a", "concept")),
            (9.0, self.page("concept-b", "concept")),
        ]
        result = wiki_search.collapse_by_provenance(scored, "concept")
        self.assertEqual(["concept-a", "concept-b"], [page["slug"] for _, page in result])

    def test_source_variant_is_hidden_when_derived_page_exists(self):
        scored = [
            (10.0, self.page("raw-source", "source")),
            (8.0, self.page("derived", "concept")),
        ]
        result = wiki_search.collapse_by_provenance(scored, "concept")
        self.assertEqual(["derived"], [page["slug"] for _, page in result])
        self.assertEqual(10.0, result[0][0])

    def test_source_score_is_transferred_without_collapsing_derived_pages(self):
        scored = [
            (100.0, self.page("raw-source", "source")),
            (2.0, self.page("concept-a", "concept")),
            (1.0, self.page("concept-b", "concept")),
        ]
        result = wiki_search.collapse_by_provenance(scored, "concept")
        self.assertEqual(
            [(100.0, "concept-a"), (1.0, "concept-b")],
            [(score, page["slug"]) for score, page in result],
        )

    def test_preferred_type_becomes_provenance_representative(self):
        scored = [
            (20.0, self.page("raw-source", "source")),
            (8.0, self.page("api-page", "api")),
            (3.0, self.page("concept-page", "concept")),
        ]
        result = wiki_search.collapse_by_provenance(scored, "concept")
        self.assertEqual(
            [(20.0, "concept-page"), (8.0, "api-page")],
            [(score, page["slug"]) for score, page in result],
        )

    def test_best_source_is_kept_when_group_has_only_sources(self):
        scored = [
            (7.0, self.page("source-a", "source")),
            (9.0, self.page("source-b", "source")),
        ]
        result = wiki_search.collapse_by_provenance(scored, None)
        self.assertEqual(["source-b"], [page["slug"] for _, page in result])

    def test_preference_does_not_copy_score_between_derived_pages(self):
        scored = [
            (9.0, self.page("api-page", "api")),
            (1.0, self.page("concept-page", "concept")),
        ]
        result = wiki_search.collapse_by_provenance(scored, "concept")
        self.assertEqual(
            [(9.0, "api-page"), (1.0, "concept-page")],
            [(score, page["slug"]) for score, page in result],
        )

    def test_lower_scoring_source_does_not_inflate_preferred_page(self):
        scored = [
            (9.0, self.page("api-page", "api")),
            (2.0, self.page("raw-source", "source")),
            (1.0, self.page("concept-page", "concept")),
        ]
        result = wiki_search.collapse_by_provenance(scored, "concept")
        self.assertEqual(
            [(9.0, "api-page"), (2.0, "concept-page")],
            [(score, page["slug"]) for score, page in result],
        )


class RetrievalBenchmarkTests(unittest.TestCase):
    def write_concept(
        self,
        root: Path,
        slug: str,
        title: str,
        body: str,
        sources: tuple[str, ...] = (),
    ) -> None:
        source_lines = "sources:\n" + "".join(f"  - {source}\n" for source in sources) if sources else ""
        (root / f"{slug}.md").write_text(
            "---\n"
            "type: concept\n"
            f'title: "{title}"\n'
            f"{source_lines}"
            "---\n\n"
            f"{body}\n",
            encoding="utf-8",
        )

    def test_bilingual_ranks_measure_direct_and_normalized_queries_separately(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_concept(root, "primer-design", "Primer design", "PCR primer pair constraints")
            self.write_concept(root, "tree-build", "Tree construction", "phylogenetic distance tree")
            cases = [{
                "id": "primer",
                "target": "primer-design",
                "en": "design PCR primer pair",
                "uk": "спроєктувати пару праймерів",
                "uk_normalized": "design PCR primer pair",
            }]

            ranks = wiki_search_benchmark.rank_targets(root, cases)

            self.assertEqual([1], ranks["en"])
            self.assertEqual([None], ranks["uk"])
            self.assertEqual([1], ranks["uk_normalized"])

    def test_comparison_uses_concept_sources_as_without_wiki_gold(self):
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory)
            wiki = repo / "wiki"
            docs = repo / "docs"
            wiki.mkdir()
            docs.mkdir()
            (docs / "primer.md").write_text("PCR primer pair constraints\n", encoding="utf-8")
            (docs / "trees.md").write_text("phylogenetic distance tree\n", encoding="utf-8")
            self.write_concept(
                wiki,
                "primer-design",
                "Primer design",
                "PCR primer pair constraints",
                ("docs/primer.md",),
            )
            cases = [{
                "id": "primer",
                "target": "primer-design",
                "en": "PCR primer pair constraints",
                "uk": "обмеження для пари праймерів",
                "uk_normalized": "PCR primer pair constraints",
            }]

            comparison = wiki_search_benchmark.rank_comparison(repo, wiki, cases)

            self.assertEqual([1], comparison["without_wiki"]["en"])
            self.assertEqual([1], comparison["with_wiki"]["en"])
            self.assertEqual([None], comparison["without_wiki"]["uk"])
            self.assertEqual([None], comparison["with_wiki"]["uk"])

    def test_summary_reports_binary_hit_at_each_cutoff(self):
        ranks = {
            "en": [1, 2, None],
            "uk": [3, 10, 11],
            "uk_normalized": [1, 1, 1],
        }

        rows = {row["field"]: row for row in wiki_search_benchmark.summarize(ranks)}

        self.assertEqual((1, 2, 2), tuple(rows["en"][f"hits_at_{k}"] for k in (1, 3, 10)))
        self.assertEqual((0, 1, 2), tuple(rows["uk"][f"hits_at_{k}"] for k in (1, 3, 10)))
        self.assertEqual((3, 3, 3), tuple(rows["uk_normalized"][f"hits_at_{k}"] for k in (1, 3, 10)))

    def test_fixed_benchmark_has_unique_complete_cases(self):
        benchmark = Path(__file__).resolve().parents[1] / "benchmarks" / "retrieval_queries.jsonl"
        cases = wiki_search_benchmark.load_cases(benchmark)

        self.assertEqual(30, len(cases))
        self.assertEqual(30, len({case["id"] for case in cases}))
        self.assertEqual(30, len({case["target"] for case in cases}))


class IndexLimitTests(unittest.TestCase):
    def lint_index(self, line_count: int, trailing_newline: bool = True) -> dict:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            indexes = root / "indexes" / "nested"
            indexes.mkdir(parents=True)
            text = "\n".join(["entry"] * line_count)
            if trailing_newline:
                text += "\n"
            (indexes / "concepts.md").write_text(text, encoding="utf-8")
            return wiki_lint.lint([], 400, 800, [], False, 5, wiki_root=root)

    def test_exactly_300_lines_is_allowed_even_with_trailing_newline(self):
        self.assertEqual([], self.lint_index(300)["oversized_indexes"])

    def test_301_lines_is_reported(self):
        findings = self.lint_index(301)["oversized_indexes"]
        self.assertEqual(1, len(findings))
        self.assertEqual(301, findings[0]["lines"])
        self.assertEqual("indexes/nested/concepts.md", findings[0]["path"].replace("\\", "/"))

    def test_root_index_is_subject_to_the_same_limit(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "index.md").write_text(
                "\n".join(["entry"] * 301) + "\n",
                encoding="utf-8",
            )
            findings = wiki_lint.lint([], 400, 800, [], False, 5, wiki_root=root)
            self.assertEqual(
                [("index.md", 301)],
                [(item["path"].replace("\\", "/"), item["lines"])
                 for item in findings["oversized_indexes"]],
            )

    def test_cli_fails_for_oversized_index(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "indexes").mkdir()
            (root / "indexes" / "concepts.md").write_text(
                "\n".join(["entry"] * 301) + "\n",
                encoding="utf-8",
            )
            result = subprocess.run(
                [sys.executable, str(SCRIPTS / "wiki_lint.py"), str(root)],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(1, result.returncode)
            self.assertIn("Index over cap", result.stdout)


class StalenessHistoryTests(unittest.TestCase):
    def test_atomic_source_and_page_refresh_is_fresh_but_later_source_change_is_stale(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            wiki = root / "wiki"
            wiki.mkdir()

            def git(*args: str) -> subprocess.CompletedProcess[str]:
                return subprocess.run(
                    ["git", *args],
                    cwd=root,
                    text=True,
                    capture_output=True,
                    check=True,
                )

            git("init", "--quiet")
            git("config", "user.email", "wiki-tests@example.invalid")
            git("config", "user.name", "Wiki Tests")
            (root / "README.md").write_text("baseline\n", encoding="utf-8")
            (wiki / "readme.md").write_text("# placeholder\n", encoding="utf-8")
            git("add", ".")
            git("commit", "--quiet", "-m", "baseline")
            baseline = git("rev-parse", "HEAD").stdout.strip()
            command = [
                sys.executable,
                str(SCRIPTS / "wiki_stale.py"),
                "wiki",
                "--repo-root",
                str(root),
            ]

            (root / "README.md").write_text("refreshed source\n", encoding="utf-8")
            (wiki / "readme.md").write_text(
                "---\n"
                "type: source\n"
                "sources:\n"
                "  - README.md\n"
                f"source_commit: {baseline}\n"
                "---\n\n"
                "# Refreshed summary\n",
                encoding="utf-8",
            )
            git("add", ".")
            staged_together = subprocess.run(
                command,
                cwd=root,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(
                0,
                staged_together.returncode,
                staged_together.stdout + staged_together.stderr,
            )
            git("commit", "--quiet", "-m", "atomic source and wiki refresh")

            fresh = subprocess.run(command, cwd=root, text=True, capture_output=True, check=False)
            self.assertEqual(0, fresh.returncode, fresh.stdout + fresh.stderr)

            (root / "README.md").write_text("new source-only change\n", encoding="utf-8")
            git("add", "README.md")
            staged_source_only = subprocess.run(
                command,
                cwd=root,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(
                1,
                staged_source_only.returncode,
                staged_source_only.stdout + staged_source_only.stderr,
            )
            self.assertIn("staged without a staged page refresh", staged_source_only.stdout)
            git("commit", "--quiet", "-m", "source only")
            stale = subprocess.run(command, cwd=root, text=True, capture_output=True, check=False)
            self.assertEqual(1, stale.returncode, stale.stdout + stale.stderr)
            self.assertIn("README.md has 2 commit(s)", stale.stdout)


class GraphLintExitTests(unittest.TestCase):
    def test_cli_fails_cleanly_for_broken_typed_edge(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "graph").mkdir()
            repository_root = Path(__file__).resolve().parents[4]
            (root / "graph" / "ontology.yaml").write_text(
                (repository_root / "wiki" / "graph" / "ontology.yaml").read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            (root / "broken.md").write_text(
                """---
type: concept
title: Broken edge
graph:
  relationships:
    - predicate: relates_to
      object: concept:missing
      source: missing-source
      evidence: deliberate test fixture
      confidence: high
      status: current
---

# Broken edge
""",
                encoding="utf-8",
            )
            result = subprocess.run(
                [sys.executable, str(SCRIPTS / "wiki_graph_lint.py"), str(root)],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(1, result.returncode)
            self.assertIn("Broken object references", result.stdout)
            self.assertIn("source: does not match any source page", result.stdout)
            self.assertNotIn("UnicodeEncodeError", result.stderr)


class WikiStatsContractTests(unittest.TestCase):
    def test_trailing_newline_does_not_inflate_reported_line_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "index.md").write_text("# Index\n\nSecond line\n", encoding="utf-8")
            (root / "page.md").write_text(
                "---\ntype: concept\ntitle: Page\n---\n\n# Page\n",
                encoding="utf-8",
            )
            result = subprocess.run(
                [sys.executable, str(SCRIPTS / "wiki_stats.py"), str(root)],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(0, result.returncode)
            self.assertRegex(result.stdout, r"Total lines:\s+6\b")
            self.assertRegex(result.stdout, r"index\.md:\s+3 lines\b")


if __name__ == "__main__":
    unittest.main()
