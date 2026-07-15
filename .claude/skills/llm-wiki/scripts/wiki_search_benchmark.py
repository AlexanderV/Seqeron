#!/usr/bin/env python3
"""Measure BM25 Hit@K on a fixed bilingual retrieval benchmark."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from wiki_markdown import configure_utf8_streams
from wiki_search import bm25_score, build_bm25, collect_pages, tokenize


configure_utf8_streams()


LANGUAGES = (
    ("en", "English (direct)"),
    ("uk", "Українська (direct)"),
    ("uk_normalized", "Українська → English normalization"),
)
SURFACES = (
    ("without_wiki", "Without the LLM Wiki"),
    ("with_wiki", "With the LLM Wiki"),
)
K_VALUES = (1, 3, 10)


def load_cases(path: Path) -> list[dict]:
    cases = []
    seen_ids = set()
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            case = json.loads(line)
        except json.JSONDecodeError as error:
            raise ValueError(f"{path}:{line_number}: invalid JSON: {error}") from error
        missing = {"id", "target", *(field for field, _ in LANGUAGES)} - case.keys()
        if missing:
            raise ValueError(f"{path}:{line_number}: missing fields: {', '.join(sorted(missing))}")
        required = ("id", "target", *(field for field, _ in LANGUAGES))
        empty = [field for field in required if not isinstance(case[field], str) or not case[field].strip()]
        if empty:
            raise ValueError(f"{path}:{line_number}: empty fields: {', '.join(empty)}")
        if case["id"] in seen_ids:
            raise ValueError(f"{path}:{line_number}: duplicate id: {case['id']}")
        seen_ids.add(case["id"])
        cases.append(case)
    if not cases:
        raise ValueError(f"{path}: benchmark contains no cases")
    return cases


def _rank_queries(
    pages: list[dict],
    cases: list[dict],
    gold_by_case: list[set[str]],
) -> dict[str, list[int | None]]:
    index = build_bm25(pages)
    ranks: dict[str, list[int | None]] = {field: [] for field, _ in LANGUAGES}
    for case, gold in zip(cases, gold_by_case, strict=True):
        for field, _ in LANGUAGES:
            query_tokens = tokenize(case[field])
            scored = [
                (bm25_score(query_tokens, page_index, index), page["slug"])
                for page_index, page in enumerate(pages)
            ]
            ranked = [
                slug
                for score, slug in sorted(scored, key=lambda item: (-item[0], item[1]))
                if score > 0
            ]
            rank = next((position for position, slug in enumerate(ranked, 1) if slug in gold), None)
            ranks[field].append(rank)
    return ranks


def _concept_pages(wiki: Path) -> list[dict]:
    pages = [page for page in collect_pages(wiki) if page["meta"].get("type") == "concept"]
    pages.sort(key=lambda page: page["slug"])
    return pages


def rank_targets(wiki: Path, cases: list[dict]) -> dict[str, list[int | None]]:
    pages = _concept_pages(wiki)
    available = {page["slug"] for page in pages}
    missing_targets = sorted({case["target"] for case in cases} - available)
    if missing_targets:
        raise ValueError(f"benchmark targets missing from concept corpus: {', '.join(missing_targets)}")
    return _rank_queries(pages, cases, [{case["target"]} for case in cases])


def _source_pages(repo_root: Path) -> list[dict]:
    paths = sorted((repo_root / "docs").rglob("*.md")) + sorted(repo_root.glob("*.md"))
    pages = []
    for path in paths:
        text = path.read_text(encoding="utf-8")
        relative = path.relative_to(repo_root).as_posix()
        pages.append({"slug": relative, "tokens": tokenize(f"{text} {relative}")})
    return pages


def rank_source_targets(repo_root: Path, wiki: Path, cases: list[dict]) -> dict[str, list[int | None]]:
    concepts = {page["slug"]: page for page in _concept_pages(wiki)}
    source_pages = _source_pages(repo_root)
    available_sources = {page["slug"] for page in source_pages}
    gold_by_case = []
    for case in cases:
        concept = concepts.get(case["target"])
        if concept is None:
            raise ValueError(f"benchmark target missing from concept corpus: {case['target']}")
        sources = concept["meta"].get("sources", []) or []
        if isinstance(sources, str):
            sources = [sources]
        gold = {str(source).replace("\\", "/") for source in sources} & available_sources
        if not gold:
            raise ValueError(f"benchmark target has no local source documents: {case['target']}")
        gold_by_case.append(gold)
    return _rank_queries(source_pages, cases, gold_by_case)


def rank_comparison(repo_root: Path, wiki: Path, cases: list[dict]) -> dict[str, dict]:
    return {
        "without_wiki": rank_source_targets(repo_root, wiki, cases),
        "with_wiki": rank_targets(wiki, cases),
    }


def summarize(ranks: dict[str, list[int | None]]) -> list[dict]:
    rows = []
    for field, label in LANGUAGES:
        values = ranks[field]
        row = {"field": field, "label": label, "queries": len(values)}
        for k in K_VALUES:
            hits = sum(rank is not None and rank <= k for rank in values)
            row[f"hits_at_{k}"] = hits
            row[f"hit_at_{k}"] = hits / len(values)
        rows.append(row)
    return rows


def summarize_comparison(comparison: dict[str, dict]) -> list[dict]:
    rows = []
    surface_labels = dict(SURFACES)
    summaries = {
        surface: {row["field"]: row for row in summarize(comparison[surface])}
        for surface, _ in SURFACES
    }
    for field, _ in LANGUAGES:
        for surface, _ in SURFACES:
            row = summaries[surface][field]
            rows.append({"surface": surface, "surface_label": surface_labels[surface], **row})
    return rows


def print_markdown(rows: list[dict]) -> None:
    by_key = {(row["field"], row["surface"]): row for row in rows}
    query_count = rows[0]["queries"]
    print(f"Each cell shows **Without → With the LLM Wiki (absolute gain)**; N = {query_count}.")
    print()
    print("| Query form | Hit@1 | Hit@3 | Hit@10 |")
    print("|---|---:|---:|---:|")
    for field, label in LANGUAGES:
        without = by_key[(field, "without_wiki")]
        with_wiki = by_key[(field, "with_wiki")]
        cells = []
        for k in K_VALUES:
            before = without[f"hit_at_{k}"]
            after = with_wiki[f"hit_at_{k}"]
            gain = (after - before) * 100
            cells.append(f"{before:.1%} → **{after:.1%}** (+{gain:.1f} pp)")
        print(f"| {label} | {' | '.join(cells)} |")


def main() -> None:
    script_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wiki", nargs="?", type=Path, default=Path("wiki"))
    parser.add_argument(
        "--repo-root",
        type=Path,
        help="Repository root for the source-document baseline (default: parent of wiki).",
    )
    parser.add_argument(
        "--queries",
        type=Path,
        default=script_root / "benchmarks" / "retrieval_queries.jsonl",
    )
    parser.add_argument("--json", action="store_true", help="Emit machine-readable results.")
    args = parser.parse_args()

    try:
        cases = load_cases(args.queries)
        repo_root = args.repo_root or args.wiki.resolve().parent
        rows = summarize_comparison(rank_comparison(repo_root, args.wiki, cases))
    except (OSError, ValueError) as error:
        print(error, file=sys.stderr)
        raise SystemExit(1) from error

    if args.json:
        print(json.dumps(rows, ensure_ascii=False, indent=2))
    else:
        print_markdown(rows)


if __name__ == "__main__":
    main()
