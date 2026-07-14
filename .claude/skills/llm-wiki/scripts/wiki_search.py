#!/usr/bin/env python3
"""
wiki_search.py — BM25 search over wiki pages with frontmatter filters.

Fallback for when index-first navigation doesn't surface the right pages.
Pure-Python implementation (no dependencies beyond stdlib) so it runs anywhere.

Usage:
    python wiki_search.py "query terms" [options]

Options:
    --wiki <dir>            Wiki directory (default: ./wiki)
    --top N                 Return top N results (default: 10)
    --type <type>           Filter by frontmatter type (source|entity|concept|synthesis|...)
    --dedup-provenance      Collapse pages with the same primary source provenance
    --prefer-type <type>    Prefer this page type within a provenance group
    --tag <tag>             Filter by tag (repeatable)
    --since YYYY-MM-DD      Only pages updated on or after this date
    --backlinks <slug>      Find pages that link to <slug>; ignores the query
    --top-linked N          Show the N most-linked-to pages (hubs); ignores the query
    --cache <path>          Persist the BM25 index to disk for faster reruns

Examples:
    python wiki_search.py "diffusion training stability" --top 5
    python wiki_search.py "alignment" --type concept --tag safety
    python wiki_search.py "alignment" --dedup-provenance --prefer-type concept
    python wiki_search.py "" --backlinks transformer
    python wiki_search.py "" --top-linked 10
"""

import argparse
import json
import math
import re
import sys
from collections import Counter, defaultdict
from datetime import date, datetime
from pathlib import Path


WIKILINK_RE = re.compile(r"\[\[([^\]|]+)(?:\|[^\]]+)?\]\]")
FRONTMATTER_RE = re.compile(r"^---\s*\n(.*?)\n---\s*\n", re.DOTALL)
TOKEN_RE = re.compile(r"[a-z0-9]+")


def parse_frontmatter(text: str) -> tuple[dict, str]:
    """Lightweight YAML-ish frontmatter parser. Returns (metadata, body)."""
    m = FRONTMATTER_RE.match(text)
    if not m:
        return {}, text
    fm_text = m.group(1)
    body = text[m.end():]
    meta = {}
    current_key = None
    for line in fm_text.split("\n"):
        if not line.strip():
            continue
        # Inline list: tags: [a, b, c]
        kv = re.match(r"^([a-zA-Z_]+):\s*(.*)$", line)
        if kv:
            key, value = kv.group(1), kv.group(2).strip()
            if value.startswith("[") and value.endswith("]"):
                items = [x.strip().strip('"').strip("'") for x in value[1:-1].split(",") if x.strip()]
                meta[key] = items
            elif value:
                meta[key] = value.strip('"').strip("'")
            else:
                meta[key] = []
                current_key = key
        elif line.startswith("  - ") and current_key:
            meta[current_key].append(line[4:].strip().strip('"').strip("'"))
    return meta, body


def tokenize(text: str) -> list[str]:
    return TOKEN_RE.findall(text.lower())


def slug_from_path(path: Path, wiki_root: Path) -> str:
    return path.stem


def extract_wikilinks(body: str) -> list[str]:
    return [m.group(1).strip() for m in WIKILINK_RE.finditer(body)]


def collect_pages(wiki_root: Path) -> list[dict]:
    """Walk the wiki and return [{path, slug, meta, body, tokens, links}]."""
    pages = []
    for md_path in wiki_root.rglob("*.md"):
        # Skip the schema, index, log, and template files
        rel = md_path.relative_to(wiki_root)
        if rel.parts[0] in {"SCHEMA.md", "index.md", "log.md"} or rel.name.startswith("."):
            continue
        if rel.parts[0] in {"indexes", "graph"}:
            continue
        try:
            text = md_path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        meta, body = parse_frontmatter(text)
        pages.append({
            "path": str(md_path),
            "rel_path": str(rel),
            "slug": slug_from_path(md_path, wiki_root),
            "meta": meta,
            "body": body,
            "tokens": tokenize(body + " " + meta.get("title", "")),
            "links": extract_wikilinks(body),
        })
    return pages


def build_bm25(pages: list[dict]) -> dict:
    """Build a BM25 index. Returns {df, avgdl, N, doc_lens, term_freqs}."""
    N = len(pages)
    df = Counter()
    doc_lens = []
    term_freqs = []
    for page in pages:
        tokens = page["tokens"]
        doc_lens.append(len(tokens))
        tf = Counter(tokens)
        term_freqs.append(tf)
        for term in tf:
            df[term] += 1
    avgdl = sum(doc_lens) / N if N else 0
    return {"N": N, "df": df, "avgdl": avgdl, "doc_lens": doc_lens, "term_freqs": term_freqs}


def bm25_score(query_tokens: list[str], doc_idx: int, idx: dict, k1: float = 1.5, b: float = 0.75) -> float:
    score = 0.0
    N = idx["N"]
    df = idx["df"]
    avgdl = idx["avgdl"]
    dl = idx["doc_lens"][doc_idx]
    tf = idx["term_freqs"][doc_idx]
    for term in query_tokens:
        if term not in df:
            continue
        idf = math.log(1 + (N - df[term] + 0.5) / (df[term] + 0.5))
        f = tf.get(term, 0)
        if f == 0:
            continue
        denom = f + k1 * (1 - b + b * (dl / avgdl if avgdl else 1))
        score += idf * (f * (k1 + 1)) / denom
    return score


def parse_date(s: str | None) -> date | None:
    if not s:
        return None
    try:
        return datetime.strptime(s[:10], "%Y-%m-%d").date()
    except (ValueError, TypeError):
        return None


def passes_filters(page: dict, args) -> bool:
    meta = page["meta"]
    if args.type and meta.get("type") != args.type:
        return False
    if args.tag:
        page_tags = set(meta.get("tags", []) or [])
        if not all(t in page_tags for t in args.tag):
            return False
    if args.since:
        since = parse_date(args.since)
        updated = parse_date(meta.get("updated"))
        if since and updated and updated < since:
            return False
        if since and not updated:
            return False
    return True


def primary_provenance(page: dict) -> str:
    """Return a stable key for source/concept variants derived from the same source."""
    meta = page["meta"]
    doc_path = meta.get("doc_path")
    if isinstance(doc_path, str) and doc_path:
        return doc_path.lower().replace("\\", "/")
    sources = meta.get("sources", []) or []
    if isinstance(sources, str):
        sources = [sources]
    if sources:
        return str(sources[0]).lower().replace("\\", "/")
    return f"page:{page['slug']}"


def collapse_by_provenance(scored_pages: list[tuple[float, dict]], prefer_type: str | None) -> list[tuple[float, dict]]:
    """Keep one result per primary source, optionally preferring a page type.

    Group order is determined by the group's best BM25 score. This preserves search
    relevance while allowing a synthesized concept page to represent its raw source.
    """
    groups = defaultdict(list)
    for score, page in scored_pages:
        groups[primary_provenance(page)].append((score, page))

    collapsed = []
    for candidates in groups.values():
        group_score = max(score for score, _ in candidates)
        preferred = [item for item in candidates if item[1]["meta"].get("type") == prefer_type]
        representative = max(preferred or candidates, key=lambda item: item[0])[1]
        collapsed.append((group_score, representative))
    collapsed.sort(key=lambda item: -item[0])
    return collapsed


def cmd_search(args, pages: list[dict]) -> None:
    filtered = [p for p in pages if passes_filters(p, args)]
    if not filtered:
        print("No pages matched the filters.", file=sys.stderr)
        return
    idx = build_bm25(filtered)
    query_tokens = tokenize(args.query)
    if not query_tokens:
        print("Empty query.", file=sys.stderr)
        return
    scored = [(bm25_score(query_tokens, i, idx), filtered[i]) for i in range(len(filtered))]
    scored = [(score, page) for score, page in scored if score > 0]
    scored.sort(key=lambda item: -item[0])
    if args.dedup_provenance:
        scored = collapse_by_provenance(scored, args.prefer_type)
    top = scored[:args.top]
    if not top:
        print("No matches.", file=sys.stderr)
        return
    print(f"Top {len(top)} results for: {args.query!r}")
    print()
    for score, page in top:
        title = page["meta"].get("title") or page["slug"]
        page_type = page["meta"].get("type", "?")
        print(f"  [{score:6.2f}] [{page_type:9}] {title}")
        print(f"             {page['rel_path']}")


def cmd_backlinks(args, pages: list[dict]) -> None:
    target = args.backlinks
    inbound = []
    for page in pages:
        if target in page["links"]:
            inbound.append(page)
    if not inbound:
        print(f"No pages link to [[{target}]].", file=sys.stderr)
        return
    print(f"Pages linking to [[{target}]] ({len(inbound)}):")
    for page in inbound:
        title = page["meta"].get("title") or page["slug"]
        print(f"  - {title}  ({page['rel_path']})")


def cmd_top_linked(args, pages: list[dict]) -> None:
    inbound_count = Counter()
    for page in pages:
        for link in page["links"]:
            inbound_count[link] += 1
    top = inbound_count.most_common(args.top_linked)
    if not top:
        print("No links found in the wiki.", file=sys.stderr)
        return
    print(f"Top {len(top)} most-linked-to pages (hubs):")
    for slug, count in top:
        # Try to find the page for the title
        match = next((p for p in pages if p["slug"] == slug), None)
        title = (match["meta"].get("title") if match else None) or slug
        marker = "" if match else "  [BROKEN LINK]"
        print(f"  {count:4d}  {title}  ({slug}){marker}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("query", nargs="?", default="", help="Query terms.")
    parser.add_argument("--wiki", type=Path, default=Path("wiki"), help="Wiki directory (default: ./wiki).")
    parser.add_argument("--top", type=int, default=10, help="Top N results (default: 10).")
    parser.add_argument("--type", help="Filter by frontmatter type.")
    parser.add_argument("--dedup-provenance", action="store_true",
                        help="Collapse results that share the same primary source provenance.")
    parser.add_argument("--prefer-type", help="Prefer this page type within a provenance group.")
    parser.add_argument("--tag", action="append", default=[], help="Filter by tag (repeatable).")
    parser.add_argument("--since", help="Only pages updated on or after YYYY-MM-DD.")
    parser.add_argument("--backlinks", help="Find pages linking to this slug.")
    parser.add_argument("--top-linked", type=int, help="Show the N most-linked-to pages.")
    parser.add_argument("--cache", type=Path, help="(reserved) Cache path for BM25 index.")
    args = parser.parse_args()

    if not args.wiki.exists():
        print(f"Wiki directory not found: {args.wiki}", file=sys.stderr)
        sys.exit(1)

    pages = collect_pages(args.wiki)
    if not pages:
        print(f"No wiki pages found under {args.wiki}", file=sys.stderr)
        sys.exit(0)

    if args.backlinks:
        cmd_backlinks(args, pages)
    elif args.top_linked:
        cmd_top_linked(args, pages)
    elif args.query:
        cmd_search(args, pages)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
