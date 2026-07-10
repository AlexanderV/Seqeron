"""Batch-map MCP tool docs to existing wiki concepts (deterministic, idempotent).

The MCP tool layer (`docs/mcp/tools/**/*.md`, ~427 docs) is a thin wrapper surface:
each tool binds one library method (its **Method ID**), and that method's algorithm
almost always already has a concept page in `wiki/concepts/`. Ingesting each tool as
its own wiki page would produce hundreds of duplicates, so instead we:

  1. Parse every *pending* tool doc (unchecked `- [ ]` rows in WIKI_INGEST_CHECKLIST.md)
     for its server, tool name, Method ID, and wrapped source file.
  2. Match each tool's Method ID to an existing concept via three deterministic bridges:
       (A) algorithm doc chain: the Method ID appears in a `docs/algorithms/**` doc that
           the backlog's "Covered via concept" table maps to a concept;
       (B) concept `sources:`: the Method ID appears in an algorithm doc that some concept
           lists in its frontmatter `sources:`;
       (C) literal: the Method ID string appears verbatim in a concept body.
     Ambiguous multi-concept hits are resolved by a conservative slug/method-name token
     overlap tie-break, else left unmatched.
  3. Idempotently append the tool name to the matched concept's `mcp_tools:` frontmatter
     list (line-based edit; preserves YAML, never duplicates, safe to re-run).
  4. Emit `<scratchpad>/mcp_catalog_data.json` (server -> [{tool, concept, wrapped_method}])
     plus the unmatched list, for the catalog-page author.

Only wiki concept frontmatter and the scratchpad JSON are written. `docs/**` is never
touched. Re-running is safe: matching is pure, and frontmatter edits are set-unions.

Usage:
    python tools/wiki-ingest/mcp_map.py [--scratch <dir>] [--dry-run]
"""
from __future__ import annotations

import json
import os
import re
import sys
import pathlib
from collections import defaultdict

ROOT = pathlib.Path(".").resolve()
CHK = pathlib.Path("WIKI_INGEST_CHECKLIST.md")
CONCEPTS_DIR = pathlib.Path("wiki/concepts")
BACKLOG = pathlib.Path("wiki/backlog.md")
ALGO_ROOT = pathlib.Path("docs/algorithms")

DEFAULT_SCRATCH = pathlib.Path(
    r"C:\Users\O990A~1.PAN\AppData\Local\Temp\claude"
    r"\d--Prototype-Bio-Seqeron\2248de47-95be-41d9-8773-de12abbfd8e2\scratchpad"
)

# A concept token is "distinctive" if it occurs in <= DF_MAX concept slug+title sets.
# Generic domain words (variant, distance, composition, kmer, gene, tree...) have a
# higher document frequency and must not, alone, justify a token-overlap match.
DF_MAX = 4

# Residual semantic mismatches the token heuristic accepts but that are wrong on
# domain grounds (right token, wrong concept — usually because no better concept
# exists yet). These fall back to the catalog's gap section instead of writing
# frontmatter. Keyed by (concept_slug, tool_name).
REJECT: set[tuple[str, str]] = {
    ("microsatellite-instability-detection", "find_microsatellites"),
    ("antibiotic-resistance-gene-detection", "find_best_match"),
    ("consensus-from-alignment", "most_frequent_kmers"),
    ("conserved-gene-clusters-common-intervals", "calculate_conservation"),
    ("contig-merge-overlap-collapse", "merge_overlapping_svs"),
    ("contig-merge-overlap-collapse", "bed_merge"),
    ("evolutionary-distance-matrix", "hamming_distance"),
    ("expression-quantification", "quantile_normalize"),
    ("genetic-code-translation", "amino_acid_composition"),
    ("k-mer-statistics", "shannon_entropy"),
    ("k-mer-statistics", "entropy_profile"),
    ("longest-repeated-substring", "longest_dinucleotide_repeat"),
    ("phylogenetic-marker-selection", "cluster_genes_by_expression"),
    ("somatic-variant-calling-tumor-normal", "normalize_variant"),
    ("homozygous-deletion-detection", "find_deletions"),
}

# Generic verbs/nouns that carry no discriminating signal for the tie-break.
STOPWORDS = {
    "find", "calculate", "compute", "detect", "get", "analyze", "analyse", "predict",
    "build", "create", "run", "estimate", "call", "count", "score", "extract", "of",
    "and", "the", "a", "from", "to", "with", "sequence", "sequences", "seq", "dna",
    "rna", "region", "regions", "site", "sites", "analysis", "analyzer", "tool",
    "tools", "data", "result", "results", "per", "all", "value", "values",
}


def read(p: pathlib.Path) -> str:
    return p.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# 1. Pending tool docs
# ---------------------------------------------------------------------------

def pending_tool_docs() -> list[str]:
    text = read(CHK)
    out = []
    for line in text.splitlines():
        m = re.match(r"- \[ \] (docs/mcp/tools/\S+\.md)\s*$", line.strip())
        if m:
            out.append(m.group(1))
    return out


def parse_tool_doc(rel: str) -> dict | None:
    p = pathlib.Path(rel)
    if not p.exists():
        return None
    text = read(p)
    server_dir = rel.split("/")[3]  # docs/mcp/tools/<server>/<name>.md

    def field(name: str) -> str | None:
        m = re.search(rf"\*\*{re.escape(name)}\*\*\s*\|\s*`?([^|`]+?)`?\s*\|", text)
        return m.group(1).strip() if m else None

    tool = field("Tool Name") or p.stem
    method_id = field("Method ID")
    server = field("Server") or server_dir

    # wrapped source .cs reference (display text preferred, else link target)
    wrapped_source = None
    ms = re.search(r"(?:Source|Implementation location):\s*\[([^\]]+)\]\(([^)]+)\)", text)
    if ms:
        wrapped_source = ms.group(1).strip()
    if not wrapped_source:
        ms2 = re.search(r"\[([^\]]*\.cs[^\]]*)\]", text)
        if ms2:
            wrapped_source = ms2.group(1).strip()

    return {
        "path": rel,
        "server_dir": server_dir,
        "server": server,
        "tool": tool,
        "method_id": method_id,
        "wrapped_source": wrapped_source,
    }


# ---------------------------------------------------------------------------
# 2. Concept + backlog + algorithm-doc indexes
# ---------------------------------------------------------------------------

def concept_frontmatter_sources(text: str) -> list[str]:
    """Extract the algorithm-doc paths listed under frontmatter `sources:`."""
    lines = text.split("\n")
    if not lines or lines[0].strip() != "---":
        return []
    end = next((i for i in range(1, len(lines)) if lines[i].strip() == "---"), None)
    if end is None:
        return []
    fm = lines[1:end]
    out, in_sources = [], False
    for l in fm:
        if re.match(r"^sources:\s*$", l):
            in_sources = True
            continue
        if in_sources:
            m = re.match(r"^\s+-\s+(\S+)", l)
            if m:
                out.append(m.group(1).replace("\\", "/"))
            elif re.match(r"^\S", l):  # next top-level key ends the block
                in_sources = False
    return out


def frontmatter_title(text: str) -> str:
    m = re.search(r'^title:\s*"?([^"\n]+?)"?\s*$', text, flags=re.M)
    return m.group(1).strip() if m else ""


def build_indexes(tools: list[dict]):
    # concept text + slug + title
    concept_text: dict[str, str] = {}
    concept_sources: dict[str, list[str]] = {}
    concept_title: dict[str, str] = {}
    for cp in sorted(CONCEPTS_DIR.glob("*.md")):
        slug = cp.stem
        t = read(cp)
        concept_text[slug] = t
        concept_sources[slug] = concept_frontmatter_sources(t)
        concept_title[slug] = frontmatter_title(t)

    # invert: algo doc -> concepts that cite it in sources
    algo_to_concepts: dict[str, set[str]] = defaultdict(set)
    for slug, srcs in concept_sources.items():
        for s in srcs:
            if s.startswith("docs/algorithms/"):
                algo_to_concepts[s].add(slug)

    # backlog covered table: algo doc -> concept slug
    backlog_map: dict[str, str] = {}
    for line in read(BACKLOG).splitlines():
        m = re.match(r"^\|\s*`(docs/algorithms/[^`]+)`\s*\|\s*\[\[([^\]|]+)(?:\|[^\]]*)?\]\]\s*\|", line)
        if m:
            backlog_map[m.group(1)] = m.group(2).strip()

    # algorithm doc texts
    algo_text: dict[str, str] = {}
    for ap in ALGO_ROOT.rglob("*.md"):
        algo_text[ap.as_posix()] = read(ap)

    # precomputed stemmed token sets per concept + document frequency
    concept_tok: dict[str, set[str]] = {}
    for slug, title in concept_title.items():
        concept_tok[slug] = {stem(t) for t in concept_tokens(slug, title)}
    df: dict[str, int] = defaultdict(int)
    for toks in concept_tok.values():
        for t in toks:
            df[t] += 1

    return concept_text, concept_title, concept_tok, df, algo_to_concepts, backlog_map, algo_text


# ---------------------------------------------------------------------------
# 3. Matching
# ---------------------------------------------------------------------------

def method_variants(method_id: str) -> list[str]:
    """`Class.A/B/C` -> ['Class.A/B/C', 'Class.A', 'Class.B', 'Class.C']."""
    if not method_id:
        return []
    variants = {method_id}
    if "." in method_id:
        cls, rest = method_id.split(".", 1)
        for part in re.split(r"[/,]", rest):
            part = part.strip()
            if part:
                variants.add(f"{cls}.{part}")
    return list(variants)


def camel_tokens(name: str) -> set[str]:
    name = name.split(".")[-1]
    parts = re.findall(r"[A-Z]+(?=[A-Z][a-z])|[A-Z]?[a-z]+|[A-Z]+|\d+", name)
    return {p.lower() for p in parts if p.lower() not in STOPWORDS and len(p) > 1}


def slug_tokens(slug: str) -> set[str]:
    return {p for p in slug.split("-") if p not in STOPWORDS and len(p) > 1}


def text_tokens(s: str) -> set[str]:
    return {p for p in re.split(r"[^a-z0-9]+", s.lower()) if p not in STOPWORDS and len(p) > 1}


def concept_tokens(slug: str, title: str) -> set[str]:
    return slug_tokens(slug) | text_tokens(title)


def tool_tokens(tool: dict) -> set[str]:
    return camel_tokens(tool["method_id"]) | text_tokens(tool["tool"])


_SUFFIXES = ("ization", "isation", "ational", "ation", "ition", "ings", "ing",
             "tions", "tion", "ements", "ement", "ments", "ment", "ness",
             "ers", "er", "ed", "es", "s")


def stem(w: str) -> str:
    for suf in _SUFFIXES:
        if w.endswith(suf) and len(w) - len(suf) >= 3:
            return w[: -len(suf)]
    return w


def teq(a: str, b: str) -> bool:
    """Token equivalence: same stem, or one stem a >=4-char prefix of the other.

    Handles stats/statistics, parse/parsing, detect/detection, align/alignment;
    short tokens (bed, orf, cai, fst) must match exactly."""
    sa, sb = stem(a), stem(b)
    if sa == sb:
        return True
    lo, hi = sorted((sa, sb), key=len)
    return len(lo) >= 4 and hi.startswith(lo)


def fuzzy_overlap(ttoks: set[str], ctoks: set[str]) -> int:
    return sum(1 for a in ttoks if any(teq(a, b) for b in ctoks))


def token_match(tool, concept_tok, df):
    """Bridge D: pick the concept whose slug+title tokens best cover the tool's tokens.

    Precision guards (all required):
      * unique winner by matched-token count (a strict margin over the runner-up);
      * at least one matched token is *distinctive* (concept-side df <= DF_MAX) —
        so a shared generic word (variant, distance, composition) never matches alone;
      * either the winner covers ALL of the tool's tokens (cov == 1.0) or it shares
        >= 2 distinctive tokens with the tool.
    """
    ttok = tool_tokens(tool)
    if not ttok:
        return None, []
    stems = {stem(t) for t in ttok}
    scored = []
    for slug, ctoks in concept_tok.items():
        matched = {a for a in stems if any(teq(a, b) for b in ctoks)}
        if not matched:
            continue
        distinctive = {a for a in matched
                       if any(teq(a, b) and df.get(b, 99) <= DF_MAX for b in ctoks)}
        scored.append((len(matched), len(distinctive), slug))
    if not scored:
        return None, []
    scored.sort(reverse=True)
    best_n, _, best_slug = scored[0]
    top5 = [s for _, _, s in scored[:5]]
    winners = [s for n, _, s in scored if n == best_n]
    if len(winners) != 1:
        return None, top5
    n_matched, n_dist, slug = scored[0]
    cov = n_matched / len(stems)
    if n_dist >= 1 and (cov == 1.0 or n_dist >= 2):
        return slug, top5
    return None, top5


def match_tool(tool, concept_text, concept_title, concept_tok, df, algo_to_concepts, backlog_map, algo_text):
    method_id = tool["method_id"]
    variants = method_variants(method_id)
    if not variants:
        return None, "no-method-id", []

    # algorithm docs whose text contains any method variant
    hit_algo = [ap for ap, txt in algo_text.items() if any(v in txt for v in variants)]

    candidates: dict[str, set[str]] = defaultdict(set)  # concept -> via-tags
    for ap in hit_algo:
        if ap in backlog_map:
            candidates[backlog_map[ap]].add("A:backlog-algo-chain")
        for slug in algo_to_concepts.get(ap, ()):
            candidates[slug].add("B:concept-sources-algo")

    # (C) literal Method ID in concept body
    for slug, txt in concept_text.items():
        if any(v in txt for v in variants):
            candidates[slug].add("C:literal-in-concept")

    def finalize(slug, via, cands):
        if slug and (slug, tool["tool"]) in REJECT:
            return None, f"rejected({via})", cands
        return slug, via, cands

    if len(candidates) == 1:
        slug = next(iter(candidates))
        return finalize(slug, "+".join(sorted(candidates[slug])), list(candidates))

    if len(candidates) > 1:
        # tie-break among the exact-bridge candidates by token overlap
        ttok = tool_tokens(tool)
        scored = sorted(
            ((fuzzy_overlap(ttok, concept_tok.get(s, set())), s) for s in candidates),
            reverse=True,
        )
        top = scored[0][0]
        winners = [s for o, s in scored if o == top]
        if top > 0 and len(winners) == 1:
            return finalize(winners[0], "+".join(sorted(candidates[winners[0]])) + "+tiebreak", list(candidates))
        # fall through to token match to try to break the tie

    # Bridge D: scored slug+title token overlap (covers the common case where the
    # method lives in an algo doc no concept cites, but a concept is named for it).
    slug, top5 = token_match(tool, concept_tok, df)
    if slug:
        return finalize(slug, "D:token-overlap", top5)
    return None, "no-candidate" if not candidates else "ambiguous", (top5 or list(candidates))


# ---------------------------------------------------------------------------
# 4. Idempotent frontmatter edit
# ---------------------------------------------------------------------------

def add_mcp_tools(path: pathlib.Path, tool_names: set[str], dry_run: bool) -> tuple[bool, list[str]]:
    text = read(path)
    lines = text.split("\n")
    if not lines or lines[0].strip() != "---":
        return False, []
    end = next((i for i in range(1, len(lines)) if lines[i].strip() == "---"), None)
    if end is None:
        return False, []
    fm = lines[1:end]

    key_idx = next((i for i, l in enumerate(fm) if re.match(r"^mcp_tools:\s*$", l)), None)
    existing: set[str] = set()
    block_end = None
    if key_idx is not None:
        j = key_idx + 1
        while j < len(fm) and re.match(r"^\s+-\s+", fm[j]):
            existing.add(fm[j].strip()[1:].strip())
            j += 1
        block_end = j

    all_tools = sorted(existing | tool_names)
    new_block = ["mcp_tools:"] + [f"  - {t}" for t in all_tools]

    if key_idx is not None:
        fm = fm[:key_idx] + new_block + fm[block_end:]
    else:
        insert_at = next((i + 1 for i, l in enumerate(fm) if re.match(r"^tags:", l)), len(fm))
        fm = fm[:insert_at] + new_block + fm[insert_at:]

    new_text = "\n".join([lines[0]] + fm + lines[end:])
    if new_text != text and not dry_run:
        path.write_text(new_text, encoding="utf-8")
    return (new_text != text), all_tools


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    dry_run = "--dry-run" in sys.argv
    scratch = DEFAULT_SCRATCH
    if "--scratch" in sys.argv:
        scratch = pathlib.Path(sys.argv[sys.argv.index("--scratch") + 1])
    scratch.mkdir(parents=True, exist_ok=True)

    pend = pending_tool_docs()
    tools = [t for t in (parse_tool_doc(r) for r in pend) if t]
    (concept_text, concept_title, concept_tok, df,
     algo_to_concepts, backlog_map, algo_text) = build_indexes(tools)

    catalog: dict[str, list[dict]] = defaultdict(list)
    unmatched: list[dict] = []
    concept_hits: dict[str, set[str]] = defaultdict(set)  # concept -> tool names

    for t in tools:
        slug, via, cands = match_tool(t, concept_text, concept_title, concept_tok, df,
                                      algo_to_concepts, backlog_map, algo_text)
        rec = {
            "tool": t["tool"],
            "concept": slug,
            "wrapped_method": t["method_id"],
            "wrapped_source": t["wrapped_source"],
            "doc": t["path"],
            "via": via,
        }
        catalog[t["server"]].append(rec)
        if slug:
            concept_hits[slug].add(t["tool"])
        else:
            unmatched.append({**rec, "candidates": cands})

    # write mcp_tools frontmatter
    touched, changed = [], 0
    for slug, names in sorted(concept_hits.items()):
        cp = CONCEPTS_DIR / f"{slug}.md"
        if not cp.exists():
            continue
        did_change, _ = add_mcp_tools(cp, names, dry_run)
        touched.append(slug)
        if did_change:
            changed += 1

    matched_n = sum(len(concept_hits[s]) for s in concept_hits)
    out = {
        "generated_from": "tools/wiki-ingest/mcp_map.py",
        "pending_total": len(pend),
        "parsed_total": len(tools),
        "matched_total": len(tools) - len(unmatched),
        "unmatched_total": len(unmatched),
        "concepts_touched": sorted(concept_hits),
        "catalog": {k: sorted(v, key=lambda r: r["tool"]) for k, v in sorted(catalog.items())},
        "unmatched": sorted(unmatched, key=lambda r: (r["server"] if "server" in r else "", r["tool"]))
        if False else sorted(unmatched, key=lambda r: r["tool"]),
    }
    (scratch / "mcp_catalog_data.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # summary
    print("=" * 64)
    print("MCP tool -> concept batch mapping")
    print("=" * 64)
    print(f"pending tool docs : {len(pend)}")
    print(f"parsed            : {len(tools)}")
    print(f"matched           : {len(tools) - len(unmatched)}")
    print(f"unmatched         : {len(unmatched)}")
    print(f"concepts touched  : {len(concept_hits)} (frontmatter changed: {changed}{' (dry-run)' if dry_run else ''})")
    print()
    print("Per-server matched / total:")
    for srv in sorted(catalog):
        rows = catalog[srv]
        m = sum(1 for r in rows if r["concept"])
        print(f"  {srv:16s} {m:3d} / {len(rows):3d}")
    print()
    if unmatched:
        print(f"Unmatched ({len(unmatched)}) — bucketed as gaps in the catalog:")
        for r in sorted(unmatched, key=lambda r: r["tool"])[:60]:
            print(f"  {r['tool']:38s} {r['wrapped_method'] or '(no method id)'}  [{r['via']}]")
        if len(unmatched) > 60:
            print(f"  ... and {len(unmatched) - 60} more (see mcp_catalog_data.json)")
    print()
    print(f"JSON -> {scratch / 'mcp_catalog_data.json'}")


if __name__ == "__main__":
    main()
