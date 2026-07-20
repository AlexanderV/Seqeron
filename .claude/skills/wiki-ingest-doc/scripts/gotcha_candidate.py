#!/usr/bin/env python3
"""
gotcha_candidate.py — deterministic sharp-edge extractor for the `wiki-ingest-doc` skill.

Step 5 of an ingest ("GOTCHA — do NOT skip") used to be a manual prose scan, which is why
`wiki/gotchas/` stayed nearly empty. This script does the mechanical 90%: for one algorithm
doc it (a) resolves the wiki concept + its `mcp_tools:` bindings, (b) extracts the doc's
trap-signalling sections (Deviations / Assumptions / Scope / Limitations / Simplified /
Corner cases) and trap phrases, (c) checks whether a gotcha already covers those tools, and
(d) emits a verdict so the subagent only has to write the prose or say "no gotcha".

The recurring, highest-value trap in this repo is **"this is a heuristic / propensity profile /
Simplified subset, NOT the trained / calibrated / full published method the name implies"** —
that pattern drove ~20 of the first gotchas, and this extractor is tuned to surface it.

stdlib only. Read-only (never writes). UTF-8-safe. Fail-loud (exit 3) if run outside the repo
root. Reuses the sibling `mcp_map.py` for concept loading + the algo-doc->concept provenance
bridges, so the tool bindings it reports match what `mcp_map.py --apply` writes.

Usage:
  python gotcha_candidate.py --wiki wiki --doc docs/algorithms/ProteinPred/Disorder_Prediction.md
  python gotcha_candidate.py --wiki wiki --doc <doc> --json
  python gotcha_candidate.py --wiki wiki --doc <doc> --check   # exit 1 if signal present but no gotcha yet

Exit codes: 0 = report emitted (or --check clean), 1 = --check: uncovered gotcha signal,
2 = no usable input, 3 = run outside repo root.
"""
import argparse
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mcp_map  # noqa: E402  (sibling script — same scripts/ dir)

# STRONG headers almost always carry a real sharp edge — a non-empty section under one of
# these sets the verdict on its own. WEAK headers (edge/corner-case tables, notes) are mostly
# benign API contracts (empty->0, null->throws); they only count when the body also hits a
# trap phrase. Both are still surfaced in the report for the agent to judge.
STRONG_HEADERS = re.compile(
    r"deviation|assumption|scope|limitation|simplif|sharp edge|"
    r"not implemented|known (issue|limitation)|out of scope",
    re.I,
)
WEAK_HEADERS = re.compile(r"corner case|caveat|edge case|\bnotes?\b", re.I)
TRAP_HEADERS = re.compile(STRONG_HEADERS.pattern + "|" + WEAK_HEADERS.pattern, re.I)

# Inline phrases that flag a trap even outside a trap section. Tuned from the first 24 gotchas:
# the dominant class is "heuristic / profile / Simplified, not the trained/calibrated/full thing".
TRAP_PHRASES = re.compile(
    r"\bheuristic\b|\bsimplified\b|propensity|presence/absence|read-tally|"
    r"not (a substitute|the full|the exact|calibrat|trained|a trained|for clinical|symmetric|left)|"
    r"does not (call|apply|use|correct|resolve|dedup|rank|filter|normali|assemble|guarantee|perform|model)|"
    r"no built-in|no (training|dispersion|correction|overlap resolution|genotype)|"
    r"forward-only|only the (first|forward|query|most specific)|"
    r"caller'?s job|callers? must not|do not (read|interpret|treat)|"
    r"silently|underestimat|overestimat|without a reference|returns only|"
    r"assumes a|assumes the|assumes that|clamp|clipped|bounded to|"
    r"only .* (bundled|shipped)|not the (published|canonical|standard)|"
    r"per-fragment|not.*bit-score|status \*?simplified",
    re.I,
)

STATUS_RE = re.compile(r"(?:Implementation Status|Status)\s*\|\s*([A-Za-z][\w /*-]*)", re.I)
SIMPLIFIED_RE = re.compile(r"status\s*\*?\s*simplified", re.I)


def _fm_list(fm_text, key):
    """Return the YAML block-list values under `key:` in a frontmatter string."""
    m = re.search(r"^" + re.escape(key) + r":\s*$([\s\S]*?)(?=^\S|\Z)", fm_text or "", re.M)
    if not m:
        return []
    return [ln.strip()[2:].strip()
            for ln in m.group(1).splitlines() if ln.strip().startswith("- ")]


def load_gotcha_tool_index(wiki_dir):
    """tool name -> gotcha slug, from every wiki/gotchas page's mcp_tools frontmatter."""
    index = {}
    for p in glob.glob(os.path.join(wiki_dir, "gotchas", "*.md")):
        slug = os.path.splitext(os.path.basename(p))[0]
        fm, _ = mcp_map.split_frontmatter(open(p, encoding="utf-8").read())
        for t in _fm_list(fm, "mcp_tools"):
            index.setdefault(t, slug)
    return index


def resolve_concept(doc_path, concept_slug, concepts, wiki_dir):
    """Find the concept for an algorithm doc: explicit slug, else provenance bridges."""
    if concept_slug:
        for c in concepts:
            if c["slug"] == concept_slug.strip().lower():
                return c
    doc_norm = doc_path.replace("\\", "/")
    prov = mcp_map.build_provenance(concepts, wiki_dir)
    slugs = set(prov["algo_to_concepts"].get(doc_norm) or set())
    if doc_norm in prov["backlog_map"]:
        slugs = slugs | {prov["backlog_map"][doc_norm]}
    slugs = {s for s in slugs if s not in mcp_map.HUB_SLUGS}
    if len(slugs) == 1:
        want = next(iter(slugs))
        for c in concepts:
            if c["slug"] == want:
                return c
    return None


def extract_sections(doc_text):
    """Return [(header, trimmed_body)] for headings whose title matches TRAP_HEADERS."""
    lines = doc_text.splitlines()
    heads = [(i, ln) for i, ln in enumerate(lines) if ln.startswith("#")]
    out = []
    for idx, (i, ln) in enumerate(heads):
        title = ln.lstrip("#").strip()
        if not TRAP_HEADERS.search(title):
            continue
        end = heads[idx + 1][0] if idx + 1 < len(heads) else len(lines)
        body = [l.rstrip() for l in lines[i + 1:end] if l.strip()]
        if not body:
            continue
        # Suppress a trivially-empty trap section ("Deviations: None") — a clean
        # true-negative that would otherwise fire the verdict.
        joined = " ".join(body).strip().lower()
        joined = re.sub(r"[*_`>\-|:.\s]+", " ", joined).strip()
        if re.fullmatch(r"(none|n a|not applicable|no deviations?( from the sources?)?( are recorded)?)", joined):
            continue
        strong = bool(STRONG_HEADERS.search(title))
        # A weak-header section only counts toward the verdict if its body hits a trap phrase.
        counts = strong or bool(TRAP_PHRASES.search(" ".join(body)))
        out.append((title, body[:12], strong, counts))
    return out


def extract_phrases(doc_text, section_spans_text):
    """Trap phrase lines that are NOT already inside a captured trap section."""
    hits, seen = [], set()
    for ln in doc_text.splitlines():
        s = ln.strip()
        if len(s) < 40 or s in section_spans_text:
            continue
        if TRAP_PHRASES.search(s):
            key = s[:80].lower()
            if key not in seen:
                seen.add(key)
                hits.append(s)
    return hits


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass
    ap = argparse.ArgumentParser(description="Extract sharp-edge (gotcha) candidates from an algorithm doc.")
    ap.add_argument("--wiki", default="wiki")
    ap.add_argument("--doc", help="algorithm/spec doc path (docs/algorithms/**)")
    ap.add_argument("--concept", help="concept slug, if the doc->concept bridge can't resolve it")
    ap.add_argument("--json", action="store_true", help="machine-readable output")
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if the doc shows gotcha signal but no gotcha covers its tools (CI gate)")
    args = ap.parse_args()

    if not args.doc:
        ap.error("give --doc <algorithm doc path>")
    if not os.path.isfile(args.doc):
        sys.stderr.write(f"[gotcha_candidate] ERROR: no such file: {args.doc}\n")
        sys.exit(2)

    concepts = mcp_map.load_concepts(args.wiki)
    if not concepts:
        sys.stderr.write(f"[gotcha_candidate] ERROR: no concepts under {args.wiki}/concepts/ — "
                         f"run from the repo root and pass the right --wiki.\n")
        sys.exit(3)

    doc_text = open(args.doc, encoding="utf-8").read()
    concept = resolve_concept(args.doc, args.concept, concepts, args.wiki)
    tools = _fm_list(concept["fm"], "mcp_tools") if concept else []

    status = ""
    ms = STATUS_RE.search(doc_text)
    if ms:
        status = ms.group(1).strip()
    simplified = bool(SIMPLIFIED_RE.search(doc_text)) or "simplif" in status.lower()

    sections = extract_sections(doc_text)
    section_text = {l for _, body, _, _ in sections for l in body}
    phrases = extract_phrases(doc_text, section_text)

    gotcha_index = load_gotcha_tool_index(args.wiki)
    covering = sorted({gotcha_index[t] for t in tools if t in gotcha_index})

    counting_sections = [s for s in sections if s[3]]
    has_signal = simplified or bool(counting_sections) or len(phrases) >= 2
    if covering:
        verdict = "COVERED"
    elif has_signal:
        verdict = "LIKELY_GOTCHA"
    else:
        verdict = "NO_SIGNAL"

    mt = re.search(r'^#\s+(.*)$', doc_text, re.M)
    doc_title = mt.group(1).strip() if mt else os.path.basename(args.doc)

    result = {
        "doc": args.doc.replace("\\", "/"),
        "doc_title": doc_title,
        "status": status,
        "simplified": simplified,
        "concept": concept["slug"] if concept else None,
        "mcp_tools": tools,
        "trap_sections": [{"header": h, "lines": body, "strong": strong, "counts": counts}
                          for h, body, strong, counts in sections],
        "trap_phrases": phrases[:8],
        "existing_gotcha": covering,
        "verdict": verdict,
    }

    if args.json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        print("=" * 72)
        print(f"GOTCHA CANDIDATE  —  {result['doc']}")
        print("=" * 72)
        print(f"title:      {doc_title}")
        print(f"status:     {status or '(none)'}{'   [SIMPLIFIED]' if simplified else ''}")
        print(f"concept:    {result['concept'] or '(unresolved — pass --concept)'}")
        print(f"mcp_tools:  {', '.join(tools) or '(none — run mcp_map.py first)'}")
        print(f"covered by: {', '.join(covering) if covering else '(no existing gotcha)'}")
        print(f"VERDICT:    {verdict}")
        if sections:
            print("\n--- trap sections (*=sets verdict; others surfaced for your judgment) ---")
            for h, body, strong, counts in sections:
                mark = "*" if counts else " "
                print(f"{mark}## {h}")
                for l in body:
                    print(f"   {l[:150]}")
        if phrases:
            print("\n--- trap phrases (outside those sections) ---")
            for p in phrases[:8]:
                print(f"   - {p[:150]}")
        if verdict == "LIKELY_GOTCHA":
            print("\nNEXT: write wiki/gotchas/<slug>.md — frontmatter mcp_tools: "
                  f"{tools or '[<tool>]'}, sources: [{result['doc']}], source_commit: $(git rev-parse HEAD);")
            print("      body = **The trap.** / **Why it bites.** / **What to rely on instead.**;")
            print(f"      then backlink from [[{result['concept'] or '<concept>'}]], add to "
                  "wiki/indexes/gotchas.md, bump wiki/index.md count.")
        elif verdict == "NO_SIGNAL":
            print("\nNEXT: no deterministic sharp-edge signal — record \"no gotcha\" in your report "
                  "(still eyeball the trigger checklist in SKILL.md).")

    if args.check:
        sys.exit(1 if verdict == "LIKELY_GOTCHA" else 0)
    sys.exit(0)


if __name__ == "__main__":
    main()
