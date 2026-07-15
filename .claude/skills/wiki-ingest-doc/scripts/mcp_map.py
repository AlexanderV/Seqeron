#!/usr/bin/env python3
"""
mcp_map.py — deterministic MCP-tool -> concept mapper for the `wiki-ingest-doc` skill.

Bundled with the skill so they travel together. Does the mechanical 90% of the
MCP-layer ingest: for each `docs/mcp/tools/**/*.md` tool doc, find the existing
`wiki/concepts/*.md` that documents the wrapped method, and idempotently append the
tool name to that concept's `mcp_tools:` YAML frontmatter list. Emits a JSON report
of matched / unmatched / ambiguous tools so a subagent only has to resolve the tail
and write the catalog page — no per-tool wiki pages, no per-tool subagent.

stdlib only. Dry-run by default; pass --apply to write.

Usage:
  # one tool doc (the per-doc skill path):
  python mcp_map.py --wiki wiki --tools docs/mcp/tools/alignment/global_align.md --apply

  # whole layer (the one-time batch path):
  python mcp_map.py --wiki wiki --all --apply

  # only tools still pending in the checklist:
  python mcp_map.py --wiki wiki --pending WIKI_INGEST_CHECKLIST.md --apply

Exit code 0 always (report is the product); check the JSON `unmatched`/`ambiguous`.
"""
import argparse
import glob
import json
import os
import re
import sys

TOOLS_ROOT = os.path.join("docs", "mcp", "tools")
ALGO_ROOT = os.path.join("docs", "algorithms")

# ---------- tool-doc parsing ----------

def _table_value(text, label):
    """Pull `**Label** | value` from the Overview table; strip backticks."""
    m = re.search(r"\*\*" + re.escape(label) + r"\*\*\s*\|\s*`?([^`|\n]+?)`?\s*(?:\||\n)", text)
    return m.group(1).strip() if m else None


def parse_tool_doc(path):
    text = open(path, encoding="utf-8").read()
    tool = _table_value(text, "Tool Name")
    if not tool:
        # fall back to the H1 or the filename
        m = re.search(r"^#\s+(\S+)", text, re.M)
        tool = m.group(1).strip() if m else os.path.splitext(os.path.basename(path))[0]
    method_id = _table_value(text, "Method ID")      # e.g. "SequenceAligner.GlobalAlign"
    server = _table_value(text, "Server")            # e.g. "Alignment"
    method_short = method_id.split(".")[-1] if method_id else None
    # wrapped .cs source reference (display text preferred, else the first .cs link)
    wrapped_source = None
    ms = re.search(r"(?:Source|Implementation location):\s*\[([^\]]+)\]\(", text)
    if not ms:
        ms = re.search(r"\[([^\]]*\.cs[^\]]*)\]", text)
    if ms:
        wrapped_source = ms.group(1).strip()
    return {
        "path": path.replace("\\", "/"),
        "tool": tool,
        "method_id": method_id,
        "method_short": method_short,
        "server": server,
        "wrapped_source": wrapped_source,
    }

# ---------- concept index ----------

# Generic verb/connector tokens that carry no identity — dropped from tool tokens.
STOP = {"find", "get", "compute", "calculate", "calc", "analyze", "analyse",
        "run", "make", "build", "to", "of", "and", "the", "with", "from", "a", "an"}

# Generic nouns kept as tokens but too weak to ANCHOR a match on their own
# (a tool whose tokens are all weak — e.g. find_best_match, sequence_identity — is
# left unmatched for the subagent rather than force-mapped to a spurious overlap).
WEAK = {"best", "match", "matches", "all", "read", "reads", "pair", "pairs",
        "sequence", "sequences", "seq", "identity", "region", "regions",
        "site", "sites", "score", "scores", "count", "counts", "value", "values",
        "result", "results", "stat", "stats", "data", "list", "set"}


def _camel_split(s):
    return re.sub(r"(?<=[a-z0-9])(?=[A-Z])", " ", s)


def tokenize(*parts):
    """Lowercase identity tokens: split on non-alnum and CamelCase; drop stopwords & <3-char."""
    text = " ".join(p for p in parts if p)
    text = _camel_split(text)
    toks = re.split(r"[^A-Za-z0-9]+", text.lower())
    return {t for t in toks if len(t) >= 3 and t not in STOP and not t.isdigit()}




def already_has_tool(fm_text, tool_name):
    """True if the concept frontmatter already lists this tool under mcp_tools:."""
    m = re.search(r"^mcp_tools:\s*$([\s\S]*?)(?=^\S|\Z)", fm_text, re.M)
    if not m:
        return False
    return any(ln.strip()[2:].strip() == tool_name
              for ln in m.group(1).splitlines() if ln.strip().startswith("- "))


# Meta / hub concept slugs that mention everything — never a match TARGET
# (they token-overlap almost any tool and would attract false positives).
HUB_SLUGS = {"algorithm-validation-evidence", "test-unit-registry", "mcp-tool-catalog",
             "validation-and-testing", "validation-findings-disposition",
             "operating-envelope-and-limitation-policy", "scientific-rigor"}

# Curated (concept_slug, tool_name) rejections — right token, WRONG concept (usually a
# cross-domain collision with no better concept yet). Ported from the canonical
# tools/wiki-ingest/mcp_map.py; complements the automatic class-family veto (catches the
# ones it misses, e.g. normalize_variant, hamming_distance). Rejected -> unmatched (gap).
REJECT = {
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


def _slug_of(path):
    return os.path.splitext(os.path.basename(path))[0]


def _fm_sources(fm):
    """Algorithm-doc paths listed under the concept frontmatter `sources:` block."""
    out, in_src = [], False
    for ln in fm.split("\n"):
        if re.match(r"^sources:\s*$", ln):
            in_src = True
            continue
        if in_src:
            m = re.match(r"^\s+-\s+(\S+)", ln)
            if m:
                out.append(m.group(1).replace("\\", "/"))
            elif re.match(r"^\S", ln):
                in_src = False
    return out


def method_variants(method_id):
    """`Class.A/B/C` -> {'Class.A/B/C', 'Class.A', 'Class.B', 'Class.C'} for provenance search."""
    if not method_id:
        return set()
    variants = {method_id}
    if "." in method_id:
        cls, rest = method_id.split(".", 1)
        for part in re.split(r"[/,]", rest):
            part = part.strip()
            if part:
                variants.add(f"{cls}.{part}")
    return variants


def load_concepts(wiki_dir):
    """List of dicts: path, slug, fm, body, sources (algo docs), ident-token set, is_hub."""
    out = []
    for p in glob.glob(os.path.join(wiki_dir, "concepts", "*.md")):
        slug = os.path.splitext(os.path.basename(p))[0]
        text = open(p, encoding="utf-8").read()
        fm, body = split_frontmatter(text)
        fm = fm or ""
        mt = re.search(r'^title:\s*"?(.*?)"?\s*$', fm, re.M)
        title = mt.group(1) if mt else ""
        out.append({"path": p, "slug": slug, "fm": fm, "body": body or text,
                    "sources": _fm_sources(fm),
                    "ident": tokenize(slug, title), "is_hub": slug in HUB_SLUGS})
    return out


def build_provenance(concepts, wiki_dir):
    """Indexes for the Method-ID provenance bridges (ported from tools/wiki-ingest/mcp_map.py):
    algo_to_concepts (bridge B: concept `sources:` cites the algo doc), backlog_map
    (bridge A: wiki/backlog.md 'covered' table maps algo doc -> concept), algo_text."""
    algo_to_concepts = {}
    for c in concepts:
        for s in c["sources"]:
            if s.startswith("docs/algorithms/"):
                algo_to_concepts.setdefault(s, set()).add(c["slug"])
    backlog_map = {}
    bp = os.path.join(wiki_dir, "backlog.md")
    if os.path.exists(bp):
        for ln in open(bp, encoding="utf-8"):
            m = re.match(r"^\|\s*`(docs/algorithms/[^`]+)`\s*\|\s*\[\[([^\]|]+)(?:\|[^\]]*)?\]\]\s*\|", ln)
            if m:
                backlog_map[m.group(1)] = m.group(2).strip()
    algo_text = {}
    for root, _, files in os.walk(ALGO_ROOT):
        for f in files:
            if f.endswith(".md"):
                ap = os.path.join(root, f).replace("\\", "/")
                algo_text[ap] = open(ap, encoding="utf-8").read()
    return {"algo_to_concepts": algo_to_concepts, "backlog_map": backlog_map, "algo_text": algo_text}


def provenance_match(tool, prov, by_slug):
    """Bridges A/B/C by Method ID -> (slug|None, via, candidate_slugs). High precision:
    a match here is trustworthy (the method's provenance chain reaches the concept)."""
    variants = method_variants(tool.get("method_id"))
    if not variants or not prov:
        return None, "", []
    cands = {}
    for ap, txt in prov["algo_text"].items():
        if any(v in txt for v in variants):
            if ap in prov["backlog_map"]:
                cands.setdefault(prov["backlog_map"][ap], set()).add("A:backlog")
            for slug in prov["algo_to_concepts"].get(ap, ()):
                cands.setdefault(slug, set()).add("B:sources")
    for c in by_slug.values():                       # (C) literal Method ID in concept body
        if not c["is_hub"] and any(v in c["body"] for v in variants):
            cands.setdefault(c["slug"], set()).add("C:literal")
    cands = {s: v for s, v in cands.items() if s in by_slug and not by_slug[s]["is_hub"]}
    if len(cands) == 1:
        s = next(iter(cands))
        return s, "+".join(sorted(cands[s])), [s]
    if len(cands) > 1:                               # tie-break by token overlap
        toks = tokenize(tool["tool"]) | (tokenize(tool["method_short"]) if tool["method_short"] else set())
        scored = sorted(((len(toks & by_slug[s]["ident"]), s) for s in cands), reverse=True)
        if scored[0][0] > 0 and (len(scored) == 1 or scored[1][0] < scored[0][0]):
            s = scored[0][1]
            return s, "+".join(sorted(cands[s])) + "+tiebreak", list(cands)
        return None, "ambiguous", list(cands)
    return None, "", []


def _covers(t, c):
    """A tool token t is covered by concept token c: equal, or one is a PREFIX of the
    other with the shorter >= 4 chars (align/alignment, cpg/cpgs, maxent/maxentscan).
    Prefix — not arbitrary substring — to avoid 'rich' inside 'enrichment' false hits."""
    if t == c:
        return True
    lo, hi = (t, c) if len(t) <= len(c) else (c, t)
    return len(lo) >= 4 and hi.startswith(lo)


def _token_match(strong_toks, ident):
    """Every STRONG tool token must be prefix-covered by the concept identity.
    Weak/generic tokens are ignored for coverage (they only fail to anchor, gated earlier)."""
    if not strong_toks:
        return False
    return all(any(_covers(t, c) for c in ident) for t in strong_toks)


def _covers_any(ident, toks):
    return any(any(_covers(t, c) for c in ident) for t in toks)


# Generic class-name suffixes that carry no family identity.
CLASS_SUFFIX = {"analyzer", "analyser", "finder", "calculator", "tool", "tools", "service",
                "util", "utils", "utility", "helper", "manager", "predictor", "detector",
                "designer", "caller", "parser", "generator", "annotator", "optimizer",
                "optimiser", "classifier", "engine", "statistics", "stats"}


def class_tokens(method_id):
    """Family tokens from the Method ID's CLASS, e.g. StructuralVariantAnalyzer -> {structural, variant}.
    Used to disambiguate / veto lexical hits so a `StructuralVariantAnalyzer` tool cannot map to a
    population-genetics 'genotype' concept."""
    if not method_id or "." not in method_id:
        return set()
    return {t for t in tokenize(method_id.split(".")[0]) if t not in CLASS_SUFFIX}


def match_concept(tool, concepts, prov=None, by_slug=None):
    """Return (concept_path|None, status, candidates). Precision-first tiers:
    (1) already in mcp_tools -> present; (2) Method-ID PROVENANCE bridge (backlog / sources /
    literal) -> matched (trustworthy, auto-appliable); (3) slug/title token + class-family veto
    -> proposed (lexical, needs review). Provenance is tried before tokens so precise wins."""
    for c in concepts:
        if already_has_tool(c["fm"], tool["tool"]):
            return c["path"], "present", [c["path"]]
    if prov is not None:
        if by_slug is None:
            by_slug = {c["slug"]: c for c in concepts}
        slug, via, _ = provenance_match(tool, prov, by_slug)
        if slug and (slug, tool["tool"]) not in REJECT:
            return by_slug[slug]["path"], "matched", [by_slug[slug]["path"]]
    tool_toks = tokenize(tool["tool"]) | (tokenize(tool["method_short"]) if tool["method_short"] else set())
    strong = tool_toks - WEAK
    if not strong:
        return None, "unmatched", []
    hits = [c for c in concepts if not c["is_hub"] and _token_match(strong, c["ident"])]
    paths = [c["path"] for c in hits]
    ctoks = class_tokens(tool.get("method_id"))
    if len(hits) > 1 and ctoks:
        refined = [c for c in hits if _covers_any(c["ident"], ctoks)]
        if len(refined) == 1:                       # class family breaks the tie
            return refined[0]["path"], "proposed", paths
    if len(hits) == 1:
        if ctoks and not _covers_any(hits[0]["ident"], ctoks):
            return None, "unmatched", paths          # lone hit contradicts the class family -> don't write
        if (hits[0]["slug"], tool["tool"]) in REJECT:
            return None, "unmatched", paths          # curated rejection
        return hits[0]["path"], "proposed", paths
    if len(hits) > 1:
        return None, "ambiguous", paths
    return None, "unmatched", []

# ---------- frontmatter editing (no yaml dep) ----------

def split_frontmatter(text):
    if not text.startswith("---\n"):
        return None, text
    end = text.find("\n---\n", 4)
    if end == -1:
        return None, text
    return text[4:end], text[end + 5:]


def append_mcp_tool(concept_path, tool_name):
    """Idempotently add `- tool_name` under `mcp_tools:` in the concept frontmatter.
    Returns 'added' | 'present' | 'error'."""
    text = open(concept_path, encoding="utf-8").read()
    fm, body = split_frontmatter(text)
    if fm is None:
        return "error"
    lines = fm.split("\n")
    key_i = next((i for i, ln in enumerate(lines) if re.match(r"^mcp_tools:\s*$", ln)), None)
    if key_i is not None:
        j = key_i + 1
        items = []
        while j < len(lines) and re.match(r"^\s*-\s+", lines[j]):
            items.append(lines[j].strip()[2:].strip())
            j += 1
        if tool_name in items:
            return "present"
        lines.insert(j, f"  - {tool_name}")
    else:
        # insert a new block before `sources:` (fallback: before source_commit, else end)
        anchor = next((i for i, ln in enumerate(lines) if ln.startswith("sources:")), None)
        if anchor is None:
            anchor = next((i for i, ln in enumerate(lines) if ln.startswith("source_commit:")), None)
        if anchor is None:
            anchor = len(lines)
        lines[anchor:anchor] = ["mcp_tools:", f"  - {tool_name}"]
    new = "---\n" + "\n".join(lines) + "\n---\n" + body
    with open(concept_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(new)
    return "added"

# ---------- tool discovery ----------

def discover(args):
    if args.tools:
        return [t.replace("\\", "/") for t in args.tools]
    if args.pending:
        rows = open(args.pending, encoding="utf-8").read().splitlines()
        pat = re.compile(r"^- \[ \] (docs/mcp/tools/.+\.md)\s*$")
        return [m.group(1) for ln in rows if (m := pat.match(ln))]
    if args.all:
        return sorted(p.replace("\\", "/")
                      for p in glob.glob(os.path.join(TOOLS_ROOT, "**", "*.md"), recursive=True))
    return []


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass
    ap = argparse.ArgumentParser(description="Map MCP tool docs to wiki concepts and append mcp_tools frontmatter.")
    ap.add_argument("--wiki", default="wiki", help="wiki root (default: wiki)")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--tools", nargs="+", help="specific tool doc path(s)")
    src.add_argument("--all", action="store_true", help="all docs/mcp/tools/**/*.md")
    src.add_argument("--pending", metavar="CHECKLIST", help="only tools still `- [ ]` in this checklist")
    ap.add_argument("--apply", action="store_true", help="write changes (default: dry-run)")
    ap.add_argument("--trust-proposed", action="store_true",
                    help="with --apply, also write the lexical 'proposed' matches "
                         "(review them first — cross-domain token collisions happen)")
    ap.add_argument("--check", action="store_true",
                    help="CI gate: exit 1 if any given tool is not yet CONFIRMED in a concept "
                         "(unmatched/ambiguous/proposed count as unreflected). Never writes.")
    ap.add_argument("--catalog", metavar="PATH",
                    help="also write a per-server catalog JSON (server -> [{tool, concept, "
                         "wrapped_method, wrapped_source, via}] + unmatched) — the input for "
                         "regenerating wiki/concepts/mcp-tool-catalog.md.")
    args = ap.parse_args()

    tool_paths = discover(args)
    concepts = load_concepts(args.wiki)
    if not concepts:
        sys.stderr.write(f"[mcp_map] ERROR: no concepts under {args.wiki}/concepts/ — "
                         f"run from the repo root and pass the right --wiki.\n")
        sys.exit(2)
    if not tool_paths:
        sys.stderr.write("[mcp_map] ERROR: no tool docs selected (bad path / empty --pending / "
                         "no docs/mcp/tools). Run from the repo root.\n")
        sys.exit(2)
    if args.check:
        args.apply = args.trust_proposed = False
    prov = build_provenance(concepts, args.wiki)
    by_slug = {c["slug"]: c for c in concepts}
    # confirmed = already in mcp_tools (0 risk). matched = Method-ID provenance bridge
    # (trustworthy — written on --apply). proposed = lexical token match (REVIEW; needs
    # --trust-proposed). ambiguous / unmatched are handed to the subagent.
    report = {"confirmed": [], "matched": [], "proposed": [], "ambiguous": [], "unmatched": [],
              "apply": bool(args.apply), "trust_proposed": bool(args.trust_proposed)}
    catalog = {}  # server -> [record]  (built when --catalog is requested)

    for tp in tool_paths:
        if not os.path.exists(tp):
            report["unmatched"].append({"tool_doc": tp, "reason": "file not found"})
            continue
        t = parse_tool_doc(tp)
        concept, status, cands = match_concept(t, concepts, prov, by_slug)
        rel = lambda p: os.path.relpath(p).replace("\\", "/")
        if args.catalog:
            catalog.setdefault(t.get("server") or "?", []).append({
                "tool": t["tool"], "concept": _slug_of(concept) if concept else None,
                "wrapped_method": t["method_id"], "wrapped_source": t.get("wrapped_source"),
                "doc": t["path"], "via": status})
        if status == "present":
            report["confirmed"].append({"tool": t["tool"], "method_id": t["method_id"],
                                        "concept": rel(concept)})
        elif status == "matched":                    # provenance bridge — trustworthy
            action = "would-add"
            if args.apply:
                action = append_mcp_tool(concept, t["tool"])
            report["matched"].append({"tool": t["tool"], "method_id": t["method_id"],
                                      "concept": rel(concept), "action": action})
        elif status == "proposed":                   # lexical — review before writing
            action = "proposed"
            if args.apply and args.trust_proposed:
                action = append_mcp_tool(concept, t["tool"])
            report["proposed"].append({"tool": t["tool"], "method_id": t["method_id"],
                                       "concept": rel(concept), "action": action})
        elif status == "ambiguous":
            report["ambiguous"].append({"tool": t["tool"], "method_id": t["method_id"],
                                        "candidates": [rel(c) for c in cands]})
        else:
            report["unmatched"].append({"tool": t["tool"], "method_id": t["method_id"],
                                        "server": t.get("server")})

    if args.catalog:
        recs = [r for v in catalog.values() for r in v]
        unmatched_recs = [r for r in recs if not r["concept"]]
        out = {
            "generated_from": ".claude/skills/wiki-ingest-doc/scripts/mcp_map.py",
            "parsed_total": len(recs),
            "matched_total": len(recs) - len(unmatched_recs),
            "unmatched_total": len(unmatched_recs),
            "catalog": {k: sorted(v, key=lambda r: r["tool"]) for k, v in sorted(catalog.items())},
            "unmatched": sorted(unmatched_recs, key=lambda r: r["tool"]),
        }
        with open(args.catalog, "w", encoding="utf-8", newline="\n") as f:
            json.dump(out, f, indent=2, ensure_ascii=False)
        sys.stderr.write(f"[mcp_map] catalog JSON -> {args.catalog} "
                         f"({out['matched_total']}/{len(recs)} mapped)\n")

    json.dump(report, sys.stdout, indent=2, ensure_ascii=False)
    sys.stdout.write("\n")
    m_added = sum(1 for x in report["matched"] if x["action"] == "added")
    p_added = sum(1 for x in report["proposed"] if x["action"] == "added")
    mode = "APPLIED" if args.apply else "dry-run"
    sys.stderr.write(
        f"[mcp_map] {len(tool_paths)} tools | confirmed {len(report['confirmed'])} | "
        f"matched {len(report['matched'])} (provenance, written {m_added}) | "
        f"proposed {len(report['proposed'])} (lexical, written {p_added}) | "
        f"ambiguous {len(report['ambiguous'])} | unmatched {len(report['unmatched'])} | {mode}\n")

    if args.check:
        unreflected = len(report["proposed"]) + len(report["ambiguous"]) + len(report["unmatched"])
        if unreflected:
            sys.stderr.write(f"[mcp_map] CHECK FAIL: {unreflected} tool(s) not yet provenance-wired.\n")
            sys.exit(1)


if __name__ == "__main__":
    main()
