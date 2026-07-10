#!/usr/bin/env python3
"""
find_concept.py — "does a concept already exist for this?" guard for the `wiki-ingest-doc` skill.

Before a subagent creates a NEW concept page, run this to avoid forking a duplicate — the
single biggest quality failure of an incremental wiki. Give it a method id, a slug, or free
text (a title / algorithm name); it reuses the sibling `mcp_map.py` token matcher and prints
the existing concept(s) that already cover it.

stdlib only. Read-only (never writes). Exit 0 if a unique concept exists (ENRICH it),
2 if none (create is OK), 1 if several candidates (a human/subagent must pick).

Usage:
  python find_concept.py --wiki wiki --method SequenceAligner.GlobalAlign
  python find_concept.py --wiki wiki --name "Codon Adaptation Index"
  python find_concept.py --wiki wiki --slug global-alignment-needleman-wunsch
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mcp_map  # noqa: E402  (sibling script — same scripts/ dir)


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass
    ap = argparse.ArgumentParser(description="Find the existing wiki concept for a method/name/slug.")
    ap.add_argument("--wiki", default="wiki")
    ap.add_argument("--method", help="a Class.Method id, e.g. SequenceAligner.GlobalAlign")
    ap.add_argument("--name", help="free text: an algorithm name / title")
    ap.add_argument("--slug", help="a candidate concept slug")
    args = ap.parse_args()
    if not (args.method or args.name or args.slug):
        ap.error("give at least one of --method / --name / --slug")

    concepts = mcp_map.load_concepts(args.wiki)
    if not concepts:
        sys.stderr.write(f"[find_concept] ERROR: no concepts under {args.wiki}/concepts/ — "
                         f"run from the repo root and pass the right --wiki.\n")
        sys.exit(3)

    # Exact slug hit is decisive.
    if args.slug:
        exact = [c for c in concepts if c["slug"] == args.slug.strip().lower()]
        if exact:
            print(json.dumps({"verdict": "exists", "concepts": [exact[0]["slug"]]}, ensure_ascii=False))
            sys.stderr.write(f"[find_concept] EXISTS (exact slug): {exact[0]['slug']} — ENRICH, don't fork.\n")
            sys.exit(0)

    pseudo = {
        "tool": (args.name or args.slug or (args.method or "").split(".")[-1] or "").lower(),
        "method_id": args.method,
        "method_short": args.method.split(".")[-1] if args.method else None,
    }
    prov = mcp_map.build_provenance(concepts, args.wiki)
    cp, status, cands = mcp_map.match_concept(pseudo, concepts, prov)
    slugs = [os.path.splitext(os.path.basename(c))[0] for c in cands]

    if status in ("present", "matched", "proposed"):
        slug = os.path.splitext(os.path.basename(cp))[0]
        print(json.dumps({"verdict": "exists", "concepts": [slug]}, ensure_ascii=False))
        sys.stderr.write(f"[find_concept] EXISTS: {slug} — ENRICH it, do not create a duplicate.\n")
        sys.exit(0)
    if status == "ambiguous":
        print(json.dumps({"verdict": "ambiguous", "concepts": slugs}, ensure_ascii=False))
        sys.stderr.write(f"[find_concept] AMBIGUOUS: {slugs} — a human/subagent must pick.\n")
        sys.exit(1)
    print(json.dumps({"verdict": "none", "concepts": []}, ensure_ascii=False))
    sys.stderr.write("[find_concept] NONE — no existing concept; creating one is OK.\n")
    sys.exit(2)


if __name__ == "__main__":
    main()
