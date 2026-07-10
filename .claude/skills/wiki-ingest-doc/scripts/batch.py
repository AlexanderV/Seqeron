#!/usr/bin/env python3
"""
batch.py — one-time backlog driver for the `wiki-ingest-doc` skill.

Ties a bundled mapper to the checklist bookkeeping for the two big one-time sweeps, so the
mechanical parts happen in ONE command instead of by hand. It NEVER runs git — it prints the
exact commit command for a human/orchestrator to run (commits stay in human hands).

  reports : report_verdict.py over ALL reports -> register CLEAN rows -> mark those docs done.
            DEFECT reports are left pending (each still needs its own ingest subagent).
  mcp     : mcp_map.py over ALL tools -> mark the CONFIRMED (already-wired) tools done.
            proposed / ambiguous / unmatched are left pending for the catalog subagent
            (this driver deliberately does NOT write the risky lexical 'proposed' matches).

stdlib only. Dry-run by default; --apply performs the writes + checklist marking.

Usage:
  python .claude/skills/wiki-ingest-doc/scripts/batch.py reports --wiki wiki --checklist WIKI_INGEST_CHECKLIST.md
  python .claude/skills/wiki-ingest-doc/scripts/batch.py reports --wiki wiki --checklist WIKI_INGEST_CHECKLIST.md --apply
  python .claude/skills/wiki-ingest-doc/scripts/batch.py mcp     --wiki wiki --checklist WIKI_INGEST_CHECKLIST.md --apply
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mcp_map          # noqa: E402
import report_verdict   # noqa: E402

MARK_DONE = os.path.join("tools", "wiki-ingest", "mark_done.py")


def _mark_done(paths):
    ok = 0
    if not os.path.exists(MARK_DONE):
        sys.stderr.write(f"[batch] WARN: {MARK_DONE} not found — skipping checklist marking.\n")
        return 0
    for p in paths:
        r = subprocess.run([sys.executable, MARK_DONE, p], capture_output=True, text=True)
        ok += (r.returncode == 0)
    return ok


def _commit_hint(msg, extra_paths=""):
    print("\n# next (run yourself — batch.py never commits):")
    print(f"git add wiki WIKI_INGEST_CHECKLIST.md {extra_paths}".rstrip())
    print(f'git commit -m "{msg}\n\nCo-Authored-By: Claude Fable 5 <noreply@anthropic.com>"')


def do_reports(args):
    paths = sorted(os.path.join(report_verdict.REPORTS_ROOT, f)
                   for f in os.listdir(report_verdict.REPORTS_ROOT) if f.endswith(".md"))
    concepts = mcp_map.load_concepts(args.wiki)
    prov = mcp_map.build_provenance(concepts, args.wiki)
    clean, defect = [], []
    for p in paths:
        v = report_verdict.classify(p)
        if v.get("method_id"):
            pseudo = {"tool": v["unit"].lower(), "method_id": v["method_id"],
                      "method_short": v["method_id"].split(".")[-1]}
            cp, st, _ = mcp_map.match_concept(pseudo, concepts, prov)
            v["concept"] = os.path.splitext(os.path.basename(cp))[0] if cp and st in ("present", "matched", "proposed") else ""
        (defect if v["defect"] else clean).append(v)
    # a clean report that already has a full wiki/sources/<unit>-report.md page (done under the
    # old policy) is already reflected — don't also add a registry row for it.
    def _has_page(v):
        return os.path.exists(os.path.join(args.wiki, "sources", v["unit"].lower() + "-report.md"))
    to_register = [v for v in clean if not _has_page(v)]
    print(f"[batch reports] {len(paths)} reports | clean {len(clean)} "
          f"({len(to_register)} to register, {len(clean) - len(to_register)} already have pages) | "
          f"defect {len(defect)} (need subagent each)")
    if args.apply:
        added, present = report_verdict.register_clean(args.wiki, to_register)
        marked = _mark_done([v["path"] for v in to_register])
        print(f"  registered {added} new rows ({present} already there); marked {marked} CLEAN reports done")
        _commit_hint("docs(wiki): batch-register CLEAN validation verdicts")
    else:
        print("  (dry-run — pass --apply to register + mark done)")
    if defect:
        print("  DEFECT reports still pending (ingest each with the skill):")
        for v in defect:
            print(f"    - {v['path']}")


def do_mcp(args):
    paths = sorted(p.replace("\\", "/")
                   for root, _, files in os.walk(os.path.join("docs", "mcp", "tools"))
                   for p in [os.path.join(root, f) for f in files if f.endswith(".md")])
    concepts = mcp_map.load_concepts(args.wiki)
    confirmed, tail = [], []
    for tp in paths:
        t = mcp_map.parse_tool_doc(tp)
        _, st, _ = mcp_map.match_concept(t, concepts)
        (confirmed if st == "present" else tail).append(tp)
    print(f"[batch mcp] {len(paths)} tools | confirmed {len(confirmed)} | tail {len(tail)} (subagent + review)")
    if args.apply:
        marked = _mark_done(confirmed)
        print(f"  marked {marked} CONFIRMED tools done (already wired into a concept)")
        print("  tail left pending — run the catalog subagent + `mcp_map.py --all` review for the rest")
        _commit_hint("docs(wiki): batch-mark confirmed MCP tools done")
    else:
        print("  (dry-run — pass --apply to mark confirmed done)")


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass
    ap = argparse.ArgumentParser(description="One-time backlog driver for wiki-ingest-doc.")
    ap.add_argument("mode", choices=["reports", "mcp"])
    ap.add_argument("--wiki", default="wiki")
    ap.add_argument("--checklist", default="WIKI_INGEST_CHECKLIST.md")
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()
    if not mcp_map.load_concepts(args.wiki):
        sys.stderr.write(f"[batch] ERROR: no concepts under {args.wiki}/concepts/ — run from the repo root.\n")
        sys.exit(2)
    (do_reports if args.mode == "reports" else do_mcp)(args)


if __name__ == "__main__":
    main()
