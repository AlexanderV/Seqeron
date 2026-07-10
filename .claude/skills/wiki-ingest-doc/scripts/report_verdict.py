#!/usr/bin/env python3
"""
report_verdict.py — classify a validation report (CLEAN vs DEFECT) and register the
cheap ones, for the `wiki-ingest-doc` skill. Bundled with the skill.

Recommendation it implements: a CLEAN per-unit report is NOT worth a full subagent or a
standalone wiki page — it is one verdict row in a shared registry. Only a report that
FOUND A DEFECT needs the full subagent (which corrects the concept + writes a page +
adds a gotcha). This script does the cheap half: parse the report, decide clean/defect,
and (with --apply) append a one-line verdict row to the wiki registry
`wiki/sources/validation-verdicts.md` for the CLEAN ones — idempotently, no subagent.

The main agent then: for `clean` -> mark_done + commit (done, no subagent);
for `defect` -> spawn the full ingest subagent.

stdlib only. Dry-run by default; --apply writes the registry rows for CLEAN reports.

Usage:
  python report_verdict.py --wiki wiki --reports docs/Validation/reports/ALIGN-SEMI-001.md
  python report_verdict.py --wiki wiki --all --apply
  python report_verdict.py --wiki wiki --pending WIKI_INGEST_CHECKLIST.md --apply
"""
import argparse
import glob
import json
import os
import re
import sys

REPORTS_ROOT = os.path.join("docs", "Validation", "reports")
REGISTRY = os.path.join("sources", "validation-verdicts.md")  # under --wiki

# A report needs the full subagent ONLY when a CODE defect was found — because then the
# concept may describe the old (buggy) behavior and must be corrected. A pure test-coverage
# / test-quality / code-echo fix does NOT touch the concept, so it is a CLEAN registry row.
# Detection looks ONLY at the Stage-B verdict field (not the whole body, which mentions
# "not a defect", "failure mode", etc. in prose).
TEST_ONLY_RE = re.compile(
    r"test[- ]quality|test[- ]coverage|coverage gap|code[- ]echo|green[- ]wash|"
    r"test[- ]data|grammar[- ]branch|branch coverage|strengthened|weak assertion", re.I)
CODE_DEFECT_RE = re.compile(
    r"code defect|off[- ]by[- ]one|rounding|probability defect|fidelity defect|"
    r"coordinate|code gap|`[A-Za-z]", re.I)  # backtick+symbol = a named code method


def _stage_b_is_code_defect(stage_b):
    sb = stage_b or ""
    if re.search(r"\bFAIL\b", sb):
        return True
    if re.search(r"no code (?:defect|change)", sb, re.I):
        return False              # report explicitly states no code defect -> CLEAN
    if re.search(r"defect|fix", sb, re.I):
        if TEST_ONLY_RE.search(sb) and not CODE_DEFECT_RE.search(sb):
            return False          # a test-only fix -> concept unaffected -> CLEAN
        if CODE_DEFECT_RE.search(sb):
            return True
        return True               # a defect/fix mentioned, not clearly test-only -> be safe
    return False


def _field(text, label):
    m = re.search(r"^-\s*\*\*" + re.escape(label) + r":\*\*\s*(.+?)\s*$", text, re.M)
    return m.group(1).strip() if m else None


def classify(path):
    text = open(path, encoding="utf-8").read()
    m = re.search(r"^#\s*Validation Report:\s*([A-Z0-9][A-Z0-9-]+)", text, re.M)
    unit = m.group(1) if m else os.path.splitext(os.path.basename(path))[0]
    stage_a = _field(text, "Stage A verdict") or ""
    stage_b = _field(text, "Stage B verdict") or ""
    end_state = _field(text, "End-state") or _field(text, "End state") or ""
    validated = _field(text, "Validated") or ""
    area = _field(text, "Area") or ""
    # area may share the "Validated" line: "2026-06-24   **Area:** Alignment"
    if not area:
        ma = re.search(r"\*\*Area:\*\*\s*([A-Za-z /_-]+)", text)
        area = ma.group(1).strip() if ma else ""
    validated = re.split(r"\s{2,}|\*\*", validated)[0].strip()

    hay = " ".join([stage_a, stage_b, end_state])
    defect = _stage_b_is_code_defect(stage_b)
    verdict = "DEFECT" if defect else ("CLEAN" if re.search(r"CLEAN|PASS", hay, re.I) else "UNKNOWN")
    tests = ""
    for pat in (r"(?:Passed|passed):\s*([\d,]+)",         # "Passed: 6535"
                r"suite\s*(?:of\s*)?([1-9][\d,]{1,})",    # "suite 6532", "suite of 17"
                r"\b([1-9]\d{1,6})\s*/\s*0\b"):           # "6535/0", "17/0" (no leading zero)
        mt = re.search(pat, text)
        if mt:
            tests = mt.group(1).replace(",", "")
            break
    # first canonical `Class.Method` on the Canonical-method line — for the concept link.
    # Skip test files (`.cs`) and assertion helpers (Is./Does./Assert./That.).
    method_id = None
    mline = re.search(r"Canonical method\(s\):\**\s*(.+)$", text, re.M)
    if mline:
        for cand in re.findall(r"`([A-Za-z_][\w]*\.[A-Za-z_][\w]*)", mline.group(1)):
            if cand.endswith(".cs") or re.match(r"(Is|Does|Assert|That|Has|Contains)\.", cand):
                continue
            method_id = cand
            break
    return {
        "unit": unit, "verdict": verdict, "defect": defect,
        "stage_a": stage_a, "stage_b": stage_b, "validated": validated,
        "area": area, "tests_passed": tests, "method_id": method_id,
        "path": path.replace("\\", "/"),
    }


def _short(v):
    a = "PASS" if re.search(r"pass", v["stage_a"], re.I) else (v["stage_a"] or "?")
    b = re.sub(r"\s*\(.*\)$", "", v["stage_b"]) or "?"
    return f"{a} / {b}"


REGISTRY_HEADER = (
    "---\n"
    "type: source\n"
    'title: "Validation verdict registry (CLEAN per-unit reports)"\n'
    "tags: [validation, governance, registry]\n"
    "sources:\n"
    "  - docs/Validation/reports/\n"
    "---\n\n"
    "# Validation verdict registry\n\n"
    "One row per **CLEAN** per-unit validation report (no defect, no page needed). "
    "Reports that found a defect get their own `wiki/sources/<unit>-report.md` and a "
    "correction to the concept instead. See [[validation-ledger]] / [[validation-protocol]] "
    "for the governance model and [[validation-and-testing]] for the campaign.\n\n"
    "| Unit | Concept | State | Stage A/B | Validated | Tests |\n"
    "|------|---------|-------|-----------|-----------|-------|\n"
)


def register_clean(wiki_dir, cleans):
    """Idempotently append CLEAN rows to the registry page. Returns (added, present)."""
    reg = os.path.join(wiki_dir, REGISTRY)
    os.makedirs(os.path.dirname(reg), exist_ok=True)
    text = open(reg, encoding="utf-8").read() if os.path.exists(reg) else REGISTRY_HEADER
    added, present = 0, 0
    new_rows = []
    for v in cleans:
        if re.search(r"^\|\s*" + re.escape(v["unit"]) + r"\s*\|", text, re.M):
            present += 1
            continue
        concept = v.get("concept") or ""
        concept_cell = f"[[{concept}]]" if concept else "?"
        new_rows.append(f"| {v['unit']} | {concept_cell} | {v['verdict']} | {_short(v)} | "
                        f"{v['validated']} | {v['tests_passed']} |")
        added += 1
    if new_rows:
        if not text.endswith("\n"):
            text += "\n"
        text += "\n".join(new_rows) + "\n"
        with open(reg, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
    return added, present


def discover(args):
    if args.reports:
        return [r.replace("\\", "/") for r in args.reports]
    if args.pending:
        rows = open(args.pending, encoding="utf-8").read().splitlines()
        pat = re.compile(r"^- \[ \] (docs/Validation/reports/.+\.md)\s*$")
        return [m.group(1) for ln in rows if (m := pat.match(ln))]
    if args.all:
        return sorted(p.replace("\\", "/") for p in glob.glob(os.path.join(REPORTS_ROOT, "*.md")))
    return []


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass
    ap = argparse.ArgumentParser(description="Classify validation reports; register CLEAN verdicts.")
    ap.add_argument("--wiki", default="wiki")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--reports", nargs="+", help="specific report path(s)")
    src.add_argument("--all", action="store_true", help="all docs/Validation/reports/*.md")
    src.add_argument("--pending", metavar="CHECKLIST", help="only reports still `- [ ]`")
    ap.add_argument("--apply", action="store_true", help="write CLEAN rows to the registry")
    ap.add_argument("--check", action="store_true",
                    help="CI gate: exit 1 if any given report is not yet reflected (a CLEAN report "
                         "not in the registry, or a DEFECT report with no wiki/sources/<unit>-report.md). "
                         "Never writes.")
    args = ap.parse_args()
    if args.check:
        args.apply = False

    paths = discover(args)
    if not paths:
        sys.stderr.write("[report_verdict] ERROR: no report docs selected (bad path / empty --pending / "
                         "no docs/Validation/reports). Run from the repo root.\n")
        sys.exit(2)

    # Resolve each report's concept via the sibling mcp_map matcher (best-effort; "" if unsure).
    try:
        import mcp_map
        _concepts = mcp_map.load_concepts(args.wiki)
        _prov = mcp_map.build_provenance(_concepts, args.wiki) if _concepts else None
    except Exception:
        _concepts, _prov = [], None

    def _resolve_concept(v):
        if not _concepts or not v.get("method_id"):
            return ""
        pseudo = {"tool": v["unit"].lower(), "method_id": v["method_id"],
                  "method_short": v["method_id"].split(".")[-1]}
        cp, status, _ = mcp_map.match_concept(pseudo, _concepts, _prov)
        return os.path.splitext(os.path.basename(cp))[0] if cp and status in ("present", "matched", "proposed") else ""

    report = {"clean": [], "defect": [], "unknown": [], "apply": bool(args.apply)}
    for p in paths:
        if not os.path.exists(p):
            report["unknown"].append({"path": p, "reason": "file not found"})
            continue
        v = classify(p)
        v["concept"] = _resolve_concept(v)
        bucket = "defect" if v["defect"] else ("clean" if v["verdict"] == "CLEAN" else "unknown")
        report[bucket].append(v)

    if args.check:
        reg_text = ""
        reg_path = os.path.join(args.wiki, REGISTRY)
        if os.path.exists(reg_path):
            reg_text = open(reg_path, encoding="utf-8").read()
        missing = []
        for v in report["clean"]:
            if not re.search(r"^\|\s*" + re.escape(v["unit"]) + r"\s*\|", reg_text, re.M):
                missing.append(v["unit"] + " (clean, not in registry)")
        for v in report["defect"]:
            page = os.path.join(args.wiki, "sources", v["unit"].lower() + "-report.md")
            if not os.path.exists(page):
                missing.append(v["unit"] + " (defect, no report page)")
        json.dump({"unreflected": missing}, sys.stdout, indent=2, ensure_ascii=False)
        sys.stdout.write("\n")
        sys.stderr.write(f"[report_verdict] CHECK: {len(missing)} unreflected of {len(paths)}\n")
        sys.exit(1 if missing else 0)

    if args.apply and report["clean"]:
        added, present = register_clean(args.wiki, report["clean"])
        report["registry"] = {"added": added, "already_present": present,
                              "path": os.path.join(args.wiki, REGISTRY).replace("\\", "/")}

    json.dump(report, sys.stdout, indent=2, ensure_ascii=False)
    sys.stdout.write("\n")
    reg = report.get("registry", {})
    sys.stderr.write(
        f"[report_verdict] {len(paths)} reports | clean {len(report['clean'])} "
        f"(registered {reg.get('added', 0)}) | DEFECT {len(report['defect'])} (need subagent) | "
        f"unknown {len(report['unknown'])} | {'APPLIED' if args.apply else 'dry-run'}\n")
    if report["defect"]:
        sys.stderr.write("  defect reports (spawn the full ingest subagent for each):\n")
        for v in report["defect"]:
            sys.stderr.write(f"    - {v['unit']}  ({v['stage_b']})\n")


if __name__ == "__main__":
    main()
