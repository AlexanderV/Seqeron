---
type: source
title: "Evidence: PROTMOTIF-FIND-001 (Protein motif search — general pattern-based PROSITE→regex engine)"
tags: [validation, protein, motif]
doc_path: docs/Evidence/PROTMOTIF-FIND-001-Evidence.md
sources:
  - docs/Evidence/PROTMOTIF-FIND-001-Evidence.md
source_commit: fb7daac68f3ea1c2203fdafb193818710566a1ee
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PROTMOTIF-FIND-001

The validation-evidence artifact for test unit **PROTMOTIF-FIND-001** — **Protein Motif
Search (Pattern-based)**: the general `ProteinMotifFinder.FindMotifByPattern` engine that
takes an **arbitrary caller-supplied PROSITE pattern** (converted to regex by
`ConvertPrositeToRegex`) and reports every hit, plus the information-content scoring
primitives `CalculateMotifScore` / `CalculateEValue`. It is one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern; the engine model,
PROSITE→regex table, overlapping-match mechanism, scoring formulae, corrected patterns and
invariants are synthesized in [[protein-motif-pattern-search]]. The fixed-dictionary
application over the same catalog is [[common-protein-motifs]] (`FindCommonMotifs`). See
[[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources (authority rank 2 — official DB / spec):**
  - **PROSITE database + User Manual** — pattern syntax (`x` = any; `[ABC]` = alternatives;
    `{ABC}` = exclusions; `(n)` / `(n,m)` = repetition; `<` / `>` = N/C-terminus; `-` =
    separator; period ends pattern) and the element-by-element **PROSITE→regex** mapping
    (`x`→`.`, `[ABC]`→`[ABC]`, `{ABC}`→`[^ABC]`, `x(n)`→`.{n}`, `x(n,m)`→`.{n,m}`, `<`→`^`,
    `>`→`$`).
  - **Official patterns verified against implementation:** PS00001, PS00004–PS00009,
    PS00016–PS00018, PS00028 (C2H2 zinc finger `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H`),
    PS00029 (leucine zipper `L-x(6)-L-x(6)-L-x(6)-L`). SKIP-FLAG=TRUE noted for the common
    false-positive-prone sites.
  - Wikipedia *Sequence motif* / *PROSITE* (rank 4) for definitions and history (Bairoch
    1988; Hulo 2007; ScanProsite De Castro 2006).
- **Documented corner cases:** overlapping matches are expected for patterns (standard regex
  misses them → solved via lookahead); case-insensitive matching; empty/null → no results;
  invalid caller regex handled gracefully.
- **Overlapping-match mechanism:** `FindMotifByPattern` uses regex lookahead `(?=(pattern))`
  to find all matches including overlaps (ScanProsite-consistent).
- **Two corrected implementation bugs (Dataset 3 & 4):**
  - **PS00007** official `[RK]-x(2)-[DE]-x(3)-Y` → `[RK].{2}[DE].{3}Y`; impl had the loosened
    `[RK].{2,3}[DE].{2,3}Y` (wrong repeat counts).
  - **PS00018 (EF-hand)** — impl used `x` instead of `{W}` at position 2 **and** dropped the
    trailing `[LIVMFYW]` element; both fixed.
- **Non-PROSITE patterns verified vs primary literature:** NLS1 `K-[KR]-x-[KR]`
  (Dingwall & Laskey 1991), NES1 `[LIVFM]-x(2,3)-[LIVFM]-x(2,3)-[LIVFM]-x-[LIVFM]`
  (la Cour 2004), SIM1 `[VIL]-x-[VIL]-[VIL]` (Hecker 2006), WW1 `P-P-x-Y` (Chen & Sudol 1995),
  SH3_1 `[RK]-x(2)-P-x(2)-P` (Mayer 2001, upgraded from a minimal PxxP core to the full Class
  I consensus).
- **Scoring (Schneider & Stephens 1990):** `CalculateMotifScore` IC = Σ log₂(20/allowed_count)
  over constrained positions (fixed pos ≈ 4.32 bits, `x` = 0 bits); `CalculateEValue`
  E = (N−L+1)·2^(−IC). An earlier BLAST (Altschul 1990) citation was **removed** — the
  E-value is a direct combinatorial probability under IC, not BLAST extreme-value statistics.
- **Test datasets:** Human Insulin P01308 N-glycosylation scan; constructed PROSITE-pattern
  verifications (PS00001 `ANFTB`→`NFTB`, `ANPSB`→no match; PS00005 `ASARK`→`SAR`; PS00006
  `ATAAD`→`TAAD`; PS00016 `ARGDA`→`RGD`; PS00017 `AXXXXGKS`→pos 0); PS00007/PS00018 correction
  fixtures.

## Deviations and assumptions

No open algorithm deviations — the change history records that **all** assumptions were
eliminated (patterns corrected, non-PROSITE patterns literature-verified, heuristic scoring
replaced by information-content scoring, overlapping-match discovery implemented). The only
standing convention is 0-based `MotifMatch.Start`/`End` vs ScanProsite's 1-based coordinates
(no correctness effect; shared with sibling ProteinMotif units).

## Recommended coverage

MUST: `FindMotifByPattern` returns correct matches and 0-based Start/End for `R-G-D`;
`ConvertPrositeToRegex` converts every PROSITE syntax element; empty/null → empty;
case-insensitive; invalid regex handled gracefully; PS00007 and PS00018 match the official
PROSITE definitions after the fixes; `FindCommonMotifs` detects PS00001 (with Proline
exclusion), PS00005, PS00006, PS00017; the `CommonMotifs` dictionary carries every official
entry with the correct pattern. SHOULD: `FindCommonMotifs` on a real protein (Human Insulin
P01308). COULD: multiple/overlapping match positions for repeated motifs.
