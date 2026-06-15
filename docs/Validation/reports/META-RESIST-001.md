# Validation Report: META-RESIST-001 — Antibiotic Resistance Gene Detection (ResFinder-style)

- **Validated:** 2026-06-15   **Area:** Metagenomics
- **Canonical method(s):** `MetagenomicsAnalyzer.FindAntibioticResistanceGenes(contigs, referenceGenes, identityThreshold, coverageThreshold)`; internal `BestUngappedMatch(contig, reference)`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS (after in-session fix)
- **End-state:** ✅ CLEAN

## Stage A — Description

### Sources opened this session (URLs + extracted numbers)

1. **Heng Li (2018), "On the definition of sequence identity"** — https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
   - BLAST identity = "the number of matching bases over the number of alignment columns."
   - Denominator = sum of M+I+D CIGAR columns; gapless ⇒ no gap columns ⇒ denominator = aligned window length.
   - Worked example: 43 matches / 50 columns = **86%**. Confirms the `matches / window` formula used by the code.
2. **Zankari et al. (2012), JAC 67(11):2640–2644 (ResFinder)** — https://academic.oup.com/jac/article/67/11/2640/707208
   - "All genes from the ResFinder database were BLASTed against the assembled genome, and the best-matching genes were given as output." (best-matching gene, single output)
   - "For a gene to be reported, it has to cover at least **2/5** of the length of the resistance gene in the database." (coverage floor relative to **reference** length)
   - "The default ID is **100%**." (web-service default is user-selectable; default 100%)
   - %ID = "the percentage of nucleotides that are identical between the best-matching resistance gene in the database and the corresponding sequence in the genome."
3. **ResFinder GitHub (genomicepidemiology/resfinder)** — https://github.com/genomicepidemiology/resfinder
   - `CGE_RESFINDER_GENE_ID` (-t) default **0.80**; `CGE_RESFINDER_GENE_COV` (-l) default **0.60**.
   - Coverage = "breadth-of coverage … proportion of a reference gene's sequence covered (0–1)."
4. **Web search corroboration (ResFinder operating point)** — multiple secondary sources (Sci Rep 2023 s41598-023-42154-6; JAC 2016 71(9):2484): "ResFinder employs … parameters set at **90% identity and 60% coverage**"; "ResFinder 4.0 uses … ≥80% identity over ≥60% of the length." 60% floor exists "so that genes lying on the edge of a contig or spread over two contigs are not missed."
5. **CARD RGI** — https://card.mcmaster.ca/analyze/rgi — "Perfect" match = 100% identical over the entire reference length; best hit by bit-score. Corroborates INV-04 and single best-hit output.

### Formula check
- **Percent identity = matches / window (gapless)** — exact match to Li (2018). ✔
- **Coverage = window / reference length** — exact match to Zankari (2012) / ResFinder README breadth-of-coverage. ✔
- **Dual-threshold reporting (identity ≥ id AND coverage ≥ cov)** — matches Zankari reporting rule. ✔
- **Best-matching gene per contig** — matches Zankari "best-matching genes" + RGI best-hit. ✔

### Edge-case semantics
Empty contig skipped, empty reference ignored, null → `ArgumentNullException`, threshold ∉ [0,1] → `ArgumentOutOfRangeException`. All defined and standard. ✔

### Findings / divergences (Stage A)
- **N-A1 (NOTE):** Default identity threshold **0.90**. The 2012 paper's *default* is 100%; the GitHub CLI (ResFinder 4.0) default is 0.80; 0.90/0.60 is the documented historical **web-service operating point** (corroborated by source 4). The spec cites this explicitly. Acceptable as a sourced, version-specific choice — recorded as a note, not a defect. The doc/comment should make the version distinction explicit (already cites Zankari 2017 web service).
- **N-A2 (DEFECT, internal inconsistency → fixed):** The doc §4.2 documented the tie-break as "ties → longer window," but the Evidence table's own M3 expected value (`%ID = 4/4 = 1.0, coverage 4/7 ≈ 0.571`) is only producible with a **shorter/higher-identity** tie-break. The documented "longer window" rule is biologically wrong (see Stage B). Description corrected this session.

Stage A is otherwise faithful to the cited primary literature. Verdict: **PASS-WITH-NOTES** (N-A1 threshold-version note; N-A2 description fixed).

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs:1084-1191`
(`FindAntibioticResistanceGenes` + `BestUngappedMatch`). Only one caller of `BestUngappedMatch`; legacy `FindResistanceGenes` (motif stub, line 1196) is a separate MCP method, out of scope.

### Defect found and fixed — tie-break direction in `BestUngappedMatch`
The original objective was "maximize match count; **on ties prefer the longer window**." A longer
window can only be reached by padding the alignment with **mismatching** flanks, which *lowers*
identity. This is not how BLAST reports an HSP (BLAST reports the best-scoring local segment, Li 2018).

Hand-computed counter-examples (independent Python re-implementation of the sliding scan):
- **M3** `contig=TTTCGTA`, `ref=CGTACGT`: true alignment is the perfect 4-base HSP `CGTA`
  (offset +3, identity 4/4 = **1.0**, coverage 4/7 ≈ 0.571 — the Evidence-table value). The old code
  instead picked offset −1 (window 6, 4 matches → identity **0.667**, coverage 6/7 ≈ 0.857),
  contradicting the unit's own Evidence table. (Both happen to be rejected, but for the wrong reason.)
- **Outcome-changing probe** `contig=GGGACGTACG`, `ref=ACGTACGTAC` (m=10): true alignment is the
  perfect 7-base suffix HSP `ACGTACG` (identity **1.0**, coverage 7/10 = 0.70 → PASSES 0.90/0.60).
  The old code picked the padded 9-wide window (7 matches → identity **0.778**, coverage 0.90),
  which **fails** the 0.90 identity threshold → the gene is **wrongly missed (false negative)**.

**Fix:** tie-break now prefers the **shorter (higher-identity)** window
(`matches == bestMatches && bestWindow != 0 && window < bestWindow`), with `bestWindow == 0` marking
the unset state. This makes the chosen alignment the highest-identity ungapped placement and never
pads with mismatches.

### Cross-verification table recomputed vs the fixed code (independent Python model)

| Case | contig / ref | identity (sourced) | coverage (sourced) | code (fixed) |
|------|--------------|--------------------|--------------------|--------------|
| M1 exact | AAACGTACGT / CGTACGT | 7/7 = 1.0 | 7/7 = 1.0 | 1.0 / 1.0 ✔ |
| M2 1-mismatch | CGTTCGT / CGTACGT | 6/7 ≈ 0.857143 | 1.0 | 0.857143 / 1.0 ✔ |
| M3 edge-4 | TTTCGTA / CGTACGT | 4/4 = 1.0 | 4/7 ≈ 0.571429 | 1.0 / 0.571429 ✔ |
| M4 partial-5 | TTCGTAC / CGTACGT | 5/5 = 1.0 | 5/7 ≈ 0.714286 | 1.0 / 0.714286 ✔ |
| M5 low-id | CGAACTT / CGTACGT | 5/7 ≈ 0.714286 | 1.0 | 0.714286 / 1.0 ✔ |
| M6 best-of-2 | CGTACGT / {CGTACGT, CGTTCGT} | 1.0 vs 6/7 | — | geneA (1.0) wins ✔ |
| C1 tie→cov | CGTACGT / {CGTACGTGG, CGTACGT} | both 1.0 | 7/9 vs 7/7 | geneFull (cov 1.0) wins ✔ |
| probe edge-perfect | GGGACGTACG / ACGTACGTAC | 7/7 = 1.0 | 7/10 = 0.70 | 1.0 / 0.70 ✔ (reported) |

All values trace to Li (2018) identity formula and Zankari (2012) coverage definition — not to code output.

### Variant/delegate consistency
Single public method + one internal helper; defaults (0.90/0.60) are public constants asserted in M7. No other variants.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** identity/coverage assertions use exact arithmetic (`6.0/7.0`, `5.0/7.0`, `4.0/7.0`, `7.0/10.0`) from the BLAST-identity / coverage definitions, not numbers read off the implementation. ✔
- **Strengthened M3:** previously asserted only `Is.Empty` (would pass against a wrong implementation and its comment misattributed the rejection reason). Now also exposes the best-match identity (1.0) and coverage (4/7) via a lowered coverage threshold, locking the sourced values. ✔
- **Added M3b regression guard:** new test `FindAntibioticResistanceGenes_EdgePerfectHsp_PreferredOverPaddedWindow` locks the outcome-changing probe (perfect 7/10 edge HSP must be REPORTED at default thresholds, not diluted to 0.778). This test FAILS against the old "longer-window" code and PASSES against the fix — a genuine guard, not green-washing. ✔
- **Coverage of branches:** M1/M2/M4 (identity+coverage paths), M3/M3b/M5 (reject-on-coverage, reject-on-identity, edge-HSP), M6 (best-of-two), C1 (identity tie → coverage), M7 (default constants), S1–S5 (null×2, threshold-range×2, empty/no-match). All Stage-A MUST/SHOULD/COULD rows covered. ✔
- **No weakening:** no assertion loosened, no tolerance widened, no skip/ignore. Exact `Within(1e-10)`/`1e-12` tolerances retained. ✔
- **Honest green:** FULL unfiltered suite = **Failed: 0, Passed: 6556** (1 pre-existing benchmark skip); `dotnet build` 0 errors, no new warnings. ✔

### Findings / defects (Stage B)
- **F-RESIST-001 (fixed):** `BestUngappedMatch` tie-break preferred the longer (mismatch-padded) window, producing lower-identity alignments and a real false-negative under default thresholds. Fixed to prefer the shorter/higher-identity window; doc §4.2 corrected; M3 strengthened; M3b guard added.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (threshold-version note N-A1; description tie-break inconsistency N-A2 fixed).
- **Stage B:** PASS (defect F-RESIST-001 completely fixed in-session; tests locked to sourced values).
- **End-state:** ✅ CLEAN — code corrected, tests added/strengthened to sourced values, full suite green.
- No remaining follow-ups. (The gapless-vs-gapped-BLAST simplification remains a documented, accepted scope limitation, not a defect.)
