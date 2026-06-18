# Validation Report: CODON-STATS-001 — Codon Usage Statistics

- **Validated:** 2026-06-15   **Area:** Codon
- **Canonical method(s):** `CodonUsageAnalyzer.GetStatistics(string|DnaSequence)`, `CalculateCai(string|DnaSequence, Dictionary<string,double>)`, `EColiOptimalCodons`, `HumanOptimalCodons`, `CodonUsageStatistics.OverallGc`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independent of the repo)

| Source | Retrieved | What it confirms |
|--------|-----------|------------------|
| Wikipedia, *Codon Adaptation Index* (WebFetch) | 2026-06-15 | `w_i = f_i / max(f_j)`; `CAI = (∏ w_i)^(1/L)` = geometric mean over L codons; range ~[0,1]. |
| seqinr `cai` reference doc (WebFetch) | 2026-06-15 | w = ratio of codon usage to most-abundant synonym; CAI = geometric mean via natural-log summation; **"Non-synonymous codons and termination codons (genetic code dependent) are excluded"**; zero-value codons floored to 0.01 (Bulmer 1988). |
| CodonW `Indices.html` (WebFetch) | 2026-06-15 | GC3s = "fraction of codons that are synonymous at the third codon position which have either a guanine or cytosine at that third codon position"; CAI excludes non-synonymous + termination codons. |
| Peden 1999 thesis PDF (WebFetch → `pdftotext`) | 2026-06-15 | §1.8.2.1.3 line 2065-2066 **verbatim**: "The index GC3s, is the frequency of G or C nucleotides present at the third position of synonymous codons (i.e. excluding Met, Trp and termination codons)." |
| Biopython v1.79 `SharpEcoliIndex` (WebFetch raw) | 2026-06-15 | Full E. coli w table; 61 sense-codon entries; stops absent; ATG=TGG=1.0. |
| Kazusa H. sapiens [gbpri] (WebFetch) | 2026-06-15 | Per-thousand frequencies; 93,487 CDS / 40,662,582 codons. |

### Formula check
- **CAI**: `w_i = f_i/max(f_j)`, `CAI = exp[(1/L) Σ ln w_i]` — matches Wikipedia + seqinr exactly.
- **Exclusions**: non-synonymous (single-codon Met ATG, Trp TGG) + termination codons excluded from CAI — confirmed verbatim by seqinr **and** CodonW.
- **GC3s**: "G/C at third position of synonymous codons, excluding Met/Trp/stop" — confirmed verbatim by Peden §1.8.2.1.3 (the exact string the spec quotes was located in the fetched PDF).
- **GC1/GC2/GC3**: per-codon-position G/C content (EMBOSS cusp "1st/2nd/3rd letter GC").
- **RSCU**: `RSCU_j = n·x_j / Σ_k x_k` (Sharp et al. 1986).

### Edge-case semantics
- Only-Met/Trp/stop → CAI 0, GC3s 0 (empty scorable/synonymous set) — sourced (seqinr/Peden).
- Empty/null string → zeroed stats / CAI 0; null DnaSequence or null reference → `ArgumentNullException` — input contract (not literature-specified, explicitly documented as such).
- Non-ACGT codon skipped; trailing partial codon ignored — EMBOSS cusp counts valid codons only.

### Independent cross-check (numbers recomputed this session)
- M2 `√(1×0.122)` = **0.3492849839314596** (Python). Matches test.
- M3 `∛(0.007×0.003×0.066)` = **0.011149474795453503** (Python). Matches test.
- E. coli w-table (TTT 0.296, CTG 1, CTC 0.037, CTA 0.007, GCT 1, GCC 0.122, GCA 0.586, GCG 0.424, CGT 1, AGG 0.002, ATA 0.003, GTT 1, GTC 0.066 …) — **every requested value matched the fetched Biopython source**.
- Human RSCU **recomputed from the fetched Kazusa per-thousand frequencies**: CTG = 6·39.6/100.2 = **2.37126**, CTC = **1.17365**, GCC = 4·27.7/69.3 = **1.59885**, GCT = **1.06205**, GTG = 4·28.1/60.7 = **1.85173**. All match the implementation's 2.3713 / 1.1737 / 1.5988 / 1.0620 / 1.8517 to 4 dp.

### Findings / divergences (→ PASS-WITH-NOTES)
1. **GC3s as percentage (×100)** vs CodonW's fraction in [0,1] — display-units only; the synonymous subset (numerator/denominator) is exactly per Peden. Documented ASM-1.
2. **Zero-w codons skipped** rather than floored to 0.01 (Bulmer 1988, used by seqinr/EMBOSS) — affects only a gene using a codon entirely absent from the reference; with the bundled tables no synonymous w is 0, so all worked examples are unaffected. Documented deviation.
3. **GC3s 6-fold subtlety**: CodonW's "synonymous at the third position" treats 6-fold families (Leu/Ser/Arg) per-block; the implementation (and Peden's own *definition*) uses the simpler "exclude Met/Trp/stop" set. The test cases use the 4-fold Ala family, so the divergence is not exercised; the validated *definition* (Peden §1.8.2.1.3) is realised faithfully. Noted, not a defect against the stated definition.

Both deviations are unit/edge-case choices, properly documented, and do not affect any sourced expected value. Stage A = **PASS-WITH-NOTES**.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs`:
- `CalculateCaiCore` (L142-184): builds `w = ref/family-max` skipping stops (`'*'`) and single-codon families (`Count==1`), then `exp(logSum/count)` over codons with `w>0`. Realises the validated formula and exclusion rule.
- `GetStatisticsCore` (L389-443): per-position GC counts; GC3s via `IsSynonymousAtThirdPosition` (degeneracy>1, excludes `'*'`) — matches Peden definition.
- `EColiOptimalCodons` (L195-213): 64 entries = 61 Biopython w + 3 stops at 0.0. Verified against fetched source.
- `HumanOptimalCodons` (L224-242): Kazusa-derived RSCU. Verified by recomputation.

### Cross-verification table recomputed vs code (full suite run)
| Case | Expected (sourced) | Test result |
|------|--------------------|-------------|
| M1 all-optimal | 1.0 | pass |
| M2 √(1×0.122) | 0.3492849839314596 | pass |
| M3 ∛(·) | 0.011149474795453503 | pass |
| M4 Met/Trp/stop only | 0.0 | pass |
| M5 ATGGCA GC3s=0, GC3=50 | 0/50 | pass |
| M6 GCCGCA GC3s=50 | 50 | pass |
| M7 CTGGTTAAA GC1/2/3 | 66.667/0/33.333 | pass |
| M8 counts/total | 4, CTG=2 | pass |
| M9 E.coli w | exact (Biopython) | pass |
| M10 human RSCU | exact (Kazusa) | pass |
| M11 RSCU(CTG)=6 | 6.0 | pass |
| S6 arbitrary CAI | 0.47706538020472955 (recomputed) | pass (strengthened) |

### Variant/delegate consistency
`string` and `DnaSequence` overloads delegate to the same `*Core`; string overloads short-circuit null/empty to zeroed/0, DnaSequence overloads throw on null per contract. Consistent.

### Test quality audit (HARD GATE)
- **Sourced, not code-echoes**: M1–M11, C1 assert exact externally-sourced/hand-derived values within 1e-10. The E. coli (M9) and human (M10) tables are asserted against values I independently recomputed from Biopython and Kazusa this session, not from the code.
- **Defects found & fixed this session:**
  1. **S6 was bounds-only** (`>=0`, `<=1`) — a check an arbitrary in-range implementation would pass. The exact value is computable (0.47706538020472955); strengthened to assert it exactly (`Within(1e-10)`) while keeping the INV-1 bounds. (`CalculateCai_ArbitrarySequence_EqualsExactGeometricMeanAndIsBetweenZeroAndOne`)
  2. **Documented "non-ACGT codon skipped" edge case (INV-5 / doc §6.1) was untested** — added `GetStatistics_NonAcgtCodon_IsSkipped` (`CTGNNNGTT` → 2 codons, no `NNN` key). Locks the documented behaviour.
- **No green-washing**: no assertion weakened, no tolerance widened, no test skipped, no expected value bent to match output.
- **Coverage**: every public method/overload exercised; all Stage-A formula paths (CAI geo-mean, exclusions, GC1/2/3, GC3s synonymous subset, RSCU) and edge cases (empty, null seq, null ref, partial codon, non-ACGT, only-Met/Trp/stop) covered.
- **Honest green**: full unfiltered suite **6528 passed, 0 failed** (1 benchmark skipped by design); `dotnet build` 0 errors; changed test file warning-free.

### Findings / defects
No algorithm defect. Two test-quality defects (weak S6 bounds-only assertion; missing non-ACGT coverage) found and **completely fixed** in session.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** (GC3s percentage units; zero-w skip vs 0.01 floor; GC3s 6-fold subtlety — all documented, none affect sourced values).
- **Stage B: PASS** (code realises the validated formulas exactly; all cross-checked values reproduced; tests now lock sourced values).
- **End-state: CLEAN** — no algorithm defect; both test-quality defects fixed; full suite green.
- **Test-quality gate: PASS** (after fixing S6 + adding non-ACGT test).
