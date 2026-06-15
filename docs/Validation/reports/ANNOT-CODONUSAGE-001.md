# Validation Report: ANNOT-CODONUSAGE-001 — Relative Synonymous Codon Usage (RSCU)

- **Validated:** 2026-06-15   **Area:** Annotation
- **Canonical method(s):** `GenomeAnnotator.GetCodonUsage(IEnumerable<string>)` and overload `GetCodonUsage(IEnumerable<string>, GeneticCode)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES (code correct; closed three test-coverage gaps)
- **End-state:** ✅ CLEAN

## Stage A — Description

### Sources opened this session (independent of the repo's own Evidence)

| Source | Retrieval | What it confirms |
|--------|-----------|------------------|
| LIRMM / Rivals "RSCU RS" methods page — https://www.lirmm.fr/~rivals/rscu/ | WebFetch | Verbatim formula `RSCU_{i,j} = n_i·x_{i,j} / Σ_{j=1..n_i} x_{i,j}`; `n_i` = "the number of codons that code amino acid i"; `x_{i,j}` = "the number of occurrences of codon j"; range "between 0 and the number of synonymous codons for that amino acid". |
| WebSearch restatement of the Sharp et al. definition | WebSearch | Same formula and variable definitions; "RSCU = ratio of observed frequency to the expected frequency given that all synonymous codons are used equally"; attributed to Sharp et al. |
| cubar `est_rscu` (CRAN) — https://rdrr.io/cran/cubar/man/est_rscu.html | WebFetch | Independent reference implementation: "observed count divided by the mean count of all synonymous codons" — algebraically identical to `n_i·x/Σx` (mean = Σx/n_i). RSCU>1 over-, <1 under-represented. |
| NCBI The Genetic Codes — Standard table 1 — https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi | WebFetch | Leu = TTA/TTG/CTT/CTC/CTA/CTG (6), Ser = 6, Arg = 6, Met = ATG only, Trp = TGG only, stops = TAA/TAG/TGA. |

### Formula check
The implementation and description use `RSCU = n_i·x_{i,j} / Σ_j x_{i,j}`. This matches the LIRMM page **verbatim** (symbols, normalisation by the family total, multiplication by family size `n_i`). The "observed / expected-under-uniform" phrasing (PMC2528880, cubar) is the same quantity: expected = Σx/n_i, so observed/expected = n_i·x/Σx. Confirmed.

### Edge-case semantics
- Single-codon amino acids (Met=ATG, Trp=TGG): `n_i=1` ⇒ RSCU = 1·x/x = 1.0 (NCBI table 1 + formula). Confirmed.
- Stop codons excluded (sense codons only): NCBI table 1 stops TAA/TAG/TGA; RSCU defined over the 61 sense codons. Confirmed.
- Range [0, n_i]: LIRMM "between 0 and the number of synonymous codons". Confirmed.
- Whole-family zero count: base RSCU has a 0 denominator (undefined). The unit reports 0.0 for each member and does **not** apply the CAI 0.5 pseudocount — correct, because the 0.5 substitution is a CAI convention (Sharp & Li 1987), not part of plain RSCU. This is a documented, sourced design choice, not a defect.

### Independent cross-check (hand computation from the formula + NCBI families)
- M1 Leu `CTTCTTCTGTTA` → CTT=2, CTG=1, TTA=1, Σ=4, n_i=6: CTT=6·2/4=**3.0**, CTG=6·1/4=**1.5**, TTA=**1.5**, others 0; Σ over family = **6.0** = n_i. ✓
- M2 Phe `TTTTTC` → TTT=1, TTC=1, n_i=2, Σ=2 → **1.0/1.0**. ✓
- M3 Met `ATGATG` → n_i=1 → **1.0**. ✓
- Trp `TGGTGGTGG` → n_i=1 → **1.0** (added test). ✓

### Findings
Stage-A description, TestSpec, and Evidence are accurate and externally corroborated. No divergence. **PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs:922–992`. Standard-code overload delegates to the `GeneticCode` overload. Counts pooled across all sequences via case-insensitive (`ToUpperInvariant`) in-frame triplet stepping (`i += 3`, stops at `Length-3` so a partial trailing codon is dropped); only A/C/G/T codons counted. Families built from `code.CodonTable` (RNA-keyed; `U→T` conversion), stop codons (`'*'`) skipped. RSCU = `nI*x/familyTotal`, with `familyTotal==0 ⇒ 0.0`. `ArgumentNullException` on null `codingSequences`/`code`. Codon table (`GeneticCode.cs:230+`) matches NCBI table 1.

### Formula realised correctly?
Yes — the code computes `(double)nI * x / familyTotal` per family, the exact `n_i·x/Σx`. Family grouping is by amino-acid char excluding `'*'`, giving exactly the 61 sense codons. Verified by running the tests: output `Count == 61`, Trp/Met = 1.0, Leu values 3.0/1.5/1.5.

### Cross-verification table recomputed vs code (via tests, `.Within(1e-10)`)
| Case | Expected (sourced) | Code output |
|------|--------------------|-------------|
| M1 CTT/CTG/TTA | 3.0 / 1.5 / 1.5; others 0 | match |
| M2 TTT/TTC | 1.0 / 1.0 | match |
| M3 ATG | 1.0 | match |
| Trp TGG | 1.0 | match |
| M4 pooled ["CTTCTT","CTGTTA"] | = M1 | match |
| M5 ATGTAA | ATG=1.0; no TAA/TAG/TGA key | match |
| M6 Σ Leu family | 6.0 | match |
| 61-codon set | Count=61, no stops | match |

### Variant/delegate consistency
S3 confirms the `GeneticCode.Standard` overload equals the default overload (delegation). The legacy `GetCodonUsage(string)` raw-count method is a different method, out of scope, left unchanged.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** all expected RSCU values (3.0, 1.5, 1.0, 6.0) trace to the LIRMM formula + NCBI families, recomputed by hand; they would fail a deliberately-wrong implementation. No tautologies.
- **No green-washing:** exact `.EqualTo(...).Within(1e-10)` assertions; no ranges/AtLeast/Contains where exact values are known; no skips/weakened tolerances.
- **Coverage gaps closed (test-only, 0 code change):**
  1. INV-02 tested Met only — added `GetCodonUsage_TryptophanSingleCodon_ReturnsOnePointZero` (TGG=1.0).
  2. INV-04/cardinality untested — added `GetCodonUsage_OutputContainsExactlyThe61SenseCodons_NoStops` (Count==61, no TAA/TAG/TGA).
  3. C2 tested only `new[]{""}` — extended to also assert the empty-enumerable branch `Array.Empty<string>()` ⇒ 61-codon all-zero map.
- **Honest green:** full unfiltered suite **6568 passed, 0 failed** (1 pre-existing skipped benchmark `MFE_Benchmark_AllScenarios`, unrelated); `dotnet build` 0 warnings / 0 errors on the changed file. Fixture grew 11 → 13 tests.

### Findings / defects
No algorithm or code defect. The implementation faithfully realises the validated RSCU formula. The only issues were three Stage-A invariants/edge branches not directly exercised by the original tests; all three are now locked with sourced assertions.

## Verdict & follow-ups
- **Stage A: PASS.** Formula verbatim-confirmed by two independent sources (LIRMM, cubar) plus the Sharp definition; families/stops confirmed against NCBI table 1.
- **Stage B: PASS-WITH-NOTES.** Code correct; closed three test-coverage gaps (Trp single-codon, 61-codon set/stop exclusion, empty-enumerable branch).
- **End-state: ✅ CLEAN.** No defect; coverage gaps fully fixed in-session; full suite green.
- No follow-ups.
