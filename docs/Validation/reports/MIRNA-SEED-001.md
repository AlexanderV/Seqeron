# Validation Report: MIRNA-SEED-001 — Seed Sequence Analysis

- **Validated:** 2026-06-12   **Area:** MiRNA
- **Canonical method(s):** `MiRnaAnalyzer.GetSeedSequence`, `MiRnaAnalyzer.CreateMiRna`, `MiRnaAnalyzer.CompareSeedRegions`
  (src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **End state:** CLEAN (no biological defect; documentation/labeling nuance noted, code values correct)

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia: MicroRNA** (cites primary refs): "nucleotides 2–7 of the miRNA (its 'seed region')" — the **canonical seed = positions 2–7, a 6-mer**.
- **Bartel (2009) Cell 136:215-233** (via PMC excerpts): "pairing to miRNA nucleotides 2–7, referred to as the miRNA seed, is critical for target recognition, with an additional pair to miRNA position 8 or an A across from miRNA position 1 often enhancing targeting efficacy." Family definition: "miRNAs with the same sequence at nucleotides 2–8."
- **TargetScan docs (7mer.html)**: 7mer-m8 = "exact match to positions 2-8 (the seed + position 8)"; 7mer-A1 = "exact match to positions 2-7 (the seed) followed by an 'A'". 8mer = positions 2-8 + A1; 6mer = positions 2-7. All positions 1-based from the 5' end of the mature miRNA.

### Definitions confirmed
- **Canonical seed = positions 2–7 (6 nt)**; the 7-nt "extended/family" region = positions 2–8 (adds m8).
- Match-type taxonomy (TargetScan / Bartel 2009): **6mer** (2–7), **7mer-A1** (2–7 + A1), **7mer-m8** (2–8), **8mer** (2–8 + A1). Effectiveness 8mer > 7mer-m8 > 7mer-A1 > 6mer.
- Family = identical at positions 2–8; any mismatch ⇒ different family.

### Independent cross-check (hand computation)
- hsa-let-7a-5p = `UGAGGUAGUAGGUUGUAUAGUU`.
  - 1-based positions 2–7 (canonical 6-mer seed) = `GAGGUA`.
  - 1-based positions 2–8 (7-mer / family region) = `GAGGUAG`.
- M-010: `GAGGUAG` vs `AGCUUAU` → 2 matches, 5 mismatches (hand-verified, matches test).
- M-011: `GAGGUAG` vs `GAGUUAG` → 6 matches, 1 mismatch (hand-verified).

### Findings / divergences (Stage A)
- The TestSpec and Evidence call **"the seed" = positions 2–8 (7 nt)**. Authoritatively, the *canonical seed* is positions **2–7 (6 nt)**; positions 2–8 is the **"family-defining" / 7-mer seed region** (Bartel 2009 uses 2–8 specifically for family identity). This is a **labeling choice, not a biological error**: the 7-nt region the code extracts is exactly the family-defining region, and the code's `IsSameFamily` (identical 2–8) is precisely Bartel's family definition. Note recorded; no correctness impact.

## Stage B — Implementation

### Code path reviewed
- `GetSeedSequence` (MiRnaAnalyzer.cs:93-99): guards null/empty/`<8 nt` → `""`; else `Substring(1, 7).ToUpperInvariant()`.
- `CreateMiRna` (:104-115): uppercases, `T→U`, extracts seed, sets `SeedStart=1, SeedEnd=7`.
- `CompareSeedRegions` (:121-147): Hamming over min length + length-diff penalty; `IsSameFamily = seed1 == seed2`; empty-seed guard returns (0,0,false).

### Formula realised correctly? (evidence)
- `Substring(1, 7)` = 0-based indices [1..7] = **1-based positions 2–8** = `GAGGUAG` for let-7a. **Correct 5' indexing**: extraction starts at position 2 (index 1), NOT position 1. The "positions 1–6 off-by-one" biological error is **absent**. The 7-nt extracted region is the family-defining seed (Bartel 2009, positions 2–8).
- `CompareSeedRegions`: Hamming distance + `IsSameFamily = identical 2–8` matches the sourced family definition.

### Cross-verification table recomputed vs code (all PASS)
| Case | Input | Expected (source) | Code result |
|------|-------|-------------------|-------------|
| let-7a seed | UGAGGUAGUAGGUUGUAUAGUU | GAGGUAG (pos 2-8) | GAGGUAG ✅ |
| miR-21 seed | UAGCUUAU… | AGCUUAU | AGCUUAU ✅ |
| let-7a vs miR-21 | GAGGUAG/AGCUUAU | 2 match / 5 mismatch | 2 / 5 ✅ |
| single mismatch | GAGGUAG/GAGUUAG | 6 / 1 | 6 / 1 ✅ |
| let-7 family | a/b/c | all GAGGUAG, same family | ✅ |
| 8-nt boundary | ABCDEFGH | BCDEFGH | BCDEFGH ✅ |
| short (<8 nt) | UAGCAUU (7) | "" | "" ✅ |

### Edge cases (all handled)
- null/empty/`<8 nt` → `""` (guard at :95). 8-nt boundary → 7-nt seed (works). U/T handled via `T→U` in `CreateMiRna` and case-normalised to upper. Empty-seed comparison → graceful (0,0,false).

### Variant/delegate consistency
- `GroupBySeedFamily` groups on `SeedSequence` (2–8) — consistent with `IsSameFamily`.
- `FindTargetSites` / `TargetSiteType` enum independently implement the full 6mer/7mer-A1/7mer-m8/8mer/offset-6mer taxonomy with correct positions (out of this unit's M/S/C scope but cross-checked: definitions align with Stage A).

### Test quality audit
- Tests assert **exact sourced values** against miRBase reference seeds (not tautologies); Hamming distances hand-verifiable; deterministic; cover null/empty/short/boundary edge cases. 29/29 parametric cases pass.

### Findings / defects (Stage B)
- **No biological defect.** Two documentation nuances (no correctness impact, not fixed to avoid churn on a green, sourced spec):
  1. XML doc comment on `GetSeedSequence` says "(positions 2-8)"; the record `SeedStart=1, SeedEnd=7` are **0-based inclusive** indices for that same region. Internally consistent but mixes 0-based fields with a 1-based doc comment.
  2. Spec/Evidence label positions 2–8 as "the seed"; canonically that is the 7-mer/family region (seed proper = 2–7). The implemented behaviour matches Bartel's family definition, so the algorithm is correct as built.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — seed/match-type biology confirmed against Bartel 2009 + TargetScan; note: "seed" labeled as 2–8 (family region) vs canonical 2–7.
- **Stage B: PASS-WITH-NOTES** — code extracts positions 2–8 with correct 5' indexing (no off-by-one), `IsSameFamily` = Bartel family definition; all cross-checks and edge cases pass.
- **End state: CLEAN** — no defect; only documentation/labeling notes. No code changed.
- Tests: SeedAnalysis filter 29/29 pass; full suite 4486 passed / 0 failed.
