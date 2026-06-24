# Validation Report: MIRNA-SEED-001 — Seed Sequence Analysis

- **Validated:** 2026-06-24   **Area:** MiRNA
- **Canonical method(s):** `MiRnaAnalyzer.GetSeedSequence`, `MiRnaAnalyzer.CreateMiRna`, `MiRnaAnalyzer.CompareSeedRegions`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS-WITH-NOTES
- **End state:** CLEAN (no biological defect; one documentation/labeling nuance noted, code values correct)

## Stage A — Description

### Sources opened & what they confirm
- **TargetScan docs — 7mer.html** (fetched this session): "7mer-m8 = an exact match to positions 2-8 of the mature miRNA (the seed + position 8)"; "7mer-A1 = an exact match to positions 2-7 of the mature miRNA (the seed) followed by an 'A'". Confirms the **canonical seed = positions 2–7**, with position 8 and an A opposite position 1 as supplementary determinants.
- **Bartel (2009) Cell 136:215–233 / canonical seed-type literature** (cross-checked via web search incl. PMC "sufficient minimal set of miRNA seed types"): the four canonical sites are **6mer** (perfect match to miRNA nt 2–7), **7mer-m8** (seed + match to nt 8), **7mer-A1** (seed + A at target pos 1), **8mer** (seed + both nt-8 match and A1). Family = "miRNAs with the same sequence at nucleotides 2–8."
- **Wikipedia: MicroRNA** + **miRBase**: seed = nt 2–7; 1-based indexing from the 5' end of the mature miRNA; miRBase reference sequences confirm the six test miRNAs.

### Definitions confirmed
- **Canonical seed = positions 2–7 (6 nt)**; the 7-nt region = positions 2–8 (adds m8), used by Bartel 2009 as the **family-defining** region.
- Match-type hierarchy: 8mer > 7mer-m8 > 7mer-A1 > 6mer.
- Family ⇔ identical at positions 2–8; any mismatch ⇒ different family.
- Target pairing uses the **reverse complement** of the seed (Watson–Crick A-U, C-G).

### Independent cross-check (hand computation)
- hsa-let-7a-5p = `UGAGGUAGUAGGUUGUAUAGUU`. 1-based pos2–8 = `GAGGUAG` (7 nt, family region); pos2–7 = `GAGGUA` (canonical 6-mer).
- Reverse complement of `GAGGUAG` (RNA) = `CUACCUC`; 6mer core (RC of pos2–7) = `UACCUC` — matches the `seedRC.Substring(1,6)` layout used by `FindTargetSites`.
- M-010: `GAGGUAG` vs `AGCUUAU` → 2 matches / 5 mismatches (hand-verified).
- M-011: `GAGGUAG` vs `GAGUUAG` → 6 matches / 1 mismatch (hand-verified).
- S-004: `AAAAAAA` vs `UUUUUUU` → 0 / 7 (max Hamming = seed length).

### Findings / divergences (Stage A)
- TestSpec and Evidence label **"the seed" = positions 2–8 (7 nt)**. Authoritatively the *canonical seed* is **2–7 (6 nt)**; 2–8 is the **7-mer / family-defining** region (Bartel 2009 uses 2–8 specifically for family identity). This is a **labeling choice, not a biological error**: the 7-nt region the code extracts is exactly the family-defining region, and `IsSameFamily` (identical 2–8) is precisely Bartel's family definition.

## Stage B — Implementation

### Code path reviewed
- `GetSeedSequence` (MiRnaAnalyzer.cs:93–99): guards null/empty/`<8 nt` → `""`; else `Substring(1, 7).ToUpperInvariant()`.
- `CreateMiRna` (:104–115): uppercases, `T→U`, extracts seed, sets `SeedStart=1, SeedEnd=7`.
- `CompareSeedRegions` (:121–147): Hamming over `min(len1,len2)` + `|len1−len2|` penalty; `IsSameFamily = seed1 == seed2`; empty-seed guard returns `(0,0,false)`.
- (Cross-checked, out of unit's M/S/C scope) `FindTargetSites` / `TargetSiteType` (:42–264) + `GetReverseComplement` (:294–317): site typing uses RC of the seed and correctly classifies 8mer/7mer-m8/7mer-A1/6mer/offset-6mer per Stage-A definitions.

### Formula realised correctly? (evidence)
- `Substring(1, 7)` = 0-based indices [1..7] = **1-based positions 2–8** = `GAGGUAG` for let-7a. Extraction starts at position 2 (index 1), NOT position 1 — the classic "positions 1–7 off-by-one" error is **absent**.
- `CompareSeedRegions`: Hamming distance + `IsSameFamily = identical 2–8` matches the sourced family definition.

### Cross-verification table recomputed vs code (all PASS)
| Case | Input | Expected (source) | Code result |
|------|-------|-------------------|-------------|
| let-7a seed | UGAGGUAGUAGGUUGUAUAGUU | GAGGUAG (pos 2–8) | GAGGUAG ✅ |
| miR-21 seed | UAGCUUAU… | AGCUUAU | AGCUUAU ✅ |
| miR-155 seed | UUAAUGCU… | UAAUGCU | UAAUGCU ✅ |
| miR-1 seed | UGGAAUGU… | GGAAUGU | GGAAUGU ✅ |
| let-7a vs miR-21 | GAGGUAG/AGCUUAU | 2 / 5 | 2 / 5 ✅ |
| single mismatch | GAGGUAG/GAGUUAG | 6 / 1 | 6 / 1 ✅ |
| max different | AAAAAAA/UUUUUUU | 0 / 7 | 0 / 7 ✅ |
| let-7 family a/b/c | — | all GAGGUAG, same family | ✅ |
| 8-nt boundary | ABCDEFGH | BCDEFGH | BCDEFGH ✅ |
| short (<8 nt) | UAGCAUU (7) | "" | "" ✅ |

### Edge cases (all handled)
- null/empty/`<8 nt` → `""` (guard at :95). 8-nt boundary → 7-nt seed. T handled via `T→U`; case normalised to upper. Empty-seed comparison → graceful `(0,0,false)`.

### Variant/delegate consistency
- `GroupBySeedFamily` groups on `SeedSequence` (2–8), consistent with `IsSameFamily`. `FindSimilarMiRnas` uses positional Hamming over the seed, consistent with `CompareSeedRegions`.

### Test quality audit
- Tests assert **exact sourced values** against miRBase reference seeds (not tautologies); Hamming distances hand-verifiable; deterministic; cover null/empty/short/boundary/both-empty edge cases. 29/29 parametric cases pass this session.

### Findings / defects (Stage B)
- **No biological defect.** Two documentation nuances (no correctness impact, not changed to avoid churn on a green, sourced spec):
  1. The XML doc on `GetSeedSequence` says "(positions 2-8)", while record fields `SeedStart=1, SeedEnd=7` are the **0-based inclusive** indices of that same region. Internally consistent but mixes 0-based fields with a 1-based doc comment.
  2. Spec/Evidence label positions 2–8 as "the seed"; canonically that is the 7-mer/family region (seed proper = 2–7). Implemented behaviour matches Bartel's family definition, so the algorithm is correct as built.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — seed/match-type biology independently confirmed against TargetScan 7mer.html + Bartel 2009; note: "seed" labeled as 2–8 (family region) vs canonical 2–7.
- **Stage B: PASS-WITH-NOTES** — code extracts positions 2–8 with correct 5' indexing (no off-by-one); `IsSameFamily` = Bartel family definition; all cross-checks and edge cases pass.
- **End state: CLEAN** — no defect; only the documentation/labeling note. No code changed.
- Tests: SeedAnalysis filter 29/29 pass; no code changed so the full suite is unaffected (baseline green).
