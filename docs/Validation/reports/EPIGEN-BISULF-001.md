# Validation Report: EPIGEN-BISULF-001 — Bisulfite Sequencing Analysis

- **Validated:** 2026-06-15   **Area:** Epigenetics
- **Canonical method(s):** `EpigeneticsAnalyzer.SimulateBisulfiteConversion(sequence, methylatedPositions)`, `EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(referenceSequence, bisulfiteReads)`, `EpigeneticsAnalyzer.GenerateMethylationProfile(sites)` (`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (retrieved, not trusted by label)

1. **Wikipedia — Bisulfite sequencing** (WebFetched). Verbatim: *"Treatment of DNA with bisulfite converts cytosine residues to uracil, but leaves 5-methylcytosine residues unaffected."* — *"All sites of unmethylated cytosines are displayed as thymines in the resulting amplified sequence of the sense strand."* — *"Only cytosines in single-stranded DNA are susceptible to attack by bisulfite, therefore denaturation of the DNA … is critical."* Confirms: unmeth C→U→T; 5mC stays C; single-strand specificity (justifies converting only the supplied strand).
2. **Frommer et al. (1992) PNAS 89:1827–1831** (WebSearch summaries of the primary paper). Confirms the conversion chemistry ("Cytosine is deaminated to uracil by bisulfite treatment, while 5mC does not react with bisulfite") and the strand-specific PCR read-out (thymine takes the place of cytosine).
3. **Krueger & Andrews (2011) Bioinformatics 27:1571–1572, Bismark** (WebFetched PMC3102221). Verbatim: *"The methylation state of positions involving cytosines is determined by comparing the read sequence with the corresponding genomic sequence."* — *"methylation calls in Bismark … discriminate between cytosines in CpG, CHG and CHH context."* Confirms: read-C = methylated, read-T = unmethylated at a reference C; three contexts.
4. **Bismark `bismark_methylation_extractor` source** (WebFetched, FelixKrueger/Bismark GitHub). Verbatim Perl: `$percent = … meth * 100 / ( meth + un )`. Confirms the percentage/level formula `meth/(meth+unmeth)` exactly.
5. **Schultz, Schmitz & Ecker (2012) Trends Genet. 28:583–585** (WebSearch of the paper's own worked example). The paper contrasts the **weighted methylation level** (sequencing-depth-weighted, i.e. read-pooled) with the unweighted mean of per-site levels, using a region of two sites: one with **90 methylated of 100 reads** and one with **1 of 2 reads**.

### Formula check
- Conversion (Frommer): unprotected C → T, protected (5mC) C → C, non-C unchanged, case preserved. Matches sources [1][2]. ✔
- Methylation level (Bismark): `level = methylated / (methylated + unmethylated)` ∈ [0,1]; C-call=meth, T-call=unmeth at a reference CpG C; non-C/T read base = not a valid call. Matches [3][4]. ✔
- Weighted per-context level (Schultz): `Σ(methylated reads)/Σ(total reads) = Σ(levelᵢ·coverageᵢ)/Σ(coverageᵢ)`. Matches [5]. ✔

### Edge-case semantics check
All Stage-A edge cases have sourced, defined behaviour: zero-coverage CpG excluded (Bismark — undefined percentage); read bases past reference end ignored; non-C/T read base at a reference C ignored (Krueger 2011); last reference base cannot start a CpG; null/empty sequence → `""`; empty site list → zero profile. ✔

### Independent cross-check (numbers retrieved this session)
- **Conversion (Frommer chemistry, hand-applied):** `ACGTCGAA` (C at idx 1,4), no protection → `ATGTTGAA`; protect {1} → `ACGTTGAA`. ✔
- **Bismark level:** one C-call + one T-call at a CpG → 1/(1+1) = 0.5, coverage 2. ✔
- **Schultz canonical example (the paper's own numbers):** sites 90/100 and 1/2 → **weighted = (90+1)/(100+2) = 91/102 = 0.8921568627…**, **unweighted mean = (0.90+0.50)/2 = 0.70**. The two differ — exactly the point the paper makes. This is a stronger anchor than the repo's invented (1.0,cov10)/(0.0,cov30)=0.25 example (which is also correct math, just not the source's numbers).

### Findings / divergences
None. The description (`docs/algorithms/Epigenetics/Bisulfite_Sequencing_Analysis.md`), TestSpec and Evidence all match the retrieved sources. The Bismark call-symbol mapping in the Evidence doc (`z/Z` CpG, `x/X` CHG, `h/H` CHH; lowercase=unmethylated) is the standard Bismark convention and correct. (A WebFetch summarisation garbled it; the doc itself is right.)

## Stage B — Implementation

### Code path reviewed
- `SimulateBisulfiteConversion` — `EpigeneticsAnalyzer.cs:379-413`
- `CalculateMethylationFromBisulfite` — `EpigeneticsAnalyzer.cs:431-489`
- `GenerateMethylationProfile` + `WeightedMethylationLevel` — `EpigeneticsAnalyzer.cs:504-562`
- `FindCpGSites` — `EpigeneticsAnalyzer.cs:148-162`

### Formula realised correctly? (evidence)
- **Conversion:** per-base loop; C/c not in `methylatedPositions` → `T`/`t` (case preserved), protected C/c kept, everything else copied. Matches Frommer exactly. ✔
- **Calling:** `FindCpGSites` enumerates reference CpG cytosines; per read, only C and T bases at a CpG key are tallied (C→meth+total, T→total); other bases ignored; `level = meth/total`; `total==0` sites skipped. Matches Bismark/Krueger exactly. The loop guard `startPos+i < Length-1` is harmless: every `FindCpGSites` key is ≤ `Length-2`, so no real CpG is ever excluded. ✔
- **Profile:** `WeightedMethylationLevel = Σ(level·coverage)/Σ(coverage)`, computed per context and globally. Matches Schultz. A coverage-0 fallback to the unweighted mean is used only for sequence-only sites (documented; does not affect any sourced read-coverage value). ✔

### Cross-verification table recomputed vs code (tests run)
| Case | Source value | Code result | Match |
|------|-------------|-------------|-------|
| Convert `ACGTCGAA`, none protected | `ATGTTGAA` | `ATGTTGAA` | ✔ |
| Convert `ACGTCGAA`, protect {1} | `ACGTTGAA` | `ACGTTGAA` | ✔ |
| Level: C+T at CpG | 0.5, cov 2 | 0.5, cov 2 | ✔ |
| Level: single T | 0.0, cov 1 | 0.0, cov 1 | ✔ |
| Level: single C | 1.0, cov 1 | 1.0, cov 1 | ✔ |
| **Schultz 90/100 + 1/2 weighted** | 91/102 = 0.8921568627 | 0.8921568627 | ✔ |
| Schultz unweighted mean (contrast) | 0.70 (must differ) | code ≠ 0.70 | ✔ |
| CHH weighted 2/2 + 0/8 | 0.2 | 0.2 | ✔ |
| MethylatedCpGSites (0.49, 0.80) ≥0.5 | 1 of 2 | 1 of 2 | ✔ |

### Variant/delegate consistency
No `*Fast`/delegate variants for these three methods. `CalculateMethylationFromBisulfite` reuses `FindCpGSites` (already validated under EPIGEN-METHYL-001); `GenerateMethylationProfile` consumes `MethylationSite` records produced by the caller.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** original 16 tests assert exact sourced values (M1 `ATGTTGAA`, M2 `ACGTTGAA`, M5 0.5/cov2, M6 0.0/cov1, M7 1.0/cov1, M8 weighted 0.25 with the explicit `≠0.5` contrast). The weighted-vs-mean contrast prevents a deliberately-wrong (unweighted) implementation from passing.
- **Coverage gap found (no algorithm defect):** three of the seven `MethylationProfile` fields — `GlobalMethylation`, `MethylatedCpGSites`, `CHHMethylation` — had no value assertion, and the weighted example used repo-invented numbers rather than Schultz's own.
- **Fix:** added 3 exact-value tests anchored to retrieved sources — `GenerateProfile_SchultzWorkedExample_ReturnsWeightedLevel` (the paper's own 90/100+1/2 → 91/102, asserts the weighted level ≠ unweighted mean 0.70, and locks `GlobalMethylation` + `MethylatedCpGSites=2`); `GenerateProfile_ChhContext_UsesWeightedLevel` (CHH weighting 0.2); `GenerateProfile_LevelBelowThreshold_NotCountedAsMethylatedSite` (0.49<0.5 cutoff → 1 methylated site of 2). No assertion weakened, no tolerance widened, no test skipped/ignored, no expected value tuned to code. Fixture 16 → 19 tests.
- **Honest green:** FULL unfiltered suite — **6539 passed, 0 failed** (1 benchmark `[Explicit]`-skipped, unrelated); `dotnet build` 0 errors; changed test file warning-free.

### Findings / defects
No algorithm defect. One test-coverage gap (untested profile fields + repo-invented weighted example), fixed in-session and locked to the Schultz paper's own numbers. Logged as A23 in the Findings Register.

## Verdict & follow-ups
- **Stage A: PASS** — description and formulas match retrieved authoritative sources verbatim.
- **Stage B: PASS** — code realises every validated formula exactly; all cross-check values reproduced.
- **End-state: ✅ CLEAN.** Test-quality gate satisfied (sourced expectations, no green-washing, full public surface and contexts now covered, full suite honest-green 6539/0).
