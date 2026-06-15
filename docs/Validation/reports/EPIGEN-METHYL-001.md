# Validation Report: EPIGEN-METHYL-001 — Methylation Site Detection, Context Classification (CpG/CHG/CHH), Methylation Profile

- **Validated:** 2026-06-15   **Area:** Epigenetics
- **Canonical method(s):** `EpigeneticsAnalyzer.GetMethylationContext(string, int)`, `EpigeneticsAnalyzer.FindMethylationSites(string)`, `EpigeneticsAnalyzer.GenerateMethylationProfile(IEnumerable<MethylationSite>)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session
- **IUPAC nucleotide ambiguity codes** — `https://www.bioinformatics.org/sms/iupac.html` (WebFetch, 2026-06-15). Confirms **H = A or C or T** (i.e. excludes G). Independent of the repo's cited Los Alamos table, agreeing with Cornish-Bowden (1985), Nucleic Acids Res. 13:3021–3030.
- **Bismark context definition / Lister 2009 contexts** — WebSearch returning PMC3102221 (Krueger & Andrews 2011) and corroborating plant-methylation literature (2026-06-15). Confirms: Bismark "discriminate[s] between cytosines in CpG, CHG and CHH context"; "**H can be either A, T or C**"; CpG and CHG are symmetric, CHH asymmetric.
- **Schultz et al. (2012) weighted methylation level** — definition reproduced in the repo Evidence and consistent with the standard formula: WML = (Σ methylated reads) / (Σ total reads) over all cytosines of a context.

### Formula / definition check
- **CpG** = C immediately followed by G. **CHG** = C, H, G. **CHH** = C, H, H, with H ∈ {A, C, T}. Matches Bismark/Lister exactly.
- **H excludes G** confirmed against an independent IUPAC table — a G in the second position makes the site CpG, never CHG/CHH.
- **Weighted methylation level** = Σ(level·coverage)/Σ(coverage), equal to the unweighted mean of per-site fractions only under equal coverage. Matches Schultz (2012).
- **Per-cytosine level** ∈ [0,1] = methylated/total. Matches convention.

### Edge-case semantics
- Non-cytosine index → no context (only Cs have a methylation context). Sourced.
- Incomplete trailing context (e.g. `CA` at end) → unclassified, because CHG vs CHH needs the third base; a terminal `CG` is still an unambiguous CpG. Sourced (Bismark context geometry).
- Non-ACGT base in any context position → undetermined context, not a classifiable site. Sourced (IUPAC + Bismark).

### Independent cross-check (numbers, hand-computed this session)
- `CGACAGCAA`: C@0→`CG`→CpG; C@3→`CAG` (H=A,then G)→CHG; C@6→`CAA` (H=A,H=A)→CHH; no other C → **3 sites** at positions 0,3,6. ✓
- WML (cov 10 each, levels 0.8 & 0.2): (0.8·10 + 0.2·10)/20 = **0.5**. ✓
- WML unequal coverage (0.8·90 + 0.2·10)/100 = **0.74** (≠ 0.5 unweighted mean — distinguishes weighted from naive mean). ✓
- Global WML over CpG/CHG/CHH (levels 0.9,0.5,0.1 cov 10): 15/30 = **0.5**. ✓

### Findings / divergences
None. The description is biologically and mathematically correct and traces to independent authoritative sources.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/EpigeneticsAnalyzer.cs`
- `GetMethylationContext` lines 179–216
- `FindMethylationSites` lines 229–251
- `GenerateMethylationProfile` lines 504–562 (`WeightedMethylationLevel` 551–562)

### Formula realised correctly?
- `GetMethylationContext`: null guards (null/empty/index OOB) → null; non-C → null; need 1 downstream base; `next=='G'`→CpG; else require `IsHBase(next)` (A/C/T) else null; need third base; `third=='G'`→CHG; else `IsHBase(third)`→CHH else null. This is exactly the validated CpG/CHG/CHH logic with H excluding G. ✓
- `FindMethylationSites`: delegates classification to `GetMethylationContext`, emits 0-based Position at each classifiable C; sequence-only level/coverage = 0 (documented representational default). ✓
- `GenerateMethylationProfile`: per-context weighted level via `WeightedMethylationLevel` = Σ(level·coverage)/Σ(coverage); empty → all-zero profile; `MethylatedCpGSites` uses descriptive 0.5 cutoff (documented assumption, count field only). Matches Schultz (2012). ✓

### Cross-verification table (recomputed vs code, via tests)
| Case | Expected (sourced) | Code result | Match |
|------|--------------------|-------------|-------|
| `GetMethylationContext("ACGT",1)` | CpG | CpG | ✓ |
| `GetMethylationContext("CAG",0)` | CHG | CHG | ✓ |
| `GetMethylationContext("CAA",0)` | CHH | CHH | ✓ |
| `GetMethylationContext("CGG",0)` | CpG (H≠G) | CpG | ✓ |
| `GetMethylationContext("CNG",0)` | None | None | ✓ |
| `GetMethylationContext("CAN",0)` | None | None | ✓ |
| `GetMethylationContext("CA",0)` | None | None | ✓ |
| `FindMethylationSites("CGACAGCAA")` | CpG@0,CHG@3,CHH@6 | same | ✓ |
| Profile WML equal cov | 0.5 | 0.5 | ✓ |
| Profile WML unequal cov | 0.74 | 0.74 | ✓ |
| Profile per-context + global | 0.9/0.5/0.1, global 0.5 | same | ✓ |
| Empty profile | all-zero | all-zero | ✓ |

### Variant/delegate consistency
`FindMethylationSites` reuses `GetMethylationContext`, so classification is identical between the two. No competing implementation. (Other methods in the class — `SimulateBisulfiteConversion`, `CalculateMethylationFromBisulfite`, DMR/chromatin — belong to other test units and are out of this unit's scope.)

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** exact equality on context enums, on positions (0,3,6), and on WML values (0.5, 0.74, 0.9/0.5/0.1, global 0.5). The unequal-coverage test (0.74) would fail a wrong unweighted-mean implementation, so it genuinely tests the formula, not the code.
- **No green-washing:** all assertions are exact (`Is.EqualTo`, `.Within(1e-10)`); no ranges where exact values are known; nothing skipped/weakened.
- **Coverage:** all Stage-A branches now exercised. Two gaps were found and fixed this session:
  1. The CHH invalid-third-base branch (`return null` when the third base is non-ACGT) was untested — added `GetMethylationContext_NonAcgtThirdBase_ReturnsNull` (`CAN` → null).
  2. `FindMethylationSites_PositionsAreZeroBasedAtCytosine` used a bare `All(...)` that is vacuously true on an empty result — strengthened with `Is.Not.Empty` guard inside `Assert.Multiple`.
- **Honest green:** full unfiltered suite `Failed: 0, Passed: 6478`; methylation fixture 22 tests pass; `dotnet build` 0 errors (4 pre-existing warnings in unrelated test files, none in changed files).

### Findings / defects
No implementation defect. Two minor test-coverage gaps fixed in-session (see audit). No code change to the algorithm was needed or made.

## Verdict & follow-ups
- **Stage A: PASS.** Description matches independent authoritative sources exactly.
- **Stage B: PASS.** Code faithfully realises the validated formulas; cross-checks reproduce hand-computed sourced values.
- **End-state: CLEAN.** No algorithm defect; two test gaps closed; build + full suite green.
