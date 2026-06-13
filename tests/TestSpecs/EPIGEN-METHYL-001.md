# Test Specification: EPIGEN-METHYL-001

**Test Unit ID:** EPIGEN-METHYL-001
**Area:** Epigenetics
**Algorithm:** Methylation Site Detection, Context Classification (CpG/CHG/CHH), Methylation Profile
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Cornish-Bowden (1985), IUPAC-IUB nucleotide nomenclature | 2 | https://doi.org/10.1093/nar/13.9.3021 (table: https://www.hiv.lanl.gov/content/sequence/HelpDocs/IUPAC.html) | 2026-06-13 |
| 2 | Krueger & Andrews (2011), Bismark | 3 | https://doi.org/10.1093/bioinformatics/btr167 (PMC3102221) | 2026-06-13 |
| 3 | Lister et al. (2009), Human DNA methylomes | 1 | https://doi.org/10.1038/nature08514 (PMID 19829295) | 2026-06-13 |
| 4 | Schultz et al. (2012), 'Leveling' the playing field | 1 | https://doi.org/10.1016/j.tig.2012.10.012 | 2026-06-13 |

### 1.2 Key Evidence Points

1. IUPAC code **H = A, C or T** ("not G") — Cornish-Bowden (1985) [1].
2. Cytosine sequence contexts: **CpG** (C followed by G), **CHG** (C, H, G), **CHH** (C, H, H); CpG/CHG symmetric, CHH asymmetric; "H can be either A, T or C" — Krueger & Andrews (2011) [2], Lister et al. (2009) [3].
3. **Weighted methylation level** = (Σ methylated reads) / (Σ total reads) over all cytosines of a context — Schultz et al. (2012) [4].
4. **Per-cytosine methylation level** = methylated reads / total reads ∈ [0,1] — Schultz (2012) [4].
5. Prevalence: IMR90 somatic = 99.98% CG; H1 ES = ~25% non-CG — Lister (2009) [3].

### 1.3 Documented Corner Cases

- **H excludes G**: a G in the H position changes the context (CpG vs CHG) — [1][2].
- **Incomplete trailing context**: a C at end with truncated downstream cannot be classed CHG vs CHH; a terminal `CG` is still unambiguous CpG — [2].
- **Non-ACGT bases** in the context positions → undetermined context, not a classifiable site — [1][2].

### 1.4 Known Failure Modes / Pitfalls

1. Treating GpC as CpG — only 5'→3' C-then-G is CpG [2].
2. Classifying `C x G` with x = G as CHG instead of CpG — H must exclude G [1].
3. Equating weighted methylation level with the unweighted mean of per-site fractions — they differ under unequal coverage [4].

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GetMethylationContext(sequence, index)` | EpigeneticsAnalyzer | Canonical | Classifies the cytosine at `index` into CpG/CHG/CHH (or none). |
| `FindMethylationSites(sequence)` | EpigeneticsAnalyzer | Canonical | Enumerates classifiable C sites with Type + 0-based Position. |
| `GenerateMethylationProfile(sites)` | EpigeneticsAnalyzer | Canonical | Aggregates per-context methylation (Schultz 2012 WML + counts). Registry name: `CalculateMethylationProfile`. |

<!-- Type values: Canonical / Delegate / Internal -->

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | A cytosine is CpG ⟺ next base is G | Yes | [2][3] |
| INV-2 | CHG/CHH require the H position(s) ∈ {A,C,T}; a G there is never H | Yes | [1] |
| INV-3 | Every site Position is 0-based and points at a `C` | Yes | [2] |
| INV-4 | Each per-site MethylationLevel ∈ [0,1] | Yes | [4] |
| INV-5 | Weighted methylation level = Σmeth/Σtotal over the context's sites; equals the mean of fractions only when coverage is equal | Yes | [4] |
| INV-6 | Context classification is case-insensitive and deterministic | Yes | [2] |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Context CpG | `GetMethylationContext("ACGT",1)` | CpG | [2][3] |
| M2 | Context CHG | `GetMethylationContext("CAG",0)` (C,A,G) | CHG | [1][2] |
| M3 | Context CHH | `GetMethylationContext("CAA",0)` (C,A,A) | CHH | [2] |
| M4 | H excludes G | `GetMethylationContext("CGG",0)` → next is G → CpG (not CHG) | CpG | [1] |
| M5 | CHG with H=C | `GetMethylationContext("CCG",0)` (C,C,G) | CHG | [1][2] |
| M6 | CHG with H=T | `GetMethylationContext("CTG",0)` (C,T,G) | CHG | [1][2] |
| M7 | Non-C index | `GetMethylationContext("ACGT",0)` (A, not C) | None | [2] |
| M8 | Incomplete trailing (CHH undecidable) | `GetMethylationContext("CA",0)` (only 2 bases, not CG) | None | [2] |
| M9 | Terminal CG still CpG | `GetMethylationContext("AACG",2)` | CpG | [2] |
| M10 | FindMethylationSites all contexts | `FindMethylationSites("CGACAGCAA")` → CpG@0, CHG@3, CHH@6 | 3 sites, exact (Type,Position) | dataset, [2] |
| M11 | FindMethylationSites positions are 0-based at C | each site.Position indexes a 'C' | true for all sites | [2] |
| M12 | Profile weighted CpG level | sites CpG(level .8,cov 10),(level .2,cov 10) → WML 0.5 | CpGMethylation = 0.5 | [4] |
| M13 | Profile counts | 2 CpG sites, one ≥0.5 → TotalCpGSites=2, MethylatedCpGSites=1 | exact counts | [4] (count cutoff = assumption) |
| M14 | Profile empty sites | `GenerateMethylationProfile([])` | all-zero profile | contract |
| M15 | Null/empty sequence | `FindMethylationSites(null/"")` | empty | contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Case-insensitive context | `GetMethylationContext("acgt",1)` | CpG | INV-6 |
| S2 | FindMethylationSites lowercase | `FindMethylationSites("cgacagcaa")` | same 3 sites as M10 | INV-6 |
| S3 | Profile per-context split | sites of CpG+CHG+CHH with distinct levels → each context averaged independently | exact CpG/CHG/CHH levels | [4] |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Non-ACGT in context | `GetMethylationContext("CNG",0)` (N not in {A,C,T}) | None | [1] |
| C2 | Index out of range | `GetMethylationContext("ACGT",10)` | None | contract |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzerTests.cs` — legacy fixture; one weak `FindMethylationSites_IdentifiesAllContexts` (no `Assert.Multiple`, no messages) and `GenerateMethylationProfile_*` tests using non-evidence inputs. Out of scope of this unit's canonical file; not modified.
- No `GetMethylationContext` method or canonical `EpigeneticsAnalyzer_Methylation_Tests.cs` existed before this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new method |
| M2 | ❌ Missing | |
| M3 | ❌ Missing | |
| M4 | ❌ Missing | H-excludes-G |
| M5 | ❌ Missing | |
| M6 | ❌ Missing | |
| M7 | ❌ Missing | |
| M8 | ❌ Missing | |
| M9 | ❌ Missing | |
| M10 | ⚠ Weak | legacy test exists but no messages / no Assert.Multiple / no count assertion |
| M11 | ❌ Missing | |
| M12 | ❌ Missing | weighted level |
| M13 | ❌ Missing | |
| M14 | ❌ Missing | |
| M15 | ❌ Missing | |
| S1 | ❌ Missing | |
| S2 | ❌ Missing | |
| S3 | ❌ Missing | |
| C1 | ❌ Missing | |
| C2 | ❌ Missing | |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/EpigeneticsAnalyzer_Methylation_Tests.cs` — all M/S/C cases for the three methods under test.
- **Remove:** nothing — legacy `EpigeneticsAnalyzerTests.cs` belongs to no Test Unit and is left untouched (out of scope).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `EpigeneticsAnalyzer_Methylation_Tests.cs` | Canonical (this unit) | 19 |
| `EpigeneticsAnalyzerTests.cs` | Legacy (untouched) | unchanged |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ | implemented | ✅ Done |
| 2 | M2 | ❌ | implemented | ✅ Done |
| 3 | M3 | ❌ | implemented | ✅ Done |
| 4 | M4 | ❌ | implemented | ✅ Done |
| 5 | M5 | ❌ | implemented | ✅ Done |
| 6 | M6 | ❌ | implemented | ✅ Done |
| 7 | M7 | ❌ | implemented | ✅ Done |
| 8 | M8 | ❌ | implemented | ✅ Done |
| 9 | M9 | ❌ | implemented | ✅ Done |
| 10 | M10 | ⚠ | rewritten from scratch (exact, Assert.Multiple, messages) | ✅ Done |
| 11 | M11 | ❌ | implemented | ✅ Done |
| 12 | M12 | ❌ | implemented | ✅ Done |
| 13 | M13 | ❌ | implemented | ✅ Done |
| 14 | M14 | ❌ | implemented | ✅ Done |
| 15 | M15 | ❌ | implemented | ✅ Done |
| 16 | S1 | ❌ | implemented | ✅ Done |
| 17 | S2 | ❌ | implemented | ✅ Done |
| 18 | S3 | ❌ | implemented | ✅ Done |
| 19 | C1 | ❌ | implemented | ✅ Done |
| 20 | C2 | ❌ | implemented | ✅ Done |

**Total items:** 20
**✅ Done:** 20 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact CpG |
| M2 | ✅ | exact CHG |
| M3 | ✅ | exact CHH |
| M4 | ✅ | CpG (H excludes G) |
| M5 | ✅ | CHG H=C |
| M6 | ✅ | CHG H=T |
| M7 | ✅ | None at non-C |
| M8 | ✅ | None on truncated |
| M9 | ✅ | terminal CpG |
| M10 | ✅ | 3 sites exact |
| M11 | ✅ | positions index 'C' |
| M12 | ✅ | WML 0.5 |
| M13 | ✅ | counts exact |
| M14 | ✅ | empty profile |
| M15 | ✅ | null/empty empty |
| S1 | ✅ | lowercase CpG |
| S2 | ✅ | lowercase sites |
| S3 | ✅ | per-context split |
| C1 | ✅ | N → None |
| C2 | ✅ | out-of-range → None |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Sequence-only sites carry MethylationLevel/Coverage = 0 (no bisulfite data) | `FindMethylationSites` default values |
| 2 | `MethylatedCpGSites` uses a descriptive 0.5 fractional cutoff (Schultz recommends a binomial test) | profile count field only; does not affect continuous levels |

---

## 7. Open Questions / Decisions

1. Registry lists the profile method as `CalculateMethylationProfile`; the implemented/established name is `GenerateMethylationProfile` (used by MCP tools and prior tests). Decision: keep `GenerateMethylationProfile` for API stability; documented here and in the Processing Registry note.
