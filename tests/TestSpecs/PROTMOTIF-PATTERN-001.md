# Test Specification: PROTMOTIF-PATTERN-001

**Test Unit ID:** PROTMOTIF-PATTERN-001
**Area:** ProteinMotif
**Algorithm:** Protein Pattern Matching Methods (FindMotifByPattern, FindMotifByProsite, ConvertPrositeToRegex, FindDomains)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | PROSITE pattern syntax (ScanProsite documentation) | 2 | https://prosite.expasy.org/scanprosite/scanprosite_doc.html | 2026-06-14 |
| 2 | de Castro et al. (2006) ScanProsite | 1 | https://doi.org/10.1093/nar/gkl124 | 2026-06-14 |
| 3 | Schneider & Stephens (1990) Sequence logos | 1 | https://doi.org/10.1093/nar/18.20.6097 | 2026-06-14 |
| 4 | PROSITE PS00001 | 5 | https://prosite.expasy.org/PS00001 | 2026-06-14 |
| 5 | PROSITE PS00005 | 5 | https://prosite.expasy.org/PS00005 | 2026-06-14 |
| 6 | PROSITE PS00016 | 5 | https://prosite.expasy.org/PS00016 | 2026-06-14 |
| 7 | PROSITE PS00017 | 5 | https://prosite.expasy.org/PS00017 | 2026-06-14 |
| 8 | PROSITE PS00029 | 5 | https://prosite.expasy.org/PS00029 | 2026-06-14 |

### 1.2 Key Evidence Points

1. PROSITE PA-line grammar: `x`→`.`, `[..]`→`[..]`, `{..}`→`[^..]`, `-` is a dropped separator, `x(n)`→`.{n}`, `x(n,m)`→`.{n,m}`, `A(n)`→`A{n}`, `<`→`^`, `>`→`$`, trailing `.` terminates the pattern. — Source 1.
2. Worked patterns: PS00001 `N-{P}-[ST]-{P}`→`N[^P][ST][^P]`; PS00005 `[ST]-x-[RK]`→`[ST].[RK]`; PS00016 `R-G-D`→`RGD`; PS00017 `[AG]-x(4)-G-K-[ST]`→`[AG].{4}GK[ST]`; PS00029 `L-x(6)-L-x(6)-L-x(6)-L`→`L.{6}L.{6}L.{6}L`. — Sources 4–8.
3. PROSITE patterns are realized as regular expressions. — Source 2.
4. Per-position information content IC = log2(N) − Σ p·log2(p); for k uniformly-allowed protein residues, IC = log2(20/k) bits; max log2(20) ≈ 4.321928094887363 bits. — Source 3.
5. The `*` Kleene star (`<{C}*>`) is a ScanProsite *query* extension, not part of the PA-line grammar. — Source 1.

### 1.3 Documented Corner Cases

- Trailing `.` ends the pattern; content after it is ignored (Source 1).
- Ranges `x(n,m)` apply only to `x`; fixed `A(n)` is allowed (Source 1).
- A fixed residue (k=1) contributes log2(20) bits; a wildcard (k=20) contributes 0 bits (Source 3).

### 1.4 Known Failure Modes / Pitfalls

1. Silently dropping unsupported metacharacters (`*`, `?`, `+`) would mis-parse a pattern — the converter must reject them. — Source 1 (extended-query metachar note).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindMotifByPattern(sequence, regexPattern, name, id)` | ProteinMotifFinder | **Canonical** | Regex match enumeration + IC scoring + E-value |
| `ConvertPrositeToRegex(prositePattern)` | ProteinMotifFinder | **Canonical** | PROSITE PA-line → .NET regex |
| `FindMotifByProsite(sequence, prositePattern, name)` | ProteinMotifFinder | **Delegate** | Composes ConvertPrositeToRegex + FindMotifByPattern |
| `FindDomains(sequence)` | ProteinMotifFinder | **Delegate** | Runs fixed signature patterns through FindMotifByPattern |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | Null/empty sequence or null/empty pattern → empty match enumeration (no exception). | Yes | Repository contract |
| INV-02 | Reported match `Sequence` equals `sequence.Substring(Start, End−Start+1)` (uppercased), with `End = Start + len − 1`. | Yes | Repository contract |
| INV-03 | `Score` = Σ over pattern positions of log2(20/k_i) where k_i = allowed residues (1 for fixed, class size for `[..]`, 20 for `x`). | Yes | Schneider & Stephens (1990) |
| INV-04 | `ConvertPrositeToRegex` maps each PROSITE atom to its regex per Source 1 (deterministic, exact). | Yes | Source 1 + PS00001/05/16/17/29 |
| INV-05 | Matching is case-insensitive; positions are 0-based. | Yes | Repository contract |
| INV-06 | Unsupported PA-line metacharacters (`*`,`?`,`+`) raise `FormatException` (reject, don't silently drop). | Yes | Source 1 |
| INV-07 | E-value ≥ 0 and equals (N−L+1)·2^(−Score). | Yes | Schneider & Stephens (1990); model def. |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Convert PS00016 literal | `ConvertPrositeToRegex("R-G-D")` | `"RGD"` | Source 6 |
| M2 | Convert PS00001 exclusion+class | `ConvertPrositeToRegex("N-{P}-[ST]-{P}")` | `"N[^P][ST][^P]"` | Source 4 |
| M3 | Convert PS00017 fixed range | `ConvertPrositeToRegex("[AG]-x(4)-G-K-[ST]")` | `"[AG].{4}GK[ST]"` | Source 7 |
| M4 | Convert PS00029 repeated x(6) | `ConvertPrositeToRegex("L-x(6)-L-x(6)-L-x(6)-L")` | `"L.{6}L.{6}L.{6}L"` | Source 8 |
| M5 | Convert anchors + trailing period | `ConvertPrositeToRegex("<M-x-K>")`/`"R-G-D."` | `"^M.K$"` / `"RGD"` | Source 1 |
| M6 | FindMotifByPattern literal position | `FindMotifByPattern("AAARGDAAA","RGD")` | 1 match, Start=3, End=5, Seq="RGD" | Sources 3,6 |
| M7 | FindMotifByPattern IC score (fixed) | Score of RGD match (3 fixed residues) | 3·log2(20) = 12.965784284662089 | Source 3 |
| M8 | FindMotifByPattern IC score (class) | Score of `[ST].[RK]` match | 2·log2(10) = 6.643856189774724 | Source 3 |
| M9 | FindMotifByProsite end-to-end PS00001 | `FindMotifByProsite("AANASAAANGTAAAA","N-{P}-[ST]-{P}")` | starts {2,8}, seqs {NASA,NGTA} | Source 4 |
| M10 | FindDomains P-loop detection | `FindDomains` on synthetic `[AG]x(4)GK[ST]` seq | 1 "Protein Kinase ATP-binding" domain at first A/G | Source 7 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Overlapping enumeration | `FindMotifByPattern("AAAA","A.A")` | matches at starts 0 and 1 | lookahead contract |
| S2 | E-value formula | E-value of RGD match | (9−3+1)·2^(−Score)=7·2^(−12.96578...) | INV-07 |
| S3 | Reject Kleene star | `ConvertPrositeToRegex("<{C}*>")` | throws FormatException | Source 1 |
| S4 | Case-insensitive | lower/mixed case same matches as upper | identical Start | INV-05 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Null/empty inputs | each method on null/empty | empty / no throw | INV-01 |
| C2 | Substring invariant | match.Sequence == substring(Start,len) | holds | INV-02 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Pre-existing (not registered to this unit) files for sibling units cover individual methods: `ProteinMotifFinder_MotifSearch_Tests.cs` (FIND-001), `ProteinMotifFinder_PrositePattern_Tests.cs` (PROSITE-001), `ProteinMotifFinder_DomainPrediction_Tests.cs` (DOMAIN-001), and legacy `ProteinMotifFinderTests.cs`.
- No canonical test file exists for PROTMOTIF-PATTERN-001 (the integrated pattern-matching umbrella unit). New file: `ProteinMotifFinder_FindMotifByPattern_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new canonical file |
| M2 | ❌ Missing | new canonical file |
| M3 | ❌ Missing | new canonical file |
| M4 | ❌ Missing | new canonical file |
| M5 | ❌ Missing | new canonical file |
| M6 | ❌ Missing | new canonical file |
| M7 | ❌ Missing | exact IC value, new |
| M8 | ❌ Missing | exact IC value, new |
| M9 | ❌ Missing | new canonical file |
| M10 | ❌ Missing | new canonical file |
| S1 | ❌ Missing | new canonical file |
| S2 | ❌ Missing | new canonical file |
| S3 | ❌ Missing | new canonical file |
| S4 | ❌ Missing | new canonical file |
| C1 | ❌ Missing | new canonical file |
| C2 | ❌ Missing | new canonical file |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindMotifByPattern_Tests.cs` — all PROTMOTIF-PATTERN-001 cases for the four registry methods.
- **Remove:** none. Pre-existing sibling-unit files (FIND/PROSITE/DOMAIN) remain for their own units; this unit does not duplicate their exact assertions — it adds the integrated, exact-IC-value contract for the pattern-matching method group registered as PATTERN-001.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ProteinMotifFinder_FindMotifByPattern_Tests.cs` | Canonical for PROTMOTIF-PATTERN-001 | 16 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented | ✅ Done |
| 11 | S1 | ❌ Missing | Implemented | ✅ Done |
| 12 | S2 | ❌ Missing | Implemented | ✅ Done |
| 13 | S3 | ❌ Missing | Implemented | ✅ Done |
| 14 | S4 | ❌ Missing | Implemented | ✅ Done |
| 15 | C1 | ❌ Missing | Implemented | ✅ Done |
| 16 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact regex assertion |
| M2 | ✅ Covered | exact regex assertion |
| M3 | ✅ Covered | exact regex assertion |
| M4 | ✅ Covered | exact regex assertion |
| M5 | ✅ Covered | exact regex assertion |
| M6 | ✅ Covered | exact Start/End/Sequence |
| M7 | ✅ Covered | exact IC `.Within(1e-10)` |
| M8 | ✅ Covered | exact IC `.Within(1e-10)` |
| M9 | ✅ Covered | exact starts + substrings |
| M10 | ✅ Covered | domain at exact position |
| S1 | ✅ Covered | overlapping starts |
| S2 | ✅ Covered | exact E-value `.Within(1e-10)` |
| S3 | ✅ Covered | FormatException |
| S4 | ✅ Covered | identical results |
| C1 | ✅ Covered | null/empty per method |
| C2 | ✅ Covered | substring invariant |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Overlapping-match enumeration is a repository contract (not PROSITE-mandated). | S1 |
| 2 | E-value uses uniform i.i.d. background (model-defined, not ScanProsite's Swiss-Prot frequencies). | S2 |

---

## 7. Open Questions / Decisions

1. Decision: `FindMotifByProsite` and `FindDomains` are Delegate-type (compose the canonical methods); they receive end-to-end smoke verification (M9, M10) rather than re-deriving every conversion atom (already covered by M1–M5 on `ConvertPrositeToRegex`).
2. Decision: existing sibling-unit test files are left untouched to avoid scope creep; PATTERN-001 owns one canonical file with exact-value contracts.
