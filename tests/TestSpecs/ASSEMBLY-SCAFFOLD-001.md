# Test Specification: ASSEMBLY-SCAFFOLD-001

**Test Unit ID:** ASSEMBLY-SCAFFOLD-001
**Area:** Assembly
**Algorithm:** Scaffolding (joining ordered contigs into scaffolds with N-gaps)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-11-25 (B10 review)

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Jackman et al. (2017). ABySS 2.0. *Genome Research* 27:768–777 (scaffold construction text) | 1 | https://genome.cshlp.org/content/27/5/768 (text via http://sjackman.ca/abyss-scaffold-paper/) | 2026-06-13 |
| 2 | NCBI AGP Specification v2.1 | 2 | https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/ | 2026-06-13 |
| 3 | Sahlin et al. (2012). Improved gap size estimation. *Bioinformatics* 28(17):2215–2222 | 1 | https://academic.oup.com/bioinformatics/article/28/17/2215/246308 | 2026-06-13 |
| 4 | Pop, Kosack & Salzberg (2004). Hierarchical Scaffolding With Bambus. *Genome Research* 14(1):149–159 | 4→primary | https://en.wikipedia.org/wiki/Scaffolding_(bioinformatics) | 2026-06-13 |
| 5 | Huson, Reinert & Myers (2002). The greedy path-merging algorithm for contig scaffolding. *J. ACM* 49(5):603–615 | 1 | https://publications.mpi-cbg.de/Huson_2002_6349.pdf | 2026-11-25 |
| 6 | ABySS `MergePaths/MergeContigs.cpp` (overlap merge for negative distances) | 1 (reference code) | https://github.com/bcgsc/abyss | 2026-11-25 |
| 7 | SSPACE-Basic (Boetzer et al. 2011) `SSPACE_Basic.pl`, `PairingAndScaffolding.pl` | 1 (reference code) | https://github.com/nsoranzo/sspace_basic | 2026-11-25 |

### 1.2 Key Evidence Points

1. Scaffold = ordered contigs along a path concatenated, "interspersed with gaps represented by a run of the character N, whose length corresponds to the estimate of the distance between those two contigs." — Jackman et al. (2017)
2. The number of `N` between two contigs equals the estimated distance. — Jackman et al. (2017)
3. Gap lengths must be positive; "Negative gaps and gap lines with zero length are not valid." — NCBI AGP Spec v2.1
4. For negative / unknown gaps, "use ... 100 as the gap size, since 100 is the GenBank/EMBL/DDBJ standard for gaps of unknown size." — NCBI AGP Spec v2.1
5. A negative estimate indicates contigs should overlap (frequently one k-mer length for de Bruijn assemblers). — Jackman et al. (2017); Sahlin et al. (2012)
6. Scaffolding orders/orients contigs into a path of distinct contigs (greedy linkage). — Pop et al. (2004) / Bambus
7. Greedy path merging: mate-edges in descending weight; an edge is selected only if its end nodes lie on different selected paths and every node keeps ≤ 2 selected edges (contig edge + one link). — Huson et al. (2002) Alg. 5.1
8. A negative distance whose −d bases agree (case-insensitive, N wildcard) is merged; otherwise ABySS emits a gap. — ABySS MergeContigs; Jackman et al. (2017)

### 1.3 Documented Corner Cases

- Zero-length gap line is invalid (AGP) → a zero estimate maps to the unknown-gap default, not "0 N".
- Negative estimate = overlap (AGP/Jackman/Sahlin) → collapsed when the −d bases agree; else unknown-gap default (100 N).
- Link into a lower-index contig, branching ends, cycles, reverse-strand links (Huson et al. 2002).

### 1.4 Known Failure Modes / Pitfalls

1. Emitting `Math.Max(1, gapSize)` for non-positive estimates (1 N) is non-conforming — the GenBank/EMBL/DDBJ unknown-gap standard is 100 N. — NCBI AGP Spec v2.1
2. Placing a contig into more than one scaffold breaks the path model. — Pop et al. (2004) / Jackman et al. (2017)

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `Scaffold(IReadOnlyList<string> contigs, IReadOnlyList<(int,int,int)> links, char gapCharacter='N')` | SequenceAssembler | Delegate | forward links, weight 1 |
| `Scaffold(IReadOnlyList<string> contigs, IReadOnlyList<ScaffoldLink> links, char gapCharacter='N')` | SequenceAssembler | Canonical | oriented, weighted links; deep evidence-based testing |
| MCP `scaffold_contigs` | AlignmentTools | Wrapper | delegates |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Between two linked contigs the scaffold contains a contiguous run of exactly `gapSize` gap characters when `gapSize > 0` | Yes | Jackman et al. (2017) — run length = distance estimate |
| INV-2 | `gapSize = 0`, or a negative estimate whose overlap does not agree, emits exactly 100 gap characters; an agreeing overlap is collapsed | Yes | NCBI AGP Spec v2.1; ABySS MergeContigs |
| INV-3 | Scaffold length for a path = Σ(contig lengths) + Σ(emitted gap lengths) − Σ(collapsed overlaps) | Yes | Jackman et al. (2017) — concatenation + gaps |
| INV-4 | Each contig appears in exactly one scaffold (no contig used twice) | Yes | Pop et al. (2004) path-of-distinct-contigs model |
| INV-5 | The original contig substrings (reverse-complemented for `-`) are preserved verbatim and in path order | Yes | Jackman et al. (2017) — "sequences ... are concatenated" |
| INV-6 | ≤ 1 selected link per contig end; no cycle; links considered in descending weight, ties by input order | Yes | Huson et al. (2002) Alg. 5.1 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Chain two links positive gaps | `["ACGT","TTGG","CCAA"]`, links `[(0,1,3),(1,2,2)]` | `"ACGTNNNTTGGNNCCAA"`, count 1 | Jackman et al. (2017) |
| M2 | Exact gap count | `["AA","TT"]`, link `(0,1,5)` | `"AANNNNNTT"`, gap run length 5 | Jackman et al. (2017) (INV-1) |
| M3 | Negative gap → 100 N | `["AAAA","TTTT"]`, link `(0,1,-5)` | `"AAAA"+"N"×100+"TTTT"`, length 108 | NCBI AGP Spec v2.1 (INV-2) |
| M4 | Zero gap → 100 N | `["AAAA","TTTT"]`, link `(0,1,0)` | gap run length 100 | NCBI AGP Spec v2.1 (zero invalid → unknown default) |
| M5 | Custom gap character | `["AA","TT"]`, link `(0,1,2)`, gapChar `'X'` | `"AAXXTT"` | Jackman et al. (2017) — "run of the character" parameterized |
| M6 | No links → one scaffold per contig | `["AAA","CCC"]`, links `[]` | `["AAA","CCC"]`, count 2 | Pop et al. (2004) length-1 paths |
| M7 | Length invariant | `["ACGT","TTGG","CCAA"]`, links `[(0,1,3),(1,2,2)]` | length = 12 + 5 = 17 | Jackman et al. (2017) (INV-3) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Contig not reused | `["AA","TT"]`, links `[(0,1,1),(1,0,1)]` | one scaffold `"AANTT"`, count 1 (index 0 already used) | INV-4 |
| S2 | Out-of-range link ignored | `["AA","TT"]`, link `(0,5,2)` | `["AA","TT"]`, count 2 | malformed-link robustness |
| S3 | Self link ignored | `["AA","TT"]`, link `(0,0,2)` | `["AA","TT"]`, count 2 | contig2==contig1 skipped |
| S4 | Null contigs throws | `Scaffold(null, links)` | `ArgumentNullException` | sibling convention |
| S5 | Null links throws | `Scaffold(contigs, null)` | `ArgumentNullException` | sibling convention |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty contigs | `Scaffold([], links)` | empty result | trivial identity |

### 4.4 B10 review additions (oracle values, Evidence table)

| ID | Test Case | Expected | Evidence |
|----|-----------|----------|----------|
| R1 | `AA,TT`, `(1,0,2)` | `TTNNAA` | Huson Alg. 5.1 (F16) |
| R2 | `AA,CC,GG`, `(2,0,1),(0,1,1)` | `GGNAANCC` | F16 |
| R3 | `AA,CC,GG`, `(1,0,1),(1,2,1)` | `CCNAA`, `GG` | ≤ 1 link per end |
| R4 | weights 1 vs 5 | heavier link wins | descending weight |
| R5 | cycle `2→0→1→2` | `GGNAANCC` | P1 ≠ P2 |
| R6 | orientation `+-`, `-+`, `--` | `AACGNNNGCAA`, `GCAANNNAACG`, `TTGCNNNAACG` | F17 |
| R7 | oriented chain / inconsistent orientation | `ACCANNAACCNNNNCATG`; `ACCANNAACC`,`CATG` | F17 |
| R8 | overlap agrees / with N / on RC / then gap | `ACGTACGGA` ×2, `GGATCCAAA`, `AAACGTTTNNNGGGG` | F18 |
| R9 | overlap disagrees | 100 N | AGP |
| R10 | ordering by lowest index | `AAA`, `GGGNNNNCCC` | documented order |
| R11 | weight 0 | `ArgumentOutOfRangeException` | weight = mate pairs |
| R12 | fuzz: oriented overlapping genome pieces | genome or its RC (1000 trials) | F17/F18 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior test file for `SequenceAssembler.Scaffold` exists. Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `Scaffold` — none found. Sibling: `SequenceAssembler_MergeContigs_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| M7 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| S5 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_Scaffold_Tests.cs` — all cases above.
- **Remove:** (none — no prior tests existed)

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceAssembler_Scaffold_Tests.cs` | Canonical | 13 |

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
| 8 | S1 | ❌ Missing | Implemented | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented | ✅ Done |
| 11 | S4 | ❌ Missing | Implemented | ✅ Done |
| 12 | S5 | ❌ Missing | Implemented | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `Scaffold_TwoPositiveGapLinks_ChainsWithExactNRuns` |
| M2 | ✅ Covered | `Scaffold_PositiveGap_EmitsExactlyThatManyGapChars` |
| M3 | ✅ Covered | `Scaffold_NegativeGap_EmitsHundredGapChars` |
| M4 | ✅ Covered | `Scaffold_ZeroGap_EmitsHundredGapChars` |
| M5 | ✅ Covered | `Scaffold_CustomGapCharacter_UsesItVerbatim` |
| M6 | ✅ Covered | `Scaffold_NoLinks_OneScaffoldPerContig` |
| M7 | ✅ Covered | `Scaffold_FollowedPath_LengthEqualsContigsPlusGaps` |
| S1 | ✅ Covered | `Scaffold_LinkToAlreadyPlacedContig_ContigNotReused` |
| S2 | ✅ Covered | `Scaffold_OutOfRangeLink_Ignored` |
| S3 | ✅ Covered | `Scaffold_SelfLink_Ignored` |
| S4 | ✅ Covered | `Scaffold_NullContigs_ThrowsArgumentNullException` |
| S5 | ✅ Covered | `Scaffold_NullLinks_ThrowsArgumentNullException` |
| C1 | ✅ Covered | `Scaffold_EmptyContigs_ReturnsEmpty` |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | A zero estimate, or a negative estimate whose overlap is not confirmed, uses the AGP unknown-gap length (100 N); ABySS/SSPACE would emit a single `n` (AGP chosen as the submission format) | M3, M4, INV-2 |

---

## 7. Open Questions / Decisions

1. Decision: This unit consumes a supplied gap estimate; the upstream maximum-likelihood distance estimator (Jackman et al. 2017; Sahlin et al. 2012) is out of scope.
2. Decision (B10): Huson's interleaving ("zipper") merge and happy/unhappy test need per-link σ, not in the input; such links are conflicts. Bundling, SSPACE ratio filter and ABySS gapped-overlap fallback not implemented (see Scaffolding.md §5.3).
