# Test Specification: ASSEMBLY-COMPARE-001

**Test Unit ID:** ASSEMBLY-COMPARE-001
**Area:** Assembly
**Algorithm:** Assembly comparison (Mash k-mer sets) and syntenic blocks (nucmer MUMs + mgaps)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect (B10 sweep U4)
**Last Updated:** 2026-11-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | Reference | Accessed |
|---|--------|---------------|-----------|----------|
| 1 | Mash (Ondov et al. 2016): J, D = −(1/k)ln(2J/(1+J)), canonical ACGT-only k-mers, boundary rules (via `KmerAnalyzer`) | 1 | doi:10.1186/s13059-016-0997-x | – |
| 2 | Containment index (Koslicki & Zabeti 2019; sourmash) via `KmerAnalyzer.ContainmentIndex` | 1 | – | – |
| 3 | MUMmer 4 `mgaps` (`src/tigr/mgaps.cc`, `include/mummer/mgaps.hh`) | 3 (reference implementation) | raw.githubusercontent.com/mummer4/mummer/master | 2026-11-25 |
| 4 | nucmer usage (defaults −l 20, −c 65, −g 90, −d 0.12, −D 5, `--mumreference`); MUMmer 3 manual | 2 | mummer.sourceforge.net/manual (snippets) | 2026-11-25 |
| 5 | Delcher et al. 2002 (NAR 30:2478), Marçais et al. 2018 (dnadiff) — model only | 1 | – | – |

### 1.2 Key Evidence Points

1. Mash compares exact (or sketched) canonical k-mer sets; windows with a non-ACGT symbol are not k-mers; identical sets → D 0, disjoint → D 1.
2. nucmer anchors are maximal unique matches (unique in the reference by default), found on both query strands; anchors are clustered by mgaps (Filter_Matches, union-find on gap ≤ maxGap and diagonal difference ≤ max(fixed, factor·gap), best-chain extraction) and kept when the summed anchor length ≥ the minimum cluster score.
3. dnadiff aligned bases / identity / breakpoints come from alignments; they are not defined for k-mer sets.

### 1.3 Documented Corner Cases

- Both sets empty: J undefined → 0, Mash D 0 (shared = union).
- Sequence shorter than k contributes no k-mer; k-mers never span contigs.
- A repeat in the reference is not an anchor (`--mumreference`); N and IUPAC codes never match.
- Blocks bridge short inversions when diagonals agree on both sides and the gap ≤ maxGap.

### 1.4 Known Failure Modes / Pitfalls

1. Averaging two containment fractions is not an identity (old `SequenceIdentity`); breakpoint counts were constants.
2. Spacing anchors every k bases and splitting on a sign flip is not a published clustering.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CompareAssemblies(asm1, asm2, name1, name2, kmerSize=21)` | GenomeAssemblyAnalyzer | Canonical (Mash model) | `KmerAnalyzer.DistinctKmers(Canonical, AcgtOnly)`, `MashDistanceFromJaccard` |
| `FindSyntenicBlocks(asm1, asm2, minBlockSize=1000, kmerSize=21, maxGap=90, diagonalFactor=0.12, diagonalDifference=5, uniqueInBoth=false)` | GenomeAssemblyAnalyzer | Canonical (nucmer front end) | suffix-tree MUMs + `MgapsClusterer` port |
| MCP `compare_assemblies`, `find_syntenic_blocks_assemblies` | ChromosomeTools | Wrapper | delegate |

---

## 3. Invariants

| ID | Invariant | Evidence |
|----|-----------|----------|
| INV-1 | 0 ≤ J, C1, C2 ≤ 1; shared ≤ min(|K1|,|K2|); C1·|K1| = C2·|K2| = shared | Source 1, 2 |
| INV-2 | CompareAssemblies is strand-insensitive (RC of an assembly gives J = 1) | Source 1 |
| INV-3 | Block coordinates are 0-based inclusive; Start ≤ End; IsInverted ⇔ match with the RC query strand | Source 3 |
| INV-4 | Result equals the oracle (brute-force MUM definition + mgaps port) | Sources 3, 4 |

---

## 4. Test Cases (`GenomeAssemblyAnalyzer_AssemblyCompare_Tests`, MCP `CompareAssembliesTests`, `FindSyntenicBlocksAssembliesTests`, MutationKillers3/4/6, `GenomeAssemblyAnalyzerTests`)

| ID | Test Case | Expected | Test |
|----|-----------|----------|------|
| M1 | 12 oracle fixtures (600/600 offline) | counts, J, D, containment, every block | `CompareAssemblies_AndSyntenicBlocks_MatchOracle` |
| M2 | ACGTACGT vs ACGTAAAA, k 4 | shared 2, J 1/3, D ln2/4 | `CompareAssemblies_MashFormula_PartialOverlap`, MCP `…PartialOverlap_ReportsMashValues` |
| M3 | Strand / pooling / non-ACGT / empty / validation | see test names | `CompareAssemblies_IsStrandInsensitive`, `…ContigsArePooled…`, `…NonAcgtWindowsAreSkipped`, `…EmptySets…`, `…Validation` |
| M4 | Planted inversion, rearrangement, unrelated | 3 / 2 / 0 blocks with strands | `FindSyntenicBlocks_PlantedInversion…`, `…Rearrangement…`, `…Unrelated…` |
| M5 | Reference repeat, N, gap > maxGap, minBlockSize | oracle values | `…ReferenceRepeatIsNotAnAnchor`, `…SingleMismatchSymbol…`, `…GapBeyondMaxGap…`, `…MinBlockSize…` |
| M6 | Identical sequences | one block over the whole sequence | MutationKillers3 `…CollinearBlock…`, `GenomeAssemblyAnalyzerTests.FindSyntenicBlocks_IdenticalSequences_FindsBlock`, MCP identical |
