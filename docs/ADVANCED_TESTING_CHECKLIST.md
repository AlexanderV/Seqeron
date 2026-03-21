# Advanced Testing Techniques — Analysis & Checklists

**Date:** 2026-03-19  
**Scope:** 79 completed test units from `ALGORITHMS_CHECKLIST_V2.md` Processing Registry  
**Library:** Seqeron.Genomics

---

## Part 1: Effectiveness Analysis

### Summary Matrix

| # | Technique | Applicability | Existing Coverage | Effort | Priority |
|---|-----------|:---:|:---:|:---:|:---:|
| 1 | **Property-Based Testing** | ★★★★★ | 22/79 files | Med | **P0** |
| 2 | **Metamorphic Testing** | ★★★★★ | 7/79 units | Med | **P0** |
| 3 | **Mutation Testing** | ★★★★☆ | 2/79 units | Low | **P1** |
| 4 | **Algebraic Testing** | ★★★★☆ | 0/79 units | Med | **P1** |
| 5 | **Snapshot / Approval Testing** | ★★★★☆ | ~40/79 units | Low | **P1** |
| 6 | **Architecture Testing** | ★★★☆☆ | 1 file, 5 rules | Low | **P2** |
| 7 | **Differential Testing** | ★★★☆☆ | 0/79 (exists for SuffixTree only) | High | **P2** |
| 8 | **Fuzzing** | ★★★☆☆ | 0/79 (exists for SuffixTree only) | Med | **P2** |
| 9 | **Characterization Testing** | ★★☆☆☆ | 0/79 explicit | Low | **P3** |
| 10 | **Combinatorial / Pairwise** | ★★☆☆☆ | 0/79 explicit | Med | **P3** |

---

### 1. Property-Based Testing (FsCheck) — ★★★★★

**Why effective:** Genomic algorithms have rich mathematical invariants (range bounds, symmetry, idempotence, complement preservation, metric axioms). FsCheck generates hundreds of random inputs per test, catching edge cases hand-written tests miss.

**Current state:** 22 property files in `Properties/` covering: GcContent, GcSkew, Sequence, EditDistance, Hamming, FASTA, FileIO, Alignment, Codon, CRISPR, K-mer, MiRNA, PatternMatching, Phylogenetic, PopGen, PrimerProbe, ProteinMotif, RepeatFinder, Restriction, RnaStructure, SequenceComposition, Splicing.

**Gap:** Property files exist for most areas but many test units within those areas lack specific property tests. Need: systematic invariant identification per algorithm.

**Recommended for ALL 79 completed units** — every algorithm has at least one expressible invariant.

---

### 2. Metamorphic Testing — ★★★★★

**Why effective:** The "oracle problem" is central to bioinformatics — for many algorithms there is no single correct answer to compare against. Metamorphic relations (MR) bypass this by relating outputs of transformed inputs. Ideal for search, alignment, scoring, and statistical algorithms.

**Current state:** `MetamorphicTests.cs` covers 7 units: PAT-IUPAC-001, PAT-PWM-001, PAT-APPROX-002, REP-STR-001, REP-TANDEM-001, REP-INV-001, REP-DIRECT-001 (18+ MR relations).

**Gap:** 72 completed units have no metamorphic tests. Many have strong natural MRs (alignment symmetry, GC preservation under complement, diversity monotonicity, etc.).

**Recommended for ~65 units** with clear metamorphic relations.

---

### 3. Mutation Testing (Stryker) — ★★★★☆

**Why effective:** Measures test suite effectiveness by verifying tests actually kill code mutations. Stryker already configured for this project. Cost: CPU-intensive but requires zero new test code — surfaces gaps in existing tests.

**Current state:** `MutationKillerTests.cs` targets MotifFinder.cs and RepeatFinder.cs. `stryker-config.json` exists. StrykerOutput/ has run results.

**Gap:** Only run against 2 source files. Should systematically run per-module and create targeted killer tests for survived mutants.

**Recommended for ALL 79 units** — config-only for most, test writes only for survivors.

---

### 4. Algebraic Testing — ★★★★☆

**Why effective:** Verifies algebraic laws: identity, associativity, commutativity, inverse, distributivity, idempotence. Genomic algorithms have many: `complement(complement(x)) = x`, `reverseComplement(reverseComplement(x)) = x`, alignment score symmetry, `parse(serialize(x)) = x`. More structured than general property tests.

**Current state:** Some algebraic properties exist implicitly in Property files (e.g., complement involution in SequenceProperties, round-trip in FastaRoundTripProperties) but not labeled or systematic.

**Gap:** No dedicated algebraic test suite. Many algorithms have untested algebraic laws.

**Recommended for ~45 units** with clear algebraic laws.

---

### 5. Snapshot / Approval Testing (Verify) — ★★★★☆

**Why effective:** Complex outputs (alignment matrices, phylogenetic trees, annotation records, parse results) are hard to assert field-by-field. Snapshot tests serialize the full output and catch any unintentional change. Already mature in this project.

**Current state:** ~20 snapshot test files in `Snapshots/` covering: Alignment, Annotation, Chromosome, Codon, CRISPR, Disorder, Epigenetics, FileIO, Metagenomics, MiRNA, MolTools, PatternMatching, Phylogenetic, Population, PrimerProbe, ProteinMotif, Repeat, Restriction, RNA, Splicing.

**Gap:** Some recently added units (ONCO-IMMUNE-001) lack snapshots. Some snapshot files may be stale.

**Recommended for ~30 remaining units** with complex structured outputs.

---

### 6. Architecture Testing (ArchUnitNET) — ★★★☆☆

**Why effective:** Prevents architectural drift as the library grows. Enforces layer boundaries (Core → Analysis → IO), naming conventions, and structural rules at IL level.

**Current state:** `Architecture/ArchitectureTests.cs` with 5 rules (Core !→ Analysis, Core !→ IO, Core !→ Alignment, IO !→ Analysis, static analyzers).

**Gap:** Rules don't cover newer modules (Metagenomics, Phylogenetics, Population, Chromosome, MolTools, Oncology). No naming convention rules for test classes.

**Recommended:** Expand to all 14 source projects, add module-specific rules. Not per-algorithm but per-module.

---

### 7. Differential Testing — ★★★☆☆

**Why effective:** Compares two independent implementations of the same algorithm. Catches subtle implementation bugs. In bioinformatics, well-known reference implementations exist (Biopython, EMBOSS, etc.).

**Current state:** `SafeVsUnsafeDifferentialTests` exists for SuffixTree (safe vs. unsafe memory paths). Zero for Genomics.

**Gap:** High effort — requires reference implementations or alternative algorithm variants. Best ROI for algorithms with multiple strategies (e.g., UPGMA vs. NJ tree construction, Nussinov vs. MFE RNA folding, NW global vs. recoded).

**Recommended for ~15 units** where dual implementations exist or can be trivially written.

---

### 8. Fuzzing — ★★★☆☆

**Why effective:** Discovers crashes, hangs, and security issues in parsers and input-processing code. Highly effective for file format parsers (FASTA, FASTQ, BED, VCF, GFF, GenBank, EMBL) and sequence validation.

**Current state:** `SuffixTreeFuzzTests` exists (header corruption). Zero for Genomics.

**Gap:** All 7 parser units and sequence validation are prime fuzz targets. Random byte sequences, malformed headers, truncated records.

**Recommended for ~15 units** — primarily FileIO parsers and input validation.

---

### 9. Characterization Testing — ★★☆☆☆

**Why effective:** Captures current behavior as a safety net before refactoring. Less useful for new code. Snapshot tests already serve this role partially.

**Current state:** Snapshot tests already function as characterization tests. No separate labeled suite.

**Gap:** Low — snapshot tests cover this need. Useful only when planning major refactors of existing algorithm internals.

**Recommended for ~10 units** as a pre-refactoring safety net (only when refactoring is planned).

---

### 10. Combinatorial / Pairwise Testing — ★★☆☆☆

**Why effective:** Reduces exponential parameter combinations to a polynomial set while maintaining coverage of pairwise interactions. Useful for algorithms with many configuration parameters.

**Current state:** None explicit. Some tests use `[TestCase]` with multiple parameter combinations but not systematic pairwise.

**Gap:** Relevant for algorithms with multiple configuration knobs: primer design (Tm method + length range + GC range + salt), codon optimization (organism + method + constraints), CRISPR (PAM type + guide length + scoring), restriction enzymes (enzyme set + circular/linear).

**Recommended for ~12 units** with ≥3 independent configuration parameters.

---



## Part 2: Checklists by Technique (separate files)

Each technique has a dedicated checklist covering ALL 86 completed algorithms.

| # | Priority | Technique | Checklist File |
|---|----------|-----------|---------------|
| 1 | P0 | Property-Based Testing (FsCheck) | [01_PROPERTY_BASED_TESTING.md](checklists/01_PROPERTY_BASED_TESTING.md) |
| 2 | P0 | Metamorphic Testing | [02_METAMORPHIC_TESTING.md](checklists/02_METAMORPHIC_TESTING.md) |
| 3 | P1 | Mutation Testing (Stryker) | [03_MUTATION_TESTING.md](checklists/03_MUTATION_TESTING.md) |
| 4 | P1 | Algebraic Testing | [04_ALGEBRAIC_TESTING.md](checklists/04_ALGEBRAIC_TESTING.md) |
| 5 | P1 | Snapshot / Approval Testing (Verify) | [05_SNAPSHOT_TESTING.md](checklists/05_SNAPSHOT_TESTING.md) |
| 6 | P2 | Architecture Testing (ArchUnitNET) | [06_ARCHITECTURE_TESTING.md](checklists/06_ARCHITECTURE_TESTING.md) |
| 7 | P2 | Differential Testing | [07_DIFFERENTIAL_TESTING.md](checklists/07_DIFFERENTIAL_TESTING.md) |
| 8 | P2 | Fuzzing | [08_FUZZING.md](checklists/08_FUZZING.md) |
| 9 | P3 | Characterization Testing | [09_CHARACTERIZATION_TESTING.md](checklists/09_CHARACTERIZATION_TESTING.md) |
| 10 | P3 | Combinatorial / Pairwise Testing | [10_COMBINATORIAL_TESTING.md](checklists/10_COMBINATORIAL_TESTING.md) |

### Execution Roadmap

1. **P0 — Property-Based:** Fill gaps in existing 22 property files; create 7 new files for uncovered areas
2. **P0 — Metamorphic:** Create area-specific MR test classes (one per area, ~15 new files)
3. **P1 — Mutation:** Run Stryker per-module, write killers for survivors
4. **P1 — Algebraic:** Add algebraic law tests (can combine with property files)
5. **P1 — Snapshot:** Add missing 31 snapshot tests to existing + new snapshot files
6. **P2 — Architecture:** Expand ArchitectureTests.cs with 14 new module rules
7. **P2 — Fuzzing:** Create fuzz test files for 7 parsers + boundary tests for all areas
8. **P2 — Differential:** Create differential test pairs for all 86 algorithms
9. **P3 — Characterization:** Create on-demand before refactoring
10. **P3 — Combinatorial:** Create pairwise suites for 15 multi-parameter algorithms
