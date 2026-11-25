# Evidence Artifact: ASSEMBLY-COMPARE-001

**Test Unit ID:** ASSEMBLY-COMPARE-001
**Algorithm:** Assembly comparison (`GenomeAssemblyAnalyzer.CompareAssemblies`, `AssemblyComparison`) and syntenic blocks (`GenomeAssemblyAnalyzer.FindSyntenicBlocks`; MCP `compare_assemblies`, `find_syntenic_blocks_assemblies`)
**Date Collected:** 2026-11-25 (B10 sweep U4)

---

## What the old code computed (and why it was replaced)

- `CompareAssemblies`: forward-strand distinct k-mers (only `N` windows skipped), `AlignedFraction1/2` = shared / distinct per assembly, `SequenceIdentity` = their arithmetic mean (no source), `Breakpoints/Inversions/Translocations` hard-coded 0. "Aligned" and "identity" are alignment concepts (dnadiff AlignedBases / AvgIdentity); a k-mer set ratio is neither. The sourced alignment-free quantities are the Mash Jaccard / distance and the containment index.
- `FindSyntenicBlocks`: every `kmerSize`-th k-mer of assembly 1 looked up in an index of assembly 2, anchors per target sequence sorted by position, a new block when a gap exceeds `minBlockSize` or the sign of the target step flips; block `End` = start of the last anchor (k-mer length missing); inversion = "second anchor moved backwards". No source for the grouping rule; repeat k-mers produced all-pairs anchors. The published model of "anchors → colinear blocks" for assemblies is MUMmer (MUMs + mgaps clustering).
- DUP_MAP §20: the Chromosome synteny functions are conceptual duplicates of gene-order synteny (`ComparativeGenomics.FindSyntenicBlocks`, MCScanX on orthologous gene pairs). They are NOT duplicates of this method: input here is raw sequence, anchors are exact unique matches. Left in place; the anchor finder is the canonical suffix-tree MUM finder (`ISuffixTreeAnalysis.FindMaximalUniqueMatches`, MUMmer 3 `-mum/-mumreference`) and the k-mer sets are `KmerAnalyzer.DistinctKmers(Canonical)`.

## Sources actually opened

### Mash (Ondov et al. 2016, Genome Biol 17:132) — via the project's canonical `KmerAnalyzer` (docs opened in MASH/KMER units; Mash `Sketch.cpp`/`CommandDistance.cpp` rules transcribed in `KmerAnalyzer.MashDistanceFromJaccard`)
J = |A∩B| / |A∪B| (eq. 1); D = −(1/k)·ln(2J/(1+J)) (eq. 4); shared = union → 0, no shared → 1; canonical k-mers, windows with non-ACGT skipped (Mash `addMinHashes`). A multi-FASTA is one genome for `mash sketch`. The paper itself was not re-opened in this session.
Containment |A∩B|/|A| (Koslicki & Zabeti 2019; sourmash `--containment`) via `KmerAnalyzer.ContainmentIndex` documentation.

### MUMmer 4 source (mummer4/mummer master), fetched with curl, neutral User-Agent
- `src/tigr/mgaps.cc` (`Filter_Matches`, printing) and `include/mummer/mgaps.hh` (`Cluster_each`, `Process_Cluster`, `ClusterMatches`, `UnionFind`, `Match_t`): read in full and ported literally (`MgapsClusterer`).
- `src/tigr/mgaps_main.cc`: standalone defaults (fixed separation 5, max separation 1000, min score 200, factor 0.05, `Use_Extents` false). nucmer's own defaults come from the nucmer usage text (WebSearch snippets of the MUMmer 3 manual and Debian `nucmer(1)`): `-l` 20, `-c` 65, `-g` 90, `--diagfactor` 0.12, `-D` 5, default `--mumreference`. `scripts/nucmer.pl` / the nucmer driver that maps these options onto `mgaps` could not be fetched, so the mapping (−c → min output score, −g → max separation, −d → separation factor, −D → fixed separation) and nucmer's use of sum-of-lengths instead of extents are **not verified from source**.
- `scripts/dnadiff.pl`, `src/tigr/delta-filter.cc`, `show-diff.cc`, `postnuc.cc` were fetched; only used to confirm that alignment-level statistics (aligned bases, AvgIdentity, breakpoints/inversions/relocations) are computed from nucmer **delta alignments** (after extension), which this library does not produce.
- Compiling the real `mgaps` as a reference binary was blocked by the sandbox (same as TRF in U3), so the oracle is a literal Python port; MUMs are checked against an independent brute-force implementation of the MUM definition (Delcher 1999 / Kurtz 2004): maximal exact match of length ≥ l whose string occurs once in the reference (and once in the query for `--mum`).
- Delcher et al. 2002 NAR 30:2478 and Marçais et al. 2018 PLoS Comput Biol 14:e1005944: cited for the model (MUM anchors → clustering → alignment; dnadiff reports); not opened.

## Decisions

1. N rule (B06 request): Mash/Jellyfish skip every window containing a non-ACGT symbol → `AcgtOnly`+`Canonical` counting, no private N-only loop. nucmer/mummer never match non-ACGT (`mummer -n`): such symbols are replaced by distinct non-matching placeholders on the two sides.
2. Strands: Mash is strand-collapsed (canonical); nucmer clusters forward and reverse-complement query matches separately.
3. Records/columns removed: `AlignedFraction*`, `Breakpoints`, `Inversions`, `Translocations`, `SequenceIdentity` (unsourced). New: `KmerSize`, `DistinctKmers1/2`, `SharedKmers`, `Containment1/2`, `JaccardIndex`, `MashDistance`.
4. Coordinates of blocks changed to 0-based **inclusive** ends of the first-anchor-start … last-anchor-end span (old End was the last anchor's start).

## Reference numbers (`tools/review/oracles/B10/compare_oracle.py`, `python3`)

- 600 random + planted cases (seed 7; 1–2 contigs of 60–220 bp; mutated copies at 0/2/8 %, reverse-complemented + inserted contigs, rotated halves, unrelated + repeats + N/n; k ∈ {4,6,8,11}, minMatch ∈ {8,10,12,15}, minBlock ∈ {15,30,65}): C# `CompareAssemblies` (counts, J, D, containment to 1e-12) and `FindSyntenicBlocks` (every coordinate and strand; 1037 blocks, 222 inverted) identical on **600/600**. 12 locked as `OracleCases` fixtures.
- Hand values: ACGTACGT vs ACGTAAAA, k 4 → canonical sets {ACGT,CGTA,GTAC} / {ACGT,CGTA,GTAA,TAAA,AAAA}, shared 2, J 1/3, D = ln2/4 = 0.173287; A^100 vs C^100, k 21 → J 0, D 1; planted inversion A+RC(B)+C+D (150-bp random segments, l 20, c 65) → [0,149]+, [150,299]−, [300,599]+; segment swap C+D+A+B → [0,299]→[300,599], [300,599]→[0,299]; a 614-bp insertion beyond maxGap 90 splits a block (gap rule).
- Observed faithful-nucmer quirk: an inversion shorter than maxGap (80 bp ≤ 90) with equal diagonals on both sides is bridged by mgaps into one forward block while still reported as its own inverted block.

## R1-WP1 addendum (2026-11-25): ASCII fold (A15-7, A11-6)

Sources opened: `Infrastructure/SequenceCase.cs` remarks (ASCII-only fold; B01 F21), `SequenceExtensions.IsValidDna` (IUPAC-IUB unambiguous DNA = Biopython `unambiguous_dna_letters`). `IsAcgt(span)` vs `IsValidDna(span)` differ only for lower-case input, which `FindRepetitiveRegions` has already upper-cased. Non-ACGT (incl. non-ASCII) symbols make a WindowMasker unit invalid exactly as `N` does (tests compare the result with the `N` variant). Note: `CalculateLocalQuality` passes the window to `SequenceComplexity.CalculateLinguisticComplexity`, which folds with `ToUpperInvariant` itself (B04 file, left as is; cross-batch note in the WP report).

