# bio-chromosome — fuller recipes, parameters, gotchas

Progressive-disclosure companion to `SKILL.md`. All tool names / Method IDs verified against
`docs/mcp/tools/chromosome/*.md`. Rigor (parse-with-a-tool, provenance, envelope, cross-check,
alpha caveat) is delegated to `bio-rigor` — not restated here.

---

## Parameter cheat-sheet (defaults from tool docs)

| Tool | Key params (default) |
|---|---|
| `analyze_centromere` | `windowSize=100000`, `minAlphaSatelliteContent=0.3` (extends while repeat ≥ 70% of threshold) |
| `analyze_telomeres` | `telomereRepeat=TTAGGG` (5′ matched vs reverse-complement `CCCTAA`), `searchLength=10000`, `minTelomereLength=500`, `criticalLength=3000`; per-window similarity ≥ 0.7 |
| `find_heterochromatin_regions` | `windowSize=100000` (step = windowSize/2), `minRepeatContent=0.5`; class by midpoint: <5%/>95% Telomeric, 45–55% Centromeric, else Constitutive |
| `predict_g_bands` | `bandSize`, `darkBandGcThreshold`, `lightBandGcThreshold` — GC below dark→gpos100, below light→gpos50, else gneg |
| `analyze_karyotype` | `expectedPloidyLevel=2` (>0); strips trailing `_N` copy suffix (`chr1_1`,`chr1_2`→`chr1`); sex chromosomes reported separately, not grouped |
| `detect_aneuploidy` | `binSize=1000000`; `logRatio=log2(meanDepth/medianDepth)`, `copyNumber=round(2^logRatio×2)` clamped [0,10] |
| `detect_ploidy` | `ploidy=round(median/expectedDiploidDepth × 2)` clamped [1,8] |
| `find_synteny_blocks` | `minGenes=3` (>0), `maxGap=10` (Mb) |
| `find_syntenic_blocks_assemblies` | `minBlockSize` (k-mer anchors) |
| `assembly_statistics` | in: `[{id,sequence}]`; length reported with and without N gaps; GC over non-N bases |
| `nx_statistics` / `nx_curve` | **descending-sorted lengths + explicit totalLength**; `threshold` in [0,100] (50→N50); inclusive "≥ x%" (QUAST) |
| `find_gaps` | `minGapLength=1` (>0); `[start,end]` inclusive; class <10 Short / <100 Medium / <1000 Long / else Scaffold |
| `extract_contigs` | `minContigLength`; emits `{id}_contig{n}` over gap-free runs |
| `compare_assemblies` | `kmerSize=21`; structural counts (breakpoints/inversions/translocations) reported as 0 by this k-mer path |

---

## Coordinate & unit notes

- **Gaps / bands / heterochromatin regions**: 1-based **inclusive** `[start,end]`; gap length =
  `end - start + 1`. Confirm per tool doc — do not assume 0-based here.
- **Synteny**: positions come from the ortholog-pair inputs you supply; `maxGap` is in **Mb**, not bp.
- **Depth binning**: bin index = `position / binSize`; copy number is depth-relative — pass a
  **genome-wide median** (`detect_aneuploidy`) or `expectedDiploidDepth` (`detect_ploidy`), never a
  single-region mean, or the ratios skew.
- **Telomeres**: the 3′ end is matched against the repeat unit; the 5′ end against its
  reverse-complement — pass the biological repeat (`TTAGGG`), the tool handles both ends.

---

## Cross-check pairings (independent paths)

- **Centromere location**: `analyze_centromere.type` (arm-ratio / Levan) vs
  `classify_chromosome_by_arm_ratio(arm_ratio(...))` — should agree; and the Centromeric region
  from `find_heterochromatin_regions` (repeat-content path) should overlap the centromere span.
- **Contiguity**: `assembly_statistics.n50` vs `au_n` (threshold-free) vs `nx_curve` — same
  assembly should trend consistently; auN ≥ N50 typically for skewed length distributions.
- **Aneuploidy**: descriptor path `analyze_karyotype` vs depth path
  `identify_whole_chromosome_aneuploidy` — a trisomy should surface in both when both inputs exist.
- **Rearrangements**: `detect_rearrangements` (gene-synteny) vs `find_syntenic_blocks_assemblies`
  (k-mer synteny orientation flags) — inversions should be flagged by both on the same event.
- **Completeness**: `assess_completeness` (marker genes) vs `estimate_completeness_from_kmers`
  (spectrum) — two orthogonal completeness estimates.

---

## Gotchas

- **GC-skew / replication origin are NOT on this server.** They are `GcSkewCalculator.*` on the
  **analysis** server: `gc_skew`, `cumulative_gc_skew`, `windowed_gc_skew`,
  `predict_replication_origin`. Route through `bio-annotation`. "Chromosome composition" *here* means
  G-bands / local GC-complexity / heterochromatin, not skew.
- **`nx_statistics` needs pre-sorted descending lengths + totalLength** — it does not sort raw
  `{id,sequence}` for you. `assembly_statistics` and `nx_curve` are the higher-level entry points.
- **`compare_assemblies` reports structural counts as 0** — it is a k-mer similarity comparison, not
  a structural aligner. For inversions/translocations use `find_syntenic_blocks_assemblies` +
  `detect_rearrangements`.
- **Descriptor karyotype ≠ depth aneuploidy.** `analyze_karyotype` counts *named* chromosome copies
  (needs `chr1_1`,`chr1_2` style descriptors); it will not detect a partial/mosaic aneuploidy that
  only shows up in read depth — use the depth tools for that.
- **Variant-level SV/CNV** (single-event breakpoints, per-locus CNV calls, effect/pathogenicity) →
  `bio-annotation`. This skill is chromosome-scale (whole-arm/whole-chromosome, synteny-scale).

---

## Envelope — STOP rule

- **CHROM-CENT-001** — `AssignSuprachromosomalFamily` (output `Sf1OrSf2Dimeric`) is **Permissive-only**
  (blocked in Strict & Moderate; default Moderate). Distinguishing **SF1 vs SF2** dimeric HORs needs
  an SF-resolved consensus-monomer reference the library does not ship (non-redistributable). If a
  task asks to *assign a suprachromosomal family / separate SF1 from SF2*, **STOP** and report the
  envelope; a caller holding an SF-resolved reference must pass it explicitly. The general centromere
  tools (`analyze_centromere`, `arm_ratio`, `classify_chromosome_by_arm_ratio`,
  `find_heterochromatin_regions`) are **not** guarded.
- Source of truth: `docs/Validation/LIMITATIONS.md` (LimitationPolicy Strict<Moderate<Permissive).

---

## C# API notes

- Both backing classes are in `Seqeron.Genomics.Chromosome`
  (`ChromosomeAnalyzer`, `GenomeAssemblyAnalyzer`). Method IDs in `tool-map.md` map 1:1 to MCP tools.
- Default LimitationPolicy is **Moderate**; the guarded SF assignment above requires a **Permissive**
  bootstrap. General chromosome/assembly tools run under the default.
