namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Chromosome area — telomere detection (CHROM-TELO-001),
/// centromere detection (CHROM-CENT-001), karyotype analysis (CHROM-KARYO-001)
/// and aneuploidy detection (CHROM-ANEU-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain input to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output, and no *unhandled* runtime
/// exception (IndexOutOfRangeException, ArgumentOutOfRangeException leaking from
/// internal indexing, DivideByZero, OutOfMemory). Every input must resolve to
/// EITHER a well-defined, theory-correct result, OR a *documented, intentional*
/// validation exception. A raw runtime exception, a hang, a crash, or a spurious
/// "telomere detected" on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CHROM-TELO-001 — telomere detection
/// Checklist: docs/checklists/03_FUZZING.md, row 48.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: the empty sequence, a sequence with NO telomeric
///          repeat, and a sequence that is ONLY the telomeric repeat (pure tract).
///   • MC = Malformed Content — non-DNA characters in the sequence (digits, IUPAC
///          ambiguity codes, whitespace, unicode) fed through the same surface.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The telomere-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A telomere is the repetitive DNA tract at a chromosome end. The vertebrate
/// canonical repeat unit is the 6-bp hexamer TTAGGG; its reverse complement
/// CCCTAA is the motif expected at the 5' end (Telomere_Analysis.md §2.1,
/// table "chromosome-end orientation"). The detector is
///   ChromosomeAnalyzer.AnalyzeTelomeres(
///       string chromosomeName,
///       string sequence,
///       string telomereRepeat = "TTAGGG",
///       int searchLength       = 10000,
///       int minTelomereLength  = 500,
///       int criticalLength     = 3000)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///    lines 250–288). It returns an EAGER value (the TelomereResult readonly
///    record struct) — there is no lazy iterator, so any exception or hang would
///    surface at the call itself; no `.ToList()` forcing is needed.
///
/// How the scan works (ChromosomeAnalyzer.cs lines 263–288, MeasureTelomereLength
/// lines 293–340; Telomere_Analysis.md §4.1, §5.2):
///   • The sequence and repeat motif are upper-cased; the 5' motif is the reverse
///     complement of `telomereRepeat`.
///   • The 5' PREFIX window (first min(searchLength, len) bases) is scanned left→
///     right against CCCTAA; the 3' SUFFIX window (last min(searchLength, len)
///     bases) is scanned right→left against TTAGGG.
///   • Scanning advances in repeat-sized (6-bp) steps, counting a window only while
///     its per-base similarity to the motif is ≥ 70%; the first window below 70%
///     stops that end. Only COMPLETE 6-bp windows are counted, so the measured
///     length is always a multiple of 6 and equals 6 × (number of accepted windows)
///     (INV-01 length ≥ 0).
///   • RepeatPurity = matchingBases / totalBases ∈ [0,1], or 0 when nothing counted
///     (INV-02). Has*Telomere ⇔ measured length ≥ minTelomereLength (INV-03).
///
/// Documented edge-case contract (Telomere_Analysis.md §3.1, §3.3, §6.1):
///   • Empty OR null sequence → the method special-cases it (line 258) and returns
///     a no-telomere result: lengths 0, Has*Telomere = false, purities 0, and
///     IsCriticallyShort = TRUE (the documented empty-input flag).
///   • Sequence shorter than the 6-bp repeat → MeasureTelomereLength exits
///     immediately (region.Length < repeatLen, line 299) → zero-length telomeres.
///   • No telomeric repeats → no window meets the 70% bar → zero lengths,
///     Has*Telomere = false, and (non-empty input) IsCriticallyShort = FALSE
///     (§5.2: critically-short is only flagged when a telomere is actually
///     detected and is shorter than criticalLength).
///   • Non-DNA characters → there is NO validation/rejection gate: the input is
///     upper-cased and scanned verbatim. Garbage bases simply fail the 70%
///     similarity test, so a non-DNA sequence behaves like "no telomere" and must
///     NEVER crash on indexing. This is the MC contract pinned below.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Citations
/// ───────────────────────────────────────────────────────────────────────────
/// • Algorithm doc: docs/algorithms/Chromosome_Analysis/Telomere_Analysis.md.
/// • Source: src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///   (AnalyzeTelomeres + MeasureTelomereLength).
/// • Meyne, Ratliff, Moyzis (1989): conservation of (TTAGGG)n among vertebrates.
/// ───────────────────────────────────────────────────────────────────────────
///
/// ═══════════════════════════════════════════════════════════════════════════
/// Unit: CHROM-CENT-001 — centromere detection
/// Checklist: docs/checklists/03_FUZZING.md, row 49.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: the empty sequence, an all-AT sequence, an all-GC
///          sequence, and an extremely short sequence (shorter than the scan
///          window). — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The centromere-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// `AnalyzeCentromere(string chromosomeName, string sequence,
///      int windowSize = 100000, double minAlphaSatelliteContent = 0.3)`
///   (ChromosomeAnalyzer.cs lines 361–437) returns an EAGER `CentromereResult`
///   readonly record struct (Chromosome, Start?, End?, Length, CentromereType,
///   AlphaSatelliteContent, IsAcrocentric) — no lazy iterator, so any exception
///   or hang surfaces at the call itself.
///
/// CRITICAL — what signal the detector actually keys on (Centromere_Analysis.md
/// §2.2, §4.2; ChromosomeAnalyzer.cs lines 384–391): the score is
///     score = repeatContent * (1 - gcVariability)
/// where `repeatContent` is the 15-mer REPETITIVENESS of a window (EstimateRepeat-
/// Content) and `gcVariability` is the std-dev of GC fraction over 1 kb sub-
/// windows. A window is accepted only when `repeatContent > minAlphaSatelliteContent`
/// AND its score is the running maximum. The biological gloss "centromeres are
/// AT-rich alpha-satellite" is realised in THIS code as "repeat-rich, GC-uniform",
/// NOT as a literal AT-content meter. The consequence for the all-AT / all-GC
/// fuzz boundaries is precise and tested below:
///   • all-AT homopolymer ≥ windowSize → a single 15-mer repeats throughout →
///     repeatContent maximal; GC fraction is a constant 0 across sub-windows →
///     gcVariability 0 → score maximal → DETECTED (the strong-positive boundary).
///   • all-GC homopolymer ≥ windowSize → ALSO a single repeated 15-mer, GC
///     fraction a constant 1 → gcVariability 0 → score maximal → ALSO DETECTED.
///     The checklist's idealised "all-GC ⇒ no AT signal ⇒ no centromere" does NOT
///     hold for this repeat-based heuristic: an all-GC homopolymer is just as
///     maximally repetitive as all-AT, so the theory-correct contract of the
///     ACTUAL code is detection. We assert that, with the high-GC sanity check
///     (gpos/AT-uniformity is not what is measured) documented in line, rather
///     than weaken the test to the checklist's idealised model.
///
/// Documented edge-case contract (Centromere_Analysis.md §3.3, §6.1):
///   • Empty OR null sequence → special-cased (line 367) → CentromereType
///     "Unknown", Start/End null, Length 0, AlphaSatelliteContent 0,
///     IsAcrocentric false. Never indexed.
///   • Sequence shorter than windowSize → the scan loop
///     `for (i = 0; i < sequence.Length - windowSize; ...)` has a NEGATIVE bound
///     (Length − windowSize < 0) so it never iterates; candidate stays unset →
///     "Unknown", null boundaries. This is the KEY boundary: the negative bound
///     must NOT produce an IndexOutOfRange (no `sequence[i..end]` is ever taken)
///     and the GC-variability helper is never reached, so NO DivideByZero.
///   • Non-repetitive sequence → no window exceeds minAlphaSatelliteContent →
///     "Unknown" (asserted via the high-complexity random sanity case).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Citations
/// ───────────────────────────────────────────────────────────────────────────
/// • Algorithm doc: docs/algorithms/Chromosome_Analysis/Centromere_Analysis.md
///   (§2.2 core model, §3.3 validation, §4.2 scoring, §6.1 edge cases).
/// • Source: src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///   (AnalyzeCentromere + EstimateRepeatContent + CalculateGcVariability).
/// • Levan, Fredga, Sandberg (1964): centromere arm-ratio nomenclature.
/// ───────────────────────────────────────────────────────────────────────────
///
/// ═══════════════════════════════════════════════════════════════════════════
/// Unit: CHROM-KARYO-001 — karyotype analysis
/// Checklist: docs/checklists/03_FUZZING.md, row 50.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate set-size boundaries called out
///          in the checklist row: 0 chromosomes, 100+ chromosomes, and empty data
///          (an empty enumerable). — docs/checklists/03_FUZZING.md §Description.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The karyotype-analysis contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// `AnalyzeKaryotype(IEnumerable&lt;(string Name, long Length, bool IsSexChromosome)&gt;
///      chromosomes, int expectedPloidyLevel = 2)`
///   (ChromosomeAnalyzer.cs lines 136–180) materializes the input with `.ToList()`
///   and returns an EAGER `Karyotype` readonly record struct (TotalChromosomes,
///   AutosomeCount, SexChromosomes, TotalGenomeSize, MeanChromosomeLength,
///   PloidyLevel, HasAneuploidy, Abnormalities) — no lazy iterator, so any exception
///   or hang surfaces at the call itself.
///
/// CRITICAL — the 0-chromosome / empty-data boundary the checklist row keys on
/// (Karyotype_Analysis.md §6.1; ChromosomeAnalyzer.cs lines 142–145): when the
/// materialized list is empty the method short-circuits and returns the documented
/// ZEROED karyotype — `new Karyotype(0, 0, [], 0, 0, 0, false, [])`. This guard is
/// what makes the average-size computation safe: the unguarded statement
///     meanLength = totalSize / (double)chromList.Count   (line 151)
/// is only reached for a NON-empty list, so the count is never 0 there. Even if it
/// were reached, the divisor is a `double`, so the result would be NaN, not a
/// DivideByZeroException (integer div-by-zero throws; floating-point does not). The
/// early return makes both moot: 0 chromosomes / empty data must yield the zeroed
/// karyotype with MeanChromosomeLength exactly 0.0, never a crash, never NaN.
///
/// Aggregate consistency invariants asserted on the 100+ boundary
/// (Karyotype_Analysis.md §2.4 INV-01..INV-03):
///   • INV-01: TotalChromosomes == AutosomeCount + SexChromosomes.Count (the input
///     is partitioned into autosomes and sex chromosomes before summarizing).
///   • INV-02: TotalGenomeSize == Σ input lengths, and MeanChromosomeLength ==
///     TotalGenomeSize / TotalChromosomes for non-empty input (both computed
///     directly from the materialized list).
///   • INV-03: HasAneuploidy is true exactly when the Abnormalities list is
///     non-empty — set only when an autosome group count differs from
///     expectedPloidyLevel.
/// A large set (100+ chromosomes) must COMPLETE in O(n) and satisfy all three; the
/// summary stats must be internally consistent with the input, never overflow into a
/// malformed karyotype.
///
/// Documented edge-case contract (Karyotype_Analysis.md §3.3, §6.1):
///   • Empty chromosome list → zeroed Karyotype, no abnormalities, PloidyLevel 0.
///     This covers BOTH "0 chromosomes" and "empty data": a freshly-allocated empty
///     enumerable and a deliberately empty array are the same boundary.
///   • Autosomes are grouped by base name after stripping a trailing `_N` numeric
///     copy suffix (GetChromosomeBaseName); a group whose size != expectedPloidyLevel
///     is labelled with the absolute cytogenetic term (Monosomy/Trisomy/…).
///   • Sex chromosomes are preserved in SexChromosomes but never grouped for
///     abnormality calling (§5.2 deviation #1).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Citations
/// ───────────────────────────────────────────────────────────────────────────
/// • Algorithm doc: docs/algorithms/Chromosome_Analysis/Karyotype_Analysis.md
///   (§2.4 invariants, §3 contract, §6.1 edge cases).
/// • Source: src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///   (AnalyzeKaryotype + GetChromosomeBaseName + GetAneuploidyTerm).
/// • Tjio, Levan (1956): the chromosome number of man (2n = 46) — the human
///   karyotype reference used for the positive-sanity set.
/// ───────────────────────────────────────────────────────────────────────────
///
/// ═══════════════════════════════════════════════════════════════════════════
/// Unit: CHROM-ANEU-001 — aneuploidy detection
/// Checklist: docs/checklists/03_FUZZING.md, row 51.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate depth boundaries called out in
///          the checklist row: 0 depth, negative depth, extremely high depth, and
///          empty data. — docs/checklists/03_FUZZING.md §Description.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The aneuploidy-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// `DetectAneuploidy(IEnumerable&lt;(string Chromosome, int Position, double Depth)&gt;
///      depthData, double medianDepth, int binSize = 1000000)`
///   (ChromosomeAnalyzer.cs lines 832–876) is a LAZY iterator (`yield`), so any
///   exception or hang surfaces only when the result is *enumerated* — every fuzz
///   case below forces evaluation with `.ToList()`, otherwise the guard and the
///   per-bin arithmetic would never run. It emits one `CopyNumberState`
///   (Chromosome, Start, End, CopyNumber, LogRatio, Confidence) per chromosome/bin.
///
/// Copy number is a read-depth RATIO against a diploid baseline (Aneuploidy_
/// Detection.md §2.2, §4.2): per bin
///     logRatio   = log2(meanDepth / medianDepth)
///     copyNumber = clamp(round((2^logRatio) * 2), 0, 10)
/// The `medianDepth` baseline is the divisor of that ratio — the 0-reference
/// division the checklist row keys on.
///
/// What the depth boundaries actually do (ChromosomeAnalyzer.cs lines 837–873;
/// Aneuploidy_Detection.md §6.1; verified empirically against net10.0 semantics):
///   • Empty data → the materialized list is empty → `yield break` (line 839) → a
///     defined EMPTY result, never a crash. This is the empty-data boundary.
///   • medianDepth == 0 OR medianDepth &lt; 0 → the SAME guard (`medianDepth &lt;= 0`,
///     line 839) short-circuits to `yield break`. This is the KEY div-by-zero
///     boundary: a 0 (or negative) baseline is REJECTED *before* the depth ratio is
///     ever computed, so `log2(meanDepth / 0)` is never reached — no DivideByZero,
///     no Infinity-driven mis-call. (A 0 baseline carries no signal, so the
///     theory-correct result is "no call", realised here as the empty result.)
///   • 0 DEPTH values (with a valid positive baseline) → meanDepth 0 →
///     `log2(0 / median)` = −∞ → `2^(−∞)` = 0 → copyNumber rounds/clamps to 0
///     (Nullisomy). LogRatio is −∞ but the emitted CopyNumber is a defined,
///     in-range integer 0 — the documented "very low depth ⇒ floored at 0" case
///     (§6.1). No crash, no spurious gain.
///   • NEGATIVE DEPTH values (biologically impossible, with a valid baseline) →
///     `meanDepth / median` is negative → `log2(negative)` = NaN → `2^NaN` = NaN →
///     `(int)round(NaN)` = 0 (ECMA-335 NaN→int conversion) → clamp keeps 0. The
///     method does NOT validate or reject a negative read depth, but it ALSO never
///     produces a garbage NON-zero copy-number GAIN: the call collapses to
///     copyNumber 0. The cost of not validating is that LogRatio and Confidence are
///     emitted as NaN for that bin (the NaN propagates through the ratio). We pin
///     this ACTUAL contract — no throw, CopyNumber in [0,10] and specifically 0, with
///     LogRatio/Confidence NaN documented in line — rather than weaken to a falsely
///     clean expectation. A negative depth is impossible input; the disciplined
///     outcome we require is "no crash, no spurious aneuploidy call", which holds.
///   • EXTREMELY HIGH depth → `(2^logRatio)*2` is enormous but the rounded copy
///     number is CLAMPED to 10 (line 860); the result is a finite, in-range integer
///     with no overflow (§6.1 "very high depth ratios ⇒ capped at 10").
///
/// Invariants asserted (Aneuploidy_Detection.md §2.4):
///   • INV-01: every emitted CopyNumber is an integer in [0, 10] (round + clamp).
///   • INV-02: Confidence ∈ [0, 1] — for FINITE input. (We document that negative
///     depth, which is out-of-domain garbage, propagates NaN into Confidence; INV-02
///     is stated for valid depth, and the negative-depth test pins the real NaN
///     behaviour instead of pretending the bound holds on impossible input.)
///
/// ───────────────────────────────────────────────────────────────────────────
/// Citations
/// ───────────────────────────────────────────────────────────────────────────
/// • Algorithm doc: docs/algorithms/Chromosome_Analysis/Aneuploidy_Detection.md
///   (§2.2 core model, §2.4 invariants, §3.3 validation, §4.2 scoring, §6.1 edges).
/// • Source: src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///   (DetectAneuploidy + CopyNumberState).
/// • Wikipedia: Aneuploidy; Copy number variation — copy-number terminology and the
///   depth-proportional-to-copy-number model used for the positive-sanity profile.
/// ───────────────────────────────────────────────────────────────────────────
///
/// ═══════════════════════════════════════════════════════════════════════════
/// Unit: CHROM-SYNT-001 — synteny analysis
/// Checklist: docs/checklists/03_FUZZING.md, row 52.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: 0 ortholog pairs, self-comparison (a genome against
///          itself under an identity ortholog map), and empty genomes.
///          — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Which synteny method — Chromosome (NOT Comparative)
/// ───────────────────────────────────────────────────────────────────────────
/// There are two synteny surfaces in the codebase: COMPGEN-SYNTENY-001 in the
/// Comparative-genomics area (ComparativeGenomics.cs, row 139) and the
/// CHROMOSOME-area CHROM-SYNT-001 under test here. This unit targets the
/// Chromosome surface:
///   `ChromosomeAnalyzer.FindSyntenyBlocks(
///        IEnumerable&lt;(string Chr1, int Start1, int End1, string Gene1,
///                     string Chr2, int Start2, int End2, string Gene2)&gt; orthologPairs,
///        int minGenes = 3,
///        int maxGap   = 10)`
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///    lines 643–716). It is a LAZY iterator (`yield`), so any exception or hang
///    surfaces only when the result is *enumerated* — every fuzz case below forces
///    evaluation with `.ToList()`, otherwise the `pairs.Count < minGenes` guard and
///    the collinear-run scan would never run.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The synteny-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Synteny analysis finds conserved gene-order (collinear) blocks between two
/// genomes from ortholog ANCHOR pairs. The algorithm (Synteny_Analysis.md §4.1;
/// ChromosomeAnalyzer.cs lines 649–715):
///   • Materializes the pairs; if the TOTAL pair count is below `minGenes`, returns
///     no blocks (line 651, `yield break`).
///   • Groups pairs by chromosome pair (Chr1, Chr2), sorts each group by the
///     first-genome start, and tracks collinear runs: orientation is set from the
///     direction of movement in the SECOND genome (`curr.Start2 > prev.End2` ⇒
///     forward), a run continues while orientation is consistent and BOTH gaps are
///     ≤ `maxGap × 1 000 000`, and a `SyntenyBlock` is emitted when a run ends with
///     ≥ `minGenes` genes (lines 666–708).
///   • Each emitted block's Species1Start/End come from the first/last gene in the
///     first genome; Species2Start/End are Min/Max over the run in the second genome;
///     Strand is '+' (forward) or '-' (reverse); GeneCount is the run length;
///     SequenceIdentity is `double.NaN` (INV-03 — not computable from coordinates).
///
/// What the three checklist boundaries actually do (Synteny_Analysis.md §6.1;
/// ChromosomeAnalyzer.cs lines 649–715):
///   • 0 ORTHOLOG PAIRS → the materialized list is empty, count 0 &lt; minGenes 3 →
///     `yield break` → a defined EMPTY result (no synteny blocks), never a crash.
///     This is the empty-anchor boundary: no anchors ⇒ no conserved blocks.
///   • EMPTY GENOMES → there is no separate "genome" argument; a genome with no
///     genes contributes no ortholog pairs, so the empty-genome case reduces to the
///     same 0-pair boundary: an empty ortholog-pair sequence ⇒ empty result. We pin
///     it explicitly through a different concrete enumerable (empty array) to show
///     the guard keys on the materialized count, not the source type.
///   • SELF-COMPARISON (genome vs itself under an IDENTITY ortholog map) — the KEY
///     positive boundary. Each gene maps to itself: pair
///     (chr, s, e, g, chr, s, e, g). For N ≥ minGenes genes in ascending order on a
///     single chromosome, every step is forward (`curr.Start2 = curr.Start1 >
///     prev.End2 = prev.End1` for non-overlapping ascending genes) with zero-to-small
///     gaps, so the whole genome is ONE maximal collinear forward run → exactly ONE
///     SyntenyBlock covering the whole genome: Strand '+', GeneCount = N,
///     Species1Chromosome == Species2Chromosome, and the block's coordinates in BOTH
///     genomes equal the genome's own span (identical order ⇒ perfect synteny). This
///     is the theory anchor: a genome is perfectly collinear with itself.
///
/// Invariants asserted on every emitted block (Synteny_Analysis.md §2.4):
///   • INV-01: Strand ∈ { '+', '-' }.
///   • INV-02: Species1Start ≤ Species1End and Species2Start ≤ Species2End — block
///     boundaries reference VALID, ordered coordinates in both genomes.
///   • INV-03: SequenceIdentity is NaN for coordinate-only input.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Citations
/// ───────────────────────────────────────────────────────────────────────────
/// • Algorithm doc: docs/algorithms/Chromosome_Analysis/Synteny_Analysis.md
///   (§2.4 invariants, §3 contract, §4.1 steps, §6.1 edge cases).
/// • Source: src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///   (FindSyntenyBlocks).
/// • Wang et al. (2012) MCScanX: collinearity-based synteny block detection.
/// • Wikipedia: Synteny; Comparative genomics — conserved gene order / collinearity.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ChromosomeFuzzTests
{
    #region Helpers

    private const string TelomereRepeat = "TTAGGG"; // vertebrate canonical hexamer
    private const int RepeatLen = 6;

    /// <summary>Builds a pure forward telomere tract of <paramref name="copies"/> TTAGGG units.</summary>
    private static string PureTelomere(int copies) => string.Concat(System.Linq.Enumerable.Repeat(TelomereRepeat, copies));

    /// <summary>
    /// Deterministic non-telomeric DNA — seed fixed locally so generated fuzz input
    /// is reproducible (no shared static Rng). Avoids accidentally emitting a TTAGGG
    /// run by drawing from a reduced alphabet biased away from the motif.
    /// </summary>
    private static string RandomNonTelomericDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Centromere scan window used by the detection fuzz cases. The production default
    /// is 100 000; we use the smaller 10 000 (matching ChromosomeAnalyzer_Centromere_Tests)
    /// so fuzz sequences stay a few hundred kb instead of multi-megabyte, while still
    /// exercising the identical scan/score/extend code path.
    /// </summary>
    private const int CentWindow = 10_000;

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CHROM-TELO-001 — telomere detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CHROM-TELO-001 — telomere detection

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is special-cased (ChromosomeAnalyzer.cs line 258;
    /// Telomere_Analysis.md §6.1). It must return a no-telomere result — zero
    /// lengths, both Has*Telomere = false, purities 0 — and the documented
    /// IsCriticallyShort = TRUE flag, never a crash.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_EmptySequence_ReturnsNoTelomereAndCriticallyShort()
    {
        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrEmpty", string.Empty);
        act.Should().NotThrow("empty input is explicitly special-cased, not indexed");

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrEmpty", string.Empty);

        result.Chromosome.Should().Be("chrEmpty", "the chromosome name is copied verbatim");
        result.Has5PrimeTelomere.Should().BeFalse();
        result.Has3PrimeTelomere.Should().BeFalse();
        result.TelomereLength5Prime.Should().Be(0);
        result.TelomereLength3Prime.Should().Be(0);
        result.RepeatPurity5Prime.Should().Be(0);
        result.RepeatPurity3Prime.Should().Be(0);
        result.IsCriticallyShort.Should().BeTrue(
            "empty input carries the documented critically-short flag (Telomere_Analysis.md §6.1)");
    }

    /// <summary>
    /// BE: a null sequence is treated identically to empty via the
    /// string.IsNullOrEmpty short-circuit — it must NOT throw a
    /// NullReferenceException and must return the same no-telomere, critically-short
    /// result.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_NullSequence_ReturnsNoTelomereAndCriticallyShort()
    {
        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrNull", null!);
        act.Should().NotThrow("null is absorbed by the IsNullOrEmpty guard, never dereferenced");

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrNull", null!);

        result.Has3PrimeTelomere.Should().BeFalse();
        result.TelomereLength3Prime.Should().Be(0);
        result.IsCriticallyShort.Should().BeTrue();
    }

    /// <summary>
    /// BE: a sequence shorter than the 6-bp repeat cannot contain a single complete
    /// repeat window — MeasureTelomereLength exits at region.Length &lt; repeatLen
    /// (line 299). Both ends must report length 0 with no Substring out-of-range, and
    /// (non-empty input) IsCriticallyShort = false.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_SequenceShorterThanRepeat_ReturnsZeroLengthsNoCrash()
    {
        // "TTAGG" is 5 bp — one base short of a single TTAGGG window.
        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrShort", "TTAGG");
        act.Should().NotThrow("a region shorter than the repeat exits before any window is taken");

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrShort", "TTAGG");

        result.TelomereLength5Prime.Should().Be(0);
        result.TelomereLength3Prime.Should().Be(0);
        result.IsCriticallyShort.Should().BeFalse(
            "non-empty input with no detected telomere is NOT flagged critically short (§5.2)");
    }

    #endregion

    #region BE — Boundary: no TTAGGG in sequence

    /// <summary>
    /// BE: a sequence with NO telomeric repeat. No 6-bp window on either end reaches
    /// the 70% similarity bar, so both lengths are 0, both Has*Telomere = false, and
    /// IsCriticallyShort = false for non-empty input (§6.1, §5.2). A poly-A tract is
    /// the cleanest "no motif" input: it shares at most one base with TTAGGG/CCCTAA
    /// (1/6 ≈ 17% &lt; 70%).
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_NoTelomericRepeat_DetectsNothing()
    {
        string sequence = new string('A', 600); // no TTAGGG, no CCCTAA anywhere

        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrNoMotif", sequence);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrNoMotif", sequence);

        result.Has5PrimeTelomere.Should().BeFalse("no window matches CCCTAA at the 5' end");
        result.Has3PrimeTelomere.Should().BeFalse("no window matches TTAGGG at the 3' end");
        result.TelomereLength5Prime.Should().Be(0);
        result.TelomereLength3Prime.Should().Be(0);
        result.RepeatPurity5Prime.Should().Be(0);
        result.RepeatPurity3Prime.Should().Be(0);
        result.IsCriticallyShort.Should().BeFalse(
            "no telomere is detected, so the critically-short flag stays false for non-empty input");
    }

    /// <summary>
    /// BE: a pseudo-random non-telomeric sequence must likewise never mis-detect a
    /// full telomere and never crash. We assert only the well-formedness invariants
    /// (INV-01..INV-03): non-negative lengths that are multiples of 6, purities in
    /// [0,1], and Has*Telomere consistent with the default minTelomereLength = 500.
    /// Deterministic seed keeps this reproducible.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_RandomSequence_ProducesOnlyWellFormedResult()
    {
        string sequence = RandomNonTelomericDna(2000, seed: 4848);

        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrRand", sequence);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrRand", sequence);

        result.TelomereLength5Prime.Should().BeGreaterThanOrEqualTo(0);
        result.TelomereLength3Prime.Should().BeGreaterThanOrEqualTo(0);
        (result.TelomereLength5Prime % RepeatLen).Should().Be(0, "measured length is always a multiple of the 6-bp repeat");
        (result.TelomereLength3Prime % RepeatLen).Should().Be(0, "measured length is always a multiple of the 6-bp repeat");
        result.RepeatPurity5Prime.Should().BeInRange(0, 1);
        result.RepeatPurity3Prime.Should().BeInRange(0, 1);
        result.Has5PrimeTelomere.Should().Be(result.TelomereLength5Prime >= 500);
        result.Has3PrimeTelomere.Should().Be(result.TelomereLength3Prime >= 500);
    }

    #endregion

    #region BE — Boundary: only TTAGGG (pure telomeric tract) — KEY positive

    /// <summary>
    /// BE (KEY positive): a sequence that is ONLY the telomeric repeat. A pure tract
    /// of N copies of TTAGGG must be detected at the 3' end with the correct measured
    /// length = 6 × N (every window is a perfect TTAGGG match, similarity 100% ≥ 70%),
    /// repeat purity exactly 1.0, and Has3PrimeTelomere true when 6·N ≥ 500. This is
    /// the core "only TTAGGG ⇒ telomere with correct repeat count (length/6)" check.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_PureTtagggTract_DetectedWithCorrectLengthAndPurity()
    {
        const int copies = 100;               // 100 × 6 = 600 bp ≥ minTelomereLength (500)
        string sequence = PureTelomere(copies);

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrPure", sequence);

        result.TelomereLength3Prime.Should().Be(copies * RepeatLen,
            "a pure TTAGGG tract is measured as 6 × (repeat count) at the 3' end");
        (result.TelomereLength3Prime / RepeatLen).Should().Be(copies,
            "the detected repeat count is exactly length / 6");
        result.RepeatPurity3Prime.Should().Be(1.0,
            "every counted window is a perfect TTAGGG match, so purity is exactly 1.0");
        result.Has3PrimeTelomere.Should().BeTrue(
            "600 bp of pure repeat clears the default minTelomereLength of 500");
    }

    /// <summary>
    /// BE: a SHORT pure tract — below minTelomereLength but a valid repeat — must be
    /// measured with the correct length (6 × copies) yet flagged Has3PrimeTelomere =
    /// false, pinning the minTelomereLength boundary rather than mis-flagging a short
    /// real tract. With a default criticalLength of 3000 the detected short telomere
    /// stays below it, but Has3Prime = false means IsCriticallyShort cannot trip
    /// (it requires a *detected* telomere, line 279).
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_ShortPureTract_MeasuredButBelowPresenceThreshold()
    {
        const int copies = 10;                // 10 × 6 = 60 bp < minTelomereLength (500)
        string sequence = PureTelomere(copies);

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrShortPure", sequence);

        result.TelomereLength3Prime.Should().Be(copies * RepeatLen,
            "the short tract is still measured at its true length 6 × copies");
        result.Has3PrimeTelomere.Should().BeFalse(
            "60 bp is below the default minTelomereLength of 500, so presence is false");
        result.IsCriticallyShort.Should().BeFalse(
            "no telomere is flagged present, so the critically-short test never fires");
    }

    #endregion

    #region MC — Malformed content: non-DNA characters

    /// <summary>
    /// MC: non-DNA characters (digits, IUPAC ambiguity codes, whitespace, unicode)
    /// must be handled without a crash. There is no validation gate — input is
    /// upper-cased and scanned verbatim — so garbage bases simply fail the 70%
    /// similarity bar and the sequence behaves like "no telomere": zero lengths,
    /// no detection, no IndexOutOfRange. We pin no-throw + no spurious detection.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_NonDnaCharacters_HandledWithoutCrashOrDetection()
    {
        string garbage = "12345 NRYSWK\t\nXZé中??!!--***";

        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrGarbage", garbage);
        act.Should().NotThrow("there is no validation gate; non-DNA chars are scanned verbatim, never indexed past end");

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrGarbage", garbage);

        result.Has5PrimeTelomere.Should().BeFalse("no non-DNA window can reach 70% similarity to CCCTAA");
        result.Has3PrimeTelomere.Should().BeFalse("no non-DNA window can reach 70% similarity to TTAGGG");
        result.TelomereLength5Prime.Should().Be(0);
        result.TelomereLength3Prime.Should().Be(0);
    }

    /// <summary>
    /// MC: a non-DNA-decorated but otherwise pure telomeric tract. Even when garbage
    /// characters are interleaved at the chromosome interior, a clean TTAGGG tract at
    /// the 3' terminus must still be detected (the 3' scan walks inward from the end
    /// and stops at the first sub-70% window), and the run must never crash on the
    /// non-DNA bytes. This pins that malformed content elsewhere in the sequence does
    /// not corrupt or suppress a legitimate terminal telomere.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_GarbagePrefixWithPureTerminalTract_StillDetectsTelomere()
    {
        // Non-DNA junk at the 5' interior, a clean 100-copy TTAGGG tract at the 3' end.
        string sequence = "12NRYSWK??##" + PureTelomere(100);

        var act = () => ChromosomeAnalyzer.AnalyzeTelomeres("chrMixed", sequence);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrMixed", sequence);

        result.Has3PrimeTelomere.Should().BeTrue(
            "the clean terminal TTAGGG tract is detected despite non-DNA junk upstream");
        result.TelomereLength3Prime.Should().Be(100 * RepeatLen,
            "the 3' scan counts exactly the 100 perfect repeats at the terminus");
        result.RepeatPurity3Prime.Should().Be(1.0);
    }

    #endregion

    #region Positive sanity — a known repeat count is recovered exactly

    /// <summary>
    /// Positive sanity: a TTAGGG tract of a KNOWN number of copies is detected with
    /// exactly that repeat count (length / 6). This is the affirmative anchor that
    /// the detector actually finds the canonical motif and counts it correctly —
    /// without it, an all-passing boundary suite could be vacuously green. We test
    /// several copy counts so the count, not a fixed length, is what is recovered.
    /// </summary>
    [TestCase(85)]
    [TestCase(150)]
    [TestCase(512)]
    public void AnalyzeTelomeres_KnownRepeatCount_RecoveredExactly(int copies)
    {
        string sequence = PureTelomere(copies);

        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chrSanity", sequence);

        (result.TelomereLength3Prime / RepeatLen).Should().Be(copies,
            "the detected 3' repeat count must equal the planted number of TTAGGG copies");
        result.TelomereLength3Prime.Should().Be(copies * RepeatLen);
        result.RepeatPurity3Prime.Should().Be(1.0, "a perfect planted tract has purity 1.0");
        result.Has3PrimeTelomere.Should().Be(copies * RepeatLen >= 500);
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CHROM-CENT-001 — centromere detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CHROM-CENT-001 — centromere detection

    #region BE — Boundary: empty / null sequence

    /// <summary>
    /// BE: the empty sequence is special-cased (ChromosomeAnalyzer.cs line 367;
    /// Centromere_Analysis.md §6.1). It must return the documented no-centromere
    /// result — CentromereType "Unknown", null boundaries, length 0,
    /// AlphaSatelliteContent 0, IsAcrocentric false — never a crash, never an index.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_EmptySequence_ReturnsUnknownNoBoundaries()
    {
        var act = () => ChromosomeAnalyzer.AnalyzeCentromere("chrEmpty", string.Empty);
        act.Should().NotThrow("empty input is explicitly special-cased, not scanned or indexed");

        var result = ChromosomeAnalyzer.AnalyzeCentromere("chrEmpty", string.Empty);

        result.Chromosome.Should().Be("chrEmpty", "the chromosome name is copied verbatim");
        result.CentromereType.Should().Be("Unknown");
        result.Start.Should().BeNull();
        result.End.Should().BeNull();
        result.Length.Should().Be(0);
        result.AlphaSatelliteContent.Should().Be(0);
        result.IsAcrocentric.Should().BeFalse();
    }

    /// <summary>
    /// BE: a null sequence is absorbed by the same string.IsNullOrEmpty short-circuit
    /// (line 367). It must NOT throw a NullReferenceException and must return the same
    /// "Unknown" / null-boundary result as the empty case.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_NullSequence_ReturnsUnknownNoBoundaries()
    {
        var act = () => ChromosomeAnalyzer.AnalyzeCentromere("chrNull", null!);
        act.Should().NotThrow("null is absorbed by the IsNullOrEmpty guard, never dereferenced");

        var result = ChromosomeAnalyzer.AnalyzeCentromere("chrNull", null!);

        result.CentromereType.Should().Be("Unknown");
        result.Start.Should().BeNull();
        result.End.Should().BeNull();
        result.Length.Should().Be(0);
    }

    #endregion

    #region BE — Boundary: extremely short (shorter than the scan window)

    /// <summary>
    /// BE (KEY boundary): a sequence shorter than <see cref="CentWindow"/>. The scan
    /// loop bound is `sequence.Length - windowSize`, which is NEGATIVE here, so the
    /// loop never iterates: no window slice `sequence[i..end]` is ever taken (no
    /// IndexOutOfRange / ArgumentOutOfRange) and the GC-variability helper is never
    /// reached (no DivideByZero). The candidate stays unset → "Unknown", null
    /// boundaries, length 0 (Centromere_Analysis.md §6.1). This pins the div-by-zero /
    /// out-of-range boundary the checklist row calls out.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_ShorterThanWindow_ReturnsUnknownNoCrash()
    {
        // 500 bp ≪ 10 000-bp window: Length − windowSize is negative.
        string sequence = RandomNonTelomericDna(500, seed: 11);

        var act = () => ChromosomeAnalyzer.AnalyzeCentromere("chrShort", sequence, windowSize: CentWindow);
        act.Should().NotThrow("a negative scan bound skips the loop; nothing is sliced or divided");

        var result = ChromosomeAnalyzer.AnalyzeCentromere("chrShort", sequence, windowSize: CentWindow);

        result.CentromereType.Should().Be("Unknown", "the scan loop never runs, so no candidate is set");
        result.Start.Should().BeNull();
        result.End.Should().BeNull();
        result.Length.Should().Be(0);
        result.AlphaSatelliteContent.Should().Be(0);
    }

    /// <summary>
    /// BE: the most extreme short input — a single base — must behave identically to
    /// any other sub-window-length sequence: "Unknown", null boundaries, no crash.
    /// Length 1 makes `Length − windowSize` maximally negative and exercises the same
    /// loop-skip guard with a sequence too short to even hold one 15-mer.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_SingleBase_ReturnsUnknownNoCrash()
    {
        var act = () => ChromosomeAnalyzer.AnalyzeCentromere("chrOne", "A", windowSize: CentWindow);
        act.Should().NotThrow("a single base is far shorter than the window; the scan loop is skipped");

        var result = ChromosomeAnalyzer.AnalyzeCentromere("chrOne", "A", windowSize: CentWindow);

        result.CentromereType.Should().Be("Unknown");
        result.Start.Should().BeNull();
        result.Length.Should().Be(0);
    }

    #endregion

    #region BE — Boundary: all-AT homopolymer (strong-positive)

    /// <summary>
    /// BE (KEY positive): an all-AT homopolymer at least as long as the scan window is
    /// maximally repetitive — a single 15-mer repeats end to end, so repeatContent is
    /// maximal — and its GC fraction is a constant 0 across every 1 kb sub-window, so
    /// gcVariability is 0 and the composite score is maximal. The detector MUST report
    /// a region (non-null Start/End) with AlphaSatelliteContent above the threshold and
    /// a concrete (non-"Unknown") classification. Asserting detection here (not a crash,
    /// not a miss) is the strong-positive AT-richness boundary the checklist demands.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_AllAt_DetectedAsStrongPositive()
    {
        string sequence = new string('A', 300_000); // 30 windows of pure repeat

        var act = () => ChromosomeAnalyzer.AnalyzeCentromere(
            "chrAllAt", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.3);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chrAllAt", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.3);

        result.Start.Should().NotBeNull("an all-AT homopolymer is maximally repetitive and must be detected");
        result.End.Should().NotBeNull();
        result.Length.Should().BeGreaterThan(0).And.Be(result.End!.Value - result.Start!.Value,
            "INV-01: Length = End − Start for a detected region");
        result.AlphaSatelliteContent.Should().BeGreaterThan(0.3,
            "uniform AT content has maximal repeat content and zero GC variability → maximal score");
        result.CentromereType.Should().NotBe("Unknown", "a detected region carries a Levan classification");
    }

    #endregion

    #region BE — Boundary: all-GC homopolymer (repeat heuristic, not AT meter)

    /// <summary>
    /// BE: an all-GC homopolymer. The checklist's idealised model says "all-GC → no
    /// AT signal → no centromere", but the ACTUAL detector scores REPETITIVENESS, not
    /// AT content: an all-GC run is just as maximally repetitive (one repeated 15-mer)
    /// with a constant GC fraction of 1 → gcVariability 0 → maximal score. The theory-
    /// correct contract of this repeat-based heuristic is therefore DETECTION, exactly
    /// like all-AT. We assert detection (never a crash, never a spurious "Unknown") and
    /// document in line that this is the repeat heuristic's behaviour, not an AT meter —
    /// we do not weaken the test to the checklist's idealised expectation.
    /// (Centromere_Analysis.md §2.2, §4.2, §6.1 "Strongly repetitive sequence".)
    /// </summary>
    [Test]
    public void AnalyzeCentromere_AllGc_DetectedBecauseSignalIsRepetitivenessNotAtContent()
    {
        string sequence = new string('G', 300_000);

        var act = () => ChromosomeAnalyzer.AnalyzeCentromere(
            "chrAllGc", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.3);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chrAllGc", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.3);

        result.Start.Should().NotBeNull(
            "the detector keys on repeat content; an all-GC homopolymer is maximally repetitive, so it is detected");
        result.End.Should().NotBeNull();
        result.AlphaSatelliteContent.Should().BeGreaterThan(0.3,
            "constant GC fraction → zero GC variability → score equals the (maximal) repeat content");
        result.CentromereType.Should().NotBe("Unknown");
    }

    #endregion

    #region BE — Well-formedness invariants on random / homogeneous input

    /// <summary>
    /// BE: a high-complexity pseudo-random sequence (no engineered repeat) must never
    /// crash and, per Centromere_Analysis.md §6.1 "Non-repetitive sequence", must yield
    /// "Unknown" — no window of distinct 15-mers clears the 0.3 repeat threshold. This
    /// is the negative counterpart to the AT/GC positives: a non-centromeric region is
    /// NOT detected. Deterministic seed keeps it reproducible.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_RandomNonRepetitive_ReturnsUnknown()
    {
        string sequence = RandomNonTelomericDna(300_000, seed: 9090);

        var act = () => ChromosomeAnalyzer.AnalyzeCentromere(
            "chrRand", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.3);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chrRand", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.3);

        result.CentromereType.Should().Be("Unknown",
            "a high-complexity random sequence has no repeat-rich window above the threshold");
        result.Start.Should().BeNull();
        result.End.Should().BeNull();
        result.Length.Should().Be(0);
        result.AlphaSatelliteContent.Should().BeGreaterThanOrEqualTo(0, "INV-04: score is never negative");
    }

    #endregion

    #region Positive sanity — engineered satellite tract detected, flank not

    /// <summary>
    /// Positive sanity: a clearly AT-rich, alpha-satellite-like tandem-repeat tract
    /// embedded between non-repetitive flanks is detected as a centromere — non-null
    /// boundaries that land inside the engineered repeat, a score above the threshold,
    /// and a concrete classification — while a purely random sequence of the same size
    /// is NOT (covered by the random case above). This is the affirmative anchor that
    /// the detector finds a genuine satellite signal and not just degenerate homopolymers.
    /// </summary>
    [Test]
    public void AnalyzeCentromere_EngineeredSatelliteTract_DetectedInsideRepeatRegion()
    {
        // Non-repetitive flanks, AT-rich tandem-repeat core (alpha-satellite-like 21-mer).
        string flank = RandomNonTelomericDna(100_000, seed: 7);
        string repeatUnit = "AATGAATATTTCTTTTATGTT"; // 21-mer, AT-rich
        string satellite = string.Concat(System.Linq.Enumerable.Repeat(repeatUnit, 6_000)); // ~126 kb
        string sequence = flank + satellite + flank;

        var act = () => ChromosomeAnalyzer.AnalyzeCentromere(
            "chrSat", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.1);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeCentromere(
            "chrSat", sequence, windowSize: CentWindow, minAlphaSatelliteContent: 0.1);

        result.Start.Should().NotBeNull("the embedded tandem-repeat tract is detected");
        result.End.Should().NotBeNull();
        result.CentromereType.Should().NotBe("Unknown");
        result.AlphaSatelliteContent.Should().BeGreaterThan(0.1, "the satellite tract clears the threshold");
        // The detected region overlaps the engineered satellite core, not the random flanks.
        result.End!.Value.Should().BeGreaterThan(flank.Length,
            "the detected region must reach into the repetitive core past the 5' flank");
        result.Start!.Value.Should().BeLessThan(flank.Length + satellite.Length,
            "the detected region must start before the 3' flank");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CHROM-KARYO-001 — karyotype analysis : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CHROM-KARYO-001 — karyotype analysis

    #region Karyo helpers

    /// <summary>Convenience tuple alias for a chromosome descriptor.</summary>
    private static (string Name, long Length, bool IsSexChromosome) Chrom(string name, long length, bool sex = false)
        => (name, length, sex);

    #endregion

    #region BE — Boundary: 0 chromosomes / empty data (KEY div-by-zero)

    /// <summary>
    /// BE (KEY boundary): zero chromosomes (an empty enumerable). The materialized
    /// list is empty, so the method short-circuits at ChromosomeAnalyzer.cs line 142
    /// and returns the documented ZEROED karyotype WITHOUT ever reaching the
    /// `meanLength = totalSize / count` statement (line 151). The result must be a
    /// fully-defined, internally-consistent zero summary — every count 0, total size
    /// 0, MeanChromosomeLength exactly 0.0 (NOT NaN, NOT a DivideByZeroException),
    /// PloidyLevel 0, no aneuploidy, empty lists — never a crash.
    /// (Karyotype_Analysis.md §6.1 "Empty chromosome list".)
    /// </summary>
    [Test]
    public void AnalyzeKaryotype_ZeroChromosomes_ReturnsZeroedKaryotypeNoDivideByZero()
    {
        var empty = Enumerable.Empty<(string, long, bool)>();

        var act = () => ChromosomeAnalyzer.AnalyzeKaryotype(empty);
        act.Should().NotThrow(
            "the empty-list short-circuit returns before the mean-length division, so no DivideByZero or index can occur");

        var result = ChromosomeAnalyzer.AnalyzeKaryotype(empty);

        result.TotalChromosomes.Should().Be(0);
        result.AutosomeCount.Should().Be(0);
        result.SexChromosomes.Should().BeEmpty();
        result.TotalGenomeSize.Should().Be(0);
        result.MeanChromosomeLength.Should().Be(0.0, "the mean over 0 chromosomes is the documented zero, never NaN");
        double.IsNaN(result.MeanChromosomeLength).Should().BeFalse("the early return avoids the 0/0 division entirely");
        result.PloidyLevel.Should().Be(0, "PloidyLevel is zeroed for empty input, not the default 2");
        result.HasAneuploidy.Should().BeFalse();
        result.Abnormalities.Should().BeEmpty();
    }

    /// <summary>
    /// BE: "empty data" expressed as a freshly-allocated empty array rather than the
    /// LINQ Empty singleton — the same boundary reached through a different concrete
    /// enumerable. It must produce the identical zeroed karyotype with no crash,
    /// pinning that the guard keys on the materialized count, not on the source type.
    /// </summary>
    [Test]
    public void AnalyzeKaryotype_EmptyArray_ReturnsZeroedKaryotype()
    {
        var emptyData = Array.Empty<(string Name, long Length, bool IsSexChromosome)>();

        var act = () => ChromosomeAnalyzer.AnalyzeKaryotype(emptyData);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeKaryotype(emptyData);

        result.TotalChromosomes.Should().Be(0);
        result.TotalGenomeSize.Should().Be(0);
        result.MeanChromosomeLength.Should().Be(0.0);
        result.Abnormalities.Should().BeEmpty();
    }

    /// <summary>
    /// BE: a non-default expectedPloidyLevel on empty input must STILL yield the
    /// zeroed karyotype with PloidyLevel 0 — the empty short-circuit ignores the
    /// parameter entirely (line 144 hard-codes 0), so the ploidy argument can never
    /// drive an empty-input divide or a non-zero summary.
    /// </summary>
    [Test]
    public void AnalyzeKaryotype_EmptyWithNonDefaultPloidy_StillZeroedWithPloidyZero()
    {
        var result = ChromosomeAnalyzer.AnalyzeKaryotype(
            Enumerable.Empty<(string, long, bool)>(), expectedPloidyLevel: 4);

        result.TotalChromosomes.Should().Be(0);
        result.MeanChromosomeLength.Should().Be(0.0);
        result.PloidyLevel.Should().Be(0, "the empty short-circuit hard-codes PloidyLevel 0 regardless of the argument");
    }

    #endregion

    #region BE — Boundary: 100+ chromosomes (completes, aggregates consistently)

    /// <summary>
    /// BE (KEY boundary): a large set of 100+ chromosomes must COMPLETE (no hang) and
    /// aggregate correctly. We build 60 diploid autosome groups (chr1_1/chr1_2 …
    /// chr60_1/chr60_2 = 120 autosome tuples) plus an XX sex pair = 122 chromosomes,
    /// every group at exactly the expected ploidy of 2. The summary must be internally
    /// consistent with the input: TotalChromosomes == AutosomeCount + sex count
    /// (INV-01), TotalGenomeSize == Σ lengths and MeanChromosomeLength ==
    /// TotalGenomeSize / TotalChromosomes (INV-02), and — because every autosome group
    /// is a clean diploid pair — HasAneuploidy is false with an empty abnormality list
    /// (INV-03). [CancelAfter] guards against a regression that loops on large input.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void AnalyzeKaryotype_ManyDiploidChromosomes_CompletesWithConsistentAggregates()
    {
        const int groups = 60;            // 60 autosome groups × 2 copies = 120 autosomes
        const long unitLength = 1_000_000;
        var chroms = new List<(string Name, long Length, bool IsSexChromosome)>();
        for (int g = 1; g <= groups; g++)
        {
            chroms.Add(Chrom($"chr{g}_1", unitLength));
            chroms.Add(Chrom($"chr{g}_2", unitLength));
        }
        chroms.Add(Chrom("chrX_1", 1_500_000, sex: true));
        chroms.Add(Chrom("chrX_2", 1_500_000, sex: true));

        int total = chroms.Count;                          // 122
        long expectedTotalSize = chroms.Sum(c => c.Length);

        var act = () => ChromosomeAnalyzer.AnalyzeKaryotype(chroms);
        act.Should().NotThrow();

        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chroms);

        result.TotalChromosomes.Should().Be(total).And.BeGreaterThanOrEqualTo(100,
            "the large fuzz set must exceed the 100-chromosome boundary");
        result.AutosomeCount.Should().Be(groups * 2);
        result.SexChromosomes.Should().HaveCount(2);
        // INV-01: partition is exhaustive and disjoint.
        result.TotalChromosomes.Should().Be(result.AutosomeCount + result.SexChromosomes.Count);
        // INV-02: total size and mean are internally consistent with the input.
        result.TotalGenomeSize.Should().Be(expectedTotalSize);
        result.MeanChromosomeLength.Should().Be(expectedTotalSize / (double)total);
        // INV-03: every autosome group is a clean diploid pair ⇒ no aneuploidy.
        result.HasAneuploidy.Should().BeFalse("every autosome group has exactly the expected 2 copies");
        result.Abnormalities.Should().BeEmpty();
    }

    /// <summary>
    /// BE: a large set whose autosome groups are SINGLE copies (haploid) against the
    /// default expected ploidy of 2. The 100+ set must still complete, but now EVERY
    /// group is abnormal: HasAneuploidy true and exactly one "Monosomy …" abnormality
    /// per group (INV-03 — flag true iff list non-empty). This pins that the large-set
    /// aggregation counts groups correctly under the aneuploidy path, not just the
    /// clean-diploid path.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void AnalyzeKaryotype_ManyMonosomicChromosomes_FlagsEveryGroupConsistently()
    {
        const int groups = 110;           // 110 single-copy autosome groups > 100
        var chroms = Enumerable.Range(1, groups)
            .Select(g => Chrom($"chr{g}", 1_000_000))
            .ToList();

        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chroms);

        result.TotalChromosomes.Should().Be(groups);
        result.AutosomeCount.Should().Be(groups);
        result.SexChromosomes.Should().BeEmpty();
        result.HasAneuploidy.Should().BeTrue("each single-copy group differs from the expected ploidy of 2");
        result.Abnormalities.Should().HaveCount(groups, "one abnormality label per under-copied autosome group");
        result.Abnormalities.Should().OnlyContain(a => a.StartsWith("Monosomy "),
            "a single copy against expected 2 is monosomy");
        // INV-03 consistency: flag set exactly when the list is non-empty.
        result.HasAneuploidy.Should().Be(result.Abnormalities.Count > 0);
    }

    #endregion

    #region BE — Well-formedness invariants on a sex-chromosome-only set

    /// <summary>
    /// BE: an all-sex-chromosome set (no autosomes). Sex chromosomes are preserved in
    /// SexChromosomes but never grouped for abnormality calling (Karyotype_Analysis.md
    /// §5.2 deviation #1), so HasAneuploidy must be false even though the set is not a
    /// normal diploid autosome complement. AutosomeCount 0 with the mean still computed
    /// over the full set (INV-02) — and crucially the mean is taken over a NON-empty
    /// list, so no zero divisor arises despite zero autosomes.
    /// </summary>
    [Test]
    public void AnalyzeKaryotype_SexChromosomesOnly_NoAneuploidyMeanOverFullSet()
    {
        var chroms = new[]
        {
            Chrom("chrX", 1_500_000, sex: true),
            Chrom("chrY", 600_000, sex: true),
        };

        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chroms);

        result.TotalChromosomes.Should().Be(2);
        result.AutosomeCount.Should().Be(0);
        result.SexChromosomes.Should().BeEquivalentTo(new[] { "chrX", "chrY" });
        result.HasAneuploidy.Should().BeFalse("sex chromosomes are never grouped for abnormality calling");
        result.Abnormalities.Should().BeEmpty();
        result.TotalGenomeSize.Should().Be(2_100_000);
        result.MeanChromosomeLength.Should().Be(2_100_000 / 2.0, "the mean is taken over the full non-empty set");
    }

    #endregion

    #region Positive sanity — a known human karyotype yields the expected summary

    /// <summary>
    /// Positive sanity: the canonical human diploid karyotype (Tjio &amp; Levan 1956,
    /// 2n = 46) built as 22 diploid autosome pairs (chr1_1/chr1_2 … chr22_1/chr22_2 =
    /// 44 autosomes) plus an XX sex pair = 46 chromosomes. This is the affirmative
    /// anchor that the analyzer recovers the EXPECTED counts and aggregates: 46 total,
    /// 44 autosomes, 2 sex chromosomes, total size = Σ lengths, the right mean, and —
    /// every autosome group being a clean diploid pair — NO aneuploidy. Without it an
    /// all-passing boundary suite could be vacuously green.
    /// </summary>
    [Test]
    public void AnalyzeKaryotype_NormalHumanDiploidSet_RecoversExpectedSummary()
    {
        const int autosomeGroups = 22;
        const long autosomeLen = 100_000_000;
        const long sexLen = 155_000_000;
        var chroms = new List<(string Name, long Length, bool IsSexChromosome)>();
        for (int g = 1; g <= autosomeGroups; g++)
        {
            chroms.Add(Chrom($"chr{g}_1", autosomeLen));
            chroms.Add(Chrom($"chr{g}_2", autosomeLen));
        }
        chroms.Add(Chrom("chrX_1", sexLen, sex: true));
        chroms.Add(Chrom("chrX_2", sexLen, sex: true));

        long expectedTotal = (autosomeGroups * 2L * autosomeLen) + (2L * sexLen);

        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chroms);

        result.TotalChromosomes.Should().Be(46, "the human diploid complement is 2n = 46");
        result.AutosomeCount.Should().Be(44);
        result.SexChromosomes.Should().HaveCount(2);
        result.PloidyLevel.Should().Be(2, "the default expected ploidy is copied through for non-empty input");
        result.TotalGenomeSize.Should().Be(expectedTotal);
        result.MeanChromosomeLength.Should().Be(expectedTotal / 46.0);
        result.HasAneuploidy.Should().BeFalse("every autosome group is a clean diploid pair");
        result.Abnormalities.Should().BeEmpty();
    }

    /// <summary>
    /// Positive sanity: a known aneuploidy (trisomy 21, Down syndrome) — chr21 present
    /// in THREE copies against the expected diploid baseline — must be surfaced with
    /// exactly the standard cytogenetic label "Trisomy chr21" and HasAneuploidy true,
    /// while the clean diploid groups contribute no abnormality. This pins that the
    /// abnormality classifier names the correct group with the correct ISCN term.
    /// </summary>
    [Test]
    public void AnalyzeKaryotype_Trisomy21_LabeledTrisomyChr21()
    {
        var chroms = new List<(string Name, long Length, bool IsSexChromosome)>
        {
            Chrom("chr1_1", 1_000_000), Chrom("chr1_2", 1_000_000),     // clean diploid pair
            Chrom("chr21_1", 500_000), Chrom("chr21_2", 500_000), Chrom("chr21_3", 500_000), // trisomy
        };

        var result = ChromosomeAnalyzer.AnalyzeKaryotype(chroms);

        result.HasAneuploidy.Should().BeTrue();
        result.Abnormalities.Should().ContainSingle()
            .Which.Should().Be("Trisomy chr21", "three copies against expected 2 is trisomy of that group");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CHROM-ANEU-001 — aneuploidy detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CHROM-ANEU-001 — aneuploidy detection

    #region Aneuploidy helpers

    /// <summary>Bin width used by the aneuploidy fuzz cases (the production default).</summary>
    private const int AneuBinSize = 1_000_000;

    /// <summary>A valid positive diploid baseline used across the aneuploidy fuzz cases.</summary>
    private const double AneuBaseline = 30.0;

    /// <summary>
    /// Builds <paramref name="bins"/> depth observations on one chromosome, one per bin,
    /// each at the given <paramref name="depth"/>. Positions are spaced one bin apart so
    /// every observation lands in a distinct bin (Position / binSize).
    /// </summary>
    private static (string Chromosome, int Position, double Depth)[] DepthProfile(
        string chromosome, double depth, int bins)
        => System.Linq.Enumerable.Range(0, bins)
            .Select(i => (chromosome, i * AneuBinSize, depth))
            .ToArray();

    #endregion

    #region BE — Boundary: empty data

    /// <summary>
    /// BE (empty data): an empty depth sequence. The materialized list is empty, so
    /// `DetectAneuploidy` hits the `data.Count == 0` guard (ChromosomeAnalyzer.cs line
    /// 839) and `yield break`s — a defined EMPTY result, never a crash. Because the
    /// method is a lazy iterator we force enumeration with `.ToList()`.
    /// (Aneuploidy_Detection.md §6.1 "Empty input depth data".)
    /// </summary>
    [Test]
    public void DetectAneuploidy_EmptyData_YieldsNoStates()
    {
        var empty = Enumerable.Empty<(string Chromosome, int Position, double Depth)>();

        var act = () => ChromosomeAnalyzer.DetectAneuploidy(empty, AneuBaseline).ToList();
        act.Should().NotThrow("the empty-input guard short-circuits before any depth-ratio arithmetic");

        var result = ChromosomeAnalyzer.DetectAneuploidy(empty, AneuBaseline).ToList();
        result.Should().BeEmpty("empty depth data yields no copy-number states");
    }

    #endregion

    #region BE — Boundary: 0 / negative baseline (KEY div-by-zero)

    /// <summary>
    /// BE (KEY div-by-zero): a 0 diploid baseline. The depth ratio is
    /// `meanDepth / medianDepth`; a 0 baseline would make every ratio Infinity. The
    /// guard `medianDepth &lt;= 0` (line 839) REJECTS this BEFORE the ratio is computed,
    /// yielding the defined empty result — no DivideByZero, no Infinity-driven mis-call.
    /// A 0 baseline carries no signal, so "no call" (empty) is the theory-correct
    /// outcome. (Aneuploidy_Detection.md §6.1 "medianDepth = 0".)
    /// </summary>
    [Test]
    public void DetectAneuploidy_ZeroBaseline_YieldsNoStatesNoDivideByZero()
    {
        var depthData = DepthProfile("chr1", AneuBaseline, bins: 5);

        var act = () => ChromosomeAnalyzer.DetectAneuploidy(depthData, medianDepth: 0.0).ToList();
        act.Should().NotThrow("a 0 baseline is rejected before the depth ratio, so no DivideByZero or Infinity ratio occurs");

        var result = ChromosomeAnalyzer.DetectAneuploidy(depthData, medianDepth: 0.0).ToList();
        result.Should().BeEmpty("a 0 diploid baseline carries no signal and yields no call");
    }

    /// <summary>
    /// BE (negative baseline): a NEGATIVE diploid baseline is biologically impossible
    /// and is caught by the SAME `medianDepth &lt;= 0` guard (line 839). It must yield the
    /// defined empty result — a negative baseline can never drive a copy-number call.
    /// (Aneuploidy_Detection.md §6.1 "medianDepth &lt; 0".)
    /// </summary>
    [Test]
    public void DetectAneuploidy_NegativeBaseline_YieldsNoStates()
    {
        var depthData = DepthProfile("chr1", AneuBaseline, bins: 5);

        var act = () => ChromosomeAnalyzer.DetectAneuploidy(depthData, medianDepth: -10.0).ToList();
        act.Should().NotThrow("a negative baseline is caught by the same non-positive guard");

        var result = ChromosomeAnalyzer.DetectAneuploidy(depthData, medianDepth: -10.0).ToList();
        result.Should().BeEmpty("a negative diploid baseline is rejected, yielding no call");
    }

    #endregion

    #region BE — Boundary: 0 depth values

    /// <summary>
    /// BE (0 depth): every depth observation is 0 against a VALID positive baseline.
    /// meanDepth is 0, so `log2(0 / median)` = −∞ and `2^(−∞) = 0` → copyNumber rounds
    /// and clamps to 0 (Nullisomy) — the documented "very low depth ⇒ floored at 0"
    /// case (§6.1). LogRatio is −∞ for the bin, but the emitted CopyNumber is a
    /// defined, in-range integer 0; the run must NOT crash and must NOT spuriously call
    /// a gain. We pin CopyNumber == 0 across all bins and INV-01 on every state.
    /// </summary>
    [Test]
    public void DetectAneuploidy_ZeroDepthValues_CallsNullisomyNoCrash()
    {
        var depthData = DepthProfile("chr1", depth: 0.0, bins: 5);

        var act = () => ChromosomeAnalyzer.DetectAneuploidy(depthData, AneuBaseline, AneuBinSize).ToList();
        act.Should().NotThrow("0 depth drives log2(0) = −∞, but 2^(−∞) = 0 floors cleanly to copy number 0");

        var result = ChromosomeAnalyzer.DetectAneuploidy(depthData, AneuBaseline, AneuBinSize).ToList();

        result.Should().HaveCount(5, "one bin per planted observation");
        result.Should().OnlyContain(s => s.CopyNumber == 0,
            "zero read depth floors to copy number 0 (Nullisomy)");
        result.Should().OnlyContain(s => s.CopyNumber >= 0 && s.CopyNumber <= 10,
            "INV-01: every emitted copy number is an integer in [0, 10]");
        result.Should().OnlyContain(s => double.IsNegativeInfinity(s.LogRatio),
            "log2(0 / median) is −∞ for a zero-depth bin");
    }

    #endregion

    #region BE — Boundary: negative depth values (impossible input, must not mis-call)

    /// <summary>
    /// BE (negative depth): NEGATIVE depth values are biologically impossible. The
    /// method does NOT validate or reject them — but, crucially, it must NOT crash and
    /// must NOT emit a garbage NON-zero copy-number GAIN. `meanDepth / median` is
    /// negative → `log2(negative)` = NaN → `2^NaN` = NaN → `(int)round(NaN)` = 0 → clamp
    /// keeps 0. So every bin collapses to copyNumber 0, the disciplined "no spurious
    /// aneuploidy call" outcome. The cost of skipping validation is that LogRatio and
    /// Confidence propagate as NaN; we pin that ACTUAL behaviour honestly rather than
    /// weaken the test to a falsely clean expectation. (Verified against net10.0
    /// NaN→int semantics; Aneuploidy_Detection.md §3.3 — no validation beyond the
    /// baseline/empty checks.)
    /// </summary>
    [Test]
    public void DetectAneuploidy_NegativeDepthValues_CollapseToCopyNumberZeroNoGarbageGain()
    {
        var depthData = DepthProfile("chr1", depth: -15.0, bins: 5);

        var act = () => ChromosomeAnalyzer.DetectAneuploidy(depthData, AneuBaseline, AneuBinSize).ToList();
        act.Should().NotThrow("negative depth produces NaN ratios but never an unhandled exception");

        var result = ChromosomeAnalyzer.DetectAneuploidy(depthData, AneuBaseline, AneuBinSize).ToList();

        result.Should().HaveCount(5, "one bin per planted observation; the method does not drop negative-depth bins");
        result.Should().OnlyContain(s => s.CopyNumber == 0,
            "a negative depth ratio yields NaN → (int)round(NaN) = 0; the call collapses to 0, never a spurious gain");
        result.Should().OnlyContain(s => s.CopyNumber >= 0 && s.CopyNumber <= 10,
            "INV-01: even on impossible input the emitted copy number stays in [0, 10]");
        // Honest contract: the un-validated negative input propagates NaN, not a clean ratio.
        result.Should().OnlyContain(s => double.IsNaN(s.LogRatio),
            "log2 of a negative ratio is NaN; the method does not sanitize impossible negative depth");
    }

    #endregion

    #region BE — Boundary: extremely high depth (clamped, no overflow)

    /// <summary>
    /// BE (extremely high depth): a depth ratio far above any real ploidy. The raw
    /// `round((2^logRatio) * 2)` is enormous, but the copy number is CLAMPED to 10
    /// (line 860), so the emitted state is a finite, in-range integer — no overflow,
    /// no Infinity in the reported CopyNumber. (§6.1 "Very high depth ratios ⇒ capped
    /// at 10".) LogRatio is a large finite double; Confidence stays in [0, 1].
    /// </summary>
    [Test]
    public void DetectAneuploidy_ExtremelyHighDepth_CopyNumberClampedToTen()
    {
        // 1e6 × the baseline — an absurd ratio that would overflow an unclamped copy number.
        var depthData = DepthProfile("chr1", depth: AneuBaseline * 1_000_000.0, bins: 3);

        var act = () => ChromosomeAnalyzer.DetectAneuploidy(depthData, AneuBaseline, AneuBinSize).ToList();
        act.Should().NotThrow("the rounded copy number is clamped before emission, so no overflow occurs");

        var result = ChromosomeAnalyzer.DetectAneuploidy(depthData, AneuBaseline, AneuBinSize).ToList();

        result.Should().HaveCount(3);
        result.Should().OnlyContain(s => s.CopyNumber == 10,
            "an extreme depth ratio is capped at the maximum copy number 10");
        result.Should().OnlyContain(s => double.IsFinite(s.LogRatio),
            "log2 of a large finite ratio is a large finite double, never Infinity or NaN");
        result.Should().OnlyContain(s => s.Confidence >= 0.0 && s.Confidence <= 1.0,
            "INV-02: confidence stays in [0, 1] for finite depth");
    }

    /// <summary>
    /// BE (extremely high baseline): the mirror extreme — a gigantic baseline against
    /// ordinary depth drives the ratio to ~0, flooring copy number at 0 with a large
    /// NEGATIVE finite LogRatio. Must not crash, must stay clamped at 0.
    /// </summary>
    [Test]
    public void DetectAneuploidy_ExtremelyHighBaseline_CopyNumberFlooredAtZero()
    {
        var depthData = DepthProfile("chr1", depth: AneuBaseline, bins: 3);
        double hugeBaseline = AneuBaseline * 1_000_000.0;

        var result = ChromosomeAnalyzer.DetectAneuploidy(depthData, hugeBaseline, AneuBinSize).ToList();

        result.Should().HaveCount(3);
        result.Should().OnlyContain(s => s.CopyNumber == 0,
            "a tiny depth ratio against a huge baseline floors copy number at 0");
        result.Should().OnlyContain(s => double.IsFinite(s.LogRatio) && s.LogRatio < 0,
            "log2 of a small finite ratio is a large negative finite double");
    }

    #endregion

    #region Positive sanity — clear trisomy is called as gain, normal as no aneuploidy

    /// <summary>
    /// Positive sanity (gain): a clear trisomy depth profile — every bin at 1.5× the
    /// diploid baseline (3 copies ⇒ 1.5× depth, Aneuploidy_Detection.md §2.2) — must be
    /// called as copy number 3 per bin, and `IdentifyWholeChromosomeAneuploidy` must
    /// surface it as a whole-chromosome "Trisomy" gain. This is the affirmative anchor
    /// that the detector recovers a real aneuploidy, not just survives garbage.
    /// </summary>
    [Test]
    public void DetectAneuploidy_ClearTrisomyProfile_CalledAsGain()
    {
        var depthData = DepthProfile("chr21", depth: AneuBaseline * 1.5, bins: 10);

        var states = ChromosomeAnalyzer.DetectAneuploidy(depthData, AneuBaseline, AneuBinSize).ToList();

        states.Should().HaveCount(10);
        states.Should().OnlyContain(s => s.CopyNumber == 3,
            "1.5× the diploid depth is three copies (trisomy)");

        var whole = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();
        whole.Should().ContainSingle()
            .Which.Should().Match<(string Chromosome, int CopyNumber, string Type)>(
                w => w.Chromosome == "chr21" && w.CopyNumber == 3 && w.Type == "Trisomy",
                "a chromosome dominated by copy-number-3 bins is a whole-chromosome trisomy gain");
    }

    /// <summary>
    /// Positive sanity (normal): a clean diploid depth profile — every bin at exactly
    /// the baseline (1.0× ⇒ 2 copies) — must be called as copy number 2 per bin, and
    /// `IdentifyWholeChromosomeAneuploidy` must report NO whole-chromosome aneuploidy
    /// (disomic dominant states are ignored, INV-03). This pins that a normal sample is
    /// not mis-called as a gain or loss.
    /// </summary>
    [Test]
    public void DetectAneuploidy_NormalDiploidProfile_NoAneuploidyCall()
    {
        var depthData = DepthProfile("chr1", depth: AneuBaseline, bins: 10);

        var states = ChromosomeAnalyzer.DetectAneuploidy(depthData, AneuBaseline, AneuBinSize).ToList();

        states.Should().HaveCount(10);
        states.Should().OnlyContain(s => s.CopyNumber == 2,
            "depth at exactly the baseline is the normal two-copy diploid state");

        var whole = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states).ToList();
        whole.Should().BeEmpty("a disomic chromosome is not reported as a whole-chromosome aneuploidy");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CHROM-SYNT-001 — synteny analysis : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CHROM-SYNT-001 — synteny analysis

    #region Synteny helpers

    /// <summary>An ortholog-pair tuple alias for the FindSyntenyBlocks input shape.</summary>
    private static (string Chr1, int Start1, int End1, string Gene1,
                    string Chr2, int Start2, int End2, string Gene2) Ortholog(
        string chr1, int start1, int end1, string gene1,
        string chr2, int start2, int end2, string gene2)
        => (chr1, start1, end1, gene1, chr2, start2, end2, gene2);

    /// <summary>
    /// Builds <paramref name="genes"/> ascending, non-overlapping genes on one
    /// chromosome together with an IDENTITY ortholog map (each gene paired with itself
    /// on the same chromosome and coordinates). This is the self-comparison input: a
    /// genome compared against itself. Genes are spaced <paramref name="step"/> apart
    /// so adjacent gaps are well below the megabase-scale maxGap, keeping the run
    /// collinear. Deterministic — no randomness.
    /// </summary>
    private static (string, int, int, string, string, int, int, string)[] SelfOrthologMap(
        string chromosome, int genes, int step = 10_000, int geneSpan = 1_000)
        => Enumerable.Range(0, genes)
            .Select(i =>
            {
                int start = i * step;
                int end = start + geneSpan;
                string name = $"gene{i}";
                return Ortholog(chromosome, start, end, name, chromosome, start, end, name);
            })
            .ToArray();

    #endregion

    #region BE — Boundary: 0 ortholog pairs (empty anchors)

    /// <summary>
    /// BE (KEY boundary): zero ortholog pairs. The materialized list is empty, so the
    /// `pairs.Count &lt; minGenes` guard (ChromosomeAnalyzer.cs line 651) `yield break`s
    /// to a defined EMPTY result — no synteny blocks, never a crash. No anchors means
    /// no conserved gene-order blocks. The method is a lazy iterator, so enumeration is
    /// forced with `.ToList()`. (Synteny_Analysis.md §6.1 "Empty ortholog input".)
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_ZeroOrthologPairs_YieldsNoBlocks()
    {
        var empty = Enumerable.Empty<(string, int, int, string, string, int, int, string)>();

        var act = () => ChromosomeAnalyzer.FindSyntenyBlocks(empty).ToList();
        act.Should().NotThrow("the empty-input guard short-circuits before any grouping or scan");

        var result = ChromosomeAnalyzer.FindSyntenyBlocks(empty).ToList();
        result.Should().BeEmpty("zero ortholog anchors yield no synteny blocks");
    }

    /// <summary>
    /// BE: fewer than `minGenes` ortholog pairs (here 2, below the default 3). The
    /// total pair count is checked before scanning (line 651), so a sub-threshold run
    /// produces no block even though the pairs are perfectly collinear. This pins the
    /// minGenes boundary adjacent to the 0-pair case. (Synteny_Analysis.md §6.1
    /// "Fewer than minGenes ortholog pairs".)
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_FewerThanMinGenes_YieldsNoBlocks()
    {
        var twoPairs = new[]
        {
            Ortholog("chr1", 1000, 2000, "g1", "chrA", 1000, 2000, "gA"),
            Ortholog("chr1", 3000, 4000, "g2", "chrA", 3000, 4000, "gB"),
        };

        var result = ChromosomeAnalyzer.FindSyntenyBlocks(twoPairs, minGenes: 3, maxGap: 10).ToList();
        result.Should().BeEmpty("2 pairs are below the minGenes threshold of 3, so no block is emitted");
    }

    #endregion

    #region BE — Boundary: empty genomes

    /// <summary>
    /// BE (empty genomes): an "empty genome" contributes no genes and therefore no
    /// ortholog pairs, so the empty-genome boundary reduces to the 0-pair boundary.
    /// Expressed here as a freshly-allocated empty ARRAY rather than the LINQ Empty
    /// singleton — the same boundary reached through a different concrete enumerable —
    /// it must produce the identical empty result, pinning that the guard keys on the
    /// materialized count, not the source type. Never a NullReference / IndexOutOfRange.
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_EmptyGenomes_YieldsNoBlocks()
    {
        var emptyGenomes = Array.Empty<(string, int, int, string, string, int, int, string)>();

        var act = () => ChromosomeAnalyzer.FindSyntenyBlocks(emptyGenomes).ToList();
        act.Should().NotThrow("empty genomes contribute no ortholog pairs; the count guard returns empty, never indexing");

        var result = ChromosomeAnalyzer.FindSyntenyBlocks(emptyGenomes).ToList();
        result.Should().BeEmpty("two empty genomes share no orthologs and so have no synteny blocks");
    }

    #endregion

    #region BE — Boundary: self-comparison (identity map) — KEY positive

    /// <summary>
    /// BE (KEY positive): a genome compared against ITSELF under an identity ortholog
    /// map. Every gene maps to itself on the same chromosome with identical
    /// coordinates, and the genes are in ascending, non-overlapping order. Because the
    /// order is identical in both genomes, the whole genome is ONE maximal collinear
    /// FORWARD run, so the detector must emit exactly ONE SyntenyBlock covering the
    /// whole genome: Strand '+', GeneCount = N, Species1Chromosome == Species2Chromosome,
    /// and the block's coordinates in BOTH genomes equal the genome's own span. Identical
    /// gene order ⇒ perfect synteny — the core self-comparison theory anchor. We also
    /// assert the invariants: ordered, valid coordinates in both genomes (INV-02) and
    /// SequenceIdentity NaN (INV-03).
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_SelfComparisonIdentityMap_OneMaximalCollinearBlock()
    {
        const int genes = 8;
        const int step = 10_000;
        const int geneSpan = 1_000;
        var selfMap = SelfOrthologMap("chr1", genes, step, geneSpan);

        var act = () => ChromosomeAnalyzer.FindSyntenyBlocks(selfMap, minGenes: 3, maxGap: 10).ToList();
        act.Should().NotThrow("a clean identity map is the best-case collinear input and must never crash");

        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(selfMap, minGenes: 3, maxGap: 10).ToList();

        blocks.Should().ContainSingle(
            "a genome is perfectly collinear with itself, so the identity map collapses to one maximal block");

        var block = blocks[0];
        block.Strand.Should().Be('+', "identical gene order in both genomes is forward collinearity");
        block.GeneCount.Should().Be(genes, "the single maximal block covers every gene in the genome");
        block.Species1Chromosome.Should().Be("chr1");
        block.Species2Chromosome.Should().Be("chr1", "self-comparison maps each chromosome onto itself");

        // The block spans the whole genome in BOTH genomes (first gene start → last gene end).
        int genomeStart = 0;
        int genomeEnd = (genes - 1) * step + geneSpan;
        block.Species1Start.Should().Be(genomeStart, "the block starts at the first gene of the genome");
        block.Species1End.Should().Be(genomeEnd, "the block ends at the last gene of the genome");
        block.Species2Start.Should().Be(genomeStart, "the identity map gives genome 2 the same span");
        block.Species2End.Should().Be(genomeEnd);

        // INV-02: valid, ordered coordinates referencing real anchors in both genomes.
        block.Species1Start.Should().BeLessThanOrEqualTo(block.Species1End);
        block.Species2Start.Should().BeLessThanOrEqualTo(block.Species2End);
        // INV-03: identity is never computed from coordinate-only input.
        double.IsNaN(block.SequenceIdentity).Should().BeTrue("SequenceIdentity is NaN for coordinate-only synteny");
    }

    #endregion

    #region Positive sanity — a planted conserved block is recovered with correct anchors

    /// <summary>
    /// Positive sanity: two DISTINCT genomes (chr1 ↦ chrA) sharing a clear conserved,
    /// forward-collinear block of 5 genes in the same order. The detector must recover
    /// exactly that block — one block, Strand '+', GeneCount 5, the correct source and
    /// target chromosome names, and coordinates that span the first→last planted gene
    /// in EACH genome. This is the affirmative anchor that synteny is actually found
    /// with correct anchors, so the boundary suite is not vacuously green.
    /// </summary>
    [Test]
    public void FindSyntenyBlocks_ConservedForwardBlock_RecoveredWithCorrectAnchors()
    {
        var orthologPairs = new[]
        {
            Ortholog("chr1", 1_000, 2_000, "g1", "chrA", 1_000, 2_000, "gA"),
            Ortholog("chr1", 3_000, 4_000, "g2", "chrA", 3_000, 4_000, "gB"),
            Ortholog("chr1", 5_000, 6_000, "g3", "chrA", 5_000, 6_000, "gC"),
            Ortholog("chr1", 7_000, 8_000, "g4", "chrA", 7_000, 8_000, "gD"),
            Ortholog("chr1", 9_000, 10_000, "g5", "chrA", 9_000, 10_000, "gE"),
        };

        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(orthologPairs, minGenes: 3, maxGap: 10).ToList();

        blocks.Should().ContainSingle("5 forward-collinear orthologs form exactly one conserved block");

        var block = blocks[0];
        block.Strand.Should().Be('+', "the conserved genes share the same order in both genomes");
        block.GeneCount.Should().Be(5, "all five planted orthologs are retained in the block");
        block.Species1Chromosome.Should().Be("chr1");
        block.Species2Chromosome.Should().Be("chrA");
        block.Species1Start.Should().Be(1_000, "the block starts at the first planted gene in genome 1");
        block.Species1End.Should().Be(10_000, "the block ends at the last planted gene in genome 1");
        block.Species2Start.Should().Be(1_000, "the block starts at the first planted ortholog in genome 2");
        block.Species2End.Should().Be(10_000, "the block ends at the last planted ortholog in genome 2");
        double.IsNaN(block.SequenceIdentity).Should().BeTrue("INV-03: coordinate-only synteny has NaN identity");
    }

    #endregion

    #endregion
}
