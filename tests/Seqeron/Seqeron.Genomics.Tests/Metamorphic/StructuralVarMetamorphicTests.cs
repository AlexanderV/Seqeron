namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Metamorphic tests for the StructuralVar area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SV-BREAKPOINT-001 — split-read breakpoint detection (StructuralVar).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 201.
///
/// API under test (StructuralVariantAnalyzer.FindBreakpoints):
///   Clusters split-read junctions within a tolerance window and reports a breakpoint per cluster,
///   with support = cluster size and quality growing with support.
///
/// Relations (derived from the clustering/support model, NOT from output):
///   • SHIFT (prepend flank shifts breakpoints): shifting every junction coordinate shifts the
///          reported breakpoint position by the same offset.
///   • MON  (more split reads ⇒ ≥ confidence): a cluster with more supporting reads has a
///          greater-or-equal quality score.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class StructuralVarMetamorphicTests
{
    private static StructuralVariantAnalyzer.SplitRead Sr(string id, int junction) =>
        new(id, "chr1", 500, junction, 20, "ACGTACGTAC");

    #region SV-BREAKPOINT-001 MON — more split reads raise the confidence

    [Test]
    [Description("MON: a breakpoint cluster with more supporting split reads has a greater-or-equal quality score.")]
    public void Breakpoints_MoreSplitReads_HigherConfidence()
    {
        double previousQuality = -1;
        int previousSupport = -1;
        foreach (int k in new[] { 2, 3, 5 })
        {
            var reads = Enumerable.Range(0, k).Select(i => Sr($"r{i}", 1000)).ToList();
            var bp = StructuralVariantAnalyzer.FindBreakpoints(reads).Single();

            bp.SupportingReads.Should().Be(k, because: "all reads cluster at the same junction");
            bp.SupportingReads.Should().BeGreaterThan(previousSupport, because: "more reads were added");
            bp.Quality.Should().BeGreaterThanOrEqualTo(previousQuality, because: "more supporting reads cannot lower the breakpoint quality");
            previousQuality = bp.Quality;
            previousSupport = bp.SupportingReads;
        }
    }

    #endregion

    #region SV-BREAKPOINT-001 SHIFT — a coordinate shift shifts the breakpoint

    [Test]
    [Description("SHIFT: shifting every split-read junction coordinate by an offset shifts the reported breakpoint position by the same offset.")]
    public void Breakpoints_CoordinateShift_ShiftsPosition()
    {
        var reads = new[] { Sr("a", 1000), Sr("b", 1001), Sr("c", 1002) };
        int originalPosition = StructuralVariantAnalyzer.FindBreakpoints(reads).Single().Position1;

        foreach (int offset in new[] { 1000, 100000 })
        {
            var shifted = reads.Select(r => Sr(r.ReadId, r.SupplementaryPosition + offset)).ToList();
            StructuralVariantAnalyzer.FindBreakpoints(shifted).Single().Position1
                .Should().Be(originalPosition + offset,
                    because: $"shifting every junction by {offset} shifts the consensus breakpoint by {offset}");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SV-CNV-001 — read-depth copy-number detection (StructuralVar).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 202.
    //
    // API under test (StructuralVariantAnalyzer.DetectCNV):
    //   Per window, CN = round(2 · 2^log2(windowDepth / referenceDepth)).
    //
    // Relations (derived from the depth-ratio model, NOT from output):
    //   • MON  (higher coverage ratio ⇒ higher CN): with a fixed reference, a window of higher read
    //          depth has a greater-or-equal integer copy number.
    //   • INV  (bin order independent): each window's CN depends only on its own depth and the fixed
    //          reference, so reordering the windows preserves the multiset of copy numbers.
    // ───────────────────────────────────────────────────────────────────────────

    private const int CnvWindow = 4;
    private const double CnvReference = 10.0;

    private static int[] DepthFromWindowDepths(params int[] windowDepths) =>
        windowDepths.SelectMany(d => Enumerable.Repeat(d, CnvWindow)).ToArray();

    #region SV-CNV-001 MON — higher coverage ratio gives higher copy number

    [Test]
    [Description("MON: with a fixed reference depth, windows of increasing read depth yield a non-decreasing integer copy number.")]
    public void Cnv_HigherDepth_HigherCopyNumber()
    {
        var depth = DepthFromWindowDepths(5, 10, 20, 40); // ratios 0.5, 1, 2, 4
        var segments = StructuralVariantAnalyzer.DetectCNV(depth, CnvWindow, CnvReference).ToList();

        segments.Should().HaveCount(4);
        int previous = -1;
        foreach (var seg in segments)
        {
            seg.CopyNumber.Should().BeGreaterThanOrEqualTo(previous, because: "a higher depth-to-reference ratio cannot give a lower copy number");
            previous = seg.CopyNumber;
        }
        segments.Last().CopyNumber.Should().BeGreaterThan(segments.First().CopyNumber, because: "the deepest window has a strictly higher copy number than the shallowest");
    }

    #endregion

    #region SV-CNV-001 INV — copy numbers are independent of bin order

    [Test]
    [Description("INV: each window's copy number depends only on its own depth and the fixed reference, so reordering the windows preserves the multiset of copy numbers.")]
    public void Cnv_BinOrder_Invariant()
    {
        int[] windowDepths = { 5, 10, 20, 40 };
        var forward = StructuralVariantAnalyzer.DetectCNV(DepthFromWindowDepths(windowDepths), CnvWindow, CnvReference)
            .Select(s => s.CopyNumber).OrderBy(x => x).ToList();
        var reordered = StructuralVariantAnalyzer.DetectCNV(DepthFromWindowDepths(windowDepths.Reverse().ToArray()), CnvWindow, CnvReference)
            .Select(s => s.CopyNumber).OrderBy(x => x).ToList();

        reordered.Should().Equal(forward, because: "the set of per-window copy numbers does not depend on window order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SV-DETECT-001 — paired-end structural-variant detection (StructuralVar).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 203.
    //
    // API under test (StructuralVariantAnalyzer.DetectSVs — the canonical SV entry point):
    //   Flags discordant read pairs (interchromosomal, span outside mean ± c·sd, or non-FR
    //   orientation — Medvedev et al. 2009 PEM signatures), clusters nearby concordant evidence,
    //   and emits one StructuralVariant per cluster meeting the minimum read-pair support. The
    //   emitted Start = min(mate-1 position) and End = max(mate-2 position) over the cluster.
    //
    // Relations (derived from the discordance + clustering model, NOT from output):
    //   • INV  (identical genomes ⇒ no SV): a read-pair set in which every pair maps concordantly
    //          (same chromosome, forward-reverse orientation, insert size within mean ± c·sd) is
    //          exactly what a sample identical to the reference produces — there is no rearrangement
    //          signature, so DetectSVs reports nothing. The positive control shows the same harness
    //          DOES call an SV once a genuine deletion-span signature is introduced, so the empty
    //          result is the model rejecting non-evidence, not a degenerate always-empty detector.
    //   • SHIFT (coordinate shift shifts SVs): discordance depends only on insert size, orientation
    //          and chromosome, and clustering only on coordinate GAPS — all invariant under a uniform
    //          translation of every mate coordinate. Hence shifting every read-pair position by an
    //          offset shifts the called SV's Start/End by exactly that offset while leaving its type,
    //          length, supporting-read count and quality unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    // BreakDancer-style library model used throughout SV-DETECT-001: mean insert 400, sd 50, c = 3
    // ⇒ a span is concordant iff it lies in [250, 550]. These are the DetectSVs defaults.
    private const int ExpectedInsert = 400;
    private const int InsertSd = 50;
    private const int ConcordantSpan = 400;   // inside [250, 550]
    private const int DeletionSpan = 1200;    // > 550 ⇒ deletion-span signature (FR mates)

    private static (string, string, int, char, string, int, char, int) Pair(
        string id, int pos1, int pos2, int insertSize, char strand1 = '+', char strand2 = '-',
        string chr1 = "chr1", string chr2 = "chr1") =>
        (id, chr1, pos1, strand1, chr2, pos2, strand2, insertSize);

    #region SV-DETECT-001 INV — a genome identical to the reference yields no structural variants

    [Test]
    [Description("INV: a fully concordant read-pair set (the signature of a sample identical to the reference) produces no SVs, whereas the same harness calls an SV once a genuine deletion-span signature is present.")]
    public void DetectSvs_IdenticalGenome_NoStructuralVariants()
    {
        // Concordant pairs: same chromosome, forward-reverse orientation, span inside mean ± 3·sd.
        // This is precisely what aligning a sample equal to the reference yields — no rearrangement.
        var concordant = new[]
        {
            Pair("c1", 1_000, 1_400, ConcordantSpan),
            Pair("c2", 1_010, 1_410, ConcordantSpan),
            Pair("c3", 5_000, 5_400, ConcordantSpan),
            Pair("c4", 5_020, 5_420, ConcordantSpan),
        };

        StructuralVariantAnalyzer.DetectSVs(concordant, ExpectedInsert, InsertSd)
            .Should().BeEmpty(because: "a genome identical to the reference produces only concordant pairs, so there is no structural-variant signature to call");

        // Positive control — replacing the spans with a deletion-span signature (everything else equal)
        // must yield at least one SV, proving the empty result above is evidence-driven, not degenerate.
        var deletionEvidence = new[]
        {
            Pair("d1", 1_000, 2_200, DeletionSpan),
            Pair("d2", 1_010, 2_210, DeletionSpan),
        };

        StructuralVariantAnalyzer.DetectSVs(deletionEvidence, ExpectedInsert, InsertSd)
            .Should().NotBeEmpty(because: "an over-long FR span supported by ≥ minSupport pairs is a deletion signature the detector must call");
    }

    #endregion

    #region SV-DETECT-001 SHIFT — translating every read-pair coordinate shifts the called SVs

    [Test]
    [Description("SHIFT: shifting every mate coordinate by an offset shifts each called SV's Start/End by the same offset while preserving its type, length, support and quality.")]
    public void DetectSvs_CoordinateShift_ShiftsStructuralVariants()
    {
        var reads = new[]
        {
            Pair("d1", 1_000, 2_200, DeletionSpan),
            Pair("d2", 1_010, 2_210, DeletionSpan),
            Pair("d3", 1_020, 2_220, DeletionSpan),
        };

        var baseline = StructuralVariantAnalyzer.DetectSVs(reads, ExpectedInsert, InsertSd).ToList();
        baseline.Should().ContainSingle(because: "the three pairs cluster into one deletion event");

        foreach (int offset in new[] { 10_000, 1_000_000 })
        {
            var shifted = reads
                .Select(r => Pair(r.Item1, r.Item3 + offset, r.Item6 + offset, r.Item8))
                .ToArray();

            var shiftedSvs = StructuralVariantAnalyzer.DetectSVs(shifted, ExpectedInsert, InsertSd).ToList();

            shiftedSvs.Should().HaveCount(baseline.Count, because: "a uniform coordinate translation does not change which pairs are discordant or how they cluster");

            for (int i = 0; i < baseline.Count; i++)
            {
                var b = baseline[i];
                var s = shiftedSvs[i];

                s.Start.Should().Be(b.Start + offset, because: $"the SV start tracks the shifted mate-1 minimum by {offset}");
                s.End.Should().Be(b.End + offset, because: $"the SV end tracks the shifted mate-2 maximum by {offset}");
                s.Type.Should().Be(b.Type, because: "the PEM signature (insert size and orientation) is unchanged by a coordinate shift");
                s.Length.Should().Be(b.Length, because: "translation preserves the span between the SV endpoints");
                s.SupportingReads.Should().Be(b.SupportingReads, because: "the same pairs cluster together after the shift");
                s.Quality.Should().Be(b.Quality, because: "quality is a function of support, which the shift preserves");
            }
        }
    }

    #endregion
}
