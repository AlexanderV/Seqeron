namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Transcriptome area (TranscriptomeAnalyzer,
/// Seqeron.Genomics.Annotation).
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of combinatorial testing.
/// Each grid cell carries a real business assertion; small grids use the exhaustive
/// <c>[Combinatorial]</c> product (a strict superset of pairwise).
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Transcriptome")]
public class TranscriptomeCombinatorialTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: TRANS-DIFF-001 — Differential expression (Transcriptome)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 198.
    // Spec: tests/TestSpecs/TRANS-DIFF-001.md (canonical FindDifferentiallyExpressed / CalculateFoldChange).
    // ADVANCED §10.
    // Dimensions: nReplicates(3) × foldChange(3) × test(2). Grid 3×3×2 = 18 (full, exhaustive).
    //
    // Model (Love 2014 DESeq2; Welch 1947; Benjamini-Hochberg 1995): log2 fold change = log2((mean2+c)/(mean1+c))
    // (positive ⇒ up in condition 2); a gene is DE iff BOTH |log2FC| ≥ threshold AND BH-adjusted Welch
    // t-test p < alpha.
    //
    // Axis mapping (documented): nReplicates → replicates per condition; foldChange → the planted
    // direction (none/up/down); test → the |log2FC| threshold (strict 1.0 vs lenient 0.0). The
    // combinatorial point: the two-criterion DE gate holds exactly (IsSignificant ⇔ both criteria),
    // and the Regulation label tracks the sign of the fold change, at every cell.
    // ═══════════════════════════════════════════════════════════════════════

    public enum FoldDirection { None, Up, Down }

    [Test, Combinatorial]
    public void TransDiff_TwoCriterionGate_AcrossReplicatesFoldChangeThreshold(
        [Values(3, 4, 5)] int nReplicates,
        [Values(FoldDirection.None, FoldDirection.Up, FoldDirection.Down)] FoldDirection fold,
        [Values(1.0, 0.0)] double log2FcThreshold)
    {
        (double Base1, double Base2) = fold switch
        {
            FoldDirection.Up => (10.0, 40.0),    // ~4× up   ⇒ log2FC ≈ +2
            FoldDirection.Down => (40.0, 10.0),  // ~4× down ⇒ log2FC ≈ −2
            _ => (10.0, 10.0),                   // no change
        };
        // Low-variance replicates (small jitter) so the t-test is well-defined.
        var c1 = Enumerable.Range(0, nReplicates).Select(i => Base1 + 0.1 * i).ToList();
        var c2 = Enumerable.Range(0, nReplicates).Select(i => Base2 + 0.1 * i).ToList();

        var de = TranscriptomeAnalyzer.FindDifferentiallyExpressed(
            new[] { ("g", (IReadOnlyList<double>)c1, (IReadOnlyList<double>)c2) }, alpha: 0.05, log2FoldChangeThreshold: log2FcThreshold)
            .Single();

        de.Log2FoldChange.Should().BeApproximately(TranscriptomeAnalyzer.CalculateFoldChange(c1, c2), 1e-9);

        string expectedReg = de.Log2FoldChange > 0 ? "Upregulated" : (de.Log2FoldChange < 0 ? "Downregulated" : "Unchanged");
        de.Regulation.Should().Be(expectedReg, "regulation tracks the sign of the fold change");

        de.IsSignificant.Should().Be(Math.Abs(de.Log2FoldChange) >= log2FcThreshold && de.AdjustedPValue < 0.05,
            "DE ⇔ |log2FC| ≥ threshold AND adjusted p < alpha");

        if (fold == FoldDirection.None)
            de.IsSignificant.Should().BeFalse("an unchanged gene is never differentially expressed");
    }

    /// <summary>
    /// Interaction witness — a strong, low-variance fold change is DE at the strict threshold; the
    /// fold-change sign distinguishes up- from down-regulation.
    /// </summary>
    [Test]
    public void TransDiff_StrongFoldChangeIsSignificant_SignSetsDirection()
    {
        var c1 = new List<double> { 10, 10.1, 10.2, 10.3 };
        var up = new List<double> { 40, 40.1, 40.2, 40.3 };

        var deUp = TranscriptomeAnalyzer.FindDifferentiallyExpressed(
            new[] { ("g", (IReadOnlyList<double>)c1, (IReadOnlyList<double>)up) }).Single();
        deUp.IsSignificant.Should().BeTrue("a ~4× change with low variance and 4 replicates is DE");
        deUp.Log2FoldChange.Should().BeGreaterThan(0);

        TranscriptomeAnalyzer.CalculateFoldChange(up, c1).Should().BeLessThan(0, "swapping conditions flips the sign");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: TRANS-EXPR-001 — Expression quantification & normalization (Transcriptome)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 199.
    // Spec: tests/TestSpecs/TRANS-EXPR-001.md (canonical CalculateTPM / CalculateFPKM / QuantileNormalize).
    // ADVANCED §10.
    // Dimensions: nReads(3) × normalization(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Wagner 2012 TPM; Mortazavi 2008 FPKM; Bolstad 2003 quantile): TPM sums to 1e6 within a
    // sample and is sequencing-depth-independent; FPKM = X·1e9/(l·N); quantile normalization makes
    // every sample share one common value distribution.
    //
    // Axis mapping (documented): normalization → which method (TPM/FPKM/Quantile); nReads → the
    // library-size scale. The combinatorial point: TPM sums to 1e6 and is invariant to read depth,
    // FPKM equals its closed form, and quantile-normalized samples are identically distributed —
    // each verified across the read-depth axis.
    // ═══════════════════════════════════════════════════════════════════════

    public enum NormMethod { Tpm, Fpkm, Quantile }

    [Test, Combinatorial]
    public void TransExpr_NormalizationInvariants_AcrossDepthAndMethod(
        [Values(1, 10, 100)] int depthScale,
        [Values(NormMethod.Tpm, NormMethod.Fpkm, NormMethod.Quantile)] NormMethod normalization)
    {
        // Three genes with lengths; counts scaled by the library-depth factor.
        var genes = new[]
        {
            ("g1", 100.0 * depthScale, 1000),
            ("g2", 300.0 * depthScale, 2000),
            ("g3", 50.0 * depthScale, 500),
        };

        switch (normalization)
        {
            case NormMethod.Tpm:
            {
                var tpm = TranscriptomeAnalyzer.CalculateTPM(genes.Select(g => (g.Item1, g.Item2, g.Item3))).ToList();
                tpm.Sum(t => t.TPM).Should().BeApproximately(1_000_000.0, 1e-3, "TPM values sum to one million");
                // Depth-independence: TPM is identical to the unscaled library's TPM.
                var baseTpm = TranscriptomeAnalyzer.CalculateTPM(genes.Select(g => (g.Item1, g.Item2 / depthScale, g.Item3))).ToList();
                for (int i = 0; i < tpm.Count; i++)
                    tpm[i].TPM.Should().BeApproximately(baseTpm[i].TPM, 1e-6, "TPM does not depend on sequencing depth");
                break;
            }
            case NormMethod.Fpkm:
            {
                double total = genes.Sum(g => g.Item2);
                foreach (var (id, count, len) in genes)
                    TranscriptomeAnalyzer.CalculateFPKM(count, len, total)
                        .Should().BeApproximately(count * 1_000_000_000.0 / (len * total), 1e-6, "FPKM = X·1e9/(l·N)");
                break;
            }
            default:
            {
                // Within-sample-distinct values (no ties) so quantile normalization yields one shared distribution.
                var samples = new IReadOnlyList<double>[]
                {
                    new double[] { 5.0 * depthScale, 2.0 * depthScale, 3.0 * depthScale, 4.0 * depthScale },
                    new double[] { 4.0 * depthScale, 1.0 * depthScale, 7.0 * depthScale, 2.0 * depthScale },
                    new double[] { 3.0 * depthScale, 9.0 * depthScale, 6.0 * depthScale, 8.0 * depthScale },
                };
                var normalized = TranscriptomeAnalyzer.QuantileNormalize(samples).Select(s => s.OrderBy(x => x).ToList()).ToList();
                for (int i = 1; i < normalized.Count; i++)
                    normalized[i].Should().Equal(normalized[0], "quantile normalization makes all samples identically distributed");
                break;
            }
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: TRANS-SPLICE-001 — Alternative splicing (PSI + AS classification) (Transcriptome)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 200.
    // Spec: tests/TestSpecs/TRANS-SPLICE-001.md (canonical CalculatePSI / DetectAlternativeSplicing).
    // ADVANCED §10.
    // Dimensions: nIsoforms(3) × junctionReads(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (rMATS Shen 2014; Wang 2008): Ψ = I/(I+S) ∈ [0,1]; DetectAlternativeSplicing classifies each
    // isoform pair of a gene into one of the five canonical AS classes (SE/RI/A5SS/A3SS/MXE).
    //
    // Axis mapping (documented): nIsoforms → isoforms per gene; junctionReads → inclusion/exclusion
    // read scenario for PSI. The combinatorial point: PSI equals I/(I+S) in [0,1], and a set of n
    // pairwise-distinct isoforms yields exactly C(n,2) classified events, each a valid AS class.
    // ═══════════════════════════════════════════════════════════════════════

    private static TranscriptomeAnalyzer.TranscriptIsoform Iso(string id, params (int Start, int End)[] exons) =>
        new(id, "gene", exons.Sum(e => e.End - e.Start + 1), exons.Length, 1.0, true, exons);

    // Pairwise-distinct isoforms: full / skipped-exon / alt-5'SS / retained-intron.
    private static readonly TranscriptomeAnalyzer.TranscriptIsoform[] Isoforms =
    {
        Iso("full", (1, 100), (200, 300), (400, 500)),
        Iso("skip", (1, 100), (400, 500)),
        Iso("a5ss", (1, 100), (200, 280), (400, 500)),
        Iso("ri",   (1, 100), (200, 500)),
    };

    private static readonly string[] ValidAsClasses =
        { "SkippedExon", "RetainedIntron", "AlternativeFivePrimeSS", "AlternativeThreePrimeSS", "MutuallyExclusiveExons" };

    [Test, Combinatorial]
    public void TransSplice_PsiAndEventCount_AcrossIsoformsAndReads(
        [Values(2, 3, 4)] int nIsoforms,
        [Values(0, 1, 2)] int readScenario)
    {
        // PSI scenario: high-inclusion / balanced / high-exclusion.
        (double I, double S) = readScenario switch
        {
            0 => (90.0, 10.0),
            1 => (50.0, 50.0),
            _ => (10.0, 90.0),
        };
        double psi = TranscriptomeAnalyzer.CalculatePSI(I, S);
        psi.Should().BeApproximately(I / (I + S), 1e-12, "Ψ = I/(I+S)");
        psi.Should().BeInRange(0.0, 1.0);

        var isoforms = Isoforms.Take(nIsoforms).ToList();
        var events = TranscriptomeAnalyzer.DetectAlternativeSplicing(isoforms).ToList();

        events.Should().HaveCount(nIsoforms * (nIsoforms - 1) / 2, "one classified event per distinct isoform pair");
        events.Should().OnlyContain(e => ValidAsClasses.Contains(e.EventType), "every event is a canonical AS class");
        events.Should().OnlyContain(e => e.GeneId == "gene" && e.Start <= e.End);
    }

    /// <summary>
    /// Interaction witnesses — the canonical AS classifications (skipped exon, alt-5'SS, retained
    /// intron) and the length-normalized rMATS PSI form.
    /// </summary>
    [Test]
    public void TransSplice_CanonicalClassifications_AndLengthNormalizedPsi()
    {
        string EventOf(TranscriptomeAnalyzer.TranscriptIsoform a, TranscriptomeAnalyzer.TranscriptIsoform b) =>
            TranscriptomeAnalyzer.DetectAlternativeSplicing(new[] { a, b }).Single().EventType;

        EventOf(Isoforms[0], Isoforms[1]).Should().Be("SkippedExon", "full vs exon-skipped ⇒ SE");
        EventOf(Isoforms[0], Isoforms[2]).Should().Be("AlternativeFivePrimeSS", "same start, different end ⇒ A5SS");
        EventOf(Isoforms[0], Isoforms[3]).Should().Be("RetainedIntron", "an exon spanning the intron ⇒ RI");

        // rMATS length normalization: equal read counts but a longer inclusion form ⇒ Ψ < 0.5.
        TranscriptomeAnalyzer.CalculatePSI(100, 100, inclusionEffectiveLength: 200, exclusionEffectiveLength: 100)
            .Should().BeApproximately((100.0 / 200) / (100.0 / 200 + 100.0 / 100), 1e-12, "length-normalized Ψ");

        TranscriptomeAnalyzer.CalculatePSI(0, 0).Should().Be(double.NaN, "no supporting reads ⇒ Ψ undefined");
    }
}
