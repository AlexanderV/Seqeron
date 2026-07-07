using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for transcriptome analysis (TranscriptomeAnalyzer): differential expression,
/// TPM quantification, and alternative-splicing PSI.
///
/// Test Units: TRANS-DIFF-001, TRANS-EXPR-001, TRANS-SPLICE-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Transcriptome")]
public class TranscriptomeProperties
{
    /// <summary>Generates a pair of non-negative replicate vectors (2..4 values each).</summary>
    private static Arbitrary<(double[] a, double[] b)> ReplicatePairArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            double[] Vec() => Enumerable.Range(0, 2 + rng.Next(3)).Select(_ => rng.NextDouble() * 100).ToArray();
            return (Vec(), Vec());
        }).ToArbitrary();

    #region TRANS-DIFF-001: S: log2FC(A,B) = −log2FC(B,A); R: p-value ∈ [0,1]; D: deterministic

    // CalculateFoldChange = log2((mean2+c)/(mean1+c)); FindDifferentiallyExpressed adds a Welch t-test
    // p-value with Benjamini-Hochberg adjustment.

    /// <summary>
    /// INV-1 (S): the log2 fold change is antisymmetric — swapping the conditions negates it.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FoldChange_IsAntisymmetric()
    {
        return Prop.ForAll(ReplicatePairArbitrary(), input =>
        {
            var (a, b) = input;
            double ab = TranscriptomeAnalyzer.CalculateFoldChange(a, b);
            double ba = TranscriptomeAnalyzer.CalculateFoldChange(b, a);
            return (Math.Abs(ab + ba) < 1e-9).Label($"log2FC(A,B)={ab} not = −log2FC(B,A)={ba}");
        });
    }

    /// <summary>
    /// INV-2 (I): equal conditions give a zero fold change.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FoldChange_EqualConditions_IsZero()
    {
        return Prop.ForAll(ReplicatePairArbitrary(), input =>
        {
            var (a, _) = input;
            return (Math.Abs(TranscriptomeAnalyzer.CalculateFoldChange(a, a)) < 1e-9)
                .Label("log2FC(A,A) must be 0");
        });
    }

    /// <summary>
    /// INV-3 (M): the fold change increases with the treatment mean — a higher-mean condition 2 gives
    /// a fold change at least as large.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FoldChange_MonotoneInTreatmentMean()
    {
        return Prop.ForAll(ReplicatePairArbitrary(), input =>
        {
            var (a, b) = input;
            var bHigher = b.Select(x => x + 10).ToArray(); // strictly higher mean
            double fc = TranscriptomeAnalyzer.CalculateFoldChange(a, b);
            double fcHigher = TranscriptomeAnalyzer.CalculateFoldChange(a, bHigher);
            return (fcHigher >= fc - 1e-9).Label($"higher treatment mean lowered FC: {fcHigher} < {fc}");
        });
    }

    /// <summary>
    /// INV-4 (R): differential-expression p-values and adjusted p-values lie in [0,1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DifferentialExpression_PValues_InUnitInterval()
    {
        var genesArb = Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            double[] Vec() => Enumerable.Range(0, 3).Select(_ => rng.NextDouble() * 50).ToArray();
            return Enumerable.Range(0, 1 + rng.Next(4))
                .Select(i => ($"g{i}", (IReadOnlyList<double>)Vec(), (IReadOnlyList<double>)Vec()))
                .ToList();
        }).ToArbitrary();

        return Prop.ForAll(genesArb, genes =>
        {
            var de = TranscriptomeAnalyzer.FindDifferentiallyExpressed(genes).ToList();
            return de.All(d => d.PValue is >= 0.0 and <= 1.0 && d.AdjustedPValue is >= 0.0 and <= 1.0)
                .Label("a differential-expression p-value fell outside [0,1]");
        });
    }

    /// <summary>
    /// INV-5 (D): Fold change is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FoldChange_IsDeterministic()
    {
        return Prop.ForAll(ReplicatePairArbitrary(), input =>
        {
            var (a, b) = input;
            return (TranscriptomeAnalyzer.CalculateFoldChange(a, b) == TranscriptomeAnalyzer.CalculateFoldChange(a, b))
                .Label("CalculateFoldChange must be deterministic");
        });
    }

    #endregion

    #region TRANS-EXPR-001: R: TPM ≥ 0; P: Σ TPM = 1e6; D: deterministic

    // CalculateTPM: TPM_i = (X_i/l_i)/Σ(X_j/l_j)·1e6, so TPM within a sample sums to one million
    // (or all zero when every count is zero).

    /// <summary>Generates 1..6 genes with non-negative counts and positive lengths.</summary>
    private static Arbitrary<(string GeneId, double RawCount, int Length)[]> GeneCountsArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int n = 1 + rng.Next(6);
            var genes = new (string, double, int)[n];
            for (int i = 0; i < n; i++)
                genes[i] = ($"g{i}", rng.Next(0, 1000), 1 + rng.Next(2000));
            return genes;
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R + P): all TPM are non-negative, and they sum to 1e6 when any count is positive (else
    /// all zero).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Tpm_NonNegative_AndSumsToMillion()
    {
        return Prop.ForAll(GeneCountsArbitrary(), genes =>
        {
            var expr = TranscriptomeAnalyzer.CalculateTPM(genes.Select(g => (g.GeneId, g.RawCount, g.Length))).ToList();
            bool nonNeg = expr.All(e => e.TPM >= 0);
            double sum = expr.Sum(e => e.TPM);
            bool anyPositive = genes.Any(g => g.RawCount > 0);
            bool sumOk = anyPositive ? Math.Abs(sum - 1_000_000.0) < 1e-3 : sum == 0;
            return (nonNeg && sumOk).Label($"TPM sum={sum}, anyPositive={anyPositive}");
        });
    }

    /// <summary>
    /// INV-2 (boundary): all-zero counts yield all-zero TPM.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Tpm_AllZeroCounts_AreZero()
    {
        var genes = new[] { ("g0", 0.0, 100), ("g1", 0.0, 200) };
        var expr = TranscriptomeAnalyzer.CalculateTPM(genes).ToList();
        Assert.That(expr.All(e => e.TPM == 0), Is.True);
    }

    /// <summary>
    /// INV-3 (D): TPM is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Tpm_IsDeterministic()
    {
        return Prop.ForAll(GeneCountsArbitrary(), genes =>
        {
            var a = TranscriptomeAnalyzer.CalculateTPM(genes.Select(g => (g.GeneId, g.RawCount, g.Length))).Select(e => e.TPM).ToList();
            var b = TranscriptomeAnalyzer.CalculateTPM(genes.Select(g => (g.GeneId, g.RawCount, g.Length))).Select(e => e.TPM).ToList();
            return a.SequenceEqual(b).Label("CalculateTPM must be deterministic");
        });
    }

    #endregion

    #region TRANS-SPLICE-001: R: PSI ∈ [0,1]; M: more inclusion → higher PSI; P: exon coordinates valid; D: deterministic

    // CalculatePSI = inclusion / (inclusion + exclusion) ∈ [0,1] (SUPPA2/rMATS); DetectAlternativeSplicing
    // emits events with valid exon coordinates.

    /// <summary>Two non-negative read counts (0..100).</summary>
    private static Arbitrary<(double inc, double exc)> PsiReadsArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            return (rng.NextDouble() * 100, rng.NextDouble() * 100);
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R): with any supporting reads the PSI lies in [0,1].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Psi_InUnitInterval()
    {
        return Prop.ForAll(PsiReadsArbitrary(), r =>
        {
            double psi = TranscriptomeAnalyzer.CalculatePSI(r.inc + 0.5, r.exc); // +0.5 guarantees support
            return (psi is >= 0.0 and <= 1.0).Label($"PSI={psi} outside [0,1]");
        });
    }

    /// <summary>
    /// INV-2 (M): increasing the inclusion reads (exclusion fixed) does not decrease the PSI.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Psi_MoreInclusion_RaisesPsi()
    {
        return Prop.ForAll(PsiReadsArbitrary(), r =>
        {
            double psi = TranscriptomeAnalyzer.CalculatePSI(r.inc, r.exc + 1);       // ensure support
            double psiMore = TranscriptomeAnalyzer.CalculatePSI(r.inc + 10, r.exc + 1);
            return (psiMore >= psi - 1e-9).Label($"more inclusion lowered PSI: {psiMore} < {psi}");
        });
    }

    /// <summary>
    /// INV-3 (P, boundary): PSI is 0 with no inclusion reads, 1 with no exclusion reads, NaN with no
    /// reads; negative counts are rejected.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Psi_BoundaryCases()
    {
        Assert.Multiple(() =>
        {
            Assert.That(TranscriptomeAnalyzer.CalculatePSI(0, 10), Is.EqualTo(0.0), "no inclusion → 0");
            Assert.That(TranscriptomeAnalyzer.CalculatePSI(10, 0), Is.EqualTo(1.0), "no exclusion → 1");
            Assert.That(double.IsNaN(TranscriptomeAnalyzer.CalculatePSI(0, 0)), Is.True, "no reads → NaN");
            Assert.Throws<ArgumentOutOfRangeException>(() => TranscriptomeAnalyzer.CalculatePSI(-1, 5));
        });
    }

    /// <summary>
    /// INV-4 (P, exon coordinates): a skipped-exon isoform pair yields events with valid (Start ≤ End)
    /// exon coordinates.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Splicing_Events_HaveValidExonCoordinates()
    {
        var full = new TranscriptomeAnalyzer.TranscriptIsoform(
            "tx1", "geneA", 300, 3, 5.0, true,
            new[] { (0, 100), (200, 300), (400, 500) });
        var skipped = new TranscriptomeAnalyzer.TranscriptIsoform(
            "tx2", "geneA", 200, 2, 5.0, true,
            new[] { (0, 100), (400, 500) });

        var events = TranscriptomeAnalyzer.DetectAlternativeSplicing(new[] { full, skipped }).ToList();
        Assert.That(events, Is.Not.Empty, "a skipped-exon event is expected");
        Assert.That(events.All(e => e.Start <= e.End), Is.True, "exon coordinates must satisfy Start ≤ End");
    }

    /// <summary>
    /// INV-5 (D): PSI is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Psi_IsDeterministic()
    {
        return Prop.ForAll(PsiReadsArbitrary(), r =>
            (TranscriptomeAnalyzer.CalculatePSI(r.inc + 1, r.exc + 1) == TranscriptomeAnalyzer.CalculatePSI(r.inc + 1, r.exc + 1))
                .Label("CalculatePSI must be deterministic"));
    }

    #endregion
}
