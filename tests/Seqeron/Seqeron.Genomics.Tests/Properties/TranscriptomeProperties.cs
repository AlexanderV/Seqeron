using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for transcriptome analysis (TranscriptomeAnalyzer): differential expression,
/// TPM quantification, and alternative-splicing PSI.
///
/// Test Units: TRANS-DIFF-001, TRANS-EXPR-001
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
}
