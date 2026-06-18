using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for transcriptome analysis (TranscriptomeAnalyzer): differential expression,
/// TPM quantification, and alternative-splicing PSI.
///
/// Test Units: TRANS-DIFF-001
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
}
