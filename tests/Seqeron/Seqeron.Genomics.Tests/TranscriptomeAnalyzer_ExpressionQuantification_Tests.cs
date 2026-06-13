// TRANS-EXPR-001 — Expression Quantification (TPM, FPKM, Quantile Normalization)
// Evidence: docs/Evidence/TRANS-EXPR-001-Evidence.md
// TestSpec: tests/TestSpecs/TRANS-EXPR-001.md
// Source: Wagner GP, Kin K, Lynch VJ (2012). Theory in Biosciences 131(4):281-285 (TPM);
//         Zhao S, Ye Z, Stanton R (2020). RNA 26(8) (TPM/RPKM formulas);
//         Pimentel H (2014), "What the FPKM?" (FPKM/TPM review);
//         Bolstad BM et al. (2003). Bioinformatics 19(2):185-193 (quantile normalization,
//         worked example via Wikipedia "Quantile normalization").

using NUnit.Framework;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class TranscriptomeAnalyzer_ExpressionQuantification_Tests
{
    private const double Tol = 1e-9;

    #region CalculateTPM

    // M1 — A(10,2000),B(20,4000),C(30,1000): RPK=(0.005,0.005,0.030), ΣRPK=0.04;
    // TPM = RPK/0.04*1e6 = (125000,125000,750000). Source: Zhao/Ye/Stanton (2020) TPM formula.
    [Test]
    public void CalculateTPM_ThreeGenes_ReturnsExactEvidenceValues()
    {
        var input = new[] { ("A", 10.0, 2000), ("B", 20.0, 4000), ("C", 30.0, 1000) };

        var result = TranscriptomeAnalyzer.CalculateTPM(input).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result[0].TPM, Is.EqualTo(125000.0).Within(Tol),
                "TPM_A = (10/2000)/0.04*1e6 = 125000 per TPM formula (Zhao/Ye/Stanton 2020).");
            Assert.That(result[1].TPM, Is.EqualTo(125000.0).Within(Tol),
                "TPM_B = (20/4000)/0.04*1e6 = 125000 (same RPK 0.005 as A).");
            Assert.That(result[2].TPM, Is.EqualTo(750000.0).Within(Tol),
                "TPM_C = (30/1000)/0.04*1e6 = 750000.");
        });
    }

    // M2 — INV-01: TPM of a non-degenerate sample sums to 10^6. Source: Wagner 2012; Zhao/Ye/Stanton 2020.
    [Test]
    public void CalculateTPM_NonDegenerateSample_SumsToOneMillion()
    {
        var input = new[] { ("A", 10.0, 2000), ("B", 20.0, 4000), ("C", 30.0, 1000) };

        double sum = TranscriptomeAnalyzer.CalculateTPM(input).Sum(g => g.TPM);

        Assert.That(sum, Is.EqualTo(1_000_000.0).Within(1e-6),
            "INV-01: TPM normalizes by Σ(X/l) then scales by 1e6, so a sample sums to exactly 1,000,000 (Wagner 2012).");
    }

    // M3 — INV-02: equal reads-per-kilobase ⇒ equal TPM. Source: Zhao/Ye/Stanton 2020 (TPM monotone in X/l).
    [Test]
    public void CalculateTPM_EqualRatePerKilobase_ProducesEqualTPM()
    {
        // A: 10/2000 = 0.005 ; B: 20/4000 = 0.005 (identical rate).
        var input = new[] { ("A", 10.0, 2000), ("B", 20.0, 4000) };

        var result = TranscriptomeAnalyzer.CalculateTPM(input).ToList();

        Assert.That(result[0].TPM, Is.EqualTo(result[1].TPM).Within(Tol),
            "INV-02: TPM depends only on X/l; equal rates ⇒ equal TPM (here both 500000).");
    }

    // S1 — empty input ⇒ empty sequence. Degenerate case.
    [Test]
    public void CalculateTPM_EmptyInput_ReturnsEmpty()
    {
        var result = TranscriptomeAnalyzer.CalculateTPM(
            Enumerable.Empty<(string, double, int)>());

        Assert.That(result, Is.Empty, "No genes to quantify ⇒ empty output.");
    }

    // C1 — ASSUMPTION-01: all-zero counts ⇒ Σ(X/l)=0 (0/0 undefined) ⇒ all TPM = 0.
    [Test]
    public void CalculateTPM_AllZeroCounts_ReturnsAllZeroTPM()
    {
        var input = new[] { ("A", 0.0, 2000), ("B", 0.0, 4000) };

        var result = TranscriptomeAnalyzer.CalculateTPM(input).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result[0].TPM, Is.EqualTo(0.0),
                "Σ(X/l)=0 makes TPM 0/0 (undefined); convention emits 0 (ASSUMPTION-01).");
            Assert.That(result[1].TPM, Is.EqualTo(0.0),
                "All genes get TPM 0 when every count is zero.");
        });
    }

    #endregion

    #region CalculateFPKM

    // M4 — FPKM = X*1e9/(l*N): 1000*1e9/(2000*1e6) = 500. Source: Zhao/Ye/Stanton 2020; Pimentel 2014.
    [Test]
    public void CalculateFPKM_SingleGene_ReturnsExactEvidenceValue()
    {
        double fpkm = TranscriptomeAnalyzer.CalculateFPKM(rawCount: 1000, length: 2000, totalReads: 1_000_000);

        Assert.That(fpkm, Is.EqualTo(500.0).Within(Tol),
            "FPKM = 1000*1e9/(2000*1e6) = 500 (RPKM formula, Zhao/Ye/Stanton 2020).");
    }

    // M5 — FPKM scales linearly with read count: doubling X doubles FPKM. Source: FPKM formula.
    [Test]
    public void CalculateFPKM_DoubledCount_DoublesFpkm()
    {
        double fpkm = TranscriptomeAnalyzer.CalculateFPKM(rawCount: 2000, length: 2000, totalReads: 1_000_000);

        Assert.That(fpkm, Is.EqualTo(1000.0).Within(Tol),
            "FPKM = 2000*1e9/(2000*1e6) = 1000; FPKM is linear in X (FPKM formula).");
    }

    // S2 — INV-03: non-positive length ⇒ FPKM undefined ⇒ 0.
    [Test]
    public void CalculateFPKM_NonPositiveLength_ReturnsZero()
    {
        double fpkm = TranscriptomeAnalyzer.CalculateFPKM(rawCount: 1000, length: 0, totalReads: 1_000_000);

        Assert.That(fpkm, Is.EqualTo(0.0),
            "INV-03: length ≤ 0 makes the per-kilobase term undefined; convention returns 0.");
    }

    // S3 — INV-03: non-positive total reads ⇒ FPKM undefined ⇒ 0.
    [Test]
    public void CalculateFPKM_NonPositiveTotalReads_ReturnsZero()
    {
        double fpkm = TranscriptomeAnalyzer.CalculateFPKM(rawCount: 1000, length: 2000, totalReads: 0);

        Assert.That(fpkm, Is.EqualTo(0.0),
            "INV-03: total mapped reads N ≤ 0 makes FPKM undefined; convention returns 0.");
    }

    #endregion

    #region QuantileNormalize

    // M6 — Wikipedia worked example (Bolstad 2003). Columns as samples:
    // s1=(5,2,3,4), s2=(4,1,4,2), s3=(3,4,6,8). Rank means r0..r3 = 2, 3, 14/3, 17/3.
    // Verbatim final matrix (rows A-D, read column-wise):
    //   s1=(5.67,2.00,3.00,4.67), s2=(5.17,2.00,5.17,3.00), s3=(2.00,3.00,4.67,5.67).
    [Test]
    public void QuantileNormalize_WikipediaExample_ReproducesEvidenceMatrix()
    {
        var samples = new[]
        {
            new[] { 5.0, 2.0, 3.0, 4.0 }, // column 1
            new[] { 4.0, 1.0, 4.0, 2.0 }, // column 2 (tied 4s at rows A and C)
            new[] { 3.0, 4.0, 6.0, 8.0 }, // column 3
        };
        double r2 = 14.0 / 3.0;   // 4.666…
        double r3 = 17.0 / 3.0;   // 5.666…

        var result = TranscriptomeAnalyzer.QuantileNormalize(samples)
            .Select(s => s.ToArray()).ToList();

        Assert.Multiple(() =>
        {
            // Column 1: 5→r3, 2→r0, 3→r1, 4→r2.
            Assert.That(result[0], Is.EqualTo(new[] { r3, 2.0, 3.0, r2 }).Within(Tol),
                "Column 1 (5,2,3,4) ⇒ (5.67,2.00,3.00,4.67) per Wikipedia/Bolstad worked example.");
            // Column 3: 3→r0, 4→r1, 6→r2, 8→r3.
            Assert.That(result[2], Is.EqualTo(new[] { 2.0, 3.0, r2, r3 }).Within(Tol),
                "Column 3 (3,4,6,8) ⇒ (2.00,3.00,4.67,5.67).");
        });
    }

    // M7 — tie rule (Bolstad 2003): the two tied 4s in column 2 each receive the average of rank
    // means r2,r3 = (14/3+17/3)/2 = 31/6 = 5.166…; verbatim final matrix shows 5.17 for both.
    [Test]
    public void QuantileNormalize_TiedValues_AssignAveragedRankMean()
    {
        var samples = new[]
        {
            new[] { 5.0, 2.0, 3.0, 4.0 },
            new[] { 4.0, 1.0, 4.0, 2.0 }, // tied 4s at positions 0 and 2
            new[] { 3.0, 4.0, 6.0, 8.0 },
        };
        double tieMean = (14.0 / 3.0 + 17.0 / 3.0) / 2.0; // 31/6 = 5.166…

        var col2 = TranscriptomeAnalyzer.QuantileNormalize(samples).ElementAt(1).ToArray();

        Assert.Multiple(() =>
        {
            Assert.That(col2[0], Is.EqualTo(tieMean).Within(Tol),
                "Tied 4 (row A) ⇒ mean of rank means r2,r3 = 31/6 = 5.17 (Bolstad tie rule).");
            Assert.That(col2[2], Is.EqualTo(tieMean).Within(Tol),
                "Tied 4 (row C) ⇒ same 31/6 = 5.17; both tied values share the averaged rank mean.");
            Assert.That(col2[1], Is.EqualTo(2.0).Within(Tol),
                "Value 1 (lowest, rank 0) ⇒ rank mean 2.00.");
            Assert.That(col2[3], Is.EqualTo(3.0).Within(Tol),
                "Value 2 (rank 1) ⇒ rank mean 3.00.");
        });
    }

    // S4 — empty input ⇒ empty sequence. Degenerate case.
    [Test]
    public void QuantileNormalize_EmptyInput_ReturnsEmpty()
    {
        var result = TranscriptomeAnalyzer.QuantileNormalize(
            Enumerable.Empty<IEnumerable<double>>());

        Assert.That(result, Is.Empty, "No samples ⇒ empty output.");
    }

    // S5 — INV-04: within a column, rank order is preserved (output non-decreasing where input increases).
    [Test]
    public void QuantileNormalize_PreservesWithinColumnRankOrder()
    {
        var samples = new[]
        {
            new[] { 10.0, 30.0, 20.0, 40.0 }, // strictly distinct values
            new[] { 5.0, 15.0, 25.0, 35.0 },
        };

        var col0 = TranscriptomeAnalyzer.QuantileNormalize(samples).First().ToArray();

        // Original column 0 ascending order of positions: 0 (10) < 2 (20) < 1 (30) < 3 (40).
        Assert.Multiple(() =>
        {
            Assert.That(col0[0], Is.LessThan(col0[2]),
                "INV-04: position of the smallest value gets the smallest rank mean.");
            Assert.That(col0[2], Is.LessThan(col0[1]),
                "INV-04: rank order is monotone — larger input ⇒ larger rank mean.");
            Assert.That(col0[1], Is.LessThan(col0[3]),
                "INV-04: the largest value gets the largest rank mean.");
        });
    }

    // C2 — identical columns: every rank mean equals the common value, so output equals input.
    [Test]
    public void QuantileNormalize_IdenticalColumns_ReturnsInputUnchanged()
    {
        var samples = new[]
        {
            new[] { 1.0, 2.0, 3.0 },
            new[] { 1.0, 2.0, 3.0 },
        };

        var result = TranscriptomeAnalyzer.QuantileNormalize(samples)
            .Select(s => s.ToArray()).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result[0], Is.EqualTo(new[] { 1.0, 2.0, 3.0 }).Within(Tol),
                "Rank mean of identical columns is the value itself ⇒ unchanged.");
            Assert.That(result[1], Is.EqualTo(new[] { 1.0, 2.0, 3.0 }).Within(Tol),
                "Second identical column is likewise unchanged.");
        });
    }

    #endregion
}
