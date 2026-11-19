namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 / ALIGN-SEMI-001 / ALIGN-MULTI-001 (audit A2-37, B09 F80): one scoring contract for the scalar
/// <see cref="ScoringMatrix"/> and the <see cref="MatrixScoring"/> model.
///
/// Contract: (1) gap scores are finite scores AT MOST 0 (a penalty) and greater than int.MinValue - a positive gap score is rejected by both
/// models with the same message (it was silently accepted by the scalar path: GlobalAlign("ACGT","AC", new ScoringMatrix(1,-1,5,5)) = 30);
/// (2) every dynamic program keeps the score of any path below 2^28 (<c>max|step| * (m + n + 2)</c>), otherwise
/// <see cref="ArgumentException"/> "Score range too large for the integer dynamic program" instead of a silently wrapped score
/// (GlobalAlign("A"x50, "A"x50, new ScoringMatrix(1_000_000_000,-1,-2,-1)) answered 1 999 999 904 for the true 5e10);
/// (3) the multiple-alignment sum-of-pairs scorers accumulate in 64 bits and refuse a total that does not fit the reported int.
/// Oracle: exact integer arithmetic (n identities x match); the accepted scores are checked against the closed forms below.
/// </summary>
[TestFixture]
public class SequenceAligner_ScoringContract_Tests
{
    private static readonly string A50 = new('A', 50);
    private const string RangeMessage = "Score range too large for the integer dynamic program";

    [Test]
    public void ScalarOverflow_IsRejectedByEveryScalarAligner_NotWrapped()
    {
        var huge = new ScoringMatrix(1_000_000_000, -1, -2, -1);

        Assert.That(Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign(A50, A50, huge))!.Message, Does.StartWith(RangeMessage));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign(A50, A50, huge, CancellationToken.None));
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlign(A50, A50, huge));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignAffine(A50, A50, huge));
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlignAffine(A50, A50, huge));
        Assert.Throws<ArgumentException>(() => SequenceAligner.SemiGlobalAlign(A50, A50, huge));
        Assert.Throws<ArgumentException>(() => SequenceAligner.SemiGlobalAlignAffine(A50, A50, huge));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignLinearSpace(A50, A50, huge));
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlignLinearSpace(A50, A50, huge));
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlignTopN(A50, A50, huge, 3, 1));
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlignAffineTopN(A50, A50, huge, 3, 1));
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlignCoOptimal(A50, A50, huge, 3));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignCoOptimal(A50, A50, huge, 3));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignAffineCoOptimal(A50, A50, huge, 3));
    }

    [Test]
    public void ScalarLargeScores_BelowTheBound_AreExact()
    {
        // 50 identities x 1e6 = 5e7 (< 2^28 = 268 435 456), max|step| * (m + n + 2) = 1e6 * 102 = 1.02e8 < 2^28
        var big = new ScoringMatrix(1_000_000, -1, -2, -1);
        Assert.That(SequenceAligner.GlobalAlign(A50, A50, big).Score, Is.EqualTo(50_000_000));
        Assert.That(SequenceAligner.LocalAlign(A50, A50, big).Score, Is.EqualTo(50_000_000));
        Assert.That(SequenceAligner.GlobalAlignAffine(A50, A50, big).Score, Is.EqualTo(50_000_000));
        Assert.That(SequenceAligner.LocalAlignAffine(A50, A50, big).Score, Is.EqualTo(50_000_000));
        Assert.That(SequenceAligner.SemiGlobalAlign(A50, A50, big).Score, Is.EqualTo(50_000_000));
    }

    [Test]
    public void ScalarRangeBound_IsTheSameAsTheMatrixPathBound()
    {
        // max|step| * (m + n + 2) < 2^28: with step 1e6 the longest accepted pair has m + n + 2 <= 268 -> m = n = 133
        var s = new ScoringMatrix(1_000_000, -1, -2, -1);
        string ok = new('A', 133), tooLong = new('A', 134);
        Assert.DoesNotThrow(() => SequenceAligner.GlobalAlign(ok, ok, s));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign(tooLong, tooLong, s));
    }

    [Test]
    public void IntMinValueGaps_AreRejected_WithArgumentException_NotIndexOutOfRange()
    {
        var min = new ScoringMatrix(1, -1, int.MinValue, int.MinValue);
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignAffine("ACGT", "ACGT", min));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign("ACGT", "ACGT", new ScoringMatrix(1, -1, 0, int.MinValue)));
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlignAffine("ACGT", "ACGT", min));
    }

    [Test]
    public void PositiveGapScore_IsRejected_ByScalarAndMatrixModelAlike()
    {
        Assert.That(Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign("ACGT", "AC", new ScoringMatrix(1, -1, 5, 5)))!.Message,
            Does.Contain("must be a finite score <= 0"));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign("ACGT", "AC", new ScoringMatrix(1, -1, 0, 1)));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignAffine("ACGT", "AC", new ScoringMatrix(1, -1, 1, 0)));
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlign("ACGT", "AC", new ScoringMatrix(1, -1, -1, 2)));
        var matrix = Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign("ACGT", "AC", new MatrixScoring(SubstitutionMatrix.EdnaFull, 1, -1)));
        Assert.That(matrix!.Message, Does.Contain("must be a finite score <= 0"));
        // zero gaps stay valid for both
        Assert.DoesNotThrow(() => SequenceAligner.GlobalAlign("ACGT", "AC", new ScoringMatrix(1, -1, 0, 0)));
        Assert.DoesNotThrow(() => SequenceAligner.GlobalAlign("ACGT", "AC", new MatrixScoring(SubstitutionMatrix.EdnaFull, 0, 0)));
    }

    [Test]
    public void IntMinValueGaps_ThroughMultipleAlign_SurfaceAsArgumentException_ForEveryK()
    {
        var min = new ScoringMatrix(1, -1, int.MinValue, int.MinValue);
        foreach (int k in new[] { 2, 3, 5 })
        {
            var seqs = Enumerable.Range(0, k).Select(i => new DnaSequence("ACGTACGTAC".Substring(0, 8 + i % 2))).ToArray();
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlign(seqs, min, new MsaOptions { GapModel = MsaGapModel.Affine }), $"k={k}");
        }
    }

    [Test]
    public void SumOfPairs_Accumulates64Bit_AndRefusesATotalBeyondInt()
    {
        // 3 rows x 2 columns of identical residues, match 2e9: 3 pairs x 2 columns x 2e9 = 1.2e10 (the old int accumulator wrapped)
        var rows = new[] { "AA", "AA", "AA" };
        var ex = Assert.Throws<ArgumentException>(() => SequenceAligner.SumOfPairsScore(rows, new ScoringMatrix(2_000_000_000, 0, 0, 0)));
        Assert.That(ex!.Message, Does.Contain("12000000000"));
        Assert.Throws<ArgumentException>(() => SequenceAligner.SumOfPairsScore(rows, new ScoringMatrix(2_000_000_000, 0, 0, 0), MsaGapModel.Affine));
        // the matrix scorer returns the exact real number (a 64-bit accumulator): 3 pairs x 2 columns x BLOSUM62 A/A = 4
        Assert.That(SequenceAligner.MatrixSumOfPairsScore(rows, new MatrixScoring(SubstitutionMatrix.Blosum62, -1, -1)), Is.EqualTo(24.0));
        // in range: unchanged
        Assert.That(SequenceAligner.SumOfPairsScore(rows, new ScoringMatrix(1_000_000, 0, 0, 0)), Is.EqualTo(6_000_000));
    }
}
