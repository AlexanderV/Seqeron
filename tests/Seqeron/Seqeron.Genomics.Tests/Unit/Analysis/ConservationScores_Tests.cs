using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// B09 WP67 (F308, R14-13): the conservation scores of Capra and Singh (2007) (<see cref="ConservationScores"/>, <see cref="MultipleAlignmentStatistics.ConservationProfile"/>, the members of <see cref="MsaColumnStatistics"/>).
/// Expected values come from an independent implementation of the formulas of the paper in Python (written from the definitions, not from the program: frequency with the pseudo count 1e-7, h = sum f ln f / ln(min(21, rows)),
/// JSD = (sum f log2(f / m) + sum b log2(b / m)) / 2 with m = (f + b) / 2 and the BLOSUM62 background of the paper, the gap penalty 1 - gap fraction) and hand derivations stated in the tests.
/// </summary>
[TestFixture]
public class ConservationScores_Tests
{
    [TestCase("AAAAA", 0.9999953455735638, 0.7990563231997274)]
    [TestCase("AAAAC", 0.6890781790246829, 0.7560746838749699)]
    [TestCase("ACDEF", -3.403541291246981e-06, 0.5511114399832872)]
    [TestCase("ARNDC", -3.4035412910249363e-06, 0.5637373791192319)]
    [TestCase("LLLLLLLLLL", 0.9999983131414484, 0.7733283727046794)]
    [TestCase("AALL", 0.49999389715757536, 0.6515500311614348)]
    public void ShannonAndJsd_EqualTheFormulasOfThePaper(string column, double shannon, double jsd)
    {
        Assert.Multiple(() =>
        {
            Assert.That(ConservationScores.ShannonEntropy(column), Is.EqualTo(shannon).Within(1e-13));
            Assert.That(ConservationScores.JsDivergence(column, ConservationScores.BlosumBackground), Is.EqualTo(jsd).Within(1e-13));
        });
    }

    [Test]
    public void GapPenalty_MultipliesByTheFractionOfNonGapRows()
    {
        // AAA--: 3 of 5 rows are residues: score x 0.6 (independent values: shannon 0.5818300761109975 / jsd 0.7990537184114778 without the penalty)
        Assert.Multiple(() =>
        {
            Assert.That(ConservationScores.WeightedGapPenalty("AAA--"), Is.EqualTo(0.6).Within(1e-15));
            Assert.That(ConservationScores.ShannonEntropy("AAA--", gapPenalty: false), Is.EqualTo(0.5818300761109975).Within(1e-13));
            Assert.That(ConservationScores.ShannonEntropy("AAA--"), Is.EqualTo(0.34909804566659847).Within(1e-13));
            Assert.That(ConservationScores.JsDivergence("AAA--", ConservationScores.BlosumBackground, gapPenalty: false), Is.EqualTo(0.7990537184114778).Within(1e-13));
            Assert.That(ConservationScores.JsDivergence("AAA--", ConservationScores.BlosumBackground), Is.EqualTo(0.47943223104688665).Within(1e-13));
        });
    }

    [Test]
    public void SumOfPairs_HandDerivedWithBlosum62()
    {
        // BLOSUM62: A/A 4, A/C 0, C/C 9, A/R -1. Column AAC: pairs (A,A) 4, (A,C) 0, (A,C) 0: sum 4 / 3 pairs = 1.3333333333333333
        // Column ACR: A/C 0, A/R -1, C/R -3: -4 / 3 = -1.3333333333333333; a gap row leaves the pairs of the others and multiplies by 2/3
        Assert.Multiple(() =>
        {
            Assert.That(ConservationScores.SumOfPairs("AAC"), Is.EqualTo(4.0 / 3).Within(1e-15));
            Assert.That(ConservationScores.SumOfPairs("ACR"), Is.EqualTo(-4.0 / 3).Within(1e-15));
            Assert.That(ConservationScores.SumOfPairs("AAC-"), Is.EqualTo(4.0 / 3 * 0.75).Within(1e-15));
            Assert.That(ConservationScores.SumOfPairs("A-"), Is.EqualTo(0.0));
            Assert.That(ConservationScores.SumOfPairs("AA", weights: [0.5, 2.0]), Is.EqualTo(4.0).Within(1e-15), "weights cancel in the mean of a conserved column");
            Assert.Throws<KeyNotFoundException>(() => ConservationScores.SumOfPairs("AJ"));
        });
    }

    [Test]
    public void WindowScore_AveragesTheNeighbours()
    {
        // lambda 0.5: position 1 of [1, 2, 3, 4, 5] with window 1 = 0.5 x mean(1, 3) + 0.5 x 2 = 2; the ends and the negative scores stay
        var w = ConservationScores.WindowScore([1, 2, 3, 4, 5], 1);
        Assert.That(w, Is.EqualTo(new[] { 1.0, 2.0, 3.0, 4.0, 5.0 }));
        var w2 = ConservationScores.WindowScore([0, 4, 0, -1, 8], 1);
        // position 1: 0.5 x mean(0, 0) + 0.5 x 4 = 2; position 2: neighbours 4 and (-1 ignored) = 0.5 x 4 + 0.5 x 0 = 2; position 3 is negative and stays
        Assert.That(w2, Is.EqualTo(new[] { 0.0, 2.0, 2.0, -1.0, 8.0 }));
    }

    [Test]
    public void SequenceWeights_HenikoffExample()
    {
        // the classical example of Henikoff and Henikoff (1994): columns of AAAA / ABBB-like data reduced to 4 rows x 4 columns
        // rows: ACDE, ACDF, ACGF, AHGF. column 1: A x4 (1 type): every row 1/(4 x 1) = .25 ... derived in the test body
        double[] w = ConservationScores.SequenceWeights(["ACDE", "ACDF", "ACGF", "AHGF"]);
        double[] expected = new double[4];
        string[] rows = ["ACDE", "ACDF", "ACGF", "AHGF"];
        for (int c = 0; c < 4; c++)
        {
            var groups = rows.GroupBy(r => r[c]).ToDictionary(g => g.Key, g => g.Count());
            for (int r = 0; r < 4; r++)
                expected[r] += 1.0 / (groups[rows[r][c]] * groups.Count);
        }

        for (int r = 0; r < 4; r++)
            expected[r] /= 4;
        Assert.That(w, Is.EqualTo(expected).Within(1e-15));
        Assert.That(w.Sum(), Is.EqualTo(1.0).Within(1e-15), "the weights of a gap-free alignment sum to 1");
    }

    [Test]
    public void ConservationProfile_MapsLikeTheScript()
    {
        // rows are upper-cased, X is a gap, B is D: column 3 'ABX' -> 'ADX' -> 'AD-'
        var profile = MultipleAlignmentStatistics.ConservationProfile(["AAcA", "AAdA", "AAbA"], ConservationMethod.JsDivergence);
        Assert.That(profile, Has.Length.EqualTo(4));
        Assert.That(profile[0], Is.EqualTo(ConservationScores.JsDivergence("AAA", ConservationScores.BlosumBackground)).Within(1e-15));
        Assert.That(profile[2], Is.EqualTo(ConservationScores.JsDivergence("CDD", ConservationScores.BlosumBackground)).Within(1e-15));
        var windowed = MultipleAlignmentStatistics.ConservationProfile(["AAcA", "AAdA", "AAbA"], ConservationMethod.ShannonEntropy, windowSize: 1);
        var plain = MultipleAlignmentStatistics.ConservationProfile(["AAcA", "AAdA", "AAbA"], ConservationMethod.ShannonEntropy);
        Assert.That(windowed[1], Is.EqualTo(0.5 * ((plain[0] + plain[2]) / 2) + 0.5 * plain[1]).Within(1e-15));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentStatistics.ConservationProfile(["AC", "A"]));
    }

    [Test]
    public void MsaColumnStatistics_ExtendedMembers()
    {
        var stats = MultipleAlignmentStatistics.Calculate(["AAC-", "AAG-", "AAC-", "ACC-"]);
        var col0 = stats.Columns[0];
        var col1 = stats.Columns[1];
        var col3 = stats.Columns[3];
        Assert.Multiple(() =>
        {
            Assert.That(col0.ShannonEntropy, Is.EqualTo(0.0));
            // column 1: A x3, C x1: -(0.75 log2 0.75 + 0.25 log2 0.25) = 0.8112781244591328
            Assert.That(col1.ShannonEntropy, Is.EqualTo(0.8112781244591328).Within(1e-15));
            Assert.That(col3.ShannonEntropy, Is.EqualTo(0.0), "an all-gap column");
            Assert.That(col0.JensenShannonDivergence(), Is.EqualTo(ConservationScores.JsDivergence("AAAA", ConservationScores.BlosumBackground)).Within(1e-15));
            Assert.That(col1.SumOfPairs(), Is.EqualTo(ConservationScores.SumOfPairs("AAAC")).Within(1e-15));
            // information content with a uniform nucleotide background: column 0 conserved = log2(1 / 0.25) = 2 bits
            Assert.That(col0.InformationContent(), Is.EqualTo(2.0).Within(1e-15));
            Assert.That(col1.RowCount, Is.EqualTo(4));
        });
    }
}
