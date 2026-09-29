namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Tests for RepeatFinder: TandemRepeatSummary (REP-TANDEM-001 D1–D4 + field conventions).
/// Microsatellite detection tests consolidated into RepeatFinder_Microsatellite_Tests.cs (REP-STR-001).
/// Inverted repeat tests consolidated into RepeatFinder_InvertedRepeat_Tests.cs (REP-INV-001).
/// </summary>
/// <remarks>
/// Expected values come from an independent Python reference (brute-force maximal primitive runs of
/// period 1–6, then aggregation: per-class counts as in MISA <c>misa.pl</c> "Distribution to different
/// repeat type classes" / Krait <c>motifTypeStatis</c> Mono…Hexa; total length = Krait <c>SUM(length)</c>;
/// coverage = union of spans / length), cross-checked against the code on 9000 random sequences.
/// </remarks>
[TestFixture]
public class RepeatFinderTests
{
    #region Tandem Repeat Summary Tests

    [Test]
    public void GetTandemRepeatSummary_MixedRepeats_CorrectSummary()
    {
        // Runs: A×6@0, AC×4@15, CGT×3@6, ATG×3@23 — disjoint, together cover all 32 bases.
        var sequence = new DnaSequence("AAAAAACGTCGTCGTACACACACATGATGATG");
        var summary = RepeatFinder.GetTandemRepeatSummary(sequence, 3);

        Assert.Multiple(() =>
        {
            Assert.That(summary.TotalRepeats, Is.EqualTo(4));
            Assert.That(summary.TotalRepeatBases, Is.EqualTo(32));
            Assert.That(summary.PercentageOfSequence, Is.EqualTo(100.0).Within(1e-12));
            Assert.That(summary.MononucleotideRepeats, Is.EqualTo(1));
            Assert.That(summary.DinucleotideRepeats, Is.EqualTo(1));
            Assert.That(summary.TrinucleotideRepeats, Is.EqualTo(2));
            Assert.That(summary.TetranucleotideRepeats, Is.EqualTo(0));
            Assert.That(summary.PentanucleotideRepeats, Is.EqualTo(0));
            Assert.That(summary.HexanucleotideRepeats, Is.EqualTo(0));
            // CGT×3@6 and ATG×3@23 both span 9 bp: tie → shorter unit, then leftmost.
            Assert.That(summary.LongestRepeat, Is.EqualTo(
                new MicrosatelliteResult(6, "CGT", 3, 9, RepeatType.Trinucleotide)));
            Assert.That(summary.MostFrequentUnit, Is.EqualTo("A"));
        });
    }

    [Test]
    public void GetTandemRepeatSummary_HigherMinRepeats_PartialCoverage()
    {
        // minRepeats = 4 keeps only A×6@0 and AC×4@15: 14 of 32 bases = 43.75 %.
        var summary = RepeatFinder.GetTandemRepeatSummary(new DnaSequence("AAAAAACGTCGTCGTACACACACATGATGATG"), 4);

        Assert.Multiple(() =>
        {
            Assert.That(summary.TotalRepeats, Is.EqualTo(2));
            Assert.That(summary.TotalRepeatBases, Is.EqualTo(14));
            Assert.That(summary.PercentageOfSequence, Is.EqualTo(43.75).Within(1e-12));
            Assert.That(summary.LongestRepeat, Is.EqualTo(
                new MicrosatelliteResult(15, "AC", 4, 8, RepeatType.Dinucleotide)));
        });
    }

    [Test]
    public void GetTandemRepeatSummary_NoRepeats_ZeroSummary()
    {
        var sequence = new DnaSequence("ACGT");
        var summary = RepeatFinder.GetTandemRepeatSummary(sequence, 3);

        Assert.Multiple(() =>
        {
            Assert.That(summary.TotalRepeats, Is.EqualTo(0));
            Assert.That(summary.TotalRepeatBases, Is.EqualTo(0));
            Assert.That(summary.PercentageOfSequence, Is.EqualTo(0.0));
            // Nothing found → no longest repeat (must be null, not a default Position-0 record).
            Assert.That(summary.LongestRepeat, Is.Null);
            Assert.That(summary.MostFrequentUnit, Is.Null);
        });
    }

    [Test]
    public void GetTandemRepeatSummary_EmptySequence_ZeroSummaryWithNullLongest()
    {
        var summary = RepeatFinder.GetTandemRepeatSummary(new DnaSequence(""), 3);

        Assert.Multiple(() =>
        {
            Assert.That(summary.TotalRepeats, Is.EqualTo(0));
            Assert.That(summary.PercentageOfSequence, Is.EqualTo(0.0));
            Assert.That(summary.LongestRepeat, Is.Null);
            Assert.That(summary.MostFrequentUnit, Is.Null);
        });
    }

    [Test]
    public void GetTandemRepeatSummary_MononucleotideCount_Correct()
    {
        var sequence = new DnaSequence("AAAAAATTTTTGGGGGCCCCC");
        var summary = RepeatFinder.GetTandemRepeatSummary(sequence, 3);

        Assert.Multiple(() =>
        {
            Assert.That(summary.MononucleotideRepeats, Is.EqualTo(4)); // A×6, T×5, G×5, C×5
            Assert.That(summary.TotalRepeats, Is.EqualTo(4));
            Assert.That(summary.TotalRepeatBases, Is.EqualTo(21));
            // Every unit occurs once: tie → the first unit in the scan order.
            Assert.That(summary.MostFrequentUnit, Is.EqualTo("A"));
            Assert.That(summary.LongestRepeat!.Value.RepeatUnit, Is.EqualTo("A"));
        });
    }

    [Test]
    public void GetTandemRepeatSummary_LongestRepeat_Identified()
    {
        var sequence = new DnaSequence("AAACAGCAGCAGCAGCAGCAGAAA"); // 6x CAG = 18bp
        var summary = RepeatFinder.GetTandemRepeatSummary(sequence, 3);

        Assert.That(summary.LongestRepeat, Is.EqualTo(
            new MicrosatelliteResult(3, "CAG", 6, 18, RepeatType.Trinucleotide)));
    }

    [Test]
    public void GetTandemRepeatSummary_PentaAndHexa_CountedAndSumToTotal()
    {
        // (AAAGA)4 (Penta E forensic motif) + CC + (TTAGGG)4 (vertebrate telomere). Runs: A×3@0, A×4@4, A×4@9,
        // A×4@14, G×3@25, G×3@31, G×3@37, G×3@43, AAAGA×4@0, TTAGGG×4@22 (MISA/Krait report all six classes).
        var summary = RepeatFinder.GetTandemRepeatSummary(
            new DnaSequence("AAAGAAAAGAAAAGAAAAGACCTTAGGGTTAGGGTTAGGGTTAGGG"), 3);

        Assert.Multiple(() =>
        {
            Assert.That(summary.TotalRepeats, Is.EqualTo(10));
            Assert.That(summary.MononucleotideRepeats, Is.EqualTo(8));
            Assert.That(summary.DinucleotideRepeats, Is.EqualTo(0));
            Assert.That(summary.TrinucleotideRepeats, Is.EqualTo(0));
            Assert.That(summary.TetranucleotideRepeats, Is.EqualTo(0));
            Assert.That(summary.PentanucleotideRepeats, Is.EqualTo(1));
            Assert.That(summary.HexanucleotideRepeats, Is.EqualTo(1));
            Assert.That(summary.MononucleotideRepeats + summary.DinucleotideRepeats + summary.TrinucleotideRepeats
                        + summary.TetranucleotideRepeats + summary.PentanucleotideRepeats
                        + summary.HexanucleotideRepeats, Is.EqualTo(summary.TotalRepeats));
            // Sum of lengths counts the homopolymers inside the penta/hexa runs again (Krait SUM(length)) ...
            Assert.That(summary.TotalRepeatBases, Is.EqualTo(71));
            // ... while coverage is the union: [0,20) ∪ [22,46) = 44 of 46 bases.
            Assert.That(summary.PercentageOfSequence, Is.EqualTo(44.0 / 46.0 * 100).Within(1e-12));
            Assert.That(summary.LongestRepeat, Is.EqualTo(
                new MicrosatelliteResult(22, "TTAGGG", 4, 24, RepeatType.Hexanucleotide)));
            // A and G each occur in 4 runs: tie → the unit whose first run comes first ("A").
            Assert.That(summary.MostFrequentUnit, Is.EqualTo("A"));
        });
    }

    [Test]
    public void GetTandemRepeatSummary_OverlappingRuns_SumExceedsCoveredBases()
    {
        // A×5@0 and AT×4@4 share base 4: 13 repeat bases over 12 covered bases (100 %).
        var summary = RepeatFinder.GetTandemRepeatSummary(new DnaSequence("AAAAATATATAT"), 3);

        Assert.Multiple(() =>
        {
            Assert.That(summary.TotalRepeats, Is.EqualTo(2));
            Assert.That(summary.TotalRepeatBases, Is.EqualTo(13));
            Assert.That(summary.PercentageOfSequence, Is.EqualTo(100.0).Within(1e-12));
            Assert.That(summary.LongestRepeat, Is.EqualTo(
                new MicrosatelliteResult(4, "AT", 4, 8, RepeatType.Dinucleotide)));
        });
    }

    [Test]
    public void GetTandemRepeatSummary_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => RepeatFinder.GetTandemRepeatSummary(null!, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.GetTandemRepeatSummary(new DnaSequence("ACGT"), 1));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.GetTandemRepeatSummary(new DnaSequence("ACGT"), 0));
    }

    #endregion
}
