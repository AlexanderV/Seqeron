using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>tandem_repeat_summary</c> MCP tool.
/// Expected values derived from the microsatellite aggregation on a single (CAG)3 STR
/// and an STR-free sequence, NOT the wrapper output.
/// </summary>
[TestFixture]
public class TandemRepeatSummaryTests
{
    [Test]
    public void TandemRepeatSummary_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.TandemRepeatSummary("CAGCAGCAG", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.TandemRepeatSummary("", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.TandemRepeatSummary(null!, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.TandemRepeatSummary("XYZ", 3));
    }

    [Test]
    public void TandemRepeatSummary_Binding_InvokesSuccessfully()
    {
        // (CAG)3 -> one trinucleotide STR spanning the whole 9-mer.
        var s = AnalysisTools.TandemRepeatSummary("CAGCAGCAG", 3);
        Assert.Multiple(() =>
        {
            Assert.That(s.TotalRepeats, Is.EqualTo(1));
            Assert.That(s.TotalRepeatBases, Is.EqualTo(9));
            Assert.That(s.PercentageOfSequence, Is.EqualTo(100.0).Within(1e-9));
            Assert.That(s.TrinucleotideRepeats, Is.EqualTo(1));
            Assert.That(s.MononucleotideRepeats, Is.EqualTo(0));
            Assert.That(s.DinucleotideRepeats, Is.EqualTo(0));
            Assert.That(s.MostFrequentUnit, Is.EqualTo("CAG"));
            Assert.That(s.LongestRepeat, Is.EqualTo(new MicrosatelliteItem(0, "CAG", 3, 9, "Trinucleotide")));
            Assert.That(s.PentanucleotideRepeats, Is.EqualTo(0));
            Assert.That(s.HexanucleotideRepeats, Is.EqualTo(0));
        });

        // No STRs.
        var none = AnalysisTools.TandemRepeatSummary("ACGT", 3);
        Assert.Multiple(() =>
        {
            Assert.That(none.TotalRepeats, Is.EqualTo(0));
            Assert.That(none.TotalRepeatBases, Is.EqualTo(0));
            Assert.That(none.PercentageOfSequence, Is.EqualTo(0.0).Within(1e-9));
            // No STR → no longest repeat (previously a default Position-0 item with a null unit).
            Assert.That(none.LongestRepeat, Is.Null);
            Assert.That(none.MostFrequentUnit, Is.Null);
        });
    }

    [Test]
    public void TandemRepeatSummary_PentaHexa_DocExample3()
    {
        // docs/mcp/tools/analysis/tandem_repeat_summary.md Example 3; values from the independent Python
        // brute-force maximal-run reference: (AAAGA)4@0, (TTAGGG)4@22 and 8 homopolymer runs; 71 repeat bases,
        // union coverage [0,20) U [22,46) = 44/46.
        var s = AnalysisTools.TandemRepeatSummary("AAAGAAAAGAAAAGAAAAGACCTTAGGGTTAGGGTTAGGGTTAGGG", 3);
        Assert.Multiple(() =>
        {
            Assert.That(s.TotalRepeats, Is.EqualTo(10));
            Assert.That(s.TotalRepeatBases, Is.EqualTo(71));
            Assert.That(s.PercentageOfSequence, Is.EqualTo(95.65217391304348).Within(1e-12));
            Assert.That(s.MononucleotideRepeats, Is.EqualTo(8));
            Assert.That(s.PentanucleotideRepeats, Is.EqualTo(1));
            Assert.That(s.HexanucleotideRepeats, Is.EqualTo(1));
            Assert.That(s.LongestRepeat, Is.EqualTo(new MicrosatelliteItem(22, "TTAGGG", 4, 24, "Hexanucleotide")));
            Assert.That(s.MostFrequentUnit, Is.EqualTo("A"));
            // perl misa.pl (definition 1-3 … 6-3) .statistics classified table: A/T 4, C/G 4, AAAAG/CTTTT 1, AACCCT/AGGGTT 1.
            Assert.That(s.CanonicalMotifCounts, Is.EquivalentTo(new Dictionary<string, int>
            {
                ["A/T"] = 4, ["C/G"] = 4, ["AAAAG/CTTTT"] = 1, ["AACCCT/AGGGTT"] = 1,
            }));
        });
    }

    [Test]
    public void TandemRepeatSummary_MisaThresholds_CanonicalMotifCounts_MatchMisaPl()
    {
        // perl misa.pl (default misa.ini) on (AC)6 / (GT)6 / (TG)6 / (A)13 separated by a 101-bp SSR-free spacer:
        // 4 SSRs (p2 (AC)6, p2 (GT)6, p2 (TG)6, p1 (A)14); .statistics "Frequency of classified repeat types": A/T 1, AC/GT 3.
        const string spacer =
            "CTAAGCCAACTGCATTGCTAGAGCGAAGTCTTCGTAATGGACCGACCGTTCTGTCCGGACTAGTGAATCGCTGTACAAGTCCGAGGCATCAAGGACTAGTA";
        string seq = "ACACACACACACT" + spacer + "GTGTGTGTGTGTA" + spacer + "TGTGTGTGTGTGA" + spacer + "AAAAAAAAAAAAAG";
        var s = AnalysisTools.TandemRepeatSummary(seq, misaThresholds: true);
        Assert.Multiple(() =>
        {
            Assert.That(s.TotalRepeats, Is.EqualTo(4));
            Assert.That(s.DinucleotideRepeats, Is.EqualTo(3));
            Assert.That(s.MononucleotideRepeats, Is.EqualTo(1));
            Assert.That(s.CanonicalMotifCounts, Is.EquivalentTo(new Dictionary<string, int> { ["A/T"] = 1, ["AC/GT"] = 3 }));
        });

        // Default thresholds (minRepeats 3) also fill the class table: (CAG)3 → AGC/CTG.
        Assert.That(AnalysisTools.TandemRepeatSummary("CAGCAGCAG", 3).CanonicalMotifCounts,
            Is.EquivalentTo(new Dictionary<string, int> { ["AGC/CTG"] = 1 }));
    }
    [Test]
    public void TandemRepeatSummary_StandardMotifLevel_MatchesKraitCounts()
    {
        // misa.pl .statistics classified table A/T 2, AC/GT 4, ACAT/ATGT 1; Krait standard motifs (level 2) A 2, AC 4, ATAC 1
        // (RepeatFinder_MisaCompound_Tests.GetCanonicalMotifFrequencies_StatSequence_MatchesMisaClassifiedTable).
        const string spacer =
            "CTAAGCCAACTGCATTGCTAGAGCGAAGTCTTCGTAATGGACCGACCGTTCTGTCCGGACTAGTGAATCGCTGTACAAGTCCGAGGCATCAAGGACTAGTA";
        string seq = "ACACACACACACT" + spacer + "CACACACACACAT" + spacer + "GTGTGTGTGTGTA" + spacer + "TGTGTGTGTGTGA" + spacer +
            "ACATACATACATACATACATG" + spacer + "AAAAAAAAAAAAG" + spacer + "TTTTTTTTTTTTTG";
        var s = AnalysisTools.TandemRepeatSummary(seq, misaThresholds: true, standardMotifLevel: 2);
        Assert.Multiple(() =>
        {
            Assert.That(s.CanonicalMotifCounts, Is.EquivalentTo(new Dictionary<string, int> { ["A/T"] = 2, ["AC/GT"] = 4, ["ACAT/ATGT"] = 1 }));
            Assert.That(s.StandardMotifCounts, Is.EquivalentTo(new Dictionary<string, int> { ["A"] = 2, ["AC"] = 4, ["ATAC"] = 1 }));
            Assert.That(AnalysisTools.TandemRepeatSummary(seq, misaThresholds: true).StandardMotifCounts, Is.Null);
            Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.TandemRepeatSummary(seq, standardMotifLevel: 5));
        });
    }

    [Test]
    public void TandemRepeatSummary_AcceptsNAndIupac_LikeFindMicrosatellites()
    {
        // WP11: tandem_repeat_summary used to reject N (DnaSequence) while find_microsatellites accepted it.
        // misa.pl .statistics on this 87-mer (1-10 2-6 3-5 4-5 5-5 6-5): size 87, 3 SSRs, unit sizes 1:1 2:2.
        const string seq = "ACACACACANCACACACACACACACAGTNNNNNNNNNNNNNNNNNNNNATATATATATATATATATRYSWKMggggggggggggggg";
        var misa = AnalysisTools.TandemRepeatSummary(seq, misaThresholds: true, misaScan: true);
        var max = AnalysisTools.TandemRepeatSummary(seq, 3);
        Assert.Multiple(() =>
        {
            Assert.That(misa.TotalRepeats, Is.EqualTo(3));
            Assert.That(misa.MononucleotideRepeats, Is.EqualTo(1));
            Assert.That(misa.DinucleotideRepeats, Is.EqualTo(2));
            Assert.That(max.TotalRepeats, Is.EqualTo(AnalysisTools.FindMicrosatellites(seq, 1, 6, 3).Items.Length));
            Assert.That(AnalysisTools.TandemRepeatSummary("AAAAAAAAAANNNNNNNNNN", misaThresholds: true).PercentageOfSequence,
                Is.EqualTo(50.0));
            Assert.Throws<ArgumentException>(() => AnalysisTools.TandemRepeatSummary("ACGU", 3));
            Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.TandemRepeatSummary("ACGT", 1));
        });
    }

    [Test]
    public void TandemRepeatSummary_MisaScan_UniformThresholds_MatchesMisaPl()
    {
        // misa.pl 1-3 2-3 3-3 4-3 5-3 6-3 on the 87-mer: 4 SSRs, unit sizes 1:1 2:3 ((AC)4an(CA)8 ... (AT)9 ... (G)15).
        const string seq = "ACACACACANCACACACACACACACAGTNNNNNNNNNNNNNNNNNNNNATATATATATATATATATRYSWKMggggggggggggggg";
        var s = AnalysisTools.TandemRepeatSummary(seq, 3, misaScan: true);
        Assert.Multiple(() =>
        {
            Assert.That(s.TotalRepeats, Is.EqualTo(4));
            Assert.That(s.MononucleotideRepeats, Is.EqualTo(1));
            Assert.That(s.DinucleotideRepeats, Is.EqualTo(3));
            Assert.That(s.CanonicalMotifCounts, Is.EquivalentTo(new Dictionary<string, int> { ["C/G"] = 1, ["AC/GT"] = 2, ["AT/AT"] = 1 }));
        });
    }

    /// <summary>
    /// misaDefinition (B04 F61): custom misa.ini definition with unit sizes 7-10. misa.pl v1.0 with
    /// <c>1-10 2-6 3-5 4-5 5-5 6-5 7-3 8-3 9-2 10-2</c> / interruptions 10 reports (ACGTTGC)3 1-21, (A)12 23-34,
    /// (TTAGGCA)3 37-57, (ACGT)6 62-85, (ATCCATGCA)2 88-105 and the row
    /// <c>c (ACGTTGC)3g(A)12cc(TTAGGCA)3ttcg(ACGT)6gg(ATCCATGCA)2 105 1 105</c>; .statistics distribution 1→1, 4→1, 7→2, 9→1.
    /// </summary>
    [Test]
    public void TandemRepeatSummary_MisaDefinition_CountsByUnitLength_MatchMisaPlStatistics()
    {
        const string seq = "ACGTTGCACGTTGCACGTTGCGAAAAAAAAAAAACCTTAGGCATTAGGCATTAGGCATTCGACGTACGTACGTACGTACGTACGTGGATCCATGCAATCCATGCAATTGCACCTTGAGACCTTGAGANNACCTTGAGAGG";
        var s = AnalysisTools.TandemRepeatSummary(seq, misaScan: true, misaDefinition: "1-10 2-6 3-5 4-5 5-5 6-5 7-3 8-3 9-2 10-2");
        var classic = AnalysisTools.TandemRepeatSummary(seq, misaThresholds: true, misaScan: true);
        Assert.Multiple(() =>
        {
            Assert.That(s.TotalRepeats, Is.EqualTo(5));
            Assert.That(s.CountsByUnitLength!.Where(kv => kv.Value > 0).OrderBy(kv => kv.Key).Select(kv => (kv.Key, kv.Value)),
                Is.EqualTo(new[] { (1, 1), (4, 1), (7, 2), (9, 1) }));
            Assert.That(s.CountsByUnitLength!.Keys.Order(), Is.EqualTo(Enumerable.Range(1, 10)));
            Assert.That(classic.CountsByUnitLength!.Keys.Order(), Is.EqualTo(Enumerable.Range(1, 6)));
            Assert.That(classic.CountsByUnitLength!.Values.Sum(), Is.EqualTo(classic.TotalRepeats));
            Assert.Throws<ArgumentException>(() => AnalysisTools.TandemRepeatSummary(seq, misaDefinition: ""));
            Assert.Throws<ArgumentException>(() => AnalysisTools.TandemRepeatSummary(seq, misaThresholds: true, misaDefinition: "7-3"));
        });
    }
}
