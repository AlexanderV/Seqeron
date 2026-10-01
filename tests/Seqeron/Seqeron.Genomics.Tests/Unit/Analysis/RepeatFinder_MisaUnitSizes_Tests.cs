namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// REP-STR-001 / REP-TANDEM-001, B04 audit WP16 (F61): MISA definitions with unit sizes above 6. misa.pl v1.0
/// (Thiel et al. 2003; raw GitHub cfljam/SSR_marker_design, perl 5.38.2) reads any number of <c>size-min</c> pairs from
/// its <c>def</c> line (lines 76–81) and searches <c>(([acgt]{size})\2{min−1,})</c> for every defined size
/// (lines 101–106); its <c>.statistics</c> "Distribution to different repeat type classes" has one row per unit size
/// with at least one SSR. Expected values below are copied from real misa.pl runs.
/// </summary>
[TestFixture]
public class RepeatFinder_MisaUnitSizes_Tests
{
    /// <summary>
    /// misa.ini <c>1-10 2-6 3-5 4-5 5-5 6-5 7-3 8-3 9-2 10-2</c>, interruptions 10. misa.pl output (seq length 140):
    /// SSRs (ACGTTGC)3 1-21, (A)12 23-34, (TTAGGCA)3 37-57, (ACGT)6 62-85, (ATCCATGCA)2 88-105 (the 8-mer
    /// ACGTACGT is rejected as redundant); one <c>.misa</c> row
    /// <c>c (ACGTTGC)3g(A)12cc(TTAGGCA)3ttcg(ACGT)6gg(ATCCATGCA)2 105 1 105</c>; <c>.statistics</c>: total 5,
    /// distribution 1→1, 4→1, 7→2, 9→1, "SSRs present in compound formation" 4.
    /// </summary>
    private const string Case7To10 =
        "ACGTTGCACGTTGCACGTTGCGAAAAAAAAAAAACCTTAGGCATTAGGCATTAGGCATTCGACGTACGTACGTACGTACGTACGTGGATCCATGCAATCCATGCAATTGCACCTTGAGACCTTGAGANNACCTTGAGAGG";

    private const string Definition7To10 = "1-10 2-6 3-5 4-5 5-5 6-5 7-3 8-3 9-2 10-2";

    [Test]
    public void MisaRegex_UnitSizes7To10_EqualMisaPlSsrsCompoundAndStatistics()
    {
        var map = RepeatFinder.ParseMisaDefinition(Definition7To10);
        var ssrs = RepeatFinder.FindMicrosatellites(Case7To10, map, MicrosatelliteScanMode.MisaRegex).OrderBy(m => m.Position);
        var compound = RepeatFinder.FindCompoundMicrosatellites(Case7To10, map, 10, MicrosatelliteScanMode.MisaRegex).Single();
        var summary = RepeatFinder.GetTandemRepeatSummary(Case7To10, map, MicrosatelliteScanMode.MisaRegex);

        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder_MisaScan_Tests.Ssrs(ssrs),
                Is.EqualTo("ACGTTGCx3@1-21;Ax12@23-34;TTAGGCAx3@37-57;ACGTx6@62-85;ATCCATGCAx2@88-105"));
            Assert.That((compound.MisaType, compound.Notation, compound.Length, compound.Start + 1, compound.End),
                Is.EqualTo(("c", "(ACGTTGC)3g(A)12cc(TTAGGCA)3ttcg(ACGT)6gg(ATCCATGCA)2", 105, 1, 105)));
            Assert.That(summary.TotalRepeats, Is.EqualTo(5));
            Assert.That(summary.CountsByUnitLength, Is.EqualTo(new Dictionary<int, int>
            {
                [1] = 1, [2] = 0, [3] = 0, [4] = 1, [5] = 0, [6] = 0, [7] = 2, [8] = 0, [9] = 1, [10] = 0,
            }));
            Assert.That(summary.CountsByUnitLength.Keys, Is.Ordered);
            Assert.That(summary.CountsByUnitLength.Where(kv => kv.Value > 0).Select(kv => (kv.Key, kv.Value)),
                Is.EqualTo(new[] { (1, 1), (4, 1), (7, 2), (9, 1) })); // misa.pl .statistics rows
            // Named class fields keep their 1-6 meaning; sizes above 6 are only in CountsByUnitLength.
            Assert.That((summary.MononucleotideRepeats, summary.TetranucleotideRepeats, summary.HexanucleotideRepeats),
                Is.EqualTo((1, 1, 0)));
            Assert.That(summary.LongestRepeat?.RepeatUnit, Is.EqualTo("ACGT"));
            Assert.That(RepeatFinder.GetTandemRepeatSummary(new DnaSequence(Case7To10.Replace("N", "")), map, MicrosatelliteScanMode.MisaRegex)
                .CountsByUnitLength[7], Is.EqualTo(2));
        });
    }

    [Test]
    public void ParseMisaDefinition_MisaIniSyntax()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.ParseMisaDefinition("1-10 2-6 3-5 4-5 5-5 6-5"), Is.EqualTo(RepeatFinder.MisaDefaultMinRepeats));
            Assert.That(RepeatFinder.ParseMisaDefinition("definition(unit_size,min_repeats):   7-5 1-10\t10-2"),
                Is.EqualTo(new Dictionary<int, int> { [1] = 10, [7] = 5, [10] = 2 }));
            Assert.That(RepeatFinder.ParseMisaDefinition("2-3, 4-2,6-2").Keys, Is.EqualTo(new[] { 2, 4, 6 }));
            foreach (string bad in new[] { "", "   ", "def", "1-10 2", "1:10", "a-3", "0-5", "3-1", "2-6 2-5", "1-10 -3", "+1-10" })
                Assert.Throws<ArgumentException>(() => RepeatFinder.ParseMisaDefinition(bad), bad);
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.ParseMisaDefinition(null!));
        });
    }

    [Test]
    public void Summary_CountsByUnitLength_UniformOverloadsHaveSizes1To6AndSumToTotal()
    {
        const string seq = "AAAAAAAAACACACACACAGGCAGCAGCAGCAGTTTATTTATTTAGATTAGATTAGATTCCAGTTCCAGTTCCAGT";
        foreach (var s in new[] { RepeatFinder.GetTandemRepeatSummary(seq, 3), RepeatFinder.GetTandemRepeatSummary(new DnaSequence(seq), 3) })
        {
            Assert.That(s.CountsByUnitLength.Keys, Is.EqualTo(new[] { 1, 2, 3, 4, 5, 6 }));
            Assert.That(s.CountsByUnitLength.Values.Sum(), Is.EqualTo(s.TotalRepeats));
            Assert.That(s.CountsByUnitLength.Values, Is.EqualTo(new[]
            {
                s.MononucleotideRepeats, s.DinucleotideRepeats, s.TrinucleotideRepeats,
                s.TetranucleotideRepeats, s.PentanucleotideRepeats, s.HexanucleotideRepeats,
            }));
        }
        Assert.That(default(TandemRepeatSummary).CountsByUnitLength, Is.Empty);
        Assert.That(RepeatFinder.GetTandemRepeatSummary("", RepeatFinder.ParseMisaDefinition("8-2 3-4"), MicrosatelliteScanMode.MisaRegex)
            .CountsByUnitLength, Is.EqualTo(new Dictionary<int, int> { [3] = 0, [8] = 0 }));
    }

    [Test]
    public void Summary_Equality_ComparesCountsByUnitLengthByContent()
    {
        var map = RepeatFinder.ParseMisaDefinition(Definition7To10);
        string acgt = Case7To10.Replace("N", "");
        var a = RepeatFinder.GetTandemRepeatSummary(acgt, map, MicrosatelliteScanMode.MisaRegex);
        var b = RepeatFinder.GetTandemRepeatSummary(new DnaSequence(acgt), map, MicrosatelliteScanMode.MisaRegex);
        var changed = a with { CountsByUnitLength = new Dictionary<int, int>(a.CountsByUnitLength) { [10] = 1 } };
        Assert.Multiple(() =>
        {
            Assert.That(a, Is.EqualTo(b));
            Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));
            Assert.That(changed, Is.Not.EqualTo(a));
            Assert.That(changed.TotalRepeats, Is.EqualTo(a.TotalRepeats));
        });
    }

    /// <summary>
    /// misa.pl's scan loop transcribed on the .NET regex engine (<see cref="RepeatFinder_MisaScan_Tests.MisaPlScan"/>),
    /// 300 random sequences rich in 7–12 bp units × 3 definitions with sizes up to 12: SSR lists equal, and the summary's
    /// per-size counts equal the per-size counts of that list (the .statistics distribution).
    /// </summary>
    [Test]
    public void MisaRegex_LongUnits_EqualTranscriptionAndSummaryCounts()
    {
        var maps = new[]
        {
            RepeatFinder.ParseMisaDefinition("1-10 2-6 3-5 4-5 5-5 6-5 7-5 8-5 9-5 10-5"),
            RepeatFinder.ParseMisaDefinition("1-3 2-2 3-2 4-2 5-2 6-2 7-2 8-2 9-2 10-2 11-2 12-2"),
            RepeatFinder.ParseMisaDefinition("7-2 8-2 9-2 10-2"),
        };
        var rng = new Random(16);
        for (int t = 0; t < 300; t++)
        {
            string seq = LongUnitRich(rng);
            var map = maps[t % maps.Length];
            var expected = RepeatFinder_MisaScan_Tests.MisaPlScan(seq, map);
            Assert.That(RepeatFinder_MisaScan_Tests.Ssrs(RepeatFinder.FindMicrosatellites(seq, map, MicrosatelliteScanMode.MisaRegex)),
                Is.EqualTo(RepeatFinder_MisaScan_Tests.Ssrs(expected)), seq);
            var counts = RepeatFinder.GetTandemRepeatSummary(seq, map, MicrosatelliteScanMode.MisaRegex).CountsByUnitLength;
            Assert.That(counts, Is.EqualTo(map.Keys.ToDictionary(p => p, p => expected.Count(m => m.RepeatUnit.Length == p))), seq);
        }
    }

    private static string LongUnitRich(Random rng)
    {
        var sb = new System.Text.StringBuilder();
        int n = rng.Next(30, 400);
        while (sb.Length < n)
        {
            if (rng.Next(2) == 0)
            {
                int p = rng.Next(1, 13);
                string unit = new(Enumerable.Range(0, p).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
                if (p >= 4 && rng.Next(6) == 0) unit = unit[..(p / 2)] + unit[..(p / 2)]; // non-primitive
                for (int c = rng.Next(2, 7); c > 0; c--) sb.Append(unit);
                if (rng.Next(3) == 0) sb.Append(unit[..rng.Next(unit.Length)]);
            }
            else
            {
                for (int c = rng.Next(1, 12); c > 0; c--) sb.Append("ACGTN"[rng.Next(rng.Next(8) == 0 ? 5 : 4)]);
            }
        }
        return sb.ToString(0, n);
    }
}
