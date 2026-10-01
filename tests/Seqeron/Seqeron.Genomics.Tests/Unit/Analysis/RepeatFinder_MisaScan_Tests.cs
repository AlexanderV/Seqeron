using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// REP-STR-001 misa.pl-parity scan (B04 completeness audit WP8, L12):
/// <see cref="MicrosatelliteScanMode.MisaRegex"/> reproduces misa.pl v1.0's SSR list (Thiel et al. 2003; raw GitHub
/// mirror cfljam/SSR_marker_design, executed with perl 5.38). Expected values are copied from real misa.pl runs
/// (an instrumented copy dumps its SSR list); the random test compares against a line-by-line transcription of
/// misa.pl's scan loop on the .NET regex engine (same leftmost-greedy semantics for this pattern).
/// </summary>
[TestFixture]
public class RepeatFinder_MisaScan_Tests
{
    /// <summary>misa.pl reports (AAAGA)9 1-45 then, resuming at 46, (AGAAA)8 46-85 (a .misa row "c (AAAGA)9(AGAAA)8 85 1 85").</summary>
    private const string OverlapCase =
        "AAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGACTTACGTAGATAAGACTTACGTAGATAAG";

    /// <summary>misa.pl's hexamer regex consumes (CTCTCT)5, rejects it, resumes at 131 and reports (TAAACT)6 131-166.</summary>
    private const string ConsumedCase =
        "GCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGGCACAAAAAAGAAAACTATCAGGAATAGAGTATAGAGTAGGGGGGGGGGGGAAACAACTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTAAACTTAAACTTAAACTTAAACTTAAACTTAAACTTAAAATTAAACT";

    internal static string Ssrs(IEnumerable<MicrosatelliteResult> r) =>
        string.Join(";", r.Select(m => $"{m.RepeatUnit}x{m.RepeatCount}@{m.Position + 1}-{m.Position + m.TotalLength}"));

    [Test]
    public void MisaRegex_OverlapAndConsumedCases_EqualMisaPlSsrList()
    {
        var map = RepeatFinder.MisaDefaultMinRepeats;
        Assert.Multiple(() =>
        {
            Assert.That(Ssrs(RepeatFinder.FindMicrosatellites(OverlapCase, map, MicrosatelliteScanMode.MisaRegex)),
                Is.EqualTo("AAAGAx9@1-45;AGAAAx8@46-85"));
            Assert.That(Ssrs(RepeatFinder.FindMicrosatellites(ConsumedCase, map, MicrosatelliteScanMode.MisaRegex)),
                Is.EqualTo("Gx12@83-94;CTx15@101-130;GCGx14@1-42;TAAACTx6@131-166"));
            // The maximal-run convention (default) reports the true run starts instead.
            Assert.That(Ssrs(RepeatFinder.FindMicrosatellites(OverlapCase, map, MicrosatelliteScanMode.MaximalRuns)),
                Is.EqualTo("AAAGAx9@1-45;AAGAAx8@45-84"));
            Assert.That(RepeatFinder.FindMicrosatellites(ConsumedCase, map, MicrosatelliteScanMode.MaximalRuns),
                Is.EqualTo(RepeatFinder.FindMicrosatellites(ConsumedCase, map)));
        });
    }

    [Test]
    public void MisaRegex_Compounds_EqualMisaPlRows()
    {
        var a = RepeatFinder.FindCompoundMicrosatellites(OverlapCase, null, 100, MicrosatelliteScanMode.MisaRegex).Single();
        var b = RepeatFinder.FindCompoundMicrosatellites(ConsumedCase, null, 100, MicrosatelliteScanMode.MisaRegex).Single();
        Assert.Multiple(() =>
        {
            Assert.That((a.MisaType, a.Notation, a.Length, a.Start + 1, a.End), Is.EqualTo(("c", "(AAAGA)9(AGAAA)8", 85, 1, 85)));
            Assert.That((b.MisaType, b.Notation, b.Length, b.Start + 1, b.End), Is.EqualTo(
                ("c", "(GCG)14ggcacaaaaaagaaaactatcaggaatagagtatagagta(G)12aaacaa(CT)15(TAAACT)6", 166, 1, 166)));
            Assert.That(RepeatFinder.FindCompoundMicrosatellites(new DnaSequence(OverlapCase), null, 100, MicrosatelliteScanMode.MisaRegex)
                .Single().Notation, Is.EqualTo(a.Notation));
        });
    }

    /// <summary>
    /// Start ties (misa.ini <c>1-3 2-2 3-2 4-2 5-2 6-2</c>, interruptions 10): (T)3 and (TTTG)5 both start at 1. misa.pl
    /// orders them by Perl's randomised hash order — PERL_HASH_SEED=1 prints <c>(TTTG)5(T)3*g(T)3g(T)3g(T)3g(T)3g(T)3</c>,
    /// PERL_HASH_SEED=2 prints <c>(T)3(TTTG)5*(T)3*g(T)3g(T)3g(T)3g(T)3</c>; here ties keep unit-length (SSR-number)
    /// order, which equals misa.pl with the tie broken by SSR number (and its seed-2 output).
    /// </summary>
    [Test]
    public void MisaRegex_StartTies_UnitLengthOrder()
    {
        var map = new Dictionary<int, int> { [1] = 3, [2] = 2, [3] = 2, [4] = 2, [5] = 2, [6] = 2 };
        const string seq = "TTTGTTTGTTTGTTTGTTTGTTT";
        var c = RepeatFinder.FindCompoundMicrosatellites(seq, map, 10, MicrosatelliteScanMode.MisaRegex).Single();
        Assert.That((c.MisaType, c.Notation, c.Start + 1, c.End), Is.EqualTo(("c*", "(T)3(TTTG)5*(T)3*g(T)3g(T)3g(T)3g(T)3", 1, 23)));
    }

    /// <summary>misa.pl's scan loop transcribed on the .NET regex engine; 400 random SSR-rich sequences × 3 definitions.</summary>
    [Test]
    public void MisaRegex_EqualsTranscriptionOfMisaPlScanLoop()
    {
        var maps = new[]
        {
            new Dictionary<int, int>(RepeatFinder.MisaDefaultMinRepeats),
            new Dictionary<int, int> { [1] = 3, [2] = 2, [3] = 2, [4] = 2, [5] = 2, [6] = 2 },
            new Dictionary<int, int> { [2] = 3, [4] = 2, [6] = 2, [7] = 2 },
        };
        var rng = new Random(12);
        for (int t = 0; t < 400; t++)
        {
            string seq = SsrRich(rng);
            var map = maps[t % maps.Length];
            Assert.That(Ssrs(RepeatFinder.FindMicrosatellites(seq, map, MicrosatelliteScanMode.MisaRegex)),
                Is.EqualTo(Ssrs(MisaPlScan(seq, map))), seq);
        }
    }

    [Test]
    public void MisaRegex_SummaryAndProgress()
    {
        var dna = new DnaSequence(ConsumedCase);
        var summary = RepeatFinder.GetTandemRepeatSummary(dna, RepeatFinder.MisaDefaultMinRepeats, MicrosatelliteScanMode.MisaRegex);
        Assert.That((summary.TotalRepeats, summary.HexanucleotideRepeats, summary.TotalRepeatBases), Is.EqualTo((4, 1, 12 + 30 + 42 + 36)));
        Assert.That(RepeatFinder.GetTandemRepeatSummary(dna, RepeatFinder.MisaDefaultMinRepeats, MicrosatelliteScanMode.MaximalRuns),
            Is.EqualTo(RepeatFinder.GetTandemRepeatSummary(dna, RepeatFinder.MisaDefaultMinRepeats)));

        var values = new List<double>();
        var progress = new SyncProgress(values.Add);
        string longSeq = string.Concat(Enumerable.Repeat("ACGTTGCAAC" + "AT" + "GGCATTACG", 2000));
        _ = RepeatFinder.FindMicrosatellites(longSeq, RepeatFinder.MisaDefaultMinRepeats, MicrosatelliteScanMode.MisaRegex,
            CancellationToken.None, progress).ToList();
        Assert.That(values, Is.Ordered);
        Assert.That(values[^1], Is.EqualTo(1.0));
        Assert.That(values.Take(values.Count - 1), Is.All.LessThan(1.0));
    }

    [Test]
    public void MisaRegex_Validation()
    {
        var map = RepeatFinder.MisaDefaultMinRepeats;
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.FindMicrosatellites(string.Empty, map, MicrosatelliteScanMode.MisaRegex), Is.Empty);
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindMicrosatellites("ACGT", map, (MicrosatelliteScanMode)9));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FindMicrosatellites((DnaSequence)null!, map, MicrosatelliteScanMode.MisaRegex));
            Assert.Throws<ArgumentException>(() => RepeatFinder.FindMicrosatellites("ACGT", new Dictionary<int, int>(), MicrosatelliteScanMode.MisaRegex));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                RepeatFinder.FindCompoundMicrosatellites("ACGT", null, -1, MicrosatelliteScanMode.MisaRegex));
            using var cts = new CancellationTokenSource();
            cts.Cancel();
            Assert.Throws<OperationCanceledException>(() =>
                RepeatFinder.FindMicrosatellites("ACACACACACACAC", map, MicrosatelliteScanMode.MisaRegex, cts.Token).ToList());
        });
    }

    private sealed class SyncProgress(Action<double> report) : IProgress<double>
    {
        public void Report(double value) => report(value);
    }

    /// <summary>
    /// misa.pl lines 101–125: for each unit size (ascending) <c>while ($seq =~ /(([acgt]{p})\2{t-1,})/ig)</c>; skip
    /// redundant motifs (<c>([ACGT]{j})\1{p/j-1}</c> for j &lt; p — a fractional count is a literal brace and never
    /// matches); <c>$end = pos($seq)</c>, <c>$start = $end - length + 1</c>. Positions here 0-based.
    /// </summary>
    internal static List<MicrosatelliteResult> MisaPlScan(string seq, IReadOnlyDictionary<int, int> map)
    {
        var result = new List<MicrosatelliteResult>();
        foreach (int p in map.Keys.OrderBy(x => x))
        {
            var search = new Regex($"(([acgt]{{{p}}})\\2{{{map[p] - 1},}})", RegexOptions.IgnoreCase | RegexOptions.CultureInvariant);
            for (var m = search.Match(seq); m.Success; m = m.NextMatch())
            {
                string motif = m.Groups[2].Value.ToUpperInvariant();
                bool redundant = false;
                for (int j = p - 1; j > 0; j--)
                {
                    if (p % j == 0 && Regex.IsMatch(motif, $"([ACGT]{{{j}}})\\1{{{p / j - 1}}}"))
                        redundant = true;
                }
                if (redundant) continue;
                int copies = m.Groups[1].Length / p;
                result.Add(new MicrosatelliteResult(m.Index, motif, copies, copies * p, RepeatType.Complex));
            }
        }
        return result;
    }

    private static string SsrRich(Random rng)
    {
        var sb = new System.Text.StringBuilder();
        int n = rng.Next(20, 300);
        while (sb.Length < n)
        {
            double r = rng.NextDouble();
            if (r < 0.45)
            {
                int p = rng.Next(1, 8);
                string unit = new(Enumerable.Range(0, p).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
                if (rng.Next(6) == 0) unit = unit[..Math.Max(1, p / 2)] + unit[..Math.Max(1, p / 2)];
                for (int c = rng.Next(2, 12); c > 0; c--) sb.Append(unit);
                if (rng.Next(3) == 0) sb.Append(unit[..rng.Next(unit.Length)]);
            }
            else if (r < 0.9)
            {
                string alphabet = rng.Next(4) switch { 0 => "ACGT", 1 => "AC", 2 => "AT", _ => "A" };
                for (int c = rng.Next(1, 15); c > 0; c--) sb.Append(alphabet[rng.Next(alphabet.Length)]);
            }
            else
            {
                sb.Append('N', rng.Next(1, 4));
            }
        }
        string s = sb.ToString(0, n);
        return rng.Next(5) == 0 ? new string(s.Select(c => rng.Next(3) == 0 ? char.ToLowerInvariant(c) : c).ToArray()) : s;
    }
}
