// 08_DIFFERENTIAL_TESTING rows 153, 169, 171, 174. Independent oracles: a Watson-Crick + wobble pairing
// table, a majority-vote consensus, a threshold-IUPAC consensus, and a brute-force minimum-Hamming best
// match.

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class MatchingRnaDifferentialTests
{
    // ---- Row 153: RNA-PAIR-001 — CanPair / GetBasePairType vs WC+wobble table ----

    [Test]
    [Category("RNA-PAIR-001")]
    public void BasePairing_MatchesWatsonCrickWobbleTable()
    {
        var wc = new HashSet<string> { "AU", "UA", "GC", "CG" };
        var wobble = new HashSet<string> { "GU", "UG" };

        foreach (char b1 in "ACGU")
        foreach (char b2 in "ACGU")
        {
            string p = $"{b1}{b2}";
            RnaSecondaryStructure.BasePairType? expected = wc.Contains(p) ? RnaSecondaryStructure.BasePairType.WatsonCrick
                                   : wobble.Contains(p) ? RnaSecondaryStructure.BasePairType.Wobble
                                   : (RnaSecondaryStructure.BasePairType?)null;
            Assert.That(RnaSecondaryStructure.CanPair(b1, b2), Is.EqualTo(expected != null), p);
            Assert.That(RnaSecondaryStructure.GetBasePairType(b1, b2), Is.EqualTo(expected), p);
        }
    }

    // ---- Row 169: MOTIF-CONS-001 — majority consensus (A<C<G<T tie-break) ----

    private static readonly char[] Alpha = { 'A', 'C', 'G', 'T' };

    private static string MajorityConsensus(string[] seqs)
    {
        int len = seqs[0].Length;
        var sb = new System.Text.StringBuilder();
        for (int col = 0; col < len; col++)
        {
            var counts = new int[4];
            foreach (var s in seqs) counts[Array.IndexOf(Alpha, s[col])]++;
            int best = 0;
            for (int b = 1; b < 4; b++) if (counts[b] > counts[best]) best = b;
            sb.Append(Alpha[best]);
        }
        return sb.ToString();
    }

    [Test]
    [Category("MOTIF-CONS-001")]
    public void CreateConsensusFromAlignment_MatchesMajorityVote()
    {
        var seqs = new[] { "ACGT", "ACGA", "ATGT" };
        Assert.That(MotifFinder.CreateConsensusFromAlignment(seqs), Is.EqualTo(MajorityConsensus(seqs)));
        // Tie -> alphabetically earliest base.
        var tie = new[] { "AC", "GC" };
        Assert.That(MotifFinder.CreateConsensusFromAlignment(tie), Is.EqualTo("AC"));
    }

    // ---- Row 171: MOTIF-GENERATE-001 — threshold-IUPAC consensus ----

    private static readonly Dictionary<string, char> Iupac = new()
    {
        ["A"] = 'A', ["C"] = 'C', ["G"] = 'G', ["T"] = 'T', ["AG"] = 'R', ["CT"] = 'Y', ["CG"] = 'S',
        ["AT"] = 'W', ["GT"] = 'K', ["AC"] = 'M', ["CGT"] = 'B', ["AGT"] = 'D', ["ACT"] = 'H', ["ACG"] = 'V',
    };

    private static string IupacConsensus(string[] seqs)
    {
        int len = seqs[0].Length, total = seqs.Length;
        double threshold = total * 0.25;
        var sb = new System.Text.StringBuilder();
        for (int col = 0; col < len; col++)
        {
            var counts = new Dictionary<char, int> { ['A'] = 0, ['C'] = 0, ['G'] = 0, ['T'] = 0 };
            foreach (var s in seqs) if (counts.TryGetValue(s[col], out int value)) counts[s[col]] = ++value;
            var present = counts.Where(kv => kv.Value > threshold).Select(kv => kv.Key).OrderBy(c => c).ToList();
            if (present.Count == 0) sb.Append(counts.MaxBy(kv => kv.Value).Key);
            else sb.Append(Iupac.GetValueOrDefault(string.Concat(present), 'N'));
        }
        return sb.ToString();
    }

    [Test]
    [Category("MOTIF-GENERATE-001")]
    public void GenerateConsensus_MatchesThresholdIupac()
    {
        foreach (var seqs in new[]
        {
            new[] { "AG", "AG", "AG", "AG" },
            new[] { "AG", "AC", "AG", "AC" },   // col1 {C,G} -> S
            new[] { "ACGT", "AGGT", "ATGT", "AAGT" },
        })
        {
            Assert.That(MotifFinder.GenerateConsensus(seqs), Is.EqualTo(IupacConsensus(seqs)),
                string.Join(",", seqs));
        }
    }

    // ---- Row 174: PAT-APPROX-003 — best match = first minimum-Hamming window ----

    [Test]
    [Category("PAT-APPROX-003")]
    [TestCase("ACGTACGT", "AGGT")]
    [TestCase("AAAAA", "AA")]
    [TestCase("ACGTACGT", "ACGT")]   // exact
    public void FindBestMatch_MatchesBruteMinHamming(string seq, string pat)
    {
        int bestDist = int.MaxValue, bestPos = -1;
        for (int i = 0; i + pat.Length <= seq.Length; i++)
        {
            int d = 0;
            for (int j = 0; j < pat.Length; j++) if (seq[i + j] != pat[j]) d++;
            if (d < bestDist) { bestDist = d; bestPos = i; }
        }
        var best = ApproximateMatcher.FindBestMatch(seq, pat);
        Assert.That(best.HasValue, Is.True);
        Assert.That(best!.Value.Position, Is.EqualTo(bestPos));
        Assert.That(best.Value.Distance, Is.EqualTo(bestDist));
    }
}
