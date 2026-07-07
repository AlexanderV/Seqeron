// 08_DIFFERENTIAL_TESTING rows 182, 223, 226. Independent oracles: manual bisulfite C->T conversion,
// a manual six-frame translation (own NCBI table-1 map + revcomp), and manual alignment statistics.

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class EpiTransAlignDifferentialTests
{
    private const double Tol = 1e-12;

    // ---- Row 182: EPIGEN-BISULF-001 — bisulfite C->T (methylated C protected) ----

    private static string BisulfiteOracle(string seq, ISet<int> methylated)
    {
        var sb = new System.Text.StringBuilder();
        for (int i = 0; i < seq.Length; i++)
        {
            char c = seq[i];
            if (c == 'C' || c == 'c')
                sb.Append(methylated.Contains(i) ? c : (c == 'C' ? 'T' : 't'));
            else sb.Append(c);
        }
        return sb.ToString();
    }

    [Test]
    [Category("EPIGEN-BISULF-001")]
    public void SimulateBisulfiteConversion_MatchesManual()
    {
        const string seq = "ACGCGTccaa";
        foreach (var meth in new[] { new HashSet<int>(), new HashSet<int> { 1 }, new HashSet<int> { 3, 6 } })
        {
            Assert.That(EpigeneticsAnalyzer.SimulateBisulfiteConversion(seq, meth),
                Is.EqualTo(BisulfiteOracle(seq, meth)), string.Join(",", meth));
        }
    }

    // ---- Row 223: TRANS-SIXFRAME-001 — six-frame translation vs manual ----

    private static readonly Dictionary<string, char> Dna = BuildDnaCode();
    private static Dictionary<string, char> BuildDnaCode()
    {
        const string bases = "TCAG";
        const string aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
        var d = new Dictionary<string, char>();
        int i = 0;
        foreach (char b1 in bases) foreach (char b2 in bases) foreach (char b3 in bases) d[$"{b1}{b2}{b3}"] = aa[i++];
        return d;
    }
    private static readonly Dictionary<char, char> Comp = new() { ['A'] = 'T', ['T'] = 'A', ['G'] = 'C', ['C'] = 'G' };
    private static string RevComp(string s) { var a = s.Select(c => Comp[c]).ToArray(); Array.Reverse(a); return new string(a); }
    private static string Tr(string s, int off)
    {
        var sb = new System.Text.StringBuilder();
        for (int i = off; i + 3 <= s.Length; i += 3) sb.Append(Dna[s.Substring(i, 3)]);
        return sb.ToString();
    }

    [Test]
    [Category("TRANS-SIXFRAME-001")]
    [TestCase("ATGAAATTTGGGCCC")]
    [TestCase("ACGTACGTACGT")]
    public void TranslateSixFrames_MatchesManual(string seq)
    {
        var rc = RevComp(seq);
        var expected = new Dictionary<int, string>
        {
            [1] = Tr(seq, 0), [2] = Tr(seq, 1), [3] = Tr(seq, 2),
            [-1] = Tr(rc, 0), [-2] = Tr(rc, 1), [-3] = Tr(rc, 2),
        };
        var actual = Translator.TranslateSixFrames(new DnaSequence(seq));
        foreach (var kv in expected)
            Assert.That(actual[kv.Key].Sequence, Is.EqualTo(kv.Value), $"frame {kv.Key}");
    }

    // ---- Row 226: ALIGN-STATS-001 — alignment statistics vs manual column count ----

    [Test]
    [Category("ALIGN-STATS-001")]
    [TestCase("ACGTACGT", "ACGAACGT")]
    [TestCase("ACGT", "ACGTACGT")]
    [TestCase("GATTACA", "GATTACA")]
    public void CalculateStatistics_MatchesManualColumnCount(string s1, string s2)
    {
        var r = SequenceAligner.GlobalAlign(new DnaSequence(s1), new DnaSequence(s2));
        string a1 = r.AlignedSequence1, a2 = r.AlignedSequence2;
        int m = 0, mm = 0, g = 0;
        for (int i = 0; i < a1.Length; i++)
        {
            if (a1[i] == '-' || a2[i] == '-') g++;
            else if (a1[i] == a2[i]) m++;
            else mm++;
        }
        int len = a1.Length;

        var st = SequenceAligner.CalculateStatistics(r);
        Assert.That(st.Matches, Is.EqualTo(m));
        Assert.That(st.Mismatches, Is.EqualTo(mm));
        Assert.That(st.Gaps, Is.EqualTo(g));
        Assert.That(st.AlignmentLength, Is.EqualTo(len));
        Assert.That(st.Identity, Is.EqualTo(100.0 * m / len).Within(Tol));
        Assert.That(st.GapPercent, Is.EqualTo(100.0 * g / len).Within(Tol));
        Assert.That(st.Similarity, Is.EqualTo(100.0 * m / len).Within(Tol)); // SimpleDna mismatch < 0 -> similar == matches
    }
}
