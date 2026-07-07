// 08_DIFFERENTIAL_TESTING rows 102, 128, 233. Independent oracles: the in-frame modulo formula, a manual
// sequence summary (each metric recomputed), and a naive per-window GC analysis.

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class SummaryFrameGcDifferentialTests
{
    private const double Tol = 1e-9;

    // ---- Row 102: ONCO-FUSION-003 — IsInFrame = (5'codingBases − 3'startPhase) % 3 == 0 ----

    [Test]
    [Category("ONCO-FUSION-003")]
    [TestCase(6, 0, true)]
    [TestCase(7, 1, true)]
    [TestCase(7, 0, false)]
    [TestCase(5, 2, true)]
    [TestCase(4, 0, false)]
    public void IsInFrame_MatchesModuloFormula(int fivePrime, int phase, bool expected)
    {
        Assert.That(OncologyAnalyzer.IsInFrame(fivePrime, phase), Is.EqualTo((fivePrime - phase) % 3 == 0));
        Assert.That(OncologyAnalyzer.IsInFrame(fivePrime, phase), Is.EqualTo(expected));
    }

    // ---- Row 128: SEQ-SUMMARY-001 — SummarizeNucleotideSequence vs per-metric recompute ----

    private static double ShannonOracle(string s)
    {
        var letters = s.ToUpperInvariant().Where(char.IsLetter).ToList();
        if (letters.Count == 0) return 0;
        return -letters.GroupBy(c => c).Select(g => (double)g.Count() / letters.Count).Sum(p => p * Math.Log2(p));
    }
    private static double TmOracle(string s)
    {
        var u = s.ToUpperInvariant();
        int a = u.Count(c => c == 'A'), t = u.Count(c => c == 'T'), g = u.Count(c => c == 'G'), c2 = u.Count(c => c == 'C');
        if (s.Length < 14) return 2.0 * (a + t) + 4.0 * (g + c2);
        int total = a + t + g + c2;
        return total == 0 ? 0 : 64.9 + 41.0 * (g + c2 - 16.4) / total;
    }

    [Test]
    [Category("SEQ-SUMMARY-001")]
    [TestCase("ACGTACGTAC")]
    [TestCase("ATGCATGCATGCATGCAT")]
    public void SummarizeNucleotideSequence_MatchesPerMetricRecompute(string seq)
    {
        var u = seq.ToUpperInvariant();
        int a = u.Count(c => c == 'A'), t = u.Count(c => c == 'T'), g = u.Count(c => c == 'G'), c2 = u.Count(c => c == 'C');
        int valid = a + t + g + c2;

        var sm = SequenceStatistics.SummarizeNucleotideSequence(seq);
        Assert.That(sm.Length, Is.EqualTo(seq.Length));
        Assert.That(sm.GcContent, Is.EqualTo(valid > 0 ? (double)(g + c2) / valid : 0).Within(Tol));
        Assert.That(sm.Entropy, Is.EqualTo(ShannonOracle(seq)).Within(Tol));
        Assert.That(sm.MeltingTemperature, Is.EqualTo(TmOracle(seq)).Within(Tol));
        Assert.That(sm.Complexity, Is.EqualTo(SequenceStatistics.CalculateLinguisticComplexity(seq)).Within(Tol));
        Assert.That(sm.Composition['A'], Is.EqualTo(a));
        Assert.That(sm.Composition['G'], Is.EqualTo(g));
    }

    // ---- Row 233: SEQ-GC-ANALYSIS-001 — windowed GC analysis vs naive per-window count ----

    [Test]
    [Category("SEQ-GC-ANALYSIS-001")]
    public void AnalyzeGcContent_WindowedMatchesNaiveCount()
    {
        const string seq = "GGGCCCAAATTTGGGCCC";
        const int win = 4, step = 2;
        var windowed = GcSkewCalculator.AnalyzeGcContent(seq, win, step).WindowedGcContent;

        var s = seq.ToUpperInvariant();
        int idx = 0;
        for (int i = 0; i + win <= s.Length; i += step, idx++)
        {
            var w = s.Substring(i, win);
            int gc = w.Count(c => c == 'G' || c == 'C');
            int total = w.Count(c => "ATGC".Contains(c));
            double expectedGc = total > 0 ? (double)gc / total * 100 : 0;

            Assert.That(windowed[idx].WindowStart, Is.EqualTo(i), $"start[{idx}]");
            Assert.That(windowed[idx].WindowEnd, Is.EqualTo(i + win - 1), $"end[{idx}]");
            Assert.That(windowed[idx].Position, Is.EqualTo(i + win / 2), $"pos[{idx}]");
            Assert.That(windowed[idx].GcContent, Is.EqualTo(expectedGc).Within(Tol), $"gc[{idx}]");
        }
        Assert.That(windowed.Count, Is.EqualTo(idx));
    }
}
