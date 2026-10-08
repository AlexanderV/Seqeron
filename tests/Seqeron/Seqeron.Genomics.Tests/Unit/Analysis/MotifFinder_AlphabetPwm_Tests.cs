using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// PAT-PWM-001 — generic-alphabet (protein) position weight matrix: <see cref="MotifFinder.CreateAlphabetPwm(IEnumerable{string}, string, double, IReadOnlyList{double}?, bool)"/>,
/// <see cref="AlphabetPositionWeightMatrix"/>, <see cref="MotifFinder.CalculateAlphabetPwmScores"/>, <see cref="MotifFinder.ScanWithAlphabetPwm"/>.
/// Oracle: Biopython 1.88 <c>motifs.create(instances, alphabet).counts.normalize(pseudocounts).log_odds(background)</c>,
/// <c>consensus</c>, <c>anticonsensus</c>, <c>max</c>, <c>min</c>, <c>mean</c>, <c>std</c>; window scores = Σ_j <c>pssm[letter][j]</c>
/// (NaN for a symbol outside the alphabet, the <c>_pwm.c</c> rule — Biopython's own <c>calculate</c> raises for non-DNA
/// alphabets). For alphabet "ACGT" the scores equal Biopython <c>calculate</c> (float32) and <c>search(both=False)</c>.
/// </summary>
[TestFixture]
public class MotifFinder_AlphabetPwm_Tests
{
    private const string Protein = "ACDEFGHIKLMNPQRSTVWY";
    private static readonly string[] Instances = { "MKVLAT", "MRVLGT", "MKILAS", "LKVLAT", "MKVIAT", "MEVLAT" };
    private const double Tol = 1e-12;

    private static AlphabetPositionWeightMatrix ProteinPwm() => MotifFinder.CreateAlphabetPwm(Instances, Protein, 0.5);

    [Test]
    public void ProteinPwm_ScalarPseudocount_EqualsBiopython()
    {
        var pwm = ProteinPwm();
        Assert.Multiple(() =>
        {
            Assert.That(pwm.Length, Is.EqualTo(6));
            Assert.That(pwm.Alphabet, Is.EqualTo(Protein));
            Assert.That(pwm['M', 0], Is.EqualTo(2.7813597135246595).Within(Tol));
            Assert.That(pwm['k', 1], Is.EqualTo(2.4918530963296748).Within(Tol), "case-insensitive lookup");
            Assert.That(pwm['A', 0], Is.EqualTo(-0.6780719051126377).Within(Tol));
            Assert.That(pwm.Consensus, Is.EqualTo("MKVLAT"));
            Assert.That(pwm.Anticonsensus, Is.EqualTo("AAAACA"), "ties → first symbol in alphabet order");
            Assert.That(pwm.MaxScore, Is.EqualTo(16.398651663952972).Within(Tol));
            Assert.That(pwm.MinScore, Is.EqualTo(-4.068431430675826).Within(Tol));
            Assert.That(pwm.Mean(), Is.EqualTo(3.809139711610947).Within(Tol));
            Assert.That(pwm.Std(), Is.EqualTo(3.8318293687371425).Within(Tol));
        });
    }

    [Test]
    public void CalculateScores_ProteinSequence_NaNForUnknownSymbols_EqualsBiopythonSums()
    {
        // 'x' and 'X' are not in the alphabet → NaN windows; lower-case 'v' is scored as 'V'.
        double[] scores = MotifFinder.CalculateAlphabetPwmScores("GGMKVLATxxMRvLGTPPMKXLAS", ProteinPwm());
        double nan = double.NaN;
        double[] expected =
        {
            -4.068431430675826, -4.068431430675826, 16.398651663952972, nan, nan, nan, nan, nan, nan, nan,
            12.939220045315675, -4.068431430675826, -4.068431430675826, -2.4834689299546704, -4.068431430675826,
            nan, nan, nan, nan,
        };
        Assert.That(scores, Has.Length.EqualTo(expected.Length));
        for (int i = 0; i < expected.Length; i++)
        {
            if (double.IsNaN(expected[i]))
                Assert.That(scores[i], Is.NaN, $"window {i}");
            else
                Assert.That(scores[i], Is.EqualTo(expected[i]).Within(Tol), $"window {i}");
        }
    }

    [Test]
    public void Scan_Protein_ForwardHitsAtThreshold_EqualsBiopythonSearchBothFalseRule()
    {
        var hits = MotifFinder.ScanWithAlphabetPwm("GGMKVLATxxMRvLGTPPMKXLAS", ProteinPwm(), 5.0).ToList();
        Assert.That(hits.Select(h => h.Position), Is.EqualTo(new[] { 2, 10 }));
        Assert.That(hits[0].MatchedSequence, Is.EqualTo("MKVLAT"));
        Assert.That(hits[1].MatchedSequence, Is.EqualTo("MRvLGT"), "window reported as given");
        Assert.That(hits[1].Score, Is.EqualTo(12.939220045315675).Within(Tol));
        Assert.That(hits[0].Pattern, Is.EqualTo("MKVLAT"));
    }

    [Test]
    public void ProteinPwm_PerSymbolPseudocountsAndBackground_EqualsBiopython()
    {
        var pseudo = Protein.Select(a => a is 'M' or 'K' ? 1.0 : 0.1).ToArray();
        var bg = Protein.Select(a => "LAGV".Contains(a) ? 2.0 : 1.0).ToArray();
        var pwm = MotifFinder.CreateAlphabetPwm(Instances, Protein, pseudo, bg);
        Assert.Multiple(() =>
        {
            Assert.That(pwm['L', 3], Is.EqualTo(2.6426779985774442).Within(Tol));
            Assert.That(pwm.MaxScore, Is.EqualTo(19.061964092904915).Within(Tol));
            Assert.That(pwm.MinScore, Is.EqualTo(-18.17848406036431).Within(Tol));
            Assert.That(pwm.Mean(bg), Is.EqualTo(10.253659399954454).Within(Tol));
            Assert.That(pwm.Std(bg), Is.EqualTo(4.8477154676833605).Within(Tol));
            Assert.That(pwm.Consensus, Is.EqualTo("MKVLAT"));
        });
    }

    [Test]
    public void ProteinPwm_ZeroPseudocount_UnseenSymbolsNegativeInfinity_EqualsBiopython()
    {
        var pwm = MotifFinder.CreateAlphabetPwm(Instances, Protein);
        Assert.Multiple(() =>
        {
            Assert.That(pwm.MinScore, Is.EqualTo(double.NegativeInfinity));
            Assert.That(pwm.MaxScore, Is.EqualTo(24.03143403943405).Within(Tol));
            Assert.That(pwm.Mean(), Is.EqualTo(21.429827293694583).Within(Tol), "−∞ cells skipped");
            Assert.That(pwm.Std(), Is.EqualTo(2.1524130640960037).Within(Tol));
            Assert.That(pwm.Anticonsensus, Is.EqualTo("AAAACA"));
        });
        Assert.That(MotifFinder.CalculateAlphabetPwmScores("WWWWWW", pwm)[0], Is.EqualTo(double.NegativeInfinity));
    }

    [Test]
    public void IgnoreUnknownSymbols_SkipsGapsAndX_EqualsBiopythonCounts()
    {
        var pwm = MotifFinder.CreateAlphabetPwm(new[] { "MK-V", "MKXV", "AK-V" }, Protein, 0.2, ignoreUnknownSymbols: true);
        Assert.Multiple(() =>
        {
            Assert.That(pwm['M', 0], Is.EqualTo(2.6520766965796927).Within(Tol));
            Assert.That(pwm['V', 3], Is.EqualTo(3.192645077942396).Within(Tol));
            Assert.That(pwm['A', 2], Is.EqualTo(-3.2034265038149186e-16).Within(Tol), "column with no counted symbol → uniform");
            Assert.That(pwm.Consensus, Is.EqualTo("MKAV"));
        });
        Assert.Throws<ArgumentException>(() => MotifFinder.CreateAlphabetPwm(new[] { "MK-V" }, Protein, 0.2));
    }

    [Test]
    public void DnaAlphabet_ScoresEqualDnaPwm()
    {
        string[] seqs = { "GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT" };
        var generic = MotifFinder.CreateAlphabetPwm(seqs, "ACGT", 0.25);
        var dna = MotifFinder.CreatePwm(seqs, 0.25);
        const string target = "CCTAGGTAAGTNacaggtcagtgg";
        double[] g = MotifFinder.CalculateAlphabetPwmScores(target, generic);
        double[] d = MotifFinder.CalculatePwmScores(target, dna);
        Assert.That(g, Has.Length.EqualTo(d.Length));
        for (int i = 0; i < g.Length; i++)
        {
            if (double.IsNaN(d[i]))
                Assert.That(g[i], Is.NaN);
            else
                Assert.That(g[i], Is.EqualTo(d[i]).Within(1e-13), $"window {i}");
        }

        Assert.That(generic.Consensus, Is.EqualTo(dna.Consensus));
        Assert.That(generic.MaxScore, Is.EqualTo(dna.MaxScore).Within(1e-13));
        Assert.That(generic.Mean(), Is.EqualTo(dna.Mean()).Within(1e-13));
        Assert.That(generic.Std(), Is.EqualTo(dna.Std()).Within(1e-13));
    }

    [Test]
    public void FromCounts_EqualsCreateFromInstances()
    {
        var counts = new double[Protein.Length, 6];
        foreach (string s in Instances)
            for (int j = 0; j < 6; j++)
                counts[Protein.IndexOf(s[j]), j]++;
        var a = AlphabetPositionWeightMatrix.FromCounts(Protein, counts, 0.5);
        var b = ProteinPwm();
        Assert.That(a.GetMatrix(), Is.EqualTo(b.GetMatrix()));
    }

    [Test]
    public void ShortSequence_EmptyScoresAndNoHits()
    {
        Assert.That(MotifFinder.CalculateAlphabetPwmScores("MKV", ProteinPwm()), Is.Empty);
        Assert.That(MotifFinder.ScanWithAlphabetPwm("", ProteinPwm(), double.NegativeInfinity), Is.Empty);
    }

    [Test]
    public void GetMatrix_ReturnsDefensiveCopy()
    {
        var pwm = ProteinPwm();
        var m = pwm.GetMatrix();
        m[0, 0] = 99;
        Assert.That(pwm.GetMatrix()[0, 0], Is.Not.EqualTo(99));
    }

    [Test]
    public void Guards()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => MotifFinder.CreateAlphabetPwm(null!, Protein));
            Assert.Throws<ArgumentNullException>(() => MotifFinder.CreateAlphabetPwm(Instances, null!));
            Assert.Throws<ArgumentException>(() => MotifFinder.CreateAlphabetPwm(Instances, ""));
            Assert.Throws<ArgumentException>(() => MotifFinder.CreateAlphabetPwm(Instances, "ACa"), "duplicate ignoring case");
            Assert.Throws<ArgumentException>(() => MotifFinder.CreateAlphabetPwm(Array.Empty<string>(), Protein));
            Assert.Throws<ArgumentException>(() => MotifFinder.CreateAlphabetPwm(new[] { "MK", "M" }, Protein));
            Assert.Throws<ArgumentException>(() => MotifFinder.CreateAlphabetPwm(new[] { "MZ" }, Protein), "Z not in alphabet");
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.CreateAlphabetPwm(Instances, Protein, -1.0));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.CreateAlphabetPwm(Instances, Protein, double.NaN));
            Assert.Throws<ArgumentException>(() => MotifFinder.CreateAlphabetPwm(Instances, Protein, new double[3]));
            Assert.Throws<ArgumentException>(() => MotifFinder.CreateAlphabetPwm(Instances, Protein, 0.5, new double[4]));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.CreateAlphabetPwm(Instances, Protein, 0.5, new double[20]));
            Assert.Throws<ArgumentException>(() => MotifFinder.CreateAlphabetPwm(new[] { "--" }, Protein, 0.0, ignoreUnknownSymbols: true),
                "zero total column (Biopython ZeroDivisionError)");
            Assert.Throws<ArgumentException>(() => new AlphabetPositionWeightMatrix("AB", new double[3, 2]));
            Assert.Throws<ArgumentException>(() => new AlphabetPositionWeightMatrix("AB", new double[,] { { double.NaN }, { 0 } }));
            Assert.Throws<ArgumentNullException>(() => MotifFinder.CalculateAlphabetPwmScores(null!, ProteinPwm()));
            Assert.Throws<ArgumentNullException>(() => MotifFinder.ScanWithAlphabetPwm("M", null!));
            Assert.Throws<ArgumentException>(() => _ = ProteinPwm()['Z', 0]);
            Assert.Throws<ArgumentOutOfRangeException>(() => _ = ProteinPwm()['M', 6]);
        });
    }
}
