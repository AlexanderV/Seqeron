// MOTIF-DISCOVER-001 (audit group E6) — RSAT dyad-analysis (spaced dyads).
// Evidence: docs/Evidence/MOTIF-DISCOVER-001-Evidence.md; TestSpec: tests/TestSpecs/MOTIF-DISCOVER-001.md
// Source: van Helden J, Rios AF, Collado-Vides J (2000) Nucleic Acids Res 28:1808-1818; rsa-tools/rsat-code master
//         10043f2 perl-scripts/dyad-analysis v1.78 (CountDyads, CalcPossibleDyads, SecondElement, CalcOccSum, SumRCDyads,
//         CalcExpFreqFromMonads, CalcDyadFrequencies, CalcExpectedOcc, CalcProba), lib/RSA.seq.lib OverlapCoeff.
// Reference values: RSAT dyad-analysis run from that clone (perl 5.38) with a dump of every pattern field at %.17g
//         before PrintResult. RSAT rounds z-scores to 2 decimals (%7.2f); the unrounded value is locked and checked
//         against RSAT's to ±0.005. Where RSAT code loses or double-counts dyads (see the tests) the documented
//         behaviour is locked and the RSAT number is quoted.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class MotifFinder_DyadAnalysis_Tests
{
    private const string T2 = "ATGCATGCATGCAAATTTGGGCCCATGCTTAGCGGATCCATGCATGCTTTAAACGTACGTAGC";

    private const double Rel = 1e-10; // percent units (1e-12 relative)

    private static DyadStatistics Get(DyadAnalysisReport r, string pattern) => r.Dyads.Single(d => d.Pattern == pattern);

    private static DyadAnalysisReport Run(DyadAnalysisOptions options, params string[] sequences)
        => MotifFinder.AnalyzeDyads(sequences.Select(s => new DnaSequence(s)), options);

    // RSAT: dyad-analysis -i t2.fa -l 3 -sp 0-2 -1str -return occ,proba,zscore,ratio,freq -lth occ 2 (default -noov)
    //   nb_possible_dyads 12288, 13 tested, monad positions 61; per spacing (valid, occ_sum, ovl_sum):
    //   0: (58, 54, 4), 1: (57, 54, 3), 2: (56, 54, 2).
    //   atgn{0}cat occ 2 overlaps 1 exp_freq 0.008062348830959418 exp_occ 0.43536683687180855 obs_freq 0.034482758620689655
    //     ratio 4.5938271604938281 var 0.42225837663093196 ov_coef 1.00390625 zscore 2.41 occ_P 0.079911679827561213
    //     occ_E 1.0388518377582958 occ_sig -0.016553612396124719;
    //   catn{1}cat (dir_rep) exp_freq 0.0067186240257995156 occ_P 0.056496561382709906 ov_coef 1.015625 zscore 2.72;
    //   atgn{2}tgc occ_P 0.10242325819904796 var 0.47364750180216564.
    [Test]
    public void SingleStrand_NoOverlap_EqualsRsat()
    {
        var r = Run(new DyadAnalysisOptions { MonadLength = 3, MinSpacing = 0, MaxSpacing = 2, Strands = OligoStrandMode.Single, MinCount = 2 }, T2);
        var d = Get(r, "ATGn{0}CAT");
        Assert.Multiple(() =>
        {
            Assert.That(r.PossibleDyads, Is.EqualTo(12288));
            Assert.That(r.TestedPatterns, Is.EqualTo(13));
            Assert.That(r.MonadOccurrences, Is.EqualTo(61));
            Assert.That(r.Spacings.Select(s => (s.PossiblePositions, s.Occurrences, s.Overlaps)),
                Is.EqualTo(new[] { (58L, 54L, 4L), (57L, 54L, 3L), (56L, 54L, 2L) }));
            Assert.That(d.Occurrences, Is.EqualTo(2));
            Assert.That(d.Overlaps, Is.EqualTo(1));
            Assert.That(d.FirstMonad, Is.EqualTo("ATG"));
            Assert.That(d.SecondMonad, Is.EqualTo("CAT"));
            Assert.That(d.Spacing, Is.Zero);
            Assert.That(d.ExpectedFrequency, Is.EqualTo(0.008062348830959418).Within(Rel).Percent);
            Assert.That(d.ExpectedOccurrences, Is.EqualTo(0.43536683687180855).Within(Rel).Percent);
            Assert.That(d.ObservedFrequency, Is.EqualTo(0.034482758620689655).Within(Rel).Percent);
            Assert.That(d.Ratio, Is.EqualTo(4.5938271604938281).Within(Rel).Percent);
            Assert.That(d.ExpectedVariance, Is.EqualTo(0.42225837663093196).Within(Rel).Percent);
            Assert.That(d.OverlapCoefficient, Is.EqualTo(1.00390625));
            Assert.That(d.ZScore, Is.EqualTo(2.41).Within(0.005));
            Assert.That(d.OccurrenceProbability, Is.EqualTo(0.079911679827561213).Within(Rel).Percent);
            Assert.That(d.OccurrenceEValue, Is.EqualTo(1.0388518377582958).Within(Rel).Percent);
            Assert.That(d.OccurrenceSignificance, Is.EqualTo(-0.016553612396124719).Within(1e-11));
            Assert.That(Get(r, "CATn{1}CAT").IsDirectRepeat, Is.True);
            Assert.That(Get(r, "CATn{1}CAT").ExpectedFrequency, Is.EqualTo(0.0067186240257995156).Within(Rel).Percent);
            Assert.That(Get(r, "CATn{1}CAT").OccurrenceProbability, Is.EqualTo(0.056496561382709906).Within(Rel).Percent);
            Assert.That(Get(r, "CATn{1}CAT").OverlapCoefficient, Is.EqualTo(1.015625));
            Assert.That(Get(r, "CATn{1}CAT").ZScore, Is.EqualTo(2.72).Within(0.005));
            Assert.That(Get(r, "ATGn{2}TGC").OccurrenceProbability, Is.EqualTo(0.10242325819904796).Within(Rel).Percent);
            Assert.That(Get(r, "ATGn{2}TGC").ExpectedVariance, Is.EqualTo(0.47364750180216564).Within(Rel).Percent);
        });
    }

    // RSAT: dyad-analysis -i t2.fa -l 2 -sp 0-6 -type dr -1str -ovlp -return occ,proba,zscore -lth occ 2
    //   nb_possible_dyads 16·7 = 112, 6 tested; atn{2}at occ 4 exp_freq 0.016649323621227886 occ_P 0.015985324783370104
    //   occ_E 0.095911948700220628 zscore 3.83; tgn{2}tg occ_P 0.038043032254936944.
    [Test]
    public void DirectRepeats_Overlapping_EqualsRsat()
    {
        var r = Run(new DyadAnalysisOptions { MonadLength = 2, MinSpacing = 0, MaxSpacing = 6, Type = DyadType.DirectRepeat,
            Strands = OligoStrandMode.Single, CountOverlapping = true, MinCount = 2 }, T2);
        var at = Get(r, "ATn{2}AT");
        Assert.Multiple(() =>
        {
            Assert.That(r.PossibleDyads, Is.EqualTo(112));
            Assert.That(r.TestedPatterns, Is.EqualTo(6));
            Assert.That(r.Dyads.All(x => x.IsDirectRepeat), Is.True);
            Assert.That(at.Occurrences, Is.EqualTo(4));
            Assert.That(at.ExpectedFrequency, Is.EqualTo(0.016649323621227886).Within(Rel).Percent);
            Assert.That(at.OccurrenceProbability, Is.EqualTo(0.015985324783370104).Within(Rel).Percent);
            Assert.That(at.OccurrenceEValue, Is.EqualTo(0.095911948700220628).Within(Rel).Percent);
            Assert.That(at.ZScore, Is.EqualTo(3.83).Within(0.005));
            Assert.That(Get(r, "TGn{2}TG").OccurrenceProbability, Is.EqualTo(0.038043032254936944).Within(Rel).Percent);
        });
    }

    // RSAT defaults (-2str -noov): dyad-analysis -i t2.fa -l 3 -sp 0-2 -return occ,proba,zscore,ratio,freq -lth occ 2
    //   nb_possible_dyads 6240; per spacing (valid, occ_sum, ovl_sum) 0: (58, 50, 8), 1: (57, 48, 9), 2: (56, 53, 3);
    //   atgn{0}cat (reverse palindrome) exp_freq 0.008062348830959418 exp_occ 0.40311744154797091 ratio 4.961333333333334
    //   occ_P 0.079911679827561213 zscore 2.46; gcan{0}tgc occ_P 0.054219681283078143.
    //   RSAT tests 9 dyads: SumRCDyads sums a pair only from its lexicographically smaller member, so the observed pairs
    //   aagn{0}cat|atgn{0}ctt and aagn{1}atg|catn{1}ctt (seen only as the larger member, twice each) are lost; the
    //   documented pair sum tests 11, so occ_E = occ_P · 11.
    [Test]
    public void BothStrands_RsatDefaults_EqualRsat_WithLostPairsRestored()
    {
        var r = Run(new DyadAnalysisOptions { MonadLength = 3, MinSpacing = 0, MaxSpacing = 2, MinCount = 2 }, T2);
        var d = Get(r, "ATGn{0}CAT");
        Assert.Multiple(() =>
        {
            Assert.That(r.PossibleDyads, Is.EqualTo(6240));
            Assert.That(r.Spacings.Select(s => (s.PossiblePositions, s.Occurrences, s.Overlaps)),
                Is.EqualTo(new[] { (58L, 50L, 8L), (57L, 48L, 9L), (56L, 53L, 3L) }));
            Assert.That(r.TestedPatterns, Is.EqualTo(11));
            Assert.That(d.IsReversePalindrome, Is.True);
            Assert.That(d.ReverseComplement, Is.EqualTo("ATGn{0}CAT"));
            Assert.That(d.ExpectedFrequency, Is.EqualTo(0.008062348830959418).Within(Rel).Percent);
            Assert.That(d.ExpectedOccurrences, Is.EqualTo(0.40311744154797091).Within(Rel).Percent);
            Assert.That(d.Ratio, Is.EqualTo(4.961333333333334).Within(Rel).Percent);
            Assert.That(d.OccurrenceProbability, Is.EqualTo(0.079911679827561213).Within(Rel).Percent);
            Assert.That(d.OccurrenceEValue, Is.EqualTo(0.079911679827561213 * 11).Within(Rel).Percent);
            Assert.That(d.ZScore, Is.EqualTo(2.46).Within(0.005));
            Assert.That(Get(r, "GCAn{0}TGC").OccurrenceProbability, Is.EqualTo(0.054219681283078143).Within(1e-12));
            var lost = Get(r, "AAGn{0}CAT");
            Assert.That(lost.ReverseComplement, Is.EqualTo("ATGn{0}CTT"));
            Assert.That(lost.Occurrences, Is.EqualTo(2));
            Assert.That(r.Dyads.Any(x => x.Pattern == "AAGn{1}ATG"), Is.True);
        });
    }

    [Test]
    public void BothStrands_PairCountsEqualSingleStrandSums_WhenOverlapping()
    {
        var one = Run(new DyadAnalysisOptions { MonadLength = 2, MaxSpacing = 3, Strands = OligoStrandMode.Single, CountOverlapping = true }, T2);
        var two = Run(new DyadAnalysisOptions { MonadLength = 2, MaxSpacing = 3, CountOverlapping = true }, T2);
        var single = one.Dyads.ToDictionary(d => d.Pattern, d => d.Occurrences);
        Assert.Multiple(() =>
        {
            foreach (var d in two.Dyads)
            {
                int expected = single.GetValueOrDefault(d.Pattern)
                               + (d.ReverseComplement != d.Pattern ? single.GetValueOrDefault(d.ReverseComplement!) : 0);
                Assert.That(d.Occurrences, Is.EqualTo(expected), d.Pattern);
            }
            Assert.That(two.Dyads.Sum(d => d.Occurrences), Is.EqualTo(one.Dyads.Sum(d => d.Occurrences)));
        });
    }

    // Brute force on t2 (independent of RSAT's SecondElement): with -type rep a dyad whose monad is its own reverse
    // complement (e.g. ATn{s}AT) is both a direct and an inverted repeat and is counted once. RSAT 1.78 lists such a
    // monad twice in SecondElement and adds its occurrences twice to occ_sum and to the tested count (40/40 random runs).
    [Test]
    public void RepeatType_PalindromicMonadCountedOnce()
    {
        var r = Run(new DyadAnalysisOptions { MonadLength = 2, MinSpacing = 2, MaxSpacing = 2, Type = DyadType.Repeat,
            Strands = OligoStrandMode.Single, CountOverlapping = true }, T2);
        int brute = 0;
        for (int p = 0; p + 6 <= T2.Length; p++)
        {
            string a = T2.Substring(p, 2), b = T2.Substring(p + 4, 2);
            if (a == b || b == DnaSequence.GetReverseComplementString(a)) brute++;
        }
        Assert.Multiple(() =>
        {
            Assert.That(r.Spacings[0].Occurrences, Is.EqualTo(brute));
            Assert.That(r.Dyads.Sum(d => d.Occurrences), Is.EqualTo(brute));
            Assert.That(r.Dyads.Select(d => d.Pattern).Distinct().Count(), Is.EqualTo(r.Dyads.Count));
        });
    }

    // RSAT without -lth occ: every (observed monad, spacing, observed monad) combination is tested.
    // dyad-analysis -i t2.fa -l 3 -sp 0-2 -1str: 3267 tested = 33 distinct monads² × 3 spacings.
    [Test]
    public void MinCountZero_TestsEveryMonadCombination_AsRsat()
    {
        var r = Run(new DyadAnalysisOptions { MonadLength = 3, MinSpacing = 0, MaxSpacing = 2, Strands = OligoStrandMode.Single, MinCount = 0 }, T2);
        Assert.That(r.TestedPatterns, Is.EqualTo(3267));
    }

    [Test]
    public void DyadFrequencyTable_UsedAndMissingFallsBackOrIsNotTested()
    {
        var table = new Dictionary<string, double> { ["atgn{0}cat"] = 0.01 };
        var r = Run(new DyadAnalysisOptions { MonadLength = 3, MaxSpacing = 0, Strands = OligoStrandMode.Single, MinCount = 2,
            Background = DyadBackgroundModel.DyadFrequencies(table) }, T2);
        var big = Run(new DyadAnalysisOptions { MonadLength = 3, MaxSpacing = 21, Strands = OligoStrandMode.Single, MinCount = 2,
            Background = DyadBackgroundModel.DyadFrequencies(table) }, T2);
        var monads = Run(new DyadAnalysisOptions { MonadLength = 3, MaxSpacing = 21, Strands = OligoStrandMode.Single, MinCount = 2 }, T2);
        Assert.Multiple(() =>
        {
            Assert.That(Get(r, "ATGn{0}CAT").ExpectedFrequency, Is.EqualTo(0.01));
            Assert.That(Get(r, "ATGn{0}CAT").ExpectedFromMonads, Is.False);
            Assert.That(r.TestedPatterns, Is.EqualTo(1)); // other observed dyads: missing value, max spacing ≤ 20 → not tested
            Assert.That(r.Dyads.Where(d => d.Pattern != "ATGn{0}CAT").All(d => double.IsNaN(d.OccurrenceProbability)), Is.True);
            // max spacing > 20: missing values fall back to the monad formula (RSAT CalcProba).
            Assert.That(Get(big, "CATn{1}CAT").ExpectedFrequency, Is.EqualTo(Get(monads, "CATn{1}CAT").ExpectedFrequency));
            Assert.That(Get(big, "CATn{1}CAT").ExpectedFromMonads, Is.True);
        });
    }

    [Test]
    public void MonadFrequencyTable_EqualToInputFrequencies_ReproducesDefault()
    {
        var counts = new Dictionary<string, double>();
        for (int p = 0; p + 3 <= T2.Length; p++)
            counts[T2.Substring(p, 3)] = counts.GetValueOrDefault(T2.Substring(p, 3)) + 1.0 / 61;
        var a = Run(new DyadAnalysisOptions { MonadLength = 3, MaxSpacing = 2, MinCount = 2 }, T2);
        var b = Run(new DyadAnalysisOptions { MonadLength = 3, MaxSpacing = 2, MinCount = 2, Background = DyadBackgroundModel.MonadFrequencies(counts) }, T2);
        Assert.That(b.Dyads.Select(d => d.ExpectedFrequency), Is.EqualTo(a.Dyads.Select(d => d.ExpectedFrequency)).Within(1e-15));
    }

    [Test]
    public void Validation()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => Run(new DyadAnalysisOptions { MonadLength = 0 }, T2), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => Run(new DyadAnalysisOptions { MinSpacing = 3, MaxSpacing = 2 }, T2), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => MotifFinder.AnalyzeDyads(null!), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => Run(new DyadAnalysisOptions { MonadLength = 2, Background = DyadBackgroundModel.MonadFrequencies(new Dictionary<string, double> { ["ACG"] = 1 }) }, T2),
                NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => DyadBackgroundModel.DyadFrequencies(new Dictionary<string, double> { ["acg"] = 1 }), NUnit.Framework.Throws.ArgumentException);
            Assert.That(Run(new DyadAnalysisOptions(), "ACG").Dyads, Is.Empty);
        });
    }
}
