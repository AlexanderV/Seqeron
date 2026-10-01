// MOTIF-DISCOVER-001 (audit group E) — RSAT oligo-analysis options: -zscore, -pseudo, -oneN/-onedeg, -lexicon,
// -calibN/-calib1, multi-sequence input; MOTIF-SHARED-001 -pseudo on the matching-sequence statistics.
// Evidence: docs/Evidence/MOTIF-DISCOVER-001-Evidence.md; TestSpec: tests/TestSpecs/MOTIF-DISCOVER-001.md
// Source: rsa-tools/rsat-code master 10043f2 — perl-scripts/oligo-analysis v1.169 (CountOligos, Degenerate,
//         CalcSubWordFrequencies, CalcExpected incl. the lexicon and -pseudo blocks, CalcOverlapCoefficient, CalcZscore,
//         CalcProba), lib/RSA.lib ReadCalibration, lib/RSA.seq.lib OverlapCoeff, lib/RSA.disco.lib NbPossibleOligos,
//         lib/RSAT/stats.pm sum_of_poisson / sum_of_negbin2.
// Reference values: RSAT oligo-analysis itself run from that clone (perl 5.38) with a dump of every pattern field at
//         %.17g before PrintResult. Two RSAT defects are bypassed in the oracle and stated per test:
//         (1) sum_of_poisson keeps $prev_value as a package global between calls, so a Poisson tail can stop after its
//             first terms (oracle: same code with a fresh lexical; the exact value equals scipy / mpmath);
//         (2) -oneN/-onedeg produce no output in RSAT 1.169 (Degenerate reads the never-filled global %IUPAC and stores
//             scalars in %patterns) — the oracle fixes only those lines, defines the IUPAC residue probabilities before
//             the overlap coefficient and sets sum_occurrences to the number of windows.
//         Negative-binomial terms are rounded to 5 significant digits by RSAT (LogToEng); exact values are locked and
//         the RSAT-printed value is quoted.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class MotifFinder_OligoAnalysisOptions_Tests
{
    // RSAT input t2.fa (one sequence, 63 nt).
    private const string T2 = "ATGCATGCATGCAAATTTGGGCCCATGCTTAGCGGATCCATGCATGCTTTAAACGTACGTAGC";

    private const double Rel = 1e-10; // percent units (1e-12 relative)

    // RSAT calibration file cal3.tab (fit-distribution layout: pattern [id] avg std var repet fitted_distrib).
    private const string Cal3 =
        "; pattern\tavg\tstd\tvar\trepet\tfitted_distrib\n" +
        "atg\t0.9\t1.2\t1.44\t1000\tnegbin\n" +
        "tgc\ttgc|gca\t1.5\t1.1\t1.21\t1000\tpoisson\n" +
        "cat\t2.0\t1.5\t2.25\t1000\tnegbin\n" +
        "gca\t0.8\t0.8\t0.64\t1000\tpoisson\n";

    private static OligoStatistics Get(OligoAnalysisReport r, string pattern) => r.Patterns.Single(p => p.Pattern == pattern);

    private static OligoAnalysisReport Run(int k, OligoAnalysisOptions options, params string[] sequences)
        => MotifFinder.AnalyzeOligos(sequences.Select(s => new DnaSequence(s)), k, options);

    #region Default options = DiscoverMotifs

    private static readonly object[] EquivalenceCases =
    {
        new object[] { "equi", OligoStrandMode.Single, true, 4, 2 },
        new object[] { "input", OligoStrandMode.Both, true, 4, 2 },
        new object[] { "input", OligoStrandMode.Both, false, 4, 1 },
        new object[] { "markov1", OligoStrandMode.Single, true, 4, 2 },
        new object[] { "markov2", OligoStrandMode.Both, false, 5, 1 },
        new object[] { "lexicon", OligoStrandMode.Single, true, 4, 2 },
    };

    private static OligoBackgroundModel Bg(string name) => name switch
    {
        "equi" => OligoBackgroundModel.Equiprobable,
        "input" => OligoBackgroundModel.BernoulliFromInput,
        "markov1" => OligoBackgroundModel.MarkovFromInput(1),
        "markov2" => OligoBackgroundModel.MarkovFromInput(2),
        _ => OligoBackgroundModel.Lexicon,
    };

    [TestCaseSource(nameof(EquivalenceCases))]
    public void AnalyzeOligos_DefaultOptions_BitIdenticalToDiscoverMotifs(string bg, OligoStrandMode strands, bool overlapping, int k, int minCount)
    {
        var a = MotifFinder.DiscoverMotifs(new DnaSequence(T2), k, minCount, Bg(bg), strands, overlapping);
        var b = Run(k, new OligoAnalysisOptions { Background = Bg(bg), Strands = strands, CountOverlapping = overlapping, MinCount = minCount }, T2);

        Assert.That(b.Patterns, Has.Count.EqualTo(a.Motifs.Count));
        Assert.Multiple(() =>
        {
            Assert.That(b.TestedPatterns, Is.EqualTo(a.TestedPatterns));
            Assert.That(b.TotalOccurrences, Is.EqualTo(a.TotalOccurrences));
            Assert.That(b.PossibleOligos, Is.EqualTo(a.PossibleOligos));
            for (int i = 0; i < a.Motifs.Count; i++)
            {
                var x = a.Motifs[i];
                var y = b.Patterns[i];
                Assert.That(y.Pattern, Is.EqualTo(x.Sequence));
                Assert.That(y.ReverseComplement, Is.EqualTo(x.ReverseComplement));
                Assert.That(y.Occurrences, Is.EqualTo(x.Count));
                Assert.That(y.Positions.Select(p => p.Position), Is.EqualTo(x.Positions));
                Assert.That(y.ExpectedFrequency, Is.EqualTo(x.ExpectedFrequency));
                Assert.That(y.ExpectedOccurrences, Is.EqualTo(x.ExpectedOccurrences));
                Assert.That(y.Ratio, Is.EqualTo(x.Ratio));
                Assert.That(y.OccurrenceProbability, Is.EqualTo(x.OccurrenceProbability));
                Assert.That(y.OccurrenceEValue, Is.EqualTo(x.OccurrenceEValue));
                Assert.That(y.OccurrenceSignificance, Is.EqualTo(x.OccurrenceSignificance));
            }
        });
    }

    #endregion

    #region -zscore

    // RSAT: oligo-analysis -i t2.fa -l 4 -1str -bg equi -return occ,proba,zscore -lth occ 2
    //   atgc exp_var 0.22613525390625 ovlp 1 zscore 12.124455829017547; tgca zscore 7.9186825333149837;
    //   tgct ovlp 1.015625 exp_var 0.23345947265625 zscore 3.6542034172140836 (period 3: 1 + q(t)q(g)q(c) = 1 + 1/64).
    [Test]
    public void ZScore_Equiprobable_EqualsRsat()
    {
        var r = Run(4, new OligoAnalysisOptions { Background = OligoBackgroundModel.Equiprobable, MinCount = 2 }, T2);
        var atgc = Get(r, "ATGC");
        var tgct = Get(r, "TGCT");
        Assert.Multiple(() =>
        {
            Assert.That(r.TestedPatterns, Is.EqualTo(10));
            Assert.That(atgc.ExpectedVariance, Is.EqualTo(0.22613525390625).Within(Rel).Percent);
            Assert.That(atgc.OverlapCoefficient, Is.EqualTo(1.0));
            Assert.That(atgc.ZScore, Is.EqualTo(12.124455829017547).Within(Rel).Percent);
            Assert.That(Get(r, "TGCA").ZScore, Is.EqualTo(7.9186825333149837).Within(Rel).Percent);
            Assert.That(tgct.OverlapCoefficient, Is.EqualTo(1.015625));
            Assert.That(tgct.ExpectedVariance, Is.EqualTo(0.23345947265625).Within(Rel).Percent);
            Assert.That(tgct.ZScore, Is.EqualTo(3.6542034172140836).Within(Rel).Percent);
            Assert.That(tgct.FittedDistribution, Is.EqualTo(OligoFittedDistribution.Binomial));
        });
    }

    // RSAT: -l 4 -2str -noov -return occ,proba,zscore -lth occ 2 (input Bernoulli): with -noov exp_var = exp_occ.
    //   atgc|gcat occ 6, overlaps 3, exp_occ = exp_var 0.46662655992102053, zscore 8.1003774086907097,
    //   occ_P 7.7298916348290109e-06; catg (palindrome) zscore 9.8684030823463527; 11 tested.
    // RSAT: -l 4 -2str -ovlp ...: atgc exp_var 0.43396550795746169 zscore 12.953679437171939.
    [Test]
    public void ZScore_BothStrands_NoOverlapAndOverlap_EqualRsat()
    {
        var noov = Run(4, new OligoAnalysisOptions { Strands = OligoStrandMode.Both, CountOverlapping = false, MinCount = 2 }, T2);
        var ovlp = Run(4, new OligoAnalysisOptions { Strands = OligoStrandMode.Both, MinCount = 2 }, T2);
        var a = Get(noov, "ATGC");
        Assert.Multiple(() =>
        {
            Assert.That(noov.TestedPatterns, Is.EqualTo(11));
            Assert.That(a.Occurrences, Is.EqualTo(6));
            Assert.That(a.Overlaps, Is.EqualTo(3));
            Assert.That(a.ExpectedVariance, Is.EqualTo(0.46662655992102053).Within(Rel).Percent);
            Assert.That(a.ExpectedVariance, Is.EqualTo(a.ExpectedOccurrences));
            Assert.That(a.ZScore, Is.EqualTo(8.1003774086907097).Within(Rel).Percent);
            Assert.That(a.OccurrenceProbability, Is.EqualTo(7.7298916348290109e-06).Within(Rel).Percent);
            Assert.That(Get(noov, "CATG").ZScore, Is.EqualTo(9.8684030823463527).Within(Rel).Percent);
            Assert.That(Get(ovlp, "ATGC").ExpectedVariance, Is.EqualTo(0.43396550795746169).Within(Rel).Percent);
            Assert.That(Get(ovlp, "ATGC").ZScore, Is.EqualTo(12.953679437171939).Within(Rel).Percent);
        });
    }

    #endregion

    #region -pseudo

    // RSAT: -l 4 -1str -markov 1 -pseudo 0.1 -return occ,proba,zscore -lth occ 2
    //   atgc exp_freq 0.030014297368480814 = 0.9·p_Markov + 0.1/256, exp_var 1.3143944969861079,
    //   occ_P 0.0091641466020121916, occ_E 0.091641466020121909, zscore 3.6626693266117103;
    //   tgct exp_freq 0.0087839988377362288, zscore 2.0792489512372634.
    [Test]
    public void Pseudo_Markov1_EqualsRsat()
    {
        var r = Run(4, new OligoAnalysisOptions { Background = OligoBackgroundModel.MarkovFromInput(1), PseudoFrequency = 0.1, MinCount = 2 }, T2);
        var atgc = Get(r, "ATGC");
        Assert.Multiple(() =>
        {
            Assert.That(r.PossibleOligos, Is.EqualTo(256));
            Assert.That(atgc.ExpectedFrequency, Is.EqualTo(0.030014297368480814).Within(Rel).Percent);
            Assert.That(atgc.ExpectedVariance, Is.EqualTo(1.3143944969861079).Within(Rel).Percent);
            Assert.That(atgc.OccurrenceProbability, Is.EqualTo(0.0091641466020121916).Within(Rel).Percent);
            Assert.That(atgc.OccurrenceEValue, Is.EqualTo(0.091641466020121909).Within(Rel).Percent);
            Assert.That(atgc.ZScore, Is.EqualTo(3.6626693266117103).Within(Rel).Percent);
            Assert.That(Get(r, "TGCT").ExpectedFrequency, Is.EqualTo(0.0087839988377362288).Within(Rel).Percent);
            Assert.That(Get(r, "TGCT").ZScore, Is.EqualTo(2.0792489512372634).Within(Rel).Percent);
        });
    }

    // RSAT: -l 4 -2str -bg equi -pseudo 0.2 -lth occ 2: NPO 136, ψ/NPO added per strand before the pair sum.
    //   atgc|gcat exp_freq 0.0091911764705882356 = 2·(0.8/256 + 0.2/136), occ_P 4.5344418714526228e-09,
    //   occ_sig 7.1673849020433105; catg (palindrome) exp_freq 0.0045955882352941178, occ_P 9.0707447416984962e-06.
    [Test]
    public void Pseudo_BothStrands_AddsPerStrandBeforePairSum_EqualsRsat()
    {
        var r = Run(4, new OligoAnalysisOptions { Background = OligoBackgroundModel.Equiprobable, Strands = OligoStrandMode.Both, PseudoFrequency = 0.2, MinCount = 2 }, T2);
        Assert.Multiple(() =>
        {
            Assert.That(r.PossibleOligos, Is.EqualTo(136));
            Assert.That(r.TestedPatterns, Is.EqualTo(15));
            Assert.That(Get(r, "ATGC").ExpectedFrequency, Is.EqualTo(0.0091911764705882356).Within(Rel).Percent);
            Assert.That(Get(r, "ATGC").OccurrenceProbability, Is.EqualTo(4.5344418714526228e-09).Within(Rel).Percent);
            Assert.That(Get(r, "ATGC").OccurrenceSignificance, Is.EqualTo(7.1673849020433105).Within(1e-11));
            Assert.That(Get(r, "CATG").ExpectedFrequency, Is.EqualTo(0.0045955882352941178).Within(Rel).Percent);
            Assert.That(Get(r, "CATG").OccurrenceProbability, Is.EqualTo(9.0707447416984962e-06).Within(Rel).Percent);
        });
    }

    // RSAT: oligo-analysis -i ms.fa -l 4 -1str -markov 1 -pseudo 0.1 -return mseq,proba -lth mseq 2
    //   acgt exp_freq 0.040055603743779117 exp_ms 1.1451050486269749 ms_P 0.073696629017800608 ms_E 18.866337028556956;
    //   tagc ms_P 0.0042797192185131767.
    [Test]
    public void FindSharedMotifs_Pseudo_EqualsRsatMseq()
    {
        var seqs = new[] { "ACGTACGTTAGC", "TTACGTAGCAAC", "GGTAGCACGTTT", "CATTTTACG" }.Select(s => new DnaSequence(s));
        var r = MotifFinder.FindSharedMotifs(seqs, 4, 2, OligoBackgroundModel.MarkovFromInput(1), OligoStrandMode.Single, 0.1);
        var acgt = r.Motifs.Single(m => m.Sequence == "ACGT");
        Assert.Multiple(() =>
        {
            Assert.That(acgt.ExpectedFrequency, Is.EqualTo(0.040055603743779117).Within(Rel).Percent);
            Assert.That(acgt.ExpectedMatchingSequences, Is.EqualTo(1.1451050486269749).Within(Rel).Percent);
            Assert.That(acgt.MatchingSequenceProbability, Is.EqualTo(0.073696629017800608).Within(Rel).Percent);
            Assert.That(acgt.MatchingSequenceEValue, Is.EqualTo(18.866337028556956).Within(Rel).Percent);
            Assert.That(r.Motifs.Single(m => m.Sequence == "TAGC").MatchingSequenceProbability,
                Is.EqualTo(0.0042797192185131767).Within(Rel).Percent);
        });
    }

    [Test]
    public void FindSharedMotifs_PseudoZero_EqualsPlainOverload()
    {
        var seqs = new[] { "ACGTACGTTAGC", "TTACGTAGCAAC", "GGTAGCACGTTT", "CATTTTACG" }.Select(s => new DnaSequence(s)).ToList();
        var a = MotifFinder.FindSharedMotifs(seqs, 3, 2, OligoBackgroundModel.BernoulliFromInput, OligoStrandMode.Both);
        var b = MotifFinder.FindSharedMotifs(seqs, 3, 2, OligoBackgroundModel.BernoulliFromInput, OligoStrandMode.Both, 0.0);
        Assert.That(b.Motifs.Select(m => (m.Sequence, m.MatchingSequenceProbability)),
            Is.EqualTo(a.Motifs.Select(m => (m.Sequence, m.MatchingSequenceProbability))));
    }

    #endregion

    #region -lexicon

    // RSAT: -l 4 -1str -lexicon -return occ,proba,zscore -lth occ 2
    //   atgc exp_freq 0.025 segments "a tgc 0.25 0.1", occ_P 0.003853224500294343, zscore 4.1736500618415135;
    //   tgct exp_freq 0.028333333333333335 segments "tgc t 0.1 0.2833333333333", ovlp 1.0152971241406616 (input
    //   composition), zscore 0.26126395937379832; gcat exp_freq 0.019444444444444445.
    [Test]
    public void Lexicon_EqualsRsat()
    {
        var r = Run(4, new OligoAnalysisOptions { Background = OligoBackgroundModel.Lexicon, MinCount = 2 }, T2);
        var atgc = Get(r, "ATGC");
        var tgct = Get(r, "TGCT");
        Assert.Multiple(() =>
        {
            Assert.That(atgc.ExpectedFrequency, Is.EqualTo(0.025).Within(Rel).Percent);
            Assert.That(atgc.LexiconSegmentation, Is.EqualTo(new LexiconSegmentation("A", "TGC", 0.25, 0.1)));
            Assert.That(atgc.OccurrenceProbability, Is.EqualTo(0.003853224500294343).Within(Rel).Percent);
            Assert.That(atgc.ZScore, Is.EqualTo(4.1736500618415135).Within(Rel).Percent);
            Assert.That(tgct.ExpectedFrequency, Is.EqualTo(0.028333333333333335).Within(Rel).Percent);
            Assert.That(tgct.LexiconSegmentation!.Value.Prefix, Is.EqualTo("TGC"));
            Assert.That(tgct.LexiconSegmentation!.Value.SuffixFrequency, Is.EqualTo(0.2833333333333).Within(1e-12));
            Assert.That(tgct.OverlapCoefficient, Is.EqualTo(1.0152971241406616).Within(Rel).Percent);
            Assert.That(tgct.ZScore, Is.EqualTo(0.26126395937379832).Within(Rel).Percent);
            Assert.That(Get(r, "GCAT").ExpectedFrequency, Is.EqualTo(0.019444444444444445).Within(Rel).Percent);
        });
    }

    // RSAT: -l 5 -1str -lexicon -pseudo 0.01 -lth occ 2: NPO 1024, 8 tested;
    //   atgca exp_freq 0.025605858701702098 occ_P 0.064536489643818168 (segments atgc|a);
    //   catgc exp_freq 0.023899452496588626 occ_P 0.013413523853406996 occ_E 0.10730819082725597.
    [Test]
    public void Lexicon_WithPseudo_EqualsRsat()
    {
        var r = Run(5, new OligoAnalysisOptions { Background = OligoBackgroundModel.Lexicon, PseudoFrequency = 0.01, MinCount = 2 }, T2);
        Assert.Multiple(() =>
        {
            Assert.That(r.TestedPatterns, Is.EqualTo(8));
            Assert.That(Get(r, "ATGCA").ExpectedFrequency, Is.EqualTo(0.025605858701702098).Within(Rel).Percent);
            Assert.That(Get(r, "ATGCA").OccurrenceProbability, Is.EqualTo(0.064536489643818168).Within(Rel).Percent);
            Assert.That(Get(r, "ATGCA").LexiconSegmentation!.Value.Prefix, Is.EqualTo("ATGC"));
            Assert.That(Get(r, "CATGC").ExpectedFrequency, Is.EqualTo(0.023899452496588626).Within(Rel).Percent);
            Assert.That(Get(r, "CATGC").OccurrenceEValue, Is.EqualTo(0.10730819082725597).Within(Rel).Percent);
        });
    }

    [Test]
    public void Lexicon_NeedsWordLengthTwo_AndBothStrandsSumsMembers()
    {
        Assert.That(() => MotifFinder.DiscoverMotifs(new DnaSequence(T2), 1, 1, OligoBackgroundModel.Lexicon),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());

        // Documented -2str rule: exp_freq(W|W') = exp_freq(W) + exp_freq(W') with single-strand lexicon tables.
        var one = Run(4, new OligoAnalysisOptions { Background = OligoBackgroundModel.Lexicon }, T2);
        var two = Run(4, new OligoAnalysisOptions { Background = OligoBackgroundModel.Lexicon, Strands = OligoStrandMode.Both }, T2);
        double expected = Get(one, "ATGC").ExpectedFrequency + Get(one, "GCAT").ExpectedFrequency;
        Assert.That(Get(two, "ATGC").ExpectedFrequency, Is.EqualTo(expected).Within(1e-13).Percent);
    }

    #endregion

    #region -oneN / -onedeg

    // RSAT 1.169: oligo-analysis -i t2.fa -l 3 -1str -onedeg -return occ,proba → "oligomers tested for significance 0",
    // no row (Degenerate bug). Oracle (Degenerate fixed, see header): -l 3 -1str -onedeg -return occ,proba,zscore
    // -lth occ 6 → NPO 528 (= 3·11·4²), 71 tested; ayg occ 8 exp_freq 0.03071422572556359 (= q(a)·(q(c)+q(t))·q(g)),
    // zscore 5.0517050514404822, occ_P 0.00054574388255023234, occ_E 0.038747815661066498.
    [Test]
    public void OneDegenerate_EqualsFixedRsatOracle()
    {
        var r = Run(3, new OligoAnalysisOptions { Degeneracy = OligoDegeneracy.OneDegenerate, MinCount = 6 }, T2);
        var ayg = Get(r, "AYG");
        Assert.Multiple(() =>
        {
            Assert.That(r.PossibleOligos, Is.EqualTo(528));
            Assert.That(r.TestedPatterns, Is.EqualTo(71));
            Assert.That(r.TotalOccurrences, Is.EqualTo(61));
            Assert.That(ayg.Occurrences, Is.EqualTo(8));
            Assert.That(ayg.Positions, Has.Count.EqualTo(8));
            Assert.That(ayg.ExpectedFrequency, Is.EqualTo(0.03071422572556359).Within(Rel).Percent);
            Assert.That(ayg.ZScore, Is.EqualTo(5.0517050514404822).Within(Rel).Percent);
            Assert.That(ayg.OccurrenceProbability, Is.EqualTo(0.00054574388255023234).Within(Rel).Percent);
            Assert.That(ayg.OccurrenceEValue, Is.EqualTo(0.038747815661066498).Within(Rel).Percent);
        });
    }

    // Oracle: -l 4 -2str -oneN -lth occ 6 → NPO 136 (RSAT formula: 4·4³ = 256 → 256 − (256 − 16)/2), 8 tested;
    //   atgn|ncat occ 11 exp_freq 0.032663859194471435 occ_P 3.5186386968671786e-06 zscore 7.6852328133423624.
    [Test]
    public void OneN_BothStrands_EqualsFixedRsatOracle()
    {
        var r = Run(4, new OligoAnalysisOptions { Degeneracy = OligoDegeneracy.OneN, Strands = OligoStrandMode.Both, MinCount = 6 }, T2);
        var atgn = Get(r, "ATGN");
        Assert.Multiple(() =>
        {
            Assert.That(r.PossibleOligos, Is.EqualTo(136));
            Assert.That(r.TestedPatterns, Is.EqualTo(8));
            Assert.That(atgn.ReverseComplement, Is.EqualTo("NCAT"));
            Assert.That(atgn.Occurrences, Is.EqualTo(11));
            Assert.That(atgn.ExpectedFrequency, Is.EqualTo(0.032663859194471435).Within(Rel).Percent);
            Assert.That(atgn.OccurrenceProbability, Is.EqualTo(3.5186386968671786e-06).Within(Rel).Percent);
            Assert.That(atgn.ZScore, Is.EqualTo(7.6852328133423624).Within(Rel).Percent);
        });
    }

    [Test]
    public void OneN_OccurrencesAreSumOfMatchingWords()
    {
        var plain = Run(3, new OligoAnalysisOptions(), T2);
        var deg = Run(3, new OligoAnalysisOptions { Degeneracy = OligoDegeneracy.OneN }, T2);
        int expected = plain.Patterns.Where(p => p.Pattern.StartsWith("AT", StringComparison.Ordinal)).Sum(p => p.Occurrences);
        Assert.That(Get(deg, "ATN").Occurrences, Is.EqualTo(expected));
    }

    #endregion

    #region -calibN / -calib1

    // RSAT: oligo-analysis -i t2.fa -l 3 -1str -calibN cal3.tab -return occ,proba,zscore -lth occ 1 → 4 tested.
    //   atg: mean 0.9 < var 1.44 → negative binomial, occ 6: exact 0.0067142352257039137 (mpmath 0.006714235225703917789;
    //        RSAT prints 0.0067141894090253264 — 5-digit LogToEng terms); zscore (6 − 0.9)/1.2 = 4.25.
    //   cat: negbin occ 5 → 0.063244578406644242 (RSAT 0.063244096728178351).
    //   gca: Poisson(0.8) occ 4 → 0.009079857800153978, zscore 4 (RSAT identical).
    //   tgc: Poisson(1.5) occ 6 → exact 0.0044559807752478468 (scipy 0.004455980775247892); RSAT 1.169 prints
    //        0.003529988861723204 because sum_of_poisson's $prev_value is a package global left by the previous call.
    //   exp_freq(tgc) = 1.5/61 = 0.024590163934426229; zscore 4.0909090909090908.
    [Test]
    public void CalibrationPerSet_PoissonAndNegativeBinomial_EqualRsatExactOracle()
    {
        var cal = OligoCalibration.Parse(Cal3);
        var r = Run(3, new OligoAnalysisOptions { Calibration = cal }, T2);
        var atg = Get(r, "ATG");
        var tgc = Get(r, "TGC");
        Assert.Multiple(() =>
        {
            Assert.That(r.TestedPatterns, Is.EqualTo(4));
            Assert.That(atg.FittedDistribution, Is.EqualTo(OligoFittedDistribution.NegativeBinomial));
            Assert.That(atg.OccurrenceProbability, Is.EqualTo(0.0067142352257039137).Within(Rel).Percent);
            Assert.That(atg.OccurrenceProbability, Is.EqualTo(0.0067141894090253264).Within(1e-3).Percent); // RSAT printed (rounded terms)
            Assert.That(atg.ZScore, Is.EqualTo(4.25).Within(Rel).Percent);
            Assert.That(atg.OccurrenceEValue, Is.EqualTo(0.026856940902815655).Within(Rel).Percent);
            Assert.That(Get(r, "CAT").OccurrenceProbability, Is.EqualTo(0.063244578406644242).Within(Rel).Percent);
            Assert.That(Get(r, "GCA").FittedDistribution, Is.EqualTo(OligoFittedDistribution.Poisson));
            Assert.That(Get(r, "GCA").OccurrenceProbability, Is.EqualTo(0.009079857800153978).Within(Rel).Percent);
            Assert.That(Get(r, "GCA").ZScore, Is.EqualTo(4.0).Within(Rel).Percent);
            Assert.That(tgc.OccurrenceProbability, Is.EqualTo(0.0044559807752478468).Within(Rel).Percent);
            Assert.That(tgc.ExpectedFrequency, Is.EqualTo(0.024590163934426229).Within(Rel).Percent);
            Assert.That(tgc.ExpectedOccurrences, Is.EqualTo(1.5));
            Assert.That(tgc.ExpectedVariance, Is.EqualTo(1.21));
            Assert.That(tgc.ZScore, Is.EqualTo(4.0909090909090908).Within(Rel).Percent);
        });
    }

    // Oracle (exact Poisson/negbin): -i two.fa (2 sequences) -l 3 -2str -noov -calib1 cal3.tab -lth occ 1 → 2 tested, NPO 32;
    //   mean/var × 2 sequences; atg|cat occ 6 overlaps 5, exp_occ 1.8, exp_var 2.88, exp_freq 2·1.8/39 = 0.092307692307692313
    //   (39 = 51 positions − 12 forbidden; RSAT copies the kept member's value to its partner before summing),
    //   negbin occ_P 0.03602153062820438, zscore 2.4748737341529163;
    //   gca|tgc Poisson(1.6) occ_P 0.0060402911115813741, zscore 3.8890872965260117.
    [Test]
    public void CalibrationPerSequence_TwoSequences_BothStrandsNoOverlap_EqualsRsatExactOracle()
    {
        var cal = OligoCalibration.Parse(Cal3, OligoCalibrationMode.PerSequence);
        var r = Run(3, new OligoAnalysisOptions { Calibration = cal, Strands = OligoStrandMode.Both, CountOverlapping = false },
            "ATGCATGCATGCAAATTTGGGCC", "CATGCTTAGCGGATCCATGCATGCTTTAAACG");
        var atg = Get(r, "ATG");
        Assert.Multiple(() =>
        {
            Assert.That(r.SequenceCount, Is.EqualTo(2));
            Assert.That(r.PossiblePositions, Is.EqualTo(51));
            Assert.That(r.TestedPatterns, Is.EqualTo(2));
            Assert.That(r.PossibleOligos, Is.EqualTo(32));
            Assert.That(atg.Occurrences, Is.EqualTo(6));
            Assert.That(atg.Overlaps, Is.EqualTo(5));
            Assert.That(atg.ExpectedOccurrences, Is.EqualTo(1.8).Within(Rel).Percent);
            Assert.That(atg.ExpectedVariance, Is.EqualTo(2.88).Within(Rel).Percent);
            Assert.That(atg.ExpectedFrequency, Is.EqualTo(0.092307692307692313).Within(Rel).Percent);
            Assert.That(atg.OccurrenceProbability, Is.EqualTo(0.03602153062820438).Within(Rel).Percent);
            Assert.That(atg.ZScore, Is.EqualTo(2.4748737341529163).Within(Rel).Percent);
            Assert.That(Get(r, "GCA").OccurrenceProbability, Is.EqualTo(0.0060402911115813741).Within(Rel).Percent);
            Assert.That(Get(r, "GCA").ZScore, Is.EqualTo(3.8890872965260117).Within(Rel).Percent);
        });
    }

    [Test]
    public void Calibration_Parse_InfersReverseComplements_AndRejectsBadInput()
    {
        var cal = OligoCalibration.Parse(Cal3);
        Assert.Multiple(() =>
        {
            Assert.That(cal.WordLength, Is.EqualTo(3));
            Assert.That(cal.Entries["CAT"], Is.EqualTo((2.0, 2.25)));     // own line kept (ATG's rc)
            Assert.That(cal.Entries["GCA"], Is.EqualTo((0.8, 0.64)));
            Assert.That(cal.Entries.ContainsKey("CGT"), Is.False);
            Assert.That(() => OligoCalibration.Parse("acg\t1.0\n"), NUnit.Framework.Throws.TypeOf<FormatException>());
            Assert.That(() => OligoCalibration.Parse("; only comments\n"), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => OligoCalibration.Parse("acg\t1\t1\t1\nac\t1\t1\t1\n"), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => Run(4, new OligoAnalysisOptions { Calibration = cal }, T2), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => Run(3, new OligoAnalysisOptions { Calibration = cal, Degeneracy = OligoDegeneracy.OneN }, T2),
                NUnit.Framework.Throws.ArgumentException);
        });
    }

    [Test]
    public void Calibration_ReverseComplementInferredFromPartner()
    {
        var cal = OligoCalibration.Parse("aac\t1.5\t1\t1\n");
        Assert.That(cal.Entries["GTT"], Is.EqualTo((1.5, 1.0)));
    }

    #endregion

    [Test]
    public void Options_Validation()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => Run(3, new OligoAnalysisOptions { PseudoFrequency = 1.5 }, T2), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => Run(0, new OligoAnalysisOptions(), T2), NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => MotifFinder.AnalyzeOligos(null!, 3), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => MotifFinder.AnalyzeOligos(new DnaSequence[] { null! }, 3), NUnit.Framework.Throws.ArgumentException);
            Assert.That(() => Run(4, new OligoAnalysisOptions { Background = OligoBackgroundModel.MarkovFromInput(3) }, T2),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>());
        });
    }

    [Test]
    public void ShortSequences_ContributeNoWindows_AndEmptyInputGivesNoPatterns()
    {
        var r = Run(4, new OligoAnalysisOptions(), "ACG", T2);
        var empty = Run(4, new OligoAnalysisOptions(), "ACG");
        Assert.Multiple(() =>
        {
            Assert.That(r.SequenceCount, Is.EqualTo(2));
            Assert.That(r.PossiblePositions, Is.EqualTo(60));
            Assert.That(r.Patterns.SelectMany(p => p.Positions).All(p => p.SequenceIndex == 1), Is.True);
            Assert.That(empty.Patterns, Is.Empty);
            Assert.That(empty.TestedPatterns, Is.Zero);
        });
    }
}
