// MOTIF-DISCOVER-001 / MOTIF-SHARED-001 — RSAT oligo-analysis significance, backgrounds and strands
// Evidence: docs/Evidence/MOTIF-DISCOVER-001-Evidence.md, docs/Evidence/MOTIF-SHARED-001-Evidence.md
// TestSpec: tests/TestSpecs/MOTIF-DISCOVER-001.md, tests/TestSpecs/MOTIF-SHARED-001.md
// Source: van Helden J, André B, Collado-Vides J (1998) J Mol Biol 281:827-842; RSAT rsa-tools/rsat-code master
//         10043f2 (perl-scripts/oligo-analysis v1.169, lib/RSA.disco.lib, lib/RSA.lib, lib/RSAT/stats.pm,
//         lib/RSAT/MarkovModel.pm).
// Reference values: RSAT oligo-analysis itself, run from a git clone with perl 5.38 (RSAT=<clone>, perl -I
//         perl-scripts/lib perl-scripts/oligo-analysis ...), with one debug line added after MultiTestCorrections
//         printing exp_freq/occ_P/occ_E/occ_sig/exp_ms/ms_P/ms_E/ms_sig with %.17g. Where RSAT code and its manual
//         disagree (-2str with a strand-asymmetric background) the value is an independent Python port of the
//         documented formula (scipy.stats.binom.sf) — stated per test.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class MotifFinder_OligoAnalysis_Tests
{
    // RSAT input t2.fa (one sequence, 63 nt).
    private const string T2 = "ATGCATGCATGCAAATTTGGGCCCATGCTTAGCGGATCCATGCATGCTTTAAACGTACGTAGC";

    // RSAT input t3.fa (35 nt) with homopolymer runs for -noov.
    private const string T3 = "AAAAAAATTTTTTGCGCGCACGTACGTAAAATTTT";

    // RSAT input ms.fa (4 sequences; the MOTIF-SHARED-001 mseq reference set).
    private static readonly string[] Ms = { "ACGTACGTTAGC", "TTACGTAGCAAC", "GGTAGCACGTTT", "CATTTTACG" };

    private const double Rel = 1e-10; // relative tolerance, in percent units for NUnit (1e-12 relative)

    private static SignificantMotif Get(OligoAnalysisResult r, string word) => r.Motifs.Single(m => m.Sequence == word);

    private static SignificantSharedMotif Get(SharedMotifAnalysisResult r, string word) => r.Motifs.Single(m => m.Sequence == word);

    private static void AssertOcc(SignificantMotif m, int occ, double expFreq, double occP, double occE, double occSig)
    {
        Assert.Multiple(() =>
        {
            Assert.That(m.Count, Is.EqualTo(occ), $"{m.Sequence} occ");
            Assert.That(m.ExpectedFrequency, Is.EqualTo(expFreq).Within(Rel).Percent, $"{m.Sequence} exp_freq");
            Assert.That(m.OccurrenceProbability, Is.EqualTo(occP).Within(Rel).Percent, $"{m.Sequence} occ_P");
            Assert.That(m.OccurrenceEValue, Is.EqualTo(occE).Within(Rel).Percent, $"{m.Sequence} occ_E");
            Assert.That(m.OccurrenceSignificance, Is.EqualTo(occSig).Within(1e-11), $"{m.Sequence} occ_sig");
        });
    }

    #region DiscoverMotifs — occ_P / occ_E / occ_sig

    // RSAT: oligo-analysis -i t2.fa -l 4 -1str -bg equi -lth occ 2 -return occ,proba
    //   nb possible positions 60, total oligo occurrences 60, oligomers tested 10.
    //   atgc occ 6 occ_P 1.4844942103072834e-07 occ_E 1.4844942103072835e-06 occ_sig 5.8284214918614703
    //   tgca occ 4 occ_P 9.5343622068379445e-05 occ_E 0.00095343622068379448 occ_sig 3.0207083534113868
    //   (scipy.stats.binom.sf(5, 60, 1/256) = 1.484494210307282e-07).
    [Test]
    public void DiscoverMotifs_Equiprobable_SingleStrand_EqualsRsatOligoAnalysis()
    {
        var r = MotifFinder.DiscoverMotifs(new DnaSequence(T2), 4, 2, OligoBackgroundModel.Equiprobable);

        Assert.Multiple(() =>
        {
            Assert.That(r.TotalOccurrences, Is.EqualTo(60));
            Assert.That(r.TestedPatterns, Is.EqualTo(10));
            Assert.That(r.PossibleOligos, Is.EqualTo(256));
        });
        AssertOcc(Get(r, "ATGC"), 6, 0.00390625, 1.4844942103072834e-07, 1.4844942103072835e-06, 5.8284214918614703);
        AssertOcc(Get(r, "TGCA"), 4, 0.00390625, 9.5343622068379445e-05, 0.00095343622068379448, 3.0207083534113868);
        Assert.That(Get(r, "ATGC").Ratio, Is.EqualTo(6 / (60 / 256.0)).Within(Rel).Percent, "ratio = occ / exp_occ");
    }

    // RSAT: -l 4 -2str (default input Bernoulli, residues summed with their complements) -lth occ 2:
    //   nb possible oligomers 136, tested 15; atgc|gcat occ 9 exp_freq 0.0077771093320170084
    //   occ_P 1.0760647330191343e-09 occ_E 1.6140970995287015e-08 occ_sig 7.7920703428970466;
    //   catg (palindrome, counted once) occ 5 exp_freq 0.0038885546660085047 occ_P 4.0637121408222015e-06
    //   occ_E 6.0955682112333022e-05 occ_sig 4.2149858044834696.
    [Test]
    public void DiscoverMotifs_InputBernoulli_BothStrands_EqualsRsatOligoAnalysis()
    {
        var r = MotifFinder.DiscoverMotifs(new DnaSequence(T2), 4, 2, OligoBackgroundModel.BernoulliFromInput,
            OligoStrandMode.Both);

        Assert.Multiple(() =>
        {
            Assert.That(r.TotalOccurrences, Is.EqualTo(60));
            Assert.That(r.TestedPatterns, Is.EqualTo(15));
            Assert.That(r.PossibleOligos, Is.EqualTo(136), "(4^4 + 4^2) / 2 reverse-complement classes");
            Assert.That(Get(r, "ATGC").ReverseComplement, Is.EqualTo("GCAT"));
            Assert.That(Get(r, "CATG").ReverseComplement, Is.EqualTo("CATG"));
        });
        AssertOcc(Get(r, "ATGC"), 9, 0.0077771093320170084, 1.0760647330191343e-09, 1.6140970995287015e-08, 7.7920703428970466);
        AssertOcc(Get(r, "CATG"), 5, 0.0038885546660085047, 4.0637121408222015e-06, 6.0955682112333022e-05, 4.2149858044834696);
    }

    // RSAT: -l 4 -1str -markov 1 (no threshold, 40 tested): atgc occ 6 exp_freq 0.032915191520534237
    //   occ_P 0.013954892377198637 occ_E 0.55819569508794542 occ_sig 0.25321351719665269.
    // RSAT: -l 4 -1str -markov 2 -lth occ 3 (4 tested): gcat occ 3 exp_freq 0.055540625279942669
    //   occ_P 0.65459731021733769 occ_E 2.6183892408693508 occ_sig -0.41803420775951439.
    [Test]
    public void DiscoverMotifs_MarkovFromInput_EqualsRsatMarkovOption()
    {
        var m1 = MotifFinder.DiscoverMotifs(new DnaSequence(T2), 4, 1, OligoBackgroundModel.MarkovFromInput(1));
        var m2 = MotifFinder.DiscoverMotifs(new DnaSequence(T2), 4, 3, OligoBackgroundModel.MarkovFromInput(2));

        Assert.Multiple(() =>
        {
            Assert.That(m1.TestedPatterns, Is.EqualTo(40));
            Assert.That(m2.TestedPatterns, Is.EqualTo(4));
        });
        AssertOcc(Get(m1, "ATGC"), 6, 0.032915191520534237, 0.013954892377198637, 0.55819569508794542, 0.25321351719665269);
        AssertOcc(Get(m2, "GCAT"), 3, 0.055540625279942669, 0.65459731021733769, 2.6183892408693508, -0.41803420775951439);
    }

    // Hand derivation (RSAT MarkovModel with pseudo-frequency 0): table f(xy) = 1..16 in AA, AC, ..., TT order.
    //   P(A) = S(A)/F = 10/136, P(T|A) = 4/10, P(G|T) = 15/58, P(C|G) = 10/42 → exp_freq(ATGC) = 600/331296.
    // RSAT -bgfile (oligos, ψ = 0.01) -l 4 -1str -lth occ 5: atgc exp_freq printed 0.00184789 (segment_proba %5g),
    //   occ_P 1.8299478950713788e-09; unrounded (Python port of MarkovModel.pm): exp_freq 0.0018478880336134456,
    //   occ_P 1.8299363778610906e-09, occ_E 3.6598727557221813e-09, occ_sig 8.436534013635193 (2 tested).
    [Test]
    public void DiscoverMotifs_MarkovFromOligoTable_EqualsRsatBackgroundFile()
    {
        var table = new Dictionary<string, double>();
        int i = 1;
        foreach (char a in "ACGT")
            foreach (char b in "ACGT")
                table[$"{a}{b}"] = i++;

        var exact = MotifFinder.DiscoverMotifs(new DnaSequence(T2), 4, 5,
            OligoBackgroundModel.MarkovFromOligoFrequencies(table, pseudoFrequency: 0.0));
        var rsat = MotifFinder.DiscoverMotifs(new DnaSequence(T2), 4, 5,
            OligoBackgroundModel.MarkovFromOligoFrequencies(table));

        Assert.That(Get(exact, "ATGC").ExpectedFrequency, Is.EqualTo(600.0 / 331296).Within(Rel).Percent);
        Assert.That(Get(rsat, "ATGC").ExpectedFrequency, Is.EqualTo(0.00184789).Within(5e-4).Percent,
            "RSAT prints segment_proba with %5g (6 significant digits)");
        AssertOcc(Get(rsat, "ATGC"), 6, 0.0018478880336134456, 1.8299363778610906e-09, 3.6598727557221813e-09, 8.436534013635193);
    }

    // RSAT -noov (windows read from the 3' end): -l 2 -1str -bg equi: aa occ 5 occ_P 0.058469994275128084
    //   occ_E 0.5846999427512809 occ_sig 0.23306694868504213 (10 tested; n = 34 windows).
    //   AAAAAAA (windows 0..5) keeps 5, 3, 1 and AAAA (27..29) keeps 29, 27 → positions {1, 3, 5, 27, 29}.
    // -2str -noov: aa|tt occ 10 occ_P 0.0070505538526756187 occ_E 0.049353876968729331 occ_sig 1.3066787258719326;
    // -2str -ovlp: aa|tt occ 17 occ_P 1.2344994901541932e-07 occ_E 8.6414964310793518e-07 occ_sig 6.0634110450805023.
    [Test]
    public void DiscoverMotifs_NonOverlapping_EqualsRsatNoov()
    {
        var single = MotifFinder.DiscoverMotifs(new DnaSequence(T3), 2, 1, OligoBackgroundModel.Equiprobable,
            OligoStrandMode.Single, countOverlapping: false);
        var both = MotifFinder.DiscoverMotifs(new DnaSequence(T3), 2, 1, OligoBackgroundModel.Equiprobable,
            OligoStrandMode.Both, countOverlapping: false);
        var bothOvlp = MotifFinder.DiscoverMotifs(new DnaSequence(T3), 2, 1, OligoBackgroundModel.Equiprobable,
            OligoStrandMode.Both);

        Assert.Multiple(() =>
        {
            Assert.That(single.TotalOccurrences, Is.EqualTo(34));
            Assert.That(single.TestedPatterns, Is.EqualTo(10));
            Assert.That(Get(single, "AA").Positions, Is.EqualTo(new[] { 1, 3, 5, 27, 29 }));
            Assert.That(both.TestedPatterns, Is.EqualTo(7));
        });
        AssertOcc(Get(single, "AA"), 5, 0.0625, 0.058469994275128084, 0.5846999427512809, 0.23306694868504213);
        AssertOcc(Get(both, "AA"), 10, 0.125, 0.0070505538526756187, 0.049353876968729331, 1.3066787258719326);
        AssertOcc(Get(bothOvlp, "AA"), 17, 0.125, 1.2344994901541932e-07, 8.6414964310793518e-07, 6.0634110450805023);
    }

    // -2str with a strand-asymmetric background: RSAT manual "exp_freq(W|Wr) = exp_freq(W) + exp_freq(Wr)".
    // (RSAT code instead copies the kept member's value, 2·exp_freq(min(W, Wr)); on t2.fa -markov 1 that gives
    // aaat|attt 0.0166535195193 = 2 × 0.0083267597597, and -markov 2 even gives exp_freq 0 for observed aaag|cttt.)
    // Python port (documented sum, Markov-1 estimated from the single strand, scipy binom.sf), -l 4 -2str -lth occ 3:
    //   atgc|gcat occ 9 exp_freq 0.06289152665530649 occ_P 0.012300624704130635 occ_E 0.04920249881652254
    //   occ_sig 1.3080128404304179 (4 tested); single-strand p(ATGC) 0.032915191520534237 + p(GCAT) = the sum.
    [Test]
    public void DiscoverMotifs_MarkovBothStrands_UsesDocumentedPairSum()
    {
        var both = MotifFinder.DiscoverMotifs(new DnaSequence(T2), 4, 3, OligoBackgroundModel.MarkovFromInput(1),
            OligoStrandMode.Both);
        var single = MotifFinder.DiscoverMotifs(new DnaSequence(T2), 4, 1, OligoBackgroundModel.MarkovFromInput(1));

        Assert.That(both.TestedPatterns, Is.EqualTo(4));
        AssertOcc(Get(both, "ATGC"), 9, 0.06289152665530649, 0.012300624704130635, 0.04920249881652254, 1.3080128404304179);
        Assert.That(Get(both, "ATGC").ExpectedFrequency,
            Is.EqualTo(Get(single, "ATGC").ExpectedFrequency + Get(single, "GCAT").ExpectedFrequency).Within(Rel).Percent);
    }

    [Test]
    public void DiscoverMotifs_EquiprobableSingleStrand_RatioEqualsLegacyEnrichment()
    {
        var seq = new DnaSequence(T2);
        var legacy = MotifFinder.DiscoverMotifs(seq, 3, 2).ToDictionary(m => m.Sequence);
        var rsat = MotifFinder.DiscoverMotifs(seq, 3, 2, OligoBackgroundModel.Equiprobable);

        Assert.That(rsat.Motifs.Select(m => m.Sequence), Is.EqualTo(legacy.Keys), "same words, same first-occurrence order");
        foreach (var m in rsat.Motifs)
        {
            Assert.That(m.Ratio, Is.EqualTo(legacy[m.Sequence].Enrichment).Within(1e-10).Percent);
            Assert.That(m.Positions, Is.EqualTo(legacy[m.Sequence].Positions));
        }
    }

    [Test]
    public void DiscoverMotifs_BothStrands_PairCountIsForwardPlusReverseComplement()
    {
        var seq = new DnaSequence(T2);
        var single = MotifFinder.DiscoverMotifs(seq, 4, 1, OligoBackgroundModel.Equiprobable).Motifs
            .ToDictionary(m => m.Sequence);
        var both = MotifFinder.DiscoverMotifs(seq, 4, 1, OligoBackgroundModel.Equiprobable, OligoStrandMode.Both);

        foreach (var m in both.Motifs)
        {
            int forward = single.TryGetValue(m.Sequence, out var f) ? f.Count : 0;
            int reverse = m.ReverseComplement != m.Sequence && single.TryGetValue(m.ReverseComplement!, out var r) ? r.Count : 0;
            Assert.That(m.Count, Is.EqualTo(forward + reverse), $"{m.Sequence}|{m.ReverseComplement}");
            Assert.That(string.CompareOrdinal(m.Sequence, m.ReverseComplement), Is.LessThanOrEqualTo(0));
        }
        Assert.That(both.Motifs.Sum(m => m.Count), Is.EqualTo(60), "each window belongs to exactly one pair");
    }

    [Test]
    public void DiscoverMotifs_Significance_InvalidArguments()
    {
        var seq = new DnaSequence(T2);
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => MotifFinder.DiscoverMotifs(null!, 4, 1, OligoBackgroundModel.Equiprobable));
            Assert.Throws<ArgumentNullException>(() => MotifFinder.DiscoverMotifs(seq, 4, 1, (OligoBackgroundModel)null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.DiscoverMotifs(seq, 0, 1, OligoBackgroundModel.Equiprobable));
            // RSAT: "Markov order cannot be higher than word length - 2".
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.DiscoverMotifs(seq, 4, 1, OligoBackgroundModel.MarkovFromInput(3)));
            Assert.Throws<ArgumentOutOfRangeException>(() => OligoBackgroundModel.MarkovFromInput(-1));
            Assert.Throws<ArgumentException>(() => OligoBackgroundModel.MarkovFromOligoFrequencies(new Dictionary<string, double>()));
            Assert.Throws<ArgumentException>(() => OligoBackgroundModel.MarkovFromOligoFrequencies(new Dictionary<string, double> { ["AC"] = 1, ["A"] = 1 }));
            Assert.Throws<ArgumentException>(() => OligoBackgroundModel.MarkovFromOligoFrequencies(new Dictionary<string, double> { ["AN"] = 1 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => OligoBackgroundModel.MarkovFromOligoFrequencies(new Dictionary<string, double> { ["AC"] = -1 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => OligoBackgroundModel.MarkovFromOligoFrequencies(new Dictionary<string, double> { ["AC"] = 1 }, 1.5));
            // Table order 3 cannot score 3-mers (RSAT segment_proba: length must exceed the order).
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.DiscoverMotifs(seq, 3, 1,
                OligoBackgroundModel.MarkovFromOligoFrequencies(new Dictionary<string, double> { ["ACGT"] = 1 })));
            // ψ = 0 and a missing transition: RSAT segment_proba "null transition" → error.
            Assert.Throws<ArgumentException>(() => MotifFinder.DiscoverMotifs(seq, 2, 1,
                OligoBackgroundModel.MarkovFromOligoFrequencies(new Dictionary<string, double> { ["AA"] = 1 }, 0.0)));
        });
    }

    [Test]
    public void DiscoverMotifs_Significance_ShortSequence_NoMotifs()
    {
        var r = MotifFinder.DiscoverMotifs(new DnaSequence("ACG"), 4, 1, OligoBackgroundModel.BernoulliFromInput);

        Assert.Multiple(() =>
        {
            Assert.That(r.Motifs, Is.Empty);
            Assert.That(r.TotalOccurrences, Is.Zero);
            Assert.That(r.TestedPatterns, Is.Zero);
        });
    }

    [Test]
    public void DiscoverMotifs_LongK_SignificanceFiniteWhereProbabilityUnderflows()
    {
        // X + X with X a 600-mer, k = 600: exp_freq = 4^-600 underflows, occ = 2 of n = 601 windows.
        // ln occ_P = ln C(601,2) + 2·ln 4^-600 + 599·ln(1 - 4^-600) + ln(1 + O(600·4^-600)) = ln 180300 - 1200·ln 4 to
        // double precision; E-value multiplier = tested patterns (600 distinct windows) → occ_sig ≈ 714.4.
        var rng = new Random(5);
        string x = new(Enumerable.Range(0, 600).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
        var r = MotifFinder.DiscoverMotifs(new DnaSequence(x + x), 600, 1, OligoBackgroundModel.Equiprobable);
        var m = Get(r, x);

        double expected = -(Math.Log(601.0 * 600 / 2) - 1200 * Math.Log(4) + Math.Log(r.TestedPatterns)) / Math.Log(10);
        Assert.Multiple(() =>
        {
            Assert.That(m.Count, Is.EqualTo(2));
            Assert.That(m.OccurrenceProbability, Is.Zero, "the p-value itself is below the double range");
            Assert.That(m.OccurrenceSignificance, Is.EqualTo(expected).Within(1e-9).Percent);
            Assert.That(r.PossibleOligos, Is.EqualTo(double.PositiveInfinity));
        });
    }

    #endregion

    #region FindSharedMotifs — ms_P / ms_E / ms_sig

    // RSAT: -l 4 -1str -return occ,mseq,proba -lth mseq 2 (input Bernoulli): nb possible positions 33, 4 sequences.
    //   acgt mseq 3 exp_freq 0.0037555250723974999 exp_ms 0.12225827585979898 ms_P 0.00011159466135309564
    //   ms_E 0.028568233306392483 ms_sig 1.5441166160753867 (ms_E = ms_P × 256).
    [Test]
    public void FindSharedMotifs_InputBernoulli_SingleStrand_EqualsRsatMseq()
    {
        var r = MotifFinder.FindSharedMotifs(Ms.Select(s => new DnaSequence(s)), 4, 2, OligoBackgroundModel.BernoulliFromInput);
        var acgt = Get(r, "ACGT");

        Assert.Multiple(() =>
        {
            Assert.That(r.SequenceCount, Is.EqualTo(4));
            Assert.That(r.PossiblePositions, Is.EqualTo(33));
            Assert.That(r.PossibleOligos, Is.EqualTo(256));
            Assert.That(r.Motifs.Select(m => m.Sequence), Is.EqualTo(
                MotifFinder.FindSharedMotifs(Ms.Select(s => new DnaSequence(s)), 4, 2).Select(m => m.Sequence)),
                "same words and order as the quorum-only overload");
            Assert.That(acgt.SequenceIndices, Is.EqualTo(new[] { 0, 1, 2 }));
            Assert.That(acgt.ExpectedFrequency, Is.EqualTo(0.0037555250723974999).Within(Rel).Percent);
            Assert.That(acgt.ExpectedMatchingSequences, Is.EqualTo(0.12225827585979898).Within(Rel).Percent);
            Assert.That(acgt.MatchingSequenceProbability, Is.EqualTo(0.00011159466135309564).Within(Rel).Percent);
            Assert.That(acgt.MatchingSequenceEValue, Is.EqualTo(0.028568233306392483).Within(Rel).Percent);
            Assert.That(acgt.MatchingSequenceSignificance, Is.EqualTo(1.5441166160753867).Within(1e-11));
        });
    }

    // RSAT: -l 4 -2str -return occ,mseq,proba -lth mseq 2: NPO 136, 7 pairs.
    //   cgta|tacg mseq 3 exp_freq 0.0076207895137936294 exp_ms 0.24464786582736009 ms_P 0.00087319490302470882
    //   ms_E 0.11875450681136041 ms_sig 0.9253498996301337;
    //   acgt (palindrome) mseq 3 exp_ms 0.12401989758949705 ms_P 0.00011644903222066987 ms_sig 1.8003252080589478.
    [Test]
    public void FindSharedMotifs_InputBernoulli_BothStrands_EqualsRsatMseq()
    {
        var r = MotifFinder.FindSharedMotifs(Ms.Select(s => new DnaSequence(s)), 4, 2,
            OligoBackgroundModel.BernoulliFromInput, OligoStrandMode.Both);
        var cgta = Get(r, "CGTA");
        var acgt = Get(r, "ACGT");

        Assert.Multiple(() =>
        {
            Assert.That(r.Motifs, Has.Count.EqualTo(7));
            Assert.That(r.PossibleOligos, Is.EqualTo(136));
            Assert.That(cgta.ReverseComplement, Is.EqualTo("TACG"));
            Assert.That(cgta.SequenceIndices, Is.EqualTo(new[] { 0, 1, 3 }), "CGTA in 0, 1 and TACG in 3");
            Assert.That(cgta.ExpectedFrequency, Is.EqualTo(0.0076207895137936294).Within(Rel).Percent);
            Assert.That(cgta.ExpectedMatchingSequences, Is.EqualTo(0.24464786582736009).Within(Rel).Percent);
            Assert.That(cgta.MatchingSequenceProbability, Is.EqualTo(0.00087319490302470882).Within(Rel).Percent);
            Assert.That(cgta.MatchingSequenceEValue, Is.EqualTo(0.11875450681136041).Within(Rel).Percent);
            Assert.That(cgta.MatchingSequenceSignificance, Is.EqualTo(0.9253498996301337).Within(1e-11));
            Assert.That(acgt.ExpectedMatchingSequences, Is.EqualTo(0.12401989758949705).Within(Rel).Percent);
            Assert.That(acgt.MatchingSequenceProbability, Is.EqualTo(0.00011644903222066987).Within(Rel).Percent);
            Assert.That(acgt.MatchingSequenceSignificance, Is.EqualTo(1.8003252080589478).Within(1e-11));
        });
    }

    // RSAT: -l 4 -1str -markov 1 -return occ,mseq,proba -lth mseq 3: tagc exp_freq 0.014423628634101493
    //   exp_ms 0.45182668021626959 ms_P 0.0052765589179613543 ms_E 1.3507990829981067 ms_sig -0.13059075713882745.
    // RSAT: -bg equi -lth mseq 3: tacg exp_ms 0.12709569604698778 ms_P 0.00012525578947139476 ms_sig 1.493962226078039.
    [Test]
    public void FindSharedMotifs_MarkovAndEquiprobable_EqualRsatMseq()
    {
        var seqs = Ms.Select(s => new DnaSequence(s)).ToList();
        var markov = Get(MotifFinder.FindSharedMotifs(seqs, 4, 3, OligoBackgroundModel.MarkovFromInput(1)), "TAGC");
        var equi = Get(MotifFinder.FindSharedMotifs(seqs, 4, 3, OligoBackgroundModel.Equiprobable), "TACG");

        Assert.Multiple(() =>
        {
            Assert.That(markov.ExpectedFrequency, Is.EqualTo(0.014423628634101493).Within(Rel).Percent);
            Assert.That(markov.ExpectedMatchingSequences, Is.EqualTo(0.45182668021626959).Within(Rel).Percent);
            Assert.That(markov.MatchingSequenceProbability, Is.EqualTo(0.0052765589179613543).Within(Rel).Percent);
            Assert.That(markov.MatchingSequenceEValue, Is.EqualTo(1.3507990829981067).Within(Rel).Percent);
            Assert.That(markov.MatchingSequenceSignificance, Is.EqualTo(-0.13059075713882745).Within(1e-11));
            Assert.That(equi.ExpectedMatchingSequences, Is.EqualTo(0.12709569604698778).Within(Rel).Percent);
            Assert.That(equi.MatchingSequenceProbability, Is.EqualTo(0.00012525578947139476).Within(Rel).Percent);
            Assert.That(equi.MatchingSequenceSignificance, Is.EqualTo(1.493962226078039).Within(1e-11));
        });
    }

    [Test]
    public void FindSharedMotifs_Significance_InvalidArguments()
    {
        var seqs = Ms.Select(s => new DnaSequence(s)).ToList();
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => MotifFinder.FindSharedMotifs(null!, 4, 2, OligoBackgroundModel.Equiprobable));
            Assert.Throws<ArgumentNullException>(() => MotifFinder.FindSharedMotifs(seqs, 4, 2, (OligoBackgroundModel)null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.FindSharedMotifs(seqs, 0, 2, OligoBackgroundModel.Equiprobable));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.FindSharedMotifs(seqs, 4, 0, OligoBackgroundModel.Equiprobable));
            Assert.Throws<ArgumentException>(() => MotifFinder.FindSharedMotifs(new[] { seqs[0], null! }, 4, 1, OligoBackgroundModel.Equiprobable));
        });
    }

    [Test]
    public void FindSharedMotifs_Significance_EmptyInput_NoMotifs()
    {
        var r = MotifFinder.FindSharedMotifs(Array.Empty<DnaSequence>(), 4, 1, OligoBackgroundModel.BernoulliFromInput);

        Assert.Multiple(() =>
        {
            Assert.That(r.Motifs, Is.Empty);
            Assert.That(r.SequenceCount, Is.Zero);
        });
    }

    #endregion
}
