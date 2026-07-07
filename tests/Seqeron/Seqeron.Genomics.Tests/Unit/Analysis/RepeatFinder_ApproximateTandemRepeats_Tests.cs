// REP-STR-001 — Approximate (Imperfect/Interrupted) Tandem-Repeat Detection (TRF model, opt-in)
// Evidence: docs/Evidence/REP-STR-001-Evidence.md
// TestSpec: tests/TestSpecs/REP-STR-001.md (§9)
// Source: Benson G (1999). Tandem repeats finder: a program to analyze DNA sequences.
//         Nucleic Acids Research 27(2):573-580. https://doi.org/10.1093/nar/27.2.573
//
// All expected statistics are hand-derived from the TRF alignment model with the recommended scoring
// (match +2, mismatch -7, indel -7) and read off the alignment columns. See TestSpec §9.3.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RepeatFinder_ApproximateTandemRepeats_Tests
{
    // Permissive minScore so the short worked-example tracts (scores 15/20/27) are reported; the
    // default-threshold behaviour (50) is asserted separately in A5/A6.
    private const int LowMinScore = 10;
    private const double Tol = 1e-9;

    #region FindApproximateTandemRepeats — perfect-alignment control (A1)

    // A1 — Perfect dinucleotide CA x5. Aligned against consensus "CA" tiled to 10 bp: 10 match columns,
    // 0 mismatch, 0 indel. score = 10*2 = 20; %matches = 100; %indels = 0; copies = 10/2 = 5.
    [Test]
    public void FindApproximateTandemRepeats_PerfectDinucleotide_Reports100PercentMatchesZeroIndels()
    {
        var results = RepeatFinder.FindApproximateTandemRepeats("CACACACACA", 1, 6, LowMinScore).ToList();

        var top = results.OrderByDescending(r => r.AlignmentScore).First();
        Assert.Multiple(() =>
        {
            Assert.That(top.Period, Is.EqualTo(2), "CA is a period-2 motif");
            Assert.That(top.Consensus, Is.EqualTo("CA"), "majority-rule consensus of a perfect CA tract");
            Assert.That(top.ConsensusSize, Is.EqualTo(2), "consensus size equals the period in this subset");
            Assert.That(top.CopyNumber, Is.EqualTo(5.0).Within(Tol), "10 aligned bases / period 2 = 5 copies");
            Assert.That(top.PercentMatches, Is.EqualTo(100.0).Within(Tol), "a perfect tract is 100% matches");
            Assert.That(top.PercentIndels, Is.EqualTo(0.0).Within(Tol), "a perfect tract has 0% indels");
            Assert.That(top.AlignmentScore, Is.EqualTo(20), "10 matches * +2 = 20 (TRF scoring)");
        });
    }

    #endregion

    #region FindApproximateTandemRepeats — one substitution (A2, A2b, A3)

    // A2 — (CAG)x6 with copy-4 first base substituted C->T at index 9: "CAGCAGCAGTAGCAGCAG".
    // Aligned against "CAG" tiled to 18 bp: 17 match, 1 mismatch, 0 indel.
    // score = 17*2 + 1*(-7) = 27; %matches = 17/18*100 = 94.444...; %indels = 0; copies = 18/3 = 6.
    [Test]
    public void FindApproximateTandemRepeats_OneSubstitutionTrinucleotide_ReportsExactPercentMatches()
    {
        var results = RepeatFinder.FindApproximateTandemRepeats("CAGCAGCAGTAGCAGCAG", 3, 3, LowMinScore).ToList();

        var top = results.OrderByDescending(r => r.AlignmentScore).First();
        Assert.Multiple(() =>
        {
            Assert.That(top.Period, Is.EqualTo(3));
            Assert.That(top.Consensus, Is.EqualTo("CAG"), "majority of 6 copies is CAG despite the single TAG");
            Assert.That(top.CopyNumber, Is.EqualTo(6.0).Within(Tol), "18 bases / period 3 = 6 copies");
            Assert.That(top.PercentMatches, Is.EqualTo(17.0 / 18.0 * 100.0).Within(Tol), "17 of 18 columns match");
            Assert.That(top.PercentIndels, Is.EqualTo(0.0).Within(Tol), "a pure substitution has no indels");
            Assert.That(top.AlignmentScore, Is.EqualTo(27), "17*2 - 7 = 27");
        });
    }

    // A2b — the SAME interrupted tract under the PERFECT detector fragments into one short CAG x3 run
    // (it breaks at the substitution and the trailing CAG x2 is below minRepeats=3). This is the gap
    // the approximate detector closes: the interrupted locus is reported as ONE repeat above, but the
    // perfect detector does not span the interruption.
    [Test]
    public void FindMicrosatellites_OnInterruptedTract_FragmentsWherePerfectDetectorBreaks()
    {
        var perfect = RepeatFinder.FindMicrosatellites("CAGCAGCAGTAGCAGCAG", 1, 6, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(perfect, Has.Count.EqualTo(1), "perfect detector reports only the leading run");
            Assert.That(perfect[0].RepeatUnit, Is.EqualTo("CAG"));
            Assert.That(perfect[0].RepeatCount, Is.EqualTo(3), "CAG x3 before the TAG interruption");
            Assert.That(perfect[0].Position, Is.EqualTo(0));
        });
    }

    // A3 — (CA)x6 with index 6 substituted to T: "CACACATACACA".
    // Aligned against "CA" tiled to 12 bp: 11 match, 1 mismatch, 0 indel.
    // score = 11*2 - 7 = 15; %matches = 11/12*100 = 91.666...; %indels = 0; copies = 12/2 = 6.
    [Test]
    public void FindApproximateTandemRepeats_OneSubstitutionDinucleotide_ReportsExactPercentMatches()
    {
        var results = RepeatFinder.FindApproximateTandemRepeats("CACACATACACA", 2, 2, LowMinScore).ToList();

        var top = results.OrderByDescending(r => r.AlignmentScore).First();
        Assert.Multiple(() =>
        {
            Assert.That(top.Period, Is.EqualTo(2));
            Assert.That(top.Consensus, Is.EqualTo("CA"));
            Assert.That(top.CopyNumber, Is.EqualTo(6.0).Within(Tol), "12 bases / period 2 = 6 copies");
            Assert.That(top.PercentMatches, Is.EqualTo(11.0 / 12.0 * 100.0).Within(Tol), "11 of 12 columns match");
            Assert.That(top.PercentIndels, Is.EqualTo(0.0).Within(Tol));
            Assert.That(top.AlignmentScore, Is.EqualTo(15), "11*2 - 7 = 15");
        });
    }

    #endregion

    #region FindApproximateTandemRepeats — one deletion / indel (A4)

    // A4 — perfect (CAG)x10 (30 bp) with one base deleted at index 15 ->
    // "CAGCAGCAGCAGCAGAGCAGCAGCAGCAG" (29 bp). Aligned against "CAG" tiled to 30 bp (10 whole copies):
    // 29 match, 0 mismatch, 1 indel (one gap to absorb the deletion); 30 columns.
    // score = 29*2 + 1*(-7) = 51; %matches = 29/30*100 = 96.666...; %indels = 1/30*100 = 3.333...;
    // copies = 29/3 = 9.666...
    [Test]
    public void FindApproximateTandemRepeats_OneDeletion_ReportsExactPercentIndels()
    {
        const string seq = "CAGCAGCAGCAGCAGAGCAGCAGCAGCAG"; // (CAG)x10 with one base deleted, 29 bp
        Assert.That(seq, Has.Length.EqualTo(29), "test fixture: 29-bp deletion tract");

        var results = RepeatFinder.FindApproximateTandemRepeats(seq, 3, 3, LowMinScore).ToList();

        var top = results.OrderByDescending(r => r.AlignmentScore).First();
        Assert.Multiple(() =>
        {
            Assert.That(top.Period, Is.EqualTo(3));
            Assert.That(top.Consensus, Is.EqualTo("CAG"));
            Assert.That(top.CopyNumber, Is.EqualTo(29.0 / 3.0).Within(Tol), "29 aligned bases / period 3");
            Assert.That(top.PercentMatches, Is.EqualTo(29.0 / 30.0 * 100.0).Within(Tol), "29 of 30 columns match");
            Assert.That(top.PercentIndels, Is.EqualTo(1.0 / 30.0 * 100.0).Within(Tol), "1 of 30 columns is an indel");
            Assert.That(top.AlignmentScore, Is.EqualTo(51), "29*2 - 7 = 51");
        });
    }

    #endregion

    #region Minscore threshold (A5, A6)

    // A5 — Benson (1999): "Only those repeats scoring at least 50 ... are reported." A perfect CA x5
    // tract scores only 20, so at the default minScore (50) it is NOT reported.
    [Test]
    public void FindApproximateTandemRepeats_DefaultMinScore_SuppressesLowScoringTract()
    {
        var results = RepeatFinder.FindApproximateTandemRepeats("CACACACACA", 1, 6).ToList(); // default minScore = 50

        Assert.That(results, Is.Empty, "CA x5 scores 20 < 50, below the Benson (1999) report threshold");
    }

    // A6 — the deletion tract scores 51 >= 50, so it IS reported at the default threshold.
    [Test]
    public void FindApproximateTandemRepeats_DefaultMinScore_ReportsTractAtOrAboveThreshold()
    {
        var results = RepeatFinder.FindApproximateTandemRepeats(
            "CAGCAGCAGCAGCAGAGCAGCAGCAGCAG", 3, 3).ToList(); // default minScore = 50

        Assert.That(results, Is.Not.Empty, "score 51 >= 50 is reportable");
        Assert.That(results.Max(r => r.AlignmentScore), Is.EqualTo(51));
        Assert.That(RepeatFinder.DefaultApproximateMinScore, Is.EqualTo(50),
            "default minimum score is the Benson (1999) recommended Minscore = 50");
    }

    #endregion

    #region Edge cases and validation (A7-A11)

    // A7 — empty sequence.
    [Test]
    public void FindApproximateTandemRepeats_EmptySequence_ReturnsEmpty()
    {
        Assert.That(RepeatFinder.FindApproximateTandemRepeats("", 1, 6, LowMinScore), Is.Empty);
    }

    // A8 — a sequence with no tandem structure yields nothing.
    [Test]
    public void FindApproximateTandemRepeats_NoRepeat_ReturnsEmpty()
    {
        Assert.That(RepeatFinder.FindApproximateTandemRepeats("ACGTGCAT", 1, 6, LowMinScore), Is.Empty);
    }

    // A9 — invalid period parameters throw.
    [Test]
    public void FindApproximateTandemRepeats_InvalidPeriodParameters_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => RepeatFinder.FindApproximateTandemRepeats("CACACACACA", 0, 6, LowMinScore).ToList(),
                "minPeriod < 1 is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => RepeatFinder.FindApproximateTandemRepeats("CACACACACA", 4, 2, LowMinScore).ToList(),
                "maxPeriod < minPeriod is invalid");
        });
    }

    // A10 — null DnaSequence throws.
    [Test]
    public void FindApproximateTandemRepeats_NullDnaSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => RepeatFinder.FindApproximateTandemRepeats((DnaSequence)null!, 1, 6, LowMinScore).ToList());
    }

    // A11 — deterministic: identical input yields identical results.
    [Test]
    public void FindApproximateTandemRepeats_RunTwice_ProducesIdenticalResults()
    {
        const string seq = "CAGCAGCAGCAGCAGAGCAGCAGCAGCAG";
        var first = RepeatFinder.FindApproximateTandemRepeats(seq, 1, 6, LowMinScore).ToList();
        var second = RepeatFinder.FindApproximateTandemRepeats(seq, 1, 6, LowMinScore).ToList();

        Assert.That(second, Is.EqualTo(first), "approximate detection is deterministic");
    }

    // DnaSequence overload agrees with the string overload (delegation smoke).
    [Test]
    public void FindApproximateTandemRepeats_DnaSequenceOverload_MatchesStringOverload()
    {
        const string seq = "CAGCAGCAGCAGCAGAGCAGCAGCAGCAG";
        var viaString = RepeatFinder.FindApproximateTandemRepeats(seq, 3, 3, LowMinScore).ToList();
        var viaDna = RepeatFinder.FindApproximateTandemRepeats(new DnaSequence(seq), 3, 3, LowMinScore).ToList();

        Assert.That(viaDna, Is.EqualTo(viaString), "DnaSequence and string overloads agree");
    }

    // A12 — single-character sequence cannot hold two contiguous copies of any period (Benson 1999:
    // "two or more contiguous copies"), so no tandem repeat is reportable.
    [Test]
    public void FindApproximateTandemRepeats_SingleCharacter_ReturnsEmpty()
    {
        Assert.That(RepeatFinder.FindApproximateTandemRepeats("A", 1, 6, LowMinScore), Is.Empty,
            "one base cannot span two copies of any period");
    }

    // A13 — a period longer than the sequence can never fit two copies (needs period*2 <= length), so
    // nothing is reported even at the permissive threshold.
    [Test]
    public void FindApproximateTandemRepeats_PeriodLongerThanSequence_ReturnsEmpty()
    {
        Assert.That(RepeatFinder.FindApproximateTandemRepeats("ACG", 6, 6, LowMinScore), Is.Empty,
            "period 6 cannot tile twice into a 3-bp sequence");
    }

    // A14 — all-N homopolymer. The detector treats 'N' as a literal base (no IUPAC ambiguity
    // special-casing): "NNNNNNNN" is a perfect period-1 tandem run of consensus "N", aligned against
    // "NNNNNNNN" -> 8 match columns, 0 mismatch, 0 indel. score = 8*2 = 16; %matches = 100; %indels = 0;
    // copies = 8/1 = 8.
    [Test]
    public void FindApproximateTandemRepeats_AllN_TreatedAsLiteralPerfectHomopolymer()
    {
        var results = RepeatFinder.FindApproximateTandemRepeats("NNNNNNNN", 1, 6, LowMinScore).ToList();

        var top = results.OrderByDescending(r => r.AlignmentScore).First();
        Assert.Multiple(() =>
        {
            Assert.That(top.Period, Is.EqualTo(1), "N homopolymer is a period-1 motif");
            Assert.That(top.Consensus, Is.EqualTo("N"), "majority-rule consensus of all-N is N (literal)");
            Assert.That(top.PercentMatches, Is.EqualTo(100.0).Within(Tol), "literal N matches N in every column");
            Assert.That(top.PercentIndels, Is.EqualTo(0.0).Within(Tol));
            Assert.That(top.AlignmentScore, Is.EqualTo(16), "8 match columns * +2 = 16");
            Assert.That(top.CopyNumber, Is.EqualTo(8.0).Within(Tol), "8 bases / period 1 = 8 copies");
        });
    }

    #endregion

    #region ComputeBernoulliStatistics — TRF Bernoulli model (B1-B10)
    //
    // Benson G (1999), Nucleic Acids Res 27(2):573-580, https://doi.org/10.1093/nar/27.2.573;
    // TRF desc/definitions (tandem.bu.edu/trf/trf.desc.html, trf.definitions.html; Benson-Genomics-Lab/TRF):
    //   "We model alignment of two tandem copies of a pattern ... by ... independent Bernoulli trials";
    //   "P(Heads), which we also call PM or matching probability, represents the average percent identity
    //   between the copies"; "PI or indel probability ... the average percentage of insertions and
    //   deletions between the copies"; statistics are "between adjacent copies ... not between the
    //   sequence and the consensus pattern"; defaults "PM = .80 and PI = .10".
    // Each Bernoulli trial = one alignment column between two ADJACENT copies (heads = match).
    // All expected PM/PI/E[matches] below are hand-derived per copy-pair; see comments.

    // B1 — Perfect CA x5: copies CA|CA|CA|CA|CA -> 4 adjacent pairs, 8 columns, all match.
    // PM = 8/8 = 1.0; PI = 0; E[matches] = PM*d = 1.0*8 = 8.
    [Test]
    public void ComputeBernoulliStatistics_PerfectTract_MatchProbabilityIsOne()
    {
        var stats = RepeatFinder.ComputeBernoulliStatistics("CACACACACA", period: 2);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Period, Is.EqualTo(2));
            Assert.That(stats.AdjacentCopyPairs, Is.EqualTo(4), "5 copies -> 4 adjacent pairs");
            Assert.That(stats.BernoulliTrials, Is.EqualTo(8), "4 pairs * 2 columns each");
            Assert.That(stats.Matches, Is.EqualTo(8), "every column matches in a perfect tract");
            Assert.That(stats.Mismatches, Is.EqualTo(0));
            Assert.That(stats.Indels, Is.EqualTo(0));
            Assert.That(stats.MatchProbability, Is.EqualTo(1.0).Within(Tol), "PM = 8/8 (perfect identity)");
            Assert.That(stats.IndelProbability, Is.EqualTo(0.0).Within(Tol), "PI = 0 (no indels)");
            Assert.That(stats.PercentMatches, Is.EqualTo(100.0).Within(Tol));
            Assert.That(stats.PercentIndels, Is.EqualTo(0.0).Within(Tol));
            Assert.That(stats.ExpectedMatches, Is.EqualTo(8.0).Within(Tol), "E[heads] = PM*d = 1.0*8");
            Assert.That(stats.MeetsExpectedMatchProbability, Is.True, "PM 1.0 >= default PM 0.80");
        });
    }

    // B2 — CAG x6 with copy 4 = TAG (one substitution): copies CAG|CAG|CAG|TAG|CAG|CAG.
    // Pairs: (CAG,CAG)3 (CAG,CAG)3 (CAG,TAG)2m+1mis (TAG,CAG)2m+1mis (CAG,CAG)3 ->
    // matches=13, mismatches=2, indels=0, trials=15. PM = 13/15; PI = 0; E[matches] = 13.
    [Test]
    public void ComputeBernoulliStatistics_OneSubstitution_EstimatesAdjacentCopyPm()
    {
        var stats = RepeatFinder.ComputeBernoulliStatistics("CAGCAGCAGTAGCAGCAG", period: 3);

        Assert.Multiple(() =>
        {
            Assert.That(stats.AdjacentCopyPairs, Is.EqualTo(5), "6 copies -> 5 adjacent pairs");
            Assert.That(stats.BernoulliTrials, Is.EqualTo(15), "5 pairs * 3 columns each");
            Assert.That(stats.Matches, Is.EqualTo(13), "the substituted base mismatches in the 2 pairs it joins");
            Assert.That(stats.Mismatches, Is.EqualTo(2));
            Assert.That(stats.Indels, Is.EqualTo(0));
            Assert.That(stats.MatchProbability, Is.EqualTo(13.0 / 15.0).Within(Tol), "PM = 13/15 between adjacent copies");
            Assert.That(stats.IndelProbability, Is.EqualTo(0.0).Within(Tol));
            Assert.That(stats.PercentMatches, Is.EqualTo(13.0 / 15.0 * 100.0).Within(Tol));
            Assert.That(stats.ExpectedMatches, Is.EqualTo(13.0).Within(Tol), "E[heads] = (13/15)*15 = 13");
        });
    }

    // B3 — CA x6 with index 6 C->T: copies CA|CA|CA|TA|CA|CA.
    // Pairs: 2 2 (CA,TA)=1m+1mis (TA,CA)=1m+1mis 2 -> matches=8, mismatches=2, trials=10.
    // PM = 8/10 = 0.80 (exactly the default threshold); PI = 0; E[matches] = 8.
    [Test]
    public void ComputeBernoulliStatistics_PmEqualToDefault_MeetsThreshold()
    {
        var stats = RepeatFinder.ComputeBernoulliStatistics("CACACATACACA", period: 2);

        Assert.Multiple(() =>
        {
            Assert.That(stats.BernoulliTrials, Is.EqualTo(10), "5 adjacent pairs * 2 columns");
            Assert.That(stats.Matches, Is.EqualTo(8));
            Assert.That(stats.Mismatches, Is.EqualTo(2));
            Assert.That(stats.MatchProbability, Is.EqualTo(0.80).Within(Tol), "PM = 8/10 = 0.80");
            Assert.That(stats.ExpectedMatches, Is.EqualTo(8.0).Within(Tol), "E[heads] = 0.80*10 = 8");
            Assert.That(stats.MeetsExpectedMatchProbability, Is.True, "PM 0.80 >= default PM 0.80 (>= is inclusive)");
        });
    }

    // B4 — a poorly conserved tract falls below the default PM = 0.80.
    // ACAC|TGTG: one adjacent pair of period 4, aligned A/T C/G A/T C/G -> 0 match, 4 mismatch.
    // PM = 0/4 = 0 < 0.80 -> MeetsExpectedMatchProbability = false.
    [Test]
    public void ComputeBernoulliStatistics_LowIdentityTract_DoesNotMeetThreshold()
    {
        var stats = RepeatFinder.ComputeBernoulliStatistics("ACACTGTG", period: 4);

        Assert.Multiple(() =>
        {
            Assert.That(stats.AdjacentCopyPairs, Is.EqualTo(1), "2 copies -> 1 adjacent pair");
            Assert.That(stats.MatchProbability, Is.EqualTo(0.0).Within(Tol), "no column matches");
            Assert.That(stats.MeetsExpectedMatchProbability, Is.False, "PM 0 < default PM 0.80");
        });
    }

    // B5 — caller-supplied PM threshold. The B2 tract has PM = 13/15 = 0.8667; it meets PM=0.80 but
    // not PM=0.90.
    [Test]
    public void ComputeBernoulliStatistics_CustomExpectedPm_FlagsAgainstThatThreshold()
    {
        var meets = RepeatFinder.ComputeBernoulliStatistics("CAGCAGCAGTAGCAGCAG", 3, expectedMatchProbability: 0.80);
        var fails = RepeatFinder.ComputeBernoulliStatistics("CAGCAGCAGTAGCAGCAG", 3, expectedMatchProbability: 0.90);

        Assert.Multiple(() =>
        {
            Assert.That(meets.MeetsExpectedMatchProbability, Is.True, "PM 0.8667 >= 0.80");
            Assert.That(fails.MeetsExpectedMatchProbability, Is.False, "PM 0.8667 < 0.90");
        });
    }

    // B6 — an indel tract reports PI > 0 (insertions/deletions between adjacent copies). Exact PI is
    // alignment-frame dependent, so we assert the qualitative Bernoulli property: with a single-base
    // deletion the indel probability is strictly positive and PM + PI account for the trials.
    [Test]
    public void ComputeBernoulliStatistics_DeletionTract_ReportsPositiveIndelProbability()
    {
        // (CAG)x10 with one base deleted (29 bp): adjacent copies must absorb the frame shift via indels.
        var stats = RepeatFinder.ComputeBernoulliStatistics("CAGCAGCAGCAGCAGAGCAGCAGCAGCAG", period: 3);

        Assert.Multiple(() =>
        {
            Assert.That(stats.IndelProbability, Is.GreaterThan(0.0), "a deletion forces indel columns between adjacent copies");
            Assert.That(stats.Matches + stats.Mismatches + stats.Indels, Is.EqualTo(stats.BernoulliTrials),
                "every Bernoulli trial is exactly one of match / mismatch / indel");
            Assert.That(stats.MatchProbability + stats.IndelProbability, Is.LessThanOrEqualTo(1.0 + Tol),
                "PM + PI <= 1 (mismatches make up the remainder)");
        });
    }

    // B7 — PM + (mismatch fraction) + PI = 1 exactly (the three Bernoulli outcomes partition the trials).
    [Test]
    public void ComputeBernoulliStatistics_Probabilities_PartitionTheTrials()
    {
        var stats = RepeatFinder.ComputeBernoulliStatistics("CAGCAGCAGTAGCAGCAG", period: 3);

        double mismatchFraction = (double)stats.Mismatches / stats.BernoulliTrials;
        Assert.That(stats.MatchProbability + mismatchFraction + stats.IndelProbability,
            Is.EqualTo(1.0).Within(Tol), "match + mismatch + indel fractions sum to 1");
    }

    // B8 — too-short tract (fewer than two copies) is rejected: a Bernoulli model of "two tandem copies"
    // is undefined.
    [Test]
    public void ComputeBernoulliStatistics_FewerThanTwoCopies_Throws()
    {
        Assert.Throws<ArgumentException>(() => RepeatFinder.ComputeBernoulliStatistics("CAG", period: 3),
            "one copy cannot model alignment of two adjacent copies");
    }

    // B9 — invalid parameters throw.
    [Test]
    public void ComputeBernoulliStatistics_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.ComputeBernoulliStatistics(null!, 2));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.ComputeBernoulliStatistics("CACA", 0));
            Assert.Throws<ArgumentOutOfRangeException>(
                () => RepeatFinder.ComputeBernoulliStatistics("CACA", 2, expectedMatchProbability: 1.5));
        });
    }

    // B10 — Benson (1999) documented defaults are exposed verbatim (PM = .80, PI = .10).
    [Test]
    public void ComputeBernoulliStatistics_DocumentedDefaults_MatchBenson1999()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.TrfDefaultMatchProbability, Is.EqualTo(0.80).Within(Tol),
                "Benson (1999) default PM = .80");
            Assert.That(RepeatFinder.TrfDefaultIndelProbability, Is.EqualTo(0.10).Within(Tol),
                "Benson (1999) default PI = .10");
        });
    }

    #endregion
}
