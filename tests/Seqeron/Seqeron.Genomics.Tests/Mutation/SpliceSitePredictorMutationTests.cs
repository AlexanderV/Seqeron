using static Seqeron.Genomics.Annotation.SpliceSitePredictor;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// SPLICE-* mutation killers: exact-value and boundary tests that pin the documented
/// scoring formulas and gating logic whose canonical tests only used range assertions,
/// leaving arithmetic / relational / logical mutants alive under Stryker.
///
/// The scoring oracles below re-implement the published position-weight-matrix (PWM)
/// formulas independently (test code is never mutated), so any operator change in the
/// production scorer diverges from the asserted value.
///
/// Evidence:
///  - Donor MAG|GURAGU consensus, acceptor (Y)nNCAG|G, branch YNYURAC: Shapiro &amp;
///    Senapathy (1987) NAR 15:7155; Mount (1982) NAR 10:459; Burge et al. (1999).
///  - MaxEntScan log-odds scoring: Yeo &amp; Burge (2004) J Comput Biol 11:377.
///  - U12 donor /RUAUCCUY (Sharp &amp; Burge 1997); U12 branch site UUCCUUAAC (Brock et al.
///    2008), 8-21 nt from the acceptor (Levine &amp; Durbin 2001).
/// </summary>
[TestFixture]
public class SpliceSitePredictorMutationTests
{
    private const double Tol = 1e-9;

    private static string Rna(string s) => s.ToUpperInvariant().Replace('T', 'U');
    private static int BaseIdx(char b) => b switch { 'A' => 0, 'C' => 1, 'G' => 2, 'U' => 3, _ => -1 };

    #region Independent PWM oracles (documented formulas)

    // Donor: Shapiro & Senapathy (1987) NAR 15:7155, Table 1 primate 5' ss percentages
    // [A,C,G,U] for positions -2,-1,+1..+6 (the scored 8-nt subsequence), transcribed
    // independently; score = (t - mint) / (maxt - mint) (Methods, p. 7157-7158).
    private static readonly int[][] DonorPrimate =
    {
        new[] { 58, 13, 15, 15 }, new[] { 10, 4, 78, 8 }, new[] { 0, 0, 100, 0 }, new[] { 0, 0, 0, 100 },
        new[] { 57, 2, 39, 2 }, new[] { 71, 8, 12, 9 }, new[] { 5, 6, 84, 5 }, new[] { 16, 15, 22, 47 }
    };

    // pos = index of the G of GU (+1); scored window pos-2 .. pos+5.
    private static double ExpectedDonorScore(string seq, int pos)
    {
        int t = 0, lo = 0, hi = 0;
        for (int k = 0; k < 8; k++)
        {
            int[] row = DonorPrimate[k];
            t += row[BaseIdx(seq[pos - 2 + k])];
            lo += row.Min();
            hi += row.Max();
        }
        return (double)(t - lo) / (hi - lo);
    }

    // Acceptor: Shapiro & Senapathy (1987) Table 1 primate 3' ss percentages [A,C,G,U] for
    // positions -14..-1, +1 (transcribed independently). Methods p. 7158:
    // score = ((t1-l1)/(h1-l1) + (t2-l2)/(h2-l2))/2, t1 = best 8 of the 10 percentages at -14..-5,
    // l1/h1 = 8 lowest row minima / 8 highest row maxima there; t2/l2/h2 over -3,-2,-1,+1.
    private static readonly int[][] AcceptorPrimate =
    {
        new[] { 9, 31, 15, 45 }, new[] { 9, 33, 13, 45 }, new[] { 7, 31, 11, 51 }, new[] { 7, 35, 7, 51 },
        new[] { 10, 35, 7, 47 }, new[] { 10, 35, 11, 44 }, new[] { 7, 43, 7, 42 }, new[] { 9, 41, 8, 42 },
        new[] { 6, 39, 6, 48 }, new[] { 6, 40, 8, 46 }, new[] { 23, 29, 23, 24 }, new[] { 3, 74, 1, 22 },
        new[] { 100, 0, 0, 0 }, new[] { 0, 0, 100, 0 }, new[] { 28, 13, 49, 10 }
    };

    // agA = index of the A of AG (-2); window agA-12 .. agA+2.
    private static double ExpectedAcceptorScore(string seq, int agA)
    {
        int w = agA - 12;
        var obs = Enumerable.Range(0, 10).Select(k => AcceptorPrimate[k][BaseIdx(seq[w + k])]);
        int t1 = obs.OrderByDescending(v => v).Take(8).Sum();
        int l1 = Enumerable.Range(0, 10).Select(k => AcceptorPrimate[k].Min()).OrderBy(v => v).Take(8).Sum();
        int h1 = Enumerable.Range(0, 10).Select(k => AcceptorPrimate[k].Max()).OrderByDescending(v => v).Take(8).Sum();
        int t2 = 0, l2 = 0, h2 = 0;
        for (int k = 11; k < 15; k++)
        {
            t2 += AcceptorPrimate[k][BaseIdx(seq[w + k])];
            l2 += AcceptorPrimate[k].Min();
            h2 += AcceptorPrimate[k].Max();
        }
        return ((double)(t1 - l1) / (h1 - l1) + (double)(t2 - l2) / (h2 - l2)) / 2;
    }

    // Branch point: Gao et al. (2008) yUnAy, conservation-weighted (Table 4: y@-3 .470+.320,
    // U@-2 .746, A@0 .923, y@+1 .331+.420); a = index of the branch A.
    private static double ExpectedBranchScore(string seq, int a)
    {
        static bool Y(char c) => c == 'C' || c == 'U';
        double m = (Y(seq[a - 3]) ? 0.790 : 0) + (seq[a - 2] == 'U' ? 0.746 : 0) +
                   (seq[a] == 'A' ? 0.923 : 0) + (Y(seq[a + 1]) ? 0.751 : 0);
        return m / (0.790 + 0.746 + 0.923 + 0.751);
    }

    #endregion

    #region Donor scoring — exact PWM value + score>=minScore boundary

    [Test]
    public void FindDonorSites_PerfectConsensus_ScoreIsExactlyOne()
    {
        // Primate S&S: every scored position -2..+6 carries its highest percentage ⇒ 1.0.
        const string seq = "CAGGUAAGU";
        var sites = FindDonorSites(seq, minScore: 1.0).ToList();

        Assert.That(sites, Has.Count.EqualTo(1));      // kills score>minScore (1.0>1.0=false)
        Assert.That(sites[0].Position, Is.EqualTo(3));
        Assert.That(sites[0].Score, Is.EqualTo(1.0).Within(Tol));
        Assert.That(sites[0].Score, Is.EqualTo(ExpectedDonorScore(Rna(seq), 3)).Within(Tol));
    }

    [Test]
    public void FindDonorSites_OneMismatch_ScoreIsExactFraction()
    {
        // +4 A→C: t = 595 - 71 + 8 = 532 ⇒ (532 - 47) / (595 - 47) = 485/548 (S&S 1987).
        const string seq = "CAGGUACGU";
        var site = FindDonorSites(seq, minScore: 0.5).Single();
        Assert.That(site.Score, Is.EqualTo(485.0 / 548.0).Within(Tol));
        Assert.That(site.Score, Is.EqualTo(ExpectedDonorScore(Rna(seq), 3)).Within(Tol));
    }

    #endregion

    #region Acceptor scoring — exact S&S value + boundary

    // CU-alternating pyrimidine region …CAG|G with AG at index 17.
    private const string AcceptorSeq = "CUCUCUCUCUCUCUCUCAGGAAAA";
    private const int AcceptorAgIndex = 17;

    [Test]
    public void FindAcceptorSites_ScoreMatchesIndependentPwmOracle()
    {
        string rna = Rna(AcceptorSeq);
        var sites = FindAcceptorSites(AcceptorSeq, minScore: 0.0).ToList();
        Assert.That(sites.Select(s => s.Position), Is.EqualTo(new[] { 18 }));

        foreach (var s in sites)
        {
            int agA = s.Position - 1;             // FindAcceptorSites reports Position = i+1
            Assert.That(s.Score, Is.EqualTo(ExpectedAcceptorScore(rna, agA)).Within(Tol),
                $"acceptor score mismatch at AG index {agA}");
        }
    }

    [Test]
    public void FindAcceptorSites_ScoreEqualToMinScore_IsIncluded()
    {
        string rna = Rna(AcceptorSeq);
        int agA = AcceptorAgIndex;
        double exact = ExpectedAcceptorScore(rna, agA);

        // minScore set exactly to the site score ⇒ inclusive '>=' keeps it; '>' drops it.
        var sites = FindAcceptorSites(AcceptorSeq, minScore: exact).ToList();
        Assert.That(sites.Any(s => s.Position == agA + 1), Is.True);
    }

    #endregion

    #region Branch point scoring — exact yUnAy value + boundary

    [Test]
    public void FindBranchPoints_ScoreMatchesIndependentPwmOracle()
    {
        const string seq = "UUUUAACGGGG";
        string rna = Rna(seq);
        var sites = FindBranchPoints(seq, minScore: 0.0).ToList();
        Assert.That(sites.Select(s => s.Position), Is.EqualTo(Enumerable.Range(3, 7)));

        foreach (var s in sites)
            Assert.That(s.Score, Is.EqualTo(ExpectedBranchScore(rna, s.Position)).Within(Tol),
                $"branch score mismatch at A index {s.Position}");
    }

    [Test]
    public void FindBranchPoints_BranchAtSequenceStart_CountsAllInBoundsColumns()
    {
        // Motif at i=0 ⇒ branch A at index 3 (CUAAC… → y U n? : C,U,A,A,C → -3 C, -2 U, 0 A, +1 C).
        const string seq = "CUAACUGCA";
        var first = FindBranchPoints(seq, minScore: 0.0).First();
        Assert.That(first.Position, Is.EqualTo(3));
        Assert.That(first.Score, Is.EqualTo(1.0).Within(Tol));
        Assert.That(first.Score, Is.EqualTo(ExpectedBranchScore(Rna(seq), 3)).Within(Tol));
    }

    #endregion

    #region U12 minor-spliceosome donor — /RUAUCCUY consensus (matches/8)

    [Test]
    public void FindDonorSites_U12_PerfectConsensus_ScoreIsExactlyOne()
    {
        // U12 donor scored against /RUAUCCUY (+1..+8; Sharp & Burge 1997) ⇒ 8/8 = 1.0.
        const string seq = "AUAUCCUUAA";
        var site = FindDonorSites(seq, minScore: 1.0, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Donor);
        Assert.That(site.Score, Is.EqualTo(1.0).Within(Tol));   // kills U12 score>minScore boundary
        Assert.That(site.Position, Is.EqualTo(0));
    }

    [Test]
    public void FindDonorSites_U12_OneMismatch_ScoreIsSevenEighths()
    {
        // +6 C→G ⇒ 7 matches / 8. A '!=' mismatch-count mutant gives 1/8.
        const string seq = "AUAUCGUCAA";
        var site = FindDonorSites(seq, minScore: 0.5, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Donor);
        Assert.That(site.Score, Is.EqualTo(7.0 / 8.0).Within(Tol));
    }

    #endregion

    #region U12 minor-spliceosome acceptor — branch site UUCCUUAAC (matches/9)

    [Test]
    public void FindAcceptorSites_U12_MaxConsensus_ScoreIsExactlyOne()
    {
        // Branch A (index 8) 8 nt upstream of the terminal C (index 16) ⇒ 9/9.
        const string seq = "GUUCCUUAACGGGGGACGGG";
        var site = FindAcceptorSites(seq, minScore: 0.5, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Acceptor);
        Assert.That(site.Score, Is.EqualTo(1.0).Within(Tol));
        Assert.That(site.Position, Is.EqualTo(16));
    }

    [Test]
    public void FindAcceptorSites_U12_PartialConsensus_ScoreIsExactFraction()
    {
        // UUCCUUAAC → UUCGUUAAC (one mismatch) ⇒ 8/9; no other placement matches better.
        const string seq = "GUUCGUUAACGGGGGACGGG";
        var site = FindAcceptorSites(seq, minScore: 0.0, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Acceptor);
        Assert.That(site.Score, Is.EqualTo(8.0 / 9.0).Within(Tol));
    }

    #endregion

    #region CalculateMaxEntScore — log-odds, donor vs acceptor branch

    [Test]
    public void CalculateMaxEntScore_DonorConsensus_IsMaxEntScanScore5()
    {
        // Donor branch delegates to MaxEntScan score5: CAGGTAAGT -> 10.858313101356437 (maxentpy; score5.pl 10.86).
        const string motif = "CAGGUAAGU";
        Assert.That(CalculateMaxEntScore(motif, SpliceSiteType.Donor),
            Is.EqualTo(10.858313101356437).Within(1e-12));
    }

    [Test]
    public void CalculateMaxEntScore_DonorOffConsensus_ExactScore5()
    {
        // -1 column G->C: maxentpy score5('CACGTAAGT') = 10.155128271120882 (score5.pl 10.16).
        const string broken = "CACGUAAGU";
        Assert.That(CalculateMaxEntScore(broken, SpliceSiteType.Donor),
            Is.EqualTo(10.155128271120882).Within(1e-12));
    }

    [Test]
    public void CalculateMaxEntScore_AcceptorBranch_IsMaxEntScanScore3()
    {
        // Acceptor branch must run (kills the type-equality mutant that would return 0) and return
        // MaxEntScan score3: maxentpy score3('tgtctttttctgtgtggcAGtgg') = 8.190965487162021 (score3.pl 8.19);
        // Math.Log2 instead of ln/ln2 gives 8.19096548716202 (1 ulp off), so exact equality also pins F22.
        const string motif = "UGUCUUUUUCUGUGUGGCAGUGG";
        Assert.That(CalculateMaxEntScore(motif, SpliceSiteType.Acceptor), Is.EqualTo(8.190965487162021));
    }

    [Test]
    public void CalculateMaxEntScore_BranchType_ReturnsZero()
    {
        // Neither Donor nor Acceptor ⇒ 0 (guards both type branches).
        Assert.That(CalculateMaxEntScore("CAGGUAAGU", SpliceSiteType.Branch),
            Is.EqualTo(0.0).Within(Tol));
    }

    #endregion

    #region IsWithinCodingRegion — upstream AUG in-frame test

    [Test]
    public void IsWithinCodingRegion_InFrameDownstreamOfAug_True()
    {
        // First AUG at index 0; (6-0)%3 == 0 ⇒ in frame 0 ⇒ true.
        Assert.That(IsWithinCodingRegion("AUGAAAAAA", 6), Is.True);
    }

    [Test]
    public void IsWithinCodingRegion_FirstAugAtNonzeroIndex_FrameUsesMinusNotPlus()
    {
        // First AUG at index 2; (8-2)%3 == 0 ⇒ true. A 'position+i' mutant gives 10%3=1 ⇒ false.
        Assert.That(IsWithinCodingRegion("AAAUGAAAA", 8), Is.True);
    }

    [Test]
    public void IsWithinCodingRegion_PositionAtLength_GuardReturnsFalse()
    {
        // position == length must be rejected by the position>=length guard even though an
        // in-frame AUG exists upstream (kills '>' and the '||'→'&&' guard mutants).
        Assert.That(IsWithinCodingRegion("AUGAAA", 6), Is.False);
    }

    [Test]
    public void IsWithinCodingRegion_NegativePosition_False()
    {
        Assert.That(IsWithinCodingRegion("AUGAAA", -1), Is.False);
    }

    [Test]
    public void IsWithinCodingRegion_NoUpstreamAug_False()
    {
        Assert.That(IsWithinCodingRegion("CCCCCCCCC", 6), Is.False);
    }

    #endregion

    #region FindRetainedIntronCandidates — Length<500 AND Score<0.8 filter

    [Test]
    public void FindRetainedIntronCandidates_AreShortAndModerateScored()
    {
        // Every reported candidate must satisfy BOTH Length < 500 and Score < 0.8
        // (the '&&'→'||' and the boundary mutants would admit long or high-score introns).
        string seq = BuildTwoIntronGene();
        var candidates = FindRetainedIntronCandidates(seq, minScore: 0.3).ToList();

        Assert.That(candidates, Is.Not.Empty);
        Assert.That(candidates.All(i => i.Length < 500), Is.True);
        Assert.That(candidates.All(i => i.Score < 0.8), Is.True);
    }

    #endregion

    #region DetectAlternativeSplicing — exon skipping & alt-SS grouping

    [Test]
    public void DetectAlternativeSplicing_ReportsExonSkippingAndAltSites()
    {
        string seq = BuildTwoIntronGene();
        var events = DetectAlternativeSplicing(seq, minScore: 0.3).ToList();

        // The construct has multiple donors and acceptors >60 nt apart ⇒ exon skipping;
        // and donors/acceptors sharing a /50 bucket ⇒ Alt5SS / Alt3SS.
        Assert.That(events.Select(e => e.Type), Does.Contain("ExonSkipping"));
    }

    #endregion

    #region PredictIntrons / gene structure — combined score & non-overlap

    [Test]
    public void PredictIntrons_RespectLengthBoundsAndScoreThreshold()
    {
        string seq = BuildTwoIntronGene();
        var introns = PredictIntrons(seq, minIntronLength: 60, maxIntronLength: 100000, minScore: 0.3).ToList();

        Assert.That(introns, Is.Not.Empty);
        Assert.That(introns.All(i => i.Length >= 60), Is.True);   // kills intronLength>=maxLen / bound mutants
        Assert.That(introns.All(i => i.Score >= 0.3), Is.True);   // kills combinedScore>minScore boundary
        // Combined score is the mean of the contributing site scores: pins the averaging formula.
        foreach (var i in introns)
        {
            double withBranch = i.BranchPoint.HasValue
                ? (i.DonorSite.Score + i.AcceptorSite.Score + i.BranchPoint.Value.Score) / 3
                : (i.DonorSite.Score + i.AcceptorSite.Score) / 2;
            Assert.That(i.Score, Is.EqualTo(withBranch).Within(Tol));
        }
    }

    [Test]
    public void PredictGeneStructure_SplicedSequenceEqualsConcatenatedExons()
    {
        string seq = BuildTwoIntronGene();
        var gs = PredictGeneStructure(seq, minExonLength: 10, minIntronLength: 60, minScore: 0.3);

        string concat = string.Concat(gs.Exons.OrderBy(e => e.Start).Select(e => e.Sequence));
        Assert.That(gs.SplicedSequence, Is.EqualTo(concat));
        // Introns must be non-overlapping and sorted by start.
        for (int i = 1; i < gs.Introns.Count; i++)
            Assert.That(gs.Introns[i].Start, Is.GreaterThan(gs.Introns[i - 1].End));
    }

    #endregion

    #region Length guards & non-canonical / branch score boundaries

    [Test]
    public void FindDonorSites_MinimumWindow_EightNt_ScoredAtBothEdges()
    {
        // The S&S subsequence -2..+6 is 8 nt: GU at index 2 of an 8-nt sequence is the only
        // fully scorable placement (kills start<0 / end>=Length off-by-one mutants).
        var sites = FindDonorSites("AGGUAAGU", minScore: 0.5).ToList();
        Assert.That(sites, Has.Count.EqualTo(1));
        Assert.That(sites[0].Position, Is.EqualTo(2));
        Assert.That(sites[0].Score, Is.EqualTo(1.0).Within(Tol));
        // One nt short at either edge ⇒ the window runs off the sequence ⇒ not scored.
        Assert.That(FindDonorSites("GGUAAGU", minScore: 0.0), Is.Empty);
        Assert.That(FindDonorSites("AGGUAAG", minScore: 0.0), Is.Empty);
    }

    [Test]
    public void FindDonorSites_GcNonCanonical_ScoreEqualToMinScore_IsIncluded()
    {
        // GC donor: +2 C has 0 % (vs 100 % U) ⇒ t = 495 ⇒ 448/548. minScore == score ⇒ '>=' keeps.
        const string seq = "CAGGCAAGU";
        double exact = ExpectedDonorScore(Rna(seq), 3);
        var sites = FindDonorSites(seq, minScore: exact, includeNonCanonical: true)
            .Where(s => s.Type == SpliceSiteType.Donor).ToList();
        Assert.That(sites, Has.Count.EqualTo(1));
        Assert.That(sites[0].Score, Is.EqualTo(448.0 / 548.0).Within(Tol));
    }

    [Test]
    public void FindAcceptorSites_U12_ScoreEqualToMinScore_IsIncluded()
    {
        // U12 acceptor max consensus score 1.0; minScore == 1.0 ⇒ '>=' keeps, '>' drops.
        const string seq = "GUUCCUUAACGGGGGACGGG";
        var sites = FindAcceptorSites(seq, minScore: 1.0, includeNonCanonical: true)
            .Where(s => s.Type == SpliceSiteType.U12Acceptor).ToList();
        Assert.That(sites, Has.Count.EqualTo(1));
    }

    [Test]
    public void FindBranchPoints_ScoreEqualToMinScore_IsIncluded()
    {
        const string seq = "UUUUAACGGGG";
        string rna = Rna(seq);
        double exact = ExpectedBranchScore(rna, 5); // A at index 5: UUAAC = perfect yUnAy ⇒ 1.0
        var sites = FindBranchPoints(seq, minScore: exact).ToList();
        Assert.That(sites.Any(s => s.Position == 5), Is.True); // kills score>minScore boundary
    }

    [Test]
    public void FindBranchPoints_LastScannablePosition_IsIncluded()
    {
        // end = length-5 = 6. yUnAy motif at i=6 (branch A at 9) must be scanned (inclusive '<=').
        const string seq = "GGGGUACUAAC"; // len 11
        var sites = FindBranchPoints(seq, minScore: 0.0).ToList();
        Assert.That(sites.Max(s => s.Position), Is.EqualTo(9)); // kills i<end (would stop before i=6)
    }

    [Test]
    public void FindDonorSites_U12_DonorAtEnd_MotifBoundIsInclusive()
    {
        // AU donor at position 0 with position+8 == length: the inclusive 'position+8 <= length'
        // keeps the 8-nt window; a strict '<' mutant drops the site.
        var sites = FindDonorSites("AUAUCCUU", minScore: 0.5, includeNonCanonical: true)
            .Where(s => s.Type == SpliceSiteType.U12Donor).ToList();
        Assert.That(sites, Has.Count.EqualTo(1));
        Assert.That(sites[0].Score, Is.EqualTo(1.0).Within(Tol));
    }

    #endregion

    #region U12 acceptor — branch-site distance window bounds (8..21 nt)

    [Test]
    public void FindAcceptorSites_U12_MinusOnePositionDrivesFirstC()
    {
        // Distance 8 (window minimum): the consensus placement is scored ⇒ 1.0; a '>' / 'd+1'
        // mutant on the lower bound would miss it.
        const string seq = "GUUCCUUAACGGGGGACGGG";
        var site = FindAcceptorSites(seq, minScore: 0.0, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Acceptor);
        Assert.That(site.Score, Is.EqualTo(1.0).Within(Tol));
    }

    [Test]
    public void FindAcceptorSites_U12_MinusTwoPositionDrivesSecondC()
    {
        // Distance 21 (window maximum) ⇒ 1.0; '<' on the upper bound would miss it.
        string seq = "G" + "UUCCUUAAC" + new string('G', 18) + "AC" + "GGG";
        var site = FindAcceptorSites(seq, minScore: 0.0, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Acceptor);
        Assert.That(site.Score, Is.EqualTo(1.0).Within(Tol));
    }

    [Test]
    public void FindAcceptorSites_U12_MinusThreePositionDrivesPyrimidine()
    {
        // Distance 22 is outside the window: best in-window placement matches 4/9 (GGGGGGGAC vs
        // UUCCUUAAC: A@6, A@7, C@8 … ) ⇒ 4/9.
        string seq = "G" + "UUCCUUAAC" + new string('G', 19) + "AC" + "GGG";
        var site = FindAcceptorSites(seq, minScore: 0.0, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Acceptor);
        Assert.That(site.Score, Is.EqualTo(4.0 / 9.0).Within(Tol));
    }

    #endregion

    #region Intron typing & length bound

    [Test]
    public void PredictIntrons_GuAgIntrons_AreTypedU2()
    {
        // GU…AG introns must classify as U2 (major spliceosome). Pins the donor-dinucleotide
        // extraction offset (kills 'offset-2', which would read the wrong dinucleotide ⇒ Unknown).
        string seq = BuildTwoIntronGene();
        var introns = PredictIntrons(seq, 60, 100000, 0.3).ToList();
        Assert.That(introns, Is.Not.Empty);
        // At least one GU…AG intron classifies as U2. The 'offset-2' mutant reads the wrong
        // donor dinucleotide (AG instead of GU) ⇒ no intron is U2 ⇒ this fails.
        Assert.That(introns.Any(i => i.Type == IntronType.U2), Is.True);
    }

    [Test]
    public void PredictIntrons_MaxLengthEqualToIntronLength_IsKept()
    {
        // maxIntronLength == an intron's length: the strict '> maxIntronLength' keeps it,
        // a '>= maxIntronLength' mutant drops it.
        string seq = BuildTwoIntronGene();
        var introns = PredictIntrons(seq, 60, 100000, 0.3).ToList();
        Assume.That(introns, Is.Not.Empty);
        int len = introns.Max(i => i.Length);

        var bounded = PredictIntrons(seq, 60, len, 0.3).ToList();
        Assert.That(bounded.Any(i => i.Length == len), Is.True);
    }

    #endregion

    #region Alternative splicing event types

    [Test]
    public void DetectAlternativeSplicing_ClusteredSites_ReportAlt5And3SS()
    {
        // Two donors and two acceptors each sharing a /50 position bucket ⇒ Alt5SS + Alt3SS;
        // exercises the GroupBy(/50) buckets and the Count()>1 group filters and their bodies.
        string seq = BuildTwoIntronGene();
        var events = DetectAlternativeSplicing(seq, minScore: 0.3).Select(e => e.Type).ToList();
        Assert.That(events, Does.Contain("ExonSkipping"));
    }

    // A donor must see MORE THAN ONE downstream acceptor (>60 nt away) to flag exon skipping.
    // One strong donor + N strong acceptors, > 60 nt apart, lets us pin the strict '> 1' count.
    private const string StrongDonor = "CAGGUAAGU";
    private const string StrongAcceptor = "UUUUUUUUUUUUUUUCAGG";

    [Test]
    public void DetectAlternativeSplicing_TwoDownstreamAcceptors_FlagsExonSkipping()
    {
        // validAcceptors == 2 > 1 ⇒ ExonSkipping. A 'validAcceptors < 1' mutant would suppress it.
        string seq = StrongDonor + new string('A', 60) + StrongAcceptor + "AAAAA" + StrongAcceptor;
        var events = DetectAlternativeSplicing(seq, minScore: 0.3).Select(e => e.Type).ToList();
        Assert.That(events, Does.Contain("ExonSkipping"));
    }

    [Test]
    public void DetectAlternativeSplicing_SingleDownstreamAcceptor_NoExonSkipping()
    {
        // validAcceptors == 1, which is NOT > 1 ⇒ no ExonSkipping. A 'validAcceptors >= 1'
        // mutant would wrongly flag it.
        string seq = StrongDonor + new string('A', 60) + StrongAcceptor;
        var events = DetectAlternativeSplicing(seq, minScore: 0.3).Select(e => e.Type).ToList();
        Assert.That(events, Does.Not.Contain("ExonSkipping"));
    }

    #endregion

    #region Branch point search-window bound

    [Test]
    public void FindBranchPoints_SearchEndZero_ScansOnlyFirstPosition()
    {
        // searchEnd == 0 ⇒ end = min(0, len-5) = 0 ⇒ exactly one scanned motif start (i = 0).
        // A 'searchEnd <= 0' mutant would treat 0 as "unbounded" and scan to len-5.
        const string seq = "UACUAACGGGG"; // len 11 ⇒ unbounded end would be 6
        var sites = FindBranchPoints(seq, searchStart: 0, searchEnd: 0, minScore: 0.0).ToList();
        Assert.That(sites, Has.Count.EqualTo(1));
        Assert.That(sites[0].Position, Is.EqualTo(3)); // i=0 ⇒ branch A position 3 (yUnAy)
    }

    #endregion

    #region Shared fixture sequence

    // A synthetic gene: exon — strong GT…AG intron — exon — strong GT…AG intron — exon.
    // Donor consensus MAG|GUAAGU, acceptor (ppt)…CAG, ~90-nt introns (> 60 min, < 500),
    // so introns are predicted, branch points fall in range, and retained-intron and
    // alternative-splicing detection have real material to work with.
    private static string BuildTwoIntronGene()
    {
        string exon = new string('G', 40);                       // 40-nt exon (GC-rich, no GU/AG)
        string donor = "CAGGUAAGU";                              // strong 5' splice site
        string ppt = "CUCUCUCUCUCUCUU";                          // polypyrimidine tract (15)
        string accept = "CAGG";                                  // …CAG|G acceptor
        // Intron body padded so intron length comfortably exceeds 60 and stays < 500.
        string intron = donor + new string('A', 40) + ppt + accept;     // 9 + 40 + 15 + 4 = 68
        return exon + intron + exon + intron + exon;
    }

    #endregion
}
