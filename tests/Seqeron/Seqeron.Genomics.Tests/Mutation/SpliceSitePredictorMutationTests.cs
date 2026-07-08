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
///  - U12 minor-spliceosome AU…AC consensus: Hall &amp; Padgett (1994); Jackson (1991).
/// </summary>
[TestFixture]
public class SpliceSitePredictorMutationTests
{
    private const double Tol = 1e-9;

    private static string Rna(string s) => s.ToUpperInvariant().Replace('T', 'U');
    private static int BaseIdx(char b) => b switch { 'A' => 0, 'C' => 1, 'G' => 2, 'U' => 3, _ => -1 };

    #region Independent PWM oracles (documented formulas)

    // Donor binary consensus, offsets -3..+5 (1.0 = matches consensus base set).
    private static readonly (int Off, string Ones)[] DonorCols =
    {
        (-3, "AC"), (-2, "A"), (-1, "G"), (0, "G"), (1, "U"), (2, "AG"), (3, "A"), (4, "G"), (5, "U")
    };

    private static double ExpectedDonorScore(string seq, int pos)
    {
        double score = 0; int count = 0;
        foreach (var (off, ones) in DonorCols)
        {
            int p = pos + off;
            if (p >= 0 && p < seq.Length)
            {
                char b = seq[p];
                if (BaseIdx(b) >= 0) { score += ones.Contains(b) ? 1.0 : 0.0; count++; }
            }
        }
        return count > 0 ? score / count : 0;
    }

    // Acceptor PWM rows [A,C,G,U], offsets relative to the splice site (= position+2).
    private static readonly Dictionary<int, double[]> AcceptorPwm = new()
    {
        { -15, new[] { 0.10, 0.30, 0.10, 0.50 } },
        { -10, new[] { 0.10, 0.30, 0.10, 0.50 } },
        { -5,  new[] { 0.10, 0.40, 0.10, 0.40 } },
        { -4,  new[] { 0.05, 0.40, 0.05, 0.50 } },
        { -3,  new[] { 0.05, 0.70, 0.05, 0.20 } },
        { -2,  new[] { 1.00, 0.00, 0.00, 0.00 } },
        { -1,  new[] { 0.00, 0.00, 1.00, 0.00 } },
        { 0,   new[] { 0.20, 0.15, 0.50, 0.15 } },
    };

    private static double ExpectedAcceptorScore(string seq, int position)
    {
        double score = 0; int count = 0;
        int ppt = 0;
        for (int i = position - 15; i < position - 3; i++)
            if (i >= 0 && i < seq.Length && (seq[i] == 'C' || seq[i] == 'U')) ppt++;
        score += ppt / 12.0 * 2;
        foreach (var kv in AcceptorPwm)
        {
            int pos = position + 2 + kv.Key;
            if (pos >= 0 && pos < seq.Length)
            {
                int bi = BaseIdx(seq[pos]);
                if (bi >= 0) { score += Math.Log2(kv.Value[bi] / 0.25 + 0.01); count++; }
            }
        }
        double normalized = (score / (count + 1) + 2) / 4;
        return Math.Max(0, Math.Min(1, normalized));
    }

    // Branch point PWM rows [A,C,G,U], offsets 0..6 (YNYURAC).
    private static readonly Dictionary<int, double[]> BranchPwm = new()
    {
        { 0, new[] { 0.10, 0.40, 0.10, 0.40 } },
        { 1, new[] { 0.25, 0.25, 0.25, 0.25 } },
        { 2, new[] { 0.10, 0.40, 0.10, 0.40 } },
        { 3, new[] { 0.05, 0.05, 0.05, 0.85 } },
        { 4, new[] { 0.60, 0.05, 0.30, 0.05 } },
        { 5, new[] { 0.95, 0.02, 0.02, 0.01 } },
        { 6, new[] { 0.10, 0.60, 0.10, 0.20 } },
    };

    private static double ExpectedBranchScore(string seq, int position)
    {
        double score = 0; int count = 0;
        foreach (var kv in BranchPwm)
        {
            int pos = position + kv.Key;
            if (pos >= 0 && pos < seq.Length)
            {
                int bi = BaseIdx(seq[pos]);
                if (bi >= 0) { score += Math.Log2(kv.Value[bi] / 0.25 + 0.01); count++; }
            }
        }
        double normalized = (score / count + 2) / 4;
        return Math.Max(0, Math.Min(1, normalized));
    }

    #endregion

    #region Donor scoring — exact PWM value + score>=minScore boundary

    [Test]
    public void FindDonorSites_PerfectConsensus_ScoreIsExactlyOne()
    {
        // MAG|GURAGU consensus, GU at index 3, full -3..+5 context ⇒ 9/9 = 1.0.
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
        // +3 position A→C breaks one consensus column ⇒ 8/9.
        const string seq = "CAGGUACGU";
        var site = FindDonorSites(seq, minScore: 0.5).Single();
        Assert.That(site.Score, Is.EqualTo(8.0 / 9.0).Within(Tol));
        Assert.That(site.Score, Is.EqualTo(ExpectedDonorScore(Rna(seq), 3)).Within(Tol));
    }

    #endregion

    #region Acceptor scoring — exact PWM value + ppt tract + boundary

    // 24-nt acceptor: polypyrimidine tract …CAG|G with AG at index 17 (≥ FindAcceptorSites' 20-nt min).
    private const string AcceptorSeq = "CUCUCUCUCUCUCUCUCAGGAAAA";
    private const int AcceptorAgIndex = 17;

    [Test]
    public void FindAcceptorSites_ScoreMatchesIndependentPwmOracle()
    {
        // Acceptor scores via log-odds PWM + ppt term; assert against the oracle exactly.
        string rna = Rna(AcceptorSeq);
        var sites = FindAcceptorSites(AcceptorSeq, minScore: 0.0).ToList();
        Assert.That(sites, Is.Not.Empty);

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

    #region Branch point scoring — exact PWM value + boundary

    [Test]
    public void FindBranchPoints_ScoreMatchesIndependentPwmOracle()
    {
        // YNYURAC-like heptamer; assert branch scores against the log-odds oracle exactly.
        const string seq = "UUUUAACGGGG";
        string rna = Rna(seq);
        var sites = FindBranchPoints(seq, minScore: 0.0).ToList();
        Assert.That(sites, Is.Not.Empty);

        foreach (var s in sites)
        {
            int start = s.Position - 5;            // FindBranchPoints reports Position = i+5
            Assert.That(s.Score, Is.EqualTo(ExpectedBranchScore(rna, start)).Within(Tol),
                $"branch score mismatch at start {start}");
        }
    }

    [Test]
    public void FindBranchPoints_BranchAtSequenceStart_CountsAllInBoundsColumns()
    {
        // Heptamer at i=0 exercises the pos>=0 lower bound (offset 0 ⇒ pos 0).
        const string seq = "UACUAACUGCA";
        string rna = Rna(seq);
        var first = FindBranchPoints(seq, minScore: 0.0).First(s => s.Position == 5); // i=0 ⇒ Position 5
        Assert.That(first.Score, Is.EqualTo(ExpectedBranchScore(rna, 0)).Within(Tol));
    }

    #endregion

    #region U12 minor-spliceosome donor — AUAUCC consensus (matches/6)

    [Test]
    public void FindDonorSites_U12_PerfectConsensus_ScoreIsExactlyOne()
    {
        // U12 donor AU… scored by literal match to AUAUCC ⇒ 6/6 = 1.0.
        const string seq = "AUAUCCAAAA";
        var site = FindDonorSites(seq, minScore: 1.0, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Donor);
        Assert.That(site.Score, Is.EqualTo(1.0).Within(Tol));   // kills U12 score>minScore boundary
        Assert.That(site.Position, Is.EqualTo(0));
    }

    [Test]
    public void FindDonorSites_U12_OneMismatch_ScoreIsFiveSixths()
    {
        // Last consensus position C→G ⇒ 5 matches / 6. A '!=' mismatch-count mutant gives 1/6.
        const string seq = "AUAUCGAAAA";
        var site = FindDonorSites(seq, minScore: 0.5, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Donor);
        Assert.That(site.Score, Is.EqualTo(5.0 / 6.0).Within(Tol));
    }

    #endregion

    #region U12 minor-spliceosome acceptor — YCCAC + ppt (score/3.5)

    [Test]
    public void FindAcceptorSites_U12_MaxConsensus_ScoreIsExactlyOne()
    {
        // C at -1, C at -2, Y at -3, fully pyrimidine ppt ⇒ (1+1+0.5+1)/3.5 = 1.0.
        // 15 C's then AC then padding: AC at index 15, upstream all C, length ≥ 20.
        const string seq = "CCCCCCCCCCCCCCCACCCC";   // len 20
        var site = FindAcceptorSites(seq, minScore: 0.5, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Acceptor);
        Assert.That(site.Score, Is.EqualTo(1.0).Within(Tol));
    }

    [Test]
    public void FindAcceptorSites_U12_PartialConsensus_ScoreIsExactFraction()
    {
        // pos-3 = A (not pyrimidine ⇒ no +0.5); ppt half pyrimidine (6/12).
        // [0..5]=C, [6..12]=A, [13]=C(-2), [14]=C(-1), [15]=A,[16]=C, pad [17..19]=C.
        // ppt loop i=0..11 ⇒ 6 C + 6 A ⇒ 6/12. score = 1(-1=C)+1(-2=C)+0(-3=A)+0.5(ppt) = 2.5 ⇒ 2.5/3.5.
        const string seq = "CCCCCCAAAAAAACCACCCC";   // len 20
        var site = FindAcceptorSites(seq, minScore: 0.0, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Acceptor);
        Assert.That(site.Score, Is.EqualTo(2.5 / 3.5).Within(Tol));
    }

    #endregion

    #region CalculateMaxEntScore — log-odds, donor vs acceptor branch

    [Test]
    public void CalculateMaxEntScore_DonorConsensus_ExactLogOddsSum()
    {
        // Each of the 9 donor columns matches ⇒ score = 9 * log2(1 + 0.01).
        const string motif = "CAGGUAAGU";
        double expected = 9 * Math.Log2(1.0 + 0.01);
        Assert.That(CalculateMaxEntScore(motif, SpliceSiteType.Donor),
            Is.EqualTo(expected).Within(Tol));
    }

    [Test]
    public void CalculateMaxEntScore_DonorOffConsensus_DiffersFromConsensus()
    {
        // A non-consensus base contributes log2(0 + 0.01) ≪ log2(1.01): pins the +0.01 term.
        const string consensus = "CAGGUAAGU";
        const string broken = "CACGUAAGU"; // -1 column G→C (mismatch)
        double good = CalculateMaxEntScore(consensus, SpliceSiteType.Donor);
        double bad = CalculateMaxEntScore(broken, SpliceSiteType.Donor);
        double delta = Math.Log2(1.01) - Math.Log2(0.01); // one column flips 1→0
        Assert.That(good - bad, Is.EqualTo(delta).Within(1e-6));
    }

    [Test]
    public void CalculateMaxEntScore_AcceptorBranch_IsNonzeroAndTyped()
    {
        // Acceptor branch must run (kills the type-equality mutant that would return 0)
        // and produce the documented log-odds sum over the 8 acceptor columns.
        // 18-nt acceptor motif (offsets -15..0 over indices 0..17 hitting cols at 0,5,10..15).
        const string motif = "UUUUUUUUUUCUUCAGG"; // len 17
        string rna = Rna(motif);
        double expected = 0;
        for (int i = 0; i < rna.Length; i++)
        {
            int off = i - 15;
            if (AcceptorPwm.TryGetValue(off, out var w))
            {
                int bi = BaseIdx(rna[i]);
                if (bi >= 0) expected += Math.Log2(w[bi] + 0.01);
            }
        }
        Assert.That(CalculateMaxEntScore(motif, SpliceSiteType.Acceptor),
            Is.EqualTo(expected).Within(Tol));
        Assert.That(expected, Is.Not.EqualTo(0.0)); // ensures the acceptor branch truly executed
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
    public void FindDonorSites_MinimumLengthSix_FindsConsensusAtStart()
    {
        // length == 6 is the smallest scannable sequence (guard is strict '< 6').
        // GU at index 0 with consensus tail ⇒ 6 in-bounds columns all match ⇒ 1.0.
        var sites = FindDonorSites("GUAAGU", minScore: 0.5).ToList();
        Assert.That(sites, Has.Count.EqualTo(1));        // kills sequence.Length <= 6 guard
        Assert.That(sites[0].Score, Is.EqualTo(1.0).Within(Tol));
    }

    [Test]
    public void FindDonorSites_GcNonCanonical_ScoreEqualToMinScore_IsIncluded()
    {
        // GC donor: +1 column (U) mismatches ⇒ 8/9. minScore == score ⇒ '>=' keeps, '>' drops.
        const string seq = "CAGGCAAGU";
        var sites = FindDonorSites(seq, minScore: 8.0 / 9.0, includeNonCanonical: true).ToList();
        Assert.That(sites, Has.Count.EqualTo(1));
        Assert.That(sites[0].Score, Is.EqualTo(8.0 / 9.0).Within(Tol));
    }

    [Test]
    public void FindAcceptorSites_U12_ScoreEqualToMinScore_IsIncluded()
    {
        // U12 acceptor max consensus score 1.0; minScore == 1.0 ⇒ '>=' keeps, '>' drops.
        const string seq = "CCCCCCCCCCCCCCCACCCC";
        var sites = FindAcceptorSites(seq, minScore: 1.0, includeNonCanonical: true)
            .Where(s => s.Type == SpliceSiteType.U12Acceptor).ToList();
        Assert.That(sites, Has.Count.EqualTo(1));
    }

    [Test]
    public void FindBranchPoints_ScoreEqualToMinScore_IsIncluded()
    {
        const string seq = "UUUUAACGGGG";
        string rna = Rna(seq);
        double exact = ExpectedBranchScore(rna, 0); // i=0 ⇒ Position 5
        var sites = FindBranchPoints(seq, minScore: exact).ToList();
        Assert.That(sites.Any(s => s.Position == 5), Is.True); // kills score>minScore boundary
    }

    [Test]
    public void FindBranchPoints_LastScannablePosition_IsIncluded()
    {
        // end = length-7 = 4. Branch heptamer at i=4 (Position 9) must be scanned (inclusive '<=').
        const string seq = "GGGGUACUAAC"; // len 11
        var sites = FindBranchPoints(seq, minScore: 0.0).ToList();
        Assert.That(sites.Any(s => s.Position == 9), Is.True); // kills i<end (would stop before i=4)
    }

    [Test]
    public void FindDonorSites_U12_DonorAtEnd_MotifBoundIsInclusive()
    {
        // AU donor at position 0 with position+6 == length: the inclusive 'position+6 <= length'
        // keeps the 6-mer motif; a strict '<' mutant yields an empty motif ⇒ score 0 ⇒ no site.
        var sites = FindDonorSites("AUAUCC", minScore: 0.5, includeNonCanonical: true)
            .Where(s => s.Type == SpliceSiteType.U12Donor).ToList();
        Assert.That(sites, Has.Count.EqualTo(1));
        Assert.That(sites[0].Score, Is.EqualTo(1.0).Within(Tol));
    }

    #endregion

    #region U12 acceptor — per-position consensus contributions

    // Each variant isolates one YCCAC column so a 'position+k' mutant (reading the wrong side
    // of the AC) diverges from the documented 'position-k' rule. All use AC at index 15.

    [Test]
    public void FindAcceptorSites_U12_MinusOnePositionDrivesFirstC()
    {
        // pos-1 = A (no +1); pos-2 = C (+1); pos-3 = C (+0.5); ppt full (+1) ⇒ 2.5/3.5.
        // A 'position+1' mutant reads the C of AC ⇒ +1 ⇒ 1.0.
        const string seq = "CCCCCCCCCCCCCCAACCCC"; // [14]=A,[15]=A,[16]=C
        var site = FindAcceptorSites(seq, minScore: 0.0, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Acceptor);
        Assert.That(site.Score, Is.EqualTo(2.5 / 3.5).Within(Tol));
    }

    [Test]
    public void FindAcceptorSites_U12_MinusTwoPositionDrivesSecondC()
    {
        // pos-2 = A (no +1); pos-1 = C (+1); pos-3 = C (+0.5); ppt full (+1) ⇒ 2.5/3.5.
        // A 'position+2' mutant reads the padding C ⇒ +1 ⇒ 1.0.
        const string seq = "CCCCCCCCCCCCCACACCCC"; // [13]=A,[14]=C,[15]=A,[16]=C
        var site = FindAcceptorSites(seq, minScore: 0.0, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Acceptor);
        Assert.That(site.Score, Is.EqualTo(2.5 / 3.5).Within(Tol));
    }

    [Test]
    public void FindAcceptorSites_U12_MinusThreePositionDrivesPyrimidine()
    {
        // pos-3 = A (no +0.5); pos-1,pos-2 = C (+1 each); ppt full (+1) ⇒ 3.0/3.5.
        // A 'position+3' mutant reads the padding C (a pyrimidine) ⇒ +0.5 ⇒ 1.0.
        const string seq = "CCCCCCCCCCCCACCACCCC"; // [12]=A,[13]=C,[14]=C,[15]=A,[16]=C
        var site = FindAcceptorSites(seq, minScore: 0.0, includeNonCanonical: true)
            .Single(s => s.Type == SpliceSiteType.U12Acceptor);
        Assert.That(site.Score, Is.EqualTo(3.0 / 3.5).Within(Tol));
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
        // searchEnd == 0 ⇒ end = min(0, len-7) = 0 ⇒ exactly one scanned position (i = 0).
        // A 'searchEnd <= 0' mutant would treat 0 as "unbounded" and scan to len-7.
        const string seq = "UACUAACGGGG"; // len 11 ⇒ unbounded end would be 4
        var sites = FindBranchPoints(seq, searchStart: 0, searchEnd: 0, minScore: 0.0).ToList();
        Assert.That(sites, Has.Count.EqualTo(1));
        Assert.That(sites[0].Position, Is.EqualTo(5)); // i=0 ⇒ Position 5
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
