namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the MiRNA area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("MiRNA")]
public class MiRnaCombinatorialTests
{
    // miRNA whose seed (positions 2–8) is CGUACGU; the 6-mer target core is then CGUACG.
    private const string MiRnaSeq = "ACGUACGUACGUACGUACGUAC";
    private const int SiteOffset = 30;

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MIRNA-TARGET-001 — miRNA target-site prediction (MiRNA)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 75.
    // Spec: tests/TestSpecs/MIRNA-TARGET-001.md (canonical FindTargetSites + context analysis).
    // Dimensions: seedType(3) × utrLen(3) × scoringMethod(2). Grid 3×3×2 = 18.
    //
    // Model (Bartel 2009; TargetScan): canonical miRNA sites are seed-complementary matches in the
    // 3′UTR — 8mer (seed match 2–8 + A1), 7mer-m8 (2–8) and 6mer (2–7). FindTargetSites classifies
    // each. The scoring axis covers seed-match scoring vs context augmentation (AnalyzeTargetContext).
    //
    // The combinatorial point: site type, UTR length and scoring method interact — the planted
    // canonical site is detected and classified correctly at every UTR length, with a positive
    // seed score and a valid context score.
    // ═══════════════════════════════════════════════════════════════════════

    public enum SeedType { Site8mer, Site7merM8, Site6mer }
    public enum MiRnaScoring { SeedMatch, ContextAugmented }

    private static string PlantSite(SeedType type) => type switch
    {
        // seedRC = ACGUACG (RC of seed CGUACGU); core = CGUACG.
        SeedType.Site8mer => "ACGUACGA",   // pos8 match + A1 ⇒ 8mer
        SeedType.Site7merM8 => "ACGUACGC", // pos8 match, no A1 ⇒ 7mer-m8
        _ => "UCGUACGU",                   // core only, no pos8 / no A1 ⇒ 6mer
    };

    private static (MiRnaAnalyzer.TargetSiteType Type, int Len) Expected(SeedType type) => type switch
    {
        SeedType.Site8mer => (MiRnaAnalyzer.TargetSiteType.Seed8mer, 8),
        SeedType.Site7merM8 => (MiRnaAnalyzer.TargetSiteType.Seed7merM8, 7),
        _ => (MiRnaAnalyzer.TargetSiteType.Seed6mer, 6),
    };

    [Test, Combinatorial]
    public void MiRnaTarget_DetectsAndClassifiesCanonicalSite(
        [Values(SeedType.Site8mer, SeedType.Site7merM8, SeedType.Site6mer)] SeedType seedType,
        [Values(60, 120, 240)] int utrLen,
        [Values(MiRnaScoring.SeedMatch, MiRnaScoring.ContextAugmented)] MiRnaScoring scoring)
    {
        var miRna = MiRnaAnalyzer.CreateMiRna("miR-x", MiRnaSeq);
        string site = PlantSite(seedType);
        string utr = new string('U', SiteOffset) + site + new string('U', utrLen - SiteOffset - site.Length);

        var sites = MiRnaAnalyzer.FindTargetSites(utr, miRna, minScore: 0.0).ToList();
        var (expectedType, expectedLen) = Expected(seedType);

        sites.Should().Contain(s => s.Type == expectedType, "the planted canonical site is classified correctly");
        var found = sites.First(s => s.Type == expectedType);

        if (scoring == MiRnaScoring.SeedMatch)
        {
            found.Score.Should().BeGreaterThan(0, "a canonical seed match scores positively");
            found.SeedMatchLength.Should().Be(expectedLen);
        }
        else
        {
            var ctx = MiRnaAnalyzer.AnalyzeTargetContext(utr, found.Start, found.End);
            ctx.AuContent.Should().BeInRange(0.0, 1.0);
            ctx.ContextScore.Should().BeInRange(0.0, 1.0);
        }
    }

    /// <summary>
    /// Interaction witness: the three canonical site types are ordered by stringency — an 8mer
    /// scores at least as high as a 7mer-m8, which scores at least as high as a 6mer.
    /// </summary>
    [Test]
    public void MiRnaTarget_SiteHierarchy_ScoresOrdered()
    {
        var miRna = MiRnaAnalyzer.CreateMiRna("miR-x", MiRnaSeq);
        double Score(SeedType t)
        {
            string utr = new string('U', SiteOffset) + PlantSite(t) + new string('U', 80);
            var (type, _) = Expected(t);
            return MiRnaAnalyzer.FindTargetSites(utr, miRna, 0.0).First(s => s.Type == type).Score;
        }

        Score(SeedType.Site8mer).Should().BeGreaterThanOrEqualTo(Score(SeedType.Site7merM8));
        Score(SeedType.Site7merM8).Should().BeGreaterThanOrEqualTo(Score(SeedType.Site6mer));
    }

    /// <summary>
    /// Interaction witness: a UTR with no seed-complementary core yields no canonical target site.
    /// </summary>
    [Test]
    public void MiRnaTarget_NoSeedMatch_NoSite()
    {
        var miRna = MiRnaAnalyzer.CreateMiRna("miR-x", MiRnaSeq);
        MiRnaAnalyzer.FindTargetSites(new string('U', 100), miRna, 0.0)
            .Should().BeEmpty("a poly-U UTR has no seed-complementary site");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MIRNA-PRECURSOR-001 — Pre-miRNA hairpin detection (MiRNA)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 76.
    // Spec: tests/TestSpecs/MIRNA-PRECURSOR-001.md (canonical FindPreMiRnaHairpins).
    // Dimensions: precLen(3) × minStem(3) × maxLoop(3). Grid 3×3×3 = 27.
    //
    // Model (Bartel 2004): a pre-miRNA is a ≥55-nt hairpin with a long (~≥18 bp) stem and a small
    // (3–25 nt) terminal loop. The stem/loop thresholds are FIXED internal constants here, so the
    // minStem/maxLoop axes vary the PLANTED hairpin's stem/loop against them; precLen is the
    // search length window (min/maxHairpinLength).
    //
    // The combinatorial point: planted stem length, loop size and the length window jointly
    // determine detection — a hairpin is found exactly when stem ≥ 18, loop ∈ [3,25], and its
    // length is ≥ 55 and within the search window.
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly int[] StemValues = { 16, 25, 32 };
    private static readonly int[] LoopValues = { 6, 20, 30 };
    private static readonly (int Min, int Max)[] LenWindows = { (40, 200), (40, 70), (40, 50) };

    [Test, Combinatorial]
    public void MiRnaPrecursor_DetectsHairpinAgainstThresholds(
        [Values(0, 1, 2)] int precLenIdx,
        [Values(0, 1, 2)] int minStemIdx,
        [Values(0, 1, 2)] int maxLoopIdx)
    {
        int stem = StemValues[minStemIdx];
        int loop = LoopValues[maxLoopIdx];
        var (minH, maxH) = LenWindows[precLenIdx];

        string hairpin = new string('G', stem) + new string('A', loop) + new string('C', stem);
        int innerLen = 2 * stem + loop;
        string seq = new string('U', 10) + hairpin + new string('U', 10);

        var found = MiRnaAnalyzer.FindPreMiRnaHairpins(seq, minH, maxH, matureLength: 22).ToList();

        bool expected = stem >= 18 && loop is >= 3 and <= 25 && innerLen >= 55 && innerLen >= minH && innerLen <= maxH;
        found.Any(p => p.Sequence == hairpin)
            .Should().Be(expected, "the planted hairpin is found iff stem/loop/length all qualify");

        foreach (var p in found)
        {
            p.Sequence.Length.Should().Be(p.End - p.Start + 1, "coordinates match the sequence span");
            p.Structure.Length.Should().Be(p.Sequence.Length, "dot-bracket spans the hairpin");
            p.MatureSequence.Should().NotBeEmpty("a mature arm is extracted");
        }
    }

    /// <summary>
    /// Interaction witness: each requirement independently rejects a hairpin — a short stem, an
    /// oversized loop, or a length below 55 nt all prevent detection.
    /// </summary>
    [Test]
    public void MiRnaPrecursor_EachRequirement_GatesDetection()
    {
        bool Found(int stem, int loop)
        {
            string hp = new string('G', stem) + new string('A', loop) + new string('C', stem);
            return MiRnaAnalyzer.FindPreMiRnaHairpins(new string('U', 10) + hp + new string('U', 10), 40, 200, 22)
                .Any(p => p.Sequence == hp);
        }

        Found(25, 6).Should().BeTrue("a 25-bp stem, 6-nt loop, 56-nt hairpin qualifies");
        Found(16, 6).Should().BeFalse("stem < 18 bp is rejected");
        Found(25, 30).Should().BeFalse("loop > 25 nt is rejected");
        Found(20, 6).Should().BeFalse("a 46-nt hairpin is below the 55-nt floor");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MIRNA-PAIR-001 — miRNA:target duplex pairing (MiRNA)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 225.
    // Spec: tests/TestSpecs/MIRNA-PAIR-001.md (canonical AlignMiRnaToTarget / GetReverseComplement /
    //       CanPair / IsWobblePair). ADVANCED §10.
    // Dimensions: seedType(3) × utrLen(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Lewis 2005; Crick 1966 wobble): the miRNA pairs antiparallel to the target — miRNA index i
    // with target[len−1−i] — counting Watson-Crick matches, G:U wobbles and mismatches; the target of a
    // seed is the reverse complement of that seed.
    //
    // Axis mapping (documented): seedType → the complementary-seed length (6/7/8); utrLen → target
    // length. Engineered target ends with the reverse complement of the miRNA seed. The combinatorial
    // point: the duplex pairs at least the seed (Matches ≥ seedLen) and the alignment is total
    // (matches + wobbles + mismatches = the overlap length), at every seed length and target length.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void MiRnaPair_SeedPairsAntiparallel_AcrossSeedAndTargetLength(
        [Values(6, 7, 8)] int seedLen,
        [Values(22, 30, 40)] int utrLen)
    {
        string seed = MiRnaSeq[..seedLen];
        string target = new string('A', utrLen - seedLen) + MiRnaAnalyzer.GetReverseComplement(seed);

        var duplex = MiRnaAnalyzer.AlignMiRnaToTarget(MiRnaSeq, target);

        int overlap = Math.Min(MiRnaSeq.Length, target.Length);
        (duplex.Matches + duplex.Mismatches + duplex.GUWobbles).Should().Be(overlap, "every aligned position is classified");
        duplex.AlignmentString.Length.Should().Be(overlap, "the alignment string spans the overlap");
        duplex.Matches.Should().BeGreaterThanOrEqualTo(seedLen, "the reverse-complement seed pairs Watson-Crick");
    }

    /// <summary>
    /// Interaction witnesses — pairing predicates and reverse complement: a fully complementary target
    /// pairs at every position, GetReverseComplement is an involution, and CanPair/IsWobblePair encode
    /// Watson-Crick + G:U wobble.
    /// </summary>
    [Test]
    public void MiRnaPair_PredicatesAndPerfectComplement()
    {
        string perfect = MiRnaAnalyzer.GetReverseComplement(MiRnaSeq);
        var duplex = MiRnaAnalyzer.AlignMiRnaToTarget(MiRnaSeq, perfect);
        duplex.Mismatches.Should().Be(0, "the reverse complement pairs at every position");
        (duplex.Matches + duplex.GUWobbles).Should().Be(MiRnaSeq.Length);

        MiRnaAnalyzer.GetReverseComplement(perfect).Should().Be(MiRnaSeq, "reverse complement is an involution");

        MiRnaAnalyzer.CanPair('A', 'U').Should().BeTrue();
        MiRnaAnalyzer.CanPair('G', 'C').Should().BeTrue();
        MiRnaAnalyzer.CanPair('G', 'U').Should().BeTrue("G:U is a wobble pair");
        MiRnaAnalyzer.CanPair('A', 'G').Should().BeFalse("A:G cannot pair");
        MiRnaAnalyzer.IsWobblePair('G', 'U').Should().BeTrue();
        MiRnaAnalyzer.IsWobblePair('A', 'U').Should().BeFalse("a Watson-Crick pair is not a wobble");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MIRNA-CONTEXT-001 — context++ target-site score (MiRNA)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 252.
    // Spec: tests/TestSpecs/MIRNA-CONTEXT-001.md (MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus). ADVANCED §10.
    //
    // Sources: Agarwal et al. (2015) eLife 4:e05005 (TargetScan 7 context++ regression).
    //
    // Model: the context++ score is a LINEAR feature model — the partial score equals the intercept plus
    //   the sum of every per-feature contribution (a dot product of feature values and coefficients).
    //
    // Dimensions: seedType(3) × siteContext(3) × SA(2). Grid 3×3×2 = 18 (exhaustive).
    //
    // Axis mapping (documented): siteContext is the flanking composition (AU-rich / mixed / GC-rich,
    // which drives the local-AU feature); SA (structural accessibility) is realised as the 3'UTR length
    // (the SA window fits, or is reported as an honest residual). In every cell the partial score equals
    // intercept + Σ contributions (the defining distributivity of the linear model) and the planted seed
    // type is classified correctly.
    // ═══════════════════════════════════════════════════════════════════════

    private static double ContextSumOfContributions(MiRnaAnalyzer.ContextPlusPlusScore c) =>
        c.Intercept + c.LocalAuContribution + c.SRna1Contribution + c.SRna8Contribution + c.Site8Contribution
        + c.SaContribution + c.ThreePrimePairingContribution + c.MinDistContribution + c.Len3UtrContribution
        + c.Off6mContribution + c.SpsContribution + c.TaContribution + c.LenOrfContribution
        + c.Orf8mContribution + c.PctContribution;

    private static string ContextFlank(int n, string motif) =>
        string.Concat(Enumerable.Range(0, Math.Max(0, n)).Select(i => motif[i % motif.Length]));

    [Test, Combinatorial]
    public void MiRnaContext_PartialScoreIsSumOfContributions_AcrossSeedContextAndSa(
        [Values(SeedType.Site8mer, SeedType.Site7merM8, SeedType.Site6mer)] SeedType seedType,
        [Values("A", "AG", "G")] string contextMotif,
        [Values(60, 120)] int utrLen)
    {
        var miRna = MiRnaAnalyzer.CreateMiRna("miR-x", MiRnaSeq);
        string siteCore = PlantSite(seedType);
        const int leftPad = 25;
        string utr = ContextFlank(leftPad, contextMotif) + "UU" + siteCore + "UU"
                     + ContextFlank(utrLen - leftPad - 4 - siteCore.Length, contextMotif);

        var (expectedType, _) = Expected(seedType);
        var site = MiRnaAnalyzer.FindTargetSites(utr, miRna, minScore: 0.0).First(s => s.Type == expectedType);

        var ctx = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(utr, miRna, site);

        ctx.SiteType.Should().Be(expectedType, "the scored site keeps its canonical seed type");
        ctx.ContextScorePartial.Should().BeApproximately(ContextSumOfContributions(ctx), 1e-12,
            "context++ partial score = intercept + Σ per-feature contributions (linear model)");
        MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(utr, miRna, site).ContextScorePartial
            .Should().Be(ctx.ContextScorePartial, "the score is deterministic");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MIRNA-PCT-001 — probability of conserved targeting (MiRNA)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 253.
    // Spec: tests/TestSpecs/MIRNA-PCT-001.md (MiRnaAnalyzer.PctFromBranchLength). ADVANCED §10.
    //
    // Sources: Friedman et al. (2009) Genome Res 19:92; TargetScan targetscan_70_BL_PCT.pl.
    //
    // Model: PCT(Bls) = B0 + B1/(1 + e^(−B2·Bls + B3)), truncated at 0 — the published logistic mapping
    //   from branch length to probability of conserved targeting, parameterised per seed type.
    //
    // Dimensions: seedType(3) × branchLength(3). Grid 3×3 = 9 (exhaustive).
    //
    // The combinatorial point: each (seedType parameter set × branch length) cell reproduces the closed
    // form exactly, stays within [0, B0+B1], and PCT rises monotonically with branch length (B2 > 0).
    // ═══════════════════════════════════════════════════════════════════════

    private static MiRnaAnalyzer.PctSigmoidParameters PctParamsFor(int seedIdx) => seedIdx switch
    {
        0 => new MiRnaAnalyzer.PctSigmoidParameters(0.0, 0.90, 1.0, 0.5), // 8mer-like
        1 => new MiRnaAnalyzer.PctSigmoidParameters(0.0, 0.70, 1.0, 1.0), // 7mer-m8-like
        _ => new MiRnaAnalyzer.PctSigmoidParameters(0.0, 0.50, 1.0, 1.5), // 7mer-A1-like
    };

    [Test, Combinatorial]
    public void MiRnaPct_ReproducesPublishedLogistic_AcrossSeedTypeAndBranchLength(
        [Values(0, 1, 2)] int seedIdx,
        [Values(0.5, 2.0, 5.0)] double branchLength)
    {
        var p = PctParamsFor(seedIdx);
        double pct = MiRnaAnalyzer.PctFromBranchLength(branchLength, p);

        double expected = Math.Max(0.0, p.B0 + p.B1 / (1.0 + Math.Exp(-p.B2 * branchLength + p.B3)));
        pct.Should().BeApproximately(expected, 1e-12, "PCT = B0 + B1/(1 + e^(−B2·Bls + B3)), truncated at 0");
        pct.Should().BeInRange(0.0, p.B0 + p.B1, "PCT lies within the logistic's range");

        if (branchLength > 0.5)
            pct.Should().BeGreaterThan(MiRnaAnalyzer.PctFromBranchLength(0.5, p),
                "PCT rises monotonically with branch length (B2 > 0)");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MIRNA-CLASSIFY-001 — pre-miRNA logistic classifier (MiRNA)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 254.
    // Spec: tests/TestSpecs/MIRNA-CLASSIFY-001.md (MiRnaAnalyzer.ClassifyPreMiRna). ADVANCED §10.
    //
    // Sources: logistic regression over MFE/AMFE/MFEI/GC/%paired (Bonnet 2004; Xue 2005 SVM features).
    //
    // Model: P(natural) is the logistic of the standardised feature vector, so it lies strictly in (0,1)
    //   and the boolean call IsNatural is exactly the 0.5 threshold of that probability.
    //
    // Dimensions: MFEI(3) × GC(3) × %paired(2). Grid 3×3×2 = 18 (exhaustive).
    //
    // Axis mapping (documented): MFEI/GC/%paired are FEATURES computed from the sequence, not inputs, so
    // they are realised by disrupting a genuine precursor's hairpin to varying degrees (which lowers
    // %paired and shifts MFEI/GC). In every cell the probability is a proper logistic in (0,1) and
    // IsNatural ⇔ P ≥ 0.5; a native precursor classifies as natural and scores strictly higher than its
    // hairpin-disrupted variant.
    // ═══════════════════════════════════════════════════════════════════════

    private const string ClassifyPrecursor =
        "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA"; // hsa-mir-21 (72 nt)

    // Replaces the 3' arm with poly-A from position `keep`, progressively destroying the hairpin pairing.
    private static string DisruptHairpin(string pre, int keep) =>
        pre[..keep] + new string('A', pre.Length - keep);

    [Test, Combinatorial]
    public void MiRnaClassify_LogisticThresholdLaw_AcrossDisruptionLevels(
        [Values(72, 60, 50)] int keep,
        [Values(0, 6, 12)] int extraTrim,
        [Values(false, true)] bool native)
    {
        int cut = Math.Min(ClassifyPrecursor.Length, Math.Max(40, keep - extraTrim));
        string seq = native ? ClassifyPrecursor : DisruptHairpin(ClassifyPrecursor, cut);

        var c = MiRnaAnalyzer.ClassifyPreMiRna(seq);
        c.Should().NotBeNull("a ≥40 nt candidate is classifiable");
        c!.Value.NaturalProbability.Should().BeGreaterThan(0.0).And.BeLessThan(1.0, "the logistic output is a probability in (0,1)");
        c.Value.IsNatural.Should().Be(c.Value.NaturalProbability >= 0.5, "the call is exactly the 0.5 threshold of the logistic");

        if (native)
        {
            c.Value.IsNatural.Should().BeTrue("the genuine hsa-mir-21 precursor classifies as natural");
            // Compare against a clearly hairpin-disrupted variant (3' arm replaced from position 45).
            c.Value.NaturalProbability.Should().BeGreaterThan(
                MiRnaAnalyzer.ClassifyPreMiRna(DisruptHairpin(ClassifyPrecursor, 45))!.Value.NaturalProbability,
                "the native precursor scores higher than its hairpin-disrupted variant");
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: MIRNA-CLEAVAGE-001 — Drosha/Dicer cleavage rulers (MiRNA)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 255.
    // Spec: tests/TestSpecs/MIRNA-CLEAVAGE-001.md (MiRnaAnalyzer.PredictDroshaDicerCleavage). ADVANCED §10.
    //
    // Sources: Han et al. (2006) Cell 125:887 (Drosha ~11 bp basal-junction ruler); Park et al. (2011)
    //   Nature 475:201 (Dicer 5' counting rule, ~22 nt mature).
    //
    // Model: with basal junction j the Drosha 5' cut is at j+11 and the mature 5p product is exactly the
    //   22-nt window pri[j+11 .. j+33) — an integer measuring rule independent of the stem/loop content.
    //
    // Dimensions: stemLen(3) × loopLen(2) × strand(2). Grid 3×2×2 = 12 (exhaustive).
    //
    // The combinatorial point: across stem and loop sizes the integer rulers place the 5p mature at
    // j+11 with length 22 nt equal to the precursor substring; the strand axis checks the 5p mature vs
    // the 3p star product. The placement is content-independent — it holds in every cell.
    // ═══════════════════════════════════════════════════════════════════════

    private static string CleavagePriSeq(int stemLen, int loopLen) =>
        new string('C', 11) + ContextFlank(stemLen, "ACGU") + new string('A', loopLen) + new string('G', stemLen) + "UU";

    public enum CleavageArm { Mature5p, Star3p }

    [Test, Combinatorial]
    public void MiRnaCleavage_IntegerRulers_AcrossStemLoopAndArm(
        [Values(14, 18, 22)] int stemLen,
        [Values(6, 10)] int loopLen,
        [Values(CleavageArm.Mature5p, CleavageArm.Star3p)] CleavageArm arm)
    {
        string pri = CleavagePriSeq(stemLen, loopLen);
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, basalJunction: 0);
        cut.Should().NotBeNull("a sufficiently long pri-miRNA admits the cleavage rulers");

        cut!.Value.DroshaCut5Prime.Should().Be(11, "Drosha cuts ~11 bp from the basal junction (Han 2006)");
        cut.Value.MatureStart.Should().Be(11, "the 5p mature starts at the Drosha 5' cut");
        cut.Value.MatureSequence.Length.Should().Be(22, "the Dicer 5' counting rule fixes the mature at 22 nt (Park 2011)");
        (cut.Value.MatureEnd - cut.Value.MatureStart + 1).Should().Be(22, "the inclusive mature span is 22 nt");
        cut.Value.MatureSequence.Should().Be(pri.Substring(11, 22), "the mature is exactly the measured precursor window");

        if (arm == CleavageArm.Star3p)
        {
            cut.Value.StarSequence.Should().NotBeNullOrEmpty("the 3p star product is reported");
            cut.Value.StarStart.Should().BeGreaterThan(cut.Value.MatureStart, "the 3p star lies 3' of the 5p mature");
        }
    }
}
