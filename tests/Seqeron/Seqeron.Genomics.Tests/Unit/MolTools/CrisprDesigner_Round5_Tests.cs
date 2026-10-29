namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Findings of audit round 5 on the CRISPR slice (B08 A5-7 .. A5-10): argument contracts that the round-4 pass left open.
/// Each test failed on the unfixed code before the fix (the auditor's reproductions were re-run here first).
/// TestSpec: TestSpecs/CRISPR-PAM-001.md, CRISPR-OFF-001.md, CRISPR-GUIDE-001.md ("Round 5 contracts").
/// </summary>
[TestFixture]
public class CrisprDesigner_Round5_Tests
{
    private const string Guide = "GACGCATAAAGATGAGACGC";

    // ------------------------------------------------------------------------------------------------
    // A5-7: FindPamSites(string) - only A/C/G/T/N reach the scan (CRISPOR's front end restricts the input to ACGTN;
    // here the string overload is that front end). Before: any other character inside the protospacer window came
    // back as a "guide", which EvaluateGuideRna / DesignGuideRnas then rejected.
    // ------------------------------------------------------------------------------------------------

    [TestCase("ACGT12345ACGTACGTACGTAGG", '1', 4)]
    [TestCase("ACGTACGT ACGTACGTACGTAGG", ' ', 8)]
    [TestCase("ACGTACGTRCGTACGTACGTAGG", 'R', 8)]
    [TestCase("ACGTACGTXCGTACGTACGTAGG", 'X', 8)]
    [TestCase("ACGTACGT-CGTACGTACGTAGG", '-', 8)]
    [TestCase("ACGTACGTUCGTACGTACGTAGG", 'U', 8)]
    [TestCase("ACGTACGT\nCGTACGTACGTAGG", '\n', 8)]
    [TestCase("ACGTACGTŁCGTACGTACGTAGG", 'Ł', 8)]
    public void FindPamSites_StringOverload_RejectsCharactersOutsideACGTN_AtTheCall(string sequence, char bad, int index)
    {
        // No ToList(): the check is eager, like the unknown-system check.
        var e = Assert.Throws<ArgumentException>(() => CrisprDesigner.FindPamSites(sequence));
        Assert.Multiple(() =>
        {
            Assert.That(e!.ParamName, Is.EqualTo("sequence"));
            Assert.That(e.Message, Does.Contain($"'{bad}'").And.Contain($"index {index}"));
        });
    }

    [Test]
    public void FindPamSites_StringOverload_RejectsABadCharacterOutsideAnyWindowToo()
    {
        // The bad character sits in the tail that no protospacer window covers; the input is still not a nucleotide string.
        Assert.Throws<ArgumentException>(() => CrisprDesigner.FindPamSites("ACGTACGTACGTACGTACGTAGGTTTT5"));
    }

    [Test]
    public void FindPamSites_StringOverload_AcceptsLowerCaseAndN_AsBefore()
    {
        const string upper = "ACGTACGTACGTACGTACGTAGGACGTNACGTACGTACGTACGTAGG";
        var a = CrisprDesigner.FindPamSites(upper).Select(s => (s.Position, s.TargetSequence, s.IsForwardStrand)).ToList();
        var b = CrisprDesigner.FindPamSites(upper.ToLowerInvariant()).Select(s => (s.Position, s.TargetSequence, s.IsForwardStrand)).ToList();

        Assert.That(b, Is.EqualTo(a), "lower case is folded");
        Assert.That(a, Is.Not.Empty);
        Assert.That(a.Select(t => t.TargetSequence), Has.None.Contain('N'), "a window holding an N is dropped (H12)");
    }

    // ------------------------------------------------------------------------------------------------
    // A5-8: CalculateMitSpecificityScore(guide, IEnumerable<PamAnnotatedOffTarget>) - the documented ArgumentException
    // for a malformed element; before: NullReferenceException / the inner scorer's parameter name / a "ZZZ" PAM that
    // silently down-weighted the hit by 0.2.
    // ------------------------------------------------------------------------------------------------

    private static PamAnnotatedOffTarget Hit(string protospacer, string pam) =>
        new(new OffTargetSite(0, protospacer, 1, new[] { 19 }, true, 0.0), pam);

    [Test]
    public void MitSpecificity_PamAnnotatedHits_NullSiteOrSequence_ThrowArgumentException_NamingTheList()
    {
        var nullSite = new PamAnnotatedOffTarget(null!, "AGG");
        var nullSequence = new PamAnnotatedOffTarget(new OffTargetSite(0, null!, 1, new[] { 19 }, true, 0.0), "AGG");

        foreach (var bad in new[] { nullSite, nullSequence })
        {
            var e = Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateMitSpecificityScore(Guide, new[] { bad }));
            Assert.That(e!.GetType(), Is.EqualTo(typeof(ArgumentException)), "not a NullReferenceException, and not ArgumentNull for a list element");
            Assert.That(e.ParamName, Is.EqualTo("offTargets"));
        }
    }

    [TestCase("ZZZ")]
    [TestCase("A1G")]
    [TestCase("N GG")]
    public void MitSpecificity_PamAnnotatedHits_RejectsAPamOutsideACGT(string pam)
    {
        var e = Assert.Throws<ArgumentException>(() =>
            CrisprDesigner.CalculateMitSpecificityScore(Guide, new[] { Hit("GACGCATAAAGATGAGACGT", pam) }));
        Assert.That(e!.ParamName, Is.EqualTo("offTargets"));
    }

    [Test]
    public void MitSpecificity_PamAnnotatedHits_BadProtospacerAlphabet_BlamesTheList_NotTheInnerScorer()
    {
        var e = Assert.Throws<ArgumentException>(() =>
            CrisprDesigner.CalculateMitSpecificityScore(Guide, new[] { Hit("NACGCATAAAGATGAGACGT", "AGG") }));
        Assert.That(e!.ParamName, Is.EqualTo("offTargets"), "was 'offTarget20', a parameter of a private helper");
    }

    [Test]
    public void MitSpecificity_PamAnnotatedHits_LowerCaseAndUracilPam_AreFolded_AsForTheProtospacer()
    {
        double upper = CrisprDesigner.CalculateMitSpecificityScore(Guide, new[] { Hit("GACGCATAAAGATGAGACGT", "AGG") });
        double lower = CrisprDesigner.CalculateMitSpecificityScore(Guide, new[] { Hit("gacgcataaagatgagacgt", "agg") });

        Assert.That(lower, Is.EqualTo(upper));
    }

    // ------------------------------------------------------------------------------------------------
    // A5-9: GuideRnaParameters.Validate - an undefined Ranking and an invalid backbone region are refused at the
    // call, with ParamName "parameters", from both DesignGuideRnas and EvaluateGuideRna.
    // ------------------------------------------------------------------------------------------------

    private static readonly DnaSequence Locus = new("ACGGATCTAGACGGATCTAGAGG" + "TTACGATCGATCGTAGCTAGCTAGTCGATCGATTGG");

    [Test]
    public void GuideRnaParameters_UndefinedRanking_IsRefused_ByBothEntryPoints()
    {
        var bad = GuideRnaParameters.Default with { Ranking = (GuideRnaRanking)99 };

        // No enumeration: DesignGuideRnas validates the parameters eagerly.
        var design = Assert.Throws<ArgumentOutOfRangeException>(() => CrisprDesigner.DesignGuideRnas(Locus, 0, Locus.Length - 1, CrisprSystemType.SpCas9, bad));
        var evaluate = Assert.Throws<ArgumentOutOfRangeException>(() => CrisprDesigner.EvaluateGuideRna("ACGGATCTAGACGGATCTAG", CrisprSystemType.SpCas9, bad));
        Assert.That(design!.ParamName, Is.EqualTo("parameters"));
        Assert.That(evaluate!.ParamName, Is.EqualTo("parameters"));
    }

    [TestCase("AGGXTAGTCCGT")]
    [TestCase("agg5tagtccgt")]
    public void GuideRnaParameters_InvalidBackboneRegion_IsRefusedAtTheCall_EvenWithoutAPamInTheRegion(string region)
    {
        var bad = GuideRnaParameters.Default with { SelfComplementarityBackboneRegions = new[] { region } };
        const string noPam = "ACACACACACACACACACACACACACACACACACAC"; // no G, so no NGG on either strand... and no CCN either

        var design = Assert.Throws<ArgumentException>(() => CrisprDesigner.DesignGuideRnas(new DnaSequence(noPam), 0, noPam.Length - 1, CrisprSystemType.SpCas9, bad));
        var evaluate = Assert.Throws<ArgumentException>(() => CrisprDesigner.EvaluateGuideRna("ACGGATCTAGACGGATCTAG", CrisprSystemType.SpCas9, bad));
        Assert.That(design!.ParamName, Is.EqualTo("parameters"), "was 'backboneRegions', and only on first enumeration with a PAM present");
        Assert.That(evaluate!.ParamName, Is.EqualTo("parameters"));
    }

    [Test]
    public void GuideRnaParameters_NullBackboneElement_IsRefused()
    {
        var bad = GuideRnaParameters.Default with { SelfComplementarityBackboneRegions = new string[] { null! } };

        Assert.Throws<ArgumentException>(() => CrisprDesigner.EvaluateGuideRna("ACGGATCTAGACGGATCTAG", CrisprSystemType.SpCas9, bad));
    }

    [Test]
    public void GuideRnaParameters_ValidBackboneAndEveryDefinedRanking_StillWork()
    {
        var withBackbone = GuideRnaParameters.Default with
        {
            SelfComplementarityBackboneRegions = new[] { CrisprDesigner.StandardSgRnaBackboneRegion },
        };

        Assert.DoesNotThrow(() => CrisprDesigner.EvaluateGuideRna("ACGGATCTAGACGGATCTAG", CrisprSystemType.SpCas9, withBackbone));
        foreach (var ranking in Enum.GetValues<GuideRnaRanking>())
        {
            Assert.DoesNotThrow(() =>
                CrisprDesigner.DesignGuideRnas(Locus, 0, Locus.Length - 1, CrisprSystemType.SpCas9, GuideRnaParameters.Default with { Ranking = ranking }).ToList(),
                ranking.ToString());
        }
    }

    // ------------------------------------------------------------------------------------------------
    // A5-10: CalculateCfdScore - the whole 3-nt PAM is validated, as CRISPOR's calcCfdScore screens the whole otSeq
    // (protospacer + PAM) with re.search('[^ATCG]', ...) (crispor.py L2172-2185). Before: the first base of a 3-nt PAM
    // was dropped unchecked, so "NGG", "XGG" and "5GG" all scored a perfect 1.
    // ------------------------------------------------------------------------------------------------

    [TestCase("NGG")]
    [TestCase("XGG")]
    [TestCase("5GG")]
    [TestCase("ngg")]
    public void CalculateCfdScore_ThreeNucleotidePam_FirstBaseMustBeACGT(string pam)
    {
        var e = Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateCfdScore(Guide, Guide, pam));
        Assert.That(e!.ParamName, Is.EqualTo("offTargetPam"));
    }

    [TestCase("AGG", 1.0)]
    [TestCase("cgg", 1.0)]
    [TestCase("TGG", 1.0)]
    [TestCase("GG", 1.0)]
    public void CalculateCfdScore_ConcretePam_IsUnchanged(string pam, double expected)
    {
        Assert.That(CrisprDesigner.CalculateCfdScore(Guide, Guide, pam), Is.EqualTo(expected));
    }

    // ------------------------------------------------------------------------------------------------
    // A5-12: boundary behaviours that survived mutation. Expected values come from the documented rules and the
    // reference's comparisons, not from the code under test.
    // ------------------------------------------------------------------------------------------------

    /// <summary>
    /// crispor.py annotateOfftargets filters alternative-PAM hits with <c>if otScore &lt; minScore: continue</c>
    /// (ALTPAMMINSCORE), so a hit whose score EQUALS the cut-off is kept. The NAG hit with one mismatch at index 5 scores
    /// (1 - 0.395) x 100 = 60.5, is down-weighted x0.2 to 12.1, and the aggregate is 100 / (100 + 12.1) x 100.
    /// </summary>
    [Test]
    public void MitSpecificity_AlternativePamHit_WhoseScoreEqualsTheCutoff_IsKept()
    {
        var genome = new DnaSequence("GACGCTTAAAGATGAGACGC" + "TAG" + "TTTTTT");
        double hit = CrisprDesigner.CalculateMitHitScore(Guide, "GACGCTTAAAGATGAGACGC");
        Assert.That(hit, Is.EqualTo(60.5).Within(1e-9), "Hsu weight at index 5 is 0.395");

        double atCutoff = CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(Guide, genome, 4, CrisprSystemType.SpCas9, hit);
        double aboveCutoff = CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(Guide, genome, 4, CrisprSystemType.SpCas9, hit + 1e-6);

        Assert.Multiple(() =>
        {
            Assert.That(atCutoff, Is.EqualTo(10000.0 / 112.1).Within(1e-9), "kept: 100 / (100 + 0.2 x 60.5) x 100");
            Assert.That(aboveCutoff, Is.EqualTo(100.0), "dropped by a cut-off just above its score");
        });
    }

    /// <summary>
    /// The repeat test is PAM-agnostic in the reference (<c>repCount</c> counts BWA alignments of the 20-mer): a guide with
    /// three PAM-less copies has no off-target hit at all, yet with a limit of two it is a repeat and scores 0, in every
    /// aggregate, not 100 as a "no hits first" ordering would give.
    /// </summary>
    [Test]
    public void MitSpecificity_RepeatCheck_PrecedesTheNoHitsShortcut()
    {
        var genome = new DnaSequence(Guide + "TTT" + "ACAC" + Guide + "TTT" + "ACAC" + Guide + "TTT");

        Assert.That(CrisprDesigner.FindOffTargets(Guide, genome, 4), Is.Empty, "no PAM, so no off-target hit");
        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateMitSpecificityScore(Guide, genome, 4, CrisprSystemType.SpCas9, null, null, 2), Is.EqualTo(0.0));
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScore(Guide, genome, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, null, null, 2), Is.EqualTo(0.0));
            Assert.That(CrisprDesigner.CalculateMitSpecificityScore(Guide, genome, 4, CrisprSystemType.SpCas9, null, null, 3), Is.EqualTo(100.0),
                "three alignments are within a limit of three");
        });
    }

    /// <summary>
    /// The library's own heuristic (docs: max(0, 100 - sum(5 per seed mismatch + 2 per distal mismatch)), hits up to
    /// four mismatches): one four-mismatch site costs 4 x 2 = 8, a five-mismatch site is not looked at.
    /// </summary>
    [Test]
    public void HeuristicSpecificity_CountsFourMismatchSites_AndNotFiveMismatchSites()
    {
        string FlipAt(params int[] indices)
        {
            var chars = Guide.ToCharArray();
            foreach (int i in indices)
                chars[i] = chars[i] == 'A' ? 'C' : 'A';
            return new string(chars);
        }

        var four = new DnaSequence(FlipAt(0, 2, 4, 6) + "TGG");
        var five = new DnaSequence(FlipAt(0, 2, 4, 6, 9) + "TGG");

        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateHeuristicSpecificityScore(Guide, four), Is.EqualTo(92.0), "4 distal mismatches x 2");
            Assert.That(CrisprDesigner.CalculateHeuristicSpecificityScore(Guide, five), Is.EqualTo(100.0), "5 mismatches are beyond the heuristic's cap of 4");
        });
    }

    /// <summary>
    /// The heuristic's seed is the PAM-proximal 12 nt (Hsu 2013: 8-12): for a PAM-after guide of 20 that is indices 8..19
    /// (index 8 costs 5, index 7 costs 2); for a PAM-before guide it is indices 0..11 (index 11 costs 5, index 12 costs 2).
    /// </summary>
    [Test]
    public void OffTargetScore_SeedIsTheTwelvePamProximalBases_ForBothGeometries()
    {
        string Flip(string s, int i) => s.Remove(i, 1).Insert(i, s[i] == 'A' ? "C" : "A");
        double Score(string guide, string site, CrisprSystemType type)
            => CrisprDesigner.FindOffTargets(guide, new DnaSequence(site), 3, type).Single(h => h.IsForwardStrand).OffTargetScore;

        const string cas12aGuide = "ATCGGATCTAGCATGCAAGTCCA";
        Assert.Multiple(() =>
        {
            Assert.That(Score(Guide, Flip(Guide, 8) + "TGG", CrisprSystemType.SpCas9), Is.EqualTo(5.0), "index 8 is the first seed base");
            Assert.That(Score(Guide, Flip(Guide, 7) + "TGG", CrisprSystemType.SpCas9), Is.EqualTo(2.0), "index 7 is distal");
            Assert.That(Score(cas12aGuide, "TTTA" + Flip(cas12aGuide, 11) + "GCGCGC", CrisprSystemType.Cas12a), Is.EqualTo(5.0), "index 11 is the last seed base");
            Assert.That(Score(cas12aGuide, "TTTA" + Flip(cas12aGuide, 12) + "GCGCGC", CrisprSystemType.Cas12a), Is.EqualTo(2.0), "index 12 is distal");
        });
    }
}
