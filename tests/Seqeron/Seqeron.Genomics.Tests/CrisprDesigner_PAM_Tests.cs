using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for CRISPR PAM site detection (CRISPR-PAM-001).
/// Covers FindPamSites and GetSystem methods.
/// 
/// Evidence sources:
/// - Wikipedia: Protospacer adjacent motif — PAM definitions, canonical NGG
/// - Wikipedia: CRISPR — System types overview
/// - Wikipedia: Cas9 — SpCas9 NGG 20 nt spacer, NAG tolerated (lower activity)
/// - Wikipedia: Cas12a — TTTV PAM (V=A,C,G), 5' PAM before target
/// - Jinek et al. (2012), Science 337:816 — SpCas9 NGG PAM, 20 nt guide
/// - Zetsche et al. (2015), Cell 163:759 — Cas12a/Cpf1 TTTV PAM, ~23 nt spacer
/// - Anders et al. (2014), Nature 513:569 — PAM recognition structural basis
/// - Ran et al. (2015), Nature 520:186 — SaCas9 NNGRRT PAM, 21 nt guide
/// - Hsu et al. (2013), Nat Biotechnol 31:827 — NAG as secondary SpCas9 PAM
/// - Liu et al. (2019), Nature 566:218 — CasX (Cas12e) TTCN PAM, 20 nt guide
/// </summary>
[TestFixture]
public class CrisprDesigner_PAM_Tests
{
    #region GetSystem Tests - CRISPR System Configuration

    [Test]
    [Description("M1-M3: GetSystem returns correct Name and PAM for all 7 systems")]
    [TestCase(CrisprSystemType.SpCas9, "SpCas9", "NGG")]            // Jinek 2012
    [TestCase(CrisprSystemType.SpCas9_NAG, "SpCas9-NAG", "NAG")]    // Hsu 2013
    [TestCase(CrisprSystemType.SaCas9, "SaCas9", "NNGRRT")]         // Ran 2015
    [TestCase(CrisprSystemType.Cas12a, "Cas12a/Cpf1", "TTTV")]      // Zetsche 2015
    [TestCase(CrisprSystemType.AsCas12a, "AsCas12a", "TTTV")]       // Zetsche 2015
    [TestCase(CrisprSystemType.LbCas12a, "LbCas12a", "TTTV")]       // Zetsche 2015
    [TestCase(CrisprSystemType.CasX, "CasX", "TTCN")]               // Liu 2019
    public void GetSystem_ReturnsCorrectNameAndPam(CrisprSystemType systemType, string expectedName, string expectedPam)
    {
        var system = CrisprDesigner.GetSystem(systemType);

        Assert.Multiple(() =>
        {
            Assert.That(system.Name, Is.EqualTo(expectedName));
            Assert.That(system.PamSequence, Is.EqualTo(expectedPam));
        });
    }

    [Test]
    [Description("M4: Each system has correct guide length per literature")]
    [TestCase(CrisprSystemType.SpCas9, 20)]      // Jinek 2012: 20 nt spacer
    [TestCase(CrisprSystemType.SpCas9_NAG, 20)]   // Same enzyme as SpCas9 (Hsu 2013)
    [TestCase(CrisprSystemType.SaCas9, 21)]        // Ran et al. 2015: 21 nt typical
    [TestCase(CrisprSystemType.Cas12a, 23)]        // Zetsche 2015: ~23 nt crRNA spacer
    [TestCase(CrisprSystemType.AsCas12a, 23)]      // Acidaminococcus sp. Cas12a (Zetsche 2015)
    [TestCase(CrisprSystemType.LbCas12a, 24)]      // Lachnospiraceae Cas12a (Zetsche 2015: 23-25 nt)
    [TestCase(CrisprSystemType.CasX, 20)]          // Liu et al. 2019: 20 nt guide
    public void GetSystem_ReturnsCorrectGuideLength(CrisprSystemType systemType, int expectedLength)
    {
        var system = CrisprDesigner.GetSystem(systemType);
        Assert.That(system.GuideLength, Is.EqualTo(expectedLength));
    }

    [Test]
    [Description("M5: Cas9 variants have PAM after target (3'), Cas12a/CasX have PAM before target (5')")]
    [TestCase(CrisprSystemType.SpCas9, true)]      // Wikipedia PAM: 3'-end PAM
    [TestCase(CrisprSystemType.SpCas9_NAG, true)]   // Same as SpCas9
    [TestCase(CrisprSystemType.SaCas9, true)]       // Ran 2015: 3' PAM
    [TestCase(CrisprSystemType.Cas12a, false)]      // Zetsche 2015: 5' T-rich PAM
    [TestCase(CrisprSystemType.AsCas12a, false)]    // Cas12a family
    [TestCase(CrisprSystemType.LbCas12a, false)]    // Cas12a family
    [TestCase(CrisprSystemType.CasX, false)]        // Liu 2019: 5' PAM
    public void GetSystem_ReturnsCorrectPamPosition(CrisprSystemType systemType, bool pamAfterTarget)
    {
        var system = CrisprDesigner.GetSystem(systemType);
        Assert.That(system.PamAfterTarget, Is.EqualTo(pamAfterTarget));
    }

    [Test]
    [Description("S4: GetSystem throws ArgumentException for unknown system type")]
    public void GetSystem_UnknownType_ThrowsArgumentException()
    {
        var invalidType = (CrisprSystemType)999;
        Assert.Throws<ArgumentException>(() => CrisprDesigner.GetSystem(invalidType));
    }

    #endregion

    #region FindPamSites - SpCas9 NGG Detection

    [Test]
    [Description("M6: FindPamSites detects NGG on forward strand with correct position")]
    public void FindPamSites_SpCas9_DetectsNGG_OnForwardStrand()
    {
        // 20bp target followed by AGG PAM at position 20
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAGG");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        var site = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == "AGG");
        Assert.That(site, Is.Not.Null, "Should find AGG PAM on forward strand");
        Assert.That(site!.Position, Is.EqualTo(20), "PAM should be at position 20");
        Assert.That(site.TargetStart, Is.EqualTo(0), "Target should start at position 0");
    }

    [Test]
    [Description("M7: NGG matches all variants - AGG, CGG, TGG, GGG (N = any nucleotide, IUPAC)")]
    [TestCase("ACGTACGTACGTACGTACGTAGG", "AGG")]
    [TestCase("ACGTACGTACGTACGTACGTCGG", "CGG")]
    [TestCase("ACGTACGTACGTACGTACGTTGG", "TGG")]
    [TestCase("ACGTACGTACGTACGTACGTGGG", "GGG")]
    public void FindPamSites_SpCas9_MatchesAllNGG_Variants(string sequenceStr, string expectedPam)
    {
        var sequence = new DnaSequence(sequenceStr);
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        var matchingSite = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == expectedPam);
        Assert.That(matchingSite, Is.Not.Null,
            $"Should find {expectedPam} PAM (N = any nucleotide per IUPAC)");
        Assert.That(matchingSite!.Position, Is.EqualTo(20),
            $"{expectedPam} PAM should be at position 20 (after 20bp target)");
    }

    [Test]
    [Description("M8: FindPamSites searches reverse strand and returns correct PAM/coordinates")]
    public void FindPamSites_SpCas9_SearchesReverseStrand()
    {
        // Forward: CCA at start, no GG anywhere → no forward hits
        // RevComp ends in ...TGG → NGG match on reverse strand
        // ForwardPos = len - revIdx - pamLen = 0
        // PamSequence stored as revcomp of TGG = CCA
        var sequence = new DnaSequence("CCAACGTACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        var reverseSites = sites.Where(s => !s.IsForwardStrand).ToList();
        Assert.That(reverseSites, Has.Count.EqualTo(1),
            "Should find exactly 1 PAM site on reverse strand");

        var site = reverseSites[0];
        Assert.Multiple(() =>
        {
            Assert.That(site.Position, Is.EqualTo(0),
                "CCA at forward pos 0 → reverse strand PAM maps to forward position 0");
            Assert.That(site.PamSequence, Is.EqualTo("CCA"),
                "PamSequence is reverse complement of TGG back to forward strand");
            Assert.That(site.TargetSequence.Length, Is.EqualTo(20),
                "Reverse strand target should be 20bp");
        });
    }

    [Test]
    [Description("S2: Reverse strand positions are correctly converted to forward-strand coordinates")]
    public void FindPamSites_SpCas9_ReverseStrand_PositionConversion()
    {
        // Build a 23-char sequence: 20bp + NGG at end = forward strand hit
        // Also build: CCN at start + 20bp = reverse strand hit
        // Forward: "CCA" + 20bp = CCAbcdefghijklmnopqrst (23bp)
        // RevComp: will have TGG at the end → NGG match
        // RevComp PAM at index 20 on revcomp → forwardPos = 23 - 20 - 3 = 0
        var sequence = new DnaSequence("CCAACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        var reverseSite = sites.FirstOrDefault(s => !s.IsForwardStrand);
        Assert.That(reverseSite, Is.Not.Null, "Should find reverse strand site");
        Assert.That(reverseSite!.Position, Is.EqualTo(0),
            "Reverse strand PAM at revcomp end should map to forward position 0");
    }

    [Test]
    [Description("M11: FindPamSites handles lowercase input (case-insensitive)")]
    public void FindPamSites_SpCas9_CaseInsensitive()
    {
        var sites = CrisprDesigner.FindPamSites("acgtacgtacgtacgtacgtagg", CrisprSystemType.SpCas9).ToList();

        Assert.That(sites.Any(s => s.PamSequence == "AGG"), Is.True,
            "Should find PAM in lowercase sequence");
    }

    [Test]
    [Description("M12: FindPamSites returns target sequence with correct length and content")]
    public void FindPamSites_SpCas9_ReturnsTargetSequence_WithCorrectContent()
    {
        // Sequence: exactly 20bp target "ACGTACGTACGTACGTACGT" + "AGG" PAM
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAGG");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        var site = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == "AGG");
        Assert.That(site, Is.Not.Null);
        Assert.That(site!.TargetSequence.Length, Is.EqualTo(20),
            "SpCas9 target should be 20bp (Jinek 2012)");
        Assert.That(site.TargetSequence, Is.EqualTo("ACGTACGTACGTACGTACGT"),
            "Target content should match the 20bp preceding the PAM");
    }

    #endregion

    #region FindPamSites - Edge Cases

    [Test]
    [Description("M9: FindPamSites returns no sites on either strand for PAM-free sequence")]
    public void FindPamSites_NoPamPresent_ReturnsEmpty()
    {
        // AC-repeat: no GG on forward, no GG on reverse (revcomp of ACAC... = GTGT... also no GG)
        var sequence = new DnaSequence("ACACACACACACACACACACACACACACACAC");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        Assert.That(sites, Is.Empty,
            "Should find no PAM sites on either strand in a sequence with no GG dinucleotide");
    }

    [Test]
    [Description("M10: FindPamSites returns empty for empty sequence")]
    public void FindPamSites_EmptySequence_ReturnsEmpty()
    {
        var sites = CrisprDesigner.FindPamSites("", CrisprSystemType.SpCas9).ToList();
        Assert.That(sites, Is.Empty);
    }

    [Test]
    [Description("S1: FindPamSites excludes forward strand sites where target would be out of bounds")]
    public void FindPamSites_TargetOutOfBounds_Excluded()
    {
        // PAM at position 0-2 means target would need to start at -20 (invalid)
        var sequence = new DnaSequence("AGGACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        var site = sites.FirstOrDefault(s => s.IsForwardStrand && s.Position == 0);
        Assert.That(site, Is.Null,
            "Should not return PAM site when target would be out of bounds");
    }

    [Test]
    [Description("S1: Cas12a excludes sites where target extends past sequence end")]
    public void FindPamSites_Cas12a_TargetOutOfBounds_Excluded()
    {
        // TTTA PAM at position 0, but only 10bp follow — need 23bp for target
        var sequence = new DnaSequence("TTTAACGTAC");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.Cas12a).ToList();

        var forwardSite = sites.FirstOrDefault(s => s.IsForwardStrand);
        Assert.That(forwardSite, Is.Null,
            "Cas12a should exclude PAM site when target would extend past sequence end");
    }

    [Test]
    [Description("S5: FindPamSites with null DnaSequence throws ArgumentNullException")]
    public void FindPamSites_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.FindPamSites((DnaSequence)null!, CrisprSystemType.SpCas9).ToList());
    }

    #endregion

    #region FindPamSites - Cas12a TTTV Detection

    [Test]
    [Description("M13: Cas12a detects TTTA with correct PAM position, target start, and content (Zetsche 2015)")]
    public void FindPamSites_Cas12a_DetectsTTTA()
    {
        // TTTA PAM (4bp) followed by 23bp target = 27bp minimum needed
        var sequence = new DnaSequence("TTTAACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.Cas12a).ToList();

        var site = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == "TTTA");
        Assert.That(site, Is.Not.Null, "Should find TTTA PAM (V includes A)");
        Assert.Multiple(() =>
        {
            Assert.That(site!.Position, Is.EqualTo(0), "TTTA PAM should be at position 0");
            Assert.That(site.TargetStart, Is.EqualTo(4),
                "Cas12a target starts after 4bp PAM (5' PAM before target)");
            Assert.That(site.TargetSequence, Is.EqualTo("ACGTACGTACGTACGTACGTACG"),
                "Target content should be the 23bp following the PAM");
        });
    }

    [Test]
    [Description("M13: Cas12a detects TTTC variant (V = C per IUPAC)")]
    public void FindPamSites_Cas12a_DetectsTTTC()
    {
        var sequence = new DnaSequence("TTTCACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.Cas12a).ToList();

        var site = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == "TTTC");
        Assert.That(site, Is.Not.Null, "Should find TTTC PAM (V includes C)");
        Assert.That(site!.Position, Is.EqualTo(0), "TTTC PAM should be at position 0");
    }

    [Test]
    [Description("M13: Cas12a detects TTTG variant (V = G per IUPAC)")]
    public void FindPamSites_Cas12a_DetectsTTTG()
    {
        var sequence = new DnaSequence("TTTGACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.Cas12a).ToList();

        var site = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == "TTTG");
        Assert.That(site, Is.Not.Null, "Should find TTTG PAM (V includes G)");
        Assert.That(site!.Position, Is.EqualTo(0), "TTTG PAM should be at position 0");
    }

    [Test]
    [Description("M14: Cas12a does NOT detect TTTT (V excludes T per IUPAC)")]
    public void FindPamSites_Cas12a_DoesNotDetectTTTT()
    {
        var sequence = new DnaSequence("TTTTACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.Cas12a).ToList();

        Assert.That(sites.Any(s => s.IsForwardStrand && s.PamSequence == "TTTT"), Is.False,
            "Should NOT find TTTT PAM (V excludes T per IUPAC)");
    }

    [Test]
    [Description("S6: AsCas12a FindPamSites produces 23bp target (Zetsche 2015)")]
    public void FindPamSites_AsCas12a_DetectsTTTV_With23bpTarget()
    {
        // Same sequence as Cas12a tests: TTTA(4bp) + 24bp available
        // AsCas12a guide=23 → target is first 23bp after PAM
        var sequence = new DnaSequence("TTTAACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.AsCas12a).ToList();

        var site = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == "TTTA");
        Assert.That(site, Is.Not.Null, "AsCas12a should detect TTTA PAM");
        Assert.Multiple(() =>
        {
            Assert.That(site!.Position, Is.EqualTo(0));
            Assert.That(site.TargetStart, Is.EqualTo(4));
            Assert.That(site.TargetSequence.Length, Is.EqualTo(23),
                "AsCas12a guide is 23 nt (Zetsche 2015)");
            Assert.That(site.TargetSequence, Is.EqualTo("ACGTACGTACGTACGTACGTACG"));
        });
    }

    [Test]
    [Description("S7: LbCas12a FindPamSites produces 24bp target (Zetsche 2015: 23-25 nt)")]
    public void FindPamSites_LbCas12a_DetectsTTTV_With24bpTarget()
    {
        // TTTA(4bp) + 24bp target = 28bp total
        var sequence = new DnaSequence("TTTAACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.LbCas12a).ToList();

        var site = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == "TTTA");
        Assert.That(site, Is.Not.Null, "LbCas12a should detect TTTA PAM");
        Assert.Multiple(() =>
        {
            Assert.That(site!.Position, Is.EqualTo(0));
            Assert.That(site.TargetStart, Is.EqualTo(4));
            Assert.That(site.TargetSequence.Length, Is.EqualTo(24),
                "LbCas12a guide is 24 nt (Zetsche 2015)");
            Assert.That(site.TargetSequence, Is.EqualTo("ACGTACGTACGTACGTACGTACGT"));
        });
    }

    #endregion

    #region FindPamSites - SaCas9 NNGRRT Detection

    [Test]
    [Description("M15: SaCas9 detects NNGRRT with both R=A — PAM AAGAAT (Ran et al. 2015)")]
    public void FindPamSites_SaCas9_DetectsNNGAAT()
    {
        // 21bp target + 6bp PAM AAGAAT = 27 chars
        // Sequence: ACGTACGTACGTACGTACGTA (21) + AAGAAT (6)
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAAAGAAT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SaCas9).ToList();

        var matchingSite = sites.FirstOrDefault(s => s.IsForwardStrand);
        Assert.That(matchingSite, Is.Not.Null,
            "Should find NNGRRT PAM with R=A");
        Assert.That(matchingSite!.PamSequence, Is.EqualTo("AAGAAT"),
            "PAM should be AAGAAT matching NNGRRT pattern (R=A at both R positions)");
    }

    [Test]
    [Description("M15: SaCas9 detects NNGRRT with R=A,G mixed — PAM AGGAGT (Ran et al. 2015)")]
    public void FindPamSites_SaCas9_DetectsNNGAGT()
    {
        // 21bp target + PAM AGGAGT: NN=AG, G=G, R=A, R=G... wait
        // NNGRRT: N=A, N=G, G=G, R=A, R=G, T=T → AGGAGT
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAAGGAGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SaCas9).ToList();

        var matchingSite = sites.FirstOrDefault(s => s.IsForwardStrand);
        Assert.That(matchingSite, Is.Not.Null,
            "Should find NNGRRT PAM with mixed R values");
        Assert.That(matchingSite!.PamSequence, Is.EqualTo("AGGAGT"),
            "PAM should be AGGAGT matching NNGRRT (R=A at pos 3, R=G at pos 4)");
    }

    [Test]
    [Description("M15: SaCas9 detects NNGRRT with both R=G — PAM AAGGHT (Ran et al. 2015)")]
    public void FindPamSites_SaCas9_DetectsNNGGGT()
    {
        // 21bp target + PAM where R=G at both positions: AAGGHT → AAGGGT
        // NNGRRT: N=A, N=A, G=G, R=G, R=G, T=T → AAGGGT
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAAAGGGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SaCas9).ToList();

        var matchingSite = sites.FirstOrDefault(s => s.IsForwardStrand);
        Assert.That(matchingSite, Is.Not.Null,
            "Should find NNGRRT PAM with R=G at both positions");
        Assert.That(matchingSite!.PamSequence, Is.EqualTo("AAGGGT"),
            "PAM should be AAGGGT matching NNGRRT (R=G,G)");
    }

    [Test]
    [Description("SaCas9 returns 21bp target sequence (Ran et al. 2015: 21 nt guide)")]
    public void FindPamSites_SaCas9_Returns21bp_TargetSequence()
    {
        // 21bp target + 6bp PAM = 27 chars
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAAAGAAT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SaCas9).ToList();

        var site = sites.FirstOrDefault(s => s.IsForwardStrand);
        Assert.That(site, Is.Not.Null, "Must find a SaCas9 PAM site");
        Assert.That(site!.TargetSequence.Length, Is.EqualTo(21),
            "SaCas9 target should be 21bp (Ran et al. 2015)");
        Assert.That(site.TargetSequence, Is.EqualTo("ACGTACGTACGTACGTACGTA"),
            "Target content should be the 21bp preceding the PAM");
    }

    [Test]
    [Description("SaCas9 rejects non-NNGRRT pattern — NNGCCT has C at R position")]
    public void FindPamSites_SaCas9_RejectsNonR_AtRPosition()
    {
        // PAM AAGCCT: N=A, N=A, G=G, R=C(invalid!), R=C(invalid!), T=T
        // C is not A or G, so this should NOT match NNGRRT
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAAAGCCT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SaCas9).ToList();

        var forwardSite = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == "AAGCCT");
        Assert.That(forwardSite, Is.Null,
            "Should NOT match NNGCCT — C is not R (R = A or G per IUPAC)");
    }

    #endregion

    #region FindPamSites - SpCas9_NAG Detection

    [Test]
    [Description("SpCas9_NAG detects NAG on forward strand (Hsu 2013, Wikipedia Cas9)")]
    public void FindPamSites_SpCas9NAG_DetectsNAG()
    {
        // 20bp target + AAG PAM at position 20
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAAG");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9_NAG).ToList();

        var site = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == "AAG");
        Assert.That(site, Is.Not.Null, "Should find AAG (N=A) matching NAG pattern");
        Assert.That(site!.Position, Is.EqualTo(20), "NAG PAM should be at position 20");
        Assert.That(site.TargetSequence.Length, Is.EqualTo(20),
            "SpCas9_NAG uses same 20 nt guide as SpCas9");
    }

    [Test]
    [Description("SpCas9_NAG matches all NAG variants — AAG, CAG, TAG, GAG (N = any nucleotide)")]
    [TestCase("ACGTACGTACGTACGTACGTAAG", "AAG")]
    [TestCase("ACGTACGTACGTACGTACGTCAG", "CAG")]
    [TestCase("ACGTACGTACGTACGTACGTTAG", "TAG")]
    [TestCase("ACGTACGTACGTACGTACGTGAG", "GAG")]
    public void FindPamSites_SpCas9NAG_MatchesAllNAG_Variants(string sequenceStr, string expectedPam)
    {
        var sequence = new DnaSequence(sequenceStr);
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9_NAG).ToList();

        var matchingSite = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == expectedPam);
        Assert.That(matchingSite, Is.Not.Null,
            $"Should find {expectedPam} PAM matching NAG pattern (N = any nucleotide per IUPAC)");
        Assert.That(matchingSite!.Position, Is.EqualTo(20),
            $"{expectedPam} PAM should be at position 20 (after 20bp target)");
    }

    #endregion

    #region FindPamSites - CasX TTCN Detection

    [Test]
    [Description("CasX detects TTCN PAM with all N variants (Liu et al. 2019)")]
    [TestCase("TTCA", "A")]
    [TestCase("TTCC", "C")]
    [TestCase("TTCG", "G")]
    [TestCase("TTCT", "T")]
    public void FindPamSites_CasX_DetectsTTCN_Variants(string pamStr, string nVariant)
    {
        // TTCN PAM (4bp) + 20bp target = 24bp
        var sequence = new DnaSequence(pamStr + "ACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.CasX).ToList();

        var site = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == pamStr);
        Assert.That(site, Is.Not.Null,
            $"Should find {pamStr} PAM (N={nVariant} per IUPAC)");
        Assert.That(site!.Position, Is.EqualTo(0), "CasX PAM before target at position 0");
        Assert.That(site.TargetStart, Is.EqualTo(4),
            "CasX target starts after 4bp PAM");
        Assert.That(site.TargetSequence.Length, Is.EqualTo(20),
            "CasX guide is 20 nt (Liu et al. 2019)");
    }

    #endregion

    #region FindPamSites - Overlapping PAM Sites

    [Test]
    [Description("S3: Overlapping PAM sites are all reported — each PAM = unique guide RNA (Wikipedia PAM)")]
    public void FindPamSites_OverlappingPamSites_AllReported()
    {
        // "AGGTGG" contains NGG at position 20 (AGG) and position 23 (TGG)
        // Need 20bp target before each PAM for SpCas9
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAGGTGG");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        var forwardSites = sites.Where(s => s.IsForwardStrand).OrderBy(s => s.Position).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(forwardSites, Has.Count.EqualTo(2),
                "Should find exactly 2 overlapping NGG sites on forward strand");
            Assert.That(forwardSites[0].PamSequence, Is.EqualTo("AGG"));
            Assert.That(forwardSites[0].Position, Is.EqualTo(20),
                "AGG PAM at first overlapping position");
            Assert.That(forwardSites[1].PamSequence, Is.EqualTo("TGG"));
            Assert.That(forwardSites[1].Position, Is.EqualTo(23),
                "TGG PAM at second overlapping position");
        });
    }

    #endregion

    #region String Overload Tests

    [Test]
    [Description("String overload produces identical results to DnaSequence overload")]
    public void FindPamSites_StringOverload_WorksIdentically()
    {
        const string sequenceStr = "ACGTACGTACGTACGTACGTAGG";
        var dnaSequence = new DnaSequence(sequenceStr);

        var stringResults = CrisprDesigner.FindPamSites(sequenceStr, CrisprSystemType.SpCas9).ToList();
        var dnaResults = CrisprDesigner.FindPamSites(dnaSequence, CrisprSystemType.SpCas9).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(stringResults.Count, Is.EqualTo(dnaResults.Count),
                "String and DnaSequence overloads should return same number of sites");

            for (int i = 0; i < stringResults.Count; i++)
            {
                Assert.That(stringResults[i].Position, Is.EqualTo(dnaResults[i].Position),
                    $"Site {i} position should match between overloads");
                Assert.That(stringResults[i].PamSequence, Is.EqualTo(dnaResults[i].PamSequence),
                    $"Site {i} PAM sequence should match between overloads");
                Assert.That(stringResults[i].IsForwardStrand, Is.EqualTo(dnaResults[i].IsForwardStrand),
                    $"Site {i} strand should match between overloads");
            }
        });
    }

    #endregion
}
