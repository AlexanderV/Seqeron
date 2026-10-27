namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// All nine CRISPR systems, member by member (audit round 4, B08 H38, H48, H50, H51): the system table including
/// <see cref="CrisprSystem.Description"/> and <see cref="CrisprSystem.CutOffset"/>, <see cref="CrisprDesigner.GetCutSite"/>
/// and <see cref="CrisprDesigner.DesignGuideRnas"/> for every system, and the bulge scan's <see cref="BulgeType.None"/>
/// convention. Before this class six systems (CasX, SaCas9, SaCas9-NNGRRN, SpCas9-NAG, SpCas9-NGA, AsCas12a) had no
/// <c>GetCutSite</c> / <c>DesignGuideRnas</c> test at all.
/// TestSpec: TestSpecs/CRISPR-PAM-001.md, CRISPR-GUIDE-001.md.
/// </summary>
[TestFixture]
public class CrisprDesigner_Systems_Tests
{
    // type, display name, PAM, guide length, PAM after target, cut offset, description (pinned verbatim)
    public static IEnumerable<TestCaseData> Systems()
    {
        yield return new TestCaseData(CrisprSystemType.SpCas9, "SpCas9", "NGG", 20, true, 3, "Streptococcus pyogenes Cas9").SetName("SpCas9");
        yield return new TestCaseData(CrisprSystemType.SpCas9_NAG, "SpCas9-NAG", "NAG", 20, true, 3, "SpCas9 with NAG PAM (lower activity)").SetName("SpCas9_NAG");
        yield return new TestCaseData(CrisprSystemType.SpCas9_NGA, "SpCas9-NGA", "NGA", 20, true, 3, "SpCas9 with NGA PAM (alternative off-target PAM of NGG, Kleinstiver 2015)").SetName("SpCas9_NGA");
        yield return new TestCaseData(CrisprSystemType.SaCas9, "SaCas9", "NNGRRT", 21, true, 3, "Staphylococcus aureus Cas9").SetName("SaCas9");
        yield return new TestCaseData(CrisprSystemType.SaCas9_NNGRRN, "SaCas9-NNGRRN", "NNGRRN", 21, true, 3, "SaCas9 with relaxed NNGRRN PAM (alternative off-target PAM of NNGRRT)").SetName("SaCas9_NNGRRN");
        yield return new TestCaseData(CrisprSystemType.Cas12a, "Cas12a/Cpf1", "TTTV", 23, false, 18, "Cas12a (Cpf1) - PAM before target").SetName("Cas12a");
        yield return new TestCaseData(CrisprSystemType.AsCas12a, "AsCas12a", "TTTV", 23, false, 18, "Acidaminococcus sp. Cas12a").SetName("AsCas12a");
        yield return new TestCaseData(CrisprSystemType.LbCas12a, "LbCas12a", "TTTV", 23, false, 18, "Lachnospiraceae bacterium Cas12a").SetName("LbCas12a");
        yield return new TestCaseData(CrisprSystemType.CasX, "CasX", "TTCN", 20, false, 13, "CasX - compact Cas protein").SetName("CasX");
    }

    /// <summary>H50: <c>Description</c> (a declared output field of the crispr_system_info tool) and the cut offset are pinned.</summary>
    [TestCaseSource(nameof(Systems))]
    public void GetSystem_EveryMember_IsPinned(
        CrisprSystemType type, string name, string pam, int guideLength, bool pamAfter, int cutOffset, string description)
    {
        var system = CrisprDesigner.GetSystem(type);

        Assert.Multiple(() =>
        {
            Assert.That(system.Name, Is.EqualTo(name));
            Assert.That(system.PamSequence, Is.EqualTo(pam));
            Assert.That(system.GuideLength, Is.EqualTo(guideLength));
            Assert.That(system.PamAfterTarget, Is.EqualTo(pamAfter));
            Assert.That(system.CutOffset, Is.EqualTo(cutOffset));
            Assert.That(system.Description, Is.EqualTo(description));
        });
    }

    [Test]
    public void GetSystem_CoversEveryEnumValue_AndRejectsAnUnknownOne()
    {
        Assert.That(Systems().Select(c => (CrisprSystemType)c.Arguments[0]!),
            Is.EquivalentTo(Enum.GetValues<CrisprSystemType>()), "a new system must be added to this table");
        var e = Assert.Throws<ArgumentException>(() => CrisprDesigner.GetSystem((CrisprSystemType)999));
        Assert.That(e!.Message, Does.Contain("Unknown CRISPR system"));
    }

    /// <summary>
    /// H51 / H38: the cut site of a PAM at forward position 40 on each strand, from the documented conventions:
    /// Cas9-type <c>Position - 3</c> / <c>Position + pamLength + 2</c>; Cas12a <c>Position + pamLength + 18</c> /
    /// <c>Position - 19</c> (Zetsche 2015 as described by CRISPOR); CasX <c>Position + pamLength + 13</c> /
    /// <c>Position - 14</c> (Liu 2019, midpoint of the non-target-strand 12-14).
    /// </summary>
    [TestCaseSource(nameof(Systems))]
    public void GetCutSite_EveryMember_FollowsItsConvention_OnBothStrands(
        CrisprSystemType type, string name, string pam, int guideLength, bool pamAfter, int cutOffset, string description)
    {
        var system = CrisprDesigner.GetSystem(type);
        string concretePam = pam.Replace('N', 'A').Replace('V', 'A').Replace('R', 'A');
        var forward = new PamSite(40, concretePam, new string('A', guideLength), 40 - guideLength, true, system);
        var reverse = new PamSite(40, concretePam, new string('A', guideLength), 40 + pam.Length, false, system);

        int expectedForward = pamAfter ? 37 : 40 + pam.Length + (type == CrisprSystemType.CasX ? 13 : 18);
        int expectedReverse = pamAfter ? 40 + pam.Length + 2 : 40 - (type == CrisprSystemType.CasX ? 14 : 19);
        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.GetCutSite(forward), Is.EqualTo(expectedForward), name + " forward");
            Assert.That(CrisprDesigner.GetCutSite(reverse), Is.EqualTo(expectedReverse), name + " reverse");
        });
    }

    /// <summary>
    /// H51 / H38: <c>DesignGuideRnas</c> decides inclusion by the cut site, so a region that ends one base before the cut
    /// must not return the guide and one that ends on it must (this also pins CasX's 13 against the former Cas12a 18).
    /// Each system gets a forward site with a GC-balanced 50 % guide and a flank on both sides.
    /// </summary>
    [TestCaseSource(nameof(Systems))]
    public void DesignGuideRnas_EveryMember_IncludesAGuideExactlyWhenItsCutSiteIsInTheRegion(
        CrisprSystemType type, string name, string pam, int guideLength, bool pamAfter, int cutOffset, string description)
    {
        var system = CrisprDesigner.GetSystem(type);
        string guide = string.Concat(Enumerable.Repeat("ACGGATCTAG", 3))[..guideLength];
        string concretePam = pam.Replace('N', 'C').Replace('V', 'C').Replace('R', 'G');
        string text = "TATATATACCATA" + (pamAfter ? guide + concretePam : concretePam + guide) + "ATTAGATCACAT";
        // keep the sequence free of an accidental PAM on either strand that would add candidates
        var sites = CrisprDesigner.FindPamSites(text, type).ToList();
        var site = sites.Single(s => s.IsForwardStrand && s.TargetSequence == guide);
        int cut = CrisprDesigner.GetCutSite(site);
        var open = GuideRnaParameters.Default with
        { MinGcContent = 0, MaxGcContent = 100, MinScore = 0, ScoreGcContent = false };

        var atCut = CrisprDesigner.DesignGuideRnas(new DnaSequence(text), cut, cut, type, open).Where(c => c.Sequence == guide).ToList();
        var oneBefore = CrisprDesigner.DesignGuideRnas(new DnaSequence(text), 0, cut - 1, type, open).Where(c => c.Sequence == guide).ToList();
        var oneAfter = CrisprDesigner.DesignGuideRnas(new DnaSequence(text), cut + 1, text.Length - 1, type, open).Where(c => c.Sequence == guide).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(atCut, Has.Count.EqualTo(1), name + ": region = the cut base");
            Assert.That(atCut[0].System, Is.EqualTo(system));
            Assert.That(atCut[0].Sequence, Has.Length.EqualTo(guideLength));
            Assert.That(atCut[0].IsForwardStrand, Is.True);
            Assert.That(atCut[0].Position, Is.EqualTo(site.TargetStart));
            Assert.That(oneBefore, Is.Empty, name + ": region ends one base before the cut");
            Assert.That(oneAfter, Is.Empty, name + ": region starts one base after the cut");
        });
    }

    /// <summary>
    /// H48: the documented convention of the bulge scan's <see cref="BulgeType.None"/> rows - no gap, so
    /// <c>BulgeSize == 0</c> and <c>BulgePosition == -1</c> - on the simplest live case: one substitution, no bulges allowed.
    /// </summary>
    [Test]
    public void FindOffTargetsWithBulges_UngappedHit_HasBulgeNone_SizeZero_PositionMinusOne()
    {
        const string guide = "GGCCGACCTGTCGCTGACGC";
        const string site = "GGCCGACCTGTCGCAGACGC"; // one substitution at index 14
        var genome = new DnaSequence("TTTTTTTTTT" + site + "AGG" + "TTTTTTTTTT");

        var hits = CrisprDesigner.FindOffTargetsWithBulges(guide, genome, 1, 0, 0).ToList();

        Assert.That(hits, Has.Count.EqualTo(1));
        var hit = hits[0];
        Assert.Multiple(() =>
        {
            Assert.That(hit.Bulge, Is.EqualTo(BulgeType.None));
            Assert.That(hit.BulgeSize, Is.Zero);
            Assert.That(hit.BulgePosition, Is.EqualTo(-1));
            Assert.That(hit.Mismatches, Is.EqualTo(1));
            Assert.That(hit.AlignedGuide, Is.EqualTo(guide));
            Assert.That(hit.AlignedTarget, Is.EqualTo(site));
            Assert.That(hit.Pam, Is.EqualTo("AGG"));
            Assert.That(hit.Position, Is.EqualTo(30));
            Assert.That(hit.TargetStart, Is.EqualTo(10));
        });
    }
}
