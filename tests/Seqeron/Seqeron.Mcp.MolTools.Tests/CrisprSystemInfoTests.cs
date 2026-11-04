using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class CrisprSystemInfoTests
{
    [Test]
    public void CrisprSystemInfo_Schema_ValidatesCorrectly()
    {
        // Enum inputs are always valid; the tool must return a non-null record for each.
        Assert.DoesNotThrow(() => MolToolsTools.crispr_system_info(CrisprSystemType.SpCas9));
        Assert.That(MolToolsTools.crispr_system_info(CrisprSystemType.Cas12a), Is.Not.Null);
    }

    [Test]
    public void CrisprSystemInfo_Binding_InvokesSuccessfully()
    {
        var sp = MolToolsTools.crispr_system_info(CrisprSystemType.SpCas9);
        var sa = MolToolsTools.crispr_system_info(CrisprSystemType.SaCas9);
        var cas12a = MolToolsTools.crispr_system_info(CrisprSystemType.Cas12a);

        Assert.Multiple(() =>
        {
            // SpCas9: NGG PAM, 20-nt guide, PAM downstream of target.
            Assert.That(sp.Name, Is.EqualTo("SpCas9"));
            Assert.That(sp.PamSequence, Is.EqualTo("NGG"));
            Assert.That(sp.GuideLength, Is.EqualTo(20));
            Assert.That(sp.PamAfterTarget, Is.True);

            // SaCas9: NNGRRT PAM, 21-nt guide.
            Assert.That(sa.PamSequence, Is.EqualTo("NNGRRT"));
            Assert.That(sa.GuideLength, Is.EqualTo(21));

            // Cas12a: TTTV PAM, 23-nt guide, PAM upstream of target.
            Assert.That(cas12a.PamSequence, Is.EqualTo("TTTV"));
            Assert.That(cas12a.GuideLength, Is.EqualTo(23));
            Assert.That(cas12a.PamAfterTarget, Is.False);

            // cutOffset: 3 bp 5' of the PAM (Cas9 type), after the 18th base (Cas12a), 13 nt after the PAM (CasX, Liu 2019).
            Assert.That(sp.CutOffset, Is.EqualTo(3));
            Assert.That(sa.CutOffset, Is.EqualTo(3));
            Assert.That(cas12a.CutOffset, Is.EqualTo(18));
            Assert.That(MolToolsTools.crispr_system_info(CrisprSystemType.CasX).CutOffset, Is.EqualTo(13));
        });
    }

    [TestCase(CrisprSystemType.SpCas9, "Streptococcus pyogenes Cas9")]
    [TestCase(CrisprSystemType.SpCas9_NGA, "SpCas9 with NGA PAM (CRISPOR's alternative off-target PAM of NGG after Tsai 2014 GUIDE-seq; the PAM of the VQR variant, Kleinstiver 2015)")]
    [TestCase(CrisprSystemType.CasX, "CasX - compact Cas protein")]
    public void CrisprSystemInfo_Description_IsReturned(CrisprSystemType type, string description)
    {
        Assert.That(MolToolsTools.crispr_system_info(type).Description, Is.EqualTo(description));
    }

    /// <summary>The documented error (docs/mcp/tools/moltools/crispr_system_info.md): an enum value outside the nine systems.</summary>
    [Test]
    public void CrisprSystemInfo_UnknownSystem_ThrowsArgumentException()
    {
        var e = Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_system_info((CrisprSystemType)999));
        Assert.That(e!.Message, Does.Contain("Unknown CRISPR system"));
    }
}
