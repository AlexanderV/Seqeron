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
        });
    }
}
