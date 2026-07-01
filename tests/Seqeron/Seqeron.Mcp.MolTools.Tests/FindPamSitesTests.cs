using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class FindPamSitesTests
{
    // 20-nt guide (no internal GG) + TGG PAM at index 20. SpCas9 = NGG, guide 20, PAM after target.
    private const string Guide = "ATATATATATATATATATAT";
    private const string Sequence = "ATATATATATATATATATATTGG";

    [Test]
    public void FindPamSites_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.find_pam_sites(Sequence));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_pam_sites(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_pam_sites(null!));
    }

    [Test]
    public void FindPamSites_Binding_InvokesSuccessfully()
    {
        var sites = MolToolsTools.find_pam_sites(Sequence, CrisprSystemType.SpCas9).Sites;

        // There must be a forward-strand PAM "TGG" at position 20 whose 20-nt target is the guide.
        var fwd = sites.Where(s => s.IsForwardStrand && s.Position == 20).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(fwd, Has.Count.EqualTo(1));
            Assert.That(fwd[0].PamSequence, Is.EqualTo("TGG"));
            Assert.That(fwd[0].TargetSequence, Is.EqualTo(Guide));
            Assert.That(fwd[0].TargetStart, Is.EqualTo(0));
        });

        // A sequence too short to contain any full guide+PAM window yields no sites.
        Assert.That(MolToolsTools.find_pam_sites("TGG", CrisprSystemType.SpCas9).Sites, Is.Empty);
    }
}
