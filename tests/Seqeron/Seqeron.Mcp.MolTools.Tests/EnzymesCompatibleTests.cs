using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class EnzymesCompatibleTests
{
    [Test]
    public void EnzymesCompatible_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.enzymes_compatible("EcoRV", "SmaI"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.enzymes_compatible("", "SmaI"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.enzymes_compatible("EcoRV", null!));

        // RESTR-FILTER-001 / F31: an unknown enzyme name is now an ERROR, not `false`.
        // A typo must not be reported as "not compatible" — that is the answer that silently
        // loses a real cloning option. find_restriction_sites already errors the same way.
        var ex = Assert.Throws<ArgumentException>(
            () => MolToolsTools.enzymes_compatible("NotAnEnzyme", "EcoRI"));
        Assert.That(ex!.Message, Does.Contain("Unknown enzyme: NotAnEnzyme"));
    }

    [Test]
    public void EnzymesCompatible_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Two blunt cutters -> compatible.
            Assert.That(MolToolsTools.enzymes_compatible("EcoRV", "SmaI").Compatible, Is.True);

            // BamHI and BglII share a 5' GATC overhang -> compatible.
            Assert.That(MolToolsTools.enzymes_compatible("BamHI", "BglII").Compatible, Is.True);

            // EcoRI (AATT overhang) vs BamHI (GATC overhang) -> not compatible.
            Assert.That(MolToolsTools.enzymes_compatible("EcoRI", "BamHI").Compatible, Is.False);

            // compatibleEnd echoes the shared end, and is null when incompatible.
            Assert.That(MolToolsTools.enzymes_compatible("EcoRV", "SmaI").CompatibleEnd, Is.EqualTo("blunt"));
            Assert.That(MolToolsTools.enzymes_compatible("BamHI", "BglII").CompatibleEnd, Is.EqualTo("GATC"));
            Assert.That(MolToolsTools.enzymes_compatible("SalI", "XhoI").CompatibleEnd, Is.EqualTo("TCGA"));
            Assert.That(MolToolsTools.enzymes_compatible("EcoRI", "BamHI").CompatibleEnd, Is.Null);

            // Same overhang LETTERS, opposite polarity -> NOT compatible. Biopython 1.88 agrees
            // (SacI % HindIII, SphI % NcoI, ApaI % NotI, FseI % AgeI are all False).
            Assert.That(MolToolsTools.enzymes_compatible("HindIII", "SacI").Compatible, Is.False,
                "5' AGCT cannot anneal to 3' AGCT");
            Assert.That(MolToolsTools.enzymes_compatible("NcoI", "SphI").Compatible, Is.False,
                "5' CATG cannot anneal to 3' CATG");
            Assert.That(MolToolsTools.enzymes_compatible("NotI", "ApaI").Compatible, Is.False,
                "5' GGCC cannot anneal to 3' GGCC");
            Assert.That(MolToolsTools.enzymes_compatible("AgeI", "FseI").Compatible, Is.False,
                "5' CCGG cannot anneal to 3' CCGG");
        });
    }
}
