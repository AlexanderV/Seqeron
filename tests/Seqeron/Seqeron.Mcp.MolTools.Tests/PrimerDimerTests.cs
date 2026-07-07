using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class PrimerDimerTests
{
    [Test]
    public void PrimerDimer_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.primer_dimer("AAAAAAAA", "AAAAAAAA"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.primer_dimer("", "AAAA"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.primer_dimer("AAAA", null!));
    }

    [Test]
    public void PrimerDimer_Binding_InvokesSuccessfully()
    {
        // primer1 = AAAAAAAA; revcomp(primer2=AAAAAAAA) = TTTTTTTT.
        // 3' window (8 bp): every A pairs with T -> 8 complementary -> dimer.
        var strong = MolToolsTools.primer_dimer("AAAAAAAA", "AAAAAAAA");
        Assert.Multiple(() =>
        {
            Assert.That(strong.ComplementaryBases, Is.EqualTo(8));
            Assert.That(strong.HasDimer, Is.True);
        });

        // primer2 = AAAACCCC -> revcomp = GGGGTTTT. Compared to AAAAAAAA:
        // A-G,A-G,A-G,A-G,A-T,A-T,A-T,A-T -> 4 complementary.
        var four = MolToolsTools.primer_dimer("AAAAAAAA", "AAAACCCC", min_complementarity: 5);
        Assert.Multiple(() =>
        {
            Assert.That(four.ComplementaryBases, Is.EqualTo(4));
            // 4 < min 5 -> not flagged.
            Assert.That(four.HasDimer, Is.False);
            // With default threshold 4, the same 4 complementary bases DO flag a dimer.
            Assert.That(MolToolsTools.primer_dimer("AAAAAAAA", "AAAACCCC").HasDimer, Is.True);
        });
    }
}
