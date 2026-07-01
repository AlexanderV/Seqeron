using System;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class LongestDinucleotideRepeatTests
{
    [Test]
    public void LongestDinucleotideRepeat_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.longest_dinucleotide_repeat("ATATAT"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.longest_dinucleotide_repeat(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.longest_dinucleotide_repeat(null!));
    }

    [Test]
    public void LongestDinucleotideRepeat_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // "ATATAT" -> AT repeated 3 times.
            Assert.That(MolToolsTools.longest_dinucleotide_repeat("ATATAT").Repeats, Is.EqualTo(3));
            // "GCGCGCGC" -> GC repeated 4 times.
            Assert.That(MolToolsTools.longest_dinucleotide_repeat("GCGCGCGC").Repeats, Is.EqualTo(4));
            // "ACGT" -> no tandem dinucleotide -> 1.
            Assert.That(MolToolsTools.longest_dinucleotide_repeat("ACGT").Repeats, Is.EqualTo(1));
            // Sequence shorter than 4 nt -> 0.
            Assert.That(MolToolsTools.longest_dinucleotide_repeat("AT").Repeats, Is.EqualTo(0));
        });
    }
}
