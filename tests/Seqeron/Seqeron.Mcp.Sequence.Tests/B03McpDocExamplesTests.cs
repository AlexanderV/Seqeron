using System.Text.Json;
using ModelContextProtocol;
using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

/// <summary>
/// FIN-B03 WP16 (B03.md F42, audit rows A29/A30): locks the numbers and wire names quoted in
/// docs/mcp/tools/sequence/{linguistic_complexity,shannon_entropy,nucleotide_composition}.md.
/// Reference values: exact-fraction Python (Σ V_k / Σ min(a^k, N−k+1)) and scipy-equivalent log2 entropy,
/// executed 2026-10-14.
/// </summary>
[TestFixture]
public class B03McpDocExamplesTests
{
    [TestCase("ATGCATGCATGC", 24.0 / 49.0)] // doc example 1
    [TestCase("AAAAAAAAAA", 2.0 / 13.0)]    // doc example 2
    [TestCase("ACGTN", 1.0)]                // a = 5 inferred → LC ≤ 1
    public void LinguisticComplexity_DocExamples(string sequence, double expected) =>
        Assert.That(SequenceTools.LinguisticComplexity(sequence).Complexity, Is.EqualTo(expected).Within(1e-15));

    [Test]
    public void LinguisticComplexity_ProteinFixedAlphabet_DocExample() =>
        Assert.That(SequenceTools.LinguisticComplexity("MKVLAAGIVGLLLAA", 15, 20).Complexity,
            Is.EqualTo(0.9).Within(1e-15));

    [TestCase("ATGCATGCAT", 1.9709505944546686)] // A3 T3 G2 C2 (doc example 1)
    [TestCase("AAAAAAAAAA", 0.0)]
    [TestCase("acgtNN-*", 2.2516291673878226)]   // N distinct letter; '-' and '*' ignored
    public void ShannonEntropy_DocAlphabet(string sequence, double expected) =>
        Assert.That(SequenceTools.ShannonEntropy(sequence).Entropy, Is.EqualTo(expected).Within(1e-12));

    [Test]
    public void NucleotideComposition_WireNames_AreCamelCase()
    {
        string json = JsonSerializer.Serialize(SequenceTools.NucleotideComposition("ATGCATGC"), McpJsonUtilities.DefaultOptions);
        using var doc = JsonDocument.Parse(json);
        var names = doc.RootElement.EnumerateObject().Select(p => p.Name).ToArray();
        Assert.That(names, Is.EquivalentTo(new[] { "length", "a", "t", "g", "c", "u", "other", "gcContent", "n" }));
    }
}
