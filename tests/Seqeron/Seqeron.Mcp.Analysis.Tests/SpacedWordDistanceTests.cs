using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>spaced_word_distance</c> MCP tool (audit round 2 WP6). Expected values: the spaced 1.2.0 program
/// (Leimeister et al. 2014; Ubuntu package) run with -r -f patterns -d JS / -d EU, equal to the Python replica and
/// scipy (jensenshannon(base=2)**2; euclidean on relative frequencies averaged over the patterns). NOT the wrapper output.
/// </summary>
[TestFixture]
public class SpacedWordDistanceTests
{
    private const string S1 = "AGGTAAGGTGGTTGAGATCTGGACTTTTGACGCCTGGAGCCCGCAGTGCTCCTCGAAAAGTAGCCATGCCTTGGGCTGCT";
    private const string S2 = "CAAAGGCCCTACCTTCTTATAGTCCTTTCAACATACAAGTATAGTTGGAAGTTCTAAGTTCAGTTTAATC";
    private static readonly string[] Patterns = { "11011", "10111", "11101" };

    [Test]
    public void SpacedWordDistance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.SpacedWordDistance(S1, S2, Patterns));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SpacedWordDistance("", S2, Patterns));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SpacedWordDistance(S1, null!, Patterns));
        Assert.Throws<ArgumentNullException>(() => AnalysisTools.SpacedWordDistance(S1, S2, null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SpacedWordDistance(S1, S2, Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SpacedWordDistance(S1, S2, new[] { "1101", "11" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SpacedWordDistance(S1, S2, new[] { "0110" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SpacedWordDistance(S1, S2, Patterns, "d2star"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.SpacedWordDistance(S1, S2, Patterns, "bogus"));
    }

    [Test]
    public void SpacedWordDistance_Binding_MatchesSpacedProgram()
    {
        Assert.Multiple(() =>
        {
            // spaced -r -d JS: 0.816322854161; -d EU: 12.4089758166 (12 significant digits).
            Assert.That(AnalysisTools.SpacedWordDistance(S1, S2, Patterns, "jensen_shannon").Distance, Is.EqualTo(0.8163228541607376).Within(1e-12));
            Assert.That(AnalysisTools.SpacedWordDistance(S1, S2, Patterns, "euclidean_counts").Distance, Is.EqualTo(12.40897581662776).Within(1e-10));
            // Paper's Euclidean on relative frequencies (default).
            Assert.That(AnalysisTools.SpacedWordDistance(S1, S2, Patterns).Distance, Is.EqualTo(0.17567404832368613).Within(1e-12));
            // Contiguous pattern = kmer_distance.
            Assert.That(AnalysisTools.SpacedWordDistance(S1, S2, new[] { "1111" }).Distance,
                Is.EqualTo(AnalysisTools.KmerDistance(S1, S2, 4).Distance));
        });
    }
}
