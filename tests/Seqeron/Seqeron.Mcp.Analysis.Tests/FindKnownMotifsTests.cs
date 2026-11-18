using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_known_motifs</c> MCP tool.
/// Expected values from the exact set-matching definition (ascending positions,
/// misses omitted, overlapping reported), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindKnownMotifsTests
{
    [Test]
    public void FindKnownMotifs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindKnownMotifs("ATGATGCC", new[] { "ATG" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindKnownMotifs("", new[] { "ATG" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindKnownMotifs(null!, new[] { "ATG" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindKnownMotifs("XYZ", new[] { "ATG" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindKnownMotifs("ATGATGCC", null!));
    }

    [Test]
    public void FindKnownMotifs_Binding_InvokesSuccessfully()
    {
        // "ATGATGCC": ATG at 0,3; CC at 6; GG absent (omitted).
        var m = AnalysisTools.FindKnownMotifs("ATGATGCC", new[] { "ATG", "CC", "GG" }).Matches;
        Assert.Multiple(() =>
        {
            Assert.That(m.ContainsKey("ATG"), Is.True);
            Assert.That(m["ATG"], Is.EqualTo(new[] { 0, 3 }));
            Assert.That(m.ContainsKey("CC"), Is.True);
            Assert.That(m["CC"], Is.EqualTo(new[] { 6 }));
            Assert.That(m.ContainsKey("GG"), Is.False);
        });

        // Overlapping AA in AAAA -> 0,1,2.
        var overlap = AnalysisTools.FindKnownMotifs("AAAA", new[] { "AA" }).Matches;
        Assert.That(overlap["AA"], Is.EqualTo(new[] { 0, 1, 2 }));
    }

    // EMBOSS fuzznuc 6.6.0 (-pmismatch 0 -complement Y) on GAATTCAAAGAATTCTTTATATAAAGTTCAAGAAC (B09 F19):
    // GAATTC '+' 1-6, 10-15 and '-' 1-6, 10-15; TATAWAW '+' 18-24, '-' 17-23 (forward-strand 1-based coordinates).
    [Test]
    public void FindKnownMotifs_BothStrandsAndIupac_MatchFuzznuc()
    {
        const string seq = "GAATTCAAAGAATTCTTTATATAAAGTTCAAGAAC";

        var eco = AnalysisTools.FindKnownMotifs(seq, new[] { "GAATTC" }, bothStrands: true);
        Assert.Multiple(() =>
        {
            Assert.That(eco.Matches["GAATTC"], Is.EqualTo(new[] { 0, 9 }));
            Assert.That(eco.ReverseMatches!["GAATTC"], Is.EqualTo(new[] { 0, 9 }));
        });

        var tata = AnalysisTools.FindKnownMotifs(seq, new[] { "TATAWAW" }, bothStrands: true, iupac: true);
        Assert.Multiple(() =>
        {
            Assert.That(tata.Matches["TATAWAW"], Is.EqualTo(new[] { 17 }));
            Assert.That(tata.ReverseMatches!["TATAWAW"], Is.EqualTo(new[] { 16 }));
        });

        // Defaults keep the exact forward-only contract: IUPAC symbols are literal, no reverse hits.
        var plain = AnalysisTools.FindKnownMotifs(seq, new[] { "TATAWAW", "GAATTC" });
        Assert.Multiple(() =>
        {
            Assert.That(plain.Matches.ContainsKey("TATAWAW"), Is.False);
            Assert.That(plain.Matches["GAATTC"], Is.EqualTo(new[] { 0, 9 }));
            Assert.That(plain.ReverseMatches, Is.Empty);
        });
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindKnownMotifs(seq, new[] { "AXG" }, iupac: true));
    }
}
