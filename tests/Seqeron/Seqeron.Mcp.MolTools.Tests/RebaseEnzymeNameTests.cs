using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Review 2026-09 B08 D4 part 1 (F59): the MCP tools that take an enzyme NAME resolve it through
/// <see cref="RestrictionAnalyzer.GetEnzyme"/>, which now covers the full vendored REBASE set (754
/// enzymes with a defined cut), so Type IIS / Type IIB enzymes outside the 39 common ones are reachable
/// by name without any MCP change. Expected values are Biopython 1.88's (BsaI.search / catalyse;
/// BaeI.search on A×20 ACAAAAGTACC T×20 = [11, 44], F56).
/// </summary>
[TestFixture]
public class RebaseEnzymeNameTests
{
    [Test]
    public void GetEnzyme_TypeIISName_IsResolvedFromRebase()
    {
        var result = MolToolsTools.get_enzyme("BsaI");

        Assert.That(result.Enzyme, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(result.Enzyme!.RecognitionSequence, Is.EqualTo("GGTCTC"));
            Assert.That((result.Enzyme.CutPositionForward, result.Enzyme.CutPositionReverse), Is.EqualTo((7, 11)));
            Assert.That(result.Enzyme.Organism, Is.EqualTo("Bacillus stearothermophilus 6-55"));
        });
    }

    [Test]
    public void GetEnzyme_EnzymeWithoutDefinedCut_ReturnsNull()
    {
        Assert.That(MolToolsTools.get_enzyme("Dpi3084I").Enzyme, Is.Null);
    }

    [Test]
    public void RestrictionDigest_TypeIISName_CutsAtBiopythonPositions()
    {
        // BsaI.search(Seq(s)) == [10, 19]; catalyse → lengths [9, 9, 22].
        var result = MolToolsTools.restriction_digest("AAGGTCTCAAAAAAAAAAAAAAAGAGACCAAAAAAAAAAA", new[] { "BsaI" });

        Assert.That(result.Fragments.Select(f => f.Length), Is.EqualTo(new[] { 9, 9, 22 }));
    }

    [Test]
    public void FindRestrictionSites_TwoCutName_ReportsBothCuts()
    {
        string s = new string('A', 20) + "ACAAAAGTACC" + new string('T', 20);
        var result = MolToolsTools.find_restriction_sites(s, new[] { "BaeI" });

        // Biopython BaeI.search == [11, 44] (1-based top-strand cuts): forward records only (non-palindromic).
        Assert.That(result.Sites.Where(x => x.IsForwardStrand).Select(x => x.CutPosition + 1).Order(),
            Is.EqualTo(new[] { 11, 44 }));
    }
}
