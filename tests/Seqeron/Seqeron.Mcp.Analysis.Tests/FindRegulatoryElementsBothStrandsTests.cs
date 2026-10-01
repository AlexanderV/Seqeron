using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// find_regulatory_elements_both_strands — delegation to <c>MotifFinder.FindRegulatoryElements(seq, bothStrands: true)</c>.
/// Expected values: Biopython nt_search on both strands, locked in MotifFinder_RegulatoryStrands_Tests.
/// </summary>
[TestFixture]
public class FindRegulatoryElementsBothStrandsTests
{
    private const string Case1 = "ATTGGTTTATAAACCGCCCATCCAATGGAAAGTCCCTGACTCAGGGCGGA";

    [Test]
    public void FindRegulatoryElementsBothStrands_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindRegulatoryElementsBothStrands(Case1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRegulatoryElementsBothStrands(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRegulatoryElementsBothStrands(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRegulatoryElementsBothStrands("ACGU"));
    }

    [Test]
    public void FindRegulatoryElementsBothStrands_MixedOrientation_EqualsBiopythonNtSearch()
    {
        var hits = AnalysisTools.FindRegulatoryElementsBothStrands(Case1).Items
            .Select(e => (e.Name, e.Position, e.Strand, e.Sequence)).ToArray();
        Assert.That(hits, Is.EqualTo(new[]
        {
            ("TATA Box", 7, "+", "TATAAA"), ("CAAT Box", 0, "-", "CCAAT"), ("CAAT Box", 21, "+", "CCAAT"),
            ("GC Box", 13, "-", "GGGCGG"), ("GC Box", 43, "+", "GGGCGG"), ("AP-1", 36, "+", "TGACTCA"),
            ("NF-κB", 26, "-", "GGGACTTTCC"),
        }));
    }

    [Test]
    public void FindRegulatoryElementsBothStrands_EqualsLibrary()
    {
        const string seq = "GGGACTTTCCATTGGCCAATTTATTTAGGCACGTGCCGCCCGGGCGG";
        var items = AnalysisTools.FindRegulatoryElementsBothStrands(seq).Items;
        var lib = global::Seqeron.Genomics.Analysis.MotifFinder.FindRegulatoryElements(
            new global::Seqeron.Genomics.Core.DnaSequence(seq), bothStrands: true).ToArray();
        Assert.That(items.Select(e => (e.Name, e.Position, e.Sequence, e.Pattern, e.Description, e.Strand)),
            Is.EqualTo(lib.Select(e => (e.Name, e.Position, e.Sequence, e.Pattern, e.Description, e.Strand.ToString()))));
    }
}
