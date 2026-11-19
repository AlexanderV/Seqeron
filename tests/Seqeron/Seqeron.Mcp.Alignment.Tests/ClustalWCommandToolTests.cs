using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>clustalw_command (audit A2-14b / c / d, B09 F104-F106): the command line of ClustalW 2.1 as a tool; the files equal what the installed clustalw 2.1 writes for the documented examples.</summary>
[TestFixture]
public class ClustalWCommandToolTests
{
    private const string Fasta = ">human\nACGTACGTTTGACCAGTAGGCATCGATGGATTACA\n>chimp\nACGTACGATTGACCAGTAGCCATCGATGGATTACA\n>mouse\nACGTAGGATTGACCGTAGCCATCGATGGTTACAC\n>rat\nACGTAGGATTGACCGTAGCCATCGAAGGTTACAC\n";

    [Test]
    public void Align_WritesTheAlignmentAndTheGuideTree()
    {
        // clustalw -INFILE=ex.fa -TYPE=DNA -ALIGN -OUTPUT=CLUSTAL -SEQNOS=ON
        var r = AlignmentTools.ClustalWCommand(["-INFILE=ex.fa", "-TYPE=DNA", "-ALIGN", "-OUTPUT=CLUSTAL", "-SEQNOS=ON"], new Dictionary<string, string> { ["ex.fa"] = Fasta });
        Assert.Multiple(() =>
        {
            Assert.That(r.Files.Keys, Is.EquivalentTo(new[] { "ex.aln", "ex.dnd" }));
            Assert.That(r.Files["ex.aln"], Is.EqualTo("CLUSTAL 2.1 multiple sequence alignment\n\n\nmouse           ACGTAGGATTGACC-GTAGCCATCGATGG-TTACAC 34\nrat             ACGTAGGATTGACC-GTAGCCATCGAAGG-TTACAC 34\n"
                + "chimp           ACGTACGATTGACCAGTAGCCATCGATGGATTACA- 35\nhuman           ACGTACGTTTGACCAGTAGGCATCGATGGATTACA- 35\n                ***** * ****** **** ****** ** ***** \n"));
            Assert.That(r.Files["ex.dnd"], Is.EqualTo("(\n(\nhuman:0.05798,\nchimp:-0.00084)\n:0.05966,\nmouse:0.00000,\nrat:0.02941);\n"));
        });
    }

    [Test]
    public void UseTree_AlignsAlongTheGivenTree_AndRefusesBadTrees()
    {
        var files = new Dictionary<string, string> { ["ex.fa"] = Fasta, ["t.dnd"] = "((human:0.05,chimp:0.05):0.02,(mouse:0.03,rat:0.03):0.02);" };
        var r = AlignmentTools.ClustalWCommand(["ex.fa", "-USETREE=t.dnd", "-TYPE=DNA", "-OUTPUT=FASTA"], files);
        Assert.That(r.Files["ex.fasta"], Does.StartWith(">"));
        files["t.dnd"] = "(human,chimp,nobody);";
        Assert.Throws<Seqeron.Genomics.Infrastructure.ClustalWTreeFileException>(() => AlignmentTools.ClustalWCommand(["ex.fa", "-USETREE=t.dnd", "-TYPE=DNA"], files));
    }

    [Test]
    public void BadArguments_AreRefused()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.ClustalWCommand([]));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ClustalWCommand(["-INFILE=ex.fa", "-BOGUS"], new Dictionary<string, string> { ["ex.fa"] = Fasta }));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ClustalWCommand(["-INFILE=missing.fa", "-ALIGN"]));
    }
}
