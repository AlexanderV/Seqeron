using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>MCP wrappers of the other programs of the MUMmer 3.23 package (B09 A3-9d / A3-9e / A3-9f): the output of the installed program on small examples.</summary>
[TestFixture]
public class MummerProgramToolsTests
{
    [Test]
    public void Mummer_WritesTheMatchesOfTheProgram()
    {
        // mummer -mum -l 8 -b -s r.fa q.fa (MUMmer 3.23)
        const string reference = ">r\nacgtacgtacgtaaccggtt\n";
        const string query = ">q\nacgtacgtaacctacgtacgtacgtaaccggtta\n";
        var r = AnalysisTools.Mummer(reference, [query], "mum", minLength: 8, both: true, showString: true);
        Assert.That(r.Text, Is.EqualTo("> q\n       1        14        20\nacgtacgtacgtaaccggtt\n> q Reverse\n       1        10        13\nacgtacgtacgta\n      12         1         9\ntaaccggtt\n"));
        Assert.That(r.Messages, Does.StartWith("# reading input file \"reference.fa\" of length 20\n"));
        // mummer -maxmatch -l 4 (every maximal match, the suffix-tree order: ref 9, 5, 1 for the first query position)
        Assert.That(AnalysisTools.Mummer(reference, [">q\nacgtacgtaacc\ntacgtacgtacgtaaccggtta\n"], "maxmatch", minLength: 4).Text,
            Is.EqualTo("> q\n       9         1         5\n       5         1        12\n       1         1         9\n       1         5         5\n       8        13         6\n       4        13        10\n       1        14        20\n       1        18         9\n       1        22         5\n"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Mummer(reference, [query], both: true, reverseOnly: true));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Mummer(reference, [query], matchMode: "x"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Mummer("", [query]));
        Assert.Throws<ArgumentException>(() => AnalysisTools.Mummer(reference, []));
        Assert.Throws<FormatException>(() => AnalysisTools.Mummer(reference, [">q\n"]));
    }
}
