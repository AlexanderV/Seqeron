using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>clustalw_profile_align / clustalw_sequences_to_profile (audit A2-14, B09 F103): the rows equal the installed clustalw 2.1 on the documented example.</summary>
[TestFixture]
public class ClustalWProfileToolsTests
{
    private static readonly string[] Profile = ["ACGTACGTTTGACCAGTAGGCATCGATGGA", "ACGTACGATTGACCAGTAGCCATCGATGGA", "ACGTAGGATTGACCGTAGCC-ATCGATGGA"];
    private static readonly string[] Others = ["TTGACCAGTAGGCATCGATCCCAAGT", "GGGTTTACCAGTAGGCATCGATGGAAACCT"];

    [Test]
    public void ProfileAlign_EqualsClustalW21()
    {
        var r = AlignmentTools.ClustalWProfileAlign(Profile, Others);
        Assert.That(r.AlignedSequences, Is.EqualTo(new[]
        {
            "ACGTACGTTTGACCAGTAGGCATCGATGGA-----",
            "ACGTACGATTGACCAGTAGCCATCGATGGA-----",
            "ACGTAGGATTGACCGTAGCC-ATCGATGGA-----",
            "-----TTGACCAGTAGGCATCGATCCCAAGT----",
            "-----GGGTTTACCAGTAGGCATCGATGGAAACCT",
        }));
    }

    [Test]
    public void SequencesToProfile_EqualsClustalW21()
    {
        var r = AlignmentTools.ClustalWSequencesToProfile(Profile, Others);
        Assert.That(r.AlignedSequences, Is.EqualTo(new[]
        {
            "ACGTACGTTTGACCAGTAGGCATCGATGGA-----",
            "ACGTACGATTGACCAGTAGCCATCGATGGA-----",
            "ACGTAGGATTGACCGTAGCC-ATCGATGGA-----",
            "--------TTGACCAGTAGGCATCGATCCCAAGT-",
            "-----GGGTTTACCAGTAGGCATCGATGGAAACCT",
        }));
    }

    [Test]
    public void BothTools_AcceptProteinAndOptions_AndRefuseBadInput()
    {
        string[] prot = ["MKVLAAGIVGLLAQ", "MKVLAGIVGLLAQP"];
        string[] more = ["MRVLAAGIIGLLQP", "KVLAAGIVALLAQ"];
        var a = AlignmentTools.ClustalWProfileAlign(prot, more, type: "protein", matrix: "blosum", noWeights: true);
        var b = AlignmentTools.ClustalWSequencesToProfile(prot, more, type: "protein", distance: "clustalw-quick", pwMatrix: "pam");
        Assert.That(a.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(prot.Concat(more)));
        Assert.That(b.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(prot.Concat(more)));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ClustalWProfileAlign(Array.Empty<string>(), Others));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ClustalWSequencesToProfile(Profile, Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ClustalWProfileAlign(Profile, Others, type: "rna"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ClustalWSequencesToProfile(Profile, Others, distance: "fast"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ClustalWProfileAlign(Profile, Others, clustering: "ward"));
    }
}
