// SPDX-License-Identifier: GPL-3.0-or-later
// Wrappers of the copyleft multiple-alignment assemblies (Kalign 3.3.5, GPL-3.0-or-later and others).
// Part of Seqeron.Mcp.MultipleAlignment, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.PartialOrderAlignment;
using Seqeron.Genomics.Kalign;
using Seqeron.Genomics.TCoffee;
using Seqeron.Genomics.ClustalW;
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.MultipleAlignment.Tools;

namespace Seqeron.Mcp.MultipleAlignment.Tests;

/// <summary>clustalw_profile_align / clustalw_sequences_to_profile (audit A2-14, B09 F103): the rows equal the installed clustalw 2.1 on the documented example.</summary>
[TestFixture]
public class ClustalWProfileToolsTests
{
    private static readonly string[] Profile = ["ACGTACGTTTGACCAGTAGGCATCGATGGA", "ACGTACGATTGACCAGTAGCCATCGATGGA", "ACGTAGGATTGACCGTAGCC-ATCGATGGA"];
    private static readonly string[] Others = ["TTGACCAGTAGGCATCGATCCCAAGT", "GGGTTTACCAGTAGGCATCGATGGAAACCT"];

    [Test]
    public void ProfileAlign_EqualsClustalW21()
    {
        var r = MultipleAlignmentTools.ClustalWProfileAlign(Profile, Others);
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
        var r = MultipleAlignmentTools.ClustalWSequencesToProfile(Profile, Others);
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
        var a = MultipleAlignmentTools.ClustalWProfileAlign(prot, more, type: "protein", matrix: "blosum", noWeights: true);
        var b = MultipleAlignmentTools.ClustalWSequencesToProfile(prot, more, type: "protein", distance: "clustalw-quick", pwMatrix: "pam");
        Assert.That(a.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(prot.Concat(more)));
        Assert.That(b.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(prot.Concat(more)));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ClustalWProfileAlign(Array.Empty<string>(), Others));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ClustalWSequencesToProfile(Profile, Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ClustalWProfileAlign(Profile, Others, type: "rna"));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ClustalWSequencesToProfile(Profile, Others, distance: "fast"));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ClustalWProfileAlign(Profile, Others, clustering: "ward"));
    }
}
