// SPDX-License-Identifier: GPL-3.0-or-later
// Wrappers of the copyleft multiple-alignment assemblies (Kalign 3.3.5, GPL-3.0-or-later and others).
// Part of Seqeron.Mcp.MultipleAlignment, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.PartialOrderAlignment;
using Seqeron.Genomics.Kalign;
using Seqeron.Genomics.TCoffee;
using Seqeron.Genomics.ClustalW;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.MultipleAlignment.Tools;

namespace Seqeron.Mcp.MultipleAlignment.Tests;

/// <summary>
/// needle_align / water_align / clustalw_guide_tree / tcoffee_blast_msa_library / tcoffee_tblastx_library (audit A3-17, B09 F141): the tools delegate to
/// <see cref="SequenceAligner"/>; the library tests lock the rows against the installed EMBOSS 6.6.0 / ClustalW 2.1 / t_coffee 13.41 on their committed oracles.
/// </summary>
[TestFixture]
public class NeedleWaterGuideTreeBlastToolTests
{

    private const string BlastpXml =
        "<BlastOutput><BlastOutput_iterations>" +
        "<Iteration><Iteration_iter-num>1</Iteration_iter-num><Iteration_hits>" +
        "<Hit><Hit_num>1</Hit_num><Hit_def>b</Hit_def><Hit_hsps><Hsp><Hsp_query-from>1</Hsp_query-from><Hsp_hit-from>1</Hsp_hit-from><Hsp_identity>4</Hsp_identity><Hsp_align-len>5</Hsp_align-len>" +
        "<Hsp_qseq>ACDEF</Hsp_qseq><Hsp_hseq>ACD-F</Hsp_hseq></Hsp></Hit_hsps></Hit></Iteration_hits></Iteration>" +
        "<Iteration><Iteration_iter-num>2</Iteration_iter-num><Iteration_hits></Iteration_hits></Iteration>" +
        "</BlastOutput_iterations></BlastOutput>";

    [Test]
    public void ClustalWGuideTree_TwoSequences_WritesHalfTheDistanceOnBothBranches()
    {
        var r = MultipleAlignmentTools.ClustalWGuideTree(["a", "b"], [[0, 0.5], [0.5, 0]]);
        Assert.Multiple(() =>
        {
            Assert.That(r.Tree, Does.StartWith("(a:0.25,b:0.25)"));
            Assert.That(r.Sequences, Is.EqualTo(2));
        });
    }

    [Test]
    public void ClustalWGuideTree_FourSequences_EqualsTheLibrary()
    {
        string[] names = ["s1", "s2", "s3", "s4"];
        double[][] d = [[0, 0.1, 0.5, 0.6], [0.1, 0, 0.55, 0.62], [0.5, 0.55, 0, 0.2], [0.6, 0.62, 0.2, 0]];
        var matrix = new double[4, 4];
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                matrix[i, j] = d[i][j];
        var r = MultipleAlignmentTools.ClustalWGuideTree(names, d);
        Assert.Multiple(() =>
        {
            Assert.That(r.Tree, Is.EqualTo(ClustalWAligner.ClustalWGuideTreeText(names, matrix)));
            Assert.That(r.Tree, Does.Contain("s1").And.Contain("s4"));
        });
    }

    [Test]
    public void ClustalWGuideTree_Validation()
    {
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ClustalWGuideTree(["only"], [[0]]));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ClustalWGuideTree(["a", "b"], [[0, 1]]));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ClustalWGuideTree(["a", "b"], [[0, 1], [1]]));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ClustalWGuideTree(null!, [[0]]));
    }

    [Test]
    public void TCoffeeBlastMsaLibrary_EqualsTheLibrary_WeightIsTheTruncatedPercentIdentity()
    {
        string[] names = ["a", "b"];
        string[] seqs = ["ACDEF", "ACDF"];
        var tool = MultipleAlignmentTools.TCoffeeBlastMsaLibrary(names, seqs, BlastpXml);
        var lib = TCoffeeAligner.TCoffeeBlastMsaLibraryEntries(names, seqs, BlastpXml);
        Assert.Multiple(() =>
        {
            Assert.That(tool.Count, Is.EqualTo(lib.Count));
            Assert.That(tool.Entries.Select(e => (e.SeqA, e.PosA, e.SeqB, e.PosB, e.Weight, e.Count)), Is.EqualTo(lib.Select(e => (e.SeqA, e.PosA, e.SeqB, e.PosB, e.Weight, e.Count))));
            // columns with a residue in both rows: A-A, C-C, D-D, F-F at 4 / 5 identity = 80
            Assert.That(tool.Entries, Is.Not.Empty);
            Assert.That(tool.Entries.All(e => e.Weight == 80 && e.SeqA < e.SeqB), Is.True);
            Assert.That(tool.Entries.Select(e => (e.PosA, e.PosB)), Is.EquivalentTo(new[] { (0, 0), (1, 1), (2, 2), (4, 3) }));
        });
    }

    [Test]
    public void TCoffeeTblastxLibrary_EqualsTheLibrary_OnAMinimalSearch()
    {
        string[] seqs = ["ACGTACGTACGT", "ACGTACGTACGT"];
        const string xml =
            "<BlastOutput><BlastOutput_iterations>" +
            "<Iteration><Iteration_iter-num>1</Iteration_iter-num><Iteration_query-len>12</Iteration_query-len><Iteration_hits>" +
            "<Hit><Hit_num>1</Hit_num><Hit_accession>1</Hit_accession><Hit_len>12</Hit_len><Hit_hsps><Hsp><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>6</Hsp_query-to><Hsp_query-frame>1</Hsp_query-frame>" +
            "<Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>6</Hsp_hit-to><Hsp_hit-frame>1</Hsp_hit-frame><Hsp_identity>2</Hsp_identity><Hsp_align-len>2</Hsp_align-len></Hsp></Hit_hsps></Hit>" +
            "</Iteration_hits></Iteration></BlastOutput_iterations></BlastOutput>";
        var tool = MultipleAlignmentTools.TCoffeeTblastxLibrary(seqs, xml);
        var lib = TCoffeeAligner.TCoffeeTblastxLibraryEntries(seqs, xml);
        Assert.Multiple(() =>
        {
            Assert.That(tool.Count, Is.EqualTo(lib.Count));
            Assert.That(tool.Entries.Select(e => (e.SeqA, e.PosA, e.SeqB, e.PosB, e.Weight)), Is.EqualTo(lib.Select(e => (e.SeqA, e.PosA, e.SeqB, e.PosB, e.Weight))));
            Assert.That(tool.Entries, Is.Not.Empty);
            Assert.That(tool.Entries.All(e => e.Weight == 100), Is.True);
        });
    }

    [Test]
    public void TCoffeeBlastLibraries_Validation()
    {
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.TCoffeeBlastMsaLibrary(["a"], [], BlastpXml));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.TCoffeeBlastMsaLibrary(["a", "b"], ["ACD"], BlastpXml));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.TCoffeeBlastMsaLibrary(["a"], ["ACD"], null!));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.TCoffeeTblastxLibrary([], "<x/>"));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.TCoffeeTblastxLibrary(["ACGT", ""], "<x/>"));
    }}
