using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// R14-10b (B09 WP67, F307): the MCP tools <c>codon_aligner</c>, <c>multiple_seq_alignment</c> and the gap-score functions / gap scores above zero of <c>pairwise_aligner</c>. Expected values: Biopython 1.88 (the doctests of Bio.Align.CodonAligner and
/// MultipleSeqAlignment, and the runs named in each test).
/// </summary>
[TestFixture]
public class LegacyAlignToolsTests
{
    [Test]
    public void PairwiseAligner_GapFunction_EqualsBiopython()
    {
        // aligner.gap_score = lambda i, l: -(2 + 0.5 * l); match 2, mismatch -1: 'Waterman-Smith-Beyer global alignment algorithm', score 11.0, one alignment
        var r = AlignmentTools.PairwiseAlignerAlign("ACGTTACGT", "ACGACGT", matchScore: 2, mismatchScore: -1, gapFunctionsJson: "{\"gap_score\": {\"kind\": \"lin\", \"p\": [2, 0.5]}}");
        Assert.Multiple(() =>
        {
            Assert.That(r.Algorithm, Is.EqualTo("Waterman-Smith-Beyer global alignment algorithm"));
            Assert.That(r.Score, Is.EqualTo(11.0));
            Assert.That(r.Count, Is.EqualTo(1));
            Assert.That(r.Alignments[0].Coordinates, Is.EqualTo(new[] { new long[] { 0, 3, 5, 9 }, [0, 3, 3, 7] }));
        });
        Assert.Throws<ArgumentException>(() => AlignmentTools.PairwiseAlignerAlign("A", "A", gapFunctionsJson: "{\"gap_score\": {\"kind\": \"nope\", \"p\": []}}"));
    }

    [Test]
    public void PairwiseAligner_PositiveGapScore_LocalEqualsBiopython()
    {
        // local, match 3, mismatch -1, gap_score +1: Smith-Waterman, score 11.0 and the three alignments Biopython lists (two of them end at cells it did not clear)
        var r = AlignmentTools.PairwiseAlignerAlign("GTTCA", "TTACC", mode: "local", matchScore: 3, mismatchScore: -1, gapScoresJson: "{\"gap_score\": 1}");
        Assert.Multiple(() =>
        {
            Assert.That(r.Algorithm, Is.EqualTo("Smith-Waterman"));
            Assert.That(r.Score, Is.EqualTo(11.0));
            Assert.That(r.Alignments.Select(a => string.Join("|", a.Coordinates.Select(c => string.Join(",", c)))), Is.EquivalentTo(new[] { "1,1,2|0,1,2", "0,1,2|1,1,2", "1,3|0,2" }));
        });
    }

    [Test]
    public void CodonAligner_FrameshiftDoctests_EqualBiopython()
    {
        var r = AlignmentTools.CodonAlign("MLGSNESRVCPKLSPS", "ATGCTGGGCTCGAACGAGTCCGTGTATGCCCTAAGCTGAGCCCGTCG");
        Assert.That(r.Score, Is.EqualTo(13.0));
        Assert.That(r.Count, Is.EqualTo(1));
        Assert.That(r.Alignments[0].Coordinates, Is.EqualTo(new[] { new long[] { 0, 7, 7, 16 }, [0, 21, 20, 47] }));
        var ambiguous = AlignmentTools.CodonAlign("FKKKKF", "TTTAAAAAAAAAAATTT", frameshiftScore: -1);
        Assert.That(ambiguous.Score, Is.EqualTo(5.0));
        Assert.That(ambiguous.Alignments.Select(a => string.Join("|", a.Coordinates.Select(c => string.Join(",", c)))),
            Is.EqualTo(new[] { "0,4,4,6|0,12,11,17", "0,3,3,6|0,9,8,17", "0,2,2,6|0,6,5,17" }));
        Assert.That(AlignmentTools.CodonAlign("MSR", "AUGUCUCGU", scoreOnly: true).Score, Is.EqualTo(3.0));
        Assert.Throws<ArgumentException>(() => AlignmentTools.CodonAlign("MSR", "AT?TCTCGT", scoreOnly: true));
    }

    [Test]
    public void CodonAligner_UndefinedReference_IsRefused()
    {
        // 'GD' against 'GGGG': the optimal alignment skips the first residue; Biopython's iterator then returns [[1], [3]] for ever (a heap overflow in PathGenerator_create_path), score() is 0.0
        Assert.That(AlignmentTools.CodonAlign("GD", "GGGG", scoreOnly: true).Score, Is.EqualTo(0.0));
        Assert.Throws<InvalidOperationException>(() => AlignmentTools.CodonAlign("GD", "GGGG"));
    }

    private const string Rows = "[{\"id\": \"Alpha\", \"seq\": \"AAAACGT\"}, {\"id\": \"Beta\", \"seq\": \"AAA-CGT\"}, {\"id\": \"Gamma\", \"seq\": \"AAAAGGT\"}, {\"id\": \"Delta\", \"seq\": \"AAAACGT\"}, {\"id\": \"Epsilon\", \"seq\": \"AAA-GGT\"}]";

    [Test]
    public void MultipleSeqAlignment_Indexing_EqualsBiopythonDoctests()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AlignmentTools.MultipleSeqAlignmentOperation("cell", Rows, row: 3, column: 4).ResultJson, Is.EqualTo("\"C\""));
            Assert.That(AlignmentTools.MultipleSeqAlignmentOperation("column", Rows, column: 4).ResultJson, Is.EqualTo("\"CCGCG\""));
            Assert.That(AlignmentTools.MultipleSeqAlignmentOperation("column", Rows, rows: "1:3", column: 4).ResultJson, Is.EqualTo("\"CG\""));
            Assert.That(AlignmentTools.MultipleSeqAlignmentOperation("slice", Rows, rows: "1:5", columns: "3:6").Text, Is.EqualTo("Alignment with 4 rows and 3 columns\n-CG Beta\nAGG Gamma\nACG Delta\n-GG Epsilon"));
            Assert.That(AlignmentTools.MultipleSeqAlignmentOperation("rows", Rows, rows: "::-1").Text, Is.EqualTo("Alignment with 5 rows and 7 columns\nAAA-GGT Epsilon\nAAAACGT Delta\nAAAAGGT Gamma\nAAA-CGT Beta\nAAAACGT Alpha"));
            Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.MultipleSeqAlignmentOperation("row", Rows, row: 9));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleSeqAlignmentOperation("append", Rows, newRecordsJson: "{\"seq\": \"AC\"}"));
        });
    }

    [Test]
    public void MultipleSeqAlignment_AddSortSubstitutions_EqualBiopythonDoctests()
    {
        var left = "[{\"id\": \"Alpha\", \"seq\": \"AAAAC\"}, {\"id\": \"Beta\", \"seq\": \"AAA-C\"}, {\"id\": \"Gamma\", \"seq\": \"AAAAG\"}]";
        var other = "{\"records\": [{\"id\": \"Alpha\", \"seq\": \"GT\"}, {\"id\": \"Beta\", \"seq\": \"GT\"}, {\"id\": \"Gamma\", \"seq\": \"GT\"}], \"annotations\": {\"tool\": \"demo\", \"name\": \"end\"}, \"column_annotations\": {\"stats\": \"CC\"}}";
        var added = AlignmentTools.MultipleSeqAlignmentOperation("add", left, annotationsJson: "{\"tool\": \"demo\", \"name\": \"start\"}", columnAnnotationsJson: "{\"stats\": \"CCCXC\"}", otherJson: other);
        var json = JsonNode.Parse(added.ResultJson)!;
        Assert.Multiple(() =>
        {
            Assert.That(added.Text, Is.EqualTo("Alignment with 3 rows and 7 columns\nAAAACGT Alpha\nAAA-CGT Beta\nAAAAGGT Gamma"));
            Assert.That(json["annotations"]!.ToJsonString(), Is.EqualTo("{\"tool\":\"demo\"}"));
            Assert.That(json["column_annotations"]!.ToJsonString(), Is.EqualTo("{\"stats\":\"CCCXCCC\"}"));
        });

        var four = "[{\"id\": \"seq1\", \"seq\": \"ACGT\"}, {\"id\": \"seq2\", \"seq\": \"A--A\"}, {\"id\": \"seq3\", \"seq\": \"ACGT\"}, {\"id\": \"seq4\", \"seq\": \"TTTC\"}]";
        var m = JsonNode.Parse(AlignmentTools.MultipleSeqAlignmentOperation("substitutions", four).ResultJson)!;
        Assert.Multiple(() =>
        {
            Assert.That(m["letters"]!.GetValue<string>(), Is.EqualTo("ACGT"));
            Assert.That(m["values"]!.ToJsonString(), Is.EqualTo("[[3,0.5,0,2.5],[0.5,1,0,2],[0,0,1,1],[2.5,2,1,1]]"));
        });

        var sorted = AlignmentTools.MultipleSeqAlignmentOperation("sort", "[{\"id\": \"Human\", \"seq\": \"ACGT\"}, {\"id\": \"Mouse\", \"seq\": \"ACGG\"}, {\"id\": \"Chicken\", \"seq\": \"ACGC\"}]", reverse: true);
        Assert.That(sorted.Text, Is.EqualTo("Alignment with 3 rows and 4 columns\nACGG Mouse\nACGT Human\nACGC Chicken"));
    }
}
