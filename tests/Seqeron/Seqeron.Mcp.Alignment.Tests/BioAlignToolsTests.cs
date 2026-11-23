using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// R14-10 (B09 F303): the MCP tools of the Bio.Align port: <c>pairwise_aligner</c> (FOGSAA, wildcard, strand), <c>alignment_counts</c> and <c>alignment_operation</c>. Expected values: Biopython 1.88 on the fitting alignment of "TTACGTTT" and "ACGT"
/// (coordinates [[0, 2, 6, 8], [0, 0, 4, 4]], match 2, mismatch -1, gap -2: score 0.0, two deletions on the left and on the right).
/// </summary>
[TestFixture]
public class BioAlignToolsTests
{
    private const string Alignment = "{\"sequences\": [\"TTACGTTT\", \"ACGT\"], \"coordinates\": [[0, 2, 6, 8], [0, 0, 4, 4]]}";

    [Test]
    public void PairwiseAligner_Fogsaa_EqualsBiopython()
    {
        var r = AlignmentTools.PairwiseAlignerAlign("TTACGTTT", "ACGT", mode: "fogsaa", matchScore: 2, mismatchScore: -1, gapScoresJson: "{\"gap_score\": -2}");
        Assert.That(r.Algorithm, Is.EqualTo("Fast Optimal Global Sequence Alignment Algorithm"));
        Assert.That(r.Score, Is.EqualTo(0.0));
        Assert.That(r.Count, Is.EqualTo(1));
        Assert.That(r.Alignments[0].Coordinates, Is.EqualTo(new[] { new long[] { 0, 2, 6, 8 }, [0, 0, 4, 4] }));
        Assert.That(r.Alignments[0].AlignedSequence1, Is.EqualTo("TTACGTTT"));
        Assert.That(r.Alignments[0].AlignedSequence2, Is.EqualTo("--ACGT--"));
        Assert.That(r.Alignments[0].Text, Is.EqualTo("target            0 TTACGTTT 8\n                  0 --||||-- 8\nquery             0 --ACGT-- 4\n"));
    }

    [Test]
    public void PairwiseAligner_WildcardAndScoreOnly()
    {
        var r = AlignmentTools.PairwiseAlignerAlign("ACNGT", "ACAGT", matchScore: 2, mismatchScore: -1, wildcard: "N", gapScoresJson: "{\"gap_score\": -2}", scoreOnly: true);
        Assert.That(r.Score, Is.EqualTo(8.0));
        Assert.That(r.Alignments, Is.Empty);
        Assert.Throws<ArgumentException>(() => AlignmentTools.PairwiseAlignerAlign("A", "A", wildcard: "NN"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.PairwiseAlignerAlign("A", "A", mode: "banded"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.PairwiseAlignerAlign("A", "A", strand: "x"));
    }

    [Test]
    public void AlignmentCounts_EqualBiopython()
    {
        var c = AlignmentTools.AlignmentCountsOf(Alignment);
        Assert.That((c.Aligned, c.Identities, c.Mismatches, c.Gaps, c.LeftGaps, c.RightGaps, c.OpenLeftDeletions, c.ExtendLeftDeletions, c.OpenRightDeletions, c.ExtendRightDeletions), Is.EqualTo((4L, 4L, 0L, 4L, 2L, 2L, 1L, 1L, 1L, 1L)));
        Assert.That(c.Score, Is.Null);
        Assert.That(c.Text, Does.StartWith("AlignmentCounts object with\n    aligned = 4:\n        identities = 4,"));
        var scored = AlignmentTools.AlignmentCountsOf(Alignment, matchScore: 2, mismatchScore: -1, gapScoresJson: "{\"gap_score\": -2}");
        Assert.That(scored.Score, Is.EqualTo(0.0));
        Assert.That(scored.GapScore, Is.EqualTo(-8.0));
        Assert.That(scored.SubstitutionScore, Is.EqualTo(8.0));
    }

    [TestCase("shape", "[2, 8]")]
    [TestCase("rows", "[\"TTACGTTT\",\"--ACGT--\"]")]
    [TestCase("aligned", "[[[2,6]],[[0,4]]]")]
    [TestCase("indices", "[[0,1,2,3,4,5,6,7],[-1,-1,0,1,2,3,-1,-1]]")]
    public void AlignmentOperation_EqualBiopython(string operation, string expected)
    {
        var r = AlignmentTools.AlignmentOperationOn(operation, Alignment);
        Assert.That(JsonNode.Parse(r.ResultJson)!.ToJsonString(), Is.EqualTo(JsonNode.Parse(expected)!.ToJsonString()));
    }

    [Test]
    public void AlignmentOperation_SliceFrequenciesSubstitutionsAndText()
    {
        var slice = JsonNode.Parse(AlignmentTools.AlignmentOperationOn("slice", Alignment, rows: ":", columns: "1:5").ResultJson)!;
        Assert.That(slice["coordinates"]!.ToJsonString(), Is.EqualTo("[[1,2,5],[0,0,3]]"));
        var freq = JsonNode.Parse(AlignmentTools.AlignmentOperationOn("frequencies", Alignment).ResultJson)!;
        Assert.That(freq["T"]!.ToJsonString(), Is.EqualTo("[1,1,0,0,0,2,1,1]"));
        Assert.That(freq["-"]!.ToJsonString(), Is.EqualTo("[1,1,0,0,0,0,1,1]"));
        var subs = JsonNode.Parse(AlignmentTools.AlignmentOperationOn("substitutions", Alignment).ResultJson)!;
        Assert.That(subs["letters"]!.GetValue<string>(), Is.EqualTo("ACGT"));
        var rc = JsonNode.Parse(AlignmentTools.AlignmentOperationOn("reverse_complement", Alignment).ResultJson)!;
        Assert.That(rc["coordinates"]!.ToJsonString(), Is.EqualTo("[[0,2,6,8],[0,0,4,4]]"));
        Assert.That(AlignmentTools.AlignmentOperationOn("text", Alignment).Text, Is.EqualTo("target            0 TTACGTTT 8\n                  0 --||||-- 8\nquery             0 --ACGT-- 4\n"));
        var parsed = JsonNode.Parse(AlignmentTools.AlignmentOperationOn("parse_printed", printedRows: ["TAGGCATACGTG", "AACG--TACGT-", "-ACGCATACTTG"]).ResultJson)!;
        Assert.That(parsed["coordinates"]!.ToJsonString(), Is.EqualTo("[[0,1,4,6,11,12],[0,1,4,4,9,9],[0,0,3,5,10,11]]"));
        Assert.That(AlignmentTools.AlignmentOperationOn("letter", Alignment, rows: "1", columns: "2").ResultJson, Is.EqualTo("\"A\""));
        Assert.Throws<ArgumentException>(() => AlignmentTools.AlignmentOperationOn("nosuch", Alignment));
        Assert.Throws<ArgumentException>(() => AlignmentTools.AlignmentOperationOn("add", Alignment));
    }

    [Test]
    public void AlignmentOperation_MapLiftsTheReadFromTheTranscriptToTheGenome()
    {
        // Biopython 1.88: genome vs transcript [[2, 5, 8], [0, 3, 6]] mapped through transcript vs read [[0, 1, 6], [0, 0, 5]] gives genome vs read [[3, 5, 8], [0, 2, 5]]
        const string transcriptOnGenome = "{\"sequences\": [\"CCAACGTTGG\", \"AACGTT\"], \"coordinates\": [[2, 5, 8], [0, 3, 6]]}";
        const string readOnTranscript = "{\"sequences\": [\"AACGTT\", \"ACGTT\"], \"coordinates\": [[0, 1, 6], [0, 0, 5]]}";
        var mapped = JsonNode.Parse(AlignmentTools.AlignmentOperationOn("map", transcriptOnGenome, readOnTranscript).ResultJson)!;
        Assert.That(mapped["coordinates"]!.ToJsonString(), Is.EqualTo("[[3,5,8],[0,2,5]]"));
        Assert.That(mapped["sequences"]![0]!["parts"]!["0"]!.GetValue<string>(), Is.EqualTo("CCAACGTTGG"));
        Assert.That(mapped["sequences"]![1]!["parts"]!["0"]!.GetValue<string>(), Is.EqualTo("ACGTT"));
    }
}
