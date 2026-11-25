using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// B09 WP67 (F308, R14-13): the MCP tools <c>alignment_summary_info</c> (Bio.Align.AlignInfo.SummaryInfo of Biopython 1.79) and <c>msa_conservation</c> (Capra and Singh 2007) and the extended columns of
/// <c>multiple_alignment_statistics</c>. Expected values: the test suite of Biopython 1.79 (<c>Tests/test_AlignInfo.py</c>) and an independent implementation of the formulas of the paper.
/// </summary>
[TestFixture]
public class AlignInfoToolsTests
{
    private const string Proteins = "[{\"id\": \"ID001\", \"seq\": \"MHQAIFIYQIGYP*LKSGYIQSIRSPEYDNW-\"}, {\"id\": \"ID002\", \"seq\": \"MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*\"}, {\"id\": \"ID003\", \"seq\": \"MHQAIFIYQIGYPYLKSGYIQSIRSPEYDNW*\"}]";
    private const string Dna = "[{\"id\": \"test1\", \"seq\": \"ACGTCGCG\"}, {\"id\": \"test2\", \"seq\": \"GGGGCCCC\"}, {\"id\": \"test3\", \"seq\": \"AAACACAC\"}]";

    [Test]
    public void SummaryInfo_Consensus_EqualsBiopythonTests()
    {
        Assert.Multiple(() =>
        {
            Assert.That(JsonNode.Parse(AlignmentSummaryInfoOperationOf("dumb_consensus", Proteins).ResultJson)!.GetValue<string>(), Is.EqualTo("MHQAIFIYQIGYXXLKSGYIQSIRSPEYDNW*"));
            Assert.That(JsonNode.Parse(AlignmentSummaryInfoOperationOf("gap_consensus", Proteins).ResultJson)!.GetValue<string>(), Is.EqualTo("MHXXIFIYQIGYXXLKSGYIQSIRSPEYXNWX"));
            Assert.That(JsonNode.Parse(AlignmentSummaryInfoOperationOf("dumb_consensus", Dna, ambiguous: "N").ResultJson)!.GetValue<string>(), Is.EqualTo("NNNNNNNN"));
        });
    }

    private static AlignmentSummaryInfoResult AlignmentSummaryInfoOperationOf(string operation, string records, string ambiguous = "X", string? charsToIgnoreJson = null, string? eFreq = null, double logBase = 2, double pseudoCount = 0, string? axisSeq = null, string? letters = null)
        => AnalysisTools.AlignmentSummaryInfoOperation(operation, records, ambiguous: ambiguous, charsToIgnoreJson: charsToIgnoreJson, eFreqTableJson: eFreq, logBase: logBase, pseudoCount: pseudoCount, axisSeq: axisSeq, letters: letters);

    [Test]
    public void SummaryInfo_PssmAndInformationContent_EqualBiopythonTests()
    {
        var pssm = AlignmentSummaryInfoOperationOf("pos_specific_score_matrix", Dna, charsToIgnoreJson: "[\"-\"]", axisSeq: "NNNNNNNN");
        Assert.That(pssm.Text, Is.EqualTo("    A   C   G   T\nN  2.0 0.0 1.0 0.0\nN  1.0 1.0 1.0 0.0\nN  1.0 0.0 2.0 0.0\nN  0.0 1.0 1.0 1.0\nN  1.0 2.0 0.0 0.0\nN  0.0 2.0 1.0 0.0\nN  1.0 2.0 0.0 0.0\nN  0.0 2.0 1.0 0.0\n"));
        var ic = JsonNode.Parse(AlignmentSummaryInfoOperationOf("information_content", Dna, charsToIgnoreJson: "[\"-\"]", eFreq: "{\"A\": 0.25, \"G\": 0.25, \"T\": 0.25, \"C\": 0.25}").ResultJson)!;
        Assert.That(ic["value"]!.GetValue<double>(), Is.EqualTo(7.32029999423075).Within(1e-6));
        Assert.That(ic["ic_vector"]!.AsArray(), Has.Count.EqualTo(8));
        var uniform = string.Join(", ", "ACDEFGHIKLMNPQRSTVWY".Select(c => $"\"{c}\": 0.05"));
        var protein = JsonNode.Parse(AlignmentSummaryInfoOperationOf("information_content", Proteins, charsToIgnoreJson: "[\"-\", \"*\"]", eFreq: "{" + uniform + "}").ResultJson)!;
        Assert.That(protein["value"]!.GetValue<double>(), Is.EqualTo(133.061475107).Within(1e-6));
    }

    [Test]
    public void SummaryInfo_ReplacementPrintAndErrors()
    {
        var weighted = "[{\"id\": \"a\", \"seq\": \"GTATC\", \"annotations\": {\"weight\": 0.5}}, {\"id\": \"b\", \"seq\": \"AT--C\", \"annotations\": {\"weight\": 0.8}}, {\"id\": \"c\", \"seq\": \"CTGTC\", \"annotations\": {\"weight\": 1.0}}]";
        var rep = JsonNode.Parse(AlignmentSummaryInfoOperationOf("replacement_dictionary", weighted, letters: "ACGT").ResultJson)!.AsArray();
        Assert.Multiple(() =>
        {
            Assert.That(rep, Has.Count.EqualTo(16));
            Assert.That(rep.First(x => x![0]!.GetValue<string>() == "G" && x[1]!.GetValue<string>() == "A")![2]!.GetValue<double>(), Is.EqualTo(0.4).Within(1e-12));
            Assert.That(rep.First(x => x![0]!.GetValue<string>() == "T" && x[1]!.GetValue<string>() == "T")![2]!.GetValue<double>(), Is.EqualTo(2.2).Within(1e-12));
            Assert.That(AlignmentSummaryInfoOperationOf("print_info_content", "[{\"id\": \"a\", \"seq\": \"AC\"}, {\"id\": \"b\", \"seq\": \"AG\"}]", eFreq: "{\"A\": 0.25, \"C\": 0.25, \"G\": 0.25}").Text, Is.EqualTo("0 A 2.000\n1 C 1.000\n"));
            Assert.That(JsonNode.Parse(AlignmentSummaryInfoOperationOf("get_column", Dna).ResultJson)!.GetValue<string>(), Is.EqualTo("AGA"));
            Assert.Throws<InvalidOperationException>(() => AlignmentSummaryInfoOperationOf("information_content", Dna), "no expected frequencies: Biopython's TypeError");
            Assert.Throws<ArgumentException>(() => AlignmentSummaryInfoOperationOf("nope", Dna));
        });
    }

    [Test]
    public void MsaConservation_EqualsTheFormulasOfThePaper()
    {
        // columns AAAAA, AAAAC, ACDEF, AAA--, ARNDC: values of the independent implementation (unit tests of ConservationScores)
        string[] rows = ["AAAAA", "AAAAC", "ACDEF", "AAA--", "ARNDC"];
        string[] sequences = [.. Enumerable.Range(0, 5).Select(i => new string([.. rows.Select(r => r[i])]))];
        // the alignment is the transpose: five rows of the five columns above
        var r = AnalysisTools.MsaConservationScores(sequences, "js_divergence", gapPenalty: true);
        Assert.Multiple(() =>
        {
            Assert.That(r.Columns, Is.EqualTo(rows));
            Assert.That(r.Scores[0], Is.EqualTo(0.7990563231997274).Within(1e-13));
            Assert.That(r.Scores[1], Is.EqualTo(0.7560746838749699).Within(1e-13));
            Assert.That(r.Scores[2], Is.EqualTo(0.5511114399832872).Within(1e-13));
            Assert.That(r.Scores[3], Is.EqualTo(0.47943223104688665).Within(1e-13), "AAA--: the divergence 0.79905 times the 3 / 5 non-gap rows");
            Assert.That(r.Scores[4], Is.EqualTo(0.5637373791192319).Within(1e-13));
            Assert.That(AnalysisTools.MsaConservationScores(sequences, "shannon_entropy").Scores[0], Is.EqualTo(0.9999953455735638).Within(1e-13));
            Assert.That(AnalysisTools.MsaConservationScores(sequences, "shannon_entropy").Scores[3], Is.EqualTo(0.34909804566659847).Within(1e-13));
            Assert.That(AnalysisTools.MsaConservationScores(sequences, "sum_of_pairs", gapPenalty: false).Scores[0], Is.EqualTo(4.0).Within(1e-13), "five A: every pair A/A = 4");
            Assert.That(AnalysisTools.MsaConservationScores(sequences, sequenceWeights: true).SequenceWeights.Sum(), Is.EqualTo(1.0).Within(1e-12));
            Assert.Throws<ArgumentException>(() => AnalysisTools.MsaConservationScores(sequences, "nope"));
        });
    }
}
