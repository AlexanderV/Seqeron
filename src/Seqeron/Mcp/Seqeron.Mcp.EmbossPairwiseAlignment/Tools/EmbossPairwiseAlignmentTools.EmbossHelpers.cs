// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the terms stated in LICENSE (this directory).
using Seqeron.Genomics.Analysis;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

public partial class EmbossPairwiseAlignmentTools
{
    // plain letters become a FASTA record; the text of a FASTA, GenBank or EMBL file is the file
    private static string DiffseqInput(string name, string text)
    {
        string t = text.TrimStart();
        if (t.StartsWith('>') || t.StartsWith("LOCUS ", StringComparison.Ordinal) || t.StartsWith("ID   ", StringComparison.Ordinal))
            return text;
        return ">" + name + "\n" + text + "\n";
    }

    private static EmbossSequence[] WithQualifiers(EmbossSequence[] sequences, string program, string? sequenceQualifiers, int input)
    {
        if (sequenceQualifiers is null)
            return sequences;
        var q = EmbossSequenceQualifiers.ParseInputs(program, sequenceQualifiers)[input];
        return [.. sequences.Select(s => s with { Qualifiers = q })];
    }

    private static EmbossSequence[] ToEmbossSequences(string[] sequences, string[]? names, string prefix, string paramName)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        if (sequences.Length == 0)
            throw new ArgumentException("At least one sequence is required.", paramName);
        var list = new EmbossSequence[sequences.Length];
        for (int i = 0; i < list.Length; i++)
        {
            RequireSequence(sequences[i], paramName);
            string name = names is not null && i < names.Length && !string.IsNullOrWhiteSpace(names[i]) ? names[i] : prefix + (i + 1).ToString(CultureInfo.InvariantCulture);
            list[i] = new EmbossSequence(name, sequences[i]);
        }

        return list;
    }

    private static EmbossProgramAlignResult ToProgramResult(EmbossProgramResult r, string? formatted = null)
    {
        var items = r.Alignments.Select(a => new PairAlignmentItem(a.First.Name, a.Second.Name, a.First.Row, a.Second.Row, a.First.Offset, a.Second.Offset, a.First.Offend, a.Second.Offend, a.Score)).ToArray();
        return new EmbossProgramAlignResult(items, items.Length, r.ErrorText, formatted, r.StandardError);
    }

    /// <summary>Description of the <c>sequenceQualifiers</c> parameter of the EMBOSS pairwise tools.</summary>
    internal const string SequenceQualifiersDescription = "EMBOSS sequence-input qualifiers as on the command line, with the input number: -sbegin1 3 -send1 -2 -sreverse1 -slower1 -supper1 -scircular1 -snucleotide1 -sprotein1 and the same with 2 for the second sequence "
        + "(begin / end: 1-based, negative = from the end; -sreverse: reverse-complement a nucleotide sequence, the range then refers to the original; -snucleotide / -sprotein force the type). Measured on EMBOSS 6.6.0: the alignment covers the range, and the formatted "
        + "alignment (default format 'pair') reports the positions of the ORIGINAL sequence, counting down for a reversed one. Omit for the whole sequences. " + EmbossSequenceQualifiers.CommandLineRules;

    private static (EmbossSequence First, EmbossSequence Second) QualifiedPair(string program, string sequence1, string sequence2, string? sequenceQualifiers, string[]? names)
    {
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        string n1 = names is { Length: > 0 } && !string.IsNullOrWhiteSpace(names[0]) ? names[0] : "a";
        string n2 = names is { Length: > 1 } && !string.IsNullOrWhiteSpace(names[1]) ? names[1] : "b";
        var q = EmbossSequenceQualifiers.ParseInputs(program, sequenceQualifiers);
        return (new EmbossSequence(n1, sequence1, "", q[0]), new EmbossSequence(n2, sequence2, "", q[1]));
    }

    private static AlignmentResultDto QualifiedPairResult(string program, EmbossProgramResult result, string? format, double gapOpen, double gapExtend, bool global, bool integerScores)
    {
        if (result.Alignments.Count == 0)
            return new AlignmentResultDto("", "", 0, global ? "Global" : "Local", 0, 0, 0, 0);
        var pair = result.Alignments[0];
        var options = new EmbossAlignOptions
        {
            Program = program,
            MatrixName = result.MatrixName,
            Matrix = SubstitutionMatrix.Get(result.MatrixName),
            Global = global,
            IntegerScores = integerScores,
            StatisticsSetByProgram = integerScores,
            GapPenalty = gapOpen,
            ExtendPenalty = gapExtend,
        };
        string f = string.IsNullOrWhiteSpace(format) ? "pair" : format.Trim();
        if (!EmbossAlignmentWriter.TryParseFormat(f.StartsWith("emboss:", StringComparison.OrdinalIgnoreCase) ? f[7..] : f, out var fmt))
            throw new ArgumentException($"Unknown EMBOSS alignment format '{format}'.", nameof(format));
        string text = EmbossAlignmentWriter.Format(fmt, [pair], options);
        int residues1 = pair.First.Row.Count(c => c != '-'), residues2 = pair.Second.Row.Count(c => c != '-');
        return new AlignmentResultDto(pair.First.Row, pair.Second.Row, (int)Math.Round(pair.Score, MidpointRounding.AwayFromZero), global ? "Global" : "Local",
            pair.First.Offset, pair.Second.Offset, pair.First.Offset + residues1, pair.Second.Offset + residues2, pair.Score, text);
    }

    /// <summary>Description of the <c>sequenceQualifiers</c> parameter of <c>merger</c>.</summary>
    internal const string MergerQualifiersDescription = "EMBOSS sequence-input qualifiers as on the command line with the input number: -sbegin1 3 -send1 -2 -sreverse1 -slower1 -supper1 -scircular1 -snucleotide1 -sprotein1 and the same with 2 for the second sequence. "
        + "Measured on EMBOSS 6.6.0: merger does not trim its sequences - -sbegin / -send have no effect on it (the whole sequences are merged and the conflict positions count them); -sreverse reverse-complements a sequence, -snucleotide / -sprotein force the type "
        + "(a protein pair is merged with BLOSUM62 unless a matrix is given). Omit to use the sequences as given. " + EmbossSequenceQualifiers.CommandLineRules;

    private static double?[][] ToJagged(double[,] matrix)
    {
        int n = matrix.GetLength(0), m = matrix.GetLength(1);
        var result = new double?[n][];
        for (int i = 0; i < n; i++)
        {
            result[i] = new double?[m];
            for (int j = 0; j < m; j++)
                result[i][j] = double.IsFinite(matrix[i, j]) ? matrix[i, j] : null;
        }
        return result;
    }

    private static global::Seqeron.Genomics.Analysis.ConsensusResidueType ParseConsensusResidueType(string? residueType) =>
        (residueType ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "auto" => global::Seqeron.Genomics.Analysis.ConsensusResidueType.Auto,
            "nucleotide" => global::Seqeron.Genomics.Analysis.ConsensusResidueType.Nucleotide,
            "protein" => global::Seqeron.Genomics.Analysis.ConsensusResidueType.Protein,
            _ => throw new ArgumentException("residueType must be 'nucleotide', 'protein' or 'auto'", nameof(residueType)),
        };

    private static AlignmentResultDto StretcherAlignCore(
        [Description("First sequence (-asequence).")] string sequence1,
        [Description("Second sequence (-bsequence).")] string sequence2,
        [Description("Optional EMBOSS data matrix (-datafile): EBLOSUM62, EPAM250, EDNAFULL, ... or the NCBI name. Omit for EDNAFULL (nucleotide) / EBLOSUM62 (protein).")] string? matrix = null,
        [Description("Gap opening penalty (-gapopen, integer >= 0): the cost of a gap of length 1. Omit for 16 (nucleotide) / 12 (protein).")] int? gapOpen = null,
        [Description("Gap extension penalty (-gapextend, integer >= 0): the cost of each further gap position. Omit for 4 (nucleotide) / 2 (protein).")] int? gapExtend = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        return ToAlignmentDto(global::Seqeron.Genomics.Emboss.PairwiseAlignment.EmbossPairwisePrograms.StretcherAlign(sequence1, sequence2, ResolveMatrix(matrix), gapOpen, gapExtend));
    }

    // the length of an input sequence as the EMBOSS reader counts it (A4-3: gap characters, digits, blanks and every other non-letter except * are dropped before the program works)
    private static int EmbossLength(string sequence) => EmbossSequenceReader.Clean(sequence).Length;
}
