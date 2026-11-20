using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A5-7: the EMBOSS sequence-input qualifiers (-sbegin -send -sreverse -slower -supper -scircular -snucleotide -sprotein) of the pairwise tools; the tools delegate to SequenceAligner.NeedleProgram / WaterProgram / StretcherProgram.
public partial class AlignmentTools
{
    /// <summary>Description of the <c>sequenceQualifiers</c> parameter of the EMBOSS pairwise tools.</summary>
    internal const string SequenceQualifiersDescription = "EMBOSS sequence-input qualifiers as on the command line, with the input number: -sbegin1 3 -send1 -2 -sreverse1 -slower1 -supper1 -scircular1 -snucleotide1 -sprotein1 and the same with 2 for the second sequence "
        + "(begin / end: 1-based, negative = from the end; -sreverse: reverse-complement a nucleotide sequence, the range then refers to the original; -snucleotide / -sprotein force the type). Measured on EMBOSS 6.6.0: the alignment covers the range, and the formatted "
        + "alignment (default format 'pair') reports the positions of the ORIGINAL sequence, counting down for a reversed one. Omit for the whole sequences.";

    private static (EmbossSequence First, EmbossSequence Second) QualifiedPair(string sequence1, string sequence2, string? sequenceQualifiers, string[]? names)
    {
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        string n1 = names is { Length: > 0 } && !string.IsNullOrWhiteSpace(names[0]) ? names[0] : "a";
        string n2 = names is { Length: > 1 } && !string.IsNullOrWhiteSpace(names[1]) ? names[1] : "b";
        return (new EmbossSequence(n1, sequence1, "", EmbossSequenceQualifiers.Parse(sequenceQualifiers, "1")), new EmbossSequence(n2, sequence2, "", EmbossSequenceQualifiers.Parse(sequenceQualifiers, "2")));
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
}
