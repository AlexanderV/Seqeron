using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A3-13: EMBOSS merger (merge two overlapping nucleotide sequences); delegates to SequenceAligner.Merger.
public partial class AlignmentTools
{
    /// <summary>Description of the <c>sequenceQualifiers</c> parameter of <c>merger</c>.</summary>
    internal const string MergerQualifiersDescription = "EMBOSS sequence-input qualifiers as on the command line with the input number: -sbegin1 3 -send1 -2 -sreverse1 -slower1 -supper1 -scircular1 -snucleotide1 -sprotein1 and the same with 2 for the second sequence. "
        + "Measured on EMBOSS 6.6.0: merger does not trim its sequences - -sbegin / -send have no effect on it (the whole sequences are merged and the conflict positions count them); -sreverse reverse-complements a sequence, -snucleotide / -sprotein force the type "
        + "(a protein pair is merged with BLOSUM62 unless a matrix is given). Omit to use the sequences as given.";

    [McpServerTool(Name = "merger", Title = "Pairwise — EMBOSS merger (merge two overlapping sequences)", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 merger: merges two overlapping nucleotide sequences into one. A global alignment with free end gaps (the program's embAlignPathCalc / embAlignWalkNWMatrix: single-precision path matrices, a gap of length k costs gapOpen + (k - 1) * gapExtend, EDNAFULL scores or the named matrix) is walked column by column: identical bases are copied, an N takes the other sequence's base, and a mismatch or gap takes the base of the sequence with the better quality around the position (2 points per A/C/G/T/U, 1 per other IUPAC code, 0 per N and off the end, windows of 5, 10, 20 on both sides, equal = sequence A) and is listed in conflicts (positions 1-based, '.' = gap). The unaligned head / tail of the sequence that sticks out is appended; mergedSequence is lower case except the bases chosen at conflicts and next to an N (upper case), as the program writes it. Verified against merger 6.6.0 on 450 random overlapping read pairs (point errors, indels, N / IUPAC codes, lower case, six gap-penalty settings): merged sequence, conflict lines and score equal on all. Note the program's own quirk: a tail of the longer sequence that overhangs after the alignment is appended twice when both sequences have unaligned residues at the end - unreachable, since the walk always ends on a matrix border. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static MergerResultDto Merger(
        [Description("First sequence (merger -asequence).")] string sequenceA,
        [Description("Second sequence (merger -bsequence).")] string sequenceB,
        [Description("Gap opening penalty (-gapopen, 0..100, default 50): a gap of length 1.")] double gapOpen = 50,
        [Description("Gap extension penalty (-gapextend, 0..10, default 5).")] double gapExtend = 5,
        [Description(MatrixParamDescription)] string? matrix = null,
        [Description(MergerQualifiersDescription)] string? sequenceQualifiers = null,
        CancellationToken cancellationToken = default)
    {
        RequireSequence(sequenceA, nameof(sequenceA));
        RequireSequence(sequenceB, nameof(sequenceB));
        var r = sequenceQualifiers is null
            ? global::Seqeron.Genomics.Alignment.SequenceAligner.Merger(sequenceA, sequenceB, ResolveMatrix(matrix), gapOpen, gapExtend, cancellationToken)
            : global::Seqeron.Genomics.Alignment.SequenceAligner.MergerProgram(
                new global::Seqeron.Genomics.Infrastructure.EmbossSequence("a", sequenceA, "", global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.Parse(sequenceQualifiers, "1")),
                new global::Seqeron.Genomics.Infrastructure.EmbossSequence("b", sequenceB, "", global::Seqeron.Genomics.Infrastructure.EmbossSequenceQualifiers.Parse(sequenceQualifiers, "2")),
                ResolveMatrix(matrix), gapOpen, gapExtend, cancellationToken);
        return new MergerResultDto(
            r.MergedSequence, r.Score, r.StartA, r.StartB, r.AlignedA, r.AlignedB,
            r.Conflicts.Select(c => new MergerConflictDto(c.PositionA, c.BaseA.ToString(), c.PositionB, c.BaseB.ToString(), c.Using.ToString())).ToArray(),
            r.Warnings.ToArray());
    }
}
