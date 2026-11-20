using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-28d: EMBOSS dotmatcher; delegates to DotMatcher.Compute.
public partial class AlignmentTools
{
    [McpServerTool(Name = "dotmatcher", Title = "Pairwise — EMBOSS dotmatcher dot plot segments", ReadOnly = true)]
    [Description("EMBOSS 6.6.0 dotmatcher (600 DNA / protein pairs equal): the dot plot of two sequences by windowed substitution-matrix scores. On every diagonal a window of windowSize residue pairs slides along; a maximal run of windows whose summed matrix score reaches threshold is one segment from the first residue of its first window to the last residue of its last window. Returns the segments in the order the program draws them, 1-based inclusive positions on both sequences (the program's own scan treats the SECOND sequence as its first: the first diagonals start on the first residue of sequence2). matrix = BLOSUM62 | PAM250 | EDNAFULL | any EMBOSS matrix name; default EDNAFULL when sequence1 consists of nucleotide codes, EBLOSUM62 otherwise; a character that is not a letter of the matrix scores 0. Sequences not longer than the window give no segments. Input sequences are read as EMBOSS reads them (A4-3, one shared reader): only letters and '*' are kept; gap characters (- . ~ ?), '#', '+', '_', digits, blanks and every other character are dropped first, and positions, lengths and scores count the cleaned sequence (e.g. needle of AC-GTTAGC.TT* and ACGTTAGCTT: score 58, length 11, as the program).")]
    public static DotMatcherResult Dotmatcher(
        [Description("First sequence (the plot's y axis in the program).")] string sequence1,
        [Description("Second sequence (the plot's x axis).")] string sequence2,
        [Description("Window size (dotmatcher -windowsize, default 10).")] int windowSize = 10,
        [Description("Score a window must reach (dotmatcher -threshold, default 23).")] int threshold = 23,
        [Description(MatrixParamDescription)] string? matrix = null)
    {
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        var segments = global::Seqeron.Genomics.Alignment.DotMatcher.Compute(sequence1, sequence2, windowSize, threshold, ResolveMatrix(matrix));
        return new DotMatcherResult(segments.Select(s => new DotMatcherSegmentDto(s.Start1, s.Start2, s.End1, s.End2)).ToArray());
    }
}
