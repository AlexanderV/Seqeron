using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The documented resource limits of the program ports (audit A3-2): every dynamic program that allocates a matrix of the product of the sequence lengths refuses an input above its limit with ONE
/// kind of <see cref="ArgumentException"/> (see <see cref="SequenceAligner.ResourceLimit"/>), computed in 64-bit arithmetic BEFORE anything is allocated (no <see cref="OutOfMemoryException"/>, no
/// 32-bit overflow, no killed process). The programs themselves have no such limit: they fail at the allocation (EMBOSS needle / water: "Died: Sequences too big. Try 'stretcher'" /
/// "Try 'matcher' or 'supermatcher'", measured with the installed 6.6.0 under a memory limit) or overflow a 32-bit index. The limits are public constants so that a caller can test an input beforehand.
/// </summary>
public static partial class SequenceAligner
{
    /// <summary><see cref="NeedleAlign"/>: cells (len1 x len2) of the three float matrices and the int compass matrix (16 bytes per cell): 2^27 cells = 2 GiB, about 11 585 x 11 585 residues.</summary>
    public const long NeedleMaxCells = 1L << 27;

    /// <summary><see cref="WaterAlign"/>: cells of the float path matrix and the int compass matrix (8 bytes per cell): 2^28 cells = 2 GiB, about 16 384 x 16 384 residues.</summary>
    public const long WaterMaxCells = 1L << 28;

    /// <summary><see cref="Supermatcher"/>: cells ((target length + 1) x <c>width</c>) of the banded float path and int compass matrices (8 bytes per cell): 2^28 cells = 2 GiB.</summary>
    public const long SupermatcherMaxPathCells = 1L << 28;

    /// <summary><see cref="MultipleAlignPoa"/> / <see cref="PoaRun"/>: cells (graph nodes x sequence length) of the score rows (8 bytes) and the two move matrices (1 byte each): 2^27 cells = 1.25 GiB.</summary>
    public const long PoaMaxMatrixCells = 1L << 27;

    /// <summary><see cref="MultipleAlignProbCons"/>: cells ((columns of profile 1 + 1) x (columns of profile 2 + 1)) of the posterior matrix and the dynamic program of the profile alignment: 2^26 cells (about 1 GiB with the DP).</summary>
    public const long ProbConsProfileMaxCells = 1L << 26;

    /// <summary><see cref="MultipleAlignMafft"/>: cells (len1 x len2) of the pair dynamic program <c>G__align11</c> (an <c>int</c> trace-back cell plus row arrays): 2^28 cells = 1 GiB for the traceback.</summary>
    public const long MafftMaxPairCells = 1L << 28;

    /// <summary>
    /// The one exception of a resource limit: "&lt;program&gt; &lt;what it keeps&gt;: &lt;size&gt; exceeds the limit of &lt;limit&gt; (&lt;constant&gt;). &lt;what the program itself does&gt;; &lt;what to use instead&gt;."
    /// </summary>
    internal static ArgumentException ResourceLimit(string program, string keeps, string size, long limit, string constant, string programNote, string alternative, string? paramName = null) =>
        new($"{program} {keeps}: {size} exceeds the limit of {limit:N0} ({constant}). {programNote}; {alternative}.", paramName);

    private static void RequireCells(long cells, long limit, string constant, string program, string keeps, string size, string programNote, string alternative)
    {
        if (cells > limit)
            throw ResourceLimit(program, keeps, size, limit, nameof(SequenceAligner) + "." + constant, programNote, alternative);
    }

    /// <summary>A needle / water / supermatcher / needleall penalty: a finite number inside the program's range (NaN passes a plain <c>x &lt; min || x &gt; max</c> test).</summary>
    internal static void RequireEmbossRange(double value, double min, double max, string parameter, string option, string program)
    {
        if (!(value >= min && value <= max))
            throw new ArgumentOutOfRangeException(parameter, value, $"{option} is a finite number between {min:0.0} and {max:0.0} in EMBOSS {program}.");
    }

    /// <summary>Integer-range contract of the stretcher (32-bit <c>int</c> program): the score of any path, <c>max(|substitution|, gapOpen + gapExtend) x (len1 + len2 + 2)</c>, must stay below 2^28 (as <see cref="PairScoring.EnsureRange"/>).</summary>
    internal static void RequireStretcherRange(SubstitutionMatrix matrix, long gapOpen, long gapExtend, long length1, long length2)
    {
        long maxAbs = Math.Max(matrix.MaxAbsScore, gapOpen + gapExtend);
        if (maxAbs * (length1 + length2 + 2) >= 1L << 28)
            throw new ArgumentException($"Score range too large for the integer dynamic program of EMBOSS stretcher (max |step| {maxAbs}, sequence lengths {length1} and {length2}; the program wraps its 32-bit int and prints \"stretcher CheckScore failed\"); use smaller gap penalties or shorter sequences.", "gapOpen");
    }
}
