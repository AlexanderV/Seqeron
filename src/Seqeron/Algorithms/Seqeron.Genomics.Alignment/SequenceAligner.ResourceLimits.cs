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

    /// <summary><see cref="MultipleAlignProbCons"/>: cells ((columns of profile 1 + 1) x (columns of profile 2 + 1)) of the posterior matrix and the dynamic program of the profile alignment: 2^26 cells (about 1 GiB with the DP).</summary>
    public const long ProbConsProfileMaxCells = 1L << 26;

    /// <summary><see cref="MultipleAlignMafft"/>: cells (len1 x len2) of the pair dynamic program <c>G__align11</c> (an <c>int</c> trace-back cell plus row arrays): 2^28 cells = 1 GiB for the traceback.</summary>
    public const long MafftMaxPairCells = 1L << 28;

    /// <summary><see cref="Minimap2Chaining.CollectAnchors"/> / <see cref="Minimap2Chaining.Chain"/>: anchors (20 bytes each, plus 36 bytes per anchor in the chaining arrays and the sort copies): 2^23 = 8 388 608 anchors, about 0.5 GiB at the peak. minimap2 has no such limit (a 100 kb poly-A query against itself is 10^10 anchors and the program dies at the allocation); it drops frequent minimizers instead (<c>mid_occ</c>, <c>mm_seed_select</c>: the registered A4-7).</summary>
    public const long MinimapMaxAnchors = 1L << 23;

    /// <summary><see cref="DotMatcher.Compute"/>: segments of the plot (16 bytes each; an alternating window score gives one segment per few cells, up to cells / 4 of them): 2^22 = 4 194 304 segments = 64 MiB. The program streams its segments to a graph file and has no such limit.</summary>
    public const long DotMatcherMaxSegments = 1L << 22;

    /// <summary><see cref="NucmerClusters.Cluster"/>: matches of ONE cluster (the chain dynamic program of <c>Process_Cluster</c> is O(c^2) per extracted chain and so is the program's: 100 000 matches are about 5 x 10^9 steps, a minute): 100 000. mgaps has no such limit; split the input by query sequence or raise the minimum match length.</summary>
    public const long NucmerMaxClusterMatches = 100_000;

    /// <summary><c>GenomicAnalyzer.MegaMerge</c> / <c>GenomicAnalyzer.SeqMatchAll</c>: maximal word matches of one sequence pair that are listed (12 bytes in the engine's list, 20 more in the result): 2^22 = 4 194 304. The programs keep every match in memory and have no such limit.</summary>
    public const long WordMatchMaxMatches = 1L << 22;

    /// <summary>The K-vertex list of the Huang-Miller SIM local aligner (the T-Coffee <c>lalign_id_pair</c> port and EMBOSS <c>matcher</c>): the K-vertex list holds at most this many candidate start points (a vertex is about 56 bytes): 2^22 = 4 194 304. A larger <c>maxAlignments</c> is allocated lazily, so <c>int.MaxValue</c> asks for "every alignment" exactly as the Waterman-Eggert order accepts it; the program calls <c>calloc(K)</c> up front and dies for such a K.</summary>
    public const long TopNSimMaxVertices = 1L << 22;

    /// <summary>
    /// The one exception of a resource limit: "&lt;program&gt; &lt;what it keeps&gt;: &lt;size&gt; exceeds the limit of &lt;limit&gt; (&lt;constant&gt;). &lt;what the program itself does&gt;; &lt;what to use instead&gt;."
    /// </summary>
    public static ArgumentException ResourceLimit(string program, string keeps, string size, long limit, string constant, string programNote, string alternative, string? paramName = null) =>
        new(string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{program} {keeps}: {size} exceeds the limit of {limit:N0} ({constant}). {programNote}; {alternative}."), paramName);

    public static void RequireCells(long cells, long limit, string constant, string program, string keeps, string size, string programNote, string alternative, string? owner = null)
    {
        if (cells > limit)
            throw ResourceLimit(program, keeps, size, limit, (owner ?? nameof(SequenceAligner)) + "." + constant, programNote, alternative);
    }

    /// <summary>A needle / water / supermatcher / needleall penalty: a finite number inside the program's range (NaN passes a plain <c>x &lt; min || x &gt; max</c> test).</summary>
    public static void RequireEmbossRange(double value, double min, double max, string parameter, string option, string program)
    {
        if (!(value >= min && value <= max))
            throw new ArgumentOutOfRangeException(parameter, value, $"{option} is a finite number between {min:0.0} and {max:0.0} in EMBOSS {program}.");
    }

    /// <summary>
    /// An ACD <c>integer</c> option as <c>acdSetInt</c> reads it (audit A6-3, measured on matcher 6.6.0): the text must be a <c>strtol</c> integer ("1.5" / "1e2" / "abc": <c>Invalid integer value</c>, the program dies), a value beyond 64 bits is invalid too, a value
    /// beyond the 32-bit range is saturated, and a value below the option's minimum is RESET to the minimum with the warning <c>integer value out of range -1 less than (reset to) 0</c> (the run goes on).
    /// </summary>
    /// <exception cref="ArgumentException"><paramref name="value"/> is not an integer the program accepts (NaN, infinite, fractional, beyond 64 bits).</exception>
    public static int EmbossInteger(double value, int minimum, string parameter, string option, System.Text.StringBuilder warnings)
    {
        if (!double.IsFinite(value) || value != Math.Floor(value) || Math.Abs(value) >= 9.2233720368547758E+18)
            throw new ArgumentException($"Invalid integer value '{value.ToString("R", System.Globalization.CultureInfo.InvariantCulture)}' for {option}: an integer is required.", parameter);
        int v = (int)Math.Clamp(value, int.MinValue, int.MaxValue);
        if (v >= minimum)
            return v;
        warnings.Append(System.Globalization.CultureInfo.InvariantCulture, $"Warning: integer value out of range {v} less than (reset to) {minimum}\n");
        return minimum;
    }

    /// <summary>Integer-range contract of the stretcher (32-bit <c>int</c> program): the score of any path, <c>max(|substitution|, gapOpen + gapExtend) x (len1 + len2 + 2)</c>, must stay below 2^28 (as <see cref="PairScoring.EnsureRange"/>).</summary>
    public static void RequireStretcherRange(SubstitutionMatrix matrix, long gapOpen, long gapExtend, long length1, long length2)
    {
        long maxAbs = Math.Max(matrix.MaxAbsScore, gapOpen + gapExtend);
        if (maxAbs * (length1 + length2 + 2) >= 1L << 28)
            throw new ArgumentException($"Score range too large for the integer dynamic program of EMBOSS stretcher (max |step| {maxAbs}, sequence lengths {length1} and {length2}; the program wraps its 32-bit int and prints \"stretcher CheckScore failed\"); use smaller gap penalties or shorter sequences.", "gapOpen");
    }
}
