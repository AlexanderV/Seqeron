namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// Scalar scoring parameters for sequence alignment (one score for every identical pair, one for every
/// non-identical pair, and two gap parameters).
/// </summary>
/// <remarks>
/// <para>
/// Gap convention. The affine aligners (<c>GlobalAlignAffine</c>, <c>LocalAlignAffine</c>,
/// <c>SemiGlobalAlignAffine</c>) charge a gap of length k as <see cref="GapOpen"/> + k·<see cref="GapExtend"/>
/// (Gotoh 1982 w(k) = v + u·k; the NCBI BLAST+ convention, "existence 5, extension 2" = 5 + 2k; Biopython
/// <c>open_gap_score</c> = <see cref="GapOpen"/> + <see cref="GapExtend"/>). The linear aligners
/// (<c>GlobalAlign</c>, <c>LocalAlign</c>, <c>SemiGlobalAlign</c>, the multiple-alignment methods) use only
/// <see cref="GapExtend"/> per gap position. EMBOSS needle/water instead charge gapopen + (k−1)·gapextend for
/// the same two numbers, so EMBOSS "-gapopen 10 -gapextend 0.5" corresponds to GapOpen = −9.5 here (not
/// representable with integers).
/// </para>
/// <para>
/// Contract (identical to <see cref="MatrixScoring"/>; the record itself validates nothing, the aligners do): <see cref="GapOpen"/> and
/// <see cref="GapExtend"/> are penalties, i.e. at most 0 (a positive gap score is rejected with <see cref="ArgumentException"/>) and greater than
/// <c>int.MinValue</c>; the scores must fit the integer dynamic programs: the score of any path through sequences of m and n residues,
/// <c>max(|Match|, |Mismatch|, |GapOpen| + |GapExtend|) · (m + n + 2)</c>, must stay below 2^28, otherwise the aligner throws
/// <see cref="ArgumentException"/> "Score range too large for the integer dynamic program" instead of returning a silently wrapped score
/// (the reported <see cref="AlignmentResult.Score"/> is a 32-bit integer). The multiple-alignment sum-of-pairs totals are accumulated in 64 bits
/// and refused when they do not fit the reported 32-bit total.
/// </para>
/// <para>
/// The presets live on <c>SequenceAligner</c> (<c>SimpleDna</c>, <c>BlastDna</c>, <c>HighIdentityDna</c>). For
/// residue-pair specific scores (IUPAC DNA, protein) use <see cref="SubstitutionMatrix"/>.
/// </para>
/// </remarks>
public sealed record ScoringMatrix(
    int Match,
    int Mismatch,
    int GapOpen,
    int GapExtend);

/// <summary>
/// Residue-pair scoring for pairwise alignment: a <see cref="SubstitutionMatrix"/> (BLOSUM62, PAM250, EDNAFULL, ...)
/// plus <b>real-valued</b> gap scores. Accepted by the <c>GlobalAlign</c>/<c>LocalAlign</c>/<c>SemiGlobalAlign</c>
/// (and <c>*Affine</c>) overloads of <c>SequenceAligner</c>; protein, IUPAC DNA and any alphabet of the matrix work.
/// </summary>
/// <remarks>
/// <para>
/// Gap convention = <see cref="ScoringMatrix"/>: a gap is a <b>score</b> (negative = penalty). The affine aligners
/// charge a gap of length k as <see cref="GapOpen"/> + k·<see cref="GapExtend"/>; the linear aligners charge
/// <see cref="GapExtend"/> per gap position (<see cref="GapOpen"/> unused). EMBOSS needle/water and Biopython
/// differ only in how the same numbers are named: EMBOSS <c>-gapopen G -gapextend E</c> (positive penalties; a gap of
/// length k costs G + (k-1)·E) is <c>GapOpen = -(G - E)</c>, <c>GapExtend = -E</c> — use
/// <see cref="Emboss(SubstitutionMatrix, double, double)"/>. Biopython <c>open_gap_score</c> = <c>GapOpen + GapExtend</c>,
/// <c>extend_gap_score</c> = <c>GapExtend</c>.
/// </para>
/// <para>
/// Gap scores are resolved exactly to 1e-6 (more decimal places are rejected); internally every score is an
/// exact multiple of 1/scale (scale = smallest integer making the gaps integral: 2 for the EMBOSS default 0.5, 10 for 0.1),
/// so ties and tracebacks are exact (no floating-point drift) and the integer dynamic programs are reused unchanged.
/// </para>
/// </remarks>
/// <param name="Matrix">Residue-pair scores; every sequence symbol must be in its alphabet (case-insensitive).</param>
/// <param name="GapOpen">Score added once per gap run (affine aligners only); at most 0.</param>
/// <param name="GapExtend">Score per gap position; at most 0.</param>
public sealed record MatrixScoring(SubstitutionMatrix Matrix, double GapOpen, double GapExtend)
{
    /// <summary>
    /// Scoring from EMBOSS needle/water penalties (positive numbers): <paramref name="gapOpenPenalty"/> is the cost of a gap
    /// of length 1, every further position costs <paramref name="gapExtendPenalty"/> (needle/water defaults: 10 and 0.5).
    /// </summary>
    public static MatrixScoring Emboss(SubstitutionMatrix matrix, double gapOpenPenalty = 10, double gapExtendPenalty = 0.5)
    {
        if (!double.IsFinite(gapOpenPenalty))
            throw new ArgumentOutOfRangeException(nameof(gapOpenPenalty), gapOpenPenalty, "The gap opening penalty must be a finite number (NaN / Infinity are rejected).");
        if (!double.IsFinite(gapExtendPenalty))
            throw new ArgumentOutOfRangeException(nameof(gapExtendPenalty), gapExtendPenalty, "The gap extension penalty must be a finite number (NaN / Infinity are rejected).");
        return new(matrix, -(gapOpenPenalty - gapExtendPenalty), -gapExtendPenalty);
    }
}

/// <summary>
/// Type of alignment performed.
/// </summary>
public enum AlignmentType
{
    /// <summary>Global (Needleman-Wunsch) alignment: both sequences end to end.</summary>
    Global,
    /// <summary>Local (Smith-Waterman) alignment: the best-scoring pair of subsequences.</summary>
    Local,
    /// <summary>Semi-global (ends-free, fitting / overlap) alignment: end gaps are not penalised.</summary>
    SemiGlobal
}

/// <summary>
/// End-gap policy of an <b>ends-free</b> (semi-global) alignment: which of the four alignment ends may overhang
/// without penalty. "Start"/"End" are the left/right ends of the alignment, "1"/"2" the sequence whose residues
/// overhang: <see cref="FreeStart1"/> = the leading residues of sequence 1 that precede the alignment are free (they are
/// shown against '-' in sequence 2), and so on. Flags not set are charged like any other gap
/// (linear d per position; affine GapOpen + k·GapExtend). No flag = global (Needleman-Wunsch) alignment.
/// </summary>
/// <remarks>
/// Equivalent settings in the reference tools (sequence 1 = Biopython target = parasail query/s1, sequence 2 = Biopython
/// query = parasail database/s2; a free overhang of sequence 1 is a gap in sequence 2, i.e. a Biopython <i>deletion</i>):
/// <see cref="FreeStart1"/> = Biopython <c>left_deletion_score = 0</c> = parasail <c>_qb</c>; <see cref="FreeEnd1"/> =
/// <c>right_deletion_score = 0</c> = <c>_qe</c>; <see cref="FreeStart2"/> = <c>left_insertion_score = 0</c> = <c>_db</c>;
/// <see cref="FreeEnd2"/> = <c>right_insertion_score = 0</c> = <c>_de</c>; <see cref="Fitting"/> = <c>sg_dx</c>
/// (<c>end_insertion_score = 0</c>); <see cref="All"/> = <c>sg</c> = EMBOSS needle's default (<c>-endweight N</c>).
/// </remarks>
[Flags]
public enum EndGaps
{
    /// <summary>All end gaps are penalised: global (Needleman-Wunsch) alignment.</summary>
    None = 0,

    /// <summary>Leading residues of sequence 1 before the alignment are free (Biopython <c>left_deletion_score = 0</c>, parasail <c>_qb</c>).</summary>
    FreeStart1 = 1,

    /// <summary>Trailing residues of sequence 1 after the alignment are free (Biopython <c>right_deletion_score = 0</c>, parasail <c>_qe</c>).</summary>
    FreeEnd1 = 2,

    /// <summary>Leading residues of sequence 2 before the alignment are free (Biopython <c>left_insertion_score = 0</c>, parasail <c>_db</c>).</summary>
    FreeStart2 = 4,

    /// <summary>Trailing residues of sequence 2 after the alignment are free (Biopython <c>right_insertion_score = 0</c>, parasail <c>_de</c>).</summary>
    FreeEnd2 = 8,

    /// <summary>Both flanks of sequence 2 are free: sequence 1 is fitted end to end into sequence 2 ("fitting"/"glocal"; the classic
    /// <c>SemiGlobalAlign</c>; Biopython <c>end_insertion_score = 0</c>, parasail <c>sg_dx</c>, Rosalind SIMS).</summary>
    Fitting = FreeStart2 | FreeEnd2,

    /// <summary>Both flanks of sequence 1 are free: sequence 2 is fitted end to end into sequence 1 (parasail <c>sg_qx</c>).</summary>
    FittingSequence2 = FreeStart1 | FreeEnd1,

    /// <summary>Overlap (suffix-prefix): a suffix of sequence 1 overlaps a prefix of sequence 2 — the leading part of sequence 1
    /// and the trailing part of sequence 2 are free (parasail <c>sg_qb_de</c>; sequence assembly "dovetail" overlap).</summary>
    OverlapSuffix1Prefix2 = FreeStart1 | FreeEnd2,

    /// <summary>Overlap (suffix-prefix), other orientation: a suffix of sequence 2 overlaps a prefix of sequence 1
    /// (parasail <c>sg_qe_db</c>).</summary>
    OverlapSuffix2Prefix1 = FreeEnd1 | FreeStart2,

    /// <summary>All four end gaps are free (parasail <c>sg</c>, Rosalind SMGB-like, EMBOSS needle's default <c>-endweight N</c>).</summary>
    All = FreeStart1 | FreeEnd1 | FreeStart2 | FreeEnd2,
}

/// <summary>
/// The score of one kind of gap in the Seqeron convention: a gap of length k scores <c>Open + k·Extend</c> in the affine aligners
/// (<see cref="Open"/> is added once per run), <c>k·Extend</c> in the linear aligners (<see cref="Open"/> is ignored there, as
/// <see cref="ScoringMatrix.GapOpen"/> is). Both at most 0 (penalties), finite.
/// </summary>
/// <param name="Open">Score added once per gap run (affine aligners only); at most 0.</param>
/// <param name="Extend">Score per gap position; at most 0.</param>
public readonly record struct GapScore(double Open, double Extend)
{
    /// <summary>A free gap (score 0 per position): an end that may overhang without penalty.</summary>
    public static GapScore Free => new(0, 0);

    /// <summary>
    /// A gap score from Biopython's convention (<c>open_gap_score</c> = score of the FIRST gap position, <c>extend_gap_score</c> = score of
    /// each further position): <c>Open = openGapScore − extendGapScore</c>, <c>Extend = extendGapScore</c>.
    /// </summary>
    public static GapScore FromBiopython(double openGapScore, double extendGapScore) => new(openGapScore - extendGapScore, extendGapScore);
}

/// <summary>
/// Gap scores that differ by gap direction and by position (Biopython <c>PairwiseAligner</c>'s twelve properties
/// <c>{target,query}_{internal,left,right}_{open,extend}_gap_score</c>): for each of the two gap directions an <b>interior</b>
/// score, a <b>left (start)</b> score (a gap run that begins the alignment) and a <b>right (end)</b> score (a gap run that ends it).
/// <list type="bullet">
/// <item>Direction 1 = a gap in sequence 2 opposite residues of sequence 1 (a Biopython <i>deletion</i>, i.e. a gap in the
/// <b>query</b>): <see cref="Interior1"/> = <c>query_internal_*</c>, <see cref="Start1"/> = <c>query_left_*</c>, <see cref="End1"/> = <c>query_right_*</c>.</item>
/// <item>Direction 2 = a gap in sequence 1 opposite residues of sequence 2 (a Biopython <i>insertion</i>, a gap in the <b>target</b>):
/// <see cref="Interior2"/> = <c>target_internal_*</c>, <see cref="Start2"/> = <c>target_left_*</c>, <see cref="End2"/> = <c>target_right_*</c>.</item>
/// </list>
/// A <c>null</c> member means "the scoring's own gap score" (<see cref="ScoringMatrix.GapOpen"/> / <see cref="ScoringMatrix.GapExtend"/>, or
/// <see cref="MatrixScoring.GapOpen"/> / <see cref="MatrixScoring.GapExtend"/>) — so a profile with every member null is the ordinary
/// global alignment. <see cref="GapScore.Free"/> in a start / end member is a free overhang: the <see cref="EndGaps"/> policies are the
/// profiles made only of free and default members (<see cref="FromEndGaps"/>). A gap run is a start run when it begins at the first
/// alignment column and an end run when it finishes at the last one (a column of the dynamic-programming matrix border); every
/// other run is interior. Local alignments never start or end in a gap: the profile has no effect on them.
/// </summary>
public sealed record GapScoreProfile(
    GapScore? Interior1 = null, GapScore? Start1 = null, GapScore? End1 = null,
    GapScore? Interior2 = null, GapScore? Start2 = null, GapScore? End2 = null)
{
    /// <summary>The profile of an <see cref="EndGaps"/> policy: <see cref="GapScore.Free"/> for each free end, the scoring's gap score otherwise.</summary>
    public static GapScoreProfile FromEndGaps(EndGaps free) => new(
        Start1: free.HasFlag(EndGaps.FreeStart1) ? GapScore.Free : null,
        End1: free.HasFlag(EndGaps.FreeEnd1) ? GapScore.Free : null,
        Start2: free.HasFlag(EndGaps.FreeStart2) ? GapScore.Free : null,
        End2: free.HasFlag(EndGaps.FreeEnd2) ? GapScore.Free : null);

    /// <summary>True when at least one start / end member is <see cref="GapScore.Free"/> (zero open and zero extend).</summary>
    public bool HasFreeEnd => IsFree(Start1) || IsFree(End1) || IsFree(Start2) || IsFree(End2);

    private static bool IsFree(GapScore? g) => g is { Open: 0, Extend: 0 };

    /// <summary>
    /// Resolves the name of a Biopython <c>PairwiseAligner</c> gap-score property (see <see cref="FromBiopython"/>) to the slots it sets: slot = direction (0 target / insertion, 1 query / deletion) * 3 + position (0 internal, 1 left, 2 right),
    /// whether it sets the open and / or the extend score, and its specificity (more specific names override more general ones).
    /// </summary>
    /// <exception cref="ArgumentException">An unknown name.</exception>
    public static (int Specificity, int[] Slots, bool SetOpen, bool SetExtend) ParseBiopythonName(string name)
    {
        ArgumentNullException.ThrowIfNull(name);
        var tokens = name.Split('_');
        if (tokens.Length < 2 || tokens[^1] != "score" || tokens.Distinct().Count() != tokens.Length)
            throw new ArgumentException($"Unknown Biopython gap-score property '{name}'.", nameof(name));
        var body = tokens[..^1];
        if (body.Any(t => t is not ("target" or "query" or "insertion" or "deletion" or "gap" or "internal" or "left" or "right" or "end" or "open" or "extend")))
            throw new ArgumentException($"Unknown Biopython gap-score property '{name}'.", nameof(name));
        bool target = body.Contains("target") || body.Contains("insertion");
        bool query = body.Contains("query") || body.Contains("deletion");
        bool gap = body.Contains("gap");
        if (!(target || query || gap) || (target && query))
            throw new ArgumentException($"Unknown Biopython gap-score property '{name}'.", nameof(name));
        int[] dirs = [0, 1];
        if (target) dirs = [0];
        else if (query) dirs = [1];
        var positions = new List<int>();
        if (body.Contains("internal")) positions.Add(0);
        if (body.Contains("left") || body.Contains("end")) positions.Add(1);
        if (body.Contains("right") || body.Contains("end")) positions.Add(2);
        if (positions.Count == 0) positions.AddRange([0, 1, 2]);
        return (body.Count(t => t is not "gap"), dirs.SelectMany(d => positions.Select(p => d * 3 + p)).Distinct().ToArray(), !body.Contains("extend"), !body.Contains("open"));
    }

    /// <summary>
    /// The profile of a Biopython <c>PairwiseAligner</c> from its gap-score properties (name → value, Biopython convention: the
    /// <c>open</c> scores are the score of the first gap position, the <c>extend</c> scores of each further one). Both name families of Biopython 1.88 are accepted:
    /// the current <c>[open_|extend_][internal_|left_|right_|end_]{insertion|deletion|gap}_score</c> (insertion = a gap in the target = sequence 1,
    /// deletion = a gap in the query) and the deprecated <c>[target_|query_][internal_|left_|right_|end_][open_|extend_]gap_score</c> — i.e. the twelve leaves
    /// <c>{open,extend}_{internal,left,right}_{insertion,deletion}_score</c> (= <c>{target,query}_{internal,left,right}_{open,extend}_gap_score</c>) and every aggregate of the
    /// property tree (<c>open_gap_score</c>, <c>end_extend_gap_score</c>, <c>left_gap_score</c>, <c>insertion_score</c>, <c>target_gap_score</c>, …; without open / extend both are set).
    /// More specific names override more general ones; equally specific names are applied in alphabetical order. An aggregate that sets only one of open / extend for a slot
    /// copies it to the other (<c>open = extend</c>, a linear gap); a slot nobody names stays <c>null</c> (the scoring's own score).
    /// </summary>
    /// <remarks>Direction note: Biopython's target is sequence 1, its query sequence 2. <c>insertion</c> / <c>target</c> = a gap in sequence 1 (<see cref="Interior2"/>, <see cref="Start2"/>, <see cref="End2"/>), <c>deletion</c> / <c>query</c> = a gap in sequence 2 (<see cref="Interior1"/>, …).</remarks>
    /// <param name="properties">Biopython property names and values.</param>
    /// <exception cref="ArgumentNullException"><paramref name="properties"/> is null.</exception>
    /// <exception cref="ArgumentException">An unknown name, or a non-finite value.</exception>
    public static GapScoreProfile FromBiopython(IReadOnlyDictionary<string, double> properties)
    {
        ArgumentNullException.ThrowIfNull(properties);
        // slot index = direction (0 target / insertion, 1 query / deletion) * 3 + position (0 internal, 1 left, 2 right)
        var open = new double?[6];
        var extend = new double?[6];
        var parsed = new List<(string Name, int Specificity, int[] Slots, bool SetOpen, bool SetExtend, double Value)>();
        foreach (var (name, value) in properties)
        {
            if (!double.IsFinite(value))
                throw new ArgumentException($"Gap score '{name}' must be finite, got {value}.", nameof(properties));
            var (specificity, slots, setOpen, setExtend) = ParseBiopythonName(name);
            parsed.Add((name, specificity, slots, setOpen, setExtend, value));
        }

        foreach (var item in parsed.OrderBy(x => x.Specificity).ThenBy(x => x.Name, StringComparer.Ordinal))
            foreach (int slot in item.Slots)
            {
                if (item.SetOpen) open[slot] = item.Value;
                if (item.SetExtend) extend[slot] = item.Value;
            }

        GapScore? Slot(int s)
        {
            if (open[s] is null && extend[s] is null) return null;
            double o = open[s] ?? extend[s]!.Value, e = extend[s] ?? open[s]!.Value;
            return GapScore.FromBiopython(o, e);
        }

        return new GapScoreProfile(
            Interior1: Slot(3), Start1: Slot(4), End1: Slot(5),
            Interior2: Slot(0), Start2: Slot(1), End2: Slot(2));
    }
}

/// <summary>
/// Result of pairwise sequence alignment.
/// </summary>
/// <param name="AlignedSequence1">Aligned (gapped, '-') segment of sequence 1.</param>
/// <param name="AlignedSequence2">Aligned (gapped, '-') segment of sequence 2.</param>
/// <param name="Score">Alignment score under the scoring model used.</param>
/// <param name="AlignmentType">Kind of alignment that produced the result.</param>
/// <param name="StartPosition1">0-based inclusive start of the aligned segment in sequence 1
/// (0 for global and for semi-global results in which sequence 1 is aligned in full; for an ends-free result
/// (<see cref="EndGaps"/>) the first index of sequence 1 consumed by the alignment path, i.e. the number of free leading
/// sequence-1 residues; −1 for a local result with no positive-scoring region).</param>
/// <param name="StartPosition2">0-based inclusive start of the aligned segment in sequence 2: 0 for global; for a
/// <b>semi-global</b> result the first index of sequence 2 that the alignment consumes, i.e. the number of
/// free leading sequence-2 residues (fitting: Biopython <c>alignment.coordinates</c>/<c>aligned</c> start of the fitted
/// interval); −1 for a local result with no positive-scoring region.</param>
/// <param name="EndPosition1">0-based inclusive end of the aligned segment in sequence 1
/// (length − 1 for global and fitting, so −1 for an empty sequence; for an ends-free result the last index of sequence 1
/// consumed by the path, the free trailing sequence-1 residues lie after it; −1 for an empty local result).</param>
/// <param name="EndPosition2">0-based inclusive end of the aligned segment in sequence 2: length − 1 for global; for a
/// <b>semi-global</b> result the last index of sequence 2 that the alignment consumes (the free trailing
/// residues lie after it; parasail <c>sg_dx</c> <c>end_ref</c>). A consumed interval that is empty
/// (a sequence aligned entirely against gaps, or empty) has <c>End = Start − 1</c>.
/// −1 for an empty local result.</param>
/// <remarks>
/// <b>One End convention (R14-1):</b> every producer (<c>GlobalAlign</c>, <c>LocalAlign</c>, the affine, linear-space, general-gap, top-N, co-optimal, SIM / matcher, EMBOSS <c>NeedleAlign</c> / <c>WaterAlign</c> /
/// <c>StretcherAlign</c> and the dual-affine aligners) reports the 0-based INCLUSIVE last consumed index; the residues of the sequence before / after the residues of a row are derived from the rows
/// (<see cref="EmbossAlignedSequence.Flanks"/>), never from an end coordinate.
/// Empty inputs: a global alignment of an empty side is the all-gap alignment (Needleman-Wunsch optimum, Biopython
/// <c>PairwiseAligner.score</c>); two empty inputs give score 0, empty rows and coordinates 0/0/−1/−1. A local
/// alignment with no positive-scoring pair is <c>("", "", 0, Local, −1, −1, −1, −1)</c>.
/// </remarks>
public sealed record AlignmentResult(
    string AlignedSequence1,
    string AlignedSequence2,
    int Score,
    AlignmentType AlignmentType,
    int StartPosition1,
    int StartPosition2,
    int EndPosition1,
    int EndPosition2)
{
    private readonly double? _exactScore;

    /// <summary>
    /// Alignment score as a real number. Equal to <see cref="Score"/> for every integer-scored aligner; the
    /// <see cref="MatrixScoring"/> aligners (substitution matrix, possibly fractional gap penalties such as
    /// EMBOSS <c>-gapextend 0.5</c>) set it to the exact score (e.g. 49.5) while <see cref="Score"/> then holds
    /// that value rounded half away from zero.
    /// </summary>
    public double ExactScore
    {
        get => _exactScore ?? Score;
        init => _exactScore = value;
    }

    /// <summary>The empty global alignment (no columns, score 0).</summary>
    public static AlignmentResult Empty => new("", "", 0, AlignmentType.Global, 0, 0, 0, 0);
}

/// <summary>
/// Statistics calculated from a pairwise alignment (EMBOSS needle/water + BLAST conventions; see
/// <c>SequenceAligner.CalculateStatistics</c>).
/// </summary>
/// <param name="Matches">Identical aligned residue pairs (EMBOSS/BLAST "identities"; Biopython <c>identities</c>).</param>
/// <param name="Mismatches">Aligned residue pairs that are not identical (Biopython <c>mismatches</c>).</param>
/// <param name="Gaps">Columns in which at least one row has a gap (EMBOSS "Gaps" count; equals the total gap length for a
/// valid pairwise alignment, Biopython <c>gaps</c>).</param>
/// <param name="AlignmentLength">Number of alignment columns, gap columns included (EMBOSS "Length"); the denominator of all percentages.</param>
/// <param name="Identity">100·Matches / AlignmentLength.</param>
/// <param name="Similarity">100·Positives / AlignmentLength (EMBOSS "Similarity", BLAST "Positives").</param>
/// <param name="GapPercent">100·Gaps / AlignmentLength.</param>
/// <param name="Positives">Aligned residue pairs (identical or not) whose substitution score is &gt; 0 (EMBOSS similarity count,
/// BLAST positives, Biopython <c>positives</c>).</param>
/// <param name="GapOpens">Number of gap runs: maximal runs of '-' in either row (BLAST tabular "gapopen", Biopython
/// <c>open_gaps</c>). A gap run in one row directly followed by a gap run in the other row counts as two.</param>
/// <param name="BlastPositives">BLAST "Positives" count: identical columns plus non-identical columns with a positive
/// substitution score (BLAST+ <c>x_FillIdentityInfo</c>: match + '+'). Differs from <paramref name="Positives"/> only for identical
/// pairs whose own score is not positive (X/X in BLOSUM62), which BLAST counts as identities.</param>
public readonly record struct AlignmentStatistics(
    int Matches,
    int Mismatches,
    int Gaps,
    int AlignmentLength,
    double Identity,
    double Similarity,
    double GapPercent,
    int Positives = 0,
    int GapOpens = 0,
    int BlastPositives = 0)
{
    /// <summary>All-zero statistics (an alignment of nothing).</summary>
    public static AlignmentStatistics Empty => new(0, 0, 0, 0, 0, 0, 0);

    /// <summary>
    /// BLAST pairwise-report percentage (NCBI BLAST+ <c>CAlignFormatUtil::GetPercentMatch</c>, used for the
    /// <c>Identities = a/L (p%)</c>, <c>Positives</c> and <c>Gaps</c> lines): an <b>integer</b>, rounded half up
    /// (<c>(int)(0.5 + 100·n/d)</c>) but never 100 unless n = d (an alignment with mismatches is shown as at most 99%).
    /// A zero-length alignment gives 0.
    /// </summary>
    public static int BlastPercent(int numerator, int denominator)
    {
        if (denominator <= 0)
            return 0;
        if (numerator == denominator)
            return 100;
        return Math.Min(99, (int)(0.5 + 100.0 * numerator / denominator));
    }

    /// <summary>BLAST "Identities" percentage, <see cref="BlastPercent"/>(<see cref="Matches"/>, <see cref="AlignmentLength"/>).</summary>
    public int BlastPercentIdentity => BlastPercent(Matches, AlignmentLength);

    /// <summary>BLAST "Positives" percentage, <see cref="BlastPercent"/>(<see cref="BlastPositives"/>, <see cref="AlignmentLength"/>).</summary>
    public int BlastPercentPositives => BlastPercent(BlastPositives, AlignmentLength);

    /// <summary>BLAST "Gaps" percentage, <see cref="BlastPercent"/>(<see cref="Gaps"/>, <see cref="AlignmentLength"/>).</summary>
    public int BlastPercentGaps => BlastPercent(Gaps, AlignmentLength);

    /// <summary>
    /// BLAST tabular (<c>-outfmt 6</c>) <c>pident</c> text: 100·Matches/AlignmentLength with three decimals (BLAST+ 2.12;
    /// "100.000" when every column matches), invariant culture.
    /// </summary>
    public string BlastTabularPidentText
    {
        get
        {
            double pident = 0.0;
            if (AlignmentLength > 0)
                pident = Matches == AlignmentLength ? 100.0 : 100.0 * Matches / AlignmentLength;
            return pident.ToString("F3", System.Globalization.CultureInfo.InvariantCulture);
        }
    }

    /// <summary>
    /// The BLAST pairwise-report statistics line: <c> Identities = 382/403 (95%), Gaps = 9/403 (2%)</c> for nucleotides, with
    /// <c>, Positives = n/L (p%)</c> after the identities for proteins (<paramref name="protein"/> = true).
    /// </summary>
    public string ToBlastString(bool protein = false)
    {
        var sb = new System.Text.StringBuilder();
        sb.Append(" Identities = ").Append(Matches).Append('/').Append(AlignmentLength).Append(" (").Append(BlastPercentIdentity).Append("%)");
        if (protein)
            sb.Append(", Positives = ").Append(BlastPositives).Append('/').Append(AlignmentLength).Append(" (").Append(BlastPercentPositives).Append("%)");
        sb.Append(", Gaps = ").Append(Gaps).Append('/').Append(AlignmentLength).Append(" (").Append(BlastPercentGaps).Append("%)");
        return sb.ToString();
    }
}

/// <summary>
/// Residue-pair substitution matrix (symbol × symbol → integer score), e.g. BLOSUM62 or the IUPAC nucleotide
/// matrix EDNAFULL. Used to decide which aligned pairs are <i>similar</i> (positive score) in
/// <c>SequenceAligner.CalculateStatistics</c> and <c>FormatAlignment</c>.
/// </summary>
/// <remarks>
/// Lookups are case-insensitive. The built-in tables are the NCBI/EMBOSS matrices
/// (<c>EBLOSUM62</c>, <c>EPAM250</c>, <c>EDNAFULL</c> from the EMBOSS 6.6 data directory; BLOSUM62 = Henikoff &amp;
/// Henikoff 1992 PNAS 89:10915, PAM250 = Dayhoff, Schwartz &amp; Orcutt 1978, EDNAFULL = Todd Lowe 1992 = NCBI NUC.4.4 +
/// U). Every cell equals Biopython 1.88 <c>substitution_matrices.load('BLOSUM62' / 'PAM250' / 'NUC.4.4')</c>
/// (EDNAFULL additionally defines U).
/// </remarks>
public sealed partial class SubstitutionMatrix
{
    private readonly int[,] _scores;
    private readonly sbyte[] _index = new sbyte[128];

    /// <summary>Matrix name (e.g. "BLOSUM62").</summary>
    public string Name { get; }

    /// <summary>Upper-case symbols in row/column order.</summary>
    public string Alphabet { get; }

    /// <summary>
    /// The power of ten the cells are stored in (R14-8): 1 for an integer matrix (BLOSUM62, EDNAFULL, ...), 10 or 100 for a matrix with fractional scores (BENNER22 ..., GONNET1992, JONES, MDM78, ...: the
    /// smallest power of ten, at most 10^6, that makes every cell an integer). The real score of a pair is <see cref="GetValue"/>; the integer accessors <see cref="TryGetScore"/> / <see cref="GetScore"/> /
    /// <see cref="ScoreAt"/> / <see cref="MaxAbsScore"/> exist only when <see cref="IsFractional"/> is false.
    /// </summary>
    public int Denominator { get; } = 1;

    /// <summary>True when a cell has a fractional score (<see cref="Denominator"/> above 1): the pairwise aligners score such a matrix exactly in integer units of 1 / scale, the integer accessors refuse it.</summary>
    public bool IsFractional => Denominator != 1;

    /// <summary>Creates a matrix from a square score table in <paramref name="alphabet"/> order.</summary>
    /// <exception cref="ArgumentException">Duplicate/non-ASCII/lower-case symbols or a non-square table.</exception>
    public SubstitutionMatrix(string name, string alphabet, int[,] scores)
    {
        ArgumentException.ThrowIfNullOrEmpty(name);
        ArgumentException.ThrowIfNullOrEmpty(alphabet);
        ArgumentNullException.ThrowIfNull(scores);
        if (alphabet.Length > sbyte.MaxValue || scores.GetLength(0) != alphabet.Length || scores.GetLength(1) != alphabet.Length)
            throw new ArgumentException("Score table must be square with one row/column per alphabet symbol.", nameof(scores));

        Array.Fill(_index, (sbyte)-1);
        for (int i = 0; i < alphabet.Length; i++)
        {
            char c = alphabet[i];
            if (c >= 128 || c != char.ToUpperInvariant(c) || _index[c] >= 0)
                throw new ArgumentException($"Alphabet symbol '{c}' is duplicated, lower-case or not ASCII.", nameof(alphabet));
            _index[c] = (sbyte)i;
        }

        Name = name;
        Alphabet = alphabet;
        _scores = (int[,])scores.Clone();
    }

    /// <summary>
    /// Creates a matrix with fractional scores (R14-8): the cells are stored exactly as integers in units of 1 / <see cref="Denominator"/>, the smallest power of ten (at most 10^6) that makes every cell integral
    /// (a table of whole numbers gives the same matrix as the integer constructor).
    /// </summary>
    /// <param name="name">Matrix name.</param>
    /// <param name="alphabet">Upper-case ASCII symbols in row / column order.</param>
    /// <param name="scores">Square table of real scores.</param>
    /// <exception cref="ArgumentException">Duplicate / lower-case / non-ASCII symbols, a non-square table, or a cell with more than 6 decimals.</exception>
    public SubstitutionMatrix(string name, string alphabet, decimal[,] scores)
        : this(name, alphabet, ScaleToIntegers(scores, out int denominator))
    {
        Denominator = denominator;      // out variable of the initializer, in scope in the body
    }

    private static int[,] ScaleToIntegers(decimal[,] scores, out int denominator)
    {
        ArgumentNullException.ThrowIfNull(scores);
        int d = 1;
        while (d <= 1_000_000 && !scores.Cast<decimal>().All(v => (v * d) == decimal.Truncate(v * d)))
            d *= 10;
        if (d > 1_000_000)
            throw new ArgumentException("A matrix score has more than 6 decimals.", nameof(scores));
        denominator = d;
        var ints = new int[scores.GetLength(0), scores.GetLength(1)];
        for (int i = 0; i < ints.GetLength(0); i++)
        {
            for (int j = 0; j < ints.GetLength(1); j++)
                ints[i, j] = checked((int)(scores[i, j] * d));
        }

        return ints;
    }

    private void RequireInteger()
    {
        if (IsFractional)
            throw new InvalidOperationException($"The matrix {Name} has fractional scores (denominator {Denominator}): use GetValue / TryGetValue, or the aligners, which score it exactly.");
    }

    /// <summary>Looks up the score of an aligned pair; false when either symbol is not in <see cref="Alphabet"/>.</summary>
    /// <exception cref="InvalidOperationException">The matrix is fractional (<see cref="IsFractional"/>): use <see cref="TryGetValue"/>.</exception>
    public bool TryGetScore(char a, char b, out int score)
    {
        RequireInteger();
        score = 0;
        a = char.ToUpperInvariant(a);
        b = char.ToUpperInvariant(b);
        if (a >= 128 || b >= 128 || _index[a] < 0 || _index[b] < 0)
            return false;
        score = _scores[_index[a], _index[b]];
        return true;
    }

    /// <summary>The real score of an aligned pair (exact for a fractional matrix, the integer for an integer one); false when either symbol is not in <see cref="Alphabet"/>.</summary>
    public bool TryGetValue(char a, char b, out decimal value)
    {
        value = 0;
        a = char.ToUpperInvariant(a);
        b = char.ToUpperInvariant(b);
        if (a >= 128 || b >= 128 || _index[a] < 0 || _index[b] < 0)
            return false;
        value = (decimal)_scores[_index[a], _index[b]] / Denominator;
        return true;
    }

    /// <summary>The real score of an aligned pair (<see cref="TryGetValue"/>).</summary>
    /// <exception cref="ArgumentException">A symbol is not in <see cref="Alphabet"/>.</exception>
    public decimal GetValue(char a, char b) =>
        TryGetValue(a, b, out decimal v) ? v : throw new ArgumentException($"Symbol pair '{a}','{b}' is not in matrix {Name} ({Alphabet}).");

    /// <summary>True when the pair scores above zero (a "similar" pair, <c>score &gt; 0</c> of the statistics and the report marks); false when it is not positive or a symbol is not in <see cref="Alphabet"/>. Works for fractional matrices.</summary>
    public bool IsPositive(char a, char b) => TryGetValue(a, b, out decimal v) && v > 0;

    /// <summary>True when <paramref name="symbol"/> (case-insensitive) is in <see cref="Alphabet"/>.</summary>
    public bool Contains(char symbol)
    {
        symbol = char.ToUpperInvariant(symbol);
        return symbol < 128 && _index[symbol] >= 0;
    }

    /// <summary>Largest absolute cell value.</summary>
    public int MaxAbsScore
    {
        get
        {
            RequireInteger();
            int max = 0;
            foreach (int v in _scores)
                max = Math.Max(max, Math.Abs(v));
            return max;
        }
    }

    /// <summary>Row/column index of an upper-case symbol, -1 when absent (alignment dynamic programs).</summary>
    public int IndexOf(char upperSymbol) => upperSymbol < 128 ? _index[upperSymbol] : -1;

    /// <summary>Raw cell by row/column index (see <see cref="IndexOf"/>).</summary>
    /// <exception cref="InvalidOperationException">The matrix is fractional: use <see cref="ScaledScoreAt"/> / <see cref="ValueAt"/>.</exception>
    public int ScoreAt(int row, int column)
    {
        RequireInteger();
        return _scores[row, column];
    }

    /// <summary>The cell in units of 1 / <see cref="Denominator"/> (the stored integer; equals <see cref="ScoreAt"/> for an integer matrix).</summary>
    public int ScaledScoreAt(int row, int column) => _scores[row, column];

    /// <summary>The real score of a cell by row / column index.</summary>
    public decimal ValueAt(int row, int column) => (decimal)_scores[row, column] / Denominator;

    /// <summary>Largest absolute real score of the matrix (also for a fractional matrix).</summary>
    public decimal MaxAbsValue
    {
        get
        {
            int max = 0;
            foreach (int v in _scores)
                max = Math.Max(max, Math.Abs(v));
            return (decimal)max / Denominator;
        }
    }

    /// <summary>Score of an aligned pair.</summary>
    /// <exception cref="ArgumentException">A symbol is not in <see cref="Alphabet"/>.</exception>
    public int GetScore(char a, char b) =>
        TryGetScore(a, b, out int s) ? s : throw new ArgumentException($"Symbol pair '{a}','{b}' is not in matrix {Name} ({Alphabet}).");

    /// <summary>BLOSUM62 (EMBOSS <c>EBLOSUM62</c> = Biopython <c>BLOSUM62</c>).</summary>
    public static SubstitutionMatrix Blosum62 => LazyBlosum62.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyBlosum62 = new(() => Build("BLOSUM62", "ARNDCQEGHILKMFPSTWYVBZX*", Blosum62Rows()));

    /// <summary>PAM250 (EMBOSS <c>EPAM250</c> = Biopython <c>PAM250</c>).</summary>
    public static SubstitutionMatrix Pam250 => LazyPam250.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyPam250 = new(() => Build("PAM250", "ARNDCQEGHILKMFPSTWYVBZX*", Pam250Rows()));

    /// <summary>IUPAC nucleotide matrix EDNAFULL (= NCBI NUC.4.4, match 5 / mismatch −4, plus U).</summary>
    public static SubstitutionMatrix EdnaFull => LazyEdnaFull.Value;
    private static readonly Lazy<SubstitutionMatrix> LazyEdnaFull = new(() => Build("EDNAFULL", "ATGCSWRYKMBVHDNU", EdnaFullRows()));


    private static SubstitutionMatrix Build(string name, string alphabet, string[] rows)
    {
        int n = alphabet.Length;
        var m = new int[n, n];
        for (int i = 0; i < n; i++)
        {
            var parts = rows[i].Split(' ');
            for (int j = 0; j < n; j++)
                m[i, j] = int.Parse(parts[j], System.Globalization.CultureInfo.InvariantCulture);
        }
        return new SubstitutionMatrix(name, alphabet, m);
    }

    private static string[] Blosum62Rows() => new[]
    {
        "4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0 -2 -1 0 -4",
        "-1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3 -1 0 -1 -4",
        "-2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3 3 0 -1 -4",
        "-2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3 4 1 -1 -4",
        "0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4",
        "-1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2 0 3 -1 -4",
        "-1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2 1 4 -1 -4",
        "0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3 -1 -2 -1 -4",
        "-2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3 0 0 -1 -4",
        "-1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3 -3 -3 -1 -4",
        "-1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1 -4 -3 -1 -4",
        "-1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2 0 1 -1 -4",
        "-1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1 -3 -1 -1 -4",
        "-2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1 -3 -3 -1 -4",
        "-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2 -2 -1 -2 -4",
        "1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2 0 0 0 -4",
        "0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0 -1 -1 0 -4",
        "-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3 -4 -3 -2 -4",
        "-2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1 -3 -2 -1 -4",
        "0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4 -3 -2 -1 -4",
        "-2 -1 3 4 -3 0 1 -1 0 -3 -4 0 -3 -3 -2 0 -1 -4 -3 -3 4 1 -1 -4",
        "-1 0 0 1 -3 3 4 -2 0 -3 -3 1 -1 -3 -1 0 -1 -3 -2 -2 1 4 -1 -4",
        "0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2 0 0 -2 -1 -1 -1 -1 -1 -4",
        "-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 1"
    };

    private static string[] Pam250Rows() => new[]
    {
        "2 -2 0 0 -2 0 0 1 -1 -1 -2 -1 -1 -3 1 1 1 -6 -3 0 0 0 0 -8",
        "-2 6 0 -1 -4 1 -1 -3 2 -2 -3 3 0 -4 0 0 -1 2 -4 -2 -1 0 -1 -8",
        "0 0 2 2 -4 1 1 0 2 -2 -3 1 -2 -3 0 1 0 -4 -2 -2 2 1 0 -8",
        "0 -1 2 4 -5 2 3 1 1 -2 -4 0 -3 -6 -1 0 0 -7 -4 -2 3 3 -1 -8",
        "-2 -4 -4 -5 12 -5 -5 -3 -3 -2 -6 -5 -5 -4 -3 0 -2 -8 0 -2 -4 -5 -3 -8",
        "0 1 1 2 -5 4 2 -1 3 -2 -2 1 -1 -5 0 -1 -1 -5 -4 -2 1 3 -1 -8",
        "0 -1 1 3 -5 2 4 0 1 -2 -3 0 -2 -5 -1 0 0 -7 -4 -2 3 3 -1 -8",
        "1 -3 0 1 -3 -1 0 5 -2 -3 -4 -2 -3 -5 0 1 0 -7 -5 -1 0 0 -1 -8",
        "-1 2 2 1 -3 3 1 -2 6 -2 -2 0 -2 -2 0 -1 -1 -3 0 -2 1 2 -1 -8",
        "-1 -2 -2 -2 -2 -2 -2 -3 -2 5 2 -2 2 1 -2 -1 0 -5 -1 4 -2 -2 -1 -8",
        "-2 -3 -3 -4 -6 -2 -3 -4 -2 2 6 -3 4 2 -3 -3 -2 -2 -1 2 -3 -3 -1 -8",
        "-1 3 1 0 -5 1 0 -2 0 -2 -3 5 0 -5 -1 0 0 -3 -4 -2 1 0 -1 -8",
        "-1 0 -2 -3 -5 -1 -2 -3 -2 2 4 0 6 0 -2 -2 -1 -4 -2 2 -2 -2 -1 -8",
        "-3 -4 -3 -6 -4 -5 -5 -5 -2 1 2 -5 0 9 -5 -3 -3 0 7 -1 -4 -5 -2 -8",
        "1 0 0 -1 -3 0 -1 0 0 -2 -3 -1 -2 -5 6 1 0 -6 -5 -1 -1 0 -1 -8",
        "1 0 1 0 0 -1 0 1 -1 -1 -3 0 -2 -3 1 2 1 -2 -3 -1 0 0 0 -8",
        "1 -1 0 0 -2 -1 0 0 -1 0 -2 0 -1 -3 0 1 3 -5 -3 0 0 -1 0 -8",
        "-6 2 -4 -7 -8 -5 -7 -7 -3 -5 -2 -3 -4 0 -6 -2 -5 17 0 -6 -5 -6 -4 -8",
        "-3 -4 -2 -4 0 -4 -4 -5 0 -1 -1 -4 -2 7 -5 -3 -3 0 10 -2 -3 -4 -2 -8",
        "0 -2 -2 -2 -2 -2 -2 -1 -2 4 2 -2 2 -1 -1 -1 0 -6 -2 4 -2 -2 -1 -8",
        "0 -1 2 3 -4 1 3 0 1 -2 -3 1 -2 -4 -1 0 0 -5 -3 -2 3 2 -1 -8",
        "0 0 1 3 -5 3 3 0 2 -2 -3 0 -2 -5 0 0 -1 -6 -4 -2 2 3 -1 -8",
        "0 -1 0 -1 -3 -1 -1 -1 -1 -1 -1 -1 -1 -2 -1 0 0 -4 -2 -1 -1 -1 -1 -8",
        "-8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 1"
    };


    private static string[] EdnaFullRows() => new[]
    {
        "5 -4 -4 -4 -4 1 1 -4 -4 1 -4 -1 -1 -1 -2 -4",
        "-4 5 -4 -4 -4 1 -4 1 1 -4 -1 -4 -1 -1 -2 5",
        "-4 -4 5 -4 1 -4 1 -4 1 -4 -1 -1 -4 -1 -2 -4",
        "-4 -4 -4 5 1 -4 -4 1 -4 1 -1 -1 -1 -4 -2 -4",
        "-4 -4 1 1 -1 -4 -2 -2 -2 -2 -1 -1 -3 -3 -1 -4",
        "1 1 -4 -4 -4 -1 -2 -2 -2 -2 -3 -3 -1 -1 -1 1",
        "1 -4 1 -4 -2 -2 -1 -4 -2 -2 -3 -1 -3 -1 -1 -4",
        "-4 1 -4 1 -2 -2 -4 -1 -2 -2 -1 -3 -1 -3 -1 1",
        "-4 1 1 -4 -2 -2 -2 -2 -1 -4 -1 -3 -3 -1 -1 1",
        "1 -4 -4 1 -2 -2 -2 -2 -4 -1 -3 -1 -1 -3 -1 -4",
        "-4 -1 -1 -1 -1 -3 -3 -1 -1 -3 -1 -2 -2 -2 -1 -1",
        "-1 -4 -1 -1 -1 -3 -1 -3 -3 -1 -2 -1 -2 -2 -1 -4",
        "-1 -1 -4 -1 -3 -1 -3 -1 -3 -1 -2 -2 -1 -2 -1 -1",
        "-1 -1 -1 -4 -3 -1 -1 -3 -1 -3 -2 -2 -2 -1 -1 -1",
        "-2 -2 -2 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2",
        "-4 5 -4 -4 -4 1 -4 1 1 -4 -1 -4 -1 -1 -2 5"
    };
}

/// <summary>
/// Result of multiple sequence alignment.
/// </summary>
public sealed record MultipleAlignmentResult(
    string[] AlignedSequences,
    string Consensus,
    int TotalScore)
{
    private readonly double? _exactTotalScore;

    /// <summary>
    /// Sum-of-pairs score as a real number. Equal to <see cref="TotalScore"/> for every integer-scored (<see cref="ScoringMatrix"/>)
    /// aligner; the <see cref="MatrixScoring"/> aligners (substitution matrix, possibly fractional gap penalties such as 0.5) set it to
    /// the exact score while <see cref="TotalScore"/> then holds that value rounded half away from zero (as
    /// <see cref="AlignmentResult.ExactScore"/> for pairwise alignments).
    /// </summary>
    public double ExactTotalScore
    {
        get => _exactTotalScore ?? TotalScore;
        init => _exactTotalScore = value;
    }

    /// <summary>The empty multiple alignment (no sequences, empty consensus, score 0).</summary>
    public static MultipleAlignmentResult Empty => new(Array.Empty<string>(), "", 0);
}

/// <summary>Gap model of the multiple aligners and of their sum-of-pairs objective (<see cref="MultipleAlignmentResult.TotalScore"/>).</summary>
public enum MsaGapModel
{
    /// <summary>
    /// Default (historical behaviour): every gap position scores <see cref="ScoringMatrix.GapExtend"/> and <see cref="ScoringMatrix.GapOpen"/>
    /// is unused; <c>TotalScore</c> is the column sum of pairs (gap–gap columns neutral).
    /// </summary>
    Linear = 0,

    /// <summary>
    /// Affine gaps (Gotoh 1982): a gap of length k costs <see cref="ScoringMatrix.GapOpen"/> + k·<see cref="ScoringMatrix.GapExtend"/>.
    /// <c>TotalScore</c> is the sum-of-pairs score with <b>natural</b> affine gap costs (Altschul 1989: every pair of rows is projected
    /// out of the alignment, gap–gap columns are skipped, and each maximal run of gap columns of the projection pays the opening once);
    /// the star methods align every sequence to the center with the exact affine pairwise aligner, the profile methods use an affine
    /// (Gotoh three-state) profile–profile dynamic program with occupancy-scaled opening/extension.
    /// </summary>
    Affine = 1,
}

/// <summary>
/// Additive options of the multiple aligners (<c>MultipleAlign</c>, <c>MultipleAlignAnchored</c>, <c>MultipleAlignProgressive</c>,
/// <c>MultipleAlignIterative</c>, <c>MultipleAlignConsistency</c>). The default instance reproduces the parameterless behaviour.
/// </summary>
public record MsaOptions
{
    /// <summary>
    /// Cancels the alignment (the same convention as <see cref="System.Threading.Tasks.ParallelOptions.CancellationToken"/>): the dynamic programs poll it at their
    /// row / iteration boundaries and the call throws <see cref="OperationCanceledException"/> without a partial result. Default: never cancelled.
    /// It does not take part in the semantics of the options (two option sets that differ only in the token align identically).
    /// </summary>
    public CancellationToken CancellationToken { get; init; }

    /// <summary>Gap model (default <see cref="MsaGapModel.Linear"/>).</summary>
    public MsaGapModel GapModel { get; init; } = MsaGapModel.Linear;

    /// <summary>
    /// <c>MultipleAlignExact</c> only: complexity guard — the largest number of dynamic-programming values (surviving lattice cells ×
    /// gap-history states) the exact sum-of-pairs program may hold (default 50 000 000 ≈ 200 MB of 4-byte scores); above it the
    /// call throws a clear <see cref="ArgumentException"/> instead of exhausting memory.
    /// </summary>
    public long ExactMaxStateCells { get; init; } = 50_000_000;

    /// <summary>
    /// <c>MultipleAlignExact</c> only: prune the k-dimensional lattice with the Carrillo–Lipman (1988) bound (default <c>true</c>;
    /// the result is identical, only the number of cells changes). <c>false</c> evaluates the whole lattice.
    /// </summary>
    public bool ExactCarrilloLipman { get; init; } = true;

    /// <summary>
    /// <c>MultipleAlignExact</c> only: a lower bound of the optimum the caller already knows (default none) - the sum-of-pairs score, under the same objective, of an
    /// alignment of the same sequences (for example the best of the guide-tree aligners of <c>Seqeron.Genomics.ClustalW</c>). The Carrillo–Lipman pruning uses the larger of it and
    /// the star alignment's score; a value above the optimum would prune the optimum away, so it must be the score of a real alignment. The optimum does not depend on it.
    /// </summary>
    public int? ExactLowerBoundHint { get; init; }

    /// <summary>
    /// Upper bound of the worker threads of the parallel pairwise stages (the center selection of the star aligners
    /// and the alignment of every sequence to the center, which run the independent pair dynamic programs in parallel):
    /// -1 (default) = the runtime's choice, 1 = sequential, n = at most n threads. 0 and values below -1 are rejected
    /// with an <see cref="ArgumentOutOfRangeException"/>. The result never depends on it (every slot is independent and written once).
    /// An exception of a pair dynamic program is rethrown with its original type, whatever the thread count.
    /// </summary>
    public int MaxDegreeOfParallelism { get; init; } = -1;
}

/// <summary>Residue alphabet of <c>SequenceAligner.MultipleAlignMuscle</c> (MUSCLE's <c>-seqtype</c>).</summary>
public enum MuscleSequenceType
{
    /// <summary>MUSCLE's guess (<c>SeqVect::GuessAlpha</c>): nucleotide when at least 95 % of the first 100 residues are A, C, G, T / U or N, else protein.</summary>
    Auto = 0,

    /// <summary>Amino acids: the 20 residues and the wildcards B, X, Z (any other letter is replaced by X).</summary>
    Protein = 1,

    /// <summary>DNA: A, C, G, T / U and the IUPAC codes (any other letter is replaced by N).</summary>
    Dna = 2,

    /// <summary>RNA: as DNA, with MUSCLE's RNA gap parameters (gap open -420, center -300).</summary>
    Rna = 3,
}

/// <summary>Profile-profile score of MUSCLE's progressive alignment (<c>-le</c>, <c>-sp</c>, <c>-sv</c>, <c>-spn</c>).</summary>
public enum MuscleScoreModel
{
    /// <summary>MUSCLE's default for the alphabet: the log-expectation score for proteins, the sum-of-pairs score <c>SPN</c> for nucleotides.</summary>
    Default = 0,

    /// <summary>Log-expectation score (Edgar 2004): <c>log(sum_ij f_i g_j p_ij / p_i p_j)</c> weighted by the column occupancies, VTML200 likelihood matrix; protein only.</summary>
    LogExpectation = 1,

    /// <summary>Sum-of-pairs score with the PAM200 matrix (<c>-sp</c>); protein only.</summary>
    SumOfPairs = 2,

    /// <summary>Sum-of-pairs score with the VTML240 matrix (<c>-sv</c>); protein only.</summary>
    SumOfPairsVtml = 3,

    /// <summary>Nucleotide sum-of-pairs score with the BLASTZ matrix (<c>-spn</c>); nucleotides only.</summary>
    SumOfPairsNucleotide = 4,
}

/// <summary>Clustering method of MUSCLE's guide trees (<c>-cluster1</c> / <c>-cluster2</c>, also <c>-cluster</c>).</summary>
public enum MuscleClusterMethod
{
    /// <summary>UPGMB (default): UPGMA with the "biased" linkage 0.1 x average + 0.9 x minimum.</summary>
    Upgmb = 0,

    /// <summary>UPGMA (average linkage).</summary>
    Upgma = 1,

    /// <summary>Single linkage (minimum).</summary>
    UpgmaMin = 2,

    /// <summary>Complete linkage (maximum).</summary>
    UpgmaMax = 3,

    /// <summary>Neighbour joining (<c>-cluster neighborjoining</c>): the unrooted tree is rooted by <see cref="MuscleRoot"/> (default: the root of the last join, "pseudo").</summary>
    NeighborJoining = 4,
}

/// <summary>Distance of MUSCLE's guide trees (<c>-distance1</c> for the unaligned sequences, <c>-distance2</c> for the draft alignment; <c>-distance</c> sets both).</summary>
public enum MuscleDistance
{
    /// <summary>MUSCLE's default: <c>kmer6_6</c> (proteins) / <c>kmer4_6</c> (nucleotides) for the first tree, <c>pctid_kimura</c> for the second.</summary>
    Default = 0,

    /// <summary>6-mers of the 6-group compressed amino-acid alphabet (A G P S T / I L M V / N D Q E B Z / R H K / F W Y / C), MAFFT's sextet5 distance.</summary>
    Kmer6_6 = 1,

    /// <summary>6-mers of the nucleotide alphabet (A, C, G, T/U, other); nucleotides only.</summary>
    Kmer4_6 = 2,

    /// <summary>3-mers of the 20-letter amino-acid alphabet, <c>fastdistkmer.cpp</c> (<c>DistKmer20_3</c>, "kmer20_3").</summary>
    Kmer20_3 = 3,

    /// <summary>4-mers of the 20-letter amino-acid alphabet, <c>fastdist.cpp</c> / <c>fastdistjones.cpp</c> (<c>FastDistKmer</c>, "kmer20_4": Jones-Taylor-Thornton style compressed 4-mers).</summary>
    Kmer20_4 = 4,

    /// <summary>3-mers of the 20-letter alphabet as bit vectors, <c>fastdistkbit.cpp</c> (<c>DistKbit20_3</c>, "kbit20_3").</summary>
    Kbit20_3 = 5,

    /// <summary>Kimura distance of the percent identity of the aligned pair (draft alignments only, "pctid_kimura").</summary>
    PctIdKimura = 6,

    /// <summary><c>-ln(max(pct identity, 0.05))</c> of the aligned pair (draft alignments only, "pctid_log").</summary>
    PctIdLog = 7,

    /// <summary>Kimura distance of the identity of a pairwise profile alignment of each pair (<c>distpwkimura.cpp</c>, "pwkimura", unaligned sequences only).</summary>
    PwKimura = 8,

    /// <summary>"pwscoredist" (<c>scoredist.cpp</c>): see the source; unaligned sequences only.</summary>
    PwScoreDist = 9,

    /// <summary>"scoredist" (<c>scoredist.cpp</c>): aligned pairs only.</summary>
    ScoreDist = 10,

    /// <summary>"edit": 1 - percent identity of the aligned pair.</summary>
    Edit = 11,
}

/// <summary>Weighting of the sequences of a profile (<c>-weight1</c> for the progressive stage, <c>-weight2</c> for the refinement; <c>-weight</c> sets both).</summary>
public enum MuscleSeqWeight
{
    /// <summary>ClustalW weights from the guide tree (default).</summary>
    ClustalW = 0,

    /// <summary>Uniform 1/n (<c>-weight none</c>).</summary>
    None = 1,

    /// <summary>Henikoff and Henikoff 1994 position-based weights (<c>-weight henikoff</c>).</summary>
    Henikoff = 2,

    /// <summary>The PSI-BLAST variant of the Henikoff weights (<c>-weight henikoffpb</c>).</summary>
    HenikoffPb = 3,

    /// <summary>Gerstein / Sonnhammer / Chothia weights (<c>-weight gsc</c>, HMMer's algorithm).</summary>
    Gsc = 4,

    /// <summary>Gotoh's three-way tree weights (<c>-weight threeway</c>).</summary>
    ThreeWay = 5,
}

/// <summary>Objective score maximised by the refinement iterations (<c>-objscore</c>).</summary>
public enum MuscleObjectiveScore
{
    /// <summary>MUSCLE's default (<c>spm</c>): XP for up to 100 sequences, SPF above.</summary>
    Spm = 0,

    /// <summary>Sum of pairs of sequences (<c>sp</c>).</summary>
    Sp = 1,

    /// <summary>Dynamic-programming score of the two profiles (<c>dp</c>).</summary>
    Dp = 2,

    /// <summary>Cross pairs, the sum of the pair scores between the two profiles (<c>xp</c>).</summary>
    Xp = 3,

    /// <summary>Sum of profile-sequence scores (<c>ps</c>, LE profile score only).</summary>
    Ps = 4,

    /// <summary>Fast approximation of the sum of pairs (<c>spf</c>).</summary>
    Spf = 5,
}

/// <summary>Treatment of terminal gaps (<c>-termgaps</c>).</summary>
public enum MuscleTermGaps
{
    /// <summary>Half the gap open penalty is charged for a terminal gap (default <c>half</c>).</summary>
    Half = 0,

    /// <summary>Terminal gaps are charged as internal gaps (<c>full</c>).</summary>
    Full = 1,

    /// <summary>Terminal gaps are charged the extension score (<c>ext</c>).</summary>
    Ext = 2,
}

/// <summary>Rooting of a neighbour-joining guide tree (<c>-root1</c> / <c>-root2</c>).</summary>
public enum MuscleRoot
{
    /// <summary>The root assigned by the clustering (default <c>pseudo</c>).</summary>
    Pseudo = 0,

    /// <summary>The unrooted tree is rooted at the middle of its longest leaf-to-leaf span (<c>midlongestspan</c>).</summary>
    MidLongestSpan = 1,

    /// <summary>The root edge minimising the average leaf distance (<c>minavgleafdist</c>).</summary>
    MinAvgLeafDist = 2,
}

/// <summary>
/// Parameters of <c>SequenceAligner.MultipleAlignMuscle</c>, the port of MUSCLE 3.8.1551 (Edgar 2004, NAR 32:1792): k-mer distance, UPGMB guide tree, ClustalW
/// tree weights, profile-profile dynamic programming with the log-expectation (protein) / sum-of-pairs (nucleotide) score, the tree-dependent refinement of iteration 2 and the
/// restricted-partition refinement of iterations 3 and above. The defaults are MUSCLE's, except <see cref="MaxIterations"/> (2 here, 8 in the program: <see cref="ProgramDefaults"/> is the option set of a plain <c>muscle -in x.fa</c>).
/// </summary>
public sealed record MuscleOptions
{
    /// <summary>
    /// The option set of a plain <c>muscle -in x.fa</c>: every MUSCLE default including <b>-maxiters 8</b> (<c>MaxIterations</c> of a new <see cref="MuscleOptions"/> is 2, the draft + tree refinement only).
    /// Verified against the installed muscle3 3.8.1551 (scripts/alignment/generate_muscle3_default_iters_oracle.py): a 40 x 300 residue family aligns differently for 8 and 16, and the program without
    /// <c>-maxiters</c> writes the rows of 8. Also the default of <c>SequenceAligner.MuscleRun</c> / <c>MuscleRefine</c> when no options are passed.
    /// </summary>
    public static MuscleOptions ProgramDefaults { get; } = new() { MaxIterations = 8 };

    /// <summary>
    /// Cancels the alignment (the same convention as <see cref="System.Threading.Tasks.ParallelOptions.CancellationToken"/>): the dynamic programs poll it at their
    /// row / iteration boundaries and the call throws <see cref="OperationCanceledException"/> without a partial result. Default: never cancelled.
    /// It does not take part in the semantics of the options (two option sets that differ only in the token align identically).
    /// </summary>
    public CancellationToken CancellationToken { get; init; }

    /// <summary>Residue alphabet (default <see cref="MuscleSequenceType.Auto"/>; the DNA / protein overloads of the method set it).</summary>
    public MuscleSequenceType SequenceType { get; init; } = MuscleSequenceType.Auto;

    /// <summary>
    /// <c>-maxiters</c>: 1 = the progressive draft only, 2 = plus the tree-dependent refinement of iteration 2 (the guide tree is rebuilt from the Kimura distances of the draft and only the
    /// changed subtrees are realigned), 3 and above = plus <c>MaxIterations - 2</c> rounds of refinement by dividing the tree into two (restricted partitions: <see cref="Anchors"/> /
    /// <see cref="ObjectiveScore"/>). Default 2 (kept for the callers of <c>MultipleAlignMuscle</c>; the program's default is <b>8</b>: <c>params.cpp</c> of muscle3 3.8.1551 has <c>g_uMaxIters = 8</c>, and the installed binary writes the rows of <c>-maxiters 8</c>, not those of 16, for a plain <c>muscle3 -in x.fa</c> although its <c>-help</c> text says "default 16" (stale). Use <see cref="ProgramDefaults"/> for the program's behaviour).
    /// </summary>
    public int MaxIterations { get; init; } = 2;

    /// <summary>Profile-profile score (<c>-le / -sp / -sv / -spn</c>; default: LE for proteins, SPN for nucleotides).</summary>
    public MuscleScoreModel ScoreModel { get; init; } = MuscleScoreModel.Default;

    /// <summary><c>-gapopen</c> (a negative score; default: -2.9 LE, -1439 SP, -300 SV, -400 SPN DNA, -420 SPN RNA).</summary>
    public double? GapOpen { get; init; }

    /// <summary><c>-gapextend</c> (a score added per extension; default 0).</summary>
    public double? GapExtend { get; init; }

    /// <summary><c>-center</c> (default -0.52 LE, 0 otherwise, -300 RNA): the offset subtracted from the profile-position score.</summary>
    public double? Center { get; init; }

    /// <summary>Distance of the first guide tree (<c>-distance1</c>; <c>-distance</c> sets this and <see cref="Distance2"/>).</summary>
    public MuscleDistance Distance1 { get; init; } = MuscleDistance.Default;

    /// <summary>Distance of the second guide tree (<c>-distance2</c>, default <see cref="MuscleDistance.PctIdKimura"/>); only the distances of aligned pairs are defined.</summary>
    public MuscleDistance Distance2 { get; init; } = MuscleDistance.Default;

    /// <summary>Clustering of the first guide tree (<c>-cluster1</c>, default UPGMB).</summary>
    public MuscleClusterMethod Cluster1 { get; init; } = MuscleClusterMethod.Upgmb;

    /// <summary>Clustering of the second guide tree of iteration 2 (<c>-cluster2</c>, default UPGMB; its distance is <see cref="Distance2"/>).</summary>
    public MuscleClusterMethod Cluster2 { get; init; } = MuscleClusterMethod.Upgmb;

    /// <summary>Rooting of the first tree (<c>-root1</c>; only changes neighbour-joining trees and the other clusterings when not <see cref="MuscleRoot.Pseudo"/>).</summary>
    public MuscleRoot Root1 { get; init; } = MuscleRoot.Pseudo;

    /// <summary>Rooting of the second tree (<c>-root2</c>).</summary>
    public MuscleRoot Root2 { get; init; } = MuscleRoot.Pseudo;

    /// <summary><c>-sueff</c> (default 0.1): the weight of the average linkage in UPGMB.</summary>
    public float Sueff { get; init; } = 0.1f;

    /// <summary>
    /// <c>-maxtrees</c> (default 1): the number of tree refinements of iteration 2 (0 switches iteration 2 off; values above 1 repeat the tree refinement only in the program's
    /// non-"low-complexity" mode, i.e. when a sequence weighting other than ClustalW is selected).
    /// </summary>
    public int MaxTrees { get; init; } = 1;

    /// <summary>Sequence weighting of the progressive stage and of iteration 2 (<c>-weight1</c>, default ClustalW).</summary>
    public MuscleSeqWeight Weight1 { get; init; } = MuscleSeqWeight.ClustalW;

    /// <summary>Sequence weighting of the refinement iterations (<c>-weight2</c>, default ClustalW).</summary>
    public MuscleSeqWeight Weight2 { get; init; } = MuscleSeqWeight.ClustalW;

    /// <summary>Objective score of the refinement iterations (<c>-objscore</c>, default <see cref="MuscleObjectiveScore.Spm"/>).</summary>
    public MuscleObjectiveScore ObjectiveScore { get; init; } = MuscleObjectiveScore.Spm;

    /// <summary>Terminal gap treatment (<c>-termgaps</c>, default half).</summary>
    public MuscleTermGaps TermGaps { get; init; } = MuscleTermGaps.Half;

    /// <summary>Use anchor columns to cut the alignment into blocks before the refinement of iterations 3+ (<c>-anchors</c>, default; <c>-noanchors</c> = false).</summary>
    public bool Anchors { get; init; } = true;

    /// <summary><c>-anchorspacing</c> (default 32): minimum spacing of the anchor columns.</summary>
    public int AnchorSpacing { get; init; } = 32;

    /// <summary><c>-smoothwindow</c> (default 7, odd): the window of the anchor score smoothing.</summary>
    public int SmoothWindow { get; init; } = 7;

    /// <summary><c>-smoothscoreceil</c> (default 3 LE, 200 SP, 90 SV, 999 SPN).</summary>
    public double? SmoothScoreCeil { get; init; }

    /// <summary><c>-minbestcolscore</c> (default 2 LE, 300 SP, 130 SV, 90 SPN): the minimum column score of an anchor column.</summary>
    public double? MinBestColScore { get; init; }

    /// <summary><c>-minsmoothscore</c> (default 1 LE, 125 SP, 40 SV, 90 SPN): the minimum smoothed score of an anchor column.</summary>
    public double? MinSmoothScore { get; init; }

    /// <summary><c>-brenner</c>: assemble the root alignment of the progressive stage from the guide tree the way <c>makerootmsab.cpp</c> does (Brenner's method).</summary>
    public bool Brenner { get; init; }

    /// <summary>
    /// <c>-diags1</c>, <c>-diags2</c>, <c>-diags</c> (diagonal optimisation) and the parameters <c>-diaglength</c> (24), <c>-diagbreak</c> (1), <c>-diagmargin</c> (5): accepted and recorded, no effect -
    /// the 3.8.1551 source compiles the diagonal dynamic program out (<c>#define VER_3_52 0</c>: <c>GlobalAlign</c> is always <c>NWSmall</c>).
    /// </summary>
    public bool Diags1 { get; init; }

    /// <inheritdoc cref="Diags1"/>
    public bool Diags2 { get; init; }

    /// <inheritdoc cref="Diags1"/>
    public int DiagLength { get; init; } = 24;

    /// <inheritdoc cref="Diags1"/>
    public int DiagBreak { get; init; } = 1;

    /// <inheritdoc cref="Diags1"/>
    public int DiagMargin { get; init; } = 5;

    /// <summary><c>-hydro</c> (5), <c>-hydrofactor</c> (1.2) and <c>-tomhydro</c>: accepted and recorded, no effect - the hydrophobicity gap modification is compiled out of 3.8.1551 (<c>#if HYDRO</c>, HYDRO is not defined).</summary>
    public int Hydro { get; init; } = 5;

    /// <inheritdoc cref="Hydro"/>
    public double HydroFactor { get; init; } = 1.2;

    /// <inheritdoc cref="Hydro"/>
    public bool TomHydro { get; init; }

    /// <summary><c>-gapopen2</c> / <c>-gapextend2</c> (double-affine scores) and <c>-gapambig</c>: accepted and recorded, no effect - <c>DOUBLE_AFFINE</c> is 0 in 3.8.1551.</summary>
    public double? GapOpen2 { get; init; }

    /// <inheritdoc cref="GapOpen2"/>
    public double? GapExtend2 { get; init; }

    /// <inheritdoc cref="GapOpen2"/>
    public double? GapAmbig { get; init; }

    /// <summary><c>-maxmb</c> and <c>-maxhours</c>: accepted, no effect on the result (memory / time budget of the program; the library does not interrupt the refinement).</summary>
    public int? MaxMb { get; init; }

    /// <inheritdoc cref="MaxMb"/>
    public double? MaxHours { get; init; }

    /// <summary>
    /// <c>-matrix</c>: the text of a user score matrix file (NCBI / BLAST style: '#' comment lines, a heading line starting with a space, one row per residue; read by <c>readmx.cpp</c>,
    /// including its quirk that the last column of a file without a trailing '*' column is not read). It replaces the score matrix of the profile score after the SP defaults are set.
    /// For DNA / RNA the program reads the file under the nucleotide alphabet (only the rows and columns A, C, G, T/U of the amino-acid matrix are read, into the 4 x 4 corner), then switches to the amino-acid
    /// alphabet and the SP score (gap open -1439, kmer6_6 distance, the methionine hack on a leading M, amino-acid wildcard counts), which this option reproduces (oracle R638 ... R645).
    /// </summary>
    public string? UserMatrix { get; init; }

    /// <summary><c>-usetree</c> / <c>-usetree_nowarn</c>: the text of a rooted binary Newick guide tree whose leaf labels are the <see cref="SequenceNames"/>; the first tree is not computed and the tree refinement of iteration 2 is skipped.</summary>
    public string? UserTree { get; init; }

    /// <summary>The FASTA labels of the sequences (default <c>s0</c>, <c>s1</c>, ...): the leaf names of <see cref="UserTree"/> and of the trees written by the program.</summary>
    public IReadOnlyList<string>? SequenceNames { get; init; }

    /// <summary><c>-stable</c> (output in input order): the 3.8.1551 program refuses it ("-stable not supported in this version of muscle"), so does the port (<see cref="NotSupportedException"/>); the library already returns the input order by default (<see cref="TreeOrder"/> = false).</summary>
    public bool Stable { get; init; }

    /// <summary><c>-clusteronly</c> (with <c>-tree1</c>): stop after the first guide tree, no alignment is computed.</summary>
    public bool ClusterOnly { get; init; }

    /// <summary><c>-computeweights</c>: stop after the progressive alignment and report the sequence weights of its rows instead of an alignment.</summary>
    public bool ComputeWeights { get; init; }

    /// <summary>
    /// <c>-refinewindow</c> (default 200), <c>-fromwindow</c> (0), <c>-towindow</c> (0 = last), <c>-windowoffset</c> (0): the windows of <c>-refinew</c> (<c>SequenceAligner.MuscleRefineWindows</c>).
    /// </summary>
    public int RefineWindow { get; init; } = 200;

    /// <inheritdoc cref="RefineWindow"/>
    public int FromWindow { get; init; }

    /// <inheritdoc cref="RefineWindow"/>
    public int ToWindow { get; init; }

    /// <inheritdoc cref="RefineWindow"/>
    public int WindowOffset { get; init; }

    /// <summary>
    /// Return the rows in the order MUSCLE writes them (the depth-first leaf order of the final guide tree, <c>-group</c>) instead of the input order (default).
    /// </summary>
    public bool TreeOrder { get; init; }
}

/// <summary>Everything one run of <c>muscle</c> writes: the alignment and the optional files (<c>-tree1</c>, <c>-tree2</c>, <c>-distmx1</c>, <c>-distmx2</c>, <c>-computeweights</c>).</summary>
public sealed record MuscleRunResult
{
    /// <summary>The alignment (rows in input order, or in output order when <see cref="MuscleOptions.TreeOrder"/> is set); empty when the run stopped early (<c>-clusteronly</c>, <c>-computeweights</c>).</summary>
    public MultipleAlignmentResult Alignment { get; init; } = MultipleAlignmentResult.Empty;

    /// <summary>The sequence names (labels) used.</summary>
    public IReadOnlyList<string> Names { get; init; } = Array.Empty<string>();

    /// <summary>The sequence ids (0-based input indexes) in the order the program writes the rows.</summary>
    public IReadOnlyList<int> OutputOrder { get; init; } = Array.Empty<int>();

    /// <summary>The first guide tree in the format of the <c>-tree1</c> file (rooted Newick with <c>%g</c> edge lengths).</summary>
    public string? Tree1 { get; init; }

    /// <summary>The guide tree after the tree refinement of iteration 2, format of <c>-tree2</c>; null when iteration 2 did not run.</summary>
    public string? Tree2 { get; init; }

    /// <summary>The distance matrix of the first tree in the format of <c>-distmx1</c>.</summary>
    public string? DistanceMatrix1 { get; init; }

    /// <summary>The distance matrix of the second tree (from the draft alignment) in the format of <c>-distmx2</c>; null when iteration 2 did not run.</summary>
    public string? DistanceMatrix2 { get; init; }

    /// <summary>The <c>-computeweights</c> file: one "name TAB weight(%.3g)" line per row of the progressive alignment.</summary>
    public string? WeightsFile { get; init; }

    /// <summary>
    /// True when the reference program's behaviour on this input is undefined: with <c>-matrix</c> on DNA / RNA the profile score loop of <c>scorepp.cpp</c> reads the never-written
    /// <c>m_uSortOrder[4]</c> of every profile column that contains all four nucleotides (a debug-built or Debian MUSCLE 3.8.1551 may crash, 8 of 60 sampled runs did). The alignment returned is the
    /// result with that memory zero (equal to the reference on every sampled run that did not crash).
    /// </summary>
    public bool ReferenceBehaviourUndefined { get; init; }
}
