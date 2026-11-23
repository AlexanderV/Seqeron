using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>The alignment mode of a <see cref="PairwiseAligner"/> (Biopython <c>aligner.mode</c>).</summary>
public enum PairwiseAlignerMode
{
    /// <summary>Global (Needleman-Wunsch, Gotoh) alignment: <c>global</c>.</summary>
    Global,
    /// <summary>Local (Smith-Waterman) alignment: <c>local</c>.</summary>
    Local,
    /// <summary>The Fast Optimal Global Sequence Alignment Algorithm (FOGSAA, Chakraborty &amp; Bandyopadhyay 2013): <c>fogsaa</c>; one optimal global alignment, found by branch and bound.</summary>
    Fogsaa,
}

/// <summary>
/// The settings of Biopython 1.88's <c>Bio.Align.PairwiseAligner</c> (R14-10): mode (<c>global</c>, <c>local</c>, <c>fogsaa</c>), match / mismatch scores, a substitution matrix, the <c>wildcard</c> letter (scores 0 against every letter
/// with scalar scoring), the twelve gap scores <c>{open,extend}_{internal,left,right}_{insertion,deletion}_score</c> with the aggregate names of Biopython (<see cref="SetGapScore"/>), gap-score functions and <c>epsilon</c>; the scoring presets
/// <c>blastn</c>, <c>megablast</c> and <c>blastp</c>. <see cref="Score"/> and <see cref="Align"/> run the library's aligners (Needleman-Wunsch / Gotoh / Smith-Waterman / Waterman-Smith-Beyer through <see cref="SequenceAligner"/>, FOGSAA through <see cref="Fogsaa"/>).
/// Target = sequence A, query = sequence B; <c>insertion</c> = a gap in the target, <c>deletion</c> = a gap in the query.
/// </summary>
public sealed partial class PairwiseAligner
{
    // slot = direction (0 insertion, 1 deletion) * 3 + position (0 internal, 1 left, 2 right); Biopython's open = the score of the first gap position
    private readonly double[] _open = [-1, -1, -1, -1, -1, -1];
    private readonly double[] _extend = [-1, -1, -1, -1, -1, -1];

    /// <summary>Creates the aligner with Biopython's defaults (global, match 1, mismatch 0, every gap score -1, epsilon 1e-6).</summary>
    public PairwiseAligner()
    {
    }

    /// <summary>Creates an aligner with a Biopython scoring preset: <c>null</c> (defaults), <c>blastn</c> (BLASTN matrix, open -7, extend -2), <c>megablast</c> (MEGABLAST matrix, -2.5 / -2.5) or <c>blastp</c> (BLASTP matrix, -12 / -1).</summary>
    /// <exception cref="ArgumentException">Unknown preset.</exception>
    public PairwiseAligner(string? scoring)
    {
        switch (scoring)
        {
            case null:
                break;
            case "blastn":
                SubstitutionMatrix = Infrastructure.SubstitutionMatrix.Get("BLASTN");
                SetGapScore("open_gap_score", -7.0);
                SetGapScore("extend_gap_score", -2.0);
                break;
            case "megablast":
                SubstitutionMatrix = Infrastructure.SubstitutionMatrix.Get("MEGABLAST");
                SetGapScore("open_gap_score", -2.5);
                SetGapScore("extend_gap_score", -2.5);
                break;
            case "blastp":
                SubstitutionMatrix = Infrastructure.SubstitutionMatrix.Get("BLASTP");
                SetGapScore("open_gap_score", -12.0);
                SetGapScore("extend_gap_score", -1.0);
                break;
            default:
                throw new ArgumentException($"Unknown scoring scheme '{scoring}'", nameof(scoring));
        }
    }

    /// <summary><c>aligner.mode</c>.</summary>
    public PairwiseAlignerMode Mode { get; set; } = PairwiseAlignerMode.Global;

    /// <summary><c>aligner.match_score</c> (ignored while a substitution matrix is set).</summary>
    public double MatchScore { get; set; } = 1.0;

    /// <summary><c>aligner.mismatch_score</c> (ignored while a substitution matrix is set).</summary>
    public double MismatchScore { get; set; }

    /// <summary><c>aligner.substitution_matrix</c>: when set, the scores of letter pairs come from it (any letter of the alphabet; case as stored).</summary>
    public SubstitutionMatrix? SubstitutionMatrix { get; set; }

    /// <summary><c>aligner.wildcard</c>: with scalar scoring this letter scores 0 against every letter (including itself); null = none. Not used with a substitution matrix.</summary>
    public char? Wildcard { get; set; }

    /// <summary><c>aligner.epsilon</c>: two scores closer than this are equal (default 1e-6).</summary>
    public double Epsilon { get; set; } = 1e-6;

    /// <summary><c>aligner.insertion_score_function</c> (<c>target_gap_function</c>): the score of a gap in the target, a function of (position in the target, gap length); when set the Waterman-Smith-Beyer algorithm is used.</summary>
    public Func<int, int, double>? InsertionScoreFunction { get; set; }

    /// <summary><c>aligner.deletion_score_function</c> (<c>query_gap_function</c>): the score of a gap in the query, a function of (position in the query, gap length).</summary>
    public Func<int, int, double>? DeletionScoreFunction { get; set; }

    private static readonly string[] SlotNames = ["internal_insertion", "left_insertion", "right_insertion", "internal_deletion", "left_deletion", "right_deletion"];

    /// <summary>
    /// Sets the gap scores a Biopython property name stands for (<c>open_gap_score</c>, <c>extend_end_deletion_score</c>, <c>target_left_gap_score</c>, <c>gap_score</c> ... every name of the property tree, current and deprecated spellings):
    /// an <c>open_</c> name sets the open scores, an <c>extend_</c> name the extend scores, a name without either both.
    /// </summary>
    /// <exception cref="ArgumentException">Unknown name or non-finite value.</exception>
    public void SetGapScore(string name, double value)
    {
        if (!double.IsFinite(value))
            throw new ArgumentException($"Gap score '{name}' must be finite, got {value}.", nameof(value));
        var (_, slots, setOpen, setExtend) = GapScoreProfile.ParseBiopythonName(name);
        foreach (int slot in slots)
        {
            if (setOpen)
                _open[slot] = value;
            if (setExtend)
                _extend[slot] = value;
            // a number replaces the gap-score function of its direction (the setters of Bio.Align)
            if (slot < 3)
                InsertionScoreFunction = null;
            else
                DeletionScoreFunction = null;
        }

        Version++;
    }

    /// <summary>
    /// Sets a gap-score function (<c>aligner.gap_score = f</c>, <c>aligner.insertion_score = f</c>, <c>aligner.deletion_score = f</c>, and the deprecated <c>target_gap_score</c> / <c>query_gap_score</c>):
    /// the score of a gap as a function of (position, length); the Waterman-Smith-Beyer algorithm is then used for that direction.
    /// </summary>
    /// <exception cref="ArgumentException">The name is not one of the three function properties.</exception>
    public void SetGapFunction(string name, Func<int, int, double> function)
    {
        ArgumentNullException.ThrowIfNull(function);
        switch (name)
        {
            case "gap_score":
                InsertionScoreFunction = function;
                DeletionScoreFunction = function;
                break;
            case "insertion_score" or "target_gap_score" or "target_gap_function":
                InsertionScoreFunction = function;
                break;
            case "deletion_score" or "query_gap_score" or "query_gap_function":
                DeletionScoreFunction = function;
                break;
            default:
                throw new ArgumentException($"'{name}' does not take a gap-score function.", nameof(name));
        }

        Version++;
    }

    /// <summary>The value of one of the twelve gap-score properties by its name (<c>open_left_insertion_score</c> ...).</summary>
    public double GetGapScore(string name)
    {
        bool open = name.StartsWith("open_", StringComparison.Ordinal);
        if (!open && !name.StartsWith("extend_", StringComparison.Ordinal))
            throw new ArgumentException($"Not a gap score leaf: '{name}'.", nameof(name));
        string body = name[(open ? 5 : 7)..];
        int slot = Array.IndexOf(SlotNames, body.Replace("_score", "", StringComparison.Ordinal));
        if (slot < 0)
            throw new ArgumentException($"Not a gap score leaf: '{name}'.", nameof(name));
        return open ? _open[slot] : _extend[slot];
    }

    internal int Version { get; private set; }

    internal double Open(int slot) => _open[slot];

    internal double Extend(int slot) => _extend[slot];

    /// <summary>True when the scoring needs the Gotoh recursion: some open score differs from its extend score.</summary>
    private bool IsAffine
    {
        get
        {
            for (int i = 0; i < 6; i++)
                if (_open[i] != _extend[i])
                    return true;
            return false;
        }
    }

    /// <summary><c>aligner.algorithm</c>: the algorithm Biopython selects (<c>Needleman-Wunsch</c>, <c>Smith-Waterman</c>, <c>Gotoh global alignment algorithm</c>, <c>Gotoh local alignment algorithm</c>, <c>Waterman-Smith-Beyer global alignment algorithm</c> / <c>local</c>, <c>Fast Optimal Global Sequence Alignment Algorithm</c>).</summary>
    public string Algorithm
    {
        get
        {
            bool local = Mode == PairwiseAlignerMode.Local;
            if (Mode == PairwiseAlignerMode.Fogsaa)
                return "Fast Optimal Global Sequence Alignment Algorithm";
            if (InsertionScoreFunction is not null || DeletionScoreFunction is not null)
                return local ? "Waterman-Smith-Beyer local alignment algorithm" : "Waterman-Smith-Beyer global alignment algorithm";
            if (IsAffine)
                return local ? "Gotoh local alignment algorithm" : "Gotoh global alignment algorithm";
            return local ? "Smith-Waterman" : "Needleman-Wunsch";
        }
    }

    /// <summary>Gap in the target (sequence A) opened in the interior: <c>open_internal_insertion_score</c>.</summary>
    public double OpenInternalInsertionScore { get => _open[0]; set => SetGapScore("open_internal_insertion_score", value); }

    /// <summary><c>extend_internal_insertion_score</c>.</summary>
    public double ExtendInternalInsertionScore { get => _extend[0]; set => SetGapScore("extend_internal_insertion_score", value); }

    /// <summary><c>open_left_insertion_score</c>.</summary>
    public double OpenLeftInsertionScore { get => _open[1]; set => SetGapScore("open_left_insertion_score", value); }

    /// <summary><c>extend_left_insertion_score</c>.</summary>
    public double ExtendLeftInsertionScore { get => _extend[1]; set => SetGapScore("extend_left_insertion_score", value); }

    /// <summary><c>open_right_insertion_score</c>.</summary>
    public double OpenRightInsertionScore { get => _open[2]; set => SetGapScore("open_right_insertion_score", value); }

    /// <summary><c>extend_right_insertion_score</c>.</summary>
    public double ExtendRightInsertionScore { get => _extend[2]; set => SetGapScore("extend_right_insertion_score", value); }

    /// <summary><c>open_internal_deletion_score</c>.</summary>
    public double OpenInternalDeletionScore { get => _open[3]; set => SetGapScore("open_internal_deletion_score", value); }

    /// <summary><c>extend_internal_deletion_score</c>.</summary>
    public double ExtendInternalDeletionScore { get => _extend[3]; set => SetGapScore("extend_internal_deletion_score", value); }

    /// <summary><c>open_left_deletion_score</c>.</summary>
    public double OpenLeftDeletionScore { get => _open[4]; set => SetGapScore("open_left_deletion_score", value); }

    /// <summary><c>extend_left_deletion_score</c>.</summary>
    public double ExtendLeftDeletionScore { get => _extend[4]; set => SetGapScore("extend_left_deletion_score", value); }

    /// <summary><c>open_right_deletion_score</c>.</summary>
    public double OpenRightDeletionScore { get => _open[5]; set => SetGapScore("open_right_deletion_score", value); }

    /// <summary><c>extend_right_deletion_score</c>.</summary>
    public double ExtendRightDeletionScore { get => _extend[5]; set => SetGapScore("extend_right_deletion_score", value); }
}
