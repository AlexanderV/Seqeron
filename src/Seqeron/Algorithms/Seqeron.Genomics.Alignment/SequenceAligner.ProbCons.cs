using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The parameters of the ProbCons 1.12 pair-HMM (<c>probcons -p file</c>, the "parameter set" of <c>Main.cc</c> <c>ReadParameters</c> / <c>PrintParameters</c>): the initial distribution of the five
/// states (a match state and two pairs of insertion states), the gap-open and gap-extension probabilities of the two pairs, and the emission tables indexed by the character code
/// (<c>emitPairs[256][256]</c> default 1e-10, <c>emitSingle[256]</c> default 1e-5; the letters of <see cref="Alphabet"/> are set for both cases).
/// </summary>
/// <remarks>
/// The tables are single precision exactly as the program holds them. <see cref="Default"/> is the program's built-in set (<c>Defaults.h</c>: emissions from BLOCKS 5.0 / BLOSUM62, transitions from
/// EM training on BAliBASE 2.0) - the same numbers T-Coffee 13.41's <c>proba_pair</c> embeds (<c>TcProbaTables</c>, verified equal) - so the pair-HMM engine of <c>TcProbaModel</c> is reused here.
/// </remarks>
public sealed class ProbConsParameters
{
    /// <summary>Number of insertion-state pairs of the Debian build (<c>-DNumInsertStates=2</c>).</summary>
    public const int InsertStatePairs = 2;

    /// <summary>Initial distribution of the five states (<c>initDistrib</c>).</summary>
    public float[] InitDistrib { get; }

    /// <summary>Gap-open probabilities (<c>gapOpen[2 x NumInsertStates]</c>).</summary>
    public float[] GapOpen { get; }

    /// <summary>Gap-extension probabilities (<c>gapExtend[2 x NumInsertStates]</c>).</summary>
    public float[] GapExtend { get; }

    /// <summary>The alphabet line of the parameter file (concatenation of its tokens).</summary>
    public string Alphabet { get; }

    /// <summary>Emission probabilities of a match, indexed by character code (256 x 256).</summary>
    public float[,] EmitPairs { get; }

    /// <summary>Emission probabilities of an insertion, indexed by character code (256).</summary>
    public float[] EmitSingle { get; }

    internal ProbConsParameters(float[] initDistrib, float[] gapOpen, float[] gapExtend, string alphabet, float[,] emitPairs, float[] emitSingle)
    {
        InitDistrib = initDistrib;
        GapOpen = gapOpen;
        GapExtend = gapExtend;
        Alphabet = alphabet;
        EmitPairs = emitPairs;
        EmitSingle = emitSingle;
    }

    private static float[,] NewPairs()
    {
        var p = new float[256, 256];
        for (int i = 0; i < 256; i++)
            for (int j = 0; j < 256; j++)
                p[i, j] = 1e-10f;
        return p;
    }

    private static float[] NewSingle()
    {
        var s = new float[256];
        Array.Fill(s, 1e-5f);
        return s;
    }

    private static void SetPair(float[,] pairs, char a, char b, float v)
    {
        int al = char.ToLowerInvariant(a) & 0xFF, au = char.ToUpperInvariant(a) & 0xFF, bl = char.ToLowerInvariant(b) & 0xFF, bu = char.ToUpperInvariant(b) & 0xFF;
        pairs[al, bl] = v;
        pairs[al, bu] = v;
        pairs[au, bl] = v;
        pairs[au, bu] = v;
        pairs[bl, al] = v;
        pairs[bl, au] = v;
        pairs[bu, al] = v;
        pairs[bu, au] = v;
    }

    internal static ProbConsParameters Blank(string alphabet) => new(new float[2 * InsertStatePairs + 1], new float[2 * InsertStatePairs], new float[2 * InsertStatePairs], alphabet, NewPairs(), NewSingle());

    /// <summary>A deep copy (the program mutates its parameter set during pre-training).</summary>
    public ProbConsParameters Clone() => new((float[])InitDistrib.Clone(), (float[])GapOpen.Clone(), (float[])GapExtend.Clone(), Alphabet, (float[,])EmitPairs.Clone(), (float[])EmitSingle.Clone());

    /// <summary>The built-in parameter set (<c>ReadParameters</c> without <c>-p</c>, <c>Defaults.h</c>).</summary>
    public static ProbConsParameters Default
    {
        get
        {
            var pairs = NewPairs();
            var single = NewSingle();
            string alphabet = TcProbaTables.ProteinAlphabet;
            int l = alphabet.Length;
            for (int i = 0; i < l; i++)
            {
                single[char.ToLowerInvariant(alphabet[i])] = TcProbaTables.EmitSingle[i];
                single[char.ToUpperInvariant(alphabet[i])] = TcProbaTables.EmitSingle[i];
                for (int j = 0; j <= i; j++)
                    SetPair(pairs, alphabet[i], alphabet[j], TcProbaTables.EmitPairs[i * l + j]);
            }

            return new ProbConsParameters((float[])TcProbaTables.InitDistrib2.Clone(), (float[])TcProbaTables.GapOpen2.Clone(), (float[])TcProbaTables.GapExtend2.Clone(), alphabet, pairs, single);
        }
    }

    /// <summary>
    /// Reads a parameter file (<c>ReadParameters</c>): three lines of numbers (initial distribution, gap open, gap extension), the alphabet line, the lower triangle of the match emissions
    /// (whitespace separated, one row per letter) and the single emissions. Missing numbers keep the value the stream extraction leaves (0 after a failed read, as in C++11).
    /// </summary>
    /// <exception cref="ArgumentException">Fewer than four lines (the program stops: "Unable to read transition parameters" / "Unable to read alphabet").</exception>
    public static ProbConsParameters Parse(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        var lines = new List<string>();
        int pos = 0;
        string? ReadLine()
        {
            if (pos >= text.Length) return null;
            int nl = text.IndexOf('\n', pos);
            string line;
            if (nl < 0)
            {
                line = text[pos..];
                pos = text.Length;
            }
            else
            {
                line = text[pos..nl];
                pos = nl + 1;
            }
            return line;
        }

        for (int i = 0; i < 3; i++)
        {
            var l = ReadLine() ?? throw new ArgumentException("Unable to read transition parameters from parameter file");
            lines.Add(l);
        }

        float[] ReadFloats(string line, int count)
        {
            var toks = line.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            var r = new float[count];
            for (int i = 0; i < count; i++)
            {
                if (i >= toks.Length || !TryParseFloat(toks[i], out r[i])) break;
            }
            return r;
        }

        var init = ReadFloats(lines[0], 2 * InsertStatePairs + 1);
        var open = ReadFloats(lines[1], 2 * InsertStatePairs);
        var ext = ReadFloats(lines[2], 2 * InsertStatePairs);
        string alphabetLine = ReadLine() ?? throw new ArgumentException("Unable to read alphabet from scoring matrix file");
        string alphabet = string.Concat(alphabetLine.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries));

        // the numbers after the alphabet line are read as one whitespace separated stream (`data >> val`)
        var rest = text[pos..].Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
        int k = 0;
        bool failed = false;
        float Next()
        {
            if (failed || k >= rest.Length || !TryParseFloat(rest[k], out float v))
            {
                failed = true;
                return 0f;
            }
            k++;
            return v;
        }

        var pairs = NewPairs();
        var single = NewSingle();
        for (int i = 0; i < alphabet.Length; i++)
            for (int j = 0; j <= i; j++)
                SetPair(pairs, alphabet[i], alphabet[j], Next());
        for (int i = 0; i < alphabet.Length; i++)
        {
            float v = Next();
            single[char.ToLowerInvariant(alphabet[i]) & 0xFF] = v;
            single[char.ToUpperInvariant(alphabet[i]) & 0xFF] = v;
        }

        return new ProbConsParameters(init, open, ext, alphabet, pairs, single);
    }

    private static bool TryParseFloat(string token, out float value)
    {
        value = 0;
        return float.TryParse(token, NumberStyles.Float, CultureInfo.InvariantCulture, out value);
    }

    /// <summary>The text <c>PrintParameters</c> writes to standard error (<c>cerr &lt;&lt; setprecision (10)</c>, i.e. <c>%.10g</c>).</summary>
    public string Describe(string message)
    {
        var sb = new StringBuilder();
        sb.Append(message).Append("\n    initDistrib[] = { ");
        for (int i = 0; i < 2 * InsertStatePairs + 1; i++) sb.Append(SequenceAligner.CFormatG10(InitDistrib[i])).Append(' ');
        sb.Append("}\n        gapOpen[] = { ");
        for (int i = 0; i < 2 * InsertStatePairs; i++) sb.Append(SequenceAligner.CFormatG10(GapOpen[i])).Append(' ');
        sb.Append("}\n      gapExtend[] = { ");
        for (int i = 0; i < 2 * InsertStatePairs; i++) sb.Append(SequenceAligner.CFormatG10(GapExtend[i])).Append(' ');
        sb.Append("}\n\n");
        return sb.ToString();
    }

    /// <summary>The parameter file the program writes after training (<c>PrintParameters</c> with a file name: <c>%.10f </c> per number).</summary>
    public string ToFileText()
    {
        var sb = new StringBuilder();
        static string F(float v)
        {
            if (!float.IsNaN(v)) return CFormat.Fixed(v, 10) + " ";
            return (float.IsNegative(v) ? "-nan" : "nan") + " ";
        }

        foreach (var v in InitDistrib.Take(2 * InsertStatePairs + 1)) sb.Append(F(v));
        sb.Append('\n');
        foreach (var v in GapOpen.Take(2 * InsertStatePairs)) sb.Append(F(v));
        sb.Append('\n');
        foreach (var v in GapExtend.Take(2 * InsertStatePairs)) sb.Append(F(v));
        sb.Append('\n');
        sb.Append(Alphabet).Append('\n');
        for (int i = 0; i < Alphabet.Length; i++)
        {
            for (int j = 0; j <= i; j++)
                sb.Append(F(EmitPairs[Alphabet[i] & 0xFF, Alphabet[j] & 0xFF]));
            sb.Append('\n');
        }

        for (int i = 0; i < Alphabet.Length; i++)
            sb.Append(F(EmitSingle[Alphabet[i] & 0xFF]));
        sb.Append('\n');
        return sb.ToString();
    }
}

/// <summary>
/// Options of ProbCons 1.12 (<c>probcons -c N -ir N -pre N -a -go X -ge Y -co C -p file</c>). The defaults are the program's: two consistency transformations, 100 iterative-refinement passes,
/// no pre-training, input order, no cosmetic gap penalties, no posterior cutoff.
/// </summary>
public sealed record ProbConsOptions
{
    /// <summary><c>-c</c> / <c>--consistency</c>: passes of the consistency transformation, 0 ... 5.</summary>
    public int ConsistencyReps { get; init; } = 2;

    /// <summary><c>-ir</c> / <c>--iterative-refinement</c>: passes of randomized-partitioning iterative refinement, 0 ... 1000. The random bipartitions come from the C library's <c>rand()</c> with its default seed 1 (glibc), as in the program.</summary>
    public int IterativeRefinementReps { get; init; } = 100;

    /// <summary><c>-pre</c> / <c>--pre-training</c>: rounds of EM pre-training of the transition probabilities on the input itself, 0 ... 20.</summary>
    public int PreTrainingReps { get; init; }

    /// <summary><c>-a</c> / <c>--alignment-order</c>: rows in alignment order instead of input order.</summary>
    public bool AlignmentOrder { get; init; }

    /// <summary><c>-go</c> / <c>--gap-open</c>: "cosmetic" gap-open penalty (&lt;= 0) of the final profile alignments (0 together with a zero extension = plain maximum expected accuracy).</summary>
    public float GapOpenPenalty { get; init; }

    /// <summary><c>-ge</c> / <c>--gap-extension</c>: "cosmetic" gap-extension penalty (&lt;= 0).</summary>
    public float GapContinuePenalty { get; init; }

    /// <summary><c>-co</c> / <c>--cutoff</c>: value in [0, 1] subtracted from every posterior probability before the profile alignments.</summary>
    public float Cutoff { get; init; }

    /// <summary><c>-p</c> / <c>--paramfile</c>: the pair-HMM parameters (null = the built-in set).</summary>
    public ProbConsParameters? Parameters { get; init; }

    /// <summary>Sequence headers (<c>&gt;header</c> lines; the first word is the name in Clustal output); null = <c>S1</c>, <c>S2</c>, ...</summary>
    public IReadOnlyList<string>? Names { get; init; }

    /// <summary><c>-v</c> / <c>--verbose</c>: the progress lines the program writes to standard error (see <see cref="ProbConsRunResult.Log"/>).</summary>
    public bool Verbose { get; init; }

    /// <summary><c>-annot</c>: also compute the per-column annotation scores of the final alignment.</summary>
    public bool Annotate { get; init; }

    /// <summary>Cancellation of the call.</summary>
    public CancellationToken CancellationToken { get; init; }
}

/// <summary>The result of <see cref="SequenceAligner.ProbConsRun"/>.</summary>
/// <param name="Alignment">Rows in output order (input order unless <see cref="ProbConsOptions.AlignmentOrder"/>), consensus and a natural BLOSUM62 sum-of-pairs score (the program reports none).</param>
/// <param name="Order">Input index of every row.</param>
/// <param name="GuideTree">The tree the program prints to standard error (<c>((a (c d)) b)</c> with the headers of the leaves).</param>
/// <param name="Annotation">The <c>-annot</c> value of every column (0 ... 200), null unless requested.</param>
/// <param name="Parameters">The parameter set that produced the alignment (after pre-training).</param>
/// <param name="Log">What the program writes to standard error after reading the input: the pre-training parameter sets, the <c>-v</c> progress lines, the tree line (<c>Alignment tree: ...</c>) and the two blank lines.</param>
public sealed record ProbConsRunResult(MultipleAlignmentResult Alignment, int[] Order, string GuideTree, int[]? Annotation, ProbConsParameters Parameters, string Log);

/// <summary>The pairwise alignment of two sequences of <see cref="SequenceAligner.ProbConsPairs"/> (<c>probcons -pairs</c> / <c>-viterbi</c>).</summary>
/// <param name="First">Input index of the first sequence.</param>
/// <param name="Second">Input index of the second sequence.</param>
/// <param name="FirstRow">The first sequence with gaps.</param>
/// <param name="SecondRow">The second sequence with gaps.</param>
public sealed record ProbConsPairAlignment(int First, int Second, string FirstRow, string SecondRow);

/// <summary>One training set of <see cref="SequenceAligner.ProbConsTrain"/>: the sequences of one input file and their headers.</summary>
/// <param name="Sequences">Residues.</param>
/// <param name="Names">Headers (null = <c>S1</c>, <c>S2</c>, ...); only the <c>-v</c> progress lines show them.</param>
public sealed record ProbConsSet(IReadOnlyList<string> Sequences, IReadOnlyList<string>? Names = null);

/// <summary>The result of <see cref="SequenceAligner.ProbConsTrain"/>: the trained parameters and, per set, the standard-error text of its pass (<c>-v</c> progress lines).</summary>
public sealed record ProbConsTrainResult(ProbConsParameters Parameters, IReadOnlyList<string> SetLogs);

/// <summary>The pairwise alignments of <see cref="SequenceAligner.ProbConsPairs"/> with the standard-error text the program writes while computing them.</summary>
public sealed record ProbConsPairsResult(IReadOnlyList<ProbConsPairAlignment> Pairs, string Log);

public static partial class SequenceAligner
{
    #region ProbCons 1.12

    // Port of ProbCons 1.12 (Do, Mahabhashyam, Brudno & Batzoglou 2005, Genome Res 15:330; Debian probcons 1.12-14, probcons_1.12.orig.tar.gz: Main.cc, ProbabilisticModel.h, SparseMatrix.h,
    // EvolutionaryTree.h, MultiSequence.h, Sequence.h, ScoreType.h, Defaults.h). The pair-HMM forward / backward recursions, the log-space `LOOKUP` / `LOG_PLUS_EQUALS` / `EXP` arithmetic and the match
    // posterior are the engine T-Coffee 13.41 embeds (`TcProbaForward`, `TcProbaBackward`, `TcLogPlusEquals`, `TcExp`): reused, not copied; the two differences of ProbCons 1.12 (the backward
    // initialisation of the last cell) are a parameter of that engine. What is ProbCons only: the character-indexed emission tables, the sparse posterior matrices (cells >= 0.01), the consistency
    // transformation (`DoRelaxation`), the "expected accuracy" distances, the UPGMA-like tree (`TreeNode::ComputeTree`), the maximum-expected-accuracy profile alignment (`BuildPosterior`,
    // `ComputeAlignment`, `ComputeAlignmentWithGapPenalties`), the randomized-partitioning iterative refinement (glibc `rand()`), the EM training (`ComputeNewParameters`), the Viterbi alignment, the
    // column annotation and the Clustal writer.

    /// <summary>
    /// ProbCons 1.12's multiple alignment of protein sequences (<c>probcons in.fa</c>): pair-HMM posterior probabilities of every pair, <c>-c</c> rounds of consistency transformation, an expected-accuracy
    /// guide tree, progressive maximum-expected-accuracy alignment and <c>-ir</c> rounds of iterative refinement. See <see cref="ProbConsRun"/> and <see cref="ProbConsOptions"/>.
    /// </summary>
    /// <param name="sequences">Residues (letters, <c>*</c>); case is kept; white space and the gap characters <c>-</c> / <c>.</c> are dropped as the program's FASTA reader does.</param>
    /// <param name="options">The program's parameters (default: the program's defaults).</param>
    /// <exception cref="ArgumentException">No sequence, an empty sequence (or one of gap characters only; the program drops such a record), a character the program rejects ("Unknown character encountered") or an option outside the program's range.</exception>
    public static MultipleAlignmentResult MultipleAlignProbCons(IEnumerable<string> sequences, ProbConsOptions? options = null)
        => ProbConsRun(sequences, options).Alignment;

    /// <summary>
    /// ProbCons 1.12 (<c>Main.cc</c> <c>DoAlign</c> / <c>ComputeFinalAlignment</c>) with the tree, the column annotation and the parameter set. The oracle is the installed Debian <c>probcons</c> 1.12-14:
    /// rows, guide tree and <c>-annot</c> values equal on every committed set (scripts/alignment/generate_probcons_oracle.py).
    /// </summary>
    /// <remarks>
    /// Faithfully reproduced behaviours (confirmed by the oracle): the emission of a character outside the 20 amino acids (B, X, Z, J, U, O, *, <c>-</c>) is the table default (1e-5 single, 1e-10
    /// pair); lower-case letters have the probabilities of the upper-case ones and keep their case; the refinement's random bipartitions use <c>rand() % 2</c> of glibc started at seed 1;
    /// <c>BuildPosterior</c>'s cutoff is subtracted over the columns <c>0 .. length-1</c> of the other sequence (not <c>1 .. length</c>), as in the source; distances are the expected accuracy
    /// divided by the shorter length, the tree merges the largest value and averages with the weight of the merged distance (<c>(d_ik + d_jk) x d_ij / 2</c>).
    /// </remarks>
    /// <exception cref="ArgumentException">See <see cref="MultipleAlignProbCons"/>.</exception>
    public static ProbConsRunResult ProbConsRun(IEnumerable<string> sequences, ProbConsOptions? options = null)
    {
        var o = options ?? new ProbConsOptions();
        var list = (sequences ?? throw new ArgumentNullException(nameof(sequences))).ToList();
        PcValidateOptions(o);
        using var cancelScope = AlignmentCancellation.Enter(o.CancellationToken);
        var engine = new ProbConsEngine(o);
        return engine.Align(list);
    }

    /// <summary>
    /// <c>probcons -pairs</c> (<paramref name="viterbi"/> = false: the maximum-expected-accuracy alignment of the (consistency transformed) posterior matrix minus <see cref="ProbConsOptions.Cutoff"/>)
    /// and <c>probcons -viterbi</c> (<paramref name="viterbi"/> = true: the Viterbi alignment of the pair-HMM; no consistency transformation): the alignment of every pair of sequences, first index &lt; second.
    /// </summary>
    public static ProbConsPairsResult ProbConsPairs(IEnumerable<string> sequences, bool viterbi, ProbConsOptions? options = null)
    {
        var o = options ?? new ProbConsOptions();
        var list = (sequences ?? throw new ArgumentNullException(nameof(sequences))).ToList();
        PcValidateOptions(o);
        using var cancelScope = AlignmentCancellation.Enter(o.CancellationToken);
        var engine = new ProbConsEngine(o);
        return engine.AllPairs(list, viterbi);
    }

    /// <summary>
    /// <c>probcons -t file [-e] a.fa b.fa ...</c>: EM training. Every set is aligned with the starting parameters (<see cref="ProbConsOptions.Parameters"/> or the defaults), the expected
    /// transition counts of every pair of every set are accumulated and averaged (over the pairs of a set, then over the sets); <paramref name="emissions"/> also re-estimates the emission tables.
    /// </summary>
    /// <returns>The trained parameter set (what the program writes to the file: <see cref="ProbConsParameters.ToFileText"/>) and the per-set progress text.</returns>
    public static ProbConsTrainResult ProbConsTrain(IEnumerable<ProbConsSet> sets, bool emissions = false, ProbConsOptions? options = null)
    {
        var o = options ?? new ProbConsOptions();
        PcValidateOptions(o);
        using var cancelScope = AlignmentCancellation.Enter(o.CancellationToken);
        var engine = new ProbConsEngine(o);
        return engine.Train(sets.ToList(), emissions);
    }

    internal static string CFormatG10(float v) => CFormat.General(v, 10);

    private static void PcValidateOptions(ProbConsOptions o)
    {
        if (o.ConsistencyReps < 0 || o.ConsistencyReps > 5)
            throw new ArgumentException("For option -c, integer must be between 0 and 5.", nameof(o));
        if (o.IterativeRefinementReps < 0 || o.IterativeRefinementReps > 1000)
            throw new ArgumentException("For option -ir, integer must be between 0 and 1000.", nameof(o));
        if (o.PreTrainingReps < 0 || o.PreTrainingReps > 20)
            throw new ArgumentException("For option -pre, integer must be between 0 and 20.", nameof(o));
        if (o.GapOpenPenalty > 0 || o.GapContinuePenalty > 0)
            throw new ArgumentException("For options -go / -ge, floating-point value must not be positive.", nameof(o));
        if (o.Cutoff < 0 || o.Cutoff > 1)
            throw new ArgumentException("For option -co, floating-point value must be between 0 and 1.", nameof(o));
    }

    /// <summary>The C library's <c>rand()</c> of glibc (TYPE_3 additive feedback generator, r[i] = r[i-3] + r[i-31], default seed 1): ProbCons never calls <c>srand</c>.</summary>
    internal sealed class GlibcRandom
    {
        private readonly int[] _state = new int[31];
        private int _f = 3, _b;

        public GlibcRandom(uint seed = 1)
        {
            // srandom_r: seed 0 means 1; the state words come from the Lehmer generator 16807 mod 2^31 - 1, then 310 outputs are discarded
            int word = seed == 0 ? 1 : unchecked((int)seed);
            _state[0] = word;
            for (int i = 1; i < 31; i++)
            {
                int hi = word / 127773, lo = word % 127773;
                word = unchecked(16807 * lo - 2836 * hi);
                if (word < 0) word += 2147483647;
                _state[i] = word;
            }

            for (int i = 0; i < 310; i++) Next();
        }

        /// <summary>The next value, 0 ... 2^31 - 1.</summary>
        public int Next()
        {
            uint val = unchecked((uint)_state[_f] + (uint)_state[_b]);
            _state[_f] = unchecked((int)val);
            int result = (int)(val >> 1);
            if (++_f >= 31) _f = 0;
            if (++_b >= 31) _b = 0;
            return result;
        }
    }

    private sealed class PcSeq
    {
        public char[] D = Array.Empty<char>(); // 1-based, D[0] = '@'
        public int Len;
        public string Header = string.Empty;
        public int SortLabel, Label;

        public static PcSeq FromResidues(string residues, string header, int label)
        {
            var d = new char[residues.Length + 1];
            d[0] = '@';
            for (int i = 0; i < residues.Length; i++) d[i + 1] = residues[i];
            return new PcSeq { D = d, Len = residues.Length, Header = header, SortLabel = label, Label = label };
        }

        public PcSeq Clone() => new() { D = (char[])D.Clone(), Len = Len, Header = Header, SortLabel = SortLabel, Label = Label };

        public PcSeq AddGaps(List<char> path, char id)
        {
            var d = new char[path.Count + 1];
            d[0] = '@';
            int src = 1;
            for (int i = 0; i < path.Count; i++)
                d[i + 1] = path[i] == 'B' || path[i] == id ? D[src++] : '-';
            return new PcSeq { D = d, Len = path.Count, Header = Header, SortLabel = SortLabel, Label = Label };
        }

        /// <summary><c>GetMapping</c>: [0] = 0, then the (1-based) column of every residue.</summary>
        public int[] Mapping()
        {
            var r = new List<int> { 0 };
            for (int i = 1; i <= Len; i++)
                if (D[i] != '-') r.Add(i);
            return r.ToArray();
        }
    }

    private sealed class PcSparse
    {
        private const float PosteriorCutoff = 0.01f;
        public int L1, L2;
        public int[] RowSize = Array.Empty<int>(), RowStart = Array.Empty<int>();
        public int[] Cols = Array.Empty<int>();
        public float[] Vals = Array.Empty<float>();

        private PcSparse() { }

        public PcSparse(int l1, int l2, float[] posterior)
        {
            L1 = l1;
            L2 = l2;
            int numCells = 0;
            for (int i = 0; i < posterior.Length; i++)
                if (posterior[i] >= PosteriorCutoff) numCells++;
            Cols = new int[numCells];
            Vals = new float[numCells];
            RowSize = new int[l1 + 1];
            RowStart = new int[l1 + 1];
            RowSize[0] = -1;
            int dp = 0;
            for (int i = 1; i <= l1; i++)
            {
                RowStart[i] = dp;
                int baseIdx = i * (l2 + 1);
                for (int j = 1; j <= l2; j++)
                {
                    float v = posterior[baseIdx + j];
                    if (v >= PosteriorCutoff)
                    {
                        Cols[dp] = j;
                        Vals[dp] = v;
                        dp++;
                    }
                }

                RowSize[i] = dp - RowStart[i];
            }
        }

        public float GetValue(int row, int col)
        {
            for (int i = 0; i < RowSize[row]; i++)
                if (Cols[RowStart[row] + i] == col) return Vals[RowStart[row] + i];
            return 0;
        }

        public PcSparse ComputeTranspose()
        {
            var ret = new PcSparse { L1 = L2, L2 = L1 };
            int numCells = Cols.Length;
            ret.Cols = new int[numCells];
            ret.Vals = new float[numCells];
            ret.RowSize = new int[L2 + 1];
            ret.RowStart = new int[L2 + 1];
            ret.RowSize[0] = -1;
            for (int i = 0; i < numCells; i++) ret.RowSize[Cols[i]]++;
            for (int i = 1; i <= L2; i++) ret.RowStart[i] = i == 1 ? 0 : ret.RowStart[i - 1] + ret.RowSize[i - 1];
            var curr = (int[])ret.RowStart.Clone();
            for (int i = 1; i <= L1; i++)
            {
                for (int j = 0; j < RowSize[i]; j++)
                {
                    int c = Cols[RowStart[i] + j];
                    ret.Cols[curr[c]] = i;
                    ret.Vals[curr[c]] = Vals[RowStart[i] + j];
                    curr[c]++;
                }
            }

            return ret;
        }

        public float[] GetPosterior()
        {
            var p = new float[(L1 + 1) * (L2 + 1)];
            for (int i = 1; i <= L1; i++)
            {
                int b = i * (L2 + 1);
                for (int j = 0; j < RowSize[i]; j++)
                    p[b + Cols[RowStart[i] + j]] = Vals[RowStart[i] + j];
            }

            return p;
        }
    }

    private sealed class PcTree
    {
        public int Label = -1;
        public PcTree? Left, Right;

        public static PcTree Compute(float[,] distMatrix)
        {
            int numSeqs = distMatrix.GetLength(0);
            var distances = (float[,])distMatrix.Clone();
            var nodes = new PcTree?[numSeqs];
            var valid = new bool[numSeqs];
            for (int i = 0; i < numSeqs; i++)
            {
                nodes[i] = new PcTree { Label = i };
                valid[i] = true;
            }

            for (int numNodesLeft = numSeqs; numNodesLeft > 1; numNodesLeft--)
            {
                float bestProb = -1;
                (int First, int Second) bestPair = (0, 0);
                for (int i = 0; i < numSeqs; i++)
                {
                    if (!valid[i]) continue;
                    for (int j = i + 1; j < numSeqs; j++)
                    {
                        if (!valid[j]) continue;
                        if (distances[i, j] > bestProb)
                        {
                            bestProb = distances[i, j];
                            bestPair = (i, j);
                        }
                    }
                }

                var parent = new PcTree { Left = nodes[bestPair.First], Right = nodes[bestPair.Second] };
                nodes[bestPair.First] = parent;
                nodes[bestPair.Second] = null;
                for (int i = 0; i < numSeqs; i++)
                {
                    if (!valid[i]) continue;
                    float v = (distances[i, bestPair.First] + distances[i, bestPair.Second]) * bestProb / 2;
                    distances[bestPair.First, i] = v;
                    distances[i, bestPair.First] = v;
                }

                valid[bestPair.Second] = false;
            }

            return nodes[0]!;
        }

        public string Print(IReadOnlyList<string> headers)
        {
            if (Label >= 0) return headers[Label];
            return "(" + Left!.Print(headers) + " " + Right!.Print(headers) + ")";
        }

        public SortedSet<int> Subtree()
        {
            if (Label != -1) return new SortedSet<int> { Label };
            var s = Left!.Subtree();
            s.UnionWith(Right!.Subtree());
            return s;
        }
    }

    /// <summary>The state of one ProbCons run: the parameter set, the pair-HMM, the random generator and the working copies of the sequences.</summary>
    private sealed class ProbConsEngine
    {
        private readonly ProbConsOptions _o;
        private readonly GlibcRandom _rand = new();
        private ProbConsParameters _p;
        private bool _alignOrder;
        private string _tree = string.Empty;
        private readonly StringBuilder _log = new();
        private const int NumInsertStates = ProbConsParameters.InsertStatePairs;
        private const int NumMatrixTypes = 1 + 2 * NumInsertStates;

        public ProbConsEngine(ProbConsOptions o)
        {
            _o = o;
            _p = (o.Parameters ?? ProbConsParameters.Default).Clone();
            _alignOrder = o.AlignmentOrder;
        }

        private void Verbose(string text)
        {
            if (_o.Verbose) _log.Append(text);
        }

        private static string Pos(List<PcSeq> seqs, int i) => (i + 1).ToString(CultureInfo.InvariantCulture) + ") " + seqs[i].Header;

        private TcProbaModel Model() => new(NumInsertStates, _p.InitDistrib, _p.GapOpen, _p.GapExtend, _p.EmitSingle, _p.EmitPairs);

        // ----- input -----

        private static List<PcSeq> Load(IReadOnlyList<string> residues, IReadOnlyList<string>? names)
        {
            RequireProgramInput(residues.Count, "ProbCons", "ERROR: No sequences read.", nameof(residues));
            var seqs = new List<PcSeq>();
            for (int i = 0; i < residues.Count; i++)
            {
                string r = residues[i] ?? throw new ArgumentException("A sequence is null.", nameof(residues));
                // The program's FASTA reader (`Sequence::Sequence(FileBuffer&, stripGaps = true)`) skips white space and drops the gap characters '-' and '.' (a '.' is read as '-'): the
                // sequences it aligns are gap-free (verified against the installed probcons: "AC-GT" is aligned as "ACGT", an all-gap record is dropped). A record that is empty after
                // that is dropped by the program; here it is an error because the result must have one row per input.
                var sb = new StringBuilder(r.Length);
                foreach (char ch in r)
                {
                    if (ch == '-' || ch == '.' || ch == ' ' || (ch >= '\t' && ch <= '\r')) continue;
                    if (!((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') || ch == '*'))
                        throw new ArgumentException($"Unknown character encountered: {ch}", nameof(residues));
                    sb.Append(ch);
                }

                if (sb.Length == 0) throw new ArgumentException($"Sequence {i + 1} has no residues (empty or only gap characters).", nameof(residues));
                r = sb.ToString();

                string header = names != null && i < names.Count ? names[i] : "S" + (i + 1).ToString(CultureInfo.InvariantCulture);
                seqs.Add(PcSeq.FromResidues(r, header, i));
            }

            return seqs;
        }

        // ----- entry points -----

        public ProbConsRunResult Align(List<string> residues)
        {
            var seqs = Load(residues, _o.Names);
            PreTrain(seqs);

            var final = DoAlign(seqs, Model(), out int[]? annotation);
            // alignment of a single sequence: the program prints the sequence itself
            var rows = final.Select(s => new string(s.D, 1, s.Len)).ToArray();
            var order = final.Select(s => s.Label).ToArray();
            var aln = PcResult(rows);
            return new ProbConsRunResult(aln, order, _tree, annotation, _p, _log.ToString());
        }

        private void PreTrain(List<PcSeq> seqs)
        {
            for (int ct = 0; ct < _o.PreTrainingReps; ct++)
            {
                var model = Model();
                _p = DoAlignTraining(seqs, model, false);
                _log.Append(_p.Describe("Recomputed parameter set:"));
            }
        }

        public ProbConsPairsResult AllPairs(List<string> residues, bool viterbi)
        {
            var seqs = Load(residues, _o.Names);
            PreTrain(seqs);
            int n = seqs.Count;
            var model = Model();
            var result = new List<ProbConsPairAlignment>();
            var sparse = new PcSparse?[n, n];
            if (!viterbi)
            {
                for (int a = 0; a < n - 1; a++)
                {
                    for (int b = a + 1; b < n; b++)
                    {
                        Verbose("Computing posterior matrix: (" + Pos(seqs, a) + " vs. (" + Pos(seqs, b) + " -- ");
                        sparse[a, b] = new PcSparse(seqs[a].Len, seqs[b].Len, Posterior(model, seqs[a], seqs[b], out _, out _));
                        Verbose("done.\n");
                    }
                }

                for (int r = 0; r < _o.ConsistencyReps; r++)
                    sparse = DoRelaxation(seqs, sparse, _o.Verbose ? _log : null);
            }

            for (int a = 0; a < n - 1; a++)
            {
                for (int b = a + 1; b < n; b++)
                {
                    Verbose("Performing pairwise alignment: (" + Pos(seqs, a) + " vs. (" + Pos(seqs, b) + " -- ");
                    List<char> path;
                    if (viterbi)
                    {
                        path = ComputeViterbiAlignment(model, seqs[a], seqs[b]).Path;
                    }
                    else
                    {
                        float[] posterior = sparse[a, b]!.GetPosterior();
                        for (int i = 0; i < posterior.Length; i++) posterior[i] -= _o.Cutoff;
                        path = ComputeAlignment(seqs[a].Len, seqs[b].Len, posterior).Path;
                    }

                    result.Add(new ProbConsPairAlignment(a, b, new string(seqs[a].AddGaps(path, 'X').D, 1, path.Count), new string(seqs[b].AddGaps(path, 'Y').D, 1, path.Count)));
                }
            }

            return new ProbConsPairsResult(result, _log.ToString());
        }

        public ProbConsTrainResult Train(List<ProbConsSet> sets, bool emissions)
        {
            var logs = new List<string>();
            var model = Model();
            var acc = _p.Clone();
            Array.Clear(acc.InitDistrib);
            Array.Clear(acc.GapOpen);
            Array.Clear(acc.GapExtend);
            if (emissions)
            {
                Array.Clear(acc.EmitPairs);
                Array.Clear(acc.EmitSingle);
            }

            foreach (var set in sets)
            {
                var seqs = Load(set.Sequences, set.Names);
                _log.Clear();
                var thisP = DoAlignTraining(seqs, model, emissions);
                logs.Add(_log.ToString());
                for (int i = 0; i < acc.InitDistrib.Length; i++) acc.InitDistrib[i] += thisP.InitDistrib[i];
                for (int i = 0; i < acc.GapOpen.Length; i++) acc.GapOpen[i] += thisP.GapOpen[i];
                for (int i = 0; i < acc.GapExtend.Length; i++) acc.GapExtend[i] += thisP.GapExtend[i];
                if (emissions)
                {
                    for (int i = 0; i < 256; i++)
                        for (int j = 0; j < 256; j++)
                            acc.EmitPairs[i, j] += thisP.EmitPairs[i, j];
                    for (int i = 0; i < 256; i++) acc.EmitSingle[i] += thisP.EmitSingle[i];
                }
            }

            int c = sets.Count;
            for (int i = 0; i < acc.InitDistrib.Length; i++) acc.InitDistrib[i] /= c;
            for (int i = 0; i < acc.GapOpen.Length; i++) acc.GapOpen[i] /= c;
            for (int i = 0; i < acc.GapExtend.Length; i++) acc.GapExtend[i] /= c;
            if (emissions)
            {
                for (int i = 0; i < 256; i++)
                    for (int j = 0; j < 256; j++)
                        acc.EmitPairs[i, j] /= c;
                for (int i = 0; i < 256; i++) acc.EmitSingle[i] /= c;
            }

            return new ProbConsTrainResult(acc, logs);
        }

        // ----- pair-HMM -----

        /// <summary>Forward / backward matrices and the match posterior of a pair (<c>ComputeForwardMatrix</c>, <c>ComputeBackwardMatrix</c>, <c>ComputePosteriorMatrix</c>).</summary>
        private static float[] Posterior(TcProbaModel m, PcSeq s1, PcSeq s2, out float[] forward, out float[] backward)
        {
            float total = ForwardBackward(m, s1, s2, out forward, out backward);
            var posterior = new float[(s1.Len + 1) * (s2.Len + 1)];
            int ij = 0;
            for (int p = 0; p < posterior.Length; p++, ij += NumMatrixTypes)
            {
                float arg = forward[ij] + backward[ij] - total;
                posterior[p] = (float)TcExp(Math.Min(0f, arg));
            }

            posterior[0] = 0;
            return posterior;
        }

        /// <summary><c>ComputeForwardMatrix</c> + <c>ComputeBackwardMatrix</c> + <c>ComputeTotalProbability</c>.</summary>
        private static float ForwardBackward(TcProbaModel m, PcSeq s1, PcSeq s2, out float[] forward, out float[] backward)
        {
            int I = s1.Len, J = s2.Len;
            RequireTcProbaFits(I, J, NumMatrixTypes, "ProbCons");
            var tIns = new float[2][,];
            tIns[0] = new float[NumMatrixTypes, I + 1];
            tIns[1] = new float[NumMatrixTypes, J + 1];
            for (int i = 1; i <= I; i++)
                for (int k = 0; k < NumMatrixTypes; k++)
                    tIns[0][k, i] = m.InsProb[s1.D[i] & 0xFF, k];
            for (int j = 1; j <= J; j++)
                for (int k = 0; k < NumMatrixTypes; k++)
                    tIns[1][k, j] = m.InsProb[s2.D[j] & 0xFF, k];
            var tMatch = new float[(I + 1) * (J + 1)];
            for (int i = 1; i <= I; i++)
                for (int j = 1; j <= J; j++)
                    tMatch[i * (J + 1) + j] = m.MatchProb[s1.D[i] & 0xFF, s2.D[j] & 0xFF];
            forward = TcProbaForward(I, J, NumMatrixTypes, NumInsertStates, m, tMatch, tIns);
            backward = TcProbaBackward(I, J, NumMatrixTypes, NumInsertStates, m, tMatch, tIns, endMatch: false);
            return TcProbaTotal(I, J, NumMatrixTypes, NumInsertStates, forward, backward);
        }

        /// <summary><c>ComputeAlignment</c>: the maximum expected accuracy alignment of a posterior matrix (three-way choice D, L, U with the program's tie order).</summary>
        private static (List<char> Path, float Total) ComputeAlignment(int l1, int l2, float[] posterior)
        {
            var oldRow = new float[l2 + 1];
            var newRow = new float[l2 + 1];
            var trace = new char[(l1 + 1) * (l2 + 1)];
            int tp = 0;
            int pp = l2 + 1;
            for (int i = 0; i <= l2; i++)
            {
                oldRow[i] = 0;
                trace[tp++] = 'L';
            }

            for (int i = 1; i <= l1; i++)
            {
                newRow[0] = 0;
                pp++;
                trace[tp++] = 'U';
                for (int j = 1; j <= l2; j++)
                {
                    float x1 = posterior[pp++] + oldRow[j - 1], x2 = newRow[j - 1], x3 = oldRow[j];
                    ChooseBestOfThree(x1, x2, x3, 'D', 'L', 'U', out newRow[j], out trace[tp++]);
                }

                (oldRow, newRow) = (newRow, oldRow);
            }

            float total = oldRow[l2];
            var alignment = new List<char>();
            int r = l1, c = l2;
            while (r != 0 || c != 0)
            {
                char ch = trace[r * (l2 + 1) + c];
                switch (ch)
                {
                    case 'L': c--; alignment.Add('Y'); break;
                    case 'U': r--; alignment.Add('X'); break;
                    case 'D': c--; r--; alignment.Add('B'); break;
                    default: throw new InvalidOperationException();
                }
            }

            alignment.Reverse();
            return (alignment, total);
        }

        private static void ChooseBestOfThree(float x1, float x2, float x3, char b1, char b2, char b3, out float x, out char b)
        {
            if (x1 >= x2)
            {
                if (x1 >= x3)
                {
                    x = x1;
                    b = b1;
                    return;
                }

                x = x3;
                b = b3;
                return;
            }

            if (x2 >= x3)
            {
                x = x2;
                b = b2;
                return;
            }

            x = x3;
            b = b3;
        }

        /// <summary><c>ComputeViterbiAlignment</c>.</summary>
        private static (List<char> Path, float Best) ComputeViterbiAlignment(TcProbaModel m, PcSeq s1, PcSeq s2)
        {
            int I = s1.Len, J = s2.Len;
            const int nmt = NumMatrixTypes;
            RequireTcProbaFits(I, J, nmt, "ProbCons");
            var viterbi = new float[nmt * (I + 1) * (J + 1)];
            Array.Fill(viterbi, TcLogZero);
            var trace = new int[nmt * (I + 1) * (J + 1)];
            Array.Fill(trace, -1);
            for (int k = 0; k < nmt; k++) viterbi[k] = m.InitialDistribution[k];
            int ij = 0, i1j = (-J - 1) * nmt, ij1 = -nmt, i1j1 = (-J - 2) * nmt;
            for (int i = 0; i <= I; i++)
            {
                int c1 = i == 0 ? '~' : s1.D[i] & 0xFF;
                for (int j = 0; j <= J; j++)
                {
                    int c2 = j == 0 ? '~' : s2.D[j] & 0xFF;
                    if (i > 0 && j > 0)
                    {
                        for (int k = 0; k < nmt; k++)
                        {
                            float newVal = viterbi[k + i1j1] + m.TransProb[k, 0] + m.MatchProb[c1, c2];
                            if (viterbi[0 + ij] < newVal)
                            {
                                viterbi[0 + ij] = newVal;
                                trace[0 + ij] = k;
                            }
                        }
                    }

                    if (i > 0)
                    {
                        for (int k = 0; k < NumInsertStates; k++)
                        {
                            float fromMatch = m.InsProb[c1, k] + viterbi[0 + i1j] + m.TransProb[0, 2 * k + 1];
                            float fromIns = m.InsProb[c1, k] + viterbi[2 * k + 1 + i1j] + m.TransProb[2 * k + 1, 2 * k + 1];
                            if (fromMatch >= fromIns)
                            {
                                viterbi[2 * k + 1 + ij] = fromMatch;
                                trace[2 * k + 1 + ij] = 0;
                            }
                            else
                            {
                                viterbi[2 * k + 1 + ij] = fromIns;
                                trace[2 * k + 1 + ij] = 2 * k + 1;
                            }
                        }
                    }

                    if (j > 0)
                    {
                        for (int k = 0; k < NumInsertStates; k++)
                        {
                            float fromMatch = m.InsProb[c2, k] + viterbi[0 + ij1] + m.TransProb[0, 2 * k + 2];
                            float fromIns = m.InsProb[c2, k] + viterbi[2 * k + 2 + ij1] + m.TransProb[2 * k + 2, 2 * k + 2];
                            if (fromMatch >= fromIns)
                            {
                                viterbi[2 * k + 2 + ij] = fromMatch;
                                trace[2 * k + 2 + ij] = 0;
                            }
                            else
                            {
                                viterbi[2 * k + 2 + ij] = fromIns;
                                trace[2 * k + 2 + ij] = 2 * k + 2;
                            }
                        }
                    }

                    ij += nmt;
                    i1j += nmt;
                    ij1 += nmt;
                    i1j1 += nmt;
                }
            }

            float bestProb = TcLogZero;
            int state = -1;
            for (int k = 0; k < nmt; k++)
            {
                float thisProb = viterbi[k + nmt * ((I + 1) * (J + 1) - 1)] + m.InitialDistribution[k];
                if (bestProb < thisProb)
                {
                    bestProb = thisProb;
                    state = k;
                }
            }

            if (state == -1) throw new InvalidOperationException("No Viterbi path.");
            var alignment = new List<char>();
            int r = I, c = J;
            while (r != 0 || c != 0)
            {
                int newState = trace[state + nmt * (r * (J + 1) + c)];
                if (state == 0) { c--; r--; alignment.Add('B'); }
                else if (state % 2 == 1) { r--; alignment.Add('X'); }
                else { c--; alignment.Add('Y'); }
                state = newState;
            }

            alignment.Reverse();
            return (alignment, bestProb);
        }

        // ----- DoAlign -----

        /// <summary>Posterior matrices of all pairs, consistency transformation, tree, final alignment (the non-training, non-pairs branch of <c>DoAlign</c>).</summary>
        private List<PcSeq> DoAlign(List<PcSeq> seqs, TcProbaModel model, out int[]? annotation)
        {
            int n = seqs.Count;
            var distances = new float[n, n];
            var sparse = new PcSparse?[n, n];
            for (int a = 0; a < n - 1; a++)
            {
                for (int b = a + 1; b < n; b++)
                {
                    Verbose("Computing posterior matrix: (" + Pos(seqs, a) + " vs. (" + Pos(seqs, b) + " -- ");
                    float[] posterior = Posterior(model, seqs[a], seqs[b], out _, out _);
                    sparse[a, b] = new PcSparse(seqs[a].Len, seqs[b].Len, posterior);
                    var al = ComputeAlignment(seqs[a].Len, seqs[b].Len, posterior);
                    float distance = al.Total / Math.Min(seqs[a].Len, seqs[b].Len);
                    distances[a, b] = distances[b, a] = distance;
                    Verbose(CFormatG10(distance) + "\n");
                }
            }

            for (int r = 0; r < _o.ConsistencyReps; r++)
                sparse = DoRelaxation(seqs, sparse, _o.Verbose ? _log : null);

            Verbose("\n");
            var tree = PcTree.Compute(distances);
            _tree = tree.Print(seqs.Select(s => s.Header).ToList());
            _log.Append("Alignment tree: ").Append(_tree).Append("\n\n");
            var alignment = ComputeFinalAlignment(tree, seqs, sparse, model);
            annotation = _o.Annotate ? WriteAnnotation(alignment, sparse) : null;
            return alignment;
        }

        /// <summary>The training branch of <c>DoAlign</c>: the parameters re-estimated from the pairs of one set (averaged over the pairs).</summary>
        private ProbConsParameters DoAlignTraining(List<PcSeq> seqs, TcProbaModel model, bool emissions)
        {
            int n = seqs.Count;
            var acc = _p.Clone();
            Array.Clear(acc.InitDistrib);
            Array.Clear(acc.GapOpen);
            Array.Clear(acc.GapExtend);
            if (emissions)
            {
                Array.Clear(acc.EmitPairs);
                Array.Clear(acc.EmitSingle);
            }

            for (int a = 0; a < n - 1; a++)
            {
                for (int b = a + 1; b < n; b++)
                {
                    Verbose("Computing posterior matrix: (" + Pos(seqs, a) + " vs. (" + Pos(seqs, b) + " -- ");
                    ForwardBackward(model, seqs[a], seqs[b], out var forward, out var backward);
                    var thisP = ComputeNewParameters(model, seqs[a], seqs[b], forward, backward, emissions);
                    for (int i = 0; i < acc.InitDistrib.Length; i++) acc.InitDistrib[i] += thisP.InitDistrib[i];
                    for (int i = 0; i < acc.GapOpen.Length; i++) acc.GapOpen[i] += thisP.GapOpen[i];
                    for (int i = 0; i < acc.GapExtend.Length; i++) acc.GapExtend[i] += thisP.GapExtend[i];
                    if (emissions)
                    {
                        for (int i = 0; i < 256; i++)
                            for (int j = 0; j < 256; j++)
                                acc.EmitPairs[i, j] += thisP.EmitPairs[i, j];
                        for (int i = 0; i < 256; i++) acc.EmitSingle[i] += thisP.EmitSingle[i];
                    }

                    Verbose("done.\n");
                }
            }

            int pairs = n * (n - 1) / 2; // integer division, as in the source (0 for one sequence: the averages become infinite / NaN)
            for (int i = 0; i < acc.InitDistrib.Length; i++) acc.InitDistrib[i] /= pairs;
            for (int i = 0; i < acc.GapOpen.Length; i++) acc.GapOpen[i] /= pairs;
            for (int i = 0; i < acc.GapExtend.Length; i++) acc.GapExtend[i] /= pairs;
            if (emissions)
            {
                for (int i = 0; i < 256; i++)
                    for (int j = 0; j < 256; j++)
                        acc.EmitPairs[i, j] /= pairs;
                for (int i = 0; i < 256; i++) acc.EmitSingle[i] /= pairs;
            }

            return acc;
        }

        /// <summary><c>ComputeNewParameters</c>: expected transition (and emission) counts of one pair -> new parameters (log-space counts, <c>exp</c> in single precision).</summary>
        private ProbConsParameters ComputeNewParameters(TcProbaModel m, PcSeq s1, PcSeq s2, float[] forward, float[] backward, bool trainEmissions)
        {
            int I = s1.Len, J = s2.Len;
            const int nmt = NumMatrixTypes;
            float totalProb = TcProbaTotal(I, J, nmt, NumInsertStates, forward, backward);
            var transCounts = new float[nmt, nmt];
            var initCounts = new float[nmt];
            var pairCounts = new float[256, 256];
            var singleCounts = new float[256];
            Array.Fill(initCounts, TcLogZero);
            Array.Fill(singleCounts, TcLogZero);
            for (int i = 0; i < nmt; i++)
                for (int j = 0; j < nmt; j++)
                    transCounts[i, j] = TcLogZero;
            for (int i = 0; i < 256; i++)
                for (int j = 0; j < 256; j++)
                    pairCounts[i, j] = TcLogZero;
            int last = nmt * ((I + 1) * (J + 1) - 1);
            int ij = 0, i1j = (-J - 1) * nmt, ij1 = -nmt, i1j1 = (-J - 2) * nmt;
            initCounts[0] = TcLogAdd(forward[0 + nmt * (1 * (J + 1) + 1)] + backward[0 + nmt * (1 * (J + 1) + 1)], forward[0 + last] + backward[0 + last]);
            for (int k = 0; k < NumInsertStates; k++)
            {
                initCounts[2 * k + 1] = TcLogAdd(forward[2 * k + 1 + nmt * (1 * (J + 1) + 0)] + backward[2 * k + 1 + nmt * (1 * (J + 1) + 0)], forward[2 * k + 1 + last] + backward[2 * k + 1 + last]);
                initCounts[2 * k + 2] = TcLogAdd(forward[2 * k + 2 + nmt * (0 * (J + 1) + 1)] + backward[2 * k + 2 + nmt * (0 * (J + 1) + 1)], forward[2 * k + 2 + last] + backward[2 * k + 2 + last]);
            }

            for (int i = 0; i <= I; i++)
            {
                int c1 = i == 0 ? '~' : char.ToUpperInvariant(s1.D[i]) & 0xFF;
                for (int j = 0; j <= J; j++)
                {
                    int c2 = j == 0 ? '~' : char.ToUpperInvariant(s2.D[j]) & 0xFF;
                    if (i > 0 && j > 0)
                    {
                        if (trainEmissions && i == 1 && j == 1)
                        {
                            TcLogPlusEquals(ref pairCounts[c1, c2], m.InitialDistribution[0] + m.MatchProb[c1, c2] + backward[0 + ij]);
                            TcLogPlusEquals(ref pairCounts[c2, c1], m.InitialDistribution[0] + m.MatchProb[c2, c1] + backward[0 + ij]);
                        }

                        for (int k = 0; k < nmt; k++)
                        {
                            TcLogPlusEquals(ref transCounts[k, 0], forward[k + i1j1] + m.TransProb[k, 0] + m.MatchProb[c1, c2] + backward[0 + ij]);
                            if ((trainEmissions && i != 1) || j != 1)
                            {
                                TcLogPlusEquals(ref pairCounts[c1, c2], forward[k + i1j1] + m.TransProb[k, 0] + m.MatchProb[c1, c2] + backward[0 + ij]);
                                TcLogPlusEquals(ref pairCounts[c2, c1], forward[k + i1j1] + m.TransProb[k, 0] + m.MatchProb[c2, c1] + backward[0 + ij]);
                            }
                        }
                    }

                    if (i > 0)
                    {
                        for (int k = 0; k < NumInsertStates; k++)
                        {
                            TcLogPlusEquals(ref transCounts[0, 2 * k + 1], forward[0 + i1j] + m.TransProb[0, 2 * k + 1] + m.InsProb[c1, k] + backward[2 * k + 1 + ij]);
                            TcLogPlusEquals(ref transCounts[2 * k + 1, 2 * k + 1], forward[2 * k + 1 + i1j] + m.TransProb[2 * k + 1, 2 * k + 1] + m.InsProb[c1, k] + backward[2 * k + 1 + ij]);
                            if (trainEmissions)
                            {
                                if (i == 1 && j == 0)
                                {
                                    TcLogPlusEquals(ref singleCounts[c1], m.InitialDistribution[2 * k + 1] + m.InsProb[c1, k] + backward[2 * k + 1 + ij]);
                                }
                                else
                                {
                                    TcLogPlusEquals(ref singleCounts[c1], forward[0 + i1j] + m.TransProb[0, 2 * k + 1] + m.InsProb[c1, k] + backward[2 * k + 1 + ij]);
                                    TcLogPlusEquals(ref singleCounts[c1], forward[2 * k + 1 + i1j] + m.TransProb[2 * k + 1, 2 * k + 1] + m.InsProb[c1, k] + backward[2 * k + 1 + ij]);
                                }
                            }
                        }
                    }

                    if (j > 0)
                    {
                        for (int k = 0; k < NumInsertStates; k++)
                        {
                            TcLogPlusEquals(ref transCounts[0, 2 * k + 2], forward[0 + ij1] + m.TransProb[0, 2 * k + 2] + m.InsProb[c2, k] + backward[2 * k + 2 + ij]);
                            TcLogPlusEquals(ref transCounts[2 * k + 2, 2 * k + 2], forward[2 * k + 2 + ij1] + m.TransProb[2 * k + 2, 2 * k + 2] + m.InsProb[c2, k] + backward[2 * k + 2 + ij]);
                            if (trainEmissions)
                            {
                                if (i == 0 && j == 1)
                                {
                                    TcLogPlusEquals(ref singleCounts[c2], m.InitialDistribution[2 * k + 2] + m.InsProb[c2, k] + backward[2 * k + 2 + ij]);
                                }
                                else
                                {
                                    TcLogPlusEquals(ref singleCounts[c2], forward[0 + ij1] + m.TransProb[0, 2 * k + 2] + m.InsProb[c2, k] + backward[2 * k + 2 + ij]);
                                    TcLogPlusEquals(ref singleCounts[c2], forward[2 * k + 2 + ij1] + m.TransProb[2 * k + 2, 2 * k + 2] + m.InsProb[c2, k] + backward[2 * k + 2 + ij]);
                                }
                            }
                        }
                    }

                    ij += nmt;
                    i1j += nmt;
                    ij1 += nmt;
                    i1j1 += nmt;
                }
            }

            for (int i = 0; i < nmt; i++)
            {
                initCounts[i] -= totalProb;
                for (int j = 0; j < nmt; j++)
                    transCounts[i, j] -= totalProb;
            }

            if (trainEmissions)
            {
                for (int i = 0; i < 256; i++)
                {
                    for (int j = 0; j < 256; j++)
                        pairCounts[i, j] -= totalProb;
                    singleCounts[i] -= totalProb;
                }
            }

            var res = ProbConsParameters.Blank(_p.Alphabet);
            float totalInitDistribCounts = 0;
            for (int i = 0; i < nmt; i++) totalInitDistribCounts += MathF.Exp(initCounts[i]);
            res.InitDistrib[0] = Math.Min(1.0f, Math.Max(0.0f, MathF.Exp(initCounts[0]) / totalInitDistribCounts));
            for (int k = 0; k < NumInsertStates; k++)
            {
                float val = (MathF.Exp(initCounts[2 * k + 1]) + MathF.Exp(initCounts[2 * k + 2])) / 2;
                res.InitDistrib[2 * k + 1] = res.InitDistrib[2 * k + 2] = Math.Min(1.0f, Math.Max(0.0f, val / totalInitDistribCounts));
            }

            float inMatchStateCounts = 0;
            for (int i = 0; i < nmt; i++) inMatchStateCounts += MathF.Exp(transCounts[0, i]);
            for (int i = 0; i < NumInsertStates; i++)
            {
                float inGapStateCounts = MathF.Exp(transCounts[2 * i + 1, 0]) + MathF.Exp(transCounts[2 * i + 1, 2 * i + 1]) + MathF.Exp(transCounts[2 * i + 2, 0]) + MathF.Exp(transCounts[2 * i + 2, 2 * i + 2]);
                res.GapOpen[2 * i] = res.GapOpen[2 * i + 1] = (MathF.Exp(transCounts[0, 2 * i + 1]) + MathF.Exp(transCounts[0, 2 * i + 2])) / (2 * inMatchStateCounts);
                res.GapExtend[2 * i] = res.GapExtend[2 * i + 1] = (MathF.Exp(transCounts[2 * i + 1, 2 * i + 1]) + MathF.Exp(transCounts[2 * i + 2, 2 * i + 2])) / inGapStateCounts;
            }

            if (trainEmissions)
            {
                float totalPairCounts = 0, totalSingleCounts = 0;
                for (int i = 0; i < 256; i++)
                {
                    for (int j = 0; j <= i; j++) totalPairCounts += MathF.Exp(pairCounts[j, i]);
                    totalSingleCounts += MathF.Exp(singleCounts[i]);
                }

                static bool IsLower(int c) => c >= 'a' && c <= 'z';
                static int Lower(int c) => c >= 'A' && c <= 'Z' ? c + 32 : c;
                for (int i = 0; i < 256; i++)
                {
                    if (IsLower(i)) continue;
                    int li = Lower(i);
                    for (int j = 0; j <= i; j++)
                    {
                        if (IsLower(j)) continue;
                        int lj = Lower(j);
                        float v = MathF.Exp(pairCounts[j, i]) / totalPairCounts;
                        res.EmitPairs[i, j] = res.EmitPairs[i, lj] = res.EmitPairs[li, j] = res.EmitPairs[li, lj] = v;
                        res.EmitPairs[j, i] = res.EmitPairs[j, li] = res.EmitPairs[lj, i] = res.EmitPairs[lj, li] = v;
                    }

                    res.EmitSingle[i] = res.EmitSingle[li] = MathF.Exp(singleCounts[i]) / totalSingleCounts;
                }
            }

            return res;
        }

        // ----- consistency transformation -----

        private static PcSparse?[,] DoRelaxation(List<PcSeq> seqs, PcSparse?[,] sparse, StringBuilder? log)
        {
            int n = seqs.Count;
            var ns = new PcSparse?[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = i + 1; j < n; j++)
                {
                    AlignmentCancellation.ThrowIfRequested();
                    int l1 = seqs[i].Len, l2 = seqs[j].Len;
                    log?.Append("Relaxing (").Append(Pos(seqs, i)).Append(" vs. (").Append(Pos(seqs, j)).Append(": ").Append(sparse[i, j]!.Cols.Length.ToString(CultureInfo.InvariantCulture)).Append(" --> ");
                    float[] posterior = sparse[i, j]!.GetPosterior();
                    for (int k = 0; k < (l1 + 1) * (l2 + 1); k++) posterior[k] += posterior[k];
                    for (int k = 0; k < n; k++)
                    {
                        if (k == i || k == j) continue;
                        if (k < i)
                        {
                            Relax1(sparse[k, i]!, sparse[k, j]!, posterior);
                        }
                        else if (k > i && k < j)
                        {
                            Relax(sparse[i, k]!, sparse[k, j]!, posterior);
                        }
                        else
                        {
                            var temp = sparse[j, k]!.ComputeTranspose();
                            Relax(sparse[i, k]!, temp, posterior);
                        }
                    }

                    for (int k = 0; k < (l1 + 1) * (l2 + 1); k++) posterior[k] /= n;

                    // mask out positions not originally in the posterior matrix
                    var xy = sparse[i, j]!;
                    for (int y = 0; y <= l2; y++) posterior[y] = 0;
                    for (int x = 1; x <= l1; x++)
                    {
                        int p = xy.RowStart[x], end = p + xy.RowSize[x];
                        int baseIdx = x * (l2 + 1);
                        int curr = 0;
                        while (p != end)
                        {
                            while (curr < xy.Cols[p])
                            {
                                posterior[baseIdx + curr] = 0;
                                curr++;
                            }

                            curr++;
                            p++;
                        }

                        while (curr <= l2)
                        {
                            posterior[baseIdx + curr] = 0;
                            curr++;
                        }
                    }

                    ns[i, j] = new PcSparse(l1, l2, posterior);
                    log?.Append(ns[i, j]!.Cols.Length.ToString(CultureInfo.InvariantCulture)).Append(" -- done.\n");
                }
            }

            return ns;
        }

        private static void Relax(PcSparse xz, PcSparse zy, float[] posterior)
        {
            int lengthX = xz.L1, lengthY = zy.L2;
            for (int i = 1; i <= lengthX; i++)
            {
                int p = xz.RowStart[i], end = p + xz.RowSize[i];
                int baseIdx = i * (lengthY + 1);
                while (p != end)
                {
                    int zp = zy.RowStart[xz.Cols[p]], zend = zp + zy.RowSize[xz.Cols[p]];
                    float xzVal = xz.Vals[p];
                    while (zp != zend)
                    {
                        posterior[baseIdx + zy.Cols[zp]] += xzVal * zy.Vals[zp];
                        zp++;
                    }

                    p++;
                }
            }
        }

        private static void Relax1(PcSparse zx, PcSparse zy, float[] posterior)
        {
            int lengthZ = zx.L1, lengthY = zy.L2;
            for (int k = 1; k <= lengthZ; k++)
            {
                int p = zx.RowStart[k], end = p + zx.RowSize[k];
                while (p != end)
                {
                    int zp = zy.RowStart[k], zend = zp + zy.RowSize[k];
                    float zxVal = zx.Vals[p];
                    int baseIdx = zx.Cols[p] * (lengthY + 1);
                    while (zp != zend)
                    {
                        posterior[baseIdx + zy.Cols[zp]] += zxVal * zy.Vals[zp];
                        zp++;
                    }

                    p++;
                }
            }
        }

        // ----- progressive alignment and refinement -----

        private List<PcSeq> ProcessTree(PcTree tree, List<PcSeq> seqs, PcSparse?[,] sparse, TcProbaModel model)
        {
            if (tree.Label == -1)
            {
                var left = ProcessTree(tree.Left!, seqs, sparse, model);
                var right = ProcessTree(tree.Right!, seqs, sparse, model);
                return AlignAlignments(left, right, sparse, model);
            }

            return new List<PcSeq> { seqs[tree.Label].Clone() };
        }

        private List<PcSeq> ComputeFinalAlignment(PcTree tree, List<PcSeq> seqs, PcSparse?[,] sparse, TcProbaModel model)
        {
            var alignment = ProcessTree(tree, seqs, sparse, model);
            var oldOrdering = new List<int>();
            if (_alignOrder)
            {
                foreach (var s in alignment) oldOrdering.Add(s.SortLabel);
                for (int i = 0; i < alignment.Count; i++) alignment[i].SortLabel = i; // SaveOrdering
                _alignOrder = false;
            }

            for (int i = 0; i < _o.IterativeRefinementReps; i++)
                alignment = DoIterativeRefinement(sparse, model, alignment);

            _log.Append('\n');
            for (int i = 0; i < oldOrdering.Count; i++) alignment[i].SortLabel = oldOrdering[i];
            return alignment;
        }

        private static List<PcSeq> Project(List<PcSeq> alignment, SortedSet<int> indices)
        {
            var idx = indices.ToArray();
            int oldLength = alignment[idx[0]].Len;
            var keep = new List<int>();
            for (int i = 1; i <= oldLength; i++)
            {
                bool found = false;
                for (int j = 0; !found && j < idx.Length; j++) found = alignment[idx[j]].D[i] != '-';
                if (found) keep.Add(i);
            }

            var res = new List<PcSeq>();
            foreach (int ix in idx)
            {
                var src = alignment[ix];
                var d = new char[keep.Count + 1];
                d[0] = '@';
                for (int c = 0; c < keep.Count; c++) d[c + 1] = src.D[keep[c]];
                res.Add(new PcSeq { D = d, Len = keep.Count, Header = src.Header, SortLabel = src.SortLabel, Label = src.Label });
            }

            return res;
        }

        private List<PcSeq> DoIterativeRefinement(PcSparse?[,] sparse, TcProbaModel model, List<PcSeq> alignment)
        {
            var groupOne = new SortedSet<int>();
            var groupTwo = new SortedSet<int>();
            for (int i = 0; i < alignment.Count; i++)
            {
                if ((_rand.Next() % 2) != 0) groupOne.Add(i);
                else groupTwo.Add(i);
            }

            if (groupOne.Count == 0 || groupTwo.Count == 0) return alignment;
            var one = Project(alignment, groupOne);
            var two = Project(alignment, groupTwo);
            return AlignAlignments(one, two, sparse, model);
        }

        /// <summary><c>BuildPosterior</c>: the posterior matrix of two profiles from the (transformed) pair matrices of their members.</summary>
        private static float[] BuildPosterior(List<PcSeq> a1, List<PcSeq> a2, PcSparse?[,] sparse, float cutoff)
        {
            int l1 = a1[0].Len, l2 = a2[0].Len;
            RequireCells(((long)l1 + 1) * ((long)l2 + 1), ProbConsProfileMaxCells, nameof(ProbConsProfileMaxCells), "ProbCons", "keeps a posterior matrix of (columns1 + 1) x (columns2 + 1) floats and the profile-alignment dynamic program",
                $"{l1:N0} x {l2:N0} alignment columns = {((long)l1 + 1) * ((long)l2 + 1):N0} cells", "The program computes the product in a 32-bit int and dies with an out-of-memory failure", "use MultipleAlignKalign or MultipleAlignMafft for sequences of this length");
            var posterior = new float[(l1 + 1) * (l2 + 1)];
            foreach (var s1 in a1)
            {
                int first = s1.Label;
                var mapping1 = s1.Mapping();
                foreach (var s2 in a2)
                {
                    int second = s2.Label;
                    var mapping2 = s2.Mapping();
                    if (first < second)
                    {
                        var matrix = sparse[first, second]!;
                        for (int ii = 1; ii <= matrix.L1; ii++)
                        {
                            int rs = matrix.RowStart[ii];
                            int baseIdx = mapping1[ii] * (l2 + 1);
                            int rowSize = matrix.RowSize[ii];
                            for (int jj = 0; jj < rowSize; jj++)
                                posterior[baseIdx + mapping2[matrix.Cols[rs + jj]]] += matrix.Vals[rs + jj];
                            for (int jj = 0; jj < matrix.L2; jj++)
                                posterior[baseIdx + mapping2[jj]] -= cutoff;
                        }
                    }
                    else
                    {
                        var matrix = sparse[second, first]!;
                        for (int jj = 1; jj <= matrix.L1; jj++)
                        {
                            int rs = matrix.RowStart[jj];
                            int baseIdx = mapping2[jj];
                            int rowSize = matrix.RowSize[jj];
                            for (int ii = 0; ii < rowSize; ii++)
                                posterior[baseIdx + mapping1[matrix.Cols[rs + ii]] * (l2 + 1)] += matrix.Vals[rs + ii];
                            for (int ii = 0; ii < matrix.L2; ii++)
                                posterior[baseIdx + mapping1[ii] * (l2 + 1)] -= cutoff;
                        }
                    }
                }
            }

            return posterior;
        }

        private List<PcSeq> AlignAlignments(List<PcSeq> a1, List<PcSeq> a2, PcSparse?[,] sparse, TcProbaModel model)
        {
            AlignmentCancellation.ThrowIfRequested();
            if (_o.Verbose)
            {
                for (int i = 0; i < a1.Count; i++) _log.Append(i == 0 ? "[" : ",").Append(a1[i].Label.ToString(CultureInfo.InvariantCulture));
                _log.Append("] vs. ");
                for (int i = 0; i < a2.Count; i++) _log.Append(i == 0 ? "[" : ",").Append(a2[i].Label.ToString(CultureInfo.InvariantCulture));
                _log.Append("]: ");
            }

            float[] posterior = BuildPosterior(a1, a2, sparse, _o.Cutoff);
            List<char> path;
            float score = 1.0f;
            if (_o.GapOpenPenalty == 0 && _o.GapContinuePenalty == 0)
                (path, score) = ComputeAlignment(a1[0].Len, a2[0].Len, posterior);
            else
                path = ComputeAlignmentWithGapPenalties(a1, a2, posterior, a1.Count, a2.Count, _o.GapOpenPenalty, _o.GapContinuePenalty);
            if (_o.Verbose)
            {
                int totLength = 0;
                foreach (var x in a1)
                    foreach (var y in a2)
                        totLength += Math.Min(x.Len, y.Len);
                _log.Append(CFormatG10(score / totLength)).Append('\n');
            }

            var result = new List<PcSeq>();
            foreach (var s in a1) result.Add(s.AddGaps(path, 'X'));
            foreach (var s in a2) result.Add(s.AddGaps(path, 'Y'));
            if (!_alignOrder) SortByLabel(result);
            return result;
        }

        /// <summary><c>MultiSequence::SortByLabel</c>: the exchange sort of the source (not stable; the labels are distinct).</summary>
        private static void SortByLabel(List<PcSeq> s)
        {
            for (int i = 0; i < s.Count - 1; i++)
                for (int j = i + 1; j < s.Count; j++)
                    if (s[i].SortLabel > s[j].SortLabel) (s[i], s[j]) = (s[j], s[i]);
        }

        /// <summary><c>ComputeAlignmentWithGapPenalties</c> (the "cosmetic" gap penalties <c>-go</c> / <c>-ge</c>).</summary>
        private static List<char> ComputeAlignmentWithGapPenalties(List<PcSeq> align1, List<PcSeq> align2, float[] posterior, int numSeqs1, int numSeqs2, float gapOpenPenalty, float gapContinuePenalty)
        {
            int l1 = align1[0].Len, l2 = align2[0].Len;
            var numActive1 = new int[l1 + 1];
            var numGapOpens1 = new int[l1 + 1];
            var numActive2 = new int[l2 + 1];
            var numGapOpens2 = new int[l2 + 1];
            foreach (var s in align1)
            {
                numActive1[0] = numGapOpens1[0] = 0;
                for (int j = 1; j <= l1; j++)
                {
                    if (s.D[j] != '-')
                    {
                        numActive1[j]++;
                        numGapOpens1[j] += (j != 1 && s.D[j - 1] != '-') ? 1 : 0;
                    }
                }
            }

            foreach (var s in align2)
            {
                numActive2[0] = numGapOpens2[0] = 0;
                for (int j = 1; j <= l2; j++)
                {
                    if (s.D[j] != '-')
                    {
                        numActive2[j]++;
                        numGapOpens2[j] += (j != 1 && s.D[j - 1] != '-') ? 1 : 0;
                    }
                }
            }

            var openingPenalty1 = new float[numSeqs1 + 1, numSeqs2 + 1];
            var continuingPenalty1 = new float[numSeqs1 + 1];
            var openingPenalty2 = new float[numSeqs2 + 1, numSeqs1 + 1];
            var continuingPenalty2 = new float[numSeqs2 + 1];
            for (int i = 0; i <= numSeqs1; i++)
                for (int j = 0; j <= numSeqs2; j++)
                    openingPenalty1[i, j] = i * (gapOpenPenalty * j + gapContinuePenalty * (numSeqs2 - j));
            for (int i = 0; i <= numSeqs1; i++)
                continuingPenalty1[i] = i * gapContinuePenalty * numSeqs2;
            for (int i = 0; i <= numSeqs2; i++)
                for (int j = 0; j <= numSeqs1; j++)
                    openingPenalty2[i, j] = i * (gapOpenPenalty * j + gapContinuePenalty * (numSeqs1 - j));
            for (int i = 0; i <= numSeqs2; i++)
                continuingPenalty2[i] = i * gapContinuePenalty * numSeqs1;

            int w = l2 + 1;
            var oldRowMatch = new float[w];
            var newRowMatch = new float[w];
            var oldRowInsertX = new float[w];
            var newRowInsertX = new float[w];
            var oldRowInsertY = new float[w];
            var newRowInsertY = new float[w];
            var trace = new char[3 * (l1 + 1) * w];
            int tp = 0;
            int pp = l2 + 1;
            for (int i = 0; i <= l2; i++)
            {
                oldRowMatch[i] = oldRowInsertX[i] = i == 0 ? 0 : TcLogZero;
                oldRowInsertY[i] = i == 0 ? 0 : oldRowInsertY[i - 1] + continuingPenalty2[numActive2[i]];
                trace[tp] = trace[tp + 1] = trace[tp + 2] = 'Y';
                tp += 3;
            }

            for (int i = 1; i <= l1; i++)
            {
                newRowMatch[0] = newRowInsertY[0] = TcLogZero;
                newRowInsertX[0] = oldRowInsertX[0] + continuingPenalty1[numActive1[i]];
                pp++;
                trace[tp] = trace[tp + 1] = trace[tp + 2] = 'X';
                tp += 3;
                for (int j = 1; j <= l2; j++)
                {
                    ChooseBestOfThree(oldRowMatch[j - 1], oldRowInsertX[j - 1], oldRowInsertY[j - 1], 'M', 'X', 'Y', out newRowMatch[j], out trace[tp++]);
                    newRowMatch[j] += posterior[pp++];
                    ChooseBestOfThree(oldRowMatch[j] + openingPenalty1[numActive1[i], numGapOpens2[j]], oldRowInsertX[j] + continuingPenalty1[numActive1[i]],
                        oldRowInsertY[j] + openingPenalty1[numActive1[i], numGapOpens2[j]], 'M', 'X', 'Y', out newRowInsertX[j], out trace[tp++]);
                    ChooseBestOfThree(newRowMatch[j - 1] + openingPenalty2[numActive2[j], numGapOpens1[i]], newRowInsertX[j - 1] + openingPenalty2[numActive2[j], numGapOpens1[i]],
                        newRowInsertY[j - 1] + continuingPenalty2[numActive2[j]], 'M', 'X', 'Y', out newRowInsertY[j], out trace[tp++]);
                }

                (oldRowMatch, newRowMatch) = (newRowMatch, oldRowMatch);
                (oldRowInsertX, newRowInsertX) = (newRowInsertX, oldRowInsertX);
                (oldRowInsertY, newRowInsertY) = (newRowInsertY, oldRowInsertY);
            }

            ChooseBestOfThree(oldRowMatch[l2], oldRowInsertX[l2], oldRowInsertY[l2], 'M', 'X', 'Y', out _, out char matrix);
            var alignment = new List<char>();
            int r = l1, c = l2;
            while (r != 0 || c != 0)
            {
                int offset = matrix == 'M' ? 0 : 2;
                if (matrix == 'X') offset = 1;
                char ch = trace[(r * w + c) * 3 + offset];
                switch (matrix)
                {
                    case 'Y': c--; alignment.Add('Y'); break;
                    case 'X': r--; alignment.Add('X'); break;
                    case 'M': c--; r--; alignment.Add('B'); break;
                    default: throw new InvalidOperationException();
                }

                matrix = ch;
            }

            alignment.Reverse();
            return alignment;
        }

        // ----- annotation -----

        private static int[] WriteAnnotation(List<PcSeq> alignment, PcSparse?[,] sparse)
        {
            int alignLength = alignment[0].Len;
            int n = alignment.Count;
            var position = new int[n];
            var lab = alignment.Select(s => s.SortLabel).ToArray();
            var scores = new int[alignLength];
            var active = new List<(int Lab, int Pos)>();
            for (int i = 1; i <= alignLength; i++)
            {
                active.Clear();
                for (int j = 0; j < n; j++)
                    if (alignment[j].D[i] != '-') active.Add((lab[j], ++position[j]));
                active.Sort();
                scores[i - 1] = ComputeScore(active, sparse);
            }

            return scores;
        }

        private static int ComputeScore(List<(int Lab, int Pos)> active, PcSparse?[,] sparse)
        {
            if (active.Count <= 1) return 0;
            float val = 0;
            for (int i = 0; i < active.Count; i++)
                for (int j = i + 1; j < active.Count; j++)
                    val += sparse[active[i].Lab, active[j].Lab]!.GetValue(active[i].Pos, active[j].Pos);
            return (int)(200 * val / (active.Count * (active.Count - 1)));
        }

        // ----- result -----

        private static MultipleAlignmentResult PcResult(string[] rows)
        {
            string consensus = BuildConsensus(rows.ToList(), rows[0].Length, anySymbol: true);
            string[] scored = rows.Select(r => r.ToUpperInvariant()).ToArray();
            try
            {
                var sp = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.2);
                double exact = MatrixSumOfPairsScore(scored, sp, MsaGapModel.Affine);
                return new MultipleAlignmentResult(rows, consensus, (int)Math.Round(exact, MidpointRounding.AwayFromZero)) { ExactTotalScore = exact };
            }
            catch (ArgumentException)
            {
                // letters the BLOSUM62 scoring does not know (B, X, Z, J, *, ...): the program reports no score, so none is reported
                return new MultipleAlignmentResult(rows, consensus, 0) { ExactTotalScore = double.NaN };
            }
        }
    }

    #endregion
}
