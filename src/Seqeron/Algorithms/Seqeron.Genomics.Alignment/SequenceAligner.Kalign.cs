using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

/// <summary>The alignment type of Kalign 3.3.5 (<c>--type</c>; <c>kalign.h</c> KALIGN_TYPE_*): selects the substitution matrix and the default gap penalties.</summary>
public enum KalignType
{
    /// <summary><c>--type dna</c>: match 5 / mismatch -4, gap open 8, extension 6, terminal extension 0 (nucleotides only).</summary>
    Dna = 0,
    /// <summary>KALIGN_TYPE_DNA_INTERNAL: as <see cref="Dna"/> with terminal gap extension 8. Library only: the command line cannot select it (<c>--type internal</c> matches <c>strstr(in, "rna")</c> first and means RNA).</summary>
    DnaInternal = 1,
    /// <summary><c>--type rna</c> and the default for nucleotides: the RIBOSUM-like 283 + {91, -114, -31, -123, ...} matrix, gap open 217, extension 39.4, terminal extension 292.6.</summary>
    Rna = 2,
    /// <summary><c>--type protein</c> and the default for proteins: CorBLOSUM66_13plus, gap open 5.5, extension 2, terminal extension 1.</summary>
    Protein = 3,
    /// <summary><c>--type divergent</c>: the GON250 matrix, gap open 55, extension 8, terminal extension 4 (proteins only).</summary>
    ProteinDivergent = 4,
    /// <summary>No <c>--type</c> option (KALIGN_TYPE_UNDEFINED): RNA matrix for nucleotides, CorBLOSUM66_13plus for proteins.</summary>
    Undefined = 5,
}

/// <summary>
/// Options of Kalign 3.3.5 (<c>kalign -i in.fa --type T --gpo X --gpe Y --tgpe Z -n N</c>, <c>kalign_run</c> of <c>lib/src/aln_wrap.c</c>).
/// Gap penalties below 0 mean "the default of the alignment type" (the program's own sentinel; <c>aln_param_init</c> applies a value only when it is &gt;= 0,
/// and, as a quirk of 3.3.5, applies <see cref="TerminalGapExtend"/> only when <see cref="GapExtend"/> is &gt;= 0 - so <c>--tgpe</c> alone has no effect and
/// <c>--gpe</c> alone sets the terminal extension to the unset -1).
/// </summary>
public sealed record KalignOptions
{
    /// <summary><c>--type</c>.</summary>
    public KalignType Type { get; init; } = KalignType.Undefined;

    /// <summary><c>--gpo</c>: gap open penalty (positive = cost), -1 = default of the type.</summary>
    public float GapOpen { get; init; } = -1f;

    /// <summary><c>--gpe</c>: gap extension penalty, -1 = default of the type.</summary>
    public float GapExtend { get; init; } = -1f;

    /// <summary><c>--tgpe</c>: terminal gap extension penalty, -1 = default of the type (see the remarks of <see cref="KalignOptions"/>).</summary>
    public float TerminalGapExtend { get; init; } = -1f;

    /// <summary>
    /// <c>-n</c> / <c>--nthreads</c> (default 4, at least 1). The program's results depend on whether it runs with one thread or more: with one thread the serial
    /// recursion (<c>recursive_aln_serial</c> / <c>do_align_serial</c>: sequence-sequence alignments are never swapped), with more the OpenMP recursion
    /// (<c>do_align</c>: the shorter sequence is always the first). The port is single threaded in both cases; only the thread count's two variants are reproduced.
    /// </summary>
    public int Threads { get; init; } = 4;

    /// <summary>Cancellation of the call.</summary>
    public CancellationToken CancellationToken { get; init; }
}

/// <summary>Result of <see cref="SequenceAligner.KalignRun"/>: the alignment and what the program decided on the way.</summary>
/// <param name="Alignment">Rows (input order, original letters, '-' for gaps) of the kept sequences, consensus and score.</param>
/// <param name="KeptIndices">Indices (into the input list) of the sequences that were aligned; sequences without residues are dropped by the program.</param>
/// <param name="Alphabet"><c>dna</c> or <c>protein</c>: the alphabet detected from the letter frequencies of the input (<c>detect_alphabet</c>).</param>
/// <param name="GuideTree">The progressive alignment steps in execution order: <c>(a, b, c)</c> = profiles a and b merged into c (indices 0 .. n-1 are the sequences, n .. 2n-2 the internal nodes), as <c>create_tasks</c> lists them.</param>
public sealed record KalignRunResult(MultipleAlignmentResult Alignment, int[] KeptIndices, string Alphabet, (int A, int B, int C)[] GuideTree);

public static partial class SequenceAligner
{
    #region Kalign 3.3.5

    // Port of Kalign 3.3.5 (Lassmann 2019, Bioinformatics 36:1928; Debian 1:3.3.5-1, kalign_3.3.5.orig.tar.gz, lib/src/*.c) for `kalign -i x.fa -o y.fa [--type T] [--gpo X --gpe Y --tgpe Z] [-n N]`:
    //   msa_op.c detect_alphabet (log-likelihood of the letter counts: DNA = acgtun, protein = acdefghiklmnpqrstvwy), kalign_essential_input_check (sequences without residues removed),
    //   alphabet.c (the DNA alphabet ACGTN, the 13-letter reduced protein alphabet for the guide tree, the 23-letter ambiguous alphabet for the alignment),
    //   sequence_distance.c + bpm.c (distance = the bit-parallel Myers edit distance of the shorter sequence (at most 255 residues) against the best window of the longer; the AVX2 variant bpm_256 is what the
    //   Debian `kalign` dispatches to on a machine with AVX2), pick_anchor.c, bisectingKmeans.c (UPGMA below 100 sequences, otherwise bisecting 2-means on the distances to 32 anchors) + euclidean_dist.c,
    //   aln_param.c (matrices and default gap penalties), aln_run.c (progressive alignment along the tree, serial and OpenMP variants), aln_setup.c (profiles, update_n), aln_controller.c (Hirschberg
    //   divide and conquer), aln_seqseq.c / aln_seqprofile.c / aln_profileprofile.c (the three-state affine dynamic programs, forward / backward / meet-up), weave_alignment.c (make_seq).
    // All arithmetic that is single precision in the C (float profiles, scores, distances) is single precision here.

    /// <summary>
    /// Kalign 3.3.5's multiple alignment of nucleotide or protein sequences (<c>kalign -i in.fa -o out.fa</c>): bit-parallel (Myers) edit distances, a UPGMA / bisecting k-means guide
    /// tree, progressive profile alignment with the affine three-state dynamic program in Hirschberg's linear space. See <see cref="KalignRun"/> for the details and
    /// <see cref="KalignOptions"/> for the parameters.
    /// </summary>
    /// <param name="sequences">Sequences (letters; characters other than letters are removed as gaps; at least two must keep a residue).</param>
    /// <param name="options">Kalign parameters (default: the program's defaults).</param>
    /// <returns>Rows in input order (sequences without residues are dropped, as the program does), the consensus and the natural affine sum-of-pairs score (not reported, 0 and NaN, for sets so large that the exact integer score overflows).</returns>
    /// <exception cref="ArgumentException">Fewer than two sequences with residues, a letter outside the alphabet of the detected type, an alignment type that does not fit the sequences (the program stops with an error) or an unusable parameter.</exception>
    public static MultipleAlignmentResult MultipleAlignKalign(IEnumerable<string> sequences, KalignOptions? options = null)
        => KalignRun(sequences, options).Alignment;

    /// <summary>
    /// Kalign 3.3.5 (<c>kalign_run</c>, <c>lib/src/aln_wrap.c</c>) with the intermediate decisions: the alignment, which input sequences were kept, the detected alphabet and the progressive
    /// alignment steps. The oracle is the installed Debian <c>kalign</c> 3.3.5 (AVX2 build): rows equal on every committed set (scripts/alignment/generate_kalign_oracle.py, kalign_oracle.tsv).
    /// </summary>
    /// <remarks>
    /// Faithfully reproduced behaviours of the program (all confirmed by the oracle): the alphabet is detected from the letter counts of the whole input (so a set of short protein sequences made of the
    /// letters A, C, G, T only is read as DNA); letters of the input keep their case in the output; every character that is not a letter is a gap and is removed before alignment; sequences
    /// with no residue are dropped; the gap-open penalty of a profile is never charged for the columns inserted into a profile (the open / close flags of <c>add_gap_info_to_path_n</c> are never set in 3.3.5, its
    /// marking loop is dead code), only the extension penalties; with one thread the shorter-first swap of <c>do_align</c> is not made.
    /// The distance estimation limits the shorter sequence to its first 255 residues (the AVX2 <c>bpm_256</c> of the Debian build; the portable build uses 1 024).
    /// </remarks>
    /// <exception cref="ArgumentException">See <see cref="MultipleAlignKalign"/>.</exception>
    public static KalignRunResult KalignRun(IEnumerable<string> sequences, KalignOptions? options = null)
    {
        var o = options ?? new KalignOptions();
        var list = new List<string>();
        foreach (var s in sequences ?? throw new ArgumentNullException(nameof(sequences)))
            list.Add(s ?? throw new ArgumentException("A sequence is null.", nameof(sequences)));
        if (o.Threads < 1)
            throw new ArgumentException("Number of threads has to be >= 1.", nameof(options));
        using var cancelScope = AlignmentCancellation.Enter(o.CancellationToken);
        var engine = new KalignEngine(list, o);
        engine.Run();
        var rows = engine.Rows;
        var sp = engine.IsDna ? MatrixScoring.Emboss(SubstitutionMatrix.EdnaFull, 10, 0.5) : MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.2);
        string consensus = BuildConsensus(rows.ToList(), rows[0].Length, anySymbol: true);
        string[] scored = rows.Select(r => r.ToUpperInvariant()).ToArray();
        MultipleAlignmentResult aln;
        try
        {
            double exact = MatrixSumOfPairsScore(scored, sp, MsaGapModel.Affine);
            aln = new MultipleAlignmentResult(rows, consensus, (int)Math.Round(exact, MidpointRounding.AwayFromZero)) { ExactTotalScore = exact };
        }
        catch (ArgumentException)
        {
            // thousands of sequences: the exact integer sum of pairs overflows its range; the program itself reports no score, so none is reported (0 / NaN)
            aln = new MultipleAlignmentResult(rows, consensus, 0) { ExactTotalScore = double.NaN };
        }
        return new KalignRunResult(aln, engine.Kept, engine.IsDna ? "dna" : "protein", engine.Tasks);
    }

    /// <summary>
    /// The bit-parallel edit distance of Kalign (<c>bpm_256</c>, <c>lib/src/bpm.c</c>, Myers 1999 J. ACM 46:395): the minimum over the end positions of <paramref name="text"/> of the edit distance
    /// between <paramref name="pattern"/> (its first 255 symbols; free start in the text) and a suffix of the text prefix - i.e. the best approximate occurrence of the pattern in the text.
    /// Symbols are the internal codes of the alphabet (0 ... 12).
    /// </summary>
    internal static int KalignBpm(byte[] text, int n, byte[] pattern, int m) => KalignEngine.Bpm256(text, n, pattern, m);

    internal sealed class KalignEngine
    {
        private const float FltMax = float.MaxValue;
        private const int Sigma = 13;

        private readonly List<string> _input;
        private readonly KalignOptions _opt;
        private int _numseq;
        private char[][] _letters = Array.Empty<char[]>();
        private byte[][] _s = Array.Empty<byte[]>();
        private int[][] _gaps = Array.Empty<int[]>();
        private int[] _len = Array.Empty<int>();
        private int[] _kept = Array.Empty<int>();
        private int[][] _sip = Array.Empty<int[]>();
        private int[] _nsip = Array.Empty<int>();
        private int[] _plen = Array.Empty<int>();
        private float[][] _profile = Array.Empty<float[]>();
        private (int A, int B, int C)[] _tasks = Array.Empty<(int, int, int)>();

        internal string[] Rows { get; private set; } = Array.Empty<string>();
        internal int[] Kept => _kept;
        internal bool IsDna { get; private set; }
        internal (int A, int B, int C)[] Tasks => _tasks;

        internal KalignEngine(List<string> input, KalignOptions options)
        {
            _input = input;
            _opt = options;
        }

        // ------------------------------------------------------------------ input (msa_io.c read_fasta, msa_op.c)

        internal void Run()
        {
            int n = _input.Count;
            var letterFreq = new int[128];
            var letters = new List<char[]>();
            for (int i = 0; i < n; i++)
            {
                var sb = new List<char>();
                string line = _input[i];
                for (int k = 0; k < line.Length; k++)
                {
                    char ch = line[k];
                    if (ch >= 128)
                        throw new ArgumentException($"Sequence {i + 1} contains the non-ASCII character U+{(int)ch:X4}; Kalign works on ASCII letters.");
                    if (ch < 32 || ch == 127)
                        break; // read_file_stdin: a line ends at its first control character
                    letterFreq[ch]++;
                    if (IsAlpha(ch))
                        sb.Add(ch);
                    // ispunct: a gap (dropped by dealign_msa); anything else (digits, blanks) is ignored
                }
                letters.Add(sb.ToArray());
            }

            // kalign_essential_input_check
            RequireProgramInput(n, "Kalign", "No alignment", "sequences");
            if (n < 2)
                throw new ArgumentException($"only {n} sequences found.");
            var keep = new List<int>();
            for (int i = 0; i < n; i++)
                if (letters[i].Length > 0)
                    keep.Add(i);
            if (keep.Count < 2)
                throw new ArgumentException($"only {keep.Count} sequences found.");
            _kept = keep.ToArray();
            _numseq = keep.Count;
            _letters = keep.Select(i => letters[i]).ToArray();
            _len = _letters.Select(l => l.Length).ToArray();
            _gaps = _len.Select(l => new int[l + 2]).ToArray();

            IsDna = DetectDna(letterFreq);

            // alphabets
            AlnParam ap = CreateParam(IsDna);
            _s = new byte[_numseq][];
            if (IsDna)
                ConvertToInternal(CreateAlphabet(0, out _));
            else
                ConvertToInternal(CreateAlphabet(1, out _));

            // guide tree
            AlignmentCancellation.ThrowIfRequested();
            BuildTreeKmeans();

            if (!IsDna)
                ConvertToInternal(CreateAlphabet(2, out _));

            // progressive alignment
            _sip = new int[2 * _numseq - 1][];
            _nsip = new int[2 * _numseq - 1];
            _plen = new int[2 * _numseq - 1];
            _profile = new float[2 * _numseq - 1][];
            for (int i = 0; i < _numseq; i++)
            {
                _sip[i] = new[] { i };
                _nsip[i] = 1;
            }
            CreateMsaTree(ap);
            Rows = MakeRows();
        }

        private static bool IsAlpha(char c) => (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');

        /// <summary>msa_op.c detect_alphabet: log-likelihood of the letter counts under a DNA and a protein letter model.</summary>
        private static bool DetectDna(int[] letterFreq)
        {
            var dna = new double[128];
            var protein = new double[128];
            for (int i = 0; i < 128; i++)
            {
                dna[i] = Math.Log(0.0001 * 1.0 / 116.0);
                protein[i] = Math.Log(0.0001 * 1.0 / 88.0);
            }
            foreach (char c in "acgtunACGTUN")
                dna[c] = Math.Log(0.9999 * 1.0 / 12.0);
            foreach (char c in "acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY")
                protein[c] = Math.Log(0.9999 * 1.0 / 40.0);
            double dnaProb = 0.0, protProb = 0.0;
            for (int i = 0; i < 128; i++)
            {
                if (letterFreq[i] != 0)
                {
                    dnaProb += dna[i] * (double)letterFreq[i];
                    protProb += protein[i] * (double)letterFreq[i];
                }
            }
            if (dnaProb == protProb)
                throw new ArgumentException("Unable to determine what alphabet to use (the letter counts fit DNA and protein equally well).");
            return dnaProb > protProb;
        }

        // ------------------------------------------------------------------ alphabet.c

        private static sbyte[] CreateAlphabet(int kind, out int length)
        {
            var t = new sbyte[128];
            for (int i = 0; i < 128; i++)
                t[i] = -1;
            void Merge(char x, char y)
            {
                int min = Math.Min(t[x], t[y]);
                t[x] = (sbyte)min;
                t[y] = (sbyte)min;
            }
            switch (kind)
            {
                case 0: // ALPHA_defDNA
                {
                    const string dnacode = "ACGTUNRYSWKMBDHV";
                    for (int i = 0; i < 16; i++)
                        t[dnacode[i]] = (sbyte)i;
                    Merge('U', 'T');
                    foreach (char c in "RYSWKMBDHV")
                        Merge('N', c);
                    break;
                }
                case 1: // ALPHA_redPROTEIN
                {
                    const string aacode = "ACDEFGHIKLMNPQRSTVWY";
                    int code = 0;
                    for (int i = 0; i < 20; i++)
                        t[aacode[i]] = (sbyte)code++;
                    t['B'] = (sbyte)code++;
                    t['Z'] = (sbyte)code++;
                    t['X'] = (sbyte)code;
                    Merge('L', 'M'); Merge('I', 'V'); Merge('K', 'R'); Merge('E', 'Q'); Merge('A', 'S'); Merge('A', 'T'); Merge('S', 'T');
                    Merge('N', 'D'); Merge('F', 'Y');
                    Merge('B', 'N'); Merge('B', 'D'); Merge('Z', 'E'); Merge('Z', 'Q');
                    break;
                }
                default: // ALPHA_ambigiousPROTEIN
                {
                    const string aacode = "ARNDCQEGHILKMFPSTWYVBZX";
                    for (int i = 0; i < 23; i++)
                        t[aacode[i]] = (sbyte)i;
                    t['U'] = 22;
                    break;
                }
            }
            // clean_and_set_to_extern
            var trans = new sbyte[32];
            for (int i = 0; i < 32; i++)
                trans[i] = -1;
            for (int i = 64; i < 96; i++)
                if (t[i] != -1)
                    trans[t[i]] = 1;
            int next = 0;
            for (int i = 0; i < 32; i++)
                if (trans[i] == 1)
                    trans[i] = (sbyte)next++;
            length = next;
            for (int i = 64; i < 96; i++)
            {
                if (t[i] != -1)
                {
                    t[i] = trans[t[i]];
                    t[i + 32] = t[i];
                }
            }
            return t;
        }

        private void ConvertToInternal(sbyte[] table)
        {
            for (int i = 0; i < _numseq; i++)
            {
                var s = _s[i] ?? new byte[_len[i]];
                for (int j = 0; j < _len[i]; j++)
                {
                    int c = table[_letters[i][j]];
                    if (c == -1)
                        throw new ArgumentException($"The character '{_letters[i][j]}' of sequence {_kept[i] + 1} is not in the alphabet of the detected sequence type (Kalign leaves its code undefined).");
                    s[j] = (byte)c;
                }
                _s[i] = s;
            }
        }

        // ------------------------------------------------------------------ aln_param.c

        private sealed class AlnParam
        {
            public float[][] Subm = new float[23][];
            public float Gpo, Gpe, Tgpe;
        }

        private AlnParam CreateParam(bool dna)
        {
            var ap = new AlnParam();
            for (int i = 0; i < 23; i++)
                ap.Subm[i] = new float[23];
            var type = _opt.Type;
            if (dna)
            {
                switch (type)
                {
                    case KalignType.Dna: SetDna(ap, 0f); break;
                    case KalignType.DnaInternal: SetDna(ap, 8f); break;
                    case KalignType.Rna: SetRna(ap); break;
                    case KalignType.Protein:
                        throw new ArgumentException("Detected DNA sequences but --type protein option was selected.");
                    default: SetRna(ap); break;
                }
            }
            else
            {
                switch (type)
                {
                    case KalignType.Protein: SetCorBlosum(ap); break;
                    case KalignType.ProteinDivergent: SetGon250(ap); break;
                    case KalignType.Dna:
                        throw new ArgumentException("Detected protein sequences but --type dna option was selected.");
                    case KalignType.DnaInternal:
                        throw new ArgumentException("Detected protein sequences but --type internal  option was selected.");
                    case KalignType.Rna:
                        throw new ArgumentException("Detected protein sequences but --type rna option was selected.");
                    default: SetCorBlosum(ap); break;
                }
            }
            if (_opt.GapOpen >= 0.0)
                ap.Gpo = _opt.GapOpen;
            if (_opt.GapExtend >= 0.0)
                ap.Gpe = _opt.GapExtend;
            if (_opt.GapExtend >= 0.0) // sic: aln_param.c tests gpe here
                ap.Tgpe = _opt.TerminalGapExtend;
            return ap;
        }

        private static void SetDna(AlnParam ap, float tgpe)
        {
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 5; j++)
                    ap.Subm[i][j] = i == j ? 5f : -4f;
            ap.Gpo = 8; ap.Gpe = 6; ap.Tgpe = tgpe;
        }

        private static void SetRna(AlnParam ap)
        {
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 5; j++)
                    ap.Subm[i][j] = 283.0f;
            double[,] add = { { 91, -114, -31, -123 }, { -114, 100, -125, -31 }, { -31, -125, 100, -114 }, { -123, -31, -114, 91 } };
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    ap.Subm[i][j] = (float)(ap.Subm[i][j] + add[i, j]);
            ap.Gpo = 217.0f; ap.Gpe = 39.4f; ap.Tgpe = 292.6f;
        }

        private static void SetCorBlosum(AlnParam ap)
        {
            for (int i = 0; i < 23; i++)
                for (int j = 0; j < 23; j++)
                    ap.Subm[i][j] = KalignTables.CorBlosum66_13plus[i * 23 + j];
            ap.Gpo = 5.5f; ap.Gpe = 2.0f; ap.Tgpe = 1.0f;
        }

        private static void SetGon250(AlnParam ap)
        {
            for (int i = 0; i < 23; i++)
                for (int j = 0; j < 23; j++)
                    ap.Subm[i][j] = KalignTables.Gon250[i * 23 + j];
            ap.Gpo = 55; ap.Gpe = 8; ap.Tgpe = 4;
        }

        // ------------------------------------------------------------------ bpm.c

        private readonly struct V256
        {
            public readonly ulong A, B, C, D;
            public V256(ulong a, ulong b, ulong c, ulong d) { A = a; B = b; C = c; D = d; }
            public static V256 operator |(V256 x, V256 y) => new(x.A | y.A, x.B | y.B, x.C | y.C, x.D | y.D);
            public static V256 operator &(V256 x, V256 y) => new(x.A & y.A, x.B & y.B, x.C & y.C, x.D & y.D);
            public static V256 operator ^(V256 x, V256 y) => new(x.A ^ y.A, x.B ^ y.B, x.C ^ y.C, x.D ^ y.D);
            public static V256 operator ~(V256 x) => new(~x.A, ~x.B, ~x.C, ~x.D);
            public static V256 Shl1(V256 x) => new(x.A << 1, (x.B << 1) | (x.A >> 63), (x.C << 1) | (x.B >> 63), (x.D << 1) | (x.C >> 63));
            public static V256 Add(V256 x, V256 y)
            {
                ulong a = x.A + y.A; ulong c = a < x.A ? 1UL : 0UL;
                ulong b = x.B + y.B; ulong c2 = b < x.B ? 1UL : 0UL; ulong b2 = b + c; if (b2 < b) c2 = 1UL;
                ulong cc = x.C + y.C; ulong c3 = cc < x.C ? 1UL : 0UL; ulong cc2 = cc + c2; if (cc2 < cc) c3 = 1UL;
                ulong d = x.D + y.D + c3;
                return new V256(a, b2, cc2, d);
            }
            public bool AnyAnd(V256 mask) => ((A & mask.A) | (B & mask.B) | (C & mask.C) | (D & mask.D)) != 0;
        }

        /// <summary>bpm.c bpm_256: Myers' bit-vector algorithm on 256-bit words (pattern at most 255 symbols); returns the smallest value of the last row of the edit-distance matrix.</summary>
        internal static int Bpm256(byte[] t, int n, byte[] p, int m)
        {
            if (m > 255)
                m = 255;
            var f = new ulong[Sigma * 4];
            for (int i = 0; i < m; i++)
                f[p[i] * 4 + (i >> 6)] |= 1UL << (i & 63);
            var bTab = new V256[Sigma];
            for (int c = 0; c < Sigma; c++)
                bTab[c] = new V256(f[c * 4], f[c * 4 + 1], f[c * 4 + 2], f[c * 4 + 3]);
            int diff = m, k = m;
            var vp = new V256(ulong.MaxValue, ulong.MaxValue, ulong.MaxValue, ulong.MaxValue);
            var vn = new V256(0, 0, 0, 0);
            int bit = m - 1;
            var mask = new V256(bit < 64 ? 1UL << bit : 0, bit >= 64 && bit < 128 ? 1UL << (bit - 64) : 0, bit >= 128 && bit < 192 ? 1UL << (bit - 128) : 0, bit >= 192 ? 1UL << (bit - 192) : 0);
            for (int i = 0; i < n; i++)
            {
                var x = bTab[t[i]] | vn;
                var xmm1 = x & vp;
                var xmm2 = V256.Add(vp, xmm1);
                xmm1 = xmm2 ^ vp;
                var d0 = xmm1 | x;
                var hn = vp & d0;
                var hp = vn | ~(vp | d0);
                x = V256.Shl1(hp);
                vn = x & d0;
                vp = V256.Shl1(hn) | ~(x | d0);
                if (hp.AnyAnd(mask))
                    diff++;
                if (hn.AnyAnd(mask))
                    diff--;
                if (diff < k)
                    k = diff;
            }
            return k;
        }

        // ------------------------------------------------------------------ sequence_distance.c

        private float CalcDistance(int a, int b)
        {
            int la = _len[a], lb = _len[b];
            return la > lb ? Bpm256(_s[a], la, _s[b], lb) : Bpm256(_s[b], lb, _s[a], la);
        }

        /// <summary>d_estimation(pair = 1): the distance matrix of the sampled sequences (the program computes every ordered pair and symmetrises with the later write: entry (a, b), a &lt; b, = distance(b, a)).</summary>
        private float[][] DistancePairs(int[] samples, int num)
        {
            var dm = new float[num][];
            for (int i = 0; i < num; i++)
                dm[i] = new float[num];
            for (int i = 0; i < num; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    AlignmentCancellation.ThrowIfRequested();
                    float dist = CalcDistance(samples[i], samples[j]);
                    int s = (_len[samples[i]] + _len[samples[j]]) / 2;
                    float add = (float)(Math.Min(10000.0, s) / 10000.0);
                    dist += add;
                    dm[i][j] = dist;
                    dm[j][i] = dist;
                }
            }
            return dm;
        }

        // ------------------------------------------------------------------ pick_anchor.c, bisectingKmeans.c, euclidean_dist.c

        private sealed class Node
        {
            public Node? Left, Right;
            public int Id = -1;
        }

        private sealed class KmeansResult
        {
            public int[] Sl = Array.Empty<int>();
            public int[] Sr = Array.Empty<int>();
            public int Nl, Nr;
            public float Score;
        }

        private float[][] _anchorDm = Array.Empty<float[]>();

        private void BuildTreeKmeans()
        {
            int numseq = _numseq;
            var samples = new int[numseq];
            for (int i = 0; i < numseq; i++)
                samples[i] = i;

            if (numseq >= 100)
            {
                int numAnchors = Math.Min(32, numseq);
                var anchors = PickAnchors(numAnchors);
                int padded = ((numAnchors + 7) / 8) * 8;
                _anchorDm = new float[numseq][];
                for (int i = 0; i < numseq; i++)
                {
                    AlignmentCancellation.ThrowIfRequested();
                    _anchorDm[i] = new float[padded];
                    for (int j = 0; j < numAnchors; j++)
                        _anchorDm[i][j] = CalcDistance(i, anchors[j]);
                }
            }
            var root = BisectingKmeans(samples, numseq);
            LabelInternal(root, numseq);
            var tasks = new List<(int, int, int)>();
            CreateTasks(root, tasks);
            tasks.Sort((x, y) => x.Item3.CompareTo(y.Item3)); // sort_tasks(TASK_ORDER_TREE): by c (distinct)
            _tasks = tasks.ToArray();
        }

        /// <summary>pick_anchor.c: the sequences sorted by length (descending, glibc qsort = merge sort with the program's non-strict comparator) sampled with a stride.</summary>
        private int[] PickAnchors(int numAnchor)
        {
            var order = new int[_numseq];
            for (int i = 0; i < _numseq; i++)
                order[i] = i;
            var tmp = new int[_numseq];
            MsortByLength(order, tmp, 0, _numseq);
            int stride = _numseq / numAnchor;
            var anchors = new int[numAnchor];
            for (int i = 0; i < numAnchor; i++)
                anchors[i] = order[i * stride];
            return anchors;
        }

        // glibc msort_with_tmp with sort_by_len as comparator: returns -1 when len(a) > len(b) and 1 otherwise (also for equal lengths, which reverses ties in a merge).
        private void MsortByLength(int[] b, int[] tmp, int start, int n)
        {
            if (n <= 1)
                return;
            int n1 = n / 2, n2 = n - n1;
            int b1 = start, b2 = start + n1;
            MsortByLength(b, tmp, b1, n1);
            MsortByLength(b, tmp, b2, n2);
            int t = start;
            while (n1 > 0 && n2 > 0)
            {
                int cmp = _len[b[b1]] > _len[b[b2]] ? -1 : 1;
                if (cmp <= 0) { tmp[t++] = b[b1++]; n1--; }
                else { tmp[t++] = b[b2++]; n2--; }
            }
            while (n1 > 0) { tmp[t++] = b[b1++]; n1--; }
            Array.Copy(tmp, start, b, start, n - n2);
        }

        private Node BisectingKmeans(int[] samples, int numSamples)
        {
            int numAnchors = Math.Min(32, _numseq);
            if (numSamples < 100)
            {
                var dm = DistancePairs(samples, numSamples);
                return Upgma(dm, samples, numSamples);
            }
            KmeansResult? best = null;
            var res = new KmeansResult?[4];
            int tries = Math.Min(40, numSamples);
            int step = numSamples / tries;
            for (int i = 0; i < tries; i += 4)
            {
                int change = 0;
                for (int j = 0; j < 4; j++)
                    res[j] = Split(samples, numAnchors, numSamples, (i + j) * step, res[j]);
                for (int j = 0; j < 4; j++)
                {
                    if (best == null)
                    {
                        change++;
                        best = res[j];
                        res[j] = null;
                    }
                    else if (best.Score > res[j]!.Score)
                    {
                        var tmp = best;
                        best = res[j];
                        res[j] = tmp;
                        change++;
                    }
                }
                if (change == 0)
                    break;
            }
            var sl = best!.Sl;
            var sr = best.Sr;
            int numL = best.Nl, numR = best.Nr;
            var n = new Node();
            n.Left = BisectingKmeans(sl, numL);
            n.Right = BisectingKmeans(sr, numR);
            return n;
        }

        private KmeansResult Split(int[] samples, int numAnchors, int numSamples, int seedPick, KmeansResult? reuse)
        {
            int numVar = ((numAnchors + 7) / 8) * 8;
            var wr = new float[numVar];
            var wl = new float[numVar];
            var cr = new float[numVar];
            var cl = new float[numVar];
            var w = new float[numVar];
            var res = reuse ?? new KmeansResult { Sl = new int[numSamples], Sr = new int[numSamples] };
            if (res.Sl.Length < numSamples)
            {
                res.Sl = new int[numSamples];
                res.Sr = new int[numSamples];
            }
            res.Score = FltMax;
            var sl = res.Sl;
            var sr = res.Sr;
            int numL, numR;
            float score;
            for (int i = 0; i < numSamples; i++)
            {
                int s = samples[i];
                for (int j = 0; j < numAnchors; j++)
                    w[j] += _anchorDm[s][j];
            }
            for (int j = 0; j < numAnchors; j++)
                w[j] /= (float)numSamples;
            {
                int s = samples[seedPick];
                for (int j = 0; j < numAnchors; j++)
                    cl[j] = _anchorDm[s][j];
            }
            for (int j = 0; j < numAnchors; j++)
                cr[j] = w[j] - (cl[j] - w[j]);
            bool differ = false;
            for (int j = 0; j < numAnchors; j++)
            {
                if (MathF.Abs(cl[j] - cr[j]) > 1.0E-6)
                {
                    differ = true;
                    break;
                }
            }
            if (!differ)
            {
                score = 0.0F;
                numL = 0;
                numR = 0;
                for (int i = 0; i < numSamples / 2; i++)
                    sl[numL++] = samples[i];
                for (int i = numSamples / 2; i < numSamples; i++)
                    sr[numR++] = samples[i];
            }
            else
            {
                int stop = 0;
                float[]? acc = null;
                while (true)
                {
                    stop++;
                    if (stop == 10000)
                        throw new InvalidOperationException("Kalign: the k-means split did not converge (the program stops with 'Failed.').");
                    numL = 0;
                    numR = 0;
                    for (int i = 0; i < numAnchors; i++)
                    {
                        wr[i] = 0.0F;
                        wl[i] = 0.0F;
                    }
                    score = 0.0f;
                    for (int i = 0; i < numSamples; i++)
                    {
                        int s = samples[i];
                        float dl = Edist256(_anchorDm[s], cl, numAnchors);
                        float dr = Edist256(_anchorDm[s], cr, numAnchors);
                        score += Math.Min(dl, dr);
                        // right when strictly closer to the right centre, or on a tie for odd positions
                        if (dr < dl || (dr == dl && (i & 1) != 0))
                        {
                            acc = wr;
                            sr[numR++] = s;
                        }
                        else
                        {
                            acc = wl;
                            sl[numL++] = s;
                        }
                        for (int j = 0; j < numAnchors; j++)
                            acc[j] += _anchorDm[s][j];
                    }
                    for (int j = 0; j < numAnchors; j++)
                    {
                        wl[j] /= (float)numL;
                        wr[j] /= (float)numR;
                    }
                    bool moved = false;
                    for (int j = 0; j < numAnchors; j++)
                    {
                        if (wl[j] != cl[j])
                        {
                            moved = true;
                            break;
                        }
                        if (wr[j] != cr[j])
                        {
                            moved = true;
                            break;
                        }
                    }
                    if (moved)
                    {
                        var t = cl; cl = wl; wl = t;
                        t = cr; cr = wr; wr = t;
                    }
                    else
                        break;
                }
            }
            res.Nl = numL;
            res.Nr = numR;
            res.Score = score;
            return res;
        }

        /// <summary>euclidean_dist.c edist_256: eight-lane single-precision sum of squares, reduced as hsum256_ps_avx does, then the square root.</summary>
        private static float Edist256(float[] a, float[] b, int len)
        {
            Span<float> r = stackalloc float[8];
            for (int i = 0; i < len; i += 8)
            {
                for (int k = 0; k < 8; k++)
                {
                    float t = a[i + k] - b[i + k];
                    float sq = t * t;
                    r[k] = r[k] + sq;
                }
            }
            // hsum256_ps_avx: low 128 + high 128, then hsum_ps_sse3
            float v0 = r[0] + r[4], v1 = r[1] + r[5], v2 = r[2] + r[6], v3 = r[3] + r[7];
            float s0 = v0 + v1; // movehdup: [v1, v1, v3, v3]; sums = v + shuf
            float s2 = v2 + v3;
            float d = s0 + s2;
            return MathF.Sqrt(d);
        }

        /// <summary>bisectingKmeans.c upgma: UPGMA on the sampled distance matrix with the program's update (average + 0.001) in single precision.</summary>
        private static Node Upgma(float[][] dm, int[] samples, int numseq)
        {
            int numprofiles = (numseq << 1) - 1;
            var asArr = new int[numseq];
            for (int i = numseq; i-- > 0;)
                asArr[i] = i + 1;
            var tree = new Node?[numseq];
            for (int i = 0; i < numseq; i++)
                tree[i] = new Node { Id = samples[i] };
            int cnode = numseq;
            int nodeA = 0, nodeB = 0;
            while (cnode != numprofiles)
            {
                float max = FltMax;
                for (int i = 0; i < numseq - 1; i++)
                {
                    if (asArr[i] != 0)
                    {
                        for (int j = i + 1; j < numseq; j++)
                        {
                            if (asArr[j] != 0 && dm[i][j] < max)
                            {
                                max = dm[i][j];
                                nodeA = i;
                                nodeB = j;
                            }
                        }
                    }
                }
                var tmp = new Node { Left = tree[nodeA], Right = tree[nodeB] };
                tree[nodeA] = tmp;
                tree[nodeB] = null;
                asArr[nodeA] = cnode + 1;
                asArr[nodeB] = 0;
                cnode++;
                for (int j = numseq; j-- > 0;)
                {
                    if (j != nodeB)
                        dm[nodeA][j] = (dm[nodeA][j] + dm[nodeB][j]) * 0.5F + 0.001F;
                }
                dm[nodeA][nodeA] = 0.0F;
                for (int j = numseq; j-- > 0;)
                    dm[j][nodeA] = dm[nodeA][j];
            }
            return tree[nodeA]!;
        }

        private static int LabelInternal(Node n, int label)
        {
            if (n.Left != null)
                label = LabelInternal(n.Left, label);
            if (n.Right != null)
                label = LabelInternal(n.Right, label);
            if (n.Id == -1)
            {
                n.Id = label;
                label++;
            }
            return label;
        }

        private static void CreateTasks(Node n, List<(int, int, int)> t)
        {
            if (n.Left != null && n.Right != null)
                t.Add((n.Left.Id, n.Right.Id, n.Id));
            if (n.Left != null)
                CreateTasks(n.Left, t);
            if (n.Right != null)
                CreateTasks(n.Right, t);
        }

        // ------------------------------------------------------------------ aln_run.c

        private struct State
        {
            public float A, Ga, Gb;
        }

        private sealed class AlnMem
        {
            public AlnParam Ap = null!;
            public State[] F = Array.Empty<State>(), B = Array.Empty<State>();
            public int[] Path = Array.Empty<int>(), TmpPath = Array.Empty<int>();
            public byte[]? Seq1, Seq2;
            public float[]? Prof1, Prof2;
            public int Sip;
            public int StartA, StartA2, StartB, EndA, EndA2, EndB, LenA, LenB;
        }

        private void CreateMsaTree(AlnParam ap)
        {
            var active = new bool[2 * _numseq - 1];
            for (int i = 0; i < _numseq; i++)
                active[i] = true;
            RecursiveAln(ap, active, _tasks.Length - 1);
        }

        private void RecursiveAln(AlnParam ap, bool[] active, int c)
        {
            var (ta, tb, tc) = _tasks[c];
            if (!active[ta] && ta >= _numseq)
                RecursiveAln(ap, active, ta - _numseq);
            if (!active[tb] && tb >= _numseq)
                RecursiveAln(ap, active, tb - _numseq);
            AlignmentCancellation.ThrowIfRequested();
            DoAlign(ap, c);
            active[tc] = true;
        }

        private void DoAlign(AlnParam ap, int taskId)
        {
            bool parallel = _opt.Threads != 1;
            var (a, b, c) = _tasks[taskId];
            var m = new AlnMem { Ap = ap };
            if (_nsip[a] == 1)
            {
                m.LenA = _len[a];
                _profile[a] = MakeProfile(ap, _s[a], m.LenA);
            }
            else
            {
                m.LenA = _plen[a];
                SetGapPenalties(_profile[a], m.LenA, _nsip[b]);
            }
            if (_nsip[b] == 1)
            {
                m.LenB = _len[b];
                _profile[b] = MakeProfile(ap, _s[b], m.LenB);
            }
            else
            {
                m.LenB = _plen[b];
                SetGapPenalties(_profile[b], m.LenB, _nsip[a]);
            }
            InitAlnMem(m);
            int lenA, lenB;
            if (_nsip[a] == 1)
            {
                if (_nsip[b] == 1)
                {
                    if (!parallel || m.LenA < m.LenB)
                    {
                        m.Seq1 = _s[a];
                        m.Seq2 = _s[b];
                        m.Prof1 = null;
                        m.Prof2 = null;
                        Runner(m, parallel);
                    }
                    else
                    {
                        lenB = m.LenB;
                        lenA = m.LenA;
                        m.EndA = lenB;
                        m.EndB = lenA;
                        m.LenA = lenB;
                        m.LenB = lenA;
                        m.Seq1 = _s[b];
                        m.Seq2 = _s[a];
                        m.Prof1 = null;
                        m.Prof2 = null;
                        Runner(m, parallel);
                        MirrorPath(m, lenA, lenB);
                        m.LenA = lenA;
                        m.LenB = lenB;
                    }
                }
                else
                {
                    lenB = m.LenB;
                    lenA = m.LenA;
                    m.EndA = lenB;
                    m.EndB = lenA;
                    m.LenA = lenB;
                    m.LenB = lenA;
                    m.Seq1 = null;
                    m.Seq2 = _s[a];
                    m.Prof1 = _profile[b];
                    m.Prof2 = null;
                    m.Sip = _nsip[b];
                    Runner(m, parallel);
                    MirrorPath(m, lenA, lenB);
                    m.LenA = lenA;
                    m.LenB = lenB;
                }
            }
            else
            {
                if (_nsip[b] == 1)
                {
                    m.Seq1 = null;
                    m.Seq2 = _s[b];
                    m.Prof1 = _profile[a];
                    m.Prof2 = null;
                    m.Sip = _nsip[a];
                    Runner(m, parallel);
                }
                else
                {
                    if (m.LenA < m.LenB)
                    {
                        m.Seq1 = null;
                        m.Seq2 = null;
                        m.Prof1 = _profile[a];
                        m.Prof2 = _profile[b];
                        Runner(m, parallel);
                    }
                    else
                    {
                        lenB = m.LenB;
                        lenA = m.LenA;
                        m.EndA = lenB;
                        m.EndB = lenA;
                        m.LenA = lenB;
                        m.LenB = lenA;
                        m.Seq1 = null;
                        m.Seq2 = null;
                        m.Prof1 = _profile[b];
                        m.Prof2 = _profile[a];
                        Runner(m, parallel);
                        MirrorPath(m, lenA, lenB);
                        m.LenA = lenA;
                        m.LenB = lenB;
                    }
                }
            }
            AddGapInfoToPath(m);
            var tmp = new float[64 * (m.Path[0] + 2)];
            if (taskId != _tasks.Length - 1)
                UpdateProfile(_profile[a], _profile[b], tmp, ap, m.Path, _nsip[a], _nsip[b]);
            _profile[a] = Array.Empty<float>();
            _profile[b] = Array.Empty<float>();
            _profile[c] = tmp;
            MakeSeq(a, b, m.Path);
            _plen[c] = m.Path[0];
            _nsip[c] = _nsip[a] + _nsip[b];
            var sipc = new int[_nsip[a] + _nsip[b]];
            int g = 0;
            for (int j = _nsip[a]; j-- > 0;)
                sipc[g++] = _sip[a][j];
            for (int j = _nsip[b]; j-- > 0;)
                sipc[g++] = _sip[b][j];
            _sip[c] = sipc;
        }

        // ------------------------------------------------------------------ aln_setup.c

        private static void InitAlnMem(AlnMem m)
        {
            m.StartA = 0;
            m.StartB = 0;
            m.EndA = m.LenA;
            m.EndB = m.LenB;
            int g = Math.Max(m.LenA, m.LenB) + 2;
            m.F = new State[g];
            m.B = new State[g];
            m.F[0] = new State { A = 0.0F, Ga = -FltMax, Gb = -FltMax };
            m.B[0] = new State { A = 0.0F, Ga = -FltMax, Gb = -FltMax };
            m.Path = new int[m.LenA + m.LenB + 2];
            m.TmpPath = new int[m.LenA + m.LenB + 2];
            for (int i = 0; i < g; i++)
                m.Path[i] = -1;
        }

        private static float[] MakeProfile(AlnParam ap, byte[] seq, int len)
        {
            var prof = new float[(len + 2) * 64];
            int p = 64 * (len + 1);
            for (int i = 0; i < 64; i++)
                prof[p + i] = 0;
            prof[p + 23 + 32] = -ap.Gpo;
            prof[p + 24 + 32] = -ap.Gpe;
            prof[p + 25 + 32] = -ap.Tgpe;
            int idx = len;
            while (idx-- > 0)
            {
                p -= 64;
                for (int j = 0; j < 64; j++)
                    prof[p + j] = 0;
                int c = seq[idx];
                prof[p + c] += 1;
                p += 32;
                for (int j = 23; j-- > 0;)
                    prof[p + j] = ap.Subm[c][j];
                prof[p + 23] = -ap.Gpo;
                prof[p + 24] = -ap.Gpe;
                prof[p + 25] = -ap.Tgpe;
                p -= 32;
            }
            p -= 64;
            for (int i = 0; i < 64; i++)
                prof[p + i] = 0;
            prof[p + 23 + 32] = -ap.Gpo;
            prof[p + 24 + 32] = -ap.Gpe;
            prof[p + 25 + 32] = -ap.Tgpe;
            return prof;
        }

        private static void SetGapPenalties(float[] prof, int len, int nsip)
        {
            int p = 64 * (len + 1);
            prof[p + 27] = prof[p + 55] * (float)nsip;
            prof[p + 28] = prof[p + 56] * (float)nsip;
            prof[p + 29] = prof[p + 57] * (float)nsip;
            int i = len + 1;
            while (i-- > 0)
            {
                p -= 64;
                prof[p + 27] = prof[p + 55] * (float)nsip;
                prof[p + 28] = prof[p + 56] * (float)nsip;
                prof[p + 29] = prof[p + 57] * (float)nsip;
            }
        }

        /// <summary>aln_setup.c add_gap_info_to_path_n. In 3.3.5 the loop that should set the gap-open / close flags (4, 8, 16) tests <c>o_path[j] != 3</c> right after <c>o_path[j] = 3</c> and never runs; only the terminal flag (32) is set.</summary>
        private static void AddGapInfoToPath(AlnMem m)
        {
            var path = m.Path;
            var o = m.TmpPath;
            int lenA = m.LenA, lenB = m.LenB;
            for (int i = 0; i < lenA + lenB + 2; i++)
                o[i] = 0;
            int j = 1;
            int b;
            if (path[1] == -1)
            {
                o[j] = 2;
                j++;
            }
            else
            {
                if (path[1] != 1)
                {
                    for (int a = 0; a < path[1] - 1; a++)
                    {
                        o[j] = 1;
                        j++;
                    }
                    o[j] = 0;
                    j++;
                }
                else
                {
                    o[j] = 0;
                    j++;
                }
            }
            b = path[1];
            for (int i = 2; i <= lenA; i++)
            {
                if (path[i] == -1)
                {
                    o[j] = 2;
                    j++;
                }
                else
                {
                    if (path[i] - 1 != b && b != -1)
                    {
                        for (int a = 0; a < path[i] - b - 1; a++)
                        {
                            o[j] = 1;
                            j++;
                        }
                        o[j] = 0;
                        j++;
                    }
                    else
                    {
                        o[j] = 0;
                        j++;
                    }
                }
                b = path[i];
            }
            if (path[lenA] < lenB && path[lenA] != -1)
            {
                for (int a = 0; a < lenB - path[lenA]; a++)
                {
                    o[j] = 1;
                    j++;
                }
            }
            o[0] = j - 1;
            o[j] = 3;
            // the "add gap info" loop (while (o_path[j] != 3)) is dead code in 3.3.5
            int ii = 1;
            while (o[ii] != 0)
            {
                o[ii] |= 32;
                ii++;
            }
            ii = o[0];
            while (o[ii] != 0)
            {
                o[ii] |= 32;
                ii--;
                if (ii < 0)
                    throw new InvalidOperationException("Kalign: an alignment without a single aligned column makes the program read before the start of its path array (undefined behaviour); not reproduced.");
            }
            (m.Path, m.TmpPath) = (m.TmpPath, m.Path);
        }

        private static void UpdateProfile(float[] profa, float[] profb, float[] newp, AlnParam ap, int[] path, int sipa, int sipb)
        {
            int ia = 0, ib = 0, io = 0;
            for (int i = 64; i-- > 0;)
                newp[i] = profa[i] + profb[i];
            ia += 64; ib += 64; io += 64;
            int c = 1;
            while (path[c] != 3)
            {
                if (path[c] == 0)
                {
                    for (int i = 64; i-- > 0;)
                        newp[io + i] = profa[ia + i] + profb[ib + i];
                    ia += 64;
                    ib += 64;
                }
                if ((path[c] & 1) != 0)
                {
                    for (int i = 64; i-- > 0;)
                        newp[io + i] = profb[ib + i];
                    ib += 64;
                    GapColumn(newp, io, path[c], ap, sipa);
                }
                if ((path[c] & 2) != 0)
                {
                    for (int i = 64; i-- > 0;)
                        newp[io + i] = profa[ia + i];
                    ia += 64;
                    GapColumn(newp, io, path[c], ap, sipb);
                }
                io += 64;
                c++;
            }
            for (int i = 64; i-- > 0;)
                newp[io + i] = profa[ia + i] + profb[ib + i];
        }

        private static void GapColumn(float[] newp, int io, int flags, AlnParam ap, int sip)
        {
            float gp;
            float fs = (float)sip;
            if ((flags & 20) == 0)
            {
                if ((flags & 32) != 0)
                {
                    newp[io + 25] += fs;
                    gp = ap.Tgpe * fs;
                }
                else
                {
                    newp[io + 24] += fs;
                    gp = ap.Gpe * fs;
                }
                for (int j = 32; j < 55; j++)
                    newp[io + j] -= gp;
            }
            else
            {
                if ((flags & 16) != 0)
                {
                    if ((flags & 32) != 0)
                    {
                        newp[io + 25] += fs;
                        gp = ap.Tgpe * fs;
                        newp[io + 23] += fs;
                        gp += ap.Gpo * fs;
                    }
                    else
                    {
                        newp[io + 23] += fs;
                        gp = ap.Gpo * fs;
                    }
                    for (int j = 32; j < 55; j++)
                        newp[io + j] -= gp;
                }
                if ((flags & 4) != 0)
                {
                    if ((flags & 32) != 0)
                    {
                        newp[io + 25] += fs;
                        gp = ap.Tgpe * fs;
                        newp[io + 23] += fs;
                        gp += ap.Gpo * fs;
                    }
                    else
                    {
                        newp[io + 23] += fs;
                        gp = ap.Gpo * fs;
                    }
                    for (int j = 32; j < 55; j++)
                        newp[io + j] -= gp;
                }
            }
        }

        private static void MirrorPath(AlnMem m, int lenA, int lenB)
        {
            var apath = m.Path;
            var opath = m.TmpPath;
            for (int i = 0; i < lenA + 2; i++)
                opath[i] = -1;
            for (int i = 1; i <= lenB; i++)
                if (apath[i] != -1)
                    opath[apath[i]] = i;
            (m.Path, m.TmpPath) = (m.TmpPath, m.Path);
        }

        // ------------------------------------------------------------------ weave_alignment.c

        private void MakeSeq(int a, int b, int[] path)
        {
            var gapA = new int[path[0] + 1];
            var gapB = new int[path[0] + 1];
            int posa = 0, posb = 0;
            int c = 1;
            while (path[c] != 3)
            {
                if (path[c] == 0)
                {
                    posa++;
                    posb++;
                }
                else if ((path[c] & 1) != 0)
                {
                    gapA[posa] += 1;
                    posb++;
                }
                else if ((path[c] & 2) != 0)
                {
                    gapB[posb] += 1;
                    posa++;
                }
                c++;
            }
            for (int i = _nsip[a]; i-- > 0;)
                UpdateGaps(_len[_sip[a][i]], _gaps[_sip[a][i]], gapA);
            for (int i = _nsip[b]; i-- > 0;)
                UpdateGaps(_len[_sip[b][i]], _gaps[_sip[b][i]], gapB);
        }

        private static void UpdateGaps(int oldLen, int[] gis, int[] newgaps)
        {
            int relPos = 0;
            for (int i = 0; i <= oldLen; i++)
            {
                int add = 0;
                for (int j = relPos; j <= relPos + gis[i]; j++)
                    if (newgaps[j] != 0)
                        add += newgaps[j];
                relPos += gis[i] + 1;
                gis[i] += add;
            }
        }

        private string[] MakeRows()
        {
            var rows = new string[_numseq];
            for (int i = 0; i < _numseq; i++)
            {
                var sb = new StringBuilder();
                for (int j = 0; j < _len[i]; j++)
                {
                    sb.Append('-', _gaps[i][j]);
                    sb.Append(_letters[i][j]);
                }
                sb.Append('-', _gaps[i][_len[i]]);
                rows[i] = sb.ToString();
            }
            return rows;
        }

        // ------------------------------------------------------------------ aln_controller.c

        private static void Runner(AlnMem m, bool parallel)
        {
            if (parallel)
                AlnRunner(m);
            else
                AlnRunnerSerial(m);
        }

        private static void AlnRunner(AlnMem m)
        {
            if (m.EndA - m.StartA < 500)
                AlnRunnerSerial(m);
            if (m.StartA >= m.EndA)
                return;
            if (m.StartB >= m.EndB)
                return;
            RunnerBody(m, serial: false);
        }

        private static void AlnRunnerSerial(AlnMem m)
        {
            if (m.StartA >= m.EndA)
                return;
            if (m.StartB >= m.EndB)
                return;
            RunnerBody(m, serial: true);
        }

        private static void RunnerBody(AlnMem m, bool serial)
        {
            AlignmentCancellation.ThrowIfRequested();
            Span<float> inputStates = stackalloc float[6];
            Span<int> oldCor = stackalloc int[5];
            inputStates[0] = m.F[0].A;
            inputStates[1] = m.F[0].Ga;
            inputStates[2] = m.F[0].Gb;
            inputStates[3] = m.B[0].A;
            inputStates[4] = m.B[0].Ga;
            inputStates[5] = m.B[0].Gb;
            int mid = ((m.EndA - m.StartA) / 2) + m.StartA;
            oldCor[0] = m.StartA;
            oldCor[1] = m.EndA;
            oldCor[2] = m.StartB;
            oldCor[3] = m.EndB;
            oldCor[4] = mid;
            m.EndA = mid;
            m.StartA2 = mid;
            m.EndA2 = oldCor[1];
            int meet, transition;
            if (m.Seq1 != null)
            {
                SeqSeqForward(m);
                SeqSeqBackward(m);
                SeqSeqMeetup(m, oldCor, out meet, out transition);
            }
            else if (m.Prof2 != null)
            {
                ProfProfForward(m);
                ProfProfBackward(m);
                ProfProfMeetup(m, oldCor, out meet, out transition);
            }
            else
            {
                SeqProfForward(m);
                SeqProfBackward(m);
                SeqProfMeetup(m, oldCor, out meet, out transition);
            }
            Continue(m, inputStates, oldCor, meet, transition, serial);
        }

        private static void Sub(AlnMem m, bool serial)
        {
            if (serial)
                AlnRunnerSerial(m);
            else
                AlnRunner(m);
        }

        private static void Continue(AlnMem m, Span<float> input, Span<int> oldCor, int meet, int transition, bool serial)
        {
            var path = m.Path;
            switch (transition)
            {
                case 1: // a -> a
                    path[oldCor[4]] = meet;
                    path[oldCor[4] + 1] = meet + 1;
                    SetF(m, input[0], input[1], input[2]);
                    SetB(m, 0.0F, -FltMax, -FltMax);
                    m.StartA = oldCor[0]; m.EndA = oldCor[4] - 1; m.StartB = oldCor[2]; m.EndB = meet - 1;
                    Sub(m, serial);
                    m.StartA = oldCor[4] + 1; m.EndA = oldCor[1]; m.StartB = meet + 1; m.EndB = oldCor[3];
                    SetF(m, 0.0F, -FltMax, -FltMax);
                    SetB(m, input[3], input[4], input[5]);
                    Sub(m, serial);
                    break;
                case 2: // a -> ga
                    path[oldCor[4]] = meet;
                    SetF(m, input[0], input[1], input[2]);
                    SetB(m, 0.0F, -FltMax, -FltMax);
                    m.StartA = oldCor[0]; m.EndA = oldCor[4] - 1; m.StartB = oldCor[2]; m.EndB = meet - 1;
                    Sub(m, serial);
                    m.StartA = oldCor[4]; m.EndA = oldCor[1]; m.StartB = meet + 1; m.EndB = oldCor[3];
                    SetF(m, -FltMax, 0.0F, -FltMax);
                    SetB(m, input[3], input[4], input[5]);
                    Sub(m, serial);
                    break;
                case 3: // a -> gb
                    path[oldCor[4]] = meet;
                    SetF(m, input[0], input[1], input[2]);
                    SetB(m, 0.0F, -FltMax, -FltMax);
                    m.StartA = oldCor[0]; m.EndA = oldCor[4] - 1; m.StartB = oldCor[2]; m.EndB = meet - 1;
                    Sub(m, serial);
                    m.StartA = oldCor[4] + 1; m.EndA = oldCor[1]; m.StartB = meet; m.EndB = oldCor[3];
                    SetF(m, -FltMax, -FltMax, 0.0F);
                    SetB(m, input[3], input[4], input[5]);
                    Sub(m, serial);
                    break;
                case 5: // ga -> a
                    path[oldCor[4] + 1] = meet + 1;
                    SetF(m, input[0], input[1], input[2]);
                    SetB(m, -FltMax, 0.0F, -FltMax);
                    m.StartA = oldCor[0]; m.EndA = oldCor[4]; m.StartB = oldCor[2]; m.EndB = meet - 1;
                    Sub(m, serial);
                    m.StartA = oldCor[4] + 1; m.EndA = oldCor[1]; m.StartB = meet + 1; m.EndB = oldCor[3];
                    SetF(m, 0.0F, -FltMax, -FltMax);
                    SetB(m, input[3], input[4], input[5]);
                    Sub(m, serial);
                    break;
                case 6: // gb -> gb (the right half always goes through aln_runner in 3.3.5)
                    SetF(m, input[0], input[1], input[2]);
                    SetB(m, -FltMax, -FltMax, 0.0F);
                    m.StartA = oldCor[0]; m.EndA = oldCor[4] - 1; m.StartB = oldCor[2]; m.EndB = meet;
                    Sub(m, serial);
                    m.StartA = oldCor[4] + 1; m.EndA = oldCor[1]; m.StartB = meet; m.EndB = oldCor[3];
                    SetF(m, -FltMax, -FltMax, 0.0F);
                    SetB(m, input[3], input[4], input[5]);
                    AlnRunner(m);
                    break;
                case 7: // gb -> a
                    path[oldCor[4] + 1] = meet + 1;
                    SetF(m, input[0], input[1], input[2]);
                    SetB(m, -FltMax, -FltMax, 0.0F);
                    m.StartA = oldCor[0]; m.EndA = oldCor[4] - 1; m.StartB = oldCor[2]; m.EndB = meet;
                    Sub(m, serial);
                    m.StartA = oldCor[4] + 1; m.EndA = oldCor[1]; m.StartB = meet + 1; m.EndB = oldCor[3];
                    SetF(m, 0.0F, -FltMax, -FltMax);
                    SetB(m, input[3], input[4], input[5]);
                    Sub(m, serial);
                    break;
            }
        }

        private static void SetF(AlnMem m, float a, float ga, float gb) => m.F[0] = new State { A = a, Ga = ga, Gb = gb };

        private static void SetB(AlnMem m, float a, float ga, float gb) => m.B[0] = new State { A = a, Ga = ga, Gb = gb };

        private static float Mx(float a, float b) => a > b ? a : b;

        private static float Mx3(float a, float b, float c) => Mx(Mx(a, b), c);

        // ------------------------------------------------------------------ aln_seqseq.c

        private static void SeqSeqForward(AlnMem m)
        {
            var s = m.F;
            var seq1 = m.Seq1!;
            var seq2 = m.Seq2!;
            int starta = m.StartA, enda = m.EndA, startb = m.StartB, endb = m.EndB;
            float pa, pga, pgb, ca, xa, xga;
            float gpo = m.Ap.Gpo, gpe = m.Ap.Gpe, tgpe = m.Ap.Tgpe;
            var subm = m.Ap.Subm;
            int j;
            s[startb] = s[0];
            if (startb != 0)
            {
                for (j = startb + 1; j < endb; j++)
                {
                    s[j].A = -FltMax;
                    s[j].Ga = Mx(s[j - 1].Ga - gpe, s[j - 1].A - gpo);
                    s[j].Gb = -FltMax;
                }
            }
            else
            {
                for (j = startb + 1; j < endb; j++)
                {
                    s[j].A = -FltMax;
                    s[j].Ga = Mx(s[j - 1].Ga, s[j - 1].A) - tgpe;
                    s[j].Gb = -FltMax;
                }
            }
            s[endb].A = -FltMax;
            s[endb].Ga = -FltMax;
            s[endb].Gb = -FltMax;
            for (int i = starta; i < enda; i++)
            {
                var subp = subm[seq1[i]];
                pa = s[startb].A;
                pga = s[startb].Ga;
                pgb = s[startb].Gb;
                s[startb].A = -FltMax;
                s[startb].Ga = -FltMax;
                xa = s[startb].A;
                xga = s[startb].Ga;
                if (startb != 0)
                    s[startb].Gb = Mx(pgb - gpe, pa - gpo);
                else
                    s[startb].Gb = Mx(pgb, pa) - tgpe;
                for (j = startb + 1; j < endb; j++)
                {
                    ca = s[j].A;
                    pa = Mx3(pa, pga - gpo, pgb - gpo);
                    pa += subp[seq2[j - 1]];
                    s[j].A = pa;
                    pga = s[j].Ga;
                    s[j].Ga = Mx(xga - gpe, xa - gpo);
                    pgb = s[j].Gb;
                    s[j].Gb = Mx(pgb - gpe, ca - gpo);
                    pa = ca;
                    xa = s[j].A;
                    xga = s[j].Ga;
                }
                ca = s[j].A;
                pa = Mx3(pa, pga - gpo, pgb - gpo);
                pa += subp[seq2[j - 1]];
                s[j].A = pa;
                s[j].Ga = -FltMax;
                if (endb != m.LenB)
                    s[j].Gb = Mx(s[j].Gb - gpe, ca - gpo);
                else
                    s[j].Gb = Mx(s[j].Gb, ca) - tgpe;
            }
        }

        private static void SeqSeqBackward(AlnMem m)
        {
            var s = m.B;
            var seq1 = m.Seq1!;
            var seq2 = m.Seq2!;
            int starta = m.StartA2, enda = m.EndA2, startb = m.StartB, endb = m.EndB;
            float pa, pga, pgb, ca, xa, xga;
            float gpo = m.Ap.Gpo, gpe = m.Ap.Gpe, tgpe = m.Ap.Tgpe;
            var subm = m.Ap.Subm;
            int j;
            s[endb] = s[0];
            if (endb != m.LenB)
            {
                for (j = endb - 1; j > startb; j--)
                {
                    s[j].A = -FltMax;
                    s[j].Ga = Mx(s[j + 1].Ga - gpe, s[j + 1].A - gpo);
                    s[j].Gb = -FltMax;
                }
            }
            else
            {
                for (j = endb - 1; j > startb; j--)
                {
                    s[j].A = -FltMax;
                    s[j].Ga = Mx(s[j + 1].Ga, s[j + 1].A) - tgpe;
                    s[j].Gb = -FltMax;
                }
            }
            s[startb].A = -FltMax;
            s[startb].Ga = -FltMax;
            s[startb].Gb = -FltMax;
            int i = enda - starta;
            while (i-- > 0)
            {
                var subp = subm[seq1[starta + i]];
                pa = s[endb].A;
                pga = s[endb].Ga;
                pgb = s[endb].Gb;
                s[endb].A = -FltMax;
                s[endb].Ga = -FltMax;
                xa = s[endb].A;
                xga = s[endb].Ga;
                if (endb != m.LenB)
                    s[endb].Gb = Mx(pgb - gpe, pa - gpo);
                else
                    s[endb].Gb = Mx(pgb, pa) - tgpe;
                for (j = endb - 1; j > startb; j--)
                {
                    ca = s[j].A;
                    pa = Mx3(pa, pga - gpo, pgb - gpo);
                    pa += subp[seq2[j]];
                    s[j].A = pa;
                    pga = s[j].Ga;
                    s[j].Ga = Mx(xga - gpe, xa - gpo);
                    pgb = s[j].Gb;
                    s[j].Gb = Mx(pgb - gpe, ca - gpo);
                    pa = ca;
                    xa = s[j].A;
                    xga = s[j].Ga;
                }
                ca = s[j].A;
                pa = Mx3(pa, pga - gpo, pgb - gpo);
                pa += subp[seq2[j]];
                s[j].A = pa;
                s[j].Ga = -FltMax;
                if (startb != 0)
                    s[j].Gb = Mx(s[j].Gb - gpe, ca - gpo);
                else
                    s[j].Gb = Mx(s[j].Gb, ca) - tgpe;
            }
        }

        private static void SeqSeqMeetup(AlnMem m, Span<int> oldCor, out int meet, out int t)
        {
            var f = m.F;
            var b = m.B;
            float gpo = m.Ap.Gpo, gpe = m.Ap.Gpe, tgpe = m.Ap.Tgpe;
            int transition = -1;
            float max = -FltMax;
            float middle = (float)(oldCor[3] - oldCor[2]) / 2.0F + (float)oldCor[2];
            float sub;
            int c = -1;
            int i;
            for (i = oldCor[2]; i < oldCor[3]; i++)
            {
                sub = MathF.Abs(middle - (float)i);
                sub /= 1000.0F;
                if (f[i].A + b[i].A - sub > max) { max = f[i].A + b[i].A - sub; transition = 1; c = i; }
                if (f[i].A + b[i].Ga - gpo - sub > max) { max = f[i].A + b[i].Ga - gpo - sub; transition = 2; c = i; }
                if (f[i].A + b[i].Gb - gpo - sub > max) { max = f[i].A + b[i].Gb - gpo - sub; transition = 3; c = i; }
                if (f[i].Ga + b[i].A - gpo - sub > max) { max = f[i].Ga + b[i].A - gpo - sub; transition = 5; c = i; }
                if (m.StartB == 0)
                {
                    if (f[i].Gb + b[i].Gb - tgpe - sub > max) { max = f[i].Gb + b[i].Gb - tgpe - sub; transition = 6; c = i; }
                }
                else
                {
                    if (f[i].Gb + b[i].Gb - gpe - sub > max) { max = f[i].Gb + b[i].Gb - gpe - sub; transition = 6; c = i; }
                }
                if (f[i].Gb + b[i].A - gpo - sub > max) { max = f[i].Gb + b[i].A - gpo - sub; transition = 7; c = i; }
            }
            i = oldCor[3];
            sub = MathF.Abs(middle - (float)i);
            sub /= 1000.0F;
            if (f[i].A + b[i].Gb - gpo - sub > max) { max = f[i].A + b[i].Gb - gpo - sub; transition = 3; c = i; }
            if (m.EndB == m.LenB)
            {
                if (f[i].Gb + b[i].Gb - tgpe - sub > max) { transition = 6; c = i; }
            }
            else
            {
                if (f[i].Gb + b[i].Gb - gpe - sub > max) { transition = 6; c = i; }
            }
            meet = c;
            t = transition;
        }

        // ------------------------------------------------------------------ aln_seqprofile.c

        private static void SeqProfForward(AlnMem m)
        {
            var s = m.F;
            var prof1 = m.Prof1!;
            var seq2 = m.Seq2!;
            int sip = m.Sip;
            float pa, pga, pgb, ca, xa, xga;
            float open = m.Ap.Gpo * (float)sip;
            float ext = m.Ap.Gpe * (float)sip;
            float text = m.Ap.Tgpe * (float)sip;
            int p = m.StartA << 6;
            int j;
            s[m.StartB] = s[0];
            if (m.StartB != 0)
            {
                for (j = m.StartB + 1; j < m.EndB; j++)
                {
                    s[j].A = -FltMax;
                    s[j].Ga = Mx(s[j - 1].Ga - ext, s[j - 1].A - open);
                    s[j].Gb = -FltMax;
                }
            }
            else
            {
                for (j = m.StartB + 1; j < m.EndB; j++)
                {
                    s[j].A = -FltMax;
                    s[j].Ga = Mx(s[j - 1].Ga, s[j - 1].A) - text;
                    s[j].Gb = -FltMax;
                }
            }
            s[m.EndB].A = -FltMax;
            s[m.EndB].Ga = -FltMax;
            s[m.EndB].Gb = -FltMax;
            for (int i = m.StartA; i < m.EndA; i++)
            {
                p += 64;
                pa = s[m.StartB].A;
                pga = s[m.StartB].Ga;
                pgb = s[m.StartB].Gb;
                s[m.StartB].A = -FltMax;
                s[m.StartB].Ga = -FltMax;
                xa = s[m.StartB].A;
                xga = s[m.StartB].Ga;
                if (m.StartB != 0)
                    s[m.StartB].Gb = Mx(pgb + prof1[p + 28], pa + prof1[p + 27]);
                else
                    s[m.StartB].Gb = Mx(pgb, pa) + prof1[p + 29];
                for (j = m.StartB + 1; j < m.EndB; j++)
                {
                    ca = s[j].A;
                    pa = Mx3(pa, pga - open, pgb + prof1[p - 37]);
                    pa += prof1[p + 32 + seq2[j - 1]];
                    s[j].A = pa;
                    pga = s[j].Ga;
                    s[j].Ga = Mx(xga - ext, xa - open);
                    pgb = s[j].Gb;
                    s[j].Gb = Mx(pgb + prof1[p + 28], ca + prof1[p + 27]);
                    pa = ca;
                    xa = s[j].A;
                    xga = s[j].Ga;
                }
                ca = s[j].A;
                pa = Mx3(pa, pga - open, pgb + prof1[p - 37]);
                pa += prof1[p + 32 + seq2[j - 1]];
                s[j].A = pa;
                s[j].Ga = -FltMax;
                if (m.EndB != m.LenB)
                    s[j].Gb = Mx(s[j].Gb + prof1[p + 28], ca + prof1[p + 27]);
                else
                    s[j].Gb = Mx(s[j].Gb, ca) + prof1[p + 29];
            }
        }

        private static void SeqProfBackward(AlnMem m)
        {
            var s = m.B;
            var prof1 = m.Prof1!;
            var seq2 = m.Seq2!;
            int sip = m.Sip;
            float pa, pga, pgb, ca, xa, xga;
            float open = m.Ap.Gpo * (float)sip;
            float ext = m.Ap.Gpe * (float)sip;
            float text = m.Ap.Tgpe * (float)sip;
            int p = (m.EndA2 + 1) << 6;
            int j;
            s[m.EndB] = s[0];
            if (m.EndB != m.LenB)
            {
                for (j = m.EndB - 1; j > m.StartB; j--)
                {
                    s[j].A = -FltMax;
                    s[j].Ga = Mx(s[j + 1].Ga - ext, s[j + 1].A - open);
                    s[j].Gb = -FltMax;
                }
            }
            else
            {
                for (j = m.EndB - 1; j > m.StartB; j--)
                {
                    s[j].A = -FltMax;
                    s[j].Ga = Mx(s[j + 1].Ga, s[j + 1].A) - text;
                    s[j].Gb = -FltMax;
                }
            }
            s[m.StartB].A = -FltMax;
            s[m.StartB].Ga = -FltMax;
            s[m.StartB].Gb = -FltMax;
            int i = m.EndA2 - m.StartA2;
            while (i-- > 0)
            {
                p -= 64;
                pa = s[m.EndB].A;
                pga = s[m.EndB].Ga;
                pgb = s[m.EndB].Gb;
                s[m.EndB].A = -FltMax;
                s[m.EndB].Ga = -FltMax;
                xa = s[m.EndB].A;
                xga = s[m.EndB].Ga;
                if (m.EndB != m.LenB)
                    s[m.EndB].Gb = Mx(pgb + prof1[p + 28], pa + prof1[p + 27]);
                else
                    s[m.EndB].Gb = Mx(pgb, pa) + prof1[p + 29];
                for (j = m.EndB - 1; j > m.StartB; j--)
                {
                    ca = s[j].A;
                    pa = Mx3(pa, pga - open, pgb + prof1[p + 91]);
                    pa += prof1[p + 32 + seq2[j]];
                    s[j].A = pa;
                    pga = s[j].Ga;
                    s[j].Ga = Mx(xga - ext, xa - open);
                    pgb = s[j].Gb;
                    s[j].Gb = Mx(pgb + prof1[p + 28], ca + prof1[p + 27]);
                    pa = ca;
                    xa = s[j].A;
                    xga = s[j].Ga;
                }
                ca = s[j].A;
                pa = Mx3(pa, pga - open, pgb + prof1[p + 91]);
                pa += prof1[p + 32 + seq2[j]];
                s[j].A = pa;
                s[j].Ga = -FltMax;
                if (m.StartB != 0)
                    s[j].Gb = Mx(s[j].Gb + prof1[p + 28], ca + prof1[p + 27]);
                else
                    s[j].Gb = Mx(s[j].Gb, ca) + prof1[p + 29];
            }
        }

        private static void SeqProfMeetup(AlnMem m, Span<int> oldCor, out int meet, out int t)
        {
            var f = m.F;
            var b = m.B;
            var prof1 = m.Prof1!;
            int sip = m.Sip;
            int transition = -1;
            float open = m.Ap.Gpo * (float)sip;
            float max = -FltMax;
            float middle = (float)(oldCor[3] - oldCor[2]) / 2.0F + (float)oldCor[2];
            float sub;
            int p = (oldCor[4] + 1) << 6;
            int c = -1;
            int i;
            for (i = oldCor[2]; i < oldCor[3]; i++)
            {
                sub = MathF.Abs(middle - (float)i);
                sub /= 1000.0F;
                if (f[i].A + b[i].A - sub > max) { max = f[i].A + b[i].A - sub; transition = 1; c = i; }
                if (f[i].A + b[i].Ga - open - sub > max) { max = f[i].A + b[i].Ga - open - sub; transition = 2; c = i; }
                if (f[i].A + b[i].Gb + prof1[p + 27] - sub > max) { max = f[i].A + b[i].Gb + prof1[p + 27] - sub; transition = 3; c = i; }
                if (f[i].Ga + b[i].A - open - sub > max) { max = f[i].Ga + b[i].A - open - sub; transition = 5; c = i; }
                if (m.StartB == 0)
                {
                    if (f[i].Gb + b[i].Gb + prof1[p + 29] - sub > max) { max = f[i].Gb + b[i].Gb + prof1[p + 29] - sub; transition = 6; c = i; }
                }
                else
                {
                    if (f[i].Gb + b[i].Gb + prof1[p + 28] - sub > max) { max = f[i].Gb + b[i].Gb + prof1[p + 28] - sub; transition = 6; c = i; }
                }
                if (f[i].Gb + b[i].A + prof1[p - 37] - sub > max) { max = f[i].Gb + b[i].A + prof1[p - 37] - sub; transition = 7; c = i; }
            }
            i = oldCor[3];
            sub = MathF.Abs(middle - (float)i);
            sub /= 1000.0F;
            if (f[i].A + b[i].Gb + prof1[p + 27] - sub > max) { max = f[i].A + b[i].Gb + prof1[p + 27] - sub; transition = 3; c = i; }
            if (m.EndB == m.LenB)
            {
                if (f[i].Gb + b[i].Gb + prof1[p + 29] - sub > max) { transition = 6; c = i; }
            }
            else
            {
                if (f[i].Gb + b[i].Gb + prof1[p + 28] - sub > max) { transition = 6; c = i; }
            }
            meet = c;
            t = transition;
        }

        // ------------------------------------------------------------------ aln_profileprofile.c

        private static float Dot(float[] prof1, int p1, float[] prof2, int p2, int[] freq, int f, float pa)
        {
            for (int c = f; c >= 0; c--)
                pa += prof1[p1 + freq[c]] * prof2[p2 + freq[c]];
            return pa;
        }

        private static int Frequencies(float[] prof1, int p1, int[] freq)
        {
            int f = 0;
            for (int j = 0; j < 23; j++)
            {
                if (prof1[p1 + j] != 0)
                {
                    freq[f] = j;
                    f++;
                }
            }
            return f - 1;
        }

        private static void ProfProfForward(AlnMem m)
        {
            var freq = new int[24];
            var prof1 = m.Prof1!;
            var prof2 = m.Prof2!;
            var s = m.F;
            float pa, pga, pgb, ca, xa, xga;
            int p1 = m.StartA << 6;
            int p2 = m.StartB << 6;
            int j;
            s[m.StartB] = s[0];
            if (m.StartB != 0)
            {
                for (j = m.StartB + 1; j < m.EndB; j++)
                {
                    p2 += 64;
                    s[j].A = -FltMax;
                    s[j].Ga = Mx(s[j - 1].Ga + prof2[p2 + 28], s[j - 1].A + prof2[p2 + 27]);
                    s[j].Gb = -FltMax;
                }
                p2 += 64;
            }
            else
            {
                for (j = m.StartB + 1; j < m.EndB; j++)
                {
                    p2 += 64;
                    s[j].A = -FltMax;
                    s[j].Ga = Mx(s[j - 1].Ga, s[j - 1].A) + prof2[p2 + 29];
                    s[j].Gb = -FltMax;
                }
                p2 += 64;
            }
            p2 -= (m.EndB - m.StartB) << 6;
            s[m.EndB].A = -FltMax;
            s[m.EndB].Ga = -FltMax;
            s[m.EndB].Gb = -FltMax;
            for (int i = m.StartA; i < m.EndA; i++)
            {
                p1 += 64;
                int f = Frequencies(prof1, p1, freq);
                pa = s[m.StartB].A;
                pga = s[m.StartB].Ga;
                pgb = s[m.StartB].Gb;
                s[m.StartB].A = -FltMax;
                s[m.StartB].Ga = -FltMax;
                xa = s[m.StartB].A;
                xga = s[m.StartB].Ga;
                if (m.StartB != 0)
                    s[m.StartB].Gb = Mx(pgb + prof1[p1 + 28], pa + prof1[p1 + 27]);
                else
                    s[m.StartB].Gb = Mx(pgb, pa) + prof1[p1 + 29];
                for (j = m.StartB + 1; j < m.EndB; j++)
                {
                    p2 += 64;
                    ca = s[j].A;
                    pa = Mx3(pa, pga + prof2[p2 - 37], pgb + prof1[p1 - 37]);
                    p2 += 32;
                    pa = Dot(prof1, p1, prof2, p2, freq, f, pa);
                    p2 -= 32;
                    s[j].A = pa;
                    pga = s[j].Ga;
                    s[j].Ga = Mx(xga + prof2[p2 + 28], xa + prof2[p2 + 27]);
                    pgb = s[j].Gb;
                    s[j].Gb = Mx(pgb + prof1[p1 + 28], ca + prof1[p1 + 27]);
                    pa = ca;
                    xa = s[j].A;
                    xga = s[j].Ga;
                }
                p2 += 64;
                ca = s[j].A;
                pa = Mx3(pa, pga + prof2[p2 - 37], pgb + prof1[p1 - 37]);
                p2 += 32;
                pa = Dot(prof1, p1, prof2, p2, freq, f, pa);
                p2 -= 32;
                s[j].A = pa;
                s[j].Ga = -FltMax;
                if (m.EndB != m.LenB)
                    s[j].Gb = Mx(s[j].Gb + prof1[p1 + 28], ca + prof1[p1 + 27]);
                else
                    s[j].Gb = Mx(s[j].Gb, ca) + prof1[p1 + 29];
                p2 -= (m.EndB - m.StartB) << 6;
            }
        }

        private static void ProfProfBackward(AlnMem m)
        {
            var freq = new int[24];
            var s = m.B;
            var prof1 = m.Prof1!;
            var prof2 = m.Prof2!;
            float pa, pga, pgb, ca, xa, xga;
            int p1 = (m.EndA2 + 1) << 6;
            int p2 = (m.EndB + 1) << 6;
            int j;
            s[m.EndB] = s[0];
            if (m.EndB != m.LenB)
            {
                for (j = m.EndB - 1; j > m.StartB; j--)
                {
                    p2 -= 64;
                    s[j].A = -FltMax;
                    s[j].Ga = Mx(s[j + 1].Ga + prof2[p2 + 28], s[j + 1].A + prof2[p2 + 27]);
                    s[j].Gb = -FltMax;
                }
                p2 -= 64;
            }
            else
            {
                for (j = m.EndB - 1; j > m.StartB; j--)
                {
                    p2 -= 64;
                    s[j].A = -FltMax;
                    s[j].Ga = Mx(s[j + 1].Ga, s[j + 1].A) + prof2[p2 + 29];
                    s[j].Gb = -FltMax;
                }
                p2 -= 64;
            }
            s[m.StartB].A = -FltMax;
            s[m.StartB].Ga = -FltMax;
            s[m.StartB].Gb = -FltMax;
            int i = m.EndA2 - m.StartA2;
            while (i-- > 0)
            {
                p1 -= 64;
                int f = Frequencies(prof1, p1, freq);
                pa = s[m.EndB].A;
                pga = s[m.EndB].Ga;
                pgb = s[m.EndB].Gb;
                s[m.EndB].A = -FltMax;
                s[m.EndB].Ga = -FltMax;
                xa = s[m.EndB].A;
                xga = s[m.EndB].Ga;
                if (m.EndB != m.LenB)
                    s[m.EndB].Gb = Mx(pgb + prof1[p1 + 28], pa + prof1[p1 + 27]);
                else
                    s[m.EndB].Gb = Mx(pgb, pa) + prof1[p1 + 29];
                p2 += (m.EndB - m.StartB) << 6;
                for (j = m.EndB - 1; j > m.StartB; j--)
                {
                    p2 -= 64;
                    ca = s[j].A;
                    pa = Mx3(pa, pga + prof2[p2 + 91], pgb + prof1[p1 + 91]);
                    p2 += 32;
                    pa = Dot(prof1, p1, prof2, p2, freq, f, pa);
                    p2 -= 32;
                    s[j].A = pa;
                    pga = s[j].Ga;
                    s[j].Ga = Mx(xga + prof2[p2 + 28], xa + prof2[p2 + 27]);
                    pgb = s[j].Gb;
                    s[j].Gb = Mx(pgb + prof1[p1 + 28], ca + prof1[p1 + 27]);
                    pa = ca;
                    xa = s[j].A;
                    xga = s[j].Ga;
                }
                p2 -= 64;
                ca = s[j].A;
                pa = Mx3(pa, pga + prof2[p2 + 91], pgb + prof1[p1 + 91]);
                p2 += 32;
                pa = Dot(prof1, p1, prof2, p2, freq, f, pa);
                p2 -= 32;
                s[j].A = pa;
                s[j].Ga = -FltMax;
                if (m.StartB != 0)
                    s[j].Gb = Mx(s[j].Gb + prof1[p1 + 28], ca + prof1[p1 + 27]);
                else
                    s[j].Gb = Mx(s[j].Gb, ca) + prof1[p1 + 29];
            }
        }

        private static void ProfProfMeetup(AlnMem m, Span<int> oldCor, out int meet, out int t)
        {
            var f = m.F;
            var b = m.B;
            var prof1 = m.Prof1!;
            var prof2 = m.Prof2!;
            int transition = -1;
            float max = -FltMax;
            float middle = (float)(oldCor[3] - oldCor[2]) / 2.0F + (float)oldCor[2];
            float sub;
            int p1 = (oldCor[4] + 1) << 6;
            int p2 = oldCor[2] << 6;
            int c = -1;
            int i;
            for (i = oldCor[2]; i < oldCor[3]; i++)
            {
                sub = MathF.Abs(middle - (float)i);
                sub /= 1000.0F;
                p2 += 64;
                if (f[i].A + b[i].A - sub > max) { max = f[i].A + b[i].A - sub; transition = 1; c = i; }
                if (f[i].A + b[i].Ga + prof2[p2 + 27] - sub > max) { max = f[i].A + b[i].Ga + prof2[p2 + 27] - sub; transition = 2; c = i; }
                if (f[i].A + b[i].Gb + prof1[p1 + 27] - sub > max) { max = f[i].A + b[i].Gb + prof1[p1 + 27] - sub; transition = 3; c = i; }
                if (f[i].Ga + b[i].A + prof2[p2 - 37] - sub > max) { max = f[i].Ga + b[i].A + prof2[p2 - 37] - sub; transition = 5; c = i; }
                if (m.StartB == 0)
                {
                    if (f[i].Gb + b[i].Gb + prof1[p1 + 29] - sub > max) { max = f[i].Gb + b[i].Gb + prof1[p1 + 29] - sub; transition = 6; c = i; }
                }
                else
                {
                    if (f[i].Gb + b[i].Gb + prof1[p1 + 28] - sub > max) { max = f[i].Gb + b[i].Gb + prof1[p1 + 28] - sub; transition = 6; c = i; }
                }
                if (f[i].Gb + b[i].A + prof1[p1 - 37] - sub > max) { max = f[i].Gb + b[i].A + prof1[p1 - 37] - sub; transition = 7; c = i; }
            }
            i = oldCor[3];
            sub = MathF.Abs(middle - (float)i);
            sub /= 1000.0F;
            if (f[i].A + b[i].Gb + prof1[p1 + 27] - sub > max) { max = f[i].A + b[i].Gb + prof1[p1 + 27] - sub; transition = 3; c = i; }
            if (m.EndB == m.LenB)
            {
                if (f[i].Gb + b[i].Gb + prof1[p1 + 29] - sub > max) { transition = 6; c = i; }
            }
            else
            {
                if (f[i].Gb + b[i].Gb + prof1[p1 + 28] - sub > max) { transition = 6; c = i; }
            }
            meet = c;
            t = transition;
        }
    }

    #endregion
}
