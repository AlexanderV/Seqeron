using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The scoring matrix file of POA 2.0 (<c>read_score_matrix</c>, <c>seq_util.c</c>): the symbols of the matrix, the integer scores and the gap penalty model
/// (<c>GAP-PENALTIES=open ext1 ext2</c>, <c>GAP-PENALTIES-X=...</c>, <c>GAP-TRUNCATION-LENGTH=</c>, <c>GAP-DECAY-LENGTH=</c>; defaults 12 / 2 / 0, truncation 16, decay 0).
/// </summary>
public sealed class PoaScoreMatrix
{
    /// <summary>The matrix symbols in file order (<c>symbol[]</c>).</summary>
    public string Symbols { get; }

    internal int[,] Score { get; }

    internal int[,] GapPenaltySet { get; }

    internal int TruncGapLength { get; }

    internal int DecayGapLength { get; }

    internal int MaxGapLength => TruncGapLength + DecayGapLength;

    internal int[] GapPenaltyX { get; }

    internal int[] GapPenaltyY { get; }

    private PoaScoreMatrix(string symbols, int[,] score, int[,] gaps, int trunc, int decay)
    {
        Symbols = symbols;
        Score = score;
        GapPenaltySet = gaps;
        TruncGapLength = trunc;
        DecayGapLength = decay;
        int max = trunc + decay;
        GapPenaltyX = new int[max + 2];
        GapPenaltyY = new int[max + 2];
        GapPenaltyX[0] = gaps[0, 0];
        GapPenaltyY[0] = gaps[1, 0];
        for (int i = 1; i < trunc; i++)
        {
            GapPenaltyX[i] = gaps[0, 1];
            GapPenaltyY[i] = gaps[1, 1];
        }

        for (int i = 0; i < decay; i++)
        {
            double decX = (gaps[0, 1] - gaps[0, 2]) / (double)(decay + 1);
            double decY = (gaps[1, 1] - gaps[1, 2]) / (double)(decay + 1);
            GapPenaltyX[i + trunc] = (int)(gaps[0, 1] - (i + 1) * decX);
            GapPenaltyY[i + trunc] = (int)(gaps[1, 1] - (i + 1) * decY);
        }

        GapPenaltyX[max] = gaps[0, 2];
        GapPenaltyY[max] = gaps[1, 2];
        GapPenaltyX[max + 1] = 0;
        GapPenaltyY[max + 1] = 0;
    }

    // C sscanf of one %d at s[pos..]: leading blanks, sign, digits.
    private static bool ScanInt(string s, ref int pos, out int value)
    {
        value = 0;
        int p = pos;
        while (p < s.Length && char.IsWhiteSpace(s[p])) p++;
        int start = p;
        if (p < s.Length && (s[p] == '+' || s[p] == '-')) p++;
        int digits = p;
        while (p < s.Length && s[p] >= '0' && s[p] <= '9') p++;
        if (p == digits) return false;
        value = (int)Math.Clamp(long.Parse(s.AsSpan(start, p - start), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture), int.MinValue, int.MaxValue);
        pos = p;
        return true;
    }

    private static bool ScanAfter(string line, string prefix, int count, int[] values)
    {
        if (!line.StartsWith(prefix, StringComparison.Ordinal)) return false;
        int pos = prefix.Length;
        for (int k = 0; k < count; k++)
        {
            if (!ScanInt(line, ref pos, out values[k])) return false;
        }

        return true;
    }

    /// <summary>Reads a matrix file (the text of the file).</summary>
    /// <exception cref="ArgumentException">A matrix line starts with a symbol that is not in the symbol line, or a score is missing (the program: "Error reading matrix file").</exception>
    public static PoaScoreMatrix Parse(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        var gaps = new int[2, 3] { { 12, 2, 0 }, { 12, 2, 0 } };
        int trunc = 16, decay = 0;
        var symbols = new StringBuilder();
        var score = new int[128, 128];
        bool foundSymbolLine = false;
        var v = new int[3];
        int pos = 0;
        while (pos < text.Length)
        {
            int nl = text.IndexOf('\n', pos);
            string line = nl < 0 ? text[pos..] : text[pos..(nl + 1)]; // fgets keeps the newline
            pos = nl < 0 ? text.Length : nl + 1;
            if (line.Length > 1023) line = line[..1023];
            if (line[0] == '#' || line[0] == '\n') continue;
            if (ScanAfter(line, "GAP-TRUNCATION-LENGTH=", 1, v))
            {
                trunc = v[0];
            }
            else if (ScanAfter(line, "GAP-DECAY-LENGTH=", 1, v))
            {
                decay = v[0];
            }
            else if (ScanAfter(line, "GAP-PENALTIES=", 3, v))
            {
                for (int k = 0; k < 3; k++) gaps[0, k] = gaps[1, k] = v[k];
            }
            else if (ScanAfter(line, "GAP-PENALTIES-X=", 3, v))
            {
                for (int k = 0; k < 3; k++) gaps[1, k] = v[k];
            }
            else if (!foundSymbolLine)
            {
                foreach (char c in line)
                    if (!char.IsWhiteSpace(c)) symbols.Append(c);
                // the score matrix is indexed by char (int score[128][128] of the 7-bit alphabet): a symbol above 127 or more than 128 symbols have no meaning in the program
                for (int k = 0; k < symbols.Length; k++)
                    if (symbols[k] > 127)
                        throw new ArgumentException($"The matrix symbol U+{(int)symbols[k]:X4} is not a 7-bit character: the score matrix of POA is indexed by char.");
                if (symbols.Length > 128)
                    throw new ArgumentException($"The matrix has {symbols.Length} symbols: the score matrix of POA is indexed by char (at most 128).");
                foundSymbolLine = true;
            }
            else
            {
                foundSymbolLine = false;
                int nsymb = symbols.Length;
                for (int isymb = nsymb - 1; isymb >= 0; isymb--)
                {
                    if (symbols[isymb] != line[0]) continue;
                    foundSymbolLine = true;
                    int j = 1;
                    for (int i = 0; i < nsymb; i++)
                    {
                        if (!ScanInt(line, ref j, out score[isymb, i]))
                            throw new ArgumentException($"Missing score value for pair {symbols[isymb]}:{symbols[i]}");
                    }

                    break;
                }

                if (!foundSymbolLine)
                    throw new ArgumentException($"Missing or unknown sequence symbol: {line[0]}");
            }
        }

        if (symbols.Length == 0) throw new ArgumentException("The matrix file has no symbol line.");
        return new PoaScoreMatrix(symbols.ToString(), score, gaps, trunc, decay);
    }

    /// <summary><c>print_score_matrix</c>: the rows and columns of <paramref name="subset"/> (the <c>-printmatrix</c> text).</summary>
    public string Print(string subset)
    {
        var sb = new StringBuilder(" ");
        foreach (char c in subset) sb.Append("  ").Append(c);
        sb.Append('\n');
        foreach (char ci in subset)
        {
            int im = Symbols.Length - 1;
            while (im >= 0 && Symbols[im] != ci) im--;
            sb.Append(ci);
            foreach (char cj in subset)
            {
                int jm = Symbols.Length - 1;
                while (jm >= 0 && Symbols[jm] != cj) jm--;
                sb.Append(Score[Math.Max(im, 0), Math.Max(jm, 0)].ToString(CultureInfo.InvariantCulture).PadLeft(3));
            }

            sb.Append('\n');
        }

        return sb.ToString();
    }
}

/// <summary>One input sequence of POA (a FASTA record): the name (the first word of the header, at most 31 characters), the title (the rest of the header) and the residues.</summary>
/// <param name="Name">Name.</param>
/// <param name="Title">Title; null or empty = <c>untitled</c>.</param>
/// <param name="Residues">Residues (white space is removed).</param>
public sealed record PoaSequence(string Name, string? Title, string Residues);

/// <summary>Letter case handling of the input (<c>-tolower</c> / <c>-toupper</c>).</summary>
public enum PoaCase
{
    /// <summary>Keep the case (the program's default).</summary>
    Keep = 0,

    /// <summary><c>-tolower</c>.</summary>
    Lower = 1,

    /// <summary><c>-toupper</c>.</summary>
    Upper = 2,
}

/// <summary>Options of POA 2.0 (<c>poa -do_global -do_progressive -read_pairscores F -fuse_all -hb -hbmin X -best -preserve_seqorder -v MATRIX</c>).</summary>
public sealed record PoaOptions
{
    /// <summary><c>-do_global</c>: global alignment (the default is local alignment with free end gaps).</summary>
    public bool Global { get; init; }

    /// <summary><c>-do_progressive</c>: progressive alignment along the order of the pairwise alignment scores (the default is iterative: every sequence is aligned to the growing partial order in input order).</summary>
    public bool Progressive { get; init; }

    /// <summary>
    /// <c>-read_pairscores</c>: the text of a tab / blank delimited file of <c>name1 name2 score</c> lines (merging along these scores). A name that is no sequence makes the program give up merging:
    /// it reports "Error generating pair scores", exits with status 0 and writes the first sequence alone.
    /// </summary>
    public string? PairScores { get; init; }

    /// <summary><c>-fuse_all</c>: fuse identical letters on align rings.</summary>
    public bool FuseAll { get; init; }

    /// <summary><c>-hb</c>: heaviest bundling (consensus sequences <c>CONSENS0</c>, ... are added to the alignment).</summary>
    public bool HeaviestBundle { get; init; }

    /// <summary><c>-hbmin</c>: the fraction of the positions of a sequence that must lie on a bundle's path (default 0.9; the program reads the value with <c>atof</c> into a float).</summary>
    public float HeaviestBundleMinimum { get; init; } = 0.9f;

    /// <summary><c>-best</c>: restrict the PIR output to the sequences of the heaviest bundle (bundle 0; every alignment column is still written).</summary>
    public bool Best { get; init; }

    /// <summary><c>-preserve_seqorder</c>: write the sequences in their input order (progressive merging reorders them).</summary>
    public bool PreserveSequenceOrder { get; init; }

    /// <summary><c>-tolower</c> / <c>-toupper</c>.</summary>
    public PoaCase Case { get; init; }

    /// <summary><c>-v</c>: write the gap penalty tables to standard output.</summary>
    public bool Verbose { get; init; }

    /// <summary><c>-printmatrix LETTERS</c>: write the rows and columns of the score matrix for these symbols to standard output.</summary>
    public string? PrintMatrix { get; init; }

    /// <summary>The <c>#</c> comment line of the FASTA file (without the <c>#</c>): the title of the partial order in the PO output; null = the title of the first sequence.</summary>
    public string? Comment { get; init; }

    /// <summary>Cancellation of the call.</summary>
    public CancellationToken CancellationToken { get; init; }
}

/// <summary>The result of <see cref="SequenceAligner.PoaRun"/>.</summary>
/// <param name="Names">Name of every row (the sequences in the program's output order, then the <c>CONSENS</c> sequences of <c>-hb</c>).</param>
/// <param name="Titles">Title of every row.</param>
/// <param name="Rows">The multiple alignment (gap character <c>-</c>), one row per name; letters are the symbols of the matrix.</param>
/// <param name="Order">For every row of the first <c>InputCount</c>, the input index of its sequence; consensus rows have -1.</param>
/// <param name="Po">The partial-order alignment in the program's PO format (<c>-po</c>).</param>
/// <param name="Pir">The alignment in the program's FASTA-PIR format (<c>-pir</c>, gap <c>.</c>, 60 columns; restricted to bundle 0 with <see cref="PoaOptions.Best"/>).</param>
/// <param name="Clustal">The alignment in the program's CLUSTAL format (<c>-clustal</c>).</param>
/// <param name="Stdout">What the program writes to standard output (<c>-v</c> gap tables, the permutation line of <c>-preserve_seqorder</c>).</param>
public sealed record PoaRunResult(string[] Names, string[] Titles, string[] Rows, int[] Order, string Po, string Pir, string Clustal, string Stdout);

public static partial class SequenceAligner
{
    #region POA 2.0

    // Port of POA V2 (Lee, Grasso & Sharlow 2002, Bioinformatics 18:452; Lee 2003, Bioinformatics 19:999; Debian poa 2.0+20060928: main.c, lpo.c, buildup_lpo.c, align_lpo_po2.c, heaviest_bundle.c,
    // lpo_format.c, seq_util.c, fasta_format.c, create_seq.c). The partial-order graph keeps the program's data structure (letters with ordered left / right links, an ordered source list, a
    // circular alignment ring), because the order of links and ring members is visible in the PO output.

    private const int PoaNoBundle = -1;

    internal sealed class PoaLetter
    {
        public readonly List<int> Left = new();
        public readonly List<int> Right = new();
        public readonly List<(int Iseq, int Ipos)> Source = new();
        public int AlignRing, RingId, Letter;
    }

    internal sealed class PoaSourceInfo
    {
        public string Name = string.Empty, Title = string.Empty;
        public int Length, Istart, Weight, BundleId = PoaNoBundle;
        public int Input = -1; // index of the sequence in the caller's list (-1: a consensus)
    }

    internal sealed class PoaLpo
    {
        public string Name = string.Empty, Title = string.Empty;
        public PoaLetter[] Letters = Array.Empty<PoaLetter>();
        public readonly List<PoaSourceInfo> Sources = new();
        public string Sequence = string.Empty; // raw residues until lpo_init
        public int Length => Letters.Length;
    }

    /// <summary>
    /// POA 2.0's multiple sequence alignment of <paramref name="sequences"/> (<c>poa -read_fasta in.fa [-do_global] [-do_progressive] ... matrix</c>): partial-order alignment
    /// (Lee et al. 2002) of every sequence to the growing partial-order graph by dynamic programming with the program's affine gap model (<see cref="PoaScoreMatrix"/>), fusion of the aligned letters
    /// into align rings, optional heaviest bundling. The default matrix is BLOSUM80 with gap penalties 12 / 6 / 6 (<c>blosum80.mat</c> of the distribution).
    /// </summary>
    /// <returns>Rows in input order (sequence residues are mapped to the symbols of the matrix; letters outside the matrix become its first symbol, as in the program).</returns>
    public static MultipleAlignmentResult MultipleAlignPoa(IEnumerable<string> sequences, PoaOptions? options = null)
    {
        var list = new List<PoaSequence>();
        int n = 1;
        foreach (var s in sequences ?? throw new ArgumentNullException(nameof(sequences)))
        {
            if (s == null || s.All(c => c == ' ' || (c >= '\t' && c <= '\r'))) throw new ArgumentException("A sequence is null or empty.", nameof(sequences));
            list.Add(new PoaSequence("S" + (n++).ToString(CultureInfo.InvariantCulture), null, s));
        }
        var run = PoaRun(list, null, options);
        var rows = new string[list.Count];
        for (int i = 0; i < run.Rows.Length; i++)
            if (run.Order[i] >= 0) rows[run.Order[i]] = run.Rows[i];
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
            return new MultipleAlignmentResult(rows, consensus, 0) { ExactTotalScore = double.NaN };
        }
    }

    /// <summary>
    /// POA 2.0 with every output of the program (<see cref="PoaRunResult"/>): the oracle is the installed Debian <c>poa</c> 2.0+20060928: rows, PO / PIR / CLUSTAL files and standard output equal on
    /// every committed run (scripts/alignment/generate_poa_oracle.py).
    /// </summary>
    /// <param name="sequences">The input records (at least one).</param>
    /// <param name="matrix">The scoring matrix (null = <see cref="PoaScoreMatrix"/> BLOSUM80 of the distribution).</param>
    /// <param name="options">The program's options.</param>
    /// <exception cref="ArgumentException">No sequence.</exception>
    public static PoaRunResult PoaRun(IReadOnlyList<PoaSequence> sequences, PoaScoreMatrix? matrix = null, PoaOptions? options = null)
    {
        RequireNoNullElements(sequences, nameof(sequences), "Sequence");
        for (int i = 0; i < sequences.Count; i++)
        {
            if (sequences[i].Name is null)
                throw new ArgumentNullException(nameof(sequences), $"Sequence {i} has no name.");
            if (sequences[i].Residues is null)
                throw new ArgumentNullException(nameof(sequences), $"Sequence {i} has no residues.");
        }

        var o = options ?? new PoaOptions();
        RequireProgramInput(sequences.Count, "POA", "Error reading sequence file", nameof(sequences));
        using var cancelScope = AlignmentCancellation.Enter(o.CancellationToken);
        var m = matrix ?? PoaScoreMatrix.Parse(PoaDefaultMatrix.Text);
        return new PoaEngine(m, o).Run(sequences);
    }

    private sealed class PoaEngine
    {
        private readonly PoaScoreMatrix _m;
        private readonly PoaOptions _o;
        private readonly StringBuilder _stdout = new();

        public PoaEngine(PoaScoreMatrix m, PoaOptions o)
        {
            _m = m;
            _o = o;
        }

        /// <summary>The characters of <paramref name="s"/> as the bytes of their UTF-8 encoding (one char 0 ... 255 per byte); ASCII text is returned unchanged.</summary>
        internal static string PoaByteView(string s)
        {
            if (s.All(c => c < 0x80)) return s;
            var bytes = Encoding.UTF8.GetBytes(s);
            return string.Create(bytes.Length, bytes, static (span, b) => { for (int i = 0; i < b.Length; i++) span[i] = (char)b[i]; });
        }

        private static string AsciiCase(string s, bool upper)
            => string.Create(s.Length, (s, upper), static (span, st) =>
            {
                for (int i = 0; i < st.s.Length; i++)
                {
                    char c = st.s[i];
                    if (st.upper && c >= 'a' && c <= 'z') c = (char)(c - 32);
                    else if (!st.upper && c >= 'A' && c <= 'Z') c = (char)(c + 32);
                    span[i] = c;
                }
            });

        private static string Trunc31(string s) => s.Length > 31 ? s[..31] : s;

        // ----- sequences -----

        private PoaLpo MakeSequence(PoaSequence s)
        {
            // The program reads BYTES (C `char`): a non-ASCII character is its UTF-8 bytes, each one residue, none of them white space, upper / lower case (C locale) touches ASCII only; the
            // bytes outside the matrix's symbols become its first symbol (`limit_residues`). Verified against the installed poa: "ACG" + U+017F + "T" aligns as 6 residues (A C G A A T).
            string residues = PoaByteView(s.Residues);
            residues = string.Concat(residues.Where(c => !(c == ' ' || (c >= '\t' && c <= '\r'))));
            if (_o.Case == PoaCase.Lower) residues = AsciiCase(residues, false);
            else if (_o.Case == PoaCase.Upper) residues = AsciiCase(residues, true);
            return new PoaLpo { Name = Trunc31(s.Name), Title = string.IsNullOrEmpty(s.Title) ? "untitled" : s.Title!, Sequence = residues };
        }

        /// <summary><c>initialize_seqs_as_lpo</c>: <c>limit_residues</c>, <c>index_symbols</c>, <c>lpo_init</c>.</summary>
        private void Initialize(PoaLpo seq, int inputIndex)
        {
            string symbols = _m.Symbols;
            int len = seq.Sequence.Length;
            var letters = new PoaLetter[len];
            for (int i = 0; i < len; i++)
            {
                char c = seq.Sequence[i];
                if (symbols.IndexOf(c) < 0) c = symbols[0]; // limit_residues
                int k = symbols.Length - 1;
                for (int j = symbols.Length - 1; j >= 0; j--)
                {
                    if (symbols[j] == c)
                    {
                        k = j;
                        break;
                    }
                }

                var l = new PoaLetter { Letter = k, AlignRing = i, RingId = i };
                if (i > 0) l.Left.Add(i - 1);
                if (i < len - 1) l.Right.Add(i + 1);
                l.Source.Add((0, i));
                letters[i] = l;
            }

            seq.Letters = letters;
            seq.Sources.Add(new PoaSourceInfo { Name = seq.Name, Title = seq.Title, Length = len, Weight = 1, BundleId = PoaNoBundle, Input = inputIndex });
        }

        // ----- alignment (align_lpo_po2.c) -----

        private struct DpScore
        {
            public int Score;
            public short GapX, GapY;
        }

        private static int[][] LeftLinks(PoaLpo lpo)
        {
            int len = lpo.Length;
            var result = new int[len][];
            for (int i = 0; i < len; i++)
            {
                var l = lpo.Letters[i];
                bool initial = false;
                foreach (var (iseq, ipos) in l.Source)
                {
                    if (iseq >= 0 && ipos == 0) initial = true;
                }

                if (initial && !(l.Left.Count > 0 && l.Left[0] == -1))
                {
                    // all initial nodes (first residue of some sequence) are left-linked to -1: free alignment to any branch in a global alignment
                    var a = new int[l.Left.Count + 1];
                    a[0] = -1;
                    for (int k = 0; k < l.Left.Count; k++) a[k + 1] = l.Left[k];
                    result[i] = a;
                }
                else if (l.Left.Count == 0)
                {
                    result[i] = new[] { -1 };
                }
                else
                {
                    result[i] = l.Left.ToArray();
                }
            }

            return result;
        }

        private static bool[] FinalNodes(PoaLpo lpo)
        {
            int len = lpo.Length;
            var r = new bool[len];
            for (int i = 0; i < len; i++)
            {
                foreach (var (iseq, ipos) in lpo.Letters[i].Source)
                {
                    if (iseq >= 0 && ipos == lpo.Sources[iseq].Length - 1) r[i] = true;
                }
            }

            return r;
        }

        /// <summary><c>align_lpo_po</c>: partial-order alignment of two partial orders; returns the score and the two reciprocal letter maps.</summary>
        private int AlignLpo(PoaLpo lx, PoaLpo ly, out int[] xToY, out int[] yToX)
        {
            int lenX = lx.Length, lenY = ly.Length;
            var seqX = lx.Letters;
            var seqY = ly.Letters;
            var xLeft = LeftLinks(lx);
            var yLeft = LeftLinks(ly);
            var finalX = FinalNodes(lx);
            var finalY = FinalNodes(ly);
            bool global = _o.Global;
            const int MinScore = -999999;
            int bestScore = MinScore, bestX = -1, bestY = -1;

            int maxGap = _m.MaxGapLength;
            int[] gpx = (int[])_m.GapPenaltyX.Clone();
            int[] gpy = (int[])_m.GapPenaltyY.Clone();
            var nextGap = new int[maxGap + 2];
            var nextPerp = new int[maxGap + 2];
            for (int i = 0; i < maxGap + 1; i++)
            {
                nextGap[i] = i < maxGap ? i + 1 : i;
                nextPerp[i] = nextGap[i];
            }

            if (!global)
            {
                gpx[maxGap + 1] = gpy[maxGap + 1] = 0;
                nextGap[maxGap + 1] = nextPerp[maxGap + 1] = maxGap + 1;
            }
            else
            {
                gpx[maxGap + 1] = gpx[0];
                gpy[maxGap + 1] = gpy[0];
                nextGap[maxGap + 1] = nextGap[0];
                nextPerp[maxGap + 1] = nextPerp[0];
            }

            // the program computes lenX * lenY in a 32-bit int (it overflows above 2^31 cells: about 46 341 x 46 341) and callocs it; here a documented limit in 64-bit arithmetic
            RequireCells((long)lenX * lenY, PoaMaxMatrixCells, nameof(PoaMaxMatrixCells), "POA", "keeps a score row per node of the new sequence (8 bytes per cell) and two move matrices (1 byte per cell)",
                string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{lenX:N0} graph nodes x {lenY:N0} residues = {(long)lenX * lenY:N0} cells"), "The program has no limit of its own: its 32-bit cell count overflows above 2^31 cells and the allocation fails before",
                "align fewer or shorter sequences (or use MultipleAlignKalign / MultipleAlignMafft for sequences of this length)");
            var moveX = new byte[lenY * lenX];
            var moveY = new byte[lenY * lenX];
            // rows are shifted by one: index k + 1 holds position k, index 0 position -1
            var initCol = new DpScore[lenY + 1];
            var rows = new DpScore[lenY + 1][];
            rows[0] = new DpScore[lenX + 1];
            var cur = rows[0];
            cur[0].Score = 0;
            cur[0].GapX = cur[0].GapY = (short)(maxGap + 1);
            for (int i = 0; i < lenX; i++)
            {
                cur[i + 1].Score = MinScore;
                foreach (int xl in xLeft[i])
                {
                    int prevGap = cur[xl + 1].GapX;
                    int tryScore = cur[xl + 1].Score - gpx[prevGap];
                    if (tryScore > cur[i + 1].Score)
                    {
                        cur[i + 1].Score = tryScore;
                        cur[i + 1].GapX = (short)nextGap[prevGap];
                        cur[i + 1].GapY = (short)nextPerp[prevGap];
                    }
                }
            }

            initCol[0] = cur[0];
            for (int i = 0; i < lenY; i++)
            {
                initCol[i + 1].Score = MinScore;
                foreach (int yl in yLeft[i])
                {
                    int prevGap = initCol[yl + 1].GapY;
                    int tryScore = initCol[yl + 1].Score - gpy[prevGap];
                    if (tryScore > initCol[i + 1].Score)
                    {
                        initCol[i + 1].Score = tryScore;
                        initCol[i + 1].GapX = (short)nextPerp[prevGap];
                        initCol[i + 1].GapY = (short)nextGap[prevGap];
                    }
                }
            }

            for (int i = 0; i < lenY; i++)
            {
                if ((i & 15) == 0) AlignmentCancellation.ThrowIfRequested();
                rows[i + 1] = new DpScore[lenX + 1];
                cur = rows[i + 1];
                cur[0] = initCol[i + 1];
                for (int j = 0; j < lenX; j++)
                {
                    int matchScore = global ? MinScore : 0;
                    int matchX = 0, matchY = 0;
                    int insertXScore = MinScore, insertYScore = MinScore;
                    int insertXx = 0, insertYy = 0, insertXGap = 0, insertYGap = 0;
                    bool possibleEnd = !global || (finalX[j] && finalY[i]);
                    int ycount = 1;
                    foreach (int yl in yLeft[i])
                    {
                        var prev = rows[yl + 1];
                        int prevGap = prev[j + 1].GapY;
                        int tryScore = prev[j + 1].Score - gpy[prevGap];
                        if (tryScore > insertYScore)
                        {
                            insertYScore = tryScore;
                            insertYy = ycount;
                            insertYGap = prevGap;
                        }

                        int xcount = 1;
                        foreach (int xl in xLeft[j])
                        {
                            tryScore = prev[xl + 1].Score;
                            if (tryScore > matchScore)
                            {
                                matchScore = tryScore;
                                matchX = xcount;
                                matchY = ycount;
                            }

                            xcount++;
                        }

                        ycount++;
                    }

                    int xc = 1;
                    foreach (int xl in xLeft[j])
                    {
                        int prevGap = cur[xl + 1].GapX;
                        int tryScore = cur[xl + 1].Score - gpx[prevGap];
                        if (tryScore > insertXScore)
                        {
                            insertXScore = tryScore;
                            insertXx = xc;
                            insertXGap = prevGap;
                        }

                        xc++;
                    }

                    matchScore += _m.Score[seqX[j].Letter, seqY[i].Letter];
                    ref var my = ref cur[j + 1];
                    int mv = i * lenX + j;
                    if (matchScore > insertYScore && matchScore > insertXScore)
                    {
                        my.Score = matchScore;
                        my.GapX = 0;
                        my.GapY = 0;
                        moveX[mv] = (byte)matchX;
                        moveY[mv] = (byte)matchY;
                    }
                    else if (insertXScore > insertYScore)
                    {
                        my.Score = insertXScore;
                        my.GapX = (short)nextGap[insertXGap];
                        my.GapY = (short)nextPerp[insertXGap];
                        moveX[mv] = (byte)insertXx;
                        moveY[mv] = 0;
                    }
                    else
                    {
                        my.Score = insertYScore;
                        my.GapX = (short)nextPerp[insertYGap];
                        my.GapY = (short)nextGap[insertYGap];
                        moveX[mv] = 0;
                        moveY[mv] = (byte)insertYy;
                    }

                    if (possibleEnd && my.Score >= bestScore && (my.Score > bestScore || (j == bestX && i < bestY) || j < bestX))
                    {
                        bestScore = my.Score;
                        bestX = j;
                        bestY = i;
                    }
                }
            }

            if (bestX >= lenX || bestY >= lenY) throw new InvalidOperationException("Bounds exceeded!");
            var xAl = new int[lenX];
            var yAl = new int[lenY];
            Array.Fill(xAl, -1);
            Array.Fill(yAl, -1);
            int bx = bestX, by = bestY;
            while (bx >= 0 && by >= 0)
            {
                int xm = moveX[by * lenX + bx], ym = moveY[by * lenX + bx];
                if (xm > 0 && ym > 0)
                {
                    xAl[bx] = by;
                    yAl[by] = bx;
                }

                if (xm == 0 && ym == 0)
                {
                    xAl[bx] = by;
                    yAl[by] = bx;
                    break;
                }

                int nbx = bx, nby = by;
                if (xm > 0) nbx = xLeft[bx][xm - 1];
                if (ym > 0) nby = yLeft[by][ym - 1];
                bx = nbx;
                by = nby;
            }

            xToY = xAl;
            yToX = yAl;
            return bestScore;
        }

        // ----- fusion (lpo.c) -----

        private static void CrosslinkRings(int a, int b, PoaLetter[] seq)
        {
            if (seq[a].RingId == seq[b].RingId) return;
            int ring;
            if (seq[a].RingId < seq[b].RingId)
            {
                ring = b;
                do seq[ring].RingId = seq[a].RingId;
                while ((ring = seq[ring].AlignRing) != b);
            }
            else
            {
                ring = a;
                do seq[ring].RingId = seq[b].RingId;
                while ((ring = seq[ring].AlignRing) != a);
            }

            ring = seq[a].AlignRing;
            seq[a].AlignRing = seq[b].AlignRing;
            seq[b].AlignRing = ring;
        }

        private static int ReindexFusion(PoaLetter[] seqX, PoaLetter[] seqY, int[] xToY, int[] yToX, int[] newX, int[] newY)
        {
            int lenX = seqX.Length, lenY = seqY.Length;
            var doFuse = new bool[lenY];
            for (int iy = lenY; iy-- > 0;)
            {
                int ix = yToX[iy];
                if (ix >= 0 && seqX[ix].Letter == seqY[iy].Letter) doFuse[iy] = true;
            }

            int newLen = 0, iY = 0, endOfRing = -1;
            for (int iX = 0; iX < lenX; iX++)
            {
                for (int iRing = iX; iRing < lenX && seqX[iRing].RingId == seqX[iX].RingId; iRing++)
                {
                    if (xToY[iRing] >= 0)
                    {
                        while (iY < xToY[iRing]) newY[iY++] = newLen++;
                        break;
                    }
                }

                if (xToY[iX] >= 0 && iY < lenY)
                {
                    for (int iRing = seqY[iY].AlignRing; iRing != iY; iRing = seqY[iRing].AlignRing)
                    {
                        if (iRing > endOfRing) endOfRing = iRing;
                    }

                    if (doFuse[iY]) newY[iY++] = newLen;
                    else newY[iY++] = newLen++;
                }

                newX[iX] = newLen++;
                while (iY <= endOfRing) newY[iY++] = newLen++;
            }

            while (iY < lenY) newY[iY++] = newLen++;
            return newLen;
        }

        private static void AddLink(List<int> list, int ipos)
        {
            if (!list.Contains(ipos)) list.Add(ipos);
        }

        /// <summary><c>fuse_lpo</c>: fuses the partial order y into x along the alignment x_to_y / y_to_x (the result replaces the letters of x).</summary>
        private static void FuseLpo(PoaLpo hx, PoaLpo hy, int[] xToY, int[] yToX)
        {
            var seqX = hx.Letters;
            var seqY = hy.Letters;
            int lenX = seqX.Length, lenY = seqY.Length;
            var newX = new int[lenX];
            var newY = new int[lenY];
            int newLen = ReindexFusion(seqX, seqY, xToY, yToX, newX, newY);
            var lpo = new PoaLetter?[newLen];
            // translate the letters of x into the new indexing
            for (int i = 0; i < lenX; i++)
            {
                var l = seqX[i];
                l.AlignRing = newX[l.AlignRing];
                l.RingId = newX[l.RingId];
                for (int k = 0; k < l.Left.Count; k++) l.Left[k] = newX[l.Left[k]];
                for (int k = 0; k < l.Right.Count; k++) l.Right[k] = newX[l.Right[k]];
                lpo[newX[i]] = l;
            }

            // the letters of y that are not fused to a letter of x start blank
            for (int iy = lenY; iy-- > 0;)
            {
                if (yToX[iy] >= 0 && newX[yToX[iy]] == newY[iy]) continue;
                int i = newY[iy];
                lpo[i] = new PoaLetter { AlignRing = i, RingId = i };
            }

            // source list of y appended to the source list of x
            var iseqNew = new int[hy.Sources.Count];
            for (int s = 0; s < hy.Sources.Count; s++)
            {
                var src = hy.Sources[s];
                hx.Sources.Add(new PoaSourceInfo { Name = src.Name, Title = src.Title, Length = src.Length, Weight = src.Weight, BundleId = src.BundleId, Input = src.Input });
                iseqNew[s] = hx.Sources.Count - 1;
            }

            var newLetters = lpo.Select(l => l!).ToArray();
            for (int iy = lenY; iy-- > 0;)
            {
                var nl = newLetters[newY[iy]];
                var ol = seqY[iy];
                nl.Letter = ol.Letter;
                foreach (var (iseq, ipos) in ol.Source) nl.Source.Add((iseqNew[iseq], ipos));
                foreach (int l in ol.Left) AddLink(nl.Left, newY[l]);
                foreach (int r in ol.Right) AddLink(nl.Right, newY[r]);
            }

            for (int iy = lenY; iy-- > 0;)
            {
                // copy the old ring of y
                for (int ipos = iy, next; (next = seqY[ipos].AlignRing) != iy; ipos = next)
                    CrosslinkRings(newY[ipos], newY[next], newLetters);
            }

            for (int ix = lenX; ix-- > 0;)
            {
                if (xToY[ix] >= 0) CrosslinkRings(newX[ix], newY[xToY[ix]], newLetters);
            }

            hx.Letters = newLetters;
        }

        private static void FuseRingIdentities(PoaLetter[] seqX, PoaLetter[] seqY, int[] alX, int[] alY)
        {
            for (int i = seqY.Length; i-- > 0;)
            {
                if (alY[i] < 0 || seqX[alY[i]].Letter == seqY[i].Letter) continue;
                for (int j = seqX[alY[i]].AlignRing; j != alY[i]; j = seqX[j].AlignRing)
                {
                    if (seqX[j].Letter == seqY[i].Letter)
                    {
                        alX[alY[i]] = -1;
                        alY[i] = j;
                        alX[j] = i;
                        break;
                    }
                }
            }
        }

        // ----- progressive / iterative merging (buildup_lpo.c) -----

        private sealed record PairScore(double Score, int I, int J);

        private static int FindSeqName(List<PoaLpo> seqs, string name)
        {
            for (int i = 0; i < seqs.Count; i++)
            {
                if (seqs[i].Name == name) return i;
                if (seqs[i].Sources.Any(s => s.Name == name)) return i;
            }

            return -1;
        }

        private List<PairScore>? ReadPairScores(List<PoaLpo> seqs)
        {
            int nseq = seqs.Count;
            var adj = new int[nseq];
            var list = new List<PairScore>();
            double minScore = 0.0;
            if (_o.PairScores != null)
            {
                var tok = _o.PairScores.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
                for (int t = 0; t + 2 < tok.Length; t += 3)
                {
                    if (!double.TryParse(tok[t + 2], NumberStyles.Float, CultureInfo.InvariantCulture, out double x)) break;
                    int i = FindSeqName(seqs, tok[t]), j = FindSeqName(seqs, tok[t + 1]);
                    if (i < 0 || j < 0) return null; // "invalid sequence pair, not found": the program gives up merging and continues with the first sequence alone (exit status 0)
                    if (i < j) (i, j) = (j, i);
                    list.Add(new PairScore(x, i, j));
                    if (x < minScore) minScore = x;
                    if (j == i - 1) adj[i] = 1;
                }
            }
            else if (_o.Progressive)
            {
                for (int i = 0; i < nseq; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        double x = AlignLpo(seqs[i], seqs[j], out _, out _);
                        list.Add(new PairScore(x, i, j));
                        if (x < minScore) minScore = x;
                        if (j == i - 1) adj[i] = 1;
                    }
                }
            }

            for (int i = 1; i < nseq; i++)
            {
                if (adj[i] == 0) list.Add(new PairScore(minScore - 1.0, i, i - 1));
            }

            // descending by score, then ascending (i, j): the comparator of qsort is a total order
            list.Sort((a, b) =>
            {
                if (a.Score > b.Score) return -1;
                if (a.Score < b.Score) return 1;
                if (a.I != b.I) return a.I.CompareTo(b.I);
                return a.J.CompareTo(b.J);
            });
            return list;
        }

        private PoaLpo Buildup(List<PoaLpo> all)
        {
            int nseq = all.Count;
            if (nseq == 1) return all[0];
            var newSeq = all[0];
            var seqCluster = new int[nseq];
            var seqIdInCluster = new int[nseq];
            var clusterSize = new int[nseq];
            var initialNseq = new int[nseq];
            int nseqTot = 0;
            for (int i = 0; i < nseq; i++)
            {
                seqCluster[i] = i;
                clusterSize[i] = all[i].Sources.Count;
                initialNseq[i] = clusterSize[i];
                nseqTot += clusterSize[i];
            }

            var score = ReadPairScores(all);
            if (score == null) return all[0];
            foreach (var sc in score)
            {
                int clusterI, clusterJ;
                if (seqCluster[sc.I] < seqCluster[sc.J])
                {
                    clusterI = seqCluster[sc.I];
                    clusterJ = seqCluster[sc.J];
                }
                else if (seqCluster[sc.J] < seqCluster[sc.I])
                {
                    clusterI = seqCluster[sc.J];
                    clusterJ = seqCluster[sc.I];
                }
                else
                {
                    continue;
                }

                newSeq = all[clusterI];
                AlignLpo(newSeq, all[clusterJ], out int[] al1, out int[] al2);
                if (_o.FuseAll) FuseRingIdentities(newSeq.Letters, all[clusterJ].Letters, al1, al2);
                FuseLpo(newSeq, all[clusterJ], al1, al2);
                all[clusterJ].Letters = Array.Empty<PoaLetter>();
                for (int i = nseq; i-- > 0;)
                {
                    if (seqCluster[i] == clusterJ)
                    {
                        seqCluster[i] = clusterI;
                        seqIdInCluster[i] += clusterSize[clusterI];
                    }
                }

                clusterSize[clusterI] += clusterSize[clusterJ];
                clusterSize[clusterJ] = 0;
            }

            if (_o.PreserveSequenceOrder)
            {
                var perm = new int[nseqTot];
                nseqTot = 0;
                for (int i = 0; i < nseq; i++)
                    for (int j = 0; j < initialNseq[i]; j++)
                        perm[seqIdInCluster[i] + j] = nseqTot++;
                foreach (int p in perm) _stdout.Append(p.ToString(CultureInfo.InvariantCulture)).Append(' ');
                _stdout.Append('\n');
                ReindexSourceSeqs(newSeq, perm);
            }

            return newSeq;
        }

        private static void ReindexSourceSeqs(PoaLpo seq, int[] perm)
        {
            int len = seq.Length, nseq = seq.Sources.Count;
            var map = new int[Math.Max(len, nseq)];
            var invmap = new int[Math.Max(len, nseq)];
            for (int i = 0; i < nseq; i++) map[i] = invmap[i] = -1;
            for (int i = 0; i < nseq; i++)
            {
                map[i] = perm[i];
                if (map[i] >= nseq || map[i] < 0 || invmap[map[i]] != -1) throw new ArgumentException("Bad argument! 'perm' must be a permutation.");
                invmap[map[i]] = i;
            }

            for (int i = 0; i < len; i++)
            {
                var src = seq.Letters[i].Source;
                for (int k = 0; k < src.Count; k++) src[k] = (map[src[k].Iseq], src[k].Ipos);
            }

            for (int i = 0; i < nseq; i++)
            {
                int j = invmap[i];
                if (j == i) continue;
                (seq.Sources[i], seq.Sources[j]) = (seq.Sources[j], seq.Sources[i]);
                map[j] = map[i];
                invmap[map[i]] = j;
                map[i] = invmap[i] = i;
            }
        }

        // ----- heaviest bundle (heaviest_bundle.c) -----

        private static int[] HeaviestBundle(PoaLpo seq, out int bestLen)
        {
            int len = seq.Length, nsource = seq.Sources.Count;
            var path = new int[len];
            var score = new int[len];
            var contains = new int[nsource];
            int bestScore = -999999, ibest = -1;
            for (int i = len; i-- > 0;)
            {
                Array.Clear(contains);
                foreach (var (iseq, ipos) in seq.Letters[i].Source)
                    if (seq.Sources[iseq].Weight > 0) contains[iseq] = ipos + 1;
                int rightScore = 0, rightOverlap = 0, bestRight = -1;
                foreach (int right in seq.Letters[i].Right)
                {
                    int myOverlap = 0;
                    foreach (var (iseq, ipos) in seq.Letters[right].Source)
                        if (contains[iseq] == ipos) myOverlap += seq.Sources[iseq].Weight;
                    if (myOverlap > rightOverlap || (myOverlap == rightOverlap && score[right] > rightScore))
                    {
                        rightOverlap = myOverlap;
                        rightScore = score[right];
                        bestRight = right;
                    }
                }

                path[i] = bestRight;
                score[i] = rightScore + rightOverlap;
                if (score[i] > bestScore)
                {
                    ibest = i;
                    bestScore = score[i];
                }
            }

            var best = new List<int>();
            for (; ibest >= 0; ibest = path[ibest]) best.Add(ibest);
            bestLen = best.Count;
            return best.ToArray();
        }

        private void GenerateBundles(PoaLpo seq, float minimumFraction)
        {
            int nbundled = 0, ibundle = 0;
            while (nbundled < seq.Sources.Count)
            {
                var path = HeaviestBundle(seq, out int pathLength);
                if (path.Length == 0 || pathLength < 10) break;
                string name = "CONSENS" + ibundle.ToString(CultureInfo.InvariantCulture);
                // assign_sequence_bundle_id
                var bundleCount = new int[seq.Sources.Count];
                for (int i = pathLength; i-- > 0;)
                    foreach (var (iseq, _) in seq.Letters[path[i]].Source) bundleCount[iseq]++;
                int count = 0;
                for (int i = seq.Sources.Count; i-- > 0;)
                {
                    if (seq.Sources[i].BundleId < 0 && seq.Sources[i].Length * minimumFraction <= bundleCount[i])
                    {
                        seq.Sources[i].BundleId = ibundle;
                        seq.Sources[i].Weight = 0;
                        count++;
                    }
                }

                string title = "consensus produced by heaviest_bundle, containing " + count.ToString(CultureInfo.InvariantCulture) + " seqs";
                // add_path_sequence
                seq.Sources.Add(new PoaSourceInfo { Name = name, Title = title, Length = pathLength, Weight = 0, BundleId = PoaNoBundle });
                int iseqNew = seq.Sources.Count - 1;
                for (int i = pathLength; i-- > 0;) seq.Letters[path[i]].Source.Add((iseqNew, i));
                seq.Sources[iseqNew].BundleId = ibundle++;
                nbundled += count;
                if (count < 1) break;
            }
        }

        // ----- output (lpo_format.c) -----

        private char Symbol(int letter) => letter < _m.Symbols.Length ? _m.Symbols[letter] : (char)letter;

        private string WritePo(PoaLpo seq)
        {
            var sb = new StringBuilder();
            sb.Append("VERSION=LPO.0.1\n");
            sb.Append("NAME=").Append(seq.Name).Append("\nTITLE=").Append(seq.Title).Append("\nLENGTH=").Append(seq.Length).Append("\nSOURCECOUNT=").Append(seq.Sources.Count).Append('\n');
            foreach (var s in seq.Sources)
                sb.Append("SOURCENAME=").Append(s.Name).Append("\nSOURCEINFO=").Append(s.Length).Append(' ').Append(s.Istart).Append(' ').Append(s.Weight).Append(' ').Append(s.BundleId).Append(' ').Append(s.Title).Append('\n');
            for (int i = 0; i < seq.Length; i++)
            {
                var l = seq.Letters[i];
                sb.Append(Symbol(l.Letter)).Append(':');
                foreach (int link in l.Left) sb.Append('L').Append(link.ToString(CultureInfo.InvariantCulture));
                foreach (var (iseq, _) in l.Source) sb.Append('S').Append(iseq.ToString(CultureInfo.InvariantCulture));
                if (l.AlignRing != i) sb.Append('A').Append(l.AlignRing.ToString(CultureInfo.InvariantCulture));
                sb.Append('\n');
            }

            return sb.ToString();
        }

        /// <summary>
        /// <c>xlate_lpo_to_al</c>: the rows (one per source sequence) over the align rings. The program also computes the ring flags of a bundle selection (<c>include_in_save</c>) but never hands
        /// them back to its callers (<c>*p_include</c> is not assigned), so <c>-best</c> selects sequences only, never columns.
        /// </summary>
        private char[][] Translate(PoaLpo seq, char gap)
        {
            int nring = 0, current = 0;
            for (int i = 0; i < seq.Length; i++)
            {
                if (seq.Letters[i].RingId != current)
                {
                    current = seq.Letters[i].RingId;
                    nring++;
                }
            }

            nring++;
            int ns = seq.Sources.Count;
            var pos = new char[ns][];
            for (int i = 0; i < ns; i++)
            {
                pos[i] = new char[nring];
                Array.Fill(pos[i], gap);
            }

            current = 0;
            int iring = 0;
            for (int i = 0; i < seq.Length; i++)
            {
                if (seq.Letters[i].RingId != current)
                {
                    current = seq.Letters[i].RingId;
                    iring++;
                }

                foreach (var (iseq, _) in seq.Letters[i].Source)
                    pos[iseq][iring] = Symbol(seq.Letters[i].Letter);
            }

            return pos;
        }

        private string WritePir(PoaLpo seq, int ibundle)
        {
            var pos = Translate(seq, '.');
            var sb = new StringBuilder();
            for (int i = 0; i < seq.Sources.Count; i++)
            {
                if (ibundle >= 0 && seq.Sources[i].BundleId != ibundle) continue;
                sb.Append('>').Append(seq.Sources[i].Name).Append(' ').Append(seq.Sources[i].Title);
                int iprint = 0;
                for (int j = 0; j < pos[i].Length; j++)
                {
                    if (iprint % 60 == 0) sb.Append('\n');
                    sb.Append(pos[i][j]);
                    iprint++;
                }

                sb.Append('\n');
            }

            return sb.ToString();
        }

        private string WriteClustal(PoaLpo seq)
        {
            var sb = new StringBuilder("CLUSTAL W (1.74) multiple sequence alignment\n\n\n");
            var pos = Translate(seq, '-');
            const int indent = 36, blockSize = 50;
            int nring = pos.Length == 0 ? 0 : pos[0].Length;
            int nparagraph = nring / blockSize;
            if (nring % blockSize != 0) nparagraph++;
            for (int ip = 0; ip < nparagraph; ip++)
            {
                for (int iseq = 0; iseq < seq.Sources.Count; iseq++)
                {
                    string name = seq.Sources[iseq].Name;
                    if (name.Length > indent - 1) sb.Append(name, 0, indent - 1).Append(' ');
                    else sb.Append(name).Append(' ', indent - name.Length);
                    for (int i = 0; i < blockSize; i++)
                    {
                        int ipos = i + ip * blockSize;
                        if (ipos >= nring) break;
                        sb.Append(pos[iseq][ipos]);
                    }

                    sb.Append('\n');
                }

                sb.Append("\n\n");
            }

            return sb.ToString();
        }

        // ----- main -----

        public PoaRunResult Run(IReadOnlyList<PoaSequence> input)
        {
            if (_o.Verbose)
            {
                var g = _m.GapPenaltySet;
                _stdout.Append("X-Gap Penalties (Open, Aff1, Aff2; LTrunc, LDecay): ").Append(g[0, 0]).Append(' ').Append(g[0, 1]).Append(' ').Append(g[0, 2]).Append(' ').Append(_m.TruncGapLength).Append(' ').Append(_m.DecayGapLength).Append('\n');
                _stdout.Append("X-Gap Penalties (0, 1, 2, ...): ");
                for (int i = 0; i <= _m.MaxGapLength; i++) _stdout.Append(_m.GapPenaltyX[i]).Append(' ');
                _stdout.Append("... \n");
                _stdout.Append("Y-Gap Penalties (Open, Aff1, Aff2; LTrunc, LDecay): ").Append(g[1, 0]).Append(' ').Append(g[1, 1]).Append(' ').Append(g[1, 2]).Append(' ').Append(_m.TruncGapLength).Append(' ').Append(_m.DecayGapLength).Append('\n');
                _stdout.Append("Y-Gap Penalties (0, 1, 2, ...): ");
                for (int i = 0; i <= _m.MaxGapLength; i++) _stdout.Append(_m.GapPenaltyY[i]).Append(' ');
                _stdout.Append("... \n");
            }

            if (_o.PrintMatrix != null) _stdout.Append(_m.Print(_o.PrintMatrix));
            var all = new List<PoaLpo>();
            for (int k = 0; k < input.Count; k++)
            {
                var lpo = MakeSequence(input[k]);
                if (lpo.Sequence.Length == 0) continue; // a record without residues is not a sequence (create_seq is only called for non-empty data)
                Initialize(lpo, k);
                all.Add(lpo);
            }

            if (all.Count == 0) throw new ArgumentException("No input sequences were provided.");
            var result = Buildup(all);
            if (_o.Comment != null) result.Title = _o.Comment;
            if (_o.HeaviestBundle) GenerateBundles(result, _o.HeaviestBundleMinimum);
            var pos = Translate(result, '-');
            var rows = pos.Select(p => new string(p)).ToArray();
            var order = result.Sources.Select(x => x.Input).ToArray();

            return new PoaRunResult(result.Sources.Select(s => s.Name).ToArray(), result.Sources.Select(s => s.Title).ToArray(), rows, order,
                WritePo(result), WritePir(result, _o.Best ? 0 : -1), WriteClustal(result), _stdout.ToString());
        }
    }

    #endregion
}
