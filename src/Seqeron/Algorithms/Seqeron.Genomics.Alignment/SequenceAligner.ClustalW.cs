using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region ClustalW pairwise stage (percent identity / guide-tree distance)

    // ClustalW 2.1 "slow" (full dynamic programming) pairwise stage, FullPairwiseAlign::pairwiseAlign
    // (src/pairwise/FullPairwiseAlign.cpp of the clustalw 2.1+lgpl-7 Debian source package, opened via
    // sources.debian.org; Larkin et al. 2007 Bioinformatics 23:2947; Thompson, Higgins & Gibson 1994 NAR 22:4673).
    //
    // What the stage really is (read from the source, confirmed against the installed `clustalw -ALIGN`):
    //   * forwardPass: an affine-gap (Gotoh) SMITH-WATERMAN matrix (scores floored at 0), the end cell is the FIRST
    //     maximum in row-major order (strict >);
    //   * reversePass: the same recurrence run backwards from the end cell without the zero floor; the start cell is the
    //     first cell (reverse row-major) whose score reaches the forward maximum;
    //   * diff: Myers & Miller (1988) linear-space divide and conquer on the region [start..end] with ClustalW's own
    //     midpoint tie rules (tb/te terminal gap-open flags, type-1/type-2 midpoint);
    //   * tracePath counts aligned pairs with identical residue codes; mmScore = 100 * count / MIN(len1, len2) in single
    //     precision; the guide-tree distance is (100 - mmScore) / 100 (= 1 - count / MIN(len1, len2)), and `-ALIGN`
    //     prints (int)mmScore as "Score:".
    // So ClustalW's pairwise percent identity is the identity of a LOCAL affine alignment, not of a global one; its
    // co-optimal choice among ties (Myers-Miller midpoint rule) can change the identity count, which is why this method is
    // a faithful port and not a call to LocalAlignAffine (whose tie rule differs).
    //
    // Parameters (UserParameters.cpp / SubMatrix.cpp::getPairwiseMatrix, DNA): pairwise matrix = IUB (swgapdnamt of
    // matrices.h) x intScale 100 (A=C=G=T match +1000, every mismatch -900; ClustalW 2.1's transition-weight overrides
    // matrix[0][4], [2][11], [2][12] address residue codes A-E, C-M, C-N, so no A/C/G/T pair is changed), gap open
    // = (int)(2 * PWGapOpen * 100) = 3000 and gap extend = (int)(PWGapExtend * 100) = 666 for the defaults 15.00 / 6.66;
    // a gap of length k costs open + k * extend (the first gap position already pays one extension).

    /// <summary>ClustalW 2.1 default DNA pairwise gap-opening penalty (<c>-PWGAPOPEN</c>).</summary>
    public const double ClustalWDnaPairwiseGapOpen = 15.0;

    /// <summary>ClustalW 2.1 default DNA pairwise gap-extension penalty (<c>-PWGAPEXT</c>).</summary>
    public const double ClustalWDnaPairwiseGapExtend = 6.66;

    // ClustalW residue codes: index in "ABCDEFGHIKLMNOPQRSTUVWXYZ-" (UserParameters::aminoAcidCodes).
    private const string ClustalWResidueCodes = ClustalWUserMatrix.ResidueCodes;

    // Packing order of the DNA matrices of matrices.h and the IUB triangle (swgapdnamt), row i holds columns 0..i.
    private const string ClustalWNucleicOrder = "ABCDGHKMNRSTUVWXY";

    private static readonly short[] ClustalWIubTriangle =
    {
        10,
        -9, 10,
        -9, 10, 10,
        10, 10, -9, 10,
        -9, 10, -9, 10, 10,
        10, 10, 10, 10, -9, 10,
        -9, 10, -9, 10, 10, 10, 10,
        10, 10, 10, 10, -9, 10, -9, 10,
        10, 10, 10, 10, 10, 10, 10, 10, 10,
        10, 10, -9, 10, 10, 10, 10, 10, 10, 10,
        -9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
        -9, 10, -9, 10, -9, 10, 10, -9, 10, -9, -9, 10,
        -9, 10, -9, 10, -9, 10, 10, -9, 10, -9, -9, 10, 10,
        10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, -9, -9, 10,
        10, 10, -9, 10, -9, 10, 10, 10, 10, 10, -9, 10, 10, 10, 10,
        10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
        -9, 10, 10, 10, -9, 10, 10, 10, 10, -9, 10, 10, 10, 10, 10, 10, 10,
    };

    private const int ClustalWIntScale = 100;

    // ---- the nucleotide matrices of ClustalW 2.1 (matrices.h swgapdnamt / clustalvdnamt, the user's -DNAMATRIX / -PWDNAMATRIX file) ----

    /// <summary><c>defaultDNAXref</c>: for each of the 17 columns of the packing order of the DNA tables of <c>matrices.h</c> its residue code.</summary>
    private static readonly int[] ClustalWDefaultDnaXref = ClustalWNucleicOrder.Select(c => ClustalWResidueCodes.IndexOf(c)).ToArray();

    /// <summary><c>clustalvdnamt</c> (ClustalW 1.6's matrix): 10 on the diagonal of A, C, G, T and U, 0 everywhere else (lower triangle in the packing order).</summary>
    private static readonly short[] ClustalWClustalvDnaTriangle = BuildClustalvDnaTriangle();

    private static short[] BuildClustalvDnaTriangle()
    {
        var t = new short[ClustalWIubTriangle.Length];
        foreach (char c in "ACGTU")
        {
            int i = ClustalWNucleicOrder.IndexOf(c);
            t[i * (i + 1) / 2 + i] = 10;
        }

        return t;
    }

    private static IEnumerable<(int Row, int Column, short Value)> ClustalWDnaTableCells(short[] triangle)
    {
        int ix = 0;
        for (int i = 0; i < ClustalWNucleicOrder.Length; i++)
            for (int j = 0; j <= i; j++)
                yield return (i, j, triangle[ix++]);
    }

    /// <summary>
    /// <c>SubMatrix::getProfileAlignMatrix</c>, DNA branch: the table (<c>swgapdnamt</c> for IUB, <c>clustalvdnamt</c> for CLUSTALW with the gap-penalty scale 0.66,
    /// <c>userDNAMat</c> with <c>DNAXref</c> and scale 1.0 for a user matrix) with its cross reference.
    /// </summary>
    private static (IReadOnlyList<(int Row, int Column, short Value)> Cells, IReadOnlyList<int> Cross, float Scale) ClustalWProfileDnaMatrixSource(
        ClustalWDnaMatrix builtIn, ClustalWUserMatrix? user)
    {
        if (user is not null)
            return (user.Cells().ToList(), user.Cross, 1.0f);
        return builtIn switch
        {
            ClustalWDnaMatrix.Iub => (ClustalWDnaTableCells(ClustalWIubTriangle).ToList(), ClustalWDefaultDnaXref, 1.0f),
            ClustalWDnaMatrix.Clustalw => (ClustalWDnaTableCells(ClustalWClustalvDnaTriangle).ToList(), ClustalWDefaultDnaXref, 0.66f),
            _ => throw new ArgumentOutOfRangeException(nameof(builtIn), builtIn, "Unknown ClustalW DNA matrix."),
        };
    }

    /// <summary>
    /// <c>matrix[Xref[a]][Xref[b]] = matrix[Xref[b]][Xref[a]] = value</c> of <c>getProfileAlignMatrix</c>. A position outside the file's header has Xref -1: the C++ then writes
    /// outside its <c>matrix[NUMRES][NUMRES]</c> array (undefined behaviour); the write is dropped here.
    /// </summary>
    private static void ClustalWSetTransition(int[,] matrix, IReadOnlyList<int> cross, int a, int b, int value)
    {
        int x = a < cross.Count ? cross[a] : -1;
        int y = b < cross.Count ? cross[b] : -1;
        if (x < 0 || y < 0)
            return;
        matrix[x, y] = value;
        matrix[y, x] = value;
    }

    /// <summary>
    /// <c>SubMatrix::getPairwiseMatrix</c>, DNA branch: <c>getMatrix</c> with intScale 100 and no shift (<c>negFlag = true</c>) on the table, then the transition weights written
    /// to the residue CODES (0,4), (2,11), (2,12) = A-E, C-M, C-N (ClustalW 2.1 does not go through the Xref here, so with ambiguity codes the C-M / C-N scores move and no A / C / G / T
    /// pair is changed).
    /// </summary>
    private static int[,] BuildClustalWDnaPairwiseMatrix(IEnumerable<(int Row, int Column, short Value)> cells, IReadOnlyList<int> cross, float transitionWeight)
    {
        var (matrix, _) = ClustalWGetMatrix(cells, cross, ClustalWIntScale);
        int transition = (int)(transitionWeight * matrix[0, 0]);
        foreach (var (a, b) in new[] { (0, 4), (2, 11), (2, 12) })
        {
            matrix[a, b] = transition;
            matrix[b, a] = transition;
        }

        return matrix;
    }

    /// <summary>The pairwise matrix and the gap-penalty scales of <c>getPairwiseMatrix</c> for a DNA setup (CLUSTALW: opening x 0.6667, extension x 0.751).</summary>
    private static (int[,] Matrix, float GapOpenScale, float GapExtendScale) ClustalWDnaPairwiseSetup(ClustalWDnaPairwiseOptions o)
    {
        if (o.UserMatrix is not null && o.UserMatrixUse == ClustalWUserDnaMatrixUse.Applied)
            return (BuildClustalWDnaPairwiseMatrix(o.UserMatrix.Cells(), o.UserMatrix.Cross, o.TransitionWeight), 1.0f, 1.0f);
        return o.Matrix switch
        {
            ClustalWDnaMatrix.Iub => (BuildClustalWDnaPairwiseMatrix(ClustalWDnaTableCells(ClustalWIubTriangle), ClustalWDefaultDnaXref, o.TransitionWeight), 1.0f, 1.0f),
            ClustalWDnaMatrix.Clustalw => (BuildClustalWDnaPairwiseMatrix(ClustalWDnaTableCells(ClustalWClustalvDnaTriangle), ClustalWDefaultDnaXref, o.TransitionWeight), 0.6667f, 0.751f),
            _ => throw new ArgumentOutOfRangeException(nameof(o), o.Matrix, "Unknown ClustalW DNA matrix."),
        };
    }

    private static readonly Lazy<(int[,] Matrix, float GapOpenScale, float GapExtendScale)> ClustalWDefaultDnaPairwiseSetup = new(() => ClustalWDnaPairwiseSetup(new ClustalWDnaPairwiseOptions()));

    /// <summary>
    /// ClustalW 2.1 pairwise alignment stage for DNA (default <c>-PWDNAMATRIX=IUB -PWGAPOPEN=15 -PWGAPEXT=6.66</c>) with the
    /// percent identity and guide-tree distance ClustalW derives from it.
    /// </summary>
    /// <remarks>
    /// Port of <c>FullPairwiseAlign::pairwiseAlign</c> (the full-DP mode of <c>clustalw -ALIGN</c>/<c>-TREE</c>; the fast k-tuple
    /// mode selected by <c>-QUICKTREE</c> is <see cref="ClustalWQuickPairwiseAlign(string, string, ClustalWQuickOptions?, int?)"/>): an affine local
    /// (Smith–Waterman–Gotoh) alignment with the IUB matrix ×100, whose identities give ClustalW's percent identity.
    /// Because ClustalW's co-optimal choice (Myers–Miller midpoint rule) is part of the result, the alignment is reproduced
    /// step by step rather than by <see cref="LocalAlignAffine(string, string, ScoringMatrix?)"/>.
    /// Verified against the installed <c>clustalw</c> 2.1 (<c>Sequences (i:j) Aligned. Score:</c>, the truncated percent
    /// identity) on 11 000+ random, low-complexity and mutated DNA pairs.
    /// </remarks>
    /// <param name="sequence1">First sequence (A, C, G, T/U and the IUB ambiguity codes B D H K M N R S V W X Y; case-insensitive).</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="gapOpen">Pairwise gap-opening penalty (positive, ClustalW <c>-PWGAPOPEN</c>, default 15).</param>
    /// <param name="gapExtend">Pairwise gap-extension penalty (positive, <c>-PWGAPEXT</c>, default 6.66).</param>
    /// <exception cref="ArgumentNullException">A sequence is null.</exception>
    /// <exception cref="ArgumentException">A symbol is not an IUB nucleotide code, or a penalty is negative.</exception>
    public static ClustalWPairwiseResult ClustalWPairwiseAlign(
        string sequence1, string sequence2,
        double gapOpen = ClustalWDnaPairwiseGapOpen, double gapExtend = ClustalWDnaPairwiseGapExtend)
        => ClustalWPairwiseAlignCore(sequence1, sequence2, gapOpen, gapExtend, 0.5f);

    /// <summary>
    /// ClustalW 2.1 pairwise alignment stage for DNA with every option of the stage: <c>-PWGAPOPEN</c>, <c>-PWGAPEXT</c>, <c>-TRANSWEIGHT</c> and the matrix
    /// <c>-PWDNAMATRIX</c> (<c>IUB</c>, <c>CLUSTALW</c> or a user file), see <see cref="ClustalWDnaPairwiseOptions"/>.
    /// </summary>
    /// <param name="sequence1">First sequence (IUB nucleotide codes, case-insensitive).</param>
    /// <param name="sequence2">Second sequence.</param>
    /// <param name="options">Pairwise options (default: IUB, 15, 6.66, transition weight 0.5).</param>
    /// <exception cref="ArgumentNullException">A sequence or <paramref name="options"/> is null.</exception>
    /// <exception cref="ArgumentException">A symbol is not an IUB nucleotide code, or a penalty is negative / NaN.</exception>
    public static ClustalWPairwiseResult ClustalWPairwiseAlign(string sequence1, string sequence2, ClustalWDnaPairwiseOptions options)
    {
        ArgumentNullException.ThrowIfNull(options);
        return ClustalWPairwiseAlignCore(sequence1, sequence2, options.GapOpen, options.GapExtend, options.TransitionWeight, options);
    }

    /// <summary><see cref="ClustalWPairwiseAlign(string, string, ClustalWDnaPairwiseOptions)"/> for <see cref="DnaSequence"/> inputs.</summary>
    public static ClustalWPairwiseResult ClustalWPairwiseAlign(DnaSequence sequence1, DnaSequence sequence2, ClustalWDnaPairwiseOptions options)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return ClustalWPairwiseAlign(sequence1.Sequence, sequence2.Sequence, options);
    }

    /// <summary><see cref="ClustalWPairwiseAlign(string, string, double, double)"/> with ClustalW's <c>-TRANSWEIGHT</c> (affects the C-M / C-N scores of ambiguity codes).</summary>
    internal static ClustalWPairwiseResult ClustalWPairwiseAlignCore(
        string sequence1, string sequence2, double gapOpen, double gapExtend, float transitionWeight, ClustalWDnaPairwiseOptions? matrixOptions = null)
    {
        if (float.IsNaN(transitionWeight) || transitionWeight < 0)
            throw new ArgumentException("TransitionWeight must be a non-negative number.");
        (int[,] Matrix, float GapOpenScale, float GapExtendScale) setup;
        if (matrixOptions is not null)
            setup = ClustalWDnaPairwiseSetup(matrixOptions);
        else if (transitionWeight == 0.5f)
            setup = ClustalWDefaultDnaPairwiseSetup.Value;
        else
            setup = ClustalWDnaPairwiseSetup(new ClustalWDnaPairwiseOptions { TransitionWeight = transitionWeight });
        return ClustalWPairwiseAlignResolved(sequence1, sequence2, gapOpen, gapExtend, setup);
    }

    /// <summary>The pairwise stage with an already built matrix / gap-scale setup (<see cref="ClustalWDnaPairwiseSetup"/>).</summary>
    internal static ClustalWPairwiseResult ClustalWPairwiseAlignResolved(
        string sequence1, string sequence2, double gapOpen, double gapExtend, (int[,] Matrix, float GapOpenScale, float GapExtendScale) setup)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        if (gapOpen < 0 || gapExtend < 0 || double.IsNaN(gapOpen) || double.IsNaN(gapExtend))
            throw new ArgumentException("ClustalW pairwise gap penalties are non-negative numbers.");

        var codes1 = EncodeClustalWDna(sequence1, nameof(sequence1));
        var codes2 = EncodeClustalWDna(sequence2, nameof(sequence2));

        // FullPairwiseAlign: an empty sequence gives distance 1.0 (percent identity 0) without any alignment.
        if (codes1.Length == 0 || codes2.Length == 0)
            return new ClustalWPairwiseResult(0, 0, 0.0, 1.0, 0, "", "", -1, -1, -1, -1);

        // (int)(2 * PWGapOpen * intScale * gapOpenScale), (int)(PWGapExtend * intScale * gapExtendScale) in single precision.
        float go = (float)gapOpen;
        float ge = (float)gapExtend;
        int gapOpenInt = (int)(2f * go * (float)ClustalWIntScale * setup.GapOpenScale);
        int gapExtendInt = (int)(ge * (float)ClustalWIntScale * setup.GapExtendScale);

        return RunClustalWPairwise(sequence1, sequence2, codes1, codes2, setup.Matrix, gapOpenInt, gapExtendInt);
    }

    /// <summary>
    /// The alignment and percent-identity part of <c>FullPairwiseAlign::pairwiseAlign</c> shared by the DNA and PROTEIN stages (residue codes, the
    /// code-indexed integer matrix and the integer gap penalties are the only differences between them).
    /// </summary>
    private static ClustalWPairwiseResult RunClustalWPairwise(
        string sequence1, string sequence2, int[] codes1, int[] codes2, int[,] matrix, int gapOpenInt, int gapExtendInt)
    {
        var state = new ClustalWPairState(codes1, codes2, matrix, gapOpenInt, gapExtendInt);
        state.ForwardPass();
        state.ReversePass();
        state.PrintPtr = 1;
        state.LastPrint = 0;
        state.Diff(state.Sb1 - 1, state.Sb2 - 1, state.Se1 - state.Sb1 + 1, state.Se2 - state.Sb2 + 1, 0, 0);

        var (rows1, rows2, identities) = state.TracePath(sequence1, sequence2);
        int minLen = Math.Min(codes1.Length, codes2.Length);

        float mm = 100.0f * identities;
        mm /= minLen;
        int reported = (int)mm;

        bool empty = state.MaxScore == 0;
        return new ClustalWPairwiseResult(
            Score: state.MaxScore,
            Identities: identities,
            PercentIdentity: 100.0 * identities / minLen,
            Distance: 1.0 - (double)identities / minLen,
            ClustalWScore: reported,
            AlignedSequence1: rows1,
            AlignedSequence2: rows2,
            Start1: empty ? -1 : state.Sb1 - 1,
            End1: empty ? -1 : state.Se1 - 1,
            Start2: empty ? -1 : state.Sb2 - 1,
            End2: empty ? -1 : state.Se2 - 1);
    }

    /// <summary><see cref="ClustalWPairwiseAlign(string, string, double, double)"/> for <see cref="DnaSequence"/> inputs.</summary>
    public static ClustalWPairwiseResult ClustalWPairwiseAlign(
        DnaSequence sequence1, DnaSequence sequence2,
        double gapOpen = ClustalWDnaPairwiseGapOpen, double gapExtend = ClustalWDnaPairwiseGapExtend)
    {
        ArgumentNullException.ThrowIfNull(sequence1);
        ArgumentNullException.ThrowIfNull(sequence2);
        return ClustalWPairwiseAlign(sequence1.Sequence, sequence2.Sequence, gapOpen, gapExtend);
    }

    /// <summary>ClustalW's guide-tree distance of two sequences (1 when one is empty), see <see cref="MsaDistanceMode.ClustalWAffine"/>.</summary>
    private static double ClustalWPairwiseDistance(string a, string b) => ClustalWPairwiseAlign(a, b).Distance;

    private static int[] EncodeClustalWDna(string sequence, string paramName)
    {
        var codes = new int[sequence.Length];
        for (int i = 0; i < sequence.Length; i++)
        {
            char c = char.ToUpperInvariant(sequence[i]);
            int code = ClustalWNucleicOrder.Contains(c) ? ClustalWResidueCodes.IndexOf(c) : -1;
            if (code < 0)
                throw new ArgumentException($"'{sequence[i]}' at position {i} is not an IUB nucleotide code known to the ClustalW DNA matrix.", paramName);
            codes[i] = code;
        }

        return codes;
    }

    /// <summary>
    /// Mutable state of one pairwise alignment, a line-by-line port of the members of <c>FullPairwiseAlign</c>
    /// (1-based residue indices as in the C++ source; HH/DD/RR/SS are the Gotoh work rows, <c>displ</c> the edit script).
    /// </summary>
    private sealed class ClustalWPairState
    {
        private readonly int[] _a;   // 1-based residue codes of sequence 1
        private readonly int[] _b;
        private readonly int[,] _matrix;
        private readonly int _gapOpen;
        private readonly int _gapExtend;
        private readonly int _n;
        private readonly int _m;
        private readonly int[] _hh;
        private readonly int[] _dd;
        private readonly int[] _rr;
        private readonly int[] _ss;
        private readonly int[] _displ;

        public int PrintPtr;
        public int LastPrint;
        public int MaxScore;
        public int Sb1, Sb2, Se1, Se2;

        public ClustalWPairState(int[] codes1, int[] codes2, int[,] matrix, int gapOpen, int gapExtend)
        {
            _n = codes1.Length;
            _m = codes2.Length;
            _a = new int[_n + 1];
            _b = new int[_m + 1];
            Array.Copy(codes1, 0, _a, 1, _n);
            Array.Copy(codes2, 0, _b, 1, _m);
            _matrix = matrix;
            _gapOpen = gapOpen;
            _gapExtend = gapExtend;
            int size = Math.Max(_n, _m) + 2;
            _hh = new int[size];
            _dd = new int[size];
            _rr = new int[size];
            _ss = new int[size];
            _displ = new int[2 * (_n + _m) + 8];
        }

        private int Score(int i, int j) => _matrix[_a[i], _b[j]];

        private int TbGap(int k, int tb) => k <= 0 ? 0 : tb + _gapExtend * k;

        private void Add(int v)
        {
            if (LastPrint < 0)
            {
                _displ[PrintPtr - 1] = v;
                _displ[PrintPtr++] = LastPrint;
            }
            else
            {
                LastPrint = _displ[PrintPtr++] = v;
            }
        }

        private void Del(int k)
        {
            if (LastPrint < 0)
                LastPrint = _displ[PrintPtr - 1] -= k;
            else
                LastPrint = _displ[PrintPtr++] = -k;
        }

        public void ForwardPass()
        {
            MaxScore = 0;
            Se1 = Se2 = 0;
            for (int i = 0; i <= _m; i++)
            {
                _hh[i] = 0;
                _dd[i] = -_gapOpen;
            }

            for (int i = 1; i <= _n; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                int hh = 0, p = 0;
                int f = -_gapOpen;
                for (int j = 1; j <= _m; j++)
                {
                    f -= _gapExtend;
                    int t = hh - _gapOpen - _gapExtend;
                    if (f < t) f = t;

                    _dd[j] -= _gapExtend;
                    t = _hh[j] - _gapOpen - _gapExtend;
                    if (_dd[j] < t) _dd[j] = t;

                    hh = p + Score(i, j);
                    if (hh < f) hh = f;
                    if (hh < _dd[j]) hh = _dd[j];
                    if (hh < 0) hh = 0;

                    p = _hh[j];
                    _hh[j] = hh;

                    if (hh > MaxScore)
                    {
                        MaxScore = hh;
                        Se1 = i;
                        Se2 = j;
                    }
                }
            }
        }

        public void ReversePass()
        {
            int cost = 0;
            Sb1 = Sb2 = 1;
            for (int i = Se2; i > 0; i--)
            {
                _hh[i] = -1;
                _dd[i] = -1;
            }

            for (int i = Se1; i > 0; i--)
            {
                AlignmentCancellation.ThrowIfRequested();
                int hh = -1, f = -1;
                int p = i == Se1 ? 0 : -1;
                for (int j = Se2; j > 0; j--)
                {
                    f -= _gapExtend;
                    int t = hh - _gapOpen - _gapExtend;
                    if (f < t) f = t;

                    _dd[j] -= _gapExtend;
                    t = _hh[j] - _gapOpen - _gapExtend;
                    if (_dd[j] < t) _dd[j] = t;

                    hh = p + Score(i, j);
                    if (hh < f) hh = f;
                    if (hh < _dd[j]) hh = _dd[j];

                    p = _hh[j];
                    _hh[j] = hh;

                    if (hh > cost)
                    {
                        cost = hh;
                        Sb1 = i;
                        Sb2 = j;
                        if (cost >= MaxScore) break;
                    }
                }

                if (cost >= MaxScore) break;
            }
        }

        /// <summary>Myers–Miller divide and conquer (<c>FullPairwiseAlign::diff</c>, the returned score is not needed); appends the edit script to <c>displ</c>.</summary>
        public void Diff(int a, int b, int m, int n, int tb, int te)
        {
            if (n <= 0)
            {
                if (m > 0) Del(m);
                return;
            }

            if (m <= 1)
            {
                if (m <= 0)
                {
                    Add(n);
                    return;
                }

                int bestH = -(tb + _gapExtend) - TbGap(n, te);
                int h = -(te + _gapExtend) - TbGap(n, tb);
                if (h > bestH) bestH = h;
                int bestJ = 0;
                for (int j = 1; j <= n; j++)
                {
                    h = Score(1 + a, j + b) - TbGap(n - j, te) - TbGap(j - 1, tb);
                    if (h > bestH)
                    {
                        bestH = h;
                        bestJ = j;
                    }
                }

                if (bestJ == 0)
                {
                    Del(1);
                    Add(n);
                }
                else
                {
                    if (bestJ > 1) Add(bestJ - 1);
                    _displ[PrintPtr++] = LastPrint = 0;
                    if (bestJ < n) Add(n - bestJ);
                }

                return;
            }

            // Divide: find the optimum midpoint (midi, midj).
            int midi = m / 2;
            _hh[0] = 0;
            int tt = -tb;
            for (int j = 1; j <= n; j++)
            {
                _hh[j] = tt = tt - _gapExtend;
                _dd[j] = tt - _gapOpen;
            }

            tt = -tb;
            for (int i = 1; i <= midi; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                int s = _hh[0];
                int hh;
                _hh[0] = hh = tt = tt - _gapExtend;
                int f = tt - _gapOpen;
                for (int j = 1; j <= n; j++)
                {
                    if ((hh = hh - _gapOpen - _gapExtend) > (f = f - _gapExtend)) f = hh;
                    int e;
                    if ((hh = _hh[j] - _gapOpen - _gapExtend) > (e = _dd[j] - _gapExtend)) e = hh;
                    hh = s + Score(i + a, j + b);
                    if (f > hh) hh = f;
                    if (e > hh) hh = e;

                    s = _hh[j];
                    _hh[j] = hh;
                    _dd[j] = e;
                }
            }

            _dd[0] = _hh[0];

            _rr[n] = 0;
            tt = -te;
            for (int j = n - 1; j >= 0; j--)
            {
                _rr[j] = tt = tt - _gapExtend;
                _ss[j] = tt - _gapOpen;
            }

            tt = -te;
            for (int i = m - 1; i >= midi; i--)
            {
                AlignmentCancellation.ThrowIfRequested();
                int s = _rr[n];
                int hh;
                _rr[n] = hh = tt = tt - _gapExtend;
                int f = tt - _gapOpen;
                for (int j = n - 1; j >= 0; j--)
                {
                    if ((hh = hh - _gapOpen - _gapExtend) > (f = f - _gapExtend)) f = hh;
                    int e;
                    if ((hh = _rr[j] - _gapOpen - _gapExtend) > (e = _ss[j] - _gapExtend)) e = hh;
                    hh = s + Score(i + 1 + a, j + 1 + b);
                    if (f > hh) hh = f;
                    if (e > hh) hh = e;

                    s = _rr[j];
                    _rr[j] = hh;
                    _ss[j] = e;
                }
            }

            _ss[n] = _rr[n];

            int midh = _hh[0] + _rr[0];
            int midj = 0;
            int type = 1;
            for (int j = 0; j <= n; j++)
            {
                int hh = _hh[j] + _rr[j];
                if (hh >= midh && (hh > midh || (_hh[j] != _dd[j] && _rr[j] == _ss[j])))
                {
                    midh = hh;
                    midj = j;
                }
            }

            for (int j = n; j >= 0; j--)
            {
                int hh = _dd[j] + _ss[j] + _gapOpen;
                if (hh > midh)
                {
                    midh = hh;
                    midj = j;
                    type = 2;
                }
            }

            // Conquer recursively around the midpoint.
            if (type == 1)
            {
                Diff(a, b, midi, midj, tb, _gapOpen);
                Diff(a + midi, b + midj, m - midi, n - midj, _gapOpen, te);
            }
            else
            {
                Diff(a, b, midi - 1, midj, tb, 0);
                Del(2);
                Diff(a + midi + 1, b + midj, m - midi - 1, n - midj, 0, te);
            }
        }

        /// <summary>
        /// <c>tracePath</c>: walks the edit script from (Sb1, Sb2), counts identical aligned code pairs and renders the
        /// two aligned rows of the local region.
        /// </summary>
        public (string Row1, string Row2, int Identities) TracePath(string original1, string original2)
        {
            int toDo = PrintPtr - 1;
            int i1 = Sb1, i2 = Sb2;
            int count = 0;
            var r1 = new System.Text.StringBuilder();
            var r2 = new System.Text.StringBuilder();
            for (int i = 1; i <= toDo; i++)
            {
                int d = _displ[i];
                if (d == 0)
                {
                    if (_a[i1] == _b[i2]) count++;
                    r1.Append(char.ToUpperInvariant(original1[i1 - 1]));
                    r2.Append(char.ToUpperInvariant(original2[i2 - 1]));
                    i1++;
                    i2++;
                }
                else if (d > 0)
                {
                    for (int k = 0; k < d; k++)
                    {
                        r1.Append('-');
                        r2.Append(char.ToUpperInvariant(original2[i2 - 1]));
                        i2++;
                    }
                }
                else
                {
                    for (int k = 0; k < -d; k++)
                    {
                        r1.Append(char.ToUpperInvariant(original1[i1 - 1]));
                        r2.Append('-');
                        i1++;
                    }
                }
            }

            return (r1.ToString(), r2.ToString(), count);
        }
    }

    #endregion
}
