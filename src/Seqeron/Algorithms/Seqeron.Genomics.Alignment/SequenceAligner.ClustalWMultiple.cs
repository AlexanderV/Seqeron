using System.Globalization;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Phylogenetics;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region ClustalW multiple-alignment stage (guide tree -> weights -> progressive profile alignment)

    // ClustalW 2.1 multiple-alignment stage for DNA: Clustal::align -> TreeInterface::getWeightsAndStepsFromDistMat ->
    // MSA::multiSeqAlign -> MyersMillerProfileAlign::profileAlign (src/Clustal.cpp, src/tree/Tree.cpp, src/tree/TreeInterface.cpp,
    // src/multipleAlign/{MSA,MyersMillerProfileAlign,ProfileBase,ProfileStandard,ProfileWithSub}.cpp, src/substitutionMatrix/SubMatrix.cpp
    // of the clustalw 2.1+lgpl-7 Debian source package, opened via sources.debian.org; Thompson, Higgins & Gibson 1994
    // NAR 22:4673; Larkin et al. 2007 Bioinformatics 23:2947).
    //
    // What the stage does (and what each step ports):
    //   1. pairwise distances (full DP or the fast k-tuple stage; A1-15 / A1-15b) -> neighbour-joining guide tree (A1-12), written
    //      with 5 decimals, re-read in single precision and re-rooted (Tree::reRoot);
    //   2. Tree::calcSeqWeights: the weight of a sequence is sum(branch length / number of leaves below the branch) along its
    //      path to the root, normalised to INT_SCALE_FACTOR (1000 on Unix builds), minimum 1 (Thompson 1994 "branch-length-derived
    //      weights"); -NOWEIGHTS gives equal weights;
    //   3. Tree::calcSimilarities: the pairwise identities used from here on are 100 - 100 * (tree path length), clamped to
    //      [0.01, 1] (they are NOT the pairwise-stage identities);
    //   4. MSA::multiSeqAlign: the internal nodes of the rooted tree in post-order are the alignment steps; sequences whose
    //      best identity to any other sequence is <= -MAXDIV (30) are DELAYED (skipped in the first pass) and then added one by one,
    //      closest to the aligned set first, with weight = tree weight * identity to the new sequence (renormalised);
    //   5. profileAlign: both groups -> profiles (residue weights from the sequence weights, substitution matrix applied to one
    //      profile, 10 * relative weight in the other); position-specific gap penalties (positions that already hold gaps get
    //      penalty * (nonGap fraction) and half the extension penalty; end gaps are free) and Myers & Miller (1988)
    //      linear-space profile alignment with the ClustalW midpoint rule (progDiff).
    //
    // DNA specifics read from the source: the residue-specific (Pascarella & Argos), hydrophilic and variable-region gap penalties
    // and the gap-separation distance (-GAPDIST) are guarded by `if (!_DNAFlag)` in ProfileBase::calcGapCoeff and therefore do not
    // exist for DNA (-NOPGAP/-NOHGAP/-NOVGAP/-GAPDIST change nothing); they are the protein extension (A1-14).
    // Profile matrix = IUB x 100 shifted to be non-negative (the lowest score, -900, is subtracted from every residue pair unless
    // -NEGATIVE), then A-G / C-T / C-U are set to int(TRANSWEIGHT * matrix[A][A]); gap open = (int)(100 * GAPOPEN) and
    // gap extension = (int)(100.0 * GAPEXT) with GAPEXT stored as a float (6.66f = 6.6599998, so 665 and not 666).
    //
    // Arithmetic that decides results is reproduced as in the C++: float (single precision) for weights, tree distances and the
    // profile cells (10 * w / sum), int truncation everywhere ClustalW truncates.

    private const int CwNumRes = 32;
    private const int CwGapPos1 = CwNumRes - 2;       // code of the gaps inserted by ClustalW
    private const int CwGapPos2 = CwNumRes - 1;
    private const int CwMaxAa = 24;                   // aminoAcidCodes.length - 2
    private const int CwGapCol = 32;                  // position of the gap-open penalty in a profile column
    private const int CwLenCol = 33;                  // position of the gap-extension penalty
    private const int CwIntScaleFactor = 1000;        // INT_SCALE_FACTOR of OS_UNIX builds (clustalw.h)
    private const float CwReducedGap = 1.0f;
    private const int ClustalWMaxHydrophilicResidues = 9;   // MAXHYDRESIDUES of clustalw.h

    /// <summary>
    /// ClustalW 2.1's multiple alignment of DNA sequences (<c>clustalw -TYPE=DNA -ALIGN</c>): the guide tree is the neighbour-joining
    /// tree of ClustalW's pairwise distances, sequences are weighted by it (Thompson et al. 1994), divergent sequences are delayed,
    /// and profiles are aligned with position-specific gap penalties and the IUB weight matrix.
    /// </summary>
    /// <remarks>
    /// A step-by-step port, not an approximation: the rows equal the output of the installed <c>clustalw</c> 2.1 (<c>-OUTORDER=INPUT</c>)
    /// on the oracle sets of <c>clustalw_msa_full_oracle.tsv</c>. DNA only (the
    /// protein-specific penalties and matrices of ClustalW's profile stage are not part of this method); the scoring matrix is used only
    /// for the reported <see cref="MultipleAlignmentResult.TotalScore"/>, which is the natural affine sum of pairs
    /// (<see cref="SumOfPairsScore"/>) under that matrix.
    /// </remarks>
    /// <param name="sequences">DNA sequences (output rows keep the input order).</param>
    /// <param name="options">ClustalW parameters (default: ClustalW's DNA defaults).</param>
    /// <param name="scoring">Scoring used for the reported total score only (default <see cref="SimpleDna"/>, affine).</param>
    public static MultipleAlignmentResult MultipleAlignClustalW(
        IEnumerable<DnaSequence> sequences, ClustalWOptions? options = null, ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var seqList = sequences.ToList();
        if (seqList.Count == 0)
            return MultipleAlignmentResult.Empty;

        var rows = MultipleAlignClustalWCore(seqList.Select(s => s.Sequence).ToArray(), options ?? new ClustalWOptions());
        return BuildClustalWResult(rows, scoring ?? SimpleDna);
    }

    /// <summary><see cref="MultipleAlignClustalW(IEnumerable{DnaSequence}, ClustalWOptions?, ScoringMatrix?)"/> for strings (A, C, G, T/U and the IUB ambiguity codes).</summary>
    public static MultipleAlignmentResult MultipleAlignClustalW(
        IEnumerable<string> sequences, ClustalWOptions? options = null, ScoringMatrix? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var seqList = sequences.ToList();
        if (seqList.Count == 0)
            return MultipleAlignmentResult.Empty;

        return BuildClustalWResult(MultipleAlignClustalWCore(seqList.ToArray(), options ?? new ClustalWOptions()), scoring ?? SimpleDna);
    }

    /// <summary>
    /// ClustalW 2.1's multiple alignment of PROTEIN sequences (<c>clustalw -TYPE=PROTEIN -ALIGN</c>): the same stage as the DNA overload (neighbour-joining guide tree
    /// of ClustalW's pairwise distances, tree weights, delayed divergent sequences, Myers–Miller profile alignment) with the protein parameters of ClustalW's profile
    /// stage — the matrix series (GONNET / BLOSUM / PAM / ID, chosen per profile pair by its percent identity), the matrix- and length-scaled gap opening, the
    /// residue-specific (Pascarella &amp; Argos 1992), hydrophilic-run and variable-region gap penalties and the gap-separation distance.
    /// </summary>
    /// <remarks>
    /// A step-by-step port: the rows equal the output of the installed <c>clustalw</c> 2.1 (<c>-OUTORDER=INPUT</c>) on the oracle sets of
    /// <c>clustalw_protein_msa_oracle.tsv</c>. Residues are the letters A-Z except J (ClustalW silently drops J and '*'; they are rejected here). The reported
    /// <see cref="MultipleAlignmentResult.TotalScore"/> is the natural affine sum of pairs under <paramref name="scoring"/> (default BLOSUM62, gap open 10 / extend 0.2
    /// in EMBOSS convention; with the default the residues O and U, which BLOSUM62 lacks, score as X).
    /// </remarks>
    /// <param name="sequences">Protein sequences (output rows keep the input order).</param>
    /// <param name="options">ClustalW PROTEIN parameters (default: ClustalW's PROTEIN defaults).</param>
    /// <param name="scoring">Scoring used for the reported total score only.</param>
    public static MultipleAlignmentResult MultipleAlignClustalW(
        IEnumerable<ProteinSequence> sequences, ClustalWProteinOptions? options = null, MatrixScoring? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var seqList = sequences.ToList();
        if (seqList.Count == 0)
            return MultipleAlignmentResult.Empty;

        return MultipleAlignClustalW(seqList.Select(s => s.Sequence), options ?? new ClustalWProteinOptions(), scoring);
    }

    /// <summary><see cref="MultipleAlignClustalW(IEnumerable{ProteinSequence}, ClustalWProteinOptions?, MatrixScoring?)"/> for strings (letters A-Z except J).</summary>
    public static MultipleAlignmentResult MultipleAlignClustalW(
        IEnumerable<string> sequences, ClustalWProteinOptions options, MatrixScoring? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentNullException.ThrowIfNull(options);
        var seqList = sequences.ToList();
        if (seqList.Count == 0)
            return MultipleAlignmentResult.Empty;

        string[] rows = MultipleAlignClustalWProteinCore(seqList.ToArray(), options);
        var list = rows.ToList();
        string consensus = BuildConsensus(list, list.Max(r => r.Length), anySymbol: true);
        var sp = scoring ?? MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.2);
        string[] scored = scoring is null ? rows.Select(r => r.Replace('O', 'X').Replace('U', 'X')).ToArray() : rows;
        double exact = rows.Length == 1 ? 0.0 : MatrixSumOfPairsScore(scored, sp, MsaGapModel.Affine);
        return new MultipleAlignmentResult(rows, consensus, (int)Math.Round(exact, MidpointRounding.AwayFromZero)) { ExactTotalScore = exact };
    }

    private static string[] MultipleAlignClustalWProteinCore(string[] raw, ClustalWProteinOptions options)
    {
        ValidateClustalWProteinMultipleOptions(options);
        int n = raw.Length;
        var codes = new int[n][];
        var upper = new string[n];
        for (int i = 0; i < n; i++)
        {
            ArgumentNullException.ThrowIfNull(raw[i]);
            codes[i] = EncodeClustalWProtein(raw[i], "sequences");
            upper[i] = raw[i].ToUpperInvariant();
        }

        if (n == 1)
            return upper;

        for (int i = 0; i < n; i++)
        {
            if (codes[i].Length == 0)
                throw new ArgumentException("ClustalW cannot align an empty sequence.", nameof(raw));
        }

        // UserParameters::setHydResidues: at most MAXHYDRESIDUES (9) leading letters, upper-cased; an empty result keeps the default
        string hyd = new string(options.HydrophilicResidues.Take(ClustalWMaxHydrophilicResidues).TakeWhile(char.IsAsciiLetter).Select(char.ToUpperInvariant).ToArray());
        var config = CwConfig.FromProtein(options) with { HydrophilicResidues = hyd.Length > 0 ? hyd : "GPSNDQEKR" };
        var state = new ClustalWMsaState(codes, config);
        return state.Run(ClustalWProteinDistanceMatrix(upper, options));
    }

    private static void ValidateClustalWProteinMultipleOptions(ClustalWProteinOptions o)
    {
        if (o.GapOpen < 0 || o.GapExtend < 0 || double.IsNaN(o.GapOpen) || double.IsNaN(o.GapExtend))
            throw new ArgumentException("ClustalW gap penalties are non-negative numbers.", nameof(o));
        if (!Enum.IsDefined(o.Matrix))
            throw new ArgumentOutOfRangeException(nameof(o), o.Matrix, "Unknown profile matrix series.");
        if (!Enum.IsDefined(o.PairwiseMode))
            throw new ArgumentOutOfRangeException(nameof(o), o.PairwiseMode, "Unknown pairwise mode.");
        if (o.GapSeparationDistance < 0)
            throw new ArgumentOutOfRangeException(nameof(o), o.GapSeparationDistance, "GapSeparationDistance is a non-negative number of columns.");
        ArgumentNullException.ThrowIfNull(o.HydrophilicResidues);
        ArgumentNullException.ThrowIfNull(o.Pairwise);
        ValidateClustalWProteinOptions(o.Pairwise);
        ValidateClustalWQuickOptions(o.Pairwise.Quick, protein: true);
    }

    /// <summary>Distances of the PROTEIN pairwise stage as ClustalW stores them (full stage in single precision, fast stage in double precision).</summary>
    private static double[,] ClustalWProteinDistanceMatrix(string[] upper, ClustalWProteinOptions options)
    {
        int n = upper.Length;
        var d = new double[n, n];
        int longest = upper.Max(s => s.Length);
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                double v = options.PairwiseMode == ClustalWPairwiseMode.Quick
                    ? ClustalWQuickPairwiseAlign(upper[i], upper[j], ClustalWSequenceType.Protein, options.Pairwise.Quick, longest).Distance
                    : ClustalWProteinSinglePrecisionDistance(upper[i], upper[j], options.Pairwise);
                d[i, j] = v;
                d[j, i] = v;
            }
        }

        return d;
    }

    private static MultipleAlignmentResult BuildClustalWResult(string[] rows, ScoringMatrix scoring)
    {
        var list = rows.ToList();
        int maxLen = list.Max(r => r.Length);
        string consensus = BuildConsensus(list, maxLen);
        int sp = list.Count == 1 ? 0 : SumOfPairsScore(rows, scoring, MsaGapModel.Affine);
        return new MultipleAlignmentResult(rows, consensus, sp);
    }

    /// <summary>
    /// The alignment rows of the ClustalW stage in input order; also exposes the intermediate guide tree artefacts for tests.
    /// </summary>
    private static string[] MultipleAlignClustalWCore(string[] raw, ClustalWOptions options)
    {
        ValidateClustalWOptions(options);
        int n = raw.Length;
        var codes = new int[n][];
        for (int i = 0; i < n; i++)
            codes[i] = EncodeClustalWDna(raw[i], "sequences");

        if (n == 1)
            return new[] { raw[0] };

        // Clustal::align removes nothing for valid input; empty sequences cannot be read by ClustalW.
        for (int i = 0; i < n; i++)
        {
            if (codes[i].Length == 0)
                throw new ArgumentException("ClustalW cannot align an empty sequence.", nameof(raw));
        }

        var state = new ClustalWMsaState(codes, CwConfig.FromDna(options));
        return state.Run(ClustalWPairwiseDistanceMatrix(raw, options));
    }

    private static void ValidateClustalWOptions(ClustalWOptions o)
    {
        if (o.GapOpen < 0 || o.GapExtend < 0 || double.IsNaN(o.GapOpen) || double.IsNaN(o.GapExtend))
            throw new ArgumentException("ClustalW gap penalties are non-negative numbers.", nameof(o));
        if (o.PairwiseGapOpen < 0 || o.PairwiseGapExtend < 0 || double.IsNaN(o.PairwiseGapOpen) || double.IsNaN(o.PairwiseGapExtend))
            throw new ArgumentException("ClustalW pairwise gap penalties are non-negative numbers.", nameof(o));
        if (float.IsNaN(o.TransitionWeight) || o.TransitionWeight < 0)
            throw new ArgumentException("TransitionWeight must be a non-negative number.", nameof(o));
        if (!Enum.IsDefined(o.PairwiseMode))
            throw new ArgumentOutOfRangeException(nameof(o), o.PairwiseMode, "Unknown pairwise mode.");
        ValidateClustalWQuickOptions(o.Quick ?? new ClustalWQuickOptions());
    }

    /// <summary>Distances of the pairwise stage as ClustalW stores them (n x n matrix in a double array).</summary>
    private static double[,] ClustalWPairwiseDistanceMatrix(string[] raw, ClustalWOptions options)
    {
        int n = raw.Length;
        var d = new double[n, n];
        int longest = raw.Max(s => s.Length);
        var quick = options.Quick ?? new ClustalWQuickOptions();
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                double v = options.PairwiseMode == ClustalWPairwiseMode.Quick
                    ? ClustalWQuickPairwiseDistance(raw[i], raw[j], quick, longest)
                    : ClustalWSinglePrecisionDistance(raw[i], raw[j], options.PairwiseGapOpen, options.PairwiseGapExtend, options.TransitionWeight);
                d[i, j] = v;
                d[j, i] = v;
            }
        }

        return d;
    }

    /// <summary>Everything the ClustalW multiple-alignment stage reads from <c>UserParameters</c>, for DNA and PROTEIN alike.</summary>
    private sealed record CwConfig(
        bool Protein,
        double GapOpen,
        double GapExtend,
        int DivergenceCutoff,
        bool NoWeights,
        bool NegativeMatrix,
        float TransitionWeight,
        ClustalWProfileMatrix Matrix,
        bool NoPreferredGapPenalties,
        bool NoHydrophilicGapPenalties,
        bool VariableRegionGapPenalties,
        int GapSeparationDistance,
        string HydrophilicResidues,
        bool UseEndGaps,
        bool EndGapPenalties,
        IReadOnlyList<int>? GapDistanceOverrides = null,
        ClustalWUserMatrix? UserMatrix = null,
        IReadOnlyList<ClustalWMatrixSeriesEntry>? UserMatrixSeries = null)
    {
        public static CwConfig FromDna(ClustalWOptions o) => new(
            false, o.GapOpen, o.GapExtend, o.DivergenceCutoff, o.NoWeights, o.NegativeMatrix, o.TransitionWeight,
            ClustalWProfileMatrix.Gonnet, false, false, false, 4, "GPSNDQEKR", false, false);

        public static CwConfig FromProtein(ClustalWProteinOptions o) => new(
            true, o.GapOpen, o.GapExtend, o.DivergenceCutoff, o.NoWeights, o.NegativeMatrix, 0.5f,
            o.Matrix, o.NoPreferredGapPenalties, o.NoHydrophilicGapPenalties, o.VariableRegionGapPenalties,
            o.GapSeparationDistance, o.HydrophilicResidues, false, false, o.UninitialisedGapDistances, o.UserMatrix, o.UserMatrixSeries);
    }

    /// <summary>Single-run state of the ClustalW multiple alignment (the members of Clustal / MSA / MyersMillerProfileAlign it needs).</summary>
    private sealed class ClustalWMsaState
    {
        private readonly int _n;
        private readonly CwConfig _o;
        private readonly int[][] _seq;          // residue codes of sequence i (0-based; replaced by the aligned row after a profile step)
        private int[,] _matrix = new int[CwNumRes, CwNumRes];
        private int _gapCoef;                   // DNA: fixed; PROTEIN: recomputed by every profileAlign (matrix average, profile lengths)
        private int _lenCoef;
        private double[,] _distMat = new double[0, 0];   // the similarities of Tree::calcSimilarities (1-based), used by the protein matrix series
        private int[] _alnWeight = Array.Empty<int>();   // Alignment::seqWeight
        private int[] _seqWeight = Array.Empty<int>();

        // per profileAlign call
        private int[][] _profile1 = Array.Empty<int[]>();
        private int[][] _profile2 = Array.Empty<int[]>();
        private int _prfLength1, _prfLength2;
        private int[] _displ = Array.Empty<int>();
        private int[] _hh = Array.Empty<int>(), _dd = Array.Empty<int>(), _rr = Array.Empty<int>(), _ss = Array.Empty<int>(), _gS = Array.Empty<int>();
        private int _printPtr, _lastPrint;
        private int _uninitialisedUses;   // profiles that used the unassigned gdist so far (ClustalWProteinOptions.UninitialisedGapDistances)
        private readonly Stack<int> _freedGapDistances = new(new[] { 0, 0 });   // gdist left in the two most recently freed profile objects (fresh heap memory: 0)

        public ClustalWMsaState(int[][] codes, CwConfig options)
        {
            _o = options;
            _n = codes.Length;
            _seq = codes;
            if (!options.Protein)
            {
                BuildMatrix();
                // gapcoef = (int)(100.0 * getGapOpen() * scale); getGapOpen() returns a float, scale (1.0 for IUB) a double.
                _gapCoef = (int)(100.0 * (double)(float)options.GapOpen * 1.0);
                _lenCoef = (int)(100.0 * (double)(float)options.GapExtend * 1.0);
            }
        }

        // SubMatrix::getProfileAlignMatrix, DNA / IUB branch (getMatrix with intScale 100, !minimise).
        private void BuildMatrix()
        {
            var raw = BuildClustalWIubMatrixRaw();
            int n = ClustalWResidueCodes.Length;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    _matrix[i, j] = raw[i, j];
            }

            // min over matrix[i][j], i = 0..maxAA, j = 1..i (sic) of the whole 32 x 32 array (zeros of unused codes included)
            int min = _matrix[0, 0];
            for (int i = 0; i <= CwMaxAa; i++)
            {
                for (int j = 1; j <= i; j++)
                    min = Math.Min(min, _matrix[i, j]);
            }

            if (!_o.NegativeMatrix && min < 0)
            {
                // every pair of residue codes of the nucleic-acid order gets -min added
                var valid = new bool[CwNumRes];
                foreach (char c in ClustalWNucleicOrder)
                    valid[ClustalWResidueCodes.IndexOf(c)] = true;
                for (int i = 0; i <= CwMaxAa; i++)
                {
                    if (!valid[i])
                        continue;
                    for (int j = 0; j <= CwMaxAa; j++)
                    {
                        if (valid[j])
                            _matrix[i, j] -= min;
                    }
                }
            }

            // gap rows / columns (grScore = ggScore = 0) are already zero.

            // transition weight: int(transitionWeight * matrix[0][0]) for A-G, C-T, C-U (defaultDNAXref of 'A','G' / 'C','T' / 'C','U')
            int transition = (int)(_o.TransitionWeight * _matrix[0, 0]);
            int a = ClustalWResidueCodes.IndexOf('A'), g = ClustalWResidueCodes.IndexOf('G'), c2 = ClustalWResidueCodes.IndexOf('C');
            int t = ClustalWResidueCodes.IndexOf('T'), u = ClustalWResidueCodes.IndexOf('U');
            foreach (var (x, y) in new[] { (a, g), (c2, t), (c2, u) })
            {
                _matrix[x, y] = transition;
                _matrix[y, x] = transition;
            }
        }

        public string[] Run(double[,] dist)
        {
            int n = _n;

            // ---- STEP 1 (pairwise distances, 0-based, computed by the caller), STEP 2: guide tree, weights, similarities, alignment steps
            var sim = new double[n + 1, n + 1];                                // 1-based, as DistMatrix
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    sim[i + 1, j + 1] = dist[i, j];
            }

            var seqWeight = new int[n + 1];
            var steps = BuildWeightsAndSteps(dist, sim, seqWeight);   // fills seqWeight (0-based, n entries) and converts sim to similarities
            _seqWeight = seqWeight;
            _alnWeight = new int[n + 1];

            // ---- STEP 3: MSA::multiSeqAlign
            _distMat = sim;
            MultiSeqAlign(sim, steps);

            var rows = new string[n];
            for (int i = 0; i < n; i++)
                rows[i] = CodesToRow(_seq[i]);
            return rows;
        }

        private static string CodesToRow(int[] row)
        {
            var sb = new System.Text.StringBuilder(row.Length);
            foreach (int c in row)
                sb.Append(c is >= 0 and <= CwMaxAa ? ClustalWResidueCodes[c] : '-');
            return sb.ToString();
        }

        // ------------------------------------------------------------------------------------------------------------
        // Tree.cpp: weights, similarities, alignment steps
        // ------------------------------------------------------------------------------------------------------------

        private List<int[]> BuildWeightsAndSteps(double[,] dist, double[,] sim, int[] seqWeight)
        {
            int n = _n;
            CwNode root;
            CwNode[] olptr;
            if (n == 2)
            {
                // UnRootedClusterTree::treeFromDistMatrix: "(a:d,b:d);" with d = distMat(1,2) / 2.0 printed with setprecision(5)
                // (general format, 5 significant digits) and read back as float; the tree is already rooted (no reRoot).
                double d = dist[0, 1] / 2.0;
                float printed = float.Parse(d.ToString("G5", CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);
                root = new CwNode();
                var l0 = new CwNode { Leaf = true, LeafIndex = 0, Parent = root, Dist = printed };
                var l1 = new CwNode { Leaf = true, LeafIndex = 1, Parent = root, Dist = printed };
                root.Left = l0;
                root.Right = l1;
                olptr = new[] { l0, l1 };
            }
            else
            {
                var unrooted = ClustalWNeighborJoining(n, dist);
                var rooted = ClustalWRootNodes(unrooted, n);
                root = rooted.Root;
                olptr = rooted.Leaves;
            }

            // Tree::orderNodes
            foreach (var leaf in olptr)
            {
                for (var p = leaf; p != null; p = p.Parent)
                    p.Order++;
            }

            CalcSeqWeights(olptr, seqWeight);
            CalcSimilarities(olptr, sim);
            return CreateSets(root, olptr);
        }

        // Tree::calcSeqWeights (distanceTree = true)
        private void CalcSeqWeights(CwNode[] olptr, int[] sweight)
        {
            int n = _n;
            if (n >= 2 && !_o.NoWeights)
            {
                var weight = new int[n];
                for (int i = 0; i < n; i++)
                {
                    float w = 0.0f;
                    for (var p = olptr[i]; p.Parent != null; p = p.Parent)
                        w += p.Dist / p.Order;
                    w = (float)((double)w * 100.0);
                    weight[i] = (int)w;
                }

                int sum = 0;
                for (int i = 0; i < n; i++)
                    sum += weight[i];
                if (sum == 0)
                {
                    for (int i = 0; i < n; i++)
                        weight[i] = 1;
                    sum = n;
                }

                for (int i = 0; i < n; i++)
                {
                    sweight[i] = weight[i] * CwIntScaleFactor / sum;
                    if (sweight[i] < 1)
                        sweight[i] = 1;
                }
            }
            else
            {
                int temp = CwIntScaleFactor / n;
                if (temp < 1)
                    temp = 1;
                for (int i = 0; i < n; i++)
                    sweight[i] = temp;
            }
        }

        // Tree::calcSimilarities: distMat := 100 - 100 * (tree path length), distances clamped to [0.01, 1]
        private void CalcSimilarities(CwNode[] olptr, double[,] sim)
        {
            int n = _n;
            var pathToRoot = new CwNode[n + 1];
            var distToNode = new float[n + 1];
            var dmat = new double[n, n];

            for (int i = 0; i < n; i++)
            {
                int depth = 0;
                float dist = 0.0f;
                for (var p = olptr[i]; p != null; p = p.Parent)
                {
                    pathToRoot[depth] = p;
                    dist += p.Dist;
                    distToNode[depth] = dist;
                    depth++;
                }

                for (int j = 0; j < i; j++)
                {
                    var p = olptr[j];
                    dist = 0.0f;
                    bool found = false;
                    int k = 0;
                    while (!found && p.Parent != null)
                    {
                        for (int q = 0; q < depth; q++)
                        {
                            if (p.Parent == pathToRoot[q])
                            {
                                found = true;
                                k = q;
                            }
                        }

                        dist += p.Dist;
                        p = p.Parent;
                    }

                    dmat[i, j] = dist + distToNode[k - 1];
                }
            }

            for (int i = 0; i < n; i++)
            {
                dmat[i, i] = 0.0;
                for (int j = 0; j < i; j++)
                {
                    if (dmat[i, j] < 0.01)
                        dmat[i, j] = 0.01;
                    if (dmat[i, j] > 1.0)
                        dmat[i, j] = 1.0;
                }
            }

            for (int i = 0; i < n; i++)
            {
                sim[i + 1, i + 1] = 0.0;
                for (int j = 0; j < i; j++)
                {
                    double value = 100.0 - dmat[i, j] * 100.0;
                    sim[i + 1, j + 1] = value;
                    sim[j + 1, i + 1] = value;
                }
            }
        }

        // Tree::createSets / groupSeqs: one step per internal node in post-order; left subtree = 1, right subtree = 2 (1-based vector)
        private List<int[]> CreateSets(CwNode root, CwNode[] olptr)
        {
            var steps = new List<int[]>();
            var leafOf = new Dictionary<CwNode, int>();
            for (int i = 0; i < olptr.Length; i++)
                leafOf[olptr[i]] = i;

            void Mark(CwNode p, int label, int[] groups)
            {
                if (p.Left == null && p.Right == null)
                {
                    groups[leafOf[p] + 1] = label;
                    return;
                }

                if (p.Left != null)
                    Mark(p.Left, label, groups);
                if (p.Right != null)
                    Mark(p.Right, label, groups);
            }

            void Visit(CwNode p)
            {
                if (p.Left == null && p.Right == null)
                    return;
                if (p.Left != null)
                    Visit(p.Left);
                if (p.Right != null)
                    Visit(p.Right);
                var groups = new int[_n + 1];
                if (p.Left != null)
                    Mark(p.Left, 1, groups);
                if (p.Right != null)
                {
                    Mark(p.Right, 2, groups);
                    steps.Add(groups);
                }
            }

            Visit(root);
            return steps;
        }

        // ------------------------------------------------------------------------------------------------------------
        // MSA.cpp multiSeqAlign
        // ------------------------------------------------------------------------------------------------------------

        private void MultiSeqAlign(double[,] distMat, List<int[]> steps)
        {
            int n = _n;
            Array.Copy(_seqWeight, _alnWeight, n + 1);   // alnPtr->addSeqWeight(seqWeight)

            var maxid = new int[n + 1];
            for (int i = 1; i <= n; i++)
            {
                maxid[i] = -1;
                for (int j = 1; j <= n; j++)
                {
                    if (j != i && maxid[i] < distMat[i, j])
                        maxid[i] = (int)distMat[i, j];
                }
            }

            var aligned = new int[n + 1];
            int ix = 0;

            foreach (var set in steps)
            {
                int entries = 0;
                for (int i = 1; i <= n; i++)
                {
                    if (set[i] != 0 && maxid[i] > _o.DivergenceCutoff)
                    {
                        entries++;
                        if (aligned[i] == 0)
                        {
                            ix++;
                            aligned[i] = 1;
                        }
                    }
                }

                if (entries > 0)
                    ProfileAlign(set, aligned);
            }

            // second pass: align the remaining (delayed, more divergent) sequences, closest to the aligned ones first
            var treeWeight = new int[n];
            for (int i = 0; i < n; i++)
                treeWeight[i] = _seqWeight[i];

            if (ix == 0)
            {
                int max = -1, iseq0 = 0;
                for (int i = 1; i <= n; i++)
                {
                    for (int j = i + 1; j <= n; j++)
                    {
                        if (max < distMat[i, j])
                        {
                            max = (int)distMat[i, j];
                            iseq0 = i;
                        }
                    }
                }

                aligned[iseq0] = 1;
                ix++;
            }

            var group = new int[n + 1];
            while (ix < n)
            {
                for (int i = 1; i <= n; i++)
                {
                    if (aligned[i] == 0)
                    {
                        maxid[i] = -1;
                        for (int j = 1; j <= n; j++)
                        {
                            if (maxid[i] < distMat[i, j] && aligned[j] != 0)
                                maxid[i] = (int)distMat[i, j];
                        }
                    }
                }

                int best = -1, iseq = 0;
                for (int i = 1; i <= n; i++)
                {
                    if (aligned[i] == 0 && maxid[i] > best)
                    {
                        best = maxid[i];
                        iseq = i;
                    }
                }

                // weight sequences with percent identity with the new sequence (tree weight x identity), renormalised
                if (!_o.NoWeights)
                {
                    for (int j = 0; j < n; j++)
                    {
                        if (aligned[j + 1] != 0)
                            _seqWeight[j] = (int)(treeWeight[j] * distMat[j + 1, iseq]);
                    }

                    int sum = 0;
                    for (int j = 0; j < n; j++)
                    {
                        if (aligned[j + 1] != 0)
                            sum += _seqWeight[j];
                    }

                    if (sum == 0)
                    {
                        for (int j = 0; j < n; j++)
                            _seqWeight[j] = 1;
                        sum = n;
                    }

                    for (int j = 0; j < n; j++)
                    {
                        if (aligned[j + 1] != 0)
                        {
                            _seqWeight[j] = _seqWeight[j] * CwIntScaleFactor / sum;
                            if (_seqWeight[j] < 1)
                                _seqWeight[j] = 1;
                        }
                    }
                }

                for (int j = 1; j <= n; j++)
                {
                    if (aligned[j] != 0)
                        group[j] = 1;
                    else if (iseq == j)
                        group[j] = 2;
                }

                Array.Copy(_seqWeight, _alnWeight, n + 1);   // alnPtr->addSeqWeight(seqWeight)
                aligned[iseq] = 1;
                ProfileAlign(group, aligned);
                ix++;
            }
        }

        // ------------------------------------------------------------------------------------------------------------
        // MyersMillerProfileAlign::profileAlign
        // ------------------------------------------------------------------------------------------------------------

        private void ProfileAlign(int[] group, int[] aligned)
        {
            int numSeq = _n;
            for (int i = 0; i < numSeq; i++)
            {
                if (aligned[i + 1] == 0)
                    group[i + 1] = 0;
            }

            int nseqs1 = 0, nseqs2 = 0;
            for (int i = 0; i < numSeq; i++)
            {
                if (group[i + 1] == 1)
                    nseqs1++;
                else if (group[i + 1] == 2)
                    nseqs2++;
            }

            if (nseqs1 == 0 || nseqs2 == 0)
                return;

            bool switchProfiles = nseqs2 > nseqs1;
            if (switchProfiles)
            {
                for (int i = 0; i < numSeq; i++)
                {
                    if (group[i + 1] == 1)
                        group[i + 1] = 2;
                    else if (group[i + 1] == 2)
                        group[i + 1] = 1;
                }
            }

            // profile 1 and profile 2 sequences (padded to the profile length with the gap code)
            var members1 = new List<int>();
            var members2 = new List<int>();
            for (int i = 0; i < numSeq; i++)
            {
                if (group[i + 1] == 1)
                    members1.Add(i);
                else if (group[i + 1] == 2)
                    members2.Add(i);
            }

            nseqs1 = members1.Count;
            nseqs2 = members2.Count;
            _prfLength1 = members1.Max(i => _seq[i].Length);
            _prfLength2 = members2.Max(i => _seq[i].Length);

            // mean percent identity between the two groups (the similarities of Tree::calcSimilarities), used by the protein matrix series
            double pcid = 0.0;
            int pairCount = 0;
            foreach (int a in members1)
            {
                foreach (int b in members2)
                {
                    pairCount++;
                    pcid += _distMat[a + 1, b + 1];
                }
            }

            pcid /= (double)(float)pairCount;

            var seqArray = new int[nseqs1 + nseqs2][];
            var weights = new int[nseqs1 + nseqs2];
            for (int k = 0; k < nseqs1; k++)
            {
                seqArray[k] = PadRow(_seq[members1[k]], _prfLength1);
                weights[k] = _alnWeight[members1[k]];
            }

            for (int k = 0; k < nseqs2; k++)
            {
                seqArray[nseqs1 + k] = PadRow(_seq[members2[k]], _prfLength2);
                weights[nseqs1 + k] = _alnWeight[members2[k]];
            }

            if (_o.Protein)
                ConfigureProteinProfile(seqArray, nseqs1, nseqs2, pcid);

            // The two profile objects (ProfileWithSub, ProfileStandard) are heap-allocated per call and freed at its end: their member gdist, which ProfileBase never
            // initialises, keeps the value the previous occupant of the chunk left (the allocator hands out the most recently freed chunk first).
            int gdist1 = _freedGapDistances.Pop(), gdist2 = _freedGapDistances.Pop();

            // profiles: with substitution matrix (profile 1), standard (profile 2)
            var gaps = new int[_prfLength1 + _prfLength2 + 3];   // alnPtr->getMaxAlnLength() + 1
            _profile1 = NewProfile(_prfLength1);
            CalcGapCoeff(_profile1, seqArray, 0, nseqs1, _prfLength1, gaps, ref gdist1);
            CalcProfileWithSub(_profile1, seqArray, 0, nseqs1, _prfLength1, gaps, weights);

            _profile2 = NewProfile(_prfLength2);
            CalcGapCoeff(_profile2, seqArray, nseqs1, nseqs1 + nseqs2, _prfLength2, gaps, ref gdist2);
            CalcStandardProfile(_profile2, seqArray, nseqs1, nseqs1 + nseqs2, _prfLength2, weights);

            // align the profiles (Myers & Miller)
            int maxAln = _prfLength1 + _prfLength2 + 2;
            _displ = new int[maxAln + 2];
            _hh = new int[maxAln + 2];
            _dd = new int[maxAln + 2];
            _rr = new int[maxAln + 2];
            _ss = new int[maxAln + 2];
            _gS = new int[maxAln + 2];
            _lastPrint = 0;
            _printPtr = 1;

            ProgDiff(0, 0, _prfLength1, _prfLength2, _profile1[0][CwGapCol], _profile1[_prfLength1][CwGapCol]);

            // progTracepath
            var alnPath1 = new int[maxAln + 2];
            var alnPath2 = new int[maxAln + 2];
            int pos = 0;
            int toDo = _printPtr - 1;
            for (int i = 1; i <= toDo; i++)
            {
                int d = _displ[i];
                if (d == 0)
                {
                    alnPath1[pos] = 2;
                    alnPath2[pos] = 2;
                    pos++;
                }
                else if (d > 0)
                {
                    for (int j = 0; j < d; j++)
                    {
                        alnPath2[pos + j] = 2;
                        alnPath1[pos + j] = 1;
                    }

                    pos += d;
                }
                else
                {
                    int k = -d;
                    for (int j = 0; j < k; j++)
                    {
                        alnPath1[pos + j] = 2;
                        alnPath2[pos + j] = 1;
                    }

                    pos += k;
                }
            }

            int alignmentLength = pos;

            // addGGaps: insert gap codes at the path's gap columns
            for (int j = 0; j < nseqs1 + nseqs2; j++)
            {
                int[] path = j < nseqs1 ? alnPath1 : alnPath2;
                int[] src = seqArray[j];
                var ta = new int[alignmentLength];
                int ixRes = 0;
                for (int i = 0; i < alignmentLength; i++)
                {
                    if (path[i] == 2)
                    {
                        ta[i] = ixRes < src.Length ? src[ixRes] : CwGapPos1;
                        ixRes++;
                    }
                    else
                    {
                        ta[i] = CwGapPos1;
                    }
                }

                seqArray[j] = ta;
            }

            // store the new rows: group 1 members first, then group 2 (profileAlign's final loops)
            int seqNum = 0;
            for (int k = 0; k < nseqs1; k++)
                _seq[members1[k]] = seqArray[seqNum++];
            for (int k = 0; k < nseqs2; k++)
                _seq[members2[k]] = seqArray[seqNum++];

            // delete profileWithSub; delete profileStandard;
            _freedGapDistances.Push(gdist1);
            _freedGapDistances.Push(gdist2);
        }

        // MyersMillerProfileAlign::profileAlign (PROTEIN part) + SubMatrix::getProfileAlignMatrix (protein branch): the matrix of the series chosen by the
        // percent identity of the two groups and the gap coefficients scaled by the matrix average score and the (logarithmic) profile lengths.
        private void ConfigureProteinProfile(int[][] seqArray, int nseqs1, int nseqs2, double pcid)
        {
            int len1 = 0;
            for (int i = 0; i < nseqs1; i++)
                len1 += CountResidues(seqArray[i], _prfLength1);
            len1 = (int)(len1 / (float)nseqs1);

            int len2 = 0;
            for (int i = nseqs1; i < nseqs1 + nseqs2; i++)
                len2 += CountResidues(seqArray[i], _prfLength2);
            len2 = (int)(len2 / (float)nseqs2);

            int minLen = Math.Min(len1, len2);
            var (cells, cross, intScale, scale) = ClustalWChooseProfileMatrix(_o.Matrix, _o.NegativeMatrix, _o.UserMatrix, _o.UserMatrixSeries, pcid, minLen);
            var (matrix, matAvgScore) = BuildClustalWProfileMatrix(cells, cross, intScale, _o.NegativeMatrix);
            _matrix = matrix;

            double logmin, logdiff;
            if (len1 == 0 || len2 == 0)
            {
                logmin = 1.0;
                logdiff = 1.0;
            }
            else
            {
                logmin = 1.0 / Math.Log10((double)minLen);
                if (len2 < len1)
                    logdiff = 1.0 + 0.5 * Math.Log10((double)((float)len2 / (float)len1));
                else if (len1 < len2)
                    logdiff = 1.0 + 0.5 * Math.Log10((double)((float)len1 / (float)len2));
                else
                    logdiff = 1.0;
                if (logdiff < 0.9)
                    logdiff = 0.9;
            }

            float gapOpen = (float)_o.GapOpen, gapExtend = (float)_o.GapExtend;
            if (_o.NegativeMatrix)
            {
                _gapCoef = ClustalWCastToInt(100.0 * gapOpen);
                _lenCoef = ClustalWCastToInt(100.0 * gapExtend);
            }
            else
            {
                if (matAvgScore <= 0)
                    _gapCoef = ClustalWCastToInt(100.0 * (float)(gapOpen + logmin));
                else
                    _gapCoef = ClustalWCastToInt(scale * matAvgScore * (float)(gapOpen / (logdiff * logmin)));   // PrfScaleValues.scale is a float: all single precision
                _lenCoef = ClustalWCastToInt(100.0 * gapExtend);
            }
        }

        // the residues (not the gap codes 30 / 31) among the first prfLength cells of a profile row
        private static int CountResidues(int[] row, int prfLength)
        {
            int n = 0;
            for (int j = 0; j < Math.Min(row.Length, prfLength); j++)
            {
                if (row[j] != CwGapPos1 && row[j] != CwGapPos2)
                    n++;
            }

            return n;
        }

        private static int[] PadRow(int[] row, int length)
        {
            var padded = new int[length];
            Array.Copy(row, padded, row.Length);
            for (int j = row.Length; j < length; j++)
                padded[j] = CwGapPos1;
            return padded;
        }

        private static int[][] NewProfile(int prfLength)
        {
            var p = new int[prfLength + 2][];
            for (int i = 0; i < p.Length; i++)
                p[i] = new int[CwLenCol + 2];
            return p;
        }

        // ProfileBase::calcGapCoeff (no gap-penalty mask: no secondary-structure input). DNA: position-specific gap penalties only; PROTEIN additionally the
        // variable-region, residue-specific (Pascarella & Argos), hydrophilic-run penalties and the gap-separation distance.
        private void CalcGapCoeff(int[][] profile, int[][] seqArray, int firstSeq, int lastSeq, int prfLength, int[] gaps, ref int gdist)
        {
            int numSeq = lastSeq - firstSeq;
            bool protein = _o.Protein;
            int pcid = 0;
            if (numSeq == 2)
                pcid = (int)ClustalWProfilePercentId(seqArray[firstSeq], seqArray[firstSeq + 1], prfLength);

            for (int j = 0; j < prfLength; j++)
                gaps[j] = 0;

            // ProfileBase.h members nVarPen / nHydPen / nPrefPen are "this penalty is switched OFF"
            bool nVarPen, nHydPen, nPrefPen;
            if (_o.VariableRegionGapPenalties && pcid > 60)
            {
                // gdist is not assigned on this path: the member of the C++ object is uninitialised memory (see ProfileAlign)
                nHydPen = nPrefPen = true;
                nVarPen = false;
                if (_o.GapDistanceOverrides is { } forced && _uninitialisedUses < forced.Count)
                    gdist = forced[_uninitialisedUses];
                _uninitialisedUses++;
            }
            else
            {
                nVarPen = true;
                nHydPen = _o.NoHydrophilicGapPenalties;
                nPrefPen = _o.NoPreferredGapPenalties;
                gdist = _o.GapSeparationDistance;
            }

            for (int i = firstSeq; i < lastSeq; i++)
            {
                // end gaps are not counted as gaps (useEndGaps == false && endGapPenalties == false)
                int isx = 0, ie = prfLength;
                if (!_o.UseEndGaps && !_o.EndGapPenalties)
                {
                    for (int j = 0; j < prfLength; j++)
                    {
                        int c = seqArray[i][j];
                        if (c < 0 || c > CwMaxAa)
                            isx++;
                        else
                            break;
                    }

                    for (int j = prfLength - 1; j >= 0; j--)
                    {
                        int c = seqArray[i][j];
                        if (c < 0 || c > CwMaxAa)
                            ie--;
                        else
                            break;
                    }
                }

                for (int j = isx; j < ie; j++)
                {
                    if (seqArray[i][j] < 0 || seqArray[i][j] > CwMaxAa)
                        gaps[j]++;
                }
            }

            int[] vWeight = Array.Empty<int>(), resWeight = Array.Empty<int>(), hydWeight = Array.Empty<int>();
            if (protein && !nVarPen)
                vWeight = CalcVariableRegionPenalties(seqArray, firstSeq, prfLength);
            if (protein && !nPrefPen)
                resWeight = CalcResidueSpecificPenalties(seqArray, firstSeq, lastSeq, prfLength);
            if (protein && !nHydPen)
                hydWeight = CalcHydrophilicPenalties(seqArray, firstSeq, lastSeq, prfLength);

            // mark the residues close to an existing gap (gapPos < 0)
            var gapPos = new int[prfLength + 2];
            int G(int idx) => idx >= 0 && idx < gaps.Length ? gaps[idx] : 0;
            if (!protein || gdist <= 0)
            {
                for (int i = 0; i < prfLength; i++)
                    gapPos[i] = gaps[i];
            }
            else
            {
                int i = 0;
                while (i < prfLength)
                {
                    if (gaps[i] <= 0)
                    {
                        gapPos[i] = gaps[i];
                        i++;
                    }
                    else
                    {
                        for (int j = -gdist + 1; j < 0; j++)
                        {
                            if (i + j >= 0 && i + j < prfLength && (gaps[i + j] == 0 || gaps[i + j] < j))
                                gapPos[i + j] = j;
                        }

                        while (G(i) > 0)
                        {
                            if (i >= prfLength)
                                break;
                            gapPos[i] = gaps[i];
                            i++;
                        }

                        int jj;
                        for (jj = 0; jj < gdist; jj++)
                        {
                            if (G(i + jj) > 0)
                                break;
                            if (i + jj >= 0 && i + jj < prfLength && (gaps[i + jj] == 0 || gaps[i + jj] < -jj))
                                gapPos[i + jj] = -jj - 1;
                        }

                        i += jj;
                    }
                }
            }

            for (int j = 0; j < prfLength; j++)
            {
                if (gapPos[j] <= 0)
                {
                    // residue-specific and hydrophilic gap penalties (PROTEIN); DNA keeps the plain penalty
                    profile[j + 1][CwGapCol] = protein ? ClustalWLocalPenalty(_gapCoef, j, nVarPen, nHydPen, nPrefPen, resWeight, hydWeight, vWeight) : _gapCoef;
                    profile[j + 1][CwLenCol] = _lenCoef;

                    // increase the gap penalty near an existing gap
                    if (gapPos[j] < 0)
                        profile[j + 1][CwGapCol] = (int)(profile[j + 1][CwGapCol] * (2.0 + 2.0 * (gdist + gapPos[j]) / gdist));
                }
                else
                {
                    float scale = ((float)(numSeq - gaps[j]) / (float)numSeq) * CwReducedGap;
                    profile[j + 1][CwGapCol] = (int)(scale * _gapCoef);
                    profile[j + 1][CwLenCol] = (int)(0.5 * _lenCoef);
                }

                // make sure no penalty is zero - even for all-gap positions
                if (profile[j + 1][CwGapCol] <= 0)
                    profile[j + 1][CwGapCol] = 1;
                if (profile[j + 1][CwLenCol] <= 0)
                    profile[j + 1][CwLenCol] = 1;
            }

            // set the penalties at the beginning and end of the profile
            if (_o.EndGapPenalties)
            {
                profile[0][CwGapCol] = _gapCoef;
                profile[0][CwLenCol] = _lenCoef;
            }
            else
            {
                profile[0][CwGapCol] = 0;
                profile[0][CwLenCol] = 0;
                profile[prfLength][CwGapCol] = 0;
                profile[prfLength][CwLenCol] = 0;
            }
        }

        // ProfileBase::percentId (note: residue codes < maxAA, i.e. Z is not counted)
        private static float ClustalWProfilePercentId(int[] s1, int[] s2, int prfLength)
        {
            int count = 0, total = 0;
            for (int i = 0; i < prfLength; i++)
            {
                if (s1[i] >= 0 && s1[i] < CwMaxAa)
                {
                    total++;
                    if (s1[i] == s2[i])
                        count++;
                }
            }

            return total == 0 ? 0f : (float)(100.0 * count / (float)total);
        }

        // ProfileBase::calcVPenalties: variable-zone weight from the identity of the first two sequences in a window of 2 * vwindow columns (vlut = identity)
        private static int[] CalcVariableRegionPenalties(int[][] aln, int firstSeq, int prfLength)
        {
            const int vwindow = 5, vll = 50;
            var weight = new int[prfLength + 2];
            for (int i = 0; i < prfLength; i++)
            {
                weight[i] = 0;
                int t = 0;
                for (int j = i - vwindow; j < i + vwindow; j++)
                {
                    if (j >= 0 && j < prfLength)
                    {
                        int ix1 = aln[firstSeq][j], ix2 = aln[firstSeq + 1][j];
                        if (ix1 < 0 || ix1 > CwMaxAa || ix2 < 0 || ix2 > CwMaxAa)
                            continue;
                        if (ClustalWResidueCodes[ix1] == ClustalWResidueCodes[ix2])
                            weight[i] += 1;
                        t++;
                    }
                }

                // now we have a weight -t < w < t
                weight[i] += t;
                weight[i] = t > 0 ? weight[i] * 100 / (2 * t) : 100;
                if (weight[i] < vll)
                    weight[i] = vll;
            }

            return weight;
        }

        // ProfileBase::calcResidueSpecificPen: Pascarella & Argos (1992) residue-specific gap modification factors (Table 1 of Thompson et al. 1994)
        private static int[] CalcResidueSpecificPenalties(int[][] aln, int firstSeq, int lastSeq, int prfLength)
        {
            const string pascarellaRes = "ACDEFGHKILMNPQRSTVYW";
            int[] pascarellaProb = { 87, 87, 104, 69, 80, 139, 100, 104, 68, 79, 71, 137, 126, 93, 128, 124, 111, 75, 100, 77 };
            int numSeq = lastSeq - firstSeq;
            var weight = new int[prfLength + 2];
            for (int i = 0; i < prfLength; i++)
            {
                weight[i] = 0;
                for (int k = firstSeq; k < lastSeq; k++)
                {
                    int ix = aln[k][i];
                    if (ix < 0 || ix > CwMaxAa)
                        continue;
                    int j = pascarellaRes.IndexOf(ClustalWResidueCodes[ix]);
                    if (j >= 0)
                        weight[i] += 180 - pascarellaProb[j];
                }

                weight[i] /= numSeq;
            }

            return weight;
        }

        // ProfileBase::calcHydrophilicPen: runs of more than three hydrophilic residues (-HGAPRESIDUES) in a sequence add 100 to every column of the run
        private int[] CalcHydrophilicPenalties(int[][] aln, int firstSeq, int lastSeq, int prfLength)
        {
            var weight = new int[prfLength + 2];
            var hyd = new int[prfLength + 2];
            string hydResidues = _o.HydrophilicResidues;
            for (int k = firstSeq; k < lastSeq; k++)
            {
                for (int i = 0; i < prfLength; i++)
                {
                    hyd[i] = 0;
                    int res = aln[k][i];
                    if (res < 0 || res > CwMaxAa)
                        continue;
                    if (hydResidues.Contains(ClustalWResidueCodes[res]))
                        hyd[i] = 1;
                }

                int p = 0;
                while (p < prfLength)
                {
                    if (hyd[p] == 0)
                    {
                        p++;
                    }
                    else
                    {
                        int s = p;
                        while (hyd[p] != 0 && p < prfLength)
                            p++;
                        int e = p;
                        if (e - s > 3)
                        {
                            for (int j = s; j < e; j++)
                                weight[j] += 100;
                        }
                    }
                }
            }

            float scale = lastSeq - firstSeq;
            for (int i = 0; i < prfLength; i++)
                weight[i] = (int)(weight[i] / scale);
            return weight;
        }

        // ProfileBase::localPenalty (PROTEIN): the gap-opening penalty of a column scaled by the variable-region, hydrophilic and residue-specific weights
        private static int ClustalWLocalPenalty(int penalty, int n, bool nVarPen, bool nHydPen, bool nPrefPen, int[] resWeight, int[] hydWeight, int[] vWeight)
        {
            bool h = false;
            float gw = 1.0f;
            if (!nVarPen)
                gw = (float)(gw * (vWeight[n] / 100.0));
            if (!nHydPen && hydWeight[n] > 0)
            {
                gw = (float)(gw * 0.5);
                h = true;
            }

            if (!nPrefPen && !h)
                gw = (float)(gw * (resWeight[n] / 100.0));
            gw *= penalty;
            return (int)gw;
        }

        // ProfileWithSub::calcProfileWithSub
        private void CalcProfileWithSub(int[][] profile, int[][] seqArray, int firstSeq, int lastSeq, int prfLength, int[] gaps, int[] seqWeight)
        {
            int numSeq = lastSeq - firstSeq;
            var weighting = new int[CwNumRes + 2][];
            for (int a = 0; a < weighting.Length; a++)
                weighting[a] = new int[prfLength + 2];

            int sum2 = 0;
            for (int s = firstSeq; s < lastSeq; s++)
                sum2 += seqWeight[s];

            for (int col = 0; col < prfLength; col++)
            {
                for (int aa = 0; aa <= CwMaxAa; aa++)
                {
                    weighting[aa][col] = 0;
                    for (int s = firstSeq; s < lastSeq; s++)
                    {
                        if (aa == seqArray[s][col])
                            weighting[aa][col] += seqWeight[s];
                    }
                }

                weighting[CwGapPos1][col] = 0;
                for (int s = firstSeq; s < lastSeq; s++)
                {
                    if (CwGapPos1 == seqArray[s][col])
                        weighting[CwGapPos1][col] += seqWeight[s];
                }

                weighting[CwGapPos2][col] = 0;
                for (int s = firstSeq; s < lastSeq; s++)
                {
                    if (CwGapPos2 == seqArray[s][col])
                        weighting[CwGapPos2][col] += seqWeight[s];
                }
            }

            for (int pos = 0; pos < prfLength; pos++)
            {
                if (gaps[pos] == numSeq)
                {
                    for (int res = 0; res <= CwMaxAa; res++)
                        profile[pos + 1][res] = _matrix[res, CwGapPos1];
                    profile[pos + 1][CwGapPos1] = _matrix[CwGapPos1, CwGapPos1];
                    profile[pos + 1][CwGapPos2] = _matrix[CwGapPos2, CwGapPos1];
                }
                else
                {
                    float scale = (float)(numSeq - gaps[pos]) / (float)numSeq;
                    for (int res = 0; res <= CwMaxAa; res++)
                    {
                        int f = 0;
                        for (int aa = 0; aa <= CwMaxAa; aa++)
                            f += weighting[aa][pos] * _matrix[aa, res];
                        f += weighting[CwGapPos1][pos] * _matrix[CwGapPos1, res];
                        f += weighting[CwGapPos2][pos] * _matrix[CwGapPos2, res];
                        profile[pos + 1][res] = (int)(((float)f / (float)sum2) * scale);
                    }

                    int f1 = 0;
                    for (int aa = 0; aa <= CwMaxAa; aa++)
                        f1 += weighting[aa][pos] * _matrix[aa, CwGapPos1];
                    f1 += weighting[CwGapPos1][pos] * _matrix[CwGapPos1, CwGapPos1];
                    f1 += weighting[CwGapPos2][pos] * _matrix[CwGapPos2, CwGapPos1];
                    profile[pos + 1][CwGapPos1] = (int)(((float)f1 / (float)sum2) * scale);

                    int f2 = 0;
                    for (int aa = 0; aa <= CwMaxAa; aa++)
                        f2 += weighting[aa][pos] * _matrix[aa, CwGapPos2];
                    f2 += weighting[CwGapPos1][pos] * _matrix[CwGapPos1, CwGapPos2];
                    f2 += weighting[CwGapPos2][pos] * _matrix[CwGapPos2, CwGapPos2];
                    profile[pos + 1][CwGapPos2] = (int)(((float)f2 / (float)sum2) * scale);
                }
            }
        }

        // ProfileStandard::calcStandardProfile
        private static void CalcStandardProfile(int[][] profile, int[][] seqArray, int firstSeq, int lastSeq, int prfLength, int[] seqWeight)
        {
            for (int r = 0; r < prfLength; r++)
            {
                int sum2 = 0;
                for (int i = firstSeq; i < lastSeq; i++)
                    sum2 += seqWeight[i];

                if (sum2 == 0)
                {
                    for (int d = 0; d <= CwMaxAa; d++)
                        profile[r + 1][d] = 0;
                    profile[r + 1][CwGapPos1] = 0;
                    profile[r + 1][CwGapPos2] = 0;
                }
                else
                {
                    for (int d = 0; d <= CwMaxAa; d++)
                    {
                        int sum1 = 0;
                        for (int i = firstSeq; i < lastSeq; i++)
                        {
                            if (d == seqArray[i][r])
                                sum1 += seqWeight[i];
                        }

                        profile[r + 1][d] = (int)(10 * (float)sum1 / (float)sum2);
                    }

                    int g1 = 0;
                    for (int i = firstSeq; i < lastSeq; i++)
                    {
                        if (CwGapPos1 == seqArray[i][r])
                            g1 += seqWeight[i];
                    }

                    profile[r + 1][CwGapPos1] = (int)(10 * (float)g1 / (float)sum2);
                    int g2 = 0;
                    for (int i = firstSeq; i < lastSeq; i++)
                    {
                        if (CwGapPos2 == seqArray[i][r])
                            g2 += seqWeight[i];
                    }

                    profile[r + 1][CwGapPos2] = (int)(10 * (float)g2 / (float)sum2);
                }
            }
        }

        // ---- penalty accessors (end gaps are never penalised: getEndGapPenalties() == false) ----

        private int PrfScore(int n, int m)
        {
            int score = 0;
            var p1 = _profile1[n];
            var p2 = _profile2[m];
            for (int ix = 0; ix <= CwMaxAa; ix++)
                score += p1[ix] * p2[ix];
            score += p1[CwGapPos1] * p2[CwGapPos1];
            score += p1[CwGapPos2] * p2[CwGapPos2];
            return score / 10;
        }

        private int OpenPenalty1(int i, int j) => !_o.EndGapPenalties && (i == 0 || i == _prfLength1) ? 0 : _profile2[j][CwGapCol] + _profile1[i][CwGapCol];

        private int ExtPenalty1(int i, int j) => !_o.EndGapPenalties && (i == 0 || i == _prfLength1) ? 0 : _profile2[j][CwLenCol];

        private int GapPenalty1(int i, int j, int k)
        {
            if (k <= 0)
                return 0;
            if (!_o.EndGapPenalties && (i == 0 || i == _prfLength1))
                return 0;
            int g = _profile2[j][CwGapCol] + _profile1[i][CwGapCol];
            int h = 0;
            for (int ix = 0; ix < k && ix + j < _prfLength2; ix++)
                h += _profile2[ix + j][CwLenCol];
            return g + h;
        }

        private int OpenPenalty2(int i, int j) => !_o.EndGapPenalties && (j == 0 || j == _prfLength2) ? 0 : _profile1[i][CwGapCol] + _profile2[j][CwGapCol];

        private int ExtPenalty2(int i, int j) => !_o.EndGapPenalties && (j == 0 || j == _prfLength2) ? 0 : _profile1[i][CwLenCol];

        private int GapPenalty2(int i, int j, int k)
        {
            if (k <= 0)
                return 0;
            if (!_o.EndGapPenalties && (j == 0 || j == _prfLength2))
                return 0;
            int g = _profile1[i][CwGapCol] + _profile2[j][CwGapCol];
            int h = 0;
            for (int ix = 0; ix < k && ix + i < _prfLength1; ix++)
                h += _profile1[ix + i][CwLenCol];
            return g + h;
        }

        private void ProgDel(int k)
        {
            if (_lastPrint < 0)
            {
                _displ[_printPtr - 1] -= k;
                _lastPrint = _displ[_printPtr - 1];
            }
            else
            {
                _displ[_printPtr++] = -k;
                _lastPrint = -k;
            }
        }

        private void ProgAdd(int k)
        {
            if (_lastPrint < 0)
            {
                _displ[_printPtr - 1] = k;
                _displ[_printPtr++] = _lastPrint;
            }
            else
            {
                _lastPrint = _displ[_printPtr++] = k;
            }
        }

        private void ProgAlign()
        {
            _displ[_printPtr++] = _lastPrint = 0;
        }

        // MyersMillerProfileAlign::progDiff (the function-level statics of the C++ are plain locals here: every use of
        // them is inside one invocation before the recursive calls)
        private void ProgDiff(int A, int B, int M, int N, int go1, int go2)
        {
            int midi, midj, type;
            int midh;
            int t, tl, g, h, hh, f, e, s;

            if (N <= 0)
            {
                if (M > 0)
                    ProgDel(M);
                return;
            }

            if (M <= 1)
            {
                if (M <= 0)
                {
                    ProgAdd(N);
                    return;
                }

                if (go1 == 0)
                    midh = -GapPenalty1(A + 1, B + 1, N);
                else
                    midh = -GapPenalty2(A + 1, B, 1) - GapPenalty1(A + 1, B + 1, N);

                midj = 0;
                for (int j = 1; j <= N; j++)
                {
                    hh = -GapPenalty1(A, B + 1, j - 1) + PrfScore(A + 1, B + j) - GapPenalty1(A + 1, B + j + 1, N - j);
                    if (hh > midh)
                    {
                        midh = hh;
                        midj = j;
                    }
                }

                if (midj == 0)
                {
                    ProgAdd(N);
                    ProgDel(1);
                }
                else
                {
                    if (midj > 1)
                        ProgAdd(midj - 1);
                    ProgAlign();
                    if (midj < N)
                        ProgAdd(N - midj);
                }

                return;
            }

            // divide sequence A in half
            midi = M / 2;

            // forward phase: HH[j], DD[j]
            _hh[0] = 0;
            t = -OpenPenalty1(A, B + 1);
            tl = -ExtPenalty1(A, B + 1);
            for (int j = 1; j <= N; j++)
            {
                t += tl;
                _hh[j] = t;
                _dd[j] = t - OpenPenalty2(A + 1, B + j);
            }

            if (go1 == 0)
                t = 0;
            else
                t = -OpenPenalty2(A + 1, B);
            tl = -ExtPenalty2(A + 1, B);
            for (int i = 1; i <= midi; i++)
            {
                s = _hh[0];
                t += tl;
                _hh[0] = hh = t;
                f = t - OpenPenalty1(A + i, B + 1);

                for (int j = 1; j <= N; j++)
                {
                    g = OpenPenalty1(A + i, B + j);
                    h = ExtPenalty1(A + i, B + j);
                    hh = hh - g - h;
                    f = f - h;
                    if (hh > f)
                        f = hh;
                    g = OpenPenalty2(A + i, B + j);
                    h = ExtPenalty2(A + i, B + j);
                    hh = _hh[j] - g - h;
                    e = _dd[j] - h;
                    if (hh > e)
                        e = hh;
                    hh = s + PrfScore(A + i, B + j);
                    if (f > hh)
                        hh = f;
                    if (e > hh)
                        hh = e;

                    s = _hh[j];
                    _hh[j] = hh;
                    _dd[j] = e;
                }
            }

            _dd[0] = _hh[0];

            // reverse phase: RR[j], SS[j]
            _rr[N] = 0;
            tl = 0;
            for (int j = N - 1; j >= 0; j--)
            {
                g = -OpenPenalty1(A + M, B + j + 1);
                tl -= ExtPenalty1(A + M, B + j + 1);
                _rr[j] = g + tl;
                _ss[j] = _rr[j] - OpenPenalty2(A + M, B + j);
                _gS[j] = OpenPenalty2(A + M, B + j);
            }

            tl = 0;
            for (int i = M - 1; i >= midi; i--)
            {
                s = _rr[N];
                if (go2 == 0)
                    g = 0;
                else
                    g = -OpenPenalty2(A + i + 1, B + N);
                tl -= ExtPenalty2(A + i + 1, B + N);
                _rr[N] = hh = g + tl;
                t = OpenPenalty1(A + i, B + N);
                f = _rr[N] - t;

                for (int j = N - 1; j >= 0; j--)
                {
                    g = OpenPenalty1(A + i, B + j + 1);
                    h = ExtPenalty1(A + i, B + j + 1);
                    hh = hh - g - h;
                    f = f - h - g + t;
                    if (hh > f)
                        f = hh;
                    t = g;
                    g = OpenPenalty2(A + i + 1, B + j);
                    h = ExtPenalty2(A + i + 1, B + j);
                    hh = _rr[j] - g - h;
                    if (i == M - 1)
                    {
                        e = _ss[j] - h;
                    }
                    else
                    {
                        e = _ss[j] - h - g + OpenPenalty2(A + i + 2, B + j);
                        _gS[j] = g;
                    }

                    if (hh > e)
                        e = hh;
                    hh = s + PrfScore(A + i + 1, B + j + 1);
                    if (f > hh)
                        hh = f;
                    if (e > hh)
                        hh = e;

                    s = _rr[j];
                    _rr[j] = hh;
                    _ss[j] = e;
                }
            }

            _ss[N] = _rr[N];
            _gS[N] = OpenPenalty2(A + midi + 1, B + N);

            // find midj such that HH[j] + RR[j] or DD[j] + SS[j] + gap is the maximum
            midh = _hh[0] + _rr[0];
            midj = 0;
            type = 1;
            for (int j = 0; j <= N; j++)
            {
                hh = _hh[j] + _rr[j];
                if (hh >= midh && (hh > midh || (_hh[j] != _dd[j] && _rr[j] == _ss[j])))
                {
                    midh = hh;
                    midj = j;
                }
            }

            for (int j = N; j >= 0; j--)
            {
                hh = _dd[j] + _ss[j] + _gS[j];
                if (hh > midh)
                {
                    midh = hh;
                    midj = j;
                    type = 2;
                }
            }

            // conquer recursively around the midpoint
            if (type == 1)
            {
                ProgDiff(A, B, midi, midj, go1, 1);
                ProgDiff(A + midi, B + midj, M - midi, N - midj, 1, go2);
            }
            else
            {
                ProgDiff(A, B, midi - 1, midj, go1, 0);
                ProgDel(2);
                ProgDiff(A + midi + 1, B + midj, M - midi - 1, N - midj, 0, go2);
            }
        }
    }

    #endregion
}
