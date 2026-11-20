namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// proba_pair: the DEFAULT primary-library method of T-Coffee 13.41 (audit A1-13b, F66; t_coffee.c L4375 `Mproba_pair //new default`).
//
// A port of `proba_pair_wise` (util_dp_suboptimal_nw.c; the pair-HMM of ProbCons 1.1 as T-Coffee embeds it) and `ProbaMatrix2CL`:
//   * the model: a match state and NumInsertStates pairs of insertion states (DNA: 1, protein: 2 = bi-phasic), initial distribution / transition
//     probabilities / emission tables of the source (single precision), log-space arithmetic with the ProbCons `LOOKUP` approximation of
//     log(1 + exp(-x)) (`LOG_PLUS_EQUALS`, `LOG_ADD`, util.c);
//   * a residue is looked up through the 20-letter "acdefghiklmnpqrstvwy" code of the dirichlet profile (`aln2prf`, `dirichlet_code2aa_lu`): a letter outside
//     those 20 (U, X, B, ...) has an all-zero profile, i.e. log-probability 0 for every emission it takes part in (a quirk of the program that is reproduced);
//     the DNA model is the program's (note: it is built with the PROTEIN initial distribution `initDistrib2Default`, first three values, and the six-letter
//     table "ACGUTN" with 1e-5 / 1e-10 for every other letter);
//   * forward and backward in float arithmetic (`forward_proba_pair_wise`, `backward_proba_pair_wise`), total probability (`ComputeTotalProbability`),
//     posterior of every cell in the match state v = EXP(min(0, f + b - total)) (double polynomial `EXP` of util.c), cells with v > 0.01 become library
//     pairs with weight (int)((float)v x 1000), the best 4 x min(I, J) of them by weight (qsort ascending then reversed: ties in decreasing (i, j) order).
// =====================================================================================================================

/// <summary>The sequence type the program decides (<c>get_string_type</c>): selects matrices, gap penalties and the pair-HMM.</summary>
internal enum TcSeqKind { Dna, Rna, Protein }

public static partial class SequenceAligner
{
    private const float TcLogZero = -200000004008175468544.0f;
    private const float TcLogUnderflowThreshold = 7.50f;

    private static float TcLookup(float x)
    {
        if (x <= 1.00f) return ((-0.009350833524763f * x + 0.130659527668286f) * x + 0.498799810682272f) * x + 0.693203116424741f;
        if (x <= 2.50f) return ((-0.014532321752540f * x + 0.139942324101744f) * x + 0.495635523139337f) * x + 0.692140569840976f;
        if (x <= 4.50f) return ((-0.004605031767994f * x + 0.063427417320019f) * x + 0.695956496475118f) * x + 0.514272634594009f;
        return ((-0.000458661602210f * x + 0.009695946122598f) * x + 0.930734667215156f) * x + 0.168037164329057f;
    }

    /// <summary><c>LOG_PLUS_EQUALS</c> (util.c).</summary>
    private static void TcLogPlusEquals(ref float x, float y)
    {
        if (x < y)
            x = (x == TcLogZero || y - x >= TcLogUnderflowThreshold) ? y : TcLookup(y - x) + x;
        else
            x = (y == TcLogZero || x - y >= TcLogUnderflowThreshold) ? x : TcLookup(x - y) + y;
    }

    /// <summary><c>LOG_ADD</c> (util.c).</summary>
    private static float TcLogAdd(float x, float y)
    {
        if (x < y) return (x == TcLogZero || y - x >= TcLogUnderflowThreshold) ? y : TcLookup(y - x) + x;
        return (y == TcLogZero || x - y >= TcLogUnderflowThreshold) ? x : TcLookup(x - y) + y;
    }

    /// <summary><c>EXP</c> (util.c): the piecewise polynomial approximation of exp(x) for x &lt;= 0 (double arithmetic).</summary>
    private static double TcExp(double x)
    {
        if (x > -2)
        {
            if (x > -0.5)
            {
                if (x > 0) return Math.Exp(x);
                return (((0.03254409303190190000 * x + 0.16280432765779600000) * x + 0.49929760485974900000) * x + 0.99995149601363700000) * x + 0.99999925508501600000;
            }

            if (x > -1)
                return (((0.01973899026052090000 * x + 0.13822379685007000000) * x + 0.48056651562365000000) * x + 0.99326940370383500000) * x + 0.99906756856399500000;
            return (((0.00940528203591384000 * x + 0.09414963667859410000) * x + 0.40825793595877300000) * x + 0.93933625499130400000) * x + 0.98369508190545300000;
        }

        if (x > -8)
        {
            if (x > -4)
                return (((0.00217245711583303000 * x + 0.03484829428350620000) * x + 0.22118199801337800000) * x + 0.67049462206469500000) * x + 0.83556950223398500000;
            return (((0.00012398771025456900 * x + 0.00349155785951272000) * x + 0.03727721426017900000) * x + 0.17974997741536900000) * x + 0.33249299994217400000;
        }

        if (x > -16)
            return (((0.00000051741713416603 * x + 0.00002721456879608080) * x + 0.00053418601865636800) * x + 0.00464101989351936000) * x + 0.01507447981459420000;
        return 0;
    }

    /// <summary>The pair-HMM of <c>proba_pair</c> for one sequence type (<c>ProbabilisticModel</c> with the program's default tables).</summary>
    private sealed class TcProbaModel
    {
        public int MatrixTypes, InsertStates;
        public float[] InitialDistribution = Array.Empty<float>();
        public float[,] TransProb = new float[0, 0];
        public readonly float[,] InsProb = new float[256, 5];
        public readonly float[,] MatchProb = new float[256, 256];

        private static readonly Lazy<TcProbaModel> DnaModel = new(() => new TcProbaModel(TcSeqKind.Dna));
        private static readonly Lazy<TcProbaModel> RnaModel = new(() => new TcProbaModel(TcSeqKind.Rna));
        private static readonly Lazy<TcProbaModel> ProteinModel = new(() => new TcProbaModel(TcSeqKind.Protein));

        private static readonly Dictionary<(TcSeqKind, string), TcProbaModel> MatrixModels = new();

        public static TcProbaModel For(TcSeqKind kind, string? matrix = null)
        {
            var m = TcResolveMatrix(matrix);
            if (m == null) return kind switch { TcSeqKind.Protein => ProteinModel.Value, TcSeqKind.Rna => RnaModel.Value, _ => DnaModel.Value };
            lock (MatrixModels)
            {
                if (!MatrixModels.TryGetValue((kind, matrix!), out var model)) MatrixModels[(kind, matrix!)] = model = new TcProbaModel(kind, m);
                return model;
            }
        }

        private TcProbaModel(TcSeqKind kind, TcMatrix? matrix = null)
        {
            // DNA: one insert state pair, 3 matrix types; RNA (`retrieve_seq_type() == "RNA"`: more than 98 % DNA letters and at least one U): bi-phasic like the
            // protein model (2 / 5) with RNAgapOpen2Default / RNAgapExtend2Default but the six-letter nucleotide emission table and the PROTEIN initial distribution
            // (the call passes initDistrib2Default, not RNAinitDistrib2Default)
            bool protein = kind == TcSeqKind.Protein;
            bool biphasic = kind != TcSeqKind.Dna;
            InsertStates = biphasic ? 2 : 1;
            MatrixTypes = biphasic ? 5 : 3;
            float[] gapOpen = TcProbaTables.DnaGapOpen, gapExtend = TcProbaTables.DnaGapExtend;
            if (protein)
            {
                gapOpen = TcProbaTables.GapOpen2;
                gapExtend = TcProbaTables.GapExtend2;
            }
            else if (kind == TcSeqKind.Rna)
            {
                gapOpen = TcProbaTables.RnaGapOpen2;
                gapExtend = TcProbaTables.RnaGapExtend2;
            }
            float[] initDistrib = TcProbaTables.InitDistrib2; // the DNA model uses the protein one as well (first MatrixTypes values)
            string alphabet = protein ? TcProbaTables.ProteinAlphabet : TcProbaTables.DnaAlphabet;
            float[] single = protein ? TcProbaTables.EmitSingle : TcProbaTables.DnaEmitSingle;
            float[] pairs = protein ? TcProbaTables.EmitPairs : TcProbaTables.DnaEmitPairs;
            int l = alphabet.Length;
            if (matrix != null) pairs = TcEmitPairsFromMatrix(matrix, alphabet, single);

            var emitSingle = new float[256];
            var emitPairs = new float[256, 256];
            for (int i = 0; i < 256; i++)
            {
                emitSingle[i] = protein ? 1 : 1e-5f;
                for (int j = 0; j < 256; j++) emitPairs[i, j] = protein ? 1 : 1e-10f;
            }

            for (int i = 0; i < l; i++)
            {
                int c1 = char.ToLowerInvariant(alphabet[i]), c1u = char.ToUpperInvariant(alphabet[i]);
                emitSingle[c1] = emitSingle[c1u] = single[i];
                for (int j = 0; j <= i; j++)
                {
                    int c2 = char.ToLowerInvariant(alphabet[j]), c2u = char.ToUpperInvariant(alphabet[j]);
                    float p = pairs[i * l + j];
                    emitPairs[c1, c2] = emitPairs[c1u, c2] = emitPairs[c1u, c2u] = emitPairs[c1, c2u] = p;
                    emitPairs[c2, c1] = emitPairs[c2u, c1] = emitPairs[c2u, c1u] = emitPairs[c2, c1u] = p;
                }
            }

            Setup(InsertStates, initDistrib, gapOpen, gapExtend, emitSingle, emitPairs);
        }

        /// <summary>
        /// The model of ProbCons 1.12 (<c>ProbabilisticModel (initDistrib, gapOpen, gapExtend, emitPairs, emitSingle)</c>, the constructor <see cref="TcProbaModel"/> shares with T-Coffee's
        /// <c>proba_pair</c>): <paramref name="insertStates"/> pairs of insertion states, emission tables indexed by the character code (0 ... 255).
        /// </summary>
        internal TcProbaModel(int insertStates, float[] initDistrib, float[] gapOpen, float[] gapExtend, float[] emitSingle, float[,] emitPairs)
        {
            InsertStates = insertStates;
            MatrixTypes = 2 * insertStates + 1;
            Setup(insertStates, initDistrib, gapOpen, gapExtend, emitSingle, emitPairs);
        }

        private void Setup(int insertStates, float[] initDistrib, float[] gapOpen, float[] gapExtend, float[] emitSingle, float[,] emitPairs)
        {
            int n = 2 * insertStates + 1;
            var transMat = new float[n, n];
            transMat[0, 0] = 1;
            for (int i = 0; i < InsertStates; i++)
            {
                transMat[0, 2 * i + 1] = gapOpen[2 * i];
                transMat[0, 2 * i + 2] = gapOpen[2 * i + 1];
                transMat[0, 0] -= gapOpen[2 * i] + gapOpen[2 * i + 1];
                transMat[2 * i + 1, 2 * i + 1] = gapExtend[2 * i];
                transMat[2 * i + 2, 2 * i + 2] = gapExtend[2 * i + 1];
                transMat[2 * i + 1, 2 * i + 2] = 0;
                transMat[2 * i + 2, 2 * i + 1] = 0;
                transMat[2 * i + 1, 0] = 1 - gapExtend[2 * i];
                transMat[2 * i + 2, 0] = 1 - gapExtend[2 * i + 1];
            }

            InitialDistribution = new float[2 * MatrixTypes + 1];
            TransProb = new float[n, n];
            for (int i = 0; i < MatrixTypes; i++)
            {
                InitialDistribution[i] = (float)Math.Log(initDistrib[i]);
                for (int j = 0; j < MatrixTypes; j++) TransProb[i, j] = (float)Math.Log(transMat[i, j]);
            }

            for (int i = 0; i < 256; i++)
            {
                for (int j = 0; j < MatrixTypes; j++) InsProb[i, j] = (float)Math.Log(emitSingle[i]);
                for (int j = 0; j < 256; j++) MatchProb[i, j] = (float)Math.Log(emitPairs[i, j]);
            }
        }
    }

    /// <summary><c>get_emitPairs(mat, alphabet, p, s)</c>: p[a][b] = s[a] s[b] exp(score(a, b) x ln 2 / 2) in single precision, divided by the sum of all entries.</summary>
    private static float[] TcEmitPairsFromMatrix(TcMatrix matrix, string alphabet, float[] single)
    {
        int l = alphabet.Length;
        float k = (float)(Math.Log(2) / 2);
        var p = new float[l * l];
        for (int a = 0; a < l; a++)
            for (int b = 0; b < l; b++)
            {
                float e = single[a] * single[b];
                int sc = matrix.Score(alphabet[a], alphabet[b]);
                p[a * l + b] = (float)(e * Math.Exp((double)sc * k));
            }

        float t = 0;
        for (int a = 0; a < l * l; a++) t += p[a];
        for (int a = 0; a < l * l; a++) p[a] /= t;
        return p;
    }

    /// <summary>
    /// Largest forward / backward matrix of the pair-HMM (<c>proba_pair</c>): <c>(len1 + 1)·(len2 + 1)·matrixTypes</c> floats (1 GiB each at the limit, forward and
    /// backward are both held: about 2.2 GB with the match table). The real <c>t_coffee</c> 13.41 has no limit of its own: <c>forward_proba_pair_wise</c>
    /// (<c>util_dp_suboptimal_nw.c</c>) computes <c>l = (len1 + 1)·(len2 + 1)·NumMatrixTypes</c> in 32-bit <c>int</c> (it wraps at 2^31 floats, about 27 000 x 27 000
    /// protein / 36 000 x 36 000 DNA residues) and <c>vcalloc</c>s it (it dies with an out-of-memory failure well before). Longer pairs get a clear
    /// <see cref="ArgumentException"/> here; use the <c>slow_pair</c> / <c>lalign_id_pair</c> primary libraries (linear-space Myers-Miller / Gotoh) for them.
    /// </summary>
    public const long TCoffeeProbaMaxMatrixFloats = 1L << 28;

    private static void RequireTcProbaFits(int length1, int length2, int matrixTypes, string program = "T-Coffee proba_pair")
    {
        long cells = ((long)length1 + 1) * ((long)length2 + 1);
        long floats = cells * matrixTypes;
        if (floats <= TCoffeeProbaMaxMatrixFloats)
            return;
        string programNote = program == "ProbCons"
            ? "ProbCons 1.12 itself computes (len1 + 1) x (len2 + 1) x numMatrixTypes in a 32-bit int (it wraps at 2^31 floats) and allocates it with new (it dies with an out-of-memory failure well before)"
            : "t_coffee itself overflows its 32-bit index at 2^31 floats and fails to allocate before";
        string alternative = program == "ProbCons"
            ? "use MultipleAlignKalign or MultipleAlignMafft for sequences of this length"
            : "use the slow_pair or lalign_id_pair primary library for sequences of this length";
        throw ResourceLimit($"The {program} pair-HMM", $"keeps a forward and a backward matrix of (len1 + 1) x (len2 + 1) x {matrixTypes} floats", $"{length1} x {length2} residues = {floats:N0} floats", TCoffeeProbaMaxMatrixFloats,
            $"SequenceAligner.{nameof(TCoffeeProbaMaxMatrixFloats)}, {TCoffeeProbaMaxMatrixFloats * 4L / (1L << 30)} GiB per matrix", programNote, alternative);
    }

    // dirichlet_code: a c d e f g h i k l m n p q r s t v w y -> 0..19; every other letter -> -1
    private const string TcDirichletLetters = "acdefghiklmnpqrstvwy";

    private static int TcDirichletCode(char c) => c < 128 ? TcDirichletLetters.IndexOf(char.ToLowerInvariant(c)) : -1;

    /// <summary>
    /// <c>proba_pair_wise</c> + <c>ProbaMatrix2CL</c> for one pair: the library pairs (1-based residues of the first and second sequence, weight) with
    /// posterior probability &gt; 0.01 of being aligned, at most 4 x min(length) of them, strongest first.
    /// </summary>
    internal static List<(int R1, int R2, int Weight)> TcProbaPairEntries(string seq1, string seq2, TcSeqKind kind, string? matrix = null)
    {
        var m = TcProbaModel.For(kind, matrix);
        int nmt = m.MatrixTypes, nis = m.InsertStates;
        int I = seq1.Length, J = seq2.Length;
        RequireTcProbaFits(I, J, nmt);
        var tIns = new float[2][,];
        tIns[0] = new float[nmt, I + 1];
        tIns[1] = new float[nmt, J + 1];
        var tMatch = new float[(I + 1) * (J + 1)];
        if (kind == TcSeqKind.Rna)
        {
            // get_tot_prob (RNA): the residue itself indexes the emission tables (case-insensitively, `c1 - 'a'`); any letter without an entry has the default 1e-5 / 1e-10
            for (int i = 1; i <= I; i++)
                for (int k = 0; k < nmt; k++)
                    tIns[0][k, i] = m.InsProb[seq1[i - 1], k];
            for (int j = 1; j <= J; j++)
                for (int k = 0; k < nmt; k++)
                    tIns[1][k, j] = m.InsProb[seq2[j - 1], k];
            for (int i = 1; i <= I; i++)
                for (int j = 1; j <= J; j++)
                    tMatch[i * (J + 1) + j] = m.MatchProb[char.ToLowerInvariant(seq1[i - 1]), char.ToLowerInvariant(seq2[j - 1])];
        }
        else
        {
            var d1 = seq1.Select(TcDirichletCode).ToArray();
            var d2 = seq2.Select(TcDirichletCode).ToArray();

            // get_tot_prob2: single sequences have one-hot profiles; a letter outside the 20 has an all-zero profile
            for (int i = 1; i <= I; i++)
                for (int k = 0; k < nmt; k++)
                    tIns[0][k, i] = d1[i - 1] >= 0 ? m.InsProb[TcDirichletLetters[d1[i - 1]], k] : 0f;
            for (int j = 1; j <= J; j++)
                for (int k = 0; k < nmt; k++)
                    tIns[1][k, j] = d2[j - 1] >= 0 ? m.InsProb[TcDirichletLetters[d2[j - 1]], k] : 0f;
            for (int i = 1; i <= I; i++)
                for (int j = 1; j <= J; j++)
                    tMatch[i * (J + 1) + j] = d1[i - 1] >= 0 && d2[j - 1] >= 0 ? m.MatchProb[TcDirichletLetters[d1[i - 1]], TcDirichletLetters[d2[j - 1]]] : 0f;
        }

        float[] fwd = TcProbaForward(I, J, nmt, nis, m, tMatch, tIns);
        float[] bwd = TcProbaBackward(I, J, nmt, nis, m, tMatch, tIns);

        float totalProb = TcProbaTotal(I, J, nmt, nis, fwd, bwd);
        const float thr = 0.01f;
        var list = new List<(int I, int J, int W)>();
        int ij = 0;
        for (int i = 0; i <= I; i++)
            for (int j = 0; j <= J; j++, ij += nmt)
            {
                if (j == 0) AlignmentCancellation.ThrowIfRequested();
                float arg = fwd[ij] + bwd[ij] - totalProb;
                double v = TcExp(arg < 0f ? arg : 0f);
                if (v > thr) list.Add((i, j, (int)((float)v * 1000f)));
            }

        // sort_int_inv: qsort ascending on the weight (stable), then the array is reversed
        var sorted = list.Select((e, idx) => (e, idx)).OrderBy(x => x.e.W).ThenBy(x => x.idx).Select(x => x.e).ToList();
        sorted.Reverse();
        int n = Math.Min(sorted.Count, 4 * Math.Min(I, J));
        return sorted.Take(n).Select(e => (e.I, e.J, e.W)).ToList();
    }

    internal static float TcProbaTotalProbability(int I, int J, int nmt, int nis, float[] forward, float[] backward) => TcProbaTotal(I, J, nmt, nis, forward, backward);

    private static float TcProbaTotal(int I, int J, int nmt, int nis, float[] forward, float[] backward)
    {
        float totalForward = TcLogZero;
        for (int k = 0; k < nmt; k++)
        {
            int idx = k + nmt * ((I + 1) * (J + 1) - 1);
            TcLogPlusEquals(ref totalForward, forward[idx] + backward[idx]);
        }

        int b0 = 0 + nmt * (1 * (J + 1) + 1);
        float totalBackward = forward[b0] + backward[b0];
        for (int k = 0; k < nis; k++)
        {
            int i1 = 2 * k + 1 + nmt * (1 * (J + 1) + 0);
            TcLogPlusEquals(ref totalBackward, forward[i1] + backward[i1]);
            int i2 = 2 * k + 2 + nmt * (0 * (J + 1) + 1);
            TcLogPlusEquals(ref totalBackward, forward[i2] + backward[i2]);
        }

        return (totalForward + totalBackward) / 2;
    }

    private static float[] TcProbaForward(int I, int J, int nmt, int nis, TcProbaModel m, float[] matchProb, float[][,] insProb)
    {
        int l = checked((int)(((long)I + 1) * (J + 1) * nmt));
        var forward = new float[l];
        Array.Fill(forward, TcLogZero);
        var transProb = m.TransProb;
        var init = m.InitialDistribution;

        forward[0 + nmt * (1 * (J + 1) + 1)] = init[0] + matchProb[J + 2];
        for (int k = 0; k < nis; k++)
        {
            forward[2 * k + 1 + nmt * (1 * (J + 1) + 0)] = init[2 * k + 1] + insProb[0][k, 1];
            forward[2 * k + 2 + nmt * (0 * (J + 1) + 1)] = init[2 * k + 2] + insProb[1][k, 1];
        }

        int ij = 0, i1j = -J - 1, ij1 = -1, i1j1 = -J - 2;
        ij *= nmt;
        i1j *= nmt;
        ij1 *= nmt;
        i1j1 *= nmt;
        int mIdx = 0;
        for (int i = 0; i <= I; i++)
            for (int j = 0; j <= J; j++, mIdx++)
            {
                if (j == 0) AlignmentCancellation.ThrowIfRequested();
                if (i > 1 || j > 1)
                {
                    if (i > 0 && j > 0)
                    {
                        forward[0 + ij] = forward[0 + i1j1] + transProb[0, 0];
                        for (int k = 1; k < nmt; k++)
                            TcLogPlusEquals(ref forward[0 + ij], forward[k + i1j1] + transProb[k, 0]);
                        forward[0 + ij] += matchProb[mIdx];
                    }

                    if (i > 0)
                        for (int k = 0; k < nis; k++)
                            forward[2 * k + 1 + ij] = insProb[0][k, i] + TcLogAdd(forward[0 + i1j] + transProb[0, 2 * k + 1], forward[2 * k + 1 + i1j] + transProb[2 * k + 1, 2 * k + 1]);

                    if (j > 0)
                        for (int k = 0; k < nis; k++)
                            forward[2 * k + 2 + ij] = insProb[1][k, j] + TcLogAdd(forward[0 + ij1] + transProb[0, 2 * k + 2], forward[2 * k + 2 + ij1] + transProb[2 * k + 2, 2 * k + 2]);
                }

                ij += nmt;
                i1j += nmt;
                ij1 += nmt;
                i1j1 += nmt;
            }

        return forward;
    }

    private static float[] TcProbaBackward(int I, int J, int nmt, int nis, TcProbaModel m, float[] matchProb, float[][,] insProb, bool endMatch = true)
    {
        int l = checked((int)(((long)I + 1) * (J + 1) * nmt));
        var backward = new float[l];
        Array.Fill(backward, TcLogZero);
        var transProb = m.TransProb;
        var init = m.InitialDistribution;

        for (int k = 0; k < nmt; k++)
            backward[nmt * ((I + 1) * (J + 1) - 1) + k] = init[k];
        // T-Coffee's copy of the pair-HMM adds the emission of the last cell to the match state's start value; ProbCons 1.12 does not (endMatch = false)
        if (endMatch) backward[nmt * ((I + 1) * (J + 1) - 1) + 0] += matchProb[(I + 1) * (J + 1) - 1];

        int ij = (I + 1) * (J + 1) - 1;
        int i1j = ij + J + 1;
        int ij1 = ij + 1;
        int i1j1 = ij + J + 2;
        ij *= nmt;
        i1j *= nmt;
        ij1 *= nmt;
        i1j1 *= nmt;

        for (int i = I; i >= 0; i--)
            for (int j = J; j >= 0; j--)
            {
                if (j == J) AlignmentCancellation.ThrowIfRequested();
                if (i < I && j < J)
                {
                    int mm = ((i + 1) * (J + 1)) + j + 1;
                    float probXY = backward[0 + i1j1] + matchProb[mm];
                    for (int k = 0; k < nmt; k++)
                        TcLogPlusEquals(ref backward[k + ij], probXY + transProb[k, 0]);
                }

                if (i < I)
                    for (int k = 0; k < nis; k++)
                    {
                        TcLogPlusEquals(ref backward[0 + ij], backward[2 * k + 1 + i1j] + insProb[0][k, i + 1] + transProb[0, 2 * k + 1]);
                        TcLogPlusEquals(ref backward[2 * k + 1 + ij], backward[2 * k + 1 + i1j] + insProb[0][k, i + 1] + transProb[2 * k + 1, 2 * k + 1]);
                    }

                if (j < J)
                    for (int k = 0; k < nis; k++)
                    {
                        TcLogPlusEquals(ref backward[0 + ij], backward[2 * k + 2 + ij1] + insProb[1][k, j + 1] + transProb[0, 2 * k + 2]);
                        TcLogPlusEquals(ref backward[2 * k + 2 + ij], backward[2 * k + 2 + ij1] + insProb[1][k, j + 1] + transProb[2 * k + 2, 2 * k + 2]);
                    }

                ij -= nmt;
                i1j -= nmt;
                ij1 -= nmt;
                i1j1 -= nmt;
            }

        return backward;
    }
}
