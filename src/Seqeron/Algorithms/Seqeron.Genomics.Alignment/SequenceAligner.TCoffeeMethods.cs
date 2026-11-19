namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// Further pair-wise methods of the T-Coffee 13.41 program (audit A2-3 / A2-4, B09): the DP routines that `-method NAME` runs on every pair of sequences.
// Each routine is a line-by-line port of the source of the Debian package t-coffee 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1 (t_coffee_source/dp_lib/*.c); the rows are
// compared with the installed `t_coffee` (scripts/alignment/generate_tcoffee_options_oracle.py, features biphasic ...).
// =====================================================================================================================

public static partial class SequenceAligner
{
    /// <summary>
    /// <c>biphasic_pair</c>: the global alignment of <c>biphasic_pair_wise</c> (util_dp_suboptimal_nw.c): five states (match; two affine gap models of each kind: opening x 2 and extension / 2, opening / 2 and
    /// extension x 2, integer divisions), ties broken by the order of the program's <c>max_int</c> calls (first listed wins), with the program's traceback that starts from the PREDECESSOR state of the best final state.
    /// Parameters as <c>slow_pair</c> (DNA: dna_idmat, -10, -1; protein: BLOSUM62, gap open (int)(10 x mean mismatch), extend -1).
    /// </summary>
    /// <param name="sequenceA">First sequence.</param>
    /// <param name="sequenceB">Second sequence.</param>
    /// <param name="protein"><c>true</c> = the protein parameters.</param>
    /// <param name="matrix"><c>-matrix</c> (protein only).</param>
    /// <param name="parameters">The <c>@EP@</c> overrides of the method (GOP, GEP, MATRIX).</param>
    public static (string RowA, string RowB) TCoffeeBiphasicPairAlign(string sequenceA, string sequenceB, bool protein = false, string? matrix = null, TCoffeeMethodParameters? parameters = null)
    {
        ArgumentNullException.ThrowIfNull(sequenceA);
        ArgumentNullException.ThrowIfNull(sequenceB);
        if (sequenceA.Length == 0 || sequenceB.Length == 0) throw new ArgumentException("Empty sequences are not supported.");
        var p = TcSlowPairParams(protein, matrix, parameters);
        return TcBiphasic(sequenceA, sequenceB, p.Matrix, p.Gop, p.Gep);
    }

    /// <summary>The move of the program's traceback: 0 = both residues, 1 = a residue of the second sequence against a gap, 2 = a residue of the first sequence against a gap (the borders first, then the state).</summary>
    private static int BiphasicStep(int x, int y, int trace)
    {
        const int M1 = 1, D1 = 2, D2 = 3, I1 = 4, I2 = 5;
        if (x == 0) return 1;
        if (y == 0) return 2;
        if (trace == M1) return 0;
        if (trace == D1 || trace == D2) return 1;
        if (trace == I1 || trace == I2) return 2;
        throw new InvalidOperationException("biphasic_pair_wise: the traceback reached an undefined state.");
    }

    private static (string, string) TcBiphasic(string a, string b, TcMatrix matrix, int gopIn, int gepIn)
    {
        var (ia, ib, _) = TcBiphasicMasks(a.Length, b.Length, (i, j) => matrix.Score(a[i], b[j]) * TcScoreK, gopIn, gepIn);
        return TcMasksToRows(a, b, ia, ib);
    }

    /// <summary><c>biphasic_pair_wise</c> on a score function of the two 0-based positions (<c>get_dp_cost</c>): the residue masks of the two rows and the score of the best final state.</summary>
    private static (List<bool> Present1, List<bool> Present2, int Score) TcBiphasicMasks(int I, int J, Func<int, int, int> cost, int gopIn, int gepIn)
    {
        const int M1 = 1, D1 = 2, D2 = 3, I1 = 4, I2 = 5, N = 6;
        var m = new int[N][,];
        var t = new int[N][,];
        for (int c = 0; c < N; c++)
        {
            m[c] = new int[I + 1, J + 1];
            t[c] = new int[I + 1, J + 1];
            for (int i = 0; i <= I; i++)
                for (int j = 0; j <= J; j++) m[c][i, j] = -999999;
        }

        int gop1 = gopIn * TcScoreK * 2, gep1 = gepIn * TcScoreK / 2;
        int gop2 = gopIn * TcScoreK / 2, gep2 = gepIn * TcScoreK * 2;

        m[M1][0, 0] = 0;
        for (int i = 1; i <= I; i++) { m[I1][i, 0] = gep1 * i; m[I2][i, 0] = gep2 * i; }
        for (int j = 1; j <= J; j++) { m[D1][0, j] = gep1 * j; m[D2][0, j] = gep2 * j; }

        for (int i = 1; i <= I; i++)
            for (int j = 1; j <= J; j++)
            {
                int sub = cost(i - 1, j - 1);
                // max_int: the first listed key wins ties (a later value must be strictly greater)
                int best = m[D1][i - 1, j - 1], key = D1;
                if (m[I1][i - 1, j - 1] > best) { best = m[I1][i - 1, j - 1]; key = I1; }
                if (m[M1][i - 1, j - 1] > best) { best = m[M1][i - 1, j - 1]; key = M1; }
                if (m[D2][i - 1, j - 1] > best) { best = m[D2][i - 1, j - 1]; key = D2; }
                if (m[I2][i - 1, j - 1] > best) { best = m[I2][i - 1, j - 1]; key = I2; }
                t[M1][i, j] = key;
                m[M1][i, j] = best + sub;

                Affine(D1, M1, m[M1][i, j - 1] + gop1, m[D1][i, j - 1], gep1);
                Affine(I1, M1, m[M1][i - 1, j] + gop1, m[I1][i - 1, j], gep1);
                Affine(D2, M1, m[M1][i, j - 1] + gop2, m[D2][i, j - 1], gep2);
                Affine(I2, M1, m[M1][i - 1, j] + gop2, m[I2][i - 1, j], gep2);

                void Affine(int state, int openKey, int open, int extend, int gep)
                {
                    bool takeExtend = extend > open;
                    t[state][i, j] = takeExtend ? state : openKey;
                    m[state][i, j] = (takeExtend ? extend : open) + gep;
                }
            }

        int bestFinal = m[M1][I, J], trace = M1;
        foreach (int st in new[] { D1, I1, D2, I2 })
            if (m[st][I, J] > bestFinal) { bestFinal = m[st][I, J]; trace = st; }

        var ia = new List<bool>(); var ib = new List<bool>();
        int x = I, y = J;
        trace = t[trace][x, y];
        while (!(x == 0 && y == 0))
        {
            int next = t[trace][x, y];
            int step = BiphasicStep(x, y, trace);
            if (step != 1) { ia.Add(true); x--; } else ia.Add(false);
            if (step != 2) { ib.Add(true); y--; } else ib.Add(false);
            trace = next;
        }

        ia.Reverse(); ib.Reverse();
        return (ia, ib, bestFinal);
    }

    /// <summary>
    /// <c>glocal_pair_wise</c> (util_dp_gotoh_nw.c L370): four states (M1, D1, I1, M2 = an unmatched pair costing <c>gep</c>, entered with <c>gop</c>) with zero-initialised borders except the listed ones, ties to the first listed key, the traceback
    /// that starts from the predecessor of the best final state, and a pair in state M2 written as two columns (the first row's residue against a gap, then the gap against the second row's residue).
    /// </summary>
    private static (List<bool> Present1, List<bool> Present2, int Score) TcGlocalMasks(int l1, int l2, Func<int, int, int> cost, int gopIn, int gepIn)
    {
        const int M1 = 1, D1 = 2, I1 = 3, M2 = 4, N = 5;
        int gop = gopIn * TcScoreK, gep = gepIn * TcScoreK, tgop = gop, unmatch = gep;
        var m = new int[N][,]; var t = new int[N][,];
        for (int c = 0; c < N; c++) { m[c] = new int[l1 + 1, l2 + 1]; t[c] = new int[l1 + 1, l2 + 1]; }
        for (int j = 1; j <= l2; j++)
        {
            m[D1][0, j] = gep * j;
            m[M1][0, j] = 2 * gep * j;
            m[M2][0, j] = 4 * gep * j;
        }

        for (int i = 1; i <= l1; i++)
        {
            m[I1][i, 0] = i * gep;
            m[M2][i, 0] = 4 * i * gep;
            m[M1][i, 0] = 2 * i * gep;
            for (int j = 1; j <= l2; j++)
            {
                int sub = cost(i - 1, j - 1);
                m[M1][i, j] = DpMax(out t[M1][i, j], (M1, m[M1][i - 1, j - 1]), (I1, m[I1][i - 1, j - 1]), (D1, m[D1][i - 1, j - 1]), (M2, m[M2][i - 1, j - 1])) + sub;
                m[D1][i, j] = DpMax(out t[D1][i, j], (M1, m[M1][i, j - 1] + gop), (D1, m[D1][i, j - 1] + gep), (M2, m[M2][i, j - 1]));
                m[I1][i, j] = DpMax(out t[I1][i, j], (M1, m[M1][i - 1, j] + gop), (I1, m[I1][i - 1, j] + gep), (M2, m[M2][i - 1, j]));
                m[M2][i, j] = DpMax(out t[M2][i, j], (M1, m[M1][i - 1, j - 1] + tgop), (I1, m[I1][i - 1, j - 1] + tgop), (D1, m[D1][i - 1, j - 1] + tgop), (M2, m[M2][i - 1, j - 1])) + unmatch;
            }
        }

        int score = DpMax(out int trace, (M1, m[M1][l1, l2]), (D1, m[D1][l1, l2]), (I1, m[I1][l1, l2]), (M2, m[M2][l1, l2]));
        int ii = l1, jj = l2;
        var al0 = new List<bool>(); var al1 = new List<bool>();
        trace = t[trace][ii, jj];
        while (!(ii == 0 && jj == 0))
        {
            int next = t[trace][ii, jj];
            // the borders first (i == 0: deletion, j == 0: insertion), then the state of the trace
            int action;
            if (ii == 0) action = D1;
            else if (jj == 0) action = I1;
            else action = trace;
            switch (action)
            {
                case M1:
                    al0.Add(true); al1.Add(true);
                    ii--; jj--;
                    break;
                case M2:
                    al0.Add(true); al1.Add(false);
                    al0.Add(false); al1.Add(true);
                    ii--; jj--;
                    break;
                case D1:
                    al0.Add(false); al1.Add(true);
                    jj--;
                    break;
                case I1:
                    al0.Add(true); al1.Add(false);
                    ii--;
                    break;
            }

            trace = next;
        }

        al0.Reverse(); al1.Reverse();
        return (al0, al1, score);
    }

    /// <summary><c>dp_max (&amp;trace, n, key, value, ...)</c>: the value of the first key, a later value replaces it only when strictly greater.</summary>
    private static int DpMax(out int trace, params (int Key, int Value)[] candidates)
    {
        int best = candidates[0].Value;
        trace = candidates[0].Key;
        for (int a = 1; a < candidates.Length; a++)
            if (candidates[a].Value > best)
            {
                best = candidates[a].Value;
                trace = candidates[a].Key;
            }

        return best;
    }

    /// <summary>The methods whose use the program refuses (see <see cref="TCoffeePairMethod.CleanSlowPair"/>, <see cref="TCoffeePairMethod.HashPair"/>, <see cref="TCoffeePairMethod.Subop1Pair"/>, <see cref="TCoffeePairMethod.Subop2Pair"/>, <see cref="TCoffeePairMethod.CfastPair"/>).</summary>
    private static void TcRefuseMethod(TCoffeePairMethod method)
    {
        switch (method)
        {
            case TCoffeePairMethod.CleanSlowPair:
                throw new TCoffeeAbortedException("ERROR: Weight Mode clean is unknown [FATAL:T-COFFEE] (seqpair2weight, util_constraints_list.c; COREDUMP in the pair job)");
            case TCoffeePairMethod.HashPair:
                throw new TCoffeeAbortedException("COREDUMP: THE FUNCTION hash_pair DOES NOT EXIST (the pair-job dispatcher of util_dp_drivers.c has no branch for 'EXECUTABLE hash_pair')");
            case TCoffeePairMethod.Subop1Pair:
                throw new TCoffeeAbortedException("ERROR: Msubop1_pair not recognised [FATAL:T-COFFEE] -- subop1_pair is not a valid method");
            case TCoffeePairMethod.Subop2Pair:
                throw new TCoffeeAbortedException("ERROR: Msubop2_pair not recognised [FATAL:T-COFFEE] -- subop2_pair is not a valid method");
            case TCoffeePairMethod.CfastPair:
                throw new TCoffeeAbortedException("ERROR: Mcfast_pair not recognised [FATAL:T-COFFEE] -- cfast_pair is not a valid method");
        }
    }

    // ----- promo_pair / blastr_pair: slow_pair on the dinucleotide-extended alphabet (EXTEND_SEQ 1, WEIGHT extend_sim) ---------------------------------------------------------------

    private static bool TcIsGapChar(char c) => "-.#*~".Contains(c);

    /// <summary>
    /// <c>extend_seq</c> (reformat.c): every residue is replaced by the code of the dinucleotide it forms with the next residue of the UNGAPPED sequence (aa -> d, ac -> e, ... tt -> y; a neighbour or itself n / x
    /// keeps the first letter / gives x; letters outside a c g t u n x give j); the last residue keeps its letter, gaps stay, the case of the first residue of the pair decides the case of the code.
    /// </summary>
    internal static string TcExtendSeq(string seq)
    {
        var buf = seq.Where(c => !TcIsGapChar(c)).ToArray();
        var ebuf = (char[])buf.Clone();
        for (int a = 0; a < buf.Length - 1; a++)
        {
            char r1 = buf[a], r2 = buf[a + 1];
            bool upper = char.IsAsciiLetterUpper(r1);
            r1 = char.ToLowerInvariant(r1);
            r2 = char.ToLowerInvariant(r2);
            if (r1 == 'u') r1 = 't';
            if (r2 == 'u') r2 = 't';
            char v;
            if (r1 == 'x' || r1 == 'n') v = 'x';
            else if (r2 == 'n' || r2 == 'x') v = r1;
            else if (r1 == 'a' && r2 == 'a') v = 'd';
            else if (r1 == 'a' && r2 == 'c') v = 'e';
            else if (r1 == 'a' && r2 == 'g') v = 'f';
            else if (r1 == 'a' && r2 == 't') v = 'h';
            else if (r1 == 'c' && r2 == 'a') v = 'i';
            else if (r1 == 'c' && r2 == 'c') v = 'k';
            else if (r1 == 'c' && r2 == 'g') v = 'l';
            else if (r1 == 'c' && r2 == 't') v = 'm';
            else if (r1 == 'g' && r2 == 'a') v = 'n';
            else if (r1 == 'g' && r2 == 'c') v = 'p';
            else if (r1 == 'g' && r2 == 'g') v = 'q';
            else if (r1 == 'g' && r2 == 't') v = 'r';
            else if (r1 == 't' && r2 == 'a') v = 's';
            else if (r1 == 't' && r2 == 'c') v = 'v';
            else if (r1 == 't' && r2 == 'g') v = 'w';
            else if (r1 == 't' && r2 == 't') v = 'y';
            else v = 'j';
            ebuf[a] = upper ? char.ToUpperInvariant(v) : v;
        }

        return TcPlaceResidues(seq, ebuf);
    }

    /// <summary><c>unextend_seq</c>: the first nucleotide of the dinucleotide code (a d e f h -> a, c i k l m -> c, g n p q r -> g, t s v w y -> t; x -> n; everything else, j included, -> j); lossy for other letters.</summary>
    internal static string TcUnextendSeq(string seq)
    {
        var buf = seq.Where(c => !TcIsGapChar(c)).ToArray();
        var ebuf = (char[])buf.Clone();
        for (int a = 0; a < buf.Length - 1; a++)
        {
            char r1 = buf[a];
            bool upper = char.IsAsciiLetterUpper(r1);
            r1 = char.ToLowerInvariant(r1);
            if (r1 == 'u') r1 = 't';
            char v;
            if (r1 == 'x') v = 'n';
            else if (r1 is 'a' or 'd' or 'e' or 'f' or 'h') v = 'a';
            else if (r1 is 'c' or 'i' or 'k' or 'l' or 'm') v = 'c';
            else if (r1 is 'g' or 'n' or 'p' or 'q' or 'r') v = 'g';
            else if (r1 is 't' or 's' or 'v' or 'w' or 'y') v = 't';
            else v = 'j';
            ebuf[a] = upper ? char.ToUpperInvariant(v) : v;
        }

        return TcPlaceResidues(seq, ebuf);
    }

    private static string TcPlaceResidues(string seq, char[] residues)
    {
        var chars = seq.ToCharArray();
        for (int a = 0, b = 0; a < chars.Length; a++)
            if (!TcIsGapChar(chars[a])) chars[a] = residues[b++];
        return new string(chars);
    }

    /// <summary>The parameters of <c>promo_pair</c> (matrix promoter_tf1, GOP -60, GEP -1) and <c>blastr_pair</c> (matrix blosumR, GOP -20, GEP 0) with the <c>@EP@</c> overrides; every sequence type.</summary>
    private static TcPairParams TcExtendedPairParams(TCoffeePairMethod method, TCoffeeMethodParameters? mp)
    {
        var d = method == TCoffeePairMethod.PromoPair
            ? new TcPairParams(TcMatrix.ByName("promoter_tf1")!, -60, -1)
            : new TcPairParams(TcMatrix.ByName("blosumR")!, -20, 0);
        return new TcPairParams(TcResolveMatrix(mp?.Matrix) ?? d.Matrix, mp?.GapOpen ?? d.Gop, mp?.GapExtend ?? d.Gep);
    }

    /// <summary>
    /// One pair job of <c>promo_pair</c> / <c>blastr_pair</c> as <c>fast_pair</c> and <c>aln2constraint_list(.., "extend_sim")</c> run it on the program's SHARED sequence array <paramref name="work"/> (the jobs of one method run
    /// one after the other in one process, <c>-n_core 1</c>): the rows are copies of the current sequences, all sequences and the rows are extended in place, the pair is aligned (Myers-Miller) on the extended letters,
    /// everything is unextended (lossy: protein letters become j for good, u becomes t), then extended once more for the weight (identity of the extended rows x 10) and the residue positions, and unextended again.
    /// </summary>
    private static List<(int S1, int R1, int S2, int R2, int W)> TcExtendedPairEntries(string[] work, int x, int y, TCoffeePairMethod method, TCoffeeMethodParameters? mp)
    {
        var p = TcExtendedPairParams(method, mp);
        string rowA = work[x], rowB = work[y];
        for (int i = 0; i < work.Length; i++) work[i] = TcExtendSeq(work[i]);
        rowA = TcExtendSeq(rowA);
        rowB = TcExtendSeq(rowB);
        var (ra, rb) = TcMyersMiller(rowA, rowB, p.Matrix, p.Gop, p.Gep);
        for (int i = 0; i < work.Length; i++) work[i] = TcUnextendSeq(work[i]);
        ra = TcUnextendSeq(ra);
        rb = TcUnextendSeq(rb);

        // aln2constraint_list: weight_mode "extend_sim" extends the alignment and the sequence set again
        ra = TcExtendSeq(ra);
        rb = TcExtendSeq(rb);
        for (int i = 0; i < work.Length; i++) work[i] = TcExtendSeq(work[i]);
        int weight = TcSeqSim(ra, rb);
        if (weight == 0) weight = 1;
        weight *= 1000 / TcMaxId;
        int[] fa = TcFixIndex(ra.Replace("-", ""), work[x]);
        int[] fb = TcFixIndex(rb.Replace("-", ""), work[y]);
        var entries = new List<(int, int, int, int, int)>();
        int na = 0, nb = 0;
        for (int c = 0; c < ra.Length; c++)
        {
            bool ga = ra[c] == '-', gb = rb[c] == '-';
            if (!ga) na++;
            if (!gb) nb++;
            if (!ga && !gb) entries.Add((x, fa[na], y, fb[nb], weight));
        }

        for (int i = 0; i < work.Length; i++) work[i] = TcUnextendSeq(work[i]);
        return entries;
    }
}
