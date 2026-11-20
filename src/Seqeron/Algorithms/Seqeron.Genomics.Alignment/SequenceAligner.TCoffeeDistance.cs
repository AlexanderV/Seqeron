namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// The distance matrix, guide tree modes and sequence weights of the T-Coffee 13.41 program (audit A1-13d, B09).
//
//   seq2distance_matrix (util_constraints_list.c) fills for every pair (a < b) an identity `id` (similarity_matrix, 0..100) and a `score` (score_similarity_matrix):
//     ktup    id = ktup_dist_mat similarity (default)                                   score = id x SCORE_K (10)
//     fast    id = get_seq_sim of the Myers-Miller alignment (idmat / blosum62mt)      score = id x 10
//     idscore id = idscore_pairseq (affine DP, own gap rule)                           score = id            (NOT scaled: a quirk, all distances are then near 100)
//     no      id = 100                                                                  score = 1000
//   make_tree: distances = normalize_array(sim_array2dist_array(score, 1000), 1000, 100) (integer division) and then
//     nj    : slow_nj_tree / fast_nj_tree on distance / 100, reroot of the trichotomy (as the default tree)
//     upgma : int_dist2upgma_tree: merge of the lowest integer distance (first on ties, rows a < b), the new row is the integer mean (x + y) / 2.
//   weight_constraint_list(seq_weight = t_coffee): compute_t_coffee_weight on `id` (float arithmetic), re_weight_constraint_list multiplies every directed entry by
//   MIN(w[s1], w[s2]) in float and truncates.
// =====================================================================================================================

public static partial class SequenceAligner
{
    private const int TcMaxIdentity = 100;

    /// <summary>Number of leading pair alignments of the distance matrix for which a gap in the first row of <c>idscore_pairseq</c> kills the program every time (measured, see <see cref="TCoffeeAbortedException.HeapDependent"/>).</summary>
    private const int TcIdscoreWildPointerCalls = 12;

    /// <summary>
    /// The state of the program's constraint list when the distance matrix is computed (t_coffee.c L5156-5172): after <c>-filter_lib</c>, the relaxation and the first
    /// <c>evaluate_constraint_list_reference</c> (so the <c>nomatch</c> penalty is already rescaled), before the sequence weights re-weight the library.
    /// </summary>
    private sealed record TcRelaxedState(TcIndex Index, TcPenalties Penalties, int Normalise, int TgMode, TcDpSettings? Dp = null);

    private sealed class TcDistances
    {
        public int[,] Sim = new int[0, 0];
        public int[,] Score = new int[0, 0];
    }

    /// <summary>
    /// <c>make_group_aa(&amp;n, mode)</c> for the first token of <c>-distance_matrix_sim_mode</c> (the text before its first '_'): the residue groups <c>is_in_same_group_aa</c> tests
    /// (two residues are similar when they share a group; case-insensitive, groups may overlap). The groups are built ONCE per process, by the first call, which is why the caller passes the program's
    /// first mode. <c>null</c> = no group at all (a mode the function does not know: nothing is similar, not even identical letters).
    /// </summary>
    internal static IReadOnlyList<string>? TcGroupsOfMode(string mode)
    {
        if (mode.Length == 0) return null; // handled by the caller: the groups of the matrix idmat
        if (mode.Contains("sim", StringComparison.Ordinal) || mode == "idmat")
            return Enumerable.Range('a', 26).Select(c => $"{(char)c}{char.ToUpperInvariant((char)c)}").ToArray();
        return mode switch
        {
            "simple" => new[] { "avilmAVILM", "dekrDEKR", "stcnqhSTCNQH", "wfyWFY", "gG", "pP" },
            "mafft" => new[] { "agjopstAGJOPST", "ilmvILMV", "bdenqzBDENQZ", "hkrHKR", "fwyFWY", "cC" },
            "clustalw" => new[] { "astaASTA", "bneqkBNEQK", "cnhqkCNHQK", "dndeqDNDEQ", "eqhrkEQHRK", "fmilvFMILV", "gmilfGMILF", "hhyHHY", "ifywIFYW", "jcJC", "kpKP" },
            "polarity" => new[] { "eqrsdnkhtEQRSDNKHT", "pP", "gG", "cC", "fywFYW", "iavlmIAVLM" },
            "vasiliky" => new[] { "rkRK", "deDE", "qhQH", "vilmVILM", "fyFY", "sS", "wW", "aA", "cC", "gG", "nN", "pP", "tT" },
            "make_all" => new[] { "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz" },
            _ => Array.Empty<string>(),
        };
    }

    /// <summary>
    /// <c>is_in_same_group_aa(r1, r2, 0, NULL, mode)</c> of the distance measure. The pair jobs run in forked children, so the first call in the main process (which builds the process-static groups) is the distance matrix
    /// itself: the groups are those of the sim mode (an empty mode = the groups of the matrix <c>idmat</c>, T and U together, J / O never equal). A1-13d had assumed that a library method makes the first call and fixes the
    /// groups of <c>idmat</c>; sets built from T and U show it does not (the sim mode <c>idmat_sim1</c> = one group per letter).
    /// </summary>
    private static bool TcSameGroupMode(char a, char b, IReadOnlyList<string>? modeGroups)
    {
        if (modeGroups == null) return TcSameGroup(a, b);
        char ua = char.ToUpperInvariant(a), ub = char.ToUpperInvariant(b);
        foreach (string g in modeGroups)
            if (g.Contains(ua) && g.Contains(ub)) return true; // is_in_set: both letter cases are in every group string
        return false;
    }

    /// <summary><c>get_seq_sim(row1, row2, "-", "idmat_sim1" | ...)</c> for the measures of <c>-distance_matrix_sim_mode</c>.</summary>
    private static int TcSeqSimMeasure(string r1, string r2, string simMode, int callIndex = 0)
    {
        string mode = simMode.TrimStart('_');
        string measure = "";
        bool hasMeasure = false;
        int us = mode.IndexOf('_');
        if (us >= 0)
        {
            measure = mode[(us + 1)..];
            mode = mode[..us];
            hasMeasure = true;
        }

        if (mode.Contains("idscore", StringComparison.Ordinal))
        {
            // get_seq_sim hands the ALIGNED rows to idscore_pairseq(s1, s2, -12, -1, blosum62mt, mode), which reads m[s1[i] - 'a'][s2[j] - 'a']: a gap in the first row indexes 52 rows before the
            // matrix (a wild pointer: the program dies with a segmentation fault, exit status 11), a gap in the second row reads a zero. The measure after '_' is not passed on (the comp mode is the part before it).
            // measured on 220 committed sets with such a gap: the program died in all 208 whose first gap lies in one of the first 12 pair alignments (a < b order); later (calls 12-22) it died in 7 of 12 sets: the memory
            // before the matrix then holds heap data that is a valid pointer in some runs
            if (r1.Contains('-'))
                throw new TCoffeeAbortedException("segmentation fault in idscore_pairseq: -distance_matrix_sim_mode '" + simMode + "' scores the aligned rows of -distance_matrix_mode fast, and a gap '-' in the first row indexes the matrix at 'a' - 45 (get_seq_sim, aln_convertion_util.c)", 11, callIndex >= TcIdscoreWildPointerCalls);
            return TcIdscorePair(r1, r2, -12, -1, TcMatrix.Blosum62, mode);
        }

        var modeGroups = TcGroupsOfMode(mode);
        int pos0 = 0, pos1 = 0, pos2 = 0, sim = 0, gap = 0;
        for (int a = 0; a < r1.Length; a++)
        {
            bool p1 = r1[a] != '-', p2 = r2[a] != '-';
            if (p1) pos1++;
            if (p2) pos2++;
            if (p1 && p2)
            {
                pos0++;
                if (TcSameGroupMode(r1[a], r2[a], modeGroups)) sim++;
            }
            else if (p1 != p2) gap++;
        }

        if (mode.Contains("cov", StringComparison.Ordinal)) return pos0 + gap == 0 ? 0 : pos0 * TcMaxIdentity / (pos0 + gap);
        if (!hasMeasure || measure == "sim1" || measure == "sim") return pos0 == 0 ? 0 : sim * TcMaxIdentity / pos0;
        if (measure == "sim2") return pos1 == 0 || pos2 == 0 ? 0 : sim * TcMaxIdentity / Math.Min(pos1, pos2);
        if (measure == "sim3") return pos1 == 0 || pos2 == 0 ? 0 : sim * TcMaxIdentity / Math.Max(pos1, pos2);
        if (measure == "gap1") return TcMaxIdentity - (r1.Length == 0 ? TcMaxIdentity : gap * TcMaxIdentity / r1.Length);
        if (measure == "logid") return TcLogId(pos0, sim);
        if (!mode.Contains("sim", StringComparison.Ordinal) || pos0 == 0) return 0;
        return sim * TcMaxIdentity / pos0;
    }

    /// <summary><c>logid_score(len = pos0, sim)</c> (note the program's argument order: its parameters are named <c>(sim, len)</c> and called with (pos0, sim)).</summary>
    private static int TcLogId(int sim, int len)
    {
        if (len == 0) return (int)(0.33 * TcMaxIdentity);
        float score = (float)sim / len;
        score = score > 0.9f ? 1.0f : (float)-Math.Log10(1.0 - score);
        score *= TcMaxIdentity;
        return (int)score;
    }

    /// <summary><c>idscore_pairseq</c>: the affine DP of util_dp_gotoh_nw.c that returns the identity directly (letters lower-cased, terminal gaps free of the opening penalty).</summary>
    private static int TcIdscorePair(string seq1, string seq2, int gop, int gep, TcMatrix matrix, string compMode)
    {
        string s1 = seq1.ToLowerInvariant(), s2 = seq2.ToLowerInvariant();
        int l1 = s1.Length, l2 = s2.Length;
        // rows: [0]/[1] score previous / current, [2]/[3] identities, [4]/[5] aligned columns
        var I = new int[6][]; var D = new int[6][]; var M = new int[6][];
        for (int r = 0; r < 6; r++) { I[r] = new int[l2 + 1]; D[r] = new int[l2 + 1]; M[r] = new int[l2 + 1]; }
        for (int j = 0; j <= l2; j++) { D[0][j] = gep * j; M[0][j] = 2 * gep * j; D[4][j] = 0; }
        for (int i = 1; i <= l1; i++)
        {
            AlignmentCancellation.ThrowIfRequested();       // A5-2: one poll per row
            I[1][0] = i * gep;
            M[1][0] = 2 * i * gep;
            for (int j = 1; j <= l2; j++)
            {
                // m[s1[i] - 'a'][s2[j] - 'a']: a gap '-' in the SECOND string indexes 52 entries before the row, which holds 0 (a gap in the first string is a wild pointer: see TcSeqSimMeasure)
                int score = char.IsAsciiLetter(s1[i - 1]) && char.IsAsciiLetter(s2[j - 1]) ? matrix.Score(char.ToUpperInvariant(s1[i - 1]), char.ToUpperInvariant(s2[j - 1])) : 0;
                int id = s1[i - 1] == s2[j - 1] ? 1 : 0;
                int igop = i == l1 || j == l2 ? 0 : gop;

                if (D[0][j] + gep > M[0][j] + igop + gep) { D[1][j] = D[0][j] + gep; D[3][j] = D[2][j]; D[5][j] = D[4][j]; }
                else { D[1][j] = M[0][j] + igop + gep; D[3][j] = M[2][j]; D[5][j] = M[4][j]; }

                if (I[1][j - 1] + gep > M[1][j - 1] + igop + gep) { I[1][j] = I[1][j - 1] + gep; I[3][j] = I[3][j - 1]; I[5][j] = I[5][j - 1]; }
                else { I[1][j] = M[1][j - 1] + igop + gep; I[3][j] = M[3][j - 1]; I[5][j] = M[5][j - 1]; }

                int match = M[0][j - 1] + score;
                if (I[1][j] > match && I[1][j] > D[1][j]) { M[1][j] = I[1][j]; M[3][j] = I[3][j]; M[5][j] = I[5][j]; }
                else if (D[1][j] > match) { M[1][j] = D[1][j]; M[3][j] = D[3][j]; M[5][j] = D[5][j]; }
                else { M[1][j] = match; M[3][j] = M[2][j - 1] + id; M[5][j] = M[4][j - 1] + 1; }
            }

            for (int r = 0; r < 6; r += 2)
            {
                (I[r], I[r + 1]) = (I[r + 1], I[r]);
                (D[r], D[r + 1]) = (D[r + 1], D[r]);
                (M[r], M[r + 1]) = (M[r + 1], M[r]);
            }
        }

        int l;
        if (compMode.Contains("sim2", StringComparison.Ordinal)) { l = Math.Min(l1, l2); return l == 0 ? 0 : M[2][l2] * 100 / l; }
        if (compMode.Contains("sim3", StringComparison.Ordinal)) { l = Math.Max(l1, l2); return l == 0 ? 0 : M[2][l2] * 100 / l; }
        if (compMode.Contains("cov", StringComparison.Ordinal)) { l = Math.Max(l1, l2); return l == 0 ? 0 : M[4][l2] * 100 / l; }
        l = M[4][l2];
        return l == 0 ? 0 : M[2][l2] * 100 / l;
    }

    /// <summary>
    /// The number of glibc <c>rand()</c> values that precede the first pair of <c>-distance_matrix_mode random</c> in the stream started by the program's <c>srand</c> (<c>vsrand(0)</c>, <c>util.c</c>): 19, the same
    /// in every run. An <c>LD_PRELOAD</c> counter on the installed t_coffee shows 20 start-up calls before the first pair for 2-7 sequences and every option of the oracle (the 15-character temporary-file name of
    /// <c>vtmpnamH</c>, the name root of <c>get_vtmpnam2_root</c>, the header tip, ...); the first of them runs before <c>srand</c> and is not part of the seeded stream, which is what the fit on the 310
    /// committed sets shows (the pairs equal the stream from its 20th value on, for every skip 19 and no other constant).
    /// </summary>
    private const int TcRandomDrawsBeforeDistanceMatrix = 19;

    /// <summary>The distance matrix of the program (<c>seq2distance_matrix</c>) for the chosen <c>-distance_matrix_mode</c>.</summary>
    private static TcDistances TcDistanceMatrix(IReadOnlyList<string> seqs, TcSeqKind kind, TCoffeeProgramOptions o, TcRelaxedState? relaxed = null)
    {
        switch (o.DistanceMatrixMode)
        {
            case TCoffeeDistanceMode.Aln:
                throw new TCoffeeAbortedException("segmentation fault in aln2sim_mat: -distance_matrix_mode aln needs the alignment that cl2distance_matrix is called without (NOALN)", 11);
            case TCoffeeDistanceMode.Ktup2:
                throw new TCoffeeAbortedException("ERROR: UNSPECIFIED UNSPECIFIED (-distance_matrix_mode ktup2: seq2distance_matrix aligns with very_fast_aln and calls itself with the mode aln)");
            case TCoffeeDistanceMode.Other:
                throw new TCoffeeAbortedException($"Error: {o.DistanceModeText} is an unknown distance_matrix_mode [FATAL:T-COFFEE]");
            case TCoffeeDistanceMode.Slow or TCoffeeDistanceMode.Cscore or TCoffeeDistanceMode.Clustalw or TCoffeeDistanceMode.VeryFast:
                if (relaxed == null) throw new InvalidOperationException($"-distance_matrix_mode {o.DistanceMatrixMode.ToString().ToLowerInvariant()} needs the relaxed library of the progressive stage.");
                return TcLibraryDistanceMatrix(seqs, kind, o, relaxed);
            case TCoffeeDistanceMode.Geometric:
                return TcGeometricDistanceMatrix(seqs);
            case TCoffeeDistanceMode.Random:
                break;      // rand() % 100 per pair in the loop below
        }

        int n = seqs.Count;
        var dm = new TcDistances { Sim = new int[n, n], Score = new int[n, n] };
        bool protein = kind == TcSeqKind.Protein;
        if (o.DistanceMatrixMode == TCoffeeDistanceMode.Ktup)
        {
            var ktup = TCoffeeKtupSimilarity(seqs, protein);
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++)
                {
                    dm.Sim[a, b] = ktup[a, b];
                    dm.Score[a, b] = a == b ? 0 : ktup[a, b] * TcScoreK;
                }

            return dm;
        }

        TcMatrix mat = protein ? TcMatrix.Blosum62 : TcMatrix.IdMat;
        int gop = protein ? (int)(mat.AverageMismatch() * 10) : -10;
        int call = 0;
        // -distance_matrix_mode random: seq2distance_matrix draws `id = rand() % 100` for every pair a < b in row order (util_constraints_list.c); the program's glibc rand is seeded by vsrand(0)
        var random = o.DistanceMatrixMode == TCoffeeDistanceMode.Random ? new GlibcRandom(o.RandomSeed) : null;
        if (random != null) for (int skip = 0; skip < TcRandomDrawsBeforeDistanceMatrix; skip++) random.Next();
        for (int a = 0; a < n; a++)
        {
            dm.Sim[a, a] = TcMaxIdentity;
            for (int b = a + 1; b < n; b++, call++)
            {
                int id, score;
                switch (o.DistanceMatrixMode)
                {
                    case TCoffeeDistanceMode.Fast:
                        {
                            var (ra, rb) = TCoffeeSlowPairAlign(seqs[a], seqs[b], protein);
                            id = TcSeqSimMeasure(ra, rb, o.DistanceMatrixSimMode, call);
                            score = id * TcScoreK;
                            break;
                        }

                    case TCoffeeDistanceMode.Idscore:
                        id = score = TcIdscorePair(seqs[a], seqs[b], gop, -1, mat, o.DistanceMatrixSimMode);
                        break;
                    case TCoffeeDistanceMode.No:
                        id = TcMaxIdentity;
                        score = id * TcScoreK;
                        break;
                    case TCoffeeDistanceMode.Random:
                        id = random!.Next() % 100;
                        score = id * TcScoreK;
                        break;
                    default:
                        throw new NotSupportedException($"distance_matrix_mode {o.DistanceMatrixMode} is not supported.");
                }

                dm.Sim[a, b] = dm.Sim[b, a] = id;
                dm.Score[a, b] = dm.Score[b, a] = score;
            }
        }

        return dm;
    }

    /// <summary>
    /// <c>seq2distance_matrix</c> for the modes that use the (relaxed) library: <c>slow</c> (the DP score of the linked pair-wise alignment of every pair divided by its length, rescaled by <c>do_normalise</c>),
    /// <c>cscore</c> (the same score without the rescaling, float; falls back to <c>idscore</c> on an empty library), <c>clustalw</c> (the DP score itself) and <c>very_fast</c> (<c>very_fast_gotoh_pair_wise</c>:
    /// the diagonals are scored with the library entries, the 20 best with a shadow of 5, then the banded DP; the identity of the aligned rows).
    /// </summary>
    private static TcDistances TcLibraryDistanceMatrix(IReadOnlyList<string> seqs, TcSeqKind kind, TCoffeeProgramOptions o, TcRelaxedState st)
    {
        int n = seqs.Count;
        var mode = o.DistanceMatrixMode;
        if (mode == TCoffeeDistanceMode.Cscore && st.Index.EdgeCount == 0)
            return TcDistanceMatrix(seqs, kind, o with { DistanceMatrixMode = TCoffeeDistanceMode.Idscore });

        var dm = new TcDistances { Sim = new int[n, n], Score = new int[n, n] };
        bool protein = kind == TcSeqKind.Protein;
        TcMatrix mat = protein ? TcMatrix.Blosum62 : TcMatrix.IdMat;
        int gop = protein ? (int)(mat.AverageMismatch() * 10) : -10;
        int call = 0;
        for (int a = 0; a < n; a++)
        {
            dm.Sim[a, a] = TcMaxIdentity;
            for (int b = a + 1; b < n; b++, call++)
            {
                int id, score;
                switch (mode)
                {
                    case TCoffeeDistanceMode.VeryFast:
                        {
                            var (ra, rb) = TcVeryFastPairRows(seqs[a], seqs[b], a, b, st, mat, gop);
                            id = TcSeqSimMeasure(ra, rb, o.DistanceMatrixSimMode, call);
                            score = id * TcScoreK;
                            break;
                        }

                    case TCoffeeDistanceMode.Slow:
                        {
                            var (ra, rb, scoreAln) = TcPairRowsOfMainDp(seqs, a, b, st);
                            id = TcSeqSimMeasure(ra, rb, o.DistanceMatrixSimMode, call);
                            float s = (int)((float)scoreAln / (ra.Length * TcScoreK));
                            if (st.Normalise != 0) s = s * TcMaxIdentity / st.Normalise;
                            if (s > TcMaxIdentity) throw new NotSupportedException("-distance_matrix_mode slow: the pair score exceeds 100 (the program then recomputes it with sub_aln2sub_aln_score).");
                            score = (int)s;
                            break;
                        }

                    case TCoffeeDistanceMode.Cscore:
                        {
                            var (ra, _, scoreAln) = TcLinkedPairRows(seqs, a, b, st.Index, st.Penalties);
                            float s = scoreAln;
                            s /= ra.Length * TcScoreK;
                            float fid = s / TcScoreK;
                            id = (int)fid;
                            score = (int)s;
                            break;
                        }

                    default: // clustalw: id = the DP score of the pair alignment
                        {
                            var (_, _, scoreAln) = TcPairRowsOfMainDp(seqs, a, b, st);
                            float fid = scoreAln;
                            id = (int)fid;
                            score = (int)(fid * TcScoreK);
                            break;
                        }
                }

                dm.Sim[a, b] = dm.Sim[b, a] = id;
                dm.Score[a, b] = dm.Score[b, a] = score;
            }
        }

        return dm;
    }

    /// <summary><c>pair_wise</c> with the main dp mode of the program (<c>-dp_mode</c>, default the linked sparse DP) for the pair (a, b) at the distance matrix stage.</summary>
    private static (string RowA, string RowB, int Score) TcPairRowsOfMainDp(IReadOnlyList<string> seqs, int a, int b, TcRelaxedState st)
    {
        if (st.Dp is null or { Mode: TCoffeeDpMode.Linked }) return TcLinkedPairRows(seqs, a, b, st.Index, st.Penalties);
        return TcDpPairRows(seqs, a, b, new TcDpContext(st.Index, null, st.Dp, st.Penalties.NoMatch, (int)Math.Min(int.MaxValue, st.Index.MaxExtValue()), seqs));
    }

    /// <summary>
    /// <c>very_fast_gotoh_pair_wise</c> of the pair (a, b) with the relaxed library: <c>evaluate_diagonals_with_clist</c> (every library entry a -> b adds the matrix cost of its residue pair to its diagonal),
    /// <c>extract_N_diag(20, window 5)</c>, the banded Gotoh DP with the matrix cost (idmat / blosum62mt x SCORE_K minus SCORE_K x nomatch), GOP / GEP of the distance mode and the terminal gap rule <c>-tg_mode</c>.
    /// </summary>
    private static (string, string) TcVeryFastPairRows(string sa, string sb, int a, int b, TcRelaxedState st, TcMatrix mat, int gop)
    {
        int l1 = sa.Length, l2 = sb.Length, nDiag = l1 + l2 - 1;
        RequireTcPairDpCells("very_fast_pair_wise (evaluate_diagonals_with_clist)", "keeps a 3-int row (about 48 bytes) per diagonal (len1 + len2 + 1 of them)", 48, (long)l1 + l2 + 1,
            string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{l1:N0} x {l2:N0} positions"));
        Func<int, int, int> cost = (i, j) => mat.Score(sa[i], sb[j]) * TcScoreK - TcScoreK * st.Penalties.NoMatch;
        var diag = new int[nDiag + 2][];
        for (int d = 0; d < nDiag + 2; d++) diag[d] = new[] { d, 0, 0 };
        for (int r1 = 1; r1 <= l1; r1++)
            foreach (var e in st.Index.Edges[a][r1])
                if (e.Seq == b)
                    diag[e.Res - r1 + l1][1] += cost(r1 - 1, e.Res - 1);
        TcSortDiagonals(diag, nDiag);
        int[] list = TcExtractNDiagonals(l1, l2, diag, 20, 5);
        var (m0, m1, _) = TcMakeFastaGotoh(l1, l2, list, cost, gop * TcScoreK, -1 * TcScoreK, st.TgMode);
        return TcMasksToRows(sa, sb, m0, m1);
    }

    /// <summary>
    /// <c>geometric</c>: every sequence is aligned (<c>fasta_pair_wise</c>, pam250mt, gop -10, gep -1, k-tuple 3, terminal gaps free, lower-case letters) with the first three sequences; the identity of those alignments is the
    /// coordinate vector of the sequence, and the similarity of a pair is <c>100 x (1 - euclidean distance / sqrt(10000 x n_coor))</c> (single precision, as in the program).
    /// </summary>
    private static TcDistances TcGeometricDistanceMatrix(IReadOnlyList<string> seqs)
    {
        int n = seqs.Count;
        int nCoor = Math.Min(3, n);
        var pam = TcMatrix.ByName("pam250mt")!;
        var g = new float[n, 3];
        for (int a = 0; a < n; a++)
            for (int b = 0; b < nCoor; b++)
            {
                string s1 = seqs[a].ToLowerInvariant(), s2 = seqs[b].ToLowerInvariant();
                var (ra, rb) = TcFastaPairRows(s1, s2, pam, -10, -1, 3, 0, 0, 0f, 2);
                g[a, b] = TcSeqSim(ra, rb);
            }

        float reference = (float)Math.Sqrt(10000 * nCoor);
        var dm = new TcDistances { Sim = new int[n, n], Score = new int[n, n] };
        for (int a = 0; a < n; a++)
        {
            dm.Sim[a, a] = TcMaxIdentity;
            for (int b = a + 1; b < n; b++)
            {
                float t = 0;
                for (int c = 0; c < nCoor; c++)
                {
                    float d = g[a, c] - g[b, c];
                    t += d * d;
                }

                float distance = (float)Math.Sqrt(t);
                float id = TcMaxIdentity * (1 - (distance / reference));
                dm.Sim[a, b] = dm.Sim[b, a] = (int)id;
                dm.Score[a, b] = dm.Score[b, a] = (int)id * TcScoreK;
            }
        }

        return dm;
    }

    /// <summary>make_tree: <c>normalize_array(sim_array2dist_array(score, 1000), 1000, 100)</c> followed by NJ (+ reroot) or UPGMA.</summary>
    private static TCoffeeTreeNode TcGuideTreeFor(IReadOnlyList<string> seqs, TcSeqKind kind, TCoffeeProgramOptions o, TcRelaxedState? relaxed = null)
    {
        int k = seqs.Count;
        if (k == 1) return TCoffeeTreeNode.Leaf(0);
        if (k == 2) return TCoffeeTreeNode.Join(TCoffeeTreeNode.Leaf(0), TCoffeeTreeNode.Leaf(1));
        switch (o.TreeMode)
        {
            case TCoffeeTreeMode.Cwph or TCoffeeTreeMode.Other:
                throw new TCoffeeAbortedException($"COREDUMP in make_tree: ERROR: {(o.TreeMode == TCoffeeTreeMode.Cwph ? "cwph" : o.TreeModeText)} is an unknown tree computation mode [FATAL:T-COFFEE] (the dispatcher knows only nj and upgma)", 11);
            case TCoffeeTreeMode.Kmeans:
                throw new TCoffeeAbortedException("-tree_mode kmeans: seq2km_tree_old builds the k-means tree and the program ends without writing an alignment", 0);
        }

        var dm = TcDistanceMatrix(seqs, kind, o, relaxed);
        var dist = new int[k, k];
        for (int a = 0; a < k; a++)
            for (int b = 0; b < k; b++)
                dist[a, b] = (TcMaxIdentity * TcScoreK - dm.Score[a, b]) * 100 / (TcMaxIdentity * TcScoreK);

        if (o.TreeMode == TCoffeeTreeMode.Upgma) return TcUpgmaTree(dist, k);
        var d = new double[k, k];
        for (int a = 0; a < k; a++)
            for (int b = 0; b < k; b++)
                d[a, b] = a == b ? 0 : (double)dist[a, b] / TcMaxIdentity;
        var unrooted = k < 100 ? TcNeighborJoining(d, k) : TcNeighborJoiningFast(d, k);
        return ToTCoffeeTree(ToGuideNode(ClustalWRootNodes(unrooted, k).Root));
    }

    /// <summary><c>int_dist2upgma_tree</c>: the pair with the smallest integer distance (first in row-major order on ties) is merged, the new row is the integer mean.</summary>
    private static TCoffeeTreeNode TcUpgmaTree(int[,] mat, int n)
    {
        var nodes = new TCoffeeTreeNode[n];
        for (int a = 0; a < n; a++) nodes[a] = TCoffeeTreeNode.Leaf(a);
        var used = new bool[n];
        TCoffeeTreeNode root = nodes[0];
        for (int left = n; left > 1; left--)
        {
            int bestA = 0, bestB = 0, best = 0;
            bool set = false;
            for (int a = 0; a < n - 1; a++)
            {
                if (used[a]) continue;
                for (int b = a + 1; b < n; b++)
                {
                    if (used[b]) continue;
                    int w = mat[a, b];
                    if (!set || w < best) { best = w; bestA = a; bestB = b; set = true; }
                }
            }

            for (int a = 0; a < n; a++)
                if (!used[a]) mat[bestA, a] = mat[a, bestA] = (mat[bestB, a] + mat[bestA, a]) / 2;
            used[bestB] = true;
            // read_tree turns the written (left, right) into (right, left): the walk of the alignment starts with the merged-in (higher index) subtree
            root = TCoffeeTreeNode.Join(nodes[bestB], nodes[bestA]);
            nodes[bestA] = root;
        }

        return root;
    }

    /// <summary><c>compute_t_coffee_weight</c> on the identities of the distance matrix (single-precision arithmetic as in the original).</summary>
    private static float[] TcTCoffeeSequenceWeights(int[,] sim, int n)
    {
        var w = new float[n];
        Array.Fill(w, 1f);
        for (int a = 0; a < n - 1; a++)
            for (int b = a + 1; b < n; b++)
            {
                float d = (float)sim[a, b] / TcMaxIdentity * 10;
                float p = (float)Math.Pow(d, 3);
                w[a] += p;
                w[b] += p;
            }

        float total = 0;
        for (int b = 0; b < n; b++)
        {
            w[b] = 2 / w[b];
            total += w[b];
        }

        for (int b = 0; b < n; b++) w[b] = w[b] * ((float)n / total);
        return w;
    }
}
