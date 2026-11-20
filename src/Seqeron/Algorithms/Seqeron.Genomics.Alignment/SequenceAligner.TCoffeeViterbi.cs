namespace Seqeron.Genomics.Alignment;

// =====================================================================================================================
// -dp_mode viterbi_pair_wise (util_dp_suboptimal_nw.c L1788, get_tot_prob2 L1184, evaluate_dirichlet.c aln2prf L651): the Viterbi path of the five-state pair-HMM of proba_pair (the PROTEIN model: 5 matrix types,
// 2 insert state pairs, ProbCons default tables, whatever the sequence type) between the residue FREQUENCY profiles of the two groups (aln2prf: the share of each of the 20 dirichlet letters among the residues of a
// column; a column without such a letter has an all-zero profile). Single precision as in the original (float tables, float products and sums).
// =====================================================================================================================

public static partial class SequenceAligner
{
    private static float[][] TcResidueProfile(TcDpContext ctx, TcGroup g)
    {
        int len = g.Columns;
        var prf = new float[len][];
        for (int c = 0; c < len; c++)
        {
            var counts = new double[20];
            double tot = 0;
            foreach (int s in g.Order)
            {
                int r = g.Rows[s][c];
                if (r <= 0) continue;
                int d = TcDirichletCode(ctx.Sequences[s][r - 1]);
                if (d < 0) continue;
                counts[d]++;
                tot++;
            }

            prf[c] = new float[20];
            if (tot > 0)
                for (int r = 0; r < 20; r++) prf[c][r] = (float)(counts[r] / tot);
        }

        return prf;
    }

    /// <summary>The limit of <c>viterbi_pair_wise</c>: the match probability (4 bytes) plus five Viterbi floats and five trace-back ints (40 bytes) per cell, 64-bit (the program computes <c>nmt x (I + 1) x (J + 1)</c> in a 32-bit int).</summary>
    internal static void TcRequireViterbiFits(int columns1, int columns2) =>
        RequireTcPairDp("viterbi_pair_wise", "keeps the match probability of every cell and five Viterbi floats and five trace-back ints per cell", 44, (long)columns1 + 1, (long)columns2 + 1, "columns");

    private static (List<bool> Present1, List<bool> Present2, int Score) TcViterbiMasks(TcDpContext ctx, TcGroup g0, TcGroup g1)
    {
        var m = TcProbaModel.For(TcSeqKind.Protein);
        const int nmt = 5, nis = 2;
        int I = g0.Columns, J = g1.Columns;
        TcRequireViterbiFits(I, J);
        var prf1 = TcResidueProfile(ctx, g0);
        var prf2 = TcResidueProfile(ctx, g1);

        // get_tot_prob2 (profile branch)
        var insProb = new float[2][,];
        insProb[0] = new float[nmt, I + 1];
        insProb[1] = new float[nmt, J + 1];
        for (int i = 1; i <= I; i++)
            for (int k = 0; k < nmt; k++)
            {
                float v = 0;
                for (int r = 0; r < 20; r++) v += prf1[i - 1][r] * m.InsProb[TcDirichletLetters[r], k];
                insProb[0][k, i] = v;
            }

        for (int j = 1; j <= J; j++)
            for (int k = 0; k < nmt; k++)
            {
                float v = 0;
                for (int r = 0; r < 20; r++) v += prf2[j - 1][r] * m.InsProb[TcDirichletLetters[r], k];
                insProb[1][k, j] = v;
            }

        var matchProb = new float[(I + 1) * (J + 1)];
        for (int i = 1; i <= I; i++)
        {
            AlignmentCancellation.ThrowIfRequested();       // A5-2: one poll per row (400 multiplications per cell)
            for (int j = 1; j <= J; j++)
            {
                float v = 0;
                for (int r1 = 0; r1 < 20; r1++)
                    for (int r2 = 0; r2 < 20; r2++)
                    {
                        float f = prf1[i - 1][r1] * prf2[j - 1][r2];
                        v += m.MatchProb[TcDirichletLetters[r1], TcDirichletLetters[r2]] * f;
                    }

                matchProb[i * (J + 1) + j] = v;
            }
        }

        int total = nmt * (I + 1) * (J + 1);
        var viterbi = new float[total];
        var traceback = new int[total];
        Array.Fill(viterbi, TcLogZero);
        Array.Fill(traceback, -1);
        for (int k = 0; k < nmt; k++) viterbi[k] = m.InitialDistribution[k];

        for (int i = 0; i <= I; i++)
        {
            AlignmentCancellation.ThrowIfRequested();       // A5-2: one poll per row
            for (int j = 0; j <= J; j++)
            {
                int ij = nmt * (i * (J + 1) + j), i1j = ij - nmt * (J + 1), ij1 = ij - nmt, i1j1 = ij - nmt * (J + 2);
                if (i > 0 && j > 0)
                {
                    for (int k = 0; k < nmt; k++)
                    {
                        float newVal = viterbi[k + i1j1] + m.TransProb[k, 0] + matchProb[i * (J + 1) + j];
                        if (viterbi[0 + ij] < newVal)
                        {
                            viterbi[0 + ij] = newVal;
                            traceback[0 + ij] = k;
                        }
                    }
                }

                if (i > 0)
                {
                    for (int k = 0; k < nis; k++)
                    {
                        float fromMatch = insProb[0][k, i] + viterbi[0 + i1j] + m.TransProb[0, 2 * k + 1];
                        float fromIns = insProb[0][k, i] + viterbi[2 * k + 1 + i1j] + m.TransProb[2 * k + 1, 2 * k + 1];
                        if (fromMatch >= fromIns) { viterbi[2 * k + 1 + ij] = fromMatch; traceback[2 * k + 1 + ij] = 0; }
                        else { viterbi[2 * k + 1 + ij] = fromIns; traceback[2 * k + 1 + ij] = 2 * k + 1; }
                    }
                }

                if (j > 0)
                {
                    for (int k = 0; k < nis; k++)
                    {
                        float fromMatch = insProb[1][k, j] + viterbi[0 + ij1] + m.TransProb[0, 2 * k + 2];
                        float fromIns = insProb[1][k, j] + viterbi[2 * k + 2 + ij1] + m.TransProb[2 * k + 2, 2 * k + 2];
                        if (fromMatch >= fromIns) { viterbi[2 * k + 2 + ij] = fromMatch; traceback[2 * k + 2 + ij] = 0; }
                        else { viterbi[2 * k + 2 + ij] = fromIns; traceback[2 * k + 2 + ij] = 2 * k + 2; }
                    }
                }
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

        var al0 = new List<bool>(); var al1 = new List<bool>();
        int rr = I, cc = J;
        while (rr != 0 || cc != 0)
        {
            int newState = traceback[state + nmt * (rr * (J + 1) + cc)];
            if (state == 0) { cc--; rr--; al0.Add(true); al1.Add(true); }
            else if (state % 2 == 1) { rr--; al0.Add(true); al1.Add(false); }
            else { cc--; al0.Add(false); al1.Add(true); }
            state = newState;
        }

        al0.Reverse(); al1.Reverse();
        return (al0, al1, (int)(bestProb * 1000f));
    }
}
