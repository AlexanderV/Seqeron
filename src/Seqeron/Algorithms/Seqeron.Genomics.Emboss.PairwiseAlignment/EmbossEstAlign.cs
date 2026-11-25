// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.PairwiseAlignment;

using Seqeron.Genomics.Alignment;

/// <summary>The result of an EST / genomic alignment (<c>EmbOEstAlign</c> of <c>embest.h</c>): the 0-based start and stop in the genomic and the EST sequence, the score and the path
/// (<c>1</c> diagonal, <c>2</c> deletion in the EST, <c>3</c> deletion in the genome, <c>0</c> / <c>-1</c> / <c>-2</c> an intron followed by its length).</summary>
internal sealed class EstAlign
{
    public int GStart, EStart, GStop, EStop, Score, Len;
    public int[] AlignPath = [];
}

/// <summary>The donor / acceptor marks of a genomic sequence (<c>embEstFindSpliceSites</c>) with the direction of the gene it was built for (<c>forward</c> = GT/AG, <c>reverse</c> = CT/AC).</summary>
internal sealed class EstSpliceSites(string sites, string name)
{
    public string Sites { get; } = sites;

    public string Name { get; } = name;
}

/// <summary>
/// The spliced-alignment engine of EMBOSS 6.6.0 <c>est2genome</c> (<c>nucleus/embest.c</c>, Richard Mott's <c>est_genome</c>, Sanger Centre), ported statement by statement: the affine-free local / global dynamic program with
/// introns (a free-length gap that costs <c>intron_penalty</c>, or <c>splice_penalty</c> between a donor and an acceptor site), the packed two-bit back-pointers, the linear-space divide and conquer
/// (<c>estAlignMidpt</c>, <c>estAlignRecursive</c>) used when the matrix would not fit into the space limit, the table of remembered intron starts (<c>estPairRemember</c>) and the three output writers
/// (<c>embEstOutBlastStyle</c>, <c>estWriteMsp</c>, <c>embEstPrintAlign</c>). Scores are 32-bit integers, the percentages are <c>float</c>s.
/// </summary>
/// <remarks>
/// The shuffle of the program (<c>embEstShuffleSeq</c>) draws the random keys but never sorts them (the sort is missing in the 6.6.0 source: measured, the sequence is unchanged), so a shuffled EST is the EST itself and the
/// random numbers have no effect on the output; they are not generated here.
/// </remarks>
internal sealed class EstEngine
{
    public const int Intron = 0, Diagonal = 1, DeleteEst = 2, DeleteGenome = 3, ForwardSplicedIntron = -1, ReverseSplicedIntron = -2;
    private const int NotASite = 1, Donor = 2, Acceptor = 4;
    private const int MinusInfinity = -10000000;

    private readonly int[,] _sim = new int[256, 256];

    // the remembered (column, row) pairs of estDoNotForget / estPairRemember
    private readonly List<(int Col, int Row)> _pairs = [];
    private bool _pairsSorted;

    /// <summary><c>embEstMatInit (match, mismatch, gap, 0, '-')</c>.</summary>
    public EstEngine(int match, int mismatch, int gap)
    {
        const int neutral = 0;
        const int pad = '-';
        for (int c1 = 0; c1 < 256; c1++)
        {
            for (int c2 = 0; c2 < 256; c2++)
            {
                if (c1 == c2)
                {
                    _sim[c1, c2] = c1 != '*' && c1 != 'n' && c1 != 'N' && c1 != '-' ? match : 0;
                }
                else if (c1 == pad || c2 == pad)
                {
                    _sim[c1, c2] = _sim[c2, c1] = -gap;
                }
                else if (c1 == 'n' || c2 == 'n' || c1 == 'N' || c2 == 'N')
                {
                    _sim[c1, c2] = _sim[c2, c1] = neutral;
                }
                else
                {
                    _sim[c1, c2] = _sim[c2, c1] = -mismatch;
                }
            }
        }

        for (int c1 = 0; c1 < 256; c1++)
        {
            int c2 = c1 >= 'A' && c1 <= 'Z' ? c1 + 32 : c1;
            _sim[c1, c2] = _sim[c1, c1];
            _sim[c2, c1] = _sim[c1, c1];
        }
    }

    /// <summary>The substitution score of two bases (<c>lsimmat</c>).</summary>
    public int Sim(char a, char b) => _sim[a & 0xFF, b & 0xFF];

    /// <summary><c>embEstFindSpliceSites</c>: <c>forward</c> marks GT donors / AG acceptors, otherwise CT donors / AC acceptors.</summary>
    public static EstSpliceSites FindSpliceSites(string genome, bool forward)
    {
        int n = genome.Length;
        var sites = new char[n];
        Array.Fill(sites, (char)NotASite);
        char Lower(int i) => char.ToLowerInvariant(genome[i]);
        char donor1 = forward ? 'g' : 'c', donor2 = 't', accept1 = 'a', accept2 = forward ? 'g' : 'c';
        for (int pos = 1; pos < n - 2; pos++)
        {
            if (Lower(pos) == donor1 && Lower(pos + 1) == donor2)
                sites[pos - 1] = (char)(sites[pos - 1] | Donor);
            if (Lower(pos) == accept1 && Lower(pos + 1) == accept2)
                sites[pos + 1] = (char)(sites[pos + 1] | Acceptor);
        }

        return new EstSpliceSites(new string(sites), forward ? "forward" : "reverse");
    }

    private static string Lower(string s) => s.ToLowerInvariant();

    /// <summary><c>ajSeqTrim</c> of the range (<paramref name="begin"/>, <paramref name="end"/>) given as 1-based positions (<c>ajCvtSposToPosStart</c>: both ends clipped into the sequence, an end before the start leaves the sequence untrimmed at the end).</summary>
    public static string Trim(string s, int begin, int end)
    {
        int n = s.Length;
        int jbegin = begin > 0 ? begin - 1 : begin;
        int jend = end > 0 ? end - 1 : end;
        int b = 1 + PositionStart(n, 0, jbegin);
        int e = 1 + PositionStart(n, b - 1, jend);
        string residues = s;
        bool trimmed = false;
        if (end != 0)
        {
            if (e >= b)
            {
                residues = residues[..e];
                trimmed = true;
            }
            else
            {
                begin = 0;
            }
        }

        if (begin != 0 && !(end != 0 && !trimmed))
            residues = residues[(b - 1)..];
        return residues;
    }

    private static int PositionStart(int length, int minimum, long position)
    {
        long j = position < 0 ? length + position : position;
        if (j >= length)
            j = length - 1;
        if (j < minimum)
            j = minimum;
        return (int)j;
    }

    private void PairInit() => PairFree();

    private void PairFree()
    {
        _pairs.Clear();
        _pairsSorted = false;
    }

    private static int PairCompare((int Col, int Row) a, (int Col, int Row) b)
    {
        int n = a.Col - b.Col;
        return n == 0 ? a.Row - b.Row : n;
    }

    private void DoNotForget(int col, int row)
    {
        _pairs.Add((col, row));
        _pairsSorted = false;
    }

    private int PairRemember(int col, int row)
    {
        if (!_pairsSorted)
        {
            // qsort of the pairs (equal pairs are interchangeable)
            _pairs.Sort(PairCompare);
            _pairsSorted = true;
        }

        var rp = (col, row);
        int left = 0, right = _pairs.Count - 1;
        int bad;
        if (_pairs.Count == 0)
            throw new InvalidOperationException("ERROR in estPairRemember(): no remembered pair");
        if (PairCompare(_pairs[left], rp) > 0 || PairCompare(_pairs[right], rp) < 0)
        {
            bad = 1;
        }
        else
        {
            bad = 0;
            while (right - left > 1)
            {
                int middle = (left + right) / 2;
                int d = PairCompare(_pairs[middle], rp);
                if (d < 0)
                    left = middle;
                else
                    right = middle;
            }

            if (PairCompare(_pairs[right], rp) < 0 || _pairs[left].Col != col || _pairs[left].Row >= row)
                bad = 2;
        }

        if (bad != 0)
            throw new InvalidOperationException("ERROR in estPairRemember() left: " + left.ToString(CultureInfo.InvariantCulture) + " col: " + col.ToString(CultureInfo.InvariantCulture) + " row: " + row.ToString(CultureInfo.InvariantCulture) + ", bad: " + bad.ToString(CultureInfo.InvariantCulture));
        return _pairs[left].Row;
    }

    /// <summary><c>embEstAlignNonRecursive</c>: the dynamic program over the whole matrix; <paramref name="backtrack"/> keeps the two-bit back-pointers and returns the path, otherwise the start of the best local alignment is carried along with the scores.</summary>
    public EstAlign AlignNonRecursive(string est, string genome, int gapPenalty, int intronPenalty, int splicePenalty, EstSpliceSites? spliceSites, bool backtrack, bool needleman, int initPath)
    {
        int glen = genome.Length, elen = est.Length;
        string ss = spliceSites?.Sites ?? string.Empty;
        int eLenPack = (elen / 4) + 1;
        var ge = new EstAlign();
        byte[][]? ppath = null;
        int[]? tempPath = null;
        int[] startLeft1 = [], startRight1 = [], startLeft2 = [], startRight2 = [], bestIntronStartLeft = [], bestIntronStartRight = [];
        if (backtrack)
        {
            ppath = new byte[glen][];
            for (int g = 0; g < glen; g++)
                ppath[g] = new byte[eLenPack];
            tempPath = new int[glen + elen];
        }
        else
        {
            startLeft1 = new int[elen + 1];
            startRight1 = new int[elen + 1];
            startLeft2 = new int[elen + 1];
            startRight2 = new int[elen + 1];
            bestIntronStartLeft = new int[elen];
            bestIntronStartRight = new int[elen];
        }

        // t1 / t2 are start1 + 1 / start2 + 1: element -1 is index 0
        int[] t1L = startLeft1, t1R = startRight1, t2L = startLeft2, t2R = startRight2;
        int[] score1 = new int[elen + 1], score2 = new int[elen + 1];
        int[] s1 = score1, s2 = score2;
        var bestIntronCoord = new int[elen + 1];
        var bestIntronScore = new int[elen + 1];
        string gseq = Lower(genome), eseq = Lower(est);
        int emax = -1, gmax = -1, maxScore = 0;
        int bestStartLeft = 0, bestStartRight = 0;

        if (!backtrack)
        {
            for (int epos = 0; epos < elen; epos++)
            {
                t1L[epos + 1] = 0;
                t1R[epos + 1] = epos;
                t2L[epos + 1] = 0;
                t2R[epos + 1] = epos;
                bestIntronStartLeft[epos] = t1L[epos + 1];
                bestIntronStartRight[epos] = t1R[epos + 1];
            }
        }

        if (needleman)
        {
            for (int epos = 0; epos < elen; epos++)
            {
                s1[epos + 1] = MinusInfinity;
                bestIntronScore[epos] = MinusInfinity;
            }
        }

        for (int gpos = 0; gpos < glen; gpos++)
        {
            (s1, s2) = (s2, s1);
            char g = gseq[gpos];
            byte[]? path = null;
            if (backtrack)
            {
                path = ppath![gpos];
            }
            else
            {
                (t1L, t2L) = (t2L, t1L);
                (t1R, t2R) = (t2R, t1R);
                t2L[0] = gpos;
                t1L[0] = gpos;
                t1R[0] = 0;
            }

            bool isAcceptor = spliceSites != null && (ss[gpos] & Acceptor) != 0;
            if (needleman)
                s1[0] = initPath == Diagonal || gpos > 0 ? MinusInfinity : 0;
            else
                s1[0] = 0;

            for (int epos = 0; epos < elen; epos++)
            {
                int diagonal = unchecked(s2[epos] + Sim(g, eseq[epos]));
                int deleteEst = unchecked(s1[epos] - gapPenalty);
                int deleteGenome = unchecked(s2[epos + 1] - gapPenalty);
                int intron = isAcceptor && (ss[bestIntronCoord[epos]] & Donor) != 0
                    ? unchecked(bestIntronScore[epos] - splicePenalty)
                    : unchecked(bestIntronScore[epos] - intronPenalty);
                int max = deleteEst > deleteGenome ? deleteEst : deleteGenome;
                if (diagonal > max)
                    max = diagonal;
                if (intron > max)
                    max = intron;

                if (needleman || max > 0)
                {
                    if (max == diagonal)
                    {
                        s1[epos + 1] = diagonal;
                        if (backtrack)
                        {
                            path![epos / 4] = (byte)(path[epos / 4] | (1 << (2 * (epos % 4))));
                        }
                        else if (t2L[epos] == -1)
                        {
                            t1L[epos + 1] = gpos;
                            t1R[epos + 1] = epos;
                        }
                        else
                        {
                            t1L[epos + 1] = t2L[epos];
                            t1R[epos + 1] = t2R[epos];
                        }
                    }
                    else if (max == deleteEst)
                    {
                        s1[epos + 1] = deleteEst;
                        if (backtrack)
                        {
                            path![epos / 4] = (byte)(path[epos / 4] | (2 << (2 * (epos % 4))));
                        }
                        else
                        {
                            t1L[epos + 1] = t1L[epos];
                            t1R[epos + 1] = t1R[epos];
                        }
                    }
                    else if (max == deleteGenome)
                    {
                        s1[epos + 1] = deleteGenome;
                        if (backtrack)
                        {
                            path![epos / 4] = (byte)(path[epos / 4] | (3 << (2 * (epos % 4))));
                        }
                        else
                        {
                            t1L[epos + 1] = t2L[epos + 1];
                            t1R[epos + 1] = t2R[epos + 1];
                        }
                    }
                    else
                    {
                        s1[epos + 1] = intron;
                        if (!backtrack)
                        {
                            t1L[epos + 1] = bestIntronStartLeft[epos];
                            t1R[epos + 1] = bestIntronStartRight[epos];
                        }
                    }
                }
                else
                {
                    s1[epos + 1] = 0;
                    if (!backtrack)
                    {
                        t1L[epos + 1] = -1;
                        t1R[epos + 1] = -1;
                    }
                }

                if (bestIntronScore[epos] < s1[epos + 1])
                {
                    if (backtrack)
                        DoNotForget(epos, bestIntronCoord[epos]);
                    bestIntronScore[epos] = s1[epos + 1];
                    bestIntronCoord[epos] = gpos;
                    if (!backtrack)
                    {
                        bestIntronStartLeft[epos] = t1L[epos + 1];
                        bestIntronStartRight[epos] = t1R[epos + 1];
                    }
                }

                if (!needleman && maxScore < s1[epos + 1])
                {
                    maxScore = s1[epos + 1];
                    emax = epos;
                    gmax = gpos;
                    if (!backtrack)
                    {
                        bestStartLeft = t1L[epos + 1];
                        bestStartRight = t1R[epos + 1];
                    }
                }
            }
        }

        if (needleman)
        {
            ge.GStop = glen - 1;
            ge.EStop = elen - 1;
            ge.Score = s1[ge.EStop + 1];
        }
        else
        {
            ge.GStop = gmax;
            ge.EStop = emax;
            ge.Score = maxScore;
        }

        if (backtrack)
        {
            int pos = 0;
            int epos = ge.EStop;
            int gpos = ge.GStop;
            int total = 0;
            int spliceType = 0;
            if (spliceSites != null)
            {
                if (spliceSites.Name == "forward")
                    spliceType = ForwardSplicedIntron;
                else if (spliceSites.Name == "reverse")
                    spliceType = ReverseSplicedIntron;
            }
            while ((needleman || total < maxScore) && epos >= 0 && gpos >= 0)
            {
                int direction = (ppath![gpos][epos / 4] & (3 << (2 * (epos % 4)))) >> (2 * (epos % 4));
                tempPath![pos++] = direction;
                if (direction == Intron)
                {
                    int gpos1;
                    if (gpos - bestIntronCoord[epos] <= 0)
                        gpos1 = PairRemember(epos, gpos);
                    else
                        gpos1 = bestIntronCoord[epos];
                    if (spliceSites != null && (ss[gpos] & Acceptor) != 0 && (ss[gpos1] & Donor) != 0)
                    {
                        total -= splicePenalty;
                        tempPath[pos - 1] = spliceType;
                    }
                    else
                    {
                        total -= intronPenalty;
                    }

                    tempPath[pos++] = gpos - gpos1;
                    gpos = gpos1;
                }
                else if (direction == Diagonal)
                {
                    total += Sim(gseq[gpos], eseq[epos]);
                    epos--;
                    gpos--;
                }
                else if (direction == DeleteEst)
                {
                    total -= gapPenalty;
                    epos--;
                }
                else
                {
                    total -= gapPenalty;
                    gpos--;
                }
            }

            gpos++;
            epos++;
            ge.GStart = gpos;
            ge.EStart = epos;
            ge.Len = pos;
            ge.AlignPath = new int[ge.Len];
            int q = 0;
            while (q < ge.Len)
            {
                if (tempPath![q] > Intron)
                {
                    ge.AlignPath[pos - q - 1] = tempPath[q];
                }
                else
                {
                    ge.AlignPath[pos - q - 2] = tempPath[q];
                    ge.AlignPath[pos - q - 1] = tempPath[q + 1];
                    q++;
                }

                q++;
            }

            PairFree();
        }
        else
        {
            ge.GStart = bestStartLeft;
            ge.EStart = bestStartRight;
        }

        return ge;
    }

    /// <summary><c>embEstAlignLinearSpace</c>: the whole matrix when it fits into <paramref name="megabytes"/>, otherwise a local alignment to find the region and then the divide and conquer.</summary>
    public EstAlign AlignLinearSpace(string est, string genome, int match, int mismatch, int gapPenalty, int intronPenalty, int splicePenalty, EstSpliceSites? spliceSites, float megabytes)
    {
        int elen = est.Length, glen = genome.Length;
        float maxArea = megabytes * (float)1.0e6;
        PairInit();
        float area = ((float)glen + (float)1.0) * ((float)elen + (float)1.0) / 4;
        if (area <= maxArea)
            return AlignNonRecursive(est, genome, gapPenalty, intronPenalty, splicePenalty, spliceSites, true, false, Diagonal);

        var ge = AlignNonRecursive(est, genome, gapPenalty, intronPenalty, splicePenalty, spliceSites, false, false, Diagonal);
        string genomeSub = Trim(genome, ge.GStart + 1, ge.GStop + 1);
        string estSub = Trim(est, ge.EStart + 1, ge.EStop + 1);
        EstSpliceSites? spliceSub = spliceSites == null ? null : new EstSpliceSites(Trim(spliceSites.Sites, ge.GStart + 1, ge.GStop + 1), spliceSites.Name);
        var rge = AlignRecursive(estSub, genomeSub, match, mismatch, gapPenalty, intronPenalty, splicePenalty, spliceSub, maxArea);
        ge.Len = rge.Len;
        ge.AlignPath = rge.AlignPath;
        return ge;
    }

    private EstAlign AlignRecursive(string est, string genome, int match, int mismatch, int gapPenalty, int intronPenalty, int splicePenalty, EstSpliceSites? spliceSites, float maxArea)
    {
        int glen = genome.Length, elen = est.Length;
        float area = ((float)glen + (float)1.0) * ((float)elen + (float)1.0) / 4;
        if (area <= maxArea)
            return AlignNonRecursive(est, genome, gapPenalty, intronPenalty, splicePenalty, spliceSites, true, true, Diagonal);

        int middle = est.Length / 2;
        AlignMidpoint(est, genome, gapPenalty, intronPenalty, splicePenalty, spliceSites, middle, out int gleft, out int gright);
        bool splitOnDel = gleft == gright;
        string leftGenome = Trim(genome, 1, gleft + 1);
        string rightGenome = Trim(genome, gright + 1, genome.Length);
        EstSpliceSites? leftSplice = null, rightSplice = null;
        if (spliceSites != null)
        {
            leftSplice = new EstSpliceSites(Trim(spliceSites.Sites, 1, gleft + 1), spliceSites.Name);
            rightSplice = new EstSpliceSites(Trim(spliceSites.Sites, gright + 1, genome.Length), spliceSites.Name);
        }

        string leftEst = Trim(est, 1, middle + 1);
        string rightEst = Trim(est, middle + 2, est.Length);
        var leftGe = AlignRecursive(leftEst, leftGenome, match, mismatch, gapPenalty, intronPenalty, splicePenalty, leftSplice, maxArea);
        var rightGe = AlignRecursive(rightEst, rightGenome, match, mismatch, gapPenalty, intronPenalty, splicePenalty, rightSplice, maxArea);
        var ge = new EstAlign { Score = leftGe.Score + rightGe.Score, GStart = 0, EStart = 0, GStop = genome.Length - 1, EStop = est.Length - 1, Len = leftGe.Len + rightGe.Len };
        ge.AlignPath = new int[ge.Len];
        int i, j;
        for (i = 0, j = 0; j < leftGe.Len; i++, j++)
            ge.AlignPath[i] = leftGe.AlignPath[j];
        if (splitOnDel)
        {
            ge.AlignPath[i++] = DeleteEst;
            for (j = 1; j < rightGe.Len; i++, j++)
                ge.AlignPath[i] = rightGe.AlignPath[j];
        }
        else
        {
            for (j = 0; j < rightGe.Len; i++, j++)
                ge.AlignPath[i] = rightGe.AlignPath[j];
        }

        return ge;
    }

    private void AlignMidpoint(string est, string genome, int gapPenalty, int intronPenalty, int splicePenalty, EstSpliceSites? spliceSites, int middle, out int gleft, out int gright)
    {
        string ss = spliceSites?.Sites ?? string.Empty;
        int elen = est.Length, glen = genome.Length;
        int[] s1 = new int[elen + 1], s2 = new int[elen + 1];                 // index epos + 1
        int[] m1L = new int[elen + 1], m1R = new int[elen + 1], m2L = new int[elen + 1], m2R = new int[elen + 1];
        var bestIntronCoord = new int[elen + 1];
        var bestIntronScore = new int[elen + 1];
        var bestIntronMidL = new int[elen + 1];
        var bestIntronMidR = new int[elen + 1];
        string gseq = Lower(genome), eseq = Lower(est);
        middle++;
        for (int epos = 0; epos < elen; epos++)
        {
            s1[epos + 1] = MinusInfinity;
            bestIntronScore[epos] = MinusInfinity;
        }

        for (int gpos = 0; gpos < glen; gpos++)
        {
            (s1, s2) = (s2, s1);
            (m1L, m2L) = (m2L, m1L);
            (m1R, m2R) = (m2R, m1R);
            char g = gseq[gpos];
            bool isAcceptor = spliceSites != null && (ss[gpos] & Acceptor) != 0;
            s1[0] = MinusInfinity;
            for (int epos = 0; epos < elen; epos++)
            {
                int diagonal = unchecked(s2[epos] + Sim(g, eseq[epos]));
                int deleteEst = unchecked(s1[epos] - gapPenalty);
                int deleteGenome = unchecked(s2[epos + 1] - gapPenalty);
                int intron = isAcceptor && (ss[bestIntronCoord[epos]] & Donor) != 0
                    ? unchecked(bestIntronScore[epos] - splicePenalty)
                    : unchecked(bestIntronScore[epos] - intronPenalty);
                int max = deleteEst > deleteGenome ? deleteEst : deleteGenome;
                if (diagonal > max)
                    max = diagonal;
                if (intron > max)
                    max = intron;
                if (max == diagonal)
                {
                    s1[epos + 1] = diagonal;
                    if (epos == middle)
                    {
                        m1L[epos + 1] = gpos - 1;
                        m1R[epos + 1] = gpos;
                    }
                    else
                    {
                        m1L[epos + 1] = m2L[epos];
                        m1R[epos + 1] = m2R[epos];
                    }
                }
                else if (max == deleteEst)
                {
                    s1[epos + 1] = deleteEst;
                    if (epos == middle)
                    {
                        m1L[epos + 1] = gpos;
                        m1R[epos + 1] = gpos;
                    }
                    else
                    {
                        m1L[epos + 1] = m1L[epos];
                        m1R[epos + 1] = m1R[epos];
                    }
                }
                else if (max == deleteGenome)
                {
                    s1[epos + 1] = deleteGenome;
                    m1L[epos + 1] = m2L[epos + 1];
                    m1R[epos + 1] = m2R[epos + 1];
                }
                else
                {
                    s1[epos + 1] = intron;
                    m1L[epos + 1] = bestIntronMidL[epos];
                    m1R[epos + 1] = bestIntronMidR[epos];
                }

                if (bestIntronScore[epos] < s1[epos + 1])
                {
                    bestIntronScore[epos] = s1[epos + 1];
                    bestIntronCoord[epos] = gpos;
                    bestIntronMidL[epos] = m1L[epos + 1];
                    bestIntronMidR[epos] = m1R[epos + 1];
                }
            }
        }

        gleft = m1L[elen];
        gright = m1R[elen];
    }

    private static char IntronMark(int way)
    {
        if (way == ForwardSplicedIntron)
            return '>';
        return way == ReverseSplicedIntron ? '<' : '?';
    }

    private static string F(double value, int decimals) => CFormat.Fixed(value, decimals);

    private static string I(long v) => v.ToString(CultureInfo.InvariantCulture);

    private static string Pad(string s, int width) => s.PadLeft(width);

    private static string PadRight(string s, int width) => s.PadRight(width);

    private static float Percent(int matches, int len) => len > 0 ? (float)((double)((float)matches / (float)len) * 100.0) : 0f;

    private static void WriteMsp(StringBuilder o, ref int matches, ref int len, ref int tsub, string genomeName, int gOffset, int gsub, int gpos, string estName, string estDesc, int estLen, int eOffset, int esub, int epos, bool reverse, bool gapped)
    {
        float percent = Percent(matches, len);
        if (percent > 0)
        {
            o.Append(gapped ? "Exon     " : "Segment  ");
            if (reverse)
                o.Append(Pad(I(tsub), 5)).Append(' ').Append(Pad(F(percent, 1), 5)).Append(' ').Append(Pad(I(gsub + 1), 5)).Append(' ').Append(Pad(I(gpos), 5)).Append(' ').Append(PadRight(genomeName, 12)).Append(' ')
                    .Append(Pad(I(eOffset + estLen - esub), 5)).Append(' ').Append(Pad(I(eOffset + estLen - epos + 1), 5)).Append(' ').Append(PadRight(estName, 12)).Append("  ").Append(estDesc).Append('\n');
            else
                o.Append(Pad(I(tsub), 5)).Append(' ').Append(Pad(F(percent, 1), 5)).Append(' ').Append(Pad(I(gOffset + gsub + 1), 5)).Append(' ').Append(Pad(I(gOffset + gpos), 5)).Append(' ').Append(PadRight(genomeName, 12)).Append(' ')
                    .Append(Pad(I(eOffset + esub + 1), 5)).Append(' ').Append(Pad(I(eOffset + epos), 5)).Append(' ').Append(PadRight(estName, 12)).Append("  ").Append(estDesc).Append('\n');
        }

        matches = len = tsub = 0;
    }

    /// <summary><c>embEstOutBlastStyle</c>: the Exon / Intron / Span lines (<paramref name="gapped"/>) or the Segment lines of the alignment.</summary>
    public void OutBlastStyle(StringBuilder o, string genome, string genomeName, int gOffset, string est, string estName, string estDesc, int eOffset, EstAlign ge, int gapPenalty, int intronPenalty, int splicePenalty, bool gapped, bool reverse)
    {
        int gsub = ge.GStart, gpos = ge.GStart, esub = ge.EStart, epos = ge.EStart;
        int tsub = 0, matches = 0, len = 0, totalMatches = 0, totalLen = 0;
        int elen = est.Length;
        int p = -1;
        while (++p < ge.Len)
        {
            if (ge.AlignPath[p] <= Intron)
            {
                WriteMsp(o, ref matches, ref len, ref tsub, genomeName, gOffset, gsub, gpos, estName, estDesc, elen, eOffset, esub, epos, reverse, gapped);
                if (gapped)
                {
                    string kind = "-Intron  ";
                    if (ge.AlignPath[p] == Intron)
                        kind = "?Intron  ";
                    else if (ge.AlignPath[p] == ForwardSplicedIntron)
                        kind = "+Intron  ";
                    int score = ge.AlignPath[p] == Intron ? -intronPenalty : -splicePenalty;
                    o.Append(kind).Append(Pad(I(score), 5)).Append(' ').Append(Pad(F(0.0, 1), 5)).Append(' ').Append(Pad(I(gOffset + gpos + 1), 5)).Append(' ').Append(Pad(I(gOffset + gpos + ge.AlignPath[p + 1]), 5)).Append(' ').Append(PadRight(genomeName, 12)).Append('\n');
                }

                gpos += ge.AlignPath[++p];
                esub = epos;
                gsub = gpos;
            }
            else if (ge.AlignPath[p] == Diagonal)
            {
                int m = Sim(genome[gpos], est[epos]);
                tsub += m;
                if (m > 0)
                {
                    matches++;
                    totalMatches++;
                }

                len++;
                totalLen++;
                gpos++;
                epos++;
            }
            else if (ge.AlignPath[p] == DeleteEst)
            {
                if (gapped)
                {
                    tsub -= gapPenalty;
                    epos++;
                    len++;
                    totalLen++;
                }
                else
                {
                    WriteMsp(o, ref matches, ref len, ref tsub, genomeName, gOffset, gsub, gpos, estName, estDesc, elen, eOffset, esub, epos, reverse, gapped);
                    epos++;
                    esub = epos;
                    gsub = gpos;
                }
            }
            else if (ge.AlignPath[p] == DeleteGenome)
            {
                if (gapped)
                {
                    tsub -= gapPenalty;
                    gpos++;
                    totalLen++;
                    len++;
                }
                else
                {
                    WriteMsp(o, ref matches, ref len, ref tsub, genomeName, gOffset, gsub, gpos, estName, estDesc, elen, eOffset, esub, epos, reverse, gapped);
                    gpos++;
                    esub = epos;
                    gsub = gpos;
                }
            }
        }

        WriteMsp(o, ref matches, ref len, ref tsub, genomeName, gOffset, gsub, gpos, estName, estDesc, elen, eOffset, esub, epos, reverse, gapped);
        if (gapped)
        {
            float percent = totalLen > 0 ? (float)((double)((float)totalMatches / (float)totalLen) * 100.0) : 0f;
            if (reverse)
                o.Append("\nSpan     ").Append(Pad(I(ge.Score), 5)).Append(' ').Append(Pad(F(percent, 1), 5)).Append(' ').Append(Pad(I(gOffset + ge.GStart + 1), 5)).Append(' ').Append(Pad(I(gOffset + ge.GStop + 1), 5)).Append(' ')
                    .Append(PadRight(genomeName, 12)).Append(' ').Append(Pad(I(eOffset + elen - ge.EStart), 5)).Append(' ').Append(Pad(I(eOffset + elen - ge.EStop), 5)).Append(' ').Append(PadRight(estName, 12)).Append("  ").Append(estDesc).Append('\n');
            else
                o.Append("\nSpan     ").Append(Pad(I(ge.Score), 5)).Append(' ').Append(Pad(F(percent, 1), 5)).Append(' ').Append(Pad(I(gOffset + ge.GStart + 1), 5)).Append(' ').Append(Pad(I(gOffset + ge.GStop + 1), 5)).Append(' ')
                    .Append(PadRight(genomeName, 12)).Append(' ').Append(Pad(I(eOffset + ge.EStart + 1), 5)).Append(' ').Append(Pad(I(eOffset + ge.EStop + 1), 5)).Append(' ').Append(PadRight(estName, 12)).Append("  ").Append(estDesc).Append('\n');
        }
    }

    /// <summary><c>embEstPrintAlign</c>: the alignment in blocks of <paramref name="width"/> columns (introns shown as their first and last bases and their length).</summary>
    public void PrintAlign(StringBuilder o, string genome, string genomeName, string est, string estName, EstAlign ge, int width)
    {
        int namelen = Math.Max(genomeName.Length, estName.Length);
        o.Append('\n');
        int cap = genome.Length + est.Length + 1;
        var gbuf = new char[cap];
        var ebuf = new char[cap];
        var sbuf = new char[cap];
        var gcoord = new int[cap];
        var ecoord = new int[cap];
        int gpos = ge.GStart, epos = ge.EStart, len = 0;
        int pos = -1;
        while (++pos < ge.Len)
        {
            int way = ge.AlignPath[pos];
            if (way == Diagonal)
            {
                gcoord[len] = gpos;
                ecoord[len] = epos;
                gbuf[len] = genome[gpos++];
                ebuf[len] = est[epos++];
                sbuf[len] = Sim(gbuf[len], ebuf[len]) > 0 ? '|' : ' ';
                len++;
            }
            else if (way == DeleteEst)
            {
                gcoord[len] = gpos;
                ecoord[len] = epos;
                gbuf[len] = '-';
                ebuf[len] = est[epos++];
                sbuf[len] = ' ';
                len++;
            }
            else if (way == DeleteGenome)
            {
                gcoord[len] = gpos;
                ecoord[len] = epos;
                gbuf[len] = genome[gpos++];
                ebuf[len] = '-';
                sbuf[len] = ' ';
                len++;
            }
            else if (way <= Intron)
            {
                int intronWidth = ge.AlignPath[pos + 1];
                int g = gpos - 1;
                int halfWidth = intronWidth > 10 ? 5 : intronWidth / 2;
                string number = " " + I(intronWidth) + " ";
                int numlen = number.Length;
                int j;
                for (j = len; j < len + halfWidth; j++)
                {
                    g++;
                    gcoord[j] = gpos - 1;
                    ecoord[j] = epos - 1;
                    gbuf[j] = char.ToLowerInvariant(genome[g]);
                    ebuf[j] = '.';
                    sbuf[j] = IntronMark(way);
                }

                len = j;
                for (j = len; j < len + numlen; j++)
                {
                    gcoord[j] = gpos - 1;
                    ecoord[j] = epos - 1;
                    gbuf[j] = '.';
                    ebuf[j] = '.';
                    sbuf[j] = number[j - len];
                }

                len = j;
                g = gpos + intronWidth - halfWidth - 1;
                for (j = len; j < len + halfWidth; j++)
                {
                    g++;
                    gcoord[j] = gpos - 1;
                    ecoord[j] = epos - 1;
                    gbuf[j] = char.ToLowerInvariant(genome[g]);
                    ebuf[j] = '.';
                    sbuf[j] = IntronMark(way);
                }

                gpos += ge.AlignPath[++pos];
                len = j;
            }
        }

        for (int i = 0; i < len; i += width)
        {
            int max = i + width > len ? len : i + width;
            o.Append(Pad(genomeName, namelen)).Append(' ').Append(Pad(I(gcoord[i] + 1), 6)).Append(' ');
            for (int j = i; j < max; j++)
                o.Append(gbuf[j]);
            o.Append(' ').Append(Pad(I(gcoord[max - 1] + 1), 6)).Append('\n');
            o.Append(' ', namelen + 8);
            for (int j = i; j < max; j++)
                o.Append(sbuf[j]);
            o.Append('\n');
            o.Append(Pad(estName, namelen)).Append(' ').Append(Pad(I(ecoord[i] + 1), 6)).Append(' ');
            for (int j = i; j < max; j++)
                o.Append(ebuf[j]);
            o.Append(' ').Append(Pad(I(ecoord[max - 1] + 1), 6)).Append("\n\n");
        }

        o.Append("\nAlignment Score: ").Append(I(ge.Score)).Append('\n');
    }
}
