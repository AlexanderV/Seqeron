// SPDX-License-Identifier: GPL-3.0-or-later
// Derived from LAST 1256 (https://gitlab.com/mcfrith/last: GappedXdropAligner.cc / .hh / Inl.hh, gaplessXdrop.cc, mcf_gap_costs.cc), Copyright 2008 - 2019 Martin C. Frith; GPL-3.0-or-later.
// Part of Seqeron.Genomics.Last, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Last;

/// <summary>The gap costs of <c>mcf::GapCosts</c> (one affine piece for deletions and one for insertions).</summary>
internal sealed class LastGapCosts
{
    public int DelOpen, DelGrow, InsOpen, InsGrow, PairCost;

    public bool IsAffine;

    public LastGapCosts(int delOpen, int delGrow, int insOpen, int insGrow, int unalignedPairCost)
    {
        if (delGrow <= 0 || insGrow <= 0)
            throw new LastProgramException("gap extension cost must be > 0");
        DelOpen = delOpen;
        DelGrow = delGrow;
        InsOpen = insOpen;
        InsGrow = insGrow;
        PairCost = unalignedPairCost > 0 ? unalignedPairCost : delOpen + delGrow + insOpen + insGrow;
        IsAffine = PairCost >= delGrow + insGrow + Math.Max(delOpen, insOpen);
    }

    /// <summary><c>GapCosts::cost</c>.</summary>
    public int Cost(long refInsertLen, long qryInsertLen)
    {
        long genCost = int.MaxValue;
        long delCost = 0;
        if (refInsertLen > 0)
        {
            delCost = int.MaxValue;
            long s = DelOpen + DelGrow * refInsertLen;
            delCost = Math.Min(delCost, s);
            if (refInsertLen >= qryInsertLen)
                genCost = Math.Min(genCost, s + (PairCost - DelGrow) * qryInsertLen);
        }

        long insCost = 0;
        if (qryInsertLen > 0)
        {
            insCost = int.MaxValue;
            long s = InsOpen + InsGrow * qryInsertLen;
            insCost = Math.Min(insCost, s);
            if (qryInsertLen >= refInsertLen)
                genCost = Math.Min(genCost, s + (PairCost - InsGrow) * refInsertLen);
        }

        return (int)Math.Min(delCost + insCost, genCost);
    }
}

/// <summary>The gapless x-drop extensions of <c>gaplessXdrop.cc</c>. Arrays are the letter-code arrays of the sequences; positions are indexes.</summary>
internal static class LastGapless
{
    public static int ForwardScore(byte[] seq1, int p1, byte[] seq2, int p2, int[][] scorer, int maxScoreDrop)
    {
        int score = 0, s = 0;
        while (true)
        {
            s += scorer[seq1[p1++]][seq2[p2++]];
            if (s < score - maxScoreDrop)
                break;
            if (s > score)
                score = s;
        }

        if (score - s < 0)
            throw new LastProgramException("score overflow in forward gapless extension");
        return score;
    }

    public static int ReverseScore(byte[] seq1, int p1, byte[] seq2, int p2, int[][] scorer, int maxScoreDrop)
    {
        int score = 0, s = 0;
        while (true)
        {
            s += scorer[seq1[--p1]][seq2[--p2]];
            if (s < score - maxScoreDrop)
                break;
            if (s > score)
                score = s;
        }

        if (score - s < 0)
            throw new LastProgramException("score overflow in reverse gapless extension");
        return score;
    }

    public static int ForwardEnd(byte[] seq1, int p1, byte[] seq2, int p2, int[][] scorer, int score)
    {
        int s = 0;
        while (s < score)
            s += scorer[seq1[p1++]][seq2[p2++]];
        return p1;
    }

    public static int ReverseEnd(byte[] seq1, int p1, byte[] seq2, int p2, int[][] scorer, int score)
    {
        int s = 0;
        while (s < score)
            s += scorer[seq1[--p1]][seq2[--p2]];
        return p1;
    }

    public static bool IsOptimal(byte[] seq1, int begin1, int end1, byte[] seq2, int begin2, int[][] scorer, int maxScoreDrop)
    {
        int score = 0, maxScore = 0;
        int p1 = begin1, p2 = begin2;
        while (p1 < end1)
        {
            score += scorer[seq1[p1++]][seq2[p2++]];
            if (score > maxScore)
                maxScore = score;
            else if (score <= 0 || p1 == end1 || score < maxScore - maxScoreDrop)
                return false;
        }

        return true;
    }

    public static int AlignmentScore(byte[] seq1, int begin1, int end1, byte[] seq2, int begin2, int[][] scorer)
    {
        int score = 0;
        for (int p1 = begin1, p2 = begin2; p1 < end1; p1++, p2++)
            score += scorer[seq1[p1]][seq2[p2]];
        return score;
    }
}

/// <summary>
/// <c>GappedXdropAligner::align</c> and <c>getNextChunk</c> for one direction of an extension (local alignment, affine or generalised-affine gap costs). The scores of the three
/// states live in antidiagonals exactly as in LAST (padding, origins), so the traceback takes the same decisions on ties.
/// </summary>
internal sealed class LastGappedAligner
{
    private const int Inf = LastScoreMatrix.Inf;
    private const int XdropPadLen = 16;
    private const int SimdLen = 4;

    private int[] _x = new int[1024], _y = new int[1024], _z = new int[1024];
    private readonly List<long> _scoreEndsAndOrigins = [];
    private int _bestAntidiagonal;
    private int _bestSeq1Position;

    private void Resize(long size)
    {
        if (_x.Length < size)
        {
            long n = Math.Max(size, _x.Length * 2L);
            Array.Resize(ref _x, (int)n);
            Array.Resize(ref _y, (int)n);
            Array.Resize(ref _z, (int)n);
        }
    }

    private void InitAntidiagonal(int antidiagonalIncludingDummies, int seq1end, int thisEnd, int numCells)
    {
        int nextEnd = thisEnd + XdropPadLen + numCells;
        int a = 2 * (antidiagonalIncludingDummies + 1);
        while (_scoreEndsAndOrigins.Count <= a)
            _scoreEndsAndOrigins.Add(0);
        _scoreEndsAndOrigins[a - 1] = nextEnd - seq1end;
        _scoreEndsAndOrigins[a] = nextEnd;
        Resize(nextEnd + SimdLen - 1);
        for (int i = 0; i < XdropPadLen; i++)
        {
            _x[thisEnd + i] = -Inf;
            _y[thisEnd + i] = -Inf;
            _z[thisEnd + i] = -Inf;
        }
    }

    private long ScoreOrigin(int antidiagonalIncludingDummies) => _scoreEndsAndOrigins[2 * antidiagonalIncludingDummies + 1];

    private long Seq1Start(int antidiagonal)
    {
        int a = 2 * (antidiagonal + 2);
        return _scoreEndsAndOrigins[a] + XdropPadLen - _scoreEndsAndOrigins[a + 1];
    }

    private int Hori(int antidiagonal, int seq1coordinate) => (int)(ScoreOrigin(antidiagonal + 1) + seq1coordinate - 1);

    private int Vert(int antidiagonal, int seq1coordinate) => (int)(ScoreOrigin(antidiagonal + 1) + seq1coordinate);

    private int Diag(int antidiagonal, int seq1coordinate) => (int)(ScoreOrigin(antidiagonal) + seq1coordinate - 1);

    private static byte At(byte[] seq, int index) => index >= 0 && index < seq.Length ? seq[index] : (byte)' ';

    private static bool IsDelimiter(int c, int[] scores) => scores[c] <= -Inf;

    /// <summary>Extends from (start1, start2) in the given direction; returns the best score (<c>align</c>). Only local alignments (globality 0) are ported.</summary>
    public int Align(byte[] seq1, int start1, byte[] seq2, int start2, bool isForward, int[][] scorer, LastGapCosts gap, int maxScoreDrop, int maxMatchScore)
    {
        int seqIncrement = isForward ? 1 : -1;
        int numCells = 1;
        int seq1end = 1;
        int diagPos = XdropPadLen - 1;
        int horiPos = XdropPadLen * 2 - 1;
        int thisPos = XdropPadLen * 2;
        int bestScore = 0;
        _scoreEndsAndOrigins.Clear();
        InitAntidiagonal(0, 0, 0, 0);
        InitAntidiagonal(1, 0, XdropPadLen, 0);
        _x[XdropPadLen - 1] = 0;
        _bestAntidiagonal = 0;
        int p1 = start1, p2 = start2;
        var pssmQueue = new List<int[]>();
        var seq2Queue = new List<byte>();
        bool isDelimiter1 = IsDelimiter(0, scorer[At(seq1, p1)]);
        bool isDelimiter2 = IsDelimiter(At(seq2, p2), scorer[0]);
        for (int i = 0; i < SimdLen; i++)
        {
            int[] x = scorer[At(seq1, p1)];
            pssmQueue.Add(x);
            p1 += seqIncrement * (IsDelimiter(0, x) ? 0 : 1);
            seq2Queue.Add(At(seq2, p2));
        }

        p2 += seqIncrement;
        int antidiagonal = 0;
        while (true)
        {
            int n = numCells - 1;
            int s1Base = pssmQueue.Count - (n + SimdLen);
            InitAntidiagonal(antidiagonal + 2, seq1end, thisPos, numCells);
            thisPos += XdropPadLen;
            int x0 = thisPos;
            int y1 = horiPos;
            int z1 = horiPos + 1;
            int x2 = diagPos;
            if (isDelimiter1 || isDelimiter2)
                maxScoreDrop = Math.Min(maxScoreDrop, n * maxMatchScore - 1);
            int minScore = bestScore - maxScoreDrop;
            if (gap.IsAffine)
            {
                for (int i = 0; i < numCells; i++)
                {
                    int s = pssmQueue[s1Base + i][seq2Queue[seq2Queue.Count - 1 - i]];
                    int x = _x[x2 + i];
                    int y = _y[y1 + i] - gap.DelGrow;
                    int z = _z[z1 + i] - gap.InsGrow;
                    int b = Math.Max(Math.Max(x, y), z);
                    if (b > bestScore)
                    {
                        bestScore = b;
                        _bestAntidiagonal = antidiagonal;
                    }

                    _x[x0 + i] = minScore > b ? -Inf : b + s;
                    _y[x0 + i] = Math.Max(b - gap.DelOpen, y);
                    _z[x0 + i] = Math.Max(b - gap.InsOpen, z);
                }
            }
            else
            {
                int y2 = diagPos, z2 = diagPos;
                for (int i = 0; i < numCells; i++)
                {
                    int x = _x[x2 + i];
                    int y = Math.Max(_y[y1 + i] - gap.DelGrow, _y[y2 + i] - gap.PairCost);
                    int z = Math.Max(_z[z1 + i] - gap.InsGrow, _z[z2 + i] - gap.PairCost);
                    int b = Math.Max(Math.Max(x, y), z);
                    if (b >= minScore)
                    {
                        if (b > bestScore)
                        {
                            bestScore = b;
                            _bestAntidiagonal = antidiagonal;
                        }

                        _x[x0 + i] = b + pssmQueue[s1Base + i][seq2Queue[seq2Queue.Count - 1 - i]];
                        _y[x0 + i] = Math.Max(b - gap.DelOpen, y);
                        _z[x0 + i] = Math.Max(b - gap.InsOpen, z);
                    }
                    else
                    {
                        _x[x0 + i] = _y[x0 + i] = _z[x0 + i] = -Inf;
                    }
                }
            }

            diagPos = horiPos;
            horiPos = thisPos - 1;
            thisPos += numCells;
            if (_x[x0 + n] > -Inf / 2)
            {
                numCells++;
                seq1end++;
                int[] x = scorer[At(seq1, p1)];
                pssmQueue.Add(x);
                p1 += seqIncrement * (IsDelimiter(0, x) ? 0 : 1);
                isDelimiter1 = IsDelimiter(0, pssmQueue[pssmQueue.Count - SimdLen]);
            }

            if (_x[x0] > -Inf / 2)
            {
                byte y = At(seq2, p2);
                seq2Queue.Add(y);
                p2 += seqIncrement;
                isDelimiter2 = IsDelimiter(y, scorer[0]);
            }
            else
            {
                numCells--;
                if (numCells == 0)
                    break;
                diagPos++;
                horiPos++;
            }

            antidiagonal++;
        }

        // calcBestSeq1position(bestScore, 2)
        long seq1beg = Seq1Start(_bestAntidiagonal + 2 - 2);
        int xi = Diag(_bestAntidiagonal, (int)seq1beg);
        int xBeg = xi;
        while (_x[xi] != bestScore)
            xi++;
        _bestSeq1Position = xi - xBeg + (int)seq1beg;
        return bestScore;
    }

    /// <summary><c>getNextChunk</c>: the next gapless block of the traceback (end1, end2, length), false at the end.</summary>
    public bool GetNextChunk(out int end1, out int end2, out int length, LastGapCosts gap)
    {
        end1 = end2 = length = 0;
        if (_bestAntidiagonal == 0)
            return false;
        end1 = _bestSeq1Position;
        end2 = _bestAntidiagonal - _bestSeq1Position;
        const int undefinedLength = -1;
        length = undefinedLength;
        int state = 0;
        while (true)
        {
            int h = Hori(_bestAntidiagonal, _bestSeq1Position);
            int v = Vert(_bestAntidiagonal, _bestSeq1Position);
            int d = Diag(_bestAntidiagonal, _bestSeq1Position);
            int x = _x[d];
            int y = _y[h] - gap.DelGrow;
            int z = _z[v] - gap.InsGrow;
            int a = _y[d] - gap.PairCost;
            int b = _z[d] - gap.PairCost;
            if (state == 1 || state == 3)
            {
                y += gap.DelOpen;
                a += gap.DelOpen;
            }

            if (state == 2 || state == 4)
            {
                z += gap.InsOpen;
                b += gap.InsOpen;
            }

            state = MaxIndex(x, y, z, a, b);
            if (length == undefinedLength && (state > 0 || _bestAntidiagonal == 0))
                length = end1 - _bestSeq1Position;
            if (length != undefinedLength && state == 0)
                return true;
            if (state < 1 || state > 2)
                _bestAntidiagonal -= 2;
            else
                _bestAntidiagonal -= 1;
            if (state != 2)
                _bestSeq1Position -= 1;
        }
    }

    /// <summary><c>maxIndex(a, b, c, d, e)</c>: the first index of the maximum, later entries must be strictly greater (the recursion of <c>GappedXdropAlignerInl.hh</c>).</summary>
    private static int MaxIndex(int a, int b, int c, int d, int e)
    {
        int[] v = [a, b, c, d, e];
        return MaxIndexRec(v, 0, 5);
    }

    private static int MaxIndexRec(int[] v, int start, int count)
    {
        // maxIndex of v[start .. start+count): last > first ? maxIndex(rest without first) + 1 : maxIndex(without last)
        if (count == 1)
            return 0;
        if (count == 2)
            return v[start + 1] > v[start] ? 1 : 0;
        return v[start + count - 1] > v[start] ? MaxIndexRec(v, start + 1, count - 1) + 1 : MaxIndexRec(v, start, count - 1);
    }
}
