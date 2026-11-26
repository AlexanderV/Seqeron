// SPDX-License-Identifier: GPL-3.0-or-later
// Derived from LAST 1256 (https://gitlab.com/mcfrith/last: SegmentPair.cc, SegmentPairPot.cc, DiagonalTable.cc, Alignment.cc, AlignmentPot.cc, AlignmentWrite.cc), Copyright 2008 - 2019 Martin C. Frith; GPL-3.0-or-later.
// Part of Seqeron.Genomics.Last, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.Last;

/// <summary>A pair of equal-length segments (a gapless alignment); a class so that the pots can mark the entries in place as LAST does. Coordinates are unsigned 32-bit as in LAST (differences wrap around).</summary>
internal sealed class LastSegmentPair
{
    public uint Start1, Start2, Size;

    public int Score;

    public LastSegmentPair(uint s1, uint s2, uint size, int score = 0)
    {
        Start1 = s1;
        Start2 = s2;
        Size = size;
        Score = score;
    }

    public uint Beg1 => Start1;

    public uint Beg2 => Start2;

    public uint End1 => Start1 + Size;

    public uint End2 => Start2 + Size;

    /// <summary>The coordinate in one sequence minus the coordinate in the other (may wrap around).</summary>
    public uint Diagonal => Start1 - Start2;

    public bool SameAs(LastSegmentPair s) => Start1 == s.Start1 && Start2 == s.Start2 && Size == s.Size && Score == s.Score;

    public LastSegmentPair Clone() => new(Start1, Start2, Size, Score);

    /// <summary><c>maxIdenticalRun</c>: shrink to the longest run of identical letters (case-insensitive).</summary>
    public void MaxIdenticalRun(byte[] seq1, byte[] seq2, byte[] map1, byte[] map2)
    {
        int s1 = (int)Beg1, s2 = (int)Beg2, e1 = (int)End1;
        int bestEnd1 = s1;
        uint bestSize = 0, runSize = 0;
        while (s1 < e1)
        {
            if (map1[seq1[s1++]] == map2[seq2[s2++]])
            {
                runSize++;
                if (runSize > bestSize)
                {
                    bestSize = runSize;
                    bestEnd1 = s1;
                }
            }
            else
            {
                runSize = 0;
            }
        }

        int bestBeg1 = bestEnd1 - (int)bestSize;
        uint offset = (uint)(bestBeg1 - (int)Beg1);
        Start1 += offset;
        Start2 += offset;
        Size = bestSize;
    }
}

/// <summary><c>DiagonalTable</c>: the furthest sequential position covered so far on each diagonal (256 bins).</summary>
internal sealed class LastDiagonalTable
{
    private const int Bins = 256;
    private readonly List<(uint First, uint Second)>[] _hits = Enumerable.Range(0, Bins).Select(_ => new List<(uint, uint)>()).ToArray();

    public bool IsCovered(uint sequentialPos, uint randomPos)
    {
        uint diagonal = sequentialPos - randomPos;
        var v = _hits[diagonal % Bins];
        for (int i = 0; i < v.Count;)
        {
            if (v[i].First >= sequentialPos)
            {
                if (v[i].Second == diagonal)
                    return true;
                i++;
            }
            else
            {
                v.RemoveAt(i);
            }
        }

        return false;
    }

    public void AddEndpoint(uint sequentialPos, uint randomPos)
    {
        uint diagonal = sequentialPos - randomPos;
        _hits[diagonal % Bins].Add((sequentialPos, diagonal));
    }
}

/// <summary><c>SegmentPairPot</c>: gapless alignments sorted by position and by score; overlaps are marked (score 0).</summary>
internal sealed class LastSegmentPairPot
{
    public List<LastSegmentPair> Items { get; set; } = [];

    private List<LastSegmentPair> _iters = [];

    public int Size => Items.Count;

    public void Add(LastSegmentPair sp) => Items.Add(sp);

    public static void Mark(LastSegmentPair s) => s.Score = 0;

    public static bool IsMarked(LastSegmentPair s) => s.Score == 0;

    public static void MarkAsGood(LastSegmentPair s) => s.Size = 0;

    public static bool IsNotMarkedAsGood(LastSegmentPair s) => s.Size != 0;

    public static int ItemCompare(LastSegmentPair x, LastSegmentPair y)
    {
        if (x.Diagonal != y.Diagonal)
            return x.Diagonal.CompareTo(y.Diagonal);
        return x.Beg1.CompareTo(y.Beg1);
    }

    private static int IterCompare(LastSegmentPair x, LastSegmentPair y)
    {
        if (x.Score != y.Score)
            return y.Score.CompareTo(x.Score);
        if (x.Start1 != y.Start1)
            return x.Start1.CompareTo(y.Start1);
        return x.Start2.CompareTo(y.Start2);
    }

    /// <summary>Sorts by position, removes duplicates, and orders the score index (<c>sort</c>).</summary>
    public void Sort()
    {
        var sorted = Items.OrderBy(x => x, Comparer<LastSegmentPair>.Create(ItemCompare)).ToList();
        var unique = new List<LastSegmentPair>();
        foreach (var s in sorted)
            if (unique.Count == 0 || !unique[^1].SameAs(s))
                unique.Add(s);
        Items = unique;
        _iters = Items.OrderBy(x => x, Comparer<LastSegmentPair>.Create(IterCompare)).ToList();
    }

    public LastSegmentPair Get(int i) => _iters[i];

    private int LowerBound(LastSegmentPair sp)
    {
        int lo = 0, hi = Items.Count;
        while (lo < hi)
        {
            int mid = lo + (hi - lo) / 2;
            if (ItemCompare(Items[mid], sp) < 0)
                lo = mid + 1;
            else
                hi = mid;
        }

        return lo;
    }

    public void MarkOverlaps(LastSegmentPair sp)
    {
        int i = LowerBound(sp);
        if (i > 0 && Items[i - 1].Diagonal == sp.Diagonal && Items[i - 1].End1 > sp.Beg1)
            Mark(Items[i - 1]);
        while (i < Items.Count && Items[i].Diagonal == sp.Diagonal && Items[i].Beg1 < sp.End1)
        {
            Mark(Items[i]);
            i++;
        }
    }

    public void MarkAllOverlaps(IEnumerable<LastSegmentPair> sps)
    {
        foreach (var sp in sps)
            MarkOverlaps(sp);
    }

    public void MarkTandemRepeats(LastSegmentPair sp, uint maxDistance)
    {
        if (Items.Count == 0)
            throw new InvalidOperationException("no items");
        var nextDiagonal = new LastSegmentPair(sp.Beg1 + 1, sp.Beg2, 0, 0);
        int n = LowerBound(nextDiagonal);
        if (n == Items.Count)
            n = 0;
        int j = n;
        do
        {
            if (Items[j].Diagonal - sp.Diagonal > maxDistance)
                break;
            if (Items[j].Beg2 >= sp.Beg2 && Items[j].End1 <= sp.End1)
                Mark(Items[j]);
            j++;
            if (j == Items.Count)
                j = 0;
        }
        while (j != n);
        var prevDiagonal = new LastSegmentPair(sp.End1 - 1, sp.End2, 0, 0);
        int p = LowerBound(prevDiagonal);
        if (p == Items.Count)
            p = 0;
        int k = p;
        do
        {
            if (k == 0)
                k = Items.Count;
            k--;
            if (sp.Diagonal - Items[k].Diagonal > maxDistance)
                break;
            if (Items[k].Beg1 >= sp.Beg1 && Items[k].End2 <= sp.End2)
                Mark(Items[k]);
        }
        while (k != p);
    }
}

/// <summary>A gapped alignment: gapless blocks, the score and the seed it grew from.</summary>
internal sealed class LastAlignment
{
    public List<LastSegmentPair> Blocks { get; private set; } = [];

    public double Score;

    public LastSegmentPair Seed = new(0, 0, 0);

    public bool Marked;

    public uint Beg1 => Blocks[0].Beg1;

    public uint Beg2 => Blocks[0].Beg2;

    public uint End1 => Blocks[^1].End1;

    public uint End2 => Blocks[^1].End2;

    private static bool IsNext(LastSegmentPair x, LastSegmentPair y) => x.End1 == y.Beg1 && x.End2 == y.Beg2;

    /// <summary><c>Alignment::makeXdrop</c> for the ordinary (untranslated, non-greedy, local) case: extend the seed to the left and to the right.</summary>
    public void MakeXdrop(LastGappedAligner aligner, byte[] seq1, byte[] seq2, int[][] scoreMatrix, int smMax, LastGapCosts gap, int maxDrop)
    {
        // score = seed.score (the probability ratios are always given in the gapped phases of lastal)
        Score = Seed.Score;
        Blocks = [];
        Extend(Blocks, aligner, seq1, seq2, (int)Seed.Beg1, (int)Seed.Beg2, false, scoreMatrix, smMax, gap, maxDrop);
        uint seedBeg1 = Seed.Beg1;
        uint seedBeg2 = Seed.Beg2;
        foreach (var b in Blocks)
        {
            b.Start1 = seedBeg1 - b.Start1 - b.Size;
            b.Start2 = seedBeg2 - b.Start2 - b.Size;
        }

        var forwardBlocks = new List<LastSegmentPair>();
        Extend(forwardBlocks, aligner, seq1, seq2, (int)Seed.End1, (int)Seed.End2, true, scoreMatrix, smMax, gap, maxDrop);
        uint seedEnd1 = Seed.End1;
        uint seedEnd2 = Seed.End2;
        foreach (var b in forwardBlocks)
        {
            b.Start1 = seedEnd1 + b.Start1;
            b.Start2 = seedEnd2 + b.Start2;
        }

        bool isMergeSeedReverse = Blocks.Count > 0 && IsNext(Blocks[^1], Seed);
        bool isMergeSeedForward = forwardBlocks.Count > 0 && IsNext(Seed, forwardBlocks[^1]);
        if (Seed.Size == 0 && !isMergeSeedReverse && !isMergeSeedForward)
        {
            Score = -LastScoreMatrix.Inf;
            return;
        }

        if (isMergeSeedReverse)
            Blocks[^1].Size += Seed.Size;
        else
            Blocks.Add(Seed.Clone());
        if (isMergeSeedForward)
        {
            Blocks[^1].Size += forwardBlocks[^1].Size;
            forwardBlocks.RemoveAt(forwardBlocks.Count - 1);
        }

        for (int i = forwardBlocks.Count - 1; i >= 0; i--)
            Blocks.Add(forwardBlocks[i]);
    }

    private void Extend(List<LastSegmentPair> chunks, LastGappedAligner aligner, byte[] seq1, byte[] seq2, int start1, int start2, bool isForward, int[][] sm, int smMax, LastGapCosts gap, int maxDrop)
    {
        if (!isForward)
        {
            start1--;
            start2--;
        }

        int extensionScore = aligner.Align(seq1, start1, seq2, start2, isForward, sm, gap, maxDrop, smMax);
        while (aligner.GetNextChunk(out int end1, out int end2, out int size, gap))
            chunks.Add(new LastSegmentPair((uint)(end1 - size), (uint)(end2 - size), (uint)size));
        Score += extensionScore;
    }

    private static int GapCost(LastSegmentPair x, LastSegmentPair y, LastGapCosts gap) => gap.Cost((int)(y.Beg1 - x.End1), (int)(y.Beg2 - x.End2));

    /// <summary><c>Alignment::isOptimal</c> for local alignments.</summary>
    public bool IsOptimal(byte[] seq1, byte[] seq2, int[][] scoreMatrix, int maxDrop, LastGapCosts gap)
    {
        int maxScore = 0;
        int runningScore = 0;
        for (int i = 0; i < Blocks.Count; i++)
        {
            var y = Blocks[i];
            if (i > 0)
            {
                var x = Blocks[i - 1];
                runningScore -= GapCost(x, y, gap);
                if (runningScore <= 0)
                    return false;
                if (runningScore < maxScore - maxDrop)
                    return false;
            }

            int s1 = (int)y.Beg1, s2 = (int)y.Beg2, e1 = (int)y.End1;
            while (s1 < e1)
            {
                runningScore += scoreMatrix[seq1[s1++]][seq2[s2++]];
                if (runningScore > maxScore)
                    maxScore = runningScore;
                else if (runningScore <= 0)
                    return false;
                else if (s1 == e1 && i + 1 == Blocks.Count)
                    return false;
                else if (runningScore < maxScore - maxDrop)
                    return false;
            }
        }

        return true;
    }

    /// <summary><c>numColumns</c> (untranslated).</summary>
    public long NumColumns()
    {
        long num = 0;
        for (int i = 0; i < Blocks.Count; i++)
        {
            var y = Blocks[i];
            if (i > 0)
            {
                var x = Blocks[i - 1];
                num += y.Beg1 - x.End1;
                uint b = x.End2, e = y.Beg2;
                if (b < e)
                    num += e - b;
            }

            num += y.Size;
        }

        return num;
    }

    private static string SignedDifference(uint x, uint y) => x >= y ? (x - y).ToString(CultureInfo.InvariantCulture) : "-" + (y - x).ToString(CultureInfo.InvariantCulture);

    /// <summary>The blocks field of the TAB format: size,gap1:gap2,size ...</summary>
    public string BlockText()
    {
        var sb = new StringBuilder();
        for (int i = 0; i < Blocks.Count; i++)
        {
            var y = Blocks[i];
            if (i > 0)
            {
                var x = Blocks[i - 1];
                if (x.Size != 0)
                    sb.Append(',');
                sb.Append(SignedDifference(y.Beg1, x.End1)).Append(':').Append(SignedDifference(y.Beg2, x.End2));
                if (y.Size != 0)
                    sb.Append(',');
            }

            if (y.Size != 0)
                sb.Append(y.Size.ToString(CultureInfo.InvariantCulture));
        }

        return sb.ToString();
    }
}
