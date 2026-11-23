using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public sealed partial class PairwiseAligner
{
    private const int HorizontalMove = 1, VerticalMove = 2, DiagonalMove = 4, StartPoint = 8;

    private sealed class FogsaaCell
    {
        public double PresentScore, Lower, Upper;
        public int Type;
        public bool Filled, IsLeftGap;
    }

    private struct FogsaaNode
    {
        public int PA, PB, TypeUpToNext, NextType;
        public double NextLower, NextUpper;
    }

    /// <summary>The warnings the last FOGSAA run issued (Biopython's <c>BiopythonWarning</c>): match below mismatch, a gap score above the mismatch score ("Algorithm may return incorrect results").</summary>
    public IReadOnlyList<string> Warnings { get; private set; } = [];

    /// <summary>
    /// The Fast Optimal Global Sequence Alignment Algorithm (FOGSAA; Chakraborty &amp; Bandyopadhyay 2013, PLoS ONE 8:e85751) as <c>_pairwisealigner.c</c> of Biopython 1.88 implements it (<c>mode = "fogsaa"</c>): a branch-and-bound search over
    /// the alignment paths with upper and lower score bounds in a priority queue; one optimal global alignment. A line-by-line port, including the choices of the C code that make it differ from the dynamic programs on ties
    /// (second-child scores use the scalar match / mismatch, the search order of equal bounds is the heap order of the C implementation).
    /// </summary>
    private (double Score, List<(char A, char B)> Columns) RunFogsaa(string a, string b, char strand, Func<char, char, double> alignScore, double match, double mismatch)
    {
        int nA = a.Length, nB = b.Length;
        double t, t2;
        int curpA = 0, curpB = 0;
        int pathend;
        var childTypes = new int[3];
        double lowerBound;
        var childLbounds = new double[3];
        var childUbounds = new double[3];
        int typeTotal = 1;
        int newType = 0, npA = 0, npB = 0;
        double newScore = 0, newLower = 0, newUpper = 0, nextLower = 0, nextUpper = 0;
        double gapOpenA = _open[0], gapOpenB = _open[3], gapExtendA = _extend[0], gapExtendB = _extend[3];
        double leftGapOpenA, leftGapOpenB, leftGapExtendA, leftGapExtendB, rightGapOpenA, rightGapOpenB, rightGapExtendA, rightGapExtendB;
        if (strand == '+')
        {
            leftGapOpenA = _open[1];
            leftGapOpenB = _open[4];
            leftGapExtendA = _extend[1];
            leftGapExtendB = _extend[4];
            rightGapOpenA = _open[2];
            rightGapOpenB = _open[5];
            rightGapExtendA = _extend[2];
            rightGapExtendB = _extend[5];
        }
        else
        {
            leftGapOpenA = _open[2];
            leftGapOpenB = _open[5];
            leftGapExtendA = _extend[2];
            leftGapExtendB = _extend[5];
            rightGapOpenA = _open[1];
            rightGapOpenB = _open[4];
            rightGapExtendA = _extend[1];
            rightGapExtendB = _extend[4];
        }

        var warnings = new List<string>();
        if (mismatch >= match)
            warnings.Add("Match score is less than mismatch score. Algorithm may return incorrect results.");
        if (_open[4] > mismatch || _open[3] > mismatch || _open[5] > mismatch || _open[1] > mismatch || _open[0] > mismatch || _open[2] > mismatch
            || _extend[4] > mismatch || _extend[3] > mismatch || _extend[5] > mismatch || _extend[1] > mismatch || _extend[0] > mismatch || _extend[2] > mismatch)
            warnings.Add("One or more gap scores are greater than mismatch score. Algorithm may return incorrect results.");
        Warnings = warnings;

        var matrix = new FogsaaCell[(nA + 1) * (nB + 1)];
        for (int i = 0; i < matrix.Length; i++)
            matrix[i] = new FogsaaCell();
        FogsaaCell M(int x, int y) => matrix[x * (nB + 1) + y];

        // FOGSAA_CALCULATE_SCORE
        void Bounds(double currScore, int currType, out double lower, out double upper, int pA, int pB)
        {
            if (nA - pA <= nB - pB)
            {
                if (pA == nA && currType == HorizontalMove)
                {
                    lower = currScore + rightGapExtendA * (nB - pB);
                    upper = currScore + rightGapExtendA * (nB - pB);
                }
                else
                {
                    lower = currScore + (nA - pA) * mismatch;
                    upper = currScore + (nA - pA) * match;
                    t = rightGapOpenA + rightGapExtendA * (nB - pB - (nA - pA) - 1);
                    t2 = gapExtendA * (nB - pB - (nA - pA));
                    if (currType == HorizontalMove && t2 > t)
                    {
                        lower += t2;
                        upper += t2;
                    }
                    else
                    {
                        lower += t;
                        upper += t;
                    }
                }
            }
            else
            {
                if (pB == nB && currType == VerticalMove)
                {
                    lower = currScore + rightGapExtendB * (nA - pA);
                    upper = currScore + rightGapExtendB * (nA - pA);
                }
                else
                {
                    lower = currScore + (nB - pB) * mismatch;
                    upper = currScore + (nB - pB) * match;
                    t = rightGapOpenB + rightGapExtendB * (nA - pA - (nB - pB) - 1);
                    t2 = gapExtendB * (nA - pA - (nB - pB));
                    if (currType == VerticalMove && t2 > t)
                    {
                        lower += t2;
                        upper += t2;
                    }
                    else
                    {
                        lower += t;
                        upper += t;
                    }
                }
            }
        }

        void Sort()
        {
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2 - i; j++)
                {
                    if (childLbounds[j] < childLbounds[j + 1] || (childLbounds[j] == childLbounds[j + 1] && childUbounds[j] < childUbounds[j + 1]))
                    {
                        (childLbounds[j], childLbounds[j + 1]) = (childLbounds[j + 1], childLbounds[j]);
                        (childTypes[j], childTypes[j + 1]) = (childTypes[j + 1], childTypes[j]);
                        (childUbounds[j], childUbounds[j + 1]) = (childUbounds[j + 1], childUbounds[j]);
                    }
                }
            }
        }

        var queue = new List<FogsaaNode>();
        bool HeapCond(int x, int y) => queue[x].NextUpper > queue[y].NextUpper || (queue[x].NextUpper == queue[y].NextUpper && queue[x].NextLower > queue[y].NextLower);

        void Insert(int pA, int pB, int typeTotalNode, int nextType, double nextLowerNode, double nextUpperNode)
        {
            queue.Add(new FogsaaNode { PA = pA, PB = pB, NextType = nextType, NextLower = nextLowerNode, TypeUpToNext = typeTotalNode, NextUpper = nextUpperNode });
            int i = queue.Count - 1;
            while (i != 0 && !HeapCond((i - 1) / 2, i))
            {
                (queue[i], queue[(i - 1) / 2]) = (queue[(i - 1) / 2], queue[i]);
                i = (i - 1) / 2;
            }
        }

        FogsaaNode Pop()
        {
            var root = queue[0];
            int size = queue.Count - 1;
            queue[0] = queue[size];
            queue.RemoveAt(size);
            int i = 0;
            while (true)
            {
                int largest = i;
                if (2 * i + 1 < size && !HeapCond(i, 2 * i + 1))
                    largest = 2 * i + 1;
                if (2 * i + 2 < size && !HeapCond(largest, 2 * i + 2))
                    largest = 2 * i + 2;
                if (largest != i)
                {
                    (queue[i], queue[largest]) = (queue[largest], queue[i]);
                    i = largest;
                }
                else
                {
                    break;
                }
            }

            return root;
        }

        var origin = M(0, 0);
        origin.PresentScore = 0;
        origin.Type = StartPoint;
        Bounds(origin.PresentScore, StartPoint, out origin.Lower, out origin.Upper, 0, 0);
        origin.IsLeftGap = true;
        lowerBound = origin.Lower;
        do
        {
            pathend = 1;
            while (curpA < nA || curpB < nB)
            {
                var curr = M(curpA, curpB);
                if (typeTotal is DiagonalMove or HorizontalMove or VerticalMove)
                {
                    // current is a 1st child
                    if (curpA <= nA - 1 && curpB <= nB - 1)
                    {
                        double p = alignScore(a[curpA], b[curpB]);
                        Bounds(curr.PresentScore + p, DiagonalMove, out childLbounds[0], out childUbounds[0], curpA + 1, curpB + 1);
                        if (curr.Type is DiagonalMove or StartPoint)
                        {
                            if (!curr.IsLeftGap)
                            {
                                Bounds(curr.PresentScore + gapOpenA, HorizontalMove, out childLbounds[1], out childUbounds[1], curpA, curpB + 1);
                                Bounds(curr.PresentScore + gapOpenB, VerticalMove, out childLbounds[2], out childUbounds[2], curpA + 1, curpB);
                            }
                            else
                            {
                                Bounds(curr.PresentScore + leftGapOpenA, HorizontalMove, out childLbounds[1], out childUbounds[1], curpA, curpB + 1);
                                Bounds(curr.PresentScore + leftGapOpenB, VerticalMove, out childLbounds[2], out childUbounds[2], curpA + 1, curpB);
                            }
                        }
                        else if (curr.Type == HorizontalMove)
                        {
                            if (!curr.IsLeftGap)
                            {
                                Bounds(curr.PresentScore + gapExtendA, HorizontalMove, out childLbounds[1], out childUbounds[1], curpA, curpB + 1);
                                Bounds(curr.PresentScore + gapOpenB, VerticalMove, out childLbounds[2], out childUbounds[2], curpA + 1, curpB);
                            }
                            else
                            {
                                Bounds(curr.PresentScore + leftGapExtendA, HorizontalMove, out childLbounds[1], out childUbounds[1], curpA, curpB + 1);
                                Bounds(curr.PresentScore + leftGapOpenB, VerticalMove, out childLbounds[2], out childUbounds[2], curpA + 1, curpB);
                            }
                        }
                        else
                        {
                            if (!curr.IsLeftGap)
                            {
                                Bounds(curr.PresentScore + gapOpenA, HorizontalMove, out childLbounds[1], out childUbounds[1], curpA, curpB + 1);
                                Bounds(curr.PresentScore + gapExtendB, VerticalMove, out childLbounds[2], out childUbounds[2], curpA + 1, curpB);
                            }
                            else
                            {
                                Bounds(curr.PresentScore + leftGapOpenA, HorizontalMove, out childLbounds[1], out childUbounds[1], curpA, curpB + 1);
                                Bounds(curr.PresentScore + leftGapExtendB, VerticalMove, out childLbounds[2], out childUbounds[2], curpA + 1, curpB);
                            }
                        }

                        childTypes[0] = DiagonalMove;
                        childTypes[1] = HorizontalMove;
                        childTypes[2] = VerticalMove;
                        Sort();
                        newType = childTypes[0];
                        if (newType == DiagonalMove)
                        {
                            npA = curpA + 1;
                            npB = curpB + 1;
                            newScore = curr.PresentScore + p;
                        }
                        else if (newType == HorizontalMove)
                        {
                            npA = curpA;
                            npB = curpB + 1;
                            if (curr.IsLeftGap)
                                newScore = curr.PresentScore + (curr.Type == HorizontalMove ? leftGapExtendA : leftGapOpenA);
                            else
                                newScore = curr.PresentScore + (curr.Type == HorizontalMove ? gapExtendA : gapOpenA);
                        }
                        else
                        {
                            npA = curpA + 1;
                            npB = curpB;
                            if (curr.IsLeftGap)
                                newScore = curr.PresentScore + (curr.Type == VerticalMove ? leftGapExtendB : leftGapOpenB);
                            else
                                newScore = curr.PresentScore + (curr.Type == VerticalMove ? gapExtendB : gapOpenB);
                        }

                        if (childUbounds[1] >= M(0, 0).Lower)
                            Insert(curpA, curpB, newType + childTypes[1], childTypes[1], childLbounds[1], childUbounds[1]);
                    }
                    else if (curpA <= nA - 1)
                    {
                        newType = VerticalMove;
                        npA = curpA + 1;
                        npB = curpB;
                        newScore = curr.PresentScore + (curr.Type == VerticalMove ? rightGapExtendB : rightGapOpenB);
                    }
                    else
                    {
                        newType = HorizontalMove;
                        npA = curpA;
                        npB = curpB + 1;
                        newScore = curr.PresentScore + (curr.Type == HorizontalMove ? rightGapExtendA : rightGapOpenA);
                    }
                }
                else if (typeTotal is DiagonalMove + HorizontalMove or DiagonalMove + VerticalMove or HorizontalMove + VerticalMove)
                {
                    // current is a 2nd child (sum of two types)
                    if (newType == DiagonalMove)
                    {
                        npA = curpA + 1;
                        npB = curpB + 1;
                        newScore = curr.PresentScore + (a[curpA] == b[curpB] ? match : mismatch);
                        if (7 - typeTotal == HorizontalMove)
                        {
                            if (curr.Type != HorizontalMove)
                                Bounds(curr.PresentScore + (curr.IsLeftGap ? leftGapOpenA : gapOpenA), HorizontalMove, out nextLower, out nextUpper, curpA, curpB + 1);
                            else
                                Bounds(curr.PresentScore + (curr.IsLeftGap ? leftGapExtendA : gapExtendA), HorizontalMove, out nextLower, out nextUpper, curpA, curpB + 1);
                        }
                        else
                        {
                            if (curr.Type != VerticalMove)
                                Bounds(curr.PresentScore + (curr.IsLeftGap ? leftGapOpenB : gapOpenB), VerticalMove, out nextLower, out nextUpper, curpA, curpB + 1);
                            else
                                Bounds(curr.PresentScore + (curr.IsLeftGap ? leftGapExtendB : gapExtendB), VerticalMove, out nextLower, out nextUpper, curpA, curpB + 1);
                        }
                    }
                    else if (newType == HorizontalMove)
                    {
                        npA = curpA;
                        npB = curpB + 1;
                        newScore = curr.PresentScore + (curr.Type == HorizontalMove ? gapExtendA : gapOpenA);
                        if (7 - typeTotal == DiagonalMove)
                        {
                            Bounds(curr.PresentScore + alignScore(a[curpA], b[curpB]), DiagonalMove, out nextLower, out nextUpper, curpA + 1, curpB + 1);
                        }
                        else if (curr.Type != VerticalMove)
                        {
                            Bounds(curr.PresentScore + (curr.IsLeftGap ? leftGapOpenB : gapOpenB), VerticalMove, out nextLower, out nextUpper, curpA, curpB + 1);
                        }
                        else
                        {
                            Bounds(curr.PresentScore + (curr.IsLeftGap ? leftGapExtendB : gapExtendB), VerticalMove, out nextLower, out nextUpper, curpA, curpB + 1);
                        }
                    }
                    else
                    {
                        npA = curpA + 1;
                        npB = curpB;
                        newScore = curr.PresentScore + (curr.Type == VerticalMove ? gapExtendB : gapOpenB);
                        if (7 - typeTotal == DiagonalMove)
                        {
                            Bounds(curr.PresentScore + alignScore(a[curpA], b[curpB]), DiagonalMove, out nextLower, out nextUpper, curpA + 1, curpB + 1);
                        }
                        else if (curr.Type != HorizontalMove)
                        {
                            Bounds(curr.PresentScore + (curr.IsLeftGap ? leftGapOpenA : gapOpenA), HorizontalMove, out nextLower, out nextUpper, curpA, curpB + 1);
                        }
                        else
                        {
                            Bounds(curr.PresentScore + (curr.IsLeftGap ? leftGapExtendA : gapExtendA), HorizontalMove, out nextLower, out nextUpper, curpA, curpB + 1);
                        }
                    }

                    if (nextUpper >= M(0, 0).Lower)
                        Insert(curpA, curpB, 7, 7 - typeTotal, nextLower, nextUpper);
                }
                else if (typeTotal == DiagonalMove + HorizontalMove + VerticalMove)
                {
                    // current is a 3rd child
                    if (newType == DiagonalMove)
                    {
                        npA = curpA + 1;
                        npB = curpB + 1;
                        newScore = curr.PresentScore + alignScore(a[curpA], b[curpB]);
                    }
                    else if (newType == HorizontalMove)
                    {
                        npA = curpA;
                        npB = curpB + 1;
                        if (curr.Type != HorizontalMove)
                            newScore = curr.PresentScore + (curr.IsLeftGap ? leftGapOpenA : gapOpenA);
                        else
                            newScore = curr.PresentScore + (curr.IsLeftGap ? leftGapExtendA : gapExtendA);
                    }
                    else
                    {
                        npA = curpA + 1;
                        npB = curpB;
                        if (curr.Type != VerticalMove)
                            newScore = curr.PresentScore + (curr.IsLeftGap ? leftGapOpenB : gapOpenB);
                        else
                            newScore = curr.PresentScore + (curr.IsLeftGap ? leftGapExtendB : gapExtendB);
                    }
                }

                // write the new node to the matrix, but skip if there's already a better path there
                var target = M(npA, npB);
                if (target.Filled && target.Type <= 4 && target.PresentScore >= newScore)
                {
                    pathend = 0;
                    break;
                }

                Bounds(newScore, newType, out newLower, out newUpper, npA, npB);
                target.PresentScore = newScore;
                target.Lower = newLower;
                target.Upper = newUpper;
                target.Type = newType;
                target.Filled = true;
                target.IsLeftGap = newType is HorizontalMove or VerticalMove && curr.IsLeftGap;

                // make the child the new current node
                curpA = npA;
                curpB = npB;
                typeTotal = 1;
                if (target.Upper < lowerBound && lowerBound - target.Upper > Epsilon)
                {
                    pathend = 0;
                    break;
                }
            }

            var here = M(curpA, curpB);
            if (here.PresentScore > lowerBound && here.PresentScore - lowerBound > Epsilon && pathend == 1)
                lowerBound = here.PresentScore;
            if (queue.Count > 0)
            {
                var root = Pop();
                curpA = root.PA;
                curpB = root.PB;
                typeTotal = root.TypeUpToNext;
                newUpper = root.NextUpper;
                newType = root.NextType;
            }
            else
            {
                break;
            }
        }
        while (lowerBound < newUpper && newUpper - lowerBound > Epsilon);

        if (lowerBound < newUpper && newUpper - lowerBound > Epsilon)
            throw new InvalidOperationException("Algorithm ended incomplete. Report this as a bug.");
        var columns = new List<(char A, char B)>();
        int x = nA, y = nB;
        while (true)
        {
            var cell = M(x, y);
            if (cell.Type is 0 or StartPoint)
                break;
            if (cell.Type == DiagonalMove)
            {
                columns.Add((a[x - 1], b[y - 1]));
                x--;
                y--;
            }
            else if (cell.Type == HorizontalMove)
            {
                columns.Add(('-', b[y - 1]));
                y--;
            }
            else
            {
                columns.Add((a[x - 1], '-'));
                x--;
            }
        }

        columns.Reverse();
        return (M(nA, nB).PresentScore, columns);
    }
}
