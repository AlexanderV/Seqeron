namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: treeOperation.c (the branch weights of the tree-dependent iterative refinement)

    /// <summary>
    /// The unrooted tree of <c>treeOperation.c</c> (<c>stopol</c>): the 2 n nodes (n leaves, then the n - 1 inner nodes of the rooted guide tree, the root merged into one branch) with their three neighbours, the lengths of the
    /// branches, the members behind each branch and the weight of each branch. <c>checkMinusLength</c> raises branch lengths below 0.001 (the program changes <c>len</c> in place).
    /// </summary>
    internal sealed class MafftUnrootedTree
    {
        private const double EfThreeWay = 1.0, MaxBw = 1.0, MinBw = 0.01, MinLen = 0.001;

        internal sealed class Node
        {
            internal readonly int[] Children = { -1, -1, -1 };
            internal readonly double[] Length = new double[3];
            internal readonly (int Step, int LorR)?[] WeightRef = new (int, int)?[3];
            internal readonly int[]?[] Members = new int[]?[3];
        }

        private readonly Node[] _node;
        private readonly int[][][] _topol; // topol[step][0 / 1] = the leaves of the two groups of a step (complete lists)
        private readonly double[][] _bw;
        private readonly int _nseq;

        /// <summary><c>treeCnv</c>: builds the unrooted tree from the rooted <paramref name="topol"/> / <paramref name="len"/> and fills the branch weights <paramref name="bw"/> (<c>calcBranchWeight</c>).</summary>
        internal MafftUnrootedTree(int nseq, int[][][] topol, double[][] len, double[][] bw)
        {
            _nseq = nseq;
            _topol = topol;
            _bw = bw;
            _node = new Node[nseq * 2];
            for (int i = 0; i < _node.Length; i++)
                _node[i] = new Node();
            if (nseq <= 2)
                return;
            var count = new int[2 * nseq];
            for (int i = 0; i < nseq - 1; i++)
                for (int j = 0; j < 2; j++)
                    if (len[i][j] < MinLen)
                        len[i][j] = MinLen;
            for (int i = nseq; i < nseq * 2; i++)
            {
                int rep = i - nseq;
                var (step, lr) = SearchParent(rep, 0, nseq - 1);
                int cc = count[step];
                _node[step].Length[cc] = len[step][lr];
                _node[step].WeightRef[cc] = (step, lr);
                _node[step].Children[cc] = i;
                _node[step].Members[cc] = topol[step][lr];
                count[step]++;
                cc = count[i];
                _node[i].Length[cc] = len[step][lr];
                _node[i].WeightRef[cc] = (step, lr);
                _node[i].Children[cc] = step;
                _node[i].Members[cc] = topol[step][lr];
                count[i]++;
            }
            for (int i = 0; i < nseq - 2; i++)
            {
                int rep = Math.Min(topol[i][0][0], topol[i][1][0]);
                var (step, lr) = SearchParent(rep, i + 1, nseq - 1);
                int cc = count[step];
                _node[step].Length[cc] = len[step][lr];
                _node[step].WeightRef[cc] = (step, lr);
                _node[step].Children[cc] = i;
                _node[step].Members[cc] = topol[step][lr];
                count[step]++;
                cc = count[i];
                _node[i].Length[cc] = len[step][lr];
                _node[i].WeightRef[cc] = (step, lr);
                _node[i].Children[cc] = step;
                _node[i].Members[cc] = NegativeMember(topol[step][lr], nseq);
                count[i]++;
            }
            int idx;
            if (_node[nseq - 2].Children[0] == nseq - 3)
                idx = 1;
            else if (_node[nseq - 2].Children[1] == nseq - 3)
                idx = 0;
            else
                throw new InvalidOperationException("ERROR in stopol ?");
            _node[nseq - 3].Length[2] = len[nseq - 2][0] + len[nseq - 2][1];
            _node[nseq - 3].WeightRef[2] = (nseq - 2, 0);
            _node[nseq - 3].Children[2] = _node[nseq - 2].Children[idx];
            int tmpint = _node[nseq - 2].Children[idx];
            _node[tmpint].Children[2] = nseq - 3;
            _node[tmpint].Length[2] = len[nseq - 2][0] + len[nseq - 2][1];
            _node[tmpint].WeightRef[2] = (nseq - 2, 0);
            CalcBranchWeight(len);
        }

        private (int Step, int LorR) SearchParent(int top, int start, int end)
        {
            for (int i = start; i < end; i++)
                for (int j = 0; j < 2; j++)
                    if (Array.IndexOf(_topol[i][j], top) >= 0)
                        return (i, j);
            throw new InvalidOperationException("Error in searchParent");
        }

        private static int[] NegativeMember(int[] query, int locnseq)
        {
            var tmp = new bool[locnseq];
            foreach (int n in query)
                tmp[n] = true;
            var res = new List<int>();
            for (int i = 0; i < locnseq; i++)
                if (!tmp[i])
                    res.Add(i);
            return res.ToArray();
        }

        private bool IsLeaf(int n) => _node[n].Children[1] == -1;

        private double SyntheticLength(int ob, int opposite)
        {
            if (IsLeaf(ob))
                return _node[ob].Length[0];
            var dirCh = new int[3];
            int count = 0, dirPa = -10;
            for (int i = 0; i < 3; i++)
            {
                if (opposite != _node[ob].Children[i])
                    dirCh[count++] = i;
                else
                    dirPa = i;
            }
            if (count != 2)
                throw new InvalidOperationException("Invalid call");
            double tmpvalue0 = SyntheticLength(_node[ob].Children[dirCh[0]], ob);
            double tmpvalue1 = SyntheticLength(_node[ob].Children[dirCh[1]], ob);
            bool nanflag = false;
            if (tmpvalue0 != 0.0)
                tmpvalue0 = 1.0 / tmpvalue0;
            else
                nanflag = true;
            if (tmpvalue1 != 0.0)
                tmpvalue1 = 1.0 / tmpvalue1;
            else
                nanflag = true;
            double value;
            if (nanflag)
            {
                value = 0.0;
            }
            else
            {
                value = tmpvalue0 + tmpvalue1;
                value = 1.0 / value;
            }
            value += _node[ob].Length[dirPa];
            return value;
        }

        private double CalcW(int ob, int op)
        {
            if (IsLeaf(ob))
                return 1.0;
            var dirCh = new int[3];
            int count = 0, dirPa = -10;
            for (int i = 0; i < 3; i++)
            {
                if (op != _node[ob].Children[i])
                    dirCh[count++] = i;
                else
                    dirPa = i;
            }
            if (count != 2)
                throw new InvalidOperationException("Invalid call of calcW");
            double a = SyntheticLength(_node[ob].Children[dirCh[0]], ob);
            double b = SyntheticLength(_node[ob].Children[dirCh[1]], ob);
            double c = SyntheticLength(_node[ob].Children[dirPa], ob);
            if (c == 0.0)
                return MaxBw;
            if (a == 0.0 || b == 0.0)
                return MinBw;
            double f = EfThreeWay;
            double s = b * c + c * a + a * b;
            double value = a * b * (c + a) * (c + b) / (c * (a + b) * f * s);
            return Math.Sqrt(value);
        }

        private void CalcBranchWeight(double[][] len)
        {
            int nseq = _nseq;
            _ = len;
            for (int i = nseq; i < nseq * 2; i++)
            {
                int rep = i - nseq;
                var (step, lr) = SearchParent(rep, 0, nseq - 1);
                if (step == nseq - 2)
                    continue;
                int topNode = step, btmNode = i;
                double topW = CalcW(topNode, btmNode);
                double btmW = CalcW(btmNode, topNode);
                _bw[step][lr] = topW * btmW;
            }
            int ii;
            for (ii = 0; ii < nseq - 3; ii++)
            {
                int rep = Math.Min(_topol[ii][0][0], _topol[ii][1][0]);
                var (step, lr) = SearchParent(rep, ii + 1, nseq - 1);
                if (step == nseq - 2)
                    continue;
                int topNode = step, btmNode = ii;
                double topW = CalcW(topNode, btmNode);
                double btmW = CalcW(btmNode, topNode);
                _bw[step][lr] = topW * btmW;
            }
            {
                int topNode = _node[nseq - 3].Children[2];
                int btmNode = ii;
                double topW = CalcW(topNode, btmNode);
                double btmW = CalcW(btmNode, topNode);
                _bw[nseq - 2][0] = topW * btmW;
                _bw[nseq - 2][1] = 1.0;
            }
        }

        private void FindNodes(int step, int lr, out int topNode, out int btmNode)
        {
            int nseq = _nseq;
            if (step == nseq - 2)
            {
                topNode = _node[nseq - 2].Children[0];
                btmNode = nseq - 3;
            }
            else
            {
                int i;
                for (i = 0; i < 3; i++)
                    if (_node[step].Members[i]![0] == _topol[step][lr][0])
                        break;
                if (i == 3)
                    throw new InvalidOperationException("Incorrect call of distFromABranch.");
                btmNode = _node[step].Children[i];
                topNode = step;
            }
        }

        private void DistFromABranchRec(double[] result, int ob, int op, double[][] len)
        {
            _ = len;
            if (IsLeaf(ob))
                return;
            var dirCh = new int[3];
            int count = 0;
            for (int i = 0; i < 3; i++)
                if (_node[ob].Children[i] != op)
                    dirCh[count++] = i;
            if (count != 2)
                throw new InvalidOperationException("Incorrect call of distFromABranch_rec");
            foreach (int n in _node[ob].Members[dirCh[0]]!)
                result[n] += _node[ob].Length[dirCh[0]];
            DistFromABranchRec(result, _node[ob].Children[dirCh[0]], ob, len);
            foreach (int n in _node[ob].Members[dirCh[1]]!)
                result[n] += _node[ob].Length[dirCh[1]];
            DistFromABranchRec(result, _node[ob].Children[dirCh[1]], ob, len);
        }

        /// <summary><c>distFromABranch</c>: the distance of every sequence from the branch (step, group).</summary>
        internal void DistFromABranch(double[] result, double[][] len, int step, int lr)
        {
            if (_nseq == 2)
            {
                result[0] = len[0][0];
                result[1] = len[0][1];
                return;
            }
            FindNodes(step, lr, out int topNode, out int btmNode);
            for (int i = 0; i < _nseq; i++)
                result[i] = 0.0;
            DistFromABranchRec(result, btmNode, topNode, len);
            DistFromABranchRec(result, topNode, btmNode, len);
        }

        private void WeightFromABranchRec(double[] result, int ob, int op)
        {
            if (IsLeaf(ob))
                return;
            var dirCh = new int[3];
            int count = 0;
            for (int i = 0; i < 3; i++)
                if (_node[ob].Children[i] != op)
                    dirCh[count++] = i;
            if (count != 2)
                throw new InvalidOperationException("Incorrect call of weightFromABranch_rec");
            var w0 = _node[ob].WeightRef[dirCh[0]]!.Value;
            foreach (int n in _node[ob].Members[dirCh[0]]!)
                result[n] *= _bw[w0.Step][w0.LorR];
            WeightFromABranchRec(result, _node[ob].Children[dirCh[0]], ob);
            var w1 = _node[ob].WeightRef[dirCh[1]]!.Value;
            foreach (int n in _node[ob].Members[dirCh[1]]!)
                result[n] *= _bw[w1.Step][w1.LorR];
            WeightFromABranchRec(result, _node[ob].Children[dirCh[1]], ob);
        }

        /// <summary><c>weightFromABranch</c>: the weight of every sequence seen from the branch (step, group): the product of the weights of the branches between it and this branch.</summary>
        internal void WeightFromABranch(double[] result, int step, int lr)
        {
            FindNodes(step, lr, out int topNode, out int btmNode);
            for (int i = 0; i < _nseq; i++)
                result[i] = 1.0;
            WeightFromABranchRec(result, btmNode, topNode);
            WeightFromABranchRec(result, topNode, btmNode);
        }
    }

    #endregion
}
