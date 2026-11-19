namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    // ClustalW 2.1 `-ITERATION=TREE | ALIGNMENT` / `-NUMITER=n` (audit A2-14, B09 F102): src/multipleAlign/Iteration.cpp (removeFirstIterate, iterationOnTreeNode), src/alignment/ObjectiveScore.cpp,
    // the parts of src/alignment/Alignment.cpp they use (cutSelectedSequencesFromAlignment, pasteSequencesIntoPosition, removeGapOnlyColsFromSelectedSeqs, countid, addSeqWeight),
    // MSA::calcPairwiseForProfileAlign / doProfileAlign and TreeInterface::getWeightsForProfileAlignNJ of clustalw 2.1+lgpl-7 (Debian).
    //
    // What the program does ("remove first" iteration): for iteration n = 1..NUMITER, for every sequence i of the finished alignment: take a copy, cut sequence i, paste it in front, strip the gaps
    // of that single sequence (the one-sequence "profile 1" has every column removed in which it has a gap) and the all-gap columns of the others (profile 2), compute the percent identities of
    // the rows (countid), do an initial profile alignment of the two profiles with equal weights (MSA::calcPairwiseForProfileAlign, which also sets distanceTree = false), recompute the identities
    // from that alignment, weight the sequences (profile 1: 1000; profile 2: the neighbour-joining tree weights of the sequences of profile 2), convert the distances to similarities and
    // re-align the two profiles (MSA::doProfileAlign, weights = max identity to the other profile x tree weight, normalised to 1000). The realignment is kept when its ObjectiveScore is LOWER than
    // the best so far (the program compares `score < bestScore`; the score adds the gap penalties as positive numbers to a BLOSUM45 (x100) sum of residue pairs). The kept alignment is the
    // best single realignment: alnPtr is only updated after the loops, so a second iteration repeats the first and stops ("Score was not improved in last iteration").
    // -ITERATION=TREE does the same on the sequences of every profile alignment of the progressive stage, in an Alignment object that has no weights (all 100) and in which every sequence has
    // the same name, so the names of the profile-2 tree leaves match the alignment's names by the LAST leaf (Tree::readTree).

    /// <summary>Alignment object of Alignment.cpp reduced to what the iteration reads: rows of residue codes (30 / 31 = gap), unique ids, the weights (numSeqs + 1 entries) and profile1NumSeqs.</summary>
    private sealed class CwAlignment
    {
        public List<int[]> Rows = [];
        public List<long> Ids = [];
        public int[] Weights = [];
        public int Profile1NumSeqs;

        public int NumSeqs => Rows.Count;

        public CwAlignment Clone() => new() { Rows = [.. Rows.Select(r => (int[])r.Clone())], Ids = [.. Ids], Weights = (int[])Weights.Clone(), Profile1NumSeqs = Profile1NumSeqs };

        public static bool IsGap(int code) => code == CwGapPos1 || code == CwGapPos2;

        // Alignment::countid(s1, s2), 1-based sequence numbers
        public float CountId(int s1, int s2)
        {
            var a = Rows[s1 - 1];
            var b = Rows[s2 - 1];
            int shorter = Math.Min(a.Length, b.Length);
            int count = 0, total = 0;
            for (int i = 0; i < shorter; i++)
            {
                int c1 = a[i], c2 = b[i];
                if (c1 >= 0 && c1 < CwMaxAa)
                {
                    total++;
                    if (c1 == c2)
                        count++;
                }
            }

            return total == 0 ? 0f : (float)(100.0 * count / total);
        }

        // Alignment::cutSelectedSequencesFromAlignment(selected = {i}) + pasteSequencesIntoPosition(.., 0) with profileNum == 0: the sequence moves to the front; all weights are reset to 100
        public void MoveToFront(int i)
        {
            var row = Rows[i - 1];
            long id = Ids[i - 1];
            Rows.RemoveAt(i - 1);
            Ids.RemoveAt(i - 1);
            Rows.Insert(0, row);
            Ids.Insert(0, id);
            Weights = new int[NumSeqs + 1];
            Array.Fill(Weights, 100);
        }

        // Alignment::removeGapOnlyColsFromSelectedSeqs for the contiguous selection first..last (1-based): columns in which every selected sequence has a gap are deleted from the selected rows
        public void RemoveGapOnlyColumns(int first, int last)
        {
            int numColumns = last - first + 1;
            for (int i = 1; i <= Rows[first - 1].Length;)
            {
                int numGaps = 0;
                for (int j = first; j <= last; j++)
                {
                    if (Rows[j - 1].Length >= i && IsGap(Rows[j - 1][i - 1]))
                        numGaps++;
                }

                if (numGaps == numColumns)
                {
                    for (int j = first; j <= last; j++)
                    {
                        var row = Rows[j - 1];
                        // for k = i + 1 .. len: row[k - 1] = row[k]; then pop_back (a row shorter than i only loses its last element, which cannot happen here: such a row prevents numGaps == numColumns)
                        if (row.Length >= i)
                            Rows[j - 1] = [.. row[..(i - 1)], .. row[i..]];
                        else if (row.Length > 0)
                            Rows[j - 1] = row[..^1];
                        if (Rows[first - 1].Length <= 0)
                            break;
                    }
                }
                else
                {
                    i++;
                }
            }
        }
    }

    // ObjectiveScore: SubMatrix::getAlnScoreMatrix = getMatrix(blosum45mt, defaultAAXref, matrix, negFlag = true, scale 100): BLOSUM45 x 100 over the residue codes, no shift, gap codes 0 (also for DNA runs)
    private static readonly Lazy<int[,]> CwObjectiveMatrix = new(() => BuildClustalWProteinMatrix(Seqeron.Genomics.Infrastructure.SubstitutionMatrix.ClustalWBlosum45, 100).Matrix);

    /// <summary>ObjectiveScore::getScore: the weighted sum over the pairs of sequences of the BLOSUM45 score of the columns between the first and last non-double-gap column plus the gap penalties (opening for the first gap column of a run in either sequence, extension afterwards), in single precision, truncated to a long.</summary>
    internal static long ClustalWObjectiveScore(IReadOnlyList<int[]> rows, IReadOnlyList<int> weights, float gapOpen, float gapExtend)
    {
        var matrix = CwObjectiveMatrix.Value;
        int size = weights.Count;
        int sumWeights = 0;
        for (int i = 0; i < size - 1; i++)
            sumWeights += weights[i];
        var norm = new float[size];
        for (int i = 0; i < size - 1; i++)
            norm[i] = (float)weights[i] / (float)sumWeights;

        int numSeqs = rows.Count;
        float scoreTotal = 0.0f;
        for (int seq1 = 1; seq1 <= numSeqs && seq1 <= size; seq1++)
        {
            float w1 = norm[seq1 - 1];
            for (int seq2 = seq1 + 1; seq2 <= numSeqs && seq2 <= size; seq2++)
            {
                float w2 = norm[seq2 - 1];
                float weight = w1 * w2;
                var a = rows[seq1 - 1];
                var b = rows[seq2 - 1];
                float pwScoreLetters = 0f, pwScoreGaps = 0f;
                if (a.Length == b.Length)
                {
                    int numCols = a.Length;
                    int colStart = 1;
                    for (int col = 1; col < numCols; col++)
                    {
                        if (!CwAlignment.IsGap(a[col - 1]) || !CwAlignment.IsGap(b[col - 1]))
                        {
                            colStart = col;
                            break;
                        }
                    }

                    int colEnd = numCols;
                    for (int col = numCols; col >= 1; col--)
                    {
                        if (!CwAlignment.IsGap(a[col - 1]) || !CwAlignment.IsGap(b[col - 1]))
                        {
                            colEnd = col;
                            break;
                        }
                    }

                    for (int col = colStart; col <= colEnd; col++)
                        pwScoreLetters += matrix[a[col - 1], b[col - 1]];

                    bool inGap1 = false, inGap2 = false;
                    for (int col = colStart; col <= colEnd; col++)
                    {
                        bool gap1 = CwAlignment.IsGap(a[col - 1]);
                        bool gap2 = CwAlignment.IsGap(b[col - 1]);
                        if (gap1 && gap2)
                            continue;
                        if (gap1)
                        {
                            if (!inGap1)
                            {
                                pwScoreGaps += gapOpen;
                                inGap1 = true;
                            }
                            else
                            {
                                pwScoreGaps += gapExtend;
                            }

                            continue;
                        }

                        if (gap2)
                        {
                            if (!inGap2)
                            {
                                pwScoreGaps += gapOpen;
                                inGap2 = true;
                            }
                            else
                            {
                                pwScoreGaps += gapExtend;
                            }

                            continue;
                        }

                        inGap1 = inGap2 = false;
                    }
                }

                float pwScore = pwScoreLetters + pwScoreGaps;
                scoreTotal += weight * pwScore;
            }
        }

        return (long)scoreTotal;
    }

    private sealed partial class ClustalWMsaState
    {
        /// <summary><c>-NEWTREE</c> / <c>-NEWTREE1</c> / <c>-NEWTREE2</c>: the run stops after the guide tree(s) are written.</summary>
        public bool GuideTreesOnly { get; init; }

        private bool _distanceTree = true;   // UserParameters::distanceTree (Tree::readTree sets it, MSA::calcPairwiseForProfileAlign clears it)
        private bool _iterating;             // doRemoveFirstIteration is NONE while removeFirstIterate runs

        /// <summary>The "FINAL score" the program prints at the end of every <c>removeFirstIterate</c> call (the ObjectiveScore of the best re-alignment), in order.</summary>
        public List<long> IterationFinalScores { get; } = [];

        // Clustal::align: Iteration::removeFirstIterate(&alignmentObj) on the finished alignment (rows = _seq, weights = the weights of the last profile alignment)
        private void RunAlignmentIteration()
        {
            var aln = new CwAlignment { Rows = [.. _seq.Select(r => (int[])r.Clone())], Ids = [.. Enumerable.Range(1, _n).Select(i => (long)i)], Weights = (int[])_alnWeight.Clone() };
            if (!RemoveFirstIterate(aln, duplicateNames: false))
                return;
            for (int i = 0; i < _n; i++)
                _seq[i] = aln.Rows[i];
        }

        // MyersMillerProfileAlign::profileAlign -> Iteration::iterationOnTreeNode: the rows of the finished profile alignment (seqArray, group 1 then group 2) are iterated in a fresh Alignment
        private void IterationOnTreeNode(int numSeqsProf1, int numSeqsProf2, int[][] seqArray)
        {
            int numSeqsInProfiles = numSeqsProf1 + numSeqsProf2;
            if (numSeqsInProfiles <= 2)
                return;
            var aln = new CwAlignment { Weights = new int[numSeqsInProfiles + 1] };
            Array.Fill(aln.Weights, 100);
            for (int j = 0; j < numSeqsInProfiles && j < seqArray.Length; j++)
            {
                aln.Rows.Add((int[])seqArray[j].Clone());
                aln.Ids.Add(j + 1);
            }

            if (!RemoveFirstIterate(aln, duplicateNames: true))
                return;
            int newPrf1Length = 0;
            for (int j = 0; j < numSeqsInProfiles; j++)
            {
                if (j < numSeqsProf1)
                    newPrf1Length = Math.Max(newPrf1Length, aln.Rows[j].Length);
            }

            for (int j = 0; j < numSeqsInProfiles; j++)
            {
                var row = new int[newPrf1Length];
                Array.Fill(row, CwGapPos2);
                Array.Copy(aln.Rows[j], row, Math.Min(aln.Rows[j].Length, newPrf1Length));
                seqArray[j] = row;
            }
        }

        // Iteration::removeFirstIterate
        private bool RemoveFirstIterate(CwAlignment aln, bool duplicateNames)
        {
            int nSeqs = aln.NumSeqs;
            if (nSeqs <= 2)
                return false;
            bool savedIterating = _iterating;
            _iterating = true;
            float gapOpen = (float)_o.GapOpen, gapExtend = (float)_o.GapExtend;
            long bestScore = ClustalWObjectiveScore(aln.Rows, aln.Weights, gapOpen, gapExtend);
            bool scoreImprovedAnyIteration = false;
            CwAlignment? bestAlignSoFar = null;

            for (int iteration = 1; iteration <= _o.NumIterations; iteration++)
            {
                bool scoreImproved = false;
                for (int i = 1; i <= nSeqs; i++)
                {
                    AlignmentCancellation.ThrowIfRequested();
                    var iterateAlign = aln.Clone();
                    iterateAlign.Profile1NumSeqs = 1;
                    iterateAlign.MoveToFront(i);
                    iterateAlign.RemoveGapOnlyColumns(1, 1);
                    iterateAlign.RemoveGapOnlyColumns(2, nSeqs);

                    var distMat = new double[nSeqs + 1, nSeqs + 1];
                    if (nSeqs - 1 >= 2)
                    {
                        for (int a = 1; a <= nSeqs; a++)
                        {
                            for (int b = a + 1; b <= nSeqs; b++)
                            {
                                double dscore = iterateAlign.CountId(a, b);
                                distMat[a, b] = distMat[b, a] = (100.0 - dscore) / 100.0;
                            }
                        }
                    }

                    var (prof1Weight, prof2Weight) = GetWeightsForProfileAlign(iterateAlign, distMat, duplicateNames);
                    DoProfileAlign(iterateAlign, distMat, prof1Weight, prof2Weight);
                    long score = ClustalWObjectiveScore(iterateAlign.Rows, iterateAlign.Weights, gapOpen, gapExtend);
                    iterateAlign.Profile1NumSeqs = 0;
                    if (score < bestScore)
                    {
                        bestScore = score;
                        bestAlignSoFar = iterateAlign;
                        scoreImproved = true;
                        scoreImprovedAnyIteration = true;
                    }
                }

                if (!scoreImproved)
                    break;
            }

            if (scoreImprovedAnyIteration)
            {
                for (int i = 1; i <= nSeqs; i++)
                {
                    int at = bestAlignSoFar!.Ids.IndexOf(aln.Ids[i - 1]);
                    aln.Rows[i - 1] = bestAlignSoFar.Rows[at];
                }
            }

            IterationFinalScores.Add(bestScore);
            _iterating = savedIterating;
            return true;
        }

        // profileAlign on the rows of an Alignment object (the state's members swapped for the call)
        private void RunProfileAlign(CwAlignment aln, double[,] distMat, int[] group, int[] aligned)
        {
            var (n, seq, alnWeight, dist) = (_n, _seq, _alnWeight, _distMat);
            try
            {
                _n = aln.NumSeqs;
                _seq = aln.Rows.ToArray();
                _alnWeight = (int[])aln.Weights.Clone();
                _distMat = distMat;
                ProfileAlign(group, aligned);
                aln.Rows = [.. _seq];
            }
            finally
            {
                (_n, _seq, _alnWeight, _distMat) = (n, seq, alnWeight, dist);
            }
        }

        // TreeInterface::getWeightsForProfileAlignNJ (-CLUSTERING=NJ): profile trees from the identities, the initial profile alignment of MSA::calcPairwiseForProfileAlign (equal weights that never reach the
        // Alignment object: the vector has numSeqs entries, the Alignment's numSeqs + 1), the identities of that alignment, the tree weights of each profile, distances -> similarities
        private (int[] Prof1, int[] Prof2) GetWeightsForProfileAlign(CwAlignment aln, double[,] distMat, bool duplicateNames)
        {
            int numSeqs = aln.NumSeqs;
            int p1 = aln.Profile1NumSeqs;
            if (_o.Upgma)
                return GetWeightsForProfileAlignUpgma(aln, distMat);

            // Clustal::profileAlign: a user tree of a profile with fewer than two sequences is ignored
            bool useTree1 = UserTree1 != null && p1 >= 2;
            bool useTree2 = UserTree2 != null && numSeqs - p1 >= 2;

            // the profile guide trees are written from the distances BEFORE the initial alignment
            CwRootedTree? tree1 = null, tree2 = null;
            if (!useTree1 && p1 >= 2)
            {
                tree1 = ProfileTree(distMat, 1, p1, out string text1);
                GuideTree1Text = text1;
            }

            if (!useTree2 && numSeqs - p1 >= 2)
            {
                tree2 = ProfileTree(distMat, p1 + 1, numSeqs - p1, out string text2);
                GuideTree2Text = text2;
            }

            if (GuideTreesOnly)
            {
                StoppedAfterTrees = true;
                return ([], []);
            }

            // MSA::calcPairwiseForProfileAlign
            _distanceTree = false;
            var group = new int[numSeqs + 1];
            for (int i = 1; i <= p1; i++)
                group[i] = 1;
            for (int i = p1 + 1; i <= numSeqs; i++)
                group[i] = 2;
            var aligned = new int[numSeqs + 1];
            Array.Fill(aligned, 1, 1, numSeqs);
            RunProfileAlign(aln, distMat, group, aligned);
            for (int i = 1; i <= numSeqs; i++)
            {
                for (int j = i + 1; j <= numSeqs; j++)
                {
                    float dscore = aln.CountId(i, j);
                    distMat[i, j] = (100.0 - (double)dscore) / 100.0;
                    distMat[j, i] = distMat[i, j];
                }
            }

            var prof1 = new int[p1];
            var prof2 = new int[numSeqs];
            if (p1 >= 2)
            {
                _distanceTree = true;   // Tree::readTree
                if (useTree1)
                    tree1 = ReadUserProfileTree(UserTree1!, 0, p1);
            }

            var w1 = TreeSeqWeights(tree1, p1, p1, duplicateNames);
            Array.Copy(w1, prof1, p1);
            if (numSeqs - p1 >= 2)
            {
                _distanceTree = true;
                if (useTree2)
                    tree2 = ReadUserProfileTree(UserTree2!, p1, numSeqs - p1);
            }

            var w2 = TreeSeqWeights(tree2, numSeqs - p1, numSeqs, duplicateNames);
            Array.Copy(w2, 0, prof2, p1, numSeqs - p1);

            // distances -> similarities
            for (int i = 1; i < numSeqs; i++)
            {
                for (int j = i + 1; j <= numSeqs; j++)
                {
                    distMat[i, j] = 100.0 - distMat[i, j] * 100.0;
                    distMat[j, i] = distMat[i, j];
                }
            }

            return (prof1, prof2);
        }

        // TreeInterface::getWeightsForProfileAlignUPGMA (-CLUSTERING=UPGMA): the profile trees and their weights come FIRST (RootedGuideTree::calcSeqWeights, tree weights only while distanceTree is true; a profile
        // of fewer than two sequences keeps weight 0), from the proper sub-matrix of each profile (UPGMAAlgorithm::generateTree uses getDistMatrix(firstSeq, numSeqs)); then calcPairwiseForProfileAlign
        // (distanceTree = false for good: no Tree::readTree in this branch) and makeSimilarityMatrix
        private (int[] Prof1, int[] Prof2) GetWeightsForProfileAlignUpgma(CwAlignment aln, double[,] distMat)
        {
            int numSeqs = aln.NumSeqs;
            int p1 = aln.Profile1NumSeqs;
            var prof1 = new int[p1];
            var prof2 = new int[numSeqs];
            bool useTree1 = UserTree1 != null && p1 >= 2;
            bool useTree2 = UserTree2 != null && numSeqs - p1 >= 2;
            if (useTree1)
            {
                // the guide-tree file is read with the NJ reader: Tree::readTree, Tree::calcSeqWeights
                var tree = ReadUserProfileTree(UserTree1!, 0, p1);
                Array.Copy(TreeSeqWeights(tree, p1, p1, false), prof1, p1);
            }
            else if (p1 >= 2)
            {
                Array.Copy(UpgmaProfileWeights(distMat, 1, p1, p1, out string text1), prof1, p1);
                GuideTree1Text = text1;
            }

            if (useTree2)
            {
                var tree = ReadUserProfileTree(UserTree2!, p1, numSeqs - p1);
                Array.Copy(TreeSeqWeights(tree, numSeqs - p1, numSeqs, false), 0, prof2, p1, numSeqs - p1);
            }
            else if (numSeqs - p1 >= 2)
            {
                Array.Copy(UpgmaProfileWeights(distMat, p1 + 1, numSeqs - p1, numSeqs, out string text2), 0, prof2, p1, numSeqs - p1);
                GuideTree2Text = text2;
            }

            if (GuideTreesOnly)
            {
                StoppedAfterTrees = true;
                return ([], []);
            }

            _distanceTree = false;
            var group = new int[numSeqs + 1];
            for (int i = 1; i <= p1; i++)
                group[i] = 1;
            for (int i = p1 + 1; i <= numSeqs; i++)
                group[i] = 2;
            var aligned = new int[numSeqs + 1];
            Array.Fill(aligned, 1, 1, numSeqs);
            RunProfileAlign(aln, distMat, group, aligned);
            for (int i = 1; i <= numSeqs; i++)
            {
                for (int j = i + 1; j <= numSeqs; j++)
                {
                    float dscore = aln.CountId(i, j);
                    distMat[i, j] = (100.0 - (double)dscore) / 100.0;
                    distMat[j, i] = distMat[i, j];
                }
            }

            for (int i = 1; i < numSeqs; i++)
            {
                for (int j = i + 1; j <= numSeqs; j++)
                {
                    distMat[i, j] = 100.0 - distMat[i, j] * 100.0;
                    distMat[j, i] = distMat[i, j];
                }
            }

            return (prof1, prof2);
        }

        // RootedGuideTree::calcSeqWeights(firstSeq, lastSeq) of the UPGMA tree of the profile sequences first .. first + count - 1
        private int[] UpgmaProfileWeights(double[,] distMat, int first, int count, int lastSeq, out string text)
        {
            var result = new int[count];
            text = "";
            if (Names != null)
            {
                var block = new double[count, count];
                for (int i = 0; i < count; i++)
                {
                    for (int j = 0; j < count; j++)
                        block[i, j] = distMat[first + i, first + j];
                }

                text = ClustalWUpgmaGuideTreeText(block, Names.Skip(first - 1).Take(count).ToList());
            }

            if (count >= 2 && _distanceTree && !_o.NoWeights)
            {
                var weights = new float[count];
                if (count >= 3)
                {
                    var sub = new double[count, count];
                    for (int i = 0; i < count; i++)
                    {
                        for (int j = 0; j < count; j++)
                            sub[i, j] = distMat[first + i, first + j];
                    }

                    var tree = ClustalWUpgmaTree(sub, count);
                    int Order(UpgmaNode t)
                    {
                        t.Order = t.Left == null && t.Right == null ? 1 : Order(t.Left!) + Order(t.Right!);
                        return t.Order;
                    }

                    void Walk(float soFar, UpgmaNode t)
                    {
                        if (t.Left == null && t.Right == null)
                        {
                            weights[t.SeqNum - 1] = soFar;
                            return;
                        }

                        float w = (float)((double)soFar + (t.Height / t.Order));
                        Walk(w, t.Left!);
                        Walk(w, t.Right!);
                    }

                    Order(tree.Root);
                    Walk(0.0f, tree.Root);
                }

                var weight = new int[count];
                for (int i = 0; i < count; i++)
                    weight[i] = (int)(weights[i] * 100);
                int sum = 0;
                for (int i = 0; i < count; i++)
                    sum += weight[i];
                if (sum == 0)
                {
                    Array.Fill(weight, 1);
                    sum = lastSeq;
                }

                for (int i = 0; i < count; i++)
                {
                    result[i] = weight[i] * CwIntScaleFactor / sum;
                    if (result[i] < 1)
                        result[i] = 1;
                }
            }
            else
            {
                int temp = CwIntScaleFactor / count;
                if (temp < 1)
                    temp = 1;
                Array.Fill(result, temp);
            }

            return result;
        }

        // generateTreeFromDistMatNJ + Tree::readTree for the profile of the sequences first .. first + count - 1 (1-based). Two sequences: the tree file holds the distance of the pair (first, first + 1)
        // halved. More: UnRootedClusterTree::treeFromDistMatrix hands the WHOLE matrix to NJTree::generateTree, which reads the entries (1 .. count, 1 .. count) - the distances of the FIRST `count`
        // sequences of the alignment, not those of the profile (a defect of ClustalW 2.1 that only shows for a profile that does not start at sequence 1) - while the leaves are named after the
        // profile's sequences (name index + first - 1); the NJ works in place and leaves its final state in that block of the matrix.
        private CwRootedTree ProfileTree(double[,] distMat, int first, int count, out string text)
        {
            text = "";
            var slice = Names == null ? null : Names.Skip(first - 1).Take(count).ToList();
            if (count == 2)
            {
                if (slice != null)
                {
                    string d2 = ClustalWPairTreeDistance(distMat[first, first + 1]).Text;
                    text = $"({slice[0]}:{d2},{slice[1]}:{d2});\n";
                }

                float printed = ClustalWPairTreeDistance(distMat[first, first + 1]).ReadBack;
                var root = new CwNode();
                var l0 = new CwNode { Leaf = true, LeafIndex = 0, Parent = root, Dist = printed };
                var l1 = new CwNode { Leaf = true, LeafIndex = 1, Parent = root, Dist = printed };
                root.Left = l0;
                root.Right = l1;
                return new CwRootedTree(root, [l0, l1], [l0, l1]);
            }

            var sub = new double[count, count];
            for (int i = 0; i < count; i++)
            {
                for (int j = 0; j < count; j++)
                    sub[i, j] = distMat[1 + i, 1 + j];
            }

            var finalState = new double[count + 1, count + 1];
            var unrooted = ClustalWNeighborJoining(count, sub, finalState);
            for (int i = 1; i <= count; i++)
            {
                for (int j = 1; j <= count; j++)
                    distMat[i, j] = finalState[i, j];
            }

            if (slice != null)
                text = ClustalWNjTreeText(slice, unrooted, 0.0);
            return ClustalWRootNodes(unrooted, count);
        }

        // Tree::readTree of a user guide tree for the sequences first .. first + count - 1 (0-based) of the alignment
        private CwRootedTree ReadUserProfileTree(string text, int first, int count)
        {
            var user = ClustalWReadTree(text, Names ?? throw new InvalidOperationException("Sequence names are required to read a guide tree."), first, count);
            _distanceTree = user.DistanceTree;
            return new CwRootedTree(user.Root, user.Leaves, user.ParsedLeaves);
        }

        // Tree::calcSeqWeights(firstSeq, lastSeq) for a group of `count` sequences read back from the tree file (null = no tree); `lastSeq` is the 1-based end of the group in the whole alignment
        private int[] TreeSeqWeights(CwRootedTree? tree, int count, int lastSeq, bool duplicateNames)
        {
            var result = new int[count];
            if (count >= 2 && _distanceTree && !_o.NoWeights && tree != null)
            {
                // Tree::orderNodes over all parsed leaves
                foreach (var leaf in tree.ParsedLeaves)
                {
                    for (var p = leaf; p != null; p = p.Parent)
                        p.Order++;
                }

                var weight = new int[count];
                for (int i = 0; i < count; i++)
                {
                    // olptr[i]: the leaf of the sequence's name; when every sequence has the same name each lookup ends at the last matching leaf of the file
                    var start = duplicateNames ? tree.ParsedLeaves[^1] : tree.Leaves[i];
                    float w = 0.0f;
                    for (var p = start; p.Parent != null; p = p.Parent)
                        w += p.Dist / p.Order;
                    w = (float)((double)w * 100.0);
                    weight[i] = (int)w;
                }

                int sum = 0;
                for (int i = 0; i < count; i++)
                    sum += weight[i];
                if (sum == 0)
                {
                    Array.Fill(weight, 1);
                    sum = lastSeq;   // `sum = i` after the loop: the loop variable ends at lastSeq
                }

                for (int i = 0; i < count; i++)
                {
                    result[i] = weight[i] * CwIntScaleFactor / sum;
                    if (result[i] < 1)
                        result[i] = 1;
                }
            }
            else
            {
                int temp = CwIntScaleFactor / count;
                if (temp < 1)
                    temp = 1;
                Array.Fill(result, temp);
            }

            return result;
        }

        // MSA::doProfileAlign
        private void DoProfileAlign(CwAlignment aln, double[,] distMat, int[] prof1Weight, int[] prof2Weight)
        {
            int numSeqs = aln.NumSeqs;
            int p1 = aln.Profile1NumSeqs;
            var seqWeight = new int[numSeqs + 1];
            for (int i = 0; i < p1; i++)
            {
                int maxid = 0;
                for (int j = p1 + 1; j <= numSeqs; j++)
                {
                    if (maxid < distMat[i + 1, j])
                        maxid = (int)distMat[i + 1, j];
                }

                seqWeight[i] = maxid * prof1Weight[i];
            }

            for (int i = p1; i < numSeqs; i++)
            {
                int maxid = -1;
                for (int j = 1; j <= p1; j++)
                {
                    if (maxid < distMat[i + 1, j])
                        maxid = (int)distMat[i + 1, j];
                }

                seqWeight[i] = maxid * prof2Weight[i];
            }

            int sum = 0;
            for (int j = 0; j < numSeqs; j++)
                sum += seqWeight[j];
            if (sum == 0)
            {
                for (int j = 0; j < numSeqs; j++)
                    seqWeight[j] = 1;
                sum = numSeqs;
            }

            for (int j = 0; j < numSeqs; j++)
            {
                seqWeight[j] = seqWeight[j] * CwIntScaleFactor / sum;
                if (seqWeight[j] < 1)
                    seqWeight[j] = 1;
            }
            aln.Weights = seqWeight;   // Alignment::addSeqWeight: the sizes agree (numSeqs + 1)
            var group = new int[numSeqs + 1];
            for (int i = 1; i <= p1; i++)
                group[i] = 1;
            for (int i = p1 + 1; i <= numSeqs; i++)
                group[i] = 2;
            var aligned = new int[numSeqs + 1];
            Array.Fill(aligned, 1, 1, numSeqs);
            RunProfileAlign(aln, distMat, group, aligned);
        }
    }
}
