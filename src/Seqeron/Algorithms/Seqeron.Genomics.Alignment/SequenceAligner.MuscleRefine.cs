namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    // MUSCLE 3.8.1551, iterations 3 and above and the non-"low complexity" pipeline (domuscle.cpp / muscle.cpp DoMuscle):
    //   refinehoriz.cpp  RefineHoriz  - restricted partitions: every edge of the guide tree splits the sequences in two, the two profiles are re-aligned and the
    //                                   new alignment is kept when the objective score improves (scorehistory.cpp stops oscillations)
    //   refinevert.cpp   RefineVert   - the same inside the blocks between the anchor columns (anchors.cpp FindAnchorCols)
    //   objscore.cpp / objscore2.cpp / spfast.cpp  the objective scores SP, DP, XP, PS, SPF (SPM = XP up to 100 sequences, SPF above)
    //   msa2.cpp SetMSAWeightsMuscle, henikoffweight*.cpp, setgscweights.cpp, threewaywt.cpp  the sequence weightings
    //   progressivealign.cpp / refinetree.cpp / realigndiffs.cpp / difftrees.cpp  the profile-less progressive alignment of the pipeline that is used when a weighting
    //                                   other than ClustalW is selected (g_bLow false)
    // All SCORE / WEIGHT / FCOUNT values are single precision, as in the C++.

    internal sealed class MuscleMsa
    {
        public int[] Ids = Array.Empty<int>();
        public char[][] Rows = Array.Empty<char[]>();
        public float[] Weights = Array.Empty<float>();

        public int SeqCount => Rows.Length;

        public int ColCount => Rows.Length == 0 ? 0 : Rows[0].Length;

        public MuscleMsa() { }

        public MuscleMsa(int[] ids, char[][] rows)
        {
            Ids = ids;
            Rows = rows;
            Weights = new float[rows.Length];
        }

        public MuscleMsa Clone() => new((int[])Ids.Clone(), Rows.Select(r => (char[])r.Clone()).ToArray());

        public void CopyFrom(MuscleMsa other)
        {
            Ids = (int[])other.Ids.Clone();
            Rows = other.Rows.Select(r => (char[])r.Clone()).ToArray();
            Weights = new float[Rows.Length];
        }

        public int SeqIndex(int id)
        {
            for (int i = 0; i < Ids.Length; i++)
                if (Ids[i] == id)
                    return i;
            throw new InvalidOperationException($"MSA::GetSeqIndex: id {id} not found");
        }

        public bool IsGapColumn(int col)
        {
            foreach (var r in Rows)
                if (r[col] != '-' && r[col] != '.')
                    return false;
            return true;
        }

        public bool ColumnHasGap(int col)
        {
            foreach (var r in Rows)
                if (r[col] == '-' || r[col] == '.')
                    return true;
            return false;
        }

        public MuscleMsa SubsetByIds(IReadOnlyList<int> ids)
        {
            var rows = new char[ids.Count][];
            var outIds = new int[ids.Count];
            for (int i = 0; i < ids.Count; i++)
            {
                outIds[i] = ids[i];
                rows[i] = (char[])Rows[SeqIndex(ids[i])].Clone();
            }

            return new MuscleMsa(outIds, rows);
        }

        public MuscleMsa SubsetByIndexes(IReadOnlyList<int> idx)
        {
            var rows = new char[idx.Count][];
            var outIds = new int[idx.Count];
            for (int i = 0; i < idx.Count; i++)
            {
                outIds[i] = Ids[idx[i]];
                rows[i] = (char[])Rows[idx[i]].Clone();
            }

            return new MuscleMsa(outIds, rows);
        }

        public MuscleMsa ColumnRange(int from, int count)
        {
            var rows = Rows.Select(r => r.AsSpan(from, count).ToArray()).ToArray();
            return new MuscleMsa((int[])Ids.Clone(), rows);
        }

        public void DeleteGappedColumns()
        {
            var keep = new List<int>();
            for (int c = 0; c < ColCount; c++)
                if (!IsGapColumn(c))
                    keep.Add(c);
            if (keep.Count == ColCount)
                return;
            Rows = Rows.Select(r => keep.Select(c => r[c]).ToArray()).ToArray();
        }
    }

    internal sealed partial class MuscleEngine
    {
        // g_SeqWeight (the method in force; domuscle.cpp switches it from Weight1 to Weight2 before the refinement) and the globals of msa2.cpp
        private MuscleSeqWeight _curWeight;
        private float[] _muscleWeights = Array.Empty<float>();   // g_MuscleWeights: CalcClustalWWeights of the tree set by SetMuscleTree, indexed by leaf INDEX but looked up by sequence ID (as the program does)
        private MuscleTree? _muscleTree;
        private int _splitNode1 = -1, _splitNode2 = -1;          // g_uTreeSplitNode1 / 2

        // parameters set by SetParameters (params.cpp)
        private MuscleTermGaps _termGaps;
        private float _smoothScoreCeil, _minBestColScore, _minSmoothScore;

        private const int LetterGap = 255;

        private int LetterEx(char c) => c == '-' || c == '.' ? LetterGap : _charToLetterEx[c];

        /// <summary>Selects the sequence weighting in force (<c>SetSeqWeightMethod</c>); used by the tests and the MSA modes.</summary>
        internal void UseWeighting(MuscleSeqWeight method) => _curWeight = method;

        /// <summary>Selects the split edge of the three-way weighting (<c>g_uTreeSplitNode1 / 2</c>).</summary>
        internal void UseSplit(int node1, int node2)
        {
            _splitNode1 = node1;
            _splitNode2 = node2;
        }

        internal void SetMuscleTree(MuscleTree tree)
        {
            _muscleTree = tree;
            if (_curWeight != MuscleSeqWeight.ClustalW)
                return;
            _muscleWeights = CalcClustalWWeights(tree);
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // msa2.cpp SetMSAWeightsMuscle

        internal void SetMsaWeights(MuscleMsa m)
        {
            int n = m.SeqCount;
            if (m.Weights.Length != n)
                m.Weights = new float[n];
            switch (_curWeight)
            {
                case MuscleSeqWeight.None:
                    if (n == 0)
                        return;
                    float u = (float)(1.0 / n);
                    for (int i = 0; i < n; i++)
                        m.Weights[i] = u;
                    return;
                case MuscleSeqWeight.Henikoff:
                    SetHenikoffWeights(m);
                    return;
                case MuscleSeqWeight.HenikoffPb:
                    SetHenikoffWeightsPb(m);
                    return;
                case MuscleSeqWeight.Gsc:
                    SetGscWeights(m);
                    return;
                case MuscleSeqWeight.ClustalW:
                    for (int i = 0; i < n; i++)
                    {
                        int id = m.Ids[i];
                        if (id >= _muscleWeights.Length)
                            throw new InvalidOperationException("SetClustalWWeightsMuscle: id out of range");
                        m.Weights[i] = _muscleWeights[id];
                    }

                    NormalizeWeights(m, 1.0f);
                    return;
                case MuscleSeqWeight.ThreeWay:
                    SetThreeWayWeights(m);
                    return;
            }

            throw new InvalidOperationException("SetMSAWeightsMuscle, invalid method");
        }

        private static void NormalizeWeights(MuscleMsa m, float desiredTotal)
        {
            float total = 0;
            for (int i = 0; i < m.SeqCount; i++)
                total += m.Weights[i];
            if (total == 0)
                return;
            float f = desiredTotal / total;
            for (int i = 0; i < m.SeqCount; i++)
                m.Weights[i] *= f;
        }

        // intmath.cpp Normalize(PROB p[], n): PROB is float
        private static void NormalizeVector(float[] p, int n)
        {
            float sum = 0;
            for (int i = 0; i < n; i++)
                sum += p[i];
            if (sum == 0.0f)
                throw new InvalidOperationException("Normalize, sum=0");
            for (int i = 0; i < n; i++)
                p[i] /= sum;
        }

        private static bool IsGapSeq(MuscleMsa m, int s)
        {
            foreach (char c in m.Rows[s])
                if (c != '-' && c != '.')
                    return false;
            return true;
        }

        // henikoffweight.cpp
        private void SetHenikoffWeights(MuscleMsa m)
        {
            int n = m.SeqCount, cols = m.ColCount;
            if (n == 0)
                return;
            if (n == 1)
            {
                m.Weights[0] = 1.0f;
                return;
            }

            if (n == 2)
            {
                m.Weights[0] = 0.5f;
                m.Weights[1] = 0.5f;
                return;
            }

            Array.Clear(m.Weights);
            var count = new uint[20];
            for (int col = 0; col < cols; col++)
            {
                Array.Clear(count);
                uint different = 0;
                for (int s = 0; s < n; s++)
                {
                    int letter = LetterEx(m.Rows[s][col]);
                    if (letter >= 20)
                        continue;
                    uint nc = count[letter] + 1;
                    count[letter] = nc;
                    if (nc == 1)
                        different++;
                }

                for (int s = 0; s < n; s++)
                {
                    int letter = LetterEx(m.Rows[s][col]);
                    if (letter >= 20)
                        continue;
                    uint denom = count[letter] * different;
                    if (denom == 0)
                        continue;
                    m.Weights[s] += (float)(1.0 / denom);
                }
            }

            for (int s = 0; s < n; s++)
                if (IsGapSeq(m, s))
                    m.Weights[s] = 0.0f;
            NormalizeVector(m.Weights, n);
        }

        // henikoffweightpb.cpp
        private void SetHenikoffWeightsPb(MuscleMsa m)
        {
            int n = m.SeqCount, cols = m.ColCount;
            if (n == 0)
                return;
            if (n == 1)
            {
                m.Weights[0] = 1.0f;
                return;
            }

            if (n == 2)
            {
                m.Weights[0] = 0.5f;
                m.Weights[1] = 0.5f;
                return;
            }

            Array.Clear(m.Weights);
            var count = new int[21];
            for (int col = 0; col < cols; col++)
            {
                Array.Clear(count);
                for (int s = 0; s < n; s++)
                    count[HenikoffBucket(m.Rows[s][col])]++;
                bool conserved = false;
                for (int l = 0; l < 21; l++)
                {
                    if (count[l] > 0)
                    {
                        conserved = count[l] == n;
                        break;
                    }
                }

                if (conserved)
                    continue;
                for (int s = 0; s < n; s++)
                    m.Weights[s] += (float)(1.0 / count[HenikoffBucket(m.Rows[s][col])]);
            }

            for (int s = 0; s < n; s++)
                if (IsGapSeq(m, s))
                    m.Weights[s] = 0.0f;
            if (m.Weights.Take(n).All(w => w == 0.0f))
                for (int s = 0; s < n; s++)
                    m.Weights[s] = 1.0f;
            NormalizeVector(m.Weights, n);
        }

        private int HenikoffBucket(char c) => c == '-' || c == '.' || _isWildcard[c] ? 20 : _charToLetterEx[c];

        // ---------------------------------------------------------------------------------------------------------------------------------
        // profile of an MSA (profilefrommsa.cpp ProfileFromMSA: the weights of the method in force, then the counts)

        internal MuscleProfPos[] ProfileFromMsa(MuscleMsa m)
        {
            SetMsaWeights(m);
            var prof = ProfileFromRows(m.Rows, m.Weights);
            for (int c = 0; c < prof.Length; c++)
                prof[c].AllGaps = m.IsGapColumn(c);
            return prof;
        }

        private static List<(char Type, int A, int B)> PathFromMsaPair(MuscleMsa a, MuscleMsa b)
        {
            var path = new List<(char, int, int)>();
            int pa = 0, pb = 0;
            for (int col = 0; col < a.ColCount; col++)
            {
                bool gapA = a.IsGapColumn(col), gapB = b.IsGapColumn(col);
                char type;
                if (!gapA && !gapB)
                {
                    type = 'M';
                    pa++;
                    pb++;
                }
                else if (gapA && !gapB)
                {
                    type = 'I';
                    pb++;
                }
                else if (!gapA)
                {
                    type = 'D';
                    pa++;
                }
                else
                {
                    continue;
                }

                path.Add((type, pa, pb));
            }

            return path;
        }

        private static bool PathEqual(List<(char Type, int A, int B)> x, List<(char Type, int A, int B)> y)
        {
            if (x.Count != y.Count)
                return false;
            for (int i = 0; i < x.Count; i++)
                if (x[i] != y[i])
                    return false;
            return true;
        }

        // aligntwomsas.cpp AlignTwoMSAs
        private (MuscleMsa Merged, List<(char Type, int A, int B)> Path) AlignTwoMsas(MuscleMsa m1, MuscleMsa m2, bool lockLeft, bool lockRight)
        {
            int la = m1.ColCount, lb = m2.ColCount;
            var pa = ProfileFromMsa(m1);
            var pb = ProfileFromMsa(m2);
            if (lockLeft)
            {
                pa[0].GapOpen = MuscleMinusInfinity;
                pb[0].GapOpen = MuscleMinusInfinity;
            }

            if (lockRight)
            {
                pa[la - 1].GapClose = MuscleMinusInfinity;
                pb[lb - 1].GapClose = MuscleMinusInfinity;
            }

            var path = NwSmall(pa, la, pb, lb);
            var rows = RowsFromPath(path, m1.Rows, m2.Rows);
            var merged = new MuscleMsa(m1.Ids.Concat(m2.Ids).ToArray(), rows);
            return (merged, path);
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // objscore2.cpp

        internal float TermGapScore(bool gap)
        {
            switch (_termGaps)
            {
                case MuscleTermGaps.Full:
                    return 0;
                case MuscleTermGaps.Half:
                    return gap ? _gapOpen / 2 : 0;
                default:
                    return gap ? _gapExtend : 0;
            }
        }

        private float ScoreSeqPairLetters(char[] r1, char[] r2)
        {
            int cols = r1.Length;
            if (cols != r2.Length)
                throw new InvalidOperationException("ScoreSeqPairLetters, different lengths");
            int colStart = 0;
            for (int c = 0; c < cols; c++)
            {
                if (!IsGap(r1[c]) || !IsGap(r2[c]))
                {
                    colStart = c;
                    break;
                }
            }

            int colEnd = cols - 1;
            for (int c = cols - 1; c >= 0; c--)
            {
                if (!IsGap(r1[c]) || !IsGap(r2[c]))
                {
                    colEnd = c;
                    break;
                }
            }

            float scoreLetters = 0;
            for (int c = colStart; c <= colEnd; c++)
            {
                int l1 = LetterEx(r1[c]);
                if (l1 >= _alphaSize)
                    continue;
                int l2 = LetterEx(r2[c]);
                if (l2 >= _alphaSize)
                    continue;
                scoreLetters += _matrix[l1, l2];
            }

            return scoreLetters;
        }

        private float ScoreSeqPairGaps(char[] r1, char[] r2)
        {
            int cols = r1.Length;
            if (cols != r2.Length)
                throw new InvalidOperationException("ScoreSeqPairGaps, different lengths");
            float scoreGaps = 0;
            bool gapping1 = false, gapping2 = false;
            int colStart = 0;
            for (int c = 0; c < cols; c++)
            {
                if (!IsGap(r1[c]) || !IsGap(r2[c]))
                {
                    colStart = c;
                    break;
                }
            }

            int colEnd = cols - 1;
            for (int c = cols - 1; c >= 0; c--)
            {
                if (!IsGap(r1[c]) || !IsGap(r2[c]))
                {
                    colEnd = c;
                    break;
                }
            }

            for (int c = colStart; c <= colEnd; c++)
            {
                bool g1 = IsGap(r1[c]), g2 = IsGap(r2[c]);
                if (g1 && g2)
                    continue;
                if (g1)
                {
                    if (!gapping1)
                    {
                        scoreGaps += c == colStart ? TermGapScore(true) : _gapOpen;
                        gapping1 = true;
                    }
                    else
                    {
                        scoreGaps += _gapExtend;
                    }

                    continue;
                }

                if (g2)
                {
                    if (!gapping2)
                    {
                        scoreGaps += c == colStart ? TermGapScore(true) : _gapOpen;
                        gapping2 = true;
                    }
                    else
                    {
                        scoreGaps += _gapExtend;
                    }

                    continue;
                }

                gapping1 = false;
                gapping2 = false;
            }

            if (gapping1 || gapping2)
            {
                scoreGaps -= _gapOpen;
                scoreGaps += TermGapScore(true);
            }

            return scoreGaps;
        }

        internal float ObjScoreSp(MuscleMsa msa)
        {
            int n = msa.SeqCount;
            float total = 0;
            for (int s1 = 0; s1 < n; s1++)
            {
                float w1 = msa.Weights[s1];
                for (int s2 = s1 + 1; s2 < n; s2++)
                {
                    float w2 = msa.Weights[s2];
                    float w = w1 * w2;
                    float letters = ScoreSeqPairLetters(msa.Rows[s1], msa.Rows[s2]);
                    float gaps = ScoreSeqPairGaps(msa.Rows[s1], msa.Rows[s2]);
                    float pair = letters + gaps;
                    total += w * pair;
                }
            }

            return total;
        }

        internal float ObjScoreXp(MuscleMsa m1, MuscleMsa m2)
        {
            if (m1.ColCount != m2.ColCount)
                throw new InvalidOperationException($"ObjScoreXP, alignment lengths differ {m1.ColCount} {m2.ColCount}");
            float total = 0;
            int pairs = 0;
            for (int s1 = 0; s1 < m1.SeqCount; s1++)
            {
                float w1 = m1.Weights[s1];
                for (int s2 = 0; s2 < m2.SeqCount; s2++)
                {
                    float w2 = m2.Weights[s2];
                    float letters = ScoreSeqPairLetters(m1.Rows[s1], m2.Rows[s2]);
                    float gaps = ScoreSeqPairGaps(m1.Rows[s1], m2.Rows[s2]);
                    float pair = letters + gaps;
                    total += w1 * w2 * pair;
                    pairs++;
                }
            }

            if (pairs == 0)
                throw new InvalidOperationException("0 == uPairCount");
            return total;
        }

        internal float ObjScoreDp(MuscleMsa m1, MuscleMsa m2)
        {
            int cols = m1.ColCount;
            if (m2.ColCount != cols)
                throw new InvalidOperationException("ObjScoreDP, must be same length");
            return ObjScoreDpProfiles(ProfileFromMsa(m1), ProfileFromMsa(m2), cols);
        }

        // objscore2.cpp ObjScoreDP_Profs
        internal float ObjScoreDpProfiles(MuscleProfPos[] pa, MuscleProfPos[] pb, int cols)
        {
            float total = 0;
            for (int c = 0; c < cols; c++)
            {
                var a = pa[c];
                var b = pb[c];
                float scoreGap = 0, scoreMatch = 0;
                if (a.AllGaps && b.AllGaps)
                {
                    scoreGap = 0;
                }
                else if (a.AllGaps)
                {
                    if (cols - 1 == c || !pa[c + 1].AllGaps)
                        scoreGap = b.GapClose;
                    if (0 == c || !pa[c - 1].AllGaps)
                        scoreGap += b.GapOpen;
                }
                else if (b.AllGaps)
                {
                    if (cols - 1 == c || !pb[c + 1].AllGaps)
                        scoreGap = a.GapClose;
                    if (0 == c || !pb[c - 1].AllGaps)
                        scoreGap += a.GapOpen;
                }
                else
                {
                    scoreMatch = ScoreProfPos2(a, b);
                }

                total += scoreMatch + scoreGap;
            }

            return total;
        }

        internal float ObjScorePs(MuscleMsa msa)
        {
            if (_pp != MuscleScoreModel.LogExpectation)
                throw new InvalidOperationException("FastScoreMSA_LASimple: LA");
            int n = msa.SeqCount, cols = msa.ColCount;
            var prof = ProfileFromMsa(msa);
            float total = 0;
            for (int s = 0; s < n; s++)
            {
                float weightSeq = msa.Weights[s];
                float scoreSeq = 0;
                var row = msa.Rows[s];
                for (int c = 0; c < cols; c++)
                {
                    var pp = prof[c];
                    if (IsGap(row[c]))
                    {
                        bool open = 0 == c || !IsGap(row[c - 1]);
                        bool close = cols - 1 == c || !IsGap(row[c + 1]);
                        if (open)
                            scoreSeq += pp.GapOpen;
                        if (close)
                            scoreSeq += pp.GapClose;
                    }
                    else if (!_isWildcard[row[c]])
                    {
                        int letter = _charToLetterEx[row[c]];
                        scoreSeq += pp.AaScores[letter];
                    }
                }

                total += weightSeq * scoreSeq;
            }

            return total;
        }

        // spfast.cpp ObjScoreSPDimer
        private float[,] GapScoreMatrix()
        {
            float t = 0.2f;
            const int LL = 0, LG = 1, GL = 2, GG = 3;
            var m = new float[4, 4];
            m[LL, LG] = _gapOpen;
            m[LG, LL] = _gapOpen;
            m[LG, GL] = _gapOpen;
            m[LG, GG] = t * _gapOpen;
            m[GL, LG] = _gapOpen;
            m[GG, LG] = t * _gapOpen;
            return m;
        }

        private static float SpFreqs(float[] freqs, int n, Func<int, int, float> mx, Func<int, float> diagRow)
        {
            float totalOffDiag = 0, totalDiag = 0;
            for (int i = 0; i < n; i++)
            {
                float fi = freqs[i];
                if (0 == fi)
                    continue;
                float diagt = fi * fi * diagRow(i);
                totalDiag += diagt;
                float sum = 0;
                for (int j = 0; j < i; j++)
                {
                    float t = freqs[j] * mx(i, j);
                    sum += t;
                }

                totalOffDiag += fi * sum;
            }

            return totalOffDiag * 2 + totalDiag;
        }

        internal float ObjScoreSpDimer(MuscleMsa msa)
        {
            var gsm = GapScoreMatrix();
            float total = 0;
            int n = msa.SeqCount, cols = msa.ColCount;
            for (int c = 0; c < cols; c++)
            {
                var freqs = new float[20];
                var gapFreqs = new float[4];
                float selfOver = 0, gapSelfOver = 0;
                for (int s = 0; s < n; s++)
                {
                    float w = msa.Weights[s];
                    bool gapThis = IsGap(msa.Rows[s][c]);
                    bool gapPrev = c != 0 && IsGap(msa.Rows[s][c - 1]);
                    int gapType = (gapThis ? 1 : 0) + 2 * (gapPrev ? 1 : 0);
                    gapFreqs[gapType] += w;
                    float gapt = w * w * gsm[gapType, gapType];
                    gapSelfOver += gapt;
                    if (gapThis)
                        continue;
                    int letter = LetterEx(msa.Rows[s][c]);
                    if (letter >= 20)
                        continue;
                    freqs[letter] += w;
                    float t = w * w * _matrix[letter, letter];
                    selfOver += t;
                }

                float spf = SpFreqs(freqs, 20, (i, j) => _matrix[i, j], i => _matrix[i, i]);
                float col = spf - selfOver;
                float spfGaps = SpFreqs(gapFreqs, 4, (i, j) => gsm[i, j], i => gsm[i, i]);
                float colGaps = spfGaps - gapSelfOver;
                total += col + colGaps;
            }

            return total / 2;
        }

        // objscore.cpp ObjScore / ObjScoreIds
        private float ObjScore(MuscleMsa msa, int[] idx1, int[] idx2)
        {
            int seqCount = msa.SeqCount;
            var os = _o.ObjectiveScore;
            if (os == MuscleObjectiveScore.Spm)
                os = seqCount <= 100 ? MuscleObjectiveScore.Xp : MuscleObjectiveScore.Spf;

            MuscleMsa? m1 = null, m2 = null;
            switch (os)
            {
                case MuscleObjectiveScore.Dp:
                case MuscleObjectiveScore.Xp:
                    m1 = msa.SubsetByIndexes(idx1);
                    m2 = msa.SubsetByIndexes(idx2);
                    SetMsaWeights(m1);
                    SetMsaWeights(m2);
                    break;
                default:
                    SetMsaWeights(msa);
                    break;
            }

            return os switch
            {
                MuscleObjectiveScore.Sp => ObjScoreSp(msa),
                MuscleObjectiveScore.Dp => ObjScoreDp(m1!, m2!),
                MuscleObjectiveScore.Xp => ObjScoreXp(m1!, m2!),
                MuscleObjectiveScore.Ps => ObjScorePs(msa),
                MuscleObjectiveScore.Spf => ObjScoreSpDimer(msa),
                _ => throw new InvalidOperationException("Invalid g_ObjScore"),
            };
        }

        private float ObjScoreIds(MuscleMsa msa, int[] ids1, int[] ids2) =>
            ObjScore(msa, ids1.Select(msa.SeqIndex).ToArray(), ids2.Select(msa.SeqIndex).ToArray());

        // ---------------------------------------------------------------------------------------------------------------------------------
        // diffobjscore.cpp GetLetterScores / ScoreColLetters, anchors.cpp FindAnchorCols

        internal float ScoreColLetters(MuscleMsa msa, int col)
        {
            int n = msa.SeqCount;
            double total = 0;
            for (int s = 0; s < n; s++)
                total += msa.Weights[s];
            if (total <= 0)
                return 0;

            var freqs = new float[20];
            float score = 0;
            for (int s = 0; s < n; s++)
            {
                int letter = LetterEx(msa.Rows[s][col]);
                if (letter >= _alphaSize)
                    continue;
                float w = msa.Weights[s];
                freqs[letter] += w;
                score -= w * w * _matrix[letter, letter];
            }

            for (int l1 = 0; l1 < _alphaSize; l1++)
            {
                float f1 = freqs[l1];
                score += f1 * f1 * _matrix[l1, l1];
                for (int l2 = l1 + 1; l2 < _alphaSize; l2++)
                {
                    float f2 = freqs[l2];
                    score += 2 * f1 * f2 * _matrix[l1, l2];
                }
            }

            score /= 2;
            return score;
        }

        private static float[] WindowSmooth(float[] score, int count, int window, double ceil)
        {
            if (window % 2 != 1)
                throw new InvalidOperationException($"WindowSmooth={window} must be odd");
            var smooth = new float[count];
            float Ceil(float x) => (float)(x > ceil ? ceil : x);
            if (count <= window)
                return smooth;
            int w2 = window / 2;
            float windowTotal = 0;
            for (int i = 0; i < window; i++)
                windowTotal += Ceil(score[i]);
            for (int i = w2; i <= count - w2 - 1; i++)
            {
                smooth[i] = windowTotal / window;
                if (i == count - w2 - 1)
                    break;
                windowTotal -= Ceil(score[i - w2]);
                windowTotal += Ceil(score[i + w2 + 1]);
            }

            return smooth;
        }

        internal List<int> FindAnchorCols(MuscleMsa msa)
        {
            int cols = msa.ColCount;
            var anchors = new List<int>();
            if (cols < 16)
                return anchors;

            var match = new float[cols];
            for (int c = 0; c < cols; c++)
                match[c] = ScoreColLetters(msa, c);
            var smooth = WindowSmooth(match, cols, _o.SmoothWindow, _smoothScoreCeil);

            var best = new List<int>();
            for (int c = 0; c < cols; c++)
            {
                if (match[c] < (double)_minBestColScore)
                    continue;
                if (smooth[c] < (double)_minSmoothScore)
                    continue;
                if (msa.ColumnHasGap(c))
                    continue;
                best.Add(c);
            }

            // MergeBestCols
            uint window = (uint)_o.AnchorSpacing;
            for (int nn = 0; nn < best.Count;)
            {
                uint bestIdx = (uint)best[nn];
                int within = 0;
                for (int i = nn + 1; i < best.Count; i++)
                {
                    uint b2 = (uint)best[i];
                    if (unchecked(b2 - bestIdx) >= window)
                        break;
                    within++;
                }

                uint anchor = bestIdx;
                if (within == 1)
                {
                    uint b2 = (uint)best[nn + 1];
                    anchor = match[bestIdx] > match[b2] ? bestIdx : b2;
                }
                else if (within > 1)
                {
                    int closestDist = (int)window;
                    uint closest = bestIdx;
                    for (int i = nn + 1; i < nn + within; i++)
                    {
                        uint colIdx = (uint)best[i];
                        int dist = unchecked((int)(colIdx - bestIdx));
                        if (dist < 0)
                            dist = -dist;
                        if (dist < closestDist)
                        {
                            closest = colIdx;
                            closestDist = dist;
                        }
                    }

                    anchor = closest;
                }

                anchors.Add((int)anchor);
                nn += within + 1;
            }

            return anchors;
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // tree helpers (phy2.cpp / phy4.cpp)

        private static void GetLeaves(MuscleTree t, int node, List<int> leaves)
        {
            if (t.IsLeaf(node))
            {
                leaves.Add(node);
                return;
            }

            GetLeaves(t, t.Left[node], leaves);
            GetLeaves(t, t.Right[node], leaves);
        }

        private static void GetLeavesExcluding(MuscleTree t, int node, int exclude, List<int> leaves)
        {
            if (node == exclude)
                return;
            if (t.IsLeaf(node))
            {
                leaves.Add(node);
                return;
            }

            GetLeavesExcluding(t, t.Left[node], exclude, leaves);
            GetLeavesExcluding(t, t.Right[node], exclude, leaves);
        }

        // Tree::GetNodeHeight (edge lengths clamped at 0, average of the two children)
        private static double NodeHeight(MuscleTree t, int node, double[] cache, bool[] has)
        {
            if (t.IsLeaf(node))
                return 0.0;
            if (has[node])
                return cache[node];
            int left = t.Left[node], right = t.Right[node];
            double dl = t.EdgeLength[left], dr = t.EdgeLength[right];
            if (double.IsNaN(dl) || double.IsNaN(dr))
                throw new InvalidOperationException("Missing edge length in tree");
            if (dl < 0)
                dl = 0;
            if (dr < 0)
                dr = 0;
            double hl = dl + NodeHeight(t, left, cache, has);
            double hr = dr + NodeHeight(t, right, cache, has);
            double h = (hl + hr) / 2;
            has[node] = true;
            cache[node] = h;
            return h;
        }

        private static int[] InternalNodesInHeightOrder(MuscleTree t)
        {
            int nodeCount = t.NodeCount;
            int internalCount = (nodeCount - 1) / 2;
            var cache = new double[nodeCount];
            var has = new bool[nodeCount];
            var idx = new int[internalCount];
            var heights = new double[internalCount];
            int k = 0;
            for (int n = 0; n < nodeCount; n++)
            {
                if (t.IsLeaf(n))
                    continue;
                idx[k] = n;
                heights[k] = NodeHeight(t, n, cache, has);
                k++;
            }

            bool done = false;
            while (!done)
            {
                done = true;
                for (int i = 0; i < internalCount - 1; i++)
                {
                    if (heights[i] > heights[i + 1])
                    {
                        (heights[i], heights[i + 1]) = (heights[i + 1], heights[i]);
                        (idx[i], idx[i + 1]) = (idx[i + 1], idx[i]);
                        done = false;
                    }
                }
            }

            return idx;
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // scorehistory.cpp

        private sealed class ScoreHistory
        {
            private readonly float[][] _score;
            private readonly bool[][] _set;
            private readonly int _iters;

            public ScoreHistory(int iters, int nodeCount)
            {
                _iters = iters;
                _score = new float[iters][];
                _set = new bool[iters][];
                for (int n = 0; n < iters; n++)
                {
                    _score[n] = new float[nodeCount * 2];
                    _set[n] = new bool[nodeCount * 2];
                }
            }

            // true = repeated (oscillating)
            public bool SetScore(int iter, int node, bool right, float score)
            {
                if (iter >= _iters)
                    throw new InvalidOperationException("ScoreHistory::SetScore-1");
                int index = node * 2 + (right ? 1 : 0);
                for (int n = 1; n < iter; n++)
                {
                    int prev = n - 1;
                    if (!_set[prev][index])
                        throw new InvalidOperationException("ScoreHistory::SetScore-3");
                    if (_score[prev][index] == score)
                        return true;
                }

                _score[iter][index] = score;
                _set[iter][index] = true;
                return false;
            }
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // refinehoriz.cpp

        // The scores are in/out: when a subset has no residue column left the program returns without setting them (uninitialised locals of RefineHeightParts that keep the value of the
        // previous edge: they live in the caller's frame and are passed by address).
        private bool TryRealign(MuscleMsa msaIn, MuscleTree tree, List<int> leaves1, List<int> leaves2, ref float scoreBefore, ref float scoreAfter, bool lockLeft, bool lockRight)
        {
            var ids1 = leaves1.Select(l => tree.LeafId[l]).ToArray();
            var ids2 = leaves2.Select(l => tree.LeafId[l]).ToArray();
            var msa1 = msaIn.SubsetByIds(ids1);
            var msa2 = msaIn.SubsetByIds(ids2);

            var pathBefore = PathFromMsaPair(msa1, msa2);
            msa1.DeleteGappedColumns();
            msa2.DeleteGappedColumns();
            if (msa1.ColCount == 0 || msa2.ColCount == 0)
                return false;

            var (realigned, pathAfter) = AlignTwoMsas(msa1, msa2, lockLeft, lockRight);
            bool anyChanges = !PathEqual(pathAfter, pathBefore);
            if (!anyChanges)
            {
                scoreBefore = 0;
                scoreAfter = 0;
                return false;
            }

            SetMsaWeights(msaIn);
            SetMsaWeights(realigned);
            scoreBefore = ObjScoreIds(msaIn, ids1, ids2);
            scoreAfter = ObjScoreIds(realigned, ids1, ids2);
            bool accept = scoreAfter > scoreBefore;
            if (accept)
                msaIn.CopyFrom(realigned);
            return accept;
        }

        private bool RefineHeightParts(MuscleMsa msa, MuscleTree tree, int[] internals, int iter, ScoreHistory history, out bool oscillating, bool lockLeft, bool lockRight, bool right)
        {
            oscillating = false;
            int seqCount = msa.SeqCount;
            int internalCount = seqCount - 1;
            bool anyAccepted = false;
            float before = 0, after = 0;
            for (int i = 0; i < internalCount; i++)
            {
                AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per realigned edge (each is a profile-profile dynamic program)
                int node = internals[i];
                int neighbor;
                if (node == tree.Root && !right)
                    continue;
                else if (right)
                    neighbor = tree.Right[node];
                else
                    neighbor = tree.Left[node];

                _splitNode1 = node;
                _splitNode2 = neighbor;

                var leaves1 = new List<int>();
                var leaves2 = new List<int>();
                GetLeaves(tree, neighbor, leaves1);
                GetLeavesExcluding(tree, tree.Root, neighbor, leaves2);

                bool accepted = TryRealign(msa, tree, leaves1, leaves2, ref before, ref after, lockLeft, lockRight);
                float scoreMax = after > before ? after : before;
                bool repeated = history.SetScore(iter, node, right, scoreMax);
                if (repeated)
                {
                    oscillating = true;
                    break;
                }

                if (accepted)
                    anyAccepted = true;
            }

            return anyAccepted;
        }

        internal bool RefineHoriz(MuscleMsa msa, MuscleTree tree, long iters, bool lockLeft, bool lockRight)
        {
            int seqCount = msa.SeqCount;
            if (seqCount < 3)
                return false;
            int internalCount = seqCount - 1;
            var forward = InternalNodesInHeightOrder(tree);
            var reverse = new int[internalCount];
            for (int n = 0; n < internalCount; n++)
                reverse[internalCount - 1 - n] = forward[n];

            int iterCount = (int)Math.Min(iters, int.MaxValue / 2);
            var history = new ScoreHistory(iterCount, 2 * seqCount - 1);
            bool anyAnyIter = false;
            for (int iter = 0; iter < iterCount; iter++)
            {
                AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per refinement iteration
                bool anyThisIter = false;
                bool reversed = iter % 2 != 0;
                var internals = reversed ? reverse : forward;
                bool oscillating = false;
                for (int i = 0; i < 2; i++)
                {
                    bool right = i == 0;
                    bool anyChanges = RefineHeightParts(msa, tree, internals, iter, history, out oscillating, lockLeft, lockRight, right);
                    if (oscillating)
                        return anyAnyIter;
                    if (anyChanges)
                    {
                        anyThisIter = true;
                        anyAnyIter = true;
                    }
                }

                if (!anyThisIter)
                    break;
            }

            return anyAnyIter;
        }

        // refinevert.cpp
        internal bool RefineVert(MuscleMsa msa, MuscleTree tree, long iters)
        {
            int colCountIn = msa.ColCount, seqCountIn = msa.SeqCount;
            if (colCountIn < 3 || seqCountIn < 3)
                return false;

            SetMsaWeights(msa);
            var anchorCols = FindAnchorCols(msa);
            int rangeCount = anchorCols.Count + 1;

            var outRows = new List<char>[seqCountIn];
            for (int s = 0; s < seqCountIn; s++)
                outRows[s] = new List<char>();
            bool anyChanges = false;
            for (int r = 0; r < rangeCount; r++)
            {
                int left = r > 0 ? anchorCols[r - 1] : 0;
                int right = r < anchorCols.Count ? anchorCols[r] : colCountIn;
                int count = right - left;
                if (count == 0)
                    continue;
                MuscleMsa range;
                if (count == 1)
                {
                    range = msa.ColumnRange(left, 1);
                    Append(outRows, msa, range);
                    continue;
                }

                range = msa.ColumnRange(left, count);
                bool lockLeft = r != 0, lockRight = rangeCount - 1 != r;
                bool anyThisBlock = RefineHoriz(range, tree, iters, lockLeft, lockRight);
                anyChanges = anyChanges || anyThisBlock;
                Append(outRows, msa, range);
            }

            if (anyChanges)
                msa.Rows = outRows.Select(l => l.ToArray()).ToArray();
            return anyChanges;
        }

        // MSAAppend: the sequences of the block are matched by id
        private static void Append(List<char>[] outRows, MuscleMsa reference, MuscleMsa block)
        {
            for (int s = 0; s < reference.SeqCount; s++)
            {
                int s2 = block.SeqIndex(reference.Ids[s]);
                outRows[s].AddRange(block.Rows[s2]);
            }
        }
    }
}
