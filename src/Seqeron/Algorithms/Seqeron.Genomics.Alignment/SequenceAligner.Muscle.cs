using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MUSCLE 3.8.1551 progressive alignment (maxiters 1 and 2)

    // Step-by-step port of MUSCLE 3.8.1551 (Edgar 2004, NAR 32:1792; muscle3.8.1551_src, Debian 3.8.1551-2 source tree) for the command line
    //   muscle -in x.fa -out y.fa -maxiters 1|2 [-seqtype ...] [-le|-sp|-sv|-spn] [-gapopen -gapextend -center] [-distance1 -cluster1 -cluster2 -sueff -maxtrees]
    // i.e. domuscle.cpp::DoMuscle with the default g_bLow pipeline:
    //   1. seq.cpp / fasta2.cpp / alpha.cpp: letters upper-cased, gaps and white space dropped, letters outside the alphabet replaced by the wildcard (X protein, N nucleotide);
    //      mhack.cpp: a leading M of a protein is treated as X during the whole alignment and restored at the end;
    //   2. fastdistmafft.cpp / fastdistnuc.cpp: k-mer distance (DistKmer6_6 protein / DistKmer4_6 nucleotide), fastclust.cpp + upgma2.cpp: UPGMB guide tree (pseudo root);
    //   3. progressivealign (ProgressiveAlignE): clwwt.cpp ClustalW tree weights, profilefrommsa.cpp leaf profiles, aligntwoprofs.cpp: nwsmall.cpp (affine DP with position-specific
    //      gap open / close, scorepp.cpp LE / SP / SPN score, termgaps.cpp half terminal gaps) + aligngivenpath.cpp (profile of the merged node from the frequencies), makerootmsa.cpp;
    //   4. maxiters 2: refinetreee.cpp: distance = Kimura of the draft MSA (msadistkimura.cpp), UPGMB, difftreese.cpp (nodes whose leaf set is unchanged are kept), realigndiffse.cpp.
    // All arithmetic that is single precision in the C++ (SCORE / FCOUNT / WEIGHT = float, dist_t = float) is single precision here; the double expressions of the C++ stay double.
    // Not ported: -diags / -diags1 / -diags2 and -anchors (finddiags.cpp, anchors.cpp), the restricted-partition refinement of iterations >= 3 (refinevert.cpp / refinehoriz.cpp with the
    // objective score of objscore2.cpp), -weight / -distance2 / -root variants, -usetree, -matrix (user matrices).

    private const float MuscleMinusInfinity = -1e37f;
    private const byte MuscleBitDm = 0x01, MuscleBitIm = 0x02, MuscleBitXm = 0x03, MuscleBitMd = 0x04, MuscleBitXd = 0x04, MuscleBitMi = 0x08, MuscleBitXi = 0x08;
    private const uint MuscleNodeChanged = uint.MaxValue;

    /// <summary>
    /// MUSCLE 3.8.1551's multiple alignment of nucleotide sequences (<c>muscle -seqtype dna -maxiters 1|2</c>), see <see cref="MultipleAlignMuscle(IEnumerable{string}, MuscleOptions?, MatrixScoring?)"/>.
    /// </summary>
    public static MultipleAlignmentResult MultipleAlignMuscle(IEnumerable<DnaSequence> sequences, MuscleOptions? options = null, MatrixScoring? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var list = sequences.ToList();
        if (list.Count == 0)
            return MultipleAlignmentResult.Empty;
        var o = (options ?? new MuscleOptions()) with { SequenceType = options?.SequenceType is MuscleSequenceType.Rna ? MuscleSequenceType.Rna : MuscleSequenceType.Dna };
        return MultipleAlignMuscle(list.Select(s => s.Sequence), o, scoring);
    }

    /// <summary>MUSCLE 3.8.1551's multiple alignment of protein sequences (<c>muscle -seqtype protein -maxiters 1|2</c>), see <see cref="MultipleAlignMuscle(IEnumerable{string}, MuscleOptions?, MatrixScoring?)"/>.</summary>
    public static MultipleAlignmentResult MultipleAlignMuscle(IEnumerable<ProteinSequence> sequences, MuscleOptions? options = null, MatrixScoring? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var list = sequences.ToList();
        if (list.Count == 0)
            return MultipleAlignmentResult.Empty;
        return MultipleAlignMuscle(list.Select(s => s.Sequence), (options ?? new MuscleOptions()) with { SequenceType = MuscleSequenceType.Protein }, scoring);
    }

    /// <summary>
    /// MUSCLE 3.8.1551's multiple alignment (Edgar 2004, NAR 32:1792), a step-by-step port of the program's default pipeline for <c>-maxiters 1</c> (the progressive draft: k-mer distance,
    /// UPGMB guide tree, ClustalW weights, profile-profile alignment with the log-expectation / sum-of-pairs score and position-specific gap penalties) and <c>-maxiters 2</c> (plus the
    /// tree-dependent refinement: rebuild the tree from the Kimura distances of the draft, realign only the changed subtrees).
    /// </summary>
    /// <remarks>
    /// Rows are returned in input order (<see cref="MuscleOptions.TreeOrder"/> = the order <c>muscle</c> writes). Letters are upper-cased, gaps ('-', '.') and white space are dropped,
    /// letters outside the alphabet are replaced by the wildcard (protein X, nucleotide N; MUSCLE warns "invalid letters found"), and MUSCLE's methionine hack applies (a leading M of a protein is aligned as X and
    /// restored in the output). The reported <see cref="MultipleAlignmentResult.TotalScore"/> is the natural affine sum of pairs under <paramref name="scoring"/> (default: BLOSUM62, open 10 / extend 0.2 for
    /// proteins; EDNAFULL, 10 / 0.5 for nucleotides), not MUSCLE's objective. Parameters not ported are listed in the region comment of the source file.
    /// The oracle of this port is the MUSCLE 3.8.1551 binary built from the same source (docs/Validation/review-2026-09/handoff/A1-10c-muscle); the installed <c>muscle</c> 5.x is a different program.
    /// </remarks>
    /// <param name="sequences">Sequences (letters; at least one, none empty).</param>
    /// <param name="options">MUSCLE parameters (default: the defaults of MUSCLE 3.8, <see cref="MuscleOptions.MaxIterations"/> 2).</param>
    /// <param name="scoring">Scoring of the reported total score only.</param>
    /// <exception cref="ArgumentException">An empty sequence, a non-letter character or an unusable parameter.</exception>
    /// <exception cref="NotSupportedException"><see cref="MuscleOptions.MaxIterations"/> above 2.</exception>
    public static MultipleAlignmentResult MultipleAlignMuscle(IEnumerable<string> sequences, MuscleOptions? options = null, MatrixScoring? scoring = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var o = options ?? new MuscleOptions();
        var list = sequences.ToList();
        if (list.Count == 0)
            return MultipleAlignmentResult.Empty;

        var engine = new MuscleEngine(o, list);
        var (rows, nucleotide) = engine.Run();
        var aligned = rows.ToList();
        string consensus = BuildConsensus(aligned, aligned.Max(r => r.Length), anySymbol: true);
        var sp = scoring ?? (nucleotide ? MatrixScoring.Emboss(SubstitutionMatrix.EdnaFull, 10, 0.5) : MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.2));
        string[] scored = nucleotide ? rows.Select(r => r.Replace('X', 'N')).ToArray() : rows;
        double exact = rows.Length == 1 ? 0.0 : MatrixSumOfPairsScore(scored, sp, MsaGapModel.Affine);
        return new MultipleAlignmentResult(rows, consensus, (int)Math.Round(exact, MidpointRounding.AwayFromZero)) { ExactTotalScore = exact };
    }

    /// <summary>The same alignment together with the intermediate artefacts the oracle checks (draft, guide tree order); internal, used by the tests.</summary>
    internal static (string[] Rows, string[] DraftRows, int[] TreeOrder, string Newick) MultipleAlignMuscleDetailed(IReadOnlyList<string> sequences, MuscleOptions options)
    {
        var engine = new MuscleEngine(options, sequences);
        engine.Run();
        return (engine.FinalRowsInputOrder, engine.DraftRowsInputOrder, engine.FinalTreeOrder, engine.GuideTreeNewick);
    }

    internal sealed class MuscleProfPos
    {
        public readonly int[] SortOrder = new int[21];
        public readonly float[] Counts = new float[20];
        public readonly float[] AaScores = new float[20];
        public float LL, LG, GL, GG, Occ, GapOpen, GapClose;
    }

    internal sealed class MuscleNode
    {
        public int[] Ids = Array.Empty<int>();       // sequence ids of the rows, in the node's row order
        public char[][] Rows = Array.Empty<char[]>();
        public MuscleProfPos[] Prof = Array.Empty<MuscleProfPos>();
        public int Length;
        public float Weight;
    }

    internal sealed class MuscleTree
    {
        public int LeafCount;
        public int Root;
        public int[] Left = Array.Empty<int>(), Right = Array.Empty<int>(), Parent = Array.Empty<int>();
        public double[] EdgeLength = Array.Empty<double>();   // length of the edge to the parent (stored from a float)
        public int[] LeafId = Array.Empty<int>();
        public int NodeCount => 2 * LeafCount - 1;
        public bool IsLeaf(int node) => node < LeafCount;

        /// <summary>Depth-first (left first) order of the nodes, children before the parent (<c>FirstDepthFirstNode</c> / <c>NextDepthFirstNode</c>).</summary>
        public List<int> PostOrder()
        {
            var order = new List<int>(NodeCount);
            var stack = new Stack<(int Node, bool Visited)>();
            stack.Push((Root, false));
            while (stack.Count > 0)
            {
                var (n, visited) = stack.Pop();
                if (IsLeaf(n) || visited)
                {
                    order.Add(n);
                    continue;
                }

                stack.Push((n, true));
                stack.Push((Right[n], false));
                stack.Push((Left[n], false));
            }

            return order;
        }
    }

    internal sealed class MuscleEngine
    {
        private readonly MuscleOptions _o;
        private readonly int _n;
        private readonly bool _nuc;
        private readonly bool _rna;
        private readonly int _alphaSize;
        private readonly int _wildcardLetterBase;
        private readonly string _wildcardChars;
        private readonly int[] _charToLetterEx = new int[256];
        private readonly bool[] _isResidue = new bool[256];
        private readonly bool[] _isWildcard = new bool[256];
        private readonly char[][] _seq;      // letters (after FixAlpha and the methionine hack)
        private readonly bool[] _mHack;

        // parameters of the profile-profile score (params.cpp SetDefaultsLE / SP / SV / SPN)
        private MuscleScoreModel _pp;
        private float _gapOpen, _gapExtend, _center;
        private bool _normalizeCounts;
        private float[,] _matrix = new float[32, 32];
        private MuscleDistance _distance1;

        public string[] FinalRowsInputOrder { get; private set; } = Array.Empty<string>();
        public string[] DraftRowsInputOrder { get; private set; } = Array.Empty<string>();
        public int[] FinalTreeOrder { get; private set; } = Array.Empty<int>();
        public string GuideTreeNewick { get; private set; } = "";

        public MuscleEngine(MuscleOptions options, IReadOnlyList<string> sequences)
        {
            _o = options;
            _n = sequences.Count;
            if (options.MaxIterations < 1)
                throw new ArgumentException("MaxIterations must be at least 1.", nameof(options));
            if (options.MaxIterations > 2)
                throw new NotSupportedException("MUSCLE's iterations 3 and above (restricted-partition refinement, refinevert.cpp / refinehoriz.cpp) are not ported; use MaxIterations 1 or 2.");
            if (options.MaxTrees < 0)
                throw new ArgumentException("MaxTrees must be non-negative.", nameof(options));
            if (!Enum.IsDefined(options.SequenceType) || !Enum.IsDefined(options.ScoreModel) || !Enum.IsDefined(options.Distance1) || !Enum.IsDefined(options.Cluster1) || !Enum.IsDefined(options.Cluster2))
                throw new ArgumentOutOfRangeException(nameof(options), "Unknown MUSCLE option value.");

            var cleaned = new string[_n];
            for (int i = 0; i < _n; i++)
            {
                ArgumentNullException.ThrowIfNull(sequences[i]);
                cleaned[i] = CleanSequence(sequences[i]);
            }


            for (int i = 0; i < cleaned.Length; i++)
            {
                if (cleaned[i].Length == 0)
                    throw new ArgumentException($"Sequence {i} is empty (MUSCLE ignores empty FASTA records; pass only non-empty sequences).", nameof(sequences));
            }

            var type = options.SequenceType;
            if (type == MuscleSequenceType.Auto)
                type = GuessAlphabet(cleaned);
            _nuc = type != MuscleSequenceType.Protein;
            _rna = type == MuscleSequenceType.Rna;
            _alphaSize = _nuc ? 4 : 20;
            _wildcardChars = _nuc ? "MRWSYKVHDBXN" : "BXZ";
            _wildcardLetterBase = _nuc ? 4 : 20;
            InitAlphabet();

            _seq = new char[_n][];
            for (int i = 0; i < _n; i++)
            {
                var s = cleaned[i].ToCharArray();
                for (int k = 0; k < s.Length; k++)
                {
                    if (!_isResidue[s[k]])
                        s[k] = _nuc ? 'N' : 'X';
                }

                _seq[i] = s;
            }

            // MHackStart: the leading M of a protein (when there are at least two sequences) is aligned as X
            _mHack = new bool[_n];
            if (!_nuc && _n > 1)
            {
                for (int i = 0; i < _n; i++)
                {
                    if (_seq[i][0] == 'M')
                    {
                        _mHack[i] = true;
                        _seq[i][0] = 'X';
                    }
                }
            }

            SetParameters();
        }

        private static string CleanSequence(string s)
        {
            var sb = new System.Text.StringBuilder(s.Length);
            foreach (char c in s)
            {
                if (char.IsWhiteSpace(c) || c == '-' || c == '.')
                    continue;
                if (c is >= 'A' and <= 'Z' or >= 'a' and <= 'z')
                    sb.Append(char.ToUpperInvariant(c));
                else
                    throw new ArgumentException($"Invalid character '{c}' in a sequence (MUSCLE ignores it with a warning; only letters, '-', '.' and white space are accepted here).");
            }

            return sb.ToString();
        }

        // SeqVect::GuessAlpha
        private static MuscleSequenceType GuessAlphabet(string[] seqs)
        {
            const int charCount = 100, minNucleoPct = 95;
            int dna = 0, rna = 0, total = 0;
            foreach (char c in string.Concat(seqs).Take(charCount))
            {
                if ("AGCTN".Contains(c))
                    dna++;
                if ("AGCUN".Contains(c))
                    rna++;
                total++;
            }

            if (total != 0 && dna * 100 / total >= minNucleoPct)
                return MuscleSequenceType.Dna;
            if (total != 0 && rna * 100 / total >= minNucleoPct)
                return MuscleSequenceType.Rna;
            return MuscleSequenceType.Protein;
        }

        private void InitAlphabet()
        {
            Array.Fill(_charToLetterEx, -1);
            string residues = _nuc ? "ACGT" : "ACDEFGHIKLMNPQRSTVWY";
            for (int i = 0; i < residues.Length; i++)
                SetLetter(residues[i], i, residue: true);
            if (_nuc)
                SetLetter('U', 3, residue: true);   // SetAlpha(DNA) falls through into SetAlphaRNA: U = T for both
            for (int i = 0; i < _wildcardChars.Length; i++)
                SetLetter(_wildcardChars[i], i + _wildcardLetterBase, residue: false);
            // protein wildcards in MUSCLE's enum order: X = 20, B = 21, Z = 22; nucleotide M R W S Y K V H D B X N = 4 .. 15
            if (!_nuc)
            {
                SetLetter('X', 20, residue: false);
                SetLetter('B', 21, residue: false);
                SetLetter('Z', 22, residue: false);
            }
        }

        private void SetLetter(char c, int letter, bool residue)
        {
            _charToLetterEx[c] = letter;
            _charToLetterEx[char.ToLowerInvariant(c)] = letter;
            _isResidue[c] = true;
            _isResidue[char.ToLowerInvariant(c)] = true;
            if (!residue)
            {
                _isWildcard[c] = true;
                _isWildcard[char.ToLowerInvariant(c)] = true;
            }
        }

        // params.cpp: SetPPScore -> SetPPDefaultParams, then the command-line overrides
        private void SetParameters()
        {
            var model = _o.ScoreModel;
            if (model == MuscleScoreModel.Default)
                model = _nuc ? MuscleScoreModel.SumOfPairsNucleotide : MuscleScoreModel.LogExpectation;
            if (_nuc && model != MuscleScoreModel.SumOfPairsNucleotide)
                throw new ArgumentException("Nucleotide sequences are aligned with the SPN score (MUSCLE switches the protein scores to SPN).");
            if (!_nuc && model == MuscleScoreModel.SumOfPairsNucleotide)
                throw new ArgumentException("The SPN score needs a nucleotide alphabet.");
            _pp = model;

            switch (model)
            {
                case MuscleScoreModel.LogExpectation:
                    Load(MuscleVtmlLa, 0f, 20);
                    _gapOpen = -2.9f;
                    _center = -0.52f;
                    _normalizeCounts = true;
                    break;
                case MuscleScoreModel.SumOfPairs:
                    Load(MusclePam200, 20f, 20);
                    _gapOpen = -1439f;
                    _center = 0f;
                    break;
                case MuscleScoreModel.SumOfPairsVtml:
                    Load(MuscleVtmlSp, 22f, 20, addInDouble: true);
                    _gapOpen = -300f;
                    _center = 0f;
                    break;
                default:
                    Load(MuscleNucSp, 60f, 4);
                    _gapOpen = _rna ? -420f : -400f;
                    _center = _rna ? -300f : 0f;
                    break;
            }

            _gapExtend = 0f;
            if (_o.GapOpen is { } go) _gapOpen = (float)go;
            if (_o.GapExtend is { } ge) _gapExtend = (float)ge;
            if (_o.Center is { } c) _center = (float)c;
            if (float.IsNaN(_gapOpen) || float.IsNaN(_gapExtend) || float.IsNaN(_center))
                throw new ArgumentException("MUSCLE gap / center parameters must be numbers.");

            _distance1 = _o.Distance1;
            if (_distance1 == MuscleDistance.Default)
                _distance1 = _nuc ? MuscleDistance.Kmer4_6 : MuscleDistance.Kmer6_6;
            if (_distance1 == MuscleDistance.Kmer4_6 && !_nuc)
                throw new ArgumentException("kmer4_6 requires a nucleotide alphabet.");
        }

        private void Load(double[] table, float center, int size, bool addInDouble = false)
        {
            _matrix = new float[32, 32];
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                {
                    double x = table[i * size + j];
                    // VTML_LA: (float) x; PAM200 / NUC_SP: (float) x + center; VTML_SP: (float) (x + center)
                    float v = (float)x;
                    if (addInDouble)
                        v = (float)(x + center);
                    else if (center != 0f)
                        v += center;
                    _matrix[i, j] = v;
                }
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        public (string[] Rows, bool Nucleotide) Run()
        {
            if (_n == 1)
            {
                var one = new[] { new string(_seq[0]) };
                FinalRowsInputOrder = DraftRowsInputOrder = one;
                FinalTreeOrder = [0];
                GuideTreeNewick = "(0);";
                return (_mHack[0] ? [RestoreM(one[0])] : one, _nuc);
            }

            var dist = ComputeKmerDistances();
            var tree = Upgma(dist, _o.Cluster1, Enumerable.Range(0, _n).ToArray());
            GuideTreeNewick = Newick(tree);

            var weights = CalcClustalWWeights(tree);
            var nodes = new MuscleNode[tree.NodeCount];
            ProgressiveAlign(tree, weights, nodes);
            var root = nodes[tree.Root];
            var rows = ReorderToDepthFirst(root, tree, nodes.Length);
            DraftRowsInputOrder = ToInputOrder(rows.ids, rows.rows);

            var msaIds = rows.ids;
            var msaRows = rows.rows;
            if (_o.MaxIterations >= 2 && _n > 2 && _o.MaxTrees != 0)
            {
                RefineTree(ref msaIds, ref msaRows, ref tree, nodes);
            }

            FinalTreeOrder = msaIds;
            FinalRowsInputOrder = ToInputOrder(msaIds, msaRows);
            var result = FinalRowsInputOrder.Select((r, i) => _mHack[i] ? RestoreM(r) : r).ToArray();
            FinalRowsInputOrder = result;
            DraftRowsInputOrder = DraftRowsInputOrder.Select((r, i) => _mHack[i] ? RestoreM(r) : r).ToArray();
            if (_o.TreeOrder)
                return (msaIds.Select(id => result[id]).ToArray(), _nuc);
            return (result, _nuc);
        }

        // MHackEnd: the first residue of the row becomes M again
        private static string RestoreM(string row)
        {
            var c = row.ToCharArray();
            for (int i = 0; i < c.Length; i++)
            {
                if (c[i] != '-')
                {
                    c[i] = 'M';
                    break;
                }
            }

            return new string(c);
        }

        private string[] ToInputOrder(int[] ids, string[] rows)
        {
            var r = new string[_n];
            for (int k = 0; k < ids.Length; k++)
                r[ids[k]] = rows[k];
            return r;
        }

        private static (int[] ids, string[] rows) ReorderToDepthFirst(MuscleNode root, MuscleTree tree, int nodeCount)
        {
            var byId = new Dictionary<int, char[]>();
            for (int k = 0; k < root.Ids.Length; k++)
                byId[root.Ids[k]] = root.Rows[k];
            var ids = new List<int>();
            foreach (int node in tree.PostOrder())
            {
                if (tree.IsLeaf(node))
                    ids.Add(tree.LeafId[node]);
            }

            return (ids.ToArray(), ids.Select(id => new string(byId[id])).ToArray());
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // fastdistmafft.cpp DistKmer6_6 / fastdistnuc.cpp DistKmer4_6

        // ResidueGroup[] of fastdistmafft.cpp, indexed by the letter number of the alphabet (A C D E F G H I K L M N P Q R S T V W Y X B Z gap): the table lists B and Z where the
        // enum has X and B, i.e. the wildcard X gets group 2, B group 2 and Z group 0 (as in the source)
        private static readonly int[] MuscleAminoGroup = { 0, 5, 2, 2, 4, 0, 3, 1, 3, 1, 1, 2, 0, 2, 3, 0, 0, 1, 4, 4, 2, 2, 0, 0 };
        private static readonly int[] MuscleNucleoGroup = { 0, 1, 2, 3, 4, 4, 4, 4 };

        internal float[,] ComputeKmerDistances()
        {
            const int tupleCount = 6 * 6 * 6 * 6 * 6 * 6;
            int[] groups = _distance1 == MuscleDistance.Kmer4_6 ? MuscleNucleoGroup : MuscleAminoGroup;
            var letters = new int[_n][];
            for (int s = 0; s < _n; s++)
            {
                letters[s] = _seq[s].Select(c => _charToLetterEx[c]).ToArray();
                if (_distance1 == MuscleDistance.Kmer4_6)
                {
                    for (int k = 0; k < letters[s].Length; k++)
                        if (letters[s][k] >= 4)
                            letters[s][k] = 4;
                }
            }

            int Tuple(int[] l, int p) =>
                groups[l[p + 5]] + groups[l[p + 4]] * 6 + groups[l[p + 3]] * 36 + groups[l[p + 2]] * 216 + groups[l[p + 1]] * 1296 + groups[l[p]] * 7776;

            var common = new uint[_n, _n];
            var count1 = new byte[tupleCount];
            var count2 = new byte[tupleCount];
            for (int s1 = 0; s1 < _n; s1++)
            {
                if (letters[s1].Length < 5)
                    continue;
                int tc1 = letters[s1].Length - 5;
                Array.Clear(count1);
                for (int p = 0; p < tc1; p++)
                    unchecked { count1[Tuple(letters[s1], p)]++; }

                for (int s2 = 0; s2 <= s1; s2++)
                {
                    if (letters[s2].Length < 5)
                        continue;
                    int tc2 = letters[s2].Length - 5;
                    Array.Clear(count2);
                    for (int p = 0; p < tc2; p++)
                        unchecked { count2[Tuple(letters[s2], p)]++; }

                    uint sum = 0;
                    for (int p = 0; p < tc2; p++)
                    {
                        int t = Tuple(letters[s2], p);
                        sum += Math.Min(count1[t], count2[t]);
                        count2[t] = 0;   // each unique tuple counted once
                    }

                    common[s1, s2] = sum;
                    common[s2, s1] = sum;
                }
            }

            var d = new float[_n, _n];
            for (int s1 = 0; s1 < _n; s1++)
            {
                double c11 = common[s1, s1];
                if (c11 == 0)
                    c11 = 1;
                for (int s2 = 0; s2 < s1; s2++)
                {
                    double c22 = common[s2, s2];
                    if (c22 == 0)
                        c22 = 1;
                    double d1 = 3.0 * (c11 - common[s1, s2]) / c11;
                    double d2 = 3.0 * (c22 - common[s1, s2]) / c22;
                    d[s1, s2] = (float)Math.Min(d1, d2);
                    d[s2, s1] = d[s1, s2];
                }
            }

            return d;
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // upgma2.cpp UPGMA2 (O(N^2), ties: lowest row)

        internal MuscleTree Upgma(float[,] dist, MuscleClusterMethod method, int[] ids)
        {
            const float bigDist = 1e29f;
            int n = ids.Length;
            int triangle = n * (n - 1) / 2;
            var g = new float[triangle];
            var nodeIndex = new int[n];
            var nearest = new int[n];
            var minDist = new float[n];
            int internalCount = n - 1;
            var uLeft = new int[internalCount];
            var uRight = new int[internalCount];
            var height = new float[internalCount];
            var leftLen = new float[internalCount];
            var rightLen = new float[internalCount];
            const int insane = 8888888;

            static int Tri(int a, int b) => a >= b ? b + a * (a - 1) / 2 : a + b * (b - 1) / 2;

            for (int i = 0; i < n; i++)
            {
                minDist[i] = bigDist;
                nodeIndex[i] = i;
                nearest[i] = insane;
            }

            for (int i = 0; i < internalCount; i++)
            {
                uLeft[i] = insane;
                uRight[i] = insane;
                leftLen[i] = bigDist;
                rightLen[i] = bigDist;
                height[i] = bigDist;
            }

            for (int i = 1; i < n; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    float d = dist[i, j];
                    g[Tri(i, j)] = d;
                    if (d < minDist[i])
                    {
                        minDist[i] = d;
                        nearest[i] = j;
                    }

                    if (d < minDist[j])
                    {
                        minDist[j] = d;
                        nearest[j] = i;
                    }
                }
            }

            float sueff = _o.Sueff;
            for (int k = 0; k < n - 1; k++)
            {
                int lmin = insane, rmin = insane;
                float dtMin = bigDist;
                for (int j = 0; j < n; j++)
                {
                    if (nodeIndex[j] == insane)
                        continue;
                    float d = minDist[j];
                    if (d < dtMin)
                    {
                        dtMin = d;
                        lmin = j;
                        rmin = nearest[j];
                    }
                }

                float dtNewMin = bigDist;
                int newNearest = insane;
                for (int j = 0; j < n; j++)
                {
                    if (j == lmin || j == rmin || nodeIndex[j] == insane)
                        continue;
                    int vL = Tri(lmin, j), vR = Tri(rmin, j);
                    float dL = g[vL], dR = g[vR];
                    float dtNew = method switch
                    {
                        MuscleClusterMethod.Upgma => (dL + dR) / 2,
                        MuscleClusterMethod.UpgmaMin => dL < dR ? dL : dR,
                        MuscleClusterMethod.UpgmaMax => dL > dR ? dL : dR,
                        _ => sueff * ((dL + dR) / 2) + (1 - sueff) * (dL < dR ? dL : dR),
                    };
                    if (nearest[j] == rmin)
                        nearest[j] = lmin;
                    g[vL] = dtNew;
                    if (dtNew < dtNewMin)
                    {
                        dtNewMin = dtNew;
                        newNearest = j;
                    }
                }

                float dLR = g[Tri(lmin, rmin)];
                float dHeightNew = dLR / 2;
                int left = nodeIndex[lmin], right = nodeIndex[rmin];
                float heightLeft = left < n ? 0 : height[left - n];
                float heightRight = right < n ? 0 : height[right - n];
                uLeft[k] = left;
                uRight[k] = right;
                leftLen[k] = dHeightNew - heightLeft;
                rightLen[k] = dHeightNew - heightRight;
                height[k] = dHeightNew;
                nodeIndex[lmin] = n + k;
                nearest[lmin] = newNearest;
                minDist[lmin] = dtNewMin;
                nodeIndex[rmin] = insane;
            }

            // Tree::Create
            var t = new MuscleTree { LeafCount = n };
            int nodes = 2 * n - 1;
            t.Left = new int[nodes];
            t.Right = new int[nodes];
            t.Parent = new int[nodes];
            t.EdgeLength = new double[nodes];
            t.LeafId = new int[nodes];
            Array.Fill(t.Left, -1);
            Array.Fill(t.Right, -1);
            Array.Fill(t.Parent, -1);
            for (int i = 0; i < n; i++)
                t.LeafId[i] = ids[i];
            for (int node = n; node < nodes; node++)
            {
                int v = node - n;
                t.Left[node] = uLeft[v];
                t.Right[node] = uRight[v];
                t.Parent[uLeft[v]] = node;
                t.Parent[uRight[v]] = node;
                t.EdgeLength[uLeft[v]] = leftLen[v];
                t.EdgeLength[uRight[v]] = rightLen[v];
            }

            t.Root = n - 2 + n;
            return t;
        }

        internal static string Newick(MuscleTree t)
        {
            string Rec(int node) => t.IsLeaf(node)
                ? t.LeafId[node].ToString()
                : $"({Rec(t.Left[node])},{Rec(t.Right[node])})";
            return t.LeafCount == 1 ? "(0);" : Rec(t.Root) + ";";
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // clwwt.cpp CalcClustalWWeights

        internal static float[] CalcClustalWWeights(MuscleTree tree)
        {
            int leafCount = tree.LeafCount;
            var w = new float[leafCount];
            if (leafCount == 1)
            {
                w[0] = 1f;
                return w;
            }

            if (leafCount == 2)
            {
                w[0] = 0.5f;
                w[1] = 0.5f;
                return w;
            }

            int nodes = tree.NodeCount;
            var leaves = new int[nodes];
            foreach (int node in tree.PostOrder())
                leaves[node] = tree.IsLeaf(node) ? 1 : leaves[tree.Left[node]] + leaves[tree.Right[node]];

            var strength = new double[nodes];
            for (int node = 0; node < nodes; node++)
            {
                if (node == tree.Root)
                {
                    strength[node] = 0.0;
                    continue;
                }

                strength[node] = tree.EdgeLength[node] / (double)leaves[node];
            }

            for (int leaf = 0; leaf < leafCount; leaf++)
            {
                // Weights[n] is indexed by the n-th leaf in node order, which is the leaf node itself (the id is looked up by the caller)
                double dWeight = 0;
                int node = leaf;
                while (node != tree.Root)
                {
                    dWeight += strength[node];
                    node = tree.Parent[node];
                }

                if (dWeight < 0.0001)
                    dWeight = 1.0;
                w[leaf] = (float)dWeight;
            }

            float sum = 0;
            for (int i = 0; i < leafCount; i++)
                sum += w[i];
            if (sum == 0)
                throw new InvalidOperationException("Normalize, sum=0");
            for (int i = 0; i < leafCount; i++)
                w[i] /= sum;
            return w;
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // progressalign (ProgressiveAlignE): leaf profiles, AlignTwoProfs along the guide tree

        private void ProgressiveAlign(MuscleTree tree, float[] weights, MuscleNode[] nodes)
        {
            foreach (int node in tree.PostOrder())
            {
                if (tree.IsLeaf(node))
                {
                    int id = tree.LeafId[node];
                    float w = weights[node];   // Weights[] of CalcClustalWWeights is indexed by leaf node; the leaf id equals the node index of the first tree
                    var n = new MuscleNode { Ids = [id], Rows = [(char[])_seq[id].Clone()], Length = _seq[id].Length, Weight = w };
                    // SetClustalWWeightsMuscle: the (single) sequence weight is normalised to 1 (f = 1 / w; w * f)
                    float wn = NormalizedSingleWeight(w);
                    n.Prof = ProfileFromRows(n.Rows, [wn]);
                    nodes[node] = n;
                }
                else
                {
                    nodes[node] = Merge(nodes[tree.Left[node]], nodes[tree.Right[node]]);
                }
            }
        }

        private static float NormalizedSingleWeight(float w)
        {
            if (w == 0f)
                return w;
            float f = 1f / w;
            return w * f;
        }

        private MuscleNode Merge(MuscleNode a, MuscleNode b)
        {
            var path = AlignTwoProfs(a, b, out var prof);
            var parent = new MuscleNode
            {
                Ids = a.Ids.Concat(b.Ids).ToArray(),
                Length = path.Count,
                Weight = a.Weight + b.Weight,
                Prof = prof,
            };
            parent.Rows = RowsFromPath(path, a, b);
            return parent;
        }

        private static char[][] RowsFromPath(List<(char Type, int A, int B)> path, MuscleNode a, MuscleNode b)
        {
            int ra = a.Rows.Length, rb = b.Rows.Length;
            var rows = new char[ra + rb][];
            for (int r = 0; r < rows.Length; r++)
                rows[r] = new char[path.Count];
            int colA = 0, colB = 0;
            for (int e = 0; e < path.Count; e++)
            {
                switch (path[e].Type)
                {
                    case 'M':
                        for (int r = 0; r < ra; r++) rows[r][e] = a.Rows[r][colA];
                        for (int r = 0; r < rb; r++) rows[ra + r][e] = b.Rows[r][colB];
                        colA++;
                        colB++;
                        break;
                    case 'D':
                        for (int r = 0; r < ra; r++) rows[r][e] = a.Rows[r][colA];
                        for (int r = 0; r < rb; r++) rows[ra + r][e] = '-';
                        colA++;
                        break;
                    default:
                        for (int r = 0; r < ra; r++) rows[r][e] = '-';
                        for (int r = 0; r < rb; r++) rows[ra + r][e] = b.Rows[r][colB];
                        colB++;
                        break;
                }
            }

            return rows;
        }

        // aligntwoprofs.cpp AlignTwoProfs: GlobalAlign (NWSmall) then AlignTwoProfsGivenPath with the weights wA / (wA + wB), wB / (wA + wB)
        private List<(char Type, int A, int B)> AlignTwoProfs(MuscleNode a, MuscleNode b, out MuscleProfPos[] profOut)
        {
            var path = NwSmall(a.Prof, a.Length, b.Prof, b.Length);
            float wA = a.Weight, wB = b.Weight;
            profOut = ProfileGivenPath(path, a.Prof, wA / (wA + wB), b.Prof, wB / (wA + wB));
            return path;
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // profilefrommsa.cpp ProfileFromMSA + msa2.cpp GetFractionalWeightedCounts

        private static bool IsGap(char c) => c == '-' || c == '.';

        internal MuscleProfPos[] ProfileFromRows(char[][] rows, float[] weights)
        {
            int seqCount = rows.Length;
            int colCount = rows[0].Length;
            var prof = new MuscleProfPos[colCount];
            for (int col = 0; col < colCount; col++)
            {
                var pp = new MuscleProfPos();
                var counts = pp.Counts;
                float wTotal = 0, fGap = 0;
                for (int s = 0; s < seqCount; s++)
                {
                    float w = weights[s];
                    char c = rows[s][col];
                    if (IsGap(c))
                    {
                        fGap += w;
                        continue;
                    }

                    if (_isWildcard[c])
                    {
                        int letter = _charToLetterEx[c];
                        if (!_nuc)
                        {
                            switch (letter)
                            {
                                case 21:   // AX_B: D or N
                                    counts[2] += w / 2;     // AX_D
                                    counts[11] += w / 2;    // AX_N
                                    break;
                                case 22:   // AX_Z: E or Q
                                    counts[3] += w / 2;     // AX_E
                                    counts[13] += w / 2;    // AX_Q
                                    break;
                                default:
                                    {
                                        float f = w / 20;
                                        for (int l = 0; l < 20; l++)
                                            counts[l] += f;
                                        break;
                                    }
                            }
                        }
                        else
                        {
                            // MSA::GetFractionalWeightedCounts compares the nucleotide letter with AX_R (14) and AX_Y (19): letter 14 (X) counts as G or A, every other wildcard as any (w / 20 each)
                            switch (letter)
                            {
                                case 14:
                                    counts[2] += w / 2;     // NX_G
                                    counts[0] += w / 2;     // NX_A
                                    break;
                                case 19:
                                    counts[1] += w / 2;
                                    counts[3] += w / 2;
                                    break;
                                default:
                                    {
                                        float f = w / 20;
                                        for (int l = 0; l < 4; l++)
                                            counts[l] += f;
                                        break;
                                    }
                            }
                        }

                        continue;
                    }

                    counts[_charToLetterEx[c]] += w;
                    wTotal += w;
                }

                pp.Occ = (float)(1.0 - fGap);
                if (_normalizeCounts && wTotal > 0)
                {
                    for (int l = 0; l < _alphaSize; l++)
                        counts[l] /= wTotal;
                }

                float startCount = 0;
                if (col == 0)
                {
                    for (int s = 0; s < seqCount; s++)
                        if (IsGap(rows[s][col]))
                            startCount += weights[s];
                }
                else
                {
                    for (int s = 0; s < seqCount; s++)
                        if (IsGap(rows[s][col]) && !IsGap(rows[s][col - 1]))
                            startCount += weights[s];
                }

                float endCount = 0;
                if (colCount - 1 == col)
                {
                    for (int s = 0; s < seqCount; s++)
                        if (IsGap(rows[s][col]))
                            endCount += weights[s];
                }
                else
                {
                    for (int s = 0; s < seqCount; s++)
                        if (IsGap(rows[s][col]) && !IsGap(rows[s][col + 1]))
                            endCount += weights[s];
                }

                float ll = 0, lg = 0, gl = 0, gg = 0;
                for (int s = 0; s < seqCount; s++)
                {
                    float w = weights[s];
                    bool here = !IsGap(rows[s][col]);
                    bool prev = col == 0 || !IsGap(rows[s][col - 1]);
                    if (here)
                    {
                        if (prev) ll += w; else gl += w;
                    }
                    else
                    {
                        if (prev) lg += w; else gg += w;
                    }
                }

                pp.LL = ll;
                pp.LG = lg;
                pp.GL = gl;
                pp.GG = gg;

                SortCounts(pp.Counts, pp.SortOrder);
                ComputeAaScores(pp);
                float sStartOcc = (float)(1.0 - startCount);
                float sEndOcc = (float)(1.0 - endCount);
                pp.GapOpen = sStartOcc * _gapOpen / 2;
                pp.GapClose = sEndOcc * _gapOpen / 2;
                prof[col] = pp;
            }

            return prof;
        }

        private void ComputeAaScores(MuscleProfPos pp)
        {
            for (int i = 0; i < _alphaSize; i++)
            {
                float sum = 0;
                for (int j = 0; j < _alphaSize; j++)
                    sum += pp.Counts[j] * _matrix[i, j];
                pp.AaScores[i] = sum;
            }
        }

        // profilefrommsa.cpp SortCounts
        private void SortCounts(float[] counts, int[] order)
        {
            for (int i = 0; i < _alphaSize; i++)
                order[i] = i;
            bool any = true;
            while (any)
            {
                any = false;
                for (int n = 0; n < _alphaSize - 1; n++)
                {
                    int i1 = order[n], i2 = order[n + 1];
                    if (counts[i1] < counts[i2])
                    {
                        order[n + 1] = i1;
                        order[n] = i2;
                        any = true;
                    }
                }
            }
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // aligngivenpath.cpp AlignTwoProfsGivenPath

        private static readonly MuscleProfPos PPStart = new() { LL = 1.0f, Occ = 1.0f };

        private MuscleProfPos[] ProfileGivenPath(List<(char Type, int A, int B)> path, MuscleProfPos[] pa, float wA, MuscleProfPos[] pb, float wB)
        {
            int edgeCount = path.Count;
            var pout = new MuscleProfPos[edgeCount];
            char prevType = 'M';
            for (int e = 0; e < edgeCount; e++)
            {
                var (type, la, lb) = path[e];
                var po = pout[e] = new MuscleProfPos();
                var ppa = la > 0 ? pa[la - 1] : PPStart;
                var ppb = lb > 0 ? pb[lb - 1] : PPStart;
                switch (type)
                {
                    case 'M':
                        SetFreqs(ppa, wA, ppb, wB, po);
                        switch (prevType)
                        {
                            case 'M':   // SetGapsMM
                                po.LL = wA * ppa.LL + wB * ppb.LL;
                                po.LG = wA * ppa.LG + wB * ppb.LG;
                                po.GL = wA * ppa.GL + wB * ppb.GL;
                                po.GG = wA * ppa.GG + wB * ppb.GG;
                                break;
                            case 'D':   // SetGapsDM
                                po.LL = wA * ppa.LL;
                                po.LG = wA * ppa.LG;
                                po.GL = wA * ppa.GL + wB * (ppb.LL + ppb.GL);
                                po.GG = wA * ppa.GG + wB * (ppb.LG + ppb.GG);
                                break;
                            default:    // SetGapsIM
                                po.LL = wB * ppb.LL;
                                po.LG = wB * ppb.LG;
                                po.GL = wB * ppb.GL + wA * (ppa.LL + ppa.GL);
                                po.GG = wB * ppb.GG + wA * (ppa.LG + ppa.GG);
                                break;
                        }

                        break;
                    case 'D':
                        SetFreqs(ppa, wA, ppb, 0f, po);
                        switch (prevType)
                        {
                            case 'M':   // SetGapsMD
                                po.LL = wA * ppa.LL;
                                po.LG = wA * ppa.LG + wB * (ppb.LL + ppb.GL);
                                po.GL = wA * ppa.GL;
                                po.GG = wA * ppa.GG + wB * (ppb.LG + ppb.GG);
                                break;
                            case 'D':   // SetGapsDD
                                po.LL = wA * ppa.LL;
                                po.LG = wA * ppa.LG;
                                po.GL = wA * ppa.GL;
                                po.GG = wA * ppa.GG + wB;
                                break;
                            default:    // SetGapsID
                                po.LL = 0;
                                po.LG = wB * ppb.GL + wB * ppb.LL;
                                po.GL = wA * ppa.GL + wA * ppa.LL;
                                po.GG = wA * (ppa.LG + ppa.GG) + wB * (ppb.LG + ppb.GG);
                                break;
                        }

                        break;
                    default:
                        SetFreqs(ppa, 0f, ppb, wB, po);
                        switch (prevType)
                        {
                            case 'M':   // SetGapsMI
                                po.LL = wB * ppb.LL;
                                po.LG = wB * ppb.LG + wA * (ppa.LL + ppa.GL);
                                po.GL = wB * ppb.GL;
                                po.GG = wB * ppb.GG + wA * (ppa.LG + ppa.GG);
                                break;
                            case 'D':   // SetGapsDI
                                po.LL = 0;
                                po.LG = wA * ppa.GL + wA * ppa.LL;
                                po.GL = wB * ppb.GL + wB * ppb.LL;
                                po.GG = wA * (ppa.LG + ppa.GG) + wB * (ppb.LG + ppb.GG);
                                break;
                            default:    // SetGapsII
                                po.LL = wB * ppb.LL;
                                po.LG = wB * ppb.LG;
                                po.GL = wB * ppb.GL;
                                po.GG = wB * ppb.GG + wA;
                                break;
                        }

                        break;
                }

                prevType = type;
            }

            // ProfScoresFromFreqs
            for (int i = 0; i < edgeCount; i++)
            {
                var pp = pout[i];
                SortCounts(pp.Counts, pp.SortOrder);
                pp.Occ = pp.LL + pp.GL;
                float fcOpen = pp.LG;
                float fcClose = i + 1 < edgeCount ? pout[i + 1].GL : pp.GG + pp.LG;
                pp.GapOpen = (float)((1.0 - fcOpen) * _gapOpen / 2.0);
                pp.GapClose = (float)((1.0 - fcClose) * _gapOpen / 2.0);
                ComputeAaScores(pp);
            }

            return pout;
        }

        private void SetFreqs(MuscleProfPos ppa, float wA, MuscleProfPos ppb, float wB, MuscleProfPos po)
        {
            if (_normalizeCounts)
            {
                float fA = ppa.Occ * wA / (wA + wB);
                float fB = ppb.Occ * wB / (wA + wB);
                float total = 0;
                for (int i = 0; i < 20; i++)
                {
                    float f = fA * ppa.Counts[i] + fB * ppb.Counts[i];
                    po.Counts[i] = f;
                    total += f;
                }

                if (total > 0)
                {
                    for (int i = 0; i < 20; i++)
                        po.Counts[i] /= total;
                }
            }
            else
            {
                for (int i = 0; i < 20; i++)
                    po.Counts[i] = wA * ppa.Counts[i] + wB * ppb.Counts[i];
            }
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // scorepp.cpp ScoreProfPos2 and termgaps.cpp SetTermGaps (TERMGAPS_Half)

        internal float ScoreProfPos2(MuscleProfPos ppa, MuscleProfPos ppb)
        {
            float score = 0;
            if (_pp == MuscleScoreModel.SumOfPairsNucleotide)
            {
                for (int n = 0; n < 4; n++)
                {
                    int letter = ppa.SortOrder[n];
                    float fc = ppa.Counts[letter];
                    if (fc == 0)
                        break;
                    score += fc * ppb.AaScores[letter];
                }

                return score - _center;
            }

            for (int n = 0; n < 20; n++)
            {
                int letter = ppa.SortOrder[n];
                float fc = ppa.Counts[letter];
                if (fc == 0)
                    break;
                score += fc * ppb.AaScores[letter];
            }

            if (_pp != MuscleScoreModel.LogExpectation)
                return score - _center;       // NS (SP) and SP (SV) are the same expression
            if (score == 0)
                return -2.5f;
            float logScore = MathF.Log(score);
            return (logScore - _center) * (ppa.Occ * ppb.Occ);
        }

        private static void SetTermGaps(MuscleProfPos[] prof, int length)
        {
            if (length == 0)
                return;
            var first = prof[0];
            var last = prof[length - 1];
            // TERMGAPS_Half: the terminal open / close penalties are zero (then multiplied by -1, which leaves a signed zero), unless they are -infinity
            if (first.GapOpen != MuscleMinusInfinity)
                first.GapOpen = 0;
            if (length > 1 && last.GapClose != MuscleMinusInfinity)
                last.GapClose = 0;
            if (first.GapOpen != MuscleMinusInfinity)
                first.GapOpen *= -1;
            if (length > 1 && last.GapClose != MuscleMinusInfinity)
                last.GapClose *= -1;
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // nwsmall.cpp NWSmall + bittraceback.cpp BitTraceBack

        internal List<(char Type, int A, int B)> NwSmall(MuscleProfPos[] pa, int lengthA, MuscleProfPos[] pb, int lengthB)
        {
            if (lengthA == 0 || lengthB == 0)
                throw new InvalidOperationException("Internal error, NWSmall: length=0");

            SetTermGaps(pa, lengthA);
            SetTermGaps(pb, lengthB);

            int prefixCountA = lengthA + 1, prefixCountB = lengthB + 1;
            float e = _gapExtend;
            var mCurr = new float[prefixCountB + 1];
            var mNext = new float[prefixCountB + 1];
            var mPrev = new float[prefixCountB + 1];
            var dRow = new float[prefixCountB + 1];
            var tb = new byte[prefixCountA][];
            for (int i = 0; i < prefixCountA; i++)
                tb[i] = new byte[prefixCountB];

            for (int j = 0; j <= lengthB; j++)
                dRow[j] = MuscleMinusInfinity;

            mPrev[0] = 0;
            for (int j = 1; j <= lengthB; j++)
                mPrev[j] = MuscleMinusInfinity;

            ResetFirst(mCurr);
            mCurr[1] = ScoreProfPos2(pa[0], pb[0]);
            // SetBitTBM(TB, 1, 1, 'M') = BIT_MM = 0: nothing to set
            for (int j = 2; j <= lengthB; j++)
            {
                mCurr[j] = ScoreProfPos2(pa[0], pb[j - 1]) + pb[0].GapOpen + (j - 2) * e + pb[j - 2].GapClose;
                tb[1][j] = (byte)((tb[1][j] & ~MuscleBitXm) | MuscleBitIm);
            }

            float iij;
            for (int i = 1; i < lengthA; i++)
            {
                var tbRow = tb[i];
                iij = MuscleMinusInfinity;
                dRow[0] = pa[0].GapOpen + (i - 1) * e;
                ResetFirst(mCurr);
                if (i == 1)
                {
                    mCurr[1] = ScoreProfPos2(pa[0], pb[0]);
                }
                else
                {
                    mCurr[1] = ScoreProfPos2(pa[i - 1], pb[0]) + pa[0].GapOpen + (i - 2) * e + pa[i - 2].GapClose;
                    tbRow[1] = (byte)((tbRow[1] & ~MuscleBitXm) | MuscleBitDm);
                }

                for (int j = 1; j < lengthB; j++)
                    mNext[j + 1] = ScoreProfPos2(pa[i], pb[j]);

                for (int j = 1; j < lengthB; j++)
                {
                    RecurseD(pa, i, j, dRow, mPrev, tbRow, e);
                    iij = RecurseI(pb, j, iij, mCurr, tbRow, e);
                    RecurseM(pa, pb, i, j, dRow, iij, mCurr, mNext, tb);
                }

                RecurseD(pa, i, lengthB, dRow, mPrev, tbRow, e);
                _ = RecurseI(pb, lengthB, iij, mCurr, tbRow, e);   // the I value of the last column of a row is not used again

                var tmp = mPrev;
                mPrev = mCurr;
                mCurr = mNext;
                mNext = tmp;
            }

            // special case for i = lengthA
            var lastRow = tb[lengthA];
            ResetFirst(mCurr);
            if (lengthA > 1)
                mCurr[1] = ScoreProfPos2(pa[lengthA - 1], pb[0]) + (lengthA - 2) * e + pa[0].GapOpen + pa[lengthA - 2].GapClose;
            else
                mCurr[1] = ScoreProfPos2(pa[lengthA - 1], pb[0]) + pa[0].GapOpen + pa[0].GapClose;
            lastRow[1] = (byte)((lastRow[1] & ~MuscleBitXm) | MuscleBitDm);

            dRow[0] = MuscleMinusInfinity;
            for (int j = 1; j <= lengthB; j++)
                RecurseD(pa, lengthA, j, dRow, mPrev, lastRow, e);

            iij = MuscleMinusInfinity;
            for (int j = 1; j <= lengthB; j++)
                iij = RecurseI(pb, j, iij, mCurr, lastRow, e);

            float mab = mCurr[lengthB];
            float dab = dRow[lengthB];
            float iab = iij;
            char edge = 'M';
            float best = mab;
            if (dab > best)
            {
                best = dab;
                edge = 'D';
            }

            if (iab > best)
                edge = 'I';

            return BitTraceBack(tb, lengthA, lengthB, edge);
        }

        private static void ResetFirst(float[] row) => row[0] = MuscleMinusInfinity;

        private static void RecurseD(MuscleProfPos[] pa, int i, int j, float[] dRow, float[] mPrev, byte[] tbRow, float e)
        {
            float dd = dRow[j] + e;
            float md = mPrev[j] + pa[i - 1].GapOpen;
            if (dd > md)
            {
                dRow[j] = dd;
            }
            else
            {
                dRow[j] = md;
                tbRow[j] = (byte)((tbRow[j] & ~MuscleBitXd) | MuscleBitMd);
            }
        }

        private static float RecurseI(MuscleProfPos[] pb, int j, float iij, float[] mCurr, byte[] tbRow, float e)
        {
            iij += e;
            float mi = mCurr[j - 1] + pb[j - 1].GapOpen;
            if (mi >= iij)
            {
                iij = mi;
                tbRow[j] = (byte)((tbRow[j] & ~MuscleBitXi) | MuscleBitMi);
            }

            return iij;
        }

        private static void RecurseM(MuscleProfPos[] pa, MuscleProfPos[] pb, int i, int j, float[] dRow, float iij, float[] mCurr, float[] mNext, byte[][] tb)
        {
            float dm = dRow[j] + pa[i - 1].GapClose;
            float im = iij + pb[j - 1].GapClose;
            float mm = mCurr[j];
            tb[i + 1][j + 1] = (byte)(tb[i + 1][j + 1] & ~MuscleBitXm);
            if (mm >= dm && mm >= im)
            {
                mNext[j + 1] += mm;
            }
            else if (dm >= mm && dm >= im)
            {
                mNext[j + 1] += dm;
                tb[i + 1][j + 1] |= MuscleBitDm;
            }
            else
            {
                mNext[j + 1] += im;
                tb[i + 1][j + 1] |= MuscleBitIm;
            }
        }

        private static List<(char Type, int A, int B)> BitTraceBack(byte[][] tb, int lengthA, int lengthB, char lastEdge)
        {
            var edges = new List<(char, int, int)>();
            int pla = lengthA, plb = lengthB;
            char type = lastEdge;
            while (true)
            {
                edges.Add((type, pla, plb));
                byte bits = tb[pla][plb];
                char next = type switch
                {
                    'M' => (bits & MuscleBitXm) switch { 0 => 'M', MuscleBitDm => 'D', MuscleBitIm => 'I', _ => throw new InvalidOperationException("Huh!?") },
                    'D' => (bits & MuscleBitXd) == MuscleBitMd ? 'M' : 'D',
                    _ => (bits & MuscleBitXi) == MuscleBitMi ? 'M' : 'I',
                };
                switch (type)
                {
                    case 'M':
                        if (pla == 0 || plb == 0)
                            throw new InvalidOperationException("BitTraceBack MA=0");
                        pla--;
                        plb--;
                        break;
                    case 'D':
                        if (pla == 0)
                            throw new InvalidOperationException("BitTraceBack DA=0");
                        pla--;
                        break;
                    default:
                        if (plb == 0)
                            throw new InvalidOperationException("BitTraceBack IB=0");
                        plb--;
                        break;
                }

                if (pla == 0 && plb == 0)
                    break;
                type = next;
            }

            edges.Reverse();
            return edges;
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // refinetreee.cpp RefineTreeE, difftreese.cpp DiffTreesE, realigndiffse.cpp RealignDiffsE

        private void RefineTree(ref int[] msaIds, ref string[] msaRows, ref MuscleTree tree, MuscleNode[] oldNodes)
        {
            int rowCount = msaIds.Length;
            // TreeFromMSA(msa, Tree2, Cluster2, DISTANCE_PctIdKimura): UPGMA2 over the rows in MSA order, leaf ids = the sequence ids
            var dist = new float[rowCount, rowCount];
            for (int i = 1; i < rowCount; i++)
                for (int j = 0; j < i; j++)
                {
                    float pctId = (float)PctIdentityPair(msaRows[i], msaRows[j]);
                    dist[i, j] = (float)KimuraDist(pctId);
                    dist[j, i] = dist[i, j];
                }

            var tree2 = Upgma(dist, _o.Cluster2, msaIds);

            // DiffTreesE(Tree2, tree)
            int nodeCount = tree2.NodeCount;
            var idToOld = new int[_n];
            for (int old = 0; old < nodeCount; old++)
                if (tree.IsLeaf(old))
                    idToOld[tree.LeafId[old]] = old;

            var newToOld = new uint[nodeCount];
            for (int nw = 0; nw < nodeCount; nw++)
                newToOld[nw] = tree2.IsLeaf(nw) ? (uint)idToOld[tree2.LeafId[nw]] : MuscleNodeChanged;

            foreach (int nw in tree2.PostOrder())
            {
                if (tree2.IsLeaf(nw))
                    continue;
                int nl = tree2.Left[nw], nr = tree2.Right[nw];
                if (newToOld[nl] == MuscleNodeChanged)
                    continue;
                if (newToOld[nr] == MuscleNodeChanged)
                    continue;
                int pl = tree.Parent[newToOld[nl]], pr = tree.Parent[newToOld[nr]];
                newToOld[nw] = pl == pr ? (uint)pl : MuscleNodeChanged;
            }

            if (newToOld[tree2.Root] != MuscleNodeChanged)
                return;   // same tree: nothing to realign

            var newNodes = new MuscleNode?[nodeCount];
            for (int nw = 0; nw < nodeCount; nw++)
            {
                if (newToOld[nw] == MuscleNodeChanged)
                    continue;
                newNodes[nw] = oldNodes[newToOld[nw]];   // MakeNode: profile, rows, length and weight of the unchanged node (the L / R swap does not change them)
            }

            foreach (int nw in tree2.PostOrder())
            {
                if (newToOld[nw] != MuscleNodeChanged)
                    continue;
                var node1 = newNodes[tree2.Left[nw]]!;
                var node2 = newNodes[tree2.Right[nw]]!;
                newNodes[nw] = Merge(node1, node2);
            }

            var root = newNodes[tree2.Root]!;
            var (ids, rows) = ReorderToDepthFirst(root, tree2, nodeCount);
            msaIds = ids;
            msaRows = rows;
            tree = tree2;
        }

        // cons.cpp MSA::GetPctIdentityPair (characters compared as they are, gaps skipped)
        private static double PctIdentityPair(string a, string b)
        {
            int pos = 0, same = 0;
            for (int i = 0; i < a.Length; i++)
            {
                char c1 = a[i], c2 = b[i];
                if (IsGap(c1) || IsGap(c2))
                    continue;
                if (c1 == c2)
                    same++;
                pos++;
            }

            return pos == 0 ? 0 : (double)same / pos;
        }

        // msadistkimura.cpp KimuraDist (ClustalW's dayhoff.h table above 75 % divergence)
        private static readonly int[] DayhoffPams = { 195, 196, 197, 198, 199, 200, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 226, 227, 228, 229, 230, 231, 232, 233, 234, 236, 237, 238, 239, 240, 241, 243, 244, 245, 246, 248, 249, 250, 252, 253, 254, 255, 257, 258, 260, 261, 262, 264, 265, 267, 268, 270, 271, 273, 274, 276, 277, 279, 281, 282, 284, 285, 287, 289, 291, 292, 294, 296, 298, 299, 301, 303, 305, 307, 309, 311, 313, 315, 317, 319, 321, 323, 325, 328, 330, 332, 335, 337, 339, 342, 344, 347, 349, 352, 354, 357, 360, 362, 365, 368, 371, 374, 377, 380, 383, 386, 389, 393, 396, 399, 403, 407, 410, 414, 418, 422, 426, 430, 434, 438, 442, 447, 451, 456, 461, 466, 471, 476, 482, 487, 493, 498, 504, 511, 517, 524, 531, 538, 545, 553, 560, 569, 577, 586, 595, 605, 615, 626, 637, 649, 661, 675, 688, 703, 719, 736, 754, 775, 796, 819, 845, 874, 907, 945, 988 };

        internal static double KimuraDist(double pctId)
        {
            double p = 1 - pctId;
            if (p < 0.75)
                return -Math.Log(1 - p - p * p / 5);
            if (p > 0.93)
                return 10.0;
            int index = (int)((p - 0.75) * 1000 + 0.5);
            if (index < 0 || index >= DayhoffPams.Length)
                throw new InvalidOperationException("Internal error in MSADistKimura::ComputeDist");
            return DayhoffPams[index] / 100.0;
        }
    }

    #endregion
}
