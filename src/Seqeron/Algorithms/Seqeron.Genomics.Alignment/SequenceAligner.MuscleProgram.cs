using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    // The other command lines of MUSCLE 3.8.1551 that run the alignment code (muscle.cpp Run):
    //   muscle -in x -out y [options]                      MuscleRun (the option set of MuscleOptions, incl. -usetree, -matrix, -tree1/-tree2, -distmx1/2, -clusteronly, -computeweights, -brenner)
    //   muscle -profile -in1 a -in2 b                      MuscleProfileAlign   (profile.cpp)
    //   muscle -refine -in x                               MuscleRefine         (refine.cpp)
    //   muscle -spscore x                                  MuscleScoreAlignment (dosp.cpp)
    //   muscle -ppscore -in1 a -in2 b                      MuscleProfileScore   (ppscore.cpp; nucleotides: the program never sets a protein score matrix in this mode)
    //   muscle -maketree -in x -out tree                   MuscleGuideTree      (maketree.cpp)
    //   muscle -refinew -in x                              MuscleRefineWindows  (refinew.cpp; nucleotides, see the remarks)
    //   muscle -sw -in1 a -in2 b                           MuscleLocalAlign     (local.cpp / sw.cpp; nucleotides)

    /// <summary>
    /// Runs MUSCLE 3.8.1551 on unaligned sequences (<c>muscle -in x.fa -out y.fa</c> with the options of <paramref name="options"/>) and returns what the program writes: the alignment and,
    /// when the options ask for them, the guide trees, distance matrices and sequence weights.
    /// </summary>
    public static MuscleRunResult MuscleRun(IEnumerable<string> sequences, MuscleOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var o = options ?? new MuscleOptions();
        var list = sequences.ToList();
        if (list.Count == 0)
            return new MuscleRunResult();
        var engine = new MuscleEngine(o, list);
        var (rows, nucleotide) = engine.Run();
        MultipleAlignmentResult alignment = MultipleAlignmentResult.Empty;
        if (!engine.Stopped)
        {
            var sp = nucleotide ? MatrixScoring.Emboss(SubstitutionMatrix.EdnaFull, 10, 0.5) : MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.2);
            string[] scored = nucleotide ? rows.Select(r => r.Replace('X', 'N')).ToArray() : rows;
            double exact = rows.Length == 1 ? 0.0 : MatrixSumOfPairsScore(scored, sp, MsaGapModel.Affine);
            alignment = new MultipleAlignmentResult(rows, BuildConsensus(rows.ToList(), rows.Max(r => r.Length), anySymbol: true), (int)Math.Round(exact, MidpointRounding.AwayFromZero)) { ExactTotalScore = exact };
        }

        return new MuscleRunResult
        {
            Alignment = alignment,
            Names = engine.Names,
            OutputOrder = engine.FinalTreeOrder,
            Tree1 = engine.Tree1Text,
            Tree2 = engine.Tree2Text,
            DistanceMatrix1 = engine.Dist1Text,
            DistanceMatrix2 = engine.Dist2Text,
            WeightsFile = engine.WeightsText,
            ReferenceBehaviourUndefined = engine.ReadsUninitialisedMemory,
        };
    }

    /// <summary>
    /// <c>muscle -profile -in1 a.fa -in2 b.fa</c>: aligns two existing alignments (profiles) by profile-profile dynamic programming; returns the rows of <paramref name="alignment1"/> followed by
    /// those of <paramref name="alignment2"/>, in the combined columns. The sequence weights of each profile are <c>-weight1</c> (default ClustalW weights of a guide tree built from the profile with
    /// <c>-cluster2</c> / <c>-distance2</c> / <c>-root1</c>).
    /// </summary>
    public static string[] MuscleProfileAlign(IReadOnlyList<string> alignment1, IReadOnlyList<string> alignment2, MuscleOptions? options = null)
    {
        var o = options ?? new MuscleOptions();
        var engine = MuscleEngine.ForMsa(o, alignment1);
        return engine.ProfileProfileMode(alignment1, alignment2);
    }

    /// <summary>
    /// <c>muscle -refine -in x.afa</c>: refines an existing alignment with <see cref="MuscleOptions.MaxIterations"/> rounds of restricted-partition refinement (the program's default is 8) on a
    /// guide tree built from the alignment (<c>-cluster2</c> / <c>-distance2</c> / <c>-root2</c>), sequence weights <c>-weight1</c>. Rows are returned in the order the program writes them.
    /// </summary>
    public static string[] MuscleRefine(IReadOnlyList<string> alignment, MuscleOptions? options = null)
    {
        var o = options ?? new MuscleOptions { MaxIterations = 8 };
        return MuscleEngine.ForMsa(o, alignment).RefineMode(alignment);
    }

    /// <summary><c>muscle -spscore x.afa</c>: the sum-of-pairs objective score (<c>-objscore sp</c> arithmetic, ClustalW weights of a guide tree built from the alignment) printed as "SP=%.4g".</summary>
    public static double MuscleScoreAlignment(IReadOnlyList<string> alignment, MuscleOptions? options = null)
    {
        var o = options ?? new MuscleOptions();
        return MuscleEngine.ForMsa(o, alignment).SpScoreMode(alignment);
    }

    /// <summary>
    /// <c>muscle -ppscore -in1 a.afa -in2 b.afa</c>: the profile-profile (dynamic-programming) score of two alignments of the same length, "Score=%.4g". Nucleotides only: for a protein the program never
    /// initialises the score matrix in this mode and crashes.
    /// </summary>
    public static double MuscleProfileScore(IReadOnlyList<string> alignment1, IReadOnlyList<string> alignment2, MuscleOptions? options = null)
    {
        var o = options ?? new MuscleOptions();
        return MuscleEngine.ForMsa(o, alignment1).ProfileScoreMode(alignment1, alignment2);
    }

    /// <summary><c>muscle -maketree -in x.afa -out t</c>: the guide tree of an alignment (<c>-cluster2</c>, <c>-distance2</c>, <c>-root2</c>) in the Newick format of the program's tree files; leaf names are <see cref="MuscleOptions.SequenceNames"/> (default s0, s1, ...).</summary>
    public static string MuscleGuideTree(IReadOnlyList<string> alignment, MuscleOptions? options = null)
    {
        var o = options ?? new MuscleOptions();
        return MuscleEngine.ForMsa(o, alignment).MakeTreeMode(alignment);
    }

    /// <summary>
    /// <c>muscle -refinew -in x.afa</c> (<see cref="MuscleOptions.RefineWindow"/>, <see cref="MuscleOptions.FromWindow"/>, <see cref="MuscleOptions.ToWindow"/>, <see cref="MuscleOptions.WindowOffset"/>): every window of
    /// columns is re-aligned from scratch with the full pipeline (<c>MUSCLE()</c>) and the windows are concatenated. Nucleotides only (for a protein the program never initialises the score matrix and
    /// crashes); a window in which a sequence has no residue is not defined by the program and is rejected.
    /// </summary>
    public static string[] MuscleRefineWindows(IReadOnlyList<string> alignment, MuscleOptions? options = null)
    {
        var o = options ?? new MuscleOptions { MaxIterations = 8 };
        return MuscleEngine.ForMsa(o, alignment).RefineWindowsMode(alignment);
    }

    /// <summary>
    /// <c>muscle -sw -in1 a.afa -in2 b.afa</c>: the local (Smith-Waterman) alignment of two profiles; the residues outside the local alignment are written in lower case and padded with '.'. Nucleotides only
    /// (for a protein the program never initialises the score matrix and crashes). The guide trees of the two profiles are built with <c>-cluster1</c> / <c>-distance1</c> / <c>-root1</c>, which for a
    /// profile of two or more sequences needs <c>-distance1 pctid_kimura</c> (the default nucleotide k-mer distance is not defined for aligned sequences and the program stops).
    /// Rows of the first profile followed by those of the second.
    /// </summary>
    public static string[] MuscleLocalAlign(IReadOnlyList<string> alignment1, IReadOnlyList<string> alignment2, MuscleOptions? options = null)
    {
        var o = options ?? new MuscleOptions();
        return MuscleEngine.ForMsa(o, alignment1).LocalAlignMode(alignment1, alignment2);
    }

    /// <summary>
    /// <c>muscle -profdb -in1 a.afa -in2 db.fa</c> (profdb.cpp: align every sequence of a database to the profile one by one). The program of 3.8.1551 never initialises the residue alphabet in
    /// this mode and crashes (segmentation fault) on every input - protein, DNA and RNA alike (5 / 5 runs of the oracle) - so the port refuses it instead of inventing an unverified behaviour.
    /// </summary>
    /// <exception cref="NotSupportedException">Always.</exception>
    public static string[] MuscleProfileDatabase(IReadOnlyList<string> alignment, IReadOnlyList<string> database, MuscleOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(alignment);
        ArgumentNullException.ThrowIfNull(database);
        throw new NotSupportedException("muscle -profdb: MUSCLE 3.8.1551 never initialises the residue alphabet in this mode and crashes on every input.");
    }

    internal sealed partial class MuscleEngine
    {
        // ---------------------------------------------------------------------------------------------------------------------------------
        // MSA modes: construction (the alphabet comes from the alignment, msa.cpp GuessAlpha; no sequences, no methionine hack)

        public static MuscleEngine ForMsa(MuscleOptions options, IReadOnlyList<string> rows)
        {
            if (rows.Count == 0)
                throw new ArgumentException("No sequences in input alignment.", nameof(rows));
            return new MuscleEngine(options, rows, rows.Count);
        }

        private MuscleEngine(MuscleOptions options, IReadOnlyList<string> msaRows, int n)
        {
            _o = options;
            _n = n;
            if (options.Stable)
                throw new NotSupportedException("-stable not supported in this version of muscle (MUSCLE 3.8.1551 refuses it).");
            if (options.MaxIterations < 1)
                throw new ArgumentException("MaxIterations must be at least 1.", nameof(options));
            if (options.MaxTrees < 0)
                throw new ArgumentException("MaxTrees must be non-negative.", nameof(options));
            var type = options.SequenceType;
            if (type == MuscleSequenceType.Auto)
                type = GuessAlphabetMsa(msaRows);
            _nuc = type != MuscleSequenceType.Protein;
            _rna = type == MuscleSequenceType.Rna;
            _gNuc = _nuc;
            _alphaSize = _nuc ? 4 : 20;
            _wildcardChars = _nuc ? "MRWSYKVHDBXN" : "BXZ";
            _wildcardLetterBase = _nuc ? 4 : 20;
            InitAlphabet();
            _seq = [];
            _mHack = new bool[n];
            _names = options.SequenceNames is { } nm ? nm.ToArray() : Enumerable.Range(0, n).Select(i => "s" + i).ToArray();
            SetParameters();
        }

        // msa.cpp MSA::GuessAlpha: the first 100 non-gap letters; RNA is tested before DNA
        private static MuscleSequenceType GuessAlphabetMsa(IReadOnlyList<string> rows)
        {
            const int charCount = 100, minNucleoPct = 95;
            int dna = 0, rna = 0, total = 0;
            foreach (char c in string.Concat(rows))
            {
                if (c == '-' || c == '.')
                    continue;
                if ("AGCTNagctn".Contains(c))
                    dna++;
                if ("AGCUNagcun".Contains(c))
                    rna++;
                total++;
                if (total >= charCount)
                    break;
            }

            if (total != 0 && rna * 100 / total >= minNucleoPct)
                return MuscleSequenceType.Rna;
            if (total != 0 && dna * 100 / total >= minNucleoPct)
                return MuscleSequenceType.Dna;
            return MuscleSequenceType.Protein;
        }

        // MSA::FixAlpha: letters that are neither residues nor gaps become the wildcard; the case of the letters is kept
        private MuscleMsa PrepareMsa(IReadOnlyList<string> rows)
        {
            if (rows.Count == 0)
                throw new ArgumentException("No sequences in input alignment.");
            int cols = rows[0].Length;
            var chars = new char[rows.Count][];
            for (int i = 0; i < rows.Count; i++)
            {
                if (rows[i].Length != cols)
                    throw new ArgumentException("The rows of an alignment must have the same length.");
                var r = rows[i].ToCharArray();
                for (int c = 0; c < r.Length; c++)
                {
                    if (r[c] == '.')
                        r[c] = '-';
                    else if (r[c] != '-' && !_isResidue[r[c]])
                        r[c] = _nuc ? 'N' : 'X';
                }

                chars[i] = r;
            }

            return new MuscleMsa(Enumerable.Range(0, rows.Count).ToArray(), chars);
        }

        private static string[] ToStrings(MuscleMsa m) => m.Rows.Select(r => new string(r)).ToArray();

        // profile.cpp / ppscore.cpp ProfileFromMSALocal: ids 0..n-1, the guide tree of the profile (when the weighting needs a tree) for the ClustalW weights
        private MuscleProfPos[] ProfileFromMsaLocal(MuscleMsa m, bool alwaysTree)
        {
            for (int i = 0; i < m.SeqCount; i++)
                m.Ids[i] = i;
            bool treeNeeded = _o.Weight2 is MuscleSeqWeight.ClustalW or MuscleSeqWeight.ThreeWay;
            if (alwaysTree || treeNeeded)
            {
                var tree = TreeFromMsaRows(m.Ids, m.Rows, _o.Root1);
                SetMuscleTree(tree);
            }

            return ProfileFromMsa(m);
        }

        public string[] ProfileProfileMode(IReadOnlyList<string> rows1, IReadOnlyList<string> rows2)
        {
            _curWeight = _o.Weight1;
            var m1 = PrepareMsa(rows1);
            var m2 = PrepareMsa(rows2);
            var p1 = ProfileFromMsaLocal(m1, alwaysTree: false);
            var p2 = ProfileFromMsaLocal(m2, alwaysTree: false);
            var path = NwSmall(p1, m1.ColCount, p2, m2.ColCount);
            return RowsFromPath(path, m1.Rows, m2.Rows).Select(r => new string(r)).ToArray();
        }

        public string[] RefineMode(IReadOnlyList<string> rows)
        {
            _curWeight = _o.Weight1;
            var msa = PrepareMsa(rows);
            var tree = TreeFromMsaRows(msa.Ids, msa.Rows);
            SetMuscleTree(tree);
            if (_o.Anchors)
                RefineVert(msa, tree, _o.MaxIterations);
            else
                RefineHoriz(msa, tree, _o.MaxIterations, false, false);
            return ToStrings(msa);
        }

        public double SpScoreMode(IReadOnlyList<string> rows)
        {
            _curWeight = _o.Weight1;
            var msa = PrepareMsa(rows);
            var tree = TreeFromMsaRows(msa.Ids, msa.Rows);
            SetMuscleTree(tree);
            SetMsaWeights(msa);
            return ObjScoreSp(msa);
        }

        public double ProfileScoreMode(IReadOnlyList<string> rows1, IReadOnlyList<string> rows2)
        {
            if (!_nuc)
                throw new NotSupportedException("muscle -ppscore on proteins: MUSCLE 3.8.1551 never initialises the score matrix in this mode (SetPPScore is only called for nucleotides) and crashes.");
            _curWeight = _o.Weight1;
            var m1 = PrepareMsa(rows1);
            var m2 = PrepareMsa(rows2);
            if (m1.ColCount != m2.ColCount)
                throw new ArgumentException("Profiles must have the same length");
            var p1 = ProfileFromMsaLocal(m1, alwaysTree: true);
            var p2 = ProfileFromMsaLocal(m2, alwaysTree: true);
            return ObjScoreDpProfiles(p1, p2, m1.ColCount);
        }

        // local.cpp Local / sw.cpp SW / tracebacksw.cpp TraceBackSW / aligngivenpathsw.cpp AlignTwoMSAsGivenPathSW
        public string[] LocalAlignMode(IReadOnlyList<string> rows1, IReadOnlyList<string> rows2)
        {
            if (!_nuc)
            {
                // local.cpp: SetPPScore is only called for nucleotides, so -distance1 is not parsed (it stays kmer6_6, the global default) and the score matrix is never set. A profile of two or more
                // sequences stops in TreeFromMSA with "DistCalcMSA: Invalid DISTANCE_1" (exit status 2, oracle R610); with one-sequence profiles the program reaches the profile score and crashes (R611, R612).
                if (rows1.Count > 1 || rows2.Count > 1)
                    throw new ArgumentException("DistCalcMSA: Invalid DISTANCE_1 - muscle -sw on proteins does not parse -distance1 (kmer6_6 is not a distance of an aligned pair).");
                throw new NotSupportedException("muscle -sw on proteins: MUSCLE 3.8.1551 never initialises the score matrix in this mode (SetPPScore is only called for nucleotides) and crashes.");
            }

            _curWeight = _o.Weight1;
            var m1 = PrepareMsa(rows1);
            var m2 = PrepareMsa(rows2);
            MuscleProfPos[] LocalProfile(MuscleMsa m)
            {
                for (int i = 0; i < m.SeqCount; i++)
                    m.Ids[i] = i;
                var distance = _o.Distance1 == MuscleDistance.Default ? MuscleDistance.Kmer4_6 : _o.Distance1;
                var tree = TreeFromMsaRows(m.Ids, m.Rows, _o.Cluster1, distance, _o.Root1, out _);
                SetMuscleTree(tree);
                return ProfileFromMsa(m);
            }

            var pa = LocalProfile(m1);
            var pb = LocalProfile(m2);
            var path = SmithWaterman(pa, m1.ColCount, pb, m2.ColCount);
            return CombineLocal(path, m1, m2);
        }

        private List<(char Type, int A, int B)> SmithWaterman(MuscleProfPos[] pa, int lengthA, MuscleProfPos[] pb, int lengthB)
        {
            int cA = lengthA + 1, cB = lengthB + 1;
            var dpm = new float[cA * cB];
            var dpd = new float[cA * cB];
            var dpi = new float[cA * cB];
            int Ix(int a, int b) => a * cB + b;
            dpm[Ix(0, 0)] = 0;
            dpd[Ix(0, 0)] = MuscleMinusInfinity;
            dpi[Ix(0, 0)] = MuscleMinusInfinity;
            foreach (var (a, b) in new[] { (1, 0), (0, 1) })
            {
                dpm[Ix(a, b)] = MuscleMinusInfinity;
                dpd[Ix(a, b)] = MuscleMinusInfinity;
                dpi[Ix(a, b)] = MuscleMinusInfinity;
            }

            for (int a = 2; a < cA; a++)
            {
                dpm[Ix(a, 0)] = MuscleMinusInfinity;
                dpd[Ix(a, 0)] = MuscleMinusInfinity;
                dpi[Ix(a, 0)] = MuscleMinusInfinity;
            }

            for (int b = 2; b < cB; b++)
            {
                dpm[Ix(0, b)] = MuscleMinusInfinity;
                dpd[Ix(0, b)] = MuscleMinusInfinity;
                dpi[Ix(0, b)] = MuscleMinusInfinity;
            }

            float scoreMax = MuscleMinusInfinity;
            int maxA = -1, maxB = -1;
            float gapCloseB = MuscleMinusInfinity;
            for (int b = 1; b < cB; b++)
            {
                var ppb = pb[b - 1];
                float gapCloseA = MuscleMinusInfinity;
                for (int a = 1; a < cA; a++)
                {
                    var ppa = pa[a - 1];
                    float scoreLL = ScoreProfPos2(ppa, ppb);
                    float mm = dpm[Ix(a - 1, b - 1)];
                    float dm = dpd[Ix(a - 1, b - 1)] + gapCloseA;
                    float im = dpi[Ix(a - 1, b - 1)] + gapCloseB;
                    float best;
                    if (mm >= dm && mm >= im)
                        best = mm;
                    else if (dm >= mm && dm >= im)
                        best = dm;
                    else
                        best = im;
                    if (best < 0)
                        best = 0;
                    best += scoreLL;
                    if (best > scoreMax)
                    {
                        scoreMax = best;
                        maxA = a;
                        maxB = b;
                    }

                    dpm[Ix(a, b)] = best;
                    float md = dpm[Ix(a - 1, b)] + pa[a - 1].GapOpen;
                    float dd = dpd[Ix(a - 1, b)];
                    dpd[Ix(a, b)] = md >= dd ? md : dd;
                    float mi = dpm[Ix(a, b - 1)] + pb[b - 1].GapOpen;
                    float ii = dpi[Ix(a, b - 1)];
                    dpi[Ix(a, b)] = mi >= ii ? mi : ii;
                    gapCloseA = ppa.GapClose;
                }

                gapCloseB = ppb.GapClose;
            }

            return TraceBackSw(pa, pb, dpm, dpd, dpi, cB, maxA, maxB);
        }

        private List<(char Type, int A, int B)> TraceBackSw(MuscleProfPos[] pa, MuscleProfPos[] pb, float[] dpm, float[] dpd, float[] dpi, int cB, int maxA, int maxB)
        {
            int Ix(int a, int b) => a * cB + b;
            static bool Eq(float a, float b) => Math.Abs(a - b) < 0.1;
            var edges = new List<(char, int, int)>();
            int la = maxA, lb = maxB;
            char edge = 'M';
            while (edge != 'S')
            {
                edges.Add((edge, la, lb));
                char prev;
                int pla = la, plb = lb;
                switch (edge)
                {
                    case 'M':
                        {
                            float score = dpm[Ix(la, lb)];
                            float scoreMatch = ScoreProfPos2(pa[la - 1], pb[lb - 1]);
                            float sm = la == 1 && lb == 1 ? scoreMatch : MuscleMinusInfinity;
                            float mm = MuscleMinusInfinity, dm = MuscleMinusInfinity, im = MuscleMinusInfinity;
                            if (la > 1 && lb > 1)
                                mm = dpm[Ix(la - 1, lb - 1)] + scoreMatch;
                            if (la > 1)
                                dm = dpd[Ix(la - 1, lb - 1)] + pa[la - 2].GapClose + scoreMatch;
                            if (lb > 1)
                                im = dpi[Ix(la - 1, lb - 1)] + pb[lb - 2].GapClose + scoreMatch;
                            if (Eq(mm, score))
                                prev = 'M';
                            else if (Eq(dm, score))
                                prev = 'D';
                            else if (Eq(im, score))
                                prev = 'I';
                            else if (Eq(sm, score) || Eq(scoreMatch, score))
                                prev = 'S';
                            else
                                throw new InvalidOperationException("TraceBack2: failed to match M");
                            pla--;
                            plb--;
                            break;
                        }

                    case 'D':
                        {
                            float score = dpd[Ix(la, lb)];
                            float md = MuscleMinusInfinity, dd = MuscleMinusInfinity, sd = MuscleMinusInfinity;
                            if (lb == 0)
                                sd = la == 1 ? pa[0].GapOpen : dpd[Ix(la - 1, 0)];
                            if (la > 1)
                            {
                                md = dpm[Ix(la - 1, lb)] + pa[la - 1].GapOpen;
                                dd = dpd[Ix(la - 1, lb)];
                            }

                            if (Eq(score, md))
                                prev = 'M';
                            else if (Eq(score, dd))
                                prev = 'D';
                            else if (Eq(score, sd))
                                prev = 'S';
                            else
                                throw new InvalidOperationException("TraceBack2: failed to match D");
                            pla--;
                            break;
                        }

                    default:
                        {
                            float score = dpi[Ix(la, lb)];
                            float mi = MuscleMinusInfinity, ii = MuscleMinusInfinity, si = MuscleMinusInfinity;
                            if (la == 0)
                                si = lb == 1 ? pb[0].GapOpen : dpi[Ix(0, lb - 1)];
                            if (lb > 1)
                            {
                                mi = dpm[Ix(la, lb - 1)] + pb[lb - 1].GapOpen;
                                ii = dpi[Ix(la, lb - 1)];
                            }

                            if (Eq(score, mi))
                                prev = 'M';
                            else if (Eq(score, ii))
                                prev = 'I';
                            else if (Eq(score, si))
                                prev = 'S';
                            else
                                throw new InvalidOperationException("TraceBack2: failed to match I");
                            plb--;
                            break;
                        }
                }

                edge = prev;
                la = pla;
                lb = plb;
            }

            edges.Reverse();
            return edges;
        }

        private static string[] CombineLocal(List<(char Type, int A, int B)> path, MuscleMsa a, MuscleMsa b)
        {
            int ra = a.SeqCount, rb = b.SeqCount;
            var rows = new List<char>[ra + rb];
            for (int i = 0; i < rows.Length; i++)
                rows[i] = [];
            int colA = 0, colB = 0;

            static char Unalign(char c) => c is '-' or '.' ? c : char.ToLowerInvariant(c);

            void Terminals(int countA, int countB)
            {
                int newCount = Math.Max(countA, countB);
                for (int n = 0; n < newCount; n++)
                {
                    for (int s = 0; s < ra; s++)
                        rows[s].Add(n < countA ? Unalign(a.Rows[s][colA + n]) : '.');
                    for (int s = 0; s < rb; s++)
                        rows[ra + s].Add(n < countB ? Unalign(b.Rows[s][colB + n]) : '.');
                }

                colA += countA;
                colB += countB;
            }

            foreach (var (type, pla, plb) in path)
            {
                int countA = pla > 0 && pla - 1 > colA ? pla - 1 - colA : 0;
                int countB = plb > 0 && plb - 1 > colB ? plb - 1 - colB : 0;
                Terminals(countA, countB);
                switch (type)
                {
                    case 'M':
                        for (int s = 0; s < ra; s++)
                            rows[s].Add(a.Rows[s][colA]);
                        for (int s = 0; s < rb; s++)
                            rows[ra + s].Add(b.Rows[s][colB]);
                        colA++;
                        colB++;
                        break;
                    case 'D':
                        for (int s = 0; s < ra; s++)
                            rows[s].Add(a.Rows[s][colA]);
                        for (int s = 0; s < rb; s++)
                            rows[ra + s].Add('-');
                        colA++;
                        break;
                    default:
                        for (int s = 0; s < ra; s++)
                            rows[s].Add('-');
                        for (int s = 0; s < rb; s++)
                            rows[ra + s].Add(b.Rows[s][colB]);
                        colB++;
                        break;
                }
            }

            Terminals(a.ColCount - colA, b.ColCount - colB);
            return rows.Select(r => new string(r.ToArray())).ToArray();
        }

        public string MakeTreeMode(IReadOnlyList<string> rows)
        {
            _curWeight = _o.Weight1;
            var msa = PrepareMsa(rows);
            // maketree.cpp never calls SetPPScore, so the command-line -distance2 is never parsed in this mode: the global default PctIdKimura is always used (oracle: R547 / R551 / R554)
            return TreeToFile(TreeFromMsaRows(msa.Ids, msa.Rows, _o.Cluster2, MuscleDistance.PctIdKimura, _o.Root2, out _));
        }

        public string[] RefineWindowsMode(IReadOnlyList<string> rows)
        {
            if (!_nuc)
                throw new NotSupportedException("muscle -refinew on proteins: MUSCLE 3.8.1551 never initialises the score matrix in this mode (SetPPScore is only called for nucleotides) and crashes.");
            var msa = PrepareMsa(rows);
            int n = msa.SeqCount, cols = msa.ColCount, window = _o.RefineWindow;
            if (window < 1)
                throw new ArgumentException("RefineWindow must be positive.");
            int windowCount = (cols + window - 1) / window;
            int to = _o.ToWindow == 0 ? windowCount - 1 : _o.ToWindow;
            var outRows = new List<char>[n];
            for (int s = 0; s < n; s++)
                outRows[s] = new List<char>();
            if (_o.WindowOffset > 0)
            {
                for (int s = 0; s < n; s++)
                    outRows[s].AddRange(msa.Rows[s].AsSpan(0, _o.WindowOffset).ToArray());
            }

            for (int w = _o.FromWindow; w <= to; w++)
            {
                int from = _o.WindowOffset + w * window;
                int last = Math.Min(from + window - 1, cols - 1);
                if (last < from)
                    throw new ArgumentException("A window lies beyond the end of the alignment (the program's behaviour is undefined).");
                var seqs = new string[n];
                for (int s = 0; s < n; s++)
                {
                    seqs[s] = new string(msa.Rows[s].AsSpan(from, last - from + 1).ToArray().Where(c => c != '-').ToArray());
                    if (seqs[s].Length == 0)
                        throw new ArgumentException($"Window {w}: sequence {s} has no residue (MUSCLE's behaviour is undefined).");
                }

                // MUSCLE(v, msaOut): the library entry of the pipeline (nucleotides: SPN score, kmer4_6 distance forced)
                var sub = new MuscleEngine(_o with { SequenceType = _rna ? MuscleSequenceType.Rna : MuscleSequenceType.Dna, Distance1 = MuscleDistance.Kmer4_6, SequenceNames = null, UserTree = null, ClusterOnly = false, ComputeWeights = false, TreeOrder = false }, seqs);
                var (win, _) = sub.Run();
                for (int s = 0; s < n; s++)
                    outRows[s].AddRange(win[s]);
            }

            return outRows.Select(l => new string(l.ToArray())).ToArray();
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // makerootmsab.cpp MakeRootMSABrenner: every sequence is aligned to the profile of the root (the sequence profile cannot be gapped), rows in input order

        private string[] BrennerRootRows(MuscleNode root)
        {
            var rootProf = root.Prof;
            int rootLength = root.Length;
            var rows = new string[_n];
            for (int s = 0; s < _n; s++)
            {
                float w = NormalizedSingleWeightFromId(s);
                var seqProf = ProfileFromRows([(char[])_seq[s].Clone()], [w]);
                foreach (var pp in seqProf)
                {
                    pp.GapOpen = MuscleMinusInfinity;
                    pp.GapClose = MuscleMinusInfinity;
                }

                var path = NwSmall(seqProf, _seq[s].Length, rootProf, rootLength);
                var dummy = new char[][] { Enumerable.Repeat('?', rootLength).ToArray() };
                var combined = RowsFromPath(path, [_seq[s]], dummy);
                rows[s] = new string(combined[0].AsSpan(0, rootLength).ToArray());
            }

            return rows;
        }

        // ProfileFromMSA(msaSeq) of a single sequence: the ClustalW weight of its id normalised to 1
        private float NormalizedSingleWeightFromId(int id)
        {
            var m = new MuscleMsa([id], [[]]);
            SetMsaWeights(m);
            return m.Weights[0];
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // phytofile.cpp Tree::ToFile (rooted), distfunc / clust SaveDF / SaveMSADist

        internal string TreeToFile(MuscleTree t)
        {
            var sb = new StringBuilder();
            void Node(int node)
            {
                bool group = !t.IsLeaf(node) || node == t.Root;
                if (group)
                    sb.Append("(\n");
                if (t.IsLeaf(node))
                {
                    sb.Append(_names[t.LeafId[node]]);
                }
                else
                {
                    Node(t.Left[node]);
                    sb.Append(",\n");
                    Node(t.Right[node]);
                }

                if (group)
                    sb.Append(')');
                if (node != t.Root && !double.IsNaN(t.EdgeLength[node]))
                    sb.Append(':').Append(FormatG(t.EdgeLength[node], 6));
                sb.Append('\n');
            }

            Node(t.Root);   // a one-leaf tree: the leaf is the root, so it is written as a group "(\ns0)\n;\n" (oracle R651)
            sb.Append(";\n");
            return sb.ToString();
        }

        // C printf("%.<precision>g")
        internal static string FormatG(double v, int precision)
        {
            if (v == 0)
                return double.IsNegative(v) ? "-0" : "0";   // printf("%g", -0.0) = "-0" (the oracle's "s2:-0")
            string s = v.ToString("G" + precision, CultureInfo.InvariantCulture);
            if (s.Contains('E'))
            {
                int e = s.IndexOf('E');
                string mant = s[..e], exp = s[(e + 1)..];
                char sign = exp[0] == '-' ? '-' : '+';
                exp = exp.TrimStart('+', '-').TrimStart('0');
                if (exp.Length < 2)
                    exp = exp.PadLeft(2, '0');
                s = mant + "e" + sign + exp;
            }

            return s;
        }

        private string DistanceMatrixText(float[,] d, int n, bool fullMatrix = false, int[]? ids = null)
        {
            var sb = new StringBuilder();
            string Name(int i)
            {
                int id = ids is null ? i : ids[i];
                string nm = id < _names.Length ? _names[id] : "s" + id;
                return (nm.Length > 10 ? nm[..10] : nm).PadLeft(10);
            }

            if (!fullMatrix)
                sb.Append(n).Append('\n');
            for (int i = 0; i < n; i++)
            {
                sb.Append(Name(i)).Append("  ");
                int upTo = fullMatrix ? n : i;
                for (int j = 0; j < upTo; j++)
                    sb.Append("  ").Append(FormatG(d[i, j], 6).PadLeft(9));
                sb.Append('\n');
            }

            return sb.ToString();
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // phyfromfile.cpp Tree::FromFile (rooted, binary): the nodes are numbered as the program creates them (a group appends its two children when its '(' is read), then renumbered to
        // the leaves-first convention of MuscleTree keeping the relative order of the leaves and of the internal nodes

        internal MuscleTree ParseUserTree(string text)
        {
            int pos = 0;
            var left = new List<int> { -1 };
            var right = new List<int> { -1 };
            var parent = new List<int> { -1 };
            var length = new List<double> { 0.0 };
            var label = new List<string?> { null };

            string? Token(out char kind)
            {
                while (pos < text.Length && char.IsWhiteSpace(text[pos]))
                    pos++;
                if (pos >= text.Length)
                    throw new ArgumentException("Tree: unexpected end of the Newick text.");
                char c = text[pos];
                kind = c;
                if (c is '(' or ')' or ':' or ';' or ',')
                {
                    pos++;
                    return c.ToString();
                }

                if (c == '[')
                {
                    while (pos < text.Length && text[pos] != ']')
                        pos++;
                    pos++;
                    return Token(out kind);
                }

                var sb = new StringBuilder();
                if (c is '\'' or '"')
                {
                    pos++;
                    while (pos < text.Length && text[pos] != c)
                        sb.Append(text[pos++]);
                    pos++;
                    kind = 's';
                    return sb.ToString();
                }

                while (pos < text.Length && !char.IsWhiteSpace(text[pos]) && "():;,".IndexOf(text[pos]) < 0)
                    sb.Append(text[pos++]);
                kind = 's';
                return sb.ToString();
            }

            void Group(int node)
            {
                string? tok = Token(out char kind);
                if (kind == 's')
                {
                    label[node] = tok;
                }
                else if (kind == '(')
                {
                    int l = left.Count, r = l + 1;
                    for (int k = 0; k < 2; k++)
                    {
                        left.Add(-1);
                        right.Add(-1);
                        parent.Add(node);
                        length.Add(double.NaN);
                        label.Add(null);
                    }

                    left[node] = l;
                    right[node] = r;
                    Group(l);
                    if (Token(out char k1) != "," || k1 != ',')
                        throw new ArgumentException("Tree: expected ',' (only binary rooted trees are supported).");
                    Group(r);
                    Token(out char k2);
                    if (k2 != ')')
                        throw new ArgumentException("Tree: expected ')' (User tree must be rooted).");
                }
                else
                {
                    throw new ArgumentException($"Tree: expected '(' or a leaf name, got '{tok}'.");
                }

                // optional edge length
                while (pos < text.Length && char.IsWhiteSpace(text[pos]))
                    pos++;
                if (pos < text.Length && text[pos] == ':')
                {
                    pos++;
                    string? len = Token(out char k3);
                    if (k3 != 's')
                        throw new ArgumentException("Tree: expected an edge length.");
                    if (node != 0)
                    {
                        if (!double.TryParse(len, NumberStyles.Float, CultureInfo.InvariantCulture, out double parsed))
                            throw new ArgumentException($"Tree: '{len}' is not a valid edge length.");
                        length[node] = parsed;
                    }
                }
            }

            Group(0);
            Token(out char end);
            if (end != ';')
                throw new ArgumentException("User tree must be rooted");
            int total = left.Count;
            var leaves = Enumerable.Range(0, total).Where(i => left[i] == -1).ToList();
            var internals = Enumerable.Range(0, total).Where(i => left[i] != -1).ToList();
            if (leaves.Count != _n)
                throw new ArgumentException("User tree does not match input sequences");
            var map = new int[total];
            for (int i = 0; i < leaves.Count; i++)
                map[leaves[i]] = i;
            for (int i = 0; i < internals.Count; i++)
                map[internals[i]] = leaves.Count + i;
            var t = new MuscleTree { LeafCount = leaves.Count, Root = map[0] };
            t.Left = new int[total];
            t.Right = new int[total];
            t.Parent = new int[total];
            t.EdgeLength = new double[total];
            t.LeafId = new int[total];
            for (int i = 0; i < total; i++)
            {
                int m = map[i];
                t.Left[m] = left[i] == -1 ? -1 : map[left[i]];
                t.Right[m] = right[i] == -1 ? -1 : map[right[i]];
                t.Parent[m] = parent[i] == -1 ? -1 : map[parent[i]];
                t.EdgeLength[m] = length[i];
                if (left[i] == -1)
                {
                    int id = Array.IndexOf(_names, label[i]);
                    if (id < 0)
                        throw new ArgumentException($"Label {label[i]} in tree does not match sequences");
                    t.LeafId[m] = id;
                }
            }

            if (t.LeafId.Take(leaves.Count).Distinct().Count() != leaves.Count)
                throw new ArgumentException("User tree does not match input sequences");
            return t;
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        // readmx.cpp ReadMx (the program reads HeadingCount - 1 columns per row: with a trailing '*' column that is every residue column, without it the last one is lost)

        internal float[,] ReadMatrixFile(string text)
        {
            var lines = text.Replace("\r", "").Split('\n');
            int li = 0;
            string line = "";
            while (true)
            {
                if (li >= lines.Length)
                    throw new ArgumentException("Premature EOF in matrix file");
                line = lines[li++];
                if (line.Length > 0 && line[0] == '#')
                    continue;
                if (line.Length > 0 && line[0] == ' ')
                    break;
                throw new ArgumentException($"Invalid line in matrix file: '{line}'");
            }

            var heading = line.Where(c => !char.IsWhiteSpace(c)).ToList();
            if (heading.Count > 0 && heading[^1] == '*')
                heading.RemoveAt(heading.Count - 1);
            int headingCount = heading.Count;
            if (headingCount < 20)
                throw new ArgumentException($"Error in matrix file: < 20 headers, line='{line}'");

            var mx = new float[32, 32];
            for (int row = 0; row < headingCount; row++)
            {
                if (li >= lines.Length)
                    throw new ArgumentException("Premature EOF in matrix file");
                line = lines[li++].Trim();
                if (line.Length > 0 && line[0] == '#')
                    continue;
                if (line.Length == 0)
                    continue;
                char c = line[0];
                if (!_isResidue[c])
                    continue;
                int rowLetter = _charToLetterEx[c];
                if (_isWildcard[c] || rowLetter >= 20)
                    continue;
                int p = 1;
                for (int col = 0; col < headingCount - 1; col++)
                {
                    if (p >= line.Length)
                        throw new ArgumentException($"Too few fields in line of matrix file: '{line}'");
                    while (p < line.Length && char.IsWhiteSpace(line[p]))
                        p++;
                    int start = p;
                    while (p < line.Length && !char.IsWhiteSpace(line[p]))
                        p++;
                    float v = AtoF(line[start..p]);
                    char header = heading[col];
                    if (_isResidue[header])
                    {
                        int colLetter = _charToLetterEx[header];
                        if (!_isWildcard[header] && colLetter < 20)
                            mx[rowLetter, colLetter] = v;
                    }

                    p += 1;
                }
            }

            return mx;
        }

        // C atof: the longest numeric prefix
        private static float AtoF(string s)
        {
            int i = 0;
            while (i < s.Length && (char.IsDigit(s[i]) || s[i] is '+' or '-' or '.' or 'e' or 'E'))
                i++;
            for (; i > 0; i--)
            {
                if (double.TryParse(s[..i], NumberStyles.Float, CultureInfo.InvariantCulture, out double d))
                    return (float)d;
            }

            return 0f;
        }
    }
}
