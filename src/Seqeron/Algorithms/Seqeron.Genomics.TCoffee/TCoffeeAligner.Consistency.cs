// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.TCoffee;

using Seqeron.Genomics.Alignment;
using static Seqeron.Genomics.Alignment.SequenceAligner;
using Seqeron.Genomics.ClustalW;
using static Seqeron.Genomics.ClustalW.ClustalWAligner;
using System.Text;
using Seqeron.Genomics.Infrastructure;

public static partial class TCoffeeAligner
{
    #region Multiple Sequence Alignment (Consistency-based — T-Coffee)

    // Consistency-based multiple sequence alignment, the T-Coffee method (Notredame, Higgins &
    // Heringa 2000), added as a FOURTH aligner alongside the star (MultipleAlign), progressive
    // (MultipleAlignProgressive) and iterative (MultipleAlignIterative) methods. All three of those
    // remain byte-for-byte unchanged. This optimises the T-Coffee CONSISTENCY objective — a distinct
    // objective class from the fixed-matrix sum-of-pairs (SP) score the other methods optimise.
    //
    // Algorithm — Notredame, Higgins & Heringa (2000) "T-Coffee: A novel method for fast and accurate
    // multiple sequence alignment", J Mol Biol 302:205-217 (DOI 10.1006/jmbi.2000.4042). Figure 1
    // pipeline: primary library -> extension -> extended library -> progressive alignment.
    //
    //   1. PRIMARY LIBRARY (p.207). For every pair of sequences (i,j) compute pairwise alignments —
    //      a GLOBAL (Needleman-Wunsch, the ClustalW library in the paper) and a LOCAL (Smith-Waterman,
    //      the Lalign library). Each aligned residue pair (residue at position p in Si aligned to the
    //      residue at position q in Sj) "receives a weight equal to percent identity within the
    //      pair-wise alignment it comes from". The global and local libraries are combined by
    //      "signal addition": "If any pair is duplicated between the two libraries, it is merged into
    //      a single entry that has a weight equal to the sum of the two weights." Pairs that never
    //      occur have weight 0.
    //
    //   2. LIBRARY EXTENSION — the triplet consistency transformation (pp.208-209). For each library
    //      pair (Si.p, Sj.q), and for every intermediate sequence Sk, if Si.p aligns to Sk.r (weight
    //      W1) and Sk.r aligns to Sj.q (weight W2), the triplet through Sk contributes min(W1, W2):
    //      "we associate that alignment with a weight equal to the minimum of W1 and W2". The extended
    //      weight is "the sum of all the weights gathered through the examination of all the triplets
    //      involving that pair" PLUS the direct primary weight — the paper's worked example: direct
    //      88, through-C min(77,100)=77, extended = 88 + 77 = 165. Uninformative triplets contribute 0,
    //      so extension never lowers a weight.
    //
    //   3. PROGRESSIVE ALIGNMENT (pp.209-210). A guide tree (the shared
    //      BuildProgressiveGuideTree: UPGMA by default; the paper uses NJ, available as MsaGuideTree.NeighborJoining —
    //      merge order only, same library and objective)
    //      directs a progressive alignment. Profiles are aligned by dynamic programming whose
    //      position-specific column score is the SUM of the EXTENDED-LIBRARY weights over all
    //      cross-profile residue pairs in the two columns (group-vs-group: "the average library scores
    //      in each column of existing alignment are taken" — here summed over the cross pairs, which
    //      is the same objective up to a per-merge constant). No fixed substitution matrix is used and
    //      "gap-opening penalties and gap-extension penalties [are] set to zero" (p.210): a gap column
    //      scores 0. "Once a gap is introduced ... it cannot be shifted later" is enforced as in the
    //      progressive aligner (whole all-gap columns inserted, existing columns never edited).
    //
    // Determinism: no RNG. Primary alignments, the guide tree, the extension sum and the DP tie-breaks
    // are all deterministic.
    //
    // Sources retrieved this session (2026-06-23):
    //   - Notredame C, Higgins DG, Heringa J (2000) J Mol Biol 302:205-217. Full text fetched/read:
    //     https://web.stanford.edu/class/gene211/pdfs/Notredame-Tcoffee.pdf (Figures 1-2, pp.207-211).
    //   - T-Coffee Technical Documentation (min() rule for triplet legs):
    //     https://tcoffee.readthedocs.io/en/latest/tcoffee_technical_documentation.html
    //   - Gotoh O (1982) J Mol Biol 162:705-708 (the DP used in the progressive phase).

    /// <summary>
    /// Performs <b>consistency-based</b> multiple sequence alignment — the T-Coffee method of
    /// Notredame, Higgins &amp; Heringa (2000) — optimising the T-Coffee consistency objective rather
    /// than the fixed-matrix sum-of-pairs score used by
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?)"/>,
    /// <see cref="MultipleAlignProgressive(IEnumerable{DnaSequence}, ScoringMatrix?)"/> and
    /// <see cref="MultipleAlignIterative(IEnumerable{DnaSequence}, ScoringMatrix?, int)"/>
    /// (all of which are unchanged).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Pipeline (Notredame et al. 2000, J Mol Biol 302:205-217):
    /// </para>
    /// <list type="number">
    /// <item>Build a <b>primary library</b> from all pairwise global (Needleman-Wunsch) and local
    /// (Smith-Waterman) alignments; each aligned residue pair is weighted by the pairwise percent
    /// identity, and global+local weights are summed (signal addition).</item>
    /// <item><b>Extend</b> the library by the triplet consistency transformation: each residue pair
    /// accumulates, over every intermediate sequence, the minimum of the two connecting weights,
    /// added to the direct weight.</item>
    /// <item><b>Progressively align</b> along a guide tree (UPGMA by default, <see cref="ConsistencyMsaOptions.GuideTree"/>) with dynamic programming whose
    /// column score is the summed extended-library weight (no substitution matrix, zero gap
    /// penalty).</item>
    /// </list>
    /// </remarks>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">
    /// Scoring matrix used only to compute the pairwise alignments that seed the library
    /// (default: <see cref="SimpleDna"/>). The progressive phase uses the library weights, not this
    /// matrix, as the substitution score.
    /// </param>
    /// <returns>Multiple alignment result (aligned rows, majority consensus, sum-of-pairs score).</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    public static MultipleAlignmentResult MultipleAlignConsistency(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring = null) =>
        MultipleAlignConsistency(sequences, scoring, new ConsistencyMsaOptions());

    /// <summary>
    /// <see cref="MultipleAlignConsistency(IEnumerable{DnaSequence}, ScoringMatrix?)"/> with additive
    /// <see cref="ConsistencyMsaOptions"/> (guide-tree distance source, <see cref="MsaDistanceMode"/>).
    /// </summary>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix seeding the pairwise library (default: <see cref="SimpleDna"/>).</param>
    /// <param name="options">Options; <c>null</c> = defaults.</param>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static MultipleAlignmentResult MultipleAlignConsistency(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring,
        ConsistencyMsaOptions? options)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        return ConsistencyAlignCore(DnaRows(sequences), scoring ?? SimpleDna, options);
    }

    /// <summary>
    /// <see cref="MultipleAlignConsistency(IEnumerable{DnaSequence}, ScoringMatrix?, ConsistencyMsaOptions?)"/> for <b>protein</b> sequences (or any alphabet of the
    /// <see cref="SubstitutionMatrix"/>): the primary library is built from the global and local pairwise alignments under the matrix (linear or
    /// affine, possibly fractional gap scores), weighted by percent identity; the extension, guide tree and progressive stages are those of the
    /// DNA overload (including <see cref="MsaConsistencyEngine.TCoffeeLinked"/>). See
    /// <see cref="MultipleAlign(IEnumerable{ProteinSequence}, MatrixScoring, ConsistencyMsaOptions?)"/> for the input contract.
    /// </summary>
    public static MultipleAlignmentResult MultipleAlignConsistency(
        IEnumerable<ProteinSequence> sequences,
        MatrixScoring scoring,
        ConsistencyMsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(ConsistencyAlignCore(rows, ps, options), ps);
    }

    /// <summary><see cref="MultipleAlignConsistency(IEnumerable{ProteinSequence}, MatrixScoring, ConsistencyMsaOptions?)"/> for raw sequence strings.</summary>
    public static MultipleAlignmentResult MultipleAlignConsistency(
        IEnumerable<string> sequences,
        MatrixScoring scoring,
        ConsistencyMsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(ConsistencyAlignCore(rows, ps, options), ps);
    }

    private static MultipleAlignmentResult ConsistencyAlignCore(string[] rawSeqs, PairScoring effectiveScoring, ConsistencyMsaOptions? options)
    {
        using var cancelScope = AlignmentCancellation.Enter(options?.CancellationToken ?? default);
        options ??= new ConsistencyMsaOptions();
        var gapModel = ResolveGapModel(options);

        var trivial = TrivialMultipleAlignment(rawSeqs);
        if (trivial != null)
            return trivial;

        int k = rawSeqs.Length;

        // Steps 1-2: primary library (global + local, signal-added) then triplet extension.
        if (options.PrimaryLocalAlignments < 1)
            throw new ArgumentOutOfRangeException(nameof(options), options.PrimaryLocalAlignments, "PrimaryLocalAlignments must be at least 1.");
        if (!Enum.IsDefined(options.SequenceWeights))
            throw new ArgumentOutOfRangeException(nameof(options), "Unknown sequence weighting.");
        var primaryLibrary = BuildPrimaryLibrary(rawSeqs, effectiveScoring, gapModel, options.PrimaryLocalAlignments);
        float[]? seqWeights = options.SequenceWeights == MsaSequenceWeights.TCoffee && k > 1 && rawSeqs.All(r => r.Length > 0)
            ? TcTCoffeeSequenceWeights(TcDistanceMatrix(rawSeqs, TcKindOf(rawSeqs, null), new TCoffeeProgramOptions()).Sim, k)
            : null;
        var library = new ExtendedLibrary(primaryLibrary, seqWeights);

        // Step 3a: guide tree (UPGMA by default, neighbour-joining on request) over the ClustalW distances (1 - identities / min(len1, len2)), reusing the progressive pipeline.
        var distance = ComputeGuideDistances(rawSeqs, effectiveScoring, options.Distance, gapModel, options.GuideTree, options.QuickPairwise, options.ClustalWProtein, options.FastDistanceCells);
        ProgressiveGuideNode root = BuildProgressiveGuideTree(k, distance, options.GuideTree);

        if (options.ConsistencyEngine == MsaConsistencyEngine.TCoffeeLinked)
            return BuildTCoffeeLinkedResult(rawSeqs, primaryLibrary, root, effectiveScoring, gapModel, seqWeights);

        // Step 3b: progressive alignment driven by the EXTENDED-LIBRARY weights.
        ConsistencyProfile aligned = AlignConsistencySubtree(root, rawSeqs, library);

        // Reproject rows back to input order.
        var orderedRows = new string[k];
        for (int r = 0; r < aligned.Rows.Count; r++)
            orderedRows[aligned.SequenceIndices[r]] = aligned.Rows[r];

        var alignedList = orderedRows.ToList();
        int maxLen = alignedList.Count == 0 ? 0 : alignedList.Max(s => s.Length);
        string consensus = BuildConsensus(alignedList, maxLen, effectiveScoring.IsMatrix);
        // TotalScore is reported in the same SP currency as the sibling aligners for comparability.
        int spScore = ComputeSumOfPairsScore(alignedList, effectiveScoring, gapModel);

        return new MultipleAlignmentResult(
            AlignedSequences: orderedRows,
            Consensus: consensus,
            TotalScore: spScore);
    }

    // ----- Library representation -------------------------------------------------------------

    // A residue-pair key (Si.posI, Sj.posJ) with i &lt; j, used to store library weights. Stored
    // canonically with the smaller sequence index first so lookups are symmetric.
    private readonly record struct ResiduePairKey(int SeqA, int PosA, int SeqB, int PosB);

    private static ResiduePairKey MakeKey(int s1, int p1, int s2, int p2) =>
        s1 < s2 || (s1 == s2 && p1 <= p2)
            ? new ResiduePairKey(s1, p1, s2, p2)
            : new ResiduePairKey(s2, p2, s1, p1);

    /// <summary>
    /// The T-Coffee <b>extended</b> library (Notredame et al. 2000, pp.207-209): the primary library
    /// (global + local pairwise alignments, weighted by percent identity and combined by signal
    /// addition) plus the triplet consistency transformation, evaluated for <b>any</b> residue pair
    /// on demand (the paper: "Extension will be carried out on each pair of residues of A and B";
    /// the reference implementation <c>residue_pair_extended_list</c> in <c>lib/dp_lib/evaluate.c</c>
    /// likewise scores every pair the dynamic program asks for, not only pairs of the primary library).
    /// </summary>
    private sealed class ExtendedLibrary
    {
        private readonly Dictionary<ResiduePairKey, double> _primary;
        private readonly Dictionary<(int Seq, int Pos), Dictionary<(int Seq, int Pos), double>> _neighbours;
        private readonly Dictionary<ResiduePairKey, double> _cache = new();
        private readonly float[]? _seqWeights;

        public ExtendedLibrary(Dictionary<ResiduePairKey, double> primary, float[]? seqWeights = null)
        {
            _primary = primary;
            _neighbours = BuildNeighbourIndex(primary);
            _seqWeights = seqWeights;
        }

        /// <summary>
        /// Inverts the primary library into a per-residue adjacency: (seq,pos) -&gt; {(otherSeq,otherPos):
        /// weight}. Both directions of each stored pair are emitted so the triplet extension can walk from
        /// either endpoint.
        /// </summary>
        private static Dictionary<(int Seq, int Pos), Dictionary<(int Seq, int Pos), double>>
            BuildNeighbourIndex(Dictionary<ResiduePairKey, double> primary)
        {
            var index = new Dictionary<(int, int), Dictionary<(int, int), double>>();

            void Add(int s1, int p1, int s2, int p2, double w)
            {
                if (!index.TryGetValue((s1, p1), out var inner))
                {
                    inner = new Dictionary<(int, int), double>();
                    index[(s1, p1)] = inner;
                }
                inner[(s2, p2)] = w;
            }

            foreach (var (key, w) in primary)
            {
                Add(key.SeqA, key.PosA, key.SeqB, key.PosB, w);
                Add(key.SeqB, key.PosB, key.SeqA, key.PosA, w);
            }

            return index;
        }

        public double Primary(ResiduePairKey key) => _primary.GetValueOrDefault(key);

        /// <summary>
        /// The column-versus-column scores of two consistency profiles as one matrix: cell (i, j) = Σ over all cross pairs of residues of <see cref="Extended"/>, i.e. exactly what <c>ColumnLibraryScore</c> adds up cell by
        /// cell, but computed from the support: a pair has a non-zero extended weight only when it is in the primary library or both residues have a common partner in a third sequence, so the candidates of a residue
        /// are its partners and the partners of its partners (B09 WP78, R14-16: no dictionary probe for each of the m · n · rows<sub>1</sub> · rows<sub>2</sub> pairs). The weights are integers (percent identities), so the sum is exact in any
        /// order. Null with sequence weights (T-Coffee <c>-seq_weight</c> multiplies by floats: the order of the additions would matter), then the cell-by-cell score is used.
        /// </summary>
        public double[,]? ColumnScoreMatrix(ConsistencyProfile p1, ConsistencyProfile p2)
        {
            if (_seqWeights is not null)
                return null;
            int m = p1.Width, n = p2.Width;
            var score = new double[m, n];
            var columnOf = new Dictionary<int, int[]>();
            for (int r2 = 0; r2 < p2.Rows.Count; r2++)
            {
                var map = p2.PositionMaps[r2];
                int size = 0;
                foreach (int pos in map)
                    size = Math.Max(size, pos + 1);
                var inverse = new int[size];
                Array.Fill(inverse, -1);
                for (int col = 0; col < map.Length; col++)
                {
                    if (map[col] >= 0)
                        inverse[map[col]] = col;
                }

                columnOf[p2.SequenceIndices[r2]] = inverse;
            }

            var candidates = new HashSet<(int Seq, int Pos)>();
            for (int r1 = 0; r1 < p1.Rows.Count; r1++)
            {
                int seqA = p1.SequenceIndices[r1];
                var map = p1.PositionMaps[r1];
                for (int col1 = 0; col1 < map.Length; col1++)
                {
                    int posA = map[col1];
                    if (posA < 0 || !_neighbours.TryGetValue((seqA, posA), out var partners))
                        continue;
                    candidates.Clear();
                    foreach (var partner in partners.Keys)
                    {
                        if (columnOf.ContainsKey(partner.Seq))
                            candidates.Add(partner);
                        if (partner.Seq == seqA || !_neighbours.TryGetValue(partner, out var second))
                            continue;
                        foreach (var q in second.Keys)
                        {
                            if (q.Seq != seqA && q.Seq != partner.Seq && columnOf.ContainsKey(q.Seq))
                                candidates.Add(q);
                        }
                    }

                    foreach (var (seqB, posB) in candidates)
                    {
                        int[] inverse = columnOf[seqB];
                        if (posB >= inverse.Length || inverse[posB] < 0)
                            continue;
                        score[col1, inverse[posB]] += Extended(MakeKey(seqA, posA, seqB, posB));
                    }
                }
            }

            return score;
        }

        /// <summary>
        /// Extended weight of (Si.p, Sj.q) = direct primary weight + Σ over every residue Sk.r of a
        /// third sequence that is in the primary library with both Si.p (W1) and Sj.q (W2) of
        /// min(W1, W2) (worked example of the paper: 88 + min(77, 100) = 165). Zero when the pair has
        /// neither direct nor triplet support.
        /// </summary>
        public double Extended(ResiduePairKey key)
        {
            if (_cache.TryGetValue(key, out double cached))
                return cached;

            double total = _primary.GetValueOrDefault(key);

            if (_neighbours.TryGetValue((key.SeqA, key.PosA), out var partnersOfA) &&
                _neighbours.TryGetValue((key.SeqB, key.PosB), out var partnersOfB))
            {
                // Intersect on (Sk, r); walk the smaller adjacency list.
                var (walk, other) = partnersOfA.Count <= partnersOfB.Count
                    ? (partnersOfA, partnersOfB)
                    : (partnersOfB, partnersOfA);

                foreach (var (partner, w1) in walk)
                {
                    if (partner.Seq == key.SeqA || partner.Seq == key.SeqB) continue; // third sequences only
                    if (other.TryGetValue(partner, out double w2))
                        total += Math.Min(w1, w2);
                }
            }

            // T-Coffee -seq_weight (re_weight_constraint_list, after the extension): the extended weight is multiplied by the smaller of the two sequence weights
            if (_seqWeights is not null)
                total *= Math.Min(_seqWeights[key.SeqA], _seqWeights[key.SeqB]);

            _cache[key] = total;
            return total;
        }
    }

    /// <summary>Primary library: global (NW) + local (SW) pairwise alignments of every sequence pair, signal-added.</summary>
    private static Dictionary<ResiduePairKey, double> BuildPrimaryLibrary(
        string[] seqs, PairScoring scoring, MsaGapModel gapModel = MsaGapModel.Linear, int localAlignments = 1)
    {
        bool affine = gapModel == MsaGapModel.Affine;
        var primary = new Dictionary<ResiduePairKey, double>();
        int k = seqs.Length;

        for (int i = 0; i < k; i++)
        {
            for (int j = i + 1; j < k; j++)
            {
                AddAlignmentToPrimary(primary, i, j, seqs[i], seqs[j],
                    affine ? AffineAlignCore(seqs[i], seqs[j], scoring, AlignmentType.Global) : GlobalAlignCore(seqs[i], seqs[j], scoring));

                if (seqs[i].Length > 0 && seqs[j].Length > 0)
                {
                    if (localAlignments > 1)
                    {
                        // Lalign library of the paper = lalign_id_pair: the Waterman-Eggert top-N non-intersecting local alignments, each weighted by its own percent identity
                        foreach (var local in LocalTopNCore(seqs[i], seqs[j], scoring, affine, localAlignments, 0))
                            AddLocalAlignmentToPrimary(primary, i, j, seqs[i], seqs[j], local);
                    }
                    else
                    {
                        var local = affine ? AffineAlignCore(seqs[i], seqs[j], scoring, AlignmentType.Local) : LocalAlignCore(seqs[i], seqs[j], scoring);
                        AddLocalAlignmentToPrimary(primary, i, j, seqs[i], seqs[j], local);
                    }
                }
            }
        }

        return primary;
    }

    /// <summary>
    /// Adds a global pairwise alignment to the primary library. Every aligned residue pair (both
    /// non-gap) gets the alignment's percent identity weight; duplicate pairs are signal-added.
    /// </summary>
    private static void AddAlignmentToPrimary(
        Dictionary<ResiduePairKey, double> primary,
        int i, int j, string si, string sj, AlignmentResult pa)
    {
        double weight = PercentIdentity(pa.AlignedSequence1, pa.AlignedSequence2);
        if (weight <= 0) return;

        int posI = 0, posJ = 0;
        string a = pa.AlignedSequence1, b = pa.AlignedSequence2;
        for (int c = 0; c < a.Length; c++)
        {
            char ca = a[c], cb = b[c];
            if (ca != '-' && cb != '-')
            {
                AddPairWeight(primary, i, posI, j, posJ, weight);
                posI++; posJ++;
            }
            else if (ca != '-') posI++;
            else if (cb != '-') posJ++;
        }
    }

    /// <summary>
    /// Adds a local pairwise alignment to the primary library. The local alignment covers a segment
    /// of each sequence starting at <c>StartPosition1</c>/<c>StartPosition2</c>; matched residue pairs
    /// in that segment are weighted by the local segment's percent identity and signal-added.
    /// </summary>
    private static void AddLocalAlignmentToPrimary(
        Dictionary<ResiduePairKey, double> primary,
        int i, int j, string si, string sj, AlignmentResult la)
    {
        double weight = PercentIdentity(la.AlignedSequence1, la.AlignedSequence2);
        if (weight <= 0) return;

        int posI = la.StartPosition1, posJ = la.StartPosition2;
        string a = la.AlignedSequence1, b = la.AlignedSequence2;
        for (int c = 0; c < a.Length; c++)
        {
            char ca = a[c], cb = b[c];
            if (ca != '-' && cb != '-')
            {
                AddPairWeight(primary, i, posI, j, posJ, weight);
                posI++; posJ++;
            }
            else if (ca != '-') posI++;
            else if (cb != '-') posJ++;
        }
    }

    private static void AddPairWeight(
        Dictionary<ResiduePairKey, double> lib, int s1, int p1, int s2, int p2, double weight)
    {
        var key = MakeKey(s1, p1, s2, p2);
        lib[key] = lib.GetValueOrDefault(key) + weight;
    }

    /// <summary>
    /// Percent identity of a pairwise alignment = ⌊100 × identical columns / aligned columns⌋, where an
    /// aligned column has a residue in both rows (gap columns are not counted). This is the library
    /// weight of T-Coffee's default <c>-weight=sim</c> (<c>get_seq_sim</c> in
    /// <c>lib/util_lib/aln_convertion_util.c</c>, mode sim1: <c>(sim*MAXID)/pos0</c>, integer division;
    /// Notredame et al. 2000 p.207: "a weight equal to percent identity within the pair-wise
    /// alignment it comes from").
    /// </summary>
    private static double PercentIdentity(string a, string b)
    {
        int aligned = 0;
        int identical = 0;
        for (int c = 0; c < a.Length; c++)
        {
            if (a[c] == '-' || b[c] == '-')
                continue;
            aligned++;
            if (a[c] == b[c])
                identical++;
        }

        return aligned == 0 ? 0 : (int)(PercentIdentityScale * identical) / aligned;
    }

    // Weight is expressed as a percentage (0-100), matching the paper's integer "Prim Weight" values.
    private const double PercentIdentityScale = 100.0;

    // ----- Consistency profile + progressive DP ----------------------------------------------

    /// <summary>
    /// A profile for consistency alignment: like <see cref="Profile"/> but every cell also records the
    /// original residue position in its source sequence (or -1 for a gap), so the library weight for a
    /// residue pair can be looked up during the position-specific DP.
    /// </summary>
    private sealed class ConsistencyProfile
    {
        public List<string> Rows { get; } = new();
        public List<int> SequenceIndices { get; } = new();
        // PositionMaps[r][col] = original residue index in sequence SequenceIndices[r], or -1 for gap.
        public List<int[]> PositionMaps { get; } = new();
        public int Width => Rows.Count == 0 ? 0 : Rows[0].Length;
    }

    private static ConsistencyProfile AlignConsistencySubtree(
        ProgressiveGuideNode node, string[] rawSeqs, ExtendedLibrary library)
    {
        if (node.IsLeaf)
        {
            int idx = node.LeafIndex;
            string seq = rawSeqs[idx];
            var leaf = new ConsistencyProfile();
            leaf.Rows.Add(seq);
            leaf.SequenceIndices.Add(idx);
            var map = new int[seq.Length];
            for (int p = 0; p < seq.Length; p++) map[p] = p;
            leaf.PositionMaps.Add(map);
            return leaf;
        }

        var left = AlignConsistencySubtree(node.Left!, rawSeqs, library);
        var right = AlignConsistencySubtree(node.Right!, rawSeqs, library);
        return AlignConsistencyProfiles(left, right, library);
    }

    /// <summary>
    /// Aligns two consistency profiles with the Needleman-Wunsch recurrence over columns, scoring a
    /// column pair by the SUM of extended-library weights over all cross-profile residue pairs (gap
    /// columns score 0 — zero gap penalty, Notredame et al. 2000 p.210). Existing columns are never
    /// edited; gaps are inserted as whole all-gap columns ("once a gap, always a gap").
    /// </summary>
    private static ConsistencyProfile AlignConsistencyProfiles(
        ConsistencyProfile p1, ConsistencyProfile p2, ExtendedLibrary library)
    {
        int m = p1.Width;
        int n = p2.Width;
        RequireProfileMatrixFits(m, n);

        var f = new double[m + 1, n + 1];
        // First row/column: aligning a column against a gap scores 0 (zero gap penalty).
        // (f[i,0] = f[0,j] = 0 for all i,j.)
        var columnScores = library.ColumnScoreMatrix(p1, p2); // null: the cell-by-cell score (sequence weights)

        for (int i = 1; i <= m; i++)
        {
            AlignmentCancellation.ThrowIfRequested();
            for (int j = 1; j <= n; j++)
            {
                double diag = f[i - 1, j - 1] + (columnScores is null ? ColumnLibraryScore(p1, i - 1, p2, j - 1, library) : columnScores[i - 1, j - 1]);
                double up = f[i - 1, j];     // gap in p2: 0
                double left = f[i, j - 1];   // gap in p1: 0
                f[i, j] = Math.Max(diag, Math.Max(up, left));
            }
        }

        var cols1 = new List<int>();
        var cols2 = new List<int>();
        int ii = m, jj = n;
        while (ii > 0 || jj > 0)
        {
            if (ii > 0 && jj > 0)
            {
                double diag = f[ii - 1, jj - 1] + (columnScores is null ? ColumnLibraryScore(p1, ii - 1, p2, jj - 1, library) : columnScores[ii - 1, jj - 1]);
                if (AreClose(f[ii, jj], diag))
                {
                    cols1.Add(ii - 1);
                    cols2.Add(jj - 1);
                    ii--; jj--;
                    continue;
                }
            }

            if (ii > 0 && (jj == 0 || AreClose(f[ii, jj], f[ii - 1, jj])))
            {
                cols1.Add(ii - 1);
                cols2.Add(-1);
                ii--;
                continue;
            }

            cols1.Add(-1);
            cols2.Add(jj - 1);
            jj--;
        }

        cols1.Reverse();
        cols2.Reverse();

        var merged = new ConsistencyProfile();
        merged.SequenceIndices.AddRange(p1.SequenceIndices);
        merged.SequenceIndices.AddRange(p2.SequenceIndices);
        AppendConsistencyRows(merged, p1, cols1);
        AppendConsistencyRows(merged, p2, cols2);
        return merged;
    }

    private static void AppendConsistencyRows(
        ConsistencyProfile dest, ConsistencyProfile source, List<int> columnMap)
    {
        for (int r = 0; r < source.Rows.Count; r++)
        {
            string row = source.Rows[r];
            int[] srcMap = source.PositionMaps[r];
            var sb = new StringBuilder(columnMap.Count);
            var newMap = new int[columnMap.Count];
            for (int c = 0; c < columnMap.Count; c++)
            {
                int col = columnMap[c];
                if (col < 0)
                {
                    sb.Append('-');
                    newMap[c] = -1;
                }
                else
                {
                    sb.Append(row[col]);
                    newMap[c] = srcMap[col];
                }
            }
            dest.Rows.Add(sb.ToString());
            dest.PositionMaps.Add(newMap);
        }
    }

    /// <summary>
    /// Column-vs-column score = sum over all cross-profile residue pairs of their extended-library
    /// weight (0 when a residue is a gap or the pair is absent from the library). This realises the
    /// position-specific library scoring of Notredame et al. (2000): the DP maximises the summed
    /// consistency weight rather than a fixed substitution matrix.
    /// </summary>
    private static double ColumnLibraryScore(
        ConsistencyProfile p1, int col1, ConsistencyProfile p2, int col2,
        ExtendedLibrary library)
    {
        double total = 0;
        for (int r1 = 0; r1 < p1.Rows.Count; r1++)
        {
            int posA = p1.PositionMaps[r1][col1];
            if (posA < 0) continue; // gap contributes 0
            int seqA = p1.SequenceIndices[r1];
            for (int r2 = 0; r2 < p2.Rows.Count; r2++)
            {
                int posB = p2.PositionMaps[r2][col2];
                if (posB < 0) continue;
                int seqB = p2.SequenceIndices[r2];
                total += library.Extended(MakeKey(seqA, posA, seqB, posB));
            }
        }
        return total;
    }

    /// <summary>
    /// Test/diagnostic accessor: builds the T-Coffee extended library for the given sequences and
    /// returns the extended weight of the residue pair (Si.posI, Sj.posJ), and that pair's primary
    /// (pre-extension) weight, so the consistency transformation can be verified directly
    /// (Notredame et al. 2000, p.209).
    /// </summary>
    internal static (double Primary, double Extended) GetLibraryWeights(
        string[] seqs, int seqI, int posI, int seqJ, int posJ, ScoringMatrix? scoring = null)
    {
        var effectiveScoring = scoring ?? SimpleDna;
        var library = new ExtendedLibrary(BuildPrimaryLibrary(seqs, effectiveScoring));
        var key = MakeKey(seqI, posI, seqJ, posJ);
        return (library.Primary(key), library.Extended(key));
    }

    /// <summary>
    /// Test/diagnostic accessor (audit A2-16): the primary library of <see cref="MultipleAlignConsistency(IEnumerable{string}, MatrixScoring, ConsistencyMsaOptions?)"/> under <paramref name="options"/>
    /// (<see cref="ConsistencyMsaOptions.PrimaryLocalAlignments"/>) as (seqA, posA, seqB, posB) -> weight, and the sequence weights of <see cref="ConsistencyMsaOptions.SequenceWeights"/> (null = none).
    /// </summary>
    internal static (Dictionary<(int SeqA, int PosA, int SeqB, int PosB), double> Primary, float[]? SeqWeights) GetPrimaryLibrary(string[] seqs, ScoringMatrix? scoring, ConsistencyMsaOptions options)
    {
        var lib = BuildPrimaryLibrary(seqs, scoring ?? SimpleDna, ResolveGapModel(options), options.PrimaryLocalAlignments);
        float[]? w = options.SequenceWeights == MsaSequenceWeights.TCoffee && seqs.Length > 1
            ? TcTCoffeeSequenceWeights(TcDistanceMatrix(seqs, TcKindOf(seqs, null), new TCoffeeProgramOptions()).Sim, seqs.Length)
            : null;
        return (lib.ToDictionary(kv => (kv.Key.SeqA, kv.Key.PosA, kv.Key.SeqB, kv.Key.PosB), kv => kv.Value), w);
    }

    /// <summary>Extended weight of a residue pair under <paramref name="options"/> (diagnostic; see <see cref="GetPrimaryLibrary"/>).</summary>
    internal static double GetExtendedWeight(string[] seqs, int seqI, int posI, int seqJ, int posJ, ScoringMatrix? scoring, ConsistencyMsaOptions options)
    {
        var (_, w) = GetPrimaryLibrary(seqs, scoring, options);
        var library = new ExtendedLibrary(BuildPrimaryLibrary(seqs, scoring ?? SimpleDna, ResolveGapModel(options), options.PrimaryLocalAlignments), w);
        return library.Extended(MakeKey(seqI, posI, seqJ, posJ));
    }

    #endregion
}
