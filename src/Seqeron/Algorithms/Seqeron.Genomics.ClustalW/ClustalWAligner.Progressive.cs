// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.ClustalW;

using Seqeron.Genomics.Alignment;
using static Seqeron.Genomics.Alignment.SequenceAligner;
using System.Text;
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Phylogenetics;

public static partial class ClustalWAligner
{
    #region Multiple Sequence Alignment (Guide-Tree Progressive / Feng-Doolittle)

    // Progressive (guide-tree) multiple sequence alignment, the Feng-Doolittle /
    // Clustal-style method, added as a SECOND aligner alongside the star MSA above.
    // The star MSA (MultipleAlign) is unchanged.
    //
    // Algorithm (Feng & Doolittle 1987; Wikipedia "Multiple sequence alignment"):
    //   1. Compute all C(k,2) pairwise global (Needleman-Wunsch) alignments and derive a
    //      distance matrix. Distance = 1 - identities / min(len1, len2): the ClustalW pairwise
    //      distance (Thompson, Higgins & Gibson 1994; pairalign.c tmat = (100 - mm_score) / 100,
    //      mm_score = 100 * count / MIN(len1, len2)), see PairwiseIdentityDistance. (Feng-Doolittle
    //      converts pairwise alignment scores to distances; identity-based distance is the
    //      standard conversion used by Clustal-style tools.)
    //   2. Build a guide tree over that distance matrix — "an efficient clustering
    //      method such as neighbor-joining or unweighted pair group method with arithmetic
    //      mean (UPGMA)" (Wikipedia "Multiple sequence alignment"): UPGMA (the classic
    //      Feng-Doolittle choice, default) or neighbour-joining (ClustalW's default), both built by the
    //      canonical PhylogeneticAnalyzer builders (ProgressiveMsaOptions.GuideTree; SequenceAligner.GuideTree.cs).
    //   3. Progressively align along the guide tree from the tips: sequence-sequence, then
    //      sequence-profile, then profile-profile, using the Needleman-Wunsch recurrence over
    //      *columns* (a profile is the set of already-aligned rows). Column-vs-column score is
    //      the average match/mismatch over all cross-profile residue pairs (sum-of-pairs profile
    //      scoring). The "once a gap, always a gap" rule (Feng & Doolittle 1987: "Once a gap is
    //      introduced ... it is preserved within all subsequent fusions") is enforced: gaps
    //      already present in a profile are never removed; merging two profiles only inserts
    //      whole new gap columns, never edits existing columns.
    //
    // Sources retrieved this session:
    //   - Feng & Doolittle (1987) "Progressive sequence alignment as a prerequisite to correct
    //     phylogenetic trees", J Mol Evol 25:351-360 (PubMed 3118049) — progressive method,
    //     guide tree, "once a gap, always a gap".
    //   - Wikipedia "Multiple sequence alignment", § Progressive alignment construction
    //     (https://en.wikipedia.org/wiki/Multiple_sequence_alignment) — guide tree built by
    //     NJ/UPGMA; sequences added sequentially along the guide tree.
    //   - Wikipedia UPGMA / Sokal & Michener (1958) — cluster averaging formula.

    /// <summary>
    /// Performs guide-tree <b>progressive</b> multiple sequence alignment (the Feng-Doolittle /
    /// Clustal-style method), as an additive alternative to the star-alignment
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?)"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Pipeline (Feng &amp; Doolittle 1987; Wikipedia "Multiple sequence alignment"):
    /// </para>
    /// <list type="number">
    /// <item>All pairwise Needleman-Wunsch alignments give a distance matrix using the ClustalW distance
    /// d = 1 − (identical aligned residue pairs) / (length of the shorter sequence)
    /// (<c>pairalign.c</c>: <c>tmat = (100 − 100·count / MIN(len1, len2)) / 100</c>).</item>
    /// <item>A guide tree is built over that matrix: UPGMA (the classic Feng-Doolittle choice, default) or ClustalW's
    /// neighbour-joining tree (<see cref="ProgressiveMsaOptions.GuideTree"/>).</item>
    /// <item>Profiles are aligned progressively from the tips with the NW recurrence over
    /// columns and sum-of-pairs profile scoring, enforcing "once a gap, always a gap".</item>
    /// </list>
    /// <para>
    /// Linear gaps by default; <see cref="MsaGapModel.Affine"/> (<see cref="ProgressiveMsaOptions.GapModel"/>) switches to a Gotoh three-state
    /// profile DP with occupancy-scaled <see cref="ScoringMatrix.GapOpen"/>/<see cref="ScoringMatrix.GapExtend"/> and the natural affine
    /// sum-of-pairs <c>TotalScore</c> (see <see cref="SumOfPairsScore"/>); exact for two sequences, heuristic for more.
    /// </para>
    /// </remarks>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <returns>Multiple alignment result (aligned rows, majority consensus, sum-of-pairs score).</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    public static MultipleAlignmentResult MultipleAlignProgressive(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring = null) =>
        MultipleAlignProgressive(sequences, scoring, new ProgressiveMsaOptions());

    /// <summary>
    /// <see cref="MultipleAlignProgressive(IEnumerable{DnaSequence}, ScoringMatrix?)"/> with additive
    /// <see cref="ProgressiveMsaOptions"/> (guide-tree distance source: the repository's linear NW alignment or ClustalW's own
    /// pairwise stage, <see cref="MsaDistanceMode"/>).
    /// </summary>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <param name="options">Options; <c>null</c> = defaults.</param>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static MultipleAlignmentResult MultipleAlignProgressive(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring,
        ProgressiveMsaOptions? options)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        return ProgressiveAlignCore(DnaRows(sequences), scoring ?? SimpleDna, options);
    }

    /// <summary>
    /// <see cref="MultipleAlignProgressive(IEnumerable{DnaSequence}, ScoringMatrix?, ProgressiveMsaOptions?)"/> for <b>protein</b> sequences (or any alphabet of the
    /// <see cref="SubstitutionMatrix"/>): the same Feng-Doolittle pipeline — identity distances of the optimal pairwise alignments under the
    /// matrix, the configured guide tree, and profile-profile alignment with the column score Σ matrix(a, b) averaged over the cross pairs
    /// (linear gaps) or the occupancy-scaled Gotoh program (<see cref="MsaGapModel.Affine"/>), real-valued (possibly fractional) gap scores.
    /// See <see cref="MultipleAlign(IEnumerable{ProteinSequence}, MatrixScoring, ProgressiveMsaOptions?)"/> for the input contract.
    /// </summary>
    public static MultipleAlignmentResult MultipleAlignProgressive(
        IEnumerable<ProteinSequence> sequences,
        MatrixScoring scoring,
        ProgressiveMsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(ProgressiveAlignCore(rows, ps, options), ps);
    }

    /// <summary><see cref="MultipleAlignProgressive(IEnumerable{ProteinSequence}, MatrixScoring, ProgressiveMsaOptions?)"/> for raw sequence strings.</summary>
    public static MultipleAlignmentResult MultipleAlignProgressive(
        IEnumerable<string> sequences,
        MatrixScoring scoring,
        ProgressiveMsaOptions? options = null)
    {
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(ProgressiveAlignCore(rows, ps, options), ps);
    }

    private static MultipleAlignmentResult ProgressiveAlignCore(string[] rawSeqs, PairScoring effectiveScoring, ProgressiveMsaOptions? options)
    {
        using var cancelScope = AlignmentCancellation.Enter(options?.CancellationToken ?? default);
        options ??= new ProgressiveMsaOptions();
        var gapModel = ResolveGapModel(options);

        var trivial = TrivialMultipleAlignment(rawSeqs);
        if (trivial != null)
            return trivial;

        int k = rawSeqs.Length;

        // Step 1: pairwise alignments -> distance matrix (ClustalW: 1 - identities / min(len1, len2)).
        var distance = ComputeGuideDistances(rawSeqs, effectiveScoring, options.Distance, gapModel, options.GuideTree, options.QuickPairwise, options.ClustalWProtein, options.FastDistanceCells);

        // Step 2: guide tree (UPGMA by default, neighbour-joining on request).
        ProgressiveGuideNode root = BuildProgressiveGuideTree(k, distance, options.GuideTree);

        // Step 3: progressive alignment along the tree (post-order: tips first).
        Profile aligned = AlignProfileSubtree(root, rawSeqs, effectiveScoring, gapModel);

        // Reproject rows back into input order (the tree visited leaves in tree order).
        var orderedRows = new string[k];
        for (int r = 0; r < aligned.Rows.Count; r++)
            orderedRows[aligned.SequenceIndices[r]] = aligned.Rows[r];

        var alignedList = orderedRows.ToList();
        int maxLen = alignedList.Count == 0 ? 0 : alignedList.Max(s => s.Length);
        string consensus = BuildConsensus(alignedList, maxLen, effectiveScoring.IsMatrix);
        int spScore = ComputeSumOfPairsScore(alignedList, effectiveScoring, gapModel);

        return new MultipleAlignmentResult(
            AlignedSequences: orderedRows,
            Consensus: consensus,
            TotalScore: spScore);
    }

    /// <summary>
    /// Pairwise distance for the guide tree, the ClustalW pairwise distance (Thompson, Higgins &amp;
    /// Gibson 1994): 1 − (identical aligned residue pairs of the optimal pairwise alignment) /
    /// (length of the shorter sequence). This is <c>tmat = (100 − mm_score) / 100</c> with
    /// <c>mm_score = 100·count / MIN(len1, len2)</c> in <c>pairalign.c</c> of ClustalW 1.83 (the
    /// reference alignment there uses affine gaps; this method uses the repository's linear
    /// Needleman-Wunsch). Two empty sequences (no comparable residue) are treated as distance 0, an
    /// empty versus a non-empty sequence as distance 1.
    /// </summary>
    private static double PairwiseIdentityDistance(string a, string b, PairScoring scoring, MsaGapModel gapModel = MsaGapModel.Linear)
    {
        if (a.Length == 0 && b.Length == 0)
            return 0.0;
        if (a.Length == 0 || b.Length == 0)
            return 1.0; // no shared residues possible

        if (gapModel != MsaGapModel.Affine)
            return 1.0 - (double)LinearIdentityCount(a, b, scoring) / Math.Min(a.Length, b.Length); // the identities of the traceback path without a matrix, B09 WP78

        var pa = AffineAlignCore(a, b, scoring, AlignmentType.Global);

        int identical = 0;
        for (int i = 0; i < pa.AlignedSequence1.Length; i++)
        {
            char c1 = pa.AlignedSequence1[i];
            if (c1 != '-' && c1 == pa.AlignedSequence2[i])
                identical++;
        }

        return 1.0 - (double)identical / Math.Min(a.Length, b.Length);
    }

    /// <summary>
    /// Symmetric guide-tree distance matrix of the progressive / iterative / consistency aligners:
    /// <see cref="MsaDistanceMode.Linear"/> = the repository's own global alignment under <paramref name="scoring"/> (linear NW, or
    /// the affine Gotoh alignment when <paramref name="gapModel"/> is <see cref="MsaGapModel.Affine"/>; <see cref="PairwiseIdentityDistance"/>); <see cref="MsaDistanceMode.ClustalWAffine"/> = ClustalW's own pairwise stage
    /// (<see cref="ClustalWPairwiseAlign(string, string, double, double)"/>, independent of <paramref name="scoring"/>);
    /// <see cref="MsaDistanceMode.ClustalWQuick"/> = ClustalW's fast k-tuple stage (<see cref="ClustalWQuickPairwiseAlign(string, string, ClustalWQuickOptions?, int?)"/>).
    /// </summary>
    public static double[,] ComputeGuideDistances(string[] rawSeqs, PairScoring scoring, MsaDistanceMode mode, MsaGapModel gapModel, MsaGuideTree guideTree = MsaGuideTree.Upgma, ClustalWQuickOptions? quickOptions = null, ClustalWProteinPairwiseOptions? proteinOptions = null, long fastDistanceCells = long.MaxValue)
    {
        if (!Enum.IsDefined(mode))
            throw new ArgumentOutOfRangeException(nameof(mode), mode, "Unknown guide-tree distance mode.");

        // Protein matrices select ClustalW's PROTEIN pairwise stages (GONNET / BLOSUM30 / PAM350 / identity, gap 10 / 0.1, k-tuple 1), nucleotide matrices
        // and the scalar scoring the DNA stages.
        bool protein = scoring.Protein;
        var proteinParameters = proteinOptions ?? new ClustalWProteinPairwiseOptions();
        var quick = protein ? proteinParameters.Quick : quickOptions ?? new ClustalWQuickOptions();
        int longest = 0;
        foreach (var s in rawSeqs)
            longest = Math.Max(longest, s.Length);

        int k = rawSeqs.Length;
        var distance = new double[k, k];

        // B09 WP78 (R14-16): above the cell limit of ProgressiveMsaOptions.FastDistanceCells the Linear distance of a pair is the k-tuple estimate (the guide tree needs a ranking of the pairs, not their optimal alignments)
        ulong[]?[]? tupleCodes = null;
        int tupleLength = scoring.Protein ? KTupleProteinLength : KTupleDnaLength;
        if (mode == MsaDistanceMode.Linear && fastDistanceCells != long.MaxValue && UseKTupleDistances(rawSeqs, fastDistanceCells))
            tupleCodes = rawSeqs.Select(s => KTupleCodes(s, tupleLength)).ToArray();

        for (int i = 0; i < k; i++)
        {
            for (int j = i + 1; j < k; j++)
            {
                double d;
                if (tupleCodes?[i] is { } codesI && tupleCodes[j] is { } codesJ)
                    d = KTupleDistance(codesI, codesJ, rawSeqs[i].Length, rawSeqs[j].Length, tupleLength);
                else if (mode == MsaDistanceMode.ClustalWQuick)
                    d = protein
                        ? ClustalWQuickPairwiseAlign(rawSeqs[i], rawSeqs[j], ClustalWSequenceType.Protein, quick, longest).Distance
                        : ClustalWQuickPairwiseDistance(rawSeqs[i], rawSeqs[j], quick, longest);
                else if (mode != MsaDistanceMode.ClustalWAffine)
                    d = PairwiseIdentityDistance(rawSeqs[i], rawSeqs[j], scoring, gapModel);
                else if (protein)
                    d = guideTree == MsaGuideTree.NeighborJoining
                        ? ClustalWProteinSinglePrecisionDistance(rawSeqs[i], rawSeqs[j], proteinParameters)
                        : ClustalWPairwiseAlign(rawSeqs[i], rawSeqs[j], proteinParameters).Distance;
                else if (guideTree == MsaGuideTree.NeighborJoining)
                    d = ClustalWSinglePrecisionDistance(rawSeqs[i], rawSeqs[j]);
                else
                    d = ClustalWPairwiseDistance(rawSeqs[i], rawSeqs[j]);
                distance[i, j] = d;
                distance[j, i] = d;
            }
        }

        return distance;
    }

    /// <summary>
    /// A profile: a block of aligned rows of equal length, each tagged with the input-sequence
    /// index it represents (so the final output can be reordered back to input order).
    /// </summary>
    private sealed class Profile
    {
        public List<string> Rows { get; } = new();
        public List<int> SequenceIndices { get; } = new();
        public int Width => Rows.Count == 0 ? 0 : Rows[0].Length;
    }

    /// <summary>
    /// Post-order traversal of the guide tree: align each subtree's profile, combining child
    /// profiles by profile-profile NW alignment. Leaves start as a single-row profile.
    /// </summary>
    private static Profile AlignProfileSubtree(
        ProgressiveGuideNode node, string[] rawSeqs, PairScoring scoring, MsaGapModel gapModel)
    {
        if (node.IsLeaf)
        {
            var leaf = new Profile();
            leaf.Rows.Add(rawSeqs[node.LeafIndex]);
            leaf.SequenceIndices.Add(node.LeafIndex);
            return leaf;
        }

        var left = AlignProfileSubtree(node.Left!, rawSeqs, scoring, gapModel);
        var right = AlignProfileSubtree(node.Right!, rawSeqs, scoring, gapModel);
        return gapModel == MsaGapModel.Affine
            ? AlignProfilesAffine(left, right, scoring)
            : AlignProfiles(left, right, scoring);
    }

    /// <summary>
    /// Largest profile-vs-profile dynamic-programming matrix (<c>(m + 1)·(n + 1)</c> doubles, 8 bytes per cell = 400 MB) of the
    /// progressive / iterative / consistency multiple aligners. The exact center-star <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?, ProgressiveMsaOptions?)"/> needs no such
    /// matrix (its pairwise stage is linear-space); a profile-profile DP has no published linear-space variant that keeps the
    /// column-pair scores of these methods, so longer inputs get a clear error instead of an <see cref="OutOfMemoryException"/>.
    /// </summary>
    public const long ProfileMatrixMaxCells = 50_000_000;

    public static void RequireProfileMatrixFits(int m, int n)
    {
        long cells = (long)(m + 1) * (n + 1);
        if (cells > ProfileMatrixMaxCells)
            throw new ArgumentException(
                string.Create(System.Globalization.CultureInfo.InvariantCulture, $"The profile-profile dynamic program needs O(m·n) memory: {m} x {n} alignment columns = {cells:N0} cells exceeds the limit of {ProfileMatrixMaxCells:N0}. ") +
                "Use MultipleAlign (exact center-star, linear-space pairwise stage) for long sequences.");
    }

    /// <summary>
    /// Aligns two profiles with the Needleman-Wunsch recurrence over columns. The column-vs-column
    /// substitution score is the average match/mismatch over all cross-profile residue pairs
    /// (sum-of-pairs profile scoring); a gap column scores GapExtend per residue averaged.
    /// Existing columns are never edited — gaps are inserted as whole new all-gap columns into one
    /// profile, enforcing "once a gap, always a gap" (Feng &amp; Doolittle 1987).
    /// </summary>
    private static Profile AlignProfiles(Profile p1, Profile p2, PairScoring scoring)
    {
        int m = p1.Width;
        int n = p2.Width;
        RequireProfileMatrixFits(m, n);

        // The column scores from the compositions of the columns (ProfileColumnTable: the same doubles as the row-by-row sums of ProfileColumnScore / ProfileColumnVsGap).
        var symbolIndex = new Dictionary<char, int>();
        foreach (var profile in new[] { p1, p2 })
            foreach (var row in profile.Rows)
                foreach (char c in row)
                    if (c != '-' && !symbolIndex.ContainsKey(c))
                        symbolIndex[c] = symbolIndex.Count;
        var symbols = new char[symbolIndex.Count];
        foreach (var (symbol, index) in symbolIndex)
            symbols[index] = symbol;
        var substitution = new int[symbols.Length][];
        for (int x = 0; x < symbols.Length; x++)
        {
            substitution[x] = new int[symbols.Length];
            for (int y = 0; y < symbols.Length; y++)
                substitution[x][y] = scoring.Sub(symbols[x], symbols[y]);
        }

        var t1 = new ProfileColumnTable(p1, symbolIndex);
        var t2 = new ProfileColumnTable(p2, symbolIndex);
        t1.Project(substitution);
        t2.Project(substitution);
        int gapExtend = scoring.GapExtend;
        int rows1 = p1.Rows.Count;
        int rows2 = p2.Rows.Count;
        var upCost = new double[m + 1];
        var leftCost = new double[n + 1];
        for (int i = 1; i <= m; i++)
            upCost[i] = ProfileColumnVsGapFast(t1, i - 1, rows2, gapExtend);
        for (int j = 1; j <= n; j++)
            leftCost[j] = ProfileColumnVsGapFast(t2, j - 1, rows1, gapExtend);

        // F(i,j): best score aligning first i columns of p1 with first j columns of p2.
        var f = new double[m + 1, n + 1];
        for (int i = 1; i <= m; i++)
            f[i, 0] = f[i - 1, 0] + upCost[i];
        for (int j = 1; j <= n; j++)
            f[0, j] = f[0, j - 1] + leftCost[j];

        for (int i = 1; i <= m; i++)
        {
            AlignmentCancellation.ThrowIfRequested();
            for (int j = 1; j <= n; j++)
            {
                double diag = f[i - 1, j - 1] + ProfileColumnScoreFast(t1, i - 1, t2, j - 1, gapExtend);
                double up = f[i - 1, j] + upCost[i];
                double left = f[i, j - 1] + leftCost[j];
                f[i, j] = Math.Max(diag, Math.Max(up, left));
            }
        }

        // Traceback, building column lists for both profiles.
        var cols1 = new List<int>(); // column index in p1, or -1 for a gap column
        var cols2 = new List<int>();
        int ii = m, jj = n;
        while (ii > 0 || jj > 0)
        {
            if (ii > 0 && jj > 0)
            {
                double diag = f[ii - 1, jj - 1] + ProfileColumnScoreFast(t1, ii - 1, t2, jj - 1, gapExtend);
                if (AreClose(f[ii, jj], diag))
                {
                    cols1.Add(ii - 1);
                    cols2.Add(jj - 1);
                    ii--; jj--;
                    continue;
                }
            }

            if (ii > 0)
            {
                double up = f[ii - 1, jj] + upCost[ii];
                if (jj == 0 || AreClose(f[ii, jj], up))
                {
                    cols1.Add(ii - 1);
                    cols2.Add(-1);
                    ii--;
                    continue;
                }
            }

            cols1.Add(-1);
            cols2.Add(jj - 1);
            jj--;
        }

        cols1.Reverse();
        cols2.Reverse();

        var merged = new Profile();
        merged.SequenceIndices.AddRange(p1.SequenceIndices);
        merged.SequenceIndices.AddRange(p2.SequenceIndices);
        AppendReprojectedRows(merged, p1, cols1);
        AppendReprojectedRows(merged, p2, cols2);
        return merged;
    }

    /// <summary>
    /// Re-emits each row of <paramref name="source"/> according to <paramref name="columnMap"/>:
    /// a non-negative entry copies that existing column (gaps in it preserved), a -1 inserts a
    /// fresh gap. Existing columns are copied verbatim ("once a gap, always a gap").
    /// </summary>
    private static void AppendReprojectedRows(Profile dest, Profile source, List<int> columnMap)
    {
        foreach (var row in source.Rows)
        {
            var sb = new StringBuilder(columnMap.Count);
            foreach (int col in columnMap)
                sb.Append(col < 0 ? '-' : row[col]);
            dest.Rows.Add(sb.ToString());
        }
    }

    private const double ProfileScoreEpsilon = 1e-9;

    public static bool AreClose(double x, double y) => Math.Abs(x - y) <= ProfileScoreEpsilon;

    #endregion

    #region Multiple Sequence Alignment (Iterative Refinement — MUSCLE tree-dependent restricted partitioning)

    // Iterative refinement of a progressive MSA, added as a THIRD aligner alongside the star
    // (MultipleAlign) and progressive (MultipleAlignProgressive) methods. Both of those remain
    // byte-for-byte unchanged. This removes the single-pass "once a gap, always a gap" limitation
    // of the progressive seed: an early gap-placement error CAN now be corrected, because the
    // alignment is repeatedly re-split and re-aligned and a strictly better arrangement replaces it.
    //
    // Algorithm — MUSCLE Stage 3, "tree-dependent restricted partitioning" (Edgar 2004):
    //   Quoting Edgar (2004), Nucleic Acids Res 32(5):1792-1797, §"Stage 3, Refinement":
    //   3.1 "An edge is chosen from TREE2 (edges are visited in order of decreasing distance from
    //        the root). TREE2 is divided into two subtrees by deleting the edge."
    //   3.2/3.3 "The profile of the multiple alignment in each subtree is computed. A new multiple
    //        alignment is produced by re-aligning the two profiles."
    //   3.4 "If the SP score is improved, the new alignment is kept, otherwise it is discarded."
    //   "Steps 3.1-3.4 are repeated until convergence or until a user-defined limit is reached."
    //
    //   The same scheme also realises the Barton-Sternberg (1987) idea of returning to an existing
    //   alignment and re-aligning a part of it against the rest; here the partition is the
    //   guide-tree edge rather than a single removed sequence, which is exactly Edgar's restricted
    //   partitioning. Re-aligning the two sub-profiles uses the EXISTING profile-profile NW
    //   (AlignProfiles) and is accepted only on a non-decreasing sum-of-pairs (SP) score, so the
    //   refined alignment is provably never worse than the progressive seed.
    //
    //   SP score: "the MSA program optimizes the sum of all of the pairs of characters at each
    //   position in the alignment (the so-called sum of pair score)" (Wikipedia "Multiple sequence
    //   alignment"). Computed by the existing ComputeSumOfPairsScore (column-based: match/mismatch
    //   from the matrix, residue-gap = GapExtend, gap-gap neutral).
    //
    // Determinism: edges are enumerated in a fixed order (decreasing distance from the root, i.e.
    // internal nodes nearest the leaves first, with a deterministic tie-break); profile splitting,
    // re-projection and the accept-on-improvement rule are deterministic. No RNG is used.
    //
    // Sources retrieved this session (2026-06-23):
    //   - Edgar RC (2004) "MUSCLE: multiple sequence alignment with high accuracy and high
    //     throughput", Nucleic Acids Res 32(5):1792-1797.
    //     https://academic.oup.com/nar/article/32/5/1792/2380623 — Stage 3 steps 3.1-3.4 quoted above.
    //   - Barton GJ, Sternberg MJ (1987) "A strategy for the rapid multiple alignment of protein
    //     sequences...", J Mol Biol 198(2):327-337. https://pubmed.ncbi.nlm.nih.gov/3430611/ —
    //     iterative refinement: re-align a part of the alignment to the rest, iterate to a final
    //     alignment.
    //   - Wikipedia "Multiple sequence alignment", §Iterative methods / sum-of-pairs score
    //     https://en.wikipedia.org/wiki/Multiple_sequence_alignment — SP-score definition; iterative
    //     methods "can return to previously calculated ... sub-MSAs ... optimizing a general
    //     objective function".

    // Default cap on full refinement passes over all edges (Edgar's "user-defined limit"). A pass
    // that makes no accepted change means convergence, so this is only an upper bound.
    private const int DefaultMaxRefinementIterations = 16;

    /// <summary>
    /// Performs <b>iterative refinement</b> of a progressive multiple sequence alignment using
    /// MUSCLE-style tree-dependent restricted partitioning (Edgar 2004, Stage 3), removing the
    /// single-pass "once a gap, always a gap" limitation of
    /// <see cref="MultipleAlignProgressive(IEnumerable{DnaSequence}, ScoringMatrix?)"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Pipeline (Edgar 2004, Nucleic Acids Res 32(5):1792-1797; Barton &amp; Sternberg 1987):
    /// </para>
    /// <list type="number">
    /// <item>Build the progressive (Feng-Doolittle; UPGMA by default, <see cref="ProgressiveMsaOptions.GuideTree"/>) seed alignment and keep its guide tree.</item>
    /// <item>Visit each internal guide-tree edge (leaves-first, deterministic order). Deleting the
    /// edge partitions the sequences into two groups.</item>
    /// <item>Project the current alignment onto each group, drop columns that became all-gap, and
    /// realign the two sub-profiles with the existing profile-profile Needleman-Wunsch.</item>
    /// <item>Accept the re-alignment only if the sum-of-pairs (SP) score does not decrease.</item>
    /// <item>Repeat full passes until no edge yields an improvement (convergence) or the iteration
    /// cap is reached.</item>
    /// </list>
    /// <para>
    /// The result is therefore guaranteed to have an SP score <b>≥</b> that of the progressive seed,
    /// and is deterministic (no RNG; fixed edge order). The star
    /// <see cref="MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?)"/> and progressive
    /// aligners are unchanged.
    /// </para>
    /// </remarks>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <param name="maxIterations">
    /// Upper bound on full refinement passes (Edgar's "user-defined limit"); must be positive.
    /// Defaults to <see cref="DefaultMaxRefinementIterations"/>. Convergence usually stops earlier.
    /// </param>
    /// <returns>Multiple alignment result (aligned rows, majority consensus, sum-of-pairs score).</returns>
    /// <exception cref="ArgumentNullException">When <paramref name="sequences"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">When <paramref name="maxIterations"/> &lt; 1.</exception>
    public static MultipleAlignmentResult MultipleAlignIterative(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring = null,
        int maxIterations = DefaultMaxRefinementIterations) =>
        MultipleAlignIterative(sequences, scoring, new ProgressiveMsaOptions(), maxIterations);

    /// <summary>
    /// <see cref="MultipleAlignIterative(IEnumerable{DnaSequence}, ScoringMatrix?, int)"/> with additive
    /// <see cref="ProgressiveMsaOptions"/> (guide-tree distance source, <see cref="MsaDistanceMode"/>).
    /// </summary>
    /// <param name="sequences">Collection of sequences to align.</param>
    /// <param name="scoring">Scoring matrix (default: <see cref="SimpleDna"/>).</param>
    /// <param name="options">Options; <c>null</c> = defaults.</param>
    /// <param name="maxIterations">Upper bound on full refinement passes (&gt;= 1).</param>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static MultipleAlignmentResult MultipleAlignIterative(
        IEnumerable<DnaSequence> sequences,
        ScoringMatrix? scoring,
        ProgressiveMsaOptions? options,
        int maxIterations = DefaultMaxRefinementIterations)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        ArgumentOutOfRangeException.ThrowIfLessThan(maxIterations, 1);
        return IterativeAlignCore(DnaRows(sequences), scoring ?? SimpleDna, options, maxIterations);
    }

    /// <summary>
    /// <see cref="MultipleAlignIterative(IEnumerable{DnaSequence}, ScoringMatrix?, ProgressiveMsaOptions?, int)"/> for <b>protein</b> sequences (or any alphabet of
    /// the <see cref="SubstitutionMatrix"/>): the same MUSCLE-style tree-dependent restricted partitioning of the progressive seed, accepted on a
    /// strictly better sum-of-pairs score under the matrix and gap model. See
    /// <see cref="MultipleAlign(IEnumerable{ProteinSequence}, MatrixScoring, ProgressiveMsaOptions?)"/> for the input contract.
    /// </summary>
    public static MultipleAlignmentResult MultipleAlignIterative(
        IEnumerable<ProteinSequence> sequences,
        MatrixScoring scoring,
        ProgressiveMsaOptions? options = null,
        int maxIterations = DefaultMaxRefinementIterations)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(maxIterations, 1);
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(IterativeAlignCore(rows, ps, options, maxIterations), ps);
    }

    /// <summary><see cref="MultipleAlignIterative(IEnumerable{ProteinSequence}, MatrixScoring, ProgressiveMsaOptions?, int)"/> for raw sequence strings.</summary>
    public static MultipleAlignmentResult MultipleAlignIterative(
        IEnumerable<string> sequences,
        MatrixScoring scoring,
        ProgressiveMsaOptions? options = null,
        int maxIterations = DefaultMaxRefinementIterations)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(maxIterations, 1);
        var (rows, ps) = PrepareMatrixMsa(sequences, scoring);
        return ScaleMsaResult(IterativeAlignCore(rows, ps, options, maxIterations), ps);
    }

    private static MultipleAlignmentResult IterativeAlignCore(string[] rawSeqs, PairScoring effectiveScoring, ProgressiveMsaOptions? options, int maxIterations)
    {
        using var cancelScope = AlignmentCancellation.Enter(options?.CancellationToken ?? default);
        options ??= new ProgressiveMsaOptions();
        var gapModel = ResolveGapModel(options);

        var trivial = TrivialMultipleAlignment(rawSeqs);
        if (trivial != null)
            return trivial;

        int k = rawSeqs.Length;

        // Step 1: progressive seed (same pipeline as MultipleAlignProgressive) + retain guide tree.
        var distance = ComputeGuideDistances(rawSeqs, effectiveScoring, options.Distance, gapModel, options.GuideTree, options.QuickPairwise, options.ClustalWProtein, options.FastDistanceCells);

        ProgressiveGuideNode root = BuildProgressiveGuideTree(k, distance, options.GuideTree);
        Profile current = AlignProfileSubtree(root, rawSeqs, effectiveScoring, gapModel);

        // Step 2-5: iterative refinement over the guide-tree edges.
        current = RefineByTreePartitioning(current, root, effectiveScoring, maxIterations, gapModel);

        // Reproject rows back into input order.
        var orderedRows = new string[k];
        for (int r = 0; r < current.Rows.Count; r++)
            orderedRows[current.SequenceIndices[r]] = current.Rows[r];

        var alignedList = orderedRows.ToList();
        int maxLen = alignedList.Count == 0 ? 0 : alignedList.Max(s => s.Length);
        string consensus = BuildConsensus(alignedList, maxLen, effectiveScoring.IsMatrix);
        int spScore = ComputeSumOfPairsScore(alignedList, effectiveScoring, gapModel);

        return new MultipleAlignmentResult(
            AlignedSequences: orderedRows,
            Consensus: consensus,
            TotalScore: spScore);
    }

    /// <summary>
    /// Repeatedly splits the current alignment at each internal guide-tree edge into two
    /// sub-profiles, realigns them with profile-profile NW, and keeps the result only when the SP
    /// score does not decrease (Edgar 2004 Stage 3, steps 3.1-3.4). Iterates full passes until a
    /// pass accepts no change (convergence) or <paramref name="maxIterations"/> is reached.
    /// </summary>
    private static Profile RefineByTreePartitioning(
        Profile current, ProgressiveGuideNode root, PairScoring scoring, int maxIterations, MsaGapModel gapModel)
    {
        // 3.1 enumerate the partitions induced by deleting each internal edge of the guide tree.
        // An internal edge is the edge above each non-root internal node; deleting it isolates that
        // node's leaf set as one group and the remaining leaves as the other. Edges are ordered by
        // decreasing distance from the root (deepest internal nodes — nearest the leaves — first),
        // with a deterministic tie-break on the sorted leaf-set, matching Edgar's visiting order.
        var partitions = EnumerateEdgePartitions(root);
        if (partitions.Count == 0)
            return current; // 2 leaves: the single edge is the root split, no internal edge to refine.

        // The SP score the existing pipeline reports for this alignment (input-order independent).
        int CurrentSpScore(Profile p) => ComputeSumOfPairsScore(p.Rows, scoring, gapModel);

        int bestScore = CurrentSpScore(current);

        for (int iter = 0; iter < maxIterations; iter++)
        {
            bool improvedThisPass = false;

            foreach (var group in partitions)
            {
                // 3.1 split the current alignment into the two leaf groups.
                var groupSet = group;
                var sideA = SplitProfile(current, idx => groupSet.Contains(idx));
                var sideB = SplitProfile(current, idx => !groupSet.Contains(idx));
                if (sideA.Rows.Count == 0 || sideB.Rows.Count == 0)
                    continue;

                // 3.2/3.3 realign the two profiles with the existing profile-profile NW.
                Profile candidate = gapModel == MsaGapModel.Affine
                    ? AlignProfilesAffine(sideA, sideB, scoring)
                    : AlignProfiles(sideA, sideB, scoring);
                int candidateScore = CurrentSpScore(candidate);

                // 3.4 keep only if the SP score does not decrease. Strict ">" would also be valid;
                // ">=" lets an equal-score rearrangement settle, but acceptance still never lowers
                // the score, so the monotonic-non-decreasing guarantee holds. To keep convergence
                // well-defined we only treat a STRICT improvement as "made progress".
                if (candidateScore > bestScore)
                {
                    current = candidate;
                    bestScore = candidateScore;
                    improvedThisPass = true;
                }
            }

            if (!improvedThisPass)
                break; // convergence: a full pass changed nothing.
        }

        return current;
    }

    /// <summary>
    /// Returns, for each internal (non-root) node of the guide tree, the set of input-sequence
    /// indices in its subtree — i.e. one side of the partition obtained by deleting that node's
    /// parent edge. Ordered by decreasing distance from the root (deepest first) with a
    /// deterministic tie-break, per Edgar (2004) "edges are visited in order of decreasing distance
    /// from the root".
    /// </summary>
    private static List<HashSet<int>> EnumerateEdgePartitions(ProgressiveGuideNode root)
    {
        var result = new List<(int Depth, List<int> Sorted, HashSet<int> Set)>();

        void Visit(ProgressiveGuideNode node, int depth, bool isRoot)
        {
            if (node.IsLeaf)
                return;

            // The root edge is not an internal edge (deleting it gives the trivial whole-vs-empty
            // split handled by the seed); only record proper internal nodes.
            if (!isRoot)
            {
                var leaves = new List<int>();
                CollectLeaves(node, leaves);
                leaves.Sort();
                result.Add((depth, leaves, new HashSet<int>(leaves)));
            }

            Visit(node.Left!, depth + 1, false);
            Visit(node.Right!, depth + 1, false);
        }

        Visit(root, 0, true);

        // Decreasing distance from the root = larger depth first; tie-break deterministically on the
        // sorted leaf list so the iteration order is fully reproducible.
        result.Sort((a, b) =>
        {
            int byDepth = b.Depth.CompareTo(a.Depth);
            if (byDepth != 0) return byDepth;
            return CompareIntLists(a.Sorted, b.Sorted);
        });

        return result.Select(r => r.Set).ToList();
    }

    private static int CompareIntLists(List<int> a, List<int> b)
    {
        int n = Math.Min(a.Count, b.Count);
        for (int i = 0; i < n; i++)
        {
            int c = a[i].CompareTo(b[i]);
            if (c != 0) return c;
        }
        return a.Count.CompareTo(b.Count);
    }

    private static void CollectLeaves(ProgressiveGuideNode node, List<int> leaves)
    {
        if (node.IsLeaf)
        {
            leaves.Add(node.LeafIndex);
            return;
        }
        CollectLeaves(node.Left!, leaves);
        CollectLeaves(node.Right!, leaves);
    }

    /// <summary>
    /// Projects the current alignment onto the rows whose input-sequence index satisfies
    /// <paramref name="keep"/>, then drops every column that is all-gap within the selected rows.
    /// The result is a valid profile of the selected subset; gaps inside retained columns are
    /// preserved (no existing column is edited — "once a gap, always a gap" within the sub-profile).
    /// </summary>
    private static Profile SplitProfile(Profile current, Func<int, bool> keep)
    {
        var rowIdx = new List<int>();
        for (int r = 0; r < current.Rows.Count; r++)
            if (keep(current.SequenceIndices[r]))
                rowIdx.Add(r);

        var sub = new Profile();
        if (rowIdx.Count == 0)
            return sub;

        int width = current.Width;
        // Identify columns that are NOT all-gap within the selected rows.
        var keptCols = new List<int>(width);
        for (int c = 0; c < width; c++)
        {
            bool allGap = true;
            foreach (int r in rowIdx)
            {
                if (current.Rows[r][c] != '-')
                {
                    allGap = false;
                    break;
                }
            }
            if (!allGap)
                keptCols.Add(c);
        }

        foreach (int r in rowIdx)
        {
            var sb = new StringBuilder(keptCols.Count);
            foreach (int c in keptCols)
                sb.Append(current.Rows[r][c]);
            sub.Rows.Add(sb.ToString());
            sub.SequenceIndices.Add(current.SequenceIndices[r]);
        }

        return sub;
    }

    #endregion
}
