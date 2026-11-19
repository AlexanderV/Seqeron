using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-36: tools for the public stage ports that had no tool (the multiple_align_* tools run only whole pipelines): ClustalW pairwise stages,
// T-Coffee pair method / lalign / primary library / guide tree / k-tuple similarity, and the sum-of-pairs score of an alignment.
// Every tool delegates to SequenceAligner; MUSCLE's guide tree of an alignment is muscle_msa_command (mode 'maketree').
public partial class AlignmentTools
{
    [McpServerTool(Name = "clustalw_pairwise_align", Title = "Pairwise — ClustalW 2.1 Pairwise Stage (full DP)", ReadOnly = true)]
    [Description("The full pairwise stage of ClustalW 2.1 (clustalw -ALIGN: FullPairwiseAlign::pairwiseAlign, the step that gives the guide-tree distances), ported step by step: an affine local (Smith-Waterman-Gotoh) alignment under ClustalW's matrix whose co-optimal choice follows the program's Myers-Miller midpoint rule, with the percent identity (identities / the shorter sequence, truncated as the program prints it) and the guide-tree distance 1 - identities / shorter length. type 'dna' (default; -PWDNAMATRIX IUB (default) or CLUSTALW, -PWGAPOPEN 15, -PWGAPEXT 6.66, -TRANSWEIGHT 0.5; IUB nucleotide codes) or 'protein' (-PWMATRIX gonnet (default, GONNET250), blosum (BLOSUM30), pam (PAM350) or id; -PWGAPOPEN 10, -PWGAPEXT 0.1). Verified against the installed clustalw 2.1 on 11 000+ DNA pairs. start / end are 0-based inclusive coordinates of the aligned region (-1 when no positive-scoring alignment exists). The multi-sequence pipeline is multiple_align_clustalw.")]
    public static ClustalWPairwiseAlignResult ClustalWPairwiseAlign(
        [Description("First sequence: DNA (A, C, G, T/U, IUB codes) or protein (letters A-Z except J), case-insensitive.")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("'dna' (default) or 'protein' (ClustalW -TYPE).")] string type = "dna",
        [Description("ClustalW -PWGAPOPEN (a penalty, positive; default 15 for DNA, 10 for protein).")] double? gapOpen = null,
        [Description("ClustalW -PWGAPEXT (default 6.66 for DNA, 0.1 for protein).")] double? gapExtend = null,
        [Description("DNA only, ClustalW -TRANSWEIGHT: weight of transitions in the pairwise matrix (default 0.5).")] float transitionWeight = 0.5f,
        [Description("type='dna': ClustalW -PWDNAMATRIX built-in, 'iub' (default) or 'clustalw'.")] string? dnaMatrix = null,
        [Description("type='protein': ClustalW -PWMATRIX, 'gonnet' (default), 'blosum', 'pam' or 'id'.")] string? pwMatrix = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        var s1 = RequireSequence(sequence1, nameof(sequence1));
        var s2 = RequireSequence(sequence2, nameof(sequence2));
        global::Seqeron.Genomics.Infrastructure.ClustalWPairwiseResult r;
        switch ((type ?? "dna").Trim().ToLowerInvariant())
        {
            case "dna":
                r = global::Seqeron.Genomics.Alignment.SequenceAligner.ClustalWPairwiseAlign(s1, s2, new global::Seqeron.Genomics.Infrastructure.ClustalWDnaPairwiseOptions
                {
                    GapOpen = gapOpen ?? 15.0,
                    GapExtend = gapExtend ?? 6.66,
                    TransitionWeight = transitionWeight,
                    Matrix = ParseClustalWDnaMatrix(dnaMatrix, nameof(dnaMatrix)),
                });
                break;
            case "protein":
                r = global::Seqeron.Genomics.Alignment.SequenceAligner.ClustalWPairwiseAlign(s1, s2, new global::Seqeron.Genomics.Infrastructure.ClustalWProteinPairwiseOptions
                {
                    GapOpen = gapOpen ?? 10.0,
                    GapExtend = gapExtend ?? 0.1,
                    Matrix = (pwMatrix ?? "gonnet").Trim().ToLowerInvariant() switch
                    {
                        "gonnet" => global::Seqeron.Genomics.Infrastructure.ClustalWPairwiseMatrix.Gonnet,
                        "blosum" => global::Seqeron.Genomics.Infrastructure.ClustalWPairwiseMatrix.Blosum,
                        "pam" => global::Seqeron.Genomics.Infrastructure.ClustalWPairwiseMatrix.Pam,
                        "id" or "identity" => global::Seqeron.Genomics.Infrastructure.ClustalWPairwiseMatrix.Identity,
                        _ => throw new ArgumentException("pwMatrix must be 'gonnet', 'blosum', 'pam' or 'id'.", nameof(pwMatrix)),
                    },
                });
                break;
            default:
                throw new ArgumentException("type must be 'dna' or 'protein'.", nameof(type));
        }

        return new ClustalWPairwiseAlignResult(r.Score, r.Identities, r.PercentIdentity, r.Distance, r.ClustalWScore,
            r.AlignedSequence1, r.AlignedSequence2, r.Start1, r.End1, r.Start2, r.End2);
    }

    [McpServerTool(Name = "clustalw_quick_pairwise_align", Title = "Pairwise — ClustalW 2.1 Fast Pairwise Stage (k-tuple)", ReadOnly = true)]
    [Description("The fast pairwise stage of ClustalW 2.1 (clustalw -QUICKTREE -ALIGN: FastPairwiseAlign::pairAlign; Wilbur & Lipman 1983 k-tuple diagonals), an exact port with the score conversion of pairwiseAlign: the best chain of k-tuple fragments (fragmentScore), the percent score fragmentScore / shorter length x 100 (percent = true; the raw fragment score otherwise) and the guide-tree distance (100 - score) / 100. type 'dna' (default; -KTUPLE 2, -WINDOW 4, -TOPDIAGS 4, -PAIRGAP 5) or 'protein' (1 / 5 / 5 / 3). Not an alignment identity. maxSequenceLength = the longest sequence of the whole input set (ClustalW bounds its fragment list by 2 x 2 x this value and reports partialAlignment beyond it; default the longer of the two). Verified against the installed clustalw 2.1.")]
    public static ClustalWQuickPairwiseAlignResult ClustalWQuickPairwiseAlign(
        [Description("First sequence: DNA (A, C, G, T/U, IUB codes) or protein (letters A-Z except J).")] string sequence1,
        [Description("Second sequence.")] string sequence2,
        [Description("'dna' (default) or 'protein'.")] string type = "dna",
        [Description("ClustalW -KTUPLE (DNA 1-4, protein 1-2; default 2 for DNA, 1 for protein).")] int? ktuple = null,
        [Description("ClustalW -WINDOW (default 4 DNA, 5 protein).")] int? window = null,
        [Description("ClustalW -TOPDIAGS (default 4 DNA, 5 protein).")] int? topDiagonals = null,
        [Description("ClustalW -PAIRGAP (default 5 DNA, 3 protein).")] int? pairGap = null,
        [Description("ClustalW -SCORE=PERCENT (true, default) or ABSOLUTE (false).")] bool percent = true,
        [Description("Longest sequence of the whole input set (>= the longer of the two); default the longer of the two.")] int? maxSequenceLength = null)
    {
        var s1 = RequireSequence(sequence1, nameof(sequence1));
        var s2 = RequireSequence(sequence2, nameof(sequence2));
        var kind = (type ?? "dna").Trim().ToLowerInvariant();
        if (kind is not ("dna" or "protein"))
            throw new ArgumentException("type must be 'dna' or 'protein'.", nameof(type));
        var defaults = kind == "dna" ? new global::Seqeron.Genomics.Infrastructure.ClustalWQuickOptions() : new global::Seqeron.Genomics.Infrastructure.ClustalWProteinPairwiseOptions().Quick;
        var options = defaults with
        {
            KTuple = ktuple ?? defaults.KTuple,
            Window = window ?? defaults.Window,
            TopDiagonals = topDiagonals ?? defaults.TopDiagonals,
            PairGap = pairGap ?? defaults.PairGap,
            Percent = percent,
        };
        var r = global::Seqeron.Genomics.Alignment.SequenceAligner.ClustalWQuickPairwiseAlign(
            s1, s2,
            kind == "dna" ? global::Seqeron.Genomics.Infrastructure.ClustalWSequenceType.Dna : global::Seqeron.Genomics.Infrastructure.ClustalWSequenceType.Protein,
            options, maxSequenceLength);
        return new ClustalWQuickPairwiseAlignResult(r.FragmentScore, r.Score, r.Distance, r.Aligned, r.PartialAlignment);
    }

    [McpServerTool(Name = "tcoffee_slow_pair_align", Title = "Pairwise — T-Coffee slow_pair (Myers-Miller global)", ReadOnly = true)]
    [Description("The 'slow_pair' pair-wise method of T-Coffee 13.41 (the global alignment behind its primary libraries): Myers-Miller linear-space global alignment with the program's matrices and penalties - DNA: dna_idmat (10 / 0), gap open -10, extend -1; protein: BLOSUM62. gapOpen / gapExtend / matrix are the method's own overrides (t_coffee -method slow_pair@EP@GOP@-30@GEP@-3@MATRIX@blosum50mt); matrix names are the T-Coffee matrix names (blosum62mt, blosum50mt, pam250mt, dna_idmat, ...). Returns the two aligned rows. Verified against the program's slow_pair libraries (see multiple_align_tcoffee).")]
    public static AlignedPairResult TCoffeeSlowPairAlign(
        [Description("First sequence (letters; gaps rejected).")] string sequenceA,
        [Description("Second sequence.")] string sequenceB,
        [Description("true = protein matrices / penalties, false (default) = DNA.")] bool protein = false,
        [Description("The method's MATRIX override (T-Coffee matrix name); default per type.")] string? matrix = null,
        [Description("The method's GOP override (gap open, an integer score, negative = penalty).")] int? gapOpen = null,
        [Description("The method's GEP override (gap extend, an integer score, negative = penalty).")] int? gapExtend = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        var parameters = gapOpen is null && gapExtend is null ? null : new global::Seqeron.Genomics.Alignment.TCoffeeMethodParameters(gapOpen, gapExtend);
        var (a, b) = global::Seqeron.Genomics.Alignment.SequenceAligner.TCoffeeSlowPairAlign(RequireSequence(sequenceA, nameof(sequenceA)), RequireSequence(sequenceB, nameof(sequenceB)), protein, matrix, parameters);
        return new AlignedPairResult(a, b);
    }

    [McpServerTool(Name = "tcoffee_lalign_top_alignments", Title = "Pairwise — T-Coffee lalign_id_pair (Huang-Miller SIM, top-N local)", ReadOnly = true)]
    [Description("The 'lalign_id_pair' pair-wise method of T-Coffee 13.41: the topAlignments (t_coffee -lalign_n_top, default 10) best NON-INTERSECTING local alignments of two sequences by Huang & Miller's SIM algorithm, with the program's scoring (DNA: match 10 / mismatch 0 style identity matrix, gaps -10 / -1; protein: BLOSUM62), best first. Each alignment has its two aligned rows, 0-based starts and score. Verified against the program's lalign libraries (see multiple_align_tcoffee).")]
    public static TCoffeeLocalAlignmentsResult TCoffeeLalignTopAlignments(
        [Description("First sequence (letters; gaps rejected).")] string sequenceA,
        [Description("Second sequence.")] string sequenceB,
        [Description("t_coffee -lalign_n_top: number of alignments to report (>= 1; default 10).")] int topAlignments = 10,
        [Description("true = protein scoring, false (default) = DNA.")] bool protein = false,
        [Description("Method GOP override (gap open score, negative).")] int? gapOpen = null,
        [Description("Method GEP override (gap extend score, negative).")] int? gapExtend = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        var parameters = gapOpen is null && gapExtend is null ? null : new global::Seqeron.Genomics.Alignment.TCoffeeMethodParameters(gapOpen, gapExtend);
        var items = global::Seqeron.Genomics.Alignment.SequenceAligner.TCoffeeLalignTopAlignments(
                RequireSequence(sequenceA, nameof(sequenceA)), RequireSequence(sequenceB, nameof(sequenceB)), topAlignments, protein, parameters)
            .Select(x => new TCoffeeLocalAlignmentItem(x.RowA, x.RowB, x.StartA, x.StartB, x.Score))
            .ToArray();
        return new TCoffeeLocalAlignmentsResult(items, items.Length);
    }

    [McpServerTool(Name = "tcoffee_primary_library", Title = "MSA — T-Coffee Primary Library", ReadOnly = true)]
    [Description("The primary library the T-Coffee 13.41 program builds for a set of sequences with the given pair-wise methods (t_coffee -seq s.fa -method m1,m2 -lib_only -out_lib x.tc_lib): residue pairs (sequence index of the first smaller than that of the second; 0-based positions) with weight 10 x percent identity of the alignment they come from and the number of times the pair was added. methods in command-line order: 'proba_pair' (default; ProbCons pair-HMM posteriors), 'slow_pair', 'lalign_id_pair', 'clustalw_pair', 'biphasic_pair' (two affine gap models), 'lalign_rs_s_pair' (SIM on every pair including each sequence against itself), 'promo_pair' / 'blastr_pair' (slow_pair on the dinucleotide-extended alphabet); the program refuses 'clean_slow_pair', 'hash_pair', 'subop1_pair', 'subop2_pair', 'cfast_pair' (TCoffeeAbortedException); with two or more methods the program's own merge loses and duplicates entries - the result is what it writes. Verified on 800+ committed libraries (see multiple_align_tcoffee). Entries are ordered by (seqA, seqB, posA, posB).")]
    public static TCoffeePrimaryLibraryResult TCoffeePrimaryLibrary(
        [Description("Sequences (at least one, none empty; letters only).")] string[] sequences,
        [Description("Pair-wise methods in command-line order ('proba_pair' (default), 'slow_pair', 'lalign_id_pair', 'clustalw_pair', 'biphasic_pair', 'lalign_rs_s_pair', 'promo_pair', 'blastr_pair').")] string[]? methods = null,
        [Description("t_coffee -lalign_n_top (default 10).")] int lalignTopAlignments = 10,
        [Description("true = protein, false = DNA, omitted = decided from the letters as the program does.")] bool? protein = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("sequences must contain at least one sequence.", nameof(sequences));
        var parsed = (methods is { Length: > 0 } ? methods : ["proba_pair"]).Select(ParseTCoffeePairMethod).ToList();
        var items = global::Seqeron.Genomics.Alignment.SequenceAligner.TCoffeePrimaryLibraryEntries(sequences, parsed, lalignTopAlignments, protein)
            .Select(e => new TCoffeeLibraryEntryItem(e.SeqA, e.PosA, e.SeqB, e.PosB, e.Weight, e.Count))
            .ToArray();
        return new TCoffeePrimaryLibraryResult(items, items.Length);
    }

    private static global::Seqeron.Genomics.Alignment.TCoffeePairMethod ParseTCoffeePairMethod(string name) => (name ?? string.Empty).Trim().ToLowerInvariant() switch
    {
        "slow_pair" => global::Seqeron.Genomics.Alignment.TCoffeePairMethod.SlowPair,
        "lalign_id_pair" => global::Seqeron.Genomics.Alignment.TCoffeePairMethod.LalignIdPair,
        "proba_pair" => global::Seqeron.Genomics.Alignment.TCoffeePairMethod.ProbaPair,
        "clustalw_pair" => global::Seqeron.Genomics.Alignment.TCoffeePairMethod.ClustalwPair,
        "biphasic_pair" => global::Seqeron.Genomics.Alignment.TCoffeePairMethod.BiphasicPair,
        "lalign_rs_s_pair" => global::Seqeron.Genomics.Alignment.TCoffeePairMethod.LalignRsSPair,
        "promo_pair" => global::Seqeron.Genomics.Alignment.TCoffeePairMethod.PromoPair,
        "blastr_pair" => global::Seqeron.Genomics.Alignment.TCoffeePairMethod.BlastrPair,
        "clean_slow_pair" => global::Seqeron.Genomics.Alignment.TCoffeePairMethod.CleanSlowPair,
        "hash_pair" => global::Seqeron.Genomics.Alignment.TCoffeePairMethod.HashPair,
        "subop1_pair" => global::Seqeron.Genomics.Alignment.TCoffeePairMethod.Subop1Pair,
        "subop2_pair" => global::Seqeron.Genomics.Alignment.TCoffeePairMethod.Subop2Pair,
        "cfast_pair" => global::Seqeron.Genomics.Alignment.TCoffeePairMethod.CfastPair,
        _ => throw new ArgumentException($"Unknown pair-wise method '{name}': expected proba_pair, slow_pair, lalign_id_pair, clustalw_pair, biphasic_pair, lalign_rs_s_pair, promo_pair, blastr_pair (or clean_slow_pair, hash_pair, subop1_pair, subop2_pair, cfast_pair, which the program refuses)."),
    };

    [McpServerTool(Name = "tcoffee_guide_tree", Title = "MSA — T-Coffee Guide Tree", ReadOnly = true)]
    [Description("The guide tree T-Coffee 13.41 builds for a set of sequences with its defaults (-tree_mode nj -distance_matrix_mode ktup): k-tuple similarity distances, the program's neighbour-joining routine and the rooting of its trifurcation into a binary tree. Returns the tree as Newick text (no branch lengths; leaf labels = names, default s0, s1, ...; left subtree = first profile of the merge) and the leaf indices left to right; the progressive merge order of multiple_align_tcoffee is the post-order of this tree. protein omitted = decided from the letters.")]
    public static TCoffeeGuideTreeResult TCoffeeGuideTree(
        [Description("Sequences (at least one, none empty; letters only).")] string[] sequences,
        [Description("true = protein, false = DNA, omitted = decided from the letters.")] bool? protein = null,
        [Description("Leaf labels (default s0, s1, ...); one per sequence.")] string[]? names = null)
    {
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("sequences must contain at least one sequence.", nameof(sequences));
        if (names is not null && names.Length != sequences.Length)
            throw new ArgumentException("names must have one label per sequence.", nameof(names));
        var tree = global::Seqeron.Genomics.Alignment.SequenceAligner.TCoffeeGuideTree(sequences, protein);
        var order = new List<int>();
        string Newick(global::Seqeron.Genomics.Alignment.TCoffeeTreeNode node)
        {
            if (node.Left is null || node.Right is null)
            {
                order.Add(node.LeafIndex);
                return names?[node.LeafIndex] ?? $"s{node.LeafIndex}";
            }

            return "(" + Newick(node.Left) + "," + Newick(node.Right) + ")";
        }

        string text = Newick(tree) + ";";
        return new TCoffeeGuideTreeResult(text, order.ToArray());
    }

    [McpServerTool(Name = "tcoffee_ktup_similarity", Title = "MSA — T-Coffee k-tuple Similarity Matrix", ReadOnly = true)]
    [Description("The pairwise k-tuple similarity matrix T-Coffee 13.41 computes for its guide tree (-distance_matrix_mode ktup: shared k-tuples of the two sequences, k = 1 for protein / 2 for DNA reduced alphabet as in the program), as the integer matrix the program stores (symmetric; the diagonal is the self similarity). Verified together with the guide tree (see tcoffee_guide_tree).")]
    public static TCoffeeKtupSimilarityResult TCoffeeKtupSimilarity(
        [Description("Sequences (at least one, none empty; letters only).")] string[] sequences,
        [Description("true = protein, false (default) = DNA.")] bool protein = false)
    {
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("sequences must contain at least one sequence.", nameof(sequences));
        var m = global::Seqeron.Genomics.Alignment.SequenceAligner.TCoffeeKtupSimilarity(sequences, protein);
        var rows = new int[m.GetLength(0)][];
        for (int i = 0; i < rows.Length; i++)
        {
            rows[i] = new int[m.GetLength(1)];
            for (int j = 0; j < rows[i].Length; j++)
                rows[i][j] = m[i, j];
        }

        return new TCoffeeKtupSimilarityResult(rows);
    }

    [McpServerTool(Name = "sum_of_pairs_score", Title = "MSA — Sum-of-Pairs Score of an Alignment", ReadOnly = true)]
    [Description("Sum-of-pairs score of a multiple alignment (equal-length rows, '-' = gap). gapModel 'linear' (default): column sum of pairs - residue pairs match / mismatch, residue against gap gapExtend, gap against gap 0 (the totalScore of the linear aligners); 'affine': the natural affine objective (Altschul 1989) - for every pair of rows the gap-gap columns are removed and every maximal gap run costs gapOpen + length x gapExtend. Scores: match / mismatch / gapOpen / gapExtend (integers; negative = penalty); with matrix (BLOSUM62, PAM250, EDNAFULL or any EMBOSS matrix) the residue pairs are scored by that table (protein / IUPAC rows, residues case-insensitive) and the gap scores may be fractional - the score is then a real number. This is the objective the exact (Carrillo-Lipman) and affine aligners optimise; MUSCLE's own SP score is muscle_msa_command mode 'spscore'.")]
    public static SumOfPairsScoreResult SumOfPairsScore(
        [Description("Equal-length aligned rows, '-' = gap.")] string[] alignedRows,
        [Description("Score for a matching residue pair (ignored with matrix).")] int match = 1,
        [Description("Score for a mismatching residue pair (ignored with matrix).")] int mismatch = -1,
        [Description("Gap-open score added once per gap run (affine model only; negative = penalty; integer unless matrix is given).")] double gapOpen = -2,
        [Description("Gap score per gap position (negative = penalty; integer unless matrix is given).")] double gapExtend = -1,
        [Description("Gap model: 'linear' (default) or 'affine'.")] string gapModel = "linear",
        [Description(MatrixParamDescription)] string? matrix = null)
    {
        if (alignedRows is null || alignedRows.Length == 0)
            throw new ArgumentException("alignedRows must contain at least one row.", nameof(alignedRows));
        var model = (gapModel ?? "linear").Trim().ToLowerInvariant() switch
        {
            "linear" => global::Seqeron.Genomics.Infrastructure.MsaGapModel.Linear,
            "affine" => global::Seqeron.Genomics.Infrastructure.MsaGapModel.Affine,
            _ => throw new ArgumentException("gapModel must be 'linear' or 'affine'.", nameof(gapModel)),
        };
        var table = ResolveMatrix(matrix);
        if (table is not null)
        {
            double score = global::Seqeron.Genomics.Alignment.SequenceAligner.MatrixSumOfPairsScore(
                alignedRows, new global::Seqeron.Genomics.Infrastructure.MatrixScoring(table, gapOpen, gapExtend), model);
            return new SumOfPairsScoreResult(score, model.ToString(), $"matrix {matrix}");
        }

        int total = global::Seqeron.Genomics.Alignment.SequenceAligner.SumOfPairsScore(
            alignedRows,
            new global::Seqeron.Genomics.Infrastructure.ScoringMatrix(match, mismatch, RequireInteger(gapOpen, nameof(gapOpen)), RequireInteger(gapExtend, nameof(gapExtend))),
            model);
        return new SumOfPairsScoreResult(total, model.ToString(), $"match {match} mismatch {mismatch}");
    }
}
