// SPDX-License-Identifier: GPL-3.0-or-later
// This file is part of Seqeron.Mcp.MultipleAlignment, licensed under the terms stated in LICENSE (this directory).
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.PartialOrderAlignment;
using Seqeron.Genomics.Kalign;
using Seqeron.Genomics.TCoffee;
using Seqeron.Genomics.ClustalW;
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Alignment.Tools;
using static Seqeron.Mcp.Alignment.Tools.AlignmentTools;

namespace Seqeron.Mcp.MultipleAlignment.Tools;

public partial class MultipleAlignmentTools
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
        global::Seqeron.Genomics.ClustalW.ClustalWPairwiseResult r;
        switch ((type ?? "dna").Trim().ToLowerInvariant())
        {
            case "dna":
                r = global::Seqeron.Genomics.ClustalW.ClustalWAligner.ClustalWPairwiseAlign(s1, s2, new global::Seqeron.Genomics.ClustalW.ClustalWDnaPairwiseOptions
                {
                    GapOpen = gapOpen ?? 15.0,
                    GapExtend = gapExtend ?? 6.66,
                    TransitionWeight = transitionWeight,
                    Matrix = ParseClustalWDnaMatrix(dnaMatrix, nameof(dnaMatrix)),
                });
                break;
            case "protein":
                r = global::Seqeron.Genomics.ClustalW.ClustalWAligner.ClustalWPairwiseAlign(s1, s2, new global::Seqeron.Genomics.ClustalW.ClustalWProteinPairwiseOptions
                {
                    GapOpen = gapOpen ?? 10.0,
                    GapExtend = gapExtend ?? 0.1,
                    Matrix = (pwMatrix ?? "gonnet").Trim().ToLowerInvariant() switch
                    {
                        "gonnet" => global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMatrix.Gonnet,
                        "blosum" => global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMatrix.Blosum,
                        "pam" => global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMatrix.Pam,
                        "id" or "identity" => global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMatrix.Identity,
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
        var defaults = kind == "dna" ? new global::Seqeron.Genomics.ClustalW.ClustalWQuickOptions() : new global::Seqeron.Genomics.ClustalW.ClustalWProteinPairwiseOptions().Quick;
        var options = defaults with
        {
            KTuple = ktuple ?? defaults.KTuple,
            Window = window ?? defaults.Window,
            TopDiagonals = topDiagonals ?? defaults.TopDiagonals,
            PairGap = pairGap ?? defaults.PairGap,
            Percent = percent,
        };
        var r = global::Seqeron.Genomics.ClustalW.ClustalWAligner.ClustalWQuickPairwiseAlign(
            s1, s2,
            kind == "dna" ? global::Seqeron.Genomics.ClustalW.ClustalWSequenceType.Dna : global::Seqeron.Genomics.ClustalW.ClustalWSequenceType.Protein,
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
        var parameters = gapOpen is null && gapExtend is null ? null : new global::Seqeron.Genomics.TCoffee.TCoffeeMethodParameters(gapOpen, gapExtend);
        var (a, b) = global::Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeSlowPairAlign(RequireSequence(sequenceA, nameof(sequenceA)), RequireSequence(sequenceB, nameof(sequenceB)), protein, matrix, parameters);
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
        var parameters = gapOpen is null && gapExtend is null ? null : new global::Seqeron.Genomics.TCoffee.TCoffeeMethodParameters(gapOpen, gapExtend);
        var items = global::Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeLalignTopAlignments(
                RequireSequence(sequenceA, nameof(sequenceA)), RequireSequence(sequenceB, nameof(sequenceB)), topAlignments, protein, parameters)
            .Select(x => new TCoffeeLocalAlignmentItem(x.RowA, x.RowB, x.StartA, x.StartB, x.Score))
            .ToArray();
        return new TCoffeeLocalAlignmentsResult(items, items.Length);
    }
    [McpServerTool(Name = "tcoffee_primary_library", Title = "MSA — T-Coffee Primary Library", ReadOnly = true)]
    [Description("The primary library the T-Coffee 13.41 program builds for a set of sequences with the given pair-wise methods (t_coffee -seq s.fa -method m1,m2 -lib_only -out_lib x.tc_lib): residue pairs (sequence index of the first smaller than that of the second; 0-based positions) with weight 10 x percent identity of the alignment they come from and the number of times the pair was added. methods in command-line order: 'proba_pair' (default; ProbCons pair-HMM posteriors), 'slow_pair', 'lalign_id_pair', 'clustalw_pair', 'biphasic_pair' (two affine gap models), 'lalign_rs_s_pair' (SIM on every pair including each sequence against itself), 'promo_pair' / 'blastr_pair' (slow_pair on the dinucleotide-extended alphabet), 'fast_pair' (FASTA heuristic of util_dp_fasta_nw.c: banded Gotoh DP on the best k-tuple diagonals), 'exon_pair' / 'exon2_pair' / 'exon3_pair' (alignment of the genomic templates given in genomicTemplates, expressed in the residues of the sequences); the program refuses 'clean_slow_pair', 'hash_pair', 'subop1_pair', 'subop2_pair', 'cfast_pair' (TCoffeeAbortedException); with two or more methods the program's own merge loses and duplicates entries - the result is what it writes. Verified on 800+ committed libraries (see multiple_align_tcoffee). Entries are ordered by (seqA, seqB, posA, posB).")]
    public static TCoffeePrimaryLibraryResult TCoffeePrimaryLibrary(
        [Description("Sequences (at least one, none empty; letters only).")] string[] sequences,
        [Description("Pair-wise methods in command-line order ('proba_pair' (default), 'slow_pair', 'lalign_id_pair', 'clustalw_pair', 'biphasic_pair', 'lalign_rs_s_pair', 'promo_pair', 'blastr_pair', 'fast_pair', 'exon_pair', 'exon2_pair', 'exon3_pair').")] string[]? methods = null,
        [Description("t_coffee -lalign_n_top (default 10).")] int lalignTopAlignments = 10,
        [Description("true = protein, false = DNA, omitted = decided from the letters as the program does.")] bool? protein = null,
        [Description("The genomic (_G_) template of every sequence for exon_pair / exon2_pair / exon3_pair (one entry per sequence, null / empty = none).")] string?[]? genomicTemplates = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("sequences must contain at least one sequence.", nameof(sequences));
        if (genomicTemplates != null && genomicTemplates.Length != sequences.Length)
            throw new ArgumentException("genomicTemplates must have one entry per sequence.", nameof(genomicTemplates));
        var parsed = (methods is { Length: > 0 } ? methods : ["proba_pair"]).Select(ParseTCoffeePairMethod).ToList();
        var items = global::Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeePrimaryLibraryEntries(sequences, parsed, lalignTopAlignments, protein, genomicTemplates?.Select(t => string.IsNullOrEmpty(t) ? null : t).ToArray())
            .Select(e => new TCoffeeLibraryEntryItem(e.SeqA, e.PosA, e.SeqB, e.PosB, e.Weight, e.Count))
            .ToArray();
        return new TCoffeePrimaryLibraryResult(items, items.Length);
    }
    [McpServerTool(Name = "tcoffee_blast_msa_library", Title = "MSA — T-Coffee blastp_msa Library (from BLAST XML)", ReadOnly = true)]
    [Description("The library t_coffee -method blastp_msa builds (-lib_only -out_lib) from the XML of the all-against-all BLAST search of the sequences (blastall -p blastp -m7, or BLAST+ -outfmt 5): for every high-scoring pair (Hsp_qseq / Hsp_hseq) the residue pairs of its gap-free alignment columns with the HSP's percent identity (integer, truncated) as the weight, an entry and its mirror, in the printed format of tcoffee_primary_library (0-based positions, sequence index of the first smaller than that of the second, self hits not printed). The BLAST search itself is the external program and is NOT run here: pass its XML in blastXml (one Iteration per sequence in input order; the hits are matched to the sequences by their Hit_def name). VERIFICATION: the library equals the one the installed t_coffee writes on every committed set with a library (70 or more oracle sets, scripts/alignment/generate_tcoffee_external_oracle.py).")]
    public static TCoffeePrimaryLibraryResult TCoffeeBlastMsaLibrary(
        [Description("Sequence names exactly as they stand in the BLAST database / query (one per sequence); the hits of the XML are matched by name.")] string[] names,
        [Description("The sequences (at least one, none empty).")] string[] sequences,
        [Description("The BLAST XML output (-m7 / -outfmt 5) of the all-against-all search; one Iteration per sequence in input order.")] string blastXml)
    {
        if (names is null || sequences is null || sequences.Length == 0)
            throw new ArgumentException("names and sequences must contain at least one sequence.", nameof(sequences));
        if (blastXml is null)
            throw new ArgumentException("blastXml is required.", nameof(blastXml));
        var items = global::Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeBlastMsaLibraryEntries(names, sequences, blastXml)
            .Select(e => new TCoffeeLibraryEntryItem(e.SeqA, e.PosA, e.SeqB, e.PosB, e.Weight, e.Count))
            .ToArray();
        return new TCoffeePrimaryLibraryResult(items, items.Length);
    }
    [McpServerTool(Name = "tcoffee_tblastx_library", Title = "MSA — T-Coffee tblastx_msa Library (from BLAST XML)", ReadOnly = true)]
    [Description("The library t_coffee -method tblastx_msa builds (tc_generic_method.pl -mode=tblastx_msa = seq2tblastx_lib) from the XML of the all-against-all tblastx search of nucleotide sequences: for every hit of a later sequence the nucleotide pairs of every high-scoring pair (the query / hit ranges of the translated frames, minus-strand frames mapped back to the plus strand, weight = the HSP's percent identity, integer truncated), an entry and its mirror, in the printed format of tcoffee_primary_library (0-based positions, first sequence index smaller than the second). The tblastx search is the external program and is NOT run here: pass its XML in blastXml (one Iteration per sequence in input order, sequences identified by their 1-based position as the program numbers them). VERIFICATION: the library equals the one the installed t_coffee writes on every committed tblastx_msa set with a library (scripts/alignment/generate_tcoffee_external_oracle.py).")]
    public static TCoffeePrimaryLibraryResult TCoffeeTblastxLibrary(
        [Description("The nucleotide sequences (at least one, none empty).")] string[] sequences,
        [Description("The tblastx XML output (-m7 / -outfmt 5) of the all-against-all search; one Iteration per sequence in input order.")] string blastXml)
    {
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("sequences must contain at least one sequence.", nameof(sequences));
        if (blastXml is null)
            throw new ArgumentException("blastXml is required.", nameof(blastXml));
        var items = global::Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeTblastxLibraryEntries(sequences, blastXml)
            .Select(e => new TCoffeeLibraryEntryItem(e.SeqA, e.PosA, e.SeqB, e.PosB, e.Weight, e.Count))
            .ToArray();
        return new TCoffeePrimaryLibraryResult(items, items.Length);
    }
    [McpServerTool(Name = "clustalw_guide_tree", Title = "MSA — ClustalW 2.1 Guide-Tree File from Distances", ReadOnly = true)]
    [Description("The guide-tree file clustalw writes for its neighbour-joining clustering (UnRootedClusterTree::treeFromDistMatrix + ClusterTreeOutput::printPhylipTree / twoWaySplit) of the sequences `names` from their pairwise distances: '(a:d,b:d);' for two sequences (d = half the distance, five significant digits, as the program writes it), else the unrooted tree as three subtrees under the first parenthesis, every node on its own lines, branch lengths with five decimals (the text of the .dnd file). The distances are ClustalW's own (the fractional differences clustalw computes in its pairwise stage; see clustalw_pairwise_align / clustalw_quick_pairwise_align). VERIFICATION: the neighbour-joining and printing routines are the ones multiple_align_clustalw runs; their guide-tree text equals the .dnd files of the installed clustalw 2.1 (committed guide-tree oracle, random DNA sets).")]
    public static ClustalWGuideTreeResult ClustalWGuideTree(
        [Description("Sequence names in file order (at least two).")] string[] names,
        [Description("ClustalW's pairwise distances: a square symmetric matrix with one row and one column per name (diagonal 0).")] double[][] distances)
    {
        if (names is null || names.Length < 2)
            throw new ArgumentException("names must contain at least two sequence names.", nameof(names));
        if (distances is null || distances.Length != names.Length || distances.Any(row => row is null || row.Length != names.Length))
            throw new ArgumentException("distances must be a square matrix with one row and one column per name.", nameof(distances));
        var matrix = new double[names.Length, names.Length];
        for (int i = 0; i < names.Length; i++)
            for (int j = 0; j < names.Length; j++)
                matrix[i, j] = distances[i][j];
        return new ClustalWGuideTreeResult(global::Seqeron.Genomics.ClustalW.ClustalWAligner.ClustalWGuideTreeText(names, matrix), names.Length);
    }
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
        var tree = global::Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeGuideTree(sequences, protein);
        var order = new List<int>();
        string Newick(global::Seqeron.Genomics.TCoffee.TCoffeeTreeNode node)
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
        var m = global::Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeKtupSimilarity(sequences, protein);
        var rows = new int[m.GetLength(0)][];
        for (int i = 0; i < rows.Length; i++)
        {
            rows[i] = new int[m.GetLength(1)];
            for (int j = 0; j < rows[i].Length; j++)
                rows[i][j] = m[i, j];
        }

        return new TCoffeeKtupSimilarityResult(rows);
    }

    private static global::Seqeron.Genomics.TCoffee.TCoffeePairMethod ParseTCoffeePairMethod(string name) => (name ?? string.Empty).Trim().ToLowerInvariant() switch
    {
        "slow_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.SlowPair,
        "lalign_id_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.LalignIdPair,
        "proba_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.ProbaPair,
        "clustalw_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.ClustalwPair,
        "biphasic_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.BiphasicPair,
        "lalign_rs_s_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.LalignRsSPair,
        "promo_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.PromoPair,
        "blastr_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.BlastrPair,
        "clean_slow_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.CleanSlowPair,
        "hash_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.HashPair,
        "subop1_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.Subop1Pair,
        "subop2_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.Subop2Pair,
        "cfast_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.CfastPair,
        "fast_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.FastPair,
        "exon_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.ExonPair,
        "exon2_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.Exon2Pair,
        "exon3_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.Exon3Pair,
        _ => throw new ArgumentException($"Unknown pair-wise method '{name}': expected proba_pair, slow_pair, lalign_id_pair, clustalw_pair, biphasic_pair, lalign_rs_s_pair, promo_pair, blastr_pair, fast_pair, exon_pair, exon2_pair, exon3_pair (or clean_slow_pair, hash_pair, subop1_pair, subop2_pair, cfast_pair, which the program refuses)."),
    };
}
