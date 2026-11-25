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
    [McpServerTool(Name = "multiple_align_progressive", Title = "MSA — Progressive Guide-Tree Alignment (Clustal-style)", ReadOnly = true)]
    [Description("Guide-tree progressive MSA (Feng & Doolittle 1987 / Clustal-style): the ClustalW pairwise distance of every pairwise alignment, d = 1 − (identical aligned residue pairs) / (length of the shorter sequence) (Thompson, Higgins & Gibson 1994; pairalign.c tmat = (100 − 100·count/MIN(len1, len2)) / 100), gives the distance matrix, a guide tree (UPGMA by default, guideTree='nj' = ClustalW's neighbour-joining tree) orders the merges, and profiles are aligned by profile–profile dynamic programming with sum-of-pairs column scores (once a gap, always a gap). Default linear gap model (gapOpen ignored); gapModel='affine' uses a Gotoh three-state profile DP (occupancy-scaled gapOpen/gapExtend, ClustalW/MUSCLE style) and reports the natural affine sum-of-pairs score. distance='clustalw' takes the guide-tree distances from ClustalW 2.1's own affine pairwise stage instead (percent identity of clustalw -ALIGN). Returns rows, majority consensus and sum-of-pairs score. PROTEIN / MATRIX: type='protein' (or matrix) uses the substitution matrix for the pairwise alignments, the profile column scores and the sum-of-pairs objective (fractional gaps allowed); rows equal an independent Python implementation on 238 sets. With type='protein', distance='clustalw' / 'clustalw-quick' use ClustalW 2.1's PROTEIN pairwise stages (GONNET250 / BLOSUM30 / PAM350 / identity matrix via pwMatrix, gap 10 / 0.1, k-tuple 1; exact ports verified on 4 204 + 3 774 pairs of clustalw -TYPE=PROTEIN -ALIGN / -QUICKTREE).")]
    public static MultipleAlignmentResultDto MultipleAlignProgressive(
        [Description("Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'.")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] double gapOpen = -2,
        [Description(GapExtendDescription)] double gapExtend = -1,
        [Description(DistanceModeDescription)] string distance = "linear",
        [Description(GapModelDescription)] string gapModel = "linear",
        [Description(GuideTreeDescription)] string guideTree = "upgma",
        [Description(KtupleDescription)] int? ktuple = null,
        [Description(WindowDescription)] int? window = null,
        [Description(TopDiagonalsDescription)] int? topDiagonals = null,
        [Description(PairGapDescription)] int? pairGap = null,
        [Description(MsaTypeDescription)] string type = "dna",
        [Description(MsaMatrixDescription)] string? matrix = null,
        [Description(PwMatrixDescription)] string? pwMatrix = null,
        [Description(PwGapOpenDescription)] double? pwGapOpen = null,
        [Description(PwGapExtendDescription)] double? pwGapExtend = null,
        [Description(FormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        [Description(FastDistanceCellsDescription)] long? fastDistanceCells = null,
        CancellationToken cancellationToken = default)
    {
        return FormatMultiple(MultipleAlignProgressiveCore(sequences, match, mismatch, gapOpen, gapExtend, distance, gapModel, guideTree, ktuple, window, topDiagonals, pairGap, type, matrix, pwMatrix, pwGapOpen, pwGapExtend, fastDistanceCells, cancellationToken), format, names);
    }
    [McpServerTool(Name = "multiple_align_iterative", Title = "MSA — Progressive + Iterative Refinement (MUSCLE-style)", ReadOnly = true)]
    [Description("Progressive MSA followed by MUSCLE-style tree-dependent restricted-partitioning refinement (Edgar 2004, stage 3): each guide-tree edge splits the alignment into two profiles that are re-aligned and kept only if the sum-of-pairs score strictly improves; repeated until convergence or maxIterations passes. The result's sum-of-pairs score is never below multiple_align_progressive's. Default linear gap model (gapOpen ignored); gapModel='affine' refines under the natural affine sum-of-pairs objective with the affine profile DP. distance='clustalw' selects ClustalW 2.1's pairwise stage for the guide tree and guideTree='nj' neighbour-joining (see multiple_align_progressive); the refinement partitions are the edges of that tree. PROTEIN / MATRIX: type='protein' (or matrix) refines under the matrix sum-of-pairs objective (see multiple_align_progressive); rows equal an independent Python implementation. With type='protein', distance='clustalw' / 'clustalw-quick' use ClustalW 2.1's PROTEIN pairwise stages (GONNET250 / BLOSUM30 / PAM350 / identity matrix via pwMatrix, gap 10 / 0.1, k-tuple 1; exact ports verified on 4 204 + 3 774 pairs of clustalw -TYPE=PROTEIN -ALIGN / -QUICKTREE).")]
    public static MultipleAlignmentResultDto MultipleAlignIterative(
        [Description("Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'.")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] double gapOpen = -2,
        [Description(GapExtendDescription)] double gapExtend = -1,
        [Description("Upper bound on full refinement passes (>= 1).")] int maxIterations = 16,
        [Description(DistanceModeDescription)] string distance = "linear",
        [Description(GapModelDescription)] string gapModel = "linear",
        [Description(GuideTreeDescription)] string guideTree = "upgma",
        [Description(KtupleDescription)] int? ktuple = null,
        [Description(WindowDescription)] int? window = null,
        [Description(TopDiagonalsDescription)] int? topDiagonals = null,
        [Description(PairGapDescription)] int? pairGap = null,
        [Description(MsaTypeDescription)] string type = "dna",
        [Description(MsaMatrixDescription)] string? matrix = null,
        [Description(PwMatrixDescription)] string? pwMatrix = null,
        [Description(PwGapOpenDescription)] double? pwGapOpen = null,
        [Description(PwGapExtendDescription)] double? pwGapExtend = null,
        [Description(FormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        [Description(FastDistanceCellsDescription)] long? fastDistanceCells = null,
        CancellationToken cancellationToken = default)
    {
        return FormatMultiple(MultipleAlignIterativeCore(sequences, match, mismatch, gapOpen, gapExtend, maxIterations, distance, gapModel, guideTree, ktuple, window, topDiagonals, pairGap, type, matrix, pwMatrix, pwGapOpen, pwGapExtend, fastDistanceCells, cancellationToken), format, names);
    }
    [McpServerTool(Name = "multiple_align_consistency", Title = "MSA — Consistency-Based Alignment (T-Coffee)", ReadOnly = true)]
    [Description("Consistency-based MSA (T-Coffee; Notredame, Higgins & Heringa 2000): a primary library of global (Needleman–Wunsch) and local (Smith–Waterman) pairwise alignments weighted by percent identity, extended by the triplet consistency transformation (direct weight + Σ min(W1, W2) over third sequences), then progressively aligned along a guide tree (UPGMA by default; guideTree='nj' = neighbour-joining as in the T-Coffee paper) with the extended-library weights as position-specific scores and zero gap penalties (T-Coffee default). The scoring parameters seed the pairwise library (gapModel='affine': Gotoh global and Smith–Waterman–Gotoh local alignments with gapOpen + k·gapExtend) and define totalScore, the sum of pairs of the result (linear column score, or natural affine costs with gapModel='affine'). distance='clustalw' takes the guide-tree distances from ClustalW 2.1's own pairwise stage. consistencyEngine='tcoffee' swaps the progressive stage for the T-Coffee 13.41 program's (library relaxation + linked pair-wise sparse DP; exact for the same library and tree). PROTEIN / MATRIX: type='protein' (or matrix) builds the primary library from the matrix global + local pair alignments (weights = percent identity) and otherwise runs the same pipeline (consistencyEngine='tcoffee' too); rows equal an independent Python implementation, and the T-Coffee progressive stage equals the t_coffee program on 160 protein sets. With type='protein', distance='clustalw' / 'clustalw-quick' use ClustalW 2.1's PROTEIN pairwise stages (GONNET250 / BLOSUM30 / PAM350 / identity matrix via pwMatrix, gap 10 / 0.1, k-tuple 1; exact ports verified on 4 204 + 3 774 pairs of clustalw -TYPE=PROTEIN -ALIGN / -QUICKTREE).")]
    public static MultipleAlignmentResultDto MultipleAlignConsistency(
        [Description("Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'.")] string[] sequences,
        [Description("Score for matching base (pairwise library and reported score).")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] double gapOpen = -2,
        [Description(GapExtendDescription)] double gapExtend = -1,
        [Description(DistanceModeDescription)] string distance = "linear",
        [Description(GapModelDescription)] string gapModel = "linear",
        [Description(GuideTreeDescription)] string guideTree = "upgma",
        [Description(KtupleDescription)] int? ktuple = null,
        [Description(WindowDescription)] int? window = null,
        [Description(TopDiagonalsDescription)] int? topDiagonals = null,
        [Description(PairGapDescription)] int? pairGap = null,
        [Description(ConsistencyEngineDescription)] string consistencyEngine = "paper",
        [Description(MsaTypeDescription)] string type = "dna",
        [Description(MsaMatrixDescription)] string? matrix = null,
        [Description(PwMatrixDescription)] string? pwMatrix = null,
        [Description(PwGapOpenDescription)] double? pwGapOpen = null,
        [Description(PwGapExtendDescription)] double? pwGapExtend = null,
        [Description("T-Coffee Lalign library: the number of non-intersecting (Waterman-Eggert) local alignments per sequence pair that enter the primary library, each weighted by its own percent identity (lalign_id_pair, -lalign_n_top; default 1 = the single best Smith-Waterman alignment of the published formula). Paper engine and TCoffeeLinked engine.")] int primaryLocalAlignments = 1,
        [Description("Sequence redundancy weights (T-Coffee -seq_weight): 'none' (default) or 'tcoffee' (compute_t_coffee_weight: w = 2 / (1 + sum (10 x identity)^3) rescaled to mean 1; an extended library weight is multiplied by the smaller of the two sequence weights).")] string sequenceWeights = "none",
        [Description(FormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        [Description(FastDistanceCellsDescription)] long? fastDistanceCells = null,
        CancellationToken cancellationToken = default)
    {
        return FormatMultiple(MultipleAlignConsistencyCore(sequences, match, mismatch, gapOpen, gapExtend, distance, gapModel, guideTree, ktuple, window, topDiagonals, pairGap, consistencyEngine, type, matrix, pwMatrix, pwGapOpen, pwGapExtend, primaryLocalAlignments, sequenceWeights, fastDistanceCells, cancellationToken), format, names);
    }
    [McpServerTool(Name = "multiple_align_clustalw", Title = "MSA — ClustalW 2.1 Multiple Alignment (exact port)", ReadOnly = true)]
    [Description("ClustalW 2.1's multiple alignment of DNA (clustalw -TYPE=DNA -ALIGN) or PROTEIN (type='protein', clustalw -TYPE=PROTEIN -ALIGN), ported step by step and verified row-for-row against the installed clustalw 2.1 (1 000 DNA sets and 1 000 protein sets, option variants included): pairwise distances (distance='clustalw' = the full affine Smith-Waterman-Gotoh stage, 'clustalw-quick' = the fast k-tuple stage of -QUICKTREE with ktuple / window / topDiagonals / pairGap), ClustalW's neighbour-joining guide tree (rooted by Tree::reRoot), sequence weights derived from the tree branch lengths (Thompson, Higgins & Gibson 1994; noWeights = equal weights), delayed divergent sequences (maxDiv, default 30: sequences with at most that percent identity to every other sequence are aligned last), and profile-profile alignment with position-specific gap penalties (positions that already hold gaps get a reduced opening penalty and half the extension penalty, end gaps are free) and Myers-Miller linear-space alignment. DNA: the IUB weight matrix shifted to non-negative scores (negative = keep it signed) with transitions weighted by transitionWeight (clustalw -TRANSWEIGHT, default 0.5); gapOpen / gapExtend = -GAPOPEN 15 / -GAPEXT 6.66. PROTEIN: the matrix series of ClustalW's profile stage (matrix = -MATRIX gonnet (default) / blosum / pam / id; the table of the series is chosen per profile pair by its percent identity, e.g. GONNET80 / 250 / 120 / 350 / 160, BLOSUM80 / 62 / 45 / 30, PAM20 / 60 / 120 / 350), a gap opening scaled by the matrix average score and the logarithm of the profile lengths (gapOpen / gapExtend = -GAPOPEN 10 / -GAPEXT 0.2), the residue-specific gap penalties of Pascarella & Argos 1992 (noPreferredGaps = -NOPGAP), the hydrophilic-run penalties (noHydrophilicGaps = -NOHGAP; hydrophilicResidues = -HGAPRESIDUES, default GPSNDQEKR, at most 9 letters), the variable-region penalties (variableRegionGaps = clustalw's -NOVGAP, sic: it switches them ON), the gap-separation distance (gapDist = -GAPDIST, default 4; 0 = off) and the pairwise stage of F62 (pwMatrix = -PWMATRIX gonnet / blosum / pam / id, pwGapOpen 10 / pwGapExtend 0.1, ktuple 1-2). User matrices (type='protein'): matrixFile = clustalw -MATRIX=file (one matrix file, or a CLUSTAL_SERIES file whose entries 'MATRIX llimit ulimit file' are looked up in matrixFiles; the entry whose percent-identity range contains the profile pair's identity is used, with the gap-penalty scale 0.5 + (pcid - llimit) / ((ulimit - llimit) * 2)) and pwMatrixFile = clustalw -PWMATRIX=file (verified on 240 + 800 clustalw -MATRIX / -PWMATRIX runs with floating-point / integer / unlabelled / tab-separated / header-subset files; only 2 sets, both -NOVGAP, hit the gdist residual below). Known residual: with variableRegionGaps and a two-sequence profile of more than 60 % identity ClustalW reads an uninitialised C++ member (gdist), so on ~0.2 % of such runs the heap decides; the port emulates the usual value. END-GAP OPTIONS (both types): useEndGaps = UserParameters::useEndGaps and endGapPenalties = UserParameters::endGapPenalties (reachable in clustalw only through its interactive menu / a patched build; verified on 340 sets x 4 modes of the patched clustalw 2.1 binary): with endGapPenalties ClustalW 2.1 writes NO output (exit status 0) when the score of a first-pass profile alignment is negative - this tool then fails with an explicit 'abandons the alignment' error. NUCLEOTIDE MATRICES (type='dna'): dnaMatrix / pwDnaMatrix = -DNAMATRIX / -PWDNAMATRIX IUB (default) or CLUSTALW; dnaMatrixFile / pwDnaMatrixFile = -DNAMATRIX=file / -PWDNAMATRIX=file (text of a matrix file; verified on 360 sets of the installed clustalw: the command line reads and validates the file but never uses it for DNA - a defect of CommandLineParser, userDnaMatrixUse='commandline' - while the interactive menu uses it, userDnaMatrixUse='applied', the default); transitionWeight = -TRANSWEIGHT also rewrites the cells of the user matrix at the letters of header positions 0 / 4, 2 / 11 and 2 / 12. Rows are returned in input order. match / mismatch / scoreGapOpen / scoreGapExtend (DNA) or scoreMatrix / scoreGapOpen / scoreGapExtend (protein; default BLOSUM62, open -9.8, extend -0.2) define only the reported totalScore (natural affine sum of pairs, Altschul 1989).")]
    public static MultipleAlignmentResultDto MultipleAlignClustalW(
        [Description("Sequences to align (at least one): DNA (A, C, G, T / IUB codes) by default, amino acids (letters A-Z except J) with type='protein'.")] string[] sequences,
        [Description("ClustalW -GAPOPEN: gap-opening penalty of the profile alignment (default 15 for DNA, 10 for protein; a penalty, positive).")] double? gapOpen = null,
        [Description("ClustalW -GAPEXT: gap-extension penalty of the profile alignment (default 6.66 for DNA, 0.2 for protein).")] double? gapExtend = null,
        [Description("Pairwise stage that gives the guide-tree distances: 'clustalw' (default; full DP, clustalw -ALIGN) or 'clustalw-quick' (fast k-tuple mode, clustalw -QUICKTREE).")] string distance = "clustalw",
        [Description("ClustalW -PWGAPOPEN: pairwise gap-opening penalty of the full pairwise stage (default 15 for DNA, 10 for protein).")] double? pwGapOpen = null,
        [Description("ClustalW -PWGAPEXT: pairwise gap-extension penalty of the full pairwise stage (default 6.66 for DNA, 0.1 for protein).")] double? pwGapExtend = null,
        [Description("DNA only, ClustalW -TRANSWEIGHT: weight of transitions (A-G, C-T/U) in the profile matrix (default 0.5; 0 = as mismatches, 1 = as matches).")] float transitionWeight = 0.5f,
        [Description("ClustalW -MAXDIV (default 30): a sequence whose highest percent identity to any other sequence is at most this value is delayed (aligned last, after the rest).")] int maxDiv = 30,
        [Description("ClustalW -NOWEIGHTS: no tree-derived sequence weights.")] bool noWeights = false,
        [Description("ClustalW -NEGATIVE: keep the substitution matrix signed instead of shifting it to non-negative scores (protein: one matrix BLOSUM40 / PAM120 / GONNET250 instead of the series).")] bool negative = false,
        [Description(KtupleDescription)] int? ktuple = null,
        [Description(WindowDescription)] int? window = null,
        [Description(TopDiagonalsDescription)] int? topDiagonals = null,
        [Description(PairGapDescription)] int? pairGap = null,
        [Description("Score for matching base in the reported totalScore only (DNA).")] int match = 1,
        [Description("Score for mismatching base in the reported totalScore only (DNA).")] int mismatch = -1,
        [Description("Gap-open score of the reported totalScore only (negative = penalty; default -2 for DNA, -9.8 for protein).")] double? scoreGapOpen = null,
        [Description("Gap-extend score of the reported totalScore only (negative = penalty; default -1 for DNA, -0.2 for protein).")] double? scoreGapExtend = null,
        [Description("'dna' (default), 'protein' (ClustalW -TYPE) or 'auto' (ClustalW's guess when -TYPE is omitted: decided by the SECOND sequence only, DNA when >= 85 % of its non-gap characters are A C G T U N; a single sequence is protein).")] string type = "dna",
        [Description("type='protein': ClustalW -MATRIX, the matrix series of the profile alignment: 'gonnet' (default), 'blosum', 'pam' or 'id'.")] string? matrix = null,
        [Description("type='protein': ClustalW -PWMATRIX, the pairwise matrix: 'gonnet' (default; GONNET250), 'blosum' (BLOSUM30), 'pam' (PAM350) or 'id'.")] string? pwMatrix = null,
        [Description("type='protein': ClustalW -NOPGAP: switch off the residue-specific (Pascarella & Argos) gap penalties.")] bool noPreferredGaps = false,
        [Description("type='protein': ClustalW -NOHGAP: switch off the hydrophilic-run gap penalties.")] bool noHydrophilicGaps = false,
        [Description("type='protein': ClustalW -NOVGAP (sic): switch the variable-region gap penalties on (clustalw's default is off).")] bool variableRegionGaps = false,
        [Description("type='protein': ClustalW -GAPDIST: gap separation distance (default 4; 0 = off): columns this close to an existing gap get a higher opening penalty.")] int gapDist = 4,
        [Description("type='protein': ClustalW -HGAPRESIDUES: residues of a hydrophilic stretch (default GPSNDQEKR; at most the first 9 letters are used).")] string hydrophilicResidues = "GPSNDQEKR",
        [Description("type='protein': substitution matrix of the reported totalScore only: BLOSUM62 (default), PAM250 or any EMBOSS protein matrix (EBLOSUM45, EPAM120, ...).")] string? scoreMatrix = null,
        [Description("type='protein': ClustalW -MATRIX=file, the TEXT of a user matrix file for the profile alignment (NCBI / EMBOSS layout: '#' comments, a header line of residue letters, one row per residue with an optional leading label, one extra '*' row, integer or floating-point values (a row whose first value contains '.' is scaled by 10); only the lower triangle is read) or of a CLUSTAL_SERIES file (lines 'MATRIX llimit ulimit file', file names resolved through matrixFiles). Replaces 'matrix'.")] string? matrixFile = null,
        [Description("type='protein': the matrix files named by a CLUSTAL_SERIES matrixFile: file name -> text of the matrix file.")] Dictionary<string, string>? matrixFiles = null,
        [Description("type='protein': ClustalW -PWMATRIX=file, the TEXT of a user matrix file for the pairwise (distance) stage (same layout as matrixFile, a single matrix). Replaces 'pwMatrix'.")] string? pwMatrixFile = null,
        [Description("End-gap option of the profile stage, UserParameters::useEndGaps (the interactive menu item 'Toggle End Gap Separation'; the command-line flag -ENDGAPS only ever writes false): gaps at the sequence ends count as gaps in the position-specific penalties. Default false.")] bool useEndGaps = false,
        [Description("End-gap option of the profile stage, UserParameters::endGapPenalties: the end gaps of a profile alignment are penalised instead of free. Default false. With it ClustalW 2.1 abandons the alignment silently (no output file) when the score of a first-pass profile alignment is negative; this tool then fails with a 'ClustalW 2.1 abandons the alignment ...' error (ClustalWAbortedException).")] bool endGapPenalties = false,
        [Description("type='dna': ClustalW -DNAMATRIX built-in: 'iub' (default) or 'clustalw' (ClustalW 1.6's identity matrix, gap penalties scaled by 0.66). Ignored when dnaMatrixFile is given.")] string? dnaMatrix = null,
        [Description("type='dna': ClustalW -PWDNAMATRIX built-in for the pairwise stage: 'iub' (default) or 'clustalw' (opening x 0.6667, extension x 0.751). Ignored when pwDnaMatrixFile is given.")] string? pwDnaMatrix = null,
        [Description("type='dna': ClustalW -DNAMATRIX=file, the TEXT of a nucleotide matrix file for the profile stage (same layout as matrixFile, a single matrix; the transition weights of getProfileAlignMatrix index the matrix with the header letters at positions 0, 2, 4, 11 and 12 - give at least 13 letters). How it is used: see userDnaMatrixUse.")] string? dnaMatrixFile = null,
        [Description("type='dna': ClustalW -PWDNAMATRIX=file, the TEXT of a nucleotide matrix file for the full pairwise stage (a single matrix, x100, never shifted). See userDnaMatrixUse.")] string? pwDnaMatrixFile = null,
        [Description("type='dna', with dnaMatrixFile / pwDnaMatrixFile: 'applied' (default; what ClustalW's interactive menu does - the matrix is used) or 'commandline' (what clustalw -DNAMATRIX=file -PWDNAMATRIX=file does: the file is read and validated but registered under the PROTEIN slot - a defect of CommandLineParser - so the DNA stages keep the IUB / built-in matrix).")] string userDnaMatrixUse = "applied",
        [Description("ClustalW -CLUSTERING: 'nj' (default; neighbour joining) or 'upgma' (the program's own rooted UPGMA tree: node-height weights and join order; verified on 600 / 600 sets (rows and guide-tree text) of the unmodified Debian source built at -O0 - the installed Debian clustalw 2.1 and an -O2 build of the same source segmentation-fault on every UPGMA run, undefined behaviour of the program).")] string clustering = "nj",
        [Description("ClustalW -ITERATION: 'none' (default), 'tree' (every profile-profile step of the progressive alignment is followed by the 'remove first' realignment of its sequences) or 'alignment' (the same realignment once on the finished alignment): each sequence in turn is cut out, the others' all-gap columns are removed, the sequence is re-aligned to the profile of the rest with tree-derived weights, and the best realignment by the program's ObjectiveScore (BLOSUM45 x 100 pair sums + gap penalties, kept when LOWER) replaces the alignment. Verified on 700 clustalw runs (rows and the printed 'FINAL score'); with clustering='upgma' the profile trees of getWeightsForProfileAlignUPGMA are ported from the source but not yet verified against a binary.")] string iteration = "none",
        [Description("ClustalW -NUMITER (default 3; positive): the maximum number of iterations of 'iteration'. The program keeps the best realignment of the first pass, so more iterations change nothing once one finds no improvement.")] int numIterations = 3,
        [Description(FormatParamDescription)] string? format = null,
        [Description(NamesParamDescription)] string[]? names = null,
        CancellationToken cancellationToken = default)
    {
        return FormatMultiple(MultipleAlignClustalWCore(sequences, gapOpen, gapExtend, distance, pwGapOpen, pwGapExtend, transitionWeight, maxDiv, noWeights, negative, ktuple, window, topDiagonals, pairGap, match, mismatch, scoreGapOpen, scoreGapExtend, type, matrix, pwMatrix, noPreferredGaps, noHydrophilicGaps, variableRegionGaps, gapDist, hydrophilicResidues, scoreMatrix, matrixFile, matrixFiles, pwMatrixFile, useEndGaps, endGapPenalties, dnaMatrix, pwDnaMatrix, dnaMatrixFile, pwDnaMatrixFile, userDnaMatrixUse, clustering, iteration, numIterations, cancellationToken), format, names);
    }
    [McpServerTool(Name = "multiple_align_tcoffee", Title = "MSA — T-Coffee 13.41 Program (primary libraries, guide tree, progressive stage)", ReadOnly = true)]
    [Description("The multiple alignment of the T-Coffee 13.41 program (t_coffee; Notredame, Higgins & Heringa 2000) reproduced end to end: its own primary library built by the program's pair-wise methods - 'proba_pair' (the program's default: ProbCons pair-HMM posteriors > 0.01, at most 4 x min(length) pairs per pair of sequences), 'slow_pair' (Myers-Miller global alignment, dna_idmat 10 / 0 with gap open -10 / extend -1 for DNA, BLOSUM62 for proteins) and 'lalign_id_pair' (Huang-Miller SIM, the lalignTopAlignments best non-intersecting local alignments), several methods merged the way the program's forked library merge does (which is not their union) - then library relaxation, the program's guide tree (k-tuple distances, neighbour joining with the trifurcating root turned into a binary tree; fast neighbour joining from 100 sequences) and the sparse linked-pair-wise progressive stage. Rows equal t_coffee -seq s.fa [-method m1,m2] -output fasta_aln on 1 500+ committed DNA and protein sets (libraries pair by pair on 800+ more), for any -n_core. The type (DNA when more than 98 % of the letters are A G C T U N R Y M K S W H B V D, otherwise protein) selects the matrices, penalties and models. RNA input (a U) uses the program's separate RNA pair-HMM (5 matrix types, 2 insert state pairs, get_tot_prob). Further command-line options of the program are accepted in 'options' (a t_coffee flag string, e.g. '-tree_mode upgma -distance_matrix_mode fast -seq_weight t_coffee -relax_lib 2 -filter_lib 20 -gapopen -50 -gapext -2 -nomatch 5 -outorder aligned'): -relax_lib N, -filter_lib T, -tree_mode nj|upgma (cwph, kmeans and unknown names are the program's refusals: TCoffeeAbortedException), -distance_matrix_mode ktup|fast|idscore|no|slow|cscore|clustalw|geometric|very_fast|random (random = rand() % 100 per pair from the glibc generator seeded by randomSeed, default 10 = the program's DEBUG_SRAND run; the five before it are computed from the RELAXED library / the FASTA DP between the relaxation and the re-weighting, as the program does; -quicktree = very_fast; aln, ktup2 and unknown names are refused as the program does; random is not reproducible: the program's rand is seeded differently on every run), -distance_matrix_sim_mode (idmat_sim1/2/3, gap1, logid), -seq_weight t_coffee|no|@w1,w2,... (weights in input order), -weight no, -gapopen, -gapext, -cosmetic_penalty, -nomatch, -do_normalise, -inorder aligned|input, -outorder input|aligned, -lalign_n_top, -method; flags without influence on the rows (-evaluate_mode, -shrink_lib, -n_core, ...) are accepted and ignored; any other flag is rejected by name. The progressive stage runs other dynamic programs with -dp_mode gotoh_pair_wise (ssec_pair_wise, gotoh_pair_wise_test) | myers_miller_pair_wise (default) | fasta_pair_wise | cfasta_pair_wise | very_fast_pair_wise on the library column scores of -extend_mode very_fast_triplet (default) | fast_triplet | slow_triplet | triplet | relative_triplet | rna0, with -max_n_pair, -tg_mode, -ndiag, -diag_threshold; the dp modes and extend modes the program refuses or dies on raise TCoffeeAbortedException, the ones that only run in the program but are not ported (viterbiL / viterbiD, co_pair_wise = Clustal Omega, linked_pair_wise_collapse, sim_pair_wise_lalign; the rna contact libraries) raise an error naming them. Every option is verified row for row against the installed t_coffee on 150+ committed sets (4 300+ sets in all). 'names' = the sequence names: the program sorts the sequences by name (strcmp order) before building the library and the tree (-inorder aligned, its default) and restores the input order on output, which is the 'name-prefix quirk' (s1 / s10 trees differ); without names the names are taken to sort in input order. -matrix NAME (38 matrices of the program: blosum30..80mt, pam120..350mt, md_*, idmat, est_idmat, ...). Further pair-wise methods: 'fast_pair' (FASTA heuristic: banded Gotoh DP on the k-tuple diagonals), 'exon_pair' / 'exon2_pair' / 'exon3_pair' (Myers-Miller alignment of the genomic _G_ templates given in genomicTemplates, matrix exon2mt, expressed in the residues of the sequences), 'biphasic_pair' (biphasic_pair_wise: two affine gap models, opening x 2 / extension / 2 and opening / 2 / extension x 2), 'lalign_rs_s_pair' (lalign_id_pair on every pair INCLUDING each sequence against itself), 'promo_pair' / 'blastr_pair' (slow_pair on the dinucleotide-extended alphabet, matrices promoter_tf1 / blosumR, gap open -60 / -20; the program's shared sequence array is changed lossily by the extension, which is reproduced); -distance_matrix_sim_mode also accepts the group names of make_group_aa (simple, mafft, clustalw, polarity, vasiliky, make_all, an unknown name = no group) and idscore (a gap in the first row of an aligned pair under -distance_matrix_mode fast makes the program segfault, exit status 11). A set on which the program itself dies - lalign_id_pair with a very small gap open (e.g. 'lalign_id_pair@EP@GOP@-3': SIM's edit-script array overflows its heap chunk), the methods clean_slow_pair / hash_pair / subop1_pair / subop2_pair / cfast_pair, the idscore segfault - raises TCoffeeAbortedException (ExitStatus, Reason). Methods that run an external program (T-Coffee's side of the interface is ported, the external program is NOT reimplemented): 'muscle_pair' (MUSCLE 3.8.1551 on every pair, the library of its alignment), 'TMalign_pair' / 'mustang_pair' / 'sap_pair' and -mode 3dcoffee (structure alignment of the PDB structures given as structureTemplates: the alignment TM-align / MUSTANG wrote is given in structureAlignerOutputs; sap_pair = TM-align because SAP is not installed; a pair without two usable structures falls back to proba_pair), 'blastp_msa' / 'tblastx_msa' (the HSPs of the all-against-all BLAST / tblastx search given as blastXml become the library, weight = truncated HSP percent identity); 'dali_pair' and 'fugue_pair' stop as in the installed program (DaliLite / fugueali are not installed, TCoffeeAbortedException), 'mafft_pair' is the MAFFT 7.505 port of the pair stage (multiple_align_mafft; 460 / 460 sets of the program equal) and -mode expresso is reported as not available (expresso downloads structures from the RCSB). More dynamic programs: -dp_mode gotoh_pair_wise_lgp (position specific gap opening / closing penalties from the gaps of the two profiles), ktup_pair_wise (k-tuple score only: the rows stay as they are, unpadded), precomputed_pair_wise (stops unless the first rows of two groups have the same length), gotoh_pair_wise_sw (Smith-Waterman form: every merge keeps only the columns of the best local alignment, rows of different length; 10 % of the committed sets are undefined behaviour of the program) test_pair_wise (a debugging routine: rows padded with gaps; a second group shorter than n columns is a segmentation fault, longer rows than the longest sequence overflow the program's row buffers), procoffee_pair_wise (sparse linked DP, older diagonal cap, groups swapped) and viterbi_pair_wise (Viterbi path of the protein pair-HMM on the residue frequency profiles); extend modes pc, mixt, g_coffee, g_coffee_quadruplets, quadruplet, test and rna / rna1-4 (no contact library = slow_triplet). totalScore is the natural affine sum of pairs (EDNAFULL / BLOSUM62), not the program's objective (0 and an empty consensus when the rows have different lengths).")]
    public static MultipleAlignmentResultDto MultipleAlignTCoffee(
        [Description("Sequences to align (at least one, none empty): letters only (the program removes gaps; here they are rejected).")] string[] sequences,
        [Description("Primary-library methods in command-line order (t_coffee -method m1,m2): 'proba_pair' (default; = plain t_coffee seqs.fa), 'slow_pair', 'lalign_id_pair', 'clustalw_pair' (the ClustalW 2.1 pair alignments), 'biphasic_pair', 'lalign_rs_s_pair', 'promo_pair', 'blastr_pair', 'fast_pair', 'exon_pair', 'exon2_pair', 'exon3_pair', each optionally with the program's parameter overrides, e.g. 'slow_pair@EP@GOP@-30@GEP@-3@MATRIX@blosum50mt'; 'clean_slow_pair', 'hash_pair', 'subop1_pair', 'subop2_pair', 'cfast_pair' are refused by the program itself and raise TCoffeeAbortedException.")] string[]? methods = null,
        [Description("t_coffee -lalign_n_top (default 10): alignments reported by lalign_id_pair for every pair of sequences.")] int lalignTopAlignments = 10,
        [Description("t_coffee -relax_lib: true = 1 (default, library relaxation), false = 0.")] bool relaxLibrary = true,
        [Description("'auto' (default; the program's decision), 'dna' or 'protein' (forces the matrices, penalties and models).")] string type = "auto",
        [Description("Further t_coffee flags, e.g. '-tree_mode upgma -relax_lib 2 -seq_weight t_coffee -gapext -2' (see the tool description for the supported list); applied after the other parameters.")] string? options = null,
        [Description("Sequence names (one per sequence, in input order): the program sorts the sequences by name before computing anything (-inorder aligned); omitted = names that sort in input order.")] string[]? names = null,
        [Description("The genomic (_G_) template sequence of every input sequence for the methods exon_pair / exon2_pair / exon3_pair (t_coffee -template_file: one entry per sequence in input order, null / empty = the sequence has no template; only pairs whose two sequences have one are aligned). Omitted = no templates.")] string?[]? genomicTemplates = null,
        [Description("The sequence of the PDB structure (a _P_ template) of every input sequence for the structure methods TMalign_pair / mustang_pair / sap_pair (one entry per sequence in input order, null / empty = no structure). A structure is attached only when it is alike enough (identity > 35 %, coverage > 50 %); a pair without two structures is aligned by proba_pair, as the program does.")] string?[]? structureTemplates = null,
        [Description("What the EXTERNAL structure aligner (TM-align / MUSTANG, not reimplemented) wrote for pairs of structures: method 'TMalign_pair', 'mustang_pair' or 'sap_pair', the input indices of the two sequences in the order the program calls it (names sorted), and the output (TMalign: its standard output, the last four lines matter; mustang: the content of results.afasta; sap_pair: the content of the file <template1>_<template2>.sap_results = the raw output of SAP 1.1.3). A missing pair = the aligner failed (no library entries). Without any 'sap_pair' item the method sap_pair uses the TMalign entries (SAP not installed: the program switches to TM-align); with at least one 'sap_pair' item SAP is installed and the program's own sap_pair glue reads the SAP files (a pair of two usable structures without a file stops the program, as it did in the oracle).")] StructureAlignerOutputInput[]? structureAlignerOutputs = null,
        [Description("The XML (blastall -p blastp -m7 / NCBI BLAST+ -outfmt 5) of the search of every sequence against a local BLAST database of all of them, one Iteration per sequence in order, for the methods blastp_msa and tblastx_msa (the BLAST search itself is the external program and is not reimplemented). The hits are matched by the names (see names; default q000, q001, ...).")] string? blastXml = null,
        [Description(FormatParamDescription)] string? format = null,
        [Description("The srand seed of -distance_matrix_mode random (rand() % 100 per pair, glibc rand). The program seeds from /dev/urandom (two runs differ) unless the environment variable DEBUG_SRAND is set, which is srand(10): the default 10 is that deterministic run (verified on 310 sets). Ignored by every other mode.")] int? randomSeed = null,
        CancellationToken cancellationToken = default)
    {
        return FormatMultiple(MultipleAlignTCoffeeCore(sequences, methods, lalignTopAlignments, relaxLibrary, type, options, names, genomicTemplates, structureTemplates, structureAlignerOutputs, blastXml, cancellationToken, randomSeed), format, names);
    }

    private static IReadOnlyDictionary<global::Seqeron.Genomics.TCoffee.TCoffeePairMethod, IReadOnlyDictionary<(int First, int Second), string>>? BuildStructureOutputs(StructureAlignerOutputInput[]? outputs)
    {
        if (outputs == null) return null;
        var map = new Dictionary<global::Seqeron.Genomics.TCoffee.TCoffeePairMethod, Dictionary<(int First, int Second), string>>();
        foreach (var o in outputs)
        {
            if (o is null)
                throw new ArgumentNullException(nameof(outputs), "structureAlignerOutputs contains a null element.");
            var method = (o.Method ?? "").Trim().ToLowerInvariant() switch
            {
                "tmalign_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.TMalignPair,
                "mustang_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.MustangPair,
                "sap_pair" => global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.SapPair,
                _ => throw new ArgumentException("structureAlignerOutputs: method must be 'TMalign_pair', 'mustang_pair' or 'sap_pair'.", nameof(outputs)),
            };
            if (!map.TryGetValue(method, out var byPair)) map[method] = byPair = new();
            byPair[(o.First, o.Second)] = o.Output ?? string.Empty;
        }

        return map.ToDictionary(kv => kv.Key, kv => (IReadOnlyDictionary<(int First, int Second), string>)kv.Value);
    }
    private static global::Seqeron.Genomics.ClustalW.ClustalWDnaMatrix ParseClustalWDnaMatrix(string? name, string parameter) => (name ?? "iub").Trim().ToLowerInvariant() switch
    {
        "iub" => global::Seqeron.Genomics.ClustalW.ClustalWDnaMatrix.Iub,
        "clustalw" => global::Seqeron.Genomics.ClustalW.ClustalWDnaMatrix.Clustalw,
        _ => throw new ArgumentException("must be 'iub' or 'clustalw' (a user matrix is given as the text of its file).", parameter),
    };
    private static global::Seqeron.Genomics.ClustalW.ClustalWUserMatrix ParseClustalWMatrix(string text, string parameter)
    {
        try
        {
            return global::Seqeron.Genomics.ClustalW.ClustalWUserMatrix.Parse(text, parameter);
        }
        catch (FormatException ex)
        {
            throw new ArgumentException(ex.Message, parameter, ex);
        }
    }

    private static MultipleAlignmentResultDto MultipleAlignProgressiveCore(
        [Description("Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'.")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] double gapOpen = -2,
        [Description(GapExtendDescription)] double gapExtend = -1,
        [Description(DistanceModeDescription)] string distance = "linear",
        [Description(GapModelDescription)] string gapModel = "linear",
        [Description(GuideTreeDescription)] string guideTree = "upgma",
        [Description(KtupleDescription)] int? ktuple = null,
        [Description(WindowDescription)] int? window = null,
        [Description(TopDiagonalsDescription)] int? topDiagonals = null,
        [Description(PairGapDescription)] int? pairGap = null,
        [Description(MsaTypeDescription)] string type = "dna",
        [Description(MsaMatrixDescription)] string? matrix = null,
        [Description(PwMatrixDescription)] string? pwMatrix = null,
        [Description(PwGapOpenDescription)] double? pwGapOpen = null,
        [Description(PwGapExtendDescription)] double? pwGapExtend = null,
        [Description(FastDistanceCellsDescription)] long? fastDistanceCells = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        return RunMsa(sequences, match, mismatch, gapOpen, gapExtend, type, matrix,
            (seqs, scoring) => global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(seqs, scoring, ParseMsaOptions(distance, gapModel, guideTree, ktuple, window, topDiagonals, pairGap, null, pwMatrix, pwGapOpen, pwGapExtend, fastDistanceCells: fastDistanceCells)),
            (seqs, scoring) => global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(seqs, scoring, ParseMsaOptions(distance, gapModel, guideTree, ktuple, window, topDiagonals, pairGap, null, pwMatrix, pwGapOpen, pwGapExtend, fastDistanceCells: fastDistanceCells)));
    }
    private static MultipleAlignmentResultDto MultipleAlignIterativeCore(
        [Description("Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'.")] string[] sequences,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] double gapOpen = -2,
        [Description(GapExtendDescription)] double gapExtend = -1,
        [Description("Upper bound on full refinement passes (>= 1).")] int maxIterations = 16,
        [Description(DistanceModeDescription)] string distance = "linear",
        [Description(GapModelDescription)] string gapModel = "linear",
        [Description(GuideTreeDescription)] string guideTree = "upgma",
        [Description(KtupleDescription)] int? ktuple = null,
        [Description(WindowDescription)] int? window = null,
        [Description(TopDiagonalsDescription)] int? topDiagonals = null,
        [Description(PairGapDescription)] int? pairGap = null,
        [Description(MsaTypeDescription)] string type = "dna",
        [Description(MsaMatrixDescription)] string? matrix = null,
        [Description(PwMatrixDescription)] string? pwMatrix = null,
        [Description(PwGapOpenDescription)] double? pwGapOpen = null,
        [Description(PwGapExtendDescription)] double? pwGapExtend = null,
        [Description(FastDistanceCellsDescription)] long? fastDistanceCells = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        if (maxIterations < 1)
            throw new ArgumentOutOfRangeException(nameof(maxIterations), "maxIterations must be >= 1.");

        var options = ParseMsaOptions(distance, gapModel, guideTree, ktuple, window, topDiagonals, pairGap, null, pwMatrix, pwGapOpen, pwGapExtend, fastDistanceCells: fastDistanceCells);
        return RunMsa(sequences, match, mismatch, gapOpen, gapExtend, type, matrix,
            (seqs, scoring) => global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignIterative(seqs, scoring, options, maxIterations),
            (seqs, scoring) => global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignIterative(seqs, scoring, options, maxIterations));
    }
    private static MultipleAlignmentResultDto MultipleAlignConsistencyCore(
        [Description("Sequences to align (at least one): DNA (A, C, G, T) by default, amino acids with type='protein'.")] string[] sequences,
        [Description("Score for matching base (pairwise library and reported score).")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description(GapOpenDescription)] double gapOpen = -2,
        [Description(GapExtendDescription)] double gapExtend = -1,
        [Description(DistanceModeDescription)] string distance = "linear",
        [Description(GapModelDescription)] string gapModel = "linear",
        [Description(GuideTreeDescription)] string guideTree = "upgma",
        [Description(KtupleDescription)] int? ktuple = null,
        [Description(WindowDescription)] int? window = null,
        [Description(TopDiagonalsDescription)] int? topDiagonals = null,
        [Description(PairGapDescription)] int? pairGap = null,
        [Description(ConsistencyEngineDescription)] string consistencyEngine = "paper",
        [Description(MsaTypeDescription)] string type = "dna",
        [Description(MsaMatrixDescription)] string? matrix = null,
        [Description(PwMatrixDescription)] string? pwMatrix = null,
        [Description(PwGapOpenDescription)] double? pwGapOpen = null,
        [Description(PwGapExtendDescription)] double? pwGapExtend = null,
        [Description("T-Coffee Lalign library: the number of non-intersecting (Waterman-Eggert) local alignments per sequence pair that enter the primary library, each weighted by its own percent identity (lalign_id_pair, -lalign_n_top; default 1 = the single best Smith-Waterman alignment of the published formula). Paper engine and TCoffeeLinked engine.")] int primaryLocalAlignments = 1,
        [Description("Sequence redundancy weights (T-Coffee -seq_weight): 'none' (default) or 'tcoffee' (compute_t_coffee_weight: w = 2 / (1 + sum (10 x identity)^3) rescaled to mean 1; an extended library weight is multiplied by the smaller of the two sequence weights).")] string sequenceWeights = "none",
        [Description(FastDistanceCellsDescription)] long? fastDistanceCells = null,
        CancellationToken cancellationToken = default)
    {
        var weightMode = (sequenceWeights ?? "none").Trim().ToLowerInvariant() switch
        {
            "none" => global::Seqeron.Genomics.TCoffee.MsaSequenceWeights.None,
            "tcoffee" or "t_coffee" => global::Seqeron.Genomics.TCoffee.MsaSequenceWeights.TCoffee,
            _ => throw new ArgumentException("sequenceWeights must be 'none' or 'tcoffee'.", nameof(sequenceWeights)),
        };
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        return RunMsa(sequences, match, mismatch, gapOpen, gapExtend, type, matrix,
            (seqs, scoring) => global::Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignConsistency(seqs, scoring, ParseMsaOptions(distance, gapModel, guideTree, ktuple, window, topDiagonals, pairGap, consistencyEngine, pwMatrix, pwGapOpen, pwGapExtend, fastDistanceCells: fastDistanceCells) with { PrimaryLocalAlignments = primaryLocalAlignments, SequenceWeights = weightMode }),
            (seqs, scoring) => global::Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignConsistency(seqs, scoring, ParseMsaOptions(distance, gapModel, guideTree, ktuple, window, topDiagonals, pairGap, consistencyEngine, pwMatrix, pwGapOpen, pwGapExtend, fastDistanceCells: fastDistanceCells) with { PrimaryLocalAlignments = primaryLocalAlignments, SequenceWeights = weightMode }));
    }
    private static MultipleAlignmentResultDto MultipleAlignClustalWCore(
        [Description("Sequences to align (at least one): DNA (A, C, G, T / IUB codes) by default, amino acids (letters A-Z except J) with type='protein'.")] string[] sequences,
        [Description("ClustalW -GAPOPEN: gap-opening penalty of the profile alignment (default 15 for DNA, 10 for protein; a penalty, positive).")] double? gapOpen = null,
        [Description("ClustalW -GAPEXT: gap-extension penalty of the profile alignment (default 6.66 for DNA, 0.2 for protein).")] double? gapExtend = null,
        [Description("Pairwise stage that gives the guide-tree distances: 'clustalw' (default; full DP, clustalw -ALIGN) or 'clustalw-quick' (fast k-tuple mode, clustalw -QUICKTREE).")] string distance = "clustalw",
        [Description("ClustalW -PWGAPOPEN: pairwise gap-opening penalty of the full pairwise stage (default 15 for DNA, 10 for protein).")] double? pwGapOpen = null,
        [Description("ClustalW -PWGAPEXT: pairwise gap-extension penalty of the full pairwise stage (default 6.66 for DNA, 0.1 for protein).")] double? pwGapExtend = null,
        [Description("DNA only, ClustalW -TRANSWEIGHT: weight of transitions (A-G, C-T/U) in the profile matrix (default 0.5; 0 = as mismatches, 1 = as matches).")] float transitionWeight = 0.5f,
        [Description("ClustalW -MAXDIV (default 30): a sequence whose highest percent identity to any other sequence is at most this value is delayed (aligned last, after the rest).")] int maxDiv = 30,
        [Description("ClustalW -NOWEIGHTS: no tree-derived sequence weights.")] bool noWeights = false,
        [Description("ClustalW -NEGATIVE: keep the substitution matrix signed instead of shifting it to non-negative scores (protein: one matrix BLOSUM40 / PAM120 / GONNET250 instead of the series).")] bool negative = false,
        [Description(KtupleDescription)] int? ktuple = null,
        [Description(WindowDescription)] int? window = null,
        [Description(TopDiagonalsDescription)] int? topDiagonals = null,
        [Description(PairGapDescription)] int? pairGap = null,
        [Description("Score for matching base in the reported totalScore only (DNA).")] int match = 1,
        [Description("Score for mismatching base in the reported totalScore only (DNA).")] int mismatch = -1,
        [Description("Gap-open score of the reported totalScore only (negative = penalty; default -2 for DNA, -9.8 for protein).")] double? scoreGapOpen = null,
        [Description("Gap-extend score of the reported totalScore only (negative = penalty; default -1 for DNA, -0.2 for protein).")] double? scoreGapExtend = null,
        [Description("'dna' (default), 'protein' (ClustalW -TYPE) or 'auto' (ClustalW's guess when -TYPE is omitted: decided by the SECOND sequence only, DNA when >= 85 % of its non-gap characters are A C G T U N; a single sequence is protein).")] string type = "dna",
        [Description("type='protein': ClustalW -MATRIX, the matrix series of the profile alignment: 'gonnet' (default), 'blosum', 'pam' or 'id'.")] string? matrix = null,
        [Description("type='protein': ClustalW -PWMATRIX, the pairwise matrix: 'gonnet' (default; GONNET250), 'blosum' (BLOSUM30), 'pam' (PAM350) or 'id'.")] string? pwMatrix = null,
        [Description("type='protein': ClustalW -NOPGAP: switch off the residue-specific (Pascarella & Argos) gap penalties.")] bool noPreferredGaps = false,
        [Description("type='protein': ClustalW -NOHGAP: switch off the hydrophilic-run gap penalties.")] bool noHydrophilicGaps = false,
        [Description("type='protein': ClustalW -NOVGAP (sic): switch the variable-region gap penalties on (clustalw's default is off).")] bool variableRegionGaps = false,
        [Description("type='protein': ClustalW -GAPDIST: gap separation distance (default 4; 0 = off): columns this close to an existing gap get a higher opening penalty.")] int gapDist = 4,
        [Description("type='protein': ClustalW -HGAPRESIDUES: residues of a hydrophilic stretch (default GPSNDQEKR; at most the first 9 letters are used).")] string hydrophilicResidues = "GPSNDQEKR",
        [Description("type='protein': substitution matrix of the reported totalScore only: BLOSUM62 (default), PAM250 or any EMBOSS protein matrix (EBLOSUM45, EPAM120, ...).")] string? scoreMatrix = null,
        [Description("type='protein': ClustalW -MATRIX=file, the TEXT of a user matrix file for the profile alignment (NCBI / EMBOSS layout: '#' comments, a header line of residue letters, one row per residue with an optional leading label, one extra '*' row, integer or floating-point values (a row whose first value contains '.' is scaled by 10); only the lower triangle is read) or of a CLUSTAL_SERIES file (lines 'MATRIX llimit ulimit file', file names resolved through matrixFiles). Replaces 'matrix'.")] string? matrixFile = null,
        [Description("type='protein': the matrix files named by a CLUSTAL_SERIES matrixFile: file name -> text of the matrix file.")] Dictionary<string, string>? matrixFiles = null,
        [Description("type='protein': ClustalW -PWMATRIX=file, the TEXT of a user matrix file for the pairwise (distance) stage (same layout as matrixFile, a single matrix). Replaces 'pwMatrix'.")] string? pwMatrixFile = null,
        [Description("End-gap option of the profile stage, UserParameters::useEndGaps (the interactive menu item 'Toggle End Gap Separation'; the command-line flag -ENDGAPS only ever writes false): gaps at the sequence ends count as gaps in the position-specific penalties. Default false.")] bool useEndGaps = false,
        [Description("End-gap option of the profile stage, UserParameters::endGapPenalties: the end gaps of a profile alignment are penalised instead of free. Default false. With it ClustalW 2.1 abandons the alignment silently (no output file) when the score of a first-pass profile alignment is negative; this tool then fails with a 'ClustalW 2.1 abandons the alignment ...' error (ClustalWAbortedException).")] bool endGapPenalties = false,
        [Description("type='dna': ClustalW -DNAMATRIX built-in: 'iub' (default) or 'clustalw' (ClustalW 1.6's identity matrix, gap penalties scaled by 0.66). Ignored when dnaMatrixFile is given.")] string? dnaMatrix = null,
        [Description("type='dna': ClustalW -PWDNAMATRIX built-in for the pairwise stage: 'iub' (default) or 'clustalw' (opening x 0.6667, extension x 0.751). Ignored when pwDnaMatrixFile is given.")] string? pwDnaMatrix = null,
        [Description("type='dna': ClustalW -DNAMATRIX=file, the TEXT of a nucleotide matrix file for the profile stage (same layout as matrixFile, a single matrix; the transition weights of getProfileAlignMatrix index the matrix with the header letters at positions 0, 2, 4, 11 and 12 - give at least 13 letters). How it is used: see userDnaMatrixUse.")] string? dnaMatrixFile = null,
        [Description("type='dna': ClustalW -PWDNAMATRIX=file, the TEXT of a nucleotide matrix file for the full pairwise stage (a single matrix, x100, never shifted). See userDnaMatrixUse.")] string? pwDnaMatrixFile = null,
        [Description("type='dna', with dnaMatrixFile / pwDnaMatrixFile: 'applied' (default; what ClustalW's interactive menu does - the matrix is used) or 'commandline' (what clustalw -DNAMATRIX=file -PWDNAMATRIX=file does: the file is read and validated but registered under the PROTEIN slot - a defect of CommandLineParser - so the DNA stages keep the IUB / built-in matrix).")] string userDnaMatrixUse = "applied",
        [Description("ClustalW -CLUSTERING: 'nj' (default; neighbour joining) or 'upgma' (the program's own rooted UPGMA tree: node-height weights and join order; verified on 600 / 600 sets (rows and guide-tree text) of the unmodified Debian source built at -O0 - the installed Debian clustalw 2.1 and an -O2 build of the same source segmentation-fault on every UPGMA run, undefined behaviour of the program).")] string clustering = "nj",
        [Description("ClustalW -ITERATION: 'none' (default), 'tree' (every profile-profile step of the progressive alignment is followed by the 'remove first' realignment of its sequences) or 'alignment' (the same realignment once on the finished alignment): each sequence in turn is cut out, the others' all-gap columns are removed, the sequence is re-aligned to the profile of the rest with tree-derived weights, and the best realignment by the program's ObjectiveScore (BLOSUM45 x 100 pair sums + gap penalties, kept when LOWER) replaces the alignment. Verified on 700 clustalw runs (rows and the printed 'FINAL score'); with clustering='upgma' the profile trees of getWeightsForProfileAlignUPGMA are ported from the source but not yet verified against a binary.")] string iteration = "none",
        [Description("ClustalW -NUMITER (default 3; positive): the maximum number of iterations of 'iteration'. The program keeps the best realignment of the first pass, so more iterations change nothing once one finds no improvement.")] int numIterations = 3,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        var iterationMode = (iteration ?? "none").Trim().ToLowerInvariant() switch
        {
            "none" => global::Seqeron.Genomics.ClustalW.ClustalWIteration.None,
            "tree" => global::Seqeron.Genomics.ClustalW.ClustalWIteration.Tree,
            "alignment" => global::Seqeron.Genomics.ClustalW.ClustalWIteration.Alignment,
            _ => throw new ArgumentException("iteration must be 'none', 'tree' or 'alignment'.", nameof(iteration)),
        };
        var mode = (distance ?? "clustalw").Trim().ToLowerInvariant() switch
        {
            "clustalw" or "clustalw_affine" or "clustalwaffine" or "full" => global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMode.Full,
            "clustalw-quick" or "clustalw_quick" or "clustalwquick" or "quick" => global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMode.Quick,
            _ => throw new ArgumentException("distance must be 'clustalw' or 'clustalw-quick'.", nameof(distance)),
        };
        var clusteringMode = (clustering ?? "nj").Trim().ToLowerInvariant() switch
        {
            "nj" or "neighborjoining" or "neighbor-joining" => global::Seqeron.Genomics.ClustalW.ClustalWClustering.NeighborJoining,
            "upgma" => global::Seqeron.Genomics.ClustalW.ClustalWClustering.Upgma,
            _ => throw new ArgumentException("clustering must be 'nj' or 'upgma'.", nameof(clustering)),
        };
        var kind = (type ?? "dna").Trim().ToLowerInvariant();
        if (kind is not ("dna" or "protein" or "auto"))
            throw new ArgumentException("type must be 'dna', 'protein' or 'auto'.", nameof(type));
        if (kind == "auto")
            kind = global::Seqeron.Genomics.ClustalW.ClustalWAligner.GuessClustalWSequenceType(sequences) == global::Seqeron.Genomics.ClustalW.ClustalWSequenceType.Dna ? "dna" : "protein";

        if (kind == "protein")
        {
            var protein = new global::Seqeron.Genomics.ClustalW.ClustalWProteinOptions
            {
                Matrix = (matrix ?? "gonnet").Trim().ToLowerInvariant() switch
                {
                    "gonnet" => global::Seqeron.Genomics.ClustalW.ClustalWProfileMatrix.Gonnet,
                    "blosum" => global::Seqeron.Genomics.ClustalW.ClustalWProfileMatrix.Blosum,
                    "pam" => global::Seqeron.Genomics.ClustalW.ClustalWProfileMatrix.Pam,
                    "id" or "identity" => global::Seqeron.Genomics.ClustalW.ClustalWProfileMatrix.Identity,
                    _ => throw new ArgumentException("matrix must be 'gonnet', 'blosum', 'pam' or 'id'.", nameof(matrix)),
                },
                GapOpen = gapOpen ?? 10.0,
                GapExtend = gapExtend ?? 0.2,
                PairwiseMode = mode,
                Pairwise = ParseMsaOptions(null, null, null, ktuple, window, topDiagonals, pairGap, null, pwMatrix, pwGapOpen, pwGapExtend).ClustalWProtein,
                DivergenceCutoff = maxDiv,
                Clustering = clusteringMode,
                Iteration = iterationMode,
                NumIterations = numIterations,
                NoWeights = noWeights,
                NegativeMatrix = negative,
                NoPreferredGapPenalties = noPreferredGaps,
                NoHydrophilicGapPenalties = noHydrophilicGaps,
                VariableRegionGapPenalties = variableRegionGaps,
                GapSeparationDistance = gapDist,
                HydrophilicResidues = hydrophilicResidues ?? "GPSNDQEKR",
                UseEndGaps = useEndGaps,
                EndGapPenalties = endGapPenalties,
            };
            if (pwMatrixFile is not null)
                protein = protein with { Pairwise = protein.Pairwise with { UserMatrix = ParseClustalWMatrix(pwMatrixFile, nameof(pwMatrixFile)) } };
            if (matrixFile is not null)
            {
                global::Seqeron.Genomics.ClustalW.ClustalWMatrixFile parsed;
                try
                {
                    parsed = global::Seqeron.Genomics.ClustalW.ClustalWMatrixFile.Parse(matrixFile,
                        name => matrixFiles is not null && matrixFiles.TryGetValue(name, out var text)
                            ? text
                            : throw new FormatException($"The series names the matrix file '{name}', which is not in matrixFiles."),
                        "matrixFile");
                }
                catch (FormatException ex)
                {
                    throw new ArgumentException(ex.Message, nameof(matrixFile), ex);
                }

                protein = protein with { UserMatrix = parsed.Matrix, UserMatrixSeries = parsed.Series };
            }
            var scoreSub = global::Seqeron.Genomics.Infrastructure.SubstitutionMatrix.TryGet(scoreMatrix ?? "BLOSUM62", out var scoreFound) && !scoreFound.IsNucleotide
                ? scoreFound
                : throw new ArgumentException("scoreMatrix must be a protein matrix name: BLOSUM62 (default), PAM250 or any other EMBOSS protein matrix (EBLOSUM45, EPAM120, ...).", nameof(scoreMatrix));
            var m = global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(
                sequences, protein, new global::Seqeron.Genomics.Infrastructure.MatrixScoring(scoreSub, scoreGapOpen ?? -9.8, scoreGapExtend ?? -0.2));
            return new MultipleAlignmentResultDto(m.AlignedSequences, m.Consensus, m.TotalScore, m.ExactTotalScore);
        }

        var options = new global::Seqeron.Genomics.ClustalW.ClustalWOptions
        {
            GapOpen = gapOpen ?? 15.0,
            GapExtend = gapExtend ?? 6.66,
            PairwiseMode = mode,
            PairwiseGapOpen = pwGapOpen ?? 15.0,
            PairwiseGapExtend = pwGapExtend ?? 6.66,
            TransitionWeight = transitionWeight,
            DivergenceCutoff = maxDiv,
            Clustering = clusteringMode,
            Iteration = iterationMode,
            NumIterations = numIterations,
            NoWeights = noWeights,
            NegativeMatrix = negative,
            UseEndGaps = useEndGaps,
            EndGapPenalties = endGapPenalties,
            DnaMatrix = ParseClustalWDnaMatrix(dnaMatrix, nameof(dnaMatrix)),
            PairwiseDnaMatrix = ParseClustalWDnaMatrix(pwDnaMatrix, nameof(pwDnaMatrix)),
            UserDnaMatrix = dnaMatrixFile is null ? null : ParseClustalWMatrix(dnaMatrixFile, nameof(dnaMatrixFile)),
            PairwiseUserDnaMatrix = pwDnaMatrixFile is null ? null : ParseClustalWMatrix(pwDnaMatrixFile, nameof(pwDnaMatrixFile)),
            UserDnaMatrixUse = (userDnaMatrixUse ?? "applied").Trim().ToLowerInvariant() switch
            {
                "applied" => global::Seqeron.Genomics.ClustalW.ClustalWUserDnaMatrixUse.Applied,
                "commandline" or "command-line" or "command_line" => global::Seqeron.Genomics.ClustalW.ClustalWUserDnaMatrixUse.CommandLine,
                _ => throw new ArgumentException("userDnaMatrixUse must be 'applied' or 'commandline'.", nameof(userDnaMatrixUse)),
            },
            Quick = new global::Seqeron.Genomics.ClustalW.ClustalWQuickOptions
            {
                KTuple = ktuple ?? 2,
                Window = window ?? 4,
                TopDiagonals = topDiagonals ?? 4,
                PairGap = pairGap ?? (ktuple.HasValue ? ktuple.Value + 4 : 5),
            },
        };
        return RunMsa(sequences, match, mismatch, RequireInteger(scoreGapOpen ?? -2, nameof(scoreGapOpen)), RequireInteger(scoreGapExtend ?? -1, nameof(scoreGapExtend)),
            (seqs, scoring) => global::Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignClustalW(seqs, options, scoring));
    }
    private static MultipleAlignmentResultDto MultipleAlignTCoffeeCore(
        [Description("Sequences to align (at least one, none empty): letters only (the program removes gaps; here they are rejected).")] string[] sequences,
        [Description("Primary-library methods in command-line order (t_coffee -method m1,m2): 'proba_pair' (default; = plain t_coffee seqs.fa), 'slow_pair', 'lalign_id_pair', 'clustalw_pair' (the ClustalW 2.1 pair alignments), 'biphasic_pair', 'lalign_rs_s_pair', 'promo_pair', 'blastr_pair', 'fast_pair', 'exon_pair', 'exon2_pair', 'exon3_pair', each optionally with the program's parameter overrides, e.g. 'slow_pair@EP@GOP@-30@GEP@-3@MATRIX@blosum50mt'; 'clean_slow_pair', 'hash_pair', 'subop1_pair', 'subop2_pair', 'cfast_pair' are refused by the program itself and raise TCoffeeAbortedException.")] string[]? methods = null,
        [Description("t_coffee -lalign_n_top (default 10): alignments reported by lalign_id_pair for every pair of sequences.")] int lalignTopAlignments = 10,
        [Description("t_coffee -relax_lib: true = 1 (default, library relaxation), false = 0.")] bool relaxLibrary = true,
        [Description("'auto' (default; the program's decision), 'dna' or 'protein' (forces the matrices, penalties and models).")] string type = "auto",
        [Description("Further t_coffee flags, e.g. '-tree_mode upgma -relax_lib 2 -seq_weight t_coffee -gapext -2' (see the tool description for the supported list); applied after the other parameters.")] string? options = null,
        [Description("Sequence names (one per sequence, in input order): the program sorts the sequences by name before computing anything (-inorder aligned); omitted = names that sort in input order.")] string[]? names = null,
        [Description("The genomic (_G_) template sequence of every input sequence for the methods exon_pair / exon2_pair / exon3_pair (t_coffee -template_file: one entry per sequence in input order, null / empty = the sequence has no template; only pairs whose two sequences have one are aligned). Omitted = no templates.")] string?[]? genomicTemplates = null,
        [Description("The sequence of the PDB structure (a _P_ template) of every input sequence for the structure methods TMalign_pair / mustang_pair / sap_pair (one entry per sequence in input order, null / empty = no structure). A structure is attached only when it is alike enough (identity > 35 %, coverage > 50 %); a pair without two structures is aligned by proba_pair, as the program does.")] string?[]? structureTemplates = null,
        [Description("What the EXTERNAL structure aligner (TM-align / MUSTANG, not reimplemented) wrote for pairs of structures: method 'TMalign_pair', 'mustang_pair' or 'sap_pair', the input indices of the two sequences in the order the program calls it (names sorted), and the output (TMalign: its standard output, the last four lines matter; mustang: the content of results.afasta; sap_pair: the content of the file <template1>_<template2>.sap_results = the raw output of SAP 1.1.3). A missing pair = the aligner failed (no library entries). Without any 'sap_pair' item the method sap_pair uses the TMalign entries (SAP not installed: the program switches to TM-align); with at least one 'sap_pair' item SAP is installed and the program's own sap_pair glue reads the SAP files (a pair of two usable structures without a file stops the program, as it did in the oracle).")] StructureAlignerOutputInput[]? structureAlignerOutputs = null,
        [Description("The XML (blastall -p blastp -m7 / NCBI BLAST+ -outfmt 5) of the search of every sequence against a local BLAST database of all of them, one Iteration per sequence in order, for the methods blastp_msa and tblastx_msa (the BLAST search itself is the external program and is not reimplemented). The hits are matched by the names (see names; default q000, q001, ...).")] string? blastXml = null,
        CancellationToken cancellationToken = default,
        int? randomSeed = null)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        ArgumentNullException.ThrowIfNull(sequences);
        if (genomicTemplates != null && genomicTemplates.Length != sequences.Length)
            throw new ArgumentException("genomicTemplates must have one entry per sequence.", nameof(genomicTemplates));
        if (structureTemplates != null && structureTemplates.Length != sequences.Length)
            throw new ArgumentException("structureTemplates must have one entry per sequence.", nameof(structureTemplates));

        if (names != null && names.Length != sequences.Length)
            throw new ArgumentException("names must have one entry per sequence.", nameof(names));
        var tcOptions = new global::Seqeron.Genomics.TCoffee.TCoffeeProgramOptions
        {
            LalignTopAlignments = lalignTopAlignments,
            RelaxLibrary = relaxLibrary,
            Names = names,
            GenomicTemplates = genomicTemplates?.Select(t => string.IsNullOrEmpty(t) ? null : t).ToArray(),
            StructureTemplates = structureTemplates?.Select(t => string.IsNullOrEmpty(t) ? null : t).ToArray(),
            StructureAlignerOutputs = BuildStructureOutputs(structureAlignerOutputs),
            BlastXml = blastXml,
            Protein = (type ?? "auto").Trim().ToLowerInvariant() switch
            {
                "auto" => null,
                "dna" => false,
                "protein" => true,
                _ => throw new ArgumentException("type must be 'auto', 'dna' or 'protein'.", nameof(type)),
            },
        };
        if (lalignTopAlignments < 1)
            throw new ArgumentException("lalignTopAlignments must be at least 1.", nameof(lalignTopAlignments));

        try
        {
            // methods go through the flag parser: 'proba_pair', 'slow_pair', 'lalign_id_pair', 'clustalw_pair' and the name@EP@GOP@-30@GEP@-3@MATRIX@blosum50mt form
            var methodList = methods is { Length: > 0 } ? methods.Select(m => (m ?? "").Trim().ToLowerInvariant()).ToArray() : new[] { "proba_pair" };
            tcOptions = tcOptions.WithCommandLine("-method " + string.Join(",", methodList));
            if (!string.IsNullOrWhiteSpace(options)) tcOptions = tcOptions.WithCommandLine(options);
            if (randomSeed.HasValue) tcOptions = tcOptions with { RandomSeed = unchecked((uint)randomSeed.Value) };
            var m2 = global::Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignTCoffeeProgram(sequences, tcOptions);
            return new MultipleAlignmentResultDto(m2.AlignedSequences, m2.Consensus, m2.TotalScore, m2.ExactTotalScore);
        }
        catch (NotSupportedException ex)
        {
            throw new ArgumentException(ex.Message, nameof(sequences), ex);
        }
    }

    private const string FastDistanceCellsDescription =
        "Work limit of the all-pairs distance stage when distance='linear': when the sum of length1 x length2 over all pairs of sequences is above this many dynamic-programming cells the distance of a pair is the k-tuple estimate " +
        "(MUSCLE's k-mer similarity as identity, d = 1 - F^(1/k), k = 8 for nucleotides, 4 for amino acids; O(length1 + length2) per pair) instead of one global alignment per pair; 0 or less = always the alignments. Default 400000000 " +
        "(the alignments cost about 7.6 ns per cell: 3 s at the limit); pairs with a sequence shorter than 2k keep the alignment distance.";
    private const string DistanceModeDescription =
        "Guide-tree distance: 'linear' (default) = 1 - identities / min(length) of the repository's global alignment under the given scores (linear NW, or Gotoh when gapModel='affine'); " +
        "'clustalw' = the same formula on ClustalW 2.1's own pairwise stage (affine local Smith-Waterman-Gotoh alignment with its IUB matrix x100 and gap penalties 15 / 6.66, " +
        "the 'slow' full-DP mode of clustalw -ALIGN; type='protein': its PROTEIN stage with the GONNET250 / BLOSUM30 / PAM350 / identity matrix (pwMatrix) and gap penalties 10 / 0.1 (pwGapOpen / pwGapExtend); " +
        "the match/mismatch/gap parameters are then used only for the profile alignment); " +
        "'clustalw-quick' = ClustalW 2.1's fast k-tuple mode (clustalw -QUICKTREE, FastPairwiseAlign.cpp: best chain of exactly matching k-mers on the top diagonals, " +
        "d = (100 - 100*fragmentScore/min(length))/100; tuned by ktuple / window / topDiagonals / pairGap).";
    private const string KtupleDescription = "distance='clustalw-quick' only: k-tuple (word) size, 1-4 for DNA (clustalw -KTUPLE, default 2), 1-2 for protein (default 1).";
    private const string WindowDescription = "distance='clustalw-quick' only: diagonals scanned around each top diagonal (clustalw -WINDOW, default 4 for DNA, 5 for protein).";
    private const string TopDiagonalsDescription = "distance='clustalw-quick' only: number of best diagonals used (clustalw -TOPDIAGS, default 4 for DNA, 5 for protein).";
    private const string PairGapDescription = "distance='clustalw-quick' only: gap penalty of the fragment chain (clustalw -PAIRGAP; DNA default 5, or ktuple + 4 when ktuple is given; protein default 3, or ktuple + 3, as the clustalw command line does).";
    private const string PwMatrixDescription = "type='protein' with distance='clustalw': ClustalW's pairwise matrix (clustalw -PWMATRIX): 'gonnet' (default; GONNET250), 'blosum' (BLOSUM30), 'pam' (PAM350) or 'id' (identity).";
    private const string PwGapOpenDescription = "type='protein' with distance='clustalw': ClustalW -PWGAPOPEN (default 10; length-dependent: 2 * matrixAvg * (open + ln min(length)) in matrix units).";
    private const string PwGapExtendDescription = "type='protein' with distance='clustalw': ClustalW -PWGAPEXT (default 0.1).";
    private const string GuideTreeDescription =
        "Guide-tree method: 'upgma' (default; Sokal & Michener 1958, the Feng-Doolittle / CLUSTAL-1 choice) or 'nj' (neighbour-joining, Saitou & Nei 1987 - " +
        "ClustalW 2.1's default tree, T-Coffee/Clustal-style; the unrooted NJ tree is rooted and ordered exactly as clustalw does: Tree::reRoot, left subtree first). " +
        "With distance='clustalw' and guideTree='nj' the tree equals the .dnd file of clustalw -ALIGN. UPGMA is built by the canonical phylogenetics builder; NJ uses ClustalW's own neighbour-joining arithmetic (exact ties as in clustalw).";
    private const string ConsistencyEngineDescription =
        "Progressive stage: 'paper' (default) = the published formula (extended-library weight summed over every cross pair of two columns, full DP, zero gap penalty); " +
        "'tcoffee' = the T-Coffee 13.41 program's own progressive stage (Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignTCoffee: library relaxation, 'linked pair-wise' sparse DP over scored column pairs, weights 10 x percent identity). " +
        "With the same primary library and guide tree 'tcoffee' reproduces t_coffee's rows exactly (800 explicit-library sets); the library here is this repository's global + local alignments, " +
        "not the program's slow_pair / lalign_id_pair / proba_pair libraries, and the tree is guideTree (not the program's k-tuple NJ tree), so a plain 't_coffee seqs.fa' run is not reproduced end to end.";

    private static global::Seqeron.Genomics.TCoffee.ConsistencyMsaOptions ParseMsaOptions(
        string? distance, string? gapModel, string? guideTree = null,
        int? ktuple = null, int? window = null, int? topDiagonals = null, int? pairGap = null, string? consistencyEngine = null,
        string? pwMatrix = null, double? pwGapOpen = null, double? pwGapExtend = null, long? fastDistanceCells = null) => new()
    {
        FastDistanceCells = fastDistanceCells ?? new global::Seqeron.Genomics.ClustalW.ProgressiveMsaOptions().FastDistanceCells,
        ClustalWProtein = new global::Seqeron.Genomics.ClustalW.ClustalWProteinPairwiseOptions
        {
            Matrix = (pwMatrix ?? "gonnet").Trim().ToLowerInvariant() switch
            {
                "gonnet" or "gon250" => global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMatrix.Gonnet,
                "blosum" or "blosum30" => global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMatrix.Blosum,
                "pam" or "pam350" => global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMatrix.Pam,
                "id" or "identity" => global::Seqeron.Genomics.ClustalW.ClustalWPairwiseMatrix.Identity,
                _ => throw new ArgumentException("pwMatrix must be 'gonnet', 'blosum', 'pam' or 'id'.", nameof(pwMatrix)),
            },
            GapOpen = pwGapOpen ?? 10.0,
            GapExtend = pwGapExtend ?? 0.1,
            Quick = new global::Seqeron.Genomics.ClustalW.ClustalWQuickOptions
            {
                KTuple = ktuple ?? 1,
                Window = window ?? 5,
                TopDiagonals = topDiagonals ?? 5,
                // clustalw's command line (PROTEIN): -KTUPLE=k sets the pair gap to k + 3, -PAIRGAP overrides; default 3
                PairGap = pairGap ?? (ktuple.HasValue ? ktuple.Value + 3 : 3),
            },
        },
        ConsistencyEngine = (consistencyEngine ?? "paper").Trim().ToLowerInvariant() switch
        {
            "paper" => global::Seqeron.Genomics.TCoffee.MsaConsistencyEngine.PaperExtension,
            "tcoffee" or "t-coffee" or "tcoffee-linked" or "tcoffee_linked" => global::Seqeron.Genomics.TCoffee.MsaConsistencyEngine.TCoffeeLinked,
            _ => throw new ArgumentException("consistencyEngine must be 'paper' or 'tcoffee'.", nameof(consistencyEngine)),
        },
        QuickPairwise = new global::Seqeron.Genomics.ClustalW.ClustalWQuickOptions
        {
            KTuple = ktuple ?? 2,
            Window = window ?? 4,
            TopDiagonals = topDiagonals ?? 4,
            // clustalw's command line: -KTUPLE=k sets the pair gap to k + 4, -PAIRGAP overrides; without either the DNA default is 5
            PairGap = pairGap ?? (ktuple.HasValue ? ktuple.Value + 4 : 5),
        },
        GuideTree = (guideTree ?? "upgma").Trim().ToLowerInvariant() switch
        {
            "upgma" => global::Seqeron.Genomics.ClustalW.MsaGuideTree.Upgma,
            "nj" or "neighbor_joining" or "neighbour_joining" or "neighborjoining" => global::Seqeron.Genomics.ClustalW.MsaGuideTree.NeighborJoining,
            _ => throw new ArgumentException("guideTree must be 'upgma' or 'nj'.", nameof(guideTree)),
        },
        Distance = (distance ?? "linear").Trim().ToLowerInvariant() switch
        {
            "linear" => global::Seqeron.Genomics.ClustalW.MsaDistanceMode.Linear,
            "clustalw" or "clustalw_affine" or "clustalwaffine" => global::Seqeron.Genomics.ClustalW.MsaDistanceMode.ClustalWAffine,
            "clustalw-quick" or "clustalw_quick" or "clustalwquick" => global::Seqeron.Genomics.ClustalW.MsaDistanceMode.ClustalWQuick,
            _ => throw new ArgumentException("distance must be 'linear', 'clustalw' or 'clustalw-quick'.", nameof(distance)),
        },
        GapModel = (gapModel ?? "linear").Trim().ToLowerInvariant() switch
        {
            "linear" => global::Seqeron.Genomics.Infrastructure.MsaGapModel.Linear,
            "affine" => global::Seqeron.Genomics.Infrastructure.MsaGapModel.Affine,
            _ => throw new ArgumentException("gapModel must be 'linear' or 'affine'.", nameof(gapModel)),
        },
    };
}
