using System.ComponentModel;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-14: ClustalW 2.1's profile alignment commands (-PROFILE1 / -PROFILE2 with -PROFILE or -SEQUENCES); both tools delegate to SequenceAligner.AlignClustalWProfiles /
// AlignClustalWSequencesToProfile.
public partial class AlignmentTools
{
    [McpServerTool(Name = "clustalw_profile_align", Title = "MSA — ClustalW 2.1 Profile Alignment (-PROFILE)", ReadOnly = true)]
    [Description("ClustalW 2.1's profile-profile alignment of two existing alignments (clustalw -PROFILE1=a -PROFILE2=b -PROFILE, Clustal::profileAlign), ported step by step and verified row-for-row against the installed clustalw 2.1 on 600 DNA and protein sets (profiles produced by clustalw itself, option variants included): the rows of both profiles are kept (gaps '-' allowed; a profile may be one sequence), the identities inside each profile weight its sequences by the neighbour-joining tree of the profile (or UPGMA, clustering='upgma'), an initial alignment with equal weights measures the identities between the profiles, and the two profiles are aligned by Myers-Miller profile-profile alignment with position-specific gap penalties (weights = identity to the other profile x tree weight). Reproduced defects of the program: the tree of a profile 2 that does not start at sequence 1 is built from the identities of the FIRST sequences of the alignment, and the identities inside the profiles are truncated to integers. Parameters as multiple_align_clustalw (type 'dna' or 'protein'; gapOpen / gapExtend, transitionWeight, noWeights, negative; protein: matrix, noPreferredGaps, noHydrophilicGaps, variableRegionGaps, gapDist, hydrophilicResidues). Rows are returned: profile 1, then profile 2, in input order.")]
    public static MultipleAlignmentResultDto ClustalWProfileAlign(
        [Description("Rows of the first alignment (at least one; letters and '-').")] string[] profile1,
        [Description("Rows of the second alignment (at least one; letters and '-').")] string[] profile2,
        [Description("'dna' (default) or 'protein' (ClustalW -TYPE).")] string type = "dna",
        [Description("ClustalW -GAPOPEN (default 15 for DNA, 10 for protein).")] double? gapOpen = null,
        [Description("ClustalW -GAPEXT (default 6.66 for DNA, 0.2 for protein).")] double? gapExtend = null,
        [Description("DNA only, ClustalW -TRANSWEIGHT (default 0.5).")] float transitionWeight = 0.5f,
        [Description("ClustalW -NOWEIGHTS: no tree-derived sequence weights.")] bool noWeights = false,
        [Description("ClustalW -NEGATIVE: keep the matrix signed (protein: BLOSUM40 / PAM120 / GONNET250 instead of the series).")] bool negative = false,
        [Description("type='protein': ClustalW -MATRIX, 'gonnet' (default), 'blosum', 'pam' or 'id'.")] string? matrix = null,
        [Description("type='protein': ClustalW -NOPGAP.")] bool noPreferredGaps = false,
        [Description("type='protein': ClustalW -NOHGAP.")] bool noHydrophilicGaps = false,
        [Description("type='protein': ClustalW -NOVGAP (sic: switches the variable-region penalties ON).")] bool variableRegionGaps = false,
        [Description("type='protein': ClustalW -GAPDIST (default 4; 0 = off).")] int gapDist = 4,
        [Description("type='protein': ClustalW -HGAPRESIDUES (default GPSNDQEKR).")] string hydrophilicResidues = "GPSNDQEKR",
        [Description("ClustalW -CLUSTERING of the profile trees: 'nj' (default) or 'upgma'.")] string clustering = "nj",
        CancellationToken cancellationToken = default)
        => RunClustalWProfileTool(false, profile1, profile2, type, gapOpen, gapExtend, null, null, null, null, null, null, null, transitionWeight, noWeights, negative, matrix, null,
            noPreferredGaps, noHydrophilicGaps, variableRegionGaps, gapDist, hydrophilicResidues, clustering, cancellationToken);

    [McpServerTool(Name = "clustalw_sequences_to_profile", Title = "MSA — ClustalW 2.1 Sequences to Profile (-SEQUENCES)", ReadOnly = true)]
    [Description("ClustalW 2.1's sequence-to-profile alignment (clustalw -PROFILE1=a -PROFILE2=b -SEQUENCES, Clustal::sequencesAlignToProfile / MSA::seqsAlignToProfile), ported step by step and verified row-for-row against the installed clustalw 2.1 on 600 DNA and protein sets: the sequences are added to the existing alignment one at a time, the one with the highest identity to the aligned set first. Distances: the identities inside the profile (countid), the full (distance='clustalw') or fast k-tuple (distance='clustalw-quick', ktuple / window / topDiagonals / pairGap) pairwise stage between every other pair of rows (computed on the rows with their gaps; pwGapOpen / pwGapExtend / pwMatrix), then the neighbour-joining tree of all rows (or UPGMA) weights them; each sequence is aligned with weights = tree weight x identity to it and Myers-Miller profile alignment with position-specific gap penalties. Other parameters as multiple_align_clustalw. Rows are returned: profile, then sequences, in input order.")]
    public static MultipleAlignmentResultDto ClustalWSequencesToProfile(
        [Description("Rows of the existing alignment (at least one; letters and '-').")] string[] profile,
        [Description("The sequences to add (at least one; gaps '-' are kept as in the program).")] string[] sequences,
        [Description("'dna' (default) or 'protein' (ClustalW -TYPE).")] string type = "dna",
        [Description("ClustalW -GAPOPEN (default 15 for DNA, 10 for protein).")] double? gapOpen = null,
        [Description("ClustalW -GAPEXT (default 6.66 for DNA, 0.2 for protein).")] double? gapExtend = null,
        [Description("Pairwise stage: 'clustalw' (default; full DP) or 'clustalw-quick' (-QUICKTREE).")] string distance = "clustalw",
        [Description("ClustalW -PWGAPOPEN (default 15 for DNA, 10 for protein).")] double? pwGapOpen = null,
        [Description("ClustalW -PWGAPEXT (default 6.66 for DNA, 0.1 for protein).")] double? pwGapExtend = null,
        [Description("type='protein': ClustalW -PWMATRIX, 'gonnet' (default), 'blosum', 'pam' or 'id'.")] string? pwMatrix = null,
        [Description(KtupleDescription)] int? ktuple = null,
        [Description(WindowDescription)] int? window = null,
        [Description(TopDiagonalsDescription)] int? topDiagonals = null,
        [Description(PairGapDescription)] int? pairGap = null,
        [Description("DNA only, ClustalW -TRANSWEIGHT (default 0.5).")] float transitionWeight = 0.5f,
        [Description("ClustalW -NOWEIGHTS: no tree-derived sequence weights.")] bool noWeights = false,
        [Description("ClustalW -NEGATIVE: keep the matrix signed.")] bool negative = false,
        [Description("type='protein': ClustalW -MATRIX, 'gonnet' (default), 'blosum', 'pam' or 'id'.")] string? matrix = null,
        [Description("type='protein': ClustalW -NOPGAP.")] bool noPreferredGaps = false,
        [Description("type='protein': ClustalW -NOHGAP.")] bool noHydrophilicGaps = false,
        [Description("type='protein': ClustalW -NOVGAP (sic: switches the variable-region penalties ON).")] bool variableRegionGaps = false,
        [Description("type='protein': ClustalW -GAPDIST (default 4; 0 = off).")] int gapDist = 4,
        [Description("type='protein': ClustalW -HGAPRESIDUES (default GPSNDQEKR).")] string hydrophilicResidues = "GPSNDQEKR",
        [Description("ClustalW -CLUSTERING of the guide tree: 'nj' (default) or 'upgma'.")] string clustering = "nj",
        CancellationToken cancellationToken = default)
        => RunClustalWProfileTool(true, profile, sequences, type, gapOpen, gapExtend, distance, pwGapOpen, pwGapExtend, ktuple, window, topDiagonals, pairGap, transitionWeight, noWeights, negative, matrix, pwMatrix,
            noPreferredGaps, noHydrophilicGaps, variableRegionGaps, gapDist, hydrophilicResidues, clustering, cancellationToken);

    private static MultipleAlignmentResultDto RunClustalWProfileTool(
        bool sequencesMode, string[] first, string[] second, string type, double? gapOpen, double? gapExtend, string? distance, double? pwGapOpen, double? pwGapExtend, int? ktuple, int? window, int? topDiagonals,
        int? pairGap, float transitionWeight, bool noWeights, bool negative, string? matrix, string? pwMatrix, bool noPreferredGaps, bool noHydrophilicGaps, bool variableRegionGaps, int gapDist,
        string hydrophilicResidues, string clustering, CancellationToken cancellationToken)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        if (first == null || first.Length == 0)
            throw new ArgumentException(sequencesMode ? "profile must contain at least one sequence." : "profile1 must contain at least one sequence.", sequencesMode ? "profile" : "profile1");
        if (second == null || second.Length == 0)
            throw new ArgumentException(sequencesMode ? "sequences must contain at least one sequence." : "profile2 must contain at least one sequence.", sequencesMode ? "sequences" : "profile2");
        var mode = (distance ?? "clustalw").Trim().ToLowerInvariant() switch
        {
            "clustalw" or "full" => global::Seqeron.Genomics.Infrastructure.ClustalWPairwiseMode.Full,
            "clustalw-quick" or "clustalw_quick" or "quick" => global::Seqeron.Genomics.Infrastructure.ClustalWPairwiseMode.Quick,
            _ => throw new ArgumentException("distance must be 'clustalw' or 'clustalw-quick'.", nameof(distance)),
        };
        var clusteringMode = (clustering ?? "nj").Trim().ToLowerInvariant() switch
        {
            "nj" or "neighborjoining" or "neighbor-joining" => global::Seqeron.Genomics.Infrastructure.ClustalWClustering.NeighborJoining,
            "upgma" => global::Seqeron.Genomics.Infrastructure.ClustalWClustering.Upgma,
            _ => throw new ArgumentException("clustering must be 'nj' or 'upgma'.", nameof(clustering)),
        };
        var kind = (type ?? "dna").Trim().ToLowerInvariant();
        if (kind is not ("dna" or "protein"))
            throw new ArgumentException("type must be 'dna' or 'protein'.", nameof(type));

        global::Seqeron.Genomics.Infrastructure.MultipleAlignmentResult m;
        if (kind == "protein")
        {
            var protein = new global::Seqeron.Genomics.Infrastructure.ClustalWProteinOptions
            {
                Matrix = (matrix ?? "gonnet").Trim().ToLowerInvariant() switch
                {
                    "gonnet" => global::Seqeron.Genomics.Infrastructure.ClustalWProfileMatrix.Gonnet,
                    "blosum" => global::Seqeron.Genomics.Infrastructure.ClustalWProfileMatrix.Blosum,
                    "pam" => global::Seqeron.Genomics.Infrastructure.ClustalWProfileMatrix.Pam,
                    "id" or "identity" => global::Seqeron.Genomics.Infrastructure.ClustalWProfileMatrix.Identity,
                    _ => throw new ArgumentException("matrix must be 'gonnet', 'blosum', 'pam' or 'id'.", nameof(matrix)),
                },
                GapOpen = gapOpen ?? 10.0,
                GapExtend = gapExtend ?? 0.2,
                PairwiseMode = mode,
                Pairwise = ParseMsaOptions(null, null, null, ktuple, window, topDiagonals, pairGap, null, pwMatrix, pwGapOpen, pwGapExtend).ClustalWProtein,
                Clustering = clusteringMode,
                NoWeights = noWeights,
                NegativeMatrix = negative,
                NoPreferredGapPenalties = noPreferredGaps,
                NoHydrophilicGapPenalties = noHydrophilicGaps,
                VariableRegionGapPenalties = variableRegionGaps,
                GapSeparationDistance = gapDist,
                HydrophilicResidues = hydrophilicResidues ?? "GPSNDQEKR",
            };
            m = sequencesMode
                ? global::Seqeron.Genomics.Alignment.SequenceAligner.AlignClustalWSequencesToProfile(first, second, protein)
                : global::Seqeron.Genomics.Alignment.SequenceAligner.AlignClustalWProfiles(first, second, protein);
        }
        else
        {
            var options = new global::Seqeron.Genomics.Infrastructure.ClustalWOptions
            {
                GapOpen = gapOpen ?? 15.0,
                GapExtend = gapExtend ?? 6.66,
                PairwiseMode = mode,
                PairwiseGapOpen = pwGapOpen ?? 15.0,
                PairwiseGapExtend = pwGapExtend ?? 6.66,
                TransitionWeight = transitionWeight,
                Clustering = clusteringMode,
                NoWeights = noWeights,
                NegativeMatrix = negative,
                Quick = new global::Seqeron.Genomics.Infrastructure.ClustalWQuickOptions
                {
                    KTuple = ktuple ?? 2,
                    Window = window ?? 4,
                    TopDiagonals = topDiagonals ?? 4,
                    PairGap = pairGap ?? (ktuple.HasValue ? ktuple.Value + 4 : 5),
                },
            };
            m = sequencesMode
                ? global::Seqeron.Genomics.Alignment.SequenceAligner.AlignClustalWSequencesToProfile(first, second, options)
                : global::Seqeron.Genomics.Alignment.SequenceAligner.AlignClustalWProfiles(first, second, options);
        }

        return new MultipleAlignmentResultDto(m.AlignedSequences, m.Consensus, m.TotalScore, m.ExactTotalScore);
    }
}
