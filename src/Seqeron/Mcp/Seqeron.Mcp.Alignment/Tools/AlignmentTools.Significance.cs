using System.ComponentModel;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Mcp.Alignment.Tools;

public partial class AlignmentTools
{
    private const string SignificanceParamDescription = "true: add the BLAST+ statistics of the score (bitScore, eValue, pValue and the Karlin-Altschul parameters) to the result. Protein matrices of the NCBI tables (BLOSUM45 / 50 / 62 / 90, PAM30 / 70 / 250 with the gap costs the tables list, e.g. BLOSUM62 existence 11 extension 1 = gapOpen -11, gapExtend -1) get Spouge's finite-size E-value exactly as blastp -comp_based_stats 0; match / mismatch with the blastn reward / penalty pairs and gap costs of the tables (e.g. 2 / -3 with gapOpen -5, gapExtend -2) the length-adjusted Karlin-Altschul E-value as blastn -task blastn. Costs outside the tables are refused (default false).";
    private const string DbLengthParamDescription = "With significance: the length of the database the subject (sequence2) belongs to, for the E-value (default: the length of sequence2, a pairwise comparison, as blastp / blastn -subject).";
    private const string DbSequencesParamDescription = "With significance on nucleotides: the number of sequences of the database (default 1).";

    /// <summary>Adds the BLAST+ significance of the score of <paramref name="d"/> (a local alignment) when <paramref name="significance"/> is true.</summary>
    private static AlignmentResultDto WithSignificance(AlignmentResultDto d, bool significance, string? matrix, int match, int mismatch, double gapOpen, double gapExtend, int length1, int length2, long? dbLength, long dbSequences)
    {
        if (!significance)
            return d;
        if (dbLength is < 1)
            throw new ArgumentOutOfRangeException(nameof(dbLength), "dbLength must be at least 1.");
        if (dbSequences < 1)
            throw new ArgumentOutOfRangeException(nameof(dbSequences), "dbSequences must be at least 1.");
        if (d.Score <= 0 || d.AlignedSequence1.Length == 0)
            return d;
        double score = d.ExactScore != 0 ? d.ExactScore : d.Score;
        var sub = ResolveMatrix(matrix);
        AlignmentSignificanceResult r;
        if (sub is not null)
        {
            r = AlignmentSignificance.ForMatrix(score, new MatrixScoring(sub, gapOpen, gapExtend), length1, length2, dbLength ?? 0);
        }
        else
        {
            int existence = RequireInteger(-gapOpen, nameof(gapOpen)), extension = RequireInteger(-gapExtend, nameof(gapExtend));
            r = AlignmentSignificance.ForNucleotide(score, new ScoringMatrix(match, mismatch, -existence, -extension), length1, dbLength ?? length2, dbSequences);
        }

        return d with { Significance = new AlignmentSignificanceDto(r.BitScore, r.EValue, r.PValue, r.Method, r.Parameters.Lambda, r.Parameters.K, r.Parameters.H, r.LengthAdjustment, r.SearchSpace) };
    }

    [McpServerTool(Name = "local_align_zscore", Title = "Pairwise — Shuffle Z-score of a Local Alignment", ReadOnly = true)]
    [Description("The shuffle Z-score of the best local alignment score of two sequences (the idea of EMBOSS prss and FASTA ssearch): sequence2 is shuffled `shuffles` times (a seeded Fisher-Yates shuffle, reproducible), every shuffle is aligned locally to sequence1 with the same linear-gap scoring as local_align, and z = (score - mean) / standard deviation of the shuffled scores. A large z (above about 5-8) means the score is far above what the composition of the sequences gives by chance. Use local_align_affine / local_align / water_align with significance=true for the analytic bit score and E-value of the Karlin-Altschul statistics.")]
    public static AlignmentZScoreResult LocalAlignZScore(
        [Description("First sequence.")] string sequence1,
        [Description("Second sequence (the one that is shuffled).")] string sequence2,
        [Description("Score for matching base.")] int match = 1,
        [Description("Score for mismatching base.")] int mismatch = -1,
        [Description("Gap score charged per gap position (<= 0; the linear model of local_align).")] double gapExtend = -1,
        [Description("Number of shuffles (>= 2; default 1000).")] int shuffles = 1000,
        [Description("Seed of the shuffles (default 1).")] int seed = 1,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(sequence1, nameof(sequence1));
        RequireSequence(sequence2, nameof(sequence2));
        if (shuffles < 2)
            throw new ArgumentOutOfRangeException(nameof(shuffles), "shuffles must be at least 2.");
        var scoring = new ScoringMatrix(match, mismatch, 0, RequireInteger(gapExtend, nameof(gapExtend)));
        var (score, mean, sd, z) = AlignmentSignificance.ShuffleZScore(sequence1, sequence2, scoring, shuffles, seed, cancellationToken);
        return new AlignmentZScoreResult(score, mean, sd, z, shuffles, seed);
    }
}
