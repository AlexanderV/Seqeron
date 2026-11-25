// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Analysis.Tools;
using static Seqeron.Mcp.Analysis.Tools.AnalysisTools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

public partial class EmbossPairwiseAlignmentTools
{
    [McpServerTool(Name = "alignment_distance_matrix", Title = "Alignment — Pairwise Distance Matrix (distmat)", ReadOnly = true)]
    [Description("Pairwise distance matrix of a multiple alignment as EMBOSS 6.6.0 distmat (every correction and option, verified on 1 100 + 600 runs of the installed program, single-precision arithmetic as in distmat.c): matches = columns where both rows carry the same non-gap symbol, g = columns where either row has a gap; 'uncorrected' D = 1 - matches/(L - g + g*gapWeight); 'jukes-cantor' -b*ln(1 - D/b) with b = 3/4 (nucleotide) or 19/20 (protein); 'kimura' = Kimura 2-parameter (nucleotides) / Kimura protein -ln(1 - D - 0.2 D^2) (proteins); 'tamura', 'tajima-nei' and 'jin-nei' (Gamma distance; calculateA / parameterA) are nucleotide methods. ambiguous counts ambiguity codes (proteins: B/Z/J/X weights; nucleotides: distmat's ajBaseAlphacharCompare with its own defects, and a gap in the first sequence of a pair makes the program loop forever = tool error). position = 123 (all codon bases, default), 12, or 1 / 2 / 3 (nucleotides); tajima-nei and protein kimura with a position other than 123 read beyond the program's buffers (undefined) = tool error. Fractions, not percent; entries with a zero denominator or a saturated logarithm are null. Rows are read as EMBOSS reads an alignment (A4-3): letters and '*' are kept, '.' '~' '-' are gaps, '?' '#' '+' are the unknown residue, every other character (digits, blanks, '_', punctuation) is dropped.")]
    public static AlignmentDistanceMatrixResult AlignmentDistanceMatrix(
        [Description("At least two aligned sequences of equal length.")] string[] sequences,
        [Description("'uncorrected' (default), 'jukes-cantor', 'kimura', 'tamura', 'tajima-nei' or 'jin-nei'.")] string method = "uncorrected",
        [Description("'auto' (default; EMBOSS decides from the first row), 'nucleotide' or 'protein'.")] string residueType = "auto",
        [Description("distmat -gapweight (default 0; uncorrected and jukes-cantor).")] double gapWeight = 0,
        [Description("distmat -ambiguous (default false).")] bool ambiguous = false,
        [Description("distmat -position: 123 (default), 12, or for nucleotides 1, 2, 3.")] int position = 123,
        [Description("distmat -calculatea: calculate the Jin-Nei parameter a instead of using parameterA (default false).")] bool calculateA = false,
        [Description("distmat -parametera: the Jin-Nei Gamma parameter a (default 1.0).")] double parameterA = 1.0,
        [Description("EMBOSS sequence-input qualifiers as on the command line: -sbegin 3 -send -2 -sreverse -slower -supper -scircular -snucleotide -sprotein (begin / end 1-based, negative = from the end; -sreverse reverse-complements a nucleotide sequence, the range then refers to the original). Measured on EMBOSS 6.6.0 per program: " + "distmat ignores -sbegin / -send; -sreverse reverse-complements every row of a nucleotide set; the set is nucleotide when its first row is (or -snucleotide); a protein set uses proteinMethod." + " " + EmbossSequenceQualifiers.CommandLineRules)] string? sequenceQualifiers = null,
        [Description("distmat -protmethod for a protein set when sequenceQualifiers is given: same names as method (default: the same as method).")] string? proteinMethod = null)
    {
        if (sequences is null || sequences.Length < 2)
            throw new ArgumentException("At least two sequences are required", nameof(sequences));
        var (m, name) = (method ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "uncorrected" => (global::Seqeron.Genomics.Analysis.MsaDistanceMethod.Uncorrected, "uncorrected"),
            "jukes-cantor" or "jc" => (global::Seqeron.Genomics.Analysis.MsaDistanceMethod.JukesCantor, "jukes-cantor"),
            "kimura" => (global::Seqeron.Genomics.Analysis.MsaDistanceMethod.Kimura, "kimura"),
            "tamura" => (global::Seqeron.Genomics.Analysis.MsaDistanceMethod.Tamura, "tamura"),
            "tajima-nei" => (global::Seqeron.Genomics.Analysis.MsaDistanceMethod.TajimaNei, "tajima-nei"),
            "jin-nei" => (global::Seqeron.Genomics.Analysis.MsaDistanceMethod.JinNei, "jin-nei"),
            _ => throw new ArgumentException("method must be 'uncorrected', 'jukes-cantor', 'kimura', 'tamura', 'tajima-nei' or 'jin-nei'", nameof(method)),
        };
        var distanceOptions = new global::Seqeron.Genomics.Emboss.PairwiseAlignment.MsaDistanceOptions
        {
            Method = m,
            ResidueType = ParseConsensusResidueType(residueType),
            GapWeight = gapWeight,
            Ambiguous = ambiguous,
            Position = position,
            CalculateA = calculateA,
            ParameterA = parameterA,
        };
        global::Seqeron.Genomics.Analysis.MsaDistanceMethod? protein = null;
        if (proteinMethod is not null)
            protein = proteinMethod.Trim().ToLowerInvariant() switch
            {
                "uncorrected" => global::Seqeron.Genomics.Analysis.MsaDistanceMethod.Uncorrected,
                "jukes-cantor" or "jc" => global::Seqeron.Genomics.Analysis.MsaDistanceMethod.JukesCantor,
                "kimura" => global::Seqeron.Genomics.Analysis.MsaDistanceMethod.Kimura,
                _ => throw new ArgumentException("proteinMethod must be 'uncorrected', 'jukes-cantor' or 'kimura'", nameof(proteinMethod)),
            };
        var d = sequenceQualifiers is null
            ? global::Seqeron.Genomics.Emboss.PairwiseAlignment.EmbossAlignmentStatistics.CalculateDistanceMatrix(sequences, distanceOptions)
            : global::Seqeron.Genomics.Emboss.PairwiseAlignment.EmbossAlignmentStatistics.CalculateDistanceMatrix(sequences, EmbossSequenceQualifiers.ParseInputs("distmat", sequenceQualifiers)[0], distanceOptions, protein);
        return new AlignmentDistanceMatrixResult(name, ToJagged(d));
    }

    [McpServerTool(Name = "multiple_alignment_statistics", Title = "Alignment — Multiple Alignment Statistics", ReadOnly = true)]
    [Description("Statistics of a gapped multiple alignment: per-sequence comparison with a reference row or the EMBOSS cons consensus (EMBOSS infoalign: seqLength, alignLength, gap runs, gap length, identities, similarities, differences, % change; verified against infoalign 6.6.0), per-column residue counts / gap fraction / most common residue / ClustalW conservation symbol ('*' identical, ':' strong, '.' weak group; ClustalW 2.1; DNA columns only '*') / Shannon entropy in bits of the non-gap residues / information content in bits (against a uniform 4-letter background for a nucleotide column, 20 letters otherwise; the per-column Jensen-Shannon divergence and sum of pairs of Capra and Singh are the tool msa_conservation), gap percent, fully conserved and gap-free columns, and the pairwise percent-identity matrix (distmat uncorrected, identity over columns where neither row has a gap; gapWeight 1 = Biopython DistanceCalculator('identity')). Rows are upper-cased, '.'/'~' read as '-'.")]
    public static MultipleAlignmentStatisticsResult MultipleAlignmentStatistics(
        [Description("At least two aligned sequences of equal length (DNA or protein, gaps '-').")] string[] sequences,
        [Description("'auto' (type of the first row as EMBOSS), 'nucleotide' (EDNAFULL) or 'protein' (EBLOSUM62; also enables ':'/'.' conservation).")] string residueType = "auto",
        [Description("infoalign -refseq: 0 = the EMBOSS cons consensus (default), k = the k-th sequence (1-based).")] int referenceSequence = 0,
        [Description("Weight of gap columns in the identity matrix denominator L - g + g*gapWeight (distmat -gapweight; default 0).")] double gapWeight = 0)
    {
        if (sequences is null || sequences.Length < 2)
            throw new ArgumentException("At least two sequences are required", nameof(sequences));
        var type = ParseConsensusResidueType(residueType);
        if (referenceSequence < 0 || referenceSequence > sequences.Length)
            throw new ArgumentOutOfRangeException(nameof(referenceSequence), "referenceSequence must be 0 (consensus) or 1..number of sequences.");

        var stats = global::Seqeron.Genomics.Analysis.MultipleAlignmentStatistics.Calculate(sequences, type);
        var info = global::Seqeron.Genomics.Emboss.PairwiseAlignment.EmbossAlignmentStatistics.CalculateSequenceStatistics(
            sequences, referenceSequence == 0 ? null : referenceSequence - 1, type);
        var identity = gapWeight == 0
            ? stats.PairwiseIdentity
            : global::Seqeron.Genomics.Analysis.MultipleAlignmentStatistics.CalculatePairwiseIdentityMatrix(sequences, gapWeight);
        return new MultipleAlignmentStatisticsResult(
            stats.SequenceCount, stats.AlignmentLength, stats.FullyConservedColumns, stats.GapFreeColumns, stats.GapPercent,
            double.IsNaN(stats.MeanPairwiseIdentity) ? null : stats.MeanPairwiseIdentity,
            global::Seqeron.Genomics.Analysis.MultipleAlignmentStatistics.ClustalConservationLine(sequences, type),
            ToJagged(identity),
            info.Select(i => new MsaSequenceItem(i.Index, i.SequenceLength, i.AlignmentLength, i.GapRuns, i.GapLength,
                i.Identities, i.Similarities, i.Differences, i.PercentChange)).ToArray(),
            stats.Columns.Select(c => new MsaColumnItem(c.Index,
                c.Counts.ToDictionary(kv => kv.Key.ToString(), kv => kv.Value), c.GapCount, c.GapFraction,
                c.MostCommonResidue.ToString(), c.MostCommonCount, c.Conservation.ToString(), c.ShannonEntropy, c.InformationContent())).ToArray());
    }
}
