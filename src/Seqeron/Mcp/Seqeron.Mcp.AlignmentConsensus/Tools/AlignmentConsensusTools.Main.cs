// SPDX-License-Identifier: GPL-3.0-only
// This file is part of Seqeron.Mcp.AlignmentConsensus, licensed under the terms stated in LICENSE (this directory).
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.AlignmentConsensus;
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Analysis.Tools;
using static Seqeron.Mcp.Analysis.Tools.AnalysisTools;

namespace Seqeron.Mcp.AlignmentConsensus.Tools;

public partial class AlignmentConsensusTools
{
    [McpServerTool(Name = "generate_decipher_consensus", Title = "Motifs — DECIPHER Consensus", ReadOnly = true)]
    [Description("Consensus with Bioconductor DECIPHER ConsensusSequence semantics (DNA, RNA or protein): at each position the least frequent characters are dropped while they represent less than 'threshold' of the sequences and the rest are encoded with an IUPAC degeneracy code (DNA/RNA) or B/Z/J/X (protein); noConsensusChar where the consensus carries less than minInformation. Gaps '-'/'.', masks '+'.")]
    public static ConsensusResult GenerateDecipherConsensus(
        [Description("Aligned sequences (unequal lengths allowed, as in DECIPHER).")] string[] sequences,
        [Description("Sequence type: 'dna' (default), 'rna' or 'protein'.")] string sequenceType = "dna",
        [Description("Fraction of sequence information that may be lost at a position, in [0, 1) (default 0.05).")] double threshold = 0.05,
        [Description("Split IUPAC degeneracy codes between their residues (default true).")] bool ambiguity = true,
        [Description("Single character from the alphabet for positions without consensus (default '+').")] string noConsensusChar = "+",
        [Description("Minimum fraction of information in (0, 1]; default 1 - threshold.")] double? minInformation = null,
        [Description("DECIPHER includeNonLetters (default false: gaps/masks are counted; true: they are left out, as DECIPHER's own examples show).")] bool includeNonLetters = false,
        [Description("Count leading/trailing gaps (default false).")] bool includeTerminalGaps = false)
    {
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("At least one sequence is required", nameof(sequences));
        if (noConsensusChar is null || noConsensusChar.Length != 1)
            throw new ArgumentException("noConsensusChar must be exactly one character", nameof(noConsensusChar));
        var type = (sequenceType ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "dna" => global::Seqeron.Genomics.AlignmentConsensus.DecipherSequenceType.Dna,
            "rna" => global::Seqeron.Genomics.AlignmentConsensus.DecipherSequenceType.Rna,
            "protein" or "aa" => global::Seqeron.Genomics.AlignmentConsensus.DecipherSequenceType.AminoAcid,
            _ => throw new ArgumentException("sequenceType must be 'dna', 'rna' or 'protein'", nameof(sequenceType)),
        };

        var consensus = global::Seqeron.Genomics.AlignmentConsensus.DecipherConsensus.GenerateDecipherConsensus(
            sequences, type, threshold, ambiguity, noConsensusChar[0], minInformation,
            includeNonLetters, includeTerminalGaps);
        return new ConsensusResult(consensus);
    }
}
