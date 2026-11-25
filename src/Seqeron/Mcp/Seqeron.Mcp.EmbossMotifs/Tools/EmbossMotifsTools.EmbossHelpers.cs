// SPDX-License-Identifier: GPL-2.0-or-later
// This file is part of Seqeron.Mcp.EmbossMotifs, licensed under the terms stated in LICENSE (this directory).
using Seqeron.Genomics.Analysis;

namespace Seqeron.Mcp.EmbossMotifs.Tools;

public partial class EmbossMotifsTools
{
    private static void RequireFuzzproInput(string sequence, string pattern)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (string.IsNullOrWhiteSpace(pattern))
            throw new ArgumentException("Pattern cannot be null or empty", nameof(pattern));
    }

    private static (ProteinPattern Parsed, IReadOnlyList<ProteinPatternHit> Hits) RunFuzzpro(string sequence, string pattern, CancellationToken cancellationToken, int mismatch = 0)
    {
        RequireFuzzproInput(sequence, pattern);
        var parsed = ProteinPattern.Parse(pattern);
        return (parsed, EmbossMotifPrograms.FindProteinPatternHits(sequence, parsed, cancellationToken, mismatch));
    }

    private static TranslatedPatternHitItem TranslatedItem(TranslatedPatternHit h) => new(h.Frame, h.ProteinStart, h.ProteinEnd, h.Translation, h.Position, h.Length, h.EmbossStart, h.EmbossEnd, h.Mismatches);

    private static (ProteinPattern Parsed, IReadOnlyList<TranslatedPatternHit> Hits) RunFuzztran(string sequence, string pattern, string frame, int table, CancellationToken cancellationToken, int mismatch = 0)
    {
        var nucleotides = RequireEmbossIupacDna(sequence, nameof(sequence));
        if (string.IsNullOrWhiteSpace(pattern))
            throw new ArgumentException("Pattern cannot be null or empty", nameof(pattern));
        var parsed = ProteinPattern.Parse(pattern);
        return (parsed, EmbossMotifPrograms.FindTranslatedPatternHits(nucleotides, pattern, frame, table, cancellationToken, mismatch));
    }

    private static global::Seqeron.Genomics.Analysis.ConsensusResidueType ParseConsensusResidueType(string? residueType) =>
        (residueType ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "auto" => global::Seqeron.Genomics.Analysis.ConsensusResidueType.Auto,
            "nucleotide" => global::Seqeron.Genomics.Analysis.ConsensusResidueType.Nucleotide,
            "protein" => global::Seqeron.Genomics.Analysis.ConsensusResidueType.Protein,
            _ => throw new ArgumentException("residueType must be 'nucleotide', 'protein' or 'auto'", nameof(residueType)),
        };
}
