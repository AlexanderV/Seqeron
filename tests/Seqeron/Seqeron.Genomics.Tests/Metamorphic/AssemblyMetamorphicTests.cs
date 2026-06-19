using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Assembly area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ASSEMBLY-CONSENSUS-001 — column-wise consensus (Assembly).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 140.
///
/// API under test (SequenceAssembler.ComputeConsensus):
///   Per alignment column, non-gap residues are tallied and the strict-majority residue is emitted
///   when it meets the threshold (Biopython dumb_consensus). The result depends only on the
///   per-column residue multiset.
///
/// Relations (derived from the column-tally rule, NOT from output):
///   • INV  (read order independent): each column's tally is a multiset, so permuting the reads
///          leaves the consensus unchanged.
///   • MON  (adding a concordant read preserves the consensus): a read equal to the current
///          consensus only reinforces each column's winning residue, so the consensus is unchanged.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class AssemblyMetamorphicTests
{
    // Reads whose every column has a clear strict majority, so the consensus is fully committed.
    private static readonly string[] AlignedReads =
    {
        "ACGTAACC",
        "ACGTAACC",
        "ACGTAACG",
        "ACGTTACC",
    };

    #region ASSEMBLY-CONSENSUS-001 INV — consensus is independent of read order

    [Test]
    [Description("INV: each column's residue tally is a multiset, so reversing (permuting) the read order yields the identical consensus.")]
    public void Consensus_ReadOrder_Invariant()
    {
        string original = SequenceAssembler.ComputeConsensus(AlignedReads);
        var permuted = AlignedReads.Reverse().ToList();

        SequenceAssembler.ComputeConsensus(permuted).Should().Be(original,
            because: "the per-column majority depends only on the residue multiset, not the read order");
    }

    #endregion

    #region ASSEMBLY-CONSENSUS-001 MON — a concordant read preserves the consensus

    [Test]
    [Description("MON: adding a read equal to the current consensus reinforces each column's winning residue without displacing it, so the consensus is unchanged.")]
    public void Consensus_AddConcordantRead_Unchanged()
    {
        string consensus = SequenceAssembler.ComputeConsensus(AlignedReads);

        var withConcordant = new List<string>(AlignedReads) { consensus };
        SequenceAssembler.ComputeConsensus(withConcordant).Should().Be(consensus,
            because: "a read matching the consensus adds a vote to each column's existing majority residue, preserving it");
    }

    #endregion
}
