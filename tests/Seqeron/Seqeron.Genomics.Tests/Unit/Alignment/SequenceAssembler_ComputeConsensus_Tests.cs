// ASSEMBLY-CONSENSUS-001 — Consensus Computation
// Evidence: docs/Evidence/ASSEMBLY-CONSENSUS-001-Evidence.md
// TestSpec: tests/TestSpecs/ASSEMBLY-CONSENSUS-001.md
// Source: Biopython Bio.Align.AlignInfo.SummaryInfo.dumb_consensus (v1.79) — decision rule
//         (unique max ∧ max_size/num_atoms ≥ threshold else ambiguous; gaps '-'/'.'  skipped;
//         tie → ambiguous; consensus length = alignment length);
//         EMBOSS cons (plurality cut-off); Wikipedia "Consensus sequence" (IUPAC N = any base).

using System;
using System.Collections.Generic;
using NUnit.Framework;

using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

[TestFixture]
public class SequenceAssembler_ComputeConsensus_Tests
{
    #region ComputeConsensus

    // M1 — Every column unanimous → that residue (frequency 1.0 ≥ threshold). (INV-02)
    [Test]
    public void ComputeConsensus_UnanimousColumns_ReturnsThatResidue()
    {
        var reads = new List<string> { "ACGT", "ACGT", "ACGT" };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        Assert.That(consensus, Is.EqualTo("ACGT"),
            "Each column has a single residue at frequency 1.0 ≥ threshold, so it is committed (Biopython dumb_consensus).");
    }

    // M2 — Single majority above threshold → majority residue. Col0 = A,A,A,T (3/4 = 0.75 ≥ 0.7). (INV-02)
    [Test]
    public void ComputeConsensus_MajorityAboveThreshold_ReturnsMajorityResidue()
    {
        var reads = new List<string> { "ACGT", "ACGT", "ACGT", "TCGT" };

        string consensus = SequenceAssembler.ComputeConsensus(reads, threshold: 0.7);

        Assert.That(consensus, Is.EqualTo("ACGT"),
            "Col0 A=3/4=0.75 ≥ 0.7 → A; cols 1-3 unanimous (Biopython threshold rule max_size/num_atoms ≥ threshold).");
    }

    // M3 — Single majority below threshold → ambiguous. Col = A,A,T (2/3 ≈ 0.667 < 0.7). (INV-02)
    [Test]
    public void ComputeConsensus_SubThresholdMajority_ReturnsAmbiguous()
    {
        var reads = new List<string> { "AC", "AC", "TC" };

        string consensus = SequenceAssembler.ComputeConsensus(reads, threshold: 0.7);

        Assert.That(consensus, Is.EqualTo("NC"),
            "Col0 A=2/3≈0.667 < 0.7 → ambiguous N; col1 C unanimous (Biopython: below threshold emits ambiguous).");
    }

    // M4 — Two residues tied for max count → ambiguous, NOT an arbitrary winner. (INV-03)
    [Test]
    public void ComputeConsensus_TiedResidues_ReturnsAmbiguous()
    {
        var reads = new List<string> { "A", "G" };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        Assert.That(consensus, Is.EqualTo("N"),
            "A=1, G=1 tie for max → len(max_atoms)>1 → ambiguous N (Biopython guard len(max_atoms)==1).");
    }

    // M5 — Gap characters '-' and '.' are skipped from the tally and never emitted. (INV-04)
    [Test]
    public void ComputeConsensus_GapCharacters_AreSkipped()
    {
        var reads = new List<string> { "A-GT", "A.GT", "ACGT" };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        Assert.That(consensus, Is.EqualTo("ACGT"),
            "Col1: '-' and '.' skipped, only C counted (1/1 ≥ threshold) → C (Biopython skips '-' and '.').");
    }

    // M6 — Ragged reads: consensus spans the longest read; shorter reads contribute nothing past their end. (INV-01)
    [Test]
    public void ComputeConsensus_RaggedReads_SpansLongestRead()
    {
        var reads = new List<string> { "ACGT", "ACG" };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        Assert.Multiple(() =>
        {
            Assert.That(consensus.Length, Is.EqualTo(4),
                "Consensus length = longest read = 4 (Biopython con_len = get_alignment_length).");
            Assert.That(consensus, Is.EqualTo("ACGT"),
                "Col3: only the length-4 read contributes T (1/1 ≥ threshold); shorter read absent there.");
        });
    }

    // M7 — All-gap column → ambiguous with no division by zero. (INV-05)
    [Test]
    public void ComputeConsensus_AllGapColumn_ReturnsAmbiguous()
    {
        var reads = new List<string> { "A-T", "A-T" };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        Assert.That(consensus, Is.EqualTo("ANT"),
            "Col1 is all gaps (num_atoms=0) → ambiguous N, no division by zero (Biopython short-circuit).");
    }

    // M8 — Empty read list → empty string. (INV-05)
    [Test]
    public void ComputeConsensus_EmptyReadList_ReturnsEmpty()
    {
        string consensus = SequenceAssembler.ComputeConsensus(new List<string>());

        Assert.That(consensus, Is.EqualTo(string.Empty),
            "No reads → no columns → empty consensus (trivial; INV-05).");
    }

    // M9 — Threshold parameter reproduces Biopython's documented 0.7 default behavior.
    [Test]
    public void ComputeConsensus_Threshold070_ReproducesBiopython()
    {
        // Col0 across both: above (4/5=0.8) and below (2/3≈0.667) the 0.7 threshold.
        var above = new List<string> { "A", "A", "A", "A", "T" }; // 4/5 = 0.8 ≥ 0.7
        var below = new List<string> { "A", "A", "T" };           // 2/3 ≈ 0.667 < 0.7

        string aboveConsensus = SequenceAssembler.ComputeConsensus(above, threshold: 0.7);
        string belowConsensus = SequenceAssembler.ComputeConsensus(below, threshold: 0.7);

        Assert.Multiple(() =>
        {
            Assert.That(aboveConsensus, Is.EqualTo("A"),
                "4/5 = 0.80 ≥ 0.7 → A (Biopython threshold formula at default 0.7).");
            Assert.That(belowConsensus, Is.EqualTo("N"),
                "2/3 ≈ 0.667 < 0.7 → ambiguous N (Biopython threshold formula at default 0.7).");
        });
    }

    // S1 — Null input throws ArgumentNullException (repository contract).
    [Test]
    public void ComputeConsensus_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAssembler.ComputeConsensus(null!),
            "Null read list violates the precondition (ArgumentNullException.ThrowIfNull).");
    }

    // S2 — Lowercase residues are normalized to uppercase before tallying.
    [Test]
    public void ComputeConsensus_LowercaseReads_AreNormalized()
    {
        var reads = new List<string> { "acgt", "ACGT" };

        string consensus = SequenceAssembler.ComputeConsensus(reads);

        Assert.That(consensus, Is.EqualTo("ACGT"),
            "char.ToUpperInvariant collapses 'a'/'A' to one residue (2/2 ≥ threshold) → uppercase output.");
    }

    // C1 — Custom ambiguous symbol ('X', Biopython default for protein) is emitted on no-consensus columns.
    [Test]
    public void ComputeConsensus_CustomAmbiguousSymbol_IsEmitted()
    {
        var reads = new List<string> { "A", "G" }; // tie → ambiguous

        string consensus = SequenceAssembler.ComputeConsensus(reads, ambiguous: 'X');

        Assert.That(consensus, Is.EqualTo("X"),
            "Tie column emits the supplied ambiguous symbol 'X' (Biopython default ambiguous='X').");
    }

    #endregion
}
