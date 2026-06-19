using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Statistics area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-COMPOSITION-001 — nucleotide composition (Statistics).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 121.
///
/// API under test (SequenceStatistics.CalculateNucleotideComposition):
///   Counts A/T/G/C/U/N and derives GC%, AT%, GC-skew and AT-skew — all functions of the
///   per-base counts only.
///
/// Relations (derived from count-based composition, NOT from output):
///   • INV  (permutation invariant): composition depends only on the multiset of bases, so any
///          reordering of the sequence yields the same composition.
///   • P    (complement swaps counts): the DNA complement maps A↔T and C↔G, so the complemented
///          sequence's A/T and C/G counts are the originals swapped.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class StatisticsMetamorphicTests
{
    private static string Complement(string dna) =>
        new(dna.Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c }).ToArray());

    #region SEQ-COMPOSITION-001 INV — composition is permutation invariant

    [Test]
    [Description("INV: nucleotide composition depends only on the multiset of bases, so reordering the sequence (here reversing it) yields the same counts and derived statistics.")]
    public void NucleotideComposition_Permutation_Invariant()
    {
        const string seq = "AACGTACGTTGGCCAATAC";
        var original = SequenceStatistics.CalculateNucleotideComposition(seq);
        var reversed = SequenceStatistics.CalculateNucleotideComposition(new string(seq.Reverse().ToArray()));

        reversed.Should().Be(original,
            because: "counts, GC%, AT% and the skews are all functions of the per-base counts, which a permutation preserves");
    }

    #endregion

    #region SEQ-COMPOSITION-001 P — the complement swaps A↔T and C↔G counts

    [Test]
    [Description("P: the DNA complement maps A↔T and C↔G, so the complemented sequence's A/T and C/G counts are the originals swapped (and the length is unchanged).")]
    public void NucleotideComposition_Complement_SwapsCounts()
    {
        const string seq = "AACGTACGTTGGCCAATAC";
        var original = SequenceStatistics.CalculateNucleotideComposition(seq);
        var complemented = SequenceStatistics.CalculateNucleotideComposition(Complement(seq));

        complemented.CountA.Should().Be(original.CountT, because: "complement turns every T into A");
        complemented.CountT.Should().Be(original.CountA, because: "complement turns every A into T");
        complemented.CountC.Should().Be(original.CountG, because: "complement turns every G into C");
        complemented.CountG.Should().Be(original.CountC, because: "complement turns every C into G");
        complemented.Length.Should().Be(original.Length, because: "complement is length-preserving");
    }

    #endregion
}
