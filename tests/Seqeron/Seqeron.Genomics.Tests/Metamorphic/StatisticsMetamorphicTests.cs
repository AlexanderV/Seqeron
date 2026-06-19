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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: SEQ-DINUC-001 — dinucleotide frequencies (Statistics).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 122.
    //
    // API under test (SequenceStatistics.CalculateDinucleotideFrequencies):
    //   f_XY = count(XY) / (N−1) over the sliding dinucleotide windows.
    //
    // Relations (derived from dinucleotide counting, NOT from output):
    //   • INV  (reverse-complement maps each dinucleotide to its revcomp): the count of XY in a
    //          sequence equals the count of revcomp(XY) in its reverse complement, so the
    //          frequency of XY equals the frequency of revcomp(XY) in the reverse complement.
    //   • SHIFT (prepend flank adds only boundary dinucleotides): the dinucleotide multiset of
    //          flank+seq is the multisets of flank and seq plus the single boundary dinucleotide.
    // ───────────────────────────────────────────────────────────────────────────

    private static string RevComp(string dna) =>
        new(dna.Reverse().Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c }).ToArray());

    private static string RevCompDinuc(string xy) => RevComp(xy);

    // Integer dinucleotide counts reconstructed from frequencies (exact for clean A/C/G/T sequences).
    private static System.Collections.Generic.Dictionary<string, int> DinucCounts(string seq)
    {
        var freq = SequenceStatistics.CalculateDinucleotideFrequencies(seq);
        int positions = seq.Length - 1;
        return freq.ToDictionary(kv => kv.Key, kv => (int)System.Math.Round(kv.Value * positions));
    }

    #region SEQ-DINUC-001 INV — reverse complement maps each dinucleotide to its revcomp

    [Test]
    [Description("INV: the count of XY in a sequence equals the count of revcomp(XY) in its reverse complement, so f_seq[XY] = f_revcomp[revcomp(XY)].")]
    public void DinucleotideFrequencies_ReverseComplement_MapsToRevcomp()
    {
        const string seq = "ACGTACGTTGGCCAATAC";
        var freqSeq = SequenceStatistics.CalculateDinucleotideFrequencies(seq);
        var freqRc = SequenceStatistics.CalculateDinucleotideFrequencies(RevComp(seq));

        foreach (var (dinuc, f) in freqSeq)
            freqRc.GetValueOrDefault(RevCompDinuc(dinuc)).Should().BeApproximately(f, 1e-12,
                because: $"dinucleotide {dinuc} maps to {RevCompDinuc(dinuc)} on the reverse-complement strand");
    }

    #endregion

    #region SEQ-DINUC-001 SHIFT — prepending a flank adds only boundary dinucleotides

    [Test]
    [Description("SHIFT: the dinucleotide multiset of flank+seq equals the multisets of the flank and the sequence plus the single boundary dinucleotide (flank-last, seq-first).")]
    public void DinucleotideFrequencies_PrependFlank_AddsOnlyBoundary()
    {
        const string seq = "ACGTACGT";

        foreach (var flank in new[] { "TT", "GCGC", "AATT" })
        {
            var combined = DinucCounts(flank + seq);

            // Expected multiset = flank dinucs ⊎ seq dinucs ⊎ {boundary}.
            var expected = new System.Collections.Generic.Dictionary<string, int>();
            void Add(string d, int n) => expected[d] = expected.GetValueOrDefault(d) + n;
            foreach (var (d, n) in DinucCounts(flank)) Add(d, n);
            foreach (var (d, n) in DinucCounts(seq)) Add(d, n);
            Add($"{flank[^1]}{seq[0]}", 1); // boundary dinucleotide

            combined.Should().BeEquivalentTo(expected,
                because: $"prepending '{flank}' adds the flank's dinucleotides and one boundary dinucleotide, preserving the sequence's own");
        }
    }

    #endregion
}
