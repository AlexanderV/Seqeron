using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Analysis area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: GENOMIC-COMMON-001 — common-region / longest-common-substring detection (Analysis).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 175.
///
/// API under test (GenomicAnalyzer.FindLongestCommonRegion / FindCommonRegions):
///   The longest common substring of two sequences, and the set of right-maximal common
///   substrings of length ≥ minLength, via a generalized suffix tree.
///
/// Relations (derived from the common-substring definition, NOT from output):
///   • INV  (input order independent): a common substring occurs in both sequences regardless of
///          which is searched, so the longest-common-region length is symmetric.
///   • SUB  (more inputs ⇒ ⊆ common): a substring common to all of {a,b,c} is common to {a,b}, so
///          requiring an additional sequence can only shrink the common set.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class GenomicAnalysisMetamorphicTests
{
    #region GENOMIC-COMMON-001 INV — longest common region length is symmetric

    [Test]
    [Description("INV: the longest common substring occurs in both sequences regardless of search order, so its length is the same for (a,b) and (b,a).")]
    public void LongestCommonRegion_Symmetric()
    {
        var a = new DnaSequence("ACGTACGTGG");
        var b = new DnaSequence("TTACGTACAA");

        var ab = GenomicAnalyzer.FindLongestCommonRegion(a, b);
        var ba = GenomicAnalyzer.FindLongestCommonRegion(b, a);

        ba.Sequence.Length.Should().Be(ab.Sequence.Length, because: "the longest common substring length does not depend on argument order");
        a.Sequence.Should().Contain(ab.Sequence, because: "a common region occurs in the first sequence");
        b.Sequence.Should().Contain(ab.Sequence, because: "a common region occurs in the second sequence");
    }

    #endregion

    #region GENOMIC-COMMON-001 SUB — requiring an extra input shrinks the common set

    [Test]
    [Description("SUB: a substring common to all of {a,b,c} must be common to {a,b}, so adding the requirement that it also occur in c can only shrink the common-region set.")]
    public void CommonRegions_MoreInputs_Subset()
    {
        const string a = "ACGTACGTGG";
        const string b = "TTACGTACAA";
        const string c = "GGACGTTTTT"; // shares "ACGT" with a/b but not the longer "ACGTAC"

        var commonAb = GenomicAnalyzer.FindCommonRegions(new DnaSequence(a), new DnaSequence(b), 3)
            .Select(r => r.Sequence).ToHashSet();
        var commonAbc = commonAb.Where(sub => c.Contains(sub)).ToHashSet();

        commonAb.Should().NotBeEmpty();
        commonAbc.IsSubsetOf(commonAb).Should().BeTrue(because: "requiring occurrence in c only removes regions");
        commonAbc.Count.Should().BeLessThan(commonAb.Count, because: "c lacks the longer (a,b) common regions such as ACGTAC");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: GENOMIC-MOTIFS-001 — known-motif search (Analysis).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 176.
    //
    // API under test (GenomicAnalyzer.FindKnownMotifs):
    //   For each query motif, returns the sorted 0-based positions of all (overlapping) occurrences.
    //
    // Relations (derived from positional exact matching, NOT from output):
    //   • SHIFT (prepend flank shifts positions): a flank with no occurrence shifts every motif
    //          position by the flank length.
    //   • INV  (deterministic): the result is a pure function with positions sorted ascending.
    // ───────────────────────────────────────────────────────────────────────────

    #region GENOMIC-MOTIFS-001 SHIFT — a prepended flank shifts the motif positions

    [Test]
    [Description("SHIFT: prepending a flank that contains no occurrence (and forms none at the junction) shifts every motif's positions by the flank length.")]
    public void KnownMotifs_PrependFlank_ShiftsPositions()
    {
        const string seq = "ACGTACGTAA";
        var motifs = new[] { "ACGT", "CGT" };
        var original = GenomicAnalyzer.FindKnownMotifs(new DnaSequence(seq), motifs);

        foreach (var flank in new[] { "TT", "GTTG" })
        {
            var shifted = GenomicAnalyzer.FindKnownMotifs(new DnaSequence(flank + seq), motifs);
            foreach (var motif in original.Keys)
                shifted[motif].Should().Equal(original[motif].Select(p => p + flank.Length),
                    because: $"the {flank.Length}-base flank relocates every '{motif}' occurrence by {flank.Length}");
        }
    }

    #endregion

    #region GENOMIC-MOTIFS-001 INV — the search is deterministic and sorted

    [Test]
    [Description("INV: FindKnownMotifs is a pure function returning ascending positions, so repeated calls give the identical result.")]
    public void KnownMotifs_Deterministic_Sorted()
    {
        const string seq = "AAAAACGTAAACGT";
        var motifs = new[] { "AA", "ACGT" };
        var first = GenomicAnalyzer.FindKnownMotifs(new DnaSequence(seq), motifs);
        var second = GenomicAnalyzer.FindKnownMotifs(new DnaSequence(seq), motifs);

        foreach (var motif in first.Keys)
        {
            first[motif].Should().BeInAscendingOrder(because: "occurrence positions are returned sorted");
            second[motif].Should().Equal(first[motif], because: "the search has no hidden state");
        }
    }

    #endregion
}
