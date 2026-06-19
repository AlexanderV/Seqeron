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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: GENOMIC-ORF-001 — open reading frame detection (Analysis).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 177.
    //
    // API under test (GenomicAnalyzer.FindOpenReadingFrames):
    //   ATG→stop ORFs over the three forward frames and the three reverse-complement frames.
    //
    // Relations (derived from the frame scan, NOT from output):
    //   • SHIFT (prepend in-frame flank shifts ORFs): a 5' flank whose length is a multiple of 3
    //          (and contains no start) keeps every forward ORF's sequence and frame, shifting its
    //          position by the flank length.
    //   • INV  (revcomp gives reverse-strand ORFs): the reverse-strand ORFs of a sequence are exactly
    //          the forward-strand ORFs of its reverse complement.
    // ───────────────────────────────────────────────────────────────────────────

    private const int OrfMinLength = 6;

    private static string RevComp(string s) => DnaSequence.GetReverseComplementString(s);

    #region GENOMIC-ORF-001 SHIFT — an in-frame flank shifts forward ORFs

    [Test]
    [Description("SHIFT: a 5' flank whose length is a multiple of 3 (and has no ATG) preserves every forward ORF's sequence and frame and shifts its position by the flank length.")]
    public void Orf_PrependInFrameFlank_ShiftsForwardOrfs()
    {
        const string seq = "CCATGAAATAACC";
        var original = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(seq), OrfMinLength)
            .Where(o => !o.IsReverseComplement).Select(o => (o.Sequence, o.Frame, o.Position)).ToHashSet();
        original.Should().NotBeEmpty();

        foreach (var flank in new[] { "GGG", "GGGGGG" }) // length multiple of 3, no ATG
        {
            var shifted = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(flank + seq), OrfMinLength)
                .Where(o => !o.IsReverseComplement).Select(o => (o.Sequence, o.Frame, o.Position)).ToHashSet();
            shifted.Should().BeEquivalentTo(original.Select(o => (o.Sequence, o.Frame, o.Position + flank.Length)),
                because: $"the in-frame {flank.Length}-base flank keeps each forward ORF's sequence and frame, shifting its position by {flank.Length}");
        }
    }

    #endregion

    #region GENOMIC-ORF-001 INV — reverse-strand ORFs equal forward ORFs of the reverse complement

    [Test]
    [Description("INV: reverse-strand ORFs are computed from the reverse complement, so they equal the forward-strand ORFs of the reverse complement (by sequence).")]
    public void Orf_ReverseStrand_EqualsForwardOrfsOfReverseComplement()
    {
        const string seq = "CCATGAAATAACCTTATTTCAT"; // forward ORF + a reverse-strand ORF

        var reverseOrfs = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(seq), OrfMinLength)
            .Where(o => o.IsReverseComplement).Select(o => o.Sequence).ToHashSet();
        var forwardOrfsOfRevComp = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(RevComp(seq)), OrfMinLength)
            .Where(o => !o.IsReverseComplement).Select(o => o.Sequence).ToHashSet();

        reverseOrfs.Should().NotBeEmpty(because: "the sequence was built to contain a reverse-strand ORF");
        forwardOrfsOfRevComp.Should().BeEquivalentTo(reverseOrfs,
            because: "scanning the reverse complement on its forward strand reproduces the original's reverse-strand ORFs");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: GENOMIC-REPEAT-001 — repeated-substring detection (Analysis).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 178.
    //
    // API under test (GenomicAnalyzer.FindRepeats):
    //   Reports every distinct substring of length ≥ minLength occurring ≥ 2 times, with positions.
    //
    // Relations (derived from the repeat definition, NOT from output):
    //   • MON  (lower minLen ⇒ superset): a lower length cut-off admits the shorter repeats too, so
    //          the repeat set grows.
    //   • SHIFT (prepend flank shifts positions): a flank containing no occurrence of a repeat shifts
    //          that repeat's positions by the flank length.
    // ───────────────────────────────────────────────────────────────────────────

    private const string RepeatSeq = "GATTACACCGATTACA"; // "GATTACA" repeated at 0 and 9

    #region GENOMIC-REPEAT-001 MON — lowering minLength yields a superset

    [Test]
    [Description("MON: a lower length cut-off admits the shorter repeats in addition to the longer ones, so the repeat set at minLength 3 is a superset of the set at minLength 7.")]
    public void Repeats_LowerMinLength_Superset()
    {
        var len7 = GenomicAnalyzer.FindRepeats(new DnaSequence(RepeatSeq), 7).Select(r => r.Sequence).ToHashSet();
        var len3 = GenomicAnalyzer.FindRepeats(new DnaSequence(RepeatSeq), 3).Select(r => r.Sequence).ToHashSet();

        len7.IsSubsetOf(len3).Should().BeTrue(because: "every length-≥7 repeat is also a length-≥3 repeat");
        len3.Count.Should().BeGreaterThan(len7.Count, because: "shorter repeats (e.g. GAT, ATTACA) appear only at the lower cut-off");
    }

    #endregion

    #region GENOMIC-REPEAT-001 SHIFT — a prepended flank shifts repeat positions

    [Test]
    [Description("SHIFT: prepending a flank that contains no occurrence of the repeat shifts that repeat's positions by the flank length.")]
    public void Repeats_PrependFlank_ShiftsPositions()
    {
        var original = GenomicAnalyzer.FindRepeats(new DnaSequence(RepeatSeq), 7).Single(r => r.Sequence == "GATTACA");

        foreach (var flank in new[] { "TT", "CCGCC" }) // contain no GATTACA
        {
            var shifted = GenomicAnalyzer.FindRepeats(new DnaSequence(flank + RepeatSeq), 7).Single(r => r.Sequence == "GATTACA");
            shifted.Positions.Should().Equal(original.Positions.Select(p => p + flank.Length),
                because: $"the {flank.Length}-base flank relocates every GATTACA occurrence by {flank.Length}");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: GENOMIC-SIMILARITY-001 — k-mer Jaccard similarity (Analysis).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 179.
    //
    // API under test (GenomicAnalyzer.CalculateSimilarity):
    //   k-mer Jaccard index |A∩B|/|A∪B|, reported as a percentage in [0, 100].
    //
    // Relations (derived from the Jaccard definition, NOT from output):
    //   • SYM  (sim(a,b)=sim(b,a)): intersection and union are symmetric in the two k-mer sets.
    //   • INV  (sim(x,x)=1): a sequence's k-mer set equals itself, giving Jaccard 1 (100 %).
    // ───────────────────────────────────────────────────────────────────────────

    #region GENOMIC-SIMILARITY-001 SYM — similarity is symmetric

    [Test]
    [Description("SYM: Jaccard intersection and union are symmetric, so CalculateSimilarity(a,b) equals CalculateSimilarity(b,a).")]
    public void Similarity_Symmetric()
    {
        var a = new DnaSequence("ACGTACGTTGGCCAATAC");
        var b = new DnaSequence("ACGTTGCAACGTGGATCC");

        GenomicAnalyzer.CalculateSimilarity(b, a).Should().BeApproximately(GenomicAnalyzer.CalculateSimilarity(a, b), 1e-12,
            because: "the Jaccard k-mer similarity does not depend on argument order");
    }

    #endregion

    #region GENOMIC-SIMILARITY-001 INV — self-similarity is 100 %

    [Test]
    [Description("INV: a sequence's k-mer set equals itself, so its Jaccard self-similarity is 1 (reported as 100 %).")]
    public void Similarity_SelfSimilarity_IsMaximal()
    {
        foreach (var seq in new[] { "ACGTACGTTGGCCAATAC", "AAAAAAAAAA", "GCGCGCGCGCGC" })
            GenomicAnalyzer.CalculateSimilarity(new DnaSequence(seq), new DnaSequence(seq))
                .Should().Be(100.0, because: $"'{seq}' shares all k-mers with itself (Jaccard 1 = 100 %)");
    }

    #endregion
}
