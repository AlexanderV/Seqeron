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
}
