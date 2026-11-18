using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class MultipleAlignTests
{
    [Test]
    public void MultipleAlign_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.MultipleAlign(new[] { "ACGT", "ACGT" }));
        // Invalid DNA in any sequence is rejected.
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlign(new[] { "ACGT", "ACXT" }));
        // Empty set yields an empty MSA (no throw).
        var empty = AlignmentTools.MultipleAlign(Array.Empty<string>());
        Assert.That(empty.AlignedSequences, Is.Empty);
    }

    [Test]
    public void MultipleAlign_Binding_InvokesSuccessfully()
    {
        // Three identical 8-mers align without gaps; consensus equals the input.
        // Sum-of-pairs score = C(3,2)=3 pairs * 8 matched columns * (+1) = 24.
        var r = AlignmentTools.MultipleAlign(new[] { "ACGTACGT", "ACGTACGT", "ACGTACGT" });
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Has.Length.EqualTo(3));
            Assert.That(r.AlignedSequences, Is.All.EqualTo("ACGTACGT"));
            Assert.That(r.Consensus, Is.EqualTo("ACGTACGT"));
            Assert.That(r.TotalScore, Is.EqualTo(24));
        });
    }

    // B09 F4: multiple_align is the exact center-star: two sequences give an optimal pairwise alignment.
    // Biopython 1.88 optimum (match 1, mismatch -1, linear gap -1) of this pair is 24; the previous
    // suffix-tree-anchored pairwise stage returned 9.
    [Test]
    public void MultipleAlign_TwoSequences_IsOptimalPairwiseAlignment()
    {
        var r = AlignmentTools.MultipleAlign(new[]
        {
            "GGAAGCCAATAATCGTGCACGAAGGTGACTGCAGGAGAGTACGAAGGT",
            "AGACAGCCAAAATCTTCAGAGTTCGACTGGAGGTGACTACGAAGAT"
        });

        Assert.That(r.TotalScore, Is.EqualTo(24));
    }

    [Test]
    public void MultipleAlignVariants_Schema_ValidatesCorrectly()
    {
        var bad = new[] { "ACGT", "ACXT" };
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignAnchored(bad));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignProgressive(bad));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignIterative(bad));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignConsistency(bad));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.MultipleAlignIterative(new[] { "ACGT", "ACGT" }, maxIterations: 0));

        Assert.That(AlignmentTools.MultipleAlignAnchored(Array.Empty<string>()).AlignedSequences, Is.Empty);
        Assert.That(AlignmentTools.MultipleAlignProgressive(Array.Empty<string>()).AlignedSequences, Is.Empty);
        Assert.That(AlignmentTools.MultipleAlignIterative(Array.Empty<string>()).AlignedSequences, Is.Empty);
        Assert.That(AlignmentTools.MultipleAlignConsistency(Array.Empty<string>()).AlignedSequences, Is.Empty);
    }

    // Hand-derived (progressive unit tests): ACGT/ACGT/AGT -> ACGT, ACGT, A-GT; SP = 3 - 1 + 3 + 3 = 8.
    [Test]
    public void MultipleAlignProgressive_Binding_MatchesHandDerivedAlignment()
    {
        var r = AlignmentTools.MultipleAlignProgressive(new[] { "ACGT", "ACGT", "AGT" });

        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "ACGT", "ACGT", "A-GT" }));
            Assert.That(r.TotalScore, Is.EqualTo(8));
        });
    }

    [Test]
    public void MultipleAlignIterative_Binding_NeverBelowProgressive()
    {
        var input = new[] { "ACGTTGCA", "ACGTGCA", "ACTTGCA", "AGTTGGCA" };
        var progressive = AlignmentTools.MultipleAlignProgressive(input);
        var iterative = AlignmentTools.MultipleAlignIterative(input);

        Assert.That(iterative.TotalScore, Is.GreaterThanOrEqualTo(progressive.TotalScore));
    }

    // Independent Python T-Coffee reference (B09 F7): CTT/ACT/CAA -> -CTT, AC-T, -CAA.
    [Test]
    public void MultipleAlignConsistency_Binding_MatchesIndependentReference()
    {
        var r = AlignmentTools.MultipleAlignConsistency(new[] { "CTT", "ACT", "CAA" });

        Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "-CTT", "AC-T", "-CAA" }));
    }

    [Test]
    public void MultipleAlignAnchored_Binding_IdenticalSequencesAreGapFree()
    {
        var r = AlignmentTools.MultipleAlignAnchored(new[] { "ACGTACGT", "ACGTACGT", "ACGTACGT" });

        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Is.All.EqualTo("ACGTACGT"));
            Assert.That(r.TotalScore, Is.EqualTo(24));
        });
    }

    // B09 F51 (A1-15): guide-tree distance source. The set below is discriminating: ClustalW 2.1's pairwise stage and the
    // repository's linear NW alignment give different distances, hence different guide trees and different rows
    // (rows of the independent Python implementation in scripts/alignment/generate_clustalw_pairwise_oracle.py).
    private static readonly string[] DistanceSet = ["CCAACCAACCAA", "ACACTCCGCCCGCACTCAGAACTTT", "TAA"];

    [Test]
    public void MultipleAlignProgressive_DistanceClustalW_UsesClustalWPairwiseStage()
    {
        var clustalw = AlignmentTools.MultipleAlignProgressive(DistanceSet, 5, -4, 0, -3, "clustalw");
        var linear = AlignmentTools.MultipleAlignProgressive(DistanceSet, 5, -4, 0, -3);
        Assert.Multiple(() =>
        {
            Assert.That(clustalw.AlignedSequences, Is.EqualTo(new[] { "-----CCAACC-AAC-C--AA----", "ACACTCCGCCCGCACTCAGAACTTT", "----------------T--AA----" }));
            Assert.That(linear.AlignedSequences, Is.EqualTo(new[] { "-----CCAACC-AAC-C--AA----", "ACACTCCGCCCGCACTCAGAACTTT", "---------------T---AA----" }));
            Assert.That(AlignmentTools.MultipleAlignProgressive(DistanceSet, 5, -4, 0, -3, "LINEAR").AlignedSequences, Is.EqualTo(linear.AlignedSequences));
        });
    }

    [Test]
    public void MultipleAlignIterativeAndConsistency_Distance_ParsedAndValidated()
    {
        Assert.DoesNotThrow(() => AlignmentTools.MultipleAlignIterative(DistanceSet, 5, -4, 0, -3, 16, "clustalw"));
        Assert.DoesNotThrow(() => AlignmentTools.MultipleAlignConsistency(DistanceSet, 5, -4, 0, -3, "clustalw"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignProgressive(DistanceSet, 5, -4, 0, -3, "upgma"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignIterative(DistanceSet, 5, -4, 0, -3, 16, "bogus"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignConsistency(DistanceSet, 5, -4, 0, -3, "bogus"));
    }
}
