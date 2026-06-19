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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ASSEMBLY-CORRECT-001 — k-mer spectrum read error correction (Assembly).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 141.
    //
    // API under test (SequenceAssembler.ErrorCorrectReads):
    //   Musket/Quake two-sided correction: a position covered only by untrusted (low-multiplicity)
    //   k-mers is substituted by the unique base that makes every covering k-mer trusted.
    //
    // Relations (derived from the k-mer-spectrum model, NOT from output):
    //   • INV  (error-free reads unchanged): identical correct reads make every k-mer trusted, so
    //          no position is corrected and the reads are returned unchanged.
    //   • MON  (more coverage ⇒ ≤ residual errors): increasing the number of correct copies raises
    //          the multiplicity of the true k-mers above the trusted cut-off, enabling the unique
    //          correction, so the residual error count of an erroneous read is non-increasing.
    // ───────────────────────────────────────────────────────────────────────────

    private const string TrueSequence = "ACGTTGCAACGTGGATCCGT";
    private const int CorrectKmerSize = 6;
    private const int TrustedCutoff = 2;

    private static string SubstituteAt(string seq, int index)
    {
        char[] arr = seq.ToCharArray();
        arr[index] = arr[index] == 'A' ? 'C' : 'A';
        return new string(arr);
    }

    private static int Hamming(string a, string b) =>
        a.Zip(b, (x, y) => x == y ? 0 : 1).Sum();

    #region ASSEMBLY-CORRECT-001 INV — error-free reads are returned unchanged

    [Test]
    [Description("INV: when all reads are identical and correct every k-mer is trusted, so no position qualifies for correction and the reads are returned unchanged.")]
    public void ErrorCorrect_ErrorFreeReads_Unchanged()
    {
        var reads = new List<string> { TrueSequence, TrueSequence, TrueSequence };

        SequenceAssembler.ErrorCorrectReads(reads, CorrectKmerSize, TrustedCutoff)
            .Should().Equal(reads, because: "every k-mer of an all-correct read set is trusted, so nothing is changed");
    }

    #endregion

    #region ASSEMBLY-CORRECT-001 MON — more coverage cannot increase residual errors

    [Test]
    [Description("MON: adding correct copies raises the true k-mers' multiplicity above the trusted cut-off, enabling the unique correction, so the residual error count of an erroneous read is non-increasing.")]
    public void ErrorCorrect_MoreCoverage_FewerOrEqualResidualErrors()
    {
        string erroneous = SubstituteAt(TrueSequence, index: 10);

        int previous = int.MaxValue;
        int lowestCoverageErrors = -1, highestCoverageErrors = -1;
        int[] coverages = { 1, 2, 3, 5 };

        foreach (int correctCopies in coverages)
        {
            var reads = Enumerable.Repeat(TrueSequence, correctCopies).Append(erroneous).ToList();
            var corrected = SequenceAssembler.ErrorCorrectReads(reads, CorrectKmerSize, TrustedCutoff);

            // The erroneous read is the last entry; residual errors = Hamming distance to the truth.
            int residual = Hamming(corrected[^1], TrueSequence);
            residual.Should().BeLessThanOrEqualTo(previous, because: $"raising correct-copy coverage to {correctCopies} cannot introduce errors");
            previous = residual;

            if (correctCopies == coverages.First()) lowestCoverageErrors = residual;
            if (correctCopies == coverages.Last()) highestCoverageErrors = residual;
        }

        lowestCoverageErrors.Should().BeGreaterThan(0, because: "at coverage 1 the true k-mers are not yet trusted, so the error cannot be corrected");
        highestCoverageErrors.Should().Be(0, because: "at high coverage the true k-mers are trusted and the substitution is uniquely corrected");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ASSEMBLY-COVER-001 — per-base coverage depth (Assembly).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 142.
    //
    // API under test (SequenceAssembler.CalculateCoverage):
    //   Each read is placed independently at its best ungapped match and increments the depth of the
    //   reference positions it spans. Depth at a position is the count of reads covering it.
    //
    // Relations (derived from independent per-read placement, NOT from output):
    //   • INV  (read order independent): each read's contribution is independent of the others, so
    //          permuting the reads yields the identical depth array.
    //   • ADD  (coverage additive over reads): the depth of a combined read set equals the
    //          element-wise sum of the depths of any partition of that set.
    // ───────────────────────────────────────────────────────────────────────────

    private const string CoverageReference = "ACGTTGCAACGTGGATCCGTACGATCGATT";
    private const int CoverageMinOverlap = 6;

    // Exact substrings of the reference, so every read places at its origin.
    private static readonly string[] CoverageReads =
    {
        CoverageReference.Substring(0, 10),
        CoverageReference.Substring(5, 12),
        CoverageReference.Substring(12, 10),
        CoverageReference.Substring(20, 10),
    };

    #region ASSEMBLY-COVER-001 INV — coverage is independent of read order

    [Test]
    [Description("INV: each read is placed independently, so reversing (permuting) the read order yields the identical per-base depth array.")]
    public void Coverage_ReadOrder_Invariant()
    {
        var original = SequenceAssembler.CalculateCoverage(CoverageReference, CoverageReads, CoverageMinOverlap);
        var permuted = SequenceAssembler.CalculateCoverage(CoverageReference, CoverageReads.Reverse().ToList(), CoverageMinOverlap);

        permuted.Should().Equal(original, because: "per-read placement does not depend on the order of the reads");
    }

    #endregion

    #region ASSEMBLY-COVER-001 ADD — coverage is additive over a read partition

    [Test]
    [Description("ADD: each read contributes independently, so the depth of the full read set equals the element-wise sum of the depths of any two-part partition.")]
    public void Coverage_Additive_OverReadPartition()
    {
        var all = SequenceAssembler.CalculateCoverage(CoverageReference, CoverageReads, CoverageMinOverlap);

        var firstHalf = CoverageReads.Take(2).ToList();
        var secondHalf = CoverageReads.Skip(2).ToList();
        var covFirst = SequenceAssembler.CalculateCoverage(CoverageReference, firstHalf, CoverageMinOverlap);
        var covSecond = SequenceAssembler.CalculateCoverage(CoverageReference, secondHalf, CoverageMinOverlap);

        var summed = covFirst.Zip(covSecond, (x, y) => x + y).ToArray();
        all.Should().Equal(summed, because: "depth is a sum of independent per-read contributions, so it is additive over any partition of the reads");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ASSEMBLY-DBG-001 — de Bruijn graph construction (Assembly).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 143.
    //
    // API under test (SequenceAssembler.BuildDeBruijnGraph):
    //   Each input k-mer becomes a directed edge from its (k-1)-mer prefix to its (k-1)-mer suffix;
    //   the graph is the out-adjacency multimap.
    //
    // Relations (derived from the k-mer edge definition, NOT from output):
    //   • INV  (read order independent): the graph is the accumulation of per-k-mer edges, so the
    //          node set and edge multiset do not depend on the order of the reads.
    //   • MON  (larger k ⇒ ≤ spurious joins): a (k-1)-mer shared by unrelated contexts creates a
    //          branching node (a spurious join); increasing k makes shared (k-1)-mers rarer, so the
    //          number of branching nodes is non-increasing.
    // ───────────────────────────────────────────────────────────────────────────

    // Number of branching nodes (≥ 2 distinct successors) — the de Bruijn graph's spurious joins.
    private static int BranchingNodes(IReadOnlyList<string> reads, int k) =>
        SequenceAssembler.BuildDeBruijnGraph(reads, k).Count(kv => kv.Value.Distinct().Count() > 1);

    // Canonical graph: each node's successor multiset sorted, so two graphs compare by content.
    private static Dictionary<string, List<string>> CanonicalGraph(IReadOnlyList<string> reads, int k) =>
        SequenceAssembler.BuildDeBruijnGraph(reads, k)
            .ToDictionary(kv => kv.Key, kv => kv.Value.OrderBy(s => s, System.StringComparer.Ordinal).ToList());

    #region ASSEMBLY-DBG-001 INV — the graph is independent of read order

    [Test]
    [Description("INV: the de Bruijn graph accumulates one edge per k-mer, so reversing (permuting) the read order yields the identical node set and edge multiset.")]
    public void DeBruijn_ReadOrder_Invariant()
    {
        var reads = new[] { "GATCGAAACCC", "AAACCCGATCC", "GATCCTTTAAA" };
        const int k = 4;

        CanonicalGraph(reads.Reverse().ToList(), k).Should().BeEquivalentTo(CanonicalGraph(reads, k),
            because: "edges are accumulated per k-mer, independent of which read (or in which order) supplied them");
    }

    #endregion

    #region ASSEMBLY-DBG-001 MON — larger k cannot increase spurious joins

    [Test]
    [Description("MON: increasing k makes shared (k-1)-mers rarer, so the number of branching nodes (spurious joins) is non-increasing; here a single 4-mer repeat is progressively resolved.")]
    public void DeBruijn_LargerK_FewerOrEqualBranchingNodes()
    {
        // The read repeats the 4-mer GATC in two different contexts (…GATCG… and …GATCC…).
        var reads = new[] { "GATCGAAACCCGATCCTTT" };

        int previous = int.MaxValue;
        int smallestK = -1, largestK = -1;
        int[] kValues = { 4, 5, 6, 7 };

        foreach (int k in kValues)
        {
            int branches = BranchingNodes(reads, k);
            branches.Should().BeLessThanOrEqualTo(previous, because: $"raising k to {k} can only resolve, never create, shared (k-1)-mer junctions");
            previous = branches;

            if (k == kValues.First()) smallestK = branches;
            if (k == kValues.Last()) largestK = branches;
        }

        smallestK.Should().BeGreaterThan(0, because: "at small k the repeated 4-mer collapses into a branching node");
        largestK.Should().Be(0, because: "once k exceeds the repeat length every (k-1)-mer is unique, so no branching remains");
    }

    #endregion
}
