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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ASSEMBLY-MERGE-001 — overlap-collapsing contig merge (Assembly).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 144.
    //
    // API under test (SequenceAssembler.MergeContigs):
    //   Merges contig1 and contig2 over a known suffix/prefix overlap, keeping the shared region
    //   once: result = contig1 + contig2[overlap:].
    //
    // Relations (derived from the merge definition, NOT from output):
    //   • INV  (merge order independent for compatible contigs): collapsing overlaps is associative
    //          along a compatible chain A→B→C, so folding left (merge(merge(A,B),C)) equals folding
    //          right (merge(A,merge(B,C))) — both reconstruct the original superstring.
    // ───────────────────────────────────────────────────────────────────────────

    #region ASSEMBLY-MERGE-001 INV — overlap-collapsing merge is associative

    [Test]
    [Description("INV: for a compatible overlap chain A→B→C, collapsing overlaps is associative, so left-folding and right-folding the merges give the same superstring (here the original sequence).")]
    public void MergeContigs_CompatibleChain_OrderIndependent()
    {
        const string genome = "ACGTACGTTTGGCCAATT";
        string a = genome.Substring(0, 10);  // ACGTACGTTT
        string b = genome.Substring(6, 10);  // GTTTGGCCAA  (overlap with A = 4: "GTTT")
        string c = genome.Substring(12, 6);  // CCAATT      (overlap with B = 4: "CCAA")
        const int overlapAb = 4;
        const int overlapBc = 4;

        string leftFold = SequenceAssembler.MergeContigs(SequenceAssembler.MergeContigs(a, b, overlapAb), c, overlapBc);
        string rightFold = SequenceAssembler.MergeContigs(a, SequenceAssembler.MergeContigs(b, c, overlapBc), overlapAb);

        leftFold.Should().Be(rightFold, because: "collapsing compatible overlaps is associative, so the merge order does not matter");
        leftFold.Should().Be(genome, because: "the compatible chain reconstructs the original superstring");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ASSEMBLY-OLC-001 — Overlap-Layout-Consensus assembly (Assembly).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 145.
    //
    // API under test (SequenceAssembler.AssembleOLC):
    //   Finds suffix/prefix overlaps, greedily chains reads by best overlap into contigs.
    //
    // Relations (derived from the overlap-chaining model, NOT from output):
    //   • INV  (read order independent): for an unambiguous tiling the best-overlap chain is
    //          determined by the read content, so reordering the reads yields the same contig set.
    //   • MON  (higher minOverlap ⇒ ≤ joins): raising the minimum overlap can only remove overlap
    //          edges, so the number of joins (reads − contigs) is non-increasing.
    // ───────────────────────────────────────────────────────────────────────────

    private const string OlcGenome = "ACGTTGCAACCGGATTCAGTCCGATACGATGCATTGAC";

    #region ASSEMBLY-OLC-001 INV — assembly is independent of read order

    [Test]
    [Description("INV: for an unambiguous tiling the greedy best-overlap chain depends only on the read content, so reordering the reads produces the same set of contigs.")]
    public void Olc_ReadOrder_Invariant()
    {
        // Three reads tiling the genome with 4-base consecutive overlaps.
        string r0 = OlcGenome.Substring(0, 14);
        string r1 = OlcGenome.Substring(10, 14);
        string r2 = OlcGenome.Substring(20, 18);
        var param = new SequenceAssembler.AssemblyParameters(MinOverlap: 4, MinContigLength: 1);

        var ordered = SequenceAssembler.AssembleOLC(new[] { r0, r1, r2 }, param).Contigs.OrderBy(c => c).ToList();
        var shuffled = SequenceAssembler.AssembleOLC(new[] { r2, r0, r1 }, param).Contigs.OrderBy(c => c).ToList();

        shuffled.Should().Equal(ordered, because: "the unambiguous best-overlap tiling does not depend on the input order of the reads");
    }

    #endregion

    #region ASSEMBLY-OLC-001 MON — higher minOverlap cannot increase the number of joins

    [Test]
    [Description("MON: raising the minimum overlap only removes overlap edges, so the number of joins (reads − contigs) is non-increasing.")]
    public void Olc_HigherMinOverlap_FewerOrEqualJoins()
    {
        // Four reads with consecutive overlaps 8, 6, 4.
        var reads = new[]
        {
            OlcGenome.Substring(0, 14),   // r0
            OlcGenome.Substring(6, 14),   // r1: overlap 8 with r0
            OlcGenome.Substring(14, 14),  // r2: overlap 6 with r1
            OlcGenome.Substring(24, 14),  // r3: overlap 4 with r2
        };

        int previous = int.MaxValue;
        int loosest = -1, strictest = -1;
        int[] minOverlaps = { 4, 5, 7, 9 };

        foreach (int minOverlap in minOverlaps)
        {
            var param = new SequenceAssembler.AssemblyParameters(MinOverlap: minOverlap, MinContigLength: 1);
            var result = SequenceAssembler.AssembleOLC(reads, param);
            int joins = result.TotalReads - result.Contigs.Count;

            joins.Should().BeLessThanOrEqualTo(previous, because: $"raising minOverlap to {minOverlap} removes overlap edges, so joins cannot increase");
            previous = joins;

            if (minOverlap == minOverlaps.First()) loosest = joins;
            if (minOverlap == minOverlaps.Last()) strictest = joins;
        }

        loosest.Should().BeGreaterThan(strictest, because: "the loosest overlap admits joins that the strictest (above every overlap length) rejects");
        strictest.Should().Be(0, because: "no overlap reaches 9, so every read stays its own contig");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ASSEMBLY-SCAFFOLD-001 — paired-end scaffolding (Assembly).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 146.
    //
    // API under test (SequenceAssembler.Scaffold):
    //   Orders contigs along paired-end link paths, inserting gap-character runs sized by the link's
    //   gap estimate (Jackman et al. 2017).
    //
    // Relations (derived from the link-path construction, NOT from output):
    //   • INV  (link order independent): when each contig has at most one forward link (an
    //          unambiguous link set) the scaffold path is fixed, so permuting the link list yields
    //          the identical scaffolds.
    // ───────────────────────────────────────────────────────────────────────────

    #region ASSEMBLY-SCAFFOLD-001 INV — scaffolds are independent of link order

    [Test]
    [Description("INV: with at most one forward link per contig the scaffold path is unambiguous, so permuting the link list produces the identical scaffolds.")]
    public void Scaffold_LinkOrder_Invariant()
    {
        var contigs = new[] { "AAAA", "CCCC", "GGGG", "TTTT" };
        var links = new[] { (0, 1, 3), (1, 2, 5), (2, 3, 2) };
        var shuffled = new[] { (2, 3, 2), (0, 1, 3), (1, 2, 5) };

        var fromOrdered = SequenceAssembler.Scaffold(contigs, links);
        var fromShuffled = SequenceAssembler.Scaffold(contigs, shuffled);

        fromShuffled.Should().Equal(fromOrdered,
            because: "an unambiguous link set determines a fixed scaffold path regardless of the order links are listed in");
    }

    #endregion
}
