using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for sequence/genome assembly algorithms (SequenceAssembler).
/// Verifies invariants from the literature each algorithm implements.
///
/// Test Units: ASSEMBLY-CONSENSUS-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Assembly")]
public class AssemblyProperties
{
    /// <summary>Generates 2..5 equal-length DNA reads (gap-free) over {A,C,G,T}.</summary>
    private static Arbitrary<string[]> AlignedReadsArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int rows = 2 + rng.Next(4);   // 2..5 reads
            int len = 3 + rng.Next(6);    // length 3..8
            const string bases = "ACGT";
            var reads = new string[rows];
            for (int r = 0; r < rows; r++)
            {
                var chars = new char[len];
                for (int i = 0; i < len; i++) chars[i] = bases[rng.Next(4)];
                reads[r] = new string(chars);
            }
            return reads;
        }).ToArbitrary();

    #region ASSEMBLY-CONSENSUS-001: R: length = longest read; P: each committed position = majority base; D: deterministic

    // ComputeConsensus follows Biopython dumb_consensus / EMBOSS cons: per column, tally non-gap
    // residues; emit the unique maximum residue when its frequency ≥ threshold, otherwise the
    // ambiguous symbol 'N'. Consensus length equals the longest read.

    /// <summary>
    /// INV-1 (R): The consensus length equals the longest input read (the full alignment width).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Consensus_Length_EqualsLongestRead()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
        {
            int longest = reads.Max(r => r.Length);
            int len = SequenceAssembler.ComputeConsensus(reads).Length;
            return (len == longest).Label($"consensus length {len} ≠ longest read {longest}");
        });
    }

    /// <summary>
    /// INV-2 (P): Every emitted character is either the ambiguous symbol or a residue actually
    /// observed in that column — the consensus never invents a base.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Consensus_EachChar_IsAmbiguousOrObserved()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
        {
            string consensus = SequenceAssembler.ComputeConsensus(reads);
            for (int pos = 0; pos < consensus.Length; pos++)
            {
                if (consensus[pos] == 'N') continue;
                var column = reads.Where(r => pos < r.Length).Select(r => char.ToUpperInvariant(r[pos])).ToHashSet();
                if (!column.Contains(consensus[pos]))
                    return false.Label($"position {pos}: '{consensus[pos]}' not present in column");
            }
            return true.Label("all columns consistent");
        });
    }

    /// <summary>
    /// INV-3 (P, majority rule): when one base occupies a strict majority of a column (count &gt; half
    /// the non-gap residues), the default-threshold (0.5) consensus must emit exactly that base.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Consensus_StrictMajorityColumn_EmitsMajorityBase()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
        {
            string consensus = SequenceAssembler.ComputeConsensus(reads);
            for (int pos = 0; pos < consensus.Length; pos++)
            {
                var column = reads.Where(r => pos < r.Length).Select(r => char.ToUpperInvariant(r[pos])).ToList();
                var grouped = column.GroupBy(c => c).Select(g => (Base: g.Key, Count: g.Count())).ToList();
                var top = grouped.OrderByDescending(g => g.Count).First();
                bool strictMajority = top.Count * 2 > column.Count;
                if (strictMajority && consensus[pos] != top.Base)
                    return false.Label($"position {pos}: majority '{top.Base}' but consensus '{consensus[pos]}'");
            }
            return true.Label("all columns consistent");
        });
    }

    /// <summary>
    /// INV-4 (D): Consensus computation is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Consensus_IsDeterministic()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
            (SequenceAssembler.ComputeConsensus(reads) == SequenceAssembler.ComputeConsensus(reads))
                .Label("ComputeConsensus must be deterministic"));
    }

    /// <summary>
    /// INV-5 (positive/boundary controls): unanimous columns reproduce the read; tied columns and
    /// the empty input behave per the dumb_consensus rule.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Consensus_GoldenAndBoundaryCases()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAssembler.ComputeConsensus(new[] { "ACGT", "ACGT", "ACGT" }), Is.EqualTo("ACGT"),
                "unanimous columns yield the agreed base");
            // Two reads disagree at every column → ties → all ambiguous.
            Assert.That(SequenceAssembler.ComputeConsensus(new[] { "AA", "CC" }), Is.EqualTo("NN"),
                "tied columns emit the ambiguous symbol");
            Assert.That(SequenceAssembler.ComputeConsensus(Array.Empty<string>()), Is.EqualTo(""),
                "empty input yields an empty consensus");
            // Gap symbols are skipped: column 1 has only 'A' (the dash is ignored) → committed.
            Assert.That(SequenceAssembler.ComputeConsensus(new[] { "A-", "AA" }), Is.EqualTo("AA"),
                "gaps are excluded from the column tally");
        });
    }

    #endregion
}
