using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Chromosome;
using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// ASSEMBLY-STATS-001 mutation killers (batch 5): exercises the suspicious-region detector's
/// High-N-content and GC-deviation branches (the canonical and earlier batches only hit the
/// low-complexity branch), the inverted syntenic-block branch, and the repetitive-region span.
/// </summary>
[TestFixture]
public class GenomeAssemblyAnalyzer_MutationKillers5_Tests
{
    private static string RandomFrom(string alphabet, int length, int seed)
    {
        var rng = new Random(seed);
        return string.Concat(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]));
    }

    [Test]
    public void FindSuspiciousRegions_HighNContent_IsFlagged()
    {
        // A 100-bp N run inside otherwise complex sequence: a 500-bp window has 100/500 = 20% N
        // (> 10%) ⇒ "High N content".
        string seq = RandomFrom("ACGT", 250, 1) + new string('N', 100) + RandomFrom("ACGT", 250, 3);
        var regions = FindSuspiciousRegions(new[] { ("s", seq) }).ToList();
        Assert.That(regions.Any(r => r.Reason.Contains("High N content")), Is.True);
    }

    [Test]
    public void FindSuspiciousRegions_GcDeviation_IsFlagged()
    {
        // GC-rich half (≈75% GC) then AT-rich half (≈25% GC): global ≈ 0.5, so 500-bp windows in
        // either half deviate by ≈0.25 (> 0.15) ⇒ "GC deviation". Both halves stay complex (4 bases).
        string seq = RandomFrom("GGGCCCAT", 1000, 1) + RandomFrom("AAATTTGC", 1000, 2);
        var regions = FindSuspiciousRegions(new[] { ("s", seq) }).ToList();
        Assert.That(regions.Any(r => r.Reason.Contains("GC deviation")), Is.True);
    }

    [Test]
    public void FindSyntenicBlocks_ReverseOrderedAnchors_AreInverted()
    {
        // seq2 carries seq1's three 21-mer blocks in reverse order, so anchors map to DECREASING
        // target positions ⇒ the block is flagged inverted.
        string seq1 = RandomFrom("ACGT", 63, 5);
        string seq2 = seq1.Substring(42, 21) + seq1.Substring(21, 21) + seq1.Substring(0, 21);
        var blocks = FindSyntenicBlocks(
            new[] { ("s1", seq1) }, new[] { ("s2", seq2) }, minBlockSize: 25, kmerSize: 21).ToList();

        Assert.That(blocks, Has.Count.EqualTo(1));
        Assert.That(blocks[0].IsInverted, Is.True);
    }

    [Test]
    public void FindRepetitiveRegions_TandemSpansWholeSequence()
    {
        // (CAG)×200 = 600 bp of high-copy 6-mers ⇒ one repetitive region spanning [0, 599].
        string seq = string.Concat(Enumerable.Repeat("CAG", 200));
        var r = FindRepetitiveRegions(new[] { ("s", seq) }, kmerSize: 6, minCopies: 3, windowSize: 60).Single();
        Assert.That(r.SequenceId, Is.EqualTo("s"));
        Assert.That(r.Start, Is.EqualTo(0));
        Assert.That(r.End, Is.EqualTo(593)); // truncated final windows drop below threshold
        Assert.That(r.Copies, Is.GreaterThanOrEqualTo(3));
    }
}
