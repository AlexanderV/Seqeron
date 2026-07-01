using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Chromosome;
using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Chromosome;

/// <summary>
/// ASSEMBLY-STATS-001 mutation killers (batch 7): leading-gap scaffold and gap-free contig edge
/// cases that pin the contig-open / gap-start boundaries of AnalyzeScaffolds and ExtractContigs.
/// </summary>
[TestFixture]
public class GenomeAssemblyAnalyzer_MutationKillers7_Tests
{
    [Test]
    public void AnalyzeScaffolds_LeadingGap_NoEmptyContigAtStart()
    {
        // 10 N then 3 A: the leading gap must NOT open a zero-length contig at index 0
        // (kills 'i >= contigStart' and 'gapStart <= 0' mutants); the only contig is the trailing AAA.
        var s = AnalyzeScaffolds(new[] { ("s", "NNNNNNNNNNAAA") }, minGapLength: 10).Single();
        Assert.That(s.Contigs.Count, Is.EqualTo(1));
        Assert.That(s.Contigs[0], Is.EqualTo(("s_contig1", 10, 12)));
        Assert.That(s.Gaps.Count, Is.EqualTo(1));
        Assert.That(s.Gaps[0].Start, Is.EqualTo(0));
        Assert.That(s.Gaps[0].End, Is.EqualTo(9)); // kills the 'i + 1' gap-end mutant
        Assert.That(s.Gaps[0].Length, Is.EqualTo(10));
    }

    [Test]
    public void ExtractContigs_GapFreeSequence_IsOneContig()
    {
        // No N at all: contigStart stays 0, so the final-contig branch must fire with the inclusive
        // 'contigStart >= 0' guard (a 'contigStart > 0' mutant would drop the whole sequence).
        Assert.That(ExtractContigs(new[] { ("s", "AAAAAA") }, minContigLength: 3).Single(),
            Is.EqualTo(("s_contig1", "AAAAAA")));
    }
}
