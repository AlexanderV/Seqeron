using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using static Seqeron.Genomics.Annotation.VariantAnnotator;
using VariantType = Seqeron.Genomics.Annotation.VariantAnnotator.VariantType;
using Variant = Seqeron.Genomics.Annotation.VariantAnnotator.Variant;
using Transcript = Seqeron.Genomics.Annotation.VariantAnnotator.Transcript;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// VARIANT-ANNOT-001 mutation killers (batch 4): the Annotate pipeline with codon-level refinement
/// (most-severe selection + RefineAnnotation against an aligned reference) and the AnnotateVariants
/// per-variant grouping.
/// </summary>
[TestFixture]
public class VariantAnnotator_MutationKillers4_Tests
{
    // Single-exon '+' transcript: exon 90..130, CDS 100..120; referenceSequence[0] == genomic 100.
    private const string RefWindow = "ATGCAAGAATTATAANAAGGG";
    private const int RefStart = 100;

    private static Transcript CodingTranscript() => new(
        "ENST_TEST", "ENSG_TEST", "GENE_TEST", "chr1",
        90, 130, '+',
        new List<(int, int)> { (90, 130) },
        new List<(int, int)> { (100, 120) },
        100, 120);

    [Test]
    public void Annotate_RefinesCodingConsequenceAgainstAlignedReference()
    {
        // GAA(Glu, codon3 @106-108) base2 A>T → GTA(Val): RefineAnnotation recovers the exact
        // amino-acid change p.E3V (AnnotateVariant alone, without sequenceStart, cannot).
        var v = new Variant("chr1", 107, "A", "T", VariantType.SNV);
        var a = Annotate(new[] { v }, new[] { CodingTranscript() }, RefWindow, RefStart).Single();

        Assert.That(a.Consequence, Is.EqualTo(ConsequenceType.MissenseVariant));
        Assert.That(a.Impact, Is.EqualTo(ImpactLevel.Moderate));
        Assert.That(a.AminoAcidChange, Is.EqualTo("p.E3V"));
    }

    [Test]
    public void Annotate_PicksMostSevereConsequenceAcrossTranscripts()
    {
        // A premature-stop transcript (HIGH) must outrank a 3'UTR overlap (MODIFIER): the lowest
        // Ensembl rank wins. CAA(Gln, codon2 @103-105) C>T → TAA(Stop).
        var stopVariant = new Variant("chr1", 103, "C", "T", VariantType.SNV);
        var utrTranscript = new Transcript(
            "ENST_UTR", "ENSG_UTR", "GENE_UTR", "chr1",
            90, 130, '+',
            new List<(int, int)> { (90, 130) },
            new List<(int, int)> { (110, 120) }, // CDS starts after 103 ⇒ 103 is 5'UTR here
            110, 120);

        var a = Annotate(new[] { stopVariant }, new[] { CodingTranscript(), utrTranscript }, RefWindow, RefStart).Single();
        Assert.That(a.Consequence, Is.EqualTo(ConsequenceType.StopGained));
        Assert.That(a.TranscriptId, Is.EqualTo("ENST_TEST"));
    }

    [Test]
    public void AnnotateVariants_GroupsAnnotationsByVariant()
    {
        var v = new Variant("chr1", 107, "A", "T", VariantType.SNV);
        var groups = AnnotateVariants(new[] { v }, new[] { CodingTranscript() }).ToList();

        Assert.That(groups, Has.Count.EqualTo(1));
        Assert.That(groups[0].Key, Is.EqualTo(v));
        Assert.That(groups[0].Any(a => a.TranscriptId == "ENST_TEST"), Is.True);
    }

    [Test]
    public void AnnotateVariants_OutOfRangeVariant_GroupsAsIntergenic()
    {
        var v = new Variant("chr1", 100000, "A", "G", VariantType.SNV);
        var group = AnnotateVariants(new[] { v }, new[] { CodingTranscript() }).Single();
        Assert.That(group.Single().Consequence, Is.EqualTo(ConsequenceType.IntergenicVariant));
    }
}
