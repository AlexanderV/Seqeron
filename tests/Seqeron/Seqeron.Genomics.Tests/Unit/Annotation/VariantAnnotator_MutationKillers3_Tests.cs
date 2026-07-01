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
/// VARIANT-ANNOT-001 mutation killers (batch 3): drives every DetermineConsequence branch through
/// AnnotateVariant on a two-exon '+' transcript — intron, 5'/3' UTR, splice acceptor/donor/region,
/// start-loss, frameshift and in-frame indels (Ensembl VEP consequence predicates).
/// </summary>
[TestFixture]
public class VariantAnnotator_MutationKillers3_Tests
{
    // exon1 100..150, intron 151..199, exon2 200..250; CDS 120..240 ('+').
    // ⇒ 5'UTR 100..119, coding 120..240, 3'UTR 241..250.
    private static Transcript TwoExon() => new(
        "ENST2", "ENSG2", "GENE2", "chr1",
        100, 250, '+',
        new List<(int, int)> { (100, 150), (200, 250) },
        new List<(int, int)> { (120, 150), (200, 240) },
        120, 240);

    // Non-null reference avoids any null path in the downstream amino-acid-change step; the
    // consequence itself is determined independently of the reference by DetermineConsequence.
    private static readonly string RefSeq = new string('A', 300);

    private static ConsequenceType Consequence(int pos, string r, string a, VariantType t) =>
        AnnotateVariant(new Variant("chr1", pos, r, a, t), new[] { TwoExon() }, RefSeq).Single().Consequence;

    [Test]
    public void IntronVariant()
        => Assert.That(Consequence(175, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.IntronVariant));

    [Test]
    public void FivePrimeUtr()
        => Assert.That(Consequence(110, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.FivePrimeUtrVariant));

    [Test]
    public void ThreePrimeUtr()
        => Assert.That(Consequence(245, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.ThreePrimeUtrVariant));

    [Test]
    public void SpliceDonor()
        // 2 bp after exon1 end (151..152).
        => Assert.That(Consequence(151, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.SpliceDonorVariant));

    [Test]
    public void SpliceAcceptor()
        // 2 bp before exon2 start (198..199).
        => Assert.That(Consequence(199, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.SpliceAcceptorVariant));

    [Test]
    public void SpliceRegion()
        // 3..8 bp into the intron after exon1 (153..158).
        => Assert.That(Consequence(155, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.SpliceRegionVariant));

    [Test]
    public void StartLost()
        // SNV at the CDS start codon (120..122).
        => Assert.That(Consequence(120, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.StartLost));

    [Test]
    public void Frameshift()
        // 1-bp coding deletion (length change not a multiple of 3).
        => Assert.That(Consequence(130, "AC", "A", VariantType.Deletion), Is.EqualTo(ConsequenceType.FrameshiftVariant));

    [Test]
    public void InframeDeletion()
        // 3-bp coding deletion (length change ≡ 0 mod 3, net loss).
        => Assert.That(Consequence(130, "ACGT", "A", VariantType.Deletion), Is.EqualTo(ConsequenceType.InframeDeletion));

    [Test]
    public void InframeInsertion()
        // 3-bp coding insertion (length change ≡ 0 mod 3, net gain).
        => Assert.That(Consequence(130, "A", "ATGC", VariantType.Insertion), Is.EqualTo(ConsequenceType.InframeInsertion));

    [Test]
    public void CodingSnv_DefaultsToMissense()
        => Assert.That(Consequence(135, "A", "G", VariantType.SNV), Is.EqualTo(ConsequenceType.MissenseVariant));

    [Test]
    public void NonCodingExon_WhenNoCds()
    {
        var t = new Transcript("ENST3", "ENSG3", "GENE3", "chr1", 100, 250, '+',
            new List<(int, int)> { (100, 150), (200, 250) },
            new List<(int, int)>(), null, null);
        var c = AnnotateVariant(new Variant("chr1", 130, "A", "G", VariantType.SNV), new[] { t }).Single().Consequence;
        Assert.That(c, Is.EqualTo(ConsequenceType.NonCodingTranscriptExonVariant));
    }
}
