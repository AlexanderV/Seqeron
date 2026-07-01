// SV-DETECT-001 — Structural Variant Detection (Paired-End Mapping signatures)
// Evidence: docs/Evidence/SV-DETECT-001-Evidence.md
// TestSpec: tests/TestSpecs/SV-DETECT-001.md
// Source: Medvedev P, Stanciu M, Brudno M (2009). Nat Methods 6(11s):S13-S20, doi:10.1038/nmeth.1374;
//         Chen K et al. (2009) BreakDancer, Nat Methods 6:677-681, doi:10.1038/nmeth.1363 (README -c/-r defaults);
//         FR proper-pair orientation (SAM FLAG 0x02 / BWA), cureffi.org.

using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.StructuralVariantAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

[TestFixture]
public class StructuralVariantAnalyzer_DetectSVs_Tests
{
    // Library parameters used throughout (Evidence dataset): mean 400, sd 50, cutoff c=3 => bounds [250, 550].
    private const int Mean = 400;
    private const int Sd = 50;
    private const double Cutoff = 3.0;

    private static ReadPairSignature Pair(
        string chr1, char strand1, string chr2, char strand2, int insertSize) =>
        new(
            ReadId: "r",
            Chromosome1: chr1,
            Position1: 1000,
            Strand1: strand1,
            Chromosome2: chr2,
            Position2: 1000 + insertSize,
            Strand2: strand2,
            InsertSize: insertSize,
            IsDiscordant: true);

    #region ClassifySV

    // M1 — Deletion: same chr, FR, span 5000 > mean + 3*sd (550). INV-03. Span larger than insert size (Medvedev 2009).
    [Test]
    public void ClassifySV_LargeSpanSameChr_ReturnsDeletion()
    {
        var pair = Pair("chr1", '+', "chr1", '-', insertSize: 5000);

        var type = ClassifySV(pair, Mean, Sd, Cutoff);

        Assert.That(type, Is.EqualTo(SVType.Deletion),
            "A mate pair spanning a deletion maps with distance greater than the insert size (Medvedev et al. 2009; BreakDancer DEL).");
    }

    // M2 — Insertion: same chr, FR, span 100 < mean - 3*sd (250). INV-04. Span smaller than insert size (Medvedev 2009).
    [Test]
    public void ClassifySV_SmallSpanSameChr_ReturnsInsertion()
    {
        var pair = Pair("chr1", '+', "chr1", '-', insertSize: 100);

        var type = ClassifySV(pair, Mean, Sd, Cutoff);

        Assert.That(type, Is.EqualTo(SVType.Insertion),
            "If the event is an insertion the mapped distance is smaller than the insert size (Medvedev et al. 2009; BreakDancer INS).");
    }

    // M3 — Inversion: same chr, same strand (FF), span within bounds. INV-02. Flipped orientation (Medvedev 2009; cureffi/BWA).
    [Test]
    public void ClassifySV_SameOrientationSameChr_ReturnsInversion()
    {
        var pair = Pair("chr1", '+', "chr1", '+', insertSize: Mean);

        var type = ClassifySV(pair, Mean, Sd, Cutoff);

        Assert.That(type, Is.EqualTo(SVType.Inversion),
            "An intra-chromosomal pair whose mates map to the same strand (FF/RR) is the basic inversion signature (Medvedev et al. 2009).");
    }

    // M4 — Translocation: mates on different chromosomes. INV-01. Linking signature across chromosomes (Medvedev 2009; CTX).
    [Test]
    public void ClassifySV_DifferentChromosomes_ReturnsTranslocation()
    {
        var pair = Pair("chr1", '+', "chr2", '-', insertSize: Mean);

        var type = ClassifySV(pair, Mean, Sd, Cutoff);

        Assert.That(type, Is.EqualTo(SVType.Translocation),
            "Mates mapping to different chromosomes form a linking/translocation signature (Medvedev et al. 2009; BreakDancer CTX).");
    }

    // M5 — Precedence: cross-chromosome AND same orientation must classify as Translocation, not Inversion (ASSUMPTION A1 / ASM-02).
    [Test]
    public void ClassifySV_DifferentChromosomesSameOrientation_ReturnsTranslocation()
    {
        var pair = Pair("chr1", '+', "chr2", '+', insertSize: Mean);

        var type = ClassifySV(pair, Mean, Sd, Cutoff);

        Assert.That(type, Is.EqualTo(SVType.Translocation),
            "Chromosome difference takes precedence over orientation: a cross-chromosome event is a translocation, not an inversion (ASM-02).");
    }

    // M9 — Duplication: same chr, RF (reverse-forward, outward-facing / everted) orientation.
    // Tandem-duplication signature per DELLY/LUMPY/Manta/SVXplorer (RF cluster ⇒ duplication).
    [Test]
    public void ClassifySV_RfEvertedOrientationSameChr_ReturnsDuplication()
    {
        var pair = Pair("chr1", '-', "chr1", '+', insertSize: Mean);

        var type = ClassifySV(pair, Mean, Sd, Cutoff);

        Assert.That(type, Is.EqualTo(SVType.Duplication),
            "An intra-chromosomal RF (everted, outward-facing) pair is the tandem-duplication signature: LUMPY/Manta/DELLY/SVXplorer all read an RF cluster as a duplication candidate, distinct from the FR deletion signature.");
    }

    #endregion

    #region FindDiscordantPairs (cutoff and orientation boundaries)

    private static (string, string, int, char, string, int, char, int) Tuple(
        string chr1, char strand1, string chr2, char strand2, int insertSize) =>
        ("r", chr1, 1000, strand1, chr2, 1000 + insertSize, strand2, insertSize);

    // M6 — Concordant FR pair (same chr, span 400 within [250,550]) is not flagged discordant. BreakDancer normal class; FR proper pair.
    [Test]
    public void FindDiscordantPairs_ConcordantFrPair_NotReturned()
    {
        var pairs = new[] { Tuple("chr1", '+', "chr1", '-', insertSize: Mean) };

        var result = FindDiscordantPairs(pairs, Mean, Sd, Cutoff).ToList();

        Assert.That(result, Is.Empty,
            "A same-chromosome FR pair with span inside mean +/- 3*sd is concordant and must not be flagged discordant (BreakDancer normal class).");
    }

    // S1 — Lower boundary inclusive: span exactly mean - 3*sd (250) is concordant (discordant iff strictly outside). INV-05.
    [Test]
    public void FindDiscordantPairs_SpanAtLowerBound_NotDiscordant()
    {
        var pairs = new[] { Tuple("chr1", '+', "chr1", '-', insertSize: Mean - (int)(Cutoff * Sd)) }; // 250

        var result = FindDiscordantPairs(pairs, Mean, Sd, Cutoff).ToList();

        Assert.That(result, Is.Empty,
            "Span exactly at the lower bound mean - c*sd (250) is inside the concordant range; a pair is discordant only when strictly outside (bounds = mean +/- c*std, BreakDancer).");
    }

    // S2 — Just below the lower bound (249) is discordant and classifies as Insertion. INV-04/INV-05.
    [Test]
    public void FindDiscordantPairs_SpanBelowLowerBound_IsDiscordant()
    {
        int span = Mean - (int)(Cutoff * Sd) - 1; // 249
        var pairs = new[] { Tuple("chr1", '+', "chr1", '-', insertSize: span) };

        var result = FindDiscordantPairs(pairs, Mean, Sd, Cutoff).ToList();

        Assert.That(result, Has.Count.EqualTo(1),
            "Span one unit below the lower bound (249 < 250) is outside mean +/- c*sd and must be flagged discordant (BreakDancer).");
        Assert.That(ClassifySV(result[0], Mean, Sd, Cutoff), Is.EqualTo(SVType.Insertion),
            "A smaller-than-expected span on the same chromosome is the insertion signature (Medvedev et al. 2009).");
    }

    // S3 — Upper boundary inclusive: span exactly mean + 3*sd (550) is concordant. INV-05.
    [Test]
    public void FindDiscordantPairs_SpanAtUpperBound_NotDiscordant()
    {
        var pairs = new[] { Tuple("chr1", '+', "chr1", '-', insertSize: Mean + (int)(Cutoff * Sd)) }; // 550

        var result = FindDiscordantPairs(pairs, Mean, Sd, Cutoff).ToList();

        Assert.That(result, Is.Empty,
            "Span exactly at the upper bound mean + c*sd (550) is inside the concordant range (inclusive bounds = mean +/- c*std, BreakDancer).");
    }

    // S4 — RF (outward-facing / everted) orientation is DISCORDANT and is the tandem-duplication
    // signature, even with an in-bounds span. Only FR (inward) is the proper short-insert orientation;
    // DELLY/LUMPY/Manta/SVXplorer all read an RF cluster as a duplication candidate. cureffi/BWA:
    // "RF, FF or RR … that's a problem." SAM proper-pair FLAG 0x02 is set only for FR.
    [Test]
    public void FindDiscordantPairs_RfOrientationWithinBounds_IsDiscordantDuplication()
    {
        var pairs = new[] { Tuple("chr1", '-', "chr1", '+', insertSize: Mean) };

        var result = FindDiscordantPairs(pairs, Mean, Sd, Cutoff).ToList();

        Assert.That(result, Has.Count.EqualTo(1),
            "An outward-facing RF pair is everted relative to the FR proper-pair orientation and must be flagged discordant (DELLY/LUMPY/Manta/SVXplorer treat an RF cluster as a duplication; cureffi/BWA: RF 'is a problem').");
        Assert.That(ClassifySV(result[0], Mean, Sd, Cutoff), Is.EqualTo(SVType.Duplication),
            "An RF (reverse-forward, everted) intra-chromosomal pair is the tandem-duplication signature (DELLY: 'paired-ends where the first and second read changed their relative order but kept the alignment strand'; SVXplorer: 'an RF cluster … as a tandem duplication').");
    }

    #endregion

    #region DetectSVs

    // M7 — Min-support gate (below): a single deletion-signature pair with minSupport=2 yields no SV. INV-06; BreakDancer -r default 2.
    [Test]
    public void DetectSVs_BelowMinSupport_EmitsNoSv()
    {
        var pairs = new[] { Tuple("chr1", '+', "chr1", '-', insertSize: 5000) };

        var svs = DetectSVs(pairs, Mean, Sd, Cutoff, clusterDistance: 500, minSupport: 2).ToList();

        Assert.That(svs, Is.Empty,
            "A cluster with fewer than the minimum supporting read pairs (default 2) is not reported (BreakDancer -r).");
    }

    // M8 — Min-support gate (meets): three clustered deletion-signature pairs yield one Deletion with SupportingReads=3. INV-06.
    [Test]
    public void DetectSVs_MeetsMinSupport_EmitsOneDeletion()
    {
        var pairs = new[]
        {
            Tuple("chr1", '+', "chr1", '-', insertSize: 5000),
            Tuple("chr1", '+', "chr1", '-', insertSize: 5000),
            Tuple("chr1", '+', "chr1", '-', insertSize: 5000),
        };

        var svs = DetectSVs(pairs, Mean, Sd, Cutoff, clusterDistance: 500, minSupport: 2).ToList();

        Assert.That(svs, Has.Count.EqualTo(1),
            "Three nearby deletion-signature pairs cluster into one SV that meets the minimum support of 2 (BreakDancer -r).");
        Assert.Multiple(() =>
        {
            Assert.That(svs[0].Type, Is.EqualTo(SVType.Deletion),
                "A larger-than-expected span on the same chromosome is the deletion signature (Medvedev et al. 2009).");
            Assert.That(svs[0].SupportingReads, Is.EqualTo(3),
                "All three clustered pairs support the call, so SupportingReads must equal 3.");
        });
    }

    // C1 — Empty input yields empty output (defined trivial behavior).
    [Test]
    public void DetectSVs_EmptyInput_ReturnsEmpty()
    {
        var pairs = Array.Empty<(string, string, int, char, string, int, char, int)>();

        var svs = DetectSVs(pairs, Mean, Sd, Cutoff).ToList();

        Assert.That(svs, Is.Empty, "No read pairs means no structural variants can be detected.");
    }

    // C2 — Null input throws ArgumentNullException (input validation contract).
    [Test]
    public void DetectSVs_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => DetectSVs(null!, Mean, Sd, Cutoff).ToList(),
            "Null read-pair input violates the precondition and must throw ArgumentNullException.");
    }

    #endregion
}
