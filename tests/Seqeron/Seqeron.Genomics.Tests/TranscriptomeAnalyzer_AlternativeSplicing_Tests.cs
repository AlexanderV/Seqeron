// TRANS-SPLICE-001 — Alternative Splicing (Event Classification + Percent Spliced In)
// Evidence: docs/Evidence/TRANS-SPLICE-001-Evidence.md
// TestSpec: tests/TestSpecs/TRANS-SPLICE-001.md
// Source: Wang ET et al. (2008). Nature 456(7221):470-476 (five AS event classes);
//         BMC Bioinformatics 13(Suppl 6):S11, PMC3330053 (PSI = I/(I+S));
//         Shen S et al. (2014). PNAS 111(51):E5593 — rMATS (length-normalized PSI);
//         Trincado JL et al. (2018). SUPPA2 (PSI read-count definition).

using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class TranscriptomeAnalyzer_AlternativeSplicing_Tests
{
    private const double Tol = 1e-10;

    private static TranscriptomeAnalyzer.TranscriptIsoform Iso(
        string transcriptId, string geneId, params (int Start, int End)[] exons)
        => new(transcriptId, geneId, Length: exons.Sum(e => e.End - e.Start + 1),
               ExonCount: exons.Length, Expression: 1.0, IsProteinCoding: true,
               Exons: exons.ToList());

    #region CalculatePSI

    // M1 — Ψ = I/(I+S) = 80/100 = 0.80. Source: PMC3330053; SUPPA2.
    [Test]
    public void CalculatePSI_InclusionEightyExclusionTwenty_ReturnsZeroPointEight()
    {
        double psi = TranscriptomeAnalyzer.CalculatePSI(80, 20);

        Assert.That(psi, Is.EqualTo(0.80).Within(Tol),
            "Ψ = I/(I+S) = 80/(80+20) = 0.80 per the PSI read-count definition (PMC3330053, SUPPA2).");
    }

    // M2 — rMATS ψ̂ = (I/l_I)/(I/l_I + S/l_S) = (80/200)/((80/200)+(20/100)) = 0.4/0.6 = 2/3.
    // Source: Shen et al. (2014) rMATS, PNAS.
    [Test]
    public void CalculatePSI_WithEffectiveLengths_UsesRMatsNormalizedValue()
    {
        double psi = TranscriptomeAnalyzer.CalculatePSI(80, 20, inclusionEffectiveLength: 200, exclusionEffectiveLength: 100);

        Assert.That(psi, Is.EqualTo(2.0 / 3.0).Within(Tol),
            "rMATS ψ̂ = (80/200)/((80/200)+(20/100)) = 0.4/0.6 = 0.6666… (Shen et al. 2014).");
    }

    // M3 — INV-03: S=0,I>0 ⇒ Ψ=1. Source: PMC3330053.
    [Test]
    public void CalculatePSI_NoExclusionReads_ReturnsOne()
    {
        double psi = TranscriptomeAnalyzer.CalculatePSI(50, 0);

        Assert.That(psi, Is.EqualTo(1.0).Within(Tol),
            "INV-03: with only inclusion reads Ψ = I/(I+0) = 1 (fully included).");
    }

    // M4 — INV-03: I=0,S>0 ⇒ Ψ=0. Source: PMC3330053.
    [Test]
    public void CalculatePSI_NoInclusionReads_ReturnsZero()
    {
        double psi = TranscriptomeAnalyzer.CalculatePSI(0, 40);

        Assert.That(psi, Is.EqualTo(0.0).Within(Tol),
            "INV-03: with only skipping reads Ψ = 0/(0+S) = 0 (fully excluded).");
    }

    // M5 — INV-02: I+S=0 ⇒ Ψ undefined (NaN). Source: PMC3330053 (0/0 undefined).
    [Test]
    public void CalculatePSI_NoReads_ReturnsNaN()
    {
        double psi = TranscriptomeAnalyzer.CalculatePSI(0, 0);

        Assert.That(double.IsNaN(psi), Is.True,
            "INV-02: 0/0 is undefined; PSI returns NaN when there are no supporting reads (PMC3330053).");
    }

    // S1 — INV-01: 0 ≤ Ψ ≤ 1 for any non-negative I,S with I+S>0. Source: part/whole ratio (PMC3330053).
    [Test]
    public void CalculatePSI_VariousNonNegativeCounts_StaysInUnitInterval()
    {
        var cases = new (double I, double S)[] { (1, 9), (7, 7), (100, 1), (3, 0), (0, 5) };

        Assert.Multiple(() =>
        {
            foreach (var (i, s) in cases)
            {
                double psi = TranscriptomeAnalyzer.CalculatePSI(i, s);
                Assert.That(psi, Is.InRange(0.0, 1.0),
                    $"INV-01: Ψ for (I={i},S={s}) must lie in [0,1] as a part/whole ratio.");
            }
        });
    }

    // S5 — negative read counts are invalid. Source: counts are non-negative (contract).
    [Test]
    public void CalculatePSI_NegativeReadCounts_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => TranscriptomeAnalyzer.CalculatePSI(-1, 5),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(),
                "Negative inclusion read count is invalid.");
            Assert.That(() => TranscriptomeAnalyzer.CalculatePSI(5, -1),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(),
                "Negative exclusion read count is invalid.");
        });
    }

    // C1 — only one effective length supplied ⇒ falls back to unnormalized I/(I+S). Source: ASSUMPTION 1.
    [Test]
    public void CalculatePSI_PartialEffectiveLength_FallsBackToUnnormalized()
    {
        double psi = TranscriptomeAnalyzer.CalculatePSI(80, 20, inclusionEffectiveLength: 200, exclusionEffectiveLength: 0);

        Assert.That(psi, Is.EqualTo(0.80).Within(Tol),
            "Normalization requires both effective lengths > 0; otherwise Ψ = I/(I+S) = 0.80 (ASSUMPTION 1).");
    }

    #endregion

    #region DetectAlternativeSplicing

    // M6 — Skipped exon: B lacks the middle exon (200,300) present in A. Source: Wang et al. (2008).
    [Test]
    public void DetectAlternativeSplicing_MiddleExonMissingInOneIsoform_ClassifiesSkippedExon()
    {
        var isoforms = new[]
        {
            Iso("A", "G1", (1, 100), (200, 300), (400, 500)),
            Iso("B", "G1", (1, 100), (400, 500)),
        };

        var events = TranscriptomeAnalyzer.DetectAlternativeSplicing(isoforms).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(events, Has.Count.EqualTo(1), "One isoform pair differing by one cassette exon yields one event.");
            Assert.That(events[0].EventType, Is.EqualTo(nameof(TranscriptomeAnalyzer.SplicingEventType.SkippedExon)),
                "An exon present in one isoform and absent from the other is a skipped/cassette exon (Wang 2008).");
        });
    }

    // M7 — Retained intron: B's single exon (1,300) spans the intron between A's (1,100) and (200,300).
    // Source: Wang et al. (2008).
    [Test]
    public void DetectAlternativeSplicing_ExonSpanningIntron_ClassifiesRetainedIntron()
    {
        var isoforms = new[]
        {
            Iso("A", "G1", (1, 100), (200, 300)),
            Iso("B", "G1", (1, 300)),
        };

        var events = TranscriptomeAnalyzer.DetectAlternativeSplicing(isoforms).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(events, Has.Count.EqualTo(1), "One difference → one event.");
            Assert.That(events[0].EventType, Is.EqualTo(nameof(TranscriptomeAnalyzer.SplicingEventType.RetainedIntron)),
                "An exon bridging the intron between two exons of the other isoform is intron retention (Wang 2008).");
        });
    }

    // M8 — Alternative 3' splice site: shared start 200, different end (300 vs 350). Source: Wang et al. (2008).
    [Test]
    public void DetectAlternativeSplicing_SharedStartDifferentEnd_ClassifiesAlternativeThreePrimeSS()
    {
        var isoforms = new[]
        {
            Iso("A", "G1", (1, 100), (200, 300)),
            Iso("B", "G1", (1, 100), (200, 350)),
        };

        var events = TranscriptomeAnalyzer.DetectAlternativeSplicing(isoforms).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(events, Has.Count.EqualTo(1), "One difference → one event.");
            Assert.That(events[0].EventType, Is.EqualTo(nameof(TranscriptomeAnalyzer.SplicingEventType.AlternativeThreePrimeSS)),
                "Unique exons sharing the 5' boundary but differing 3' end → A3SS (Wang 2008).");
        });
    }

    // M9 — Alternative 5' splice site: shared end 300, different start (1 vs ... ) on the donor exon.
    // A=(1,100); B=(1,150) share start 1 differ end → A3SS; use shared END to force A5SS:
    // A=(200,300) and B=(150,300) share end 300, differ start. Source: Wang et al. (2008).
    [Test]
    public void DetectAlternativeSplicing_SharedEndDifferentStart_ClassifiesAlternativeFivePrimeSS()
    {
        var isoforms = new[]
        {
            Iso("A", "G1", (1, 100), (200, 300)),
            Iso("B", "G1", (1, 100), (150, 300)),
        };

        var events = TranscriptomeAnalyzer.DetectAlternativeSplicing(isoforms).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(events, Has.Count.EqualTo(1), "One difference → one event.");
            Assert.That(events[0].EventType, Is.EqualTo(nameof(TranscriptomeAnalyzer.SplicingEventType.AlternativeFivePrimeSS)),
                "Unique exons sharing the 3' boundary but differing 5' start → A5SS (Wang 2008).");
        });
    }

    // M10 — Mutually exclusive exons: each isoform uses exactly one of two non-overlapping middle exons.
    // Source: Wang et al. (2008).
    [Test]
    public void DetectAlternativeSplicing_OneOfTwoAlternativeExons_ClassifiesMutuallyExclusive()
    {
        var isoforms = new[]
        {
            Iso("A", "G1", (1, 100), (200, 300), (500, 600)),
            Iso("B", "G1", (1, 100), (350, 400), (500, 600)),
        };

        var events = TranscriptomeAnalyzer.DetectAlternativeSplicing(isoforms).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(events, Has.Count.EqualTo(1), "One pair, one MXE difference → one event.");
            Assert.That(events[0].EventType, Is.EqualTo(nameof(TranscriptomeAnalyzer.SplicingEventType.MutuallyExclusiveExons)),
                "Each isoform carries exactly one of two non-overlapping alternative exons → MXE (Wang 2008).");
        });
    }

    // S2 — INV-05: a gene with a single isoform cannot define an event. Source: Wang et al. (2008).
    [Test]
    public void DetectAlternativeSplicing_SingleIsoformPerGene_ReturnsNoEvents()
    {
        var isoforms = new[] { Iso("A", "G1", (1, 100), (200, 300)) };

        var events = TranscriptomeAnalyzer.DetectAlternativeSplicing(isoforms).ToList();

        Assert.That(events, Is.Empty,
            "INV-05: an AS event requires two isoforms of the same gene; one isoform yields nothing (Wang 2008).");
    }

    // S3 — identical isoforms have no structural difference → no event.
    [Test]
    public void DetectAlternativeSplicing_IdenticalIsoforms_ReturnsNoEvents()
    {
        var isoforms = new[]
        {
            Iso("A", "G1", (1, 100), (200, 300)),
            Iso("B", "G1", (1, 100), (200, 300)),
        };

        var events = TranscriptomeAnalyzer.DetectAlternativeSplicing(isoforms).ToList();

        Assert.That(events, Is.Empty,
            "Structurally identical isoforms differ in no exon, so no alternative-splicing event is reported.");
    }

    // S4 — null and empty inputs yield an empty result (tolerant input contract).
    [Test]
    public void DetectAlternativeSplicing_NullOrEmptyInput_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(TranscriptomeAnalyzer.DetectAlternativeSplicing(null!), Is.Empty,
                "Null isoform input is treated as an empty sequence.");
            Assert.That(TranscriptomeAnalyzer.DetectAlternativeSplicing(
                Array.Empty<TranscriptomeAnalyzer.TranscriptIsoform>()), Is.Empty,
                "Empty isoform input yields no events.");
        });
    }

    // INV-05 — the detected event references the gene of the compared isoforms.
    [Test]
    public void DetectAlternativeSplicing_DetectedEvent_CarriesGeneId()
    {
        var isoforms = new[]
        {
            Iso("A", "GENEX", (1, 100), (200, 300), (400, 500)),
            Iso("B", "GENEX", (1, 100), (400, 500)),
        };

        var events = TranscriptomeAnalyzer.DetectAlternativeSplicing(isoforms).ToList();

        Assert.That(events[0].GeneId, Is.EqualTo("GENEX"),
            "INV-05: the event is attributed to the gene whose isoforms were compared.");
    }

    #endregion
}
