using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using static Seqeron.Genomics.Annotation.VariantAnnotator;
using Variant = Seqeron.Genomics.Annotation.VariantAnnotator.Variant;
using VariantType = Seqeron.Genomics.Annotation.VariantAnnotator.VariantType;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the canonical Variant-annotation unit — VARIANT-ANNOT-001.
/// The unit under test is the VEP-like functional-consequence engine in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs:
/// <see cref="VariantAnnotator.PredictFunctionalImpact"/> (codon-translation
/// consequence + IMPACT + HGVS change), its batch most-severe wrapper
/// <see cref="VariantAnnotator.Annotate"/>, the coordinate-routing
/// <see cref="VariantAnnotator.AnnotateVariant"/>, and the supporting
/// <see cref="VariantAnnotator.GetImpactLevel"/> / <see cref="VariantAnnotator.GetConsequenceRank"/>
/// mappings. This is the surface that maps a variant (chrom, 1-based position,
/// forward-strand REF/ALT) to a Sequence Ontology consequence term + IMPACT
/// against a reference window and a transcript model.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / MALFORMED inputs and asserts the code
/// NEVER fails in an undisciplined way: no hang, no fabricated output, no
/// *unhandled* runtime fault — specifically NO IndexOutOfRange / ArgumentOutOfRange
/// leaking from an internal Substring on an out-of-bounds coordinate, NO
/// NullReferenceException on an empty / absent input, and NO fabricated /
/// undefined consequence term. Every input must resolve to EITHER a well-defined,
/// theory-correct consequence OR a *documented, intentional* validation outcome
/// (an <see cref="ArgumentException"/> for an empty reference window; an
/// <see cref="ArgumentNullException"/> for null variant/transcript enumerables).
/// A silently-wrong or fabricated consequence on garbage input is just as much a
/// bug as a crash. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: VARIANT-ANNOT-001 — Variant annotation (functional consequence)
/// Checklist: docs/checklists/03_FUZZING.md, row 186.
/// Fuzz strategy exercised for THIS unit:
///   • MC = Malformed Content (невалідний контент). Targets (checklist row 186):
///     "out-of-bounds, unknown consequence, empty".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// Scope vs row 91 (OncologyAnnotationFuzzTests.cs, ONCO-ANNOT-001): the SAME
/// canonical annotation unit. Row 91 focused its MC sweep on "ref≠genome" and
/// "unknown chrom". This file is COMPLEMENTARY and goes deeper on the three
/// row-186 MC targets without duplicating row 91's exact assertions:
///   • OUT-OF-BOUNDS — not just OOB position, but reference-WINDOW alignment
///     boundaries: codon partly/wholly before the window (negative offset),
///     codon straddling the window end, sequenceStart shifted so the in-CDS
///     codon falls out of range. The documented guard is ExtractCodon → null
///     (refinement skipped, coarse coding term retained), never an OOB throw.
///   • UNKNOWN CONSEQUENCE — every reported term must be a DEFINED member of the
///     SO ConsequenceType enum with a known severity rank (GetConsequenceRank ≠
///     int.MaxValue) and a term/IMPACT pair consistent with GetImpactLevel; an
///     input that overlaps NO feature defaults to IntergenicVariant/MODIFIER
///     (the documented default consequence), never a fabricated term; an
///     untranslatable/ambiguous codon is reported CodingSequenceVariant, NOT an
///     invented SO term and NOT synonymous (§3.3, §6.1 / INV-02).
///   • EMPTY — empty/null referenceSequence ⇒ ArgumentException; empty ALT
///     allele ⇒ no NRE / no empty-Substring crash (verifies the row-91
///     TranslateOrUnknown fix still holds); empty TRANSCRIPT set ⇒ intergenic;
///     empty VARIANT set ⇒ empty result; null enumerables ⇒ ArgumentNullException.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// docs/algorithms/Variants/Variant_Annotation.md (VARIANT-ANNOT-001):
///   • A variant is located relative to a transcript (upstream/downstream/intron/
///     UTR/splice/coding); a coding substitution is refined by translating the
///     affected REFERENCE and alternate codons with the NCBI Standard code, then
///     applying the VEP VariationEffect.pm predicates (§2.2, §4.1, INV-02).
///   • Codon extraction is a GUARDED Substring (ExtractCodon) returning null —
///     refinement skipped, coarse coding term retained — when the codon falls
///     outside the supplied window (§5.2, ASM-02).
///   • A variant overlapping NO transcript ⇒ IntergenicVariant / MODIFIER
///     (AnnotateVariant intergenic fallback; documented default consequence).
///   • Validation (§3.3, §6.1): PredictFunctionalImpact throws ArgumentException
///     on null/empty referenceSequence; Annotate throws ArgumentNullException on
///     null variants/annotations; an IUPAC-ambiguous codon is excluded from
///     synonymous and reported CodingSequenceVariant.
///   • INV-01: IMPACT(term) equals the Constants.pm class for that term.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Test reference fixture (hand-computed; mirrors the canonical FunctionalImpact tests)
/// ───────────────────────────────────────────────────────────────────────────
/// Single-exon '+'-strand transcript. Exon (90..130) wraps a CDS (100..120) so a
/// coding SNV sits inside the exon and the CDS starts at genomic 100.
/// referenceSequence[0] == genomic 100 (sequenceStart = 100).
///   RefWindow = "ATGCAAGAATTATAANAAGGG"  (21 bases, genomic 100..120)
///   codon1 ATG(100-102 Met/M start)  codon2 CAA(103-105 Gln/Q)
///   codon3 GAA(106-108 Glu/E)        codon4 TTA(109-111 Leu/L)
///   codon5 TAA(112-114 Stop)         codon6 NAA(115-117 ambiguous)
///   codon7 GGG(118-120 Gly/G)
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class VariantAnnotationFuzzTests
{
    private const string RefWindow = "ATGCAAGAATTATAANAAGGG"; // genomic 100..120
    private const int RefStart = 100;

    private static Transcript CodingTranscript() => new(
        "ENST_VA", "ENSG_VA", "VA_TEST", "chr1",
        90, 130, '+',
        new List<(int, int)> { (90, 130) },
        new List<(int, int)> { (100, 120) },
        100, 120);

    private static Variant Snv(string chrom, int position, string refAllele, string altAllele) =>
        new(chrom, position, refAllele, altAllele, VariantType.SNV);

    // A well-formed result: the consequence is a DEFINED SO term with a KNOWN
    // severity rank, and the IMPACT matches GetImpactLevel for that term (INV-01).
    // "Known rank" is the no-fabrication anchor: an undefined/invented term would
    // have rank int.MaxValue (GetConsequenceRank fallback).
    private static void AssertWellFormed(FunctionalImpact fi)
    {
        Enum.IsDefined(typeof(ConsequenceType), fi.Consequence).Should().BeTrue(
            "every annotation must resolve to a DEFINED SO consequence term, never a fabricated one");
        GetConsequenceRank(fi.Consequence).Should().NotBe(int.MaxValue,
            "a reported term must carry a known Constants.pm severity rank (no fabricated/unranked term)");
        Enum.IsDefined(typeof(ImpactLevel), fi.Impact).Should().BeTrue();
        fi.Impact.Should().Be(GetImpactLevel(fi.Consequence),
            "IMPACT must equal the Constants.pm class for the reported term (INV-01)");
    }

    private static void AssertWellFormed(VariantAnnotation a)
    {
        Enum.IsDefined(typeof(ConsequenceType), a.Consequence).Should().BeTrue(
            "every annotation must resolve to a DEFINED SO consequence term, never a fabricated one");
        GetConsequenceRank(a.Consequence).Should().NotBe(int.MaxValue,
            "a reported term must carry a known Constants.pm severity rank (no fabricated/unranked term)");
        a.Impact.Should().Be(GetImpactLevel(a.Consequence),
            "IMPACT must equal the Constants.pm class for the reported term (INV-01)");
    }

    // ════════════════════════════════════════════════════════════════════════
    #region VARIANT-ANNOT-001 — variant annotation (PredictFunctionalImpact / Annotate)
    // ════════════════════════════════════════════════════════════════════════

    // ───────────────────────────── POSITIVE SANITY ──────────────────────────
    // The business anchors: the engine must be CORRECT on good input before we
    // assert it stays DISCIPLINED on garbage. Three documented contracts:
    //   (a) a well-formed in-feature coding SNV → the hand-computed consequence;
    //   (b) an out-of-feature (intergenic) position → the documented default term;
    //   (c) the matching-contig coding call still works (guards not over-broad).

    // GAA(Glu/E codon3 @106-108) base2 A>T -> GTA(Val/V): missense, MODERATE, p.E3V.
    // Src: VariationEffect.pm missense predicate + NCBI gc.prt.
    [Test]
    public void Positive_WellFormedInFeatureSnv_AnnotatesHandComputedMissense()
    {
        var fi = PredictFunctionalImpact(Snv("chr1", 107, "A", "T"), CodingTranscript(), RefWindow, RefStart);

        AssertWellFormed(fi);
        fi.Consequence.Should().Be(ConsequenceType.MissenseVariant,
            "GAA(Glu) -> GTA(Val) changes the amino acid with preserved length");
        fi.Impact.Should().Be(ImpactLevel.Moderate);
        fi.AminoAcidChange.Should().Be("p.E3V");
    }

    [Test]
    public void Positive_IntergenicPosition_AnnotatesDocumentedDefaultConsequence()
    {
        // A position far beyond the transcript (>5kb upstream / >500bp downstream)
        // overlaps no feature → the DOCUMENTED default consequence is intergenic.
        var result = Annotate(new[] { Snv("chr1", 9_000_000, "A", "T") }, new[] { CodingTranscript() },
            RefWindow, RefStart).ToList();

        result.Should().HaveCount(1);
        AssertWellFormed(result[0]);
        result[0].Consequence.Should().Be(ConsequenceType.IntergenicVariant,
            "a variant overlapping no transcript falls through to the documented intergenic default");
        result[0].Impact.Should().Be(ImpactLevel.Modifier);
    }

    [Test]
    public void Positive_MatchingContigCodingCall_StillWorks()
    {
        var result = Annotate(new[] { Snv("chr1", 107, "A", "T") }, new[] { CodingTranscript() },
            RefWindow, RefStart).ToList();

        result.Should().HaveCount(1);
        AssertWellFormed(result[0]);
        result[0].Consequence.Should().Be(ConsequenceType.MissenseVariant,
            "the documented guards must not swallow a genuine in-feature coding SNV");
    }

    // ════════════════════════════════════════════════════════════════════════
    // MC TARGET 1 — OUT-OF-BOUNDS
    // The codon read is a guarded Substring (ExtractCodon: offset<0 OR offset+3 >
    // length ⇒ null). When the codon falls outside the window the refinement is
    // skipped and the coarse coding term is retained — NEVER an OOB throw. We
    // probe every boundary of that guard.
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void Mc_OutOfBounds_CodonStraddlesWindowEnd_DoesNotThrow_RetainsCodingTerm()
    {
        // CDS extends to 200 but RefWindow only covers 100..120. A coding SNV at
        // genomic 150 has codon offset 50, well past the 21-base window end.
        var t = new Transcript("ENST_VA", "ENSG_VA", "VA_TEST", "chr1",
            90, 250, '+',
            new List<(int, int)> { (90, 250) },
            new List<(int, int)> { (100, 200) },
            100, 200);
        var variant = Snv("chr1", 150, "A", "T");

        FunctionalImpact fi = default;
        var act = () => fi = PredictFunctionalImpact(variant, t, RefWindow, RefStart);

        act.Should().NotThrow("a codon past the window end must be a guarded skip, never an OOB throw");
        AssertWellFormed(fi);
        IsCoding(fi.Consequence).Should().BeTrue(
            "a coding SNV with an out-of-window codon retains the coarse coding term (refinement skipped)");
    }

    [Test]
    public void Mc_OutOfBounds_SequenceStartShiftedSoCodonIsBeforeWindow_DoesNotThrow()
    {
        // The window's declared start is shifted to 130, so the genuinely-in-CDS
        // variant at 107 maps to a NEGATIVE offset (107-130 = -23). ExtractCodon's
        // offset<0 guard must catch it: a guarded skip, never a negative-index throw.
        var variant = Snv("chr1", 107, "A", "T");

        FunctionalImpact fi = default;
        var act = () => fi = PredictFunctionalImpact(variant, CodingTranscript(), RefWindow, sequenceStart: 130);

        act.Should().NotThrow("a codon before the window (negative offset) must be a guarded skip");
        AssertWellFormed(fi);
        IsCoding(fi.Consequence).Should().BeTrue(
            "negative-offset codon retains the coarse coding term, no refinement, no throw");
    }

    [Test]
    public void Mc_OutOfBounds_CodonEndsExactlyAtWindowEnd_IsInRange_NoThrow()
    {
        // Boundary: codon7 GGG spans genomic 118..120; offset 18, 18+3==21==length.
        // ExtractCodon's "offset+3 > length" guard is strict, so this is IN range
        // and DOES refine — the guard must not be off-by-one and reject the last codon.
        var variant = Snv("chr1", 118, "G", "A"); // GGG(Gly) -> AGG(Arg): missense p.G7R

        FunctionalImpact fi = default;
        var act = () => fi = PredictFunctionalImpact(variant, CodingTranscript(), RefWindow, RefStart);

        act.Should().NotThrow();
        AssertWellFormed(fi);
        fi.Consequence.Should().Be(ConsequenceType.MissenseVariant,
            "the last codon exactly filling the window must still refine (guard not off-by-one)");
        fi.AminoAcidChange.Should().Be("p.G7R");
    }

    [Test]
    public void Mc_OutOfBounds_ExtremeAndNonPositivePositions_DoNotThrow_AreNotCoding()
    {
        var t = CodingTranscript();
        foreach (int pos in new[] { int.MinValue + 4, -1_000_000, -1, 0, int.MaxValue - 3 })
        {
            FunctionalImpact fi = default;
            var act = () => fi = PredictFunctionalImpact(Snv("chr1", pos, "A", "T"), t, RefWindow, RefStart);

            act.Should().NotThrow($"an extreme/non-positive position ({pos}) must not overflow or OOB-throw");
            AssertWellFormed(fi);
            IsCoding(fi.Consequence).Should().BeFalse(
                $"position {pos} lies far outside the transcript and cannot be a coding consequence");
        }
    }

    // ════════════════════════════════════════════════════════════════════════
    // MC TARGET 2 — UNKNOWN CONSEQUENCE
    // No input may produce a fabricated/undefined SO term. An input that maps to
    // no known feature must default to the DOCUMENTED IntergenicVariant/MODIFIER;
    // an untranslatable/ambiguous codon must be CodingSequenceVariant (NOT a made-up
    // term, NOT synonymous). Every result carries a known severity rank.
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void Mc_UnknownConsequence_AmbiguousCodon_IsCodingSequenceVariant_NotSynonymous()
    {
        // codon6 NAA(@115-117) is IUPAC-ambiguous → translates to 'X'. A substitution
        // here cannot be synonymous (X excluded) and must NOT fabricate a term: the
        // documented outcome is CodingSequenceVariant (§3.3, §6.1, INV-02).
        var fi = PredictFunctionalImpact(Snv("chr1", 116, "A", "T"), CodingTranscript(), RefWindow, RefStart);

        AssertWellFormed(fi);
        fi.Consequence.Should().Be(ConsequenceType.CodingSequenceVariant,
            "an ambiguous (X) codon is excluded from synonymous and reported coding_sequence_variant");
        fi.Consequence.Should().NotBe(ConsequenceType.SynonymousVariant);
    }

    [Test]
    public void Mc_UnknownConsequence_OutOfAlphabetAltCodon_IsCodingSequenceVariant_NoCrash()
    {
        // An ALT carrying an out-of-alphabet base makes the alternate codon
        // untranslatable. The row-91 TranslateOrUnknown fix maps it to 'X' rather
        // than throwing, so it must surface as coding_sequence_variant — a DEFINED
        // term, never a fabricated one and never a crash.
        var rng = new Random(186_2001);
        const string junk = "XZ?#@$1*qK";

        for (int i = 0; i < 100; i++)
        {
            char bad = junk[rng.Next(junk.Length)];
            FunctionalImpact fi = default;
            var act = () => fi = PredictFunctionalImpact(
                Snv("chr1", 107, "A", bad.ToString()), CodingTranscript(), RefWindow, RefStart);

            act.Should().NotThrow($"out-of-alphabet ALT base '{bad}' must not crash the translator");
            AssertWellFormed(fi);
            fi.Consequence.Should().Be(ConsequenceType.CodingSequenceVariant,
                $"untranslatable alt codon (ALT '{bad}') → coding_sequence_variant, not a fabricated term");
        }
    }

    [Test]
    public void Mc_UnknownConsequence_NoFeatureOverlap_DefaultsToIntergenic_AcrossPositions()
    {
        // Sweep positions that overlap nothing on the transcript's chromosome. Each
        // must default to the DOCUMENTED intergenic term — never an undefined term,
        // never a coding call invented from thin air.
        var rng = new Random(186_2002);
        var transcripts = new[] { CodingTranscript() }; // chr1, 90..130

        for (int i = 0; i < 300; i++)
        {
            // Pick positions well outside the [Start-5000, End+500] influence window.
            int pos = rng.Next(0, 2) == 0 ? rng.Next(-50_000, -6_000) : rng.Next(2_000_000, 9_000_000);
            var result = Annotate(new[] { Snv("chr1", pos, "A", "T") }, transcripts, RefWindow, RefStart).ToList();

            result.Should().HaveCount(1);
            AssertWellFormed(result[0]);
            result[0].Consequence.Should().Be(ConsequenceType.IntergenicVariant,
                $"position {pos} overlaps no feature → documented intergenic default");
            result[0].Impact.Should().Be(ImpactLevel.Modifier);
        }
    }

    [Test]
    public void Mc_UnknownConsequence_AllResultsCarryKnownSeverityRank_AndDefinedImpact()
    {
        // No matter the (degenerate) input, the reported term is always a defined SO
        // member with a known Constants.pm rank and a self-consistent IMPACT. This is
        // the strongest no-fabrication guard across a broad input grid.
        var rng = new Random(186_2003);
        var t = CodingTranscript();
        string[] alleles = { "A", "C", "G", "T", "N", "", "AT", "n", "x", "  ", "*" };

        for (int i = 0; i < 1000; i++)
        {
            int pos = rng.Next(80, 140); // around/inside the transcript
            string refA = alleles[rng.Next(alleles.Length)];
            string altA = alleles[rng.Next(alleles.Length)];
            var variant = new Variant("chr1", pos, refA, altA, VariantType.SNV);

            FunctionalImpact fi = default;
            var act = () => fi = PredictFunctionalImpact(variant, t, RefWindow, RefStart);

            act.Should().NotThrow($"(pos={pos}, ref='{refA}', alt='{altA}') must resolve cleanly");
            AssertWellFormed(fi);
        }
    }

    // ════════════════════════════════════════════════════════════════════════
    // MC TARGET 3 — EMPTY
    // Empty/null reference window ⇒ ArgumentException (documented). Empty ALT ⇒ no
    // NRE / empty-Substring crash. Empty transcript set ⇒ intergenic. Empty variant
    // set ⇒ empty output. Null enumerables ⇒ ArgumentNullException.
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void Mc_Empty_NullOrEmptyReferenceWindow_ThrowsArgumentException()
    {
        var actEmpty = () => PredictFunctionalImpact(Snv("chr1", 107, "A", "T"), CodingTranscript(), "", RefStart);
        actEmpty.Should().Throw<ArgumentException>("a null/empty reference window is a documented ArgumentException (§3.3)");

        var actNull = () => PredictFunctionalImpact(Snv("chr1", 107, "A", "T"), CodingTranscript(), null!, RefStart);
        actNull.Should().Throw<ArgumentException>("a null reference window is a documented ArgumentException (§3.3)");
    }

    [Test]
    public void Mc_Empty_AltAllele_DoesNotThrow_NoFabricatedHighImpact()
    {
        // An empty/blank ALT must not NullReference or empty-Substring crash (the
        // row-91 fix). With no ALT base the alternate codon == reference codon, so a
        // missing allele can NEVER synthesise a premature stop or a real AA change.
        var t = CodingTranscript();
        foreach (string alt in new[] { "", " ", "   " })
        {
            FunctionalImpact fi = default;
            var act = () => fi = PredictFunctionalImpact(
                new Variant("chr1", 107, "A", alt, VariantType.SNV), t, RefWindow, RefStart);

            act.Should().NotThrow($"empty/blank ALT ('{alt}') must not NRE or empty-Substring crash");
            AssertWellFormed(fi);
            fi.Consequence.Should().NotBe(ConsequenceType.StopGained,
                "an absent ALT cannot fabricate a premature stop");
            fi.Consequence.Should().NotBe(ConsequenceType.StopLost,
                "an absent ALT cannot fabricate a stop-loss");
        }
    }

    [Test]
    public void Mc_Empty_TranscriptSet_AnnotatesIntergenic_NoThrow()
    {
        // No transcript to overlap → the documented intergenic fallback in
        // AnnotateVariant, never a crash on the empty relevant-transcript list.
        var result = Annotate(new[] { Snv("chr1", 107, "A", "T") },
            Array.Empty<Transcript>(), RefWindow, RefStart).ToList();

        result.Should().HaveCount(1, "an empty transcript set still yields one annotation per variant");
        AssertWellFormed(result[0]);
        result[0].Consequence.Should().Be(ConsequenceType.IntergenicVariant,
            "with no transcripts every variant is intergenic");
        result[0].Impact.Should().Be(ImpactLevel.Modifier);
    }

    [Test]
    public void Mc_Empty_VariantSet_YieldsEmptyResult_NoThrow()
    {
        List<VariantAnnotation> result = null!;
        var act = () => result = Annotate(Array.Empty<Variant>(), new[] { CodingTranscript() },
            RefWindow, RefStart).ToList();

        act.Should().NotThrow("an empty variant set must enumerate to an empty result, not throw");
        result.Should().BeEmpty();
    }

    [Test]
    public void Mc_Empty_NullVariantsOrTranscripts_ThrowArgumentNull()
    {
        var act1 = () => Annotate(null!, new[] { CodingTranscript() }, RefWindow, RefStart).ToList();
        act1.Should().Throw<ArgumentNullException>("null variants is a documented ArgumentNullException (§3.3)");

        var act2 = () => Annotate(new[] { Snv("chr1", 107, "A", "T") }, null!, RefWindow, RefStart).ToList();
        act2.Should().Throw<ArgumentNullException>("null annotations is a documented ArgumentNullException (§3.3)");
    }

    // ════════════════════════════════════════════════════════════════════════
    // MC COMBINED — all three targets at once
    // A randomised stew of OOB positions, empty/junk alleles and feature/no-feature
    // overlaps must NEVER throw and must always yield exactly one well-formed
    // annotation per variant whose term is defined, ranked and IMPACT-consistent.
    // [CancelAfter] guards against any accidental hang.
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    [CancelAfter(30_000)]
    public void Mc_RandomisedMalformedSweep_NeverCrashes_AlwaysWellFormed()
    {
        var rng = new Random(186_0042);
        var transcripts = new[] { CodingTranscript() };
        // Positions span far-negative, OOB, in-feature and far-downstream.
        int[] posBuckets = { -50_000, -1, 0, 95, 107, 116, 125, 5_000, 9_000_000 };
        string[] alleles = { "A", "C", "G", "T", "N", "", "AT", "n", "x", "  ", "*", "?" };

        for (int i = 0; i < 2000; i++)
        {
            int pos = rng.Next(0, 2) == 0
                ? posBuckets[rng.Next(posBuckets.Length)]
                : rng.Next(-60_000, 9_000_000);
            string refA = alleles[rng.Next(alleles.Length)];
            string altA = alleles[rng.Next(alleles.Length)];
            var variant = new Variant("chr1", pos, refA, altA, VariantType.SNV);

            List<VariantAnnotation> result = null!;
            var act = () => result = Annotate(new[] { variant }, transcripts, RefWindow, RefStart).ToList();

            act.Should().NotThrow(
                $"malformed record (pos={pos}, ref='{refA}', alt='{altA}') must resolve cleanly");
            result.Should().HaveCount(1);
            AssertWellFormed(result[0]);
        }
    }

    #endregion

    // Local mirror of the engine's internal coding-consequence predicate (the
    // source method is private). Kept in lockstep with VariantAnnotator
    // IsCodingConsequence so the fuzz assertions reason about the same set.
    private static bool IsCoding(ConsequenceType c) => c switch
    {
        ConsequenceType.StopGained => true,
        ConsequenceType.FrameshiftVariant => true,
        ConsequenceType.StopLost => true,
        ConsequenceType.StartLost => true,
        ConsequenceType.InframeInsertion => true,
        ConsequenceType.InframeDeletion => true,
        ConsequenceType.MissenseVariant => true,
        ConsequenceType.SynonymousVariant => true,
        ConsequenceType.ProteinAlteringVariant => true,
        _ => false
    };
}
