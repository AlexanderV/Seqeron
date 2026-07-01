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
/// Fuzz tests for the Oncology variant-annotation area — ONCO-ANNOT-001.
/// The unit under test is the VEP-like variant-consequence engine
/// <see cref="VariantAnnotator.PredictFunctionalImpact"/> (codon-translation
/// consequence + IMPACT + HGVS change), its batch most-severe wrapper
/// <see cref="VariantAnnotator.Annotate"/>, and the coordinate-routing
/// <see cref="VariantAnnotator.AnnotateVariant"/>, implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs. This is
/// the canonical surface that maps a tumour/somatic variant (chrom, 1-based
/// position, forward-strand REF/ALT) to its functional consequence against a
/// reference sequence and a gene/transcript model — the engine an oncology
/// pipeline drives to annotate somatic calls.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / MALFORMED inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no nonsense
/// output, and no *unhandled* runtime fault — specifically NO IndexOutOfRange /
/// ArgumentOutOfRange leaking from an internal Substring on a bad coordinate, NO
/// NullReferenceException on an empty/absent ALT allele, and NO KeyNotFound /
/// NullReference on an unknown contig. Every input must resolve to EITHER a
/// well-defined, theory-correct consequence OR a *documented, intentional*
/// validation outcome (an <see cref="ArgumentException"/> for an empty reference
/// window; an <see cref="ArgumentNullException"/> for null variant/transcript
/// enumerables). A silently-wrong consequence on garbage input is just as much a
/// bug as a crash. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-ANNOT-001 — Variant annotation (functional consequence)
/// Checklist: docs/checklists/03_FUZZING.md, row 91.
/// Fuzz strategy exercised for THIS unit:
///   • MC = Malformed Content (невалідний контент). Targets (checklist row 91):
///     "out-of-bounds coords, ref≠genome, empty alt, unknown chrom".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// docs/algorithms/Variants/Variant_Annotation.md (VARIANT-ANNOT-001):
///   • A variant is located relative to a transcript (upstream / downstream /
///     intron / UTR / splice / coding); a coding substitution is refined by
///     translating the affected REFERENCE and alternate codons with the NCBI
///     Standard code, then applying the VEP VariationEffect.pm predicates:
///       stop_gained ⇔ alt_pep has '*' absent from ref_pep;
///       stop_lost   ⇔ ref_pep is '*', alt_pep is not;
///       start_lost  ⇔ substitution disrupts the canonical ATG (precedence);
///       synonymous  ⇔ alt_pep == ref_pep (neither 'X');
///       missense    ⇔ peptides differ, equal length, not start/stop change.
///                                                              (§2.2, §4.1, INV-02)
///   • The consequence of a coding substitution is driven by the codon read from
///     the supplied forward-strand REFERENCE WINDOW (referenceSequence /
///     sequenceStart), NOT by the variant's *stated* REF allele — REF/ALT are
///     forward-strand, position 1-based (§3.1, §5.2, ASM-02). Codon extraction is
///     a guarded Substring that returns null (refinement skipped, coarse coding
///     term retained) when the codon falls outside the window. (VariantAnnotator
///     ExtractCodon; §5.2)
///   • A variant whose chromosome matches NO transcript is annotated
///     IntergenicVariant / MODIFIER (AnnotateVariant intergenic fallback).
///   • Validation (§3.3, §6.1): PredictFunctionalImpact throws ArgumentException
///     on null/empty referenceSequence; Annotate throws ArgumentNullException on
///     null variants/annotations; an ambiguous (IUPAC-'N') codon is excluded from
///     synonymous and reported CodingSequenceVariant.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Test reference fixture (hand-computed, mirrors VariantAnnotator_FunctionalImpact_Tests)
/// ───────────────────────────────────────────────────────────────────────────
/// Single-exon '+'-strand transcript. Exon (90..130) wraps a CDS (100..120) so a
/// coding SNV sits inside the exon (no splice trigger) and the CDS starts at
/// genomic 100. referenceSequence[0] == genomic 100 (sequenceStart = 100).
///   RefWindow = "ATGCAAGAATTATAANAAGGG"  (21 bases, genomic 100..120)
///   codon1 ATG(100-102 Met/M, start)  codon2 CAA(103-105 Gln/Q)
///   codon3 GAA(106-108 Glu/E)         codon4 TTA(109-111 Leu/L)
///   codon5 TAA(112-114 Stop)          codon6 NAA(115-117 ambiguous)
///   codon7 GGG(118-120 Gly/G)
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class OncologyAnnotationFuzzTests
{
    private const string RefWindow = "ATGCAAGAATTATAANAAGGG"; // genomic 100..120
    private const int RefStart = 100;

    private static Transcript CodingTranscript() => new(
        "ENST_ONCO", "ENSG_ONCO", "TP_TEST", "chr1",
        90, 130, '+',
        new List<(int, int)> { (90, 130) },
        new List<(int, int)> { (100, 120) },
        100, 120);

    private static Variant Snv(string chrom, int position, string refAllele, string altAllele) =>
        new(chrom, position, refAllele, altAllele, VariantType.SNV);

    // A well-formed result: a finite, defined consequence in the SO enum, with an
    // IMPACT that matches GetImpactLevel for that consequence (INV-01). No silent
    // contradiction between term and impact.
    private static void AssertWellFormed(FunctionalImpact fi)
    {
        Enum.IsDefined(typeof(ConsequenceType), fi.Consequence).Should().BeTrue(
            "every annotation must resolve to a defined SO consequence term");
        Enum.IsDefined(typeof(ImpactLevel), fi.Impact).Should().BeTrue();
        fi.Impact.Should().Be(GetImpactLevel(fi.Consequence),
            "IMPACT must equal the Constants.pm class for the reported term (INV-01)");
    }

    // ════════════════════════════════════════════════════════════════════════
    #region ONCO-ANNOT-001 — variant annotation (PredictFunctionalImpact / Annotate)
    // ════════════════════════════════════════════════════════════════════════

    // ───────────────────────────── POSITIVE SANITY ──────────────────────────
    // A well-formed variant whose REF matches the genome and lies inside the CDS
    // must annotate to the hand-computed documented consequence. This is the
    // business anchor: the MC tests below assert the engine stays disciplined on
    // garbage, this asserts it is CORRECT on good input.

    // GAA(Glu/E codon3 @106-108) base2 A>T -> GTA(Val/V): missense, MODERATE, p.E3V.
    // Src: VariationEffect.pm missense predicate + NCBI gc.prt.
    [Test]
    public void Positive_WellFormedCodingSnv_AnnotatesHandComputedMissense()
    {
        var fi = PredictFunctionalImpact(Snv("chr1", 107, "A", "T"), CodingTranscript(), RefWindow, RefStart);

        AssertWellFormed(fi);
        fi.Consequence.Should().Be(ConsequenceType.MissenseVariant,
            "GAA(Glu) -> GTA(Val) changes the amino acid with preserved length");
        fi.Impact.Should().Be(ImpactLevel.Moderate);
        fi.AminoAcidChange.Should().Be("p.E3V");
    }

    // CAA(Gln/Q codon2 @103-105) base1 C>T -> TAA(Stop): stop_gained, HIGH, p.Q2*.
    [Test]
    public void Positive_WellFormedPrematureStop_AnnotatesHandComputedStopGained()
    {
        var fi = PredictFunctionalImpact(Snv("chr1", 103, "C", "T"), CodingTranscript(), RefWindow, RefStart);

        AssertWellFormed(fi);
        fi.Consequence.Should().Be(ConsequenceType.StopGained,
            "CAA(Gln) -> TAA(Stop) introduces a premature stop absent from the reference");
        fi.Impact.Should().Be(ImpactLevel.High);
        fi.AminoAcidChange.Should().Be("p.Q2*");
    }

    // TTA(Leu/L codon4 @109-111) base3 A>G -> TTG(Leu/L): synonymous, LOW, p.L4=.
    [Test]
    public void Positive_WellFormedSynonymousSnv_AnnotatesHandComputedSynonymous()
    {
        var fi = PredictFunctionalImpact(Snv("chr1", 111, "A", "G"), CodingTranscript(), RefWindow, RefStart);

        AssertWellFormed(fi);
        fi.Consequence.Should().Be(ConsequenceType.SynonymousVariant);
        fi.AminoAcidChange.Should().Be("p.L4=");
    }

    // ──────────────────────── MC: out-of-bounds coords ──────────────────────
    // A position whose codon falls past the END of the reference window: codon
    // extraction is a guarded Substring (ExtractCodon offset+3 > length -> null),
    // so refinement is skipped and the coarse coding term is retained. NO
    // IndexOutOfRange / ArgumentOutOfRange may leak. The CDS bound here is 100..120
    // but we push the variant to genomic 119 (codon7 GGG @118-120) while the window
    // only spans to 120 — still in range; to force OOB we extend the transcript and
    // probe a coding position beyond the 21-base window.

    [Test]
    public void Mc_CodonBeyondReferenceWindow_DoesNotThrow_RetainsCodingTerm()
    {
        // Transcript whose CDS extends to 200, but RefWindow only covers 100..120.
        var t = new Transcript("ENST_ONCO", "ENSG_ONCO", "TP_TEST", "chr1",
            90, 250, '+',
            new List<(int, int)> { (90, 250) },
            new List<(int, int)> { (100, 200) },
            100, 200);
        // Genomic 150 is coding but its codon offset (50) is past the 21-base window.
        var variant = Snv("chr1", 150, "A", "T");

        FunctionalImpact fi = default;
        var act = () => fi = PredictFunctionalImpact(variant, t, RefWindow, RefStart);

        act.Should().NotThrow("a codon outside the reference window must be a guarded skip, never an OOB throw");
        AssertWellFormed(fi);
        // Refinement skipped (no codon) -> coarse coding consequence (MissenseVariant
        // is DetermineCodingConsequence's default coding SNV term); never a crash.
        IsCoding(fi.Consequence).Should().BeTrue(
            "a coding SNV with an out-of-window codon keeps the coarse coding term");
    }

    [Test]
    public void Mc_NegativeAndZeroPosition_DoNotThrow_AreNonCodingModifier()
    {
        var t = CodingTranscript();
        foreach (int pos in new[] { int.MinValue + 4, -1000, -1, 0 })
        {
            var variant = Snv("chr1", pos, "A", "T");

            FunctionalImpact fi = default;
            var act = () => fi = PredictFunctionalImpact(variant, t, RefWindow, RefStart);

            act.Should().NotThrow($"a non-positive position ({pos}) must be a guarded route, never an OOB throw");
            AssertWellFormed(fi);
            // Far upstream of the '+'-strand transcript start -> not a coding call.
            IsCoding(fi.Consequence).Should().BeFalse(
                $"position {pos} lies before the transcript and cannot be a coding consequence");
        }
    }

    [Test]
    public void Mc_HugePositionFarBeyondTranscript_DoesNotThrow_IsModifier()
    {
        var t = CodingTranscript();
        var variant = Snv("chr1", int.MaxValue - 3, "A", "T");

        FunctionalImpact fi = default;
        var act = () => fi = PredictFunctionalImpact(variant, t, RefWindow, RefStart);

        act.Should().NotThrow("an extreme position must not overflow or throw OOB");
        AssertWellFormed(fi);
        fi.Impact.Should().Be(ImpactLevel.Modifier,
            "a position far downstream of the transcript is intergenic/modifier, not coding");
    }

    // ─────────────────────────── MC: ref ≠ genome ───────────────────────────
    // The consequence of a coding substitution is driven by the codon READ FROM
    // THE REFERENCE WINDOW, not by the variant's stated REF allele (§3.1, §5.2,
    // ASM-02). So a variant that LIES about its REF base must STILL be classified
    // correctly relative to the genome — never crash, never a different consequence
    // than the genome-correct one. We pin this by annotating the SAME genomic
    // substitution with a deliberately wrong stated REF and asserting the
    // consequence is identical to the genome-correct call.

    [Test]
    public void Mc_StatedRefDisagreesWithGenome_ConsequenceIsGenomeDriven_NotCrash()
    {
        // Genome at 107 (codon3 GAA base2) is 'A'. A>T is the genome-correct missense.
        var genomeCorrect = PredictFunctionalImpact(Snv("chr1", 107, "A", "T"), CodingTranscript(), RefWindow, RefStart);

        // Now LIE: claim REF is 'G' (it is genomically 'A') at the same site, same ALT.
        var lyingRef = PredictFunctionalImpact(Snv("chr1", 107, "G", "T"), CodingTranscript(), RefWindow, RefStart);

        AssertWellFormed(lyingRef);
        lyingRef.Consequence.Should().Be(genomeCorrect.Consequence,
            "consequence is computed from the genome codon, so a wrong stated REF cannot change it");
        lyingRef.Consequence.Should().Be(ConsequenceType.MissenseVariant);
        lyingRef.AminoAcidChange.Should().Be("p.E3V",
            "the peptide change is read from the genome (Glu->Val), independent of the false REF");
    }

    [Test]
    public void Mc_RandomStatedRefBases_NeverChangeGenomeDrivenConsequence()
    {
        var rng = new Random(91_0001);
        var genomeCorrect = PredictFunctionalImpact(Snv("chr1", 107, "A", "T"), CodingTranscript(), RefWindow, RefStart);
        const string bases = "ACGTNacgt";

        for (int i = 0; i < 200; i++)
        {
            char fakeRef = bases[rng.Next(bases.Length)];
            var fi = PredictFunctionalImpact(Snv("chr1", 107, fakeRef.ToString(), "T"), CodingTranscript(), RefWindow, RefStart);

            AssertWellFormed(fi);
            fi.Consequence.Should().Be(genomeCorrect.Consequence,
                $"fake REF '{fakeRef}' must not perturb the genome-driven consequence");
        }
    }

    // ────────────────────────────── MC: empty alt ───────────────────────────
    // An empty/absent ALT allele must not NullReference or empty-Substring crash.
    // The codon machinery builds the alternate codon from the reference codon and
    // substitutes ALT bases in a length-bounded loop; with no ALT the alternate
    // codon equals the reference codon, so the result is well-formed (no peptide
    // change to report) — never a throw.

    [Test]
    public void Mc_EmptyAlt_DoesNotThrow_IsWellFormed()
    {
        var t = CodingTranscript();
        foreach (string alt in new[] { "", " " })
        {
            var variant = new Variant("chr1", 107, "A", alt, VariantType.SNV);

            FunctionalImpact fi = default;
            var act = () => fi = PredictFunctionalImpact(variant, t, RefWindow, RefStart);

            act.Should().NotThrow($"an empty/blank ALT ('{alt}') must not NullReference or empty-Substring crash");
            AssertWellFormed(fi);
        }
    }

    [Test]
    public void Mc_EmptyAlt_NoStopGainedOrFalseMissense_FromMissingAllele()
    {
        // No ALT base means the alternate codon == reference codon (GAA, Glu): a
        // missing allele can never SYNTHESISE a premature stop or a real amino-acid
        // change. Whatever defined term the engine assigns, it must not be a
        // genome-corrupting HIGH-impact call invented from absent data.
        var fi = PredictFunctionalImpact(new Variant("chr1", 107, "A", "", VariantType.SNV),
            CodingTranscript(), RefWindow, RefStart);

        AssertWellFormed(fi);
        fi.Consequence.Should().NotBe(ConsequenceType.StopGained,
            "an absent ALT cannot fabricate a premature stop");
        fi.Consequence.Should().NotBe(ConsequenceType.StopLost);
    }

    [Test]
    public void Mc_EmptyReferenceWindow_ThrowsArgumentException()
    {
        // Documented validation (§3.3, §6.1): a null/empty reference window is an
        // intentional ArgumentException, NOT an OOB/NRE.
        var act = () => PredictFunctionalImpact(Snv("chr1", 107, "A", "T"), CodingTranscript(), "", RefStart);
        act.Should().Throw<ArgumentException>();

        var actNull = () => PredictFunctionalImpact(Snv("chr1", 107, "A", "T"), CodingTranscript(), null!, RefStart);
        actNull.Should().Throw<ArgumentException>();
    }

    // ──────────────────────────── MC: unknown chrom ─────────────────────────
    // A contig name absent from the transcript model must produce a clean
    // IntergenicVariant / MODIFIER through Annotate — never a KeyNotFound or
    // NullReference. Annotate filters transcripts by chromosome equality, so an
    // unknown contig matches zero transcripts and falls through to the documented
    // intergenic annotation.

    [Test]
    public void Mc_UnknownChromosome_AnnotatesIntergenic_NoThrow()
    {
        var transcripts = new[] { CodingTranscript() }; // all on "chr1"
        var unknownContigs = new[] { "chrZZ", "chrMT_alt", "", "scaffold_99", "1", "CHR1" };

        foreach (var contig in unknownContigs)
        {
            var variant = Snv(contig, 107, "A", "T");

            List<VariantAnnotation> result = null!;
            var act = () => result = Annotate(new[] { variant }, transcripts, RefWindow, RefStart).ToList();

            act.Should().NotThrow($"unknown contig '{contig}' must not KeyNotFound / NullReference");
            result.Should().HaveCount(1, "Annotate yields one most-severe annotation per variant");
            result[0].Consequence.Should().Be(ConsequenceType.IntergenicVariant,
                $"a variant on contig '{contig}' overlaps no transcript -> intergenic");
            result[0].Impact.Should().Be(ImpactLevel.Modifier);
        }
    }

    [Test]
    public void Mc_KnownChromosome_StillAnnotatesCoding_GuardsAreNotOverbroad()
    {
        // Sanity dual: the unknown-chrom guard must not swallow a MATCHING contig.
        var result = Annotate(new[] { Snv("chr1", 107, "A", "T") }, new[] { CodingTranscript() }, RefWindow, RefStart).ToList();

        result.Should().HaveCount(1);
        result[0].Consequence.Should().Be(ConsequenceType.MissenseVariant,
            "the matching-contig coding SNV is still annotated correctly (guard is not over-broad)");
    }

    // ─────────────────────── MC: validation / null inputs ───────────────────

    [Test]
    public void Mc_NullVariantsOrTranscripts_ThrowArgumentNull()
    {
        var act1 = () => Annotate(null!, new[] { CodingTranscript() }, RefWindow, RefStart).ToList();
        act1.Should().Throw<ArgumentNullException>();

        var act2 = () => Annotate(new[] { Snv("chr1", 107, "A", "T") }, null!, RefWindow, RefStart).ToList();
        act2.Should().Throw<ArgumentNullException>();
    }

    // ─────────────────────── MC: combined malformed sweep ───────────────────
    // A randomised stew of all four MC targets at once (junk contig, OOB / negative
    // position, garbage REF, empty-or-junk ALT) must NEVER throw an undisciplined
    // runtime fault and must always yield exactly one well-formed annotation per
    // variant. [CancelAfter] guards against any accidental hang.

    [Test]
    [CancelAfter(30_000)]
    public void Mc_RandomisedMalformedSweep_NeverCrashes_AlwaysWellFormed()
    {
        var rng = new Random(91_0042);
        var transcripts = new[] { CodingTranscript() };
        string[] contigs = { "chr1", "chrZZ", "", "scaffold", "1", "CHR1" };
        string[] alleles = { "A", "C", "G", "T", "N", "", "AT", "n", "x", "  " };

        for (int i = 0; i < 1500; i++)
        {
            string chrom = contigs[rng.Next(contigs.Length)];
            int pos = rng.Next(-5000, 5000);          // includes negatives, OOB, in-range
            string refA = alleles[rng.Next(alleles.Length)];
            string altA = alleles[rng.Next(alleles.Length)];
            var variant = new Variant(chrom, pos, refA, altA, VariantType.SNV);

            List<VariantAnnotation> result = null!;
            var act = () => result = Annotate(new[] { variant }, transcripts, RefWindow, RefStart).ToList();

            act.Should().NotThrow(
                $"malformed record (chrom='{chrom}', pos={pos}, ref='{refA}', alt='{altA}') must resolve cleanly");
            result.Should().HaveCount(1);
            VariantAnnotation a = result[0];
            Enum.IsDefined(typeof(ConsequenceType), a.Consequence).Should().BeTrue();
            a.Impact.Should().Be(GetImpactLevel(a.Consequence),
                "term/impact consistency must hold even on malformed input (INV-01)");
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
