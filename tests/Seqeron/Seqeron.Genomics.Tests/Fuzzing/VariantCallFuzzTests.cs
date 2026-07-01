using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Variants area — VARIANT-CALL-001 (Variant Detection).
/// The unit under test is the deterministic, alignment-based variant caller
/// <see cref="VariantCaller.CallVariants"/> together with its column-scan core
/// <see cref="VariantCaller.CallVariantsFromAlignment"/>, the SNP/indel projections
/// (<see cref="VariantCaller.FindSnps"/>, <see cref="VariantCaller.FindSnpsDirect"/>,
/// <see cref="VariantCaller.FindInsertions"/>, <see cref="VariantCaller.FindDeletions"/>),
/// the point-mutation classifier <see cref="VariantCaller.ClassifyMutation"/>, the
/// aggregate <see cref="VariantCaller.CalculateTiTvRatio"/>, and the statistics roll-up
/// <see cref="VariantCaller.CalculateStatistics"/>; implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no nonsense output, and no *unhandled* runtime exception
/// (DivideByZero / NullReference / Overflow / NaN). Every input must resolve to
/// EITHER a well-defined, theory-correct value OR a *documented, intentional*
/// outcome (here, an ArgumentNullException for null inputs and an
/// ArgumentException for unequal aligned lengths). — docs/ADVANCED_TESTING_CHECKLIST.md
/// §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: VARIANT-CALL-001 — Variant Detection (Variants)
/// Checklist: docs/checklists/03_FUZZING.md, row 187.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// MAPPING of the generic checklist BE targets ("zero depth, tumor=normal,
/// all-N") onto THIS unit's contract. VARIANT-CALL-001 is an *alignment-based*
/// caller — it has no pileup, no read depth, no allele fraction, no quality and
/// no tumor/normal channels (that is a DIFFERENT unit, ONCO-SOMATIC-001, row 87).
/// The targets translate to the boundary degeneracies of the alignment-scan
/// contract:
///   • "zero depth"   → zero-information sites: empty aligned strings and empty
///       DnaSequences ⇒ EMPTY result (no call), and the only division in the
///       unit — VariantDensity = variants/refLen — must be guarded at refLen == 0
///       to return 0, NOT throw DivideByZero nor yield NaN (§5; CalculateStatistics).
///       The other aggregate, Ti/Tv with #Tv == 0, must return the documented 0
///       sentinel, NOT +∞/NaN (INV-06, §5.4).
///   • "tumor=normal" → the alignment analogue of "two identical evidence
///       channels": identical reference == query ⇒ ZERO variants, i.e. no false
///       positive call (INV-01, §2.4, §6.1).
///   • "all-N"        → a non-informative alphabet. DnaSequence validation rejects
///       'N' (A/C/G/T only), so the typed entry points throw a *documented*
///       ArgumentException; the raw-string entry point CallVariantsFromAlignment
///       accepts arbitrary chars, so all-N-vs-all-N ⇒ no variants and ref-vs-all-N
///       ⇒ SNPs that classify (purine/pyrimidine) WITHOUT crash or NaN.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (docs/algorithms/Variants/Variant_Detection.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • Column rule: ref-gap ⇒ Insertion; query-gap ⇒ Deletion; mismatch ⇒ SNP;
///       else match (no variant).                                       (§2.2, INV-04)
///   • Identical reference and query ⇒ zero variants.                   (INV-01, §6.1)
///   • Every emitted SNP has distinct single-base REF and ALT.          (INV-02)
///   • Every variant's 0-based Position ∈ [0, reference.Length].        (INV-03)
///   • ClassifyMutation = Transition iff {ref,alt}⊆{A,G} or ⊆{C,T};
///       else Transversion (for a SNP); Other for non-SNP.             (INV-05, §4.2)
///   • CalculateTiTvRatio = #Ti / #Tv, or 0 when #Tv == 0.             (INV-06, §5.4)
///   • Null reference/query ⇒ ArgumentNullException;
///       empty aligned input ⇒ empty result;
///       unequal aligned lengths ⇒ ArgumentException.                  (§3.3, §6.1)
///   • Classification is case-insensitive.                             (§3.3, §6.1)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class VariantCallFuzzTests
{
    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented per-variant contract on EVERY emitted variant, no
    // matter how degenerate the input: the 0-based Position must be a finite,
    // non-negative reference coordinate (INV-03), the type must come from the
    // documented set, and a SNP must have distinct single-base REF/ALT (INV-02)
    // while an indel must carry the "-" gap sentinel on exactly one side (§2.2).
    // This is what stops a fuzz test from rubber-stamping nonsense output green.
    private static void AssertWellFormedVariant(Variant v, int refLen)
    {
        v.Position.Should().BeGreaterThanOrEqualTo(0, "Position is a 0-based reference coordinate (INV-03)");
        v.Position.Should().BeLessThanOrEqualTo(refLen, "Position ∈ [0, reference.Length] (INV-03)");
        v.QueryPosition.Should().BeGreaterThanOrEqualTo(0, "QueryPosition is a 0-based query coordinate");
        v.Type.Should().BeOneOf(VariantType.SNP, VariantType.Insertion, VariantType.Deletion);

        switch (v.Type)
        {
            case VariantType.SNP:
                v.ReferenceAllele.Should().HaveLength(1).And.NotBe("-", "a SNP REF is a single base (INV-02)");
                v.AlternateAllele.Should().HaveLength(1).And.NotBe("-", "a SNP ALT is a single base (INV-02)");
                v.ReferenceAllele.Should().NotBe(v.AlternateAllele, "a SNP column is a mismatch (INV-02)");
                break;
            case VariantType.Insertion:
                v.ReferenceAllele.Should().Be("-", "an insertion has a ref-side gap (§2.2)");
                v.AlternateAllele.Should().NotBe("-", "an insertion carries the inserted query base (§2.2)");
                break;
            case VariantType.Deletion:
                v.ReferenceAllele.Should().NotBe("-", "a deletion carries the deleted reference base (§2.2)");
                v.AlternateAllele.Should().Be("-", "a deletion has a query-side gap (§2.2)");
                break;
        }
    }

    private static void AssertAllWellFormed(IReadOnlyList<Variant> variants, int refLen)
    {
        foreach (var v in variants)
            AssertWellFormedVariant(v, refLen);
    }

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-CALL-001 — Variant Detection (positive sanity)
    // ═════════════════════════════════════════════════════════════════════════

    // ── POSITIVE sanity: a clear alt-allele site IS called, with hand-computed
    //    position / alleles / classification (docs §7.1 worked example). ───────
    [Test]
    public void CallVariants_ClearSnp_HandComputedCall()
    {
        // Docs §7.1: ref ATGC vs query ATTC ⇒ one SNP G→T at 0-based position 2,
        // classified Transversion (purine G ↔ pyrimidine T).
        var variants = VariantCaller.CallVariants(new DnaSequence("ATGC"), new DnaSequence("ATTC")).ToList();

        variants.Should().HaveCount(1, "exactly one base differs ⇒ exactly one SNP");
        var v = variants[0];
        v.Type.Should().Be(VariantType.SNP);
        v.Position.Should().Be(2);
        v.ReferenceAllele.Should().Be("G");
        v.AlternateAllele.Should().Be("T");
        VariantCaller.ClassifyMutation(v).Should().Be(MutationType.Transversion, "G↔T is purine↔pyrimidine (INV-05)");
        AssertWellFormedVariant(v, refLen: 4);
    }

    [Test]
    public void CallVariantsFromAlignment_ClearIndels_HandComputedCalls()
    {
        // Insertion: ref "A-C" vs query "ATC" ⇒ Insertion of T at refPos 1.
        var ins = VariantCaller.CallVariantsFromAlignment("A-C", "ATC").ToList();
        ins.Should().ContainSingle().Which.Type.Should().Be(VariantType.Insertion);
        ins[0].AlternateAllele.Should().Be("T");

        // Deletion: ref "ATC" vs query "A-C" ⇒ Deletion of T at refPos 1.
        var del = VariantCaller.CallVariantsFromAlignment("ATC", "A-C").ToList();
        del.Should().ContainSingle().Which.Type.Should().Be(VariantType.Deletion);
        del[0].ReferenceAllele.Should().Be("T");
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-CALL-001 — BE: "zero depth" (zero-information / empty sites)
    // ═════════════════════════════════════════════════════════════════════════

    // Empty aligned strings ⇒ empty result (no call), the alignment analogue of a
    // site with zero reads. §3.3, §6.1.
    [Test]
    public void CallVariantsFromAlignment_BothEmpty_YieldsNoCall()
    {
        VariantCaller.CallVariantsFromAlignment("", "").Should().BeEmpty("empty aligned input ⇒ empty result (§6.1)");
    }

    [TestCase("", "ACGT")]
    [TestCase("ACGT", "")]
    [TestCase(null, "ACGT")]
    [TestCase("ACGT", null)]
    public void CallVariantsFromAlignment_OneSideEmptyOrNull_YieldsNoCall(string? aRef, string? aQuery)
    {
        // string.IsNullOrEmpty(...) ⇒ yield break: no call, no NullReference (§6.1).
        VariantCaller.CallVariantsFromAlignment(aRef!, aQuery!).Should().BeEmpty();
    }

    // Empty DnaSequences ⇒ empty variant set, and — critically — the only division
    // in the unit (VariantDensity = variants/refLen) must be GUARDED at refLen == 0
    // to return 0, NOT throw DivideByZero nor produce NaN. §5; CalculateStatistics.
    [Test]
    public void CalculateStatistics_EmptyReference_NoDivideByZero_DensityIsZero()
    {
        var stats = VariantCaller.CalculateStatistics(new DnaSequence(""), new DnaSequence(""));

        stats.TotalVariants.Should().Be(0, "no positions ⇒ no variants");
        stats.ReferenceLength.Should().Be(0);
        double.IsNaN(stats.VariantDensity).Should().BeFalse("density must never be NaN at zero length");
        double.IsInfinity(stats.VariantDensity).Should().BeFalse("density must never be ±∞ at zero length");
        stats.VariantDensity.Should().Be(0.0, "documented zero-length guard returns 0 (no DivideByZero)");
        double.IsNaN(stats.TiTvRatio).Should().BeFalse("Ti/Tv must never be NaN");
    }

    // Ti/Tv with no transversions (#Tv == 0) must return the documented 0 sentinel,
    // NOT +∞ or NaN. INV-06, §5.4.
    [Test]
    public void CalculateTiTvRatio_ZeroTransversions_ReturnsZeroSentinel_NotNaNOrInfinity()
    {
        // A→G is a transition; with no transversions the denominator is 0.
        var onlyTransitions = new[]
        {
            new Variant(0, "A", "G", VariantType.SNP, 0),
            new Variant(1, "C", "T", VariantType.SNP, 1),
        };

        double ratio = VariantCaller.CalculateTiTvRatio(onlyTransitions);

        double.IsNaN(ratio).Should().BeFalse("undefined ratio is mapped to a sentinel, not NaN (INV-06)");
        double.IsInfinity(ratio).Should().BeFalse("undefined ratio is mapped to a sentinel, not +∞ (INV-06)");
        ratio.Should().Be(0.0, "#Tv == 0 ⇒ documented 0 sentinel (§5.4)");
    }

    [Test]
    public void CalculateTiTvRatio_EmptyVariantSet_ReturnsZero()
    {
        VariantCaller.CalculateTiTvRatio(Array.Empty<Variant>()).Should().Be(0.0);
    }

    // Fuzz: many empty / single-column / no-difference alignments never crash and
    // never over-call. BE sweep around the empty/zero boundary.
    [Test]
    [CancelAfter(20_000)]
    public void CallVariantsFromAlignment_DegenerateColumns_NeverCrash([Random(1, 1_000_000, 25)] int seed)
    {
        var rng = new Random(seed);
        const string bases = "ACGT-";

        for (int t = 0; t < 40; t++)
        {
            int len = rng.Next(0, 6); // includes 0 (empty) and 1 (single column)
            var refChars = new char[len];
            var qryChars = new char[len];
            for (int i = 0; i < len; i++)
            {
                refChars[i] = bases[rng.Next(bases.Length)];
                qryChars[i] = bases[rng.Next(bases.Length)];
            }

            var aRef = new string(refChars);
            var aQuery = new string(qryChars);

            var act = () => VariantCaller.CallVariantsFromAlignment(aRef, aQuery).ToList();
            act.Should().NotThrow("equal-length aligned strings are always valid input");

            var variants = act();
            AssertAllWellFormed(variants, refLen: aRef.Count(c => c != '-'));
        }
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-CALL-001 — BE: "tumor=normal" (identical channels ⇒ no call)
    // ═════════════════════════════════════════════════════════════════════════

    // The alignment analogue of identical tumor and normal evidence: identical
    // reference == query MUST yield zero variants — no false-positive call.
    // INV-01, §2.4, §6.1.
    [Test]
    public void CallVariants_IdenticalReferenceAndQuery_NoCall()
    {
        var seq = new DnaSequence("ATGCATGCATGC");
        VariantCaller.CallVariants(seq, seq).Should().BeEmpty("identical sequences differ nowhere (INV-01)");
    }

    [Test]
    public void CallVariantsFromAlignment_IdenticalGappedStrings_NoCall()
    {
        VariantCaller.CallVariantsFromAlignment("AT-GC", "AT-GC")
            .Should().BeEmpty("identical aligned strings (incl. shared gap columns) differ nowhere (INV-01)");
    }

    [Test]
    public void CalculateStatistics_IdenticalSequences_ZeroEverything()
    {
        var seq = new DnaSequence("GATTACAGATTACA");
        var stats = VariantCaller.CalculateStatistics(seq, seq);

        stats.TotalVariants.Should().Be(0);
        stats.Snps.Should().Be(0);
        stats.Insertions.Should().Be(0);
        stats.Deletions.Should().Be(0);
        stats.VariantDensity.Should().Be(0.0, "no variants ⇒ density 0");
        stats.TiTvRatio.Should().Be(0.0, "no SNPs ⇒ Ti/Tv 0");
    }

    // Fuzz: ANY randomly generated sequence compared against ITSELF yields no
    // variants — the "identical channels" invariant must hold universally.
    [Test]
    [CancelAfter(20_000)]
    public void CallVariants_SequenceAgainstItself_NeverCalls([Random(1, 1_000_000, 25)] int seed)
    {
        var rng = new Random(seed);
        const string bases = "ACGT";

        for (int t = 0; t < 30; t++)
        {
            int len = rng.Next(1, 40);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = bases[rng.Next(bases.Length)];

            var seq = new DnaSequence(new string(chars));
            VariantCaller.CallVariants(seq, seq).Should().BeEmpty(
                "a sequence compared against itself has zero differences (INV-01)");
        }
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-CALL-001 — BE: "all-N" (non-informative alphabet)
    // ═════════════════════════════════════════════════════════════════════════

    // The typed entry points validate the alphabet (A/C/G/T only); 'N' is a
    // DOCUMENTED ArgumentException, not a silent miscall or crash.
    [Test]
    public void DnaSequence_AllN_IsRejected_DocumentedArgumentException()
    {
        var act = () => new DnaSequence("NNNNN");
        act.Should().Throw<ArgumentException>("DnaSequence accepts only A/C/G/T; 'N' is invalid input");
    }

    // The raw-string entry point accepts arbitrary chars. All-N vs all-N ⇒ NO
    // variants (every column matches N==N) and no crash / NaN.
    [Test]
    public void CallVariantsFromAlignment_AllNVsAllN_NoCall()
    {
        VariantCaller.CallVariantsFromAlignment("NNNNNN", "NNNNNN")
            .Should().BeEmpty("N == N is a match column ⇒ no variant, no crash");
    }

    // Ref (informative) vs all-N query ⇒ every column is a mismatch ⇒ all SNPs;
    // they must be well-formed and classify without crash / NaN even though N is
    // neither purine nor pyrimidine.
    [Test]
    public void CallVariantsFromAlignment_RefVsAllN_AllSnps_ClassifyWithoutCrash()
    {
        const string aRef = "ACGT";
        const string aQuery = "NNNN";

        var variants = VariantCaller.CallVariantsFromAlignment(aRef, aQuery).ToList();

        variants.Should().HaveCount(4, "each base differs from N ⇒ one SNP per column");
        variants.Should().OnlyContain(v => v.Type == VariantType.SNP);
        AssertAllWellFormed(variants, refLen: 4);

        foreach (var v in variants)
        {
            // N is non-purine ⇒ treated as the pyrimidine branch; A/G (purine) ⇒
            // Transversion, C/T (non-purine) ⇒ Transition. Must NOT crash / NaN.
            var act = () => VariantCaller.ClassifyMutation(v);
            act.Should().NotThrow("classification reads only the first base; N must not crash it");
            act().Should().BeOneOf(MutationType.Transition, MutationType.Transversion);
        }

        double ratio = VariantCaller.CalculateTiTvRatio(variants);
        double.IsNaN(ratio).Should().BeFalse("Ti/Tv over N-bearing SNPs must never be NaN");
        double.IsInfinity(ratio).Should().BeFalse("Ti/Tv over N-bearing SNPs must never be ±∞");
    }

    // Fuzz: aligned strings drawn from an N-heavy / arbitrary alphabet (incl. gap,
    // lowercase, digits) never crash the column scan, the classifier, or Ti/Tv.
    [Test]
    [CancelAfter(20_000)]
    public void CallVariantsFromAlignment_NoisyAlphabet_NeverCrash([Random(1, 1_000_000, 25)] int seed)
    {
        var rng = new Random(seed);
        const string alphabet = "ACGTNNNn-acgt0 *"; // N-heavy, gaps, lowercase, digit, space, junk

        for (int t = 0; t < 40; t++)
        {
            int len = rng.Next(0, 12);
            var refChars = new char[len];
            var qryChars = new char[len];
            for (int i = 0; i < len; i++)
            {
                refChars[i] = alphabet[rng.Next(alphabet.Length)];
                qryChars[i] = alphabet[rng.Next(alphabet.Length)];
            }

            var aRef = new string(refChars);
            var aQuery = new string(qryChars);

            List<Variant> variants = null!;
            var call = () => variants = VariantCaller.CallVariantsFromAlignment(aRef, aQuery).ToList();
            call.Should().NotThrow("equal-length arbitrary-content aligned strings must never crash the scan");

            AssertAllWellFormed(variants, refLen: aRef.Count(c => c != '-'));

            foreach (var v in variants.Where(v => v.Type == VariantType.SNP))
            {
                var classify = () => VariantCaller.ClassifyMutation(v);
                classify.Should().NotThrow("classification must never crash on a non-DNA SNP base");
            }

            double ratio = VariantCaller.CalculateTiTvRatio(variants);
            double.IsNaN(ratio).Should().BeFalse("Ti/Tv must never be NaN");
            double.IsInfinity(ratio).Should().BeFalse("Ti/Tv must never be ±∞");
        }
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-CALL-001 — BE: documented input-validation boundaries
    // ═════════════════════════════════════════════════════════════════════════

    [Test]
    public void CallVariants_NullReferenceOrQuery_ThrowsArgumentNullException()
    {
        var nullRef = () => VariantCaller.CallVariants(null!, new DnaSequence("ACGT")).ToList();
        var nullQry = () => VariantCaller.CallVariants(new DnaSequence("ACGT"), null!).ToList();
        nullRef.Should().Throw<ArgumentNullException>("null reference is invalid input (§3.3)");
        nullQry.Should().Throw<ArgumentNullException>("null query is invalid input (§3.3)");
    }

    [Test]
    public void CallVariantsFromAlignment_UnequalAlignedLengths_ThrowsArgumentException()
    {
        var act = () => VariantCaller.CallVariantsFromAlignment("ACGT", "ACG").ToList();
        act.Should().Throw<ArgumentException>("aligned columns must align ⇒ equal lengths required (§6.1)");
    }

    [Test]
    public void CalculateTiTvRatio_Null_ThrowsArgumentNullException()
    {
        var act = () => VariantCaller.CalculateTiTvRatio(null!);
        act.Should().Throw<ArgumentNullException>();
    }

    #endregion
}
