using System;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for MIRNA-CLASSIFY-001 — the opt-in <b>trained</b> natural-vs-background
/// pre-miRNA hairpin classifier (<see cref="MiRnaAnalyzer.ClassifyPreMiRna"/>): "is this
/// sequence a genuine (natural) pre-miRNA hairpin, or background/shuffled?".
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and asserts
/// the code NEVER fails in an undisciplined way: no hang, no unhandled runtime exception
/// (IndexOutOfRange from stem/loop indexing of a degenerate dot-bracket, DivideByZero or
/// NaN from MFEI = AMFE/(G+C)% when GC = 0, an Exp overflow in the logistic link), and no
/// out-of-contract output (a feature outside its domain, or a classification with a NaN
/// probability). Every input must resolve to EITHER a well-defined, theory-correct result
/// (including the null/"not classifiable" outcome) OR a documented validation exception. — §8.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MIRNA-CLASSIFY-001 — pre-miRNA hairpin classifier
/// Checklist: docs/checklists/03_FUZZING.md, row 254.
/// Strategies (per §Description): MC = Malformed Content, BE = Boundary Exploitation.
/// Fuzz targets for row 254: "non-hairpin seq, empty seq, non-ACGU, extreme GC, very short".
/// Source doc: docs/algorithms/MiRNA/Pre_miRNA_Detection.md (§2.2.1 features + §5.3 the
///   logistic-regression classifier + §7.1 the hsa-mir-21 worked example).
/// Source: src/.../Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs
///   • ClassifyPreMiRna(string, double threshold = 0.5, int minLoopSize = 3)   (~line 2463)
///   • ExtractPreMiRnaFeatures(string, int minLoopSize = 3)                      (~line 2421)
///   • ScorePreMiRnaFeatures(PreMiRnaFeatures)                                   (~line 2480)
///   • CalculateMfeIndex(double freeEnergy, int length, double gcPercent)        (~line 1971)
///   • record struct PreMiRnaFeatures(FreeEnergy,Amfe,Mfei,GcContent,
///       PairedFraction,StemBasePairs,LoopSize,Length)                          (~line 2358)
///   • record struct PreMiRnaClassification(Features,NaturalProbability,IsNatural) (~line 2377)
///
/// ───────────────────────────────────────────────────────────────────────────
/// The classifier contract under test (independently derived from the doc, NOT read off code)
/// ───────────────────────────────────────────────────────────────────────────
/// ClassifyPreMiRna folds the candidate ONCE with the validated Zuker–Stiegler MFE engine
/// (RNA-STRUCT-001, Turner 2004) and reads published structure/sequence features from the
/// REAL MFE dot-bracket + base composition (Pre_miRNA_Detection.md §2.2.1, §5.3):
///   • FreeEnergy = ΔG° (kcal/mol) of the MFE structure — ≤ 0 (Bonnet 2004) [6].
///   • AMFE  = 100·|ΔG°|/length                       (Zhang 2006) [7].
///   • MFEI  = AMFE / (G+C)%                            (Zhang 2006) [7]  ← the GC=0 hazard.
///   • GcContent     ∈ [0,1] (a FRACTION).
///   • PairedFraction = 2·#pairs / length ∈ [0,1].
///   • StemBasePairs, LoopSize, Length ≥ 0.
/// It then standardises [FreeEnergy, Amfe, Mfei, GcContent, PairedFraction] by the bundled
/// training mean/std and applies P(natural) = sigmoid(b0 + Σ b_j·z_j) — a logistic-regression
/// model trained on public-domain miRBase positives vs Altschul-Erickson (1985) di-shuffled
/// negatives [13][14]. IsNatural ⇔ P ≥ threshold (default 0.5, the Bayes cutoff).
///
/// Pinned invariants (theory, re-derived here):
///   INV-CL-1  (NULLABLE CONTRACT): ClassifyPreMiRna / ExtractPreMiRnaFeatures return null for
///             null/empty input (and for a fold that degenerates to length 0) — never throw,
///             never an IndexOutOfRange off an empty dot-bracket. The empty-sequence boundary.
///   INV-CL-2  (FINITE FEATURES, all domains): for ANY non-empty input — non-hairpin, homopolymer,
///             extreme-GC, very short, non-ACGU junk — EVERY reported feature is finite (no NaN,
///             no ±Inf): FreeEnergy ≤ 0, GcContent ∈ [0,1], PairedFraction ∈ [0,1], MFEI / AMFE
///             finite and ≥ 0, StemBasePairs/LoopSize/Length ≥ 0, span widths consistent.
///   INV-CL-3  (NO-NaN ON EXTREME GC — the headline hazard): MFEI = AMFE/(G+C)%. A homopolymer
///             with no G/C (all-A, all-U) has (G+C)% = 0; a naive MFEI would be 0/0 = NaN or
///             DivByZero. CalculateMfeIndex guards (length≤0 || gcPercent≤0) → 0, so MFEI is a
///             finite 0, GcContent = 0, and the classification probability is finite. All-G /
///             all-C (GC = 100%) is the opposite extreme — finite, no overflow.
///   INV-CL-4  (PROBABILITY ∈ [0,1], no overflow escape): P(natural) = sigmoid(z) ∈ (0,1) for any
///             finite z; the std==0 guard in ScorePreMiRnaFeatures prevents a /0, so P is finite
///             and the IsNatural call is P ≥ threshold — never NaN, never escaping [0,1].
///   INV-CL-5  (GENUINE HAIRPIN POSITIVE vs HOMOPOLYMER/SHUFFLE NEGATIVE — the discriminative
///             pin, NOT a no-op): the real hsa-mir-21 precursor (Pre_miRNA_Detection.md §7.1)
///             classifies IsNatural == true with high P; an all-base homopolymer (no structure)
///             and a dinucleotide-shuffled (composition-matched) version classify IsNatural ==
///             false. A model that called everything natural — or a test that passed against one —
///             would be invalid.
///
/// If a source-derived assertion and the code disagree, the CODE is wrong (fixed minimally per
/// the doc): a NaN MFEI on a GC=0 homopolymer, an IndexOutOfRange from degenerate folding, a NaN
/// probability, or a positive classification of a homopolymer would each be a REAL bug.
///
/// LimitationPolicy: these are pure static methods; the assembly module-initializer
/// (_LimitationPolicyTestBootstrap) already runs under LimitationMode.Permissive (see
/// MiRnaPctFuzzTests). RNA folding is O(n^3) so the folding tests carry [CancelAfter].
///
/// All inputs are fixed / deterministically generated; the random helper uses a LOCALLY seeded
/// `new Random(seed)` (no shared static Rng), so every fuzz input is reproducible.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MiRnaClassifyFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomRna(int length, int seed)
    {
        const string bases = "ACGU";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal feature-vector invariants (INV-CL-2/3): every reported feature is
    /// finite and inside its theoretical domain — the core "no out-of-contract output" bar.
    /// </summary>
    private static void AssertFeaturesWellFormed(MiRnaAnalyzer.PreMiRnaFeatures f)
    {
        double.IsNaN(f.FreeEnergy).Should().BeFalse("MFE (ΔG°) is finite, never NaN (INV-CL-2)");
        double.IsInfinity(f.FreeEnergy).Should().BeFalse("MFE (ΔG°) is finite, never ±Inf (INV-CL-2)");
        f.FreeEnergy.Should().BeLessThanOrEqualTo(1e-9, "the MFE of an optimal structure is ≤ 0 (Bonnet 2004)");

        double.IsNaN(f.Amfe).Should().BeFalse("AMFE = 100·|ΔG°|/length is finite (INV-CL-2)");
        double.IsInfinity(f.Amfe).Should().BeFalse("AMFE is finite (INV-CL-2)");
        f.Amfe.Should().BeGreaterThanOrEqualTo(0.0, "AMFE uses |ΔG°| ⇒ ≥ 0 (Zhang 2006)");

        // THE headline hazard: MFEI = AMFE/(G+C)% must never NaN/Inf on extreme GC (GC=0).
        double.IsNaN(f.Mfei).Should().BeFalse("MFEI = AMFE/(G+C)% must never be NaN on extreme GC (INV-CL-3)");
        double.IsInfinity(f.Mfei).Should().BeFalse("MFEI must never be ±Inf on extreme GC (INV-CL-3)");
        f.Mfei.Should().BeGreaterThanOrEqualTo(0.0, "MFEI is a non-negative ratio of non-negative quantities");

        f.GcContent.Should().BeInRange(0.0, 1.0, "G+C content is a fraction in [0,1] (INV-CL-2)");
        f.PairedFraction.Should().BeInRange(0.0, 1.0, "paired fraction = 2·#pairs/length ∈ [0,1] (INV-CL-2)");

        f.StemBasePairs.Should().BeGreaterThanOrEqualTo(0, "base-pair count is ≥ 0 (INV-CL-2)");
        f.LoopSize.Should().BeGreaterThanOrEqualTo(0, "loop size is ≥ 0 (INV-CL-2)");
        f.Length.Should().BeGreaterThan(0, "a classifiable (non-empty) candidate has positive length (INV-CL-2)");
        (2 * f.StemBasePairs).Should().BeLessThanOrEqualTo(f.Length,
            "paired bases (2·#pairs) cannot exceed the sequence length (INV-CL-2)");
    }

    /// <summary>Asserts the classification probability is a finite member of [0,1] (INV-CL-4).</summary>
    private static void AssertProbabilityWellFormed(MiRnaAnalyzer.PreMiRnaClassification c)
    {
        double.IsNaN(c.NaturalProbability).Should().BeFalse("P(natural) is finite, never NaN (INV-CL-4)");
        double.IsInfinity(c.NaturalProbability).Should().BeFalse("P(natural) is finite, never ±Inf (INV-CL-4)");
        c.NaturalProbability.Should().BeInRange(0.0, 1.0, "P(natural) = sigmoid(z) ∈ [0,1] (INV-CL-4)");
        c.IsNatural.Should().Be(c.NaturalProbability >= 0.5,
            "IsNatural ⇔ P(natural) ≥ the default 0.5 threshold (INV-CL-4)");
    }

    // hsa-mir-21 (MI0000077) precursor — the doc's §7.1 worked example natural positive.
    private const string HsaMir21 =
        "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA";

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  MIRNA-CLASSIFY-001 — pre-miRNA hairpin classifier : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region MIRNA-CLASSIFY-001 — BE: empty / null sequence (the nullable-contract boundary)

    /// <summary>
    /// BE — "empty seq": the empty and null sequence. ExtractPreMiRnaFeatures /
    /// ClassifyPreMiRna short-circuit IsNullOrEmpty → null BEFORE folding, so they never
    /// touch the MFE engine or index an empty dot-bracket (INV-CL-1, §3.3). The result is a
    /// well-defined null ("not classifiable"), never an exception.
    /// </summary>
    [Test]
    public void Classify_EmptyOrNull_ReturnsNullNeverThrows()
    {
        var actEmpty = () => ClassifyPreMiRna("");
        actEmpty.Should().NotThrow("an empty sequence is a documented null result, not an error (INV-CL-1)")
                .Subject.Should().BeNull("an empty sequence has no foldable structure to classify");

        var actNull = () => ClassifyPreMiRna(null!);
        actNull.Should().NotThrow("a null sequence is treated as empty, not an error (INV-CL-1)")
               .Subject.Should().BeNull("a null sequence yields no classification");

        // The feature-extraction surface shares the same nullable contract.
        ExtractPreMiRnaFeatures("").Should().BeNull("empty → no features (INV-CL-1)");
        ((object?)ExtractPreMiRnaFeatures(null!)).Should().BeNull("null → no features (INV-CL-1)");
    }

    #endregion

    #region MIRNA-CLASSIFY-001 — BE: very short sequence (below any hairpin)

    /// <summary>
    /// BE — "very short": a sequence too short to fold a hairpin. A 1..12-nt input either folds
    /// to a trivial all-unpaired structure (ΔG° = 0, 0 pairs, MFEI 0) or, for length 0 after
    /// normalisation, returns null — but it must NEVER crash on stem/loop indexing of a tiny
    /// dot-bracket, NEVER produce a NaN MFEI, and NEVER be called a natural pre-miRNA (no stem
    /// ⇒ background). We sweep every short length and assert the well-formed, negative outcome.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void Classify_VeryShort_FiniteFeatures_NotNatural()
    {
        foreach (int len in new[] { 1, 2, 3, 4, 5, 8, 10, 12 })
        {
            string tiny = new string('A', len);
            var act = () => ClassifyPreMiRna(tiny);
            var result = act.Should().NotThrow(
                $"a {len}-nt sequence must not crash the classifier (no IndexOutOfRange off a tiny fold)").Subject;

            // All-A folds to no pairs; the call is non-null with finite, in-domain features.
            result.Should().NotBeNull($"a {len}-nt non-empty sequence still folds to a (trivial) structure");
            var c = result!.Value;
            AssertFeaturesWellFormed(c.Features);
            AssertProbabilityWellFormed(c);
            c.Features.StemBasePairs.Should().Be(0, "a short homopolymer has no base pairs ⇒ no stem");
            c.IsNatural.Should().BeFalse("a tiny structureless sequence is NOT a natural pre-miRNA (INV-CL-5)");
        }
    }

    #endregion

    #region MIRNA-CLASSIFY-001 — BE: non-hairpin sequence (no dominant stem-loop)

    /// <summary>
    /// BE — "non-hairpin seq": a full-length sequence that does NOT fold into a genuine pre-miRNA
    /// hairpin. Random sequences (Bonnet 2004: random RNA folds to markedly LESS stable structures
    /// than real pre-miRNAs) and a deliberately non-pairing block sequence must classify as
    /// background (IsNatural == false) with finite, in-domain features — never a crash, never a NaN,
    /// never a spurious positive. We sweep fixed-seed random RNAs across pre-miRNA-scale lengths.
    /// </summary>
    [Test]
    [CancelAfter(120_000)]
    public void Classify_NonHairpin_FiniteFeatures_BackgroundCall()
    {
        // Hand-built non-hairpin: a block A-run then C-run cannot fold a dominant hairpin stem.
        string block = new string('A', 35) + new string('C', 35); // 70 nt, low structure
        var blockResult = ClassifyPreMiRna(block);
        blockResult.Should().NotBeNull("a 70-nt sequence folds to some (possibly trivial) structure");
        AssertFeaturesWellFormed(blockResult!.Value.Features);
        AssertProbabilityWellFormed(blockResult.Value);

        foreach (int seed in new[] { 1, 7, 42, 2026 })
        {
            foreach (int len in new[] { 60, 80 })
            {
                string raw = RandomRna(len, seed);
                var act = () => ClassifyPreMiRna(raw);
                var result = act.Should().NotThrow(
                    $"a random {len}-nt sequence must not crash the classifier (seed {seed})").Subject;
                result.Should().NotBeNull();
                var c = result!.Value;
                AssertFeaturesWellFormed(c.Features);
                AssertProbabilityWellFormed(c);
                // Random sequences are not natural pre-miRNAs (Bonnet 2004); the trained model
                // calls them background. (Asserted as the dominant outcome of the discriminative pin.)
                c.IsNatural.Should().BeFalse(
                    $"a random non-pre-miRNA sequence classifies as background (seed {seed}, len {len}; INV-CL-5)");
            }
        }
    }

    #endregion

    #region MIRNA-CLASSIFY-001 — BE/MC: extreme GC (the MFEI = AMFE/(G+C)% no-NaN hazard)

    /// <summary>
    /// BE+MC — "extreme GC": the headline numerical hazard. MFEI = AMFE/(G+C)%, so a homopolymer
    /// with NO G or C (all-A, all-U) drives the denominator (G+C)% to 0 — a naive MFEI would be
    /// 0/0 = NaN or a DivideByZero. CalculateMfeIndex guards (gcPercent ≤ 0) → 0, so MFEI is a
    /// finite 0 and GcContent is 0. The opposite extreme (all-G, all-C ⇒ GC = 100%) must likewise
    /// be finite, with no overflow. None of these self-pairing-free homopolymers is a natural
    /// pre-miRNA. We assert finite features, a finite probability, and the background call at both
    /// GC extremes — directly exercising the GC=0 denominator edge (INV-CL-3).
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void Classify_ExtremeGc_NoNaNMfei_FiniteProbability()
    {
        // Direct unit pin on the MFEI helper: GC% = 0 ⇒ guarded to 0, never 0/0 = NaN.
        double mfeiGcZero = CalculateMfeIndex(-10.0, 60, 0.0);
        double.IsNaN(mfeiGcZero).Should().BeFalse("MFEI(GC%=0) must be the guarded 0, not 0/0 = NaN (INV-CL-3)");
        mfeiGcZero.Should().Be(0.0, "the (G+C)%=0 denominator edge is guarded to 0, not a division");
        // GC% = 100 is the opposite extreme — a normal finite division.
        double mfeiGcFull = CalculateMfeIndex(-30.0, 60, 100.0);
        double.IsNaN(mfeiGcFull).Should().BeFalse("MFEI(GC%=100) is a finite division (INV-CL-3)");
        mfeiGcFull.Should().BeApproximately((100.0 * 30.0 / 60.0) / 100.0, 1e-12,
            "MFEI = (100·|ΔG°|/length)/(G+C)% — hand-traced at GC% = 100");
        // Length = 0 edge is also guarded (no DivByZero on the AMFE 100·|ΔG°|/length term).
        CalculateMfeIndex(-5.0, 0, 50.0).Should().Be(0.0, "length 0 is guarded to MFEI 0, not /0");

        // End-to-end: GC-extreme homopolymers through the full classifier.
        foreach (char b in new[] { 'A', 'U', 'G', 'C' })
        {
            string homo = new string(b, 70);
            var act = () => ClassifyPreMiRna(homo);
            var result = act.Should().NotThrow(
                $"an all-'{b}' homopolymer (extreme GC) must not crash the classifier").Subject;
            result.Should().NotBeNull();
            var c = result!.Value;

            AssertFeaturesWellFormed(c.Features);          // includes the no-NaN-MFEI assertion
            AssertProbabilityWellFormed(c);

            // All-A / all-U have GC = 0; all-G / all-C have GC = 1 — both extremes, both finite.
            if (b == 'A' || b == 'U')
            {
                c.Features.GcContent.Should().Be(0.0, $"an all-'{b}' homopolymer has zero G+C content (GC=0 edge)");
                c.Features.Mfei.Should().Be(0.0, "MFEI is the guarded 0 when (G+C)% = 0 (INV-CL-3)");
            }
            else
            {
                c.Features.GcContent.Should().Be(1.0, $"an all-'{b}' homopolymer is 100% G+C");
            }

            // No homopolymer self-pairs into a hairpin ⇒ background, never a spurious positive.
            c.IsNatural.Should().BeFalse($"an all-'{b}' homopolymer is NOT a natural pre-miRNA (INV-CL-5)");
        }
    }

    #endregion

    #region MIRNA-CLASSIFY-001 — MC: non-ACGU characters

    /// <summary>
    /// MC — "non-ACGU": DNA 'T', 'N' ambiguity codes, digits, punctuation, whitespace, mixed case.
    /// Per the doc the input is normalised (uppercase, T→U) before folding; characters outside the
    /// pairing alphabet simply never pair, so junk content fails to fold a hairpin rather than
    /// crashing the classifier. We verify: (a) a DNA-spelled real precursor (T→U) classifies the
    /// SAME as its RNA spelling (T→U equivalence), and (b) arbitrary 'N'/digit/punctuation junk is
    /// folded and classified with finite features and a finite probability — never a throw, never a
    /// NaN, never a malformed feature.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void Classify_NonAcguCharacters_NormalizedOrRejected_NeverCrash()
    {
        // (a) DNA spelling of hsa-mir-21 (+ lowercase) folds identically to the RNA spelling.
        string dnaMir21 = HsaMir21.Replace('U', 'T').ToLowerInvariant();
        var rna = ClassifyPreMiRna(HsaMir21);
        var dna = ClassifyPreMiRna(dnaMir21);
        rna.Should().NotBeNull();
        dna.Should().NotBeNull("a DNA-spelled precursor folds via T→U normalisation, not a crash");
        dna!.Value.NaturalProbability.Should().BeApproximately(rna!.Value.NaturalProbability, 1e-9,
            "T→U + uppercasing make the DNA spelling equivalent to the RNA spelling");
        dna.Value.IsNatural.Should().Be(rna.Value.IsNatural,
            "the DNA spelling yields the same natural/background call as the RNA spelling");

        // (b) Arbitrary junk content: 'N', digits, punctuation, whitespace — never pair, never crash.
        foreach (string junk in new[]
                 {
                     new string('N', 70),                                  // ambiguity homopolymer
                     new string('N', 30) + "12345!@#$%^&*() \t" + new string('N', 30),
                     "GCGC" + new string('N', 50) + "GCGC",                 // structure-breaking N core
                 })
        {
            var act = () => ClassifyPreMiRna(junk);
            var result = act.Should().NotThrow(
                $"junk content (len {junk.Length}) must not crash the classifier").Subject;

            if (result is not null)
            {
                var c = result.Value;
                AssertFeaturesWellFormed(c.Features);
                AssertProbabilityWellFormed(c);
            }
        }
    }

    #endregion

    #region MIRNA-CLASSIFY-001 — discriminative pin: genuine hairpin POSITIVE vs shuffle NEGATIVE

    /// <summary>
    /// The discriminative pin (INV-CL-5): the classifier is NOT a no-op that calls everything
    /// natural (or everything background). The real hsa-mir-21 (MI0000077) precursor — the doc's
    /// §7.1 worked example — folds to a markedly stable hairpin (Bonnet 2004) and MUST classify
    /// IsNatural == true with high P(natural); a dinucleotide-shuffled, COMPOSITION-MATCHED version
    /// (Altschul-Erickson 1985 background convention) destroys the structure and MUST classify
    /// IsNatural == false. A test that passed against an everything-positive model would be invalid;
    /// this asserts both arms of the discrimination on independently chosen inputs.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void Classify_GenuineHairpinPositive_vs_ShuffleBackgroundNegative()
    {
        // POSITIVE: the genuine miRBase precursor.
        var natural = ClassifyPreMiRna(HsaMir21);
        natural.Should().NotBeNull("the genuine hsa-mir-21 precursor is classifiable");
        var nc = natural!.Value;
        AssertFeaturesWellFormed(nc.Features);
        AssertProbabilityWellFormed(nc);
        nc.IsNatural.Should().BeTrue("the real hsa-mir-21 precursor classifies as a natural pre-miRNA (INV-CL-5)");
        nc.NaturalProbability.Should().BeGreaterThan(0.5,
            "a genuine precursor scores above the 0.5 Bayes cutoff");
        // It really folds to a stable, well-paired hairpin (the discriminative structural signal).
        nc.Features.FreeEnergy.Should().BeLessThan(0.0, "a real precursor folds to a stable (ΔG° < 0) hairpin");
        nc.Features.StemBasePairs.Should().BeGreaterThan(0, "a real precursor has a paired stem");
        nc.Features.Mfei.Should().BeGreaterThan(0.0, "a structured precursor has a positive, finite MFEI");

        // NEGATIVE: a composition-matched dinucleotide shuffle scores as background.
        // Several fixed seeds: the structural signal is destroyed regardless of the shuffle draw.
        foreach (int seed in new[] { 999, 7, 20260626 })
        {
            string shuffled = DinucleotideShuffle(HsaMir21, new Random(seed));
            var background = ClassifyPreMiRna(shuffled);
            background.Should().NotBeNull($"the shuffled sequence is still classifiable (seed {seed})");
            var bc = background!.Value;
            AssertFeaturesWellFormed(bc.Features);
            AssertProbabilityWellFormed(bc);
            bc.IsNatural.Should().BeFalse(
                $"a dinucleotide-shuffled, composition-matched sequence is background (seed {seed}; INV-CL-5)");
            bc.NaturalProbability.Should().BeLessThan(nc.NaturalProbability,
                $"the shuffle scores strictly below the genuine precursor (seed {seed})");
        }
    }

    /// <summary>
    /// ScorePreMiRnaFeatures consistency + std-guard: the standalone scorer reproduces the
    /// classifier's probability on the same extracted features (INV-CL-4), and a degenerate feature
    /// vector (all-zero, the std==0 guard path) yields a finite probability in [0,1], never a /0.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void ScorePreMiRnaFeatures_MatchesClassifier_AndHandlesDegenerateVector()
    {
        var features = ExtractPreMiRnaFeatures(HsaMir21);
        features.Should().NotBeNull();
        var call = ClassifyPreMiRna(HsaMir21);
        call.Should().NotBeNull();

        double standalone = ScorePreMiRnaFeatures(features!.Value);
        standalone.Should().BeApproximately(call!.Value.NaturalProbability, 1e-12,
            "ScorePreMiRnaFeatures reproduces the classifier probability on the same features (INV-CL-4)");
        standalone.Should().BeInRange(0.0, 1.0, "the standalone score is a probability in [0,1]");

        // A degenerate all-zero feature vector exercises the std==0 guard path; must stay finite.
        var zero = new MiRnaAnalyzer.PreMiRnaFeatures(0, 0, 0, 0, 0, 0, 0, 0);
        double zeroP = ScorePreMiRnaFeatures(zero);
        double.IsNaN(zeroP).Should().BeFalse("a degenerate feature vector must not produce NaN (INV-CL-4)");
        double.IsInfinity(zeroP).Should().BeFalse("a degenerate feature vector must not produce ±Inf (INV-CL-4)");
        zeroP.Should().BeInRange(0.0, 1.0, "the score stays a probability even for a degenerate vector");
    }

    #endregion
}
