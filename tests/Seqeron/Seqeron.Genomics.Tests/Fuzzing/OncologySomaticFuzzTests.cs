using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology somatic-mutation-calling area — ONCO-SOMATIC-001.
/// The unit under test is the deterministic tumor/normal allele-fraction classifier
/// <see cref="OncologyAnalyzer.CallSomaticMutations"/> (and its single-variant entry
/// point <see cref="OncologyAnalyzer.Classify"/>, the <see cref="OncologyAnalyzer.FilterGermlineVariants"/>
/// projection, and the underlying VAF helper <see cref="OncologyAnalyzer.CalculateVAF"/>),
/// implemented in src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no nonsense output, and no *unhandled* runtime exception
/// (DivideByZero / NullReference / Overflow). Every input must resolve to EITHER
/// a well-defined, theory-correct value OR a *documented, intentional* outcome
/// (here, an <see cref="ArgumentOutOfRangeException"/> for malformed read counts).
/// For somatic allele-fraction calling the headline hazards are:
///   • a DivideByZeroException computing VAF = alt/depth when total depth == 0
///     (the documented contract is VAF = 0 at zero coverage, INV-05 — NOT a throw,
///     NOT a NaN);
///   • a VAF that escapes [0, 1] — in particular alt > depth must NOT yield a
///     VAF > 1 (malformed counts are a documented ArgumentOutOfRangeException);
///   • a negative read count silently producing a negative or out-of-range VAF;
///   • a SomaticScore escaping its documented [0, 1] = max(0, f_t − f_n) range;
///   • identical tumor == normal evidence being mis-called Somatic (it must be
///     Germline or NotDetected — there is no tumor-vs-normal separation);
///   • a crash / NaN on empty read evidence or non-informative all-N alleles.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-SOMATIC-001 — Somatic mutation calling (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 87.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///   • MC = Malformed Content — невалідний контент (alt &gt; depth, all-N bases).
///     Targets (checklist row 87): "0 depth, alt>depth, tumor=normal, empty reads,
///     all-N bases".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Somatic_Mutation_Calling.md (docs/algorithms/Oncology/Somatic_Mutation_Calling.md):
///   • VAF f = altReads / totalReads, in [0, 1]                          (§2.2)
///   • totalReads == 0 ⇒ VAF = 0 (uncovered site = allele absent)        (INV-05, §6.1)
///   • Read counts must satisfy 0 ≤ alt ≤ total, else
///       ArgumentOutOfRangeException                                     (§3.3)
///   • Decision rule:
///       f_t &lt; τ_t                         → NotDetected               (INV-02, §4.2)
///       f_t ≥ τ_t ∧ f_n ≤ τ_n             → Somatic                    (INV-01, §4.2)
///       f_t ≥ τ_t ∧ f_n &gt; τ_n             → Germline                   (§4.2)
///       defaults τ_t = 0.05, τ_n = 0.01                                 (§2.2, §3.1)
///   • SomaticScore = max(0, f_t − f_n) when Somatic, else 0; in [0, 1]  (INV-03, §3.2)
///   • Tumor-only (normal total = 0) ⇒ f_n = 0 ⇒ Somatic if f_t ≥ τ_t    (§6.1)
///   • f_n ≥ f_t ⇒ score = 0                                             (§6.1)
///   • FilterGermlineVariants returns EXACTLY the Somatic subset, in
///       input order                                                     (INV-04, §3.2)
///   • Output preserves input order and count                            (§3.3)
///   • The classifier reads ONLY the integer read counts — allele
///       strings (e.g. all-"N") are carried through untouched and never
///       affect the decision.                                           (record VariantObservation, §3)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologySomaticFuzzTests
{
    private const double TumorTau = DefaultTumorVafThreshold;   // 0.05
    private const double NormalTau = DefaultNormalVafThreshold; // 0.01

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented numeric contract on EVERY call result: both VAFs must
    // be FINITE and inside [0, 1] (no DivideByZero NaN, no alt>depth overshoot
    // leaking through), the SomaticScore must sit inside the documented
    // [0, 1] = max(0, f_t − f_n) range, and the score must be exactly 0 for any
    // non-Somatic status. This is what stops a fuzz test from rubber-stamping a
    // NaN / out-of-range result green.
    private static void AssertWellFormedCall(SomaticCall c)
    {
        double.IsNaN(c.TumorVaf).Should().BeFalse("tumor VAF must never be NaN");
        double.IsNaN(c.NormalVaf).Should().BeFalse("normal VAF must never be NaN");
        double.IsNaN(c.SomaticScore).Should().BeFalse("somatic score must never be NaN");

        double.IsInfinity(c.TumorVaf).Should().BeFalse("tumor VAF must be finite");
        double.IsInfinity(c.NormalVaf).Should().BeFalse("normal VAF must be finite");
        double.IsInfinity(c.SomaticScore).Should().BeFalse("somatic score must be finite");

        c.TumorVaf.Should().BeInRange(0.0, 1.0, "VAF = alt/total ∈ [0, 1] (§2.2)");
        c.NormalVaf.Should().BeInRange(0.0, 1.0, "VAF = alt/total ∈ [0, 1] (§2.2)");

        // INV-03: SomaticScore ∈ [0, 1], = max(0, f_t − f_n).
        c.SomaticScore.Should().BeInRange(0.0, 1.0, "score = max(0, f_t − f_n) ∈ [0, 1]");

        if (c.Status == SomaticStatus.Somatic)
        {
            c.SomaticScore.Should().BeApproximately(Math.Max(0.0, c.TumorVaf - c.NormalVaf), 1e-12);
        }
        else
        {
            c.SomaticScore.Should().Be(0.0, "score is 0 for non-Somatic calls (§3.2)");
        }
    }

    private static VariantObservation Obs(
        int tAlt, int tTotal, int nAlt, int nTotal,
        string refA = "A", string altA = "T") =>
        new("chr1", 100, refA, altA, tAlt, tTotal, nAlt, nTotal);

    #region ONCO-SOMATIC-001 — Somatic mutation calling (positive sanity)

    // ── POSITIVE sanity: a clear somatic case is CALLED, with hand-computed VAF/score ──
    [Test]
    public void Classify_ClearSomatic_HandComputedVafAndScore()
    {
        // Docs §7.1 worked example: f_t = 25/100 = 0.25 ≥ τ_t, f_n = 0/100 = 0 ≤ τ_n
        // → Somatic; score = 0.25 − 0.00 = 0.25.
        var call = Classify(Obs(25, 100, 0, 100));

        AssertWellFormedCall(call);
        call.TumorVaf.Should().BeApproximately(0.25, 1e-12);
        call.NormalVaf.Should().Be(0.0);
        call.Status.Should().Be(SomaticStatus.Somatic);
        call.SomaticScore.Should().BeApproximately(0.25, 1e-12);
    }

    // ── POSITIVE sanity: a germline case (allele present in BOTH) is NOT somatic ──
    [Test]
    public void Classify_ClearGermline_NotCalledSomatic()
    {
        // Docs §7.1 worked example: f_t = 48/100 = 0.48 ≥ τ_t but f_n = 50/100 = 0.50 > τ_n
        // → Germline; score 0.
        var call = Classify(Obs(48, 100, 50, 100));

        AssertWellFormedCall(call);
        call.TumorVaf.Should().BeApproximately(0.48, 1e-12);
        call.NormalVaf.Should().BeApproximately(0.50, 1e-12);
        call.Status.Should().Be(SomaticStatus.Germline);
        call.SomaticScore.Should().Be(0.0);
    }

    [Test]
    public void Classify_SubThresholdTumor_IsNotDetected()
    {
        // f_t = 4/100 = 0.04 < τ_t (0.05) ⇒ NotDetected (INV-02, §4.2), even with
        // a perfectly clean normal.
        var call = Classify(Obs(4, 100, 0, 100));

        AssertWellFormedCall(call);
        call.TumorVaf.Should().BeApproximately(0.04, 1e-12);
        call.Status.Should().Be(SomaticStatus.NotDetected);
    }

    #endregion

    #region ONCO-SOMATIC-001 — BE: zero depth (no DivideByZero, VAF = 0, INV-05)

    [Test]
    public void Classify_TumorDepthZero_VafZero_NotDetected_NoDivideByZero()
    {
        // Tumor depth 0: VAF = 0 (INV-05), so f_t = 0 < τ_t ⇒ NotDetected.
        // The hazard is alt/depth = 0/0 → DivideByZero / NaN; the documented
        // contract is VAF = 0.
        var act = () => Classify(Obs(0, 0, 0, 100));

        act.Should().NotThrow();
        var call = act();
        AssertWellFormedCall(call);
        call.TumorVaf.Should().Be(0.0);
        call.Status.Should().Be(SomaticStatus.NotDetected);
    }

    [Test]
    public void Classify_NormalDepthZero_TumorOnly_IsSomatic()
    {
        // Tumor-only mode (normal total = 0 ⇒ f_n = 0, ℓ_n = 1 analogue): a present
        // tumor allele is Somatic (§6.1). f_t = 30/100 = 0.30 ≥ τ_t, f_n = 0/0 = 0.
        var call = Classify(Obs(30, 100, 0, 0));

        AssertWellFormedCall(call);
        call.NormalVaf.Should().Be(0.0, "normal total 0 ⇒ VAF 0, not NaN (INV-05)");
        call.TumorVaf.Should().BeApproximately(0.30, 1e-12);
        call.Status.Should().Be(SomaticStatus.Somatic);
        call.SomaticScore.Should().BeApproximately(0.30, 1e-12);
    }

    [Test]
    public void Classify_BothDepthsZero_NoThrow_NotDetected()
    {
        // Both samples uncovered: f_t = f_n = 0; f_t < τ_t ⇒ NotDetected, no crash.
        var act = () => Classify(Obs(0, 0, 0, 0));

        act.Should().NotThrow();
        var call = act();
        AssertWellFormedCall(call);
        call.TumorVaf.Should().Be(0.0);
        call.NormalVaf.Should().Be(0.0);
        call.Status.Should().Be(SomaticStatus.NotDetected);
    }

    [Test]
    public void CalculateVAF_ZeroDepth_ReturnsZero_NoDivideByZero()
    {
        // Direct VAF helper at zero coverage (INV-05). alt is also 0 (alt ≤ total).
        var act = () => CalculateVAF(0, 0);

        act.Should().NotThrow();
        act().Should().Be(0.0);
    }

    #endregion

    #region ONCO-SOMATIC-001 — MC: alt > depth (documented reject, no VAF > 1)

    [Test]
    public void CalculateVAF_AltExceedsDepth_ThrowsArgumentOutOfRange_NeverReturnsAboveOne()
    {
        // Malformed counts: alt (50) > total (10). Documented contract (§3.3):
        // ArgumentOutOfRangeException — the result must NEVER be a VAF of 5.0.
        var act = () => CalculateVAF(50, 10);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void Classify_TumorAltExceedsTotal_ThrowsArgumentOutOfRange()
    {
        // alt > total in the tumor must be rejected, not silently yield VAF > 1.
        var act = () => Classify(Obs(101, 100, 0, 100));

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void Classify_NormalAltExceedsTotal_ThrowsArgumentOutOfRange()
    {
        var act = () => Classify(Obs(20, 100, 60, 50));

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void Classify_NegativeAltReads_ThrowsArgumentOutOfRange()
    {
        // Negative counts (-1 boundary, BE) must be rejected, not produce a
        // negative VAF.
        var act = () => Classify(Obs(-1, 100, 0, 100));

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void Classify_NegativeTotalReads_ThrowsArgumentOutOfRange()
    {
        var act = () => Classify(Obs(0, -5, 0, 100));

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void Classify_AltEqualsTotal_VafExactlyOne_NoOvershoot()
    {
        // Boundary alt == total: VAF must be exactly 1.0, never above (BE upper edge).
        var call = Classify(Obs(80, 80, 0, 100));

        AssertWellFormedCall(call);
        call.TumorVaf.Should().Be(1.0);
        call.Status.Should().Be(SomaticStatus.Somatic);
    }

    [Test]
    public void CalculateVAF_RandomMalformedCounts_AlwaysGuardedOrInRange()
    {
        // Fuzz: random non-negative count pairs. When alt ≤ total the VAF must be
        // a finite value in [0, 1]; when alt > total the call must THROW — never
        // a nonsense VAF > 1.
        var rng = new Random(8701);
        for (int i = 0; i < 4000; i++)
        {
            int total = rng.Next(0, 5000);
            int alt = rng.Next(0, 5000);

            if (alt <= total)
            {
                double vaf = CalculateVAF(alt, total);
                double.IsNaN(vaf).Should().BeFalse();
                vaf.Should().BeInRange(0.0, 1.0);
                if (total == 0) vaf.Should().Be(0.0);
            }
            else
            {
                var act = () => CalculateVAF(alt, total);
                act.Should().Throw<ArgumentOutOfRangeException>();
            }
        }
    }

    #endregion

    #region ONCO-SOMATIC-001 — MC: tumor == normal (identical ⇒ never Somatic)

    [Test]
    public void Classify_IdenticalTumorAndNormal_PresentAllele_IsGermlineNotSomatic()
    {
        // Identical evidence in tumor and normal: f_t = f_n = 0.40. With the allele
        // clearly present in BOTH (f_n > τ_n) there is no tumor-vs-normal separation
        // ⇒ Germline, NOT Somatic.
        var call = Classify(Obs(40, 100, 40, 100));

        AssertWellFormedCall(call);
        call.TumorVaf.Should().BeApproximately(0.40, 1e-12);
        call.NormalVaf.Should().BeApproximately(0.40, 1e-12);
        call.Status.Should().Be(SomaticStatus.Germline);
        call.Status.Should().NotBe(SomaticStatus.Somatic);
        call.SomaticScore.Should().Be(0.0, "f_n ≥ f_t ⇒ no separation (§6.1)");
    }

    [Test]
    public void Classify_IdenticalTumorAndNormal_AbsentAllele_IsNotDetected()
    {
        // Identical low evidence: f_t = f_n = 0.02 < τ_t ⇒ NotDetected (LoD gate
        // applied first, INV-02). Still not Somatic.
        var call = Classify(Obs(2, 100, 2, 100));

        AssertWellFormedCall(call);
        call.Status.Should().Be(SomaticStatus.NotDetected);
        call.Status.Should().NotBe(SomaticStatus.Somatic);
    }

    [Test]
    public void Classify_RandomIdenticalTumorEqualsNormal_NeverSomatic()
    {
        // Fuzz: any variant with tumor evidence == normal evidence can never be a
        // somatic call — there is no separation f_t ≠ f_n to exploit.
        var rng = new Random(424242);
        for (int i = 0; i < 3000; i++)
        {
            int total = rng.Next(1, 4000);
            int alt = rng.Next(0, total + 1);

            var call = Classify(Obs(alt, total, alt, total));

            AssertWellFormedCall(call);
            call.TumorVaf.Should().Be(call.NormalVaf);
            call.Status.Should().NotBe(SomaticStatus.Somatic,
                "identical tumor==normal evidence has no somatic separation");
        }
    }

    #endregion

    #region ONCO-SOMATIC-001 — BE: empty reads / empty input

    [Test]
    public void CallSomaticMutations_EmptyInput_ReturnsEmpty_NoThrow()
    {
        // Empty pileup / no variants ⇒ empty result (§6.1), no exception.
        var act = () => CallSomaticMutations(Array.Empty<VariantObservation>());

        act.Should().NotThrow();
        act().Should().BeEmpty();
    }

    [Test]
    public void CallSomaticMutations_NullInput_ThrowsArgumentNull()
    {
        var act = () => CallSomaticMutations(null!);

        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void FilterGermlineVariants_EmptyInput_ReturnsEmpty()
    {
        FilterGermlineVariants(Array.Empty<VariantObservation>()).Should().BeEmpty();
    }

    [Test]
    public void CallSomaticMutations_AllZeroDepthVariants_NoCallsCrash()
    {
        // A batch of "empty read" variants (zero coverage in both samples): every
        // call must be NotDetected with VAF 0, no DivideByZero anywhere.
        var variants = Enumerable.Range(0, 50)
            .Select(_ => Obs(0, 0, 0, 0))
            .ToArray();

        var calls = CallSomaticMutations(variants);

        calls.Should().HaveCount(50, "output preserves input count (§3.3)");
        foreach (var c in calls)
        {
            AssertWellFormedCall(c);
            c.Status.Should().Be(SomaticStatus.NotDetected);
        }
    }

    #endregion

    #region ONCO-SOMATIC-001 — MC: all-N bases (non-informative alleles carried, decision unaffected)

    [Test]
    public void Classify_AllNBases_DecisionDependsOnlyOnReadCounts()
    {
        // All-N ref/alt alleles carry no nucleotide information, but the classifier
        // reads ONLY the integer read counts. A clear somatic signal must still be
        // called Somatic, and the (uninformative) alleles must be carried through
        // untouched.
        var call = Classify(Obs(30, 100, 0, 100, refA: "N", altA: "N"));

        AssertWellFormedCall(call);
        call.Variant.ReferenceAllele.Should().Be("N");
        call.Variant.AlternateAllele.Should().Be("N");
        call.TumorVaf.Should().BeApproximately(0.30, 1e-12);
        call.Status.Should().Be(SomaticStatus.Somatic);
    }

    [Test]
    public void Classify_EmptyAndNullAlleleStrings_NoCrash()
    {
        // Degenerate allele strings (empty, "NNNN", null) must not crash — the
        // decision is purely numeric.
        var act1 = () => Classify(Obs(20, 100, 0, 100, refA: "", altA: ""));
        var act2 = () => Classify(Obs(20, 100, 0, 100, refA: "NNNN", altA: "NNNN"));
        var act3 = () => Classify(Obs(20, 100, 0, 100, refA: null!, altA: null!));

        act1.Should().NotThrow();
        act2.Should().NotThrow();
        act3.Should().NotThrow();
        AssertWellFormedCall(act1());
        AssertWellFormedCall(act2());
        AssertWellFormedCall(act3());
    }

    #endregion

    #region ONCO-SOMATIC-001 — invariants under broad random fuzzing

    [Test]
    [CancelAfter(30000)]
    public void CallSomaticMutations_BroadRandomFuzz_AlwaysWellFormedAndConsistent()
    {
        // Broad fuzz over well-formed (0 ≤ alt ≤ total) random batches mixing the
        // full range of depths (including 0) and arbitrary all-N / empty alleles.
        // Asserts the full documented contract on every call AND the decision-rule
        // consistency (INV-01, INV-02, §4.2), independent of the classifier's
        // internal branch order.
        var rng = new Random(13579);
        string[] alleleChoices = { "A", "N", "NN", "", "ACGT" };

        for (int iter = 0; iter < 400; iter++)
        {
            int n = rng.Next(0, 12);
            var variants = new VariantObservation[n];
            for (int i = 0; i < n; i++)
            {
                int tTotal = rng.Next(0, 2000);
                int tAlt = tTotal == 0 ? 0 : rng.Next(0, tTotal + 1);
                int nTotal = rng.Next(0, 2000);
                int nAlt = nTotal == 0 ? 0 : rng.Next(0, nTotal + 1);
                variants[i] = Obs(
                    tAlt, tTotal, nAlt, nTotal,
                    refA: alleleChoices[rng.Next(alleleChoices.Length)],
                    altA: alleleChoices[rng.Next(alleleChoices.Length)]);
            }

            var calls = CallSomaticMutations(variants);
            calls.Should().HaveCount(n, "output preserves input count (§3.3)");

            for (int i = 0; i < n; i++)
            {
                var c = calls[i];
                AssertWellFormedCall(c);

                double ft = c.TumorVaf;
                double fn = c.NormalVaf;

                // Re-derive the documented decision rule and compare (§4.2).
                SomaticStatus expected =
                    ft < TumorTau ? SomaticStatus.NotDetected
                    : fn <= NormalTau ? SomaticStatus.Somatic
                    : SomaticStatus.Germline;
                c.Status.Should().Be(expected);

                // INV-01 ⇔ characterization of Somatic.
                if (c.Status == SomaticStatus.Somatic)
                {
                    ft.Should().BeGreaterThanOrEqualTo(TumorTau);
                    fn.Should().BeLessThanOrEqualTo(NormalTau);
                }
            }

            // INV-04: FilterGermlineVariants ≡ the Somatic subset, in input order.
            var filtered = FilterGermlineVariants(variants);
            var expectedSomatic = calls.Where(c => c.Status == SomaticStatus.Somatic).ToList();
            filtered.Should().HaveCount(expectedSomatic.Count);
            for (int j = 0; j < filtered.Count; j++)
            {
                filtered[j].Status.Should().Be(SomaticStatus.Somatic);
                filtered[j].Variant.Should().Be(expectedSomatic[j].Variant);
            }
        }
    }

    [Test]
    [CancelAfter(20000)]
    public void Classify_RandomThresholds_NeverViolatesScoreBounds()
    {
        // Fuzz the thresholds themselves across their valid [0, 1] domain together
        // with random well-formed evidence: the score must stay in [0, 1] and equal
        // max(0, f_t − f_n) for any Somatic call regardless of the cutoffs used.
        var rng = new Random(271828);
        for (int i = 0; i < 5000; i++)
        {
            int tTotal = rng.Next(1, 1000);
            int tAlt = rng.Next(0, tTotal + 1);
            int nTotal = rng.Next(1, 1000);
            int nAlt = rng.Next(0, nTotal + 1);
            double tauT = rng.NextDouble();
            double tauN = rng.NextDouble();

            var call = Classify(Obs(tAlt, tTotal, nAlt, nTotal), tauT, tauN);

            AssertWellFormedCall(call);
            if (call.Status == SomaticStatus.Somatic)
            {
                call.SomaticScore.Should()
                    .BeApproximately(Math.Max(0.0, call.TumorVaf - call.NormalVaf), 1e-12);
            }
        }
    }

    [Test]
    public void Classify_OutOfRangeThreshold_ThrowsArgumentOutOfRange()
    {
        // Thresholds outside [0, 1] (NaN, > 1, < 0) are a documented rejection (§3.3).
        var actHigh = () => Classify(Obs(30, 100, 0, 100), tumorVafThreshold: 1.5);
        var actLow = () => Classify(Obs(30, 100, 0, 100), normalVafThreshold: -0.1);
        var actNaN = () => Classify(Obs(30, 100, 0, 100), tumorVafThreshold: double.NaN);

        actHigh.Should().Throw<ArgumentOutOfRangeException>();
        actLow.Should().Throw<ArgumentOutOfRangeException>();
        actNaN.Should().Throw<ArgumentOutOfRangeException>();
    }

    #endregion
}
