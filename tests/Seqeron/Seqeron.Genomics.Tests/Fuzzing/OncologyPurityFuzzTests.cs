using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology tumour-purity-estimation area — ONCO-PURITY-001.
/// The unit under test is tumour purity ρ estimation, implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs as three
/// public entry points:
///   • <see cref="OncologyAnalyzer.EstimatePurityFromVaf(double)"/> —
///       single clonal heterozygous copy-neutral diploid SNV, ρ = 2·VAF;
///   • <see cref="OncologyAnalyzer.EstimatePurityFromVAF(IEnumerable{VariantObservation})"/> —
///       median of ρ = 2·VAF over read-count variants;
///   • <see cref="OncologyAnalyzer.EstimatePurity(IEnumerable{PurityVariant})"/> —
///       median of the allele-specific inversion ρ = 2v/[m + v(2 − n_tot)].
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no nonsense
/// output, and no *unhandled* runtime exception (DivideByZero / Overflow / NaN).
/// Every input must resolve to EITHER a well-defined, theory-correct purity OR a
/// *documented, intentional* outcome (an <see cref="ArgumentNullException"/> for
/// a null collection, <see cref="ArgumentException"/> for an empty collection,
/// <see cref="ArgumentOutOfRangeException"/> for a VAF/(m,n_tot) state that
/// cannot map to a purity in [0, 1]).
/// For tumour purity the headline hazards are:
///   • a purity that ESCAPES [0, 1] — INV-01 says purity is a fraction of cells
///     and must stay in [0, 1]; in particular VAF > 0.5 under the diploid model
///     implies ρ > 1 and is a documented throw (§3.3, §6.1), never a leaked 1.4;
///   • a DivideByZero / NaN when the denominator m + v(2 − n_tot) reaches 0 in
///     the allele-specific inversion (documented throw, §3.3), or 0/0 on an
///     EMPTY collection (documented ArgumentException, §6.1) — NOT a NaN median;
///   • a SINGLE-variant collection crashing the median aggregation
///     (variance-of-one / empty-after-filter) — the median of one value is that
///     value, and the single-VAF overload is O(1) closed form (§4.3);
///   • a NaN VAF silently flowing through 2·VAF into a NaN purity (documented
///     ArgumentOutOfRangeException, §3.3).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-PURITY-001 — Tumour purity estimation (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 106.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 106): "all-VAF=0, all-VAF=1, single variant".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Mapping of the BE targets onto the documented contract:
///   • "all-VAF=0"  ⇒ no variant evidence ⇒ ρ = 2·0 = 0 (INV-04, §6.1); the
///       median of all-zero per-variant purities is 0, never NaN.
///   • "all-VAF=1"  ⇒ fully clonal/homozygous upper boundary. Under the diploid
///       model VAF = 1 implies ρ = 2 > 1 and is a documented throw (§3.3, §6.1) —
///       the upper boundary that PRODUCES ρ = 1 is VAF = 0.5 (ρ = 2·0.5 = 1,
///       §6.1). Under the allele-specific model VAF = 1 IS reachable as a fully
///       clonal LOH locus (m = 2, n_tot = 2 ⇒ ρ = 1, clamped to ≤ 1).
///   • "single variant" ⇒ the median of a one-element collection equals that
///       element; no variance-of-one / divide-by-count-minus-one crash.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Tumor_Purity_Estimation.md (docs/algorithms/Oncology/Tumor_Purity_Estimation.md):
///   • ρ = 2·VAF for a clonal het copy-neutral diploid SNV (m=1, n_tot=2)
///       (§2.2, INV-02; CNAqc: purity 60% ⇔ VAF 30%).
///   • π = 2v/[m + v(2 − n_tot)] — exact inverse of v = mπ/[2(1−π)+π·n_tot]
///       (§2.2, INV-03).
///   • 0 ≤ purity ≤ 1 (INV-01); inputs yielding ρ outside [0,1] are rejected.
///   • VAF = 0 ⇒ purity = 0; estimate is monotone non-decreasing in VAF
///       (§2.4 INV-04, §6.1).
///   • VAF = 0.5 (diploid) ⇒ purity = 1.0; VAF > 0.5 ⇒ ArgumentOutOfRangeException
///       (§6.1).
///   • Collection overloads aggregate per-variant purities by their MEDIAN
///       (lower-mid average for even counts) (§4.1, §5.2).
///   • Null collection ⇒ ArgumentNullException; empty collection ⇒
///       ArgumentException (purity undefined) (§3.3, §6.1).
///   • Allele-specific: m < 1 / n_tot < 1 / non-positive denominator /
///       ρ ∉ [0,1] ⇒ ArgumentOutOfRangeException (§3.3, §6.1).
///   • Worked example: EstimatePurityFromVaf(0.30) = 0.60; the 2:1 segment
///       PurityVariant(2/3, m=2, n_tot=3) ⇒ ρ = 1.0 (§7.1).
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyPurityFuzzTests
{
    // ── Well-formed-purity assertion helper ──────────────────────────────────
    // Pins the documented numeric contract on EVERY accepted purity: finite,
    // never NaN/Infinity (no DivideByZero leak, no 0/0 median) and inside [0, 1]
    // (INV-01). This is what stops a fuzz test from rubber-stamping a NaN or an
    // out-of-range purity such as 1.4.
    private static void AssertWellFormedPurity(double purity)
    {
        double.IsNaN(purity).Should().BeFalse("purity must never be NaN (no 0/0 in the inversion or median)");
        double.IsInfinity(purity).Should().BeFalse("purity must be finite");
        purity.Should().BeInRange(0.0, 1.0, "purity is a fraction of cells ⇒ ρ ∈ [0, 1] (INV-01)");
    }

    // A clonal heterozygous copy-neutral diploid SNV at the given VAF, expressed
    // as a read-count VariantObservation (alt = round(vaf·depth), total = depth)
    // for the collection overload. Used to build all-VAF=0 / all-VAF=0.5 fixtures.
    private static VariantObservation HetDiploidVariant(double vaf, int depth)
    {
        int alt = (int)Math.Round(vaf * depth);
        return new VariantObservation(
            Chromosome: "chr1",
            Position: 1,
            ReferenceAllele: "A",
            AlternateAllele: "T",
            TumorAltReads: alt,
            TumorTotalReads: depth,
            NormalAltReads: 0,
            NormalTotalReads: depth);
    }

    #region ONCO-PURITY-001 — Positive sanity (documented formula on hand-built examples)

    [Test]
    public void EstimatePurityFromVaf_DocWorkedExample_RhoIsTwiceVaf()
    {
        // Docs §2.2 / §7.1: a real purity of 60% corresponds to VAF 30%
        // (ρ = 2·VAF). This pins the headline closed form, not just "in range".
        EstimatePurityFromVaf(0.30).Should().BeApproximately(0.60, 1e-12);

        // CNAqc 55–65% purity band ⇔ 27.5–32.5% VAF band (§2.2).
        EstimatePurityFromVaf(0.275).Should().BeApproximately(0.55, 1e-12);
        EstimatePurityFromVaf(0.325).Should().BeApproximately(0.65, 1e-12);
    }

    [Test]
    public void EstimatePurityFromVAF_ClusteredHetVafs_MedianIsTwiceVaf()
    {
        // A set of clonal heterozygous somatic VAFs clustered at ~0.30 must yield
        // purity ≈ 0.60 via the documented median-of-(2·VAF) relationship (§4.1).
        var variants = new[]
        {
            HetDiploidVariant(0.29, 1000),
            HetDiploidVariant(0.30, 1000),
            HetDiploidVariant(0.31, 1000),
        };

        double purity = EstimatePurityFromVAF(variants);

        AssertWellFormedPurity(purity);
        // Median VAF is 0.30 (alt 300/1000) ⇒ ρ = 0.60, exact for these counts.
        purity.Should().BeApproximately(0.60, 1e-9);
    }

    [Test]
    public void EstimatePurity_AlleleSpecific_DocLohExample_RhoIsOne()
    {
        // Docs §7.1 numerical walk-through: a 2:1 / LOH segment with n_tot = 3,
        // m = 2, VAF = 2/3 corresponds to fully pure tumour (ρ = 1.0):
        // π = 2·(2/3)/[2 + (2/3)(2−3)] = (4/3)/(4/3) = 1.0.
        double rho = EstimatePurity(new[]
        {
            new PurityVariant(2.0 / 3.0, Multiplicity: 2, TumorTotalCopyNumber: 3),
        });

        AssertWellFormedPurity(rho);
        rho.Should().BeApproximately(1.0, 1e-12);
    }

    [Test]
    public void EstimatePurity_AlleleSpecific_InvertsTheForwardModel()
    {
        // INV-03: the inversion recovers the purity that generated the VAF.
        // Build a forward VAF v = mπ/[2(1−π) + π·n_tot] for chosen (π, m, n_tot),
        // then assert EstimatePurity recovers π exactly.
        var rng = new Random(106_001);
        for (int i = 0; i < 200; i++)
        {
            double pi = 0.05 + rng.NextDouble() * 0.9;     // π ∈ (0.05, 0.95)
            int m = 1 + rng.Next(3);                       // m ∈ {1,2,3}
            int nTot = Math.Max(m, 1 + rng.Next(4));       // n_tot ≥ m, ≥ 1
            double v = m * pi / (2.0 * (1.0 - pi) + pi * nTot);
            if (v < 0.0 || v > 1.0)
            {
                continue;
            }

            double recovered = EstimatePurity(new[]
            {
                new PurityVariant(v, Multiplicity: m, TumorTotalCopyNumber: nTot),
            });

            AssertWellFormedPurity(recovered);
            recovered.Should().BeApproximately(pi, 1e-9,
                "the inversion is the exact algebraic inverse of the forward model (INV-03)");
        }
    }

    #endregion

    #region ONCO-PURITY-001 — BE: all-VAF=0 (no variant evidence ⇒ purity 0, no DivideByZero/NaN)

    [Test]
    [CancelAfter(10_000)]
    public void EstimatePurityFromVaf_ZeroVaf_ReturnsZero()
    {
        // VAF = 0 ⇒ ρ = 2·0 = 0 (INV-04, §6.1). The lower boundary; must be an
        // exact 0, not a NaN or a tiny negative from float noise.
        double purity = EstimatePurityFromVaf(0.0);

        AssertWellFormedPurity(purity);
        purity.Should().Be(0.0);
    }

    [Test]
    [CancelAfter(10_000)]
    public void EstimatePurityFromVAF_AllVafZero_MedianIsZero_NoNaN()
    {
        // Every variant has zero alt support ⇒ every per-variant ρ = 0 ⇒ the
        // median is 0. The hazard is a 0/0 NaN slipping through the median; the
        // documented result is purity 0 (no variant evidence, INV-04).
        var rng = new Random(106_002);
        int count = 1 + rng.Next(12);
        var variants = Enumerable.Range(0, count)
            .Select(_ => HetDiploidVariant(0.0, 100 + rng.Next(900)))
            .ToArray();

        double purity = EstimatePurityFromVAF(variants);

        AssertWellFormedPurity(purity);
        purity.Should().Be(0.0, "no variant evidence ⇒ purity 0 (INV-04)");
    }

    [Test]
    [CancelAfter(10_000)]
    public void EstimatePurity_AlleleSpecific_AllVafZero_MedianIsZero()
    {
        // Allele-specific overload with VAF = 0 across mixed (m, n_tot) states:
        // ρ = 2·0/[m + 0] = 0 for every state ⇒ median 0, no DivideByZero.
        var rng = new Random(106_003);
        int count = 1 + rng.Next(10);
        var variants = Enumerable.Range(0, count)
            .Select(_ => new PurityVariant(0.0, Multiplicity: 1 + rng.Next(3), TumorTotalCopyNumber: 1 + rng.Next(4)))
            .ToArray();

        double purity = EstimatePurity(variants);

        AssertWellFormedPurity(purity);
        purity.Should().Be(0.0);
    }

    #endregion

    #region ONCO-PURITY-001 — BE: all-VAF=1 (fully clonal/homozygous upper boundary, no >1 leak)

    [Test]
    [CancelAfter(10_000)]
    public void EstimatePurityFromVaf_VafHalf_ReturnsExactlyOne()
    {
        // The diploid upper boundary that PRODUCES ρ = 1 is VAF = 0.5
        // (ρ = 2·0.5 = 1, §6.1). Must be exactly 1.0, not 0.9999… or a >1 leak.
        double purity = EstimatePurityFromVaf(0.5);

        AssertWellFormedPurity(purity);
        purity.Should().Be(1.0);
    }

    [Test]
    [CancelAfter(10_000)]
    public void EstimatePurityFromVaf_VafOne_DiploidModel_Throws_NoPurityAboveOne()
    {
        // VAF = 1 under the diploid het model implies ρ = 2 > 1, which is
        // impossible (INV-01). Documented contract: ArgumentOutOfRangeException
        // (§3.3, §6.1) — the result must NEVER be a leaked purity of 2.0.
        var act = () => EstimatePurityFromVaf(1.0);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    [CancelAfter(10_000)]
    public void EstimatePurityFromVAF_AllVafOne_DiploidModel_Throws()
    {
        // A collection where every variant is homozygous (alt == total ⇒ VAF 1)
        // hits the same diploid > 0.5 guard on the first variant ⇒ documented
        // ArgumentOutOfRangeException, never a median above 1.
        var variants = new[]
        {
            HetDiploidVariant(1.0, 500),
            HetDiploidVariant(1.0, 500),
        };

        var act = () => EstimatePurityFromVAF(variants);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    [CancelAfter(10_000)]
    public void EstimatePurity_AlleleSpecific_AllVafOne_FullyClonalLoh_ClampedToOne()
    {
        // "all-VAF=1" IS reachable under the allele-specific model as a fully
        // clonal LOH locus: m = 2, n_tot = 2 ⇒
        // π = 2·1/[2 + 1·(2−2)] = 2/2 = 1.0. The upper bound, clamped to ≤ 1 —
        // no >1 leak even when every variant is fully clonal/homozygous.
        var variants = new[]
        {
            new PurityVariant(1.0, Multiplicity: 2, TumorTotalCopyNumber: 2),
            new PurityVariant(1.0, Multiplicity: 2, TumorTotalCopyNumber: 2),
            new PurityVariant(1.0, Multiplicity: 2, TumorTotalCopyNumber: 2),
        };

        double purity = EstimatePurity(variants);

        AssertWellFormedPurity(purity);
        purity.Should().Be(1.0, "fully clonal LOH (m=2, n_tot=2) at VAF 1 ⇒ ρ = 1 (INV-01 upper bound)");
    }

    [Test]
    [CancelAfter(10_000)]
    public void EstimatePurity_AlleleSpecific_VafOne_OnAmplifiedSegment_ThrowsNotLeakAboveOne()
    {
        // VAF = 1 on an amplified diploid-multiplicity-1 segment (m = 1,
        // n_tot = 4): π = 2·1/[1 + 1·(2−4)] = 2/(−1) = −2 ⇒ non-positive
        // denominator / ρ ∉ [0,1] ⇒ documented ArgumentOutOfRangeException
        // (§3.3) — NOT a negative or >1 purity leak.
        var act = () => EstimatePurity(new[]
        {
            new PurityVariant(1.0, Multiplicity: 1, TumorTotalCopyNumber: 4),
        });

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    #endregion

    #region ONCO-PURITY-001 — BE: single variant (median-of-one, no variance-of-one crash)

    [Test]
    [CancelAfter(10_000)]
    public void EstimatePurityFromVAF_SingleVariant_MedianEqualsThatVariant()
    {
        // A one-element collection must not crash the median aggregation
        // (no divide-by-(n−1), no empty-after-filter). The median of one value is
        // that value ⇒ ρ = 2·VAF for the single variant.
        var rng = new Random(106_004);
        for (int i = 0; i < 100; i++)
        {
            double vaf = rng.NextDouble() * 0.5;          // diploid-valid VAF ∈ [0, 0.5]
            int depth = 100 + rng.Next(900);
            int alt = (int)Math.Round(vaf * depth);
            double expected = 2.0 * ((double)alt / depth); // ρ = 2·(alt/total)

            var single = new[] { HetDiploidVariant(vaf, depth) };

            double purity = EstimatePurityFromVAF(single);

            AssertWellFormedPurity(purity);
            purity.Should().BeApproximately(expected, 1e-9,
                "median of a single per-variant purity is that purity");
        }
    }

    [Test]
    [CancelAfter(10_000)]
    public void EstimatePurity_AlleleSpecific_SingleVariant_NoCrash()
    {
        // Single-element allele-specific collection across random valid states:
        // never throws on valid (m, n_tot, VAF) that map into [0,1]; result is
        // well-formed and equals the single per-variant inversion.
        var rng = new Random(106_005);
        for (int i = 0; i < 150; i++)
        {
            double pi = 0.1 + rng.NextDouble() * 0.8;
            int m = 1 + rng.Next(2);
            int nTot = Math.Max(m, 1 + rng.Next(3));
            double v = m * pi / (2.0 * (1.0 - pi) + pi * nTot);
            if (v < 0.0 || v > 1.0)
            {
                continue;
            }

            var single = new[] { new PurityVariant(v, Multiplicity: m, TumorTotalCopyNumber: nTot) };

            double purity = EstimatePurity(single);

            AssertWellFormedPurity(purity);
            purity.Should().BeApproximately(pi, 1e-9);
        }
    }

    #endregion

    #region ONCO-PURITY-001 — BE: empty / null collections (documented throws, no 0/0 median)

    [Test]
    public void EstimatePurityFromVAF_EmptyCollection_ThrowsArgumentException()
    {
        // Median of an empty set is 0/0 (NaN); the documented contract is a
        // throw (purity undefined, §3.3 §6.1) — NOT a NaN median.
        var act = () => EstimatePurityFromVAF(Array.Empty<VariantObservation>());

        act.Should().Throw<ArgumentException>();
    }

    [Test]
    public void EstimatePurity_AlleleSpecific_EmptyCollection_ThrowsArgumentException()
    {
        var act = () => EstimatePurity(Array.Empty<PurityVariant>());

        act.Should().Throw<ArgumentException>();
    }

    [Test]
    public void EstimatePurityFromVAF_NullCollection_ThrowsArgumentNullException()
    {
        var act = () => EstimatePurityFromVAF(null!);

        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    public void EstimatePurity_AlleleSpecific_NullCollection_ThrowsArgumentNullException()
    {
        var act = () => EstimatePurity(null!);

        act.Should().Throw<ArgumentNullException>();
    }

    #endregion

    #region ONCO-PURITY-001 — BE: malformed scalars (NaN / negative / out-of-range VAF, bad m,n_tot)

    [Test]
    public void EstimatePurityFromVaf_NaNVaf_ThrowsArgumentOutOfRange_NoNaNPurity()
    {
        // A NaN VAF must NOT flow through 2·VAF into a NaN purity; documented
        // ArgumentOutOfRangeException (§3.3).
        var act = () => EstimatePurityFromVaf(double.NaN);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void EstimatePurityFromVaf_NegativeVaf_ThrowsArgumentOutOfRange()
    {
        var act = () => EstimatePurityFromVaf(-0.01);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void EstimatePurityFromVaf_VafAboveOne_ThrowsArgumentOutOfRange()
    {
        var act = () => EstimatePurityFromVaf(1.5);

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void EstimatePurity_AlleleSpecific_MultiplicityBelowOne_ThrowsArgumentOutOfRange()
    {
        var act = () => EstimatePurity(new[]
        {
            new PurityVariant(0.3, Multiplicity: 0, TumorTotalCopyNumber: 2),
        });

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    [Test]
    public void EstimatePurity_AlleleSpecific_CopyNumberBelowOne_ThrowsArgumentOutOfRange()
    {
        var act = () => EstimatePurity(new[]
        {
            new PurityVariant(0.3, Multiplicity: 1, TumorTotalCopyNumber: 0),
        });

        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    #endregion

    #region ONCO-PURITY-001 — Invariant sweep (purity stays in [0,1]; monotone non-decreasing in VAF)

    [Test]
    [CancelAfter(20_000)]
    public void EstimatePurityFromVaf_RandomValidVafs_AlwaysInRange_AndMonotone()
    {
        // INV-01 / INV-04 sweep: across many valid diploid VAFs the purity stays
        // finite and in [0,1], and is monotone non-decreasing in VAF.
        var rng = new Random(106_006);
        for (int i = 0; i < 500; i++)
        {
            double a = rng.NextDouble() * 0.5;
            double b = rng.NextDouble() * 0.5;
            double lo = Math.Min(a, b);
            double hi = Math.Max(a, b);

            double pLo = EstimatePurityFromVaf(lo);
            double pHi = EstimatePurityFromVaf(hi);

            AssertWellFormedPurity(pLo);
            AssertWellFormedPurity(pHi);
            pHi.Should().BeGreaterThanOrEqualTo(pLo - 1e-12,
                "purity is monotone non-decreasing in VAF for fixed m, n_tot (INV-04)");
        }
    }

    #endregion
}
