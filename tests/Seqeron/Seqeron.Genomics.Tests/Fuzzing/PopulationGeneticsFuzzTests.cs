using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Population;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the PopGen area — allele frequency.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (DivideByZeroException, OverflowException,
/// NullReferenceException, …) and no invalid *numeric* result (NaN/±Infinity, a
/// frequency outside its mathematical range, a pair that fails to sum to 1).
/// Every input must result in EITHER a well-defined, theory-correct value, OR a
/// *documented, intentional* validation exception
/// (ArgumentException / ArgumentOutOfRangeException). A raw runtime exception, a
/// silent NaN/garbage frequency, or a hang on garbage input is a bug, not a
/// passing test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: POP-FREQ-001 — allele frequency (PopGen)
/// Checklist: docs/checklists/03_FUZZING.md, row 43.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — the degenerate population boundaries that
///           stress the frequency arithmetic: 0 samples (the DivideByZero
///           boundary — total allele count == 0), 1 individual / 1 allele
///           (the monomorphic, fully-fixed locus), "all the same allele"
///           (every sampled chromosome the same → frequency 1.0 for it, 0 for
///           the other), and negative counts (-1, the biologically impossible
///           input that must be rejected or carried as a *defined* — not a
///           silently-garbage — result).
/// — docs/checklists/03_FUZZING.md §Description ("BE = граничні значення: 0, -1,
///   MaxInt, empty").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The allele-frequency contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Allele frequency is the relative abundance of an allele at one diploid
/// biallelic locus: each individual contributes two chromosome copies, so for
/// genotype counts n_AA (homozygous major), n_Aa (heterozygous), n_aa
/// (homozygous minor):
///     p = (2·n_AA + n_Aa) / (2·(n_AA + n_Aa + n_aa))
///     q = (2·n_aa + n_Aa) / (2·(n_AA + n_Aa + n_aa))
/// and the minor allele frequency MAF = min(f, 1 − f).
///   — docs/algorithms/Population_Genetics/Allele_Frequency.md §2.2, §2.4
///     (INV-01: p + q = 1 for any non-empty sample — THE key invariant;
///      INV-02: 0 ≤ p,q ≤ 1; INV-03: 0 ≤ MAF ≤ 0.5; INV-04: major + minor copies
///      equal the total allele count). Sources: Gillespie (2004) [1];
///      Wikipedia "Allele frequency" [2]; "Genotype frequency" [4].
///
/// POP-FREQ-001 spans TWO documented entry points with DIFFERENT, intentional
/// validation contracts. Fuzzing pins both, and the boundary between them, so
/// neither can silently drift:
///
/// (1) The STRICT genotype-count surface —
///     PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(int, int, int)
///     (src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs
///      lines 138–159). Documented validation (Allele_Frequency.md §3.1, §3.3,
///      §6.1):
///       • ANY negative genotype count (−1, int.MinValue, …) →
///         a *documented, intentional* ArgumentOutOfRangeException — the
///         biologically-impossible-input gate, NOT a silent negative frequency
///         and NOT a crash. This is the KEY negative-counts fuzz target.
///       • 0 samples (all three counts 0) → total allele count == 0 → the
///         explicit `totalAlleles == 0` guard returns the DEFINED (0, 0); there
///         is NO division, so no DivideByZeroException and no NaN. This is the
///         KEY 0-samples fuzz concern.
///       • 1 individual, all the same allele (homozygous major) → (1.0, 0.0);
///         all homozygous minor → (0.0, 1.0); all heterozygous → (0.5, 0.5).
///       • for ANY non-empty, non-negative count triple the result obeys
///         p + q = 1 exactly (INV-01), both in [0, 1] (INV-02), even at scale
///         (large counts must not overflow the int allele-copy arithmetic into a
///         wrong/negative total). The total allele copies = 2·(sum); large
///         inputs are kept under the int boundary so 2·sum cannot overflow — we
///         pin the in-range arithmetic, not an overflowing one.
///
/// (2) The LENIENT genotype-code surface —
///     PopulationGeneticsAnalyzer.CalculateMAF(IEnumerable&lt;int&gt;)
///     (PopulationGeneticsAnalyzer.cs lines 164–177). LENIENT by documented
///     design (Allele_Frequency.md §3.3, §5.2): it assumes the 0/1/2 diploid
///     genotype encoding (0 = hom ref, 1 = het, 2 = hom alt) and does NOT
///     validate each element. Documented behavior:
///       • empty genotype list → 0 (explicit short-circuit; no division by zero).
///       • monomorphic list (all 0 → no alt copies; or all 2 → all alt copies) →
///         MAF 0: a fixed locus has no minor allele (§6.1).
///       • well-formed 0/1/2 input → MAF = min(f, 1 − f) ∈ [0, 0.5] (INV-03).
///       • because the surface does NOT validate, an OUT-OF-ENCODING code (a
///         negative genotype, or a code &gt; 2) is summed verbatim. Fuzzing pins
///         the *actual defined* consequence — the method stays total: it does
///         NOT throw, does NOT divide by zero (the list is non-empty so the
///         denominator 2·count &gt; 0), and does NOT produce NaN/Infinity. It may
///         return a value OUTSIDE [0, 0.5] for such malformed input — a
///         documented consequence of the "caller must normalize" contract
///         (§5.3 "Intentionally simplified"), NOT a crash. We pin that this
///         surface never crashes and never NaNs on the malformed-code boundary,
///         so the strict/lenient split is explicit and cannot drift into a
///         silent DivideByZero.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class PopulationGeneticsFuzzTests
{
    #region Helpers

    private const double Tolerance = 1e-9;

    /// <summary>Deterministic RNG — seed fixed so generated fuzz inputs are reproducible.</summary>
    private static readonly Random Rng = new(20260620);

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  POP-FREQ-001 — allele frequency : fuzz targets
    //  Strict surface:  CalculateAlleleFrequencies(int, int, int)
    //  Lenient surface: CalculateMAF(IEnumerable<int>)
    // ═══════════════════════════════════════════════════════════════════

    #region POP-FREQ-001 — allele frequency

    #region Positive sanity — known genotype set yields the expected frequencies summing to 1

    /// <summary>
    /// Positive control: the worked example from the algorithm doc
    /// (Allele_Frequency.md §7.1, the four-o'clock plant) — 49 AA, 42 Aa, 9 aa
    /// → 200 allele copies, 140 of A, 60 of a → p(A) = 0.70, q(a) = 0.30. The
    /// frequencies must be exactly those values, must each lie in [0, 1]
    /// (INV-02), and must sum to 1 (INV-01). The same population, expressed as a
    /// 0/1/2 genotype list, must give MAF = min(0.70, 0.30) = 0.30 (INV-03).
    /// This anchors the fuzz battery: before pinning that garbage does no harm,
    /// confirm a known-good input gives the textbook-correct answer.
    /// </summary>
    [Test]
    public void AlleleFrequencies_KnownGenotypeSet_MatchesWorkedExampleAndSumsToOne()
    {
        var (p, q) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(
            homozygousMajor: 49, heterozygous: 42, homozygousMinor: 9);

        p.Should().BeApproximately(0.70, Tolerance,
            because: "p(A) = (2·49 + 42)/200 = 140/200 = 0.70 (Allele_Frequency.md §7.1)");
        q.Should().BeApproximately(0.30, Tolerance,
            because: "q(a) = (2·9 + 42)/200 = 60/200 = 0.30 (Allele_Frequency.md §7.1)");

        p.Should().BeInRange(0.0, 1.0, "each allele frequency is a proportion (INV-02)");
        q.Should().BeInRange(0.0, 1.0, "each allele frequency is a proportion (INV-02)");
        (p + q).Should().BeApproximately(1.0, Tolerance,
            because: "every sampled chromosome is one allele or the other → p + q = 1 (INV-01)");

        // The same population as a 0/1/2 genotype list: 49×0, 42×1, 9×2.
        var genotypes = Enumerable.Repeat(0, 49)
            .Concat(Enumerable.Repeat(1, 42))
            .Concat(Enumerable.Repeat(2, 9));
        double maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

        maf.Should().BeApproximately(0.30, Tolerance,
            because: "MAF = min(p, q) = min(0.70, 0.30) = 0.30 (INV-03)");
        maf.Should().BeInRange(0.0, 0.5, "MAF is the smaller of two complementary frequencies (INV-03)");
    }

    #endregion

    #region BE — 0 samples (the DivideByZero boundary: total allele count == 0)

    /// <summary>
    /// BE / DivideByZero boundary: with 0 samples every genotype count is 0, so
    /// the total allele count 2·(0+0+0) is 0. The contract must NOT divide by
    /// zero and must NOT return NaN — the explicit `totalAlleles == 0` guard
    /// returns the DEFINED (0, 0) (Allele_Frequency.md §3.3, §6.1;
    /// PopulationGeneticsAnalyzer.cs lines 150–153). This is the single most
    /// important fuzz concern for this unit: an unguarded denominator here is a
    /// DivideByZeroException crash in production.
    /// </summary>
    [Test]
    public void AlleleFrequencies_ZeroSamples_ReturnsDefinedZeroPairWithNoDivideByZero()
    {
        (double p, double q) result = default;
        var act = () => result = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(0, 0, 0);

        act.Should().NotThrow<DivideByZeroException>(
            "the zero-total-allele boundary is guarded; the denominator is never divided by zero");
        act.Should().NotThrow("0 samples is a defined boundary input (returns (0,0)), not an error");

        result.p.Should().Be(0.0, "no alleles are present → major frequency is the defined 0");
        result.q.Should().Be(0.0, "no alleles are present → minor frequency is the defined 0");
        double.IsNaN(result.p).Should().BeFalse("the guard avoids a 0/0 NaN");
        double.IsNaN(result.q).Should().BeFalse("the guard avoids a 0/0 NaN");
    }

    /// <summary>
    /// BE / DivideByZero boundary, lenient surface: an EMPTY genotype list is the
    /// 0-samples boundary for CalculateMAF. The denominator would be 2·0 = 0, so
    /// the contract short-circuits to the DEFINED 0 BEFORE any division
    /// (Allele_Frequency.md §3.3, §6.1; PopulationGeneticsAnalyzer.cs lines
    /// 168–169) — no DivideByZeroException, no NaN.
    /// </summary>
    [Test]
    public void Maf_EmptyGenotypeList_ReturnsDefinedZeroWithNoDivideByZero()
    {
        double maf = double.NaN;
        var act = () => maf = PopulationGeneticsAnalyzer.CalculateMAF(Array.Empty<int>());

        act.Should().NotThrow<DivideByZeroException>(
            "the empty-list boundary short-circuits before the 2·count denominator is used");
        act.Should().NotThrow("an empty genotype list is a defined boundary input, not an error");

        maf.Should().Be(0.0, "an empty sample has no minor allele frequency (defined 0)");
        double.IsNaN(maf).Should().BeFalse("the empty short-circuit avoids a 0/0 NaN");
    }

    #endregion

    #region BE — 1 allele / 1 individual (monomorphic, fully-fixed locus)

    /// <summary>
    /// BE: the minimal non-empty population — a single individual. The locus is
    /// monomorphic (fully fixed): one homozygous-major individual fixes the major
    /// allele at frequency 1.0 and the minor at 0.0; one homozygous-minor fixes
    /// the reverse; one heterozygote is the balanced 0.5/0.5. In every case the
    /// pair must sum to 1 (INV-01) and each frequency stays in [0, 1] (INV-02).
    /// This pins that the length-1 boundary neither off-by-ones the 2·count
    /// allele accounting nor escapes the [0,1] range.
    /// </summary>
    [TestCase(1, 0, 0, 1.0, 0.0, TestName = "AlleleFrequencies_SingleHomozygousMajor_IsFixedMajor")]
    [TestCase(0, 0, 1, 0.0, 1.0, TestName = "AlleleFrequencies_SingleHomozygousMinor_IsFixedMinor")]
    [TestCase(0, 1, 0, 0.5, 0.5, TestName = "AlleleFrequencies_SingleHeterozygote_IsBalanced")]
    public void AlleleFrequencies_SingleIndividual_IsMonomorphicOrBalanced(
        int hMajor, int het, int hMinor, double expectedP, double expectedQ)
    {
        var (p, q) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(hMajor, het, hMinor);

        p.Should().BeApproximately(expectedP, Tolerance,
            because: "the major-allele frequency of a single defined genotype is exact");
        q.Should().BeApproximately(expectedQ, Tolerance,
            because: "the minor-allele frequency of a single defined genotype is exact");

        p.Should().BeInRange(0.0, 1.0, "a frequency is a proportion (INV-02)");
        q.Should().BeInRange(0.0, 1.0, "a frequency is a proportion (INV-02)");
        (p + q).Should().BeApproximately(1.0, Tolerance, "p + q = 1 for any non-empty sample (INV-01)");
    }

    /// <summary>
    /// BE: monomorphic locus on the lenient MAF surface. A genotype list that is
    /// all homozygous reference (all 0) carries no alternate copies; a list that
    /// is all homozygous alternate (all 2) carries only alternate copies. Either
    /// way the locus is FIXED, so MAF = 0 — a fixed locus has no minor allele
    /// (Allele_Frequency.md §6.1 "Monomorphic genotype list → 0"). Verified for a
    /// single individual and a longer fixed run; MAF must stay in [0, 0.5].
    /// </summary>
    [TestCase(0, 1, TestName = "Maf_AllHomozygousReference_Single_IsZero")]
    [TestCase(0, 25, TestName = "Maf_AllHomozygousReference_Run_IsZero")]
    [TestCase(2, 1, TestName = "Maf_AllHomozygousAlternate_Single_IsZero")]
    [TestCase(2, 25, TestName = "Maf_AllHomozygousAlternate_Run_IsZero")]
    public void Maf_MonomorphicGenotypeList_IsZero(int code, int count)
    {
        double maf = PopulationGeneticsAnalyzer.CalculateMAF(Enumerable.Repeat(code, count));

        maf.Should().BeApproximately(0.0, Tolerance,
            because: "a fixed (monomorphic) locus has no minor allele → MAF = 0 (INV-03, §6.1)");
        maf.Should().BeInRange(0.0, 0.5, "MAF can never escape [0, 0.5] for valid 0/1/2 input (INV-03)");
    }

    #endregion

    #region BE — all the same allele (every chromosome copy identical → frequency 1.0)

    /// <summary>
    /// BE "all same allele": when every sampled chromosome carries the same
    /// allele, that allele's frequency is exactly 1.0 and the other is exactly
    /// 0.0. A population of N homozygous-major individuals fixes the major allele;
    /// N homozygous-minor fixes the minor. Verified across a fixed-seed sweep of
    /// population sizes so the result is independent of N, stays in [0, 1]
    /// (INV-02), and the pair sums to 1 (INV-01). This is the population-scale
    /// analogue of the homopolymer boundary.
    /// </summary>
    [Test]
    public void AlleleFrequencies_AllSameAllele_IsFrequencyOneForItZeroForOther()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int n = Rng.Next(1, 100_000);

            var (pMaj, qMaj) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(n, 0, 0);
            pMaj.Should().BeApproximately(1.0, Tolerance,
                because: $"all {n} individuals homozygous-major → major frequency 1.0 regardless of N");
            qMaj.Should().BeApproximately(0.0, Tolerance, "the absent allele has frequency 0");
            (pMaj + qMaj).Should().BeApproximately(1.0, Tolerance, "p + q = 1 (INV-01)");
            pMaj.Should().BeInRange(0.0, 1.0, "INV-02"); qMaj.Should().BeInRange(0.0, 1.0, "INV-02");

            var (pMin, qMin) = PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(0, 0, n);
            pMin.Should().BeApproximately(0.0, Tolerance, "the absent allele has frequency 0");
            qMin.Should().BeApproximately(1.0, Tolerance,
                because: $"all {n} individuals homozygous-minor → minor frequency 1.0 regardless of N");
            (pMin + qMin).Should().BeApproximately(1.0, Tolerance, "p + q = 1 (INV-01)");
            pMin.Should().BeInRange(0.0, 1.0, "INV-02"); qMin.Should().BeInRange(0.0, 1.0, "INV-02");
        }
    }

    /// <summary>
    /// BE "all same allele", lenient surface: a population that is entirely
    /// homozygous alternate (all genotype code 2) carries the alternate allele at
    /// frequency 1.0; the alternate-allele frequency f = 1 makes
    /// MAF = min(1, 0) = 0 — the locus is fixed, so there is no minor allele. The
    /// "all same" boundary must therefore yield MAF = 0, never a value outside
    /// [0, 0.5].
    /// </summary>
    [Test]
    public void Maf_AllSameAlternateAllele_IsZero()
    {
        for (int trial = 0; trial < 16; trial++)
        {
            int n = Rng.Next(1, 50_000);

            double maf = PopulationGeneticsAnalyzer.CalculateMAF(Enumerable.Repeat(2, n));

            maf.Should().BeApproximately(0.0, Tolerance,
                because: $"all {n} genotypes hom-alt → alt frequency 1.0 → MAF = min(1, 0) = 0");
            maf.Should().BeInRange(0.0, 0.5, "MAF ∈ [0, 0.5] (INV-03)");
        }
    }

    #endregion

    #region BE — negative counts (biologically impossible input)

    /// <summary>
    /// BE "−1": a negative genotype count is biologically impossible. The STRICT
    /// surface must reject it with the *documented, intentional*
    /// ArgumentOutOfRangeException — NOT a silent negative/garbage frequency and
    /// NOT a DivideByZero (Allele_Frequency.md §3.1, §3.3;
    /// PopulationGeneticsAnalyzer.cs lines 143–148). Each of the three count
    /// arguments is guarded independently, and the extreme int.MinValue is
    /// covered to pin that no overflow slips a negative past the guard.
    /// </summary>
    [TestCase(-1, 0, 0, TestName = "AlleleFrequencies_NegativeHomozygousMajor_Throws")]
    [TestCase(0, -1, 0, TestName = "AlleleFrequencies_NegativeHeterozygous_Throws")]
    [TestCase(0, 0, -1, TestName = "AlleleFrequencies_NegativeHomozygousMinor_Throws")]
    [TestCase(-5, 3, 2, TestName = "AlleleFrequencies_NegativeMajorAmongPositives_Throws")]
    [TestCase(int.MinValue, 0, 0, TestName = "AlleleFrequencies_IntMinValueMajor_Throws")]
    [TestCase(10, -10, 10, TestName = "AlleleFrequencies_NegativeHetAmongPositives_Throws")]
    public void AlleleFrequencies_NegativeCount_ThrowsDocumentedArgumentOutOfRange(
        int hMajor, int het, int hMinor)
    {
        var act = () => PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(hMajor, het, hMinor);

        act.Should().Throw<ArgumentOutOfRangeException>(
            "a negative genotype count is biologically impossible and is rejected at the validation " +
            "gate, never miscounted into a negative frequency and never crashed");
    }

    /// <summary>
    /// BE "−1", lenient surface: CalculateMAF does NOT validate its genotype codes
    /// (Allele_Frequency.md §3.3, §5.2/§5.3 "caller must normalize"). An
    /// out-of-encoding code (a negative genotype, or a code &gt; 2) is summed
    /// verbatim. The contract guarantee we pin here is the fuzz-safety one: the
    /// method must stay TOTAL — it must NOT throw, must NOT divide by zero (the
    /// list is non-empty, so the denominator 2·count &gt; 0), and must NEVER
    /// produce NaN or ±Infinity. The returned value MAY fall outside [0, 0.5] for
    /// such malformed input — that is the documented consequence of feeding
    /// out-of-encoding codes, not a crash. This makes the strict/lenient split
    /// explicit: the malformed input the strict surface would reject is here
    /// carried as a defined-but-out-of-range number, never a DivideByZero.
    /// </summary>
    [Test]
    public void Maf_OutOfEncodingGenotypeCodes_StaysFiniteAndNeverThrows()
    {
        // A deterministic sweep including negatives and codes > 2, mixed with the
        // valid 0/1/2 encoding. Every list is non-empty (denominator > 0).
        for (int trial = 0; trial < 64; trial++)
        {
            int length = Rng.Next(1, 50);
            var genotypes = new int[length];
            for (int i = 0; i < length; i++)
                genotypes[i] = Rng.Next(-3, 6); // spans -3..-1 (impossible), 0/1/2 (valid), 3..5 (>2)

            double maf = double.NaN;
            var act = () => maf = PopulationGeneticsAnalyzer.CalculateMAF(genotypes);

            act.Should().NotThrow(
                "the lenient MAF surface is total over int genotype codes — it never throws on garbage");
            double.IsNaN(maf).Should().BeFalse(
                "a non-empty list has a non-zero denominator, so the result is never a 0/0 NaN");
            double.IsInfinity(maf).Should().BeFalse(
                "the arithmetic on bounded int sums and a non-zero denominator never yields ±Infinity");
        }
    }

    /// <summary>
    /// BE "−1", explicit defined-result pin for the lenient surface: a single
    /// genotype code of −1 sums to altAlleles = −1 over totalAlleles = 2, giving
    /// altFreq = −0.5 and MAF = min(−0.5, 1.5) = −0.5. We pin this EXACT value to
    /// document that the surface does NOT silently coerce the impossible input
    /// into [0, 0.5] — it carries the malformed code through deterministically
    /// (the "caller must normalize" contract). The point is that the behavior is
    /// *defined and stable*, not that −0.5 is a valid frequency.
    /// </summary>
    [Test]
    public void Maf_SingleNegativeOneCode_ReturnsDefinedOutOfRangeValue_NotACrash()
    {
        double maf = double.NaN;
        var act = () => maf = PopulationGeneticsAnalyzer.CalculateMAF(new[] { -1 });

        act.Should().NotThrow("an out-of-encoding code is carried, not validated, by the lenient surface");
        maf.Should().BeApproximately(-0.5, Tolerance,
            because: "alt copies = -1 over 2 total → altFreq -0.5 → MAF min(-0.5, 1.5) = -0.5; " +
                     "a defined (if biologically meaningless) value, documenting the no-validation contract");
    }

    #endregion

    #endregion
}
