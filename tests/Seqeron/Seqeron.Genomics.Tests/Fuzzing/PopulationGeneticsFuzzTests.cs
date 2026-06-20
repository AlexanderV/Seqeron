using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Population;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the PopGen area — allele frequency (POP-FREQ-001),
/// nucleotide diversity π (POP-DIV-001) and Hardy-Weinberg equilibrium
/// (POP-HW-001).
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
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: POP-DIV-001 — nucleotide diversity π (PopGen)
/// Checklist: docs/checklists/03_FUZZING.md, row 44.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — the degenerate sample boundaries that stress
///           the pairwise-diversity arithmetic: 0 sequences and 1 sequence (the
///           n &lt; 2 boundary — no unordered pair exists, so the number of pairs
///           C(n,2) = n(n−1)/2 is 0 → an unguarded denominator would
///           DivideByZero/NaN), "all sequences identical" (no variation → every
///           pairwise difference d_ij is 0 → π must be exactly 0 — THE key
///           identity), and very long sequences (the O(n²·L) cost boundary — must
///           still complete and never hang).
/// — docs/checklists/03_FUZZING.md §Description ("BE = граничні значення: 0, -1,
///   MaxInt, empty").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The nucleotide-diversity contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Nucleotide diversity π is the average number of nucleotide differences PER
/// SITE between all unordered pairs of aligned sequences:
///     π = ( Σ_{i&lt;j} d_ij ) / ( C(n,2) · L )
/// where d_ij is the number of differing positions between sequences i and j,
/// n is the number of sequences, and L is the (common, aligned) sequence length.
///   — docs/algorithms/Population_Genetics/Diversity_Statistics.md §2.2, §2.4
///     (INV-01: π ≥ 0 — it is a normalized count of differences; INV-02: π = 0
///      for identical sequences — every d_ij is 0). Sources: Nei &amp; Li (1979)
///     [1]; Tajima (1989) [3]; Wikipedia "Tajima's D" [5].
///
/// Entry point under test —
///   PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(
///       IEnumerable&lt;IReadOnlyList&lt;char&gt;&gt;)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs
///    lines 205–234). Documented validation/edge behavior
///   (Diversity_Statistics.md §3.3, §6.1):
///     • n &lt; 2 sequences (0 or 1) → the explicit `seqList.Count &lt; 2` guard
///       returns the DEFINED 0 BEFORE any pair loop; there is NO C(n,2) division
///       (comparisons would be 0), so NO DivideByZeroException and NO NaN. This is
///       the KEY 0-/1-sequence fuzz concern — the C(n,2) = 0 boundary.
///     • all sequences identical → every d_ij = 0 → totalDiff = 0 → π = 0 exactly,
///       independent of how many sequences or how long (INV-02).
///     • two sequences differing at every site → π = 1.0 (each of the L sites is
///       one mismatch in the single pair).
///     • very long sequences (n small, L large) → the O(n²·L) scan completes
///       within a bounded time; pinned under `[CancelAfter]` so a hang is a
///       failure, not an infinite wait.
///     • π is ALWAYS ≥ 0 and finite for valid aligned input (INV-01) — a count of
///       mismatches over a positive denominator can never be negative, NaN, or
///       ±Infinity.
///   The implementation does NOT validate equal length (§5.4): callers must supply
///   aligned, equal-length sequences. For n ≥ 2 the denominator is C(n,2)·L; with
///   a positive common L this is &gt; 0, so the division is well-defined.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: POP-HW-001 — Hardy-Weinberg equilibrium (PopGen)
/// Checklist: docs/checklists/03_FUZZING.md, row 45.
/// Fuzz strategies exercised for THIS unit:
///   • BE  = Boundary Exploitation — the degenerate genotype-count boundaries
///           that stress the allele-frequency arithmetic and the chi-square
///           goodness-of-fit: 0 genotypes (n = 0 → the DivideByZero boundary in
///           BOTH the p/q estimate 1/(2n) AND every expected-count term, whose
///           expected category is 0); a single genotype (the minimal non-empty
///           sample); an all-heterozygous sample (the classic HWE-violation case
///           — observed het excess vs the p=q=0.5 expectation); and an
///           all-homozygous sample (a het DEFICIT, the inbreeding-pattern
///           departure, and the monomorphic fully-fixed sub-case where an
///           expected category is exactly 0).
/// — docs/checklists/03_FUZZING.md §Description ("BE = граничні значення: 0, -1,
///   MaxInt, empty").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The Hardy-Weinberg contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// For a biallelic locus the test estimates allele frequencies from observed
/// genotype counts (n_AA, n_Aa, n_aa), n = n_AA + n_Aa + n_aa:
///     p = (2·n_AA + n_Aa) / (2n),   q = 1 − p
/// computes the HWE-expected counts E_AA = p²n, E_Aa = 2pqn, E_aa = q²n, and a
/// chi-square goodness-of-fit statistic with df = 1:
///     χ² = Σ_i (O_i − E_i)² / E_i   (only over categories with E_i &gt; 0)
/// returning a HardyWeinbergResult(p/q-derived expected counts, χ², p-value,
/// InEquilibrium = PValue ≥ significanceLevel).
///   — docs/algorithms/Population_Genetics/Hardy_Weinberg_Test.md §2.2, §2.4
///     (INV-01: p + q = 1; INV-03: expected counts sum to n; INV-04: χ² ≥ 0;
///      INV-05: p-value ∈ [0, 1]). Sources: Hardy (1908) [1]; Weinberg (1908) [2];
///      Emigh (1980) [3]; Wikipedia "Hardy-Weinberg principle" [5].
///
/// Entry point under test —
///   PopulationGeneticsAnalyzer.TestHardyWeinberg(string, int, int, int, double)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs
///    lines 424–471). Documented validation/edge behavior
///   (Hardy_Weinberg_Test.md §3.3, §5.2, §6.1):
///     • n = 0 (0 genotypes) → the explicit `n == 0` guard returns the DEFINED
///       result (χ² = 0, PValue = 1, InEquilibrium = true) BEFORE any 1/(2n)
///       division; there is NO DivideByZeroException and NO NaN. This is the KEY
///       0-genotypes fuzz concern (Hardy_Weinberg_Test.md §3.3, §6.1).
///     • an expected category equal to 0 (e.g. the monomorphic fully-fixed
///       sample, where E_Aa = E_aa = 0) → that χ² term is SKIPPED by the explicit
///       `expected > 0` guard, so the 0-expected division is never performed —
///       NO DivideByZero, NO NaN (Hardy_Weinberg_Test.md §5.2, §6.1). This is the
///       KEY 0-expected fuzz concern.
///     • a single genotype (e.g. (1,0,0)) → a DEFINED result: monomorphic, so
///       observed equals expected → χ² = 0, InEquilibrium = true (§6.1).
///     • all heterozygous (0, N, 0) → p = q = 0.5, expected (0.25N, 0.5N, 0.25N);
///       the observed het EXCESS gives χ² = N exactly (a meaningful, scaling,
///       nonzero departure — the classic HWE violation, §6.1 / unit test
///       HW-M9 pins N=100 → χ² = 100).
///     • all homozygous mixed (N, 0, N) → p = q = 0.5 again, expected
///       (0.25·2N, 0.5·2N, 0.25·2N); the observed het DEFICIT (inbreeding pattern)
///       gives χ² = 2N exactly (a meaningful nonzero departure).
///     • all homozygous fixed (N, 0, 0) or (0, 0, N) → p = 1 or 0, expected equals
///       observed → χ² = 0, InEquilibrium = true (the monomorphic case where two
///       expected categories are 0 and their χ² terms are skipped, §6.1).
///     • χ² is ALWAYS ≥ 0 and finite for any non-negative count triple (INV-04),
///       and the p-value stays in [0, 1] (INV-05).
///   The method does NOT validate negative genotype counts or significanceLevel
///   range (§5.4 accepted deviation #1); the fuzz targets here are all
///   non-negative, so we pin the defined behavior of the documented boundary.
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

    // ═══════════════════════════════════════════════════════════════════
    //  POP-DIV-001 — nucleotide diversity π : fuzz targets
    //  Surface: CalculateNucleotideDiversity(IEnumerable<IReadOnlyList<char>>)
    // ═══════════════════════════════════════════════════════════════════

    #region POP-DIV-001 — nucleotide diversity

    #region Helpers — sequence builders

    /// <summary>Materializes a string into the IReadOnlyList&lt;char&gt; the diversity API expects.</summary>
    private static IReadOnlyList<char> Seq(string s) => s.ToCharArray();

    /// <summary>A deterministic random aligned ACGT sequence of the given length.</summary>
    private static IReadOnlyList<char> RandomSeq(int length)
    {
        const string alphabet = "ACGT";
        var buffer = new char[length];
        for (int i = 0; i < length; i++)
            buffer[i] = alphabet[Rng.Next(alphabet.Length)];
        return buffer;
    }

    #endregion

    #region Positive sanity — known sample yields the documented π and π ≥ 0

    /// <summary>
    /// Positive control: the worked example pinned in the algorithm doc and the
    /// unit tests (Diversity_Statistics.md §7.1; the Wikipedia Tajima's D dataset)
    /// — five aligned length-20 sequences with 20 total pairwise differences over
    /// C(5,2) = 10 unordered pairs → k̂ = 20/10 = 2.0 → π = k̂ / L = 2.0/20 = 0.1.
    /// The result must be EXACTLY 0.1 (per the documented formula
    /// π = Σ_{i&lt;j} d_ij / (C(n,2)·L)), must be ≥ 0 (INV-01), and must be finite.
    /// This anchors the fuzz battery: before pinning that boundary input does no
    /// harm, confirm a known-good sample gives the textbook-correct π.
    /// </summary>
    [Test]
    public void NucleotideDiversity_KnownSample_MatchesDocumentedPiAndIsNonNegative()
    {
        var sequences = new List<IReadOnlyList<char>>
        {
            Seq("00000000000000000000"), // Y
            Seq("00100000000010000010"), // A
            Seq("00000000000010000010"), // B
            Seq("00000010000000000010"), // C
            Seq("00000010000010000010"), // D
        };

        double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        pi.Should().BeApproximately(0.1, Tolerance,
            because: "20 pairwise differences over C(5,2)=10 pairs → k̂=2.0 → π=2.0/20=0.1 " +
                     "(Diversity_Statistics.md §7.1, the documented formula)");
        pi.Should().BeGreaterThanOrEqualTo(0.0, "π is a normalized count of differences (INV-01)");
        double.IsNaN(pi).Should().BeFalse("a known-good sample never yields NaN");
        double.IsInfinity(pi).Should().BeFalse("a known-good sample never yields ±Infinity");

        // A second independent anchor: two sequences differing at EVERY site → π = 1.0
        // (each of the L sites contributes one mismatch in the single unordered pair).
        double piAllDiff = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(
            new List<IReadOnlyList<char>> { Seq("AAAA"), Seq("CCCC") });
        piAllDiff.Should().BeApproximately(1.0, Tolerance,
            because: "every site differs in the only pair → π = 4/(1·4) = 1.0 (§6.1)");
    }

    #endregion

    #region BE — 0 sequences (the C(n,2) = 0 division boundary)

    /// <summary>
    /// BE / DivideByZero boundary: with 0 sequences the number of unordered pairs
    /// C(0,2) = 0 and there is no length to normalize by. The contract must NOT
    /// divide by zero and must NOT return NaN — the explicit `seqList.Count &lt; 2`
    /// guard returns the DEFINED 0 BEFORE the pair loop ever runs
    /// (Diversity_Statistics.md §3.3, §6.1; PopulationGeneticsAnalyzer.cs lines
    /// 210–211). This is the single most important fuzz concern for this unit: an
    /// unguarded C(n,2)·L denominator here is a 0/0 NaN (or DivideByZero) in
    /// production.
    /// </summary>
    [Test]
    public void NucleotideDiversity_ZeroSequences_ReturnsDefinedZeroWithNoDivideByZero()
    {
        double pi = double.NaN;
        var act = () => pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(
            Array.Empty<IReadOnlyList<char>>());

        act.Should().NotThrow<DivideByZeroException>(
            "the empty-sample boundary short-circuits before the C(n,2)·L denominator is used");
        act.Should().NotThrow("0 sequences is a defined boundary input (returns 0), not an error");

        pi.Should().Be(0.0, "no sequences → no pairwise differences → π is the defined 0 (§6.1)");
        double.IsNaN(pi).Should().BeFalse("the n<2 guard avoids a 0/0 NaN");
    }

    #endregion

    #region BE — 1 sequence (no unordered pair exists → π = 0)

    /// <summary>
    /// BE: a single sequence — the minimal-but-still-degenerate sample. There is no
    /// unordered pair (C(1,2) = 0), so pairwise diversity is undefined for n &lt; 2;
    /// the contract returns the DEFINED 0 (Diversity_Statistics.md §3.3, §6.1).
    /// Verified across a fixed-seed sweep of lengths so the single-sequence boundary
    /// returns 0 regardless of the sequence content or length, with no
    /// DivideByZero and no NaN.
    /// </summary>
    [Test]
    public void NucleotideDiversity_SingleSequence_ReturnsDefinedZeroForAnyContent()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int length = Rng.Next(1, 256);
            var single = new List<IReadOnlyList<char>> { RandomSeq(length) };

            double pi = double.NaN;
            var act = () => pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(single);

            act.Should().NotThrow<DivideByZeroException>(
                "a single sequence has no pair, so the C(n,2) denominator is never reached");
            act.Should().NotThrow("one sequence is a defined boundary input, not an error");
            pi.Should().Be(0.0, $"no unordered pair exists for n=1 (length {length}) → defined π = 0");
            double.IsNaN(pi).Should().BeFalse("the n<2 guard avoids a 0/0 NaN");
        }
    }

    #endregion

    #region BE — all sequences identical (no variation → π = 0, the key identity)

    /// <summary>
    /// BE / KEY identity (INV-02): when every sequence in the sample is identical,
    /// every pairwise difference count d_ij is 0, so π = 0 EXACTLY — there is no
    /// variation to measure (Diversity_Statistics.md §2.4 INV-02, §6.1). Verified
    /// across a fixed-seed sweep of sample sizes and lengths so the identity holds
    /// independent of n and L. This is the diversity analogue of the "all same
    /// allele" boundary: a monomorphic alignment must read as zero diversity, never
    /// a spurious non-zero value, never NaN.
    /// </summary>
    [Test]
    public void NucleotideDiversity_AllIdenticalSequences_IsExactlyZero()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int n = Rng.Next(2, 40);          // at least a pair, so a real division happens
            int length = Rng.Next(1, 200);
            var template = RandomSeq(length);

            // n copies of the SAME content (a fresh array per copy, same characters).
            var sequences = Enumerable.Range(0, n)
                .Select(_ => (IReadOnlyList<char>)template.ToArray())
                .ToList();

            double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

            pi.Should().Be(0.0,
                because: $"every one of the C({n},2) pairs is identical (length {length}) → all d_ij = 0 → π = 0 (INV-02)");
            pi.Should().BeGreaterThanOrEqualTo(0.0, "π ≥ 0 always (INV-01)");
            double.IsNaN(pi).Should().BeFalse("an identical, non-empty sample divides 0 by a positive denominator, not 0/0");
        }
    }

    #endregion

    #region BE — random aligned samples are always non-negative and finite

    /// <summary>
    /// BE / INV-01 sweep: across a fixed-seed battery of random, equal-length
    /// aligned samples (varying n ≥ 2 and L ≥ 1) π must ALWAYS be a well-defined,
    /// finite, non-negative per-site rate — a count of mismatches over a positive
    /// C(n,2)·L denominator can never be negative, NaN, or ±Infinity, and since π
    /// is a per-site average difference rate it can never exceed 1.0. This pins the
    /// total-function contract over the random-but-aligned interior, complementing
    /// the degenerate-boundary tests.
    /// </summary>
    [Test]
    public void NucleotideDiversity_RandomAlignedSamples_AlwaysInUnitIntervalAndFinite()
    {
        for (int trial = 0; trial < 64; trial++)
        {
            int n = Rng.Next(2, 25);
            int length = Rng.Next(1, 120);
            var sequences = Enumerable.Range(0, n)
                .Select(_ => RandomSeq(length))
                .ToList();

            double pi = double.NaN;
            var act = () => pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

            act.Should().NotThrow("a random equal-length aligned sample is always a defined input");
            double.IsNaN(pi).Should().BeFalse("a positive C(n,2)·L denominator never yields a 0/0 NaN");
            double.IsInfinity(pi).Should().BeFalse("bounded mismatch counts over a positive denominator never give ±Infinity");
            pi.Should().BeGreaterThanOrEqualTo(0.0, "π is a normalized count of differences (INV-01)");
            pi.Should().BeLessThanOrEqualTo(1.0, "π is a per-site average mismatch rate, capped at one mismatch per site");
        }
    }

    #endregion

    #region BE — very long sequences (the O(n²·L) cost boundary must complete, not hang)

    /// <summary>
    /// BE / hang boundary: nucleotide diversity is O(n²·L). With a small sample of
    /// VERY LONG aligned sequences the pairwise scan must still COMPLETE within a
    /// bounded time and return a finite, in-range π — it must never hang
    /// (Diversity_Statistics.md §4.3). We keep n small (a handful of sequences) and
    /// L large (hundreds of thousands of sites) so the L dimension is stressed
    /// without an n² blow-up, and pin the run under `[CancelAfter]` so a hang is a
    /// hard test failure rather than an infinite wait. The constructed sample
    /// differs at exactly one known site, giving a tiny but exactly-pinned π that
    /// also proves the long scan is arithmetically correct, not just fast.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void NucleotideDiversity_VeryLongSequences_CompletesWithFinitePi()
    {
        const int length = 500_000;
        var a = new char[length];
        var b = new char[length];
        Array.Fill(a, 'A');
        Array.Fill(b, 'A');
        b[length / 2] = 'C'; // exactly one differing site across the single pair

        var sequences = new List<IReadOnlyList<char>> { a, b };

        double pi = double.NaN;
        var act = () => pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        act.Should().NotThrow("the long-sequence O(n²·L) scan must complete without crashing");
        double.IsNaN(pi).Should().BeFalse("a positive denominator never yields NaN even at scale");
        double.IsInfinity(pi).Should().BeFalse("the bounded difference count never overflows into ±Infinity");
        pi.Should().BeApproximately(1.0 / length, Tolerance,
            because: "exactly one of the L=500000 sites differs in the single pair → π = 1/L");
        pi.Should().BeGreaterThanOrEqualTo(0.0, "π ≥ 0 always (INV-01)");
    }

    /// <summary>
    /// BE / hang boundary, identical-at-scale: a small sample of very long but
    /// IDENTICAL sequences must complete and read as exactly π = 0 — the INV-02
    /// identity must survive the long O(n²·L) scan, proving the zero-diversity
    /// short-circuit is by content (all d_ij = 0), not an accident of small inputs.
    /// Pinned under `[CancelAfter]` so the long scan cannot silently hang.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void NucleotideDiversity_VeryLongIdenticalSequences_CompletesAtExactlyZero()
    {
        const int length = 400_000;
        var template = new char[length];
        Array.Fill(template, 'G');

        var sequences = new List<IReadOnlyList<char>>
        {
            (char[])template.Clone(),
            (char[])template.Clone(),
            (char[])template.Clone(),
        };

        double pi = double.NaN;
        var act = () => pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(sequences);

        act.Should().NotThrow("the long identical-sequence scan must complete without crashing");
        pi.Should().Be(0.0, "every one of the C(3,2)=3 pairs is identical at scale → π = 0 (INV-02)");
        double.IsNaN(pi).Should().BeFalse("an identical sample divides 0 by a positive denominator, not 0/0");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  POP-HW-001 — Hardy-Weinberg equilibrium : fuzz targets
    //  Surface: TestHardyWeinberg(string, int, int, int, double)
    // ═══════════════════════════════════════════════════════════════════

    #region POP-HW-001 — Hardy-Weinberg equilibrium

    #region Positive sanity — known genotype set yields the documented chi-square and HWE verdict

    /// <summary>
    /// Positive control: Ford's scarlet tiger moth dataset pinned in the algorithm
    /// doc and the unit tests (Hardy_Weinberg_Test.md §7.1) — 1469 AA, 138 Aa,
    /// 5 aa, n = 1612 → expected counts ≈ (1467.40, 141.21, 3.40) and χ² ≈ 0.8309,
    /// which is below the df=1 α=0.05 critical value 3.841, so the sample is NOT
    /// rejected (InEquilibrium = true). The result must match the documented χ²,
    /// the expected counts must sum to n (INV-03), p + q must equal 1 (INV-01),
    /// χ² must be ≥ 0 (INV-04) and the p-value must lie in [0, 1] (INV-05). This
    /// anchors the fuzz battery: before pinning that boundary input does no harm,
    /// confirm a known-good genotype set gives the textbook-correct verdict.
    /// </summary>
    [Test]
    public void HardyWeinberg_KnownGenotypeSet_MatchesDocumentedChiSquareAndVerdict()
    {
        var r = PopulationGeneticsAnalyzer.TestHardyWeinberg(
            "rs-moth", observedAA: 1469, observedAa: 138, observedaa: 5);

        r.ChiSquare.Should().BeApproximately(0.8309, 0.01,
            because: "Ford's moth data → χ² ≈ 0.8309 (Hardy_Weinberg_Test.md §7.1)");
        r.InEquilibrium.Should().BeTrue("χ² ≈ 0.83 < 3.841 → not rejected at α=0.05 (§7.1)");

        r.ExpectedAA.Should().BeApproximately(1467.40, 0.05, "E_AA = p²n (§7.1)");
        r.ExpectedAa.Should().BeApproximately(141.21, 0.05, "E_Aa = 2pqn (§7.1)");
        r.Expectedaa.Should().BeApproximately(3.40, 0.05, "E_aa = q²n (§7.1)");

        (r.ExpectedAA + r.ExpectedAa + r.Expectedaa).Should().BeApproximately(1612.0, Tolerance,
            because: "expected counts sum to n (INV-03)");

        double p = (2.0 * 1469 + 138) / (2.0 * 1612);
        double q = 1 - p;
        (p + q).Should().BeApproximately(1.0, Tolerance, "p + q = 1 (INV-01)");

        r.ChiSquare.Should().BeGreaterThanOrEqualTo(0.0, "χ² ≥ 0 (INV-04)");
        r.PValue.Should().BeInRange(0.0, 1.0, "the chi-square tail probability is a probability (INV-05)");
        double.IsNaN(r.ChiSquare).Should().BeFalse("a known-good sample never yields NaN");
    }

    #endregion

    #region BE — 0 genotypes (the DivideByZero boundary: n = 0)

    /// <summary>
    /// BE / DivideByZero boundary: with 0 genotypes n = 0, so BOTH the allele-
    /// frequency estimate p = (…)/(2n) AND every expected-count term would divide
    /// by zero. The contract must NOT divide by zero and must NOT return NaN — the
    /// explicit `n == 0` guard returns the DEFINED result (χ² = 0, PValue = 1,
    /// InEquilibrium = true) BEFORE any division (Hardy_Weinberg_Test.md §3.3,
    /// §6.1; PopulationGeneticsAnalyzer.cs lines 433–436). This is the single most
    /// important fuzz concern for this unit: an unguarded 1/(2n) here is a
    /// DivideByZero/NaN in production.
    /// </summary>
    [Test]
    public void HardyWeinberg_ZeroGenotypes_ReturnsDefinedEquilibriumWithNoDivideByZero()
    {
        PopulationGeneticsAnalyzer.HardyWeinbergResult r = default;
        var act = () => r = PopulationGeneticsAnalyzer.TestHardyWeinberg("rs-empty", 0, 0, 0);

        act.Should().NotThrow<DivideByZeroException>(
            "the n=0 boundary is guarded; neither 1/(2n) nor an expected-count term is divided by zero");
        act.Should().NotThrow("0 genotypes is a defined boundary input (no evidence against the null), not an error");

        r.ChiSquare.Should().Be(0.0, "no data → no deviation → χ² is the defined 0 (§6.1)");
        r.PValue.Should().Be(1.0, "no data → no evidence against equilibrium → PValue = 1 (§6.1)");
        r.InEquilibrium.Should().BeTrue("no data is treated as no evidence against the null model (§6.1)");
        double.IsNaN(r.ChiSquare).Should().BeFalse("the n=0 guard avoids a 0/0 NaN");
        double.IsNaN(r.PValue).Should().BeFalse("the n=0 guard avoids a 0/0 NaN");
    }

    #endregion

    #region BE — single genotype (the minimal non-empty sample)

    /// <summary>
    /// BE: the minimal non-empty sample — a single individual. Each is monomorphic
    /// for one genotype, so observed equals expected and χ² = 0, InEquilibrium =
    /// true (Hardy_Weinberg_Test.md §6.1 "Single monomorphic sample (1,0,0)").
    /// For the single homozygotes (1,0,0)/(0,0,1) the locus is fixed → two expected
    /// categories are exactly 0 and their χ² terms are skipped (the 0-expected
    /// boundary); for the single heterozygote (0,1,0) p = q = 0.5 but with a single
    /// observation the het-only sample is the all-het case at N=1 → χ² = 1, but
    /// χ²=1 on df=1 has a p-value ≈ 0.317 ≥ 0.05, so a single heterozygote is too
    /// small a sample to be a significant departure → InEquilibrium = true. We pin
    /// each defined result and confirm no DivideByZero and no NaN.
    /// </summary>
    [TestCase(1, 0, 0, 0.0, true, TestName = "HardyWeinberg_SingleHomozygousMajor_IsDefinedEquilibrium")]
    [TestCase(0, 0, 1, 0.0, true, TestName = "HardyWeinberg_SingleHomozygousMinor_IsDefinedEquilibrium")]
    [TestCase(0, 1, 0, 1.0, true, TestName = "HardyWeinberg_SingleHeterozygote_IsDefinedHetExcess")]
    public void HardyWeinberg_SingleGenotype_IsDefinedResult(
        int aa, int ab, int bb, double expectedChi, bool expectedInEq)
    {
        PopulationGeneticsAnalyzer.HardyWeinbergResult r = default;
        var act = () => r = PopulationGeneticsAnalyzer.TestHardyWeinberg("rs-single", aa, ab, bb);

        act.Should().NotThrow("a single genotype is a defined minimal sample, not an error");
        r.ChiSquare.Should().BeApproximately(expectedChi, Tolerance,
            because: "a single genotype gives an exact, documented χ² (§6.1)");
        r.InEquilibrium.Should().Be(expectedInEq, "the single-genotype verdict is defined and stable (§6.1)");

        r.ChiSquare.Should().BeGreaterThanOrEqualTo(0.0, "χ² ≥ 0 (INV-04)");
        double.IsNaN(r.ChiSquare).Should().BeFalse("a defined single-genotype sample never yields NaN");
        r.PValue.Should().BeInRange(0.0, 1.0, "the p-value is a probability (INV-05)");

        (r.ExpectedAA + r.ExpectedAa + r.Expectedaa).Should().BeApproximately(1.0, Tolerance,
            because: "the single expected counts sum to n = 1 (INV-03)");
    }

    #endregion

    #region BE — all heterozygous (the classic HWE-violation case: het excess → χ² = N)

    /// <summary>
    /// BE / KEY HWE-violation case: an all-heterozygous sample (0, N, 0) gives
    /// p = q = 0.5, so the expected counts are (0.25N, 0.5N, 0.25N) — but the
    /// observed are (0, N, 0), a pure het EXCESS. The chi-square then expands to
    /// 0.25N + 0.5N + 0.25N = N EXACTLY, a meaningful, nonzero, scaling departure
    /// from equilibrium (Hardy_Weinberg_Test.md §6.1; unit test HW-M9 pins N=100 →
    /// χ² = 100). Verified across a fixed-seed sweep of N so χ² = N holds at scale,
    /// the departure is always significant (InEquilibrium = false for N large
    /// enough to clear the 3.841 critical value), expected counts sum to n
    /// (INV-03), and χ² stays finite and ≥ 0 (INV-04) — never NaN.
    /// </summary>
    [Test]
    public void HardyWeinberg_AllHeterozygous_HasHetExcessChiSquareEqualToN()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int n = Rng.Next(4, 100_000); // ≥ 4 so χ² = N > 3.841 → a real, significant departure

            var r = PopulationGeneticsAnalyzer.TestHardyWeinberg("rs-allhet", 0, n, 0);

            r.ChiSquare.Should().BeApproximately(n, n * 1e-9 + Tolerance,
                because: $"all {n} heterozygous → p=q=0.5, het excess → χ² = N = {n} exactly (§6.1)");
            r.ChiSquare.Should().BeGreaterThan(3.841,
                "the het excess is a significant departure from HWE at α=0.05 (df=1 critical value 3.841)");
            r.InEquilibrium.Should().BeFalse(
                "a pure het-excess sample is rejected as out of equilibrium — the classic HWE violation (§6.1)");

            r.ExpectedAA.Should().BeApproximately(0.25 * n, Tolerance, "E_AA = p²n = 0.25N");
            r.ExpectedAa.Should().BeApproximately(0.5 * n, Tolerance, "E_Aa = 2pqn = 0.5N");
            r.Expectedaa.Should().BeApproximately(0.25 * n, Tolerance, "E_aa = q²n = 0.25N");
            (r.ExpectedAA + r.ExpectedAa + r.Expectedaa).Should().BeApproximately(n, Tolerance,
                because: "expected counts sum to n (INV-03)");

            r.ChiSquare.Should().BeGreaterThanOrEqualTo(0.0, "χ² ≥ 0 (INV-04)");
            double.IsNaN(r.ChiSquare).Should().BeFalse("a non-empty sample never yields a 0/0 NaN");
            double.IsInfinity(r.ChiSquare).Should().BeFalse("the bounded arithmetic never yields ±Infinity");
            r.PValue.Should().BeInRange(0.0, 1.0, "the p-value is a probability (INV-05)");
        }
    }

    #endregion

    #region BE — all homozygous (het deficit, and the monomorphic 0-expected boundary)

    /// <summary>
    /// BE "all homozygous", het-DEFICIT case: a sample split between the two
    /// homozygotes with NO heterozygotes (N, 0, N) gives p = q = 0.5 over a sample
    /// of 2N, so the expected counts are (0.5N, N, 0.5N) — but the observed are
    /// (N, 0, N), a pure het DEFICIT (the inbreeding pattern). The chi-square
    /// expands to 0.5N + N + 0.5N = 2N EXACTLY, a meaningful, nonzero departure
    /// from equilibrium. Verified across a fixed-seed sweep of N so χ² = 2N holds
    /// at scale, the departure is significant, expected counts sum to 2N (INV-03),
    /// and χ² stays finite and ≥ 0 — never NaN.
    /// </summary>
    [Test]
    public void HardyWeinberg_AllHomozygousMixed_HasHetDeficitChiSquareEqualToTwoN()
    {
        for (int trial = 0; trial < 32; trial++)
        {
            int n = Rng.Next(2, 50_000); // χ² = 2N ≥ 4 > 3.841 → a real, significant departure

            var r = PopulationGeneticsAnalyzer.TestHardyWeinberg("rs-homdef", n, 0, n);

            r.ChiSquare.Should().BeApproximately(2.0 * n, n * 1e-9 + Tolerance,
                because: $"(N,0,N) → p=q=0.5 over 2N, het deficit → χ² = 2N = {2 * n} exactly");
            r.ChiSquare.Should().BeGreaterThan(3.841,
                "the het deficit is a significant departure from HWE at α=0.05");
            r.InEquilibrium.Should().BeFalse(
                "a pure het-deficit (inbreeding) sample is rejected as out of equilibrium");

            (r.ExpectedAA + r.ExpectedAa + r.Expectedaa).Should().BeApproximately(2.0 * n, Tolerance,
                because: "expected counts sum to n = 2N (INV-03)");
            r.ChiSquare.Should().BeGreaterThanOrEqualTo(0.0, "χ² ≥ 0 (INV-04)");
            double.IsNaN(r.ChiSquare).Should().BeFalse("a non-empty sample never yields a 0/0 NaN");
            double.IsInfinity(r.ChiSquare).Should().BeFalse("the bounded arithmetic never yields ±Infinity");
            r.PValue.Should().BeInRange(0.0, 1.0, "the p-value is a probability (INV-05)");
        }
    }

    /// <summary>
    /// BE "all homozygous", MONOMORPHIC fully-fixed case — the KEY 0-expected
    /// boundary: when the whole sample is one homozygote (N, 0, 0) or (0, 0, N) the
    /// locus is FIXED (p = 1 or p = 0), so two of the three expected categories are
    /// EXACTLY 0. Their χ² terms (which would divide by the 0 expected count) are
    /// SKIPPED by the explicit `expected > 0` guard, so there is NO DivideByZero
    /// and NO NaN, and because observed equals expected the statistic is χ² = 0 →
    /// InEquilibrium = true (Hardy_Weinberg_Test.md §5.2, §6.1;
    /// PopulationGeneticsAnalyzer.cs lines 450–455). This pins the 0-expected
    /// division boundary the prompt flags as the chi-square danger spot.
    /// </summary>
    [Test]
    public void HardyWeinberg_MonomorphicFixed_SkipsZeroExpectedTermWithNoDivideByZero()
    {
        for (int trial = 0; trial < 16; trial++)
        {
            int n = Rng.Next(1, 100_000);

            foreach (var (aa, ab, bb, label) in new[]
                     {
                         (n, 0, 0, "fixed major"),
                         (0, 0, n, "fixed minor"),
                     })
            {
                PopulationGeneticsAnalyzer.HardyWeinbergResult r = default;
                var act = () => r = PopulationGeneticsAnalyzer.TestHardyWeinberg("rs-fixed", aa, ab, bb);

                act.Should().NotThrow<DivideByZeroException>(
                    $"the 0-expected categories of a {label} locus are skipped, not divided by");
                act.Should().NotThrow($"a {label} monomorphic sample of {n} is a defined input");

                r.ChiSquare.Should().BeApproximately(0.0, Tolerance,
                    because: $"a {label} locus has observed == expected → χ² = 0 (§6.1)");
                r.InEquilibrium.Should().BeTrue($"the {label} monomorphic counts match expectation exactly (§6.1)");

                double.IsNaN(r.ChiSquare).Should().BeFalse(
                    "the expected>0 guard skips the 0-expected term, so no 0/0 NaN arises");
                double.IsInfinity(r.ChiSquare).Should().BeFalse(
                    "skipping the 0-expected term keeps the statistic finite");
                r.PValue.Should().BeInRange(0.0, 1.0, "the p-value is a probability (INV-05)");
            }
        }
    }

    #endregion

    #endregion
}
