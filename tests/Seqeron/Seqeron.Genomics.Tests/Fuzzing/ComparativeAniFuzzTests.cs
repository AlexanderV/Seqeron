using System;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Comparative-Genomics area — Average Nucleotide Identity
/// (COMPGEN-ANI-001), the ANIb genome-to-genome relatedness estimator.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output (NaN / ±Infinity / a value outside
/// the documented range), and no *unhandled* runtime exception
/// (DivideByZeroException on zero fragments, IndexOutOfRangeException when the
/// sequence is shorter than a fragment, OverflowException). Every input must resolve
/// to EITHER a well-defined, theory-correct result, OR a *documented, intentional*
/// validation exception (ArgumentOutOfRangeException). A raw runtime exception, a
/// hang, a NaN, or an ANI outside [0, 1] is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: COMPGEN-ANI-001 — Average Nucleotide Identity (ANIb)
/// Checklist: docs/checklists/03_FUZZING.md, row 131.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate genome boundaries called out in
///          the checklist row: identical genomes, no shared k-mers (no conserved
///          fragment), the empty genome, and the single-base genome.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The ANI contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// ANI (Goris et al. 2007, ANIb) measures whole-genome relatedness as the mean
/// nucleotide identity of the conserved (alignable) fragments. The query genome is
/// cut into consecutive, NON-overlapping fragments of length L (default 1020); each
/// fragment is placed (ungapped, full-length) against the reference; the per-fragment
/// identity is recalculated over the WHOLE fragment length:
///   id(fᵢ) = (matching bases of fᵢ at its best placement) / L
///   cov(fᵢ) = alignable length / L   (1.0 for any full-length placement, else 0)
/// A fragment qualifies iff id(fᵢ) &gt; minIdentity (0.30) AND cov(fᵢ) ≥
/// minAlignableFraction (0.70). Then
///   ANI(Q, R) = mean over qualifying fᵢ of id(fᵢ)
/// (Average_Nucleotide_Identity.md §2.2; the cut-off clause is the verbatim Goris
/// rule §2.3). The API entry under test is
///   ComparativeGenomics.CalculateANI(string query, string reference,
///       int fragmentLength = 1020, double minIdentity = 0.30,
///       double minAlignableFraction = 0.70)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs
///    lines 1069–1097), with the per-fragment best placement in the private
///   BestUngappedFragmentMatch (lines 1108–1135).
///
/// THE DOCUMENTED INVARIANTS (Average_Nucleotide_Identity.md §2.4):
///   • INV-01: 0 ≤ ANI ≤ 1 — every id(fᵢ) ∈ [0,1] and the mean preserves the bound.
///   • INV-02: identical genomes → ANI = 1.0 — every fragment is a perfect substring.
///   • INV-03: only fragments with id &gt; 0.30 AND cov ≥ 0.70 contribute (Goris cut-off).
///   • INV-04: fragmentation is consecutive, non-overlapping; a trailing partial
///             fragment (&lt; L) is dropped.
///   • INV-05: no qualifying fragment / empty / null input → 0 (mean over empty set).
/// Every positive-result test pins the documented metric (INV-01 range, INV-02
/// self-identity, the exact averaging formula on a hand-built example); this is the
/// load-bearing correctness check that distinguishes a true ANI from a miscount.
///
/// Documented parameter / validation contract (Average_Nucleotide_Identity.md §3.3,
/// §6.1; ComparativeGenomics.cs lines 1076–1079):
///   • Null or empty query/reference → returns 0 (NOT an exception) — the
///     `string.IsNullOrEmpty` guard runs FIRST (lines 1076–1077), so the empty
///     short-circuit wins even when fragmentLength is itself degenerate.
///   • fragmentLength ≤ 0 → ArgumentOutOfRangeException (lines 1078–1079).
///   • query shorter than fragmentLength → 0 (no full fragment fits; the
///     `start + fragmentLength <= query.Length` loop bound never admits a window).
///   • reference shorter than fragmentLength → 0 (no fragment can align ≥ 70 %, so
///     cov = 0 fails the alignable cut-off; INV-03).
/// Input is uppercased (ToUpperInvariant) before comparison (§3.3, case-insensitive).
///
/// The four BE checklist targets map to these documented behaviours:
///   • identical genomes → ANI = 1.0 exactly (INV-02): the canonical self-identity.
///   • no shared k-mers  → ANI = 0 (INV-05): no fragment reaches id &gt; 0.30, so the
///                          qualifying set is empty; the mean-over-empty-set is reported
///                          as 0. There is NO log transform in ANIb, so no log(0)/−∞
///                          risk — the only zero-set hazard is the DivideByZero of an
///                          empty average, which the `Count > 0 ? … : 0` guard prevents.
///   • empty genome      → 0 (the IsNullOrEmpty guard; no DivideByZero on zero fragments).
///   • single base       → 0 (a 1-base genome is shorter than any fragment ≥ 1 unless
///                          fragmentLength = 1; we probe both the shorter-than-L case →
///                          0 and the fragmentLength = 1 case → a well-formed identity,
///                          with NO IndexOutOfRange on the length-1 span extraction).
/// A positive-sanity test pins the documented worked example (§7.1): identical
/// genomes → 1.0; one substitution in a 4-base fragment → (1+1+1+0.75)/4 = 0.9375;
/// and an unrelated genome with no conserved fragment → 0.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ComparativeAniFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// A well-formed ANI result is a finite number inside the documented closed range
    /// [0, 1] (Average_Nucleotide_Identity.md §2.4 INV-01, §3.2). NaN, ±Infinity, a
    /// negative value or a value above 1 is undisciplined output, not a passing test.
    /// </summary>
    private static void AssertWellFormedAni(double ani)
    {
        double.IsNaN(ani).Should().BeFalse("ANI must never be NaN (no 0/0 division of an empty fragment set)");
        double.IsInfinity(ani).Should().BeFalse("ANI must never be ±Infinity (ANIb has no log/division transform that can diverge)");
        ani.Should().BeInRange(0.0, 1.0, "ANI is a mean of per-fragment identities, each in [0,1], so the mean stays in [0,1] (INV-01)");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  COMPGEN-ANI-001 — Average Nucleotide Identity : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region COMPGEN-ANI-001 — Average Nucleotide Identity

    #region BE — Boundary: identical genomes (self-identity ceiling)

    /// <summary>
    /// BE: identical genomes are the upper boundary — the canonical self-identity. Every
    /// query fragment is a perfect substring of the reference, so id(fᵢ) = 1 for each and
    /// the mean is exactly 1.0 (INV-02, Average_Nucleotide_Identity.md §2.4). We pin
    /// EXACT equality (not merely ≈ 1): the ceiling must not drift to 0.999 or saturate
    /// above 1. The Min(..., MaxIdentity) clamp (ComparativeGenomics.cs line 1133)
    /// guarantees it never exceeds 1.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CalculateANI_IdenticalGenomes_IsExactlyOne()
    {
        string genome = RandomDna(600, seed: 1311);

        double ani = ComparativeGenomics.CalculateANI(genome, genome, fragmentLength: 100);

        AssertWellFormedAni(ani);
        ani.Should().Be(1.0, "every fragment of a genome is a perfect substring of itself, so each id = 1 and the mean is exactly 1 (INV-02)");
    }

    /// <summary>
    /// BE: identical genomes self-identity must hold at the off-by-one fragmentation edge
    /// where the genome length is an exact multiple of L (no trailing partial) AND where
    /// it is not (a trailing partial fragment is dropped, INV-04). In both cases every
    /// FULL fragment is perfect, so ANI = 1.0 — the dropped tail cannot pull it down.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CalculateANI_IdenticalGenomes_OneRegardlessOfTrailingPartialFragment()
    {
        // Length 200 = exact 2 × 100 fragments (no tail).
        string exact = RandomDna(200, seed: 4242);
        ComparativeGenomics.CalculateANI(exact, exact, fragmentLength: 100)
            .Should().Be(1.0, "two full perfect fragments with no trailing partial → mean 1.0");

        // Length 250 = 2 × 100 full fragments + a dropped 50-base tail (INV-04).
        string withTail = RandomDna(250, seed: 4243);
        ComparativeGenomics.CalculateANI(withTail, withTail, fragmentLength: 100)
            .Should().Be(1.0, "the trailing partial fragment is dropped, so only the two perfect full fragments are averaged → 1.0 (INV-04)");
    }

    #endregion

    #region BE — Boundary: no shared k-mers (no conserved fragment → ANI floor)

    /// <summary>
    /// BE: two genomes with NO conserved fragment are the lower boundary. We build a
    /// poly-A query and a poly-C reference: every position mismatches, so each fragment's
    /// best identity is 0 &lt; minIdentity (0.30), no fragment qualifies, the qualifying set
    /// is empty, and ANI is reported as 0 (INV-05). The hazard this target probes is the
    /// empty-average DivideByZero — averaging an empty list would be 0/0 = NaN; the
    /// `Count > 0 ? Average() : 0` guard (ComparativeGenomics.cs line 1096) returns a
    /// clean 0 instead. ANIb has NO log transform, so there is no log(0) = −∞ hazard
    /// here. We pin exactly 0 and a finite, in-range result.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CalculateANI_NoSharedKmers_IsZeroNoDivideByZeroNoNaN()
    {
        string query = new string('A', 500);
        string reference = new string('C', 500);

        double ani = ComparativeGenomics.CalculateANI(query, reference, fragmentLength: 100);

        AssertWellFormedAni(ani);
        ani.Should().Be(0.0, "every fragment mismatches at every base → id = 0 < 0.30 → no qualifying fragment → mean-over-empty-set reported as 0 (INV-05); the Count>0 guard prevents the 0/0 NaN");
    }

    /// <summary>
    /// BE: a borderline divergence just BELOW the 30 % identity cut-off must still floor
    /// to 0 — pinning that INV-03 is a strict `id &gt; 0.30`, not `id ≥ 0.30 - ε`. We make a
    /// query whose every fragment matches the reference in only ~25 % of positions, so no
    /// fragment qualifies and ANI = 0. This guards the cut-off from silently admitting
    /// sub-threshold noise.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CalculateANI_AllFragmentsBelowIdentityCutoff_IsZero()
    {
        // Reference is all 'A'. Query fragment "ACCC" repeated: best placement matches
        // only the single 'A' per 4 bases → id = 0.25 < 0.30 (cut-off), so excluded.
        string reference = new string('A', 400);
        string query = string.Concat(System.Linq.Enumerable.Repeat("ACCC", 100)); // 400 bases

        double ani = ComparativeGenomics.CalculateANI(query, reference, fragmentLength: 4);

        AssertWellFormedAni(ani);
        ani.Should().Be(0.0, "every fragment scores id = 0.25, which is NOT > 0.30, so no fragment qualifies (INV-03) → 0");
    }

    #endregion

    #region BE — Boundary: empty genome (validation guard, zero-fragment safety)

    /// <summary>
    /// BE: the empty genome is the degenerate floor of the size axis. The
    /// `string.IsNullOrEmpty` guard (ComparativeGenomics.cs lines 1076–1077) runs FIRST
    /// and returns 0 — there is nothing to fragment, so there can be no DivideByZero on a
    /// zero-fragment average and no out-of-range span. We probe empty/null on EITHER side
    /// and the empty-vs-empty pair, all of which must short-circuit to 0 (INV-05, §6.1).
    /// </summary>
    [Test]
    public void CalculateANI_EmptyOrNullGenome_IsZeroNoThrow()
    {
        string real = RandomDna(200, seed: 99);

        FluentActions.Invoking(() => ComparativeGenomics.CalculateANI("", real, fragmentLength: 100))
            .Should().NotThrow("an empty query short-circuits to 0 before any fragmentation");
        ComparativeGenomics.CalculateANI("", real, fragmentLength: 100)
            .Should().Be(0.0, "empty query → nothing to fragment → 0 (§6.1)");

        ComparativeGenomics.CalculateANI(real, "", fragmentLength: 100)
            .Should().Be(0.0, "empty reference → no fragment can align → 0 (§6.1)");

        ComparativeGenomics.CalculateANI(null!, real, fragmentLength: 100)
            .Should().Be(0.0, "null query is handled by the IsNullOrEmpty guard, not an NRE (§6.1)");
        ComparativeGenomics.CalculateANI(real, null!, fragmentLength: 100)
            .Should().Be(0.0, "null reference is handled by the IsNullOrEmpty guard, not an NRE (§6.1)");

        ComparativeGenomics.CalculateANI("", "", fragmentLength: 100)
            .Should().Be(0.0, "empty vs empty → 0, never a 0/0 NaN");
    }

    /// <summary>
    /// BE: a degenerate fragmentLength ≤ 0 is the boundary of the fragment-size axis and
    /// is the ONE input the contract rejects with an exception (§3.3,
    /// ComparativeGenomics.cs lines 1078–1079) — a non-positive fragment length is
    /// meaningless. We pin that it throws ArgumentOutOfRangeException carrying the
    /// documented "fragmentLength" parameter name, for both 0 and a negative value, so the
    /// floor cannot drift into defining a zero/negative-length fragment (which would make
    /// the loop bound `start + 0 <= length` run forever — a hang).
    /// </summary>
    [Test]
    public void CalculateANI_NonPositiveFragmentLength_ThrowsArgumentOutOfRange()
    {
        string real = RandomDna(50, seed: 7);

        FluentActions.Invoking(() => ComparativeGenomics.CalculateANI(real, real, fragmentLength: 0))
            .Should().Throw<ArgumentOutOfRangeException>("a fragment length of 0 is meaningless and would make the fragmentation loop never advance")
            .Which.ParamName.Should().Be("fragmentLength");

        FluentActions.Invoking(() => ComparativeGenomics.CalculateANI(real, real, fragmentLength: -10))
            .Should().Throw<ArgumentOutOfRangeException>("a negative fragment length is nonsensical; the contract rejects all fragmentLength <= 0")
            .Which.ParamName.Should().Be("fragmentLength");
    }

    /// <summary>
    /// BE: the empty/null short-circuit must WIN over a degenerate fragmentLength. The
    /// IsNullOrEmpty guard runs before the fragmentLength validation (ComparativeGenomics.cs
    /// lines 1076 then 1078), so an empty genome returns 0 even when paired with a
    /// non-positive fragmentLength — it does NOT throw. This pins the documented guard
    /// ORDER so the two boundary cases cannot be reordered into a spurious throw.
    /// </summary>
    [Test]
    public void CalculateANI_EmptyGenome_WinsOverDegenerateFragmentLength()
    {
        ComparativeGenomics.CalculateANI("", "", fragmentLength: 0)
            .Should().Be(0.0, "the empty short-circuit precedes fragmentLength validation, so this returns 0, not a throw");
        ComparativeGenomics.CalculateANI("", "ACGT", fragmentLength: -5)
            .Should().Be(0.0, "empty query wins over the negative fragmentLength (guard order, §3.3)");
    }

    #endregion

    #region BE — Boundary: single base (sequence shorter than the fragment)

    /// <summary>
    /// BE: a single-base genome is the smallest non-empty input and the classic
    /// shorter-than-fragment boundary. With the default fragment length (1020) no
    /// length-1020 window fits a 1-base query, so the fragmentation loop never runs and
    /// ANI = 0 (no qualifying fragment, INV-05; §6.1 "query shorter than fragmentLength →
    /// 0"). The hazard this probes is an IndexOutOfRange / negative-length span on
    /// `query.AsSpan(start, fragmentLength)` — pinning 0 confirms the loop bound
    /// `start + fragmentLength <= query.Length` excludes the window so no over-long span
    /// is ever taken.
    /// </summary>
    [Test]
    public void CalculateANI_SingleBaseGenome_ShorterThanFragment_IsZeroNoIndexOutOfRange()
    {
        FluentActions.Invoking(() => ComparativeGenomics.CalculateANI("A", "ACGTACGT"))
            .Should().NotThrow("a 1-base query shorter than the default 1020 fragment must not take an over-long span");

        double ani = ComparativeGenomics.CalculateANI("A", "ACGTACGT");
        AssertWellFormedAni(ani);
        ani.Should().Be(0.0, "no length-1020 window fits a 1-base query → no fragment → 0 (§6.1)");

        // Single base on BOTH sides, default fragment length → still no fitting window.
        ComparativeGenomics.CalculateANI("A", "A").Should().Be(0.0, "both 1-base, far shorter than the default fragment → 0");
    }

    /// <summary>
    /// BE: at fragmentLength = 1 a single base IS a full fragment, exercising the
    /// length-1 span extraction directly. A matching single base scores id = 1 (the lone
    /// fragment qualifies, cov = 1.0) → ANI = 1.0; a mismatching single base scores
    /// id = 0 &lt; 0.30 → no qualifying fragment → 0. Both must be finite and in range with
    /// no IndexOutOfRange on the minimal span. This pins the lower fragment-size edge.
    /// </summary>
    [Test]
    public void CalculateANI_SingleBase_FragmentLengthOne_IsWellFormed()
    {
        double match = ComparativeGenomics.CalculateANI("A", "A", fragmentLength: 1);
        AssertWellFormedAni(match);
        match.Should().Be(1.0, "a single matching base at fragmentLength 1 is one perfect fragment → ANI 1.0");

        double mismatch = ComparativeGenomics.CalculateANI("A", "C", fragmentLength: 1);
        AssertWellFormedAni(mismatch);
        mismatch.Should().Be(0.0, "a single mismatching base scores id = 0, which is not > 0.30 → no qualifying fragment → 0 (INV-03)");

        // A single base found inside a longer reference still qualifies (best placement matches).
        ComparativeGenomics.CalculateANI("G", "ACGT", fragmentLength: 1)
            .Should().Be(1.0, "the lone 'G' fragment finds a perfect placement in the reference → id = 1 → ANI 1.0");
    }

    #endregion

    #region Positive sanity — the documented metric on hand-built examples

    /// <summary>
    /// Positive sanity (Average_Nucleotide_Identity.md §7.1 worked example): pins the EXACT
    /// averaging formula, not just "green". Identical 16-base genomes at fragmentLength 4 →
    /// four perfect fragments → ANI 1.0. One substitution in the last fragment ("TTTT" →
    /// "TTTA") makes its identity 3/4 = 0.75 (> 0.30, qualifies), so ANI =
    /// (1 + 1 + 1 + 0.75) / 4 = 0.9375 EXACTLY. This is the load-bearing check that the
    /// per-fragment recalculated identity and the mean are computed correctly.
    /// </summary>
    [Test]
    public void CalculateANI_WorkedExample_MatchesDocumentedValuesExactly()
    {
        double identical = ComparativeGenomics.CalculateANI("AAAACCCCGGGGTTTT", "AAAACCCCGGGGTTTT", fragmentLength: 4);
        AssertWellFormedAni(identical);
        identical.Should().Be(1.0, "four perfect fragments → mean identity 1.0 (§7.1)");

        double oneSubstitution = ComparativeGenomics.CalculateANI("AAAACCCCGGGGTTTA", "AAAACCCCGGGGTTTT", fragmentLength: 4);
        AssertWellFormedAni(oneSubstitution);
        oneSubstitution.Should().BeApproximately(0.9375, 1e-12,
            "three perfect fragments (id 1) + one with a single substitution (id 0.75) → (1+1+1+0.75)/4 = 0.9375 (§7.1)");
    }

    /// <summary>
    /// Positive sanity — the documented similarity gradient. A closely-related genome
    /// (a few scattered substitutions) yields a HIGH ANI strictly below 1.0; an unrelated
    /// genome (no conserved fragment) yields 0; the identical genome yields exactly 1.0.
    /// This pins that ANI orders by divergence as the theory requires (more mutations →
    /// lower ANI), the core business meaning of the metric.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void CalculateANI_OrdersByDivergence_HighForRelatedZeroForUnrelated()
    {
        string reference = RandomDna(400, seed: 20260620);
        char[] mutated = reference.ToCharArray();
        // Introduce a handful of substitutions (~3 %) spread across the genome.
        var rng = new Random(555);
        for (int i = 0; i < 12; i++)
        {
            int pos = rng.Next(reference.Length);
            char cur = mutated[pos];
            char repl = cur == 'A' ? 'C' : 'A';
            mutated[pos] = repl;
        }
        string related = new string(mutated);

        double aniIdentical = ComparativeGenomics.CalculateANI(reference, reference, fragmentLength: 100);
        double aniRelated = ComparativeGenomics.CalculateANI(related, reference, fragmentLength: 100);
        double aniUnrelated = ComparativeGenomics.CalculateANI(new string('A', 400), new string('C', 400), fragmentLength: 100);

        AssertWellFormedAni(aniIdentical);
        AssertWellFormedAni(aniRelated);
        AssertWellFormedAni(aniUnrelated);

        aniIdentical.Should().Be(1.0, "identical → exactly 1.0 (INV-02)");
        aniRelated.Should().BeGreaterThan(0.90, "a genome with only ~3 % divergence is highly related → high ANI");
        aniRelated.Should().BeLessThan(1.0, "a divergent genome is NOT identical → ANI strictly below 1.0");
        aniUnrelated.Should().Be(0.0, "no conserved fragment → 0 (INV-05)");
        aniRelated.Should().BeGreaterThan(aniUnrelated, "more conserved sequence must yield a higher ANI than no conservation");
    }

    /// <summary>
    /// Positive sanity — case-insensitivity (§3.3). Lowercase input is uppercased before
    /// comparison, so a lowercase genome compared against its uppercase twin is still a
    /// perfect self-identity (ANI 1.0). Pinning this guards the documented ToUpperInvariant
    /// normalization from silently treating case as divergence.
    /// </summary>
    [Test]
    public void CalculateANI_LowercaseInput_IsUppercasedBeforeComparison()
    {
        string genome = RandomDna(200, seed: 31337);

        double ani = ComparativeGenomics.CalculateANI(genome.ToLowerInvariant(), genome, fragmentLength: 50);

        AssertWellFormedAni(ani);
        ani.Should().Be(1.0, "case is normalized via ToUpperInvariant, so a lowercase genome is identical to its uppercase twin (§3.3)");
    }

    /// <summary>
    /// Fuzz sweep — across many random genomes, fragment lengths and divergence levels,
    /// CalculateANI must ALWAYS return a finite, in-range [0,1] value and never throw.
    /// This is the broad undisciplined-failure net behind the targeted boundary tests.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void CalculateANI_RandomGenomes_AlwaysFiniteAndInRange()
    {
        for (int seed = 0; seed < 40; seed++)
        {
            var rng = new Random(seed);
            int len = rng.Next(1, 300);
            string query = RandomDna(len, seed * 7 + 1);
            string reference = RandomDna(rng.Next(1, 300), seed * 7 + 2);
            int fragLen = rng.Next(1, 60);

            double ani = 0;
            FluentActions.Invoking(() => ani = ComparativeGenomics.CalculateANI(query, reference, fragmentLength: fragLen))
                .Should().NotThrow($"random genomes must never crash ANI (seed {seed}, fragLen {fragLen})");

            AssertWellFormedAni(ani);
        }
    }

    #endregion

    #endregion
}
