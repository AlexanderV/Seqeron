using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the K-mer area — K-MER GENERATION
/// (KMER-GENERATE-001): enumerating the complete k-mer universe Σ^k.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain parameter values to a
/// unit and asserts that the code NEVER fails in an undisciplined way: no hang or
/// infinite loop, no state corruption, no nonsense output, and no *unhandled*
/// runtime exception (IndexOutOfRangeException, OutOfMemoryException, a leaking
/// NullReferenceException). Every input must resolve to EITHER a well-defined,
/// theory-correct result, OR a *documented, intentional* validation exception
/// (ArgumentException / ArgumentOutOfRangeException). A raw runtime exception, a
/// hang, an OOM, a wrong count (off the |Σ|^k formula), or a duplicate k-mer is a
/// bug, not a passing test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: KMER-GENERATE-001 — k-mer generation (k-fold Cartesian product Σ^k)
/// Checklist: docs/checklists/03_FUZZING.md, row 159.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: k = 0 (and the k ≤ 0 floor: -1, int.MinValue), a
///          "large k" stressing the n^k explosion, and a non-DNA / empty alphabet.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The k-mer generation contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// K-mer Generation enumerates EVERY possible length-k string over an alphabet —
/// the complete k-mer universe, independent of any sequence. For an alphabet of
/// size n there are exactly n^k such k-mers, formed as the k-fold Cartesian
/// product Σ^k (K-mer_Generation.md §2.2; Wikipedia "K-mer"; itertools.product).
/// The API entry under test is
///   KmerAnalyzer.GenerateAllKmers(int k, string alphabet = "ACGT")
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs lines
///   299–325), backed by the private recursive prefix-extension
///   GenerateKmersRecursive (lines 310–325) which realises the Cartesian product.
///
/// The DOCUMENTED invariants (K-mer_Generation.md §2.4):
///   • INV-01  output COUNT = n^k  (4^k for the default DNA alphabet).
///   • INV-02  all emitted k-mers are DISTINCT (the output is a set of size n^k)
///             — PROVIDED the alphabet has no duplicate characters (§6.2: the
///             alphabet is NOT de-duplicated, so a repeated character yields
///             repeated k-mers; supply a clean alphabet for distinct output).
///   • INV-03  every k-mer has length EXACTLY k and contains ONLY alphabet chars.
///   • INV-04  for a SORTED alphabet the emission order is lexicographic, with the
///             RIGHTMOST position advancing fastest (odometer ordering). The default
///             "ACGT" is already sorted, so DNA k-mers come AAA, AAC, …, TTT.
///
/// The DOCUMENTED parameter contract (K-mer_Generation.md §3.1, §3.3, §6.1;
/// KmerAnalyzer.cs lines 301–305):
///   • k must be > 0; k ≤ 0 → ArgumentOutOfRangeException(nameof(k)). NOTE: this
///     unit does NOT treat k = 0 as "the single empty k-mer" — the contract is an
///     EXCEPTION ("a k-mer length must be positive"), so the BE k = 0 target here
///     asserts the documented throw, not a count of 1.
///   • alphabet must be non-empty; null or empty → ArgumentException(nameof(
///     alphabet)). An empty alphabet has no symbols ⇒ no k-mers, surfaced as a
///     guard, never an OOM or a silent empty enumerable.
///   • the alphabet is used VERBATIM — case-sensitive, no normalisation, no sort,
///     no de-dup. A custom / non-DNA alphabet generates over the supplied symbols.
///
/// Laziness (KmerAnalyzer.cs lines 310–325, §5.2): enumeration is `yield`-driven,
/// so the n^k universe is STREAMED — a caller can read the head without
/// materialising all n^k strings. This is the guard exploited by the "large k"
/// BE target: we confirm the count FORMULA and a couple of leading k-mers for a
/// large k WITHOUT ever materialising 4^k strings (no OOM), and we Take a bounded
/// prefix under [CancelAfter] so an accidental eager-materialisation regression
/// would surface as a timeout rather than hanging the suite.
///
/// The three checklist targets map to these documented behaviours:
///   • k = 0 (and k ≤ 0)   → ArgumentOutOfRangeException; no crash, no empty list.
///   • large k             → still exactly 4^k k-mers (formula), streamed lazily;
///                           leading k-mers correct; no OOM, no hang.
///   • non-DNA / empty alpha→ generates over the supplied alphabet (count = n^k,
///                           length-k, only-alphabet chars); empty alphabet throws
///                           ArgumentException; no crash.
/// A POSITIVE sanity test pins generate(k=1) == the alphabet itself and
/// generate(k=2) over DNA == all 16 dinucleotides exactly (distinct, length 2,
/// lexicographic), cross-checked against count = 4^k.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class KmerGenerateFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomAlphabet(int distinctCount, int seed)
    {
        // A pool of DISTINCT non-DNA-shaped symbols (letters + digits + punctuation),
        // so the produced alphabet has no repeats (INV-02 precondition holds).
        const string pool = "BDEFHIJKLMNOPQRSUVWXYZ0123456789-_.*#@";
        var rng = new Random(seed);
        return new string(pool.OrderBy(_ => rng.Next()).Take(distinctCount).ToArray());
    }

    /// <summary>
    /// Asserts the result is a WELL-FORMED k-mer universe for the given k and
    /// alphabet: there are EXACTLY |alphabet|^k k-mers (INV-01), every k-mer is
    /// DISTINCT (INV-02, valid because callers below pass duplicate-free alphabets),
    /// every k-mer has length EXACTLY k and contains ONLY alphabet characters
    /// (INV-03). This is the load-bearing structural oracle reused across the fuzz
    /// cases (K-mer_Generation.md §2.4).
    /// </summary>
    private static void AssertWellFormedUniverse(IReadOnlyList<string> kmers, int k, string alphabet)
    {
        var alphabetSet = new HashSet<char>(alphabet);

        // INV-01: count = |Σ|^k exactly.
        long expectedCount = checked((long)Math.Pow(alphabet.Length, k));
        kmers.Count.Should().Be((int)expectedCount,
            $"INV-01: the k-mer universe over a size-{alphabet.Length} alphabet at k={k} has exactly {alphabet.Length}^{k} members");

        // INV-02: all distinct (alphabet is duplicate-free in these cases).
        kmers.Distinct().Should().HaveCount(kmers.Count,
            "INV-02: every k-mer is a unique length-k tuple over a duplicate-free alphabet");

        // INV-03: every k-mer has length exactly k and only alphabet chars.
        kmers.Should().OnlyContain(km => km.Length == k,
            "INV-03: every generated k-mer has length exactly k");
        kmers.Should().OnlyContain(km => km.All(c => alphabetSet.Contains(c)),
            "INV-03: every k-mer character is drawn from the supplied alphabet");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  KMER-GENERATE-001 — k-mer generation (Σ^k) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region KMER-GENERATE-001 — k-mer generation

    #region Positive sanity — the documented Σ^k contract on known cases

    /// <summary>
    /// POSITIVE / sanity: generate(k=1) over the default DNA alphabet is the alphabet
    /// itself — n^1 = n single-character k-mers, one per symbol, in alphabet order
    /// (K-mer_Generation.md §6.1 "k = 1 → one k-mer per alphabet symbol").
    /// </summary>
    [Test]
    public void GenerateAllKmers_KEquals1_IsTheAlphabetItself()
    {
        var kmers = KmerAnalyzer.GenerateAllKmers(1).ToList();

        kmers.Should().Equal(new[] { "A", "C", "G", "T" },
            "k = 1 yields exactly the alphabet's symbols, each as a length-1 k-mer, in order");
        AssertWellFormedUniverse(kmers, 1, "ACGT");
    }

    /// <summary>
    /// POSITIVE / sanity: generate(k=2) over DNA is ALL 16 (= 4^2) dinucleotides,
    /// distinct, length 2, in lexicographic odometer order AA, AC, …, TT — the
    /// worked example from K-mer_Generation.md §7.1 (INV-01–INV-04).
    /// </summary>
    [Test]
    public void GenerateAllKmers_KEquals2_AreAll16DinucleotidesLexicographic()
    {
        var kmers = KmerAnalyzer.GenerateAllKmers(2).ToList();

        kmers.Should().Equal(
            "AA", "AC", "AG", "AT",
            "CA", "CC", "CG", "CT",
            "GA", "GC", "GG", "GT",
            "TA", "TC", "TG", "TT");
        AssertWellFormedUniverse(kmers, 2, "ACGT");
    }

    /// <summary>
    /// POSITIVE: the count formula holds across a range of small DNA k — the universe
    /// is exactly 4^k for k = 1..6 (4, 16, 64, 256, 1024, 4096), distinct, length-k,
    /// only-ACGT, lexicographically ordered for the sorted default alphabet (INV-01–04).
    /// </summary>
    [Test]
    [TestCase(1, 4)]
    [TestCase(2, 16)]
    [TestCase(3, 64)]
    [TestCase(4, 256)]
    [TestCase(5, 1024)]
    [TestCase(6, 4096)]
    public void GenerateAllKmers_SmallDnaK_HasExactly4PowK(int k, int expectedCount)
    {
        var kmers = KmerAnalyzer.GenerateAllKmers(k).ToList();

        kmers.Should().HaveCount(expectedCount, $"INV-01: 4^{k} = {expectedCount}");
        AssertWellFormedUniverse(kmers, k, "ACGT");
        kmers.Should().BeInAscendingOrder(StringComparer.Ordinal,
            "INV-04: the sorted default alphabet 'ACGT' produces lexicographic output");
    }

    #endregion

    #region BE — Boundary: k = 0 and the k ≤ 0 floor

    /// <summary>
    /// BE / HEADLINE: k = 0 is NOT "the single empty k-mer" for this unit — the
    /// documented contract is an EXCEPTION ("a k-mer length must be positive",
    /// K-mer_Generation.md §3.3, §6.1; KmerAnalyzer.cs lines 301–302). We assert the
    /// documented ArgumentOutOfRangeException on the param k, never an empty/1-element
    /// enumerable and never a crash.
    ///
    /// NOTE on laziness: GenerateAllKmers validates EAGERLY (the throw is in the public
    /// method body, before the iterator is created at line 307), so the exception
    /// surfaces on the call itself — no need to enumerate to trigger it.
    /// </summary>
    [Test]
    public void GenerateAllKmers_KEqualsZero_ThrowsArgumentOutOfRange()
    {
        Action act = () => KmerAnalyzer.GenerateAllKmers(0);

        act.Should().Throw<ArgumentOutOfRangeException>(
                "k = 0 is rejected: a k-mer length must be positive (§3.3)")
            .Which.ParamName.Should().Be("k");
    }

    /// <summary>
    /// BE: the rest of the k ≤ 0 floor — negative k and int.MinValue — all throw the
    /// same documented ArgumentOutOfRangeException(nameof(k)); no overflow, no crash.
    /// </summary>
    [Test]
    [TestCase(-1)]
    [TestCase(-42)]
    [TestCase(int.MinValue)]
    public void GenerateAllKmers_NegativeK_ThrowsArgumentOutOfRange(int k)
    {
        Action act = () => KmerAnalyzer.GenerateAllKmers(k);

        act.Should().Throw<ArgumentOutOfRangeException>(
                "k ≤ 0 is rejected regardless of magnitude")
            .Which.ParamName.Should().Be("k");
    }

    #endregion

    #region BE — Boundary: large k (the n^k explosion — count formula + lazy streaming)

    /// <summary>
    /// BE / HEADLINE: "large k" stresses the n^k explosion. We never MATERIALISE the
    /// full universe — instead we (a) confirm the count FORMULA holds at a moderately
    /// large k that is still safe to enumerate (k = 8 over DNA → 4^8 = 65 536 k-mers),
    /// and (b) exploit the documented LAZY streaming (§5.2): for a genuinely huge
    /// k = 20 (4^20 ≈ 1.1×10^12, far beyond memory) we read only a bounded PREFIX of
    /// the enumerable and check the leading k-mers, so an accidental eager-
    /// materialisation regression would OOM/timeout instead of passing. [CancelAfter]
    /// guards against a hang.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void GenerateAllKmers_LargeK_FormulaHoldsAndStreamsLazilyWithoutOom()
    {
        // (a) moderately large k that is safe to fully enumerate: 4^8 = 65 536.
        var k8 = KmerAnalyzer.GenerateAllKmers(8).ToList();
        AssertWellFormedUniverse(k8, 8, "ACGT");
        k8.Should().HaveCount(65_536, "INV-01: 4^8 = 65 536");
        k8.First().Should().Be("AAAAAAAA", "lexicographic first k-mer is the all-min-symbol string");
        k8.Last().Should().Be("TTTTTTTT", "lexicographic last k-mer is the all-max-symbol string");

        // (b) genuinely huge k = 20 (4^20 ≈ 1.1e12) — enumerable MUST be lazy, so a
        // bounded prefix returns instantly without trying to build 4^20 strings.
        const int hugeK = 20;
        var prefix = KmerAnalyzer.GenerateAllKmers(hugeK).Take(5).ToList();

        prefix.Should().HaveCount(5, "Take(5) on a lazy enumerable yields just 5 k-mers, never the full universe");
        prefix.Should().OnlyContain(km => km.Length == hugeK,
            "INV-03: even at the huge k each streamed k-mer has length exactly k");
        prefix[0].Should().Be(new string('A', hugeK),
            "the FIRST k-mer in odometer order is the all-'A' string (rightmost varies fastest)");
        prefix[1].Should().Be(new string('A', hugeK - 1) + "C",
            "INV-04: the rightmost position advances first → A…AC follows A…AA");
    }

    #endregion

    #region BE — Boundary: non-DNA alphabet and empty alphabet

    /// <summary>
    /// BE / HEADLINE: a custom NON-DNA alphabet generates over the SUPPLIED symbols —
    /// the count is |Σ|^k, not 4^k. Protein-flavoured 20-symbol alphabet at k = 2 →
    /// 400 distinct di-residues; the alphabet is used verbatim (§3.1, §5.2). We pin
    /// the count, distinctness, length and only-alphabet membership.
    /// </summary>
    [Test]
    public void GenerateAllKmers_NonDnaProteinAlphabet_GeneratesOverSuppliedSymbols()
    {
        const string protein = "ACDEFGHIKLMNPQRSTVWY"; // 20 standard amino acids, sorted.
        var kmers = KmerAnalyzer.GenerateAllKmers(2, protein).ToList();

        AssertWellFormedUniverse(kmers, 2, protein);
        kmers.Should().HaveCount(400, "INV-01: 20^2 = 400 di-residues");
        kmers.Should().BeInAscendingOrder(StringComparer.Ordinal,
            "INV-04: the sorted protein alphabet yields lexicographic output");
    }

    /// <summary>
    /// BE: a randomly-built duplicate-free NON-DNA alphabet of varying size still
    /// satisfies the universe contract count = n^k over the supplied symbols, for a
    /// spread of (alphabetSize, k). Inputs are LOCALLY seeded for reproducibility.
    /// </summary>
    [Test]
    [TestCase(1, 5, 11)]   // 1^5 = 1   (single-symbol alphabet, any k → homopolymer)
    [TestCase(2, 6, 22)]   // 2^6 = 64
    [TestCase(3, 4, 33)]   // 3^4 = 81
    [TestCase(5, 3, 44)]   // 5^3 = 125
    [TestCase(7, 2, 55)]   // 7^2 = 49
    public void GenerateAllKmers_RandomNonDnaAlphabet_CountIsAlphabetSizePowK(
        int alphabetSize, int k, int seed)
    {
        string alphabet = RandomAlphabet(alphabetSize, seed);
        alphabet.Distinct().Should().HaveCount(alphabetSize, "the test alphabet is duplicate-free by construction");

        var kmers = KmerAnalyzer.GenerateAllKmers(k, alphabet).ToList();

        AssertWellFormedUniverse(kmers, k, alphabet);
    }

    /// <summary>
    /// BE / edge: a single-letter alphabet at any k yields EXACTLY ONE k-mer — the
    /// homopolymer (1^k = 1, K-mer_Generation.md §6.1). No crash, no empty result.
    /// </summary>
    [Test]
    [TestCase(1)]
    [TestCase(3)]
    [TestCase(10)]
    public void GenerateAllKmers_SingleLetterAlphabet_YieldsExactlyOneHomopolymer(int k)
    {
        var kmers = KmerAnalyzer.GenerateAllKmers(k, "X").ToList();

        kmers.Should().ContainSingle("1^k = 1 — a single-symbol alphabet has exactly one length-k string")
            .Which.Should().Be(new string('X', k));
        AssertWellFormedUniverse(kmers, k, "X");
    }

    /// <summary>
    /// BE / HEADLINE: an EMPTY alphabet has no symbols ⇒ no k-mers; the documented
    /// guard throws ArgumentException(nameof(alphabet)) eagerly — never an OOM, a hang,
    /// or a silent empty enumerable (K-mer_Generation.md §3.3, §6.1; KmerAnalyzer.cs
    /// lines 304–305).
    /// </summary>
    [Test]
    public void GenerateAllKmers_EmptyAlphabet_ThrowsArgumentException()
    {
        Action act = () => KmerAnalyzer.GenerateAllKmers(3, "");

        act.Should().Throw<ArgumentException>("an empty alphabet has no symbols to combine")
            .Which.ParamName.Should().Be("alphabet");
    }

    /// <summary>
    /// BE: a null alphabet is the other half of the documented guard — same
    /// ArgumentException(nameof(alphabet)), surfaced eagerly, no NullReferenceException
    /// leaking from internal iteration (§3.3; KmerAnalyzer.cs line 304).
    /// </summary>
    [Test]
    public void GenerateAllKmers_NullAlphabet_ThrowsArgumentException()
    {
        Action act = () => KmerAnalyzer.GenerateAllKmers(3, null!);

        act.Should().Throw<ArgumentException>("a null alphabet has no symbols to combine")
            .Which.ParamName.Should().Be("alphabet");
    }

    /// <summary>
    /// BE / documented quirk: the alphabet is NOT de-duplicated (§5.2, §6.2), so a
    /// repeated character DELIBERATELY produces repeated k-mers — count is still
    /// alphabet.Length^k literally (2^2 = 4 for "AA"), but with duplicates. We pin this
    /// documented behaviour explicitly so a future "silent de-dup" change would be
    /// caught: INV-02 (distinctness) holds ONLY for duplicate-free alphabets.
    /// </summary>
    [Test]
    public void GenerateAllKmers_AlphabetWithDuplicateChar_ProducesRepeatedKmersPerContract()
    {
        // "AA".Length == 2, so 2^2 = 4 k-mers are emitted — but only "AA" is distinct.
        var kmers = KmerAnalyzer.GenerateAllKmers(2, "AA").ToList();

        kmers.Should().HaveCount(4, "the alphabet is used verbatim: |\"AA\"|^2 = 2^2 = 4 emissions");
        kmers.Should().AllBe("AA", "every emission is 'AA' — duplicate alphabet chars yield repeated k-mers (§6.2)");
        kmers.Distinct().Should().ContainSingle(
            "INV-02 distinctness is NOT guaranteed for a duplicated alphabet — this is documented, not a bug");
    }

    #endregion

    #endregion
}
