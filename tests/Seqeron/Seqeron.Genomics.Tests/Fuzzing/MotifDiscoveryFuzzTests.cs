using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Matching area — Overrepresented k-mer Motif Discovery
/// (MOTIF-DISCOVER-001), the single-sequence de novo motif discoverer
/// <see cref="MotifFinder.DiscoverMotifs(DnaSequence, int, int)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain parameter values to a
/// unit and asserts the code NEVER fails in an undisciplined way: no hang/infinite
/// loop (the single O(N·k) window scan must always terminate), no state corruption,
/// no nonsense output (a reported motif whose length is not k, a Count that does not
/// equal the true overlapping-occurrence count, a Count below the minCount threshold,
/// a duplicate motif string, a non-positive or mis-computed Enrichment), and no
/// *unhandled* runtime exception — in particular NO IndexOutOfRange / negative-length
/// Substring when k &gt; N, and NO DivideByZero / Infinity / NaN Enrichment on the
/// no-window boundary. Every input must resolve to EITHER a well-defined,
/// theory-correct result OR a *documented, intentional* validation exception
/// (ArgumentNullException for a null sequence; ArgumentOutOfRangeException for k &lt; 1
/// — contract §3.3, §6.1). A raw runtime exception, a hang, a false motif below the
/// threshold, or a wrong Count/Enrichment is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MOTIF-DISCOVER-001 — Overrepresented k-mer Motif Discovery
/// Checklist: docs/checklists/03_FUZZING.md, row 170.
/// Algorithm doc: docs/algorithms/Motif_Discovery/Overrepresented_Kmer_Discovery.md
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the row:
///          – k = 1: the smallest VALID k → single-residue motifs counted over every
///            position; no crash, every reported motif has length 1 and is a single
///            nucleotide (§3.1 k ≥ 1; §2.2 windows = N − 1 + 1 = N).
///          – k &gt; len: k larger than the sequence → ZERO length-k windows → an EMPTY
///            result, with NO IndexOutOfRange / negative-length Substring and NO
///            DivideByZero / Infinity / NaN Enrichment (§3.3, §6.1 "k &gt; N → empty").
///          – no recurrence: a sequence in which no k-mer occurs ≥ minCount → an EMPTY
///            result, with NO false motif emitted below the threshold (INV-03).
///          – k &lt; 1 (k = 0, k = −1, int.MinValue): the documented
///            ArgumentOutOfRangeException, never an empty/garbage result (§3.3, §6.1).
/// — docs/checklists/03_FUZZING.md §Description (BE = граничні значення: 0, -1, MaxInt, empty).
///
/// Note on Malformed Content: the unit accepts a <see cref="DnaSequence"/>, which is
/// uppercased and validated to the {A,C,G,T} alphabet at construction, so out-of-domain
/// residues cannot reach this method; this is therefore a pure boundary (BE) row over
/// the integer parameters k / minCount and the sequence length, per the checklist.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Overrepresented_Kmer_Discovery.md §2.4, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
/// For a sequence of length N, every length-k window start 0 … N−k is enumerated;
/// each distinct k-mer's Count is the number of overlapping windows it occupies and
/// its Positions are those 0-based starts (INV-01). A k-mer is reported iff its
/// Count ≥ minCount (INV-03). Its Enrichment is the observed/expected ratio
/// Count / ((N − k + 1) / 4^k) under the i.i.d. uniform background (INV-02), which is
/// strictly positive whenever any k-mer exists (INV-04). Validation (§3.3): null
/// sequence → ArgumentNullException; k &lt; 1 → ArgumentOutOfRangeException; k &gt; N →
/// empty result (no windows).
///   MotifFinder.DiscoverMotifs(DnaSequence, int k = 6, int minCount = 2)
///       → IEnumerable&lt;DiscoveredMotif(Sequence, Count, Positions, Enrichment)&gt;
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MotifDiscoveryFuzzTests
{
    private const int DnaAlphabetSize = 4;
    private static readonly char[] Alphabet = { 'A', 'C', 'G', 'T' };

    #region Helpers

    /// <summary>A random ACGT string of the given length.</summary>
    private static string RandomDna(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(Alphabet[rng.Next(Alphabet.Length)]);
        return sb.ToString();
    }

    /// <summary>
    /// Independent oracle implementing the documented decision rule verbatim
    /// (Overrepresented_Kmer_Discovery.md §4.1): slide a length-k window over the
    /// sequence, tally overlapping occurrences with their 0-based starts, compute
    /// E = (N − k + 1) / 4^k, and return every k-mer with Count ≥ minCount and its
    /// enrichment Count / E. Built from the spec, not from the unit's code.
    /// </summary>
    private static Dictionary<string, (int Count, double Enrichment, List<int> Positions)> Oracle(
        string seq, int k, int minCount)
    {
        var positions = new Dictionary<string, List<int>>();
        for (int i = 0; i <= seq.Length - k; i++)
        {
            string kmer = seq.Substring(i, k);
            if (!positions.TryGetValue(kmer, out var list))
                positions[kmer] = list = new List<int>();
            list.Add(i);
        }

        double expected = (seq.Length - k + 1.0) / Math.Pow(DnaAlphabetSize, k);

        var result = new Dictionary<string, (int, double, List<int>)>();
        foreach (var (kmer, pos) in positions)
            if (pos.Count >= minCount)
                result[kmer] = (pos.Count, pos.Count / expected, pos);
        return result;
    }

    /// <summary>
    /// Asserts a discovery result is WELL-FORMED per the documented contract,
    /// independent of the (possibly degenerate) input:
    ///   • every reported motif string has length exactly k (window length);
    ///   • motif strings are DISTINCT (one record per k-mer);
    ///   • Count ≥ minCount (INV-03 — the threshold is honoured, no false motif);
    ///   • Count == Positions.Count and each position is a valid window start whose
    ///     substring equals the motif (INV-01);
    ///   • Enrichment is finite, strictly positive, and equals Count / E (INV-02/04).
    /// </summary>
    private static void AssertWellFormed(IReadOnlyList<DiscoveredMotif> motifs, string seq, int k, int minCount)
    {
        motifs.Select(m => m.Sequence).Should().OnlyHaveUniqueItems("one record per distinct k-mer");

        double expected = (seq.Length - k + 1.0) / Math.Pow(DnaAlphabetSize, k);

        foreach (var m in motifs)
        {
            m.Sequence.Should().HaveLength(k, "every reported motif is a length-k window (INV-01)");
            m.Count.Should().BeGreaterThanOrEqualTo(minCount, "INV-03: only k-mers meeting the threshold are reported");
            m.Positions.Should().HaveCount(m.Count, "INV-01: Count == number of occurrence positions");

            foreach (int p in m.Positions)
            {
                p.Should().BeInRange(0, seq.Length - k, "positions are valid 0-based window starts (INV-01)");
                seq.Substring(p, k).Should().Be(m.Sequence, "INV-01: each position is a genuine occurrence");
            }

            double.IsFinite(m.Enrichment).Should().BeTrue("INV-04: enrichment is finite (no DivideByZero/NaN/Infinity)");
            m.Enrichment.Should().BeGreaterThan(0, "INV-04: enrichment > 0 for every reported motif");
            m.Enrichment.Should().BeApproximately(m.Count / expected, 1e-9 * Math.Max(1, m.Count / expected),
                "INV-02: enrichment == Count / ((N − k + 1) / 4^k)");
        }
    }

    #endregion

    #region MOTIF-DISCOVER-001 — Overrepresented k-mer Motif Discovery (BE: k=1, k>len, no recurrence)

    #region Positive sanity — hand-computed documented discovery

    // Documented worked example (§7.1): "ATGCATGCATGC" (N=12), k=4 → window count 9,
    // E = 9/4^4 = 9/256; "ATGC" occurs at 0,4,8 (Count 3), enrichment 3/(9/256) = 768/9.
    [Test]
    public void DiscoverMotifs_DocumentedWorkedExample_ReportsAtgcWithCountAndEnrichment()
    {
        var motifs = MotifFinder.DiscoverMotifs(new DnaSequence("ATGCATGCATGC"), k: 4, minCount: 2).ToList();

        var atgc = motifs.Single(m => m.Sequence == "ATGC");
        atgc.Count.Should().Be(3, "ATGC occupies windows 0,4,8 (§7.1)");
        atgc.Positions.Should().Equal(new[] { 0, 4, 8 }, "0-based window starts of the recurring motif (INV-01)");
        atgc.Enrichment.Should().BeApproximately(768.0 / 9.0, 1e-9, "3 / (9/256) = 768/9 (§7.1, INV-02)");

        AssertWellFormed(motifs, "ATGCATGCATGC", 4, 2);
    }

    // A clearly recurring motif ("ACGT" repeated four times) is reported with its true
    // overlapping count; the documented POSITIVE recurrence case.
    [Test]
    public void DiscoverMotifs_ClearlyRecurringMotif_IsReportedWithDocumentedCount()
    {
        const string seq = "ACGTACGTACGTACGT"; // N=16, four tandem copies of ACGT
        var motifs = MotifFinder.DiscoverMotifs(new DnaSequence(seq), k: 4, minCount: 2).ToList();

        var acgt = motifs.Single(m => m.Sequence == "ACGT");
        acgt.Count.Should().Be(4, "ACGT starts at 0,4,8,12");
        acgt.Positions.Should().Equal(0, 4, 8, 12);
        acgt.Enrichment.Should().BeGreaterThan(1.0, "a tandem-repeated motif is overrepresented vs chance");

        AssertWellFormed(motifs, seq, 4, 2);
    }

    // Overlapping occurrences ARE counted (§2.3 ASM-02): "AAAA", k=2 → "AA" at 0,1,2 (Count 3).
    [Test]
    public void DiscoverMotifs_OverlappingOccurrences_AreCounted()
    {
        var motifs = MotifFinder.DiscoverMotifs(new DnaSequence("AAAA"), k: 2, minCount: 2).ToList();

        var aa = motifs.Single(m => m.Sequence == "AA");
        aa.Count.Should().Be(3, "overlapping windows 0,1,2 all yield AA (ASM-02)");
        aa.Positions.Should().Equal(0, 1, 2);
        AssertWellFormed(motifs, "AAAA", 2, 2);
    }

    #endregion

    #region BE — Boundary: k = 1 (smallest valid k → single-residue motifs, no crash)

    // k=1 over a sequence with a repeated residue → single-nucleotide motifs whose
    // counts match the residue frequencies; every reported motif has length 1.
    [Test]
    public void DiscoverMotifs_K1_ReportsSingleResidueMotifs()
    {
        // A×4, C×2, G×1, T×1 over N=8.
        const string seq = "AAAACCGT";
        var motifs = MotifFinder.DiscoverMotifs(new DnaSequence(seq), k: 1, minCount: 2).ToList();

        motifs.Should().OnlyContain(m => m.Sequence.Length == 1, "k=1 ⇒ single-residue motifs");
        motifs.Single(m => m.Sequence == "A").Count.Should().Be(4);
        motifs.Single(m => m.Sequence == "C").Count.Should().Be(2);
        motifs.Should().NotContain(m => m.Sequence == "G", "G occurs once, below minCount=2 (INV-03)");
        motifs.Should().NotContain(m => m.Sequence == "T", "T occurs once, below minCount=2 (INV-03)");

        AssertWellFormed(motifs, seq, 1, 2);
    }

    // Fuzz: k=1 on random sequences never throws and is always well-formed.
    [Test]
    [CancelAfter(30_000)]
    public void DiscoverMotifs_K1_RandomSequences_NeverThrows_MatchesOracle()
    {
        var rng = new Random(170_001);
        for (int trial = 0; trial < 800; trial++)
        {
            int n = rng.Next(0, 60);
            int minCount = rng.Next(1, 5);
            string seq = RandomDna(rng, n);

            var motifs = MotifFinder.DiscoverMotifs(new DnaSequence(seq), k: 1, minCount: minCount).ToList();

            AssertWellFormed(motifs, seq, 1, minCount);
            AssertMatchesOracle(motifs, seq, 1, minCount);
        }
    }

    #endregion

    #region BE — Boundary: k > len (no windows → empty, no IndexOutOfRange / negative-length Substring)

    // k strictly larger than N → zero windows → empty result, no negative-length Substring.
    [Test]
    public void DiscoverMotifs_KGreaterThanLength_ReturnsEmpty()
    {
        MotifFinder.DiscoverMotifs(new DnaSequence("ACGT"), k: 5, minCount: 1)
            .Should().BeEmpty("k > N has no length-k windows (§6.1)");
    }

    // k == N + 1 boundary (k = N still has one window; k = N+1 has none).
    [Test]
    public void DiscoverMotifs_KEqualsLengthPlusOne_ReturnsEmpty_WhileKEqualsLengthHasOneWindow()
    {
        var single = MotifFinder.DiscoverMotifs(new DnaSequence("ACGTACGT"), k: 8, minCount: 1).ToList();
        single.Should().ContainSingle("k == N yields exactly one window (the whole sequence)");
        single[0].Sequence.Should().Be("ACGTACGT");
        single[0].Count.Should().Be(1);

        MotifFinder.DiscoverMotifs(new DnaSequence("ACGTACGT"), k: 9, minCount: 1)
            .Should().BeEmpty("k = N + 1 ⇒ no windows (§6.1)");
    }

    // Empty sequence with any k ≥ 1 → no windows → empty result, no crash.
    [Test]
    public void DiscoverMotifs_EmptySequence_ReturnsEmpty()
    {
        MotifFinder.DiscoverMotifs(new DnaSequence(""), k: 1, minCount: 1)
            .Should().BeEmpty("an empty sequence has no windows for any k ≥ 1");
        MotifFinder.DiscoverMotifs(new DnaSequence(""), k: 6, minCount: 2)
            .Should().BeEmpty("an empty sequence has no windows for any k ≥ 1");
    }

    // Fuzz: k strictly greater than N never throws and always returns empty —
    // no IndexOutOfRange, no negative-length Substring, no NaN/Infinity enrichment.
    [Test]
    [CancelAfter(30_000)]
    public void DiscoverMotifs_KGreaterThanLength_RandomSequences_AlwaysEmpty_NeverThrows()
    {
        var rng = new Random(170_002);
        for (int trial = 0; trial < 800; trial++)
        {
            int n = rng.Next(0, 40);
            string seq = RandomDna(rng, n);
            int k = n + rng.Next(1, 20); // strictly > n

            Action act = () => MotifFinder.DiscoverMotifs(new DnaSequence(seq), k: k, minCount: 1).ToList();

            act.Should().NotThrow("k > N is a documented empty-result boundary, not an error (§6.1)");
            MotifFinder.DiscoverMotifs(new DnaSequence(seq), k: k, minCount: 1)
                .Should().BeEmpty("no length-k windows exist when k > N");
        }
    }

    #endregion

    #region BE — Boundary: no recurrence (nothing meets the threshold → empty, no false motif)

    // A sequence in which every k-mer is unique → with minCount ≥ 2 nothing recurs → empty.
    [Test]
    public void DiscoverMotifs_NoRecurrence_AllKmersUnique_ReturnsEmpty()
    {
        // Every 3-mer of "ACGTAG" : ACG,CGT,GTA,TAG — all distinct.
        MotifFinder.DiscoverMotifs(new DnaSequence("ACGTAG"), k: 3, minCount: 2)
            .Should().BeEmpty("no 3-mer recurs ⇒ no motif meets minCount=2 (INV-03)");
    }

    // minCount above every possible count → empty even though windows exist (no false motif).
    [Test]
    public void DiscoverMotifs_MinCountAboveAnyPossibleCount_ReturnsEmpty()
    {
        // "AAAA", k=1: A occurs 4 times; minCount=5 is unattainable.
        MotifFinder.DiscoverMotifs(new DnaSequence("AAAA"), k: 1, minCount: 5)
            .Should().BeEmpty("no k-mer can reach an unattainable threshold (INV-03)");
    }

    // minCount == int.MaxValue → unattainable → empty, no overflow/crash.
    [Test]
    public void DiscoverMotifs_MinCountMaxValue_ReturnsEmpty()
    {
        MotifFinder.DiscoverMotifs(new DnaSequence("ACGTACGTACGT"), k: 4, minCount: int.MaxValue)
            .Should().BeEmpty("int.MaxValue threshold is unattainable, no overflow");
    }

    // Fuzz: when every k-mer is distinct, minCount ≥ 2 yields nothing — no false motif ever.
    [Test]
    [CancelAfter(30_000)]
    public void DiscoverMotifs_NoRecurrence_NeverEmitsFalseMotif()
    {
        var rng = new Random(170_003);
        for (int trial = 0; trial < 600; trial++)
        {
            int n = rng.Next(0, 50);
            int k = rng.Next(1, 12);
            int minCount = rng.Next(2, 6);
            string seq = RandomDna(rng, n);

            var motifs = MotifFinder.DiscoverMotifs(new DnaSequence(seq), k: k, minCount: minCount).ToList();

            // Whatever is reported must genuinely meet the threshold; nothing fabricated.
            AssertWellFormed(motifs, seq, k, minCount);
            AssertMatchesOracle(motifs, seq, k, minCount);
        }
    }

    #endregion

    #region BE — Boundary: k < 1 guard (k = 0, k = -1, int.MinValue → ArgumentOutOfRangeException)

    [TestCase(0)]
    [TestCase(-1)]
    [TestCase(int.MinValue)]
    public void DiscoverMotifs_KBelowOne_ThrowsArgumentOutOfRange(int k)
    {
        Action act = () => MotifFinder.DiscoverMotifs(new DnaSequence("ACGTACGT"), k: k, minCount: 2).ToList();

        act.Should().Throw<ArgumentOutOfRangeException>(
            "k < 1 is the documented validation contract (§3.3, §6.1)")
            .Which.ParamName.Should().Be("k");
    }

    // Null sequence → documented ArgumentNullException.
    [Test]
    public void DiscoverMotifs_NullSequence_Throws()
    {
        Action act = () => MotifFinder.DiscoverMotifs(null!, k: 4, minCount: 2).ToList();

        act.Should().Throw<ArgumentNullException>("null sequence is the documented validation contract (§3.3)");
    }

    #endregion

    #region BE — Broad fuzz: random k / minCount / length never crash, match the documented rule

    [Test]
    [CancelAfter(60_000)]
    public void DiscoverMotifs_RandomInputs_NeverThrows_MatchesOracle()
    {
        var rng = new Random(170_004);
        for (int trial = 0; trial < 2000; trial++)
        {
            int n = rng.Next(0, 80);
            int k = rng.Next(1, 14);
            int minCount = rng.Next(1, 6);
            string seq = RandomDna(rng, n);

            var motifs = MotifFinder.DiscoverMotifs(new DnaSequence(seq), k: k, minCount: minCount).ToList();

            AssertWellFormed(motifs, seq, k, minCount);
            AssertMatchesOracle(motifs, seq, k, minCount);
        }
    }

    // Homopolymer (§6.1): a single k-mer dominates with high enrichment, no crash.
    [Test]
    public void DiscoverMotifs_Homopolymer_SingleDominantMotif()
    {
        const string seq = "AAAAAAAAAAAA"; // N=12
        var motifs = MotifFinder.DiscoverMotifs(new DnaSequence(seq), k: 4, minCount: 2).ToList();

        motifs.Should().ContainSingle("only one distinct 4-mer (AAAA) exists in a homopolymer");
        motifs[0].Sequence.Should().Be("AAAA");
        motifs[0].Count.Should().Be(9, "9 overlapping windows of AAAA in N=12");
        motifs[0].Enrichment.Should().BeGreaterThan(1.0, "the lone k-mer is heavily overrepresented (§6.1)");
        AssertWellFormed(motifs, seq, 4, 2);
    }

    #endregion

    #endregion

    #region Oracle cross-check helper

    /// <summary>
    /// Asserts the unit's reported motif SET matches the independent oracle exactly
    /// (same k-mers, same counts, same positions, same enrichment) — the strict
    /// theory-correct cross-check used by the fuzz loops.
    /// </summary>
    private static void AssertMatchesOracle(IReadOnlyList<DiscoveredMotif> motifs, string seq, int k, int minCount)
    {
        var oracle = Oracle(seq, k, minCount);

        motifs.Select(m => m.Sequence).Should().BeEquivalentTo(oracle.Keys,
            "the reported motif set equals the documented Count ≥ minCount set");

        foreach (var m in motifs)
        {
            var (count, enrichment, positions) = oracle[m.Sequence];
            m.Count.Should().Be(count, "Count matches the documented overlapping-window count (INV-01)");
            m.Positions.Should().Equal(positions, "Positions match the documented 0-based window starts (INV-01)");
            m.Enrichment.Should().BeApproximately(enrichment, 1e-9 * Math.Max(1, enrichment),
                "Enrichment matches the documented O/E ratio (INV-02)");
        }
    }

    #endregion
}
