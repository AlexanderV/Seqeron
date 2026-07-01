using System;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Assembly area — Contig Merging (ASSEMBLY-MERGE-001), the
/// suffix–prefix overlap-collapse primitive
/// <see cref="SequenceAssembler.MergeContigs(string, string, int)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to the unit and
/// asserts the code NEVER fails in an undisciplined way: no hang, no
/// IndexOutOfRange / ArgumentOutOfRange in the suffix-prefix join
/// (<c>contig2.Substring(overlapLength)</c>), no nonsense output (a result that
/// does not start with <c>contig1</c>, a wrong length, a deleted character), no
/// non-deterministic output, and no *unhandled* runtime exception — in particular
/// the join must NOT throw when <c>overlapLength</c> is 0, negative, <c>int.MaxValue</c>,
/// or larger than the shorter contig. Every input must resolve to EITHER a
/// well-defined, theory-correct superstring OR a *documented, intentional*
/// validation exception (<see cref="ArgumentNullException"/> for a null contig —
/// contract §3.3, §6.1). A raw runtime exception, a hang, a wrong-length result,
/// a dropped <c>contig1</c> prefix, or an order-dependent / non-deterministic
/// result is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ASSEMBLY-MERGE-001 — Contig Merging (Suffix–Prefix Overlap Collapse)
/// Checklist: docs/checklists/03_FUZZING.md, row 144.
/// Algorithm doc: docs/algorithms/Extended_Assembly/Contig_Merging.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row (graничні значення: 0, -1, MaxInt, empty):
///          – NO OVERLAP: <c>overlapLength = 0</c> (and negative / out-of-range) →
///            the two contigs are NOT joined at a phantom overlap; they are simply
///            CONCATENATED <c>c1 + c2</c>, no false join below the threshold, no
///            crash (INV-02, INV-03, §6.1).
///          – FULL CONTAINMENT: one contig is wholly a prefix of the other; at the
///            boundary overlap <c>l = min(|c1|,|c2|) = |contig2|</c> the contained
///            contig is ABSORBED — <c>c2[l..]</c> is empty so the result is exactly
///            <c>contig1</c>, no duplication (INV-01 at the boundary, INV-04).
///          – IDENTICAL CONTIGS: <c>c1 == c2</c> merged at the full overlap
///            <c>l = |c1| = |c2|</c> de-duplicates to a SINGLE copy (<c>= c1</c>),
///            idempotently and deterministically — no infinite/runaway growth, no
///            doubled sequence (INV-01, INV-04).
/// — docs/checklists/03_FUZZING.md §Description (BE = Boundary Exploitation:
///   граничні значення 0, -1, MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Contig_Merging.md §2.4, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
/// An overlap is a length-<c>l</c> suffix of <c>c1</c> that equals a length-<c>l</c>
/// prefix of <c>c2</c>; merging keeps one copy:
///   merge(c1, c2, l) = c1 + c2[l..],  |result| = |c1| + |c2| − l   (INV-01)
/// A valid overlap is bounded: <c>0 &lt; l ≤ min(|c1|,|c2|)</c>. Otherwise (l = 0,
/// l &lt; 0, or l &gt; min(|c1|,|c2|)) the result is plain concatenation c1 + c2
/// (INV-02, INV-03). The result ALWAYS starts with <c>c1</c> verbatim and ends with
/// the non-overlapped tail <c>c2[l..]</c> (INV-04); the merge deletes nothing from
/// <c>c1</c> and nothing from <c>c2[l..]</c>. The method is case-sensitive and
/// alphabet-agnostic. Null <c>c1</c>/<c>c2</c> → ArgumentNullException (§3.3, §6.1).
///   SequenceAssembler.MergeContigs(string contig1, string contig2, int overlapLength)
///   → string
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AssemblyMergeFuzzTests
{
    // Documented no-overlap sentinel (Contig_Merging.md §5.2; SequenceAssembler NoOverlap = 0).
    private const int NoOverlap = 0;

    #region Helpers

    private static readonly char[] DnaAlphabet = { 'A', 'C', 'G', 'T' };
    // A wider alphabet to exercise the "alphabet-agnostic, case-preserved" contract (§3.3).
    private static readonly char[] WideAlphabet =
        { 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'N', '-', '.', 'x', '7', ' ' };

    private static string RandomString(Random rng, int length, char[] alphabet)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(alphabet[rng.Next(alphabet.Length)]);
        return sb.ToString();
    }

    /// <summary>
    /// Independent oracle implementing the documented merge rule verbatim
    /// (Contig_Merging.md §4.1): a valid overlap <c>0 &lt; l ≤ min(|c1|,|c2|)</c>
    /// collapses to <c>c1 + c2[l..]</c>; any other <c>l</c> concatenates <c>c1 + c2</c>.
    /// Used to cross-check the unit under test on fuzzed inputs without re-using its
    /// implementation.
    /// </summary>
    private static string Oracle(string c1, string c2, int l)
    {
        if (l <= NoOverlap || l > Math.Min(c1.Length, c2.Length))
            return c1 + c2;
        return c1 + c2.Substring(l);
    }

    /// <summary>
    /// Asserts a merged superstring is WELL-FORMED per the documented contract
    /// regardless of the (possibly degenerate) overlap length:
    ///   INV-04 the result ALWAYS starts with <c>c1</c> verbatim (no character of
    ///          <c>c1</c> is ever deleted);
    ///   INV-01/INV-02/INV-03 the length is exactly <c>|c1| + |c2| − effectiveL</c>,
    ///          where <c>effectiveL</c> is <c>l</c> for a valid overlap else 0
    ///          (a valid overlap removes EXACTLY <c>l</c> prefix chars from <c>c2</c>,
    ///          an invalid one removes none);
    ///   INV-04 the tail after the <c>c1</c> prefix equals the non-overlapped
    ///          remainder <c>c2[effectiveL..]</c>.
    /// </summary>
    private static void AssertWellFormed(string result, string c1, string c2, int l)
    {
        int effectiveL = (l <= NoOverlap || l > Math.Min(c1.Length, c2.Length)) ? 0 : l;

        result.Should().StartWith(c1, "INV-04: the result always starts with contig1 verbatim");
        result.Should().HaveLength(c1.Length + c2.Length - effectiveL,
            "INV-01/02/03: |result| = |c1| + |c2| − (valid overlap, else 0)");

        // The portion after the c1 prefix is exactly the non-overlapped tail of c2.
        string tail = result.Substring(c1.Length);
        tail.Should().Be(c2.Substring(effectiveL),
            "INV-04: the suffix after c1 is exactly the non-overlapped tail c2[l..]");
    }

    #endregion

    #region ASSEMBLY-MERGE-001 — Contig Merging (BE: no overlap, full containment, identical contigs)

    #region Positive sanity — hand-computed documented merges

    // Worked example from the doc (§7.1): suffix "AA" of BAA == prefix "AA" of AAB, overlap 2.
    [Test]
    public void MergeContigs_DocWorkedExample_BAA_AAB_Overlap2_BAAB()
    {
        string merged = SequenceAssembler.MergeContigs("BAA", "AAB", 2);

        merged.Should().Be("BAAB", "doc §7.1: c1 + c2[2..] keeps one copy of the overlap 'AA'");
        AssertWellFormed(merged, "BAA", "AAB", 2);
    }

    // A clear suffix-prefix overlap collapses to length |c1| + |c2| − l (INV-01).
    [Test]
    public void MergeContigs_ClearSuffixPrefixOverlap_JoinsAtOverlapLength()
    {
        // c1 = "ACGTACGT", c2 = "ACGTTTTT", overlap 4 ("ACGT").
        string merged = SequenceAssembler.MergeContigs("ACGTACGT", "ACGTTTTT", 4);

        merged.Should().Be("ACGTACGTTTTT", "c1 + c2[4..] = ACGTACGT + TTTT");
        merged.Should().HaveLength(8 + 8 - 4, "INV-01: length = |c1| + |c2| − l");
        AssertWellFormed(merged, "ACGTACGT", "ACGTTTTT", 4);
    }

    // The overlap is computed by FindOverlap and fed to MergeContigs — the documented
    // overlap-discovery → collapse pipeline (§2.5, §5.3): two reads with a real overlap
    // merge into the joined superstring; non-overlapping reads stay separate (no false join).
    [Test]
    public void MergeContigs_WithFindOverlap_RealOverlapJoins_NoOverlapStaysSeparate()
    {
        // Real 25-char suffix/prefix overlap (≥ default minOverlap 20).
        const string shared = "ACGTACGTACGTACGTACGTACGTA"; // 25
        string c1 = "TTTTT" + shared;        // suffix == shared
        string c2 = shared + "GGGGG";        // prefix == shared

        var ov = SequenceAssembler.FindOverlap(c1, c2);
        ov.Should().NotBeNull("a 25-char identical suffix/prefix is a valid overlap");
        string merged = SequenceAssembler.MergeContigs(c1, c2, ov!.Value.length);
        merged.Should().Be("TTTTT" + shared + "GGGGG", "the shared region appears exactly once");
        AssertWellFormed(merged, c1, c2, ov.Value.length);

        // Two clearly non-overlapping contigs → FindOverlap reports nothing → caller
        // concatenates (overlap 0); no false join below the threshold.
        string a = "AAAAAAAAAAAAAAAAAAAAAAAAA"; // 25 A
        string b = "GGGGGGGGGGGGGGGGGGGGGGGGG"; // 25 G
        SequenceAssembler.FindOverlap(a, b).Should().BeNull("no suffix/prefix match → no overlap");
        SequenceAssembler.MergeContigs(a, b, NoOverlap)
            .Should().Be(a + b, "no usable overlap → plain concatenation (INV-02)");
    }

    #endregion

    #region BE — Boundary: NO OVERLAP (l = 0 / negative / out-of-range ⇒ concatenate, no false join)

    // overlapLength = 0 → plain concatenation, no overlap collapsed (INV-02, §6.1).
    [Test]
    public void MergeContigs_ZeroOverlap_PlainConcatenation()
    {
        string merged = SequenceAssembler.MergeContigs("ACGT", "TGCA", NoOverlap);

        merged.Should().Be("ACGTTGCA", "overlap 0 ⇒ concatenate, nothing collapsed (INV-02)");
        AssertWellFormed(merged, "ACGT", "TGCA", NoOverlap);
    }

    // overlapLength = -1 (and other negatives) → concatenation, NOT a Substring crash (§6.1, BE: -1).
    [Test]
    public void MergeContigs_NegativeOverlap_ConcatenatesNoCrash()
    {
        foreach (int l in new[] { -1, -7, int.MinValue })
        {
            Action act = () => SequenceAssembler.MergeContigs("ACGT", "TGCA", l);
            act.Should().NotThrow("a negative overlap is not valid → concatenate, no crash (INV-03)");

            SequenceAssembler.MergeContigs("ACGT", "TGCA", l)
                .Should().Be("ACGTTGCA", "non-positive overlap ⇒ plain concatenation (INV-03)");
        }
    }

    // overlapLength > min(|c1|,|c2|), incl. int.MaxValue → concatenation, NO IndexOutOfRange in the
    // Substring join (§6.1, BE: MaxInt; the overlap is bounded by the shorter contig).
    [Test]
    public void MergeContigs_OverlapExceedsShorterContig_ConcatenatesNoIndexOutOfRange()
    {
        foreach (int l in new[] { 5, 100, int.MaxValue })
        {
            Action act = () => SequenceAssembler.MergeContigs("ACGT", "TG", l);
            act.Should().NotThrow("overlap > min(|c1|,|c2|) ⇒ concatenate, no Substring overflow (INV-03)");

            SequenceAssembler.MergeContigs("ACGT", "TG", l)
                .Should().Be("ACGTTG", "out-of-range overlap ⇒ plain concatenation (INV-03)");
        }
    }

    // Fuzz: random contigs with overlapLength = 0 ALWAYS yield exact concatenation, never throw.
    [Test]
    public void MergeContigs_RandomContigsZeroOverlap_AlwaysConcatenate()
    {
        var rng = new Random(144_001);
        for (int trial = 0; trial < 500; trial++)
        {
            string c1 = RandomString(rng, rng.Next(0, 40), WideAlphabet);
            string c2 = RandomString(rng, rng.Next(0, 40), WideAlphabet);

            string merged = SequenceAssembler.MergeContigs(c1, c2, NoOverlap);

            merged.Should().Be(c1 + c2, "overlap 0 ⇒ plain concatenation regardless of content (INV-02)");
            AssertWellFormed(merged, c1, c2, NoOverlap);
        }
    }

    #endregion

    #region BE — Boundary: FULL CONTAINMENT (contained contig absorbed at l = min, no duplication)

    // contig2 is wholly a prefix of contig1; at the boundary overlap l = |c2| = min, the contained
    // contig is fully absorbed: c2[l..] is empty, so the result is exactly c1 (no duplication).
    [Test]
    public void MergeContigs_ContainedContigAbsorbed_ResultIsContainer()
    {
        // c2 = "ACGT" is a prefix of c1 = "ACGTACGTACGT"; overlap = |c2| = 4 = min.
        string merged = SequenceAssembler.MergeContigs("ACGTACGTACGT", "ACGT", 4);

        merged.Should().Be("ACGTACGTACGT", "the contained contig is absorbed: c2[4..] is empty → c1");
        merged.Should().HaveLength(12, "INV-01 at the boundary: 12 + 4 − 4 = 12, no duplication");
        AssertWellFormed(merged, "ACGTACGTACGT", "ACGT", 4);
    }

    // Boundary overlap l = min(|c1|,|c2|) exactly: the entire shorter prefix collapses (§6.1 row
    // "overlapLength = min" — full collapse of the shorter prefix).
    [Test]
    public void MergeContigs_OverlapEqualsMin_FullCollapseOfShorterPrefix()
    {
        // c1 = "GGGGAC" (6), c2 = "ACTT" (4), overlap l = min = 4 ("ACTT" is the whole of c2,
        // collapsing the c2 prefix entirely → result = c1).
        string merged = SequenceAssembler.MergeContigs("GGGGAC", "ACTT", 4);

        merged.Should().Be("GGGGAC", "l = min collapses the whole shorter contig's prefix (INV-01)");
        AssertWellFormed(merged, "GGGGAC", "ACTT", 4);
    }

    // Empty contig is the degenerate full-containment case: "" is contained in anything.
    // (§6.1: empty c1 or c2 with l = 0 ⇒ the other contig.)
    [Test]
    public void MergeContigs_EmptyContig_YieldsTheOther()
    {
        SequenceAssembler.MergeContigs("", "ACGT", NoOverlap)
            .Should().Be("ACGT", "empty c1 ⇒ c2 (§6.1)");
        SequenceAssembler.MergeContigs("ACGT", "", NoOverlap)
            .Should().Be("ACGT", "empty c2 ⇒ c1 (§6.1)");
        SequenceAssembler.MergeContigs("", "", NoOverlap)
            .Should().BeEmpty("two empty contigs ⇒ empty result");

        // An empty contig has min length 0, so any positive overlap is out-of-range ⇒ concatenate.
        SequenceAssembler.MergeContigs("", "ACGT", 3)
            .Should().Be("ACGT", "overlap > min(0,4)=0 ⇒ concatenate (INV-03)");
    }

    // Fuzz: a contig fully contained as the prefix of a longer one, merged at l = |contained|,
    // ALWAYS yields exactly the container — no duplication of the contained sequence.
    [Test]
    public void MergeContigs_RandomFullContainment_AlwaysReturnsContainer()
    {
        var rng = new Random(144_002);
        for (int trial = 0; trial < 500; trial++)
        {
            string container = RandomString(rng, rng.Next(1, 50), DnaAlphabet);
            int containedLen = rng.Next(1, container.Length + 1);
            string contained = container.Substring(0, containedLen); // genuine prefix of container

            // contained ⊆ container at the front; overlap = |contained| = min(|container|,|contained|).
            string merged = SequenceAssembler.MergeContigs(container, contained, containedLen);

            merged.Should().Be(container, "a contained prefix is absorbed → just the container");
            AssertWellFormed(merged, container, contained, containedLen);
        }
    }

    #endregion

    #region BE — Boundary: IDENTICAL CONTIGS (dedup to a single copy, idempotent, deterministic)

    // Two identical contigs merged at the full overlap l = |c| de-duplicate to a SINGLE copy (= c),
    // NOT a doubled sequence; idempotent.
    [Test]
    public void MergeContigs_IdenticalContigsFullOverlap_DedupToSingleCopy()
    {
        const string c = "ACGTACGTAC";
        string merged = SequenceAssembler.MergeContigs(c, c, c.Length);

        merged.Should().Be(c, "identical contigs at full overlap collapse to one copy (no doubling)");
        merged.Should().HaveLength(c.Length, "INV-01: |c| + |c| − |c| = |c|");
        AssertWellFormed(merged, c, c, c.Length);
    }

    // Determinism: merging the SAME identical pair repeatedly yields the IDENTICAL result every time
    // (no non-determinism, no runaway growth / infinite merge).
    [Test]
    [CancelAfter(10_000)]
    public void MergeContigs_IdenticalContigs_DeterministicAndStable()
    {
        const string c = "GATTACAGATTACA";
        string first = SequenceAssembler.MergeContigs(c, c, c.Length);
        for (int i = 0; i < 1000; i++)
        {
            SequenceAssembler.MergeContigs(c, c, c.Length)
                .Should().Be(first, "the merge is a deterministic pure function (no state, no growth)");
        }
        first.Should().Be(c, "identical full-overlap merge is idempotent (= c)");
    }

    // Identical contigs at overlap 0 → the documented concatenation (a doubled sequence): proves the
    // dedup above is driven by the overlap length, not silent de-duplication.
    [Test]
    public void MergeContigs_IdenticalContigsZeroOverlap_ConcatenatesToDouble()
    {
        const string c = "ACGT";
        SequenceAssembler.MergeContigs(c, c, NoOverlap)
            .Should().Be("ACGTACGT", "overlap 0 ⇒ concatenate even for identical contigs (INV-02)");
    }

    // Fuzz: identical random contigs at full overlap ALWAYS dedup to one copy; never crash, never grow.
    [Test]
    public void MergeContigs_RandomIdenticalContigsFullOverlap_AlwaysSingleCopy()
    {
        var rng = new Random(144_003);
        for (int trial = 0; trial < 500; trial++)
        {
            string c = RandomString(rng, rng.Next(1, 60), WideAlphabet);

            string merged = SequenceAssembler.MergeContigs(c, c, c.Length);

            merged.Should().Be(c, "identical contigs at full overlap collapse to a single copy");
            AssertWellFormed(merged, c, c, c.Length);
        }
    }

    #endregion

    #region BE — Validation and broad fuzz

    // Null contigs → documented ArgumentNullException (§3.3, §6.1).
    [Test]
    public void MergeContigs_NullContigs_Throw()
    {
        ((Action)(() => SequenceAssembler.MergeContigs(null!, "ACGT", 0)))
            .Should().Throw<ArgumentNullException>("null contig1 is the documented validation contract (§3.3)");
        ((Action)(() => SequenceAssembler.MergeContigs("ACGT", null!, 0)))
            .Should().Throw<ArgumentNullException>("null contig2 is the documented validation contract (§3.3)");
    }

    // Broad fuzz: random contigs × random overlap (incl. negative, 0, in-range, far out-of-range,
    // int.MaxValue) NEVER throw an undocumented exception and ALWAYS match the documented rule
    // and the well-formedness invariants.
    [Test]
    [CancelAfter(30_000)]
    public void MergeContigs_RandomContigsAndOverlap_NeverThrows_MatchesOracle()
    {
        var rng = new Random(144_004);
        for (int trial = 0; trial < 5000; trial++)
        {
            string c1 = RandomString(rng, rng.Next(0, 40), WideAlphabet);
            string c2 = RandomString(rng, rng.Next(0, 40), WideAlphabet);

            // A spread of overlap lengths covering every contract branch.
            int l = rng.Next(6) switch
            {
                0 => NoOverlap,                              // no overlap
                1 => -rng.Next(1, 100),                      // negative
                2 => rng.Next(0, Math.Max(1, Math.Min(c1.Length, c2.Length) + 1)), // in/near range
                3 => Math.Min(c1.Length, c2.Length) + rng.Next(1, 50), // just out of range
                4 => int.MaxValue,                           // BE: MaxInt
                _ => rng.Next(int.MinValue, int.MaxValue),   // anything
            };

            string merged = SequenceAssembler.MergeContigs(c1, c2, l);

            merged.Should().Be(Oracle(c1, c2, l),
                "the unit matches the documented merge rule under fuzzed contigs/overlap");
            AssertWellFormed(merged, c1, c2, l);
        }
    }

    // Case sensitivity / alphabet-agnostic: the merge manipulates characters verbatim, performs
    // NO T↔U normalization and NO case folding (§3.3).
    [Test]
    public void MergeContigs_CaseSensitiveAlphabetAgnostic_NoNormalization()
    {
        // "acgt" suffix of c1 does NOT case-fold to match "ACGT"; with l = 0 it just concatenates.
        SequenceAssembler.MergeContigs("xxacgt", "ACGTyy", NoOverlap)
            .Should().Be("xxacgtACGTyy", "no case folding, no T↔U normalization (§3.3)");

        // A genuine same-case overlap collapses; the result preserves the exact characters.
        SequenceAssembler.MergeContigs("ggUUU", "UUUcc", 3)
            .Should().Be("ggUUUcc", "exact 'UUU' overlap collapses; U is treated as a literal char");
    }

    #endregion

    #endregion
}
