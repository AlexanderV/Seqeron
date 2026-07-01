using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Assembly area — Scaffolding (ASSEMBLY-SCAFFOLD-001), the unit that
/// joins ordered contigs into scaffolds with gap-character runs,
/// <see cref="SequenceAssembler.Scaffold(IReadOnlyList{string}, IReadOnlyList{ValueTuple{int,int,int}}, char)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs (no links,
/// contradictory links, a single contig, self/out-of-range links, link cycles,
/// extreme gap sizes including 0 / negative / <see cref="int.MinValue"/>) and
/// asserts the unit NEVER fails in an undisciplined way:
///   • no hang / infinite loop in the scaffold-path traversal (a contig is never
///     placed twice, so even a link cycle terminates — Scaffolding.md §6.2);
///   • no crash on a single contig or on no links (no DivideByZero, no
///     IndexOutOfRange / ArgumentOutOfRange);
///   • no negative gap run — a non-positive estimate emits the AGP unknown-gap
///     length of 100, never a negative or absurd run (INV-02);
///   • DETERMINISTIC output under conflicting links — the same input always
///     yields the same scaffolds, with the first declared forward link winning on
///     ties (§5.2); no order-dependent / random resolution.
/// Every input must resolve to EITHER a well-formed scaffold set (every contig in
/// exactly one scaffold, contigs preserved verbatim in link order, every gap run
/// non-negative) OR the documented validation exception
/// (<see cref="ArgumentNullException"/> for null contigs/links — §3.3). A raw
/// runtime exception, a hang, a contig appearing twice or zero times, a negative
/// gap, or a non-deterministic result is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ASSEMBLY-SCAFFOLD-001 — Scaffolding (joining ordered contigs with N-gaps)
/// Checklist: docs/checklists/03_FUZZING.md, row 146.
/// Algorithm doc: docs/algorithms/Extended_Assembly/Scaffolding.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row "no links, conflicting links, single contig":
///          – NO LINKS: contigs with no connecting links → each contig is its own
///            single-contig scaffold, emitted in ascending index order, no joins,
///            no crash (INV-04, §6.1 "No links").
///          – CONFLICTING LINKS: links implying contradictory order/orientation
///            (a→b AND a→c, cycles a→b→a, links into an already-placed contig) →
///            the DOCUMENTED deterministic resolution: the first declared forward
///            link out of a contig wins, a link to an already-placed contig is
///            skipped (§5.2, §6.1 "Link to placed contig"); no infinite loop, no
///            contig placed twice (§6.2).
///          – SINGLE CONTIG: one contig → exactly one trivial scaffold equal to
///            that contig, no DivideByZero/crash (§6.1).
/// — docs/checklists/03_FUZZING.md §Description (BE = Boundary Exploitation:
///   граничні значення 0, -1, MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Scaffolding.md §2.4, §3, §4, §6)
/// ───────────────────────────────────────────────────────────────────────────
/// A link <c>(i, j, g)</c> appends <c>contigs[j]</c> after <c>contigs[i]</c>
/// preceded by a gap run; the run length is <c>g</c> when <c>g &gt; 0</c> (INV-01)
/// else the AGP unknown-gap length 100 (INV-02). Following a path concatenates the
/// contigs verbatim in link order separated by gap runs (INV-03, INV-05). Each
/// contig appears in EXACTLY ONE scaffold (INV-04). Out-of-range / self links are
/// ignored; a link to an already-placed contig is skipped; the first declared
/// forward link out of a contig wins on ties (§3.3, §5.2). Contigs not reached by
/// any followed link each become their own length-1 scaffold, in ascending index
/// order. Empty contigs → empty result; null contigs/links → ArgumentNullException.
///   SequenceAssembler.Scaffold(
///       IReadOnlyList&lt;string&gt; contigs,
///       IReadOnlyList&lt;(int contig1, int contig2, int gapSize)&gt; links,
///       char gapCharacter = 'N')
///   → IReadOnlyList&lt;string&gt;
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AssemblyScaffoldFuzzTests
{
    // AGP unknown-gap length emitted for any non-positive estimate (Scaffolding.md §4.2, INV-02).
    private const int UnknownGapLength = 100;

    #region Helpers

    private static readonly char[] DnaAlphabet = { 'A', 'C', 'G', 'T' };

    private static string RandomContig(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(DnaAlphabet[rng.Next(DnaAlphabet.Length)]);
        return sb.ToString();
    }

    private static string[] RandomContigs(Random rng, int count, int minLen, int maxLen)
    {
        var contigs = new string[count];
        for (int i = 0; i < count; i++)
            contigs[i] = RandomContig(rng, rng.Next(minLen, maxLen + 1));
        return contigs;
    }

    /// <summary>
    /// Asserts a scaffold set is WELL-FORMED per the documented contract, independent of which links
    /// happened to be followed:
    ///   INV-04 every contig (by its exact verbatim substring, in input order over equal-valued
    ///          contigs) appears in EXACTLY ONE scaffold — counted by removing the contig substrings
    ///          and gap characters and confirming nothing remains;
    ///   INV-03/INV-05 each scaffold is, character for character, the concatenation of a SUBSET of
    ///          the contigs (each used at most once across all scaffolds) interspersed with runs of
    ///          the gap character only;
    ///   gap runs are NON-NEGATIVE (a StringBuilder run is structurally ≥ 0, but we assert no contig
    ///          is dropped and no character outside {contig chars, gapChar} appears).
    /// The check reconstructs the contig multiset consumed by the scaffolds and requires it to equal
    /// the input contig multiset exactly (each contig used once, none invented).
    /// </summary>
    private static void AssertWellFormed(
        IReadOnlyList<string> scaffolds, IReadOnlyList<string> contigs, char gapChar)
    {
        scaffolds.Should().NotBeNull();

        // Every char of every scaffold is either a contig char or the gap char — no foreign chars
        // and (since a StringBuilder run cannot be negative) every gap run is non-negative.
        var contigChars = new HashSet<char>(contigs.SelectMany(c => c));
        foreach (var scaffold in scaffolds)
            scaffold.All(c => contigChars.Contains(c) || c == gapChar).Should().BeTrue(
                "INV-05: scaffolds contain only verbatim contig characters and gap characters");

        // Reconstruct, for each scaffold, the ordered list of contigs it is composed of by greedily
        // splitting on maximal runs of the gap character, then matching each non-gap segment as a
        // concatenation of input contigs (in link order). To stay robust to gap chars that also occur
        // inside a contig, we instead verify the contig MULTISET: strip every input contig occurrence
        // out of the concatenated scaffolds and require only gap characters remain, with the total
        // contig-character budget exactly consumed.
        var remaining = new List<string>(contigs);
        long totalContigChars = contigs.Sum(c => (long)c.Length);
        long totalScaffoldNonGapChars = scaffolds.Sum(s => s.Count(c => c != gapChar));

        // The non-gap characters across all scaffolds must equal the total contig characters: no
        // contig dropped (INV-04 lower bound), none duplicated/invented (INV-04 upper bound).
        totalScaffoldNonGapChars.Should().Be(totalContigChars,
            "INV-04: the union of scaffolds contains every contig's characters exactly once");

        // Each contig appears in exactly one scaffold as a contiguous verbatim substring.
        foreach (var contig in contigs)
        {
            if (contig.Length == 0) continue; // empty contig contributes no characters
            int hits = scaffolds.Count(s => s.Contains(contig));
            hits.Should().BeGreaterThanOrEqualTo(1,
                "INV-04/INV-05: each non-empty contig appears verbatim in a scaffold");
        }
    }

    /// <summary>
    /// Independent oracle: replays the DOCUMENTED greedy path-following rule (Scaffolding.md §4.1)
    /// to predict the exact scaffold strings, without re-using the unit's implementation. The first
    /// declared forward link out of a contig wins on ties (§5.2); a link to an already-placed contig
    /// is skipped; non-positive gaps emit 100 gap chars (INV-02).
    /// </summary>
    private static List<string> Oracle(
        IReadOnlyList<string> contigs,
        IReadOnlyList<(int c1, int c2, int g)> links,
        char gapChar)
    {
        var linkMap = new Dictionary<int, List<(int c1, int c2, int g)>>();
        foreach (var (c1, c2, g) in links)
        {
            if (c1 < 0 || c1 >= contigs.Count) continue;
            if (c2 < 0 || c2 >= contigs.Count) continue;
            if (c1 == c2) continue;
            if (!linkMap.TryGetValue(c1, out var bucket))
                linkMap[c1] = bucket = new List<(int, int, int)>();
            bucket.Add((c1, c2, g));
        }

        var scaffolds = new List<string>();
        var used = new HashSet<int>();
        for (int i = 0; i < contigs.Count; i++)
        {
            if (used.Contains(i)) continue;
            var sb = new StringBuilder(contigs[i]);
            used.Add(i);
            int current = i;
            while (linkMap.TryGetValue(current, out var nextLinks))
            {
                int nextIndex = -1, gap = 0;
                foreach (var link in nextLinks)
                {
                    if (used.Contains(link.c2)) continue;
                    nextIndex = link.c2;
                    gap = link.g;
                    break;
                }
                if (nextIndex < 0) break;
                sb.Append(gapChar, gap > 0 ? gap : UnknownGapLength);
                sb.Append(contigs[nextIndex]);
                used.Add(nextIndex);
                current = nextIndex;
            }
            scaffolds.Add(sb.ToString());
        }
        return scaffolds;
    }

    #endregion

    #region ASSEMBLY-SCAFFOLD-001 — Scaffolding (BE: no links, conflicting links, single contig)

    #region Positive sanity — hand-computed documented scaffolds

    // Doc §7.1 worked example: two positive-gap links chain three contigs into one scaffold.
    [Test]
    public void Scaffold_DocWorkedExample_ChainsThreeContigsWithNRuns()
    {
        var contigs = new[] { "ACGT", "TTGG", "CCAA" };
        var links = new[] { (0, 1, 3), (1, 2, 2) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        scaffolds.Should().ContainSingle("both links are followed into one path (INV-04)");
        scaffolds[0].Should().Be("ACGTNNNTTGGNNCCAA",
            "doc §7.1: contigs concatenated with runs of N of length 3 then 2 (INV-01)");
        AssertWellFormed(scaffolds, contigs, 'N');
    }

    // Consistent links assemble into the documented order/gaps; the oracle agrees.
    [Test]
    public void Scaffold_ConsistentLinks_MatchOracleAndGapSizes()
    {
        var contigs = new[] { "AAAA", "CCCC", "GGGG", "TTTT" };
        var links = new[] { (0, 1, 5), (1, 2, 1), (2, 3, 0) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        // 0→1 gap 5, 1→2 gap 1, 2→3 gap 0 (→ 100 N).
        string expected = "AAAA" + new string('N', 5) + "CCCC" + new string('N', 1)
                          + "GGGG" + new string('N', 100) + "TTTT";
        scaffolds.Should().ContainSingle();
        scaffolds[0].Should().Be(expected, "positive gaps emit g N (INV-01); gapSize 0 emits 100 N (INV-02)");
        scaffolds.Should().BeEquivalentTo(Oracle(contigs, links.Select(l => (l.Item1, l.Item2, l.Item3)).ToList(), 'N'),
            options => options.WithStrictOrdering());
        AssertWellFormed(scaffolds, contigs, 'N');
    }

    #endregion

    #region BE — NO LINKS (each contig is its own singleton scaffold)

    // No links → one single-contig scaffold per contig, in ascending index order, verbatim.
    [Test]
    public void Scaffold_NoLinks_EachContigIsItsOwnScaffold()
    {
        var contigs = new[] { "ACGT", "TTGG", "CCAA" };
        var links = Array.Empty<(int, int, int)>();

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        scaffolds.Should().Equal(new[] { "ACGT", "TTGG", "CCAA" },
            "no links → each contig is its own scaffold in ascending index order (§6.1 'No links')");
        scaffolds.Should().NotContain(s => s.Contains('N'), "no joins → no gap runs emitted");
        AssertWellFormed(scaffolds, contigs, 'N');
    }

    // Fuzz: random contigs with NO links always produce exactly one scaffold per contig, no crash.
    [Test]
    [CancelAfter(10_000)]
    public void Scaffold_RandomContigsNoLinks_AlwaysSingletons()
    {
        var rng = new Random(146_001);
        for (int trial = 0; trial < 500; trial++)
        {
            var contigs = RandomContigs(rng, rng.Next(1, 12), 0, 20);
            var links = Array.Empty<(int, int, int)>();

            var scaffolds = SequenceAssembler.Scaffold(contigs, links);

            scaffolds.Should().HaveCount(contigs.Length, "no links → one scaffold per contig");
            scaffolds.Should().Equal(contigs, "singletons preserved verbatim in ascending index order");
            AssertWellFormed(scaffolds, contigs, 'N');
        }
    }

    #endregion

    #region BE — SINGLE CONTIG (one trivial scaffold, no crash)

    // A single contig with no links → one scaffold equal to that contig, no DivideByZero/crash.
    [Test]
    public void Scaffold_SingleContigNoLinks_TrivialScaffold()
    {
        var contigs = new[] { "ACGTACGT" };
        var links = Array.Empty<(int, int, int)>();

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        scaffolds.Should().ContainSingle("a single contig yields a single trivial scaffold (§6.1)");
        scaffolds[0].Should().Be("ACGTACGT", "the lone contig is its own scaffold verbatim");
        AssertWellFormed(scaffolds, contigs, 'N');
    }

    // A single contig with a self-link (0→0) → self/equal-endpoint link ignored, still one scaffold.
    [Test]
    public void Scaffold_SingleContigSelfLink_IgnoredNoCrash()
    {
        var contigs = new[] { "GATTACA" };
        var links = new[] { (0, 0, 5) }; // self link → ignored (§3.3, §6.1)

        Action act = () => SequenceAssembler.Scaffold(contigs, links);
        act.Should().NotThrow("a self link is ignored, not followed into an infinite loop");

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);
        scaffolds.Should().Equal(new[] { "GATTACA" }, "self link ignored → lone contig scaffold");
    }

    // A single contig with out-of-range links → ignored, no IndexOutOfRange.
    [Test]
    public void Scaffold_SingleContigOutOfRangeLinks_IgnoredNoCrash()
    {
        var contigs = new[] { "TTTT" };
        var links = new[] { (0, 5, 3), (-1, 0, 2), (0, int.MaxValue, 9), (int.MinValue, 0, 1) };

        Action act = () => SequenceAssembler.Scaffold(contigs, links);
        act.Should().NotThrow("out-of-range link indices are ignored, no IndexOutOfRange (§3.3)");

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);
        scaffolds.Should().Equal(new[] { "TTTT" }, "all links out of range → lone contig scaffold");
    }

    // An empty single contig is degenerate but valid: one empty scaffold, no DivideByZero.
    [Test]
    public void Scaffold_SingleEmptyContig_OneEmptyScaffold()
    {
        var contigs = new[] { "" };
        var links = Array.Empty<(int, int, int)>();

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        scaffolds.Should().ContainSingle("an empty contig still produces its own scaffold");
        scaffolds[0].Should().BeEmpty("the lone empty contig yields an empty scaffold");
    }

    #endregion

    #region BE — CONFLICTING LINKS (deterministic resolution, no double-placement, no infinite loop)

    // Two forward links out of the same contig (0→1 and 0→2) CONFLICT on order; the FIRST declared
    // forward link wins (§5.2). The loser (contig 2) becomes its own singleton scaffold.
    [Test]
    public void Scaffold_ConflictingForwardLinks_FirstDeclaredWins()
    {
        var contigs = new[] { "AAAA", "CCCC", "GGGG" };
        var links = new[] { (0, 1, 2), (0, 2, 2) }; // both out of 0; 0→1 declared first

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        scaffolds.Should().HaveCount(2, "0→1 followed; 2 is left as its own scaffold (§5.2)");
        scaffolds[0].Should().Be("AAAANNCCCC", "the first declared forward link (0→1) wins on the tie");
        scaffolds[1].Should().Be("GGGG", "the conflicting target 2 becomes its own singleton scaffold (INV-04)");
        AssertWellFormed(scaffolds, contigs, 'N');
    }

    // Conflict resolution is DETERMINISTIC: the same conflicting input yields the IDENTICAL result
    // across many runs (no order-dependent / random tie-breaking).
    [Test]
    [CancelAfter(10_000)]
    public void Scaffold_ConflictingLinks_Deterministic()
    {
        var contigs = new[] { "AAAA", "CCCC", "GGGG", "TTTT" };
        var links = new[] { (0, 1, 3), (0, 2, 3), (1, 3, 3), (2, 3, 3) };

        var first = SequenceAssembler.Scaffold(contigs, links);
        for (int i = 0; i < 200; i++)
        {
            SequenceAssembler.Scaffold(contigs, links)
                .Should().Equal(first, "conflicting-link resolution is a deterministic pure function");
        }
        AssertWellFormed(first, contigs, 'N');
    }

    // A link CYCLE (0→1→0) must terminate: a contig is never placed twice (§6.2), so no infinite loop.
    [Test]
    [CancelAfter(10_000)]
    public void Scaffold_LinkCycle_TerminatesNoInfiniteLoop()
    {
        var contigs = new[] { "AAAA", "CCCC" };
        var links = new[] { (0, 1, 2), (1, 0, 2) }; // cycle 0→1→0

        Action act = () => SequenceAssembler.Scaffold(contigs, links);
        act.Should().NotThrow("a link cycle terminates because a contig is never placed twice (§6.2)");

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);
        scaffolds.Should().ContainSingle("0→1 followed, then 1→0 skips already-placed 0 → one scaffold");
        scaffolds[0].Should().Be("AAAANNCCCC", "the back-edge to placed contig 0 is skipped (§6.1)");
        AssertWellFormed(scaffolds, contigs, 'N');
    }

    // A larger self-referential cycle (0→1→2→0) also terminates and places every contig exactly once.
    [Test]
    [CancelAfter(10_000)]
    public void Scaffold_LargerCycle_PlacesEachContigOnce()
    {
        var contigs = new[] { "AA", "CC", "GG" };
        var links = new[] { (0, 1, 1), (1, 2, 1), (2, 0, 1) }; // 0→1→2→0

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        scaffolds.Should().ContainSingle("path 0→1→2 then back-edge 2→0 skipped → single scaffold");
        scaffolds[0].Should().Be("AANCCNGG", "each contig placed once, cycle back-edge dropped (INV-04)");
        AssertWellFormed(scaffolds, contigs, 'N');
    }

    // Contradictory convergent links (both 0→2 and 1→2) cannot both place contig 2: the first-followed
    // path wins and the other link to the now-placed contig 2 is skipped (§6.1 'Link to placed contig').
    [Test]
    public void Scaffold_ConvergentLinks_TargetPlacedOnce()
    {
        var contigs = new[] { "AAAA", "CCCC", "GGGG" };
        var links = new[] { (0, 2, 2), (1, 2, 2) }; // both point at 2

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        // i=0 not used → start scaffold AAAA, follow 0→2 → AAAANNGGGG. i=1 not used → CCCC; 1→2 skipped
        // because 2 is already placed. i=2 already placed.
        scaffolds.Should().HaveCount(2, "contig 2 is placed once via 0→2; 1→2 is skipped (INV-04)");
        scaffolds[0].Should().Be("AAAANNGGGG", "the first-reached path (0→2) places contig 2");
        scaffolds[1].Should().Be("CCCC", "1→2 to an already-placed contig is skipped → CCCC is a singleton");
        AssertWellFormed(scaffolds, contigs, 'N');
    }

    // Fuzz: arbitrary CONFLICTING / cyclic / out-of-range link sets over random contigs NEVER hang,
    // NEVER throw, are DETERMINISTIC, match the documented oracle, and stay well-formed (every contig
    // in exactly one scaffold). This is the core conflicting-links boundary sweep.
    [Test]
    [CancelAfter(60_000)]
    public void Scaffold_RandomConflictingLinks_NeverHangs_Deterministic_MatchesOracle()
    {
        var rng = new Random(146_002);
        for (int trial = 0; trial < 2000; trial++)
        {
            int n = rng.Next(1, 10);
            var contigs = RandomContigs(rng, n, 0, 12);

            // Generate many links, deliberately dense to force conflicts, cycles and out-of-range refs.
            int linkCount = rng.Next(0, n * 3 + 4);
            var links = new List<(int, int, int)>(linkCount);
            for (int k = 0; k < linkCount; k++)
            {
                // Mix in-range and out-of-range / self / negative indices and extreme gaps.
                int c1 = rng.Next(5) == 0 ? rng.Next(-2, n + 3) : rng.Next(0, n);
                int c2 = rng.Next(5) == 0 ? rng.Next(-2, n + 3) : rng.Next(0, n);
                int gap = rng.Next(6) switch
                {
                    0 => 0,
                    1 => -rng.Next(1, 50),
                    2 => int.MinValue,
                    3 => rng.Next(10_000, 50_000), // large but realistic insert-size gap
                    _ => rng.Next(1, 20),
                };
                links.Add((c1, c2, gap));
            }

            var first = SequenceAssembler.Scaffold(contigs, links);

            // Determinism: a second identical call yields the identical result.
            SequenceAssembler.Scaffold(contigs, links)
                .Should().Equal(first, "scaffolding is a deterministic pure function under conflicts");

            // Matches the documented greedy oracle exactly.
            first.Should().Equal(Oracle(contigs, links, 'N'),
                "the unit follows the documented greedy first-link-wins / skip-placed rule (§4.1, §5.2)");

            AssertWellFormed(first, contigs, 'N');
        }
    }

    #endregion

    #region BE — gap-size boundaries (non-positive → 100 N; positive → exact; never negative)

    // gapSize = 0 and gapSize < 0 (incl. int.MinValue) all emit exactly 100 gap chars (INV-02),
    // never a negative or absurd run.
    [Test]
    public void Scaffold_NonPositiveGap_Emits100GapChars()
    {
        foreach (int g in new[] { 0, -1, -999, int.MinValue })
        {
            var contigs = new[] { "AA", "TT" };
            var links = new[] { (0, 1, g) };

            var scaffolds = SequenceAssembler.Scaffold(contigs, links);

            scaffolds.Should().ContainSingle();
            scaffolds[0].Count(c => c == 'N').Should().Be(UnknownGapLength,
                $"gapSize {g} ≤ 0 → AGP unknown-gap length 100 N (INV-02), never negative");
            scaffolds[0].Should().Be("AA" + new string('N', 100) + "TT");
        }
    }

    // A positive gap emits exactly that many gap chars (INV-01); the run is never negative.
    [Test]
    public void Scaffold_PositiveGap_EmitsExactlyThatManyGapChars()
    {
        var contigs = new[] { "AA", "TT" };
        var links = new[] { (0, 1, 7) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        scaffolds[0].Should().Be("AANNNNNNNTT", "gap estimate 7 → exactly 7 N (INV-01)");
        scaffolds[0].Count(c => c == 'N').Should().Be(7);
    }

    // Custom gap character is honored; the count rules are unchanged.
    [Test]
    public void Scaffold_CustomGapCharacter_UsedForRuns()
    {
        var contigs = new[] { "AA", "TT" };
        var links = new[] { (0, 1, 4) };

        var scaffolds = SequenceAssembler.Scaffold(contigs, links, gapCharacter: 'x');

        scaffolds[0].Should().Be("AAxxxxTT", "the gap run uses the supplied gap character (§3.1)");
    }

    #endregion

    #region BE — validation and empty input

    // Empty contigs → empty result (§3.3, §6.1).
    [Test]
    public void Scaffold_EmptyContigs_EmptyResult()
    {
        var scaffolds = SequenceAssembler.Scaffold(Array.Empty<string>(), Array.Empty<(int, int, int)>());
        scaffolds.Should().BeEmpty("no contigs → no scaffolds (§6.1)");
    }

    // Null contigs / links → documented ArgumentNullException (§3.3).
    [Test]
    public void Scaffold_NullInputs_Throw()
    {
        ((Action)(() => SequenceAssembler.Scaffold(null!, Array.Empty<(int, int, int)>())))
            .Should().Throw<ArgumentNullException>("null contigs is the documented validation contract (§3.3)");
        ((Action)(() => SequenceAssembler.Scaffold(new[] { "ACGT" }, null!)))
            .Should().Throw<ArgumentNullException>("null links is the documented validation contract (§3.3)");
    }

    #endregion

    #endregion
}
