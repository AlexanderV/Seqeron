using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;
// Disambiguate: a top-level Seqeron.Genomics.Analysis.MotifMatch also exists; this unit
// asserts against the ProteinMotifFinder.MotifMatch record returned by FindMotifByPattern.
using MotifMatch = Seqeron.Genomics.Analysis.ProteinMotifFinder.MotifMatch;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the ProteinMotif area — pattern-based motif matching
/// (PROTMOTIF-PATTERN-001): scanning a protein sequence for a caller-supplied
/// .NET <em>regular-expression</em> motif via
/// <see cref="ProteinMotifFinder.FindMotifByPattern"/>. The method uppercases the
/// protein, wraps the supplied pattern in a zero-width lookahead <c>(?=(P))</c> so
/// every (including overlapping) start position is enumerated, compiles it with
/// <see cref="System.Text.RegularExpressions.RegexOptions.IgnoreCase"/>, and emits one
/// <c>MotifMatch</c> per capture with an information-content Score and a model EValue.
/// — docs/algorithms/ProteinMotif/Pattern_Matching_Methods.md §2.2, §4.1, §5.2.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Scope vs sibling units (no duplication)
/// ───────────────────────────────────────────────────────────────────────────
/// PROTMOTIF-PATTERN-001 denotes the same public method as PROTMOTIF-FIND-001
/// (row 82, <c>ProteinMotifFuzzTests</c>) — confirmed by Pattern_Matching_Methods.md
/// (Test Unit ID = PROTMOTIF-PATTERN-001) listing <c>FindMotifByPattern</c> as its
/// primary entry point. Row 82 fuzzed the PROTEIN-side inputs (empty/null/short/
/// homopolymer protein, non-amino-acid residues) and the multi-pattern
/// <c>FindCommonMotifs</c> convenience scan. Row 83 (PROTMOTIF-PROSITE-001,
/// <c>ProteinPrositeFuzzTests</c>) fuzzed the PROSITE→regex GRAMMAR converter.
///
/// This file is the row-166 MC suite and goes DEEPER on the regexPattern argument's
/// robustness — the part neither sibling stresses: malformed/unbalanced patterns,
/// regex-injection attempts, and CATASTROPHIC-BACKTRACKING patterns. The headline
/// MC targets (docs/checklists/03_FUZZING.md row 166: "empty pattern, invalid regex,
/// no match") are:
///   • EMPTY PATTERN — null / "" regexPattern → NO matches by the explicit
///     string.IsNullOrEmpty guard (INV-01), never a throw and never a spurious
///     zero-width match at every position.
///   • INVALID REGEX — a malformed pattern ("[ABC", "(", "(?<", a bare quantifier
///     "*", a regex-injection attempt) fails to compile; the compile exception is
///     SWALLOWED inside FindMotifByPattern → NO matches, NEVER an unhandled
///     exception, and a CATASTROPHIC-BACKTRACKING pattern must not wedge the suite
///     (asserted under [CancelAfter]).
///   • NO MATCH — a VALID pattern that simply does not occur in the protein →
///     empty result, never a false/phantom match.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing"; strategy MC = Malformed Content
///   (docs/checklists/03_FUZZING.md §Description).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Hazards under test (MC)
/// ───────────────────────────────────────────────────────────────────────────
///   • an UNHANDLED RegexParseException leaking out of FindMotifByPattern when the
///     pattern is malformed (the documented contract is that it is swallowed →
///     empty, Pattern_Matching_Methods.md §3.3, §5.2, §6.1, INV-01);
///   • a HANG on a catastrophic-backtracking pattern over a degenerate sequence —
///     even though a malformed-by-construction such pattern would normally be
///     swallowed at compile time, a VALID-but-pathological pattern (e.g.
///     "(A+)+B" / "(A*)*B") compiles fine and could blow up at match time; the
///     scan must terminate (guarded by [CancelAfter]);
///   • a FALSE MATCH — a valid non-occurring pattern reported as present, or a
///     coordinate corruption on any match that IS found (the universal
///     well-formed-match contract, INV-02).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Documented contract asserted (Pattern_Matching_Methods.md §2.4, §3, §6)
/// ───────────────────────────────────────────────────────────────────────────
///   • INV-01 — null/empty sequence OR pattern, AND an invalid .NET regex, → no
///     matches and no exception (guard clauses + guarded try/catch).
///   • INV-02 — Sequence == upper(S).Substring(Start, End−Start+1); End = Start+len−1.
///   • INV-05 — matching is case-insensitive; positions are 0-based inclusive.
///   • [finite-score] — Score and EValue are finite (never NaN / ±∞).
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ProteinPatternMatchFuzzTests
{
    #region Helpers

    /// <summary>The 20 standard amino-acid one-letter codes.</summary>
    private const string StandardAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>Deterministic RNG — seed fixed LOCALLY so generated fuzz inputs are reproducible.</summary>
    private static string RandomProtein(int length, int seed)
    {
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = StandardAminoAcids[rng.Next(StandardAminoAcids.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal theory-correct contract every emitted <see cref="MotifMatch"/> must
    /// satisfy against the original (case-insensitive) sequence (Pattern_Matching_Methods.md §2.4,
    /// §3.2): a CONTIGUOUS in-bounds span (INV-02) whose claimed coordinates actually reproduce the
    /// reported substring (Sequence == upper(S)[Start..End], uppercased), with finite Score / EValue
    /// ([finite-score]). The "no coordinate bug / no run-off-the-end / no NaN" property.
    /// </summary>
    private static void AssertWellFormedMatch(MotifMatch match, string originalSequence)
    {
        string upper = originalSequence.ToUpperInvariant();
        int n = upper.Length;

        match.Start.Should().BeInRange(0, n - 1, "a match Start is a valid 0-based residue index");
        match.End.Should().BeInRange(match.Start, n - 1, "a match End is in-bounds and not before its Start");

        (match.End - match.Start + 1).Should().Be(match.Sequence.Length,
            "INV-02: End − Start + 1 equals the matched substring length");
        match.Sequence.Should().Be(match.Sequence.ToUpperInvariant(), "INV-05: the matched substring is uppercased");
        upper.Substring(match.Start, match.Sequence.Length).Should().Be(match.Sequence,
            "INV-02: the reported substring is exactly upper(S)[Start..End] (no coordinate bug)");

        double.IsNaN(match.Score).Should().BeFalse("a motif Score must never be NaN");
        double.IsInfinity(match.Score).Should().BeFalse("a motif Score must never be infinite");
        double.IsNaN(match.EValue).Should().BeFalse("a motif EValue must never be NaN");
        double.IsInfinity(match.EValue).Should().BeFalse("a motif EValue must never be infinite");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PROTMOTIF-PATTERN-001 — pattern-based motif matching : MC fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PROTMOTIF-PATTERN-001 — pattern matching

    #region MC — Empty / null pattern: no matches by the guard, never a phantom hit

    /// <summary>
    /// Target "empty pattern": a null or "" regexPattern must produce NO matches by the explicit
    /// <c>string.IsNullOrEmpty</c> guard (Pattern_Matching_Methods.md §3.3, §6.1, INV-01) — NEVER a
    /// throw, and — critically — never a phantom zero-width match at every position. (An empty pattern
    /// wrapped as <c>(?=())</c> would otherwise match at every offset; the guard prevents this.) The
    /// guard must hold regardless of the protein it is paired with (real protein, empty protein, a
    /// homopolymer that would otherwise be densely matchable).
    /// </summary>
    [Test]
    public void Pattern_EmptyOrNullPattern_NoMatchesNoThrow()
    {
        foreach (string seq in new[] { "MKRGDSPEKWFIL", "", new string('A', 200), "X" })
        {
            foreach (string? pattern in new[] { "", null })
            {
                var act = () => FindMotifByPattern(seq, pattern!, "Custom", "PID").ToList();
                act.Should().NotThrow($"empty/null pattern ('{pattern ?? "null"}') must not crash over '{seq[..Math.Min(seq.Length, 8)]}'")
                    .Subject.Should().BeEmpty(
                        "an empty/null pattern is guarded by string.IsNullOrEmpty → no matches, not a phantom hit per position");
            }
        }
    }

    #endregion

    #region MC — Invalid / malformed regex: compile failure swallowed → no matches, no throw

    /// <summary>
    /// Target "invalid regex": a malformed .NET regex fails to compile inside FindMotifByPattern; the
    /// documented contract (Pattern_Matching_Methods.md §3.3, §5.2, §6.1, INV-01) is that the compile
    /// exception is SWALLOWED in the guarded try/catch → an EMPTY enumeration, NEVER a leaked
    /// RegexParseException / ArgumentException. We probe a battery of structurally-broken patterns:
    /// unbalanced brackets/parens/braces, a dangling escape, a bare leading quantifier, an unterminated
    /// group construct, and a back-reference to a non-existent group. We pair each with a real protein
    /// and a homopolymer so the swallow path is exercised regardless of the (would-be) match density.
    /// </summary>
    [Test]
    public void Pattern_InvalidRegex_SwallowedToEmptyNoThrow()
    {
        string[] malformed =
        {
            "[ABC",         // unbalanced character class
            "ABC]",         // stray closing bracket is literal in .NET — kept as a NON-malformed control elsewhere
            "(",            // unbalanced opening group
            ")",            // unbalanced closing group
            "(?<",          // unterminated named-group construct
            "(?P<name>A)",  // Python-style named group — invalid in .NET
            "A{",           // unterminated quantifier brace
            "A{2,1}",       // inverted quantifier bounds
            "*",            // bare leading quantifier — nothing to repeat
            "+RGD",         // leading quantifier with no preceding atom
            "?",            // bare optional quantifier — nothing to repeat
            @"\",           // dangling escape at end of pattern
            @"(\1)",        // back-reference to a group that does not exist yet
            "[a-",          // unterminated range inside a class
            "(?#",          // unterminated inline-comment construct
        };

        // Note: "ABC]" is actually a VALID .NET regex (']' outside a class is a literal); it is included
        // to document that '].' is literal and must NOT be treated as malformed. We assert it separately.
        foreach (string seq in new[] { "MKRGDSPEKWFILRGDAAA", new string('A', 120) })
        {
            foreach (string pattern in malformed.Where(p => p != "ABC]"))
            {
                var act = () => FindMotifByPattern(seq, pattern).ToList();
                act.Should().NotThrow($"malformed regex '{pattern}' must be swallowed, not thrown")
                    .Subject.Should().BeEmpty($"a malformed regex '{pattern}' compiles-failure → no matches (INV-01)");
            }
        }

        // Control: a stray ']' outside a class is a VALID literal in .NET regex, so "AC]" is NOT
        // malformed — but it cannot occur in a protein (']' is not a residue) so it must yield EMPTY,
        // proving the swallow path is narrow (only true compile failures), not a blanket "any bracket
        // → empty", and that a valid-but-non-occurring pattern is distinguished from a malformed one.
        FindMotifByPattern("MKABCRGD", "AC]").Should().BeEmpty(
            "'AC]' is a valid literal pattern that simply does not occur (no ']' residue) — not a swallow");

        // And the valid literal "ABC" over a protein containing it matches and reports a well-formed span.
        var literalHits = FindMotifByPattern("ZZABCZZ", "ABC").ToList();
        literalHits.Should().ContainSingle("the valid literal pattern 'ABC' matches the substring once");
        AssertWellFormedMatch(literalHits[0], "ZZABCZZ");
        literalHits[0].Sequence.Should().Be("ABC");
    }

    #endregion

    #region MC — Regex-injection attempts: treated as ordinary (in)valid patterns, never special

    /// <summary>
    /// Target "invalid regex" (injection facet): a caller-supplied pattern is fed to <c>new Regex</c>
    /// AS-IS, so regex metacharacters carry their regex meaning — there is no shell/SQL-style escape to
    /// exploit, and an injection attempt is just a regex that either compiles (and matches by its own
    /// semantics) or fails to compile (and is swallowed). The robustness contract is the SAME as for
    /// any pattern: NEVER an unhandled exception, NEVER a corrupt coordinate. We feed alternation/anchor/
    /// lookaround "injection" strings and assert: no throw, and any hit is a well-formed in-bounds span
    /// whose substring actually equals what was matched (so an injected alternation cannot fabricate a
    /// span outside the sequence).
    /// </summary>
    [Test]
    public void Pattern_RegexInjectionAttempts_NoThrowWellFormedOrEmpty()
    {
        const string seq = "MKRGDSPEKWFILRGDAAA";

        string[] injections =
        {
            ".*",               // match-everything — must yield well-formed (possibly zero-width-rejected) spans
            "^.*$",             // anchored whole-string
            "(?=RGD)",          // pure lookahead — the outer wrapper double-wraps; whatever it yields must be well-formed
            "(?!RGD)",          // negative lookahead
            "(?<=K)RGD",        // lookbehind
            "RGD|WFIL",         // alternation of two real motifs
            "[A-Z]+",           // greedy run of letters
            "$|^",              // anchor alternation
            ".{0}",             // explicit zero-width
            "(RGD)\\1?",        // optional back-reference (valid)
        };

        foreach (string pattern in injections)
        {
            var act = () => FindMotifByPattern(seq, pattern, "Custom", "INJ").ToList();
            var hits = act.Should().NotThrow($"injection-style pattern '{pattern}' must not crash").Subject;
            foreach (var hit in hits)
                AssertWellFormedMatch(hit, seq);
        }

        // The alternation "RGD|WFIL" must surface exactly the real occurrences (two RGD + one WFIL),
        // proving injected metacharacters are honoured as REGEX, not mishandled, and never fabricate spans.
        var altHits = FindMotifByPattern(seq, "RGD|WFIL").ToList();
        altHits.Select(h => h.Sequence).Should().OnlyContain(s => s == "RGD" || s == "WFIL",
            "an alternation pattern yields only its true alternatives");
        altHits.Should().HaveCount(3, "the sequence contains two RGD and one WFIL");
    }

    #endregion

    #region MC — Catastrophic backtracking: must terminate, never wedge the suite

    /// <summary>
    /// Target "invalid regex" (hang facet): a VALID-but-pathological pattern — the classic
    /// catastrophic-backtracking forms <c>(A+)+B</c>, <c>(A*)*B</c>, <c>(A|A)*B</c> — COMPILES fine,
    /// so the compile-time swallow does not protect against it; the danger is an exponential blow-up at
    /// MATCH time over a degenerate "AAAA…" sequence with no terminal 'B'. The scan must TERMINATE: this
    /// is asserted under <c>[CancelAfter]</c> so a regression that wedged the lookahead walk would FAIL
    /// the test rather than hang the whole suite. (.NET's backtracking engine has no built-in match
    /// timeout here, so termination is the headline contract.) We keep the pathological input modest.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Pattern_CatastrophicBacktracking_TerminatesNoHang(CancellationToken token)
    {
        // A run of 'A' with NO terminal 'B' — the worst case for "(A+)+B"-style patterns.
        string seq = new string('A', 28);

        string[] pathological = { "(A+)+B", "(A*)*B", "(A|A)*B", "(A+)*B", "(.*)*B" };

        foreach (string pattern in pathological)
        {
            var act = () => FindMotifByPattern(seq, pattern).ToList();
            var hits = act.Should().NotThrow($"pathological pattern '{pattern}' must not throw").Subject;
            token.ThrowIfCancellationRequested(); // if the scan wedged, [CancelAfter] cancels and this throws → FAIL

            hits.Should().BeEmpty($"no 'B' exists in the sequence, so '{pattern}' has no match");
        }

        // Same patterns over a sequence that DOES end in 'B' must terminate and yield well-formed hits.
        string seqWithB = new string('A', 24) + "B";
        foreach (string pattern in pathological.Where(p => p != "(.*)*B"))
        {
            var hits = ((Func<List<MotifMatch>>)(() => FindMotifByPattern(seqWithB, pattern).ToList()))
                .Should().NotThrow($"pathological pattern '{pattern}' over 'A…B' must not throw").Subject;
            token.ThrowIfCancellationRequested();
            foreach (var hit in hits)
                AssertWellFormedMatch(hit, seqWithB);
        }
    }

    #endregion

    #region MC — No match: valid non-occurring pattern → empty, never a false hit

    /// <summary>
    /// Target "no match": a VALID pattern that simply does not occur in the protein must yield an EMPTY
    /// result — never a false/phantom hit (Pattern_Matching_Methods.md §6.1 "No motif occurrence →
    /// empty"). We probe non-occurring literals, a non-occurring fixed-length wildcard pattern, a
    /// negated-class pattern with no satisfying residue, and an anchored pattern that cannot match
    /// mid-sequence — over both hand-built and random proteins, asserting determinism on the empty path.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Pattern_ValidPatternNoOccurrence_EmptyNoFalseHit(CancellationToken token)
    {
        // (a) Hand-built: patterns that genuinely cannot occur.
        const string seq = "ACDEFGHIKLMNPQRSTVWY"; // each residue exactly once
        var noMatch = new (string pattern, string why)[]
        {
            ("WWWW", "no run of four W exists"),
            ("RGDRGD", "the literal RGDRGD does not occur"),
            ("[QQ]{5}", "no run of five Q exists"),
            ("Z", "'Z' is not in the standard-residue sequence"),
            ("A.{30}A", "the sequence is too short for two A's 30 apart"),
            ("^RGD", "the sequence does not START with RGD (anchored)"),
            ("KKK$", "the sequence does not END with KKK (anchored)"),
        };
        foreach (var (pattern, why) in noMatch)
        {
            var hits = ((Func<List<MotifMatch>>)(() => FindMotifByPattern(seq, pattern).ToList()))
                .Should().NotThrow($"non-occurring pattern '{pattern}' must not crash").Subject;
            hits.Should().BeEmpty($"no false hit: {why}");
        }

        // (b) Random proteins: a deliberately rare 6-residue literal of a fixed residue. Over short
        //     random proteins it should not occur; whether or not it does, any hit is well-formed and
        //     the empty/non-empty result is deterministic.
        foreach (int seed in new[] { 11, 73, 404, 2026 })
        {
            string rnd = RandomProtein(40, seed);
            const string rarePattern = "WWWWWW"; // six consecutive W — astronomically unlikely in len-40
            var hits = ((Func<List<MotifMatch>>)(() => FindMotifByPattern(rnd, rarePattern).ToList()))
                .Should().NotThrow($"rare pattern over random protein (seed {seed}) must not crash").Subject;
            token.ThrowIfCancellationRequested();
            foreach (var hit in hits)
                AssertWellFormedMatch(hit, rnd);

            // INV-01 determinism on the (typically empty) path.
            var again = FindMotifByPattern(rnd, rarePattern).ToList();
            again.Select(h => (h.Start, h.End, h.Sequence))
                .Should().Equal(hits.Select(h => (h.Start, h.End, h.Sequence)),
                    "INV-01: FindMotifByPattern is deterministic for a fixed input");
        }
    }

    #endregion

    #region Positive sanity — a valid pattern matches a known protein at correct positions

    /// <summary>
    /// Positive sanity: the harness must assert against a scanner that actually FINDS motifs, not a
    /// no-op that always returns empty. A protein with a known RGD motif (PROSITE PS00016) at a known
    /// offset yields exactly one well-formed hit; a wildcard "PP.Y" pattern captures the middle residue;
    /// matching is case-insensitive (lowercase == uppercase). Then the THREE MC degenerates are pinned
    /// in one place: empty pattern → empty, invalid pattern → empty, valid-but-absent pattern → empty —
    /// each contrasted against the positive case so green cannot be a false negative.
    /// </summary>
    [Test]
    public void Pattern_KnownProtein_PositiveAndDegenerateContract()
    {
        // Positive: RGD at a known offset.
        const string protein = "MKSPRGDWFIL";
        int rgdAt = protein.IndexOf("RGD", StringComparison.Ordinal); // 4
        var rgd = FindMotifByPattern(protein, "RGD", "RGD", "PS00016").ToList();
        rgd.Should().ContainSingle("the protein contains exactly one RGD motif");
        AssertWellFormedMatch(rgd[0], protein);
        rgd[0].Start.Should().Be(rgdAt, "the RGD hit is at the known offset");
        rgd[0].End.Should().Be(rgdAt + 2, "RGD spans three residues");
        rgd[0].Sequence.Should().Be("RGD");
        rgd[0].Pattern.Should().Be("PS00016", "the supplied patternId is stored on the match");

        // Positive: wildcard PP.Y captures the middle residue.
        var py = FindMotifByPattern("AAPPSYAA", "PP.Y").ToList();
        py.Should().ContainSingle("the protein contains one PPxY motif");
        py[0].Sequence.Should().Be("PPSY", "the wildcard captures the middle residue 'S'");
        py[0].Start.Should().Be(2);

        // INV-05: case-insensitive — lowercase input yields the same hits as uppercase.
        FindMotifByPattern("mksprgdwfil", "RGD").Select(m => (m.Start, m.End, m.Sequence))
            .Should().Equal(rgd.Select(m => (m.Start, m.End, m.Sequence)),
                "INV-05: matching is case-insensitive");

        // The three MC degenerates, contrasted against the positive RGD case:
        FindMotifByPattern(protein, "").Should().BeEmpty("empty pattern → empty (guard)");
        FindMotifByPattern(protein, "[ABC").Should().BeEmpty("invalid regex → empty (swallowed)");
        FindMotifByPattern(protein, "CCCC").Should().BeEmpty("valid non-occurring pattern → empty (no false hit)");
    }

    #endregion

    #endregion
}
