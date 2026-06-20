using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the RnaStructure DOT-BRACKET notation unit — RNA-DOTBRACKET-001.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output, and no *unhandled* runtime
/// exception. Every input must resolve to EITHER a well-defined, theory-correct
/// result, OR a *documented, intentional* validation outcome.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: RNA-DOTBRACKET-001 — dot-bracket / extended WUSS notation
/// Checklist: docs/checklists/03_FUZZING.md, row 149.
/// Algorithm doc: docs/algorithms/RnaStructure/Dot_Bracket_Notation.md
///   (Test Unit ID RNA-DOTBRACKET-001, §3 Contract, §5.3/§5.4, §6 Edge Cases).
/// Methods under test (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/
///   RnaSecondaryStructure.cs):
///   • ParseDotBracket(string) → IEnumerable&lt;(int Position1, int Position2)&gt;
///       decodes the notation into 0-based (5', 3') base-pair index tuples
///       (Dot_Bracket_Notation.md §3.2, §4.1).
///   • ValidateDotBracket(string) → bool
///       true iff the string is well-formed: every closer matches an earlier
///       unmatched opener of the SAME family, and no opener is left unclosed
///       (Dot_Bracket_Notation.md §2.4 INV-04, §3.2, §4.1 step 5).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzz strategy for THIS unit
/// ───────────────────────────────────────────────────────────────────────────
///   • MC = Malformed Content (docs/checklists/03_FUZZING.md §Description).
///     Targets (row 149): "unbalanced brackets, illegal chars, empty, length
///     mismatch". Mapped to the documented contract:
///
///     – UNBALANCED brackets: more openers than closers ("(((") leaves openers
///       unclosed; more closers than openers ("...)" / ")(") meets a closer with
///       an empty/absent same-family stack; mismatched families ("(]") never
///       match up. ValidateDotBracket MUST return false for all of these (never
///       throw). ParseDotBracket is best-effort (Dot_Bracket_Notation.md §5.4
///       deviation #1): a stray closer is silently DROPPED — it must NEVER pop an
///       empty stack, throw, or emit a corrupt pair (e.g. a pair with the wrong
///       endpoint or a closer index for which there was no opener).
///
///     – ILLEGAL chars: any character that is neither a recognized bracket nor a
///       letter (digits, punctuation, whitespace, null byte, unicode) is treated
///       as UNPAIRED and skipped (Dot_Bracket_Notation.md §3.3, §5.3 last bullet).
///       Validation/parse must not crash and must not invent a pair from junk.
///
///     – EMPTY: null and "" are a documented VALID, pair-free structure —
///       ValidateDotBracket(null) == ValidateDotBracket("") == true and
///       ParseDotBracket(null)/("") yield no pairs (Dot_Bracket_Notation.md §5.4
///       deviation #2, §6.1, §3.3). Never an exception.
///
///     – LENGTH mismatch: this API takes ONLY the notation string — there is no
///       paired sequence argument, so positions are 0-based indices INTO THE
///       NOTATION STRING ITSELF (Dot_Bracket_Notation.md §3.1, §2.1). The
///       "length mismatch" failure mode (IndexOutOfRange indexing an external
///       sequence by a structure position) therefore CANNOT occur here: we pin
///       that every returned index is in [0, len) of the notation string for
///       arbitrary, even pathological, inputs — the parser never reaches past its
///       own string and never returns an out-of-range index.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Theory-correct contract asserted on every parse result
/// (Dot_Bracket_Notation.md §2.4 INV-01..INV-04)
/// ───────────────────────────────────────────────────────────────────────────
///   • INV-01 — every returned pair has Position1 &lt; Position2 (opener pushed
///     before its closer is read; the closer index is the larger).
///   • In-range — both endpoints index inside the notation string [0, len).
///   • Partner-uniqueness — each position appears in AT MOST ONE returned pair
///     (one stack entry per opener, popped exactly once).
///   • INV-02 — both endpoints belong to the same family / letter case-pair:
///     notation[i] and notation[j] are a matching open/close (per-family stacks).
///   • INV-03 — for a string ValidateDotBracket accepts, ParseDotBracket returns
///     exactly one pair per opening symbol (balanced ⇒ every opener matched), and
///     the decoded pairs are properly NESTED within each family.
///
/// All inputs are deterministic (locally-seeded new Random(seed)); the long /
/// nesting-depth cases are [CancelAfter]-guarded against an accidental hang. The
/// algorithms are single-pass O(n) (Dot_Bracket_Notation.md §4.3), so the guard
/// is purely a regression net.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RnaDotBracketFuzzTests
{
    #region Helpers

    private const string OpenSymbols = "([{<";
    private const string CloseSymbols = ")]}>";

    private static readonly IReadOnlyDictionary<char, char> CloseToOpen = new Dictionary<char, char>
    {
        [')'] = '(', [']'] = '[', ['}'] = '{', ['>'] = '<',
    };

    /// <summary>Deterministic RNG — seed fixed locally so fuzz inputs are reproducible.</summary>
    private static string RandomNotation(int length, int seed)
    {
        // Mix openers, closers, dots, letters (pseudoknot) and junk so the parser
        // meets every code path: matched, mismatched, stray closer, illegal char.
        const string alphabet = "()[]{}<>...AaBb.,-:_# 9";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Builds a well-formed, properly-nested dot-bracket string of a single family:
    /// <c>depth</c> openers, <c>loop</c> dots, then <c>depth</c> closers — e.g.
    /// "(((...)))". Returns the expected (i,j) pairs, outer→inner.
    /// </summary>
    private static (string Notation, List<(int, int)> Pairs) NestedHairpin(int depth, int loop, char open = '(')
    {
        char close = OpeningToClosingChar(open);
        string notation = new string(open, depth) + new string('.', loop) + new string(close, depth);
        var pairs = new List<(int, int)>();
        int n = notation.Length;
        for (int k = 0; k < depth; k++)
            pairs.Add((k, n - 1 - k)); // (0, n-1), (1, n-2), ...
        return (notation, pairs);
    }

    private static char OpeningToClosingChar(char open) => open switch
    {
        '(' => ')', '[' => ']', '{' => '}', '<' => '>',
        _ => throw new ArgumentException("not an opener", nameof(open)),
    };

    /// <summary>
    /// Asserts the UNIVERSAL invariants every parse result must satisfy for ANY
    /// input string, well-formed or malformed (Dot_Bracket_Notation.md §2.4):
    /// INV-01 (i&lt;j), in-range endpoints into THIS string (no length-mismatch
    /// index escape), partner-uniqueness, and INV-02 (same-family matched
    /// open/close at the two endpoints). This is the "no corrupt pair list" guard.
    /// </summary>
    private static void AssertWellFormedPairs(string notation, IReadOnlyList<(int Position1, int Position2)> pairs)
    {
        int n = notation.Length;
        var seen = new HashSet<int>();

        foreach (var (p1, p2) in pairs)
        {
            // INV-01: opener before closer.
            p1.Should().BeLessThan(p2, "INV-01: a returned pair has Position1 < Position2");

            // In-range: positions index INTO the notation string itself — never past it.
            p1.Should().BeInRange(0, n - 1, "Position1 must index inside the notation string (no out-of-range escape)");
            p2.Should().BeInRange(0, n - 1, "Position2 must index inside the notation string (no out-of-range escape)");

            // Partner-uniqueness: each position is in at most one pair.
            seen.Add(p1).Should().BeTrue($"position {p1} must appear in at most one base pair");
            seen.Add(p2).Should().BeTrue($"position {p2} must appear in at most one base pair");

            // INV-02: the two endpoints are a matching open/close of the SAME family
            // (bracket families) OR an uppercase/lowercase letter pair (pseudoknot).
            char a = notation[p1];
            char b = notation[p2];
            bool bracketFamily = CloseToOpen.TryGetValue(b, out char opener) && opener == a;
            bool letterFamily = char.IsLetter(a) && char.IsLetter(b)
                && char.IsUpper(a) && char.IsLower(b)
                && char.ToUpperInvariant(b) == a;
            (bracketFamily || letterFamily).Should().BeTrue(
                $"INV-02: pair endpoints '{a}'@{p1} and '{b}'@{p2} must be a matching same-family open/close");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  RNA-DOTBRACKET-001 — dot-bracket parsing / validation : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region RNA-DOTBRACKET-001 — dot-bracket notation

    #region MC — Malformed Content: unbalanced brackets

    /// <summary>
    /// MC — unbalanced (more OPENERS than closers). "(((" leaves openers unclosed,
    /// so ValidateDotBracket must return false (INV-04: a non-empty stack after the
    /// scan ⇒ invalid). ParseDotBracket is best-effort: an unmatched OPENER is never
    /// emitted as a pair (it has no closer), so the result is empty here — and never
    /// a crash (Dot_Bracket_Notation.md §4.1 step 5, §5.4 deviation #1, §6.1).
    /// </summary>
    [Test]
    public void Unbalanced_MoreOpenersThanClosers_InvalidAndNoCorruptPairs()
    {
        foreach (string s in new[] { "(((", "(", "((((....)))", "<<<", "[[[..]" })
        {
            RnaSecondaryStructure.ValidateDotBracket(s).Should().BeFalse(
                $"'{s}' leaves an opener unclosed — INV-04 rejects an unbalanced string");

            var pairs = ((Func<List<(int, int)>>)(() => RnaSecondaryStructure.ParseDotBracket(s).ToList()))
                .Should().NotThrow("parse is best-effort on malformed input, never a crash").Subject;
            AssertWellFormedPairs(s, pairs);
            // Only fully-matched openers may surface; a leftover opener never becomes a pair.
            pairs.Count.Should().BeLessThanOrEqualTo(
                Math.Min(s.Count(c => OpenSymbols.Contains(c)), s.Count(c => CloseSymbols.Contains(c))),
                "no pair can exceed the number of matched open/close symbols");
        }
    }

    /// <summary>
    /// MC — unbalanced (more CLOSERS than openers, the empty-stack-pop danger).
    /// "...)" and ")(" and "())" present a closer with no open partner. The KEY
    /// guard: the parse must NOT pop an empty stack (InvalidOperationException) —
    /// the stray closer is silently DROPPED (Dot_Bracket_Notation.md §5.4 #1, §6.1
    /// row "())"), and ValidateDotBracket rejects the string (a closer with an empty
    /// same-family stack ⇒ false). We pin both: no throw, and any pair that DOES
    /// surface is well-formed (never a corrupt (closer,?) tuple).
    /// </summary>
    [Test]
    public void Unbalanced_MoreClosersThanOpeners_DoesNotPopEmptyStack()
    {
        foreach (string s in new[] { ")", "...)", ")(", "())", ")))", ".)]}>", "()())" })
        {
            RnaSecondaryStructure.ValidateDotBracket(s).Should().BeFalse(
                $"'{s}' has a closer with no earlier opener — INV-04 rejects it");

            var pairs = ((Func<List<(int, int)>>)(() => RnaSecondaryStructure.ParseDotBracket(s).ToList()))
                .Should().NotThrow($"a stray closer in '{s}' must be dropped, never an empty-stack pop").Subject;
            AssertWellFormedPairs(s, pairs);
        }

        // Targeted contract example from the doc (§6.1): "())" parses to {(0,1)} only.
        var resolved = RnaSecondaryStructure.ParseDotBracket("())").ToList();
        resolved.Should().Equal(new[] { (0, 1) }, "best-effort: the matched pair survives, the stray ')' is dropped");
    }

    /// <summary>
    /// MC — unbalanced via MISMATCHED families: "(]" is count-balanced but the
    /// closer ']' belongs to a different family than the opener '(' — partners must
    /// match up (Dot_Bracket_Notation.md §2.4 INV-04, §6.1 row "(]"). Validation
    /// must reject it; parse drops the unmatched ']' AND leaves '(' unclosed, so no
    /// cross-family pair is ever fabricated. We sweep all cross-family combinations.
    /// </summary>
    [Test]
    public void Unbalanced_MismatchedFamilies_RejectedAndNeverCrossPaired()
    {
        for (int o = 0; o < OpenSymbols.Length; o++)
        {
            for (int c = 0; c < CloseSymbols.Length; c++)
            {
                if (o == c) continue; // same family is the matched case, skipped here.
                string s = $"{OpenSymbols[o]}{CloseSymbols[c]}";
                RnaSecondaryStructure.ValidateDotBracket(s).Should().BeFalse(
                    $"'{s}' mixes families — partners must match up (INV-04)");

                var pairs = RnaSecondaryStructure.ParseDotBracket(s).ToList();
                AssertWellFormedPairs(s, pairs);
                pairs.Should().BeEmpty(
                    $"'{s}': the closer of a different family pops nothing — no cross-family pair is invented");
            }
        }
    }

    #endregion

    #region MC — Malformed Content: illegal characters

    /// <summary>
    /// MC — illegal chars. Any character that is neither a recognized bracket nor a
    /// letter (digits, punctuation, whitespace, null byte, unicode) is UNPAIRED and
    /// skipped (Dot_Bracket_Notation.md §3.3, §5.3). A pure-junk string is therefore
    /// a VALID, pair-free structure (every char is unpaired ⇒ balanced ⇒ valid), and
    /// parse yields no pairs — never a crash or an invented pair from junk.
    /// </summary>
    [Test]
    public void IllegalChars_AreTreatedAsUnpaired_NoCrashNoInventedPairs()
    {
        // NOTE: these are NON-letter, NON-bracket characters only. The WUSS alphabet
        // is A–Z; char.IsLetter is true for ASCII letters (which the unit treats as
        // PAIRING symbols), so unicode LETTERS are intentionally excluded here and
        // covered separately below (Dot_Bracket_Notation.md §6.2). Emoji such as 🧬
        // are surrogate pairs (char.IsLetter == false) → unpaired junk, kept here.
        foreach (string junk in new[]
        {
            "12345", "!@#$%^&*", "   \t\n ", "\0\0\0", "..,-:_..", "🧬🧬", "999.999"
        })
        {
            RnaSecondaryStructure.ValidateDotBracket(junk).Should().BeTrue(
                $"'{junk}' contains no brackets/letters — every char is unpaired, so it is balanced (valid)");

            var pairs = ((Func<List<(int, int)>>)(() => RnaSecondaryStructure.ParseDotBracket(junk).ToList()))
                .Should().NotThrow($"illegal characters in '{junk}' must be skipped, not crash").Subject;
            pairs.Should().BeEmpty($"no bracket/letter ⇒ no base pair can be decoded from '{junk}'");
        }
    }

    /// <summary>
    /// MC — non-ASCII LETTERS. The unit recognizes openers/closers via
    /// char.IsLetter / case (invariant culture), so non-ASCII letters outside the
    /// WUSS A–Z convention are treated as pairing symbols — a documented limitation
    /// (Dot_Bracket_Notation.md §6.2). The contract that still MUST hold even on
    /// such out-of-convention input: no crash, no empty-stack pop, and any pair that
    /// surfaces is still well-formed (uppercase opener before its lowercase closer,
    /// in range). A lone unmatched non-ASCII letter is simply unbalanced (invalid),
    /// never an exception.
    /// </summary>
    [Test]
    public void IllegalChars_NonAsciiLetters_HandledNotCrashed()
    {
        foreach (string s in new[] { "αβγδ", "Σσ", "ΑΒ..βα", "Ωω" })
        {
            var pairs = ((Func<List<(int, int)>>)(() => RnaSecondaryStructure.ParseDotBracket(s).ToList()))
                .Should().NotThrow($"non-ASCII letters in '{s}' must not crash the parse").Subject;
            AssertWellFormedPairs(s, pairs);

            ((Func<bool>)(() => RnaSecondaryStructure.ValidateDotBracket(s)))
                .Should().NotThrow($"validation must not crash on non-ASCII letters in '{s}'");
        }

        // A matched uppercase→lowercase non-ASCII pair behaves like a letter pair.
        RnaSecondaryStructure.ValidateDotBracket("ΑΒ..βα").Should().BeTrue(
            "Α..α and Β..β nest as valid uppercase/lowercase letter pairs (invariant-culture case)");
    }

    /// <summary>
    /// MC — illegal chars INTERLEAVED with valid brackets. Junk between real
    /// brackets must be ignored without disturbing the bracket matching: "(1.2)"
    /// pairs the '(' at 0 with the ')' at 4 (positions are still indices into the
    /// FULL string, junk included), and "(a#b)" likewise. The junk neither shifts
    /// nor breaks the decoded indices (Dot_Bracket_Notation.md §3.3, INV-01/INV-02).
    /// </summary>
    [Test]
    public void IllegalChars_InterleavedWithBrackets_AreSkippedNotShifted()
    {
        // "(1.2)" — junk at 1,2,3; '(' at 0 pairs with ')' at 4.
        var p1 = RnaSecondaryStructure.ParseDotBracket("(1.2)").ToList();
        AssertWellFormedPairs("(1.2)", p1);
        p1.Should().Equal(new[] { (0, 4) }, "junk between brackets is skipped; indices stay positions in the full string");
        RnaSecondaryStructure.ValidateDotBracket("(1.2)").Should().BeTrue("brackets balance; junk is unpaired");

        // "((# .))" — junk '#' and space ignored; outer/inner brackets still nest.
        var p2 = RnaSecondaryStructure.ParseDotBracket("((# .))").ToList();
        AssertWellFormedPairs("((# .))", p2);
        p2.Should().BeEquivalentTo(new[] { (1, 5), (0, 6) }, "nested () survive embedded junk, inner-first popped");
        RnaSecondaryStructure.ValidateDotBracket("((# .))").Should().BeTrue();
    }

    #endregion

    #region MC — Malformed Content: empty

    /// <summary>
    /// MC — empty. null and "" are a documented VALID, pair-free structure:
    /// ValidateDotBracket(null) == ValidateDotBracket("") == true and
    /// ParseDotBracket(null)/("") yield no pairs (Dot_Bracket_Notation.md §5.4 #2,
    /// §6.1, §3.3). Never an exception, never an index into an empty string.
    /// </summary>
    [Test]
    public void Empty_NullAndEmpty_AreValidAndYieldNoPairs()
    {
        RnaSecondaryStructure.ValidateDotBracket("").Should().BeTrue("empty is an unambiguously balanced structure");
        RnaSecondaryStructure.ValidateDotBracket(null!).Should().BeTrue("null is treated as an empty (valid) structure");

        var empty = ((Func<List<(int, int)>>)(() => RnaSecondaryStructure.ParseDotBracket("").ToList()))
            .Should().NotThrow("empty input is a documented no-op").Subject;
        empty.Should().BeEmpty("the empty structure has no base pairs");

        var nul = ((Func<List<(int, int)>>)(() => RnaSecondaryStructure.ParseDotBracket(null!).ToList()))
            .Should().NotThrow("null input is a documented no-op").Subject;
        nul.Should().BeEmpty("a null structure has no base pairs");

        // An all-unpaired string is likewise valid with no pairs (§6.1 ".....").
        RnaSecondaryStructure.ValidateDotBracket(".....").Should().BeTrue("dots are unpaired ⇒ balanced");
        RnaSecondaryStructure.ParseDotBracket(".....").Should().BeEmpty("no brackets ⇒ no pairs");
    }

    #endregion

    #region MC — Malformed Content: length / index integrity (no out-of-range escape)

    /// <summary>
    /// MC — "length mismatch" mapped to this single-string API. There is NO paired
    /// sequence argument (Dot_Bracket_Notation.md §3.1), so a structure position can
    /// never index PAST an external sequence: positions are 0-based indices into the
    /// notation string itself (§2.1). The failure mode to exclude is an out-of-range
    /// index leaking from the parser. We feed adversarial mixes (stray closers,
    /// mismatched families, junk, deep nesting) and pin that EVERY returned index is
    /// strictly within [0, len) of the input string — the parser never reaches past
    /// its own bounds and never returns a corrupt index.
    /// </summary>
    [Test]
    public void LengthMismatch_EveryReturnedIndexIsWithinTheStringBounds()
    {
        foreach (string s in new[]
        {
            "([)]", ")(][}{><", "(((...)))junk", "(a)(b)cc", "<[({.})]>", "..)((..", "(]([)]"
        })
        {
            var pairs = ((Func<List<(int, int)>>)(() => RnaSecondaryStructure.ParseDotBracket(s).ToList()))
                .Should().NotThrow($"'{s}' must parse without throwing").Subject;

            // The core "no length-mismatch index escape" assertion.
            foreach (var (p1, p2) in pairs)
            {
                p1.Should().BeInRange(0, s.Length - 1, $"index {p1} must stay inside '{s}' (len {s.Length})");
                p2.Should().BeInRange(0, s.Length - 1, $"index {p2} must stay inside '{s}' (len {s.Length})");
            }
            AssertWellFormedPairs(s, pairs);
        }
    }

    #endregion

    #region MC — Malformed Content: random / pathological strings stay well-formed

    /// <summary>
    /// MC — random fuzz strings over a mixed alphabet (brackets, closers, dots,
    /// letters, junk). Across fixed seeds and lengths the parser must NEVER crash or
    /// hang, and EVERY emitted pair must satisfy the universal contract: INV-01
    /// (i&lt;j), in-range, partner-unique, INV-02 (same-family endpoints). Validation
    /// must also never throw. [CancelAfter] guards the O(n) scan against a hang
    /// regression. This is the broad "no corrupt pair list on garbage" net.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void RandomNotation_AlwaysParsesToWellFormedPairs(CancellationToken token)
    {
        foreach (int seed in new[] { 1, 7, 42, 1000, 20260620 })
        {
            foreach (int len in new[] { 0, 1, 8, 50, 500 })
            {
                string s = RandomNotation(len, seed);

                var pairs = ((Func<List<(int, int)>>)(() => RnaSecondaryStructure.ParseDotBracket(s).ToList()))
                    .Should().NotThrow($"random notation must not crash (seed {seed}, len {len})").Subject;
                token.ThrowIfCancellationRequested();

                AssertWellFormedPairs(s, pairs);

                ((Func<bool>)(() => RnaSecondaryStructure.ValidateDotBracket(s)))
                    .Should().NotThrow($"validation must not crash on random notation (seed {seed}, len {len})");
                token.ThrowIfCancellationRequested();
            }
        }
    }

    /// <summary>
    /// MC / OVF — deep nesting. A very deeply-nested single-family hairpin
    /// "((( ... )))" (1000 deep) stresses the stack without a stack overflow or a
    /// hang. It is well-formed, so it must VALIDATE and parse to exactly <c>depth</c>
    /// nested pairs (0,n-1),(1,n-2),... — the well-formed positive at scale.
    /// [CancelAfter] guards the linear scan.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void DeepNesting_ValidatesAndParsesToNestedPairs(CancellationToken token)
    {
        const int depth = 1000;
        var (notation, expected) = NestedHairpin(depth, loop: 3);

        RnaSecondaryStructure.ValidateDotBracket(notation).Should().BeTrue("a deep but balanced hairpin is well-formed");
        token.ThrowIfCancellationRequested();

        var pairs = RnaSecondaryStructure.ParseDotBracket(notation).ToList();
        token.ThrowIfCancellationRequested();

        AssertWellFormedPairs(notation, pairs);
        pairs.Should().HaveCount(depth, "INV-03: one pair per opening symbol when balanced");
        // Inner-first popping: the order is (depth-1, depth+loop), ... but as a set
        // it must equal the nested outer→inner pairing.
        pairs.Should().BeEquivalentTo(expected, "a fully nested hairpin decodes to (k, n-1-k) pairs");
    }

    #endregion

    #region Positive sanity — well-formed notation decodes to the documented pairs

    /// <summary>
    /// Positive sanity: the canonical worked example (Dot_Bracket_Notation.md §7.1).
    /// "((((....))))" decodes to {(0,11),(1,10),(2,9),(3,8)} — a properly nested
    /// 4-bp stem closing a 4-nt loop — and validates. "..." decodes to no pairs.
    /// This proves the fuzz harness asserts against a parser that actually FINDS
    /// structure, not a no-op.
    /// </summary>
    [Test]
    public void Positive_NestedHairpin_DecodesToDocumentedPairs()
    {
        var pairs = RnaSecondaryStructure.ParseDotBracket("((((....))))").ToList();
        AssertWellFormedPairs("((((....))))", pairs);
        pairs.Should().BeEquivalentTo(new[] { (0, 11), (1, 10), (2, 9), (3, 8) },
            "the canonical nested hairpin (Dot_Bracket_Notation.md §7.1)");
        RnaSecondaryStructure.ValidateDotBracket("((((....))))").Should().BeTrue();

        // INV-03 cross-check: pair count == number of opening symbols.
        pairs.Should().HaveCount("((((....))))".Count(c => c == '('),
            "INV-03: balanced ⇒ one pair per opener");

        // All-unpaired decodes to nothing and is valid.
        RnaSecondaryStructure.ParseDotBracket("...").Should().BeEmpty("dots are unpaired");
        RnaSecondaryStructure.ValidateDotBracket("...").Should().BeTrue();
    }

    /// <summary>
    /// Positive sanity: a PSEUDOKNOT with crossing families. "([)]" is valid because
    /// the () and [] families are independent pairing systems matched on their own
    /// stacks, so the helices may cross (Dot_Bracket_Notation.md §6.1 row "([)]",
    /// §7.1). It decodes to {(0,2),(1,3)} — note the pairs CROSS, which a single-stack
    /// parser could not produce. The letter-pair pseudoknot "AA..aa" likewise decodes.
    /// </summary>
    [Test]
    public void Positive_Pseudoknot_CrossingFamiliesDecodeIndependently()
    {
        var bracket = RnaSecondaryStructure.ParseDotBracket("([)]").ToList();
        AssertWellFormedPairs("([)]", bracket);
        bracket.Should().BeEquivalentTo(new[] { (0, 2), (1, 3) },
            "crossing () and [] families are matched independently (a true pseudoknot)");
        RnaSecondaryStructure.ValidateDotBracket("([)]").Should().BeTrue(
            "([)] is well-formed: each family balances on its own stack");

        // Letter-pair pseudoknot: uppercase opens, matching lowercase closes.
        var letters = RnaSecondaryStructure.ParseDotBracket("AA..aa").ToList();
        AssertWellFormedPairs("AA..aa", letters);
        letters.Should().BeEquivalentTo(new[] { (1, 4), (0, 5) },
            "A-letter pair: inner-first popping gives (1,4) then (0,5)");
        RnaSecondaryStructure.ValidateDotBracket("AA..aa").Should().BeTrue();

        // Extended families nested together also validate and decode (doc §2.4 example).
        var nestedFamilies = RnaSecondaryStructure.ParseDotBracket("<[({.})]>").ToList();
        AssertWellFormedPairs("<[({.})]>", nestedFamilies);
        RnaSecondaryStructure.ValidateDotBracket("<[({.})]>").Should().BeTrue(
            "properly nested distinct families form a valid structure");
        nestedFamilies.Should().HaveCount(4, "four nested family pairs");
    }

    #endregion

    #endregion
}
