using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;
// Disambiguate: a top-level Seqeron.Genomics.Analysis.MotifMatch also exists; this unit
// asserts against the ProteinMotifFinder.MotifMatch record returned by FindMotifByProsite.
using MotifMatch = Seqeron.Genomics.Analysis.ProteinMotifFinder.MotifMatch;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the ProteinMotif area — PROSITE pattern matching
/// (PROTMOTIF-PROSITE-001): translating a PROSITE PA-line pattern into an
/// equivalent .NET regular expression and searching a protein with it. The two
/// public entry points under test are the converter
/// <see cref="ProteinMotifFinder.ConvertPrositeToRegex"/> and the end-to-end
/// PROSITE scan <see cref="ProteinMotifFinder.FindMotifByProsite"/> (which converts
/// the pattern, then delegates the converted regex to the overlap-aware
/// <c>FindMotifByPattern</c>). The sibling unit PROTMOTIF-FIND-001 (row 82, the
/// general regex motif scanner) is covered separately in
/// <c>ProteinMotifFuzzTests</c>; this file focuses on the PROSITE GRAMMAR and the
/// converter's "reject, don't silently drop" policy.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, malformed and adversarial inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no NaN, and no *unhandled* runtime exception
/// (IndexOutOfRange / NullReference / ArgumentOutOfRange). Every input must resolve
/// to EITHER a well-defined, theory-correct result OR a *documented, intentional*
/// outcome. For a PROSITE→regex converter the headline hazards are:
///   • a NullReferenceException on a null pattern (guarded by IsNullOrEmpty);
///   • an IndexOutOfRange when an unterminated bracket/paren/brace runs off the end
///     (the parser must fall back gracefully, never index past the string);
///   • a REGEX-INJECTION escape: a PROSITE pattern is CONVERTED, not used verbatim,
///     so a metacharacter smuggled inside the pattern (`*`, `?`, `+`, `\`, `|`, `(`,
///     `$`, `.`, …) must EITHER convert to a well-defined regex token, OR be rejected
///     with a documented FormatException, OR be swallowed by the delegated matcher —
///     but must NEVER crash with an unexpected exception type, hang, or corrupt the
///     reported match coordinates;
///   • a catastrophic-backtracking / runaway scan on a long or pathological pattern
///     (kept [CancelAfter]-guarded so a regression FAILS rather than wedges the suite).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROTMOTIF-PROSITE-001 — PROSITE pattern matching (convert + scan)
/// Checklist: docs/checklists/03_FUZZING.md, row 83.
/// Fuzz strategies exercised for THIS unit:
///   • MC = Malformed Content — invalid / truncated PROSITE syntax:
///       – empty / null pattern → ConvertPrositeToRegex returns "" and
///         FindMotifByProsite yields NO matches (explicit IsNullOrEmpty guard;
///         PROSITE_Pattern_Matching.md §3.3, §6.1).
///       – unterminated `[`, `{`, `(` (no closing `]`/`}`/`)`): the parser's
///         IndexOf-based scan finds no terminator and falls back to emitting the
///         lone metacharacter / advancing, never running off the end (no
///         IndexOutOfRange); the converted regex may then be invalid and the
///         delegated matcher swallows it (PROSITE_Pattern_Matching.md §3.3).
///       – trailing content after the `.` terminator is ignored
///         (INV-02 PROSITE_Pattern_Matching.md §2.4; §6.1 "Trailing period").
///   • INJ = Injection — regex metacharacters smuggled through the converter:
///       – the ScanProsite extended Kleene star `*` (and stray `?`/`+`) is NOT part
///         of the PA-line grammar and MUST raise FormatException — "reject, don't
///         silently drop" (Pattern_Matching_Methods.md INV-06, §5.2; the converter's
///         own source comment). This is the headline injection contract: a pattern
///         author cannot smuggle an unbounded quantifier past the grammar.
///       – metacharacters that ARE legal PROSITE atoms — `<` `>` `[` `]` `{` `}`
///         `(` `)` `-` `.` `x` — convert to their DEFINED regex tokens (`^`, `$`,
///         classes, `[^…]`, quantifiers, wildcard) per the §4.2 grammar table.
///       – metacharacters with NO PROSITE meaning that ALSO are not in the reject
///         list (`|`, `\`, `$`, `^`, raw `.` mid-class) are exercised end-to-end to
///         confirm they never crash with an unexpected exception type or corrupt
///         coordinates — whatever regex they produce, every emitted match is a
///         well-formed in-bounds span.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes MC, INJ);
///   targets: "Invalid PROSITE pattern syntax, regex injection, empty pattern".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The PROSITE-matching contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// ConvertPrositeToRegex maps each PA-line atom to a regex token deterministically
/// (PROSITE_Pattern_Matching.md §4.2; Pattern_Matching_Methods.md §2.2 INV-04):
///   A (letter) → A         x → .              [ABC] → [ABC]      {ABC} → [^ABC]
///   - (sep)    → dropped    x(n) → .{n}        x(n,m) → .{n,m}    A(n) → A{n}
///   <          → ^          >    → $           [G>]  → (?:G|$)    .    → ends parse
/// FindMotifByProsite converts the PA line, then delegates the converted regex to
/// FindMotifByPattern, which uppercases the protein, wraps the regex in a zero-width
/// lookahead "(?=(…))" so OVERLAPPING occurrences are discovered, and emits one
/// MotifMatch per capture: { Start (incl. 0-based), End (incl. 0-based), Sequence
/// (uppercased matched substring), MotifName, Pattern (the ORIGINAL PROSITE string),
/// Score (information content, bits), EValue }.
///   — docs/algorithms/ProteinMotif/PROSITE_Pattern_Matching.md §2.2, §3, §4;
///     docs/algorithms/ProteinMotif/Pattern_Matching_Methods.md §2.2, §3, §4.
///
/// Methods under test (src/.../Seqeron.Genomics.Analysis/ProteinMotifFinder.cs):
///   string ConvertPrositeToRegex(string prositePattern)
///       — PROSITE PA-line → .NET regex; "" for null/empty; FormatException on `*`/`?`/`+`.
///   IEnumerable&lt;MotifMatch&gt; FindMotifByProsite(string proteinSequence,
///       string prositePattern, string motifName = "Custom")
///       — convert, then delegate to the overlap-aware FindMotifByPattern.
///
/// Theory-correct invariants asserted:
///   • INV-conv — ConvertPrositeToRegex maps each atom to its documented regex token
///     (PROSITE_Pattern_Matching.md §2.4 INV-01, §4.2; Pattern_Matching_Methods.md INV-04).
///   • INV-term — a `.` terminates parsing; trailing content is ignored
///     (PROSITE_Pattern_Matching.md §2.4 INV-02).
///   • INV-reject — `*`/`?`/`+` raise FormatException, never a silent drop
///     (Pattern_Matching_Methods.md INV-06).
///   • INV-span — every emitted match is a CONTIGUOUS in-bounds span whose
///     coordinates reproduce its substring: S[Start..End] == Sequence (uppercased),
///     0 ≤ Start ≤ End ≤ n−1 (PROSITE_Pattern_Matching.md §5.2; INV-02 of
///     Pattern_Matching_Methods.md). Pattern carries the ORIGINAL PROSITE string.
///   • INV-finite — Score and EValue are finite (never NaN / ±∞).
///   • INV-case — matching is case-insensitive (input uppercased + IgnoreCase).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Complexity / hang-safety
/// ───────────────────────────────────────────────────────────────────────────
/// Conversion is a single O(p) left-to-right pass and matching is one lookahead
/// regex walk (PROSITE_Pattern_Matching.md §4.3). The adversarial / long-pattern
/// and homopolymer targets maximise both pattern cost and overlapping-hit count;
/// they are kept modest and [CancelAfter]-guarded so a regression that turned the
/// converter into a hang, or produced a catastrophically-backtracking regex, would
/// FAIL rather than wedge the suite.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ProteinPrositeFuzzTests
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
    /// satisfy against the original (case-insensitive) sequence and the original PROSITE pattern
    /// (PROSITE_Pattern_Matching.md §3.2, §5.2; Pattern_Matching_Methods.md INV-02): the match is a
    /// CONTIGUOUS in-bounds subsequence (INV-span) whose claimed coordinates actually reproduce the
    /// reported (uppercased) substring, whose Score / EValue are finite (INV-finite), and whose
    /// <c>Pattern</c> is the ORIGINAL PROSITE string (not the converted regex). This is the headline
    /// "no coordinate bug, no run-off-the-end, no NaN" property.
    /// </summary>
    private static void AssertWellFormedMatch(MotifMatch match, string originalSequence, string prositePattern)
    {
        string upper = originalSequence.ToUpperInvariant();
        int n = upper.Length;

        // INV-span — in-bounds, non-empty, contiguous span.
        match.Start.Should().BeInRange(0, n - 1, "a match Start is a valid 0-based residue index");
        match.End.Should().BeInRange(match.Start, n - 1, "a match End is in-bounds and not before its Start");
        (match.End - match.Start + 1).Should().Be(match.Sequence.Length,
            "End − Start + 1 equals the matched substring length");
        match.Sequence.Should().Be(match.Sequence.ToUpperInvariant(), "the matched substring is uppercased");
        upper.Substring(match.Start, match.Sequence.Length).Should().Be(match.Sequence,
            "INV-span: the reported substring is exactly S[Start..End] of the uppercased input (no coordinate bug)");

        // FindMotifByProsite stores the ORIGINAL PROSITE pattern in Pattern (§5.2).
        match.Pattern.Should().Be(prositePattern,
            "MotifMatch.Pattern carries the original PROSITE pattern string, not the converted regex");

        // INV-finite — score and E-value are finite.
        double.IsNaN(match.Score).Should().BeFalse("a motif Score must never be NaN");
        double.IsInfinity(match.Score).Should().BeFalse("a motif Score must never be infinite");
        double.IsNaN(match.EValue).Should().BeFalse("a motif EValue must never be NaN");
        double.IsInfinity(match.EValue).Should().BeFalse("a motif EValue must never be infinite");
    }

    /// <summary>Asserts the produced regex is itself a compilable .NET regex (a well-defined conversion).</summary>
    private static void AssertCompilableRegex(string regex)
    {
        var compile = () => _ = new Regex(regex);
        compile.Should().NotThrow($"a converted PROSITE pattern must yield a compilable regex (\"{regex}\")");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PROTMOTIF-PROSITE-001 — PROSITE pattern matching : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PROTMOTIF-PROSITE-001 — PROSITE pattern matching

    #region Positive sanity — known PROSITE patterns convert + match correctly

    /// <summary>
    /// Positive sanity (NOT a rubber stamp): the converter must reproduce the DOCUMENTED §4.2 grammar
    /// table verbatim. We pin every atom against the exact regex the docs promise — literal residues,
    /// the `x` wildcard, `x(n)` / `x(n,m)` repetition, `A(n)` fixed counts on a letter, character
    /// classes `[…]`, exclusion classes `{…}` → `[^…]`, `-` separators dropped, anchors `<`→`^` /
    /// `>`→`$`, the rare `[G>]`→`(?:G|$)` C-terminal bracket, and the `.` terminator — including the
    /// canonical PROSITE entries (N-glycosylation PS00001 → <c>N[^P][ST][^P]</c>, CK2 PS00006 →
    /// <c>[ST].{2}[DE]</c>, zinc-finger PS00028 with `x(2,4)`/`x(3,5)`).
    /// (PROSITE_Pattern_Matching.md §4.2, §2.4 INV-01/02/03; Pattern_Matching_Methods.md §2.2 INV-04.)
    /// </summary>
    [Test]
    public void ConvertPrositeToRegex_DocumentedGrammar_MapsEachAtomVerbatim()
    {
        var cases = new (string Prosite, string ExpectedRegex)[]
        {
            // Literal residues and separators.
            ("R-G-D", "RGD"),                                   // PS00016 cell-attachment
            ("P-P-x-Y", "PP.Y"),                                // WW PY motif, single wildcard
            // Canonical PROSITE entries from CommonMotifs.
            ("N-{P}-[ST]-{P}", "N[^P][ST][^P]"),                // PS00001 N-glycosylation
            ("[ST]-x-[RK]", "[ST].[RK]"),                       // PS00005 PKC
            ("[ST]-x(2)-[DE]", "[ST].{2}[DE]"),                 // PS00006 CK2 — x(n)
            ("[RK](2)-x-[ST]", "[RK]{2}.[ST]"),                 // PS00004 — (n) on a class
            ("[AG]-x(4)-G-K-[ST]", "[AG].{4}GK[ST]"),           // PS00017 P-loop
            ("C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H",      // PS00028 zinc finger — x(n,m)
                "C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H"),
            ("L-x(6)-L-x(6)-L-x(6)-L", "L.{6}L.{6}L.{6}L"),     // PS00029 leucine zipper
            // Anchors.
            ("<M", "^M"),                                       // N-terminus
            ("R-G-D>", "RGD$"),                                 // C-terminus
            ("<M-x-[AG]>", "^M.[AG]$"),                         // both anchors
            // Fixed count on a literal residue: A(n) → A{n} (valid per §2.2).
            ("A(3)", "A{3}"),
            // Rare [G>] C-terminal bracket → (?:G|$)  (PS00267 / PS00539; INV-03).
            ("F-[GSTV]-P-R-L-[G>]", "F[GSTV]PRL(?:G|$)"),
            // Trailing '.' terminator: parsing stops, suffix ignored (INV-term, §6.1).
            ("R-G-D.this-is-ignored", "RGD"),
        };

        foreach (var (prosite, expected) in cases)
        {
            string regex = ConvertPrositeToRegex(prosite);
            regex.Should().Be(expected,
                $"PROSITE \"{prosite}\" must convert to the documented regex per §4.2 grammar table");
            AssertCompilableRegex(regex);
        }
    }

    /// <summary>
    /// Positive sanity end-to-end: a known PROSITE pattern run through FindMotifByProsite must locate
    /// its occurrence at the KNOWN offset with the correct span and overlap behavior, not a no-op.
    /// We use the N-glycosylation site N-{P}-[ST]-{P} (PS00001 → <c>N[^P][ST][^P]</c>) at a hand-built
    /// offset, verify the exclusion class genuinely EXCLUDES proline (an N-P-S-A site must NOT match),
    /// confirm case-insensitivity, and confirm overlap-aware enumeration via a wildcard pattern in a
    /// homopolymer. (PROSITE_Pattern_Matching.md §5.2; Pattern_Matching_Methods.md §7.1, INV-02/05.)
    /// </summary>
    [Test]
    public void FindMotifByProsite_KnownPattern_FoundAtCorrectPositionWithExclusions()
    {
        // N-glycosylation N-{P}-[ST]-{P}: "AA NQS A AA" has a valid site N-Q-S-A at offset 2.
        const string glyco = "AANQSAAA";
        var hits = FindMotifByProsite(glyco, "N-{P}-[ST]-{P}", "ASN_GLYCOSYLATION").ToList();
        hits.Should().ContainSingle("the protein contains exactly one N-glycosylation site");
        AssertWellFormedMatch(hits[0], glyco, "N-{P}-[ST]-{P}");
        hits[0].Start.Should().Be(2, "the N-glycosylation site begins at the N at offset 2");
        hits[0].End.Should().Be(5, "the four-residue site spans [2,5]");
        hits[0].Sequence.Should().Be("NQSA", "the matched substring is exactly the N-Q-S-A site");

        // The {P} exclusion genuinely excludes proline: N-P-S-A must NOT match (position 2 is {P}).
        FindMotifByProsite("AANPSAAA", "N-{P}-[ST]-{P}").Should().BeEmpty(
            "the {P} exclusion forbids proline at position 2, so N-P-S-A is not an N-glycosylation site");
        // …and the third position requires S or T: N-Q-A-A (no S/T) must NOT match.
        FindMotifByProsite("AANQAAAA", "N-{P}-[ST]-{P}").Should().BeEmpty(
            "[ST] requires Ser or Thr at position 3, so N-Q-A-A does not match");

        // Case-insensitivity: a lowercase protein yields the same hit (INV-case).
        var lower = FindMotifByProsite("aanqsaaa", "N-{P}-[ST]-{P}").ToList();
        lower.Select(m => (m.Start, m.End, m.Sequence))
            .Should().Equal(hits.Select(m => (m.Start, m.End, m.Sequence)),
                "PROSITE matching is case-insensitive (input uppercased + IgnoreCase)");

        // Overlap-aware enumeration: x(2) over a homopolymer reports every overlapping start.
        const string poly = "AAAAA";
        var overlap = FindMotifByProsite(poly, "x-x", "TwoAny").ToList(); // → ".." width-2 wildcard
        overlap.Should().HaveCount(poly.Length - 1,
            "the converted '..' matches every overlapping 2-residue window (lookahead-based scan)");
        foreach (var h in overlap)
        {
            AssertWellFormedMatch(h, poly, "x-x");
            h.Sequence.Should().Be("AA", "each overlapping width-2 hit is 'AA' in the homopolymer");
        }
    }

    #endregion

    #region MC — Empty / null pattern: "" regex, no matches, no throw

    /// <summary>
    /// Target "empty pattern": ConvertPrositeToRegex returns the EMPTY STRING for a null or empty
    /// pattern, and FindMotifByProsite then yields NO matches — by the explicit IsNullOrEmpty guard,
    /// never a NullReferenceException (PROSITE_Pattern_Matching.md §3.3, §6.1 "Empty or null pattern";
    /// Pattern_Matching_Methods.md INV-01). A whitespace-only or all-separator pattern likewise
    /// collapses to an empty/no-op regex and yields nothing. A null/empty PROTEIN with any pattern
    /// also yields nothing.
    /// </summary>
    [Test]
    public void Prosite_EmptyOrNullPattern_EmptyRegexNoMatchesNoThrow()
    {
        // ConvertPrositeToRegex: null / empty → "".
        foreach (string? pattern in new[] { "", null })
        {
            var convert = () => ConvertPrositeToRegex(pattern!);
            convert.Should().NotThrow($"empty/null pattern ('{pattern ?? "null"}') must not crash the converter")
                .Subject.Should().BeEmpty("an empty/null PROSITE pattern converts to the empty regex string");

            // FindMotifByProsite over a real protein with an empty/null pattern → no matches, no throw.
            var find = () => FindMotifByProsite("MKRGDSPEKWFIL", pattern!).ToList();
            find.Should().NotThrow($"empty/null pattern ('{pattern ?? "null"}') must not crash FindMotifByProsite")
                .Subject.Should().BeEmpty("an empty/null PROSITE pattern yields no matches");
        }

        // All-separator / whitespace-only patterns collapse to a no-op and yield nothing.
        foreach (string pattern in new[] { "-", "---", "- - -" })
        {
            var find = () => FindMotifByProsite("MKRGDSPEKWFIL", pattern).ToList();
            // Note: a space is not a PROSITE atom and would be rejected; restrict to separator-only here.
            if (pattern.Contains(' '))
                continue;
            find.Should().NotThrow($"separator-only pattern ('{pattern}') must not crash")
                .Subject.Should().BeEmpty("a pattern of only '-' separators has no atoms and matches nothing");
            ConvertPrositeToRegex(pattern).Should().BeEmpty("'-' separators are dropped, leaving the empty regex");
        }

        // Null / empty PROTEIN with a real pattern → no matches, no throw.
        foreach (string? seq in new[] { "", null })
        {
            var find = () => FindMotifByProsite(seq!, "R-G-D").ToList();
            find.Should().NotThrow($"empty/null protein ('{seq ?? "null"}') must not crash FindMotifByProsite")
                .Subject.Should().BeEmpty("an empty/null protein yields no PROSITE matches");
        }
    }

    #endregion

    #region MC — Invalid / truncated PROSITE syntax: graceful, no IndexOutOfRange

    /// <summary>
    /// Target "invalid PROSITE pattern syntax": malformed or TRUNCATED patterns — an unterminated `[`,
    /// `{`, or `(` with no closing delimiter — must resolve to a DISCIPLINED outcome, NEVER an
    /// IndexOutOfRangeException / NullReferenceException from the IndexOf-based parser running off the
    /// end. Two disciplined outcomes are documented and observed:
    ///   • the parser falls back to emitting the lone metacharacter / advancing one position when no
    ///     terminator is found (ProteinMotifFinder.cs branches for `[`, `{`, `(`), so conversion
    ///     returns SOME string (possibly an invalid .NET regex, which the delegated matcher then
    ///     SWALLOWS, yielding no hits — Pattern_Matching_Methods.md §6.1); OR
    ///   • the truncated tail exposes a stray non-grammar character (e.g. a dangling digit after a
    ///     half-open `x(2`), which is an UNSUPPORTED construct and is REJECTED with the documented
    ///     <see cref="FormatException"/> (reject-don't-drop, INV-06) — never a silent drop.
    /// This test pins that the ONLY exception type permitted to escape is FormatException, and that
    /// any hits produced are well-formed (PROSITE_Pattern_Matching.md §3.3; Pattern_Matching_Methods.md §6.1).
    /// </summary>
    [Test]
    public void ConvertPrositeToRegex_TruncatedDelimiters_GracefulNoCrash()
    {
        const string protein = "MKRGDSTPEKWFILACDE";
        foreach (string pattern in new[]
                 {
                     "[ABC",            // unterminated class
                     "{ABC",            // unterminated exclusion
                     "x(2",             // unterminated x repetition (dangling digit → FormatException)
                     "A(3",             // unterminated letter repetition (dangling digit → FormatException)
                     "[ST]-x(",         // unterminated trailing repetition
                     "[",               // lone open bracket
                     "{",               // lone open brace
                     "(",               // lone open paren
                     "[]",              // empty class
                     "{}",              // empty exclusion
                     "()",              // empty paren / zero-length repetition
                     "[<G",             // truncated N-terminal bracket
                     "[G>",             // truncated C-terminal bracket
                     "R-[ST]-x(2,",     // unterminated range (dangling digits/comma → FormatException)
                 })
        {
            // Conversion: disciplined outcome only — a string, OR the documented FormatException.
            // It must NEVER be an IndexOutOfRange / NullReference (running off the end of the string).
            string? regex = null;
            try
            {
                regex = ConvertPrositeToRegex(pattern);
            }
            catch (FormatException)
            {
                // documented reject-don't-drop outcome for a stray non-grammar char in the truncated tail
            }
            // The crucial assertion: no OTHER exception type escaped (the try above only catches
            // FormatException, so an IndexOutOfRange/NullReference would have failed the test already).

            // End-to-end: whatever regex results, FindMotifByProsite must not throw an UNHANDLED
            // exception — a valid regex yields well-formed hits, an invalid one yields nothing, and a
            // rejected pattern surfaces the SAME FormatException.
            List<MotifMatch> hits = new();
            var find = () =>
            {
                try { hits = FindMotifByProsite(protein, pattern, "Truncated").ToList(); }
                catch (FormatException) { /* documented rejection — acceptable */ }
            };
            find.Should().NotThrow(
                $"truncated PROSITE \"{pattern}\" must not crash FindMotifByProsite with an unhandled exception type");
            foreach (var h in hits)
                AssertWellFormedMatch(h, protein, pattern);
        }
    }

    /// <summary>
    /// Target "invalid PROSITE pattern syntax" (terminator semantics): the `.` ENDS parsing
    /// (PROSITE_Pattern_Matching.md §2.4 INV-02, §6.1 "Trailing period"; Pattern_Matching_Methods.md
    /// §2.2). Everything after the first `.` — including otherwise-illegal junk or even a `*` — is
    /// IGNORED, so a pattern whose only "invalid" content sits past the terminator still converts
    /// cleanly. A leading `.` produces the empty regex (immediate termination).
    /// </summary>
    [Test]
    public void ConvertPrositeToRegex_PeriodTerminator_IgnoresTrailingContent()
    {
        ConvertPrositeToRegex("R-G-D.").Should().Be("RGD", "a trailing '.' terminates after RGD");
        ConvertPrositeToRegex("R-G-D.x-x-x").Should().Be("RGD", "content after '.' is ignored");
        ConvertPrositeToRegex("R-G-D.*?+|\\garbage").Should().Be("RGD",
            "even otherwise-illegal junk after the '.' terminator is ignored (parsing already stopped)");
        ConvertPrositeToRegex(".R-G-D").Should().BeEmpty("a leading '.' terminates immediately → empty regex");

        // End-to-end the terminated pattern still finds RGD.
        var hits = FindMotifByProsite("AARGDAA", "R-G-D.ignored").ToList();
        hits.Should().ContainSingle("the pattern up to '.' is R-G-D, which occurs once");
        hits[0].Sequence.Should().Be("RGD");
        hits[0].Pattern.Should().Be("R-G-D.ignored", "MotifMatch.Pattern keeps the ORIGINAL (untruncated) pattern");
    }

    #endregion

    #region INJ — Kleene star / unsupported quantifiers: reject, don't silently drop

    /// <summary>
    /// Target "regex injection" (headline): the ScanProsite extended Kleene star `*` (e.g. <c>&lt;{C}*&gt;</c>),
    /// and the stray quantifiers `?` and `+`, are NOT part of the PA-line grammar and MUST be REJECTED
    /// with a <see cref="FormatException"/> — the converter's explicit "reject, don't silently drop"
    /// policy (Pattern_Matching_Methods.md INV-06, §5.2; PROSITE_Pattern_Matching.md §2.2; the
    /// converter's own source comment). This is the core anti-injection guarantee: a pattern author
    /// cannot smuggle an UNBOUNDED quantifier past the grammar to be silently dropped (mis-parsing the
    /// pattern) or to inflate the regex into a catastrophic-backtracking bomb. The exception MUST be
    /// the documented FormatException — never an IndexOutOfRange / NullReference / regex parse error —
    /// and FindMotifByProsite must propagate the SAME FormatException (it converts before matching).
    /// </summary>
    [Test]
    public void ConvertPrositeToRegex_UnsupportedQuantifiers_ThrowFormatException()
    {
        foreach (string pattern in new[]
                 {
                     "<{C}*>",          // canonical ScanProsite Kleene-star query (docs example)
                     "C*",              // star on a literal
                     "[ST]*",           // star on a class
                     "x*",              // star on the wildcard
                     "R-G-D*",          // star at the tail
                     "C?",              // optional quantifier (unsupported)
                     "C+",              // one-or-more quantifier (unsupported)
                     "R-G-D-*-Y",       // star mid-pattern
                 })
        {
            var convert = () => ConvertPrositeToRegex(pattern);
            convert.Should().ThrowExactly<FormatException>(
                    $"the unsupported quantifier in PROSITE \"{pattern}\" must be rejected, not silently dropped")
                .WithMessage("*Unsupported PROSITE construct*",
                    "the FormatException names the offending construct (reject-don't-drop policy)");

            // FindMotifByProsite converts first, so it propagates the SAME FormatException.
            var find = () => FindMotifByProsite("MKRGDSTPEKWFIL", pattern).ToList();
            find.Should().ThrowExactly<FormatException>(
                "FindMotifByProsite converts before matching, so the rejection surfaces end-to-end");
        }

        // Contrast: a star AFTER the '.' terminator is NOT reached, so it does NOT throw.
        var safe = () => ConvertPrositeToRegex("R-G-D.*");
        safe.Should().NotThrow("a '*' after the '.' terminator is never parsed, so no rejection occurs")
            .Subject.Should().Be("RGD");
    }

    #endregion

    #region INJ — Legal PROSITE metacharacters: convert to defined regex tokens

    /// <summary>
    /// Target "regex injection" (constrained metacharacters): the characters that ARE legal PROSITE
    /// atoms — anchors `<` `>`, brackets `[` `]`, braces `{` `}`, parens `(` `)`, the `-` separator,
    /// the `x` wildcard and the `.` terminator — convert to their DEFINED regex tokens, NOT to a
    /// verbatim injection. Because the pattern is CONVERTED rather than used as-is, the only `.`/`^`/`$`
    /// in the output come from the grammar (`x`→`.`, `<`→`^`, `>`→`$`), so an author cannot inject a
    /// raw anchor or wildcard except through its grammatical meaning. We verify each anchor matches at
    /// the correct boundary end-to-end and that `[G>]` matches EITHER a literal G OR the C-terminus.
    /// (PROSITE_Pattern_Matching.md §2.4 INV-03, §4.2, §6.1; Pattern_Matching_Methods.md INV-04.)
    /// </summary>
    [Test]
    public void Prosite_LegalMetacharacters_ConvertToDefinedTokens()
    {
        // N-terminus anchor: <M matches only at position 0.
        ConvertPrositeToRegex("<M").Should().Be("^M");
        FindMotifByProsite("MKKKK", "<M").Should().ContainSingle("'<M' anchors M to the N-terminus")
            .Which.Start.Should().Be(0);
        FindMotifByProsite("KMKKK", "<M").Should().BeEmpty("'<M' cannot match an internal M");

        // C-terminus anchor: D> matches only at the last position.
        ConvertPrositeToRegex("R-G-D>").Should().Be("RGD$");
        FindMotifByProsite("AARGD", "R-G-D>").Should().ContainSingle("'R-G-D>' anchors RGD to the C-terminus")
            .Which.End.Should().Be(4);
        FindMotifByProsite("RGDAA", "R-G-D>").Should().BeEmpty("'R-G-D>' cannot match RGD that is not C-terminal");

        // [G>]: literal G OR end-of-sequence (INV-03).
        ConvertPrositeToRegex("[G>]").Should().Be("(?:G|$)");
        // …matches a literal G in the middle.
        FindMotifByProsite("AAGAA", "A-[G>]").Should().Contain(m => m.Sequence == "AG",
            "'[G>]' matches a literal G inside the sequence");
        // …and matches end-of-sequence after a non-G residue (the '$' branch), reported as a
        // zero-extra-width span ending exactly at the terminus.
        var endHit = FindMotifByProsite("AAAAK", "K-[G>]").ToList();
        endHit.Should().NotBeEmpty("'[G>]' matches end-of-sequence after the K at the C-terminus");
        foreach (var h in endHit)
            AssertWellFormedMatch(h, "AAAAK", "K-[G>]");
        // …but the '$' branch cannot match mid-sequence with no G (§6.1).
        FindMotifByProsite("AAKAA", "K-[G>]").Should().BeEmpty(
            "'[G>]' after a mid-sequence K with no following G cannot match (the '$' branch needs the terminus)");
    }

    #endregion

    #region INJ — Foreign metacharacters with no PROSITE meaning: no unexpected crash

    /// <summary>
    /// Target "regex injection" (foreign metacharacters): characters that are NEITHER legal PROSITE
    /// atoms NOR in the explicit reject list (`*`/`?`/`+`) — namely `|`, `\`, `$`, `^`, `@`, `#`, a
    /// digit, a comma outside a range — are fed through the converter and end-to-end. The fuzz bar:
    /// whatever happens, it must be a DISCIPLINED outcome — either a clean conversion, OR the
    /// documented FormatException (these are not letters and fall through the converter's letter /
    /// grammar branches to the reject clause), OR (if the produced regex is invalid) a swallowed
    /// compile failure yielding no matches — but NEVER an IndexOutOfRange / NullReference / hang, and
    /// any emitted match is a well-formed in-bounds span. This proves the converter constrains
    /// injection: a smuggled metacharacter cannot crash the unit or corrupt coordinates.
    /// (Pattern_Matching_Methods.md INV-06, §6.1; PROSITE_Pattern_Matching.md §3.3.)
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Prosite_ForeignMetacharacters_DisciplinedOutcomeNoCrash(CancellationToken token)
    {
        const string protein = "MKRGDSTPEKWFILACDEYHQN";

        foreach (string pattern in new[]
                 {
                     "R-G|D",           // alternation metacharacter
                     "R-G\\D",          // backslash escape
                     "R-G$D",           // dollar mid-pattern (raw)
                     "R-G^D",           // caret mid-pattern (raw)
                     "R-G@D",           // at-sign (no meaning)
                     "R-G#D",           // hash (no meaning)
                     "R-G2D",           // a DIGIT as a position element
                     "R-[ST],x",        // stray comma outside a range
                     " RGD",       // embedded null byte
                     "RGD‮",       // unicode right-to-left override
                     "R-G-Д",           // non-ASCII (Cyrillic) "letter"
                 })
        {
            // Conversion: either succeeds (returns a string) or throws ONLY the documented FormatException.
            string? regex = null;
            FormatException? rejected = null;
            try
            {
                regex = ConvertPrositeToRegex(pattern);
            }
            catch (FormatException fx)
            {
                rejected = fx; // documented reject-don't-drop outcome
            }

            (regex is not null || rejected is not null).Should().BeTrue(
                $"foreign-metacharacter PROSITE \"{pattern}\" must resolve to a conversion OR a FormatException");

            // End-to-end must mirror the converter: no unhandled exception type, only FormatException
            // is permitted to escape; any emitted hit is a well-formed in-bounds span.
            List<MotifMatch> hits = new();
            var find = () =>
            {
                try { hits = FindMotifByProsite(protein, pattern, "Foreign").ToList(); }
                catch (FormatException) { /* documented rejection — acceptable */ }
            };
            find.Should().NotThrow(
                $"foreign-metacharacter PROSITE \"{pattern}\" must not crash with an unhandled exception type");
            token.ThrowIfCancellationRequested();

            foreach (var h in hits)
                AssertWellFormedMatch(h, protein, pattern);
        }
    }

    #endregion

    #region INJ/OVF — Long & adversarial patterns: no hang, no backtracking blow-up

    /// <summary>
    /// Target "regex injection" (resource exhaustion): long and pathological PROSITE patterns must
    /// convert in a single O(p) pass and produce a regex that scans without catastrophic backtracking
    /// (PROSITE_Pattern_Matching.md §4.3). We build a long alternating literal/wildcard chain, a long
    /// repetition chain (`x(2)` repeated), and a deeply nested-looking class chain, run them against a
    /// long homopolymer and a long random protein under [CancelAfter], and require termination with
    /// only well-formed spans. A regression that turned the converter into a hang or produced a
    /// backtracking bomb would FAIL here rather than wedge the suite.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Prosite_LongAndAdversarialPatterns_TerminateWellFormed(CancellationToken token)
    {
        // (a) Long literal/wildcard chain: A-x-A-x-… (200 atoms).
        string chain = string.Join("-", Enumerable.Range(0, 100).Select(i => i % 2 == 0 ? "A" : "x"));
        string chainRegex = ConvertPrositeToRegex(chain);
        AssertCompilableRegex(chainRegex);
        token.ThrowIfCancellationRequested();

        // (b) Long repetition chain: x(2)-x(2)-… (50 reps → .{2} fifty times).
        string reps = string.Join("-", Enumerable.Repeat("x(2)", 50));
        string repsRegex = ConvertPrositeToRegex(reps);
        repsRegex.Should().Be(string.Concat(Enumerable.Repeat(".{2}", 50)),
            "each x(2) maps to .{2}; the chain is a straightforward concatenation (no blow-up)");
        AssertCompilableRegex(repsRegex);
        token.ThrowIfCancellationRequested();

        // (c) Long class chain: [AG]-[ST]-… (60 classes).
        string classes = string.Join("-", Enumerable.Repeat("[AG]", 60));
        AssertCompilableRegex(ConvertPrositeToRegex(classes));
        token.ThrowIfCancellationRequested();

        // Run each against a long homopolymer and a long random protein — must terminate, no hang.
        string[] proteins = { new string('A', 500), RandomProtein(500, 83) };
        foreach (string protein in proteins)
        {
            foreach (string pattern in new[] { chain, reps, classes })
            {
                List<MotifMatch> hits = null!;
                var find = () => hits = FindMotifByProsite(protein, pattern, "Long").ToList();
                find.Should().NotThrow($"long PROSITE pattern must not crash over a length-{protein.Length} protein");
                token.ThrowIfCancellationRequested();
                foreach (var h in hits)
                    AssertWellFormedMatch(h, protein, pattern);
            }
        }
    }

    #endregion

    #region INJ — Random PROSITE patterns: convert is total, scan is disciplined

    /// <summary>
    /// Property fuzz: across fixed seeds, RANDOMLY ASSEMBLED PROSITE patterns (built only from LEGAL
    /// PA-line atoms — letters, `x`, `x(n)`, `[…]`, `{…}`, `<`, `>`, `-`) must ALWAYS convert to a
    /// COMPILABLE regex (the converter is total over the legal grammar — INV-conv), and running them
    /// against random proteins must never hang or emit a malformed match (INV-span / INV-finite).
    /// Determinism (INV) is pinned by re-running. This exercises the converter and scanner over a wide
    /// swath of grammatically-valid inputs, not just hand-built ones.
    /// (PROSITE_Pattern_Matching.md §2.4, §4.3; Pattern_Matching_Methods.md INV-01..INV-07.)
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Prosite_RandomLegalPatterns_AlwaysCompileAndScanWellFormed(CancellationToken token)
    {
        foreach (int seed in new[] { 11, 83, 271, 2026 })
        {
            var rng = new Random(seed);
            for (int trial = 0; trial < 40; trial++)
            {
                string pattern = BuildRandomLegalProsite(rng);

                // INV-conv — conversion is total over the legal grammar and yields a compilable regex.
                string regex = null!;
                var convert = () => regex = ConvertPrositeToRegex(pattern);
                convert.Should().NotThrow($"a grammatically-legal PROSITE pattern must convert (seed {seed}: \"{pattern}\")");
                AssertCompilableRegex(regex);

                // INV-span / INV-finite — scanning random proteins yields only well-formed spans.
                string protein = RandomProtein(rng.Next(1, 80), seed * 31 + trial);
                List<MotifMatch> hits = null!;
                var find = () => hits = FindMotifByProsite(protein, pattern, "Rand").ToList();
                find.Should().NotThrow($"scanning a random protein with \"{pattern}\" must not crash");
                token.ThrowIfCancellationRequested();
                foreach (var h in hits)
                    AssertWellFormedMatch(h, protein, pattern);

                // Determinism: re-converting and re-scanning yields identical results.
                ConvertPrositeToRegex(pattern).Should().Be(regex, "ConvertPrositeToRegex is deterministic");
                FindMotifByProsite(protein, pattern, "Rand")
                    .Select(h => (h.Start, h.End, h.Sequence))
                    .Should().Equal(hits.Select(h => (h.Start, h.End, h.Sequence)),
                        "FindMotifByProsite is deterministic for a fixed input");
            }
        }
    }

    /// <summary>
    /// Assembles a random but GRAMMATICALLY-LEGAL PROSITE PA-line from a fixed RNG: a `-`-separated
    /// chain of atoms drawn from { literal residue, x, x(n), x(n,m), [class], {exclusion} }, optionally
    /// anchored with `&lt;` / `&gt;`. Uses only constructs the §4.2 grammar table supports, so conversion
    /// must always succeed (no `*`/`?`/`+`).
    /// </summary>
    private static string BuildRandomLegalProsite(Random rng)
    {
        int atomCount = rng.Next(1, 7);
        var atoms = new List<string>();
        for (int i = 0; i < atomCount; i++)
        {
            switch (rng.Next(6))
            {
                case 0: // literal residue
                    atoms.Add(StandardAminoAcids[rng.Next(StandardAminoAcids.Length)].ToString());
                    break;
                case 1: // wildcard
                    atoms.Add("x");
                    break;
                case 2: // x(n)
                    atoms.Add($"x({rng.Next(1, 5)})");
                    break;
                case 3: // x(n,m)
                    int n = rng.Next(1, 4);
                    atoms.Add($"x({n},{n + rng.Next(1, 3)})");
                    break;
                case 4: // [class]
                    atoms.Add("[" + new string(Enumerable.Range(0, rng.Next(1, 4))
                        .Select(_ => StandardAminoAcids[rng.Next(StandardAminoAcids.Length)]).ToArray()) + "]");
                    break;
                default: // {exclusion}
                    atoms.Add("{" + new string(Enumerable.Range(0, rng.Next(1, 4))
                        .Select(_ => StandardAminoAcids[rng.Next(StandardAminoAcids.Length)]).ToArray()) + "}");
                    break;
            }
        }

        string body = string.Join("-", atoms);
        if (rng.Next(4) == 0) body = "<" + body;     // optional N-terminus anchor
        if (rng.Next(4) == 0) body += ">";           // optional C-terminus anchor
        return body;
    }

    #endregion

    #endregion
}
