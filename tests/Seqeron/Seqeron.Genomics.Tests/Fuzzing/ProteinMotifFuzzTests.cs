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

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the ProteinMotif area — general protein motif FINDING
/// (PROTMOTIF-FIND-001): locating short sequence motifs / patterns in a protein
/// string by regular-expression scan. The two public entry points under test are
/// the single-pattern scan <see cref="ProteinMotifFinder.FindMotifByPattern"/> and
/// the multi-pattern convenience scan <see cref="ProteinMotifFinder.FindCommonMotifs"/>
/// (which iterates the fixed in-source <c>CommonMotifs</c> dictionary and delegates
/// each entry to <c>FindMotifByPattern</c>). Sibling units PROTMOTIF-PROSITE-001
/// (row 83, PROSITE pattern syntax / regex injection) and PROTMOTIF-DOMAIN-001
/// (row 84, domain finding) are covered separately; this file focuses on the
/// GENERAL motif-finding contract.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output, and no *unhandled* runtime
/// exception (IndexOutOfRange / NullReference / ArgumentOutOfRange). Every input
/// must resolve to EITHER a well-defined, theory-correct result OR a *documented,
/// intentional* outcome. For a regex-based motif scanner that uppercases the input,
/// wraps the supplied pattern in a zero-width lookahead and walks every match span,
/// the headline hazards are:
///   • a NullReferenceException when the protein (or pattern) is null;
///   • an IndexOutOfRangeException when the protein is SHORTER than the motif, so
///     no match span exists (the scan must simply yield nothing, never run off the end);
///   • a runaway / quadratic scan on a homopolymer ("AAAA…") where a pattern can
///     match at every adjacent offset (the lookahead reports every overlapping hit
///     but must terminate in linear-in-matches time);
///   • a mis-reported match whose [Start..End] span does NOT actually equal the
///     substring it claims (a coordinate bug).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROTMOTIF-FIND-001 — protein motif finding (regex pattern scan)
/// Checklist: docs/checklists/03_FUZZING.md, row 82.
/// Fuzz strategies exercised for THIS unit:
///   • MC = Malformed Content — out-of-alphabet / junk residues:
///       – non-amino-acid characters (digits, punctuation, whitespace, lowercase,
///         the unknown placeholder 'X', the extended IUPAC codes B/Z/J/O/U): the
///         scanner does NOT validate the amino-acid alphabet beyond regex evaluation
///         (Motif_Search.md §3.3), so junk residues simply fail to match an
///         amino-acid pattern and yield NO hits — never a crash. Lowercase input is
///         uppercased first, so it scores identically to uppercase (§5.2, §6.1).
///   • BE = Boundary Exploitation — the length / composition corners that could
///     crash or hang:
///       – empty protein: "" / null → NO matches by the explicit guard
///         (string.IsNullOrEmpty), never a NullReference (§3.3, §6.1).
///       – extremely short: a protein SHORTER than the motif (1 residue vs a 4-residue
///         pattern, or a 0-length protein) → NO match, the regex simply finds nothing,
///         no IndexOutOfRange on the absent span (§6.1 "No motif occurrence").
///       – all-same char (homopolymer "AAAA…"): a defined outcome — a pattern present
///         in the homopolymer is matched at EVERY offset it occurs (overlapping hits
///         via the lookahead wrapper), a pattern absent from it yields nothing, and the
///         scan terminates (no quadratic hang) (§4.2, §5.2 "overlapping occurrences").
/// — docs/checklists/03_FUZZING.md §Description (strategy codes MC, BE);
///   targets: "Empty protein, non-amino acid chars, extremely short, all same char".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The motif-finding contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Given a protein sequence S and a regex motif pattern P, report every contiguous
/// subsequence of S that satisfies P. The implementation uppercases S, compiles
/// "(?=(P))" (a zero-width lookahead capturing the body) so OVERLAPPING occurrences
/// are discovered, and emits one MotifMatch per capture: { Start (incl. 0-based),
/// End (incl. 0-based), Sequence (the uppercased matched substring), MotifName,
/// Pattern, Score (information content), EValue }. An invalid regex is swallowed and
/// yields nothing rather than throwing (§4.1, §5.2).
///   — docs/algorithms/ProteinMotif/Motif_Search.md §2.2, §3.1–§3.3, §4.1–§4.2.
///
/// Methods under test (src/.../Seqeron.Genomics.Analysis/ProteinMotifFinder.cs):
///   IEnumerable&lt;MotifMatch&gt; FindMotifByPattern(string proteinSequence, string regexPattern,
///       string motifName = "Custom", string patternId = "")  — single-pattern scan.
///   IEnumerable&lt;MotifMatch&gt; FindCommonMotifs(string proteinSequence)  — scans the fixed
///       CommonMotifs dictionary, delegating each entry to FindMotifByPattern.
///   — Motif_Search.md §5.1.
///
/// Documented input handling (Motif_Search.md §3.3, §6.1):
///   • null / "" protein → NO matches (explicit string.IsNullOrEmpty guard), no throw.
///   • null / "" / invalid regexPattern → FindMotifByPattern yields NO matches (guarded
///     or the compile failure is swallowed).
///   • Matching is case-INSENSITIVE: the input is uppercased and the regex is compiled
///     with RegexOptions.IgnoreCase; lowercase == uppercase.
///   • The amino-acid alphabet is NOT validated beyond regex evaluation — non-AA characters
///     (digits, punctuation, X, B/Z/J/O/U) simply fail to satisfy an amino-acid pattern.
///   • Coordinates are inclusive 0-based; Sequence is the exact uppercased substring at
///     [Start..End].
///
/// Theory-correct invariants asserted (Motif_Search.md §2.4):
///   • INV-01 — deterministic: re-scanning the same sequence + pattern yields identical hits.
///   • INV-02 — every reported match is a CONTIGUOUS subsequence: S[Start..End] equals
///     Sequence, with 0 ≤ Start ≤ End ≤ n−1 (the headline "no run-off-the-end / no
///     coordinate bug" property).
///   • [span-shape] — End − Start + 1 == Sequence.Length, and Sequence is uppercased.
///   • [finite-score] — Score and EValue are finite (never NaN / ±∞).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Complexity / hang-safety
/// ───────────────────────────────────────────────────────────────────────────
/// The scan is a single regex walk over n residues returning k matches (Motif_Search.md
/// §4.3). The homopolymer and long-junk targets maximise the number of overlapping hits
/// (a short pattern can match at almost every offset of "AAAA…"); they are kept modest and
/// [CancelAfter]-guarded so a regression that turned the lookahead scan into a hang or a
/// super-linear blow-up would FAIL rather than wedge the suite.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ProteinMotifFuzzTests
{
    #region Helpers

    /// <summary>The 20 standard amino-acid one-letter codes.</summary>
    private const string StandardAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
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
    /// satisfy against the original (case-insensitive) sequence (Motif_Search.md §2.4, §3.2):
    /// the match is a CONTIGUOUS in-bounds subsequence (INV-02) whose claimed coordinates actually
    /// reproduce the reported substring ([span-shape]: S[Start..End] == Sequence, uppercased), and
    /// whose Score / EValue are finite ([finite-score]). This is the headline "no coordinate bug,
    /// no run-off-the-end, no NaN" property.
    /// </summary>
    private static void AssertWellFormedMatch(MotifMatch match, string originalSequence)
    {
        string upper = originalSequence.ToUpperInvariant();
        int n = upper.Length;

        // INV-02 — in-bounds, non-empty, contiguous span.
        match.Start.Should().BeInRange(0, n - 1, "a match Start is a valid 0-based residue index");
        match.End.Should().BeInRange(match.Start, n - 1, "a match End is in-bounds and not before its Start");

        // [span-shape] — the span length equals the reported substring length, which is uppercased
        // and equals the actual substring of the (uppercased) input at [Start..End].
        (match.End - match.Start + 1).Should().Be(match.Sequence.Length,
            "End − Start + 1 equals the matched substring length");
        match.Sequence.Should().Be(match.Sequence.ToUpperInvariant(), "the matched substring is uppercased");
        upper.Substring(match.Start, match.Sequence.Length).Should().Be(match.Sequence,
            "INV-02: the reported substring is exactly S[Start..End] of the uppercased input (no coordinate bug)");

        // [finite-score] — score and E-value are finite.
        double.IsNaN(match.Score).Should().BeFalse("a motif Score must never be NaN");
        double.IsInfinity(match.Score).Should().BeFalse("a motif Score must never be infinite");
        double.IsNaN(match.EValue).Should().BeFalse("a motif EValue must never be NaN");
        double.IsInfinity(match.EValue).Should().BeFalse("a motif EValue must never be infinite");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PROTMOTIF-FIND-001 — protein motif finding : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PROTMOTIF-FIND-001 — protein motif finding

    #region BE — Empty / null protein: no matches, no NullReference

    /// <summary>
    /// Target "Empty protein": "" and null must produce NO matches — by the explicit
    /// string.IsNullOrEmpty guard, NEVER a NullReferenceException — for BOTH the single-pattern
    /// scan and the CommonMotifs convenience scan. This is the headline no-crash contract on the
    /// degenerate empty input (Motif_Search.md §3.3, §6.1 "Null or empty proteinSequence →
    /// returns no matches"). A null/empty/invalid regex on a non-empty protein likewise yields
    /// nothing rather than throwing (§3.3, §6.1).
    /// </summary>
    [Test]
    public void FindMotif_EmptyOrNullProtein_NoMatchesNoThrow()
    {
        foreach (string? seq in new[] { "", null })
        {
            var byPattern = () => FindMotifByPattern(seq!, "RGD").ToList();
            byPattern.Should().NotThrow($"empty/null protein ('{seq ?? "null"}') must not crash FindMotifByPattern")
                .Subject.Should().BeEmpty("empty/null protein yields no motif matches");

            var common = () => FindCommonMotifs(seq!).ToList();
            common.Should().NotThrow($"empty/null protein ('{seq ?? "null"}') must not crash FindCommonMotifs")
                .Subject.Should().BeEmpty("empty/null protein yields no common-motif matches");
        }

        // A null / empty / invalid regex on a real protein → no matches, no throw (§3.3).
        foreach (string? pattern in new[] { "", null })
        {
            var act = () => FindMotifByPattern("MKRGDSPEKWFIL", pattern!).ToList();
            act.Should().NotThrow($"null/empty regex ('{pattern ?? "null"}') must not crash")
                .Subject.Should().BeEmpty("a null/empty regex yields no matches");
        }

        // Invalid regex (unbalanced bracket) is swallowed → no matches, no throw (§5.2).
        var invalid = () => FindMotifByPattern("MKRGDSPEKWFIL", "[ABC").ToList();
        invalid.Should().NotThrow("an invalid regex must be swallowed, not thrown")
            .Subject.Should().BeEmpty("an invalid regex yields no matches");
    }

    #endregion

    #region BE — Extremely short protein: shorter than the motif → no match, no crash

    /// <summary>
    /// Target "extremely short": a protein STRICTLY SHORTER than the motif it is scanned for has
    /// no place for the pattern to fit, so the scan must yield NO match and NEVER run off the end
    /// (no IndexOutOfRange on the absent span) — Motif_Search.md §6.1 "No motif occurrence in the
    /// sequence → empty collection". We probe a 1-residue protein against a 4-residue motif and a
    /// 3-residue motif, plus every single standard residue scanned for the multi-residue RGD /
    /// PPxY motifs (all too short). A single residue scanned for a 1-residue motif that DOES match
    /// is the positive boundary: exactly one hit at position 0 covering the whole 1-char string.
    /// </summary>
    [Test]
    public void FindMotif_ProteinShorterThanMotif_NoMatchNoCrash()
    {
        // (a) 1-residue protein vs multi-residue motifs → no match, no crash.
        foreach (char aa in StandardAminoAcids)
        {
            string seq = aa.ToString();
            foreach (string pattern in new[] { "PP.Y", "RGD", "[RK]{2}.[ST]" })
            {
                var act = () => FindMotifByPattern(seq, pattern).ToList();
                act.Should().NotThrow($"a 1-residue protein ('{aa}') vs motif '{pattern}' must not crash")
                    .Subject.Should().BeEmpty($"a 1-residue protein cannot contain the multi-residue motif '{pattern}'");
            }
        }

        // (b) A protein shorter than a multi-residue PROSITE motif yields nothing (FindCommonMotifs too).
        foreach (string seq in new[] { "M", "MK", "RGD"[..2] /* "RG" */ })
        {
            var act = () => FindCommonMotifs(seq).ToList();
            var hits = act.Should().NotThrow($"a short protein ('{seq}') must not crash FindCommonMotifs").Subject;
            foreach (var hit in hits)
                AssertWellFormedMatch(hit, seq); // any hit (e.g. a 1-wide pattern) must still be well-formed
        }

        // (c) Positive boundary: a 1-residue protein scanned for a 1-residue motif that matches →
        //     exactly one hit covering the whole string.
        var single = FindMotifByPattern("R", "R").ToList();
        single.Should().ContainSingle("a 1-residue protein contains the 1-residue motif exactly once");
        AssertWellFormedMatch(single[0], "R");
        single[0].Start.Should().Be(0, "the 1-residue match starts at residue 0");
        single[0].End.Should().Be(0, "the 1-residue match ends at residue 0");
        single[0].Sequence.Should().Be("R", "the matched substring is the whole 1-char protein");

        // And a 1-residue protein scanned for a DIFFERENT 1-residue motif → no match.
        FindMotifByPattern("R", "K").Should().BeEmpty("the 1-residue protein 'R' does not contain motif 'K'");
    }

    #endregion

    #region MC — Non-amino-acid characters: handled (no match), never a crash

    /// <summary>
    /// Target "non-amino acid chars": the scanner does NOT validate the amino-acid alphabet beyond
    /// regex evaluation (Motif_Search.md §3.3), so out-of-alphabet residues — digits, punctuation,
    /// whitespace, the unknown placeholder 'X', and the extended IUPAC codes B/Z/J/O/U — simply fail
    /// to satisfy an amino-acid pattern and yield NO hits, never a crash. Any hit that IS emitted
    /// (e.g. the '.' wildcard 'x' position of a PROSITE pattern legitimately spanning a junk char)
    /// must still be a well-formed, in-bounds, contiguous span. Lowercase input is uppercased first,
    /// so it matches identically to its uppercase form (§5.2, §6.1).
    /// </summary>
    [Test]
    public void FindMotif_NonAminoAcidChars_HandledNoCrash()
    {
        // (a) Pure junk scanned for an amino-acid motif → no match, no crash.
        foreach (string junk in new[]
                 {
                     "1234567890",          // digits
                     "!@#$%^&*()",          // punctuation
                     "   \t  \n  ",         // whitespace only
                     "BZJOUBZJOU",          // extended IUPAC ambiguity codes
                     "XXXXXXXXXX",          // unknown placeholder
                 })
        {
            foreach (string pattern in new[] { "RGD", "PP.Y", "N[^P][ST][^P]" })
            {
                var act = () => FindMotifByPattern(junk, pattern).ToList();
                var hits = act.Should().NotThrow(
                    $"junk protein ('{junk.Trim()}') vs motif '{pattern}' must not crash").Subject;
                hits.Should().BeEmpty($"out-of-alphabet residues cannot satisfy the amino-acid motif '{pattern}'");
            }

            // FindCommonMotifs over pure junk: never throws; any emitted hit is well-formed.
            var commonHits = ((Func<List<MotifMatch>>)(() => FindCommonMotifs(junk).ToList()))
                .Should().NotThrow($"junk protein ('{junk.Trim()}') must not crash FindCommonMotifs").Subject;
            foreach (var hit in commonHits)
                AssertWellFormedMatch(hit, junk);
        }

        // (b) Junk embedded among real residues: the scan must still locate the genuine motif and
        //     report a well-formed span; the junk neither crashes the scan nor corrupts coordinates.
        foreach (string seq in new[]
                 {
                     "MK1RGD2SP",           // RGD flanked by digits
                     "PP X Y PPSY",         // junk-separated, plus a real PP.Y near the end
                     "rgdRGDrgd",           // lowercase + uppercase RGD (case-insensitive)
                 })
        {
            var hits = ((Func<List<MotifMatch>>)(() => FindMotifByPattern(seq, "RGD", "RGD").ToList()))
                .Should().NotThrow($"motif scan over '{seq}' with embedded junk must not crash").Subject;
            foreach (var hit in hits)
            {
                AssertWellFormedMatch(hit, seq);
                hit.Sequence.Should().Be("RGD", "an RGD hit is exactly the RGD substring");
            }
        }

        // (c) Case-insensitivity: lowercase protein matches identically to uppercase (§5.2, §6.1).
        var lower = FindMotifByPattern("mkrgdsp", "RGD").ToList();
        var upper = FindMotifByPattern("MKRGDSP", "RGD").ToList();
        lower.Select(m => (m.Start, m.End, m.Sequence))
            .Should().Equal(upper.Select(m => (m.Start, m.End, m.Sequence)),
                "matching is case-insensitive: lowercase input yields the same hits as uppercase");
    }

    #endregion

    #region BE — All-same char (homopolymer): defined outcome, overlapping hits, no hang

    /// <summary>
    /// Target "all same char": a homopolymer "AAAA…" has a DEFINED motif-finding outcome — a pattern
    /// present in the homopolymer is matched at EVERY offset it occurs (overlapping hits via the
    /// lookahead wrapper, Motif_Search.md §4.2, §5.2), a pattern absent from it yields nothing, and
    /// the scan must TERMINATE in time linear in the number of matches (no quadratic hang on the
    /// dense overlapping-match case). For a homopolymer of length n scanned for a single-residue
    /// motif equal to that residue, the lookahead reports exactly n overlapping hits (one per
    /// position); for a motif of width w of that residue it reports n − w + 1 hits. Every hit must
    /// be a well-formed, contiguous span made entirely of the homopolymer residue.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindMotif_Homopolymer_DefinedOverlappingHitsNoHang(CancellationToken token)
    {
        foreach (char aa in StandardAminoAcids)
        {
            const int n = 100;
            string seq = new string(aa, n);

            // (a) Single-residue motif equal to the homopolymer residue → one hit per position.
            var singleHits = ((Func<List<MotifMatch>>)(() => FindMotifByPattern(seq, aa.ToString()).ToList()))
                .Should().NotThrow($"a homopolymer of '{aa}' scanned for '{aa}' must not crash").Subject;
            token.ThrowIfCancellationRequested();

            singleHits.Should().HaveCount(n,
                $"a length-{n} homopolymer of '{aa}' contains the 1-residue motif at every position (overlapping)");
            for (int i = 0; i < n; i++)
            {
                AssertWellFormedMatch(singleHits[i], seq);
                singleHits[i].Start.Should().Be(i, "overlapping hits are reported at every consecutive offset");
                singleHits[i].Sequence.Should().Be(aa.ToString(), "each hit is the single homopolymer residue");
            }

            // (b) Width-4 motif of that residue → n − 4 + 1 overlapping hits, all the same substring.
            string quad = new string(aa, 4);
            string quadPattern = $"{aa}{{4}}";
            var quadHits = ((Func<List<MotifMatch>>)(() => FindMotifByPattern(seq, quadPattern).ToList()))
                .Should().NotThrow($"a homopolymer of '{aa}' scanned for '{quadPattern}' must not crash").Subject;
            token.ThrowIfCancellationRequested();

            quadHits.Should().HaveCount(n - 3,
                $"a length-{n} homopolymer contains the width-4 motif at {n - 3} overlapping offsets");
            foreach (var hit in quadHits)
            {
                AssertWellFormedMatch(hit, seq);
                hit.Sequence.Should().Be(quad, "each width-4 hit is the 4-residue homopolymer substring");
            }

            // (c) A motif that REQUIRES a different residue is absent → no hits.
            char other = aa == 'A' ? 'C' : 'A';
            FindMotifByPattern(seq, other.ToString()).Should().BeEmpty(
                $"a homopolymer of '{aa}' does not contain the foreign single-residue motif '{other}'");
        }

        // FindCommonMotifs over a long homopolymer must terminate and emit only well-formed hits.
        foreach (char aa in new[] { 'A', 'L', 'G', 'P' })
        {
            string seq = new string(aa, 200);
            var hits = ((Func<List<MotifMatch>>)(() => FindCommonMotifs(seq).ToList()))
                .Should().NotThrow($"FindCommonMotifs over a 200-long '{aa}' homopolymer must not hang/crash").Subject;
            token.ThrowIfCancellationRequested();
            foreach (var hit in hits)
                AssertWellFormedMatch(hit, seq);
        }
    }

    #endregion

    #region Positive sanity — a known motif is found at the right position

    /// <summary>
    /// Positive sanity: the harness must assert against a scanner that actually FINDS motifs at the
    /// correct coordinates, not a no-op. A protein containing a single, unambiguous RGD cell-attachment
    /// motif (PROSITE PS00016) at a KNOWN offset must yield exactly one hit whose Start/End/Sequence
    /// pin that occurrence; a protein with two well-separated RGD motifs must yield exactly two hits at
    /// the two known offsets; and the WW-domain PY motif "PP.Y" (a wildcard-bearing pattern) must be
    /// located with the wildcard residue correctly captured. FindCommonMotifs must additionally surface
    /// the RGD hit (its accession is in CommonMotifs). This pins the headline "a found motif actually
    /// occurs at the reported position" contract (Motif_Search.md §2.4 INV-02, §4.1).
    /// </summary>
    [Test]
    public void FindMotif_KnownMotif_FoundAtCorrectPosition()
    {
        // Single RGD at offset 4: "MKSP RGD WFIL".
        const string oneRgd = "MKSPRGDWFIL";
        int expectedStart = oneRgd.IndexOf("RGD", StringComparison.Ordinal); // 4
        expectedStart.Should().Be(4, "sanity-check the hand-built offset");

        var rgdHits = FindMotifByPattern(oneRgd, "RGD", "RGD", "PS00016").ToList();
        rgdHits.Should().ContainSingle("the protein contains exactly one RGD motif");
        AssertWellFormedMatch(rgdHits[0], oneRgd);
        rgdHits[0].Start.Should().Be(expectedStart, "the RGD hit is reported at the known offset 4");
        rgdHits[0].End.Should().Be(expectedStart + 2, "RGD spans three residues [4,6]");
        rgdHits[0].Sequence.Should().Be("RGD", "the matched substring is exactly RGD");

        // Two well-separated RGD motifs.
        const string twoRgd = "RGDKKKKKKRGD";
        var twoHits = FindMotifByPattern(twoRgd, "RGD").ToList();
        twoHits.Should().HaveCount(2, "two well-separated RGD motifs yield two hits");
        twoHits.Select(h => h.Start).Should().Equal(new[] { 0, 9 }, "the two RGD hits are at the two known offsets");
        foreach (var hit in twoHits)
        {
            AssertWellFormedMatch(hit, twoRgd);
            hit.Sequence.Should().Be("RGD");
        }

        // Wildcard-bearing PY motif "PP.Y" (PROSITE WW1 / PPxY): the '.' captures the middle residue.
        const string pyProtein = "AAPPSYAA"; // PP S Y at offset 2
        var pyHits = FindMotifByPattern(pyProtein, "PP.Y", "WW_BINDING").ToList();
        pyHits.Should().ContainSingle("the protein contains exactly one PPxY motif");
        AssertWellFormedMatch(pyHits[0], pyProtein);
        pyHits[0].Start.Should().Be(2, "PPSY begins at offset 2");
        pyHits[0].Sequence.Should().Be("PPSY", "the wildcard captures the middle residue 'S'");

        // FindCommonMotifs must also surface the RGD occurrence (PS00016 is in CommonMotifs).
        var common = FindCommonMotifs(oneRgd).ToList();
        common.Should().Contain(m => m.Sequence == "RGD" && m.Start == expectedStart,
            "FindCommonMotifs surfaces the RGD motif from the CommonMotifs dictionary at the correct offset");
        foreach (var hit in common)
            AssertWellFormedMatch(hit, oneRgd);
    }

    /// <summary>
    /// Positive sanity over RANDOM proteins: across fixed seeds and lengths the scan must never crash,
    /// hang, or emit a malformed match, and every emitted MotifMatch must satisfy the full contract
    /// (in-bounds contiguous span whose coordinates reproduce its substring, finite score / E-value).
    /// INV-01 determinism is pinned by re-scanning the same input and requiring identical hits. This
    /// pins span-correctness and termination on arbitrary sequences, not just hand-built motifs.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindMotif_RandomProtein_AlwaysWellFormedAndDeterministic(CancellationToken token)
    {
        foreach (int seed in new[] { 7, 31, 137, 2026 })
        {
            foreach (int len in new[] { 1, 5, 21, 60, 250 })
            {
                string seq = RandomProtein(len, seed);

                var act = () => FindCommonMotifs(seq).ToList();
                var hits = act.Should().NotThrow($"random protein must not crash (seed {seed}, len {len})").Subject;
                token.ThrowIfCancellationRequested();

                foreach (var hit in hits)
                    AssertWellFormedMatch(hit, seq);

                // INV-01 — deterministic: the same input yields identical hits.
                var again = FindCommonMotifs(seq).ToList();
                again.Select(h => (h.MotifName, h.Start, h.End, h.Sequence))
                    .Should().Equal(hits.Select(h => (h.MotifName, h.Start, h.End, h.Sequence)),
                        "INV-01: FindCommonMotifs is deterministic for a fixed input");
            }
        }
    }

    #endregion

    #endregion
}
