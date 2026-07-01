using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Matching area — Regulatory Element Scan (MOTIF-REGULATORY-001),
/// the fixed-library consensus scanner
/// <see cref="MotifFinder.FindRegulatoryElements(DnaSequence)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and adversarial inputs to a unit and asserts
/// the code NEVER fails in an undisciplined way: no NullReference / IndexOutOfRange,
/// no hang, no FALSE HIT (a reported element whose matched substring does not actually
/// occur at the reported offset or does not IUPAC-match its consensus), and no MISSED
/// occurrence (every offset that satisfies the consensus, including OVERLAPPING ones,
/// is reported). The unit is specification-driven: it scans a FIXED LIBRARY of 12
/// published regulatory consensus strings (§2.2) and reports EVERY 0-based start
/// `i` with `0 <= i <= n-m` where the window IUPAC-matches the consensus (INV-01..03).
/// The boundaries of such a scanner are exactly: EMPTY (empty sequence → empty result,
/// no offset satisfies `0 <= i <= n-m`, no NullReference; a longer-than-sequence pattern
/// → no IndexOutOfRange because the loop never executes), NO ELEMENT (a sequence
/// containing none of the 12 consensus strings → empty result, no false hit), and
/// OVERLAPPING (a sequence where one element recurs at overlapping offsets, or two
/// different library elements overlap → ALL occurrences reported with correct
/// coordinates, none missed or duplicated). Every input must resolve to EITHER a
/// well-defined, theory-correct result OR the single documented validation exception
/// (ArgumentNullException for a null sequence — §3.3, §6.1). A raw runtime exception,
/// a hang, a false hit, a missed overlap, or an off-by-one coordinate is a bug, not a
/// passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MOTIF-REGULATORY-001 — Regulatory Element Scan
/// Checklist: docs/checklists/03_FUZZING.md, row 172.
/// Algorithm doc: docs/algorithms/Motif_Discovery/Regulatory_Elements.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the boundaries called out in the row
///     ("empty, no element, overlapping"):
///       – EMPTY: an empty sequence → empty result with NO NullReference (§6.1); a
///         library pattern (up to 13 bp, the Kozak string) longer than the whole
///         sequence → that pattern contributes no hit and NO IndexOutOfRange because
///         the scan loop bound `i <= n-m` is negative and never runs (§4.1).
///       – NO ELEMENT: a sequence built from a single base, or random DNA verified to
///         contain none of the 12 consensus strings → empty result, no false hit
///         (INV-02/INV-03 the contrapositive: nothing reported that does not match).
///       – OVERLAPPING: a sequence where one element recurs at OVERLAPPING offsets
///         (e.g. the E-box `CANNTG` family in a homopolymer-ish run, or `AATAAA`
///         poly(A) tiled so windows overlap) and where two DIFFERENT library elements
///         share bases → EVERY satisfying offset is reported exactly once, with the
///         correct 0-based position, no missed or duplicated overlap (INV-03
///         exhaustiveness).
/// — docs/checklists/03_FUZZING.md §Description (BE = граничні значення 0/-1/MaxInt/empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Regulatory_Elements.md §2.2, §2.4, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
/// For a sequence S of length n and a library pattern P of length m, report every
/// start index i with 0 <= i <= n-m such that for all j, S[i+j] is in the IUPAC base
/// set of P[j] (plain bases match themselves; N matches A/C/G/T). The 12 library
/// entries (5'→3') and their reported Name are:
///   TATA Box       TATAAA          | CAAT Box       CCAAT
///   GC Box         GGGCGG          | -10 Box        TATAAT
///   -35 Box        TTGACA          | Kozak          GCCGCCACCATGG
///   Shine-Dalgarno AGGAGG          | Poly(A) Signal AATAAA
///   E-box          CANNTG (IUPAC)  | AP-1           TGACTCA
///   NF-κB          GGGACTTTCC      | CREB           TGACGTCA
/// Each hit reports Name, 0-based Position, the matched Sequence (length = m),
/// the Pattern, and a Description (INV-01). Each matched Sequence IUPAC-matches its
/// Pattern (INV-02). The scan is exhaustive over every offset (INV-03). Results are
/// grouped library-entry by library-entry, increasing position within each entry
/// (§3.2, §5.2). Null sequence → ArgumentNullException; empty sequence → empty result
/// (§3.3, §6.1). Only the given strand is scanned (§6.2).
///   MotifFinder.FindRegulatoryElements(DnaSequence) → IEnumerable&lt;RegulatoryElement&gt;
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MotifRegulatoryFuzzTests
{
    private static readonly char[] Alphabet = { 'A', 'C', 'G', 'T' };

    /// <summary>
    /// The fixed library (§2.2) as (Name, Pattern) pairs, copied verbatim from the
    /// algorithm doc table. This is the independent oracle's reference catalog; the
    /// test owns it so a drift in the source library is caught, not masked.
    /// </summary>
    private static readonly (string Name, string Pattern)[] Library =
    {
        ("TATA Box", "TATAAA"),
        ("CAAT Box", "CCAAT"),
        ("GC Box", "GGGCGG"),
        ("-10 Box", "TATAAT"),
        ("-35 Box", "TTGACA"),
        ("Kozak", "GCCGCCACCATGG"),
        ("Shine-Dalgarno", "AGGAGG"),
        ("Poly(A) Signal", "AATAAA"),
        ("E-box", "CANNTG"),
        ("AP-1", "TGACTCA"),
        ("NF-κB", "GGGACTTTCC"),
        ("CREB", "TGACGTCA"),
    };

    #region Helpers

    /// <summary>
    /// True if base <paramref name="seqChar"/> is in the IUPAC set of consensus
    /// symbol <paramref name="patChar"/> (the library uses only plain bases and N).
    /// </summary>
    private static bool IupacMatch(char patChar, char seqChar) => patChar switch
    {
        'A' or 'C' or 'G' or 'T' => patChar == seqChar,
        'N' => seqChar is 'A' or 'C' or 'G' or 'T',
        _ => throw new InvalidOperationException($"Library pattern uses unexpected symbol '{patChar}'."),
    };

    /// <summary>
    /// Independent oracle implementing the documented scan verbatim
    /// (Regulatory_Elements.md §2.2, §4.1) WITHOUT re-using the unit: for each library
    /// entry, report every 0-based start i with 0 &lt;= i &lt;= n-m whose window
    /// IUPAC-matches the consensus. Returns (Name, Position) tuples grouped entry by
    /// entry, increasing position — the documented yield order (§3.2, §5.2).
    /// </summary>
    private static List<(string Name, int Position)> Oracle(string seq)
    {
        var hits = new List<(string, int)>();
        if (string.IsNullOrEmpty(seq)) return hits;

        foreach (var (name, pattern) in Library)
        {
            int m = pattern.Length;
            for (int i = 0; i <= seq.Length - m; i++) // negative bound when m > n: no iterations
            {
                bool ok = true;
                for (int j = 0; j < m && ok; j++)
                    ok = IupacMatch(pattern[j], seq[i + j]);
                if (ok) hits.Add((name, i));
            }
        }
        return hits;
    }

    /// <summary>
    /// Asserts a RESULT SET is WELL-FORMED per the documented contract regardless of the
    /// (possibly degenerate) input: INV-01 each hit's matched Sequence has length ==
    /// its Pattern and equals S[Position..Position+m); coordinates are in bounds and
    /// 0-based; INV-02 the matched Sequence IUPAC-matches its Pattern; Name/Pattern are
    /// drawn from the fixed library (no fabricated entry).
    /// </summary>
    private static void AssertWellFormed(IReadOnlyList<RegulatoryElement> hits, string seq)
    {
        var byName = Library.ToDictionary(e => e.Name, e => e.Pattern);
        foreach (var h in hits)
        {
            byName.Should().ContainKey(h.Name, "every reported element is from the fixed library (§2.2)");
            h.Pattern.Should().Be(byName[h.Name], "INV-04: the reported pattern is the library consensus");

            h.Sequence.Should().HaveLength(h.Pattern.Length, "INV-01: matched length == pattern length");
            h.Position.Should().BeInRange(0, seq.Length - h.Pattern.Length,
                "INV-01: 0-based start in [0, n-m], so the window is fully in bounds");
            seq.Substring(h.Position, h.Pattern.Length).Should().Be(h.Sequence,
                "INV-01: the matched Sequence is exactly S[Position..Position+m)");

            for (int j = 0; j < h.Pattern.Length; j++)
                IupacMatch(h.Pattern[j], h.Sequence[j]).Should().BeTrue(
                    "INV-02: each matched base is in the IUPAC set of the consensus symbol");
        }
    }

    private static string RandomDna(Random rng, int len)
    {
        var sb = new StringBuilder(len);
        for (int i = 0; i < len; i++) sb.Append(Alphabet[rng.Next(Alphabet.Length)]);
        return sb.ToString();
    }

    private static List<(string Name, int Position)> Actual(IEnumerable<RegulatoryElement> hits) =>
        hits.Select(h => (h.Name, h.Position)).ToList();

    #endregion

    #region MOTIF-REGULATORY-001 — Regulatory Element Scan

    // ─── BE: positive sanity — known elements at known offsets ──────────────────

    [Test]
    [Category("Fuzzing")]
    public void Sanity_TataBox_LocatedAtKnownOffset_WithCorrectCoordinates()
    {
        // "GGG" + "TATAAA" + "GGG" → the only library hit is the TATA Box at offset 3
        // (the doc's worked example, §7.1). GGG GGG is not GGGCGG (GC box), so no GC box.
        var seq = new DnaSequence("GGGTATAAAGGG");

        var hits = MotifFinder.FindRegulatoryElements(seq).ToList();

        AssertWellFormed(hits, seq.Sequence);
        hits.Should().ContainSingle();
        hits[0].Name.Should().Be("TATA Box");
        hits[0].Position.Should().Be(3);
        hits[0].Sequence.Should().Be("TATAAA");
        hits[0].Pattern.Should().Be("TATAAA");
    }

    [Test]
    [Category("Fuzzing")]
    public void Sanity_EBoxDegenerate_MatchedViaIupacN()
    {
        // CACGTG and CAGCTG are both members of the E-box family CANNTG (§6.1 example).
        var seq = new DnaSequence("AAACACGTGAAA");

        var hits = MotifFinder.FindRegulatoryElements(seq).ToList();

        AssertWellFormed(hits, seq.Sequence);
        hits.Should().Contain(h => h.Name == "E-box" && h.Position == 3 && h.Sequence == "CACGTG");
    }

    // ─── BE: empty / null ───────────────────────────────────────────────────────

    [Test]
    [Category("Fuzzing")]
    public void Be_EmptySequence_YieldsNoElements_NoCrash()
    {
        var seq = new DnaSequence(string.Empty); // DnaSequence normalises "" / null → empty

        var hits = MotifFinder.FindRegulatoryElements(seq).ToList();

        hits.Should().BeEmpty("no offset satisfies 0 <= i <= n-m when n = 0 (§6.1)");
    }

    [Test]
    [Category("Fuzzing")]
    public void Be_NullSequence_ThrowsArgumentNullException()
    {
        // §3.3 / §6.1: null sequence → ArgumentNullException. The method is an iterator,
        // so the guard fires on enumeration.
        Action act = () => MotifFinder.FindRegulatoryElements(null!).ToList();

        act.Should().Throw<ArgumentNullException>();
    }

    [Test]
    [Category("Fuzzing")]
    [CancelAfter(30000)]
    public void Be_SequenceShorterThanEveryPattern_YieldsNothing_NoIndexOutOfRange()
    {
        // Shortest library pattern is 5 bp (CCAAT); longest is 13 bp (Kozak). A 1..4 bp
        // sequence is shorter than all of them → the scan bound i <= n-m is negative for
        // every entry → no iterations, no IndexOutOfRange.
        var rng = new Random(172_001);
        for (int len = 1; len <= 4; len++)
        {
            for (int t = 0; t < 50; t++)
            {
                var seq = new DnaSequence(RandomDna(rng, len));
                Action act = () =>
                {
                    var hits = MotifFinder.FindRegulatoryElements(seq).ToList();
                    hits.Should().BeEmpty();
                    AssertWellFormed(hits, seq.Sequence);
                };
                act.Should().NotThrow($"a {len}-bp sequence is shorter than every library pattern");
            }
        }
    }

    [Test]
    [Category("Fuzzing")]
    public void Be_KozakLongestPattern_OneBaseShort_YieldsNoKozak_NoCrash()
    {
        // Kozak = GCCGCCACCATGG (13 bp). A 12-bp prefix is one base short: the scan bound
        // for that entry is i <= -1 → no iteration → no Kozak hit, no IndexOutOfRange.
        var seq = new DnaSequence("GCCGCCACCATG"); // 12 bp

        var hits = MotifFinder.FindRegulatoryElements(seq).ToList();

        AssertWellFormed(hits, seq.Sequence);
        hits.Should().NotContain(h => h.Name == "Kozak");
    }

    // ─── BE: no element — no false hit ──────────────────────────────────────────

    [Test]
    [Category("Fuzzing")]
    public void Be_Homopolymer_HasNoLibraryElement_YieldsEmpty()
    {
        // A run of a single base cannot contain any of the 12 mixed-base consensus
        // strings → empty result, no false hit.
        foreach (char b in Alphabet)
        {
            var seq = new DnaSequence(new string(b, 60));
            var hits = MotifFinder.FindRegulatoryElements(seq).ToList();
            hits.Should().BeEmpty($"a poly-{b} run contains none of the library consensus strings");
        }
    }

    [Test]
    [Category("Fuzzing")]
    [CancelAfter(60000)]
    public void Be_RandomDna_NoFalseHits_MatchesIndependentOracle()
    {
        // Differential against the independent oracle on random DNA: whatever the unit
        // reports, it must be EXACTLY the set of offsets that genuinely IUPAC-match a
        // library consensus — no false hit, no missed hit, identical ordering.
        var rng = new Random(172_002);
        for (int t = 0; t < 600; t++)
        {
            int len = rng.Next(0, 80);
            string raw = RandomDna(rng, len);
            var seq = new DnaSequence(raw);

            var hits = MotifFinder.FindRegulatoryElements(seq).ToList();

            AssertWellFormed(hits, seq.Sequence);
            Actual(hits).Should().Equal(Oracle(seq.Sequence),
                "the unit reports exactly the documented matching offsets, in the documented order");
        }
    }

    // ─── BE: overlapping occurrences ────────────────────────────────────────────

    [Test]
    [Category("Fuzzing")]
    public void Be_PolyA_OverlappingOccurrences_AllReported()
    {
        // "AATAAAATAAA": the poly(A) signal AATAAA (6 bp) occurs at offset 0 and again
        // at offset 5 (overlapping the trailing AAA of the first window region). The scan
        // is exhaustive over every offset (INV-03), so BOTH must be reported.
        var seq = new DnaSequence("AATAAAATAAA");

        var hits = MotifFinder.FindRegulatoryElements(seq).ToList();
        AssertWellFormed(hits, seq.Sequence);

        var polyA = hits.Where(h => h.Name == "Poly(A) Signal").Select(h => h.Position).ToList();
        polyA.Should().Equal(Oracle(seq.Sequence).Where(o => o.Name == "Poly(A) Signal").Select(o => o.Position),
            "every overlapping poly(A) occurrence is reported, none missed or duplicated");
        polyA.Should().Contain(0).And.Contain(5);
    }

    [Test]
    [Category("Fuzzing")]
    public void Be_EBoxRun_OverlappingDegenerateOccurrences_AllReported()
    {
        // The E-box CANNTG = C, two free bases (N), then TG. "CACATGTG" has two E-box
        // windows whose 6-bp spans OVERLAP:
        //   off0: C A C A T G → CANNTG (N=C,A) then TG at 4,5 → MATCH, window [0..5]
        //   off2: C A T G T G → CANNTG (N=T,G) then TG at 6,7 → MATCH, window [2..7]
        // Windows [0..5] and [2..7] share bases 2..5, so this is a genuine overlap. The
        // exhaustive scan (INV-03) must report BOTH; the unit must equal the oracle set.
        var seq = new DnaSequence("CACATGTG");

        var hits = MotifFinder.FindRegulatoryElements(seq).ToList();
        AssertWellFormed(hits, seq.Sequence);

        var ebox = hits.Where(h => h.Name == "E-box").Select(h => h.Position).ToList();
        var oracleEbox = Oracle(seq.Sequence).Where(o => o.Name == "E-box").Select(o => o.Position).ToList();
        ebox.Should().Equal(oracleEbox,
            "the unit reports exactly the oracle's set of E-box offsets, overlaps included");
        ebox.Should().Equal(new[] { 0, 2 }, "both overlapping E-box windows are reported at offsets 0 and 2");
    }

    [Test]
    [Category("Fuzzing")]
    public void Be_TwoDifferentElementsOverlap_BothReported()
    {
        // -10 Box (TATAAT) and TATA Box (TATAAA) are distinct library elements that can
        // sit at adjacent, overlapping windows. In "GTATAATATAAA":
        //   -10 Box TATAAT at offset 1, window [1..6]
        //   TATA Box TATAAA at offset 6, window [6..11]
        // The windows share base 6, so the two DIFFERENT elements overlap; both must be
        // reported at their correct 0-based offsets (grouped per entry, §3.2/§5.2).
        var seq = new DnaSequence("GTATAATATAAA");

        var hits = MotifFinder.FindRegulatoryElements(seq).ToList();
        AssertWellFormed(hits, seq.Sequence);

        Actual(hits).Should().Equal(Oracle(seq.Sequence),
            "two different overlapping library elements are both reported at their correct offsets");
        hits.Should().Contain(h => h.Name == "-10 Box");
        hits.Should().Contain(h => h.Name == "TATA Box");
    }

    [Test]
    [Category("Fuzzing")]
    [CancelAfter(60000)]
    public void Be_SeededElementsInRandomFlanks_AllPlantedHitsRecovered()
    {
        // Plant a known element at a known offset inside random non-matching flanks and
        // assert it is recovered at the right coordinates, while the overall result still
        // equals the oracle (no spurious flank hits beyond what genuinely matches).
        var rng = new Random(172_003);
        for (int t = 0; t < 300; t++)
        {
            var (name, pattern) = Library[rng.Next(Library.Length)];
            // Concretise an E-box CANNTG by replacing N with random bases so it is valid DNA.
            string concrete = string.Concat(pattern.Select(c => c == 'N' ? Alphabet[rng.Next(4)] : c));

            string left = RandomDna(rng, rng.Next(0, 12));
            string right = RandomDna(rng, rng.Next(0, 12));
            string raw = left + concrete + right;
            var seq = new DnaSequence(raw);

            var hits = MotifFinder.FindRegulatoryElements(seq).ToList();

            AssertWellFormed(hits, seq.Sequence);
            // The planted element must appear at its planted offset (possibly among others).
            hits.Should().Contain(h => h.Name == name && h.Position == left.Length && h.Sequence == concrete,
                "the planted element is recovered at its planted 0-based offset");
            // And the full result is exactly the oracle's — no false hit, no missed hit.
            Actual(hits).Should().Equal(Oracle(seq.Sequence));
        }
    }

    #endregion
}
