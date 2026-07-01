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
/// Fuzz tests for the Analysis area — Known Motif Search (GENOMIC-MOTIFS-001), the classical exact
/// set-matching facade
/// <see cref="GenomicAnalyzer.FindKnownMotifs(DnaSequence, IEnumerable{string})"/>: for each supplied
/// (caller-curated) query motif that occurs in the subject DNA, the 0-based start positions of
/// <b>all</b> of its occurrences, returned sorted ascending and keyed by the upper-cased motif.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain motif-set / sequence inputs to the unit and
/// asserts the code NEVER fails in an undisciplined way: no NullReference/crash on an empty (empty
/// string) subject, no IndexOutOfRange when a motif is LONGER than the subject (the canonical
/// off-by-one trap — a length-5 motif against a length-3 subject must yield no hit, not throw), no
/// hang/infinite loop (O(n) suffix-tree build + O(|m|+occ) per query must always terminate), and no
/// nonsense output: never a FALSE hit (a reported position where the subject does NOT spell the
/// motif), never a MISSED overlapping occurrence, never an off-by-one coordinate, never an unsorted
/// or duplicated position list, never a key that isn't the upper-cased motif, never a motif with
/// zero occurrences present in the result. The only declared validation exception is a null motif
/// ENUMERABLE → ArgumentNullException (§3.3); an empty subject, an empty/whitespace motif element,
/// and a non-matching motif are all accepted and yield a well-defined (possibly empty) result. A raw
/// runtime exception, a hang, a false hit, a missed overlap, a wrong coordinate, or a
/// non-deterministic output is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: GENOMIC-MOTIFS-001 — Known Motif Search (exact multi-motif set matching)
/// Checklist: docs/checklists/03_FUZZING.md, row 176 (BE — "empty, no motif, overlapping").
/// Algorithm doc: docs/algorithms/Motif_Analysis/Known_Motif_Search.md
/// Distinct from row 175 GENOMIC-COMMON-001 (two-sequence longest-common-substring facade) and the
/// MotifFinder / ProteinMotifFinder rows: THIS unit is the supplied-pattern exact set-matching
/// facade on GenomicAnalyzer — the motif source is the CALLER's curated set (no fixed catalog, no
/// IUPAC degeneracy, no reverse-complement strand search — §6.2).
///
/// Fuzz strategy exercised for THIS unit (BE = Boundary Exploitation — граничні значення: 0, -1,
/// MaxInt, empty — docs/checklists/03_FUZZING.md §Description), mapped to the row's three targets:
///   • EMPTY — the degenerate subject/motif-set boundaries: empty (empty-string) subject → empty
///     dictionary, NO NullReference (§6.1); empty motif SET → empty dictionary; empty/whitespace
///     motif ELEMENTS skipped (the empty string is not a motif — deviation #1); a null motif
///     ENUMERABLE → ArgumentNullException (the ONLY declared validation exception, §3.3).
///   • NO MOTIF — a subject containing NO occurrence of the queried motif(s) (e.g. all-A subject vs
///     "C…" motif, or a motif LONGER than the subject) → that motif is OMITTED from the result, NO
///     false hit fabricated, NO IndexOutOfRange (INV-04, §6.1).
///   • OVERLAPPING — a motif occurring at OVERLAPPING offsets (the headline correctness property):
///     ALL overlapping occurrences reported with correct 0-based coordinates — e.g. "AAA" in
///     "AAAAA" → {0,1,2}, mirroring Biopython count_overlap (INV-02, §6.1, worked example §7.1).
///
/// Note on Malformed Content / Injection: the subject is a <see cref="DnaSequence"/> (uppercased and
/// validated to {A,C,G,T} at construction, so out-of-domain residues / null bytes / unicode cannot
/// reach the subject); a motif STRING is upper-cased but NOT validated — a non-ACGT motif simply
/// fails to match (returns no positions) rather than raising (§6.2). This is therefore a pure
/// boundary (BE) row over the subject shape (empty / non-matching / overlapping) and the motif set
/// (empty / whitespace / over-long / non-ACGT / duplicate / mixed-case), exactly as the row says.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Known_Motif_Search.md §2.4, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
///   INV-01 every reported position p for motif m satisfies T[p..p+|m|-1] == m (0-based exact);
///   INV-02 ALL occurrences reported, including overlapping ones (e.g. AAA in AAAAA → {0,1,2});
///   INV-03 each motif's positions are sorted strictly ascending and DISTINCT (a set);
///   INV-04 a motif with zero occurrences is OMITTED from the result;
///   INV-05 result keys are the UPPER-CASED motif strings; duplicate upper-cased keys collapse;
///   empty/whitespace motif elements are skipped; empty subject → empty result; null motifs → throw.
///   GenomicAnalyzer.FindKnownMotifs(DnaSequence, IEnumerable&lt;string&gt;) → Dictionary&lt;string, IReadOnlyList&lt;int&gt;&gt;
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class GenomicMotifsFuzzTests
{
    private static readonly char[] Alphabet = { 'A', 'C', 'G', 'T' };

    #region Helpers

    /// <summary>A random ACGT string of the given length (length 0 ⇒ empty string).</summary>
    private static string RandomDna(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(Alphabet[rng.Next(Alphabet.Length)]);
        return sb.ToString();
    }

    /// <summary>
    /// Independent O(n·|m|) naive oracle for the exact (overlap-aware) start positions of a single
    /// motif in a subject — the textbook exact-matching set { i : T[i..i+|m|-1] == P } [Gusfield].
    /// Built from the spec (a sliding-window scan over EVERY start, including overlapping ones), not
    /// from the unit. An empty motif yields the empty set (the empty string is not a motif).
    /// </summary>
    private static List<int> OracleOccurrences(string subject, string motif)
    {
        var positions = new List<int>();
        if (string.IsNullOrEmpty(motif) || motif.Length > subject.Length)
            return positions; // over-long / empty motif ⇒ no occurrence (no IndexOutOfRange)
        for (int i = 0; i + motif.Length <= subject.Length; i++)
        {
            if (subject.AsSpan(i, motif.Length).SequenceEqual(motif.AsSpan()))
                positions.Add(i); // overlap-aware: every start that matches, not just disjoint ones
        }
        return positions; // already ascending and distinct by construction
    }

    /// <summary>
    /// Asserts a FindKnownMotifs result is WELL-FORMED per the documented contract, cross-checked
    /// against the independent naive oracle over the SAME (upper-cased) subject:
    ///   • each KEY is a non-empty UPPER-CASED string (INV-05) and equals its own ToUpperInvariant;
    ///   • each position list is non-empty (zero-occurrence motifs are OMITTED — INV-04), sorted
    ///     strictly ascending and DISTINCT (INV-03), and every position is an in-bounds 0-based
    ///     start (no IndexOutOfRange);
    ///   • the SUBSTRING at each reported position EXACTLY equals the key motif (INV-01 — no false
    ///     hit, no off-by-one);
    ///   • the reported position set is EXACTLY the oracle's overlap-aware occurrence set for that
    ///     motif (INV-02 — all overlapping occurrences, none missed, none fabricated).
    /// </summary>
    private static void AssertWellFormed(
        Dictionary<string, IReadOnlyList<int>> result, string subject)
    {
        result.Should().NotBeNull();
        foreach (var (key, positions) in result)
        {
            key.Should().NotBeNullOrWhiteSpace("empty/whitespace motifs are skipped, never keyed");
            key.Should().Be(key.ToUpperInvariant(), "result keys are the upper-cased motif (INV-05)");

            positions.Should().NotBeNull();
            positions.Should().NotBeEmpty("a motif with zero occurrences is omitted (INV-04)");
            positions.Should().BeInAscendingOrder("positions are sorted ascending (INV-03)");
            positions.Should().OnlyHaveUniqueItems("positions form a distinct set (INV-03)");

            foreach (int p in positions)
            {
                p.Should().BeInRange(0, subject.Length - key.Length,
                    "every reported position is an in-bounds 0-based start (INV-01)");
                subject.Substring(p, key.Length).Should().Be(key,
                    "the subject spells the motif at the reported position — no false hit / off-by-one (INV-01)");
            }

            // INV-02: exactly the overlap-aware occurrence set — none missed, none fabricated.
            positions.Should().Equal(OracleOccurrences(subject, key),
                "all (overlapping) occurrences reported, exactly matching the naive oracle (INV-02)");
        }
    }

    #endregion

    #region GENOMIC-MOTIFS-001 — Known Motif Search (FindKnownMotifs) — BE (Boundary Exploitation)

    // ── EMPTY — degenerate subject / motif-set / motif-element boundaries ──────────────────────

    [Test]
    public void FindKnownMotifs_EmptySubject_ReturnsEmpty_NoNullReference()
    {
        var seq = new DnaSequence(string.Empty);

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "ACGT", "A", "GAATTC" });

        result.Should().BeEmpty("no motif occurs in an empty subject (§6.1) — and no NullReference");
        AssertWellFormed(result, seq.Sequence);
    }

    [Test]
    public void FindKnownMotifs_EmptyMotifSet_ReturnsEmpty()
    {
        var seq = new DnaSequence("ACGTACGTACGT");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, Array.Empty<string>());

        result.Should().BeEmpty("an empty motif set has nothing to search (§6.1)");
    }

    [Test]
    public void FindKnownMotifs_EmptyAndWhitespaceMotifs_AreSkipped()
    {
        var seq = new DnaSequence("ACGTACGT");

        // The empty string is not a motif; whitespace-only motifs are skipped (deviation #1).
        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "", "   ", "\t", "ACGT" });

        result.Keys.Should().NotContain("", "the empty string is not a motif (deviation #1)");
        result.Keys.Should().NotContain("   ");
        result.Should().ContainKey("ACGT", "the one real motif is searched");
        AssertWellFormed(result, seq.Sequence);
    }

    [Test]
    public void FindKnownMotifs_AllEmptyOrWhitespaceMotifs_ReturnsEmpty()
    {
        var seq = new DnaSequence("ACGTACGT");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "", " ", "\t", "\n", "  " });

        result.Should().BeEmpty("every motif element is empty/whitespace and skipped (deviation #1)");
    }

    [Test]
    public void FindKnownMotifs_NullMotifEnumerable_Throws()
    {
        var seq = new DnaSequence("ACGTACGT");

        Action act = () => GenomicAnalyzer.FindKnownMotifs(seq, null!);

        act.Should().Throw<ArgumentNullException>("null motifs is the only declared validation exception (§3.3)");
    }

    // ── NO MOTIF — subject contains no occurrence; motif longer than subject ───────────────────

    [Test]
    public void FindKnownMotifs_MotifAbsent_IsOmitted_NoFalseHit()
    {
        var seq = new DnaSequence("AAAAAAAAAA"); // all A

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "C", "GG", "ACGT", "T" });

        result.Should().BeEmpty("none of these motifs occurs in an all-A subject (INV-04) — no false hit");
    }

    [Test]
    public void FindKnownMotifs_MotifLongerThanSubject_NoIndexOutOfRange()
    {
        var seq = new DnaSequence("ACG"); // length 3

        // A motif longer than the subject is the canonical off-by-one / IndexOutOfRange trap.
        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "ACGTACGT", "ACGT", "ACG" });

        result.Should().NotContainKey("ACGTACGT", "over-long motif cannot occur — no IndexOutOfRange");
        result.Should().NotContainKey("ACGT");
        result.Should().ContainKey("ACG", "the full-length motif occurs once at position 0");
        result["ACG"].Should().Equal(0);
        AssertWellFormed(result, seq.Sequence);
    }

    [Test]
    public void FindKnownMotifs_NonAcgtMotif_DoesNotMatch_NoThrow()
    {
        var seq = new DnaSequence("ACGTACGT");

        // Motif strings are not validated; a non-ACGT motif simply fails to match (§6.2).
        Action act = () =>
        {
            var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "ACNT", "ACGT" });
            result.Should().NotContainKey("ACNT", "a non-ACGT motif cannot match the subject (§6.2)");
            result.Should().ContainKey("ACGT");
            AssertWellFormed(result, seq.Sequence);
        };

        act.Should().NotThrow("a non-ACGT motif fails to match rather than raising (§6.2)");
    }

    // ── OVERLAPPING — the headline correctness property ────────────────────────────────────────

    [Test]
    public void FindKnownMotifs_OverlappingOccurrences_AllReported()
    {
        var seq = new DnaSequence("AAAAA"); // worked example §7.1

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "AAA" });

        result.Should().ContainKey("AAA");
        result["AAA"].Should().Equal(new[] { 0, 1, 2 },
            "AAA in AAAAA has three overlapping occurrences {0,1,2} (INV-02, Biopython count_overlap)");
        AssertWellFormed(result, seq.Sequence);
    }

    [Test]
    public void FindKnownMotifs_PositiveSanity_KnownMotifsAtKnownOffsets()
    {
        // EcoRI sites at 0 and 9; a periodic AT motif overlapping; a non-matching motif omitted.
        var seq = new DnaSequence("GAATTCAAAGAATTC");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "GAATTC", "AAA", "CCCC" });

        result["GAATTC"].Should().Equal(new[] { 0, 9 }, "the two EcoRI sites (worked example §7.1)");
        result["AAA"].Should().Equal(new[] { 6 }, "the single AAA run at offset 6 (positions 6 only — len-3 run)");
        result.Should().NotContainKey("CCCC", "a non-matching motif is omitted (INV-04)");
        AssertWellFormed(result, seq.Sequence);
    }

    [Test]
    public void FindKnownMotifs_OverlappingPeriodicMotif_AllStartsReported()
    {
        // "ATAT" occurs at overlapping offsets 0 and 2 in "ATATAT".
        var seq = new DnaSequence("ATATAT");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "ATAT", "AT" });

        result["ATAT"].Should().Equal(new[] { 0, 2 }, "overlapping periodic occurrences (INV-02)");
        result["AT"].Should().Equal(new[] { 0, 2, 4 }, "every AT start (INV-02)");
        AssertWellFormed(result, seq.Sequence);
    }

    // ── Case-insensitivity & duplicate-key collapse ───────────────────────────────────────────

    [Test]
    public void FindKnownMotifs_LowerCaseAndDuplicateMotifs_CollapseToUpperCasedKey()
    {
        var seq = new DnaSequence("ACGTACGT");

        var result = GenomicAnalyzer.FindKnownMotifs(seq, new[] { "acgt", "ACGT", "AcGt" });

        result.Keys.Should().ContainSingle().Which.Should().Be("ACGT",
            "all three normalize to the same upper-cased key (INV-05, duplicate collapse)");
        result["ACGT"].Should().Equal(0, 4);
        AssertWellFormed(result, seq.Sequence);
    }

    // ── Randomized fuzz sweep: every result well-formed against the independent oracle ─────────

    [Test]
    [CancelAfter(30_000)]
    public void FindKnownMotifs_RandomSubjectsAndMotifs_AlwaysWellFormed()
    {
        var rng = new Random(176_001); // locally seeded for reproducibility

        for (int iter = 0; iter < 400; iter++)
        {
            // Subject length 0..40 — includes the empty boundary.
            string subject = RandomDna(rng, rng.Next(0, 41));
            var seq = new DnaSequence(subject);

            // 1..5 motifs of length 0..6 (length 0 ⇒ empty, skipped; >subject ⇒ no hit).
            int motifCount = rng.Next(1, 6);
            var motifs = new List<string>(motifCount);
            for (int k = 0; k < motifCount; k++)
                motifs.Add(RandomDna(rng, rng.Next(0, 7)));

            Dictionary<string, IReadOnlyList<int>> result = null!;
            Action act = () => result = GenomicAnalyzer.FindKnownMotifs(seq, motifs);
            act.Should().NotThrow("no input shape may crash the unit (subject='{0}', motifs=[{1}])",
                subject, string.Join(",", motifs));

            AssertWellFormed(result, seq.Sequence);

            // Cross-check: every NON-empty motif that the oracle finds MUST be present with the
            // exact occurrence set; every motif the oracle finds nothing for MUST be absent.
            foreach (var m in motifs.Where(s => !string.IsNullOrWhiteSpace(s))
                                    .Select(s => s.ToUpperInvariant()).Distinct())
            {
                var expected = OracleOccurrences(subject, m);
                if (expected.Count == 0)
                    result.Should().NotContainKey(m, "zero-occurrence motif omitted (INV-04)");
                else
                    result[m].Should().Equal(expected, "exact overlap-aware occurrence set (INV-01/02)");
            }
        }
    }

    [Test]
    public void FindKnownMotifs_Deterministic_AcrossRepeatedCalls()
    {
        var seq = new DnaSequence("AAAAGAATTCAAAA");
        var motifs = new[] { "AAAA", "GAATTC", "AA" };

        var first = GenomicAnalyzer.FindKnownMotifs(seq, motifs);
        var second = GenomicAnalyzer.FindKnownMotifs(seq, motifs);

        second.Keys.Should().BeEquivalentTo(first.Keys);
        foreach (var key in first.Keys)
            second[key].Should().Equal(first[key], "the result is deterministic across calls");
    }

    #endregion
}
