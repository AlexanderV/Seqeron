using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;
// Disambiguate: a top-level Seqeron.Genomics.Analysis.MotifMatch also exists; this unit
// asserts against the ProteinMotifFinder.MotifMatch record returned by FindCommonMotifs.
using MotifMatch = Seqeron.Genomics.Analysis.ProteinMotifFinder.MotifMatch;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the ProteinMotif area — COMMON MOTIF FINDING
/// (PROTMOTIF-COMMON-001): scanning a single protein sequence against the curated
/// in-source <c>CommonMotifs</c> PROSITE-style library and reporting every occurrence
/// of every library pattern. The single public entry point under test is
/// <see cref="ProteinMotifFinder.FindCommonMotifs(string)"/>, which iterates
/// <c>CommonMotifs.Values</c> and delegates each entry's regex to
/// <c>FindMotifByPattern</c>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Scope vs sibling unit PROTMOTIF-FIND-001 (row 82)
/// ───────────────────────────────────────────────────────────────────────────
/// Cross-referencing docs/algorithms/ProteinMotif/Common_Motif_Finding.md (Test Unit ID
/// PROTMOTIF-COMMON-001) and the source confirms there is NO multi-sequence "common
/// across a collection" method — the ONLY common-motif entry point is the single-sequence
/// <c>FindCommonMotifs(string)</c> scanning the fixed library. Row 82 (PROTMOTIF-FIND-001)
/// already exercised the GENERAL motif-finding contract via the same method on the empty /
/// non-amino-acid / extremely-short / homopolymer (all-same-char) facets. This file is
/// COMPLEMENTARY: it owns the row-164 BE facets row 82 did not assert —
///   • "single sequence": the documented single-input contract — one sequence in, its
///     library motifs are the trivially-"common" set; the per-pattern grouping / ordering
///     contract (§5.2) and the "common across a COLLECTION" set-intersection semantics that
///     any multi-sequence analysis would build on (computed here at the test level, since no
///     such method exists), including the single-element collection degenerate case.
///   • "no common motif": a sequence (and a collection) sharing NO library motif → the EMPTY
///     common set, no FALSE common motif, no crash, no DivideByZero in any frequency fraction.
///   • "identical inputs": all sequences identical → every library motif of that sequence is
///     "common", the result is DETERMINISTIC (INV-04) and is NOT double-counted.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption, no nonsense
/// output, and no *unhandled* runtime exception (NullReference / IndexOutOfRange /
/// DivideByZero / ArgumentOutOfRange). Every input must resolve to EITHER a well-defined,
/// theory-correct result OR a *documented, intentional* outcome. For a library scan that
/// uppercases the input and runs every curated PROSITE regex, the headline hazards are:
///   • a NullReferenceException on a null sequence (or, for a "common across sequences"
///     analysis, a null/empty collection or a single-element collection);
///   • a FALSE common motif — reporting a motif as present/common when no window of the
///     sequence actually satisfies its pattern (a correctness bug);
///   • non-deterministic output across identical re-scans (would break INV-04);
///   • a DivideByZero in any conservation/frequency fraction over an empty collection;
///   • a mis-reported match whose [Start..End] span does NOT equal the substring it claims
///     (a coordinate bug — INV-01/INV-02).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROTMOTIF-COMMON-001 — common motif finding (curated PROSITE library scan)
/// Checklist: docs/checklists/03_FUZZING.md, row 164.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the count / composition corners that could crash, raise
///     a false common motif, or behave non-deterministically:
///       – "single sequence": one sequence in → its library motifs are the (trivially)
///         common set; a single-element collection's common set equals that one sequence's
///         motif set; no crash on the degenerate count-of-one (Common_Motif_Finding.md §4.1).
///       – "no common motif": a sequence with NO library hit → empty result; a collection of
///         sequences sharing NO motif → empty intersection, no FALSE common motif, no
///         DivideByZero (§6.1 "sequence with no motif → empty result").
///       – "identical inputs": identical sequences → every motif of that sequence is common,
///         deterministic (INV-04), not double-counted (§2.4 INV-04, §5.2).
/// — docs/checklists/03_FUZZING.md §Description (strategy code BE);
///   targets: "single sequence, no common motif, identical inputs".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The common-motif-finding contract under test (Common_Motif_Finding.md)
/// ───────────────────────────────────────────────────────────────────────────
/// Given a protein sequence S, FindCommonMotifs uppercases S and, for every pattern in the
/// fixed CommonMotifs library, reports every occurrence as a MotifMatch { Start (incl. 0-based),
/// End (incl. 0-based), Sequence (the uppercased matched substring), MotifName (PROSITE entry
/// name), Pattern (PROSITE accession), Score, EValue }. Matches are grouped by pattern in
/// dictionary-iteration order; within a pattern they are in increasing start order (§5.2).
/// null / empty input yields an empty result (§3.3, §6.1).
///   — Common_Motif_Finding.md §3.1–§3.3, §4.1–§4.2, §5.1–§5.2.
///
/// Theory-correct invariants asserted (Common_Motif_Finding.md §2.4):
///   • INV-01 — each reported match's substring equals the residues at its coordinates.
///   • INV-02 — 0 ≤ Start ≤ End < length for every match.
///   • INV-03 — overlapping occurrences are all reported.
///   • INV-04 — the scan is deterministic.
/// Plus the "common across a collection" set semantics this unit owns:
///   • [common-single] — the common set over a one-sequence collection equals that sequence's
///     own motif-accession set.
///   • [common-none]   — the common set over sequences sharing no accession is EMPTY (no false
///     common motif).
///   • [common-identical] — the common set over identical sequences equals each sequence's
///     motif-accession set (every motif is common), with no double-counting.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Complexity / hang-safety
/// ───────────────────────────────────────────────────────────────────────────
/// The scan is O(p · n): a constant number p of bounded-width PROSITE regex walks over n
/// residues (Common_Motif_Finding.md §4.3). The dense / long-input targets are kept modest and
/// [CancelAfter]-guarded so a regression that turned a library scan into a hang or super-linear
/// blow-up would FAIL rather than wedge the suite.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ProteinCommonMotifFuzzTests
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
    /// satisfy against the original (case-insensitive) sequence (Common_Motif_Finding.md §2.4):
    /// the match is a CONTIGUOUS in-bounds subsequence (INV-02) whose claimed coordinates actually
    /// reproduce the reported substring (INV-01: S[Start..End] == Sequence, uppercased), and whose
    /// Score / EValue are finite. This is the headline "no coordinate bug, no run-off-the-end,
    /// no NaN" property.
    /// </summary>
    private static void AssertWellFormedMatch(MotifMatch match, string originalSequence)
    {
        string upper = originalSequence.ToUpperInvariant();
        int n = upper.Length;

        // INV-02 — in-bounds, non-empty, contiguous span.
        match.Start.Should().BeInRange(0, n - 1, "a match Start is a valid 0-based residue index");
        match.End.Should().BeInRange(match.Start, n - 1, "a match End is in-bounds and not before its Start");

        // INV-01 — the span length equals the reported substring length, which is uppercased and
        // equals the actual substring of the (uppercased) input at [Start..End].
        (match.End - match.Start + 1).Should().Be(match.Sequence.Length,
            "End − Start + 1 equals the matched substring length");
        match.Sequence.Should().Be(match.Sequence.ToUpperInvariant(), "the matched substring is uppercased");
        upper.Substring(match.Start, match.Sequence.Length).Should().Be(match.Sequence,
            "INV-01: the reported substring is exactly S[Start..End] of the uppercased input (no coordinate bug)");

        // Finite score / E-value (no NaN / ±∞ → guards against a DivideByZero leaking into scoring).
        double.IsNaN(match.Score).Should().BeFalse("a motif Score must never be NaN");
        double.IsInfinity(match.Score).Should().BeFalse("a motif Score must never be infinite");
        double.IsNaN(match.EValue).Should().BeFalse("a motif EValue must never be NaN");
        double.IsInfinity(match.EValue).Should().BeFalse("a motif EValue must never be infinite");
    }

    /// <summary>
    /// The set of library motif ACCESSIONS that occur in <paramref name="sequence"/>, computed by
    /// <see cref="ProteinMotifFinder.FindCommonMotifs(string)"/>. This is the per-sequence "which
    /// common motifs are present" projection on which a "common across sequences" analysis is built.
    /// </summary>
    private static HashSet<string> PresentMotifAccessions(string sequence) =>
        FindCommonMotifs(sequence).Select(m => m.Pattern).ToHashSet();

    /// <summary>
    /// The motif accessions COMMON to every sequence in <paramref name="sequences"/> — the
    /// intersection of each sequence's present-motif set. By definition an empty collection has no
    /// common motif; a single-element collection's common set equals that sequence's own set. There
    /// is NO division anywhere (the set-intersection model is DivideByZero-free by construction).
    /// </summary>
    private static HashSet<string> CommonMotifAccessions(IReadOnlyList<string> sequences)
    {
        if (sequences.Count == 0)
            return new HashSet<string>();

        HashSet<string> common = PresentMotifAccessions(sequences[0]);
        for (int i = 1; i < sequences.Count; i++)
            common.IntersectWith(PresentMotifAccessions(sequences[i]));
        return common;
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PROTMOTIF-COMMON-001 — common motif finding : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PROTMOTIF-COMMON-001 — common motif finding

    #region BE — Single sequence: documented single-input contract, no crash on count-of-one

    /// <summary>
    /// Target "single sequence": the ONLY common-motif entry point takes a single sequence and
    /// reports its library motifs (Common_Motif_Finding.md §4.1) — those motifs are, trivially,
    /// the set "common" to a collection of just that one sequence. A single sequence (including
    /// the degenerate empty / null sequence) must never crash, must surface every present library
    /// motif as a well-formed match, and the "common across a one-element collection" set must
    /// equal that single sequence's own present-motif set ([common-single]).
    /// </summary>
    [Test]
    public void Common_SingleSequence_PresentMotifsAreTriviallyCommon()
    {
        // A sequence carrying several distinct library motifs: an RGD (PS00016), an
        // N-glycosylation site NFTA (PS00001) and a P-loop ATP/GTP-A "A....GKS" (PS00017).
        const string seq = "AARGDKKNFTAKKAVVVVGKSKK";

        var present = PresentMotifAccessions(seq);
        present.Should().Contain("PS00016", "the RGD cell-attachment motif occurs in the sequence");
        present.Should().Contain("PS00001", "the N-glycosylation site NFTA occurs in the sequence");

        // Every emitted match is well-formed and its claimed accession genuinely occurs.
        var hits = FindCommonMotifs(seq).ToList();
        hits.Should().NotBeEmpty("the constructed sequence carries multiple library motifs");
        foreach (var hit in hits)
        {
            AssertWellFormedMatch(hit, seq);
            present.Should().Contain(hit.Pattern, "every reported hit's accession is a present motif");
        }

        // [common-single]: the common set over the one-element collection equals that sequence's set.
        CommonMotifAccessions(new[] { seq })
            .Should().BeEquivalentTo(present,
                "[common-single]: a single-sequence collection's common motifs are exactly that sequence's motifs");

        // Degenerate single sequences (empty / null / whitespace-only): no crash, no motif, and the
        // one-element common set is empty — never a NullReference on the count-of-one boundary.
        foreach (string? degenerate in new[] { "", " ", null })
        {
            ((Func<List<MotifMatch>>)(() => FindCommonMotifs(degenerate!).ToList()))
                .Should().NotThrow($"a degenerate single sequence ('{degenerate ?? "null"}') must not crash")
                .Subject.Should().BeEmpty("a degenerate single sequence carries no library motif");

            CommonMotifAccessions(new[] { degenerate ?? string.Empty })
                .Should().BeEmpty("a one-element collection of a degenerate sequence has no common motif");
        }
    }

    #endregion

    #region BE — No common motif: empty result / empty intersection, no false common motif

    /// <summary>
    /// Target "no common motif": a sequence containing NO library motif yields an EMPTY result
    /// (Common_Motif_Finding.md §6.1 "sequence with no motif → empty result"), and a COLLECTION of
    /// sequences that share no library motif yields an EMPTY common set — never a FALSE common motif
    /// and never a DivideByZero in any frequency fraction over the (set-based) analysis. We use a
    /// homopolymer "WWWW…" (W never appears at any core position of the curated patterns) as the
    /// no-motif sequence, and a pair {RGD-only, P-loop-only} whose present-motif sets are disjoint.
    /// </summary>
    [Test]
    public void Common_NoCommonMotif_EmptyResultNoFalsePositive()
    {
        // (a) A sequence with no library motif → empty result, no false common motif.
        const string noMotif = "WWWWWWWWWWWWWWWW";
        var noHits = ((Func<List<MotifMatch>>)(() => FindCommonMotifs(noMotif).ToList()))
            .Should().NotThrow("a no-motif sequence must not crash").Subject;
        noHits.Should().BeEmpty("[common-none]: a sequence with no library motif yields no match");
        PresentMotifAccessions(noMotif).Should().BeEmpty("no library accession is present in a no-motif sequence");

        // (b) Two sequences whose present-motif sets are DISJOINT → empty common set.
        const string rgdOnly = "WWWRGDWWWWWW";   // carries RGD (PS00016) but not the P-loop
        const string ploopOnly = "WWWAVVVVGKSWWW"; // carries P-loop "[AG]x(4)GK[ST]" but not RGD

        var rgdSet = PresentMotifAccessions(rgdOnly);
        var ploopSet = PresentMotifAccessions(ploopOnly);
        rgdSet.Should().Contain("PS00016", "the first sequence carries the RGD motif");
        ploopSet.Should().Contain("PS00017", "the second sequence carries the P-loop motif");

        // Construct the disjointness premise robustly: the two sequences must genuinely share no
        // accession for the "no common motif" target to be meaningful.
        rgdSet.Overlaps(ploopSet).Should().BeFalse("the two constructed sequences share no library motif");

        CommonMotifAccessions(new[] { rgdOnly, ploopOnly })
            .Should().BeEmpty("[common-none]: sequences sharing no motif have an empty common set");

        // (c) A collection that includes a no-motif sequence can never have a common motif.
        CommonMotifAccessions(new[] { rgdOnly, noMotif })
            .Should().BeEmpty("[common-none]: any collection containing a no-motif sequence has no common motif");

        // (d) The empty collection is the DivideByZero corner of any frequency-based analysis: it
        //     must yield the empty common set, not throw.
        ((Func<HashSet<string>>)(() => CommonMotifAccessions(Array.Empty<string>())))
            .Should().NotThrow("the empty collection must not divide by zero")
            .Subject.Should().BeEmpty("an empty collection has no common motif");
    }

    #endregion

    #region BE — Identical inputs: every motif common, deterministic, no double-count

    /// <summary>
    /// Target "identical inputs": when all sequences are IDENTICAL, every library motif present in
    /// that sequence is "common" to the collection, the per-sequence scan is DETERMINISTIC (INV-04)
    /// and the common set equals the single sequence's motif set without double-counting
    /// ([common-identical]). We re-scan the same sequence and require byte-identical hit lists
    /// (INV-04), then form a collection of N copies and require its common set to equal the
    /// single-scan set exactly — adding more identical copies must not grow OR shrink it.
    /// </summary>
    [Test]
    public void Common_IdenticalInputs_AllMotifsCommonDeterministicNoDoubleCount()
    {
        const string seq = "AARGDKKNFTAKKAVVVVGKSKKSPRKK";

        // INV-04 — determinism: identical re-scans yield identical hit lists.
        var first = FindCommonMotifs(seq).ToList();
        var second = FindCommonMotifs(seq).ToList();
        first.Should().NotBeEmpty("the constructed sequence carries library motifs");
        second.Select(h => (h.MotifName, h.Pattern, h.Start, h.End, h.Sequence))
            .Should().Equal(first.Select(h => (h.MotifName, h.Pattern, h.Start, h.End, h.Sequence)),
                "INV-04: FindCommonMotifs is deterministic for a fixed input");

        var single = PresentMotifAccessions(seq);

        // [common-identical]: a collection of N identical copies has the SAME common set as one copy —
        // every motif is common, none is double-counted or lost as copies are added.
        foreach (int copies in new[] { 1, 2, 5, 20 })
        {
            var collection = Enumerable.Repeat(seq, copies).ToArray();
            CommonMotifAccessions(collection)
                .Should().BeEquivalentTo(single,
                    $"[common-identical]: {copies} identical copies share exactly the one sequence's motif set");
        }

        // No double-count within a single scan: the per-(accession,Start) hit is unique (a given
        // library pattern reports each occurrence once, not once per identical input or per re-scan).
        first.Select(h => (h.Pattern, h.Start)).Should().OnlyHaveUniqueItems(
            "each library motif occurrence is reported exactly once (no double-count)");
    }

    #endregion

    #region BE — Random inputs: never crash/hang, well-formed, deterministic, set-laws hold

    /// <summary>
    /// Positive sanity over RANDOM single sequences and collections: across fixed seeds/lengths the
    /// scan must never crash, hang, or emit a malformed match, every MotifMatch must satisfy the full
    /// contract (INV-01/INV-02, finite score), the scan must be deterministic (INV-04), and the
    /// "common across a collection" set must obey its laws — common(self) = present(self), common over
    /// a set containing a no-motif member is empty, and intersection is order-independent. This pins
    /// the targets on arbitrary sequences, not just hand-built motifs.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Common_RandomInputs_WellFormedDeterministicSetLaws(CancellationToken token)
    {
        foreach (int seed in new[] { 11, 53, 211, 2026 })
        {
            foreach (int len in new[] { 1, 5, 21, 80, 300 })
            {
                string seq = RandomProtein(len, seed);

                var hits = ((Func<List<MotifMatch>>)(() => FindCommonMotifs(seq).ToList()))
                    .Should().NotThrow($"random protein must not crash (seed {seed}, len {len})").Subject;
                token.ThrowIfCancellationRequested();

                foreach (var hit in hits)
                    AssertWellFormedMatch(hit, seq);

                // INV-04 — determinism.
                FindCommonMotifs(seq).Select(h => (h.MotifName, h.Pattern, h.Start, h.End, h.Sequence))
                    .Should().Equal(hits.Select(h => (h.MotifName, h.Pattern, h.Start, h.End, h.Sequence)),
                        "INV-04: FindCommonMotifs is deterministic for a fixed input");

                // [common-single]: the common set over a one-element collection equals present(seq).
                var present = PresentMotifAccessions(seq);
                CommonMotifAccessions(new[] { seq }).Should().BeEquivalentTo(present,
                    "[common-single]: a one-element collection's common set is the sequence's own motif set");

                // [common-none]: adding a guaranteed no-motif member empties the common set.
                CommonMotifAccessions(new[] { seq, "WWWWWWWWWWWWWWWW" }).Should().BeEmpty(
                    "[common-none]: any collection containing a no-motif sequence has no common motif");
            }
        }

        // Intersection is order-independent (a "common" set cannot depend on collection order).
        string a = RandomProtein(120, 7);
        string b = RandomProtein(120, 7); // identical to a (same seed) → guarantees a non-trivial common set is possible
        string c = "WWWRGDWWWRGDWWW";
        CommonMotifAccessions(new[] { a, b, c })
            .Should().BeEquivalentTo(CommonMotifAccessions(new[] { c, b, a }),
                "the common set is independent of collection ordering");
    }

    #endregion

    #endregion
}
