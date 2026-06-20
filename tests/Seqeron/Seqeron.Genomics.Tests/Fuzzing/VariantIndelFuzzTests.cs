using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Variants area — VARIANT-INDEL-001 (Indel Detection).
/// The unit under test is the alignment-based indel detector:
/// <see cref="VariantCaller.FindInsertions"/>, <see cref="VariantCaller.FindDeletions"/>
/// and their union <see cref="VariantCaller.FindIndels"/>, which filter the
/// deterministic column-scan caller
/// <see cref="VariantCaller.CallVariantsFromAlignment"/> (fed by
/// <see cref="SequenceAligner.GlobalAlign"/> through <c>CallVariants</c>) to one
/// indel class; implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no nonsense output, and no *unhandled* runtime exception
/// (IndexOutOfRange / NullReference / ArgumentOutOfRange / negative-length
/// substring). Every input must resolve to EITHER a well-defined, theory-correct
/// value OR a *documented, intentional* outcome (here, an ArgumentNullException
/// for null inputs; an empty result for empty inputs). —
/// docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: VARIANT-INDEL-001 — Indel Detection (Variants)
/// Checklist: docs/checklists/03_FUZZING.md, row 188.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// MAPPING of the generic checklist BE targets ("length-0 indel, indel at edge,
/// empty") onto THIS unit's documented contract
/// (docs/algorithms/Variants/Indel_Detection.md):
///   • "length-0 indel" → a degenerate "indel" with NO actual length change. By
///       the column rule an indel column requires a gap on exactly ONE side
///       (§2.2). A gap–gap column, a substitution column (mismatch, length
///       preserving), and an all-match alignment all carry NO length change ⇒
///       NO indel must be emitted — never a zero-length / phantom indel record
///       (INV-01, §6.1). Each emitted indel is exactly one base (one gap column),
///       so a "length-0" record is impossible by construction; the tests pin it.
///   • "indel at edge" → an insertion / deletion at the very FIRST or very LAST
///       reference column. The single-base allele slice at the boundary must not
///       throw IndexOutOfRange / produce a negative-length substring, and the
///       reported 0-based Position must be the correct edge coordinate ∈
///       [0, reference.Length] (INV-06, §3.2).
///   • "empty" → empty reference and/or empty query. CallVariants ⇒ empty aligned
///       strings ⇒ documented EMPTY result (no crash); null ⇒ documented
///       ArgumentNullException (§3.3, §6.1).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (docs/algorithms/Variants/Indel_Detection.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • Column rule: ref-gap opposite a query base ⇒ Insertion (the query base);
///       query-gap opposite a reference base ⇒ Deletion (the reference base);
///       mismatch ⇒ SNP (NOT an indel); shared match ⇒ no variant.   (§2.2, §4.1)
///   • Directional length (VCF): insertion ⇒ ALT longer than REF, encoded
///       in-memory as ReferenceAllele == "-" and a one-base AlternateAllele;
///       deletion ⇒ REF longer than ALT, encoded as AlternateAllele == "-" and
///       a one-base ReferenceAllele.                                  (§2.2, INV-03/04)
///   • A multi-base indel is reported as k consecutive single-base columns
///       (one event per base).                                        (INV-05, §6.1)
///   • FindInsertions ⇒ only Type==Insertion; FindDeletions ⇒ only
///       Type==Deletion; FindIndels ⇒ the union, no SNPs.             (INV-02)
///   • Every reported indel Position ∈ [0, reference.Length].         (INV-06)
///   • Identical sequences ⇒ no insertions and no deletions.          (INV-01, §6.1)
///   • Null reference/query ⇒ ArgumentNullException; empty sequences ⇒ empty
///       result.                                                       (§3.3, §6.1)
///   • Indels are NOT left-aligned / parsimony-normalized (ASM-02): exact
///       Position is asserted only where the global alignment is provably unique
///       (no internal repeat permitting an equal-score left shift); counts /
///       types / alleles are asserted generally.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class VariantIndelFuzzTests
{
    private const string GapAllele = "-";

    // ── Well-formed-indel assertion helper ───────────────────────────────────
    // Pins the documented per-indel contract on EVERY emitted indel, no matter how
    // degenerate the input. This is what stops a fuzz test from rubber-stamping
    // nonsense output green — a phantom length-0 indel, an out-of-range position,
    // or a type/allele inconsistent with the gap side would all fail here:
    //   • Type is Insertion XOR Deletion (never SNP) (INV-02).
    //   • Exactly ONE side is the "-" gap sentinel and the OTHER side is a single,
    //     non-gap base — i.e. a real length change, never a length-0 record
    //     (insertion ⇒ ALT longer than REF; deletion ⇒ REF longer than ALT) (INV-03/04).
    //   • Position is a finite, non-negative reference coordinate ∈ [0, refLen]
    //     (INV-06); QueryPosition ≥ 0.
    private static void AssertWellFormedIndel(Variant v, int refLen)
    {
        v.Type.Should().BeOneOf(new[] { VariantType.Insertion, VariantType.Deletion },
            "this unit reports only indels (INV-02)");

        v.Position.Should().BeGreaterThanOrEqualTo(0, "Position is a 0-based reference coordinate (INV-06)");
        v.Position.Should().BeLessThanOrEqualTo(refLen, "Position ∈ [0, reference.Length] (INV-06)");
        v.QueryPosition.Should().BeGreaterThanOrEqualTo(0, "QueryPosition is a 0-based query coordinate");

        if (v.Type == VariantType.Insertion)
        {
            // Insertion ⇒ ALT longer than REF: ref-side gap, one inserted base (§2.2, INV-03).
            v.ReferenceAllele.Should().Be(GapAllele, "an insertion has a ref-side gap (§2.2)");
            v.AlternateAllele.Should().HaveLength(1).And.NotBe(GapAllele,
                "an insertion carries exactly one inserted base — never a length-0 record (INV-03)");
            v.AlternateAllele.Length.Should().BeGreaterThan(v.ReferenceAllele.TrimEnd('-').Length,
                "insertion ⇒ ALT longer than REF (real length change, not length-0)");
        }
        else
        {
            // Deletion ⇒ REF longer than ALT: query-side gap, one deleted base (§2.2, INV-04).
            v.AlternateAllele.Should().Be(GapAllele, "a deletion has a query-side gap (§2.2)");
            v.ReferenceAllele.Should().HaveLength(1).And.NotBe(GapAllele,
                "a deletion carries exactly one deleted base — never a length-0 record (INV-04)");
            v.ReferenceAllele.Length.Should().BeGreaterThan(v.AlternateAllele.TrimEnd('-').Length,
                "deletion ⇒ REF longer than ALT (real length change, not length-0)");
        }
    }

    private static void AssertAllWellFormedIndels(IReadOnlyList<Variant> indels, int refLen)
    {
        foreach (var v in indels)
            AssertWellFormedIndel(v, refLen);
    }

    private static int RefLen(string alignedRef) => alignedRef.Count(c => c != '-');

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-INDEL-001 — Indel Detection (positive sanity)
    // ═════════════════════════════════════════════════════════════════════════

    // ── POSITIVE sanity: a known insertion and a known deletion are each reported
    //    with hand-computed type / length / coordinates (docs §7.1, Evidence
    //    dataset "Alignment-derived indel columns"). Inputs are repeat-free so the
    //    optimal global alignment is unique ⇒ Position is deterministic (ASM-02). ──
    [Test]
    public void FindInsertions_KnownSingleBaseInsertion_HandComputedCall()
    {
        // docs §7.1: ref ATGCAT vs query ATGTCAT ⇒ one T inserted after ref index 2.
        var reference = new DnaSequence("ATGCAT");
        var query = new DnaSequence("ATGTCAT");

        var insertions = VariantCaller.FindInsertions(reference, query).ToList();

        insertions.Should().ContainSingle("exactly one base is inserted ⇒ exactly one insertion (INV-05)");
        var v = insertions[0];
        v.Type.Should().Be(VariantType.Insertion);
        v.ReferenceAllele.Should().Be("-", "insertion ⇒ ref-side gap sentinel (INV-03)");
        v.AlternateAllele.Should().Be("T", "the inserted query base");
        v.Position.Should().Be(3, "the insertion sits at ref column 3 in the unique alignment (docs §7.1)");
        AssertWellFormedIndel(v, refLen: reference.Length);

        // It is an insertion, so FindDeletions reports nothing for this pair.
        VariantCaller.FindDeletions(reference, query).Should().BeEmpty("the pair carries no deletion");
    }

    [Test]
    public void FindDeletions_KnownSingleBaseDeletion_HandComputedCall()
    {
        // Mirror of the insertion case: ref ATGTCAT vs query ATGCAT ⇒ one T deleted.
        var reference = new DnaSequence("ATGTCAT");
        var query = new DnaSequence("ATGCAT");

        var deletions = VariantCaller.FindDeletions(reference, query).ToList();

        deletions.Should().ContainSingle("exactly one base is deleted ⇒ exactly one deletion (INV-05)");
        var v = deletions[0];
        v.Type.Should().Be(VariantType.Deletion);
        v.AlternateAllele.Should().Be("-", "deletion ⇒ query-side gap sentinel (INV-04)");
        v.ReferenceAllele.Should().Be("T", "the deleted reference base");
        v.Position.Should().Be(3, "the deletion sits at ref column 3 in the unique alignment");
        AssertWellFormedIndel(v, refLen: reference.Length);

        VariantCaller.FindInsertions(reference, query).Should().BeEmpty("the pair carries no insertion");
    }

    [Test]
    public void FindIndels_IsUnionOfInsertionsAndDeletions_NoSnps()
    {
        // ref ATCGAT vs query AGCGT : C/G mismatch at index 1 (SNP) AND a deletion
        // of the trailing 'A'. The indel projection must EXCLUDE the SNP.
        var reference = new DnaSequence("ATCGATA");
        var query = new DnaSequence("ATCGAT"); // last base 'A' deleted

        var indels = VariantCaller.FindIndels(reference, query).ToList();
        var insertions = VariantCaller.FindInsertions(reference, query).ToList();
        var deletions = VariantCaller.FindDeletions(reference, query).ToList();

        indels.Should().OnlyContain(v => v.Type == VariantType.Insertion || v.Type == VariantType.Deletion,
            "FindIndels is the union of insertions and deletions — no SNPs (INV-02)");
        indels.Count.Should().Be(insertions.Count + deletions.Count,
            "the union equals insertions ∪ deletions (disjoint classes)");
        deletions.Should().ContainSingle().Which.ReferenceAllele.Should().Be("A");
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-INDEL-001 — BE: "length-0 indel" (no length change ⇒ no record)
    // ═════════════════════════════════════════════════════════════════════════

    // A substitution column is length-PRESERVING ⇒ it is a SNP, never an indel.
    // A "length-0 indel" must never be emitted (INV-01, §6.1).
    [Test]
    public void FindIndels_SubstitutionOnly_NoIndelEmitted()
    {
        // Equal-length, one mismatch ⇒ a SNP, no length change.
        var reference = new DnaSequence("ACGTACGT");
        var query = new DnaSequence("ACGAACGT"); // index 3 T→A

        VariantCaller.FindInsertions(reference, query).Should().BeEmpty("a substitution is length-preserving (§6.1)");
        VariantCaller.FindDeletions(reference, query).Should().BeEmpty("a substitution is length-preserving (§6.1)");
        VariantCaller.FindIndels(reference, query).Should().BeEmpty("no length change ⇒ no indel (INV-01)");
    }

    // Identical sequences ⇒ no length change anywhere ⇒ zero indels (INV-01).
    [Test]
    public void FindIndels_IdenticalSequences_NoIndelEmitted()
    {
        var seq = new DnaSequence("GATTACAGATTACA");
        VariantCaller.FindInsertions(seq, seq).Should().BeEmpty("identical sequences differ nowhere (INV-01)");
        VariantCaller.FindDeletions(seq, seq).Should().BeEmpty("identical sequences differ nowhere (INV-01)");
        VariantCaller.FindIndels(seq, seq).Should().BeEmpty("identical sequences differ nowhere (INV-01)");
    }

    // A gap–gap column (both sides gap) is a degenerate "length-0" alignment
    // column: it changes neither coordinate and must yield NO indel record.
    // Exercised at the raw column-scan level (CallVariantsFromAlignment), the
    // shared core of the indel projections (§5.1).
    [Test]
    public void GapGapColumns_AreNotIndels_NoLengthZeroRecord()
    {
        var calls = VariantCaller.CallVariantsFromAlignment("A-C-G", "A-C-G").ToList();
        calls.Should().BeEmpty("shared gap–gap columns are length-0 non-events ⇒ no indel (INV-01)");

        // A gap–gap column interleaved with a real insertion: only the real event
        // (the single inserted base) is reported, never the empty gap–gap column.
        var mixed = VariantCaller.CallVariantsFromAlignment("A--C", "A-TC").ToList();
        mixed.Should().ContainSingle("only the one ref-gap-vs-base column is an indel");
        var v = mixed[0];
        v.Type.Should().Be(VariantType.Insertion);
        v.AlternateAllele.Should().Be("T");
        AssertWellFormedIndel(v, refLen: RefLen("A--C"));
    }

    // Every emitted indel is, by construction, exactly one base on one side and a
    // gap on the other — a real length change. Fuzz a wide range of pairs and pin
    // that NO indel is ever a length-0 phantom record.
    [Test]
    [CancelAfter(30_000)]
    public void FindIndels_RandomPairs_NeverEmitLengthZeroIndel([Random(1, 1_000_000, 25)] int seed)
    {
        var rng = new Random(seed);
        const string bases = "ACGT";

        for (int t = 0; t < 30; t++)
        {
            var reference = new DnaSequence(RandomDna(rng, bases, 1, 24));
            var query = new DnaSequence(RandomDna(rng, bases, 0, 24));

            List<Variant> indels = null!;
            var call = () => indels = VariantCaller.FindIndels(reference, query).ToList();
            call.Should().NotThrow("any A/C/G/T pair is valid input — must never crash the indel scan");

            AssertAllWellFormedIndels(indels, refLen: reference.Length);

            // Reinforce the "no length-0 indel" target explicitly: every record
            // changes length on exactly one side.
            indels.Should().OnlyContain(
                v => (v.ReferenceAllele == GapAllele) ^ (v.AlternateAllele == GapAllele),
                "every indel has a gap on exactly one side — never a length-0 record");
        }
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-INDEL-001 — BE: "indel at edge" (first / last column)
    // ═════════════════════════════════════════════════════════════════════════

    // Insertion at the VERY FIRST column: the single-base allele slice at the left
    // boundary must not throw IndexOutOfRange / negative-length substring, and the
    // reported Position must be the correct edge coordinate. Exercised on the raw
    // column scan so the alignment is fixed and the edge column is explicit.
    [Test]
    public void Insertion_AtFirstColumn_CorrectCoordinatesNoIndexError()
    {
        // ref "-ACGT" vs query "TACGT" ⇒ a leading insertion of T before ref pos 0.
        List<Variant> calls = null!;
        var act = () => calls = VariantCaller.CallVariantsFromAlignment("-ACGT", "TACGT").ToList();
        act.Should().NotThrow("a leading gap column must not crash the scan");

        calls.Should().ContainSingle();
        var v = calls[0];
        v.Type.Should().Be(VariantType.Insertion);
        v.AlternateAllele.Should().Be("T");
        v.Position.Should().Be(0, "an insertion before the first reference base sits at ref column 0 (edge)");
        AssertWellFormedIndel(v, refLen: RefLen("-ACGT"));
    }

    // Deletion at the VERY LAST column: the right-boundary slice must not throw.
    [Test]
    public void Deletion_AtLastColumn_CorrectCoordinatesNoIndexError()
    {
        // ref "ACGTA" vs query "ACGT-" ⇒ deletion of the trailing 'A' at ref pos 4.
        List<Variant> calls = null!;
        var act = () => calls = VariantCaller.CallVariantsFromAlignment("ACGTA", "ACGT-").ToList();
        act.Should().NotThrow("a trailing gap column must not crash the scan");

        calls.Should().ContainSingle();
        var v = calls[0];
        v.Type.Should().Be(VariantType.Deletion);
        v.ReferenceAllele.Should().Be("A");
        v.Position.Should().Be(4, "the deletion of the last reference base sits at ref column 4 (edge)");
        AssertWellFormedIndel(v, refLen: RefLen("ACGTA"));
    }

    // Edge indels through the TYPED entry point (global alignment, repeat-free
    // inputs so the placement is unique). A leading extra query base and a
    // trailing extra query base each surface as an insertion without any boundary
    // crash; positions stay within [0, reference.Length].
    [Test]
    public void FindInsertions_LeadingExtraQueryBase_EdgeHandled()
    {
        // query has an extra 'T' prepended to a repeat-free reference.
        var reference = new DnaSequence("ACGTAC");
        var query = new DnaSequence("TACGTAC");

        List<Variant> insertions = null!;
        var act = () => insertions = VariantCaller.FindInsertions(reference, query).ToList();
        act.Should().NotThrow("a leading inserted base is an edge event, not a crash");

        insertions.Should().ContainSingle("exactly one extra base ⇒ one insertion (INV-05)");
        insertions[0].AlternateAllele.Should().Be("T");
        AssertAllWellFormedIndels(insertions, refLen: reference.Length);
    }

    [Test]
    public void FindDeletions_TrailingMissingReferenceBase_EdgeHandled()
    {
        // query is the reference with its final base dropped (repeat-free ref).
        var reference = new DnaSequence("ACGTACG");
        var query = new DnaSequence("ACGTAC"); // trailing 'G' deleted

        List<Variant> deletions = null!;
        var act = () => deletions = VariantCaller.FindDeletions(reference, query).ToList();
        act.Should().NotThrow("a trailing deleted base is an edge event, not a crash");

        deletions.Should().ContainSingle("exactly one missing base ⇒ one deletion (INV-05)");
        deletions[0].ReferenceAllele.Should().Be("G");
        deletions[0].Position.Should().BeLessThanOrEqualTo(reference.Length, "edge position stays in-bounds (INV-06)");
        AssertAllWellFormedIndels(deletions, refLen: reference.Length);
    }

    // Fuzz the edges directly on the raw column scan: leading / trailing gap runs
    // on either side of random aligned bodies. The boundary slices must never
    // throw and every indel position must stay in [0, refLen].
    [Test]
    [CancelAfter(30_000)]
    public void CallVariantsFromAlignment_EdgeGapRuns_NeverCrash([Random(1, 1_000_000, 25)] int seed)
    {
        var rng = new Random(seed);

        for (int t = 0; t < 40; t++)
        {
            // Build a random equal-length aligned pair with deliberate gap runs at
            // the two ends so indels land on the first / last columns.
            int body = rng.Next(0, 8);
            int lead = rng.Next(0, 4);
            int trail = rng.Next(0, 4);

            var (aRef, aQuery) = BuildEdgeGapAlignment(rng, lead, body, trail);

            List<Variant> calls = null!;
            var act = () => calls = VariantCaller.CallVariantsFromAlignment(aRef, aQuery).ToList();
            act.Should().NotThrow("equal-length aligned strings with edge gaps must never crash the scan");

            AssertAllWellFormedIndels(
                calls.Where(v => v.Type != VariantType.SNP).ToList(),
                refLen: RefLen(aRef));
        }
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region VARIANT-INDEL-001 — BE: "empty" (empty / null reference or query)
    // ═════════════════════════════════════════════════════════════════════════

    [Test]
    public void FindIndels_BothEmptySequences_NoIndelNoCrash()
    {
        var empty = new DnaSequence("");
        VariantCaller.FindInsertions(empty, empty).Should().BeEmpty("empty vs empty ⇒ no indel (§6.1)");
        VariantCaller.FindDeletions(empty, empty).Should().BeEmpty("empty vs empty ⇒ no indel (§6.1)");
        VariantCaller.FindIndels(empty, empty).Should().BeEmpty("empty vs empty ⇒ no indel (§6.1)");
    }

    // Empty reference vs non-empty query: every query base is "inserted" relative
    // to an empty reference ⇒ all insertions, all at Position 0, no crash.
    [Test]
    public void FindInsertions_EmptyReferenceNonEmptyQuery_AllInsertionsNoCrash()
    {
        var reference = new DnaSequence("");
        var query = new DnaSequence("ACGT");

        List<Variant> insertions = null!;
        var act = () => insertions = VariantCaller.FindInsertions(reference, query).ToList();
        act.Should().NotThrow("empty reference is documented valid input ⇒ empty/guarded result, no crash (§6.1)");

        insertions.Should().OnlyContain(v => v.Type == VariantType.Insertion,
            "every query base is an insertion against an empty reference");
        VariantCaller.FindDeletions(reference, query).Should().BeEmpty("nothing to delete from an empty reference");
        AssertAllWellFormedIndels(insertions, refLen: reference.Length); // refLen == 0 ⇒ Position must be 0
    }

    // Non-empty reference vs empty query: every reference base is "deleted" ⇒ all
    // deletions, no crash, positions in [0, reference.Length].
    [Test]
    public void FindDeletions_NonEmptyReferenceEmptyQuery_AllDeletionsNoCrash()
    {
        var reference = new DnaSequence("ACGT");
        var query = new DnaSequence("");

        List<Variant> deletions = null!;
        var act = () => deletions = VariantCaller.FindDeletions(reference, query).ToList();
        act.Should().NotThrow("empty query is documented valid input ⇒ empty/guarded result, no crash (§6.1)");

        deletions.Should().OnlyContain(v => v.Type == VariantType.Deletion,
            "every reference base is a deletion against an empty query");
        VariantCaller.FindInsertions(reference, query).Should().BeEmpty("nothing to insert into an empty query");
        AssertAllWellFormedIndels(deletions, refLen: reference.Length);
    }

    [Test]
    public void FindIndels_NullReferenceOrQuery_ThrowsArgumentNullException()
    {
        var nullRef = () => VariantCaller.FindIndels(null!, new DnaSequence("ACGT")).ToList();
        var nullQry = () => VariantCaller.FindIndels(new DnaSequence("ACGT"), null!).ToList();
        nullRef.Should().Throw<ArgumentNullException>("null reference is invalid input (§3.3)");
        nullQry.Should().Throw<ArgumentNullException>("null query is invalid input (§3.3)");
    }

    // Fuzz the empty boundary: empty against random queries (and vice versa) must
    // never crash, and the result must be a one-sided indel set (all insertions
    // for empty-ref, all deletions for empty-query).
    [Test]
    [CancelAfter(30_000)]
    public void FindIndels_EmptyAgainstRandom_NeverCrash([Random(1, 1_000_000, 25)] int seed)
    {
        var rng = new Random(seed);
        const string bases = "ACGT";
        var empty = new DnaSequence("");

        for (int t = 0; t < 25; t++)
        {
            var other = new DnaSequence(RandomDna(rng, bases, 0, 18));

            List<Variant> emptyRef = null!;
            var actA = () => emptyRef = VariantCaller.FindIndels(empty, other).ToList();
            actA.Should().NotThrow("empty reference vs random query must never crash");
            emptyRef.Should().OnlyContain(v => v.Type == VariantType.Insertion,
                "an empty reference can only be the target of insertions");
            AssertAllWellFormedIndels(emptyRef, refLen: 0);

            List<Variant> emptyQry = null!;
            var actB = () => emptyQry = VariantCaller.FindIndels(other, empty).ToList();
            actB.Should().NotThrow("random reference vs empty query must never crash");
            emptyQry.Should().OnlyContain(v => v.Type == VariantType.Deletion,
                "an empty query can only witness deletions");
            AssertAllWellFormedIndels(emptyQry, refLen: other.Length);
        }
    }

    #endregion

    // ═════════════════════════════════════════════════════════════════════════
    #region Helpers
    // ═════════════════════════════════════════════════════════════════════════

    private static string RandomDna(Random rng, string bases, int minLen, int maxLen)
    {
        int len = rng.Next(minLen, maxLen + 1);
        var chars = new char[len];
        for (int i = 0; i < len; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    // Builds an equal-length aligned pair with `lead` gap columns at the start and
    // `trail` gap columns at the end (gap on a random side per column ⇒ leading /
    // trailing insertions or deletions), wrapping a random `body` of match /
    // mismatch / gap columns. Returns (alignedRef, alignedQuery) of equal length.
    private static (string aRef, string aQuery) BuildEdgeGapAlignment(Random rng, int lead, int body, int trail)
    {
        const string bases = "ACGT";
        var r = new System.Text.StringBuilder();
        var q = new System.Text.StringBuilder();

        void AppendEdgeGap()
        {
            // gap on exactly one side ⇒ a real edge indel column
            if (rng.Next(2) == 0) { r.Append('-'); q.Append(bases[rng.Next(4)]); }
            else { r.Append(bases[rng.Next(4)]); q.Append('-'); }
        }

        for (int i = 0; i < lead; i++) AppendEdgeGap();

        for (int i = 0; i < body; i++)
        {
            switch (rng.Next(4))
            {
                case 0: { char c = bases[rng.Next(4)]; r.Append(c); q.Append(c); break; }      // match
                case 1: { r.Append(bases[rng.Next(4)]); q.Append(bases[rng.Next(4)]); break; } // (maybe) mismatch
                case 2: { r.Append('-'); q.Append(bases[rng.Next(4)]); break; }                // insertion
                default: { r.Append(bases[rng.Next(4)]); q.Append('-'); break; }               // deletion
            }
        }

        for (int i = 0; i < trail; i++) AppendEdgeGap();

        return (r.ToString(), q.ToString());
    }

    #endregion
}
