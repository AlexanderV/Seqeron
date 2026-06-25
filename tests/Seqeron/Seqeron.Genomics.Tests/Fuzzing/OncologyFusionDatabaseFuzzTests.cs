using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology known-fusion-DATABASE-LOOKUP area — ONCO-FUSION-002.
/// The unit under test is the directional HGNC-designation lookup
/// <see cref="OncologyAnalyzer.MatchKnownFusions"/> (and its designation builder
/// <see cref="OncologyAnalyzer.GetFusionAnnotation"/>), implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// This file is scoped to FUSION-002 DB-LOOKUP / ANNOTATION ONLY. Fusion DETECTION
/// (<c>DetectFusions</c>, ONCO-FUSION-001) and breakpoint/protein analysis
/// (<c>AnalyzeBreakpoint</c>/<c>PredictFusionProtein</c>, ONCO-FUSION-003) are
/// separate checklist rows and are NOT exercised here.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state
/// corruption, no nonsense output, and no *unhandled* runtime exception
/// (DivideByZero / NullReference / Overflow). Every input must resolve to EITHER
/// a well-defined, theory-correct value OR a *documented, intentional* outcome
/// (here, <see cref="ArgumentNullException"/> for a null known-fusion map and
/// <see cref="ArgumentException"/> for a null/empty/whitespace partner symbol).
/// For a directional database lookup the headline hazards are:
///   • a NullReferenceException when the supplied known-fusion collection is EMPTY
///     — an empty set is a legal "nothing is known" DB; every query must resolve
///     to a documented IsKnown=false, never crash on an empty enumeration;
///   • a FALSE POSITIVE near-miss: a partner that differs from a stored entry by a
///     single character (or is a similar-but-distinct gene) MUST NOT match — the
///     contract is an EXACT directional key match, never fuzzy/substring matching;
///   • an ORIENTATION bug: the designation is directional (A::B ≠ B::A, INV-02/03),
///     so a reciprocal query against a DB holding only the forward key must NOT
///     match — swapping partners must change the answer;
///   • case folding: case-varied symbols (eml4/ALK) DO match a stored EML4::ALK
///     (INV-04), regardless of whether the supplied dictionary's own comparer is
///     case-sensitive (the method falls back to a case-insensitive scan).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-FUSION-002 — Known fusion database lookup (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 101.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення. Targets (checklist row 101):
///     "empty DB, exact match, near-miss partner".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Known_Fusion_Database_Lookup.md (docs/algorithms/Oncology/Known_Fusion_Database_Lookup.md):
///   • Designation = gene5p + "::" + gene3p, 5' partner ALWAYS first    (§2.2, INV-01)
///   • Directional: designation(A,B) ≠ designation(B,A) for A≠B; a match requires
///     the directional key 5'::3' — the reciprocal key does NOT match   (§2.2, INV-02/03, §6.1)
///   • A fusion is *known* IFF its exact directional key is a member of the
///     caller-supplied set; the annotation is the caller-supplied value (§2.2, §3.2)
///   • Empty set ⇒ no query is known (no crash on empty enumeration)    (§2.2, §6.1)
///   • Symbol matching is case-insensitive (ordinal-ignore-case); the returned
///     Designation preserves the INPUT case verbatim                    (§3.3, INV-04, §5.2)
///   • Null/empty/whitespace partner symbol ⇒ ArgumentException         (§3.3, §6.1)
///   • Null known-fusion map ⇒ ArgumentNullException                    (§3.3, §6.1)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyFusionDatabaseFuzzTests
{
    // ── FusionCall builder (only the partner symbols matter for the key) ──────
    // Support counts / reading frame are irrelevant to the DB lookup; they are
    // fixed to documented-valid values so the FusionCall is itself well-formed.
    private static FusionCall MakeCall(string g5, string g3) =>
        new(g5, g3, JunctionReads: 5, DiscordantMates: 4, TotalSupport: 9,
            ReadingFrame: FusionReadingFrame.InFrame);

    // Builds a case-insensitive known-fusion DB (the recommended comparer, §3.1).
    private static Dictionary<string, string> KnownDb(
        params (string Designation, string Annotation)[] entries)
    {
        var db = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        foreach (var (designation, annotation) in entries)
        {
            db[designation] = annotation;
        }

        return db;
    }

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented structural contract on EVERY KnownFusionMatch: the
    // Designation is the exact directional gene5p::gene3p key (INV-01) preserving
    // input case, and the IsKnown / Annotation pair is internally consistent
    // (a known match carries the caller's annotation; an unknown one carries
    // null). This is what stops a fuzz test from rubber-stamping an inconsistent
    // (IsKnown=false but Annotation!=null, or wrong designation) result green.
    private static void AssertWellFormedMatch(KnownFusionMatch m, string g5, string g3)
    {
        m.Designation.Should().Be(
            g5 + FusionDesignationSeparator + g3,
            "Designation = gene5p::gene3p, input case preserved verbatim (INV-01, §5.2)");

        if (!m.IsKnown)
        {
            m.Annotation.Should().BeNull("an unknown fusion carries no annotation (§3.2)");
        }
    }

    #region ONCO-FUSION-002 — positive sanity (a real recurrent fusion IS known; a fabricated pair is NOT)

    // The worked example from §7.1: a DB holding the real recurrent drivers
    // BCR::ABL1 and EML4::ALK annotates EML4::ALK as known with the caller's text,
    // while a fabricated pair (FOO::BAR) that is not in the DB is unknown. A fuzz
    // suite that never asserts a TRUE positive AND a TRUE negative proves nothing,
    // so this pins the documented behaviour end-to-end.
    [Test]
    public void MatchKnownFusions_DocumentedWorkedExample_KnownIsAnnotated_FabricatedIsNot()
    {
        var db = KnownDb(
            ("BCR::ABL1", "Chronic myeloid leukemia driver"),
            ("EML4::ALK", "NSCLC driver, ALK TKI target"));

        var known = MatchKnownFusions(MakeCall("EML4", "ALK"), db);
        known.IsKnown.Should().BeTrue("EML4::ALK is an exact directional member of the DB (§7.1)");
        known.Designation.Should().Be("EML4::ALK");
        known.Annotation.Should().Be("NSCLC driver, ALK TKI target", "the caller's annotation is returned (§3.2)");
        AssertWellFormedMatch(known, "EML4", "ALK");

        var fabricated = MatchKnownFusions(MakeCall("FOO", "BAR"), db);
        fabricated.IsKnown.Should().BeFalse("FOO::BAR is not a member of the DB (no false positive)");
        fabricated.Annotation.Should().BeNull();
        AssertWellFormedMatch(fabricated, "FOO", "BAR");
    }

    // GetFusionAnnotation is the directional designation builder (INV-01): 5' first,
    // joined by the HGNC double colon, input case preserved.
    [TestCase("BCR", "ABL1", "BCR::ABL1")]
    [TestCase("EML4", "ALK", "EML4::ALK")]
    [TestCase("TMPRSS2", "ERG", "TMPRSS2::ERG")]
    public void GetFusionAnnotation_BuildsDirectionalHgncDesignation(string g5, string g3, string expected)
    {
        GetFusionAnnotation(g5, g3).Should().Be(expected, "gene5p::gene3p, 5' first (INV-01, §2.2)");
    }

    #endregion

    #region ONCO-FUSION-002 / BE — empty DB (boundary: empty collection)

    // BE target "empty DB": an empty known-fusion set is the lower boundary of DB
    // size. Every query must resolve to a documented IsKnown=false — never a
    // NullReference / crash on enumerating an empty collection (§6.1).
    [Test]
    public void MatchKnownFusions_EmptyDb_EveryQueryIsUnknown_NoCrash()
    {
        var emptyDb = KnownDb(); // zero entries

        var m = MatchKnownFusions(MakeCall("EML4", "ALK"), emptyDb);

        m.IsKnown.Should().BeFalse("an empty DB knows nothing (§6.1)");
        m.Annotation.Should().BeNull();
        AssertWellFormedMatch(m, "EML4", "ALK");
    }

    // The empty boundary must also hold for a CASE-SENSITIVE empty dictionary
    // (a different probe path: the fallback scan iterates an empty collection).
    [Test]
    public void MatchKnownFusions_EmptyCaseSensitiveDb_IsUnknown_NoCrash()
    {
        var emptyDb = new Dictionary<string, string>(StringComparer.Ordinal); // empty, case-sensitive

        var m = MatchKnownFusions(MakeCall("BCR", "ABL1"), emptyDb);

        m.IsKnown.Should().BeFalse("an empty DB knows nothing, whatever its comparer");
    }

    // Fuzz: random partner pairs against an empty DB must ALWAYS be unknown and
    // never throw, whatever the (well-formed, non-empty) symbols are.
    [Test]
    [CancelAfter(10000)]
    public void MatchKnownFusions_FuzzAgainstEmptyDb_AlwaysUnknown([Random(1, 100000, 30)] int seed)
    {
        var rng = new Random(seed);
        var emptyDb = KnownDb();

        string g5 = $"G{rng.Next(0, 100000)}";
        string g3 = $"H{rng.Next(0, 100000)}";

        var m = MatchKnownFusions(MakeCall(g5, g3), emptyDb);

        m.IsKnown.Should().BeFalse("empty DB ⇒ nothing known, for any partner symbols (seed {0})", seed);
        AssertWellFormedMatch(m, g5, g3);
    }

    #endregion

    #region ONCO-FUSION-002 / BE — exact match (boundary: query equals a stored key)

    // BE target "exact match": a query whose directional key is EXACTLY a stored
    // entry is known and carries that entry's annotation (§2.2, §3.2).
    [Test]
    public void MatchKnownFusions_ExactDirectionalKey_IsKnown_WithThatAnnotation()
    {
        var db = KnownDb(("BCR::ABL1", "CML driver"));

        var m = MatchKnownFusions(MakeCall("BCR", "ABL1"), db);

        m.IsKnown.Should().BeTrue("the exact directional key is in the DB");
        m.Annotation.Should().Be("CML driver");
        AssertWellFormedMatch(m, "BCR", "ABL1");
    }

    // INV-04: case-varied symbols match a stored entry case-insensitively, and the
    // returned Designation preserves the INPUT case verbatim (§5.2).
    [TestCase("eml4", "alk")]
    [TestCase("EML4", "ALK")]
    [TestCase("Eml4", "Alk")]
    [TestCase("eMl4", "aLk")]
    public void MatchKnownFusions_CaseVariedQuery_MatchesStoredEntry_PreservesInputCase(string g5, string g3)
    {
        var db = KnownDb(("EML4::ALK", "NSCLC driver"));

        var m = MatchKnownFusions(MakeCall(g5, g3), db);

        m.IsKnown.Should().BeTrue("matching is case-insensitive (INV-04)");
        m.Annotation.Should().Be("NSCLC driver");
        m.Designation.Should().Be(g5 + "::" + g3, "the designation preserves the input case verbatim (§5.2)");
    }

    // Case-insensitive matching must ALSO hold when the supplied dictionary uses a
    // CASE-SENSITIVE comparer — the method falls back to a case-insensitive scan so
    // callers are not silently case-trapped (§5.2). This is the fallback-path test.
    [Test]
    public void MatchKnownFusions_CaseSensitiveDb_StillMatchesCaseInsensitively_ViaFallbackScan()
    {
        var caseSensitiveDb = new Dictionary<string, string>(StringComparer.Ordinal)
        {
            ["EML4::ALK"] = "NSCLC driver",
        };

        var m = MatchKnownFusions(MakeCall("eml4", "alk"), caseSensitiveDb);

        m.IsKnown.Should().BeTrue("the fallback case-insensitive scan finds the differently-cased key (§5.2)");
        m.Annotation.Should().Be("NSCLC driver");
    }

    // Fuzz: for every entry placed in a DB, querying its EXACT partners must return
    // known with that exact annotation — the stored key is always re-findable.
    [Test]
    [CancelAfter(10000)]
    public void MatchKnownFusions_FuzzExactReKey_AlwaysFindsItsOwnEntry([Random(1, 100000, 30)] int seed)
    {
        var rng = new Random(seed);
        int n = rng.Next(1, 25);
        var pairs = new List<(string G5, string G3, string Ann)>(n);
        var db = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);

        for (int i = 0; i < n; i++)
        {
            // Disjoint symbol spaces (A* vs B*) so 5' and 3' can never collide.
            string g5 = $"A{rng.Next(0, 1000)}_{i}";
            string g3 = $"B{rng.Next(0, 1000)}_{i}";
            string ann = $"ann-{rng.Next(0, 1000000)}";
            db[g5 + "::" + g3] = ann;
            pairs.Add((g5, g3, ann));
        }

        foreach (var (g5, g3, ann) in pairs)
        {
            var m = MatchKnownFusions(MakeCall(g5, g3), db);
            m.IsKnown.Should().BeTrue("the exact stored key must re-find its entry (seed {0})", seed);
            m.Annotation.Should().Be(ann, "the matched annotation is the stored one (seed {0})", seed);
            AssertWellFormedMatch(m, g5, g3);
        }
    }

    #endregion

    #region ONCO-FUSION-002 / BE — near-miss partner (boundary: one character off / reciprocal — NO false positive)

    // BE target "near-miss partner": a query where ONE partner differs by a single
    // character from the stored entry must NOT match — matching is exact, never
    // fuzzy. This is the headline false-positive hazard for a DB lookup.
    [TestCase("EML4", "ALKK", TestName = "3' partner off by one char (ALK→ALKK) → not known")]
    [TestCase("EML", "ALK", TestName = "5' partner off by one char (EML4→EML) → not known")]
    [TestCase("EML4", "ALX", TestName = "3' partner one letter changed (ALK→ALX) → not known")]
    [TestCase("FML4", "ALK", TestName = "5' partner one letter changed (EML4→FML4) → not known")]
    public void MatchKnownFusions_NearMissPartner_IsNotKnown_NoFalsePositive(string g5, string g3)
    {
        var db = KnownDb(("EML4::ALK", "NSCLC driver"));

        var m = MatchKnownFusions(MakeCall(g5, g3), db);

        m.IsKnown.Should().BeFalse("a near-miss partner is a DIFFERENT fusion; matching is exact (no false positive)");
        m.Annotation.Should().BeNull();
        AssertWellFormedMatch(m, g5, g3);
    }

    // A query that is a STRICT PREFIX or SUPERSTRING of a stored partner must not
    // match — guards against any accidental substring/StartsWith matching.
    [TestCase("EML", "AL")]
    [TestCase("EML4X", "ALKX")]
    [TestCase("E", "A")]
    public void MatchKnownFusions_SubstringOrSuperstringPartners_AreNotKnown(string g5, string g3)
    {
        var db = KnownDb(("EML4::ALK", "NSCLC driver"));

        var m = MatchKnownFusions(MakeCall(g5, g3), db);

        m.IsKnown.Should().BeFalse("exact-key matching, never substring/prefix matching (no false positive)");
    }

    // INV-02/03: the reciprocal partner order is a DIFFERENT fusion. A DB holding
    // only the forward key A::B must NOT match a reciprocal query B::A.
    [Test]
    public void MatchKnownFusions_ReciprocalQuery_AgainstForwardOnlyDb_IsNotKnown()
    {
        var db = KnownDb(("EML4::ALK", "NSCLC driver"));

        var forward = MatchKnownFusions(MakeCall("EML4", "ALK"), db);
        var reciprocal = MatchKnownFusions(MakeCall("ALK", "EML4"), db);

        forward.IsKnown.Should().BeTrue("the forward directional key is in the DB");
        reciprocal.IsKnown.Should().BeFalse("A::B ≠ B::A; the reciprocal is a different fusion (INV-02/03)");
        reciprocal.Designation.Should().Be("ALK::EML4", "the reciprocal designation reflects the queried order");
    }

    // The "boundary" near-miss right at the colon: a query whose CONCATENATION
    // happens to spell a stored key but with the split in the wrong place must not
    // match (designation is built from the two distinct fields, not a raw string).
    // e.g. stored "AB::C" must NOT be matched by partners ("A","B::C") which build
    // "A::B::C" — a different designation.
    [Test]
    public void MatchKnownFusions_ColonInPartner_BuildsDistinctDesignation_DoesNotCollide()
    {
        var db = KnownDb(("AB::C", "stored"));

        var m = MatchKnownFusions(MakeCall("A", "B::C"), db);

        m.Designation.Should().Be("A::B::C", "the designation is gene5p::gene3p verbatim (INV-01)");
        m.IsKnown.Should().BeFalse("A::B::C ≠ AB::C; no spurious collision across the separator");
    }

    // Broad fuzz: a DB of distinct exact entries, queried with random near-miss
    // mutations of those keys. A mutated partner (different symbol) must NEVER be
    // reported as known — the false-positive guard under many random mutations.
    [Test]
    [CancelAfter(15000)]
    public void MatchKnownFusions_FuzzNearMissMutations_NeverFalsePositive([Random(1, 100000, 40)] int seed)
    {
        var rng = new Random(seed);
        int n = rng.Next(1, 20);
        var db = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        var realPairs = new List<(string G5, string G3)>(n);

        for (int i = 0; i < n; i++)
        {
            string g5 = $"GENE5_{i}";
            string g3 = $"GENE3_{i}";
            db[g5 + "::" + g3] = $"ann{i}";
            realPairs.Add((g5, g3));
        }

        // Query a near-miss: pick a real pair, mutate exactly ONE partner so the
        // designation cannot be any other real key (disjoint indexed namespaces).
        var (baseG5, baseG3) = realPairs[rng.Next(realPairs.Count)];
        bool mutate5 = rng.Next(2) == 0;
        string qG5 = mutate5 ? baseG5 + "X" + rng.Next(0, 9) : baseG5;
        string qG3 = mutate5 ? baseG3 : baseG3 + "X" + rng.Next(0, 9);

        var m = MatchKnownFusions(MakeCall(qG5, qG3), db);

        m.IsKnown.Should().BeFalse(
            "a near-miss partner ({0}::{1}) is not an exact key; matching must not be fuzzy (seed {2})",
            qG5, qG3, seed);
        AssertWellFormedMatch(m, qG5, qG3);
    }

    #endregion

    #region ONCO-FUSION-002 / BE — degenerate inputs (null map, null/empty/whitespace partners — documented throws)

    // A null known-fusion map is the documented ArgumentNullException (§3.3, §6.1).
    [Test]
    public void MatchKnownFusions_NullDb_ThrowsArgumentNullException()
    {
        Action act = () => MatchKnownFusions(MakeCall("EML4", "ALK"), null!);

        act.Should().Throw<ArgumentNullException>("a null known-fusion map is the documented throw (§3.3)");
    }

    // Null / empty / whitespace partner symbols are the documented ArgumentException
    // (§3.3, §6.1) — for BOTH the lookup and the designation builder.
    [TestCase(null, "ALK")]
    [TestCase("EML4", null)]
    [TestCase("", "ALK")]
    [TestCase("EML4", "")]
    [TestCase("   ", "ALK")]
    [TestCase("EML4", "\t")]
    public void MatchKnownFusions_NullEmptyOrWhitespacePartner_ThrowsArgumentException(string? g5, string? g3)
    {
        var db = KnownDb(("EML4::ALK", "NSCLC driver"));

        Action act = () => MatchKnownFusions(MakeCall(g5!, g3!), db);

        act.Should().Throw<ArgumentException>("null/empty/whitespace partner symbols are rejected (§3.3)");
    }

    [TestCase(null, "ALK")]
    [TestCase("EML4", null)]
    [TestCase("", "")]
    [TestCase("  ", "ALK")]
    public void GetFusionAnnotation_NullEmptyOrWhitespacePartner_ThrowsArgumentException(string? g5, string? g3)
    {
        Action act = () => GetFusionAnnotation(g5!, g3!);

        act.Should().Throw<ArgumentException>("the designation builder rejects empty partners (§3.3)");
    }

    // Null map must throw BEFORE any partner work — i.e. a null map with otherwise
    // valid partners still throws ArgumentNullException (validation ordering).
    [Test]
    public void MatchKnownFusions_NullDb_TakesPrecedence_OverValidPartners()
    {
        Action act = () => MatchKnownFusions(MakeCall("BCR", "ABL1"), null!);

        act.Should().Throw<ArgumentNullException>();
    }

    #endregion

    #region ONCO-FUSION-002 / BE — mixed-DB invariants under fuzz

    // Broad fuzz over a populated DB queried with a 50/50 mix of REAL keys and
    // fabricated pairs. Whatever the mix, the result must (a) never throw, (b) be
    // known IFF the directional key is genuinely a member, (c) carry the correct
    // annotation on a hit and null on a miss, and (d) always report the
    // verbatim-case designation.
    [Test]
    [CancelAfter(15000)]
    public void MatchKnownFusions_FuzzMixedQueries_KnownIffMember([Random(1, 100000, 40)] int seed)
    {
        var rng = new Random(seed);
        int n = rng.Next(1, 20);
        var db = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        var members = new List<(string G5, string G3, string Ann)>(n);

        for (int i = 0; i < n; i++)
        {
            string g5 = $"M5_{i}";
            string g3 = $"M3_{i}";
            string ann = $"a{i}";
            db[g5 + "::" + g3] = ann;
            members.Add((g5, g3, ann));
        }

        for (int q = 0; q < 30; q++)
        {
            string g5, g3;
            bool expectKnown;
            string? expectAnn;

            if (rng.Next(2) == 0)
            {
                // Query a genuine member (possibly case-varied).
                var pick = members[rng.Next(members.Count)];
                bool lower = rng.Next(2) == 0;
                g5 = lower ? pick.G5.ToLowerInvariant() : pick.G5;
                g3 = lower ? pick.G3.ToLowerInvariant() : pick.G3;
                expectKnown = true;
                expectAnn = pick.Ann;
            }
            else
            {
                // Query a fabricated pair from a disjoint namespace (never a member).
                g5 = $"X5_{rng.Next(0, 100000)}";
                g3 = $"X3_{rng.Next(0, 100000)}";
                expectKnown = false;
                expectAnn = null;
            }

            var m = MatchKnownFusions(MakeCall(g5, g3), db);

            m.IsKnown.Should().Be(expectKnown,
                "known IFF the directional key is a member ({0}::{1}, seed {2})", g5, g3, seed);
            m.Annotation.Should().Be(expectAnn, "annotation matches membership (seed {0})", seed);
            AssertWellFormedMatch(m, g5, g3);
        }
    }

    #endregion
}
