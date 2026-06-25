// PARSE-EMBL-001 — EMBL/INSDC feature-location remote-aware sequence assembly
// Evidence: docs/Evidence/PARSE-EMBL-001-Evidence.md (Enhancement 2026-06-26 section)
// TestSpec: tests/TestSpecs/PARSE-EMBL-001.md (RESV-01..RESV-17)
// Source: INSDC. The DDBJ/ENA/GenBank Feature Table Definition (v11.x), §3.4 Location / §3.5
//         Operators. https://www.insdc.org/submitting-standards/feature-table/
//         (DDBJ mirror: https://www.ddbj.nig.ac.jp/ddbj/feature-table-e.html; accessed 2026-06-26)
//
// Assembly rules (all hand-derived in the TestSpec / Evidence):
//   - Base 1 = first base (5' end); span n..m selects bases n..m, 1-based inclusive.
//   - join(a,b,...) / order(a,b,...): concatenate elements in listed order.
//   - complement(loc): reverse-complement of the enclosed span, so
//       complement(join(a,b)) == join(complement(b),complement(a))  (segment order reversed).
//   - remote accession[.version]:span -> caller-supplied resolver returns the remote
//       sequence; the library slices it 1-based inclusive.
//   - < / > partials -> slice the stated number verbatim (only available coordinate).
//   - missing/null resolver -> remote element contributes empty string; local segments intact.

using NUnit.Framework;
using Seqeron.Genomics.IO;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class FeatureLocationHelper_ResolveLocationSequence_Tests
{
    // Local entry used across cases: positions 1..10 = A C G T A C G T A C.
    private const string Local = "ACGTACGTAC";

    // Resolver returning known remote sequences. Each base is traceable.
    //   J00194 (15 bp): positions 1..15 = G G G G G C C C C C A A A A A
    //   X (4 bp):       T T G G
    //   Y (4 bp):       A A A C
    private static FeatureLocationHelper.RemoteSequenceResolver StubResolver =>
        (acc, _) => acc switch
        {
            "J00194" => "GGGGGCCCCCAAAAA",
            "X" => "TTGG",
            "Y" => "AAAC",
            _ => null
        };

    #region ResolveLocationSequence — local / mixed / remote ordering

    // RESV-01 — join(1..3,7..9): 1..3 = "ACG", 7..9 = "GTA" -> "ACGGTA". Remote resolver
    // must never be invoked for a location with no remote element (§3.5 join).
    [Test]
    public void ResolveLocationSequence_LocalOnlyJoin_ResolverNotInvoked()
    {
        bool invoked = false;
        FeatureLocationHelper.RemoteSequenceResolver tracking = (a, v) => { invoked = true; return null; };

        var result = FeatureLocationHelper.ResolveLocationSequence("join(1..3,7..9)", Local, tracking);

        Assert.Multiple(() =>
        {
            Assert.That(result, Is.EqualTo("ACGGTA"),
                "join(1..3,7..9) concatenates bases 1-3 (ACG) then 7-9 (GTA) in listed order (INSDC §3.5 join).");
            Assert.That(invoked, Is.False,
                "A local-only location must not invoke the caller's remote resolver.");
        });
    }

    // RESV-02 — J00194.1:5..14 on remote "GGGGGCCCCCAAAAA": bases 5..14 (1-based inclusive)
    // = G(5) C C C C C(10) A A A A(14) -> "GCCCCCAAAA" (§3.4(e), 1-based inclusive remote slice).
    [Test]
    public void ResolveLocationSequence_RemoteOnly_SlicesOneBasedInclusive()
    {
        var result = FeatureLocationHelper.ResolveLocationSequence("J00194.1:5..14", "", StubResolver);

        Assert.That(result, Is.EqualTo("GCCCCCAAAA"),
            "Remote span 5..14 selects bases 5-14 of the resolver's 'GGGGGCCCCCAAAAA' (INSDC §3.4(e), 1-based inclusive).");
    }

    // RESV-03 — join(1..10,J00194.1:5..14): local 1..10 = "ACGTACGTAC", remote 5..14 =
    // "GCCCCCAAAA" -> "ACGTACGTACGCCCCCAAAA" (interleaved order preserved).
    [Test]
    public void ResolveLocationSequence_MixedLocalAndRemote_PreservesOrder()
    {
        var result = FeatureLocationHelper.ResolveLocationSequence(
            "join(1..10,J00194.1:5..14)", Local, StubResolver);

        Assert.That(result, Is.EqualTo("ACGTACGTACGCCCCCAAAA"),
            "Mixed join concatenates local 1..10 (ACGTACGTAC) then remote J00194 5..14 (GCCCCCAAAA) in listed order.");
    }

    #endregion

    #region ResolveLocationSequence — complement semantics (§3.5)

    // RESV-04 — complement(Y.1:1..4): remote Y 1..4 = "AAAC"; reverse-complement("AAAC") =
    // complement each (TTTG) then reverse -> "GTTT" (§3.5 complement reads complement 5'->3').
    [Test]
    public void ResolveLocationSequence_ComplementOfRemote_ReverseComplements()
    {
        var result = FeatureLocationHelper.ResolveLocationSequence("complement(Y.1:1..4)", "", StubResolver);

        Assert.That(result, Is.EqualTo("GTTT"),
            "complement(Y.1:1..4) reverse-complements remote 'AAAC' -> 'GTTT' (INSDC §3.5 complement).");
    }

    // RESV-05 — complement(join(1..5,X.1:1..4)) on local "ACGTA", remote X = "TTGG":
    // inner join = "ACGTA"+"TTGG" = "ACGTATTGG"; reverse-complement("ACGTATTGG"):
    // complement -> "TGCATAACC", reverse -> "CCAATACGT".
    [Test]
    public void ResolveLocationSequence_ComplementOfJoin_ReversesOrderAndComplements()
    {
        var result = FeatureLocationHelper.ResolveLocationSequence(
            "complement(join(1..5,X.1:1..4))", "ACGTA", StubResolver);

        Assert.That(result, Is.EqualTo("CCAATACGT"),
            "complement(join(1..5,X.1:1..4)) joins ACGTA+TTGG then reverse-complements the whole -> CCAATACGT (INSDC §3.5).");
    }

    // RESV-06 — the two forms must be identical (INSDC §3.5 equivalence):
    // complement(join(1..5,X.1:1..4)) == join(complement(X.1:1..4),complement(1..5)).
    [Test]
    public void ResolveLocationSequence_ComplementOfJoin_EqualsJoinOfComplements()
    {
        var complementOfJoin = FeatureLocationHelper.ResolveLocationSequence(
            "complement(join(1..5,X.1:1..4))", "ACGTA", StubResolver);
        var joinOfComplements = FeatureLocationHelper.ResolveLocationSequence(
            "join(complement(X.1:1..4),complement(1..5))", "ACGTA", StubResolver);

        Assert.Multiple(() =>
        {
            Assert.That(complementOfJoin, Is.EqualTo(joinOfComplements),
                "INSDC §3.5: complement(join(a,b)) and join(complement(b),complement(a)) produce identical results.");
            Assert.That(complementOfJoin, Is.EqualTo("CCAATACGT"),
                "Both forms equal the hand-derived CCAATACGT.");
        });
    }

    // RESV-07 — the INSDC §3.5 spec equivalence example, on a constructed local sequence:
    //   complement(join(2691..4571,4918..5163)) == join(complement(4918..5163),complement(2691..4571)).
    // Local = 2690 'A' + 1881 'C' (positions 2691..4571) + 346 'A' (4572..4917) + 246 'G' (4918..5163).
    // join = 1881*C + 246*G; reverse-complement = 246*C + 1881*G.
    [Test]
    public void ResolveLocationSequence_SpecEquivalenceExample_LargeJoin()
    {
        string local = new string('A', 2690) + new string('C', 1881)
                       + new string('A', 346) + new string('G', 246);

        var form1 = FeatureLocationHelper.ResolveLocationSequence(
            "complement(join(2691..4571,4918..5163))", local, StubResolver);
        var form2 = FeatureLocationHelper.ResolveLocationSequence(
            "join(complement(4918..5163),complement(2691..4571))", local, StubResolver);
        string expected = new string('C', 246) + new string('G', 1881);

        Assert.Multiple(() =>
        {
            Assert.That(form1, Is.EqualTo(form2),
                "INSDC §3.5 spec example: the two location forms are equivalent.");
            Assert.That(form1, Is.EqualTo(expected),
                "Reverse-complement of (1881 C + 246 G) is (246 C + 1881 G).");
        });
    }

    #endregion

    #region ResolveLocationSequence — missing / null resolver

    // RESV-08 — join(1..3,J00194.1:5..14) with a null resolver: remote contributes nothing;
    // local 1..3 = "ACG" is still assembled -> "ACG".
    [Test]
    public void ResolveLocationSequence_NullResolver_RemoteSegmentEmpty()
    {
        var result = FeatureLocationHelper.ResolveLocationSequence(
            "join(1..3,J00194.1:5..14)", Local, null);

        Assert.That(result, Is.EqualTo("ACG"),
            "With no resolver, the remote element contributes nothing; local 1..3 (ACG) is still assembled.");
    }

    // RESV-09 — resolver returns null for the requested accession: same defined behaviour.
    [Test]
    public void ResolveLocationSequence_ResolverReturnsNull_RemoteSegmentEmpty()
    {
        FeatureLocationHelper.RemoteSequenceResolver nullReturning = (_, __) => null;

        var result = FeatureLocationHelper.ResolveLocationSequence(
            "join(1..3,J00194.1:5..14)", Local, nullReturning);

        Assert.That(result, Is.EqualTo("ACG"),
            "A resolver returning null leaves the remote element empty; local 1..3 (ACG) is unaffected.");
    }

    #endregion

    #region ResolveLocationSequence — partials, single base, order

    // RESV-10 — <1..5: the '<' partial marker is ignored for slicing; bases 1..5 = "ACGTA".
    [Test]
    public void ResolveLocationSequence_PartialStartMarker_SlicesStatedNumber()
    {
        var result = FeatureLocationHelper.ResolveLocationSequence("<1..5", Local, StubResolver);

        Assert.That(result, Is.EqualTo("ACGTA"),
            "The '<' partial marker uses the stated number 1 (only available coordinate); bases 1..5 = ACGTA.");
    }

    // RESV-11 — join(1..3,>7..10): local 1..3 = "ACG", >7..10 slices stated 7..10 = "GTAC"
    // (positions 7,8,9,10 = G T A C) -> "ACGGTAC".
    [Test]
    public void ResolveLocationSequence_PartialEndMarker_SlicesStatedNumber()
    {
        var result = FeatureLocationHelper.ResolveLocationSequence("join(1..3,>7..10)", Local, StubResolver);

        Assert.That(result, Is.EqualTo("ACGGTAC"),
            "The '>' partial marker uses the stated bound; 1..3 (ACG) + 7..10 (GTAC) = ACGGTAC.");
    }

    // RESV-12 — single base number "2" selects base 2 = "C" (§3.4(a)).
    [Test]
    public void ResolveLocationSequence_SingleBase_SelectsOneBase()
    {
        var result = FeatureLocationHelper.ResolveLocationSequence("2", Local, StubResolver);

        Assert.That(result, Is.EqualTo("C"),
            "A single base number n selects base n; base 2 of ACGTACGTAC is C (INSDC §3.4(a)).");
    }

    // RESV-13 — order(1..3,7..9) assembles the same ordered concatenation as join -> "ACGGTA".
    [Test]
    public void ResolveLocationSequence_OrderOperator_AssemblesLikeJoin()
    {
        var result = FeatureLocationHelper.ResolveLocationSequence("order(1..3,7..9)", Local, StubResolver);

        Assert.That(result, Is.EqualTo("ACGGTA"),
            "order(...) assembles the same ordered concatenation as join for sequence extraction (INSDC §3.5).");
    }

    #endregion

    #region ResolveLocationSequence — defensive + overload + resolver arguments

    // RESV-14 — empty / null raw location returns the empty string (contract parity with ExtractSequence).
    [Test]
    public void ResolveLocationSequence_EmptyOrNullLocation_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(FeatureLocationHelper.ResolveLocationSequence("", Local, StubResolver),
                Is.EqualTo(string.Empty), "Empty location yields the empty string.");
            Assert.That(FeatureLocationHelper.ResolveLocationSequence(null!, Local, StubResolver),
                Is.EqualTo(string.Empty), "Null location yields the empty string (no throw).");
        });
    }

    // RESV-15 — the EmblParser.Location overload uses the parsed RawLocation and reproduces RESV-03.
    [Test]
    public void ResolveLocationSequence_EmblLocationOverload_UsesRawLocation()
    {
        var location = EmblParser.ParseLocation("join(1..10,J00194.1:5..14)");

        var result = FeatureLocationHelper.ResolveLocationSequence(location, Local, StubResolver);

        Assert.That(result, Is.EqualTo("ACGTACGTACGCCCCCAAAA"),
            "The Location overload assembles from RawLocation, matching the raw-string overload (RESV-03).");
    }

    // RESV-16 — the resolver receives the accession AND the parsed version (J00194.1 -> version 1).
    [Test]
    public void ResolveLocationSequence_RemoteVersionPassedToResolver()
    {
        string? seenAcc = null;
        int? seenVer = -1;
        FeatureLocationHelper.RemoteSequenceResolver capturing = (acc, ver) =>
        {
            seenAcc = acc; seenVer = ver;
            return "GGGGGCCCCCAAAAA";
        };

        FeatureLocationHelper.ResolveLocationSequence("J00194.1:5..14", "", capturing);

        Assert.Multiple(() =>
        {
            Assert.That(seenAcc, Is.EqualTo("J00194"),
                "The resolver receives the bare accession (no version digits).");
            Assert.That(seenVer, Is.EqualTo(1),
                "The resolver receives the parsed sequence version (.1 -> 1).");
        });
    }

    // RESV-17 — a remote reference with no version passes version == null to the resolver.
    [Test]
    public void ResolveLocationSequence_RemoteNoVersion_PassesNullVersion()
    {
        int? seenVer = -1;
        FeatureLocationHelper.RemoteSequenceResolver capturing = (acc, ver) =>
        {
            seenVer = ver;
            return "GGGGGCCCCCAAAAA";
        };

        var result = FeatureLocationHelper.ResolveLocationSequence("J00194:5..14", "", capturing);

        Assert.Multiple(() =>
        {
            Assert.That(seenVer, Is.Null,
                "A remote reference without an explicit version passes null to the resolver (INSDC §3.4(e)).");
            Assert.That(result, Is.EqualTo("GCCCCCAAAA"),
                "Slicing is unaffected by the missing version: bases 5..14 = GCCCCCAAAA.");
        });
    }

    #endregion
}
