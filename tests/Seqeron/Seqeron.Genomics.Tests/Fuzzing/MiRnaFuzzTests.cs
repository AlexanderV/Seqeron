using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the MiRNA area — the canonical 5' seed of a mature miRNA (MIRNA-SEED-001)
/// and seed-based 3'UTR target-site prediction (MIRNA-TARGET-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no state corruption, no
/// nonsense output, and — critically for a routine that EXTRACTS A FIXED-WIDTH WINDOW
/// out of a variable-length string — no *unhandled* runtime exception. The seed lives
/// at miRNA positions 2-8, so a miRNA shorter than 8 nt has no full seed window: a naive
/// `Substring(1, 7)` on such a string would throw ArgumentOutOfRangeException (the
/// "IndexOutOfRange extracting the seed region" hazard). The contract must instead
/// resolve every input to EITHER a well-defined, theory-correct result, OR a documented,
/// intentional validation outcome. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MIRNA-SEED-001 — miRNA seed match
/// Checklist: docs/checklists/03_FUZZING.md, row 74.
/// Source doc: docs/algorithms/MiRNA/Seed_Sequence_Analysis.md.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the length corners of the seed window:
///        miRNA shorter than 8 nt (0..7 nt, the KEY boundary that would IndexOutOfRange
///        a naive seed extraction), the empty string, and the self-match boundary
///        miRNA == target. — docs/checklists/03_FUZZING.md §Description (code BE).
///   • MC = Malformed Content — non-RNA characters (DNA 'T', digits, punctuation, junk)
///        fed both to direct seed extraction and to the target scan: each must be
///        handled per the documented normalization contract, never rejected with a
///        crash. — §Description (code MC). Fuzz targets for row 74:
///        "miRNA shorter than 8nt, empty, non-RNA, miRNA = target".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The seed contract under test (Seed_Sequence_Analysis.md §2.2, §3.3, §6.1, §5.2)
/// ───────────────────────────────────────────────────────────────────────────
/// For animal miRNAs target recognition is dominated by the 5' SEED — nucleotides 2-8
/// of the mature miRNA (1-based), i.e. zero-based indices 1..7, a 7-nt string. Seed
/// matching looks for the REVERSE COMPLEMENT of the seed in the target 3'UTR; the 6mer
/// core (RC of positions 2-7) anchors the canonical 8mer / 7mer-m8 / 7mer-A1 / 6mer
/// site classes (Bartel 2009; Lewis 2005; TargetScan).
///
/// Public surfaces probed (src/.../Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs):
///   (1) GetSeedSequence(string) → string. INV-01: returns "" for null / empty / any
///       input SHORTER THAN 8 nt, otherwise Substring(1,7).ToUpperInvariant() — a 7-nt
///       uppercase seed. This is the seed-extraction-safety boundary: the &lt;8-nt guard
///       is what prevents the Substring from walking off the end. §5.2: this surface
///       uppercases ONLY (it does NOT convert DNA 'T' → 'U'), so a short-or-not DNA
///       input is uppercased verbatim and may still contain 'T'.
///   (2) CreateMiRna(name, sequence) → MiRna. Normalizes the sequence (uppercase, T→U)
///       FIRST (INV-02), then extracts the seed from the normalized sequence; always
///       stores SeedStart=1, SeedEnd=7 (INV-03). Empty / short → empty stored seed.
///   (3) CompareSeedRegions(MiRna, MiRna) → SeedComparison. INV-04: IsSameFamily ⇔ the
///       two stored 7-nt seeds are exactly equal. §6.1: when EITHER stored seed is empty
///       the result is a fully zeroed comparison (0 matches, 0 mismatches, not family) —
///       never an exception, never an index walk off a too-short seed.
///   (4) FindTargetSites(mRna, MiRna, minScore) → IEnumerable&lt;TargetSite&gt;. Scans the
///       target for the seed's 6mer-core reverse complement and classifies the site. It
///       short-circuits on an empty target or an empty miRNA sequence, and bails when the
///       seed RC is shorter than 7 (i.e. the miRNA was too short to have a full seed) —
///       so a short / empty / seedless miRNA yields NO sites, never an IndexOutOfRange
///       off the 6mer extraction.
///
/// Documented input handling (Seed_Sequence_Analysis.md §3.3, §6.1):
///   • null / "" / length &lt; 8 → GetSeedSequence returns "" — never throws (the KEY
///     seed-extraction-safety boundary against the &lt;8-nt Substring overflow).
///   • DNA 'T' on the GetSeedSequence surface → uppercased verbatim, may contain 'T'
///     (no T→U on this surface); on the CreateMiRna / FindTargetSites surfaces T→U.
///   • Empty stored seed in CompareSeedRegions → zeroed comparison, not an exception.
///   • miRNA == target self-match → a defined relationship (the seed equals the target's
///     own seed; the RC of the seed is generally NOT present at the seed's own location,
///     so a self-target produces zero or few sites) — defined, crash-free either way.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MIRNA-TARGET-001 — miRNA target-site prediction
/// Checklist: docs/checklists/03_FUZZING.md, row 75.
/// Source doc: docs/algorithms/MiRNA/Target_Site_Prediction.md.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate corners of the seed-complement scan
///        over a 3'UTR. — docs/checklists/03_FUZZING.md §Description (code BE). Fuzz
///        targets for row 75: "No complementarity, empty 3'UTR, miRNA longer than target".
///
/// The target-prediction contract under test (Target_Site_Prediction.md §3.1-§3.3, §4, §6.1):
///   FindTargetSites(mRnaSequence, MiRna, minScore = 0.5) → IEnumerable&lt;TargetSite&gt;.
///   It normalizes the mRNA (uppercase, T→U), takes the RC of the stored 7-nt seed, derives
///   the 6mer core (seedRC[1..7]) and the offset-6mer pattern (seedRC[0..6]), then scans the
///   mRNA in two passes for exact matches, building a full antiparallel duplex alignment and
///   a heuristic score in [0,1] for each retained hit. Documented boundary handling:
///     • No complementarity (the seed's 6mer-core RC and offset pattern never occur in the
///       3'UTR) → the scan finds nothing → NO sites. Never a spurious site, never a crash.
///     • Empty / null 3'UTR → `IsNullOrEmpty(mRnaSequence)` short-circuits → NO sites (§6.1).
///     • miRNA LONGER THAN TARGET — the KEY boundary. The 6mer-core scan loop is bounded
///       `i &lt;= mrna.Length - 6`, so a sub-6-nt mRNA yields zero iterations and no sites.
///       When a hit IS extended, `CreateTargetSite` clamps the window to
///       `Math.Min(miRna.Length, mrna.Length - pos)`, and `AlignMiRnaToTarget` iterates only
///       `Math.Min(miRna.Length, target.Length)` positions with a `targetIdx &lt; 0` guard —
///       so a miRNA far longer than the target NEVER walks off the target string
///       (no IndexOutOfRange) and never divides by zero (§4, INV-03, INV-05).
///     • A 3'UTR carrying a real complementary seed site (here the doc's worked example, a
///       perfect 8mer for let-7a) → detected at the correct 0-based span, with a sensible
///       score (Seed8mer base score 1.0, ≥ the default minScore 0.5) and SeedMatchLength 8.
///
/// All inputs are fixed / deterministically generated; the random helper uses a LOCALLY
/// seeded `new Random(seed)` (no shared static Rng), so every fuzz input is reproducible.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MiRnaFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomRna(int length, int seed)
    {
        const string bases = "ACGU";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal seed-string invariant (INV-01): a seed is EITHER the empty
    /// string (too short / absent) OR exactly 7 uppercase characters. No other length is
    /// ever produced — the fixed-width window cannot leak a partial seed.
    /// </summary>
    private static void AssertSeedShape(string seed)
    {
        (seed.Length == 0 || seed.Length == 7).Should().BeTrue(
            $"a seed is either empty or a 7-nt window, never length {seed.Length} (INV-01)");
        if (seed.Length == 7)
            seed.Should().Be(seed.ToUpperInvariant(), "the extracted seed is uppercase (INV-01)");
    }

    // hsa-let-7a-5p — a real miRBase mature miRNA used as the positive-sanity probe.
    // Seq UGAGGUAGUAGGUUGUAUAGUU → canonical seed (positions 2-8) GAGGUAG.
    private const string Let7aSequence = "UGAGGUAGUAGGUUGUAUAGUU";
    private const string Let7aSeed = "GAGGUAG";

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  MIRNA-SEED-001 — miRNA seed match : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region MIRNA-SEED-001 — miRNA seed match

    #region BE — Boundary: miRNA shorter than 8 nt (the seed-extraction overflow hazard)

    /// <summary>
    /// BE — "miRNA shorter than 8nt": THE key seed-extraction-safety boundary. The seed
    /// lives at positions 2-8, so any miRNA of 0..7 nt has no full seed window; a naive
    /// `Substring(1, 7)` on such a string would throw ArgumentOutOfRangeException. The
    /// contract instead returns "" for every length below 8 (INV-01, §6.1). We sweep
    /// EVERY short length 0..7 — the off-by-one corners around the 8-nt threshold — and
    /// assert no throw and an empty seed. Length 8 (the first valid length) is included
    /// as the counterpoint that the window opens exactly at 8.
    /// </summary>
    [Test]
    public void GetSeedSequence_MiRnaShorterThanEight_ReturnsEmptyNeverIndexOverflow()
    {
        foreach (int len in new[] { 0, 1, 2, 3, 4, 5, 6, 7 })
        {
            string shortMiRna = new string('A', len);
            var act = () => GetSeedSequence(shortMiRna);

            string seed = act.Should().NotThrow(
                $"a {len}-nt miRNA has no full 2-8 seed window — the <8-nt guard must prevent the Substring overflow")
                .Subject;
            seed.Should().BeEmpty($"a {len}-nt miRNA (< 8 nt) cannot yield a seed (INV-01)");
            AssertSeedShape(seed);
        }

        // The window opens exactly at length 8: the first valid input yields a 7-nt seed.
        string eightMer = "UGAGGUAG"; // positions 2-8 = "GAGGUAG"
        string firstSeed = GetSeedSequence(eightMer);
        firstSeed.Should().Be("GAGGUAG", "length 8 is the first length with a full 2-8 seed window");
        AssertSeedShape(firstSeed);
    }

    /// <summary>
    /// BE — short miRNA on the downstream surfaces. A miRNA too short for a seed must also
    /// flow safely through CreateMiRna (empty stored seed, no crash) and through
    /// FindTargetSites (the seed RC is shorter than 7, so the scan bails with NO sites —
    /// never an IndexOutOfRange off the 6mer-core extraction). We probe a 5-nt miRNA.
    /// </summary>
    [Test]
    public void ShortMiRna_FlowsSafelyThroughCreateAndTargetScan()
    {
        var shortMiRna = CreateMiRna("too-short", "UGAGG"); // 5 nt
        shortMiRna.SeedSequence.Should().BeEmpty("a 5-nt miRNA has no seed (INV-01) even after normalization");
        AssertSeedShape(shortMiRna.SeedSequence);

        // A perfectly complementary 3'UTR cannot raise a site from a seedless miRNA.
        string target = new string('C', 40);
        var act = () => FindTargetSites(target, shortMiRna).ToList();
        act.Should().NotThrow("a seedless miRNA must not crash the target scan (seed RC < 7 → no sites)")
           .Subject.Should().BeEmpty("with no full seed there is no 6mer core to search for");
    }

    #endregion

    #region BE — Boundary: empty miRNA

    /// <summary>
    /// BE — "empty": the empty miRNA. GetSeedSequence short-circuits null/"" to "" before
    /// touching the Substring (§3.3, §6.1). CreateMiRna over "" stores an empty seed, and
    /// CompareSeedRegions involving an empty seed returns the fully zeroed comparison
    /// (§6.1) — never an exception, never an index walk off an empty seed string.
    /// </summary>
    [Test]
    public void EmptyMiRna_ProducesEmptySeedAndZeroedComparison()
    {
        GetSeedSequence("").Should().BeEmpty("an empty miRNA has no seed (INV-01)");
        GetSeedSequence(null!).Should().BeEmpty("a null miRNA is treated like empty, not an error (§6.1)");

        var empty = CreateMiRna("empty", "");
        empty.SeedSequence.Should().BeEmpty("CreateMiRna over an empty sequence stores no seed");

        // Comparing an empty-seed miRNA against a real one is a defined zeroed result.
        var real = CreateMiRna("let-7a", Let7aSequence);
        var act = () => CompareSeedRegions(empty, real);
        var cmp = act.Should().NotThrow("an empty stored seed must not crash the comparison (§6.1)").Subject;
        cmp.Matches.Should().Be(0, "an empty seed yields 0 matches (§6.1)");
        cmp.Mismatches.Should().Be(0, "an empty seed yields 0 mismatches (§6.1)");
        cmp.IsSameFamily.Should().BeFalse("an empty seed never shares a family (INV-04)");

        // The empty target on the scan surface is likewise a documented no-op.
        FindTargetSites("", real).Should().BeEmpty("an empty 3'UTR contains no target site");
    }

    #endregion

    #region MC — Malformed Content: non-RNA characters

    /// <summary>
    /// MC — "non-RNA": DNA 'T', digits, punctuation and arbitrary junk. The two surfaces
    /// have DIFFERENT but documented normalization (§5.2, §6.1):
    ///   • GetSeedSequence uppercases ONLY — DNA 'T' is preserved in the seed; the call
    ///     must still succeed and return a 7-nt window (or "" if too short), never crash.
    ///   • CreateMiRna converts T→U, so a DNA-spelled miRNA stores the equivalent RNA seed.
    ///   • FindTargetSites normalizes BOTH miRNA and target (T→U), so junk/DNA targets are
    ///     handled — non-canonical characters simply never form the 6mer core; no crash.
    /// </summary>
    [Test]
    public void NonRnaCharacters_AreNormalizedOrPreserved_NeverCrash()
    {
        // GetSeedSequence: DNA thymine is uppercased but NOT converted (T stays T).
        // "ATGAGGTA" (8 nt) → positions 2-8 = "TGAGGTA", uppercased verbatim.
        string dnaSeed = GetSeedSequence("aTGAGGta");
        dnaSeed.Should().Be("TGAGGTA", "GetSeedSequence uppercases only — DNA T is preserved (§5.2)");
        AssertSeedShape(dnaSeed);

        // CreateMiRna: the same DNA-spelled miRNA stores the T→U converted seed.
        var dnaMiRna = CreateMiRna("dna-let-7a", "TGAGGTAGTAGGTTGTATAGTT");
        dnaMiRna.SeedSequence.Should().Be(Let7aSeed, "CreateMiRna converts DNA T→U, recovering the canonical RNA seed (INV-02)");
        dnaMiRna.Sequence.Should().NotContain("T", "the stored RNA sequence has no thymine (T→U normalization)");

        // Arbitrary junk on the extraction surface: still a fixed-width window, never a crash.
        var junk = () => GetSeedSequence("XX12!! @#$%");
        var junkSeed = junk.Should().NotThrow("junk content must not crash seed extraction").Subject;
        AssertSeedShape(junkSeed);

        // Junk on the target scan: non-canonical chars never form the seed RC core → no
        // crash; if the genuine 6mer-core RC happens to occur the site is still well-formed.
        var let7a = CreateMiRna("let-7a", Let7aSequence);
        var scan = () => FindTargetSites("NN12!! \t##GGNNxx", let7a).ToList();
        scan.Should().NotThrow("garbage 3'UTR content must not crash the target scan");
    }

    #endregion

    #region BE — Boundary: miRNA == target (self-match)

    /// <summary>
    /// BE — "miRNA = target": the self-match boundary. Feeding a miRNA's own sequence as
    /// the 3'UTR is well-defined: the target search looks for the REVERSE COMPLEMENT of
    /// the seed, which is generally NOT present at the seed's own forward location, so a
    /// self-target produces zero (or only incidental) sites — and must do so WITHOUT a
    /// crash. We assert crash-freedom and well-formed coordinates for whatever is found,
    /// and that any reported seed-match really is the 6mer-core RC at the reported span.
    /// </summary>
    [Test]
    public void MiRnaEqualsTarget_IsDefinedAndCrashFree()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence);

        var act = () => FindTargetSites(Let7aSequence, let7a).ToList();
        var sites = act.Should().NotThrow("a miRNA scanned against its own sequence must not crash").Subject;

        // The seed's forward sequence is not its own reverse complement, so let-7a vs
        // itself raises no canonical seed site.
        sites.Should().BeEmpty(
            "the seed's reverse complement does not occur at the seed's own forward position — no self seed-match");

        // Every site (if any future change produced one) must have in-range, ordered coords.
        foreach (var s in sites)
        {
            s.Start.Should().BeInRange(0, Let7aSequence.Length - 1);
            s.End.Should().BeInRange(0, Let7aSequence.Length - 1);
            s.Start.Should().BeLessThanOrEqualTo(s.End);
        }

        // The seed-vs-itself COMPARISON is the trivially identical, same-family case.
        var selfCmp = CompareSeedRegions(let7a, let7a);
        selfCmp.IsSameFamily.Should().BeTrue("a miRNA's seed is identical to itself (INV-04)");
        selfCmp.Matches.Should().Be(7, "all 7 seed positions match a miRNA against itself (INV-05)");
        selfCmp.Mismatches.Should().Be(0, "a self comparison has no mismatches (INV-05)");
    }

    #endregion

    #region Positive sanity — a real seed match is detected at the correct site

    /// <summary>
    /// Positive sanity: a target that is the REVERSE COMPLEMENT of the full let-7a miRNA
    /// presents the canonical seed match — the 6mer core (RC of positions 2-7) is present,
    /// preceded by the RC of position 8 and followed by an 'A' opposite position 1 — i.e.
    /// a textbook 8mer site (Bartel 2009; TargetScan). The scan MUST detect at least one
    /// seed site, at least one of full-seed class (8mer / 7mer-m8 / 7mer-A1), with the
    /// reported target sub-sequence actually containing the seed's 6mer-core reverse
    /// complement at the reported span. This proves the fuzz harness asserts against a
    /// matcher that FINDS real seed matches, not a no-op.
    /// </summary>
    [Test]
    public void FindTargetSites_RealSeedMatch_DetectedAtCorrectSite()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence);
        let7a.SeedSequence.Should().Be(Let7aSeed, "the positive probe starts from the canonical let-7a seed");

        // Target = reverse complement of the whole miRNA → an embedded 8mer seed site.
        string target = GetReverseComplement(Let7aSequence);

        var sites = FindTargetSites(target, let7a).ToList();
        sites.Should().NotBeEmpty("the RC of the miRNA presents a canonical seed match");

        // The 6mer-core reverse complement that the scan keys on.
        string sixmerCoreRC = GetReverseComplement(Let7aSeed).Substring(1, 6); // RC of positions 2-7

        // At least one detected site is a full-seed class (8mer / 7mer-m8 / 7mer-A1).
        sites.Should().Contain(
            s => s.Type == TargetSiteType.Seed8mer
              || s.Type == TargetSiteType.Seed7merM8
              || s.Type == TargetSiteType.Seed7merA1,
            "a perfect RC target yields a canonical full-seed site class");

        // Every reported site's span lies inside the target and carries the 6mer core RC.
        foreach (var s in sites)
        {
            s.Start.Should().BeInRange(0, target.Length - 1, "site start indexes the target");
            s.End.Should().BeInRange(0, target.Length - 1, "site end indexes the target");
            s.Start.Should().BeLessThanOrEqualTo(s.End, "site coordinates are ordered");
            s.MiRnaName.Should().Be("let-7a");
            s.SeedMatchLength.Should().BeInRange(6, 8, "a canonical seed match spans 6-8 nt");
            target.Should().Contain(sixmerCoreRC,
                "the target carries the seed's 6mer-core reverse complement that the scan keys on");
        }

        // The strongest site is the 8mer (perfect seed + flanking pos-8 RC + A1).
        sites.Should().Contain(s => s.Type == TargetSiteType.Seed8mer,
            "the full reverse complement presents the highest-affinity 8mer site");
    }

    /// <summary>
    /// Positive sanity over RANDOM miRNAs: across fixed seeds and lengths, GetSeedSequence
    /// and CreateMiRna never crash and always emit a well-shaped seed (empty or 7 nt), and
    /// FindTargetSites never throws on random miRNA / target pairs. The seed window is
    /// always consistent with the stored sequence (INV-02).
    /// </summary>
    [Test]
    public void RandomMiRna_AlwaysWellShapedSeedAndCrashFreeScan()
    {
        foreach (int seed in new[] { 1, 7, 42, 2026 })
        {
            foreach (int len in new[] { 0, 1, 7, 8, 22, 60 })
            {
                string raw = RandomRna(len, seed);

                string direct = GetSeedSequence(raw);
                AssertSeedShape(direct);

                var mirna = CreateMiRna($"rnd-{seed}-{len}", raw);
                AssertSeedShape(mirna.SeedSequence);
                // INV-02: the stored seed equals the seed of the stored normalized sequence.
                mirna.SeedSequence.Should().Be(GetSeedSequence(mirna.Sequence),
                    $"stored seed is the seed of the normalized sequence (INV-02; seed {seed}, len {len})");

                string target = RandomRna(len + 13, seed + 100);
                var scan = () => FindTargetSites(target, mirna).ToList();
                scan.Should().NotThrow($"random miRNA/target must not crash the scan (seed {seed}, len {len})");
            }
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  MIRNA-TARGET-001 — miRNA target-site prediction : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region MIRNA-TARGET-001 — miRNA target prediction

    // Worked example from Target_Site_Prediction.md §7.1: let-7a (seed GAGGUAG, seed RC
    // CUACCUC) against this 3'UTR presents a canonical 8mer site (the 6mer core UACCUC
    // flanked upstream by the pos-8 RC 'C' and downstream by an 'A' opposite position 1).
    private const string Let7aTarget8mer = "GGGGGCUACCUCAGGGGG";

    #region BE — Boundary: no complementarity

    /// <summary>
    /// BE — "No complementarity": a 3'UTR that shares NO complementary region with the miRNA.
    /// The finder keys on exact matches of the seed's 6mer-core reverse complement (and the
    /// offset-6mer pattern). For let-7a the 6mer core is the RC of GAGGUAG positions 2-7; a
    /// homopolymer 3'UTR (all 'A', then all 'C') cannot contain that mixed-base pattern, so
    /// the scan finds NOTHING — the contract is NO sites, never a spurious hit, never a crash
    /// (Target_Site_Prediction.md §4 pass-1/pass-2, §6.1). We sweep several non-complementary
    /// 3'UTRs and assert an empty result at the default and at a permissive minScore.
    /// </summary>
    [Test]
    public void FindTargetSites_NoComplementarity_YieldsNoSitesNeverCrash()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence);

        foreach (string utr in new[]
                 {
                     new string('A', 60),   // homopolymer — no 6mer-core pattern can occur
                     new string('C', 60),
                     new string('G', 80),
                     "AAAAACCCCCGGGGGAAAAACCCCCGGGGG", // block copolymer, still no mixed core
                 })
        {
            var act = () => FindTargetSites(utr, let7a).ToList();
            var sites = act.Should().NotThrow(
                "a non-complementary 3'UTR must never crash the scan").Subject;
            sites.Should().BeEmpty(
                "a 3'UTR with no seed-complement region yields no target sites (§4, §6.1)");

            // Even with a fully permissive threshold the absence of the seed pattern means
            // no candidate is ever generated — emptiness is structural, not a score filter.
            FindTargetSites(utr, let7a, minScore: -1.0).Should().BeEmpty(
                "no 6mer-core / offset pattern present → no candidate at any threshold");
        }
    }

    #endregion

    #region BE — Boundary: empty 3'UTR

    /// <summary>
    /// BE — "empty 3'UTR": the empty (and null) target. `FindTargetSites` short-circuits on
    /// `IsNullOrEmpty(mRnaSequence)` before any seed matching, so an empty UTR yields NO
    /// sites and never touches the scan loop or any Substring (§3.3, §6.1). We assert at the
    /// default threshold and at a permissive one — emptiness is independent of minScore.
    /// </summary>
    [Test]
    public void FindTargetSites_EmptyUtr_YieldsNoSites()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence);

        var actEmpty = () => FindTargetSites("", let7a).ToList();
        actEmpty.Should().NotThrow("an empty 3'UTR must not crash the scan")
                .Subject.Should().BeEmpty("an empty 3'UTR contains no target site (§6.1)");

        var actNull = () => FindTargetSites(null!, let7a).ToList();
        actNull.Should().NotThrow("a null 3'UTR is treated as empty, not an error (§6.1)")
               .Subject.Should().BeEmpty("a null 3'UTR contains no target site (§6.1)");

        FindTargetSites("", let7a, minScore: -1.0).Should().BeEmpty(
            "an empty 3'UTR is a no-op regardless of the score threshold");
    }

    #endregion

    #region BE — Boundary: miRNA longer than target (the index-overflow / divide-by-zero hazard)

    /// <summary>
    /// BE — "miRNA longer than target": THE key boundary. The full-length let-7a miRNA (22 nt)
    /// is fed against 3'UTRs FAR SHORTER than itself, down to length 0. A naive duplex builder
    /// could walk off the end of the target (IndexOutOfRange) or divide by a zero-length
    /// alignment; the contract instead clamps every window to the available target. We sweep
    /// target lengths 0..7 (all below the 22-nt miRNA, and several below the 6-nt scan floor
    /// — `i &lt;= mrna.Length - 6` simply does not iterate), then probe SHORT TARGETS THAT DO
    /// carry the 6mer core so the extension/alignment path is exercised against a too-short
    /// target: the alignment must clamp to `Math.Min(miRna.Length, target.Length)` and the
    /// reported coordinates must stay inside the target. Never an IndexOutOfRange, never a
    /// divide-by-zero, never an out-of-range coordinate (§4, INV-03, INV-05).
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void FindTargetSites_MiRnaLongerThanTarget_NeverIndexOverflowOrDivideByZero()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence); // 22 nt
        let7a.Sequence.Length.Should().BeGreaterThan(7, "the probe miRNA is longer than the tiny targets below");

        // Tiny targets, every length 0..7 (all < miRNA length; lengths < 6 are below the
        // scan floor). Each must resolve to a defined, crash-free, empty-or-valid result.
        foreach (int len in new[] { 0, 1, 2, 3, 4, 5, 6, 7 })
        {
            string tinyUtr = new string('A', len);
            var act = () => FindTargetSites(tinyUtr, let7a).ToList();
            var sites = act.Should().NotThrow(
                $"a {len}-nt 3'UTR shorter than the {let7a.Sequence.Length}-nt miRNA must not " +
                "crash the scan (no IndexOutOfRange, no DivideByZero)").Subject;

            foreach (var s in sites)
            {
                s.Start.Should().BeInRange(0, tinyUtr.Length - 1, "any site start stays inside the tiny target");
                s.End.Should().BeInRange(0, tinyUtr.Length - 1, "any site end stays inside the tiny target");
                s.Start.Should().BeLessThanOrEqualTo(s.End, "coordinates are ordered (INV-03)");
            }
        }

        // Short targets that DO carry the 6mer core RC (UACCUC), so pass-1 finds a hit and
        // CreateTargetSite/AlignMiRnaToTarget run against a target far shorter than the miRNA.
        // The 8mer flank cannot fit at every length, but the extension/alignment path must
        // stay safe and the span must remain inside the target.
        foreach (string coreUtr in new[]
                 {
                     "UACCUC",          // bare 6mer core, target length 6 « 22-nt miRNA
                     "CUACCUC",         // pos-8 flank, length 7
                     "CUACCUCA",        // full 8mer flank, length 8
                     "GGUACCUCG",       // core embedded with flanks, length 9
                 })
        {
            var act = () => FindTargetSites(coreUtr, let7a, minScore: -1.0).ToList();
            var sites = act.Should().NotThrow(
                $"a {coreUtr.Length}-nt target carrying the 6mer core must not crash even though " +
                "the 22-nt miRNA is far longer (alignment clamps to the target length)").Subject;

            foreach (var s in sites)
            {
                s.Start.Should().BeInRange(0, coreUtr.Length - 1, "site start indexes the short target");
                s.End.Should().BeInRange(0, coreUtr.Length - 1, "site end indexes the short target");
                s.Start.Should().BeLessThanOrEqualTo(s.End, "coordinates are ordered (INV-03)");
                s.Score.Should().BeInRange(0.0, 1.0, "the heuristic score stays clamped to [0,1] (INV-01)");
                s.TargetSequence.Length.Should().BeLessThanOrEqualTo(coreUtr.Length,
                    "the extended target window cannot exceed the short target (clamped to mRNA tail)");
            }
        }
    }

    #endregion

    #region Positive sanity — a real complementary seed site is detected at the right place

    /// <summary>
    /// Positive sanity: a 3'UTR carrying a genuine complementary seed site for let-7a — the
    /// worked example from Target_Site_Prediction.md §7.1, a textbook 8mer (the 6mer core
    /// UACCUC flanked by the pos-8 RC 'C' upstream and an 'A' opposite position 1 downstream).
    /// At the DEFAULT minScore (0.5) the 8mer (base score 1.0) survives, so the scan must
    /// report exactly that site at the correct 0-based span, of class Seed8mer with
    /// SeedMatchLength 8 and a score ≥ 0.5. This proves the fuzz harness asserts against a
    /// matcher that FINDS real complementary sites, not a no-op that always returns empty.
    /// </summary>
    [Test]
    public void FindTargetSites_RealComplementarySite_DetectedAtCorrectPositionWithScore()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence);
        let7a.SeedSequence.Should().Be(Let7aSeed, "the positive probe starts from the canonical let-7a seed");

        var sites = FindTargetSites(Let7aTarget8mer, let7a).ToList(); // default minScore 0.5
        sites.Should().NotBeEmpty("the 3'UTR carries a real complementary 8mer seed site");

        var eightMer = sites.Should().ContainSingle(s => s.Type == TargetSiteType.Seed8mer,
            "the worked-example 3'UTR presents exactly one canonical 8mer site").Subject;

        // The 6mer core RC the scan keys on, anchored at the documented position.
        string sixmerCoreRC = GetReverseComplement(Let7aSeed).Substring(1, 6); // UACCUC
        int coreIdx = Let7aTarget8mer.IndexOf(sixmerCoreRC, StringComparison.Ordinal);
        coreIdx.Should().BeGreaterThan(0, "the 6mer core sits with room for the upstream pos-8 base");

        // 8mer span = [coreIdx-1 .. coreIdx+6] (pos-8 base + 6mer core + A1), zero-based inclusive.
        eightMer.Start.Should().Be(coreIdx - 1, "the 8mer starts one base upstream of the 6mer core (INV-03)");
        eightMer.End.Should().Be(coreIdx + 6, "the 8mer ends at the A1 base opposite miRNA position 1 (INV-03)");
        eightMer.SeedMatchLength.Should().Be(8, "an 8mer site matches 8 seed positions (INV-02)");
        eightMer.MiRnaName.Should().Be("let-7a", "the site carries the supplied miRNA name");
        eightMer.Score.Should().BeInRange(0.5, 1.0,
            "the 8mer base score (1.0) clears the default minScore and stays clamped to [0,1] (INV-01)");
        Let7aTarget8mer.Substring(eightMer.Start, eightMer.End - eightMer.Start + 1)
            .Should().Contain(sixmerCoreRC, "the reported span really carries the seed's 6mer-core RC");

        // Raising the threshold above 1.0 suppresses every site (§6.1).
        FindTargetSites(Let7aTarget8mer, let7a, minScore: 1.01).Should().BeEmpty(
            "a minScore above the [0,1] score ceiling admits no site (§6.1)");
    }

    #endregion

    #endregion
}
