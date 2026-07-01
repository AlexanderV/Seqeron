using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

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
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MIRNA-PRECURSOR-001 — miRNA precursor (pre-miRNA hairpin) detection
/// Checklist: docs/checklists/03_FUZZING.md, row 76.
/// Source doc: docs/algorithms/MiRNA/Pre_miRNA_Detection.md.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the length / structure corners of the precursor
///        scan: candidates SHORTER THAN the minimum hairpin window (0 .. minHairpinLength-1,
///        including empty and null), and sequences that have NO hairpin (no uninterrupted
///        end-to-end complementary stem). — docs/checklists/03_FUZZING.md §Description (code BE).
///   • MC = Malformed Content — homopolymer (all-one-base) runs that are self-identical
///        and therefore NEVER self-complementary, plus non-RNA / junk characters (DNA 'T'
///        which is normalized to 'U', digits, punctuation, 'N' ambiguity codes). Each must
///        be handled per the documented contract, never crash. — §Description (code MC).
///        Fuzz targets for row 76: "Too short, no hairpin, all one base, non-RNA chars".
///
/// The precursor contract under test (Pre_miRNA_Detection.md §2.2, §3.1-§3.3, §6.1):
///   FindPreMiRnaHairpins(sequence, minHairpinLength=55, maxHairpinLength=120,
///     matureLength=22) → IEnumerable&lt;PreMiRna&gt;. It normalizes the input to uppercase
///   RNA (T→U), slides a start across the sequence, and for every candidate window of length
///   [minHairpinLength, maxHairpinLength] requires UNINTERRUPTED complementary pairing
///   (Watson-Crick + G:U wobble) from both ends inward. A candidate is accepted only when the
///   uninterrupted stem is ≥ 18 bp AND the residual loop (n − 2·stem) is in [3, 25] nt. Each
///   accepted hairpin yields a PreMiRna with a balanced dot-bracket Structure of the same
///   length as Sequence (INV-01), stem ≥ 18 (INV-02), loop in [3,25] (INV-03), zero-based
///   inclusive Start/End in the scanned input (INV-04), and equal-length 5' Mature / 3' Star
///   arms of length min(matureLength, stem) (INV-05). Documented boundary handling (§3.3, §6.1):
///     • null / "" / length &lt; minHairpinLength → NO candidates; never throws (the early
///       `IsNullOrEmpty || Length &lt; minHairpinLength` guard prevents any Substring overflow).
///     • No hairpin (random / non-pairing sequence with no uninterrupted complementary stem)
///       → the stem count stops at the first mirrored mismatch, falls below 18 bp → rejected
///       → NO precursor. Never a spurious call, never a crash.
///     • All-one-base homopolymer (e.g. all 'A') → a base is never complementary to itself
///       (A:A, C:C, G:G never pair; only U:U is also a non-pair) → stem = 0 → NO hairpin →
///       NO precursor, at any window length.
///     • Non-RNA chars: DNA 'T' is normalized to 'U' (a DNA-spelled hairpin folds like its RNA
///       form); digits / punctuation / 'N' never pair → simply fail the stem test → NO crash.
///     • A constructed perfect synthetic hairpin (≥ 18-bp uninterrupted stem + 3-25 nt loop,
///       55-120 nt) → ACCEPTED as a precursor with a balanced structure and finite energy.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MIRNA-PAIR-001 — miRNA–target base-pairing (antiparallel duplex)
/// Checklist: docs/checklists/03_FUZZING.md, row 225.
/// Source doc: docs/algorithms/MiRNA/MiRNA_Target_Pairing.md.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate corners of the antiparallel
///        duplex pairing primitive. — docs/checklists/03_FUZZING.md §Description (code BE).
///        Fuzz targets for row 225: "no complementarity, perfect match, short miRNA".
///
/// The pairing contract under test (MiRNA_Target_Pairing.md §2.2, §3.3, §6.1, §7.1):
///   AlignMiRnaToTarget(miRnaSequence, targetSequence) → MiRnaDuplex. Both inputs are
///   normalised (uppercase, T→U). The duplex is ungapped and ANTIPARALLEL: miRNA index i
///   pairs with target index len(target)−1−i, over the overlap min(len(m),len(t)) (§2.2,
///   §3.3, §4.1). Every overlap position is classified into EXACTLY ONE class (INV-05):
///     • Watson-Crick match {A-U,U-A,G-C,C-G} → alignment symbol `|`, Matches++ (§2.2[3]).
///     • G:U wobble {G-U,U-G} → alignment symbol `:`, GUWobbles++ — counted SEPARATELY
///       from canonical matches (§2.2[4], §6.1; INV-02). Wobble ⊆ pairable (INV-03).
///     • Mismatch otherwise → space, Mismatches++ (A pairs only with U, C only with G).
///   Gaps is always 0 (ungapped). FreeEnergy is a Turner-2004 nearest-neighbor stacking
///   sum over runs of consecutive paired positions; a fully Watson-Crick duplex has
///   FreeEnergy ≤ 0, and an all-mismatch duplex has FreeEnergy = 0 (no paired stacks),
///   per INV-06. Supporting predicates: CanPair (INV-01), IsWobblePair (INV-02),
///   GetReverseComplement (reverse + complement A↔U/G↔C, T→A, length-preserving, INV-04).
///   Documented BE boundary handling (§3.3, §6.1):
///     • No complementarity — a target sharing no complementary base at any antiparallel
///       position yields a fully unpaired alignment: Matches 0, GUWobbles 0, Mismatches
///       = overlap, AlignmentString all spaces, FreeEnergy 0 — defined, never a crash.
///     • Perfect match — a target that is the antiparallel Watson-Crick complement of the
///       miRNA yields Matches = overlap, GUWobbles 0, Mismatches 0, AlignmentString all
///       `|`, FreeEnergy ≤ 0 (the doc's "AAAA"/"UUUU" worked example, §7.1).
///     • Short miRNA — a miRNA SHORTER than the target (down to empty / null): the overlap
///       is the shorter length, the antiparallel index never walks off either string
///       (no IndexOutOfRange), and an empty/null input returns the empty MiRnaDuplex
///       (all counts 0, empty strings) — never an exception (§3.3, §6.1).
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

    // ═══════════════════════════════════════════════════════════════════
    //  MIRNA-PRECURSOR-001 — miRNA precursor (pre-miRNA hairpin) detection
    // ═══════════════════════════════════════════════════════════════════

    #region MIRNA-PRECURSOR-001 — miRNA precursor detection

    #region Helpers — synthetic hairpin construction

    /// <summary>
    /// Builds a perfect synthetic pre-miRNA hairpin: a 5' stem, a non-pairing loop, and a
    /// 3' stem that is the positional reverse complement of the 5' stem, so that mirrored
    /// positions from the two ends inward pair uninterruptedly (the implementation's
    /// strongest-case topology — Pre_miRNA_Detection.md §6.1, ASM-01). The loop is a
    /// homopolymer 'A' run, whose first base (A) does not pair with the mirrored 3'-stem
    /// base, so the uninterrupted stem stops exactly at the stem boundary.
    /// </summary>
    private static string BuildHairpin(string stem5, int loopLength)
    {
        // 3' stem = reverse complement of 5' stem so position j pairs with position n-1-j.
        string stem3 = GetReverseComplement(stem5); // RC reverses + complements → positional mirror
        return stem5 + new string('A', loopLength) + stem3;
    }

    /// <summary>
    /// Asserts the universal precursor invariants on a reported PreMiRna (INV-01..INV-05):
    /// balanced dot-bracket structure of equal length, stem ≥ 18 bp, loop in [3,25] nt,
    /// in-range ordered coordinates, equal-length mature/star arms, and a finite energy.
    /// </summary>
    private static void AssertPreMiRnaShape(MiRnaAnalyzer.PreMiRna p, int scannedLength)
    {
        p.Structure.Length.Should().Be(p.Sequence.Length,
            "the dot-bracket structure has the same length as the candidate sequence (INV-01)");
        int open = p.Structure.Count(c => c == '(');
        int close = p.Structure.Count(c => c == ')');
        open.Should().Be(close, "the dot-bracket stem is balanced — equal '(' and ')' (INV-01)");
        open.Should().BeGreaterThanOrEqualTo(18, "an accepted precursor has stem ≥ 18 bp (INV-02)");

        int loop = p.Structure.Count(c => c == '.');
        loop.Should().BeInRange(3, 25, "an accepted precursor has loop length in [3,25] nt (INV-03)");

        p.Start.Should().BeInRange(0, scannedLength - 1, "Start is a zero-based index into the scanned input (INV-04)");
        p.End.Should().BeInRange(0, scannedLength - 1, "End is a zero-based index into the scanned input (INV-04)");
        p.Start.Should().BeLessThanOrEqualTo(p.End, "coordinates are ordered (INV-04)");
        (p.End - p.Start + 1).Should().Be(p.Sequence.Length, "the span width equals the candidate length (INV-04)");

        p.MatureSequence.Length.Should().Be(p.StarSequence.Length,
            "mature and star arms are extracted with the same bounded length (INV-05)");

        double.IsNaN(p.FreeEnergy).Should().BeFalse("the hairpin free energy is finite, never NaN");
        double.IsInfinity(p.FreeEnergy).Should().BeFalse("the hairpin free energy is finite, never infinite");
    }

    #endregion

    #region BE — Boundary: candidate too short (below the minimum hairpin window)

    /// <summary>
    /// BE — "Too short": a sequence shorter than the minimum hairpin window cannot contain a
    /// precursor. `FindPreMiRnaHairpins` guards with `IsNullOrEmpty(sequence) || Length &lt;
    /// minHairpinLength` BEFORE any window Substring (§3.3, §6.1), so every too-short input —
    /// including the empty string and null — yields NO candidates and never throws. We sweep a
    /// representative ladder of short lengths up to one below the default 55-nt floor, then add
    /// the empty and null corners. The first valid scan length (= minHairpinLength) is exercised
    /// in the positive-sanity test as the counterpoint.
    /// </summary>
    [Test]
    public void FindPreMiRnaHairpins_TooShort_YieldsNoCandidatesNeverThrows()
    {
        foreach (int len in new[] { 0, 1, 2, 10, 30, 53, 54 }) // all < default minHairpinLength (55)
        {
            string tooShort = new string('A', len);
            var act = () => FindPreMiRnaHairpins(tooShort).ToList();
            var hits = act.Should().NotThrow(
                $"a {len}-nt input is shorter than the 55-nt hairpin floor — the guard prevents any window overflow")
                .Subject;
            hits.Should().BeEmpty($"a {len}-nt input cannot contain a 55-120 nt precursor (§3.3, §6.1)");
        }

        FindPreMiRnaHairpins("").Should().BeEmpty("an empty sequence yields no precursor (§6.1)");

        var actNull = () => FindPreMiRnaHairpins(null!).ToList();
        actNull.Should().NotThrow("a null sequence is treated as empty, not an error (§6.1)")
               .Subject.Should().BeEmpty("a null sequence yields no precursor (§6.1)");

        // Even a candidate of EXACTLY the right length but below the stem/loop thresholds is
        // rejected — too-short here also means too-short a stem. A 55-nt non-pairing run yields none.
        FindPreMiRnaHairpins(new string('A', 55)).Should().BeEmpty(
            "a 55-nt homopolymer reaches the window floor but has no stem → no precursor");
    }

    #endregion

    #region BE — Boundary: no hairpin (no uninterrupted complementary stem)

    /// <summary>
    /// BE — "No hairpin": a sequence whose ends are NOT complementary has no uninterrupted
    /// stem. The stem counter stops at the first mirrored mismatch, so a random / deliberately
    /// non-pairing sequence falls below the 18-bp stem floor and is rejected → NO precursor
    /// (§4 step 3-4, §6.1). We probe (a) several fixed-seed random sequences long enough to be
    /// scanned, and (b) a hand-built sequence whose first and last bases cannot pair (A vs C),
    /// forcing stem = 0. None may produce a candidate, none may crash.
    /// </summary>
    [Test]
    public void FindPreMiRnaHairpins_NoHairpin_YieldsNoPrecursorNeverCrash()
    {
        // Hand-built: every mirrored pair is A-vs-C (A:C never pairs) → stem stops at 0.
        string noStem = new string('A', 40) + new string('C', 40); // 80 nt, ends A…C never pair
        var act0 = () => FindPreMiRnaHairpins(noStem).ToList();
        act0.Should().NotThrow("a sequence with no complementary ends must not crash the scan")
            .Subject.Should().BeEmpty("no uninterrupted stem → no precursor (§4, §6.1)");

        // Random non-structured sequences: overwhelmingly no uninterrupted ≥18-bp stem.
        foreach (int seed in new[] { 1, 7, 42, 2026 })
        {
            foreach (int len in new[] { 60, 90, 120 })
            {
                string raw = RandomRna(len, seed);
                var act = () => FindPreMiRnaHairpins(raw).ToList();
                var hits = act.Should().NotThrow(
                    $"a random {len}-nt sequence must not crash the precursor scan (seed {seed})").Subject;

                // Whatever the scan returns, every candidate must still satisfy the contract:
                // a real ≥18-bp uninterrupted stem with a 3-25 nt loop. Random hits are rare but
                // when they occur they must be well-formed, never malformed.
                foreach (var p in hits)
                    AssertPreMiRnaShape(p, raw.Length);
            }
        }
    }

    #endregion

    #region MC — Malformed Content: all-one-base homopolymer

    /// <summary>
    /// MC — "all one base": a homopolymer run has no self-complementarity — A:A, C:C, G:G
    /// never pair (and U:U does not pair either), so a single-base sequence has stem = 0 at
    /// every window and can never form a hairpin → NO precursor (§2.2 pairing set, §6.1). We
    /// sweep all four RNA bases (and DNA 'T', which normalizes to 'U') across lengths spanning
    /// the scan window, asserting an always-empty, always crash-free result.
    /// </summary>
    [Test]
    public void FindPreMiRnaHairpins_Homopolymer_NeverFormsHairpin()
    {
        foreach (char b in new[] { 'A', 'C', 'G', 'U', 'T' }) // T normalizes to U → still homopolymer
        {
            foreach (int len in new[] { 55, 60, 80, 120, 150 })
            {
                string homo = new string(b, len);
                var act = () => FindPreMiRnaHairpins(homo).ToList();
                act.Should().NotThrow(
                    $"a {len}-nt all-'{b}' homopolymer must not crash the scan")
                    .Subject.Should().BeEmpty(
                    $"a homopolymer base never pairs with itself → no hairpin → no precursor (base '{b}', len {len})");
            }
        }
    }

    #endregion

    #region MC — Malformed Content: non-RNA characters

    /// <summary>
    /// MC — "non-RNA chars": DNA 'T', digits, punctuation, whitespace, and 'N' ambiguity codes.
    /// Per §3.3/§6.1 the input is normalized to uppercase RNA (T→U) before scanning, and any
    /// character outside {A,U,G,C} simply never satisfies the pairing test (CanPair → false on
    /// 'N'/digits/punctuation), so junk content fails the stem requirement rather than crashing.
    /// We verify: (a) a DNA-spelled valid hairpin still folds (T→U equivalence), (b) arbitrary
    /// junk is scanned without a throw and yields no malformed candidate.
    /// </summary>
    [Test]
    public void FindPreMiRnaHairpins_NonRnaCharacters_NormalizedOrRejected_NeverCrash()
    {
        // (a) The same valid hairpin spelled with DNA 'T' folds identically (T→U normalization).
        string rnaStem = "GCGCGCGCGCAUAUGCGCGC";        // 20 nt, GC-rich for a clean stem
        string rnaHairpin = BuildHairpin(rnaStem, 16);  // 56 nt, stem 20, loop 16 → accepted
        string dnaHairpin = rnaHairpin.Replace('U', 'T').ToLowerInvariant(); // DNA + lowercase junk

        var rnaHits = FindPreMiRnaHairpins(rnaHairpin).ToList();
        var dnaHits = FindPreMiRnaHairpins(dnaHairpin).ToList();
        rnaHits.Should().NotBeEmpty("a perfect synthetic RNA hairpin is accepted as a precursor");
        dnaHits.Should().NotBeEmpty("the same hairpin spelled with DNA 'T' folds identically (T→U, §6.1)");
        dnaHits.Should().HaveCount(rnaHits.Count,
            "T→U + uppercasing makes the DNA spelling equivalent to the RNA spelling");

        // (b) Arbitrary junk content: digits, punctuation, whitespace, 'N' — never pairs, never crashes.
        string junk = new string('N', 30) + "12!! \t@#$%^&*()" + new string('N', 30); // ≥55 nt, non-pairing
        var actJunk = () => FindPreMiRnaHairpins(junk).ToList();
        var junkHits = actJunk.Should().NotThrow("garbage content must not crash the precursor scan").Subject;
        junkHits.Should().BeEmpty("non-canonical characters never pair → no stem → no precursor");

        // (c) A valid hairpin with junk flanking context is still found and well-formed.
        string flanked = "12!!##" + rnaHairpin + "$$%%NN";
        var flankedHits = FindPreMiRnaHairpins(flanked).ToList();
        flankedHits.Should().NotBeEmpty("a valid hairpin embedded in junk flanks is still detected");
        foreach (var p in flankedHits)
            AssertPreMiRnaShape(p, flanked.Length);
    }

    #endregion

    #region Positive sanity — a constructed pre-miRNA hairpin is recognized as a precursor

    /// <summary>
    /// Positive sanity: a constructed perfect pre-miRNA hairpin — a 20-bp uninterrupted GC-rich
    /// stem closing a 16-nt loop (56 nt, inside the 55-120 nt window) — MUST be recognized as a
    /// precursor (§6.1 "perfect synthetic hairpin … Accepted"). The reported candidate carries a
    /// balanced dot-bracket structure with a 20-bp stem (INV-01/02), a 16-nt loop (INV-03),
    /// in-range ordered coordinates (INV-04), equal-length 5' mature / 3' star arms (INV-05), and
    /// a finite Turner free energy. This proves the fuzz harness asserts against a detector that
    /// FINDS real hairpins, not a no-op that always returns empty.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void FindPreMiRnaHairpins_ConstructedHairpin_RecognizedAsPrecursor()
    {
        string stem5 = "GCGCGCGCGCAUAUGCGCGC"; // 20 nt
        const int loopLen = 16;
        string hairpin = BuildHairpin(stem5, loopLen); // 20 + 16 + 20 = 56 nt
        hairpin.Length.Should().Be(56, "the constructed hairpin sits inside the 55-120 nt precursor window");

        var sites = FindPreMiRnaHairpins(hairpin).ToList();
        sites.Should().NotBeEmpty("a perfect synthetic hairpin is accepted as a precursor (§6.1)");

        // Every reported candidate satisfies the full precursor contract.
        foreach (var p in sites)
            AssertPreMiRnaShape(p, hairpin.Length);

        // The strongest candidate is the full-length window: a 20-bp stem closing a 16-nt loop.
        var full = sites.Should().Contain(p => p.Sequence == hairpin,
            "the full constructed window is itself a valid precursor candidate").Subject;
        full.Structure.Count(c => c == '(').Should().Be(20, "the constructed stem is 20 bp");
        full.Structure.Count(c => c == '.').Should().Be(loopLen, "the constructed loop is 16 nt");
        full.MatureSequence.Length.Should().Be(20, "mature arm = min(matureLength 22, stem 20) (INV-05)");
        full.MatureSequence.Should().Be(stem5, "the 5' mature arm is the 5' stem of the constructed hairpin");
        full.Structure.Should().StartWith("((((((((((", "the structure opens with the 5' stem parentheses (INV-01)");
        full.Structure.Should().EndWith("))))))))))", "the structure closes with the 3' stem parentheses (INV-01)");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  MIRNA-PAIR-001 — miRNA–target base-pairing (antiparallel duplex) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region MIRNA-PAIR-001 — miRNA-target pairing

    #region Helpers — duplex invariant assertions

    /// <summary>
    /// Asserts the universal duplex invariants of MiRNA_Target_Pairing.md on a reported
    /// MiRnaDuplex over a non-empty overlap of width <paramref name="overlap"/>:
    /// INV-05 — Matches + Mismatches + GUWobbles == overlap, Gaps == 0, AlignmentString
    /// has exactly <paramref name="overlap"/> symbols drawn from { '|', ':', ' ' }, and the
    /// per-symbol counts agree with the totals. FreeEnergy is finite (never NaN/Infinity).
    /// </summary>
    private static void AssertDuplexShape(MiRnaAnalyzer.MiRnaDuplex d, int overlap)
    {
        d.Gaps.Should().Be(0, "the aligner is ungapped — Gaps is always 0 (§3.2)");
        d.AlignmentString.Length.Should().Be(overlap,
            "the alignment has exactly one symbol per overlap position (INV-05)");
        d.AlignmentString.ToCharArray().Should().OnlyContain(c => c == '|' || c == ':' || c == ' ',
            "every alignment symbol is Watson-Crick '|', wobble ':' or mismatch space (§3.2)");

        (d.Matches + d.Mismatches + d.GUWobbles).Should().Be(overlap,
            "every overlap position is classified into exactly one class with no gaps (INV-05)");
        d.Matches.Should().Be(d.AlignmentString.Count(c => c == '|'),
            "Matches counts the Watson-Crick '|' symbols (INV-05)");
        d.GUWobbles.Should().Be(d.AlignmentString.Count(c => c == ':'),
            "GUWobbles counts the wobble ':' symbols (INV-05)");
        d.Mismatches.Should().Be(d.AlignmentString.Count(c => c == ' '),
            "Mismatches counts the space symbols (INV-05)");

        double.IsNaN(d.FreeEnergy).Should().BeFalse("the duplex free energy is finite, never NaN");
        double.IsInfinity(d.FreeEnergy).Should().BeFalse("the duplex free energy is finite, never infinite");
    }

    #endregion

    #region Predicate sanity — CanPair / IsWobblePair encode the documented pairing set (INV-01..INV-04)

    /// <summary>
    /// Positive sanity over the pairing PRIMITIVES, derived independently from the spec
    /// (MiRNA_Target_Pairing.md §2.2 / INV-01..INV-04; Watson-Crick A-U,G-C [3] + Crick
    /// 1966 G:U wobble [4]). CanPair is true ⇔ {A-U,U-A,G-C,C-G,G-U,U-G}; IsWobblePair is
    /// true ⇔ {G-U,U-G} and the four Watson-Crick pairs are NOT wobble (INV-02); wobble ⊆
    /// pairable (INV-03). T is treated as U (case-insensitive). Every other ordered base
    /// pair (the full 6×6 alphabet incl. 'N') is exhaustively asserted non-pairing, so a
    /// matcher that mis-classified A:G or A:C as a pair would FAIL — the test pins the rule,
    /// not the code's table. GetReverseComplement is the doc's let-7a seed example (INV-04).
    /// </summary>
    [Test]
    public void PairingPredicates_EncodeDocumentedWatsonCrickAndWobbleRule()
    {
        // The complete set of valid pairs, derived from the spec — NOT read off the code.
        var validPairs = new HashSet<(char, char)>
        {
            ('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G'), ('G', 'U'), ('U', 'G'),
        };
        var wobblePairs = new HashSet<(char, char)> { ('G', 'U'), ('U', 'G') };

        // 'T' must behave exactly like 'U' (T→U normalization, §3.1). Sweep the alphabet.
        const string alphabet = "ACGUTN";
        foreach (char a in alphabet)
        {
            foreach (char b in alphabet)
            {
                char an = a == 'T' ? 'U' : a;
                char bn = b == 'T' ? 'U' : b;

                bool expectPair = validPairs.Contains((an, bn));
                bool expectWobble = wobblePairs.Contains((an, bn));

                MiRnaAnalyzer.CanPair(a, b).Should().Be(expectPair,
                    $"CanPair('{a}','{b}') ⇔ Watson-Crick or G:U wobble (INV-01)");
                MiRnaAnalyzer.IsWobblePair(a, b).Should().Be(expectWobble,
                    $"IsWobblePair('{a}','{b}') ⇔ {{G-U,U-G}} only (INV-02)");

                if (expectWobble)
                    MiRnaAnalyzer.CanPair(a, b).Should().BeTrue("wobble ⊆ pairable (INV-03)");
            }
        }

        // A:G and A:C are explicitly NON-pairing (A pairs only with U) — the key biology.
        MiRnaAnalyzer.CanPair('A', 'G').Should().BeFalse("A pairs only with U (§2.2[3])");
        MiRnaAnalyzer.CanPair('A', 'C').Should().BeFalse("A pairs only with U (§2.2[3])");
        MiRnaAnalyzer.CanPair('C', 'U').Should().BeFalse("C pairs only with G, no C:U wobble (§2.2[3])");

        // Watson-Crick pairs are emphatically NOT wobble (INV-02).
        foreach (var (a, b) in new[] { ('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G') })
            MiRnaAnalyzer.IsWobblePair(a, b).Should().BeFalse(
                $"the Watson-Crick pair {a}:{b} is not a wobble (INV-02)");

        // GetReverseComplement: reverse + complement, length-preserving (INV-04) — the doc's
        // let-7a-5p seed example, derived by hand: reverse(GAGGUAG)=GAUGGAG, complement→CUACCUC.
        MiRnaAnalyzer.GetReverseComplement(Let7aSeed).Should().Be("CUACCUC",
            "RC of the let-7a seed GAGGUAG is CUACCUC (§7.1, INV-04)");
        MiRnaAnalyzer.GetReverseComplement("").Should().BeEmpty("RC of empty is empty (§6.1)");
        MiRnaAnalyzer.GetReverseComplement(null!).Should().BeEmpty("RC of null is empty (§6.1)");
    }

    #endregion

    #region BE — Boundary: no complementarity (fully unpaired duplex)

    /// <summary>
    /// BE — "no complementarity": a target that pairs with NO miRNA position. The aligner
    /// pairs miRNA[i] antiparallel with target[len−1−i]; if no such pair is in the pairing
    /// set the whole alignment is spaces: Matches 0, GUWobbles 0, Mismatches = overlap,
    /// AlignmentString all spaces, and FreeEnergy = 0 (no consecutive paired stack to sum,
    /// INV-06 lower bound). We use a poly-A miRNA against a poly-A target — antiparallel A:A
    /// never pairs (A pairs only with U) — at several overlap widths, plus a hand-built
    /// G-vs-G case, and assert a fully unpaired, crash-free, defined result (§2.2, §6.1).
    /// </summary>
    [Test]
    public void AlignMiRnaToTarget_NoComplementarity_FullyUnpairedNeverCrash()
    {
        foreach (int n in new[] { 1, 4, 8, 22 })
        {
            string miRna = new string('A', n);
            string target = new string('A', n); // antiparallel A:A — never pairs
            var act = () => MiRnaAnalyzer.AlignMiRnaToTarget(miRna, target);
            var d = act.Should().NotThrow("a non-complementary duplex must not crash").Subject;

            AssertDuplexShape(d, n);
            d.Matches.Should().Be(0, "A:A is not a Watson-Crick pair (§2.2[3])");
            d.GUWobbles.Should().Be(0, "A:A is not a G:U wobble (§2.2[4])");
            d.Mismatches.Should().Be(n, "every antiparallel A:A position is a mismatch (INV-05)");
            d.AlignmentString.Should().Be(new string(' ', n), "a fully unpaired alignment is all spaces");
            d.FreeEnergy.Should().Be(0.0, "no paired stack contributes energy — ΔG = 0 (INV-06 lower bound)");
        }

        // Hand-built G-vs-G: antiparallel G:G never pairs either.
        var gg = MiRnaAnalyzer.AlignMiRnaToTarget("GGGG", "GGGG");
        AssertDuplexShape(gg, 4);
        gg.Mismatches.Should().Be(4, "antiparallel G:G is a mismatch at every position");
        gg.FreeEnergy.Should().Be(0.0, "no paired stack → ΔG = 0");
    }

    #endregion

    #region BE — Boundary: perfect match (fully Watson-Crick-paired duplex)

    /// <summary>
    /// BE — "perfect match": the doc's §7.1 worked example. "AAAA" vs "UUUU" — antiparallel
    /// A:U at all four positions → 4 Watson-Crick matches, 0 wobbles, 0 mismatches,
    /// AlignmentString "||||", and FreeEnergy ≤ 0 (a fully WC duplex is stabilising, INV-06).
    /// We assert that exact hand-derived result, then generalise: a target that is the
    /// antiparallel Watson-Crick complement of an arbitrary miRNA (its reverse complement)
    /// must be FULLY paired — Matches = length, 0 mismatches/wobbles, all-'|' alignment,
    /// ΔG ≤ 0. This pins the perfect-complementarity contract on a value derived from the
    /// spec, not from the code.
    /// </summary>
    [Test]
    public void AlignMiRnaToTarget_PerfectMatch_AllWatsonCrickStabilisingDuplex()
    {
        // Doc §7.1 worked example, derived independently: A:U is Watson-Crick at all 4 pos.
        var d = MiRnaAnalyzer.AlignMiRnaToTarget("AAAA", "UUUU");
        AssertDuplexShape(d, 4);
        d.Matches.Should().Be(4, "antiparallel A:U is Watson-Crick at all four positions (§7.1)");
        d.GUWobbles.Should().Be(0, "no G:U pair in an A/U duplex (§7.1)");
        d.Mismatches.Should().Be(0, "a perfectly complementary duplex has no mismatch (§7.1)");
        d.AlignmentString.Should().Be("||||", "a perfect Watson-Crick duplex aligns as all '|' (§7.1)");
        d.FreeEnergy.Should().BeLessThanOrEqualTo(0.0,
            "a fully Watson-Crick-paired duplex is stabilising — ΔG ≤ 0 (INV-06)");

        // Generalisation: the reverse complement of any miRNA is its perfect antiparallel
        // Watson-Crick target. RC(GAGGUAG)=CUACCUC; let-7a's RC target is fully paired.
        // (Note: RC complements with WC only — never introduces G:U — so wobbles stay 0.)
        foreach (string miRna in new[] { "GAGGUAG", "ACGCACGC", Let7aSequence })
        {
            string perfectTarget = MiRnaAnalyzer.GetReverseComplement(miRna);
            var p = MiRnaAnalyzer.AlignMiRnaToTarget(miRna, perfectTarget);
            int overlap = Math.Min(miRna.Length, perfectTarget.Length);

            AssertDuplexShape(p, overlap);
            p.Matches.Should().Be(overlap,
                $"the reverse complement of '{miRna}' pairs Watson-Crick at every position (INV-04)");
            p.GUWobbles.Should().Be(0, "reverse-complement pairing is pure Watson-Crick, no wobble");
            p.Mismatches.Should().Be(0, "a perfect antiparallel complement has no mismatch");
            p.AlignmentString.Should().Be(new string('|', overlap), "perfect complement → all '|'");
            p.FreeEnergy.Should().BeLessThanOrEqualTo(0.0,
                "a fully Watson-Crick duplex is stabilising — ΔG ≤ 0 (INV-06)");
        }
    }

    /// <summary>
    /// Positive sanity for the G:U wobble class, derived from the spec §7.1: "GGGG" vs
    /// "UUUU" — antiparallel G:U is a wobble at all four positions → GUWobbles 4, Matches 0,
    /// Mismatches 0, AlignmentString "::::" (wobble counted SEPARATELY from Watson-Crick,
    /// INV-02). This guards the matcher against silently folding wobble into match. We do
    /// NOT assert a sign on the wobble ΔG: the Turner wobble-stack table legitimately
    /// contains positive entries, so INV-06's ΔG ≤ 0 guarantee is scoped to Watson-Crick.
    /// </summary>
    [Test]
    public void AlignMiRnaToTarget_GuWobbleDuplex_CountedSeparatelyFromMatches()
    {
        var d = MiRnaAnalyzer.AlignMiRnaToTarget("GGGG", "UUUU");
        AssertDuplexShape(d, 4);
        d.GUWobbles.Should().Be(4, "antiparallel G:U is a wobble at all four positions (§7.1)");
        d.Matches.Should().Be(0, "a G:U wobble is NOT a Watson-Crick match (INV-02, §6.1)");
        d.Mismatches.Should().Be(0, "every position pairs (as a wobble) — no mismatch (§7.1)");
        d.AlignmentString.Should().Be("::::", "a fully wobble duplex aligns as all ':' (§7.1)");
    }

    #endregion

    #region BE — Boundary: short miRNA (overlap clamping, empty/null, no index overflow)

    /// <summary>
    /// BE — "short miRNA": a miRNA SHORTER than the target. The duplex overlaps only the
    /// shorter length (INV-05); the antiparallel index target[len−1−i] runs i in 0..min−1
    /// so it never walks off either string (the `targetIdx &lt; 0` guard + the min() overlap
    /// bound — §3.3, §4.1). We sweep miRNA lengths 0..8 against a fixed 22-nt target and
    /// assert: no throw, overlap == min length, well-formed counts. The empty / null miRNA
    /// (and empty / null target) corner returns the documented EMPTY MiRnaDuplex — all
    /// counts 0, empty strings — never an exception (§3.3, §6.1).
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void AlignMiRnaToTarget_ShortMiRna_OverlapClampedNeverIndexOverflow()
    {
        string longTarget = new string('U', 22); // a 22-nt target, longer than the short miRNAs

        for (int n = 0; n <= 8; n++)
        {
            string shortMiRna = new string('A', n); // poly-A → A:U Watson-Crick against poly-U
            var act = () => MiRnaAnalyzer.AlignMiRnaToTarget(shortMiRna, longTarget);
            var d = act.Should().NotThrow(
                $"a {n}-nt miRNA against a 22-nt target must not crash (overlap clamps to {n})").Subject;

            if (n == 0)
            {
                // Empty miRNA → the documented empty duplex (§3.3, §6.1).
                d.Matches.Should().Be(0);
                d.Mismatches.Should().Be(0);
                d.GUWobbles.Should().Be(0);
                d.Gaps.Should().Be(0);
                d.AlignmentString.Should().BeEmpty("an empty miRNA yields the empty duplex (§6.1)");
            }
            else
            {
                AssertDuplexShape(d, n);
                // poly-A miRNA vs poly-U target: every antiparallel A:U is Watson-Crick.
                d.Matches.Should().Be(n, $"all {n} antiparallel A:U positions are Watson-Crick");
                d.Mismatches.Should().Be(0, "poly-A vs poly-U is fully complementary");
            }
        }

        // Empty / null on either surface → the documented empty MiRnaDuplex, never a throw.
        foreach (var (m, t) in new[] { ("", longTarget), (longTarget, ""), ((string?)null, longTarget)!, (longTarget, (string?)null)! })
        {
            var act = () => MiRnaAnalyzer.AlignMiRnaToTarget(m!, t!);
            var d = act.Should().NotThrow("empty/null input returns the empty duplex, never an error (§3.3)").Subject;
            d.Matches.Should().Be(0);
            d.Mismatches.Should().Be(0);
            d.GUWobbles.Should().Be(0);
            d.Gaps.Should().Be(0);
            d.AlignmentString.Should().BeEmpty("the empty duplex has an empty alignment (§6.1)");
            d.FreeEnergy.Should().Be(0.0, "the empty duplex has zero free energy (§6.1)");
        }
    }

    #endregion

    #region Randomized boundary sweep — crash-free, invariant-true over fixed-seed random pairs

    /// <summary>
    /// Randomized boundary sweep: over fixed seeds and a range of miRNA/target lengths
    /// (including the unequal-length and short-miRNA boundaries), AlignMiRnaToTarget never
    /// throws and every reported duplex satisfies the universal invariants — overlap ==
    /// min(len), class counts partition the overlap (INV-05), Gaps 0, finite ΔG. The RNG is
    /// LOCALLY seeded (no shared static Rng) so every input is reproducible.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void AlignMiRnaToTarget_RandomPairs_AlwaysInvariantTrueNeverCrash()
    {
        foreach (int seed in new[] { 1, 7, 42, 2026 })
        {
            foreach (int mLen in new[] { 0, 1, 5, 8, 22, 40 })
            {
                foreach (int tLen in new[] { 0, 1, 6, 22, 35 })
                {
                    string miRna = RandomRna(mLen, seed);
                    string target = RandomRna(tLen, seed + 500);

                    var act = () => MiRnaAnalyzer.AlignMiRnaToTarget(miRna, target);
                    var d = act.Should().NotThrow(
                        $"random duplex must not crash (seed {seed}, mLen {mLen}, tLen {tLen})").Subject;

                    int overlap = (mLen == 0 || tLen == 0) ? 0 : Math.Min(mLen, tLen);
                    AssertDuplexShape(d, overlap);
                }
            }
        }
    }

    #endregion

    #endregion
}
