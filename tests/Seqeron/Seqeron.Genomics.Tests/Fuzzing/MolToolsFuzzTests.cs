using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the MolTools area — CRISPR PAM (protospacer adjacent motif)
/// site finding (CRISPR-PAM-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds malformed, out-of-domain and boundary inputs to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang or infinite
/// loop, no state corruption, no nonsense output (spurious or out-of-bounds
/// sites), and no *unhandled* runtime exception (IndexOutOfRangeException,
/// NullReferenceException, ArgumentOutOfRangeException leaking from internal
/// indexing). Every input must resolve to EITHER a well-defined, theory-correct
/// result, OR a *documented, intentional* validation exception
/// (ArgumentException / ArgumentNullException). A raw runtime exception, a hang,
/// or a blow-up of bogus "PAM sites" on garbage input is a bug, not a passing
/// test. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CRISPR-PAM-001 — PAM site finding
/// Checklist: docs/checklists/03_FUZZING.md, row 18.
/// Fuzz strategies exercised for THIS unit:
///   • MC = Malformed Content — an invalid CRISPR system selector (an undefined
///          CrisprSystemType enum value, i.e. a "PAM specification" the API does
///          not know) and non-DNA junk characters in the sequence.
///   • BE = Boundary Exploitation — the empty sequence and the degenerate
///          "PAM/guide span longer than the sequence" case (no scan window fits).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes); row 18 targets:
///   "Invalid PAM sequences, non-DNA characters, empty seq, PAM longer than seq".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The PAM-site-finding contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A PAM is a short motif (e.g. SpCas9 NGG, where N = any base) immediately
/// adjacent to a CRISPR protospacer; Cas nucleases require it to bind and cleave
/// (PAM_Site_Detection.md §2.1). The detector scans BOTH strands for the
/// system-specific PAM motif (specified with IUPAC ambiguity codes), extracts the
/// guide-length target on the side the system dictates (PAM-after-target for Cas9,
/// PAM-before-target for Cas12a), and yields a PamSite ONLY when that target fits
/// within sequence bounds (PAM_Site_Detection.md §2.2, §2.4 INV-01..INV-03).
///
/// API entry: CrisprDesigner.FindPamSites(...)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs lines 40–60)
/// in two overloads — typed DnaSequence and raw string:
///   • FindPamSites(DnaSequence, CrisprSystemType = SpCas9): null → ArgumentNullException
///     (ThrowIfNull, line 44). The DnaSequence itself is a *validated* type whose
///     constructor rejects any non-A/C/G/T character with ArgumentException, so a
///     non-DNA sequence can never even reach this overload.
///   • FindPamSites(string, CrisprSystemType = SpCas9): null OR empty → the empty
///     result via `yield break` (lines 55–56); otherwise the input is upper-cased
///     and scanned RAW (no A/C/G/T validation). This is the surface that can be fed
///     malformed content directly, so it is the primary MC target here.
///
/// IMPORTANT — what "invalid PAM" means for THIS API (the central design fact the
/// checklist row probes): the PAM motif is NOT a free-form caller string. It is
/// selected by the CrisprSystemType enum and resolved by GetSystem(...) into one of
/// a fixed set of IUPAC-valid motifs (NGG, NAG, NNGRRT, TTTV, TTCN — §4.2). So an
/// "invalid PAM specification" is an UNDEFINED enum value, which GetSystem rejects
/// with ArgumentException (CrisprDesigner.cs line 27). Because every built-in PAM
/// motif contains only valid IUPAC codes, the IUPAC matcher's own
/// ArgumentOutOfRangeException guard (IupacHelper.MatchesIupac, the `_ => throw`
/// arm) is never reached from the PAM side — the motif is always well-formed; it is
/// the *sequence* base that is tested against the (valid) PAM code. We pin both
/// facets: an undefined system throws, and the IUPAC matcher rejects truly invalid
/// codes (so the validity of the fixed motifs is load-bearing, not incidental).
///
/// Non-DNA characters in the SEQUENCE (MC) — pinned on both surfaces:
///   • Typed surface: rejected up front by the DnaSequence constructor
///     (ArgumentException, "Invalid nucleotide") — junk never reaches the scanner.
///   • Raw-string surface: tolerated, NOT crashing. A junk character is tested as
///     the `nucleotide` argument of MatchesIupac(seqChar, pamChar); it matches no
///     IUPAC code, so it can never satisfy a PAM position and never invents a site.
///     The reverse-strand pass complements via GetComplementBase, whose final arm
///     passes any non-IUPAC char THROUGH unchanged (SequenceExtensions.cs line 156),
///     so no exception and no out-of-range indexing — at worst, fewer or no sites.
///     We pin that junk in the sequence yields no bogus sites and no crash, and that
///     interleaving junk into an otherwise-valid PAM window suppresses that site.
///
/// Empty sequence (BE): the raw-string overload short-circuits null/empty to the
/// empty enumerable (yield break); the typed overload over an empty DnaSequence has
/// a forward-scan bound `i ≤ len − pamLen` that is negative, so the loop never runs.
/// No PAM sites, no division, no indexing past the end (PAM_Site_Detection.md §6.1).
///
/// PAM(+guide span) longer than the sequence (BE): when the sequence is shorter than
/// the PAM pattern, the forward and reverse scan bounds `i ≤ len − pamLen` are
/// negative, so neither loop body runs and no Substring is taken past the end. Even
/// when the PAM motif itself fits but the guide-length target cannot (sequence
/// shorter than PAM+guide), the in-loop `targetStart >= 0 && targetEnd < seq.Length`
/// guard suppresses the site (INV-02). Either way: an empty result, never an
/// IndexOutOfRangeException (PAM_Site_Detection.md §6.1, "Sequence shorter than PAM
/// plus guide → no valid sites").
///
/// Documented invariants pinned on every positive site (PAM_Site_Detection.md §2.4):
/// INV-01 the matched PAM satisfies the system motif under IUPAC matching;
/// INV-02 the extracted target fits within sequence bounds; INV-03 both strands are
/// searched. The raw-string overload's scan body lives inside a `yield` iterator, so
/// every test forces enumeration (`.ToList()`) — the documented short-circuit and
/// any hang surface only on enumeration, and any hang would manifest as a
/// non-terminating materialization.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MolToolsFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CRISPR-PAM-001 — PAM site finding : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CRISPR-PAM-001 — PAM site finding

    #region MC — Invalid PAM specification (undefined CRISPR system)

    /// <summary>
    /// MC: an "invalid PAM" for this API is an undefined CRISPR system selector —
    /// the PAM motif is chosen by the CrisprSystemType enum, not a free-form string
    /// (CrisprDesigner.cs §PAM Definitions). An out-of-range enum value resolves to
    /// the `_ => throw new ArgumentException` arm of GetSystem (line 27), so the
    /// detector REJECTS the bad specification rather than scanning with a garbage
    /// motif or crashing. Pinned on both surfaces; the raw-string overload resolves
    /// the system inside its iterator, so we force enumeration to surface the throw.
    /// </summary>
    [Test]
    public void FindPamSites_UndefinedCrisprSystem_ThrowsArgumentException()
    {
        const CrisprSystemType bogus = (CrisprSystemType)9999;

        var typed = () => CrisprDesigner.FindPamSites(new DnaSequence("ACGTACGTACGTACGTACGTACGTAGG"), bogus).ToList();
        var raw = () => CrisprDesigner.FindPamSites("ACGTACGTACGTACGTACGTACGTAGG", bogus).ToList();

        typed.Should().Throw<ArgumentException>(
            "an undefined CrisprSystemType is an invalid PAM specification; GetSystem rejects it rather than scanning a garbage motif");
        raw.Should().Throw<ArgumentException>(
            "the raw-string overload resolves the same GetSystem mapping and rejects the undefined system on enumeration");
    }

    /// <summary>
    /// MC: the fixed PAM motifs (NGG, NAG, NNGRRT, TTTV, TTCN) are load-bearing —
    /// they are valid IUPAC patterns, which is WHY scanning never trips the IUPAC
    /// matcher's own validity guard. This pins that guard directly: a genuinely
    /// invalid IUPAC code IS rejected by IupacHelper.MatchesIupac (the `_ => throw`
    /// arm), proving the no-crash guarantee for PAM scanning rests on the motifs
    /// being well-formed, not on the matcher silently swallowing junk.
    /// </summary>
    [Test]
    public void MatchesIupac_InvalidPamCode_ThrowsArgumentOutOfRange()
    {
        var act = () => IupacHelper.MatchesIupac('A', 'Z');

        act.Should().Throw<ArgumentOutOfRangeException>(
            "'Z' is not one of the 15 IUPAC codes; the matcher rejects an invalid PAM code rather than matching arbitrarily");
    }

    #endregion

    #region MC — Non-DNA characters in the sequence

    /// <summary>
    /// MC: non-DNA junk in the sequence fed to the TYPED overload is rejected up
    /// front — the DnaSequence constructor validates A/C/G/T and throws
    /// ArgumentException on the first offending character (DnaSequence.cs
    /// ValidateSequence), so junk never reaches the PAM scanner at all.
    /// </summary>
    [Test]
    public void FindPamSites_NonDnaSequence_TypedOverload_RejectedAtConstruction()
    {
        var act = () => new DnaSequence("ACGT$#@!XYZ123");

        act.Should().Throw<ArgumentException>(
            "the validated DnaSequence type rejects any non-A/C/G/T character, so the typed PAM scanner only ever sees clean DNA");
    }

    /// <summary>
    /// MC: non-DNA junk fed to the RAW-STRING overload must NOT crash and must NOT
    /// invent sites. Each sequence character is tested as the `nucleotide` argument
    /// of MatchesIupac(seqChar, pamChar); a junk char matches no IUPAC code, so it
    /// can never satisfy a PAM position. The reverse-strand pass complements via
    /// GetComplementBase, whose fall-through arm passes non-IUPAC chars through
    /// unchanged — no exception, no out-of-range indexing. Pure-junk input yields no
    /// PAM sites. We force enumeration so the in-iterator scan actually runs.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPamSites_NonDnaSequence_RawOverload_NoSitesNoCrash()
    {
        var act = () => CrisprDesigner.FindPamSites("$#@!XYZ123 ￿qwerty", CrisprSystemType.SpCas9).ToList();

        act.Should().NotThrow(
            "non-DNA characters are tested against the PAM via IUPAC matching and complemented via a pass-through arm; neither path indexes out of range");

        CrisprDesigner.FindPamSites("$#@!XYZ123 ￿qwerty", CrisprSystemType.SpCas9)
            .Should().BeEmpty(
                "a character that is not A/C/G/T matches no PAM position, so pure junk can never produce a spurious PAM site");
    }

    /// <summary>
    /// MC: junk INTERLEAVED into an otherwise-valid PAM window must suppress that
    /// site, not crash. "…AGG" at the right offset is a real forward SpCas9 site;
    /// replacing the final 'G' with '#' breaks the NGG match (the '#' satisfies no
    /// IUPAC code), so no forward site is reported there — proving junk neither
    /// fabricates nor silently "rounds to" a valid PAM.
    /// </summary>
    [Test]
    public void FindPamSites_JunkInsidePamWindow_SuppressesThatSite()
    {
        // 20 A's (guide) + "AG#" — the NGG window at index 20 is broken by '#'.
        string broken = new string('A', 20) + "AG#";

        var sites = CrisprDesigner.FindPamSites(broken, CrisprSystemType.SpCas9).ToList();

        sites.Should().NotContain(s => s.IsForwardStrand && s.Position == 20,
            "the '#' in the third PAM position fails IUPAC matching, so the forward NGG site at index 20 is correctly NOT reported");
        sites.Should().OnlyContain(s => s.PamSequence.All(c => c == 'A' || c == 'C' || c == 'G' || c == 'T'),
            "any site that IS reported still has a clean A/C/G/T PAM string — junk never leaks into a reported PAM");
    }

    #endregion

    #region BE — Empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. The raw-string overload
    /// short-circuits null/empty to the empty enumerable (yield break, lines 55–56);
    /// the typed overload over an empty DnaSequence has a negative forward-scan bound
    /// `i ≤ 0 − pamLen`, so the loop never runs. Neither path divides, indexes past
    /// the end, or hangs (PAM_Site_Detection.md §6.1). Pinned for both surfaces plus
    /// the raw null input.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPamSites_EmptySequence_IsEmptyAndDoesNotThrow()
    {
        var typed = () => CrisprDesigner.FindPamSites(new DnaSequence(string.Empty), CrisprSystemType.SpCas9).ToList();
        var rawEmpty = () => CrisprDesigner.FindPamSites(string.Empty, CrisprSystemType.SpCas9).ToList();
        var rawNull = () => CrisprDesigner.FindPamSites((string)null!, CrisprSystemType.SpCas9).ToList();

        typed.Should().NotThrow("an empty sequence has no scan window; the forward bound is negative so the loop never runs");
        rawEmpty.Should().NotThrow("the raw-string overload short-circuits empty input to an empty result");
        rawNull.Should().NotThrow("the raw-string overload treats null input as empty, not as an error");

        CrisprDesigner.FindPamSites(new DnaSequence(string.Empty), CrisprSystemType.SpCas9).Should().BeEmpty();
        CrisprDesigner.FindPamSites(string.Empty, CrisprSystemType.SpCas9).Should().BeEmpty();
        CrisprDesigner.FindPamSites((string)null!, CrisprSystemType.SpCas9).Should().BeEmpty();
    }

    /// <summary>
    /// BE/INJ: a null DnaSequence is the boundary of "no typed input". The typed
    /// overload guards it with ArgumentNullException (ThrowIfNull, line 44), raised
    /// eagerly at the call — never a NullReferenceException dereferencing
    /// `sequence.Sequence` (PAM_Site_Detection.md §6.1).
    /// </summary>
    [Test]
    public void FindPamSites_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => CrisprDesigner.FindPamSites((DnaSequence)null!, CrisprSystemType.SpCas9);

        act.Should().Throw<ArgumentNullException>(
            "the typed overload null-guards its sequence; null is rejected, never dereferenced into a NullReferenceException");
    }

    #endregion

    #region BE — PAM (and guide span) longer than the sequence

    /// <summary>
    /// BE: a sequence SHORTER than the PAM pattern itself is the degenerate "PAM
    /// longer than seq" case. SaCas9's PAM is 6 nt (NNGRRT); on a 3-nt sequence the
    /// scan bounds `i ≤ len − pamLen` are negative on both strands, so neither loop
    /// runs and no Substring is taken past the end — an empty result, never an
    /// IndexOutOfRangeException (PAM_Site_Detection.md §6.1). Pinned on both surfaces.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPamSites_PamLongerThanSequence_IsEmptyAndDoesNotThrow()
    {
        var typed = () => CrisprDesigner.FindPamSites(new DnaSequence("ACG"), CrisprSystemType.SaCas9).ToList();
        var raw = () => CrisprDesigner.FindPamSites("ACG", CrisprSystemType.SaCas9).ToList();

        typed.Should().NotThrow(
            "the 6-nt NNGRRT PAM cannot fit in a 3-nt sequence; the scan bound is negative so no Substring is taken past the end");
        raw.Should().NotThrow("the raw-string overload is equally guarded against indexing past the sequence end");

        CrisprDesigner.FindPamSites(new DnaSequence("ACG"), CrisprSystemType.SaCas9).Should().BeEmpty(
            "a PAM longer than the whole sequence yields no sites, not a crash");
        CrisprDesigner.FindPamSites("ACG", CrisprSystemType.SaCas9).Should().BeEmpty();
    }

    /// <summary>
    /// BE: the PAM motif FITS but the guide-length target does not — the sequence is
    /// long enough to hold an NGG but too short for the 20-nt protospacer upstream of
    /// it. The in-loop bounds guard `targetStart >= 0 && targetEnd < seq.Length`
    /// (line 87) suppresses the site (INV-02) rather than slicing a negative-start
    /// Substring. "AGG" alone matches NGG at index 0 but targetStart = 0 − 20 = −20,
    /// so no site is yielded. Pinned: no crash, empty result.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPamSites_PamFitsButGuideDoesNot_IsEmptyAndDoesNotThrow()
    {
        var typed = () => CrisprDesigner.FindPamSites(new DnaSequence("AGG"), CrisprSystemType.SpCas9).ToList();
        var raw = () => CrisprDesigner.FindPamSites("AGG", CrisprSystemType.SpCas9).ToList();

        typed.Should().NotThrow(
            "the NGG matches but the 20-nt target would start at -20; the bounds guard suppresses it instead of slicing past the start");
        raw.Should().NotThrow();

        CrisprDesigner.FindPamSites(new DnaSequence("AGG"), CrisprSystemType.SpCas9).Should().BeEmpty(
            "no protospacer fits upstream of an NGG at index 0, so the site is correctly not yielded");
        CrisprDesigner.FindPamSites("AGG", CrisprSystemType.SpCas9).Should().BeEmpty();
    }

    #endregion

    #region Positive sanity — known PAM sites are found at the correct positions

    /// <summary>
    /// Positive sanity: alongside the degenerate probes, a textbook SpCas9 forward
    /// site must be detected at the CORRECT position with the CORRECT PAM and target,
    /// so the boundary hardening never silently breaks the core function. The
    /// sequence is 20 A's (the protospacer) followed by "AGG": the only GG pair sits
    /// at indices 21–22, so NGG matches at PAM position 20 (N = seq[20] = 'A'). The
    /// forward site therefore has Position 20, PamSequence "AGG", and a 20-nt all-A
    /// target starting at index 0 (INV-01, INV-02).
    /// </summary>
    [Test]
    public void FindPamSites_KnownForwardNggSite_DetectedAtCorrectPosition()
    {
        string seq = new string('A', 20) + "AGG"; // length 23

        var forward = CrisprDesigner.FindPamSites(seq, CrisprSystemType.SpCas9)
            .Where(s => s.IsForwardStrand)
            .ToList();

        var site = forward.Should().ContainSingle(s => s.Position == 20).Subject;
        site.PamSequence.Should().Be("AGG", "the NGG window at index 20 is 'AGG' (N = the upstream 'A')");
        site.PamSequence[1].Should().Be('G');
        site.PamSequence[2].Should().Be('G', "INV-01: the matched PAM satisfies the NGG motif under IUPAC matching");
        site.TargetSequence.Should().Be(new string('A', 20), "the 20-nt protospacer immediately upstream of the PAM is all-A");
        site.TargetStart.Should().Be(0, "INV-02: the target fits within bounds, starting at the first base");
        site.System.PamSequence.Should().Be("NGG");
    }

    /// <summary>
    /// Positive sanity / INV-03: both strands are searched. A forward NGG on one
    /// strand appears as a reverse-complement CCN motif on the other; a construct
    /// carrying a valid reverse-strand site must therefore produce at least one site
    /// with IsForwardStrand == false. "CCT" + 20 A's contains "CC" at the 5' end,
    /// which is the reverse complement of an "GG" PAM read on the opposite strand, so
    /// the reverse-strand scan must report a non-forward site. Pinned: the detector
    /// does not silently search only the forward strand.
    /// </summary>
    [Test]
    public void FindPamSites_ReverseStrandSite_IsAlsoSearched()
    {
        // Reverse complement of ("CCT" + 20×A) is (20×T + "AGG"): a forward NGG on the
        // reverse-complement strand → a reverse-strand site on the original.
        string seq = "CCT" + new string('A', 20);

        var sites = CrisprDesigner.FindPamSites(seq, CrisprSystemType.SpCas9).ToList();

        sites.Should().Contain(s => !s.IsForwardStrand,
            "INV-03: the reverse strand is scanned, so a reverse-complement (CCN) PAM yields a non-forward site");
        sites.Should().OnlyContain(s =>
                s.TargetStart >= 0 &&
                s.PamSequence.Length == 3,
            "every yielded site is well-formed: a 3-nt PAM and a non-negative target start");
    }

    /// <summary>
    /// Positive sanity: a Cas12a (Cpf1) site exercises the PAM-BEFORE-target branch
    /// (the opposite orientation from Cas9), so the boundary work does not break the
    /// alternate code path. Cas12a PAM is TTTV (V = A/C/G), guide length 23, PAM 5'
    /// of the target. "TTTA" + 23 valid bases places a TTTV at index 0 with a 23-nt
    /// target immediately downstream. We pin that at least one Cas12a site is found
    /// with the PAM upstream of (before) its target start.
    /// </summary>
    [Test]
    public void FindPamSites_Cas12aPamBeforeTarget_DetectedWithCorrectOrientation()
    {
        string seq = "TTTA" + new string('C', 23); // TTTV at 0, 23-nt target downstream

        var forward = CrisprDesigner.FindPamSites(seq, CrisprSystemType.Cas12a)
            .Where(s => s.IsForwardStrand)
            .ToList();

        var site = forward.Should().ContainSingle(s => s.Position == 0).Subject;
        site.PamSequence.Should().Be("TTTA", "TTTV matches 'TTTA' (V = A) at the 5' end");
        site.System.PamAfterTarget.Should().BeFalse("Cas12a is a PAM-before-target system");
        site.TargetStart.Should().Be(4, "the 23-nt target lies immediately DOWNSTREAM of the 4-nt PAM");
        site.TargetSequence.Should().Be(new string('C', 23), "the protospacer is the 23 C's following the PAM");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and
    /// produce only well-formed sites — no out-of-bounds targets, no malformed PAM
    /// strings, no hang — so the degenerate-boundary guards do not corrupt the scan
    /// on ordinary input. Every site must satisfy INV-01 (PAM matches NGG under
    /// IUPAC) and INV-02 (target within bounds) regardless of random content.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void FindPamSites_RandomSequence_ProducesOnlyWellFormedSites()
    {
        string seq = RandomDna(2000, seed: 18_001);

        var sites = CrisprDesigner.FindPamSites(seq, CrisprSystemType.SpCas9).ToList();

        sites.Should().NotBeEmpty("a 2000-nt random sequence is overwhelmingly likely to contain NGG sites with room for a guide");
        // INV-01 is strand-aware: a forward site's PamSequence matches NGG (pos1=G, pos2=G);
        // a reverse site's PamSequence is reverse-complemented back to forward orientation, so it
        // reads as CCN (pos0=C, pos1=C) — the reverse complement of an NGG read on the other strand.
        sites.Should().OnlyContain(s =>
                s.PamSequence.Length == 3 &&
                (s.IsForwardStrand
                    ? (s.PamSequence[1] == 'G' && s.PamSequence[2] == 'G')
                    : (s.PamSequence[0] == 'C' && s.PamSequence[1] == 'C')) &&  // INV-01
                s.TargetSequence.Length == 20 &&                                // guide length
                s.TargetStart >= 0,                                             // INV-02
            "every site on random input is well-formed: an NGG PAM (forward) or its CCN reverse-complement (reverse) with a 20-nt in-bounds target");
    }

    #endregion

    #endregion
}
