using System;
using System.Globalization;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;
using SignalPeptide = Seqeron.Genomics.Analysis.ProteinMotifFinder.SignalPeptide;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the ProteinMotif area — SIGNAL-PEPTIDE cleavage-site prediction
/// (PROTMOTIF-SP-001) via
/// <see cref="ProteinMotifFinder.PredictSignalPeptide(string, bool, double)"/>: the
/// von Heijne (1986) position-specific weight-matrix method, a faithful re-implementation
/// of the EMBOSS <c>sigcleave</c> reference program. Each candidate cleavage site is
/// scored by summing log-odds residue weights over a fixed 15-residue window
/// (positions −13..+2 relative to the cleavage site); the global argmax is reported as
/// the predicted mature-protein start, and a configurable weight threshold (default 3.5)
/// flags whether that site is a *likely* signal peptide.
///   — docs/algorithms/ProteinMotif/Signal_Peptide_Prediction.md §1, §2.2.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang or infinite loop, no
/// IndexOutOfRange / NullReference, no NaN/±∞ score, and no fabricated coordinate.
/// Every input must resolve to EITHER a well-defined, theory-correct result OR a
/// *documented, intentional* outcome (null). For a fixed-width weight-matrix scan that
/// slides a 15-column window across the sequence, the headline hazards are:
///   • a NullReferenceException when the sequence is null (the explicit
///     string.IsNullOrEmpty guard must short-circuit to null — §3.3, §6.1);
///   • an IndexOutOfRangeException on a SHORT sequence or at the window EDGES, where
///     window columns map to sequence offsets i−13..i+1 that fall outside [0, n) — the
///     scan must skip out-of-range columns, never run off either end (§4.1 step 3);
///   • a wrong null-guard BOUNDARY: 14 residues → null, exactly 15 residues → non-null
///     (one full window scored) — the documented min length is one full window
///     (§3.3, §6.1, §5.4 ASM-1);
///   • a NaN/±∞ Score from the log-odds arithmetic on an all-hydrophobic homopolymer or
///     a sequence dominated by non-standard residues;
///   • a CleavagePosition out of range, or a SignalSequence / WindowSequence that does
///     not match the position it claims (a coordinate bug).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROTMOTIF-SP-001 — signal-peptide cleavage-site prediction (von Heijne matrix)
/// Checklist: docs/checklists/03_FUZZING.md, row 167.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the length / composition corners that could crash,
///     hang, NaN out, or fabricate a site. Targets (row 167):
///       – "no signal": a globular protein with no signal-peptide architecture → the
///         best site still scores BELOW the 3.5 threshold, so the result is non-null
///         (the matrix always returns a best site for an in-window sequence — §6.2,
///         INV-05) but IsLikelySignalPeptide is false: NO false positive.
///       – "very short": a protein SHORTER than one full 15-residue window → null per
///         the documented guard (§3.3, §6.1, §5.4); no IndexOutOfRange on the window.
///         The exact boundary (14 → null, 15 → non-null) is pinned.
///       – "all-hydrophobic": an all-Leu / all-Ile homopolymer resembles a signal
///         peptide's hydrophobic h-region, yet the conserved small (−3,−1) residues are
///         absent (driving the 1e-10 penalty columns); the documented Score must be
///         FINITE and the flag well-defined, the window scoring must not run off the
///         ends, and the score must not be ±∞ despite the log-odds penalties.
/// — docs/checklists/03_FUZZING.md §Description (strategy code BE).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The signal-peptide contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Given a protein S (≥ 15 residues, case-insensitive), report the best-scoring cleavage
/// site as a SignalPeptide { CleavagePosition (1-based mature start), Score (von Heijne
/// weight), SignalSequence (residues 1..CleavagePosition−1), WindowSequence (the −13..+2
/// window, up to 15 residues), IsLikelySignalPeptide (Score ≥ minWeight) }.
///   • S(i) = Σ over 15 window columns of W(residue, column), where
///     W(a,p) = ln(C(a,p)/E(a)) with the 1e-10 pseudocount at columns −3/−1 and 1.0
///     elsewhere; out-of-range columns and non-standard residues contribute 0 (INV-01,
///     §2.2, §3.3).
///   • Returned site = global argmax of S(i) (INV-02).
///   • Cleavage is between CleavagePosition−1 and CleavagePosition, 1-based (INV-03).
///   • Output is case-independent (INV-04).
///   • IsLikelySignalPeptide ⇔ Score ≥ minWeight, default 3.5 (INV-05).
///   • null when the sequence is null/empty or shorter than 15 residues (§3.3, §6.1).
///   — docs/algorithms/ProteinMotif/Signal_Peptide_Prediction.md §2.2, §2.4, §3.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ProteinSignalPeptideFuzzTests
{
    #region Helpers

    // The 20 standard residues recognised by the von Heijne matrix (Signal_Peptide_Prediction.md
    // §3.3; ProteinMotifFinder MatrixResidues). Anything else contributes 0 to the score.
    private const string StandardResidues = "ACDEFGHIKLMNPQRSTVWY";

    // One full scoring window. Below this length PredictSignalPeptide returns null (§6.1).
    private const int WindowWidth = 15;

    // The EMBOSS sigcleave worked example: UniProt P17644 (ACH2_DROME), 576 aa.
    // Maximum score 13.739 at mature-protein start residue 42 — the canonical positive
    // signal peptide (Signal_Peptide_Prediction.md §7.1).
    private const string Ach2Drome =
        "MAPGCCTTRPRPIALLAHIWRHCKPLCLLLVLLLLCETVQANPDAKRLYDDLLSNYNRLI" +
        "RPVSNNTDTVLVKLGLRLSQLIDLNLKDQILTTNVWLEHEWQDHKFKWDPSEYGGVTELY" +
        "VPSEHIWLPDIVLYNNADGEYVVTTMTKAILHYTGKVVWTPPAIFKSSCEIDVRYFPFDQ" +
        "QTCFMKFGSWTYDGDQIDLKHISQKNDKDNKVEIGIDLREYYPSVEWDILGVPAERHEKY" +
        "YPCCAEPYPDIFFNITLRRKTLFYTVNLIIPCVGISYLSVLVFYLPADSGEKIALCISIL" +
        "LSQTMFFLLISEIIPSTSLALPLLGKYLLFTMLLVGLSVVITIIILNIHYRKPSTHKMRP" +
        "WIRSFFIKRLPKLLLMRVPKDLLRDLAANKINYGLKFSKTKFGQALMDEMQMNSGGSSPD" +
        "SLRRMQGRVGAGGCNGMHVTTATNRFSGLVGALGGGLSTLSGYNGLPSVLSGLDDSLSDV" +
        "AARKKYPFELEKAIHNVMFIQHHMQRQDEFNAEDQDWGFVAMVMDRLFLWLFMIASLVGT" +
        "FVILGEAPSLYDDTKAIDVQLSDVAKQIYNLTEKKN";

    private static string RandomStandardProtein(int length, int seed)
    {
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = StandardResidues[rng.Next(StandardResidues.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Asserts the documented well-formedness invariants of a non-null result, independent of
    /// whether the input is a "real" signal peptide (Signal_Peptide_Prediction.md §2.4, §3.2):
    /// Score finite (never NaN/±∞), CleavagePosition in 1..length, the SignalSequence equal to
    /// the upper-cased residues 1..CleavagePosition−1, and the WindowSequence equal to the
    /// upper-cased −13..+2 window clipped to the sequence.
    /// </summary>
    private static void AssertWellFormed(SignalPeptide sp, string originalSequence)
    {
        string upper = originalSequence.ToUpperInvariant();
        int len = upper.Length;

        double.IsNaN(sp.Score).Should().BeFalse("the von Heijne score must never be NaN (INV-01)");
        double.IsInfinity(sp.Score).Should().BeFalse(
            "the log-odds sum over 15 finite-weight columns must never be ±∞ (INV-01)");

        sp.CleavagePosition.Should().BeInRange(1, len,
            "the 1-based mature start must lie within the sequence (INV-03, §3.2)");

        int bestSite = sp.CleavagePosition - 1; // 0-based +1 residue
        sp.SignalSequence.Should().Be(upper.Substring(0, bestSite),
            "SignalSequence = residues 1..CleavagePosition−1 of the upper-cased input (§3.2, INV-04)");

        int windowStart = Math.Max(0, bestSite - 13);          // FirstColumnOffset = −13
        int windowEnd = Math.Min(len - 1, bestSite + 1);       // position +2, inclusive
        sp.WindowSequence.Should().Be(upper.Substring(windowStart, windowEnd - windowStart + 1),
            "WindowSequence = the −13..+2 scoring window clipped to the sequence (§3.2)");
        sp.WindowSequence.Length.Should().BeInRange(1, WindowWidth,
            "the window is at most 15 columns wide (§3.2)");

        sp.IsLikelySignalPeptide.Should().Be(sp.Score >= ProteinMotifFinder.DefaultMinWeight,
            "with the default minWeight the flag must equal Score ≥ 3.5 (INV-05)");
    }

    #endregion

    #region PROTMOTIF-SP-001 — signal peptide cleavage-site prediction

    #region BE — Null / empty / very short: documented null, no NullReference, no IndexOutOfRange

    // The explicit string.IsNullOrEmpty guard must short-circuit before any window indexing.
    [Test]
    public void PredictSignalPeptide_NullOrEmpty_ReturnsNullNoThrow()
    {
        PredictSignalPeptide(null!).Should().BeNull("null sequence → null (§3.3, §6.1)");
        PredictSignalPeptide(string.Empty).Should().BeNull("empty sequence → null (§3.3, §6.1)");
    }

    // "very short" target: every length 1..14 must return null (one full window is required),
    // and exactly 15 must return a non-null result — the documented null-guard boundary.
    // This is the IndexOutOfRange trap: window columns map to i−13..i+1, so a naive scan over
    // a short sequence would index out of range.
    [Test]
    public void PredictSignalPeptide_ShorterThanOneWindow_ReturnsNull()
    {
        for (int len = 1; len < WindowWidth; len++)
        {
            string seq = RandomStandardProtein(len, seed: 4200 + len);
            PredictSignalPeptide(seq).Should().BeNull(
                $"a {len}-residue sequence is shorter than one 15-residue window → null (§6.1, §5.4)");
        }
    }

    // Pin the exact boundary: 14 → null, 15 → non-null (no off-by-one in the guard).
    [Test]
    public void PredictSignalPeptide_NullGuardBoundary_Is14VersusExactly15()
    {
        string fourteen = RandomStandardProtein(14, seed: 991);
        string fifteen = fourteen + StandardResidues[0]; // extend to exactly one full window

        PredictSignalPeptide(fourteen).Should().BeNull("14 residues < one window → null");

        var sp = PredictSignalPeptide(fifteen);
        sp.Should().NotBeNull("exactly 15 residues fits one full window → non-null (§6.1)");
        AssertWellFormed(sp!.Value, fifteen);
    }

    #endregion

    #region BE — All-hydrophobic homopolymers: finite score, well-formed, no crash, no NaN/Inf

    // "all-hydrophobic" target: all-Leu / all-Ile / all-Val etc. A hydrophobic stretch
    // resembles a signal peptide's h-region, but the conserved small (−3,−1) residues are
    // absent, so those columns hit the 1e-10 penalty pseudocount. The score must stay FINITE
    // (a sum of finite log-odds), the result well-formed, and the window must not run off the
    // ends despite every position scoring identically.
    [Test]
    [CancelAfter(30000)]
    public void PredictSignalPeptide_AllHydrophobicHomopolymer_FiniteScoreWellFormed(CancellationToken token)
    {
        foreach (char residue in new[] { 'L', 'I', 'V', 'F', 'A', 'M', 'W' })
        {
            foreach (int len in new[] { WindowWidth, 15, 30, 100, 500 })
            {
                token.ThrowIfCancellationRequested();
                string seq = new string(residue, len);

                foreach (bool prokaryote in new[] { false, true })
                {
                    var sp = PredictSignalPeptide(seq, prokaryote);
                    sp.Should().NotBeNull(
                        $"an all-{residue} sequence of {len} residues fits ≥ one window → non-null (§6.2)");
                    AssertWellFormed(sp!.Value, seq);
                }
            }
        }
    }

    // Even an all-hydrophobic homopolymer is a heuristic best site; a Leu/Ile h-region without
    // the conserved small-residue cleavage motif typically does NOT clear the 3.5 threshold.
    // We do not over-constrain the exact score, only that the flag is consistent with it.
    [Test]
    public void PredictSignalPeptide_AllLeucine_FlagConsistentWithScore()
    {
        string allLeu = new string('L', 60);
        var sp = PredictSignalPeptide(allLeu);

        sp.Should().NotBeNull();
        AssertWellFormed(sp!.Value, allLeu);
        // The flag is purely Score ≥ minWeight; pin that contract rather than the heuristic verdict.
        sp.Value.IsLikelySignalPeptide.Should().Be(sp.Value.Score >= 3.5);
    }

    #endregion

    #region BE — "No signal" globular protein: best site returned but not flagged (no false positive)

    // A signal-peptide-free interior fragment must still return a best site (the matrix always
    // does for an in-window sequence — §6.2, INV-05), but it must NOT be flagged as a likely
    // signal peptide at the default threshold: no false positive.
    [Test]
    public void PredictSignalPeptide_GlobularInterior_NotFlaggedNoFalsePositive()
    {
        // A stretch of the ACH2_DROME MATURE region (well past the residue-42 cleavage site):
        // ordinary globular protein with no N-terminal signal architecture.
        string mature = Ach2Drome.Substring(120, 80);
        var sp = PredictSignalPeptide(mature);

        sp.Should().NotBeNull("any in-window sequence yields a best site (§6.2)");
        AssertWellFormed(sp!.Value, mature);
        sp.Value.IsLikelySignalPeptide.Should().BeFalse(
            "a mature-region fragment has no signal-peptide character → best score below 3.5, no false positive");
    }

    // Random standard-residue junk: never crashes, always well-formed, and deterministic across
    // repeated calls. A guard against runaway scans and coordinate fabrication.
    [Test]
    [CancelAfter(30000)]
    public void PredictSignalPeptide_RandomStandardJunk_AlwaysWellFormedAndDeterministic(CancellationToken token)
    {
        for (int seed = 0; seed < 200; seed++)
        {
            token.ThrowIfCancellationRequested();
            int len = 15 + (seed % 90);
            string seq = RandomStandardProtein(len, seed);

            var first = PredictSignalPeptide(seq);
            var second = PredictSignalPeptide(seq);

            first.Should().NotBeNull($"a {len}-residue sequence fits ≥ one window (seed {seed})");
            second.Should().Be(first, "scoring is pure / deterministic for a fixed input");
            AssertWellFormed(first!.Value, seq);
        }
    }

    #endregion

    #region BE — Non-standard residues: contribute 0, no exception (§3.3, §6.1)

    // Residues outside the 20-letter alphabet (X, B, Z, *, digits, punctuation) are not in the
    // matrix and must contribute 0 — never throw. An all-non-standard window must still produce
    // a finite, well-formed result (every column contributes 0 → Score 0).
    [Test]
    public void PredictSignalPeptide_NonStandardResidues_ContributeZeroNoThrow()
    {
        string allX = new string('X', 30);
        var spX = PredictSignalPeptide(allX);
        spX.Should().NotBeNull("an all-X sequence still fits a window; every column contributes 0 (§3.3)");
        AssertWellFormed(spX!.Value, allX);
        spX.Value.Score.Should().Be(0.0, "all 15 columns map to unmapped residues → Σ 0 (§3.3)");
        spX.Value.IsLikelySignalPeptide.Should().BeFalse("Score 0 < 3.5 → not flagged");

        // Mixed junk alphabet: must not throw and must stay well-formed.
        string junk = "XB*ZJUO12 #@-XB*ZJUO12 #@-XB*Z";
        var spJunk = PredictSignalPeptide(junk);
        spJunk.Should().NotBeNull();
        AssertWellFormed(spJunk!.Value, junk);
        double.IsNaN(spJunk.Value.Score).Should().BeFalse();
        double.IsInfinity(spJunk.Value.Score).Should().BeFalse();
    }

    #endregion

    #region BE — Case independence (INV-04)

    [Test]
    public void PredictSignalPeptide_CaseInsensitive_IdenticalResult()
    {
        var upper = PredictSignalPeptide(Ach2Drome.ToUpperInvariant());
        var lower = PredictSignalPeptide(Ach2Drome.ToLowerInvariant());

        upper.Should().NotBeNull();
        lower.Should().Be(upper, "input is upper-cased before scoring → case-independent (INV-04)");
    }

    #endregion

    #region Positive sanity — a genuine signal peptide is found, flagged, with a plausible cleavage site

    // The EMBOSS sigcleave worked example: ACH2_DROME scores 13.739 at mature-start residue 42,
    // well above 3.5 → flagged. This anchors the fuzz suite to the documented business contract:
    // a real signal peptide IS detected at the right place, so the BE "no false positive" tests
    // above are not vacuously green (§7.1, INV-02, INV-05).
    [Test]
    public void PredictSignalPeptide_KnownSignalPeptide_FlaggedAtDocumentedSite()
    {
        var sp = PredictSignalPeptide(Ach2Drome);

        sp.Should().NotBeNull("ACH2_DROME is a documented signal peptide");
        AssertWellFormed(sp!.Value, Ach2Drome);

        sp.Value.CleavagePosition.Should().Be(42,
            "EMBOSS sigcleave reports the mature protein starting at residue 42 (cleavage 41|42) (§7.1)");
        sp.Value.Score.Should().BeApproximately(13.7390400704164, 1e-9,
            "Σ log-odds over −13..+2 equals the EMBOSS maximum score 13.739 (INV-01, §7.1)");
        sp.Value.IsLikelySignalPeptide.Should().BeTrue("13.739 ≥ 3.5 → likely signal peptide (INV-05)");
        sp.Value.WindowSequence.Should().Be("LLVLLLLCETVQANP",
            "the 15-residue −13..+2 window at the best site (§7.1)");
        sp.Value.SignalSequence.Should().Be(Ach2Drome.Substring(0, 41),
            "the signal peptide is residues 1..41 (§3.2)");
    }

    #endregion

    #endregion
}
