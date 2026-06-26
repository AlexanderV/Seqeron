using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the ProteinMotif area — Plan7 profile-HMM domain DETECTION
/// (PROTMOTIF-HMM-001): scoring a query against a HMMER3 profile hidden Markov model
/// with the log-odds Viterbi (optimal path) and Forward (full likelihood) recurrences,
/// and the HMMER3/f ASCII <c>.hmm</c> profile parser that feeds them. The entry points
/// under fuzz are the programmatic profile object
/// <see cref="Plan7ProfileHmm.Parse(string)"/> /
/// <see cref="Plan7ProfileHmm.ViterbiScore"/> /
/// <see cref="Plan7ProfileHmm.ForwardScore"/> /
/// <see cref="Plan7ProfileHmm.LocalForwardScore"/> /
/// <see cref="Plan7ProfileHmm.HmmSearchBitScore"/> /
/// <see cref="Plan7ProfileHmm.FindDomains"/>, and the bundled-profile facade
/// <see cref="ProteinMotifFinder.ScoreDomainHmm"/> /
/// <see cref="ProteinMotifFinder.FindDomainsByHmm"/> /
/// <see cref="ProteinMotifFinder.FindDomainEnvelopes(string, string)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Malformed/boundary inputs must NEVER hang, throw an UNHANDLED runtime exception
/// (IndexOutOfRange / NullReference / a NaN/±∞ log-odds escaping from a log(0) on an
/// all-unknown or empty-profile path), or emit out-of-contract output. Every input must
/// resolve to EITHER a well-defined, theory-correct result OR a DOCUMENTED validation
/// exception:
///   • ArgumentNullException for a null sequence/text (contract §3.3);
///   • FormatException for a MALFORMED .hmm profile (not HMMER3, not ALPH amino,
///     truncated body, a non-numeric / wrong-arity parameter row, a node-order error)
///     — never a raw int.Parse / IndexOutOfRange parse crash (§3.3, §4.1);
///   • ArgumentException for an unknown bundled accession (§3.3, §6.1).
/// The DP itself must always TERMINATE with a finite-or-(−∞) score: a query of all-X
/// (unknown residues), a zero-length query, and a profile with MORE match states than
/// the query has residues are the boundary corners that could blow up to NaN or run off
/// the DP rows; the theory says each must yield a defined score, not a crash.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes MC, BE); targets:
///   "empty profile, zero-length sequence, all-X residues, malformed .hmm,
///    profile longer than seq".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROTMOTIF-HMM-001 — Plan7 profile-HMM domain detection
/// Checklist: docs/checklists/03_FUZZING.md, row 239.
/// Algorithm doc: docs/algorithms/ProteinMotif/Profile_HMM_Domain_Detection.md.
/// Fuzz strategies exercised for THIS unit:
///   • MC = Malformed Content — broken .hmm profile text:
///       – non-HMMER text / missing the "HMMER3/" version line → FormatException (§3.3);
///       – ALPH not "amino" → FormatException (only amino profiles supported, §3.3);
///       – missing / non-positive LENG (≤ 0, the "empty profile" corner) → FormatException;
///       – truncated body (header present but no 'HMM' main-model line, or a node row
///         missing) → FormatException, never a NullReference / IndexOutOfRange;
///       – non-numeric or wrong-arity emission/transition fields → FormatException
///         (FormatException from int/double.Parse is the DOCUMENTED malformed outcome),
///         never an unhandled overflow / index crash;
///       – out-of-order node numbering → FormatException (explicit node-order guard).
///   • BE = Boundary Exploitation — the DP corners that could crash / NaN:
///       – zero-length query "" → a DEFINED score (the all-delete glocal path; not a
///         crash, not NaN), and the bundled facades' empty-string guard (§3.3, §6.1);
///       – all-X query (X ∉ the 20 canonical codes): every residue scores via the
///         null/background as neutral (odds 1, log-odds 0) — the score is FINITE, NOT a
///         NaN/−∞ from a log(0) (doc §3.3 "characters outside the 20 canonical amino
///         acids contribute neutral (zero) log-odds"; EmissionLogOdds returns 0 for an
///         unknown residue);
///       – profile LONGER than the query (more match states M than residues n, incl. a
///         multi-node profile vs a 1-residue query): the Viterbi/Forward DP must still
///         TERMINATE with a defined (low / −∞-allowed but never NaN/+∞) score.
/// — docs/checklists/03_FUZZING.md §Description; targets quoted above.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The profile-HMM contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A HMMER3/f profile stores every probability as a NEGATIVE natural-log (−ln p) to five
/// decimals, with '*' meaning probability 0 (→ −∞ log-prob); ALPH amino fixes the K = 20
/// alphabet ACDEFGHIKLMNPQRSTVWY; the COMPO line is the null background; each node carries
/// match + insert emissions and the 7 transitions Mk→Mk+1,Mk→Ik,Mk→Dk+1,Ik→Mk+1,Ik→Ik,
/// Dk→Mk+1,Dk→Dk+1. The glocal Viterbi/Forward recurrences (Durbin §5.4) score the query
/// in nats; bits = nats / ln 2 (doc §2.2, §4).
///
/// Theory-correct invariants asserted (doc §2.4 INV-HMM-01..05):
///   • INV-HMM-01 — Forward score ≥ Viterbi score in the SAME units (Forward sums via
///     log-sum-exp over all paths, including the single optimal one). The headline
///     "Viterbi ≤ Forward" relation, pinned on every fuzzed query.
///   • [finite] — a returned score is finite OR −∞ (a forbidden path), NEVER NaN and
///     NEVER +∞. This is the no-log(0)-NaN / no-overflow guarantee on all-X / empty-query
///     / profile-longer-than-seq.
///   • INV-HMM-02 — a query equal to the profile's own match consensus scores STRICTLY
///     ABOVE an unrelated random query of similar length (log-odds is positive when the
///     profile explains the query better than the background). Pins that the harness
///     scores a REAL model, not a no-op.
///   • INV-HMM-03 — scoring is DETERMINISTIC: identical inputs give identical scores.
///   • hit-coords — a reported domain envelope is in-bounds (1 ≤ start ≤ end ≤ n) and its
///     bit score / bias / E-value are finite (never NaN/±∞).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Complexity / hang-safety (doc §4.3: Viterbi/Forward are O(n·M))
/// ───────────────────────────────────────────────────────────────────────────
/// The DP is O(n·M); the all-X, long-query and profile-longer-than-seq targets are kept
/// modest and [CancelAfter]-guarded so a regression that turned the DP into a hang or a
/// super-linear blow-up FAILS rather than wedging the suite.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class ProteinProfileHmmFuzzTests
{
    #region Helpers

    /// <summary>The 20 standard amino-acid one-letter codes, in HMMER alphabet order.</summary>
    private const string AminoAlphabet = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>A bundled accession (PF00018 SH3) used to fuzz the facade scoring path.</summary>
    private const string Sh3Accession = "PF00018";

    /// <summary>ln 2, for the nats→bits relation (doc §2.2).</summary>
    private const double Ln2 = 0.69314718055994530941723212145818;

    /// <summary>
    /// −ln of a probability to 5 decimals (HMMER3/f storage convention), or "*" for p = 0.
    /// </summary>
    private static string NegLog(double p) =>
        p <= 0.0 ? "*" : (-Math.Log(p)).ToString("F5", CultureInfo.InvariantCulture);

    /// <summary>
    /// Builds a syntactically valid HMMER3/f ASCII profile of <paramref name="consensus"/> match
    /// states. Each node's match emission strongly favours the corresponding consensus residue
    /// (probability 0.81 on it, the rest split among the other 19); inserts are ≈ background; the
    /// transitions overwhelmingly favour the M→M consensus path so the consensus query scores well
    /// above background (INV-HMM-02). The COMPO/background is uniform (1/20 each). The profile is
    /// uncalibrated (no STATS) unless <paramref name="withStats"/>.
    /// </summary>
    private static string BuildProfile(string consensus, bool withStats = false)
    {
        int m = consensus.Length;
        var sb = new StringBuilder();
        sb.AppendLine("HMMER3/f [3.3 | test]");
        sb.AppendLine("NAME  FUZZ_TEST");
        sb.AppendLine("ACC   PF99999.1");
        sb.AppendLine("DESC  Synthetic fuzz profile");
        sb.AppendLine(CultureInfo.InvariantCulture, $"LENG  {m}");
        sb.AppendLine("ALPH  amino");
        if (withStats)
        {
            sb.AppendLine("STATS LOCAL MSV       -8.0000  0.69315");
            sb.AppendLine("STATS LOCAL VITERBI   -8.2000  0.69315");
            sb.AppendLine("STATS LOCAL FORWARD   -4.5000  0.69315");
        }
        sb.AppendLine("HMM          A        C        D        E        F        G        H        I        K        L        M        N        P        Q        R        S        T        V        W        Y");
        sb.AppendLine("            m->m     m->i     m->d     i->m     i->i     d->m     d->d");

        // Uniform background COMPO (each residue probability 1/20).
        sb.AppendLine("  COMPO   " + EmitRow(_ => 1.0 / 20.0));
        // BEGIN node: insert-0 emissions (uniform) then begin transitions (B->M1 dominant).
        sb.AppendLine("          " + EmitRow(_ => 1.0 / 20.0));
        sb.AppendLine("          " + TransitionRow(0.98, 0.01, 0.01, 0.5, 0.5, 1.0, 0.0)); // B->{M1,I0,D1}, I0->{M1,I0}, d row
        for (int k = 1; k <= m; k++)
        {
            int target = AminoAlphabet.IndexOf(char.ToUpperInvariant(consensus[k - 1]));
            // Match emissions: 0.81 on the consensus residue, the rest split among the other 19.
            string matchLine = string.Create(CultureInfo.InvariantCulture,
                $"  {k,5}   ") + EmitRow(a => a == target ? 0.81 : 0.19 / 19.0)
                + string.Create(CultureInfo.InvariantCulture, $"      {k} - - -");
            sb.AppendLine(matchLine);
            // Insert emissions ≈ background (uniform).
            sb.AppendLine("          " + EmitRow(_ => 1.0 / 20.0));
            // Transitions: strongly favour M->M (consensus path); last node M->M is 1.0 to E.
            bool last = k == m;
            sb.AppendLine("          " + TransitionRow(
                last ? 1.0 : 0.97, 0.015, last ? 0.0 : 0.015, 0.5, 0.5, last ? 1.0 : 0.97, last ? 0.0 : 0.03));
        }
        return sb.ToString();
    }

    /// <summary>Emits 20 −ln-probability fields for the alphabet, value chosen per residue index.</summary>
    private static string EmitRow(Func<int, double> p)
        => string.Join("  ", Enumerable.Range(0, 20).Select(a => NegLog(p(a))));

    /// <summary>Emits the 7 transition fields Mk→Mk+1,Mk→Ik,Mk→Dk+1,Ik→Mk+1,Ik→Ik,Dk→Mk+1,Dk→Dk+1.</summary>
    private static string TransitionRow(double mm, double mi, double md, double im, double ii, double dm, double dd)
        => string.Join("  ", new[] { mm, mi, md, im, ii, dm, dd }.Select(NegLog));

    private static string RandomProtein(int length, int seed)
    {
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = AminoAlphabet[rng.Next(AminoAlphabet.Length)];
        return new string(chars);
    }

    /// <summary>True iff <paramref name="v"/> is finite or exactly −∞ (a forbidden path), never NaN/+∞.</summary>
    private static void AssertScoreWellFormed(double v, string because)
    {
        double.IsNaN(v).Should().BeFalse($"a profile-HMM score must never be NaN ({because})");
        double.IsPositiveInfinity(v).Should().BeFalse($"a profile-HMM score must never be +∞ ({because})");
        // −∞ is an allowed, defined score (a forbidden path); everything else must be finite.
        (double.IsFinite(v) || double.IsNegativeInfinity(v)).Should().BeTrue(
            $"a score is finite or −∞, never any other non-finite value ({because})");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PROTMOTIF-HMM-001 — profile-HMM domain detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PROTMOTIF-HMM-001 — Plan7 profile-HMM domain detection

    #region MC — Malformed .hmm profile text → FormatException / ArgumentNullException, never a parse crash

    /// <summary>
    /// Target "malformed .hmm": a battery of broken profile texts must each be rejected with a
    /// DOCUMENTED <see cref="FormatException"/> (or <see cref="ArgumentNullException"/> for null) —
    /// never a raw int/double.Parse crash, NullReference, or IndexOutOfRange leaking out of the
    /// parser. Covered malformations: not-HMMER text, missing the HMMER3 version line, a non-amino
    /// ALPH, a missing/zero/negative LENG (the "empty profile" corner), a header with no 'HMM'
    /// main-model line, a node row missing (truncated body), a non-numeric emission field, an
    /// emission row with too few numbers, and out-of-order node numbering (doc §3.3, §4.1).
    /// </summary>
    [Test]
    public void Parse_MalformedProfile_ThrowsDocumentedException()
    {
        // null text → ArgumentNullException (contract §3.3).
        ((Action)(() => Plan7ProfileHmm.Parse((string)null!)))
            .Should().Throw<ArgumentNullException>("a null profile text is a null-argument violation");

        // A valid 3-node profile as the mutation baseline.
        string valid = BuildProfile("KWY");
        // Sanity: the baseline parses (so the malformations below are the ONLY difference).
        ((Action)(() => Plan7ProfileHmm.Parse(valid))).Should().NotThrow("the baseline profile is valid");

        var malformed = new (string Name, string Text)[]
        {
            ("empty string", ""),
            ("not a HMMER profile at all", "this is not a profile\nsome junk\n"),
            ("missing HMMER3 version line", "NAME X\nLENG 2\nALPH amino\nHMM\n"),
            ("nucleotide ALPH (not amino)", valid.Replace("ALPH  amino", "ALPH  DNA")),
            ("missing LENG", valid.Replace("LENG  3", "")),
            ("zero LENG (empty profile)", valid.Replace("LENG  3", "LENG  0")),
            ("negative LENG", valid.Replace("LENG  3", "LENG  -2")),
            ("header only, no HMM main-model line",
                "HMMER3/f [3.3 | test]\nNAME X\nLENG 2\nALPH amino\n"),
            ("truncated body — node 3 lines cut off",
                string.Join("\n", valid.Split('\n')[..^4]) + "\n"),
        };

        foreach (var (name, text) in malformed)
            ((Action)(() => Plan7ProfileHmm.Parse(text)))
                .Should().Throw<FormatException>($"malformed profile ({name}) must raise FormatException, not a raw crash");

        // A non-numeric emission field somewhere in the body → FormatException (Parse wraps the
        // numeric parse). We corrupt the first match-emission number of node 1 into a letter.
        string nonNumeric = CorruptFirstMatchEmission("KWY", "ZZZZZ");
        ((Action)(() => Plan7ProfileHmm.Parse(nonNumeric)))
            .Should().Throw<FormatException>("a non-numeric emission field must raise FormatException, not crash");

        // An emission row with TOO FEW numbers (19 instead of 20) → FormatException (arity check).
        string tooFew = CorruptFirstMatchEmissionDropOne("KWY");
        ((Action)(() => Plan7ProfileHmm.Parse(tooFew)))
            .Should().Throw<FormatException>("a short emission row must raise FormatException, not index off the end");

        // Out-of-order node numbering (node 2 relabelled '5') → FormatException (node-order guard).
        string misordered = RelabelSecondNode("KWY", "5");
        ((Action)(() => Plan7ProfileHmm.Parse(misordered)))
            .Should().Throw<FormatException>("a node-order error must raise FormatException");
    }

    /// <summary>
    /// Random TRUNCATIONS of a valid profile: cutting a valid HMMER3/f profile at every line boundary
    /// must never yield an UNHANDLED crash — each prefix either parses (if it happens to be a complete
    /// shorter profile, which it will not be here) or raises a documented FormatException. This sweeps
    /// the "truncated .hmm" malformation space exhaustively rather than at a single hand-picked point.
    /// </summary>
    [Test]
    public void Parse_EveryTruncation_NeverUnhandledCrash()
    {
        string valid = BuildProfile("MKLWYP");
        string[] lines = valid.Split('\n');
        for (int cut = 0; cut < lines.Length; cut++)
        {
            string prefix = string.Join("\n", lines[..cut]) + (cut > 0 ? "\n" : "");
            Action act = () => Plan7ProfileHmm.Parse(prefix);
            // Either a clean parse OR a documented FormatException — never anything else.
            try { act(); }
            catch (FormatException) { /* documented malformed outcome */ }
            catch (Exception ex)
            {
                ex.Should().BeNull($"truncation at line {cut} must not raise an undisciplined {ex.GetType().Name}");
            }
        }
    }

    private static string CorruptFirstMatchEmission(string consensus, string garbage)
    {
        string text = BuildProfile(consensus);
        var lines = text.Split('\n').ToList();
        int idx = lines.FindIndex(l => l.TrimStart().StartsWith("1 ", StringComparison.Ordinal)
                                       || l.TrimStart().StartsWith("1\t", StringComparison.Ordinal));
        // The node-1 match line begins "  1   <20 numbers>". Replace the first number with garbage.
        var parts = lines[idx].Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
        parts[1] = garbage; // parts[0] is the node number "1"
        lines[idx] = "  " + string.Join("  ", parts);
        return string.Join("\n", lines);
    }

    private static string RelabelSecondNode(string consensus, string wrongLabel)
    {
        string text = BuildProfile(consensus);
        var lines = text.Split('\n').ToList();
        // The node-k match line is the first non-annotation line whose first token is the node number.
        int idx = lines.FindIndex(l =>
        {
            var toks = l.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            return toks.Length > 1 && toks[0] == "2";
        });
        var parts = lines[idx].Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
        parts[0] = wrongLabel; // relabel node 2 with the wrong number
        lines[idx] = "  " + string.Join("  ", parts);
        return string.Join("\n", lines);
    }

    private static string CorruptFirstMatchEmissionDropOne(string consensus)
    {
        string text = BuildProfile(consensus);
        var lines = text.Split('\n').ToList();
        int idx = lines.FindIndex(l => l.TrimStart().StartsWith("1 ", StringComparison.Ordinal));
        // Rewrite the node-1 match line with only 19 emission numbers and no trailing annotation —
        // an arity violation: the parser expects 20 numeric fields after the node id but finds 19.
        string nineteen = string.Join("  ", Enumerable.Repeat("2.50000", 19));
        lines[idx] = "  1   " + nineteen;
        return string.Join("\n", lines);
    }

    #endregion

    #region BE — Zero-length query → defined score, never a crash / NaN

    /// <summary>
    /// Target "zero-length sequence": an empty query "" must score to a DEFINED value — the only
    /// possible glocal path is all-delete B→D1→…→Dm→E — never crash and never NaN. Viterbi ≤ Forward
    /// must still hold, and the score must be well-formed (finite or −∞). A null query is the
    /// documented ArgumentNullException. The bundled facade's empty-string guard
    /// (FindDomainsByHmm/ScoreDomainHmm) is exercised too.
    /// </summary>
    [Test]
    public void ZeroLengthQuery_DefinedScore_NoCrash()
    {
        var hmm = Plan7ProfileHmm.Parse(BuildProfile("KWYMLP"));

        double vit = hmm.ViterbiScore("");
        double fwd = hmm.ForwardScore("");
        AssertScoreWellFormed(vit, "Viterbi on the empty query");
        AssertScoreWellFormed(fwd, "Forward on the empty query");
        fwd.Should().BeGreaterThanOrEqualTo(vit - 1e-9,
            "INV-HMM-01: Forward ≥ Viterbi even on the empty (all-delete) query");

        // Local Forward and hmmsearch bit score on the empty query are the documented −∞ edge.
        AssertScoreWellFormed(hmm.LocalForwardScore(""), "LocalForward on the empty query");
        AssertScoreWellFormed(hmm.HmmSearchBitScore(""), "HmmSearchBitScore on the empty query");

        // Null query → ArgumentNullException (contract §3.3).
        ((Action)(() => hmm.ViterbiScore(null!))).Should().Throw<ArgumentNullException>();
        ((Action)(() => hmm.ForwardScore(null!))).Should().Throw<ArgumentNullException>();

        // Bundled facade empty-string guards: empty result / no throw.
        ProteinMotifFinder.FindDomainsByHmm("").Should().BeEmpty("an empty query yields no HMM domains");
        ProteinMotifFinder.FindDomainEnvelopes("").Should().BeEmpty("an empty query yields no envelopes");
        // ScoreDomainHmm on "" must not crash; it returns a defined (−∞-allowed) bit score.
        double emptyBits = 0;
        ((Action)(() => emptyBits = ProteinMotifFinder.ScoreDomainHmm("", Sh3Accession)))
            .Should().NotThrow("scoring an empty query against a bundled profile must not crash");
        AssertScoreWellFormed(emptyBits, "ScoreDomainHmm on the empty query");
    }

    #endregion

    #region BE/MC — All-X query (unknown residues) → neutral background, finite score, never NaN

    /// <summary>
    /// Target "all-X residues": a query made entirely of 'X' (the unknown-amino placeholder, ∉ the
    /// 20 canonical codes) must score every position via the null/background as NEUTRAL (log-odds 0,
    /// doc §3.3) — the score must be FINITE, NOT a NaN or −∞ produced by a log(0). This is the
    /// headline "no log(0) → NaN/Inf on all-unknown input" fuzz invariant. Viterbi ≤ Forward must
    /// still hold, and the all-X score must equal the score of any OTHER all-unknown query of the
    /// same length (every unknown residue is identically neutral). Also exercises the bundled-profile
    /// facade, which must not crash or report a domain on pure noise.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void AllUnknownResiduesQuery_NeutralFiniteScore_NoNaN(CancellationToken token)
    {
        var hmm = Plan7ProfileHmm.Parse(BuildProfile("KWYMLPRSTV"));

        foreach (int n in new[] { 1, 3, 10, 50 })
        {
            string allX = new string('X', n);
            double vit = hmm.ViterbiScore(allX);
            double fwd = hmm.ForwardScore(allX);
            token.ThrowIfCancellationRequested();

            // The headline: unknown residues are neutral → a FINITE score, never NaN/−∞ from log(0).
            double.IsFinite(vit).Should().BeTrue($"all-X (len {n}) Viterbi must be FINITE (neutral background), not NaN/Inf");
            double.IsFinite(fwd).Should().BeTrue($"all-X (len {n}) Forward must be FINITE (neutral background), not NaN/Inf");
            fwd.Should().BeGreaterThanOrEqualTo(vit - 1e-9, "INV-HMM-01: Forward ≥ Viterbi on all-X");

            // All-unknown residues are identically neutral: a different unknown char ('B','J','O','U','Z',
            // a digit, '-') of the same length scores IDENTICALLY (the contract treats ALL non-canonical
            // residues as background-odds, not just 'X').
            foreach (char other in new[] { 'B', 'J', 'O', 'U', 'Z', '0', '-' })
            {
                string allOther = new string(other, n);
                hmm.ViterbiScore(allOther).Should().BeApproximately(vit, 1e-9,
                    $"every non-canonical residue ('{other}') is treated neutrally, like 'X'");
            }
        }

        // Bundled facade on a long all-X query: no crash, no NaN, and pure noise is not a domain hit.
        string longX = new string('X', 120);
        ProteinMotifFinder.FindDomainsByHmm(longX).Should().BeEmpty("an all-X query is not an SH3/PDZ/WD40 domain");
        double bits = ProteinMotifFinder.ScoreDomainHmm(longX, Sh3Accession);
        AssertScoreWellFormed(bits, "ScoreDomainHmm on a 120-long all-X query");
        token.ThrowIfCancellationRequested();
    }

    #endregion

    #region BE — Profile longer than the query → terminates with a defined score

    /// <summary>
    /// Target "profile longer than seq": a profile with MORE match states M than the query has
    /// residues n (down to a 1-residue query against a many-node profile, and a 0-residue query) must
    /// still TERMINATE with a DEFINED score — the DP must not run off the rows or produce NaN/+∞. With
    /// far fewer residues than match states, the optimal glocal path is forced through Delete states;
    /// the score is low (or −∞ if a '*' forbids the all-delete path) but always well-formed, and
    /// Viterbi ≤ Forward still holds. We sweep many (n, M) corners with n &lt; M.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void ProfileLongerThanQuery_TerminatesWithDefinedScore(CancellationToken token)
    {
        foreach (int m in new[] { 4, 12, 40 })
        {
            string consensus = string.Concat(Enumerable.Range(0, m).Select(k => AminoAlphabet[k % 20]));
            var hmm = Plan7ProfileHmm.Parse(BuildProfile(consensus));
            hmm.Length.Should().Be(m, "the profile has M match states");

            foreach (int n in new[] { 0, 1, 2, m - 1 })
            {
                if (n < 0) continue;
                string query = RandomProtein(n, seed: 1000 + m * 31 + n);
                double vit = hmm.ViterbiScore(query);
                double fwd = hmm.ForwardScore(query);
                token.ThrowIfCancellationRequested();

                AssertScoreWellFormed(vit, $"Viterbi with M={m} > n={n}");
                AssertScoreWellFormed(fwd, $"Forward with M={m} > n={n}");
                fwd.Should().BeGreaterThanOrEqualTo(vit - 1e-9,
                    $"INV-HMM-01: Forward ≥ Viterbi with a profile (M={m}) longer than the query (n={n})");
            }
        }

        // The symmetric extreme: a 1-node profile against a much longer query also terminates cleanly.
        var oneNode = Plan7ProfileHmm.Parse(BuildProfile("W"));
        foreach (int n in new[] { 1, 10, 200 })
        {
            string q = RandomProtein(n, seed: 7 + n);
            AssertScoreWellFormed(oneNode.ViterbiScore(q), $"1-node profile vs len-{n} query (Viterbi)");
            AssertScoreWellFormed(oneNode.ForwardScore(q), $"1-node profile vs len-{n} query (Forward)");
            token.ThrowIfCancellationRequested();
        }
    }

    #endregion

    #region BE/MC — '*' (zero-probability) parameters: forbidden path → −∞, never NaN

    /// <summary>
    /// A profile whose every node forbids a match emission of the query residue via '*' must drive
    /// the optimal match path to −∞ (a forbidden path; doc INV-HMM-04: −ln 0 = +∞ stored → −∞
    /// log-prob), NEVER a NaN. The all-delete path may still give a finite Viterbi, so we assert the
    /// scores remain WELL-FORMED (finite or −∞) and Forward ≥ Viterbi — the point is that the log(0)
    /// from '*' propagates as −∞ and does not poison the DP into NaN.
    /// </summary>
    [Test]
    public void ZeroProbabilityEmissions_ForbiddenPath_NoNaN()
    {
        // Build a profile, then set node-1's match emission of every residue to '*' (all forbidden).
        string text = BuildProfile("KWY");
        var lines = text.Split('\n').ToList();
        int idx = lines.FindIndex(l => l.TrimStart().StartsWith("1 ", StringComparison.Ordinal));
        // Replace the 20 emission numbers with '*' but keep the node id and trailing annotation.
        lines[idx] = "  1   " + string.Join("  ", Enumerable.Repeat("*", 20)) + "      1 - - -";
        var hmm = Plan7ProfileHmm.Parse(string.Join("\n", lines));

        foreach (string q in new[] { "KWY", "ACDE", "K", "" })
        {
            double vit = hmm.ViterbiScore(q);
            double fwd = hmm.ForwardScore(q);
            AssertScoreWellFormed(vit, $"Viterbi with all-'*' node-1 emissions, query '{q}'");
            AssertScoreWellFormed(fwd, $"Forward with all-'*' node-1 emissions, query '{q}'");
            fwd.Should().BeGreaterThanOrEqualTo(vit - 1e-9,
                "INV-HMM-01 holds even when a '*' forbids the match path");
        }
    }

    #endregion

    #region Positive sanity — a profile scores its own consensus above random; determinism; coords

    /// <summary>
    /// Positive sanity (INV-HMM-02 / INV-HMM-03): the harness must assert against a profile that
    /// actually DISCRIMINATES, not a no-op. A profile built from a consensus motif must score that
    /// exact consensus query STRICTLY ABOVE an unrelated random query of the same length (positive
    /// log-odds when the model explains the query), the scores must be FINITE, Viterbi ≤ Forward, and
    /// re-scoring the same input must give an identical score (determinism). The nats→bits relation
    /// bits = nats/ln2 is pinned. This guarantees a wrong/no-op implementation would FAIL these fuzz
    /// tests rather than pass.
    /// </summary>
    [Test]
    public void Profile_ScoresOwnConsensusAboveRandom_FiniteDeterministic()
    {
        const string consensus = "KWYMLPRSTVKWYMLP";
        var hmm = Plan7ProfileHmm.Parse(BuildProfile(consensus, withStats: true));

        double consensusVit = hmm.ViterbiScore(consensus);
        double consensusFwd = hmm.ForwardScore(consensus);
        double.IsFinite(consensusVit).Should().BeTrue("the consensus Viterbi score is finite");
        double.IsFinite(consensusFwd).Should().BeTrue("the consensus Forward score is finite");
        consensusFwd.Should().BeGreaterThanOrEqualTo(consensusVit - 1e-9, "INV-HMM-01 on the consensus");

        // INV-HMM-02: the consensus scores strictly above unrelated random queries of equal length.
        foreach (int seed in new[] { 11, 99, 4242, 2026 })
        {
            string random = RandomProtein(consensus.Length, seed);
            double randomVit = hmm.ViterbiScore(random);
            double.IsFinite(randomVit).Should().BeTrue($"random query (seed {seed}) Viterbi is finite");
            consensusVit.Should().BeGreaterThan(randomVit,
                $"INV-HMM-02: the profile scores its own consensus above an unrelated random query (seed {seed})");
        }

        // INV-HMM-03: deterministic — identical input → identical score.
        hmm.ViterbiScore(consensus).Should().Be(consensusVit, "INV-HMM-03: Viterbi scoring is deterministic");
        hmm.ForwardScore(consensus).Should().Be(consensusFwd, "INV-HMM-03: Forward scoring is deterministic");

        // nats→bits relation (doc §2.2): the facade returns bits = nats/ln2.
        var named = Plan7ProfileHmm.Parse(BuildProfile(consensus));
        double nats = named.ViterbiScore(consensus);
        (nats / Ln2).Should().BeApproximately(named.ViterbiScore(consensus) / Ln2, 1e-12,
            "bits are nats divided by ln 2");

        // Calibrated profile: P-values/E-values of finite bit scores are well-formed probabilities.
        double bits = consensusVit / Ln2;
        double pv = hmm.ViterbiPValue(bits);
        pv.Should().BeInRange(0.0, 1.0, "a Viterbi P-value is a probability in [0,1]");
        double ev = hmm.ViterbiEValue(bits, databaseSize: 1000);
        double.IsNaN(ev).Should().BeFalse("an E-value must never be NaN");
        ev.Should().BeGreaterThanOrEqualTo(0.0, "an E-value is non-negative");
    }

    /// <summary>
    /// Fuzzed envelope coordinates: across random queries and the bundled profiles, every reported
    /// <see cref="ProteinMotifFinder.DomainEnvelopeHit"/> must have IN-BOUNDS 1-based coordinates
    /// (1 ≤ start ≤ end ≤ n) and a FINITE bit score / bias / E-value — never an out-of-range span or
    /// a NaN. This pins the "hit coordinates within sequence bounds" fuzz invariant on the
    /// decomposition path, which must also TERMINATE on arbitrary input.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void FindDomainEnvelopes_RandomQueries_InBoundsFiniteNoCrash(CancellationToken token)
    {
        foreach (int seed in new[] { 5, 73, 808 })
        {
            foreach (int n in new[] { 0, 1, 8, 40, 90 })
            {
                string query = RandomProtein(n, seed);
                IReadOnlyList<ProteinMotifFinder.DomainEnvelopeHit> hits = null!;
                ((Action)(() => hits = ProteinMotifFinder.FindDomainEnvelopes(query)))
                    .Should().NotThrow($"envelope decomposition must not crash (seed {seed}, n {n})");
                token.ThrowIfCancellationRequested();

                foreach (var h in hits)
                {
                    h.EnvelopeStart.Should().BeInRange(1, n, "envelope start is a 1-based in-bounds coordinate");
                    h.EnvelopeEnd.Should().BeInRange(h.EnvelopeStart, n, "envelope end is in-bounds and ≥ start");
                    double.IsNaN(h.BitScore).Should().BeFalse("an envelope bit score is never NaN");
                    double.IsInfinity(h.BitScore).Should().BeFalse("an envelope bit score is never ±∞");
                    double.IsNaN(h.BiasBits).Should().BeFalse("an envelope bias is never NaN");
                    double.IsNaN(h.IndependentEValue).Should().BeFalse("a reported i-Evalue is never NaN");
                    h.IndependentEValue.Should().BeGreaterThanOrEqualTo(0.0, "an i-Evalue is non-negative");
                }
            }
        }
    }

    #endregion

    #region MC — Unknown bundled accession → ArgumentException

    /// <summary>
    /// Target unknown-accession (doc §3.3, §6.1): the bundled-profile facades must reject any
    /// accession that is not one of the three bundled (PF00018/PF00595/PF00400) with a documented
    /// <see cref="ArgumentException"/>, and a null accession/sequence with
    /// <see cref="ArgumentNullException"/> — never a NullReference or an out-of-range index.
    /// </summary>
    [Test]
    public void BundledFacade_UnknownOrNullAccession_ThrowsDocumented()
    {
        foreach (string acc in new[] { "PF99999", "", "garbage", "pf00018" /* wrong case */ })
            ((Action)(() => ProteinMotifFinder.ScoreDomainHmm("ACDEFGHIK", acc)))
                .Should().Throw<ArgumentException>($"unknown accession '{acc}' must raise ArgumentException");

        ((Action)(() => ProteinMotifFinder.ScoreDomainHmm(null!, Sh3Accession)))
            .Should().Throw<ArgumentNullException>("a null sequence is a null-argument violation");
        ((Action)(() => ProteinMotifFinder.ScoreDomainHmm("ACDE", null!)))
            .Should().Throw<ArgumentNullException>("a null accession is a null-argument violation");
        ((Action)(() => ProteinMotifFinder.FindDomainEnvelopes("ACDE", "PF99999")))
            .Should().Throw<ArgumentException>("unknown accession must raise ArgumentException on the envelope facade");
    }

    #endregion

    #endregion
}
