using System.Text;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Statistics-area isoelectric-point (pI) unit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (KeyNotFoundException, DivideByZero,
/// NullReferenceException, NaN/Infinity result, …). Every input must result in
/// EITHER a well-defined, theory-correct value, OR a *documented, intentional*
/// validation result. A raw runtime exception, a NaN, a pI outside [0, 14], or a
/// non-terminating bisection on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-PI-001 — Isoelectric Point (pI) Calculation (Statistics)
/// Checklist: docs/checklists/03_FUZZING.md, row 125.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate composition boundaries where the
///     net-charge curve barely (or never) crosses zero inside the search window:
///       – no charged residues  → only the two termini titrate; pI = 6.10 (INV-04),
///         the bisection still converges, no NaN, no infinite loop.
///       – all-acidic           → net charge ≤ 0 across most of [0,14]; pI is driven
///         to the LOW end and must stay clamped ≥ 0; bisection terminates.
///       – all-basic            → net charge ≥ 0 across most of [0,14]; pI is driven
///         to the HIGH end and must stay clamped ≤ 14; bisection terminates.
///       – empty / null         → documented input-guard sentinel 7.0 (pI undefined
///         for a zero-length protein, ASM-03); no bisection over an empty count map,
///         no NaN, no DivideByZero.
/// — docs/checklists/03_FUZZING.md §Description (strategy code BE).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The isoelectric-point contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// API entry (SequenceStatistics.cs, src/.../Seqeron.Genomics.Analysis):
///   • SequenceStatistics.CalculateIsoelectricPoint(string) — protein pI in [0,14]
///
/// Documented behaviour (Isoelectric_Point.md, Test Unit ID SEQ-PI-001):
///   • §2.2 / §4: pI is the pH at which the Henderson–Hasselbalch net charge crosses
///     zero, located by bisection over the standard window [0, 14] using the EMBOSS
///     Epk.dat pKa scale; rounded to 2 decimals.
///   • §2.4 INV-01 / §6.1: 0 ≤ pI ≤ 14 — the bisection is confined to [0, 14], so the
///     result can NEVER escape the window (all-acidic clamps near the low end,
///     all-basic near the high end).
///   • §2.4 INV-02 / §6.1: pI is composition-only — permuting residues leaves pI
///     unchanged (charge is summed over counts, not positions).
///   • §2.4 INV-03: net charge is monotonically non-increasing in pH (root unique).
///   • §2.4 INV-04 / §6.1: a termini-only sequence (no ionizable side chain) →
///     pI = (8.6 + 3.6)/2 = 6.10.
///   • §3.3 / §6.1: null/empty → sentinel 7.0; non-ionizable residues are ignored,
///     case-insensitive, no exception for any string input.
///   • §6.1 worked edges: "DDDD" → 3.23 (low), "KKKK" → 11.27 (high), "A" → 6.10.
///
/// The pKa table pinned below is an independent oracle copy of the documented EMBOSS
/// Epk.dat constants (Isoelectric_Point.md §2.2 / §4.2). The pI returned by the unit
/// is cross-checked against this oracle by confirming the net charge ≈ 0 there.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SequenceIsoelectricPointFuzzTests
{
    #region Helpers — documented bounds, sentinel, and independent net-charge oracle

    /// <summary>Bisection window lower bound (Isoelectric_Point.md §2.2 / §4).</summary>
    private const double MinPh = 0.0;

    /// <summary>Bisection window upper bound (Isoelectric_Point.md §2.2 / §4).</summary>
    private const double MaxPh = 14.0;

    /// <summary>Documented input-guard sentinel for null/empty (§6.1, ASM-03).</summary>
    private const double NeutralSentinel = 7.0;

    /// <summary>Bisection pH resolution: result is rounded to 2 decimals (§4.1).</summary>
    private const double PiResolution = 0.01;

    /// <summary>The seven ionizable side-chain residues — the ONLY residues that move pI.</summary>
    private const string IonizableResidues = "DECYHKR";

    /// <summary>Non-ionizable standard residues — present only to dilute / pad sequences.</summary>
    private const string NeutralResidues = "AGILMFPSTWVNQ";

    /// <summary>N-terminus pKa (basic), EMBOSS Epk.dat (§2.2).</summary>
    private const double NTermPka = 8.6;

    /// <summary>C-terminus pKa (acidic), EMBOSS Epk.dat (§2.2).</summary>
    private const double CTermPka = 3.6;

    /// <summary>Independent oracle copy of the documented EMBOSS pKa table:
    /// (pKa, +1 = basic / −1 = acidic). Isoelectric_Point.md §2.2 / §4.2.</summary>
    private static readonly IReadOnlyDictionary<char, (double pKa, int sign)> Pka =
        new Dictionary<char, (double, int)>
        {
            { 'D', (3.9, -1) },  { 'E', (4.1, -1) },  { 'C', (8.5, -1) },
            { 'Y', (10.1, -1) }, { 'H', (6.5, 1) },   { 'K', (10.8, 1) },
            { 'R', (12.5, 1) }
        };

    /// <summary>Independent Henderson–Hasselbalch net-charge oracle at a given pH that
    /// mirrors the documented formula exactly: basic +1/(1+10^(pH−pKa)),
    /// acidic −1/(1+10^(pKa−pH)), both termini counted once (§2.2).</summary>
    private static double NetChargeOracle(string seq, double pH)
    {
        // Both termini always present (a non-empty protein has both ends).
        double charge = 1.0 / (1.0 + Math.Pow(10, pH - NTermPka));
        charge -= 1.0 / (1.0 + Math.Pow(10, CTermPka - pH));

        foreach (char c in seq.ToUpperInvariant())
        {
            if (!Pka.TryGetValue(c, out var p)) continue;
            if (p.sign > 0)
                charge += 1.0 / (1.0 + Math.Pow(10, pH - p.pKa));
            else
                charge -= 1.0 / (1.0 + Math.Pow(10, p.pKa - pH));
        }

        return charge;
    }

    /// <summary>Independent pI oracle: the documented bisection (§4) over [0,14] for the pH
    /// where the oracle net charge crosses zero, using the same precision and 2-decimal
    /// rounding as the unit. Lets a non-empty input's pI be checked against an INDEPENDENT
    /// root finder (theory-correct equivalence) instead of a slope-dependent absolute-charge
    /// threshold — the residual net charge at a finite-resolution root scales with the local
    /// curve slope, which is large for charge-dense sequences (a real property, not a bug).</summary>
    private static double PiOracle(string seq)
    {
        double low = MinPh, high = MaxPh, pH = NeutralSentinel;
        while (high - low > PiResolution)
        {
            pH = (low + high) / 2.0;
            if (NetChargeOracle(seq, pH) > 0) low = pH; else high = pH;
        }
        return Math.Round(pH, 2);
    }

    /// <summary>The universal well-formedness contract for ANY input: pI is a finite number
    /// (never NaN / ±Infinity) and lies within the documented closed window [0, 14] (INV-01).
    /// The bisection is confined to [0, 14], so no input may make the result escape it.</summary>
    private static void AssertWellFormed(double pi)
    {
        double.IsFinite(pi).Should().BeTrue(
            "pI must be finite — no NaN, no ±Infinity (§8 fuzzing contract)");
        pi.Should().BeInRange(MinPh, MaxPh,
            "INV-01: pI is confined to the bisection window [0, 14]");
    }

    /// <summary>Random string of arbitrary BMP code points (control chars, the null byte,
    /// lone surrogate halves, unicode letters/digits) — fuzz fodder with no pKa entry.</summary>
    private static string RandomBmpChars(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (char)rng.Next(0x0000, 0x10000);
        return new string(chars);
    }

    /// <summary>Random sequence over a given alphabet.</summary>
    private static string RandomOver(Random rng, string alphabet, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) sb.Append(alphabet[rng.Next(alphabet.Length)]);
        return sb.ToString();
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-PI-001 — isoelectric point : fuzz targets (BE)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — hand-checkable documented pI values

    /// <summary>
    /// Positive baseline (not a boundary): the documented worked edge values must reproduce
    /// at the bisection resolution, and the net charge at the RETURNED pI must be ≈ 0 — the
    /// defining property of the isoelectric point. Confirms the suite asserts the real
    /// BUSINESS contract (net charge crosses zero), not merely a non-throwing call.
    ///   • "A"    → 6.10  (termini-only midpoint, INV-04 / §7.1)
    ///   • "DDDD" → 3.23  (acidic, §6.1)
    ///   • "KKKK" → 11.27 (basic, §6.1)
    /// — Isoelectric_Point.md §6.1 / §7.1.
    /// </summary>
    [Test]
    public void Pi_DocumentedEdgeValues_MatchAndNetChargeIsZeroThere()
    {
        var cases = new (string seq, double expected)[]
        {
            ("A", 6.10), ("DDDD", 3.23), ("KKKK", 11.27)
        };

        foreach (var (seq, expected) in cases)
        {
            double pi = SequenceStatistics.CalculateIsoelectricPoint(seq);

            pi.Should().BeApproximately(expected, 0.02,
                $"'{seq}' has documented pI {expected} (§6.1/§7.1)");

            // Defining property: net charge at the returned pI must be ≈ 0.
            Math.Abs(NetChargeOracle(seq, pi)).Should().BeLessThan(0.05,
                $"net charge at the returned pI of '{seq}' must be ≈ 0 (pI definition)");

            AssertWellFormed(pi);
        }
    }

    /// <summary>
    /// Positive baseline: an acidic protein (D/E-rich) has pI &lt; 7 and a basic protein
    /// (K/R-rich) has pI &gt; 7 — the qualitative direction the model must always honour.
    /// — Isoelectric_Point.md §2.1 / §6.1.
    /// </summary>
    [Test]
    public void Pi_AcidicBelowSeven_BasicAboveSeven()
    {
        SequenceStatistics.CalculateIsoelectricPoint("DDEEDDEE")
            .Should().BeLessThan(7.0, "an acidic protein focuses below neutral pH");

        SequenceStatistics.CalculateIsoelectricPoint("KKRRKKRR")
            .Should().BeGreaterThan(7.0, "a basic protein focuses above neutral pH");
    }

    /// <summary>
    /// Positive baseline: a balanced peptide with one acidic (Asp) and one basic (Lys)
    /// residue plus the termini has a near-neutral pI; the net charge at the returned pI
    /// must be ≈ 0. Hand-checkable: the single +/− side chains roughly cancel and the
    /// termini set a mid-range crossing. — Isoelectric_Point.md §2.2 (charge cancellation).
    /// </summary>
    [Test]
    public void Pi_BalancedAcidBasePeptide_IsNearNeutral_NetChargeZero()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint("AGKDAG");

        pi.Should().BeInRange(5.0, 9.0, "one Lys + one Asp + termini ⇒ near-neutral pI");
        Math.Abs(NetChargeOracle("AGKDAG", pi)).Should().BeLessThan(0.05,
            "net charge at the returned pI must be ≈ 0 (pI definition)");
        AssertWellFormed(pi);
    }

    /// <summary>
    /// Positive baseline: INV-02 — pI is composition-only, so permuting residues leaves it
    /// unchanged. A canonical sequence and many random permutations must yield identical pI.
    /// — Isoelectric_Point.md §2.4 INV-02.
    /// </summary>
    [Test]
    public void Pi_IsPermutationInvariant()
    {
        const string canonical = "DEKRHCYAGDEKR";
        double expected = SequenceStatistics.CalculateIsoelectricPoint(canonical);

        var rng = new Random(125_010);
        for (int i = 0; i < 200; i++)
        {
            string shuffled = new string(canonical.OrderBy(_ => rng.Next()).ToArray());
            SequenceStatistics.CalculateIsoelectricPoint(shuffled)
                .Should().Be(expected, "INV-02: pI depends on composition, not order");
        }
    }

    #endregion

    #region BE — Boundary: no charged residues (termini-only crossing → 6.10)

    /// <summary>
    /// BE: a sequence with NO ionizable side chain is the boundary where the net-charge curve
    /// is driven ONLY by the two termini. The documented pI is exactly the pKa midpoint
    /// (8.6 + 3.6)/2 = 6.10 (INV-04) and is independent of how many neutral residues are
    /// present. The bisection must still converge (finite, terminates) — no NaN, no hang.
    /// — Isoelectric_Point.md §2.4 INV-04 / §7.1.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void Pi_NoChargedResidues_IsTerminiMidpoint_610()
    {
        foreach (string seq in new[] { "A", "AG", "GGGGG", "AGILMFPSTWVNQ", new string('G', 5000) })
        {
            double pi = SequenceStatistics.CalculateIsoelectricPoint(seq);
            pi.Should().BeApproximately(6.10, 0.02,
                $"'{(seq.Length > 20 ? string.Concat(seq.AsSpan(0, 20), "…") : seq)}' has only termini ⇒ pI = (8.6+3.6)/2 = 6.10 (INV-04)");
            Math.Abs(NetChargeOracle(seq, pi)).Should().BeLessThan(0.05,
                "net charge at the termini-only pI must be ≈ 0");
            AssertWellFormed(pi);
        }
    }

    /// <summary>
    /// BE/INV-04: pI for a no-charged-residue peptide is invariant to its neutral content —
    /// any random string of purely neutral residues yields 6.10 regardless of length or
    /// composition. Random fuzz over the neutral alphabet, all must converge to the midpoint.
    /// — Isoelectric_Point.md §2.4 INV-04.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void Pi_RandomNeutralOnly_AlwaysTerminiMidpoint()
    {
        var rng = new Random(125_020);
        for (int i = 0; i < 1000; i++)
        {
            int len = rng.Next(1, 100);
            string seq = RandomOver(rng, NeutralResidues, len);
            SequenceStatistics.CalculateIsoelectricPoint(seq)
                .Should().BeApproximately(6.10, 0.02,
                    "any all-neutral sequence ⇒ termini-only midpoint 6.10 (INV-04)");
        }
    }

    #endregion

    #region BE — Boundary: all-acidic (pI driven to the low end, clamped ≥ 0)

    /// <summary>
    /// BE: an all-acidic sequence pushes the net-charge curve negative across nearly all of
    /// [0, 14], driving pI toward the LOW boundary. The result must be a LOW pH (well below
    /// neutral), must stay clamped at or above the documented lower bound 0 (INV-01), and the
    /// bisection must terminate. Tested for D-only, E-only, mixed acidic, and at scale.
    /// — Isoelectric_Point.md §6.1 (acidic-only → low pI) / INV-01.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void Pi_AllAcidic_IsLow_ClampedAndConverges()
    {
        foreach (string seq in new[] { "DDDD", "EEEE", "DEDEDEDE", "DDDDEEEECCCC",
                                       new string('D', 2000), new string('E', 2000) })
        {
            double pi = SequenceStatistics.CalculateIsoelectricPoint(seq);

            AssertWellFormed(pi);                 // INV-01: pI ≥ 0, finite, no hang
            pi.Should().BeLessThan(7.0,
                $"all-acidic '{(seq.Length > 12 ? string.Concat(seq.AsSpan(0, 12), "…") : seq)}' ⇒ pI well below neutral");
            pi.Should().BeGreaterThanOrEqualTo(MinPh,
                "INV-01: pI is clamped to the lower bound 0, never escaping below");
        }
    }

    /// <summary>
    /// BE: adding more acidic residues only ever pushes pI DOWN (or holds it) — never up.
    /// Monotone direction check that the low-end boundary is approached from the correct side
    /// and the curve stays well-formed at every step. — Isoelectric_Point.md §2.4 INV-03.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void Pi_MoreAcidicResidues_DoNotRaisePi()
    {
        double previous = SequenceStatistics.CalculateIsoelectricPoint("A");
        for (int n = 1; n <= 50; n++)
        {
            double pi = SequenceStatistics.CalculateIsoelectricPoint("A" + new string('D', n));
            pi.Should().BeLessThanOrEqualTo(previous + PiResolution,
                "adding acidic residues monotonically lowers (never raises) pI");
            AssertWellFormed(pi);
            previous = pi;
        }
    }

    #endregion

    #region BE — Boundary: all-basic (pI driven to the high end, clamped ≤ 14)

    /// <summary>
    /// BE: an all-basic sequence pushes the net-charge curve positive across nearly all of
    /// [0, 14], driving pI toward the HIGH boundary. The result must be a HIGH pH (well above
    /// neutral), must stay clamped at or below the documented upper bound 14 (INV-01), and the
    /// bisection must terminate. Tested for K-only, R-only, mixed basic, and at scale.
    /// — Isoelectric_Point.md §6.1 (basic-only → high pI) / INV-01.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void Pi_AllBasic_IsHigh_ClampedAndConverges()
    {
        foreach (string seq in new[] { "KKKK", "RRRR", "KRKRKRKR", "KKKKRRRRHHHH",
                                       new string('K', 2000), new string('R', 2000) })
        {
            double pi = SequenceStatistics.CalculateIsoelectricPoint(seq);

            AssertWellFormed(pi);                 // INV-01: pI ≤ 14, finite, no hang
            pi.Should().BeGreaterThan(7.0,
                $"all-basic '{(seq.Length > 12 ? string.Concat(seq.AsSpan(0, 12), "…") : seq)}' ⇒ pI well above neutral");
            pi.Should().BeLessThanOrEqualTo(MaxPh,
                "INV-01: pI is clamped to the upper bound 14, never escaping above");
        }
    }

    /// <summary>
    /// BE: adding more basic residues only ever pushes pI UP (or holds it) — never down.
    /// Monotone direction check that the high-end boundary is approached from the correct side
    /// and the curve stays well-formed at every step. — Isoelectric_Point.md §2.4 INV-03.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void Pi_MoreBasicResidues_DoNotLowerPi()
    {
        double previous = SequenceStatistics.CalculateIsoelectricPoint("A");
        for (int n = 1; n <= 50; n++)
        {
            double pi = SequenceStatistics.CalculateIsoelectricPoint("A" + new string('K', n));
            pi.Should().BeGreaterThanOrEqualTo(previous - PiResolution,
                "adding basic residues monotonically raises (never lowers) pI");
            AssertWellFormed(pi);
            previous = pi;
        }
    }

    #endregion

    #region BE — Boundary: empty / null (documented sentinel 7.0, no bisection over emptiness)

    /// <summary>
    /// BE: the empty string is the zero-length boundary — pI is undefined for a protein with
    /// no residues (ASM-03), so the documented result is the neutral input-guard sentinel 7.0,
    /// reached WITHOUT bisecting over an empty count map (no NaN, no DivideByZero, no hang).
    /// — Isoelectric_Point.md §3.3 / §6.1 (null/empty → 7.0).
    /// </summary>
    [Test]
    public void Pi_EmptyString_IsNeutralSentinel()
    {
        double pi = SequenceStatistics.CalculateIsoelectricPoint(string.Empty);
        pi.Should().Be(NeutralSentinel, "empty ⇒ documented sentinel 7.0 (ASM-03)");
        AssertWellFormed(pi);
    }

    /// <summary>
    /// BE: null is treated identically to empty (IsNullOrEmpty short-circuit) — sentinel 7.0,
    /// no NullReferenceException. — Isoelectric_Point.md §3.3 (null/empty → 7.0).
    /// </summary>
    [Test]
    public void Pi_Null_IsNeutralSentinel_NoThrow()
    {
        var act = () => SequenceStatistics.CalculateIsoelectricPoint(null!);
        act.Should().NotThrow("null is documented as 'no sequence', not an error");
        act().Should().Be(NeutralSentinel, "null ⇒ documented sentinel 7.0");
    }

    /// <summary>
    /// BE: a sequence consisting ONLY of non-ionizable / unrecognized characters has zero
    /// ionizable side chains, so it reduces to the termini-only case (pI 6.10) WITHOUT a
    /// KeyNotFound on any unknown char. This is the "no charged residues among junk" boundary,
    /// distinct from the empty-sequence sentinel. — Isoelectric_Point.md §3.3 (unknown ignored).
    /// </summary>
    [TestCase("X")]
    [TestCase("BJZ")]
    [TestCase("---***")]
    [TestCase("1234567890")]
    [TestCase("AGILMFPSTWV")]
    public void Pi_OnlyNonIonizableOrJunk_ReducesToTerminiMidpoint(string seq)
    {
        var act = () => SequenceStatistics.CalculateIsoelectricPoint(seq);
        act.Should().NotThrow($"'{seq}' has no ionizable side chain but must not throw");
        act().Should().BeApproximately(6.10, 0.02,
            "no ionizable side chain ⇒ termini-only midpoint 6.10 (INV-04)");
    }

    #endregion

    #region BE / fuzz — random garbage: never throws, always well-formed, net charge ≈ 0

    /// <summary>
    /// Fuzz: a large batch of arbitrary BMP strings (control chars, the null byte, lone
    /// surrogate halves, unicode letters/digits, occasionally seeded with real ionizable
    /// residues) must NEVER throw and ALWAYS yield a well-formed pI (finite, within [0,14]).
    /// Core fuzz guarantee: no KeyNotFound, no NaN, no out-of-range pI, no non-terminating
    /// bisection on garbage of any shape or length (incl. 0). For non-empty inputs the
    /// returned pI is cross-checked against the independent bisection oracle root.
    /// — docs/ADVANCED_TESTING_CHECKLIST.md §8 + Isoelectric_Point.md INV-01.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Pi_RandomGarbageStrings_NeverThrow_WellFormed()
    {
        var rng = new Random(125_001);

        for (int iteration = 0; iteration < 3000; iteration++)
        {
            int len = rng.Next(0, 200);
            string input = RandomBmpChars(rng, len);

            double pi = NeutralSentinel;
            var act = () => pi = SequenceStatistics.CalculateIsoelectricPoint(input);

            act.Should().NotThrow($"garbage (len {len}) must never crash pI");
            AssertWellFormed(pi);

            if (input.Length > 0)
                pi.Should().BeApproximately(PiOracle(input), 2 * PiResolution,
                    "the unit's pI must match the independent bisection root for any non-empty input");
        }
    }

    /// <summary>
    /// Fuzz: randomly built sequences over the real amino-acid alphabet (ionizable + neutral
    /// mixed, n ≥ 1, padded with random junk) must stay within [0,14] and have net charge ≈ 0
    /// at the returned pI over many shapes. Confirms the zero-crossing contract holds for
    /// arbitrary realistic compositions, not just the hand-picked edges.
    /// — Isoelectric_Point.md §2.2 / INV-01.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Pi_RandomRealSequences_NetChargeZero_WithinBounds()
    {
        var rng = new Random(125_002);
        string alphabet = IonizableResidues + NeutralResidues;
        const string junk = "XBZJ-*0123 .?\t\n";

        for (int iteration = 0; iteration < 2000; iteration++)
        {
            int residues = rng.Next(1, 150);
            var sb = new StringBuilder();
            for (int i = 0; i < residues; i++)
            {
                sb.Append(alphabet[rng.Next(alphabet.Length)]);
                int pad = rng.Next(0, 3);
                for (int p = 0; p < pad; p++) sb.Append(junk[rng.Next(junk.Length)]);
            }
            string seq = sb.ToString();

            double pi = SequenceStatistics.CalculateIsoelectricPoint(seq);

            AssertWellFormed(pi);
            pi.Should().BeApproximately(PiOracle(seq), 2 * PiResolution,
                "the unit's pI must match the independent bisection root for arbitrary real compositions");
        }
    }

    #endregion
}
