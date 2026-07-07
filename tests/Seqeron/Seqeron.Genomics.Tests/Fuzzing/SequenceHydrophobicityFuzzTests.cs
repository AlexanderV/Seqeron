namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Statistics-area protein hydrophobicity (Kyte-Doolittle GRAVY) unit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (DivideByZeroException,
/// KeyNotFoundException, NullReferenceException, OverflowException, NaN result, …).
/// Every input must result in EITHER a well-defined, theory-correct value, OR a
/// *documented, intentional* validation exception. A raw runtime exception, a NaN,
/// or a hang on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-HYDRO-001 — Hydrophobicity Analysis / GRAVY (Statistics)
/// Checklist: docs/checklists/03_FUZZING.md, row 123.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — empty / single residue, and the degenerate
///          zero-recognized-residue boundary (the GRAVY denominator = 0).
///   • MC = Malformed Content — non-amino-acid characters, all-X, IUPAC ambiguity
///          codes, gaps, stop codons, unicode/control junk that has no scale entry.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes BE, MC).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The GRAVY contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// API entry: SequenceStatistics.CalculateHydrophobicity(string)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs
///    lines 372–396), returning double.
///
/// Documented behaviour (Hydrophobicity_Analysis.md, Test Unit ID SEQ-HYDRO-001):
///   • §2.2 / §3.2: GRAVY = (Σ kd(s_i)) / n — the SUM of Kyte-Doolittle hydropathy
///     values of RECOGNIZED residues divided by the count of RECOGNIZED residues
///     (not raw string length). Biopython gravy() = total_gravy / length.
///   • §3.3 / §5.4: only the 20 standard residues are in the scale [1]. ANY other
///     character (IUPAC ambiguity B/Z/X/J, gaps '-', stop '*', digits, punctuation,
///     unicode junk) is SKIPPED — not added to the sum and NOT counted toward the
///     denominator. No KeyError/KeyNotFound is raised (the source uses TryGetValue).
///   • §6.1 / INV-05: null / empty sequence → GRAVY 0, no exception.
///   • §6.1: when NO residue is recognized (all-X, all-junk) the recognized count is
///     0, so the documented result is GRAVY 0 — there is NO DivideByZero (the source
///     guards `count > 0 ? sum/count : 0`).
///   • §2.4 INV-01: GRAVY of a single recognized residue equals that residue's kd.
///   • §2.4 INV-04: GRAVY is case-insensitive (input is ToUpperInvariant'd).
///
/// Kyte-Doolittle scale pinned in this fixture (Hydrophobicity_Analysis.md §4.2,
/// Biopython kd) — the materially behaviour-defining lookup table:
///   A 1.8  R −4.5  N −3.5  D −3.5  C 2.5  Q −3.5  E −3.5  G −0.4  H −3.2  I 4.5
///   L 3.8  K −3.9  M 1.9   F 2.8   P −1.6 S −0.8  T −0.7  W −0.9  Y −1.3  V 4.2
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SequenceHydrophobicityFuzzTests
{
    #region Helpers

    private const double Tolerance = 1e-9;

    /// <summary>The 20 standard one-letter residues — the ONLY characters with a kd entry.</summary>
    private const string StandardResidues = "ARNDCEQGHILKMFPSTWYV";

    /// <summary>The Kyte-Doolittle hydropathy scale, an independent oracle copy of the
    /// documented table (Hydrophobicity_Analysis.md §4.2). Used to cross-check GRAVY.</summary>
    private static readonly IReadOnlyDictionary<char, double> Kd = new Dictionary<char, double>
    {
        { 'A', 1.8 },  { 'R', -4.5 }, { 'N', -3.5 }, { 'D', -3.5 },
        { 'C', 2.5 },  { 'E', -3.5 }, { 'Q', -3.5 }, { 'G', -0.4 },
        { 'H', -3.2 }, { 'I', 4.5 },  { 'L', 3.8 },  { 'K', -3.9 },
        { 'M', 1.9 },  { 'F', 2.8 },  { 'P', -1.6 }, { 'S', -0.8 },
        { 'T', -0.7 }, { 'W', -0.9 }, { 'Y', -1.3 }, { 'V', 4.2 }
    };

    /// <summary>Independent GRAVY oracle: mean kd over recognized (uppercased) residues,
    /// 0 when none recognized — mirrors the documented contract exactly.</summary>
    private static double GravyOracle(string seq)
    {
        if (string.IsNullOrEmpty(seq)) return 0.0;
        double sum = 0;
        int count = 0;
        foreach (char c in seq.ToUpperInvariant())
            if (Kd.TryGetValue(c, out double v)) { sum += v; count++; }
        return count > 0 ? sum / count : 0.0;
    }

    /// <summary>The universal well-formedness contract that must hold for ANY input: the
    /// GRAVY value is a finite number (never NaN / ±Infinity) and lies within the
    /// physically possible kd envelope [−4.5, +4.5] (the most hydrophilic and most
    /// hydrophobic residue values), since a mean of values in [−4.5, 4.5] stays in
    /// that interval. 0 (the empty / no-recognized-residue sentinel) is inside it.</summary>
    private static void AssertWellFormed(double gravy)
    {
        double.IsFinite(gravy).Should().BeTrue(
            "GRAVY must be a finite number — no NaN, no ±Infinity, no DivideByZero (INV-05/contract)");
        gravy.Should().BeInRange(-4.5 - Tolerance, 4.5 + Tolerance,
            "GRAVY is a mean of kd values, all of which lie in [−4.5, +4.5]");
    }

    /// <summary>Random string of arbitrary BMP code points (0x0000–0xFFFF): control
    /// chars, the null byte, lone surrogate halves, unicode letters/digits — none of
    /// which (except the 20 standard residues) has a kd entry. Fuzz fodder.</summary>
    private static string RandomBmpChars(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (char)rng.Next(0x0000, 0x10000);
        return new string(chars);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-HYDRO-001 — protein hydrophobicity (GRAVY) : fuzz targets (BE, MC)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — hand-computed exact GRAVY

    /// <summary>
    /// Positive baseline (not a boundary): the documented worked example must reproduce
    /// EXACTLY. "FLIV" → (2.8 + 3.8 + 4.5 + 4.2)/4 = 15.3/4 = 3.825, and the hydrophilic
    /// "RKDE" → (−4.5 + −3.9 + −3.5 + −3.5)/4 = −15.4/4 = −3.85. Confirms the suite
    /// asserts the real BUSINESS contract (mean of kd), not just non-throwing.
    /// — Hydrophobicity_Analysis.md §7.1 (FLIV ⇒ 3.825) and §7.1 walk-through (RKDE ⇒ −3.85).
    /// </summary>
    [Test]
    public void Gravy_DocumentedWorkedExamples_MatchHandComputedExactly()
    {
        SequenceStatistics.CalculateHydrophobicity("FLIV")
            .Should().BeApproximately(3.825, Tolerance,
                "(2.8+3.8+4.5+4.2)/4 = 15.3/4 = 3.825 (§7.1)");

        SequenceStatistics.CalculateHydrophobicity("RKDE")
            .Should().BeApproximately(-3.85, Tolerance,
                "(−4.5−3.9−3.5−3.5)/4 = −15.4/4 = −3.85 (§7.1 walk-through)");
    }

    /// <summary>
    /// Positive baseline: a purely hydrophobic peptide (all-Ile, the most hydrophobic
    /// residue) yields a high POSITIVE GRAVY equal to its kd, while a purely hydrophilic
    /// peptide (all-Arg, the most hydrophilic residue) yields a large NEGATIVE GRAVY.
    /// Pins the documented sign convention (positive = hydrophobic, negative = hydrophilic)
    /// and the kd extremes.
    /// — Hydrophobicity_Analysis.md §2.1 (Ile +4.5, Arg −4.5); §3.2 (sign convention).
    /// </summary>
    [Test]
    public void Gravy_HydrophobicVsHydrophilicHomopeptides_HaveDocumentedExtremes()
    {
        SequenceStatistics.CalculateHydrophobicity("IIIIIIIIII")
            .Should().BeApproximately(4.5, Tolerance,
                "an all-Ile peptide is maximally hydrophobic (kd = +4.5)");

        SequenceStatistics.CalculateHydrophobicity("RRRRRRRRRR")
            .Should().BeApproximately(-4.5, Tolerance,
                "an all-Arg peptide is maximally hydrophilic (kd = −4.5)");
    }

    /// <summary>
    /// Positive baseline: GRAVY is case-insensitive (INV-04) — the lowercase form of a
    /// mixed peptide must equal the uppercase GRAVY exactly. Guards against a missing
    /// ToUpperInvariant that would skip lowercase residues and corrupt the denominator.
    /// — Hydrophobicity_Analysis.md §2.4 INV-04 / §3.3 (input uppercased before lookup).
    /// </summary>
    [Test]
    public void Gravy_LowercaseInput_EqualsUppercaseGravy()
    {
        double upper = SequenceStatistics.CalculateHydrophobicity("ACDEFGHIKLMNPQRSTVWY");
        double lower = SequenceStatistics.CalculateHydrophobicity("acdefghiklmnpqrstvwy");

        lower.Should().BeApproximately(upper, Tolerance, "INV-04: GRAVY is case-insensitive");
        AssertWellFormed(upper);
    }

    #endregion

    #region BE — Boundary: empty / null (GRAVY 0, no DivideByZero)

    /// <summary>
    /// BE: the empty string is the lower size boundary — zero residues, so the GRAVY
    /// denominator (recognized-residue count) is 0. The documented result is GRAVY 0,
    /// reached via the `count > 0 ? sum/count : 0` guard — NO DivideByZero, no NaN.
    /// — Hydrophobicity_Analysis.md §6.1 / INV-05 (empty → GRAVY 0).
    /// </summary>
    [Test]
    public void Gravy_EmptyString_IsZero_NoDivideByZero()
    {
        var act = () => SequenceStatistics.CalculateHydrophobicity(string.Empty);

        act.Should().NotThrow("the empty string is a defined boundary, not an error");

        double gravy = act();
        gravy.Should().Be(0.0, "no residues to average ⇒ GRAVY 0 (INV-05)");
        AssertWellFormed(gravy);
    }

    /// <summary>
    /// BE: null is treated identically to empty (IsNullOrEmpty short-circuit,
    /// SequenceStatistics.cs line 380) — GRAVY 0, no NullReferenceException.
    /// — Hydrophobicity_Analysis.md §3.3 (null/empty → GRAVY 0).
    /// </summary>
    [Test]
    public void Gravy_Null_IsZero_NoThrow()
    {
        var act = () => SequenceStatistics.CalculateHydrophobicity(null!);

        act.Should().NotThrow("null is documented as 'no sequence', not an error");
        act().Should().Be(0.0);
    }

    #endregion

    #region BE — Boundary: single residue (GRAVY = that residue's kd)

    /// <summary>
    /// BE/INV-01: a single recognized residue is the n=1 boundary — GRAVY = sum/count
    /// with count 1, so it must equal that residue's kd value EXACTLY. Verified for ALL
    /// 20 standard residues against the oracle scale; guards off-by-one in the
    /// denominator and any wrong scale value.
    /// — Hydrophobicity_Analysis.md §2.4 INV-01 (single residue ⇒ its kd).
    /// </summary>
    [Test]
    public void Gravy_SingleStandardResidue_EqualsItsKdValue()
    {
        foreach (char r in StandardResidues)
        {
            double gravy = SequenceStatistics.CalculateHydrophobicity(r.ToString());

            gravy.Should().BeApproximately(Kd[r], Tolerance,
                $"GRAVY of the single residue '{r}' equals its kd ({Kd[r]}) (INV-01)");
            AssertWellFormed(gravy);
        }
    }

    /// <summary>
    /// BE/INV-01/INV-04: a single LOWERCASE recognized residue still yields its kd
    /// (uppercased before lookup) — n=1 boundary combined with case-insensitivity.
    /// </summary>
    [TestCase('i')]
    [TestCase('r')]
    [TestCase('a')]
    public void Gravy_SingleLowercaseResidue_EqualsItsKdValue(char r)
    {
        char upper = char.ToUpperInvariant(r);
        SequenceStatistics.CalculateHydrophobicity(r.ToString())
            .Should().BeApproximately(Kd[upper], Tolerance,
                $"single lowercase '{r}' ⇒ kd of '{upper}' (INV-01/INV-04)");
    }

    /// <summary>
    /// BE/MC: a single NON-amino-acid character (digit, gap, stop, ambiguity, junk) is
    /// the n=1 boundary where the ONLY residue is unrecognized — recognized count 0, so
    /// GRAVY is the documented 0 sentinel, with NO DivideByZero and NO KeyNotFound.
    /// — Hydrophobicity_Analysis.md §6.1 (no recognized residue ⇒ GRAVY 0).
    /// </summary>
    [TestCase("X")]
    [TestCase("B")]
    [TestCase("Z")]
    [TestCase("J")]
    [TestCase("-")]
    [TestCase("*")]
    [TestCase("1")]
    [TestCase("?")]
    public void Gravy_SingleUnrecognizedChar_IsZero_NoThrow(string seq)
    {
        var act = () => SequenceStatistics.CalculateHydrophobicity(seq);

        act.Should().NotThrow($"'{seq}' has no kd entry but must not throw (TryGetValue)");
        act().Should().Be(0.0, "the only residue is unrecognized ⇒ recognized count 0 ⇒ GRAVY 0");
    }

    #endregion

    #region MC — Malformed Content: all-X / all-unknown (no DivideByZero)

    /// <summary>
    /// MC: every residue is the unknown placeholder 'X' (no kd entry). EVERY residue is
    /// skipped, the recognized count is 0, and the documented result is GRAVY 0 — the
    /// `count > 0` guard prevents a DivideByZero by the zero denominator. Holds at the
    /// n=1 boundary, at scale, and case-insensitively for lowercase 'x'.
    /// — Hydrophobicity_Analysis.md §6.1 / §5.4 (X skipped; no recognized residue ⇒ 0).
    /// </summary>
    [TestCase("X")]
    [TestCase("XX")]
    [TestCase("XXXXXXXXXX")]
    [TestCase("xxxx")]
    public void Gravy_AllX_IsZero_NoDivideByZero(string seq)
    {
        var act = () => SequenceStatistics.CalculateHydrophobicity(seq);

        act.Should().NotThrow("an all-X input must never throw (zero recognized residues)");
        double gravy = act();
        gravy.Should().Be(0.0, "X is undefined in the kd scale ⇒ all skipped ⇒ GRAVY 0");
        AssertWellFormed(gravy);
    }

    /// <summary>
    /// MC: a long all-X homopolymer — the all-unknown boundary at scale. Still GRAVY 0,
    /// still no DivideByZero, and bounded in time (hang guard).
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void Gravy_LongAllX_IsZero_NoCrash()
    {
        double gravy = SequenceStatistics.CalculateHydrophobicity(new string('X', 5000));

        gravy.Should().Be(0.0, "no X is recognized ⇒ GRAVY 0 at scale");
        AssertWellFormed(gravy);
    }

    /// <summary>
    /// MC: a mixed string of ALL the common non-standard codes (ambiguity B/Z/X/J, gap,
    /// stop, digits, punctuation) — none recognized, so the recognized count is 0 and
    /// GRAVY is 0 with no KeyNotFound and no DivideByZero.
    /// — Hydrophobicity_Analysis.md §3.3 (any non-standard char skipped, no exception).
    /// </summary>
    [Test]
    public void Gravy_AllNonStandardCodes_IsZero_NoThrow()
    {
        var act = () => SequenceStatistics.CalculateHydrophobicity("BZXJ-*0123.,;:!@#");

        act.Should().NotThrow("no non-standard character has a kd entry; none may throw");
        act().Should().Be(0.0, "zero recognized residues ⇒ GRAVY 0");
    }

    #endregion

    #region MC — Non-amino-acid chars interleaved: counted by recognized residues only

    /// <summary>
    /// MC: GRAVY divides by the count of RECOGNIZED residues, NOT raw length. Interleaving
    /// junk between recognized residues must NOT change the GRAVY — the denominator
    /// excludes the junk. "I" and "IXX**1-I" must both reduce to the recognized residues:
    /// the latter has two I's ⇒ (4.5+4.5)/2 = 4.5, identical to a clean "II". Guards
    /// against a length-based denominator that would dilute the mean with skipped chars.
    /// — Hydrophobicity_Analysis.md §5.2 (divides by recognized count, not string length).
    /// </summary>
    [Test]
    public void Gravy_JunkInterleaved_DenominatorExcludesUnrecognized()
    {
        double clean = SequenceStatistics.CalculateHydrophobicity("II");
        double dirty = SequenceStatistics.CalculateHydrophobicity("IXX**1-I");

        clean.Should().BeApproximately(4.5, Tolerance, "(4.5+4.5)/2 = 4.5");
        dirty.Should().BeApproximately(clean, Tolerance,
            "junk is skipped — only the two recognized I residues are averaged");

        // A documented mixed example: "AX-RX" reduces to A,R ⇒ (1.8 + −4.5)/2 = −1.35.
        SequenceStatistics.CalculateHydrophobicity("AX-RX")
            .Should().BeApproximately((1.8 - 4.5) / 2.0, Tolerance,
                "recognized residues A,R only ⇒ (1.8−4.5)/2 = −1.35");
    }

    /// <summary>
    /// MC: a recognized peptide PADDED with arbitrary unrecognized characters anywhere
    /// (leading, trailing, interior) yields the SAME GRAVY as the bare peptide, over many
    /// random shapes — the recognized-only denominator is robust to junk position/volume.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Gravy_RecognizedPeptidePaddedWithJunk_EqualsBareGravy()
    {
        var rng = new Random(20260620);

        for (int iteration = 0; iteration < 1000; iteration++)
        {
            // Build a non-empty bare peptide from standard residues.
            int pepLen = rng.Next(1, 30);
            var pep = new char[pepLen];
            for (int i = 0; i < pepLen; i++)
                pep[i] = StandardResidues[rng.Next(StandardResidues.Length)];
            string bare = new string(pep);

            // Splice in random unrecognized junk between residues.
            const string junk = "XBZJ-*0123 .?\t\n";
            var sb = new System.Text.StringBuilder();
            foreach (char c in bare)
            {
                int pad = rng.Next(0, 4);
                for (int p = 0; p < pad; p++) sb.Append(junk[rng.Next(junk.Length)]);
                sb.Append(c);
            }
            int tail = rng.Next(0, 4);
            for (int p = 0; p < tail; p++) sb.Append(junk[rng.Next(junk.Length)]);
            string padded = sb.ToString();

            double bareGravy = SequenceStatistics.CalculateHydrophobicity(bare);
            double paddedGravy = SequenceStatistics.CalculateHydrophobicity(padded);

            paddedGravy.Should().BeApproximately(bareGravy, Tolerance,
                "unrecognized padding is excluded from both sum and denominator");
            bareGravy.Should().BeApproximately(GravyOracle(bare), Tolerance,
                "bare GRAVY equals the independent oracle mean");
            AssertWellFormed(paddedGravy);
        }
    }

    #endregion

    #region MC / BE — Random garbage: never throws, always well-formed

    /// <summary>
    /// MC/BE: a large batch of arbitrary BMP strings (control chars, the null byte, lone
    /// surrogate halves, unicode letters/digits, occasionally seeded with real residues)
    /// must NEVER throw and must ALWAYS yield a well-formed GRAVY (finite, in [−4.5,4.5])
    /// that matches the independent oracle. Core fuzz guarantee: no DivideByZero, no
    /// KeyNotFound, no NaN, no overflow on garbage of any shape or length (incl. 0).
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Gravy_RandomGarbageStrings_NeverThrow_MatchOracle()
    {
        var rng = new Random(987654321);

        for (int iteration = 0; iteration < 3000; iteration++)
        {
            int len = rng.Next(0, 200);
            string input = RandomBmpChars(rng, len);

            double gravy = 0;
            var act = () => gravy = SequenceStatistics.CalculateHydrophobicity(input);

            act.Should().NotThrow($"garbage input (len {len}) must never crash GRAVY");
            gravy.Should().BeApproximately(GravyOracle(input), Tolerance,
                "GRAVY matches the independent recognized-mean oracle");
            AssertWellFormed(gravy);
        }
    }

    /// <summary>
    /// MC/BE: randomly built sequences over the FULL standard residue alphabet (n ≥ 1)
    /// must equal the independent oracle mean exactly over many shapes — cross-checks the
    /// summation and recognized-count denominator against a simple oracle, and confirms
    /// every result is well-formed.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void Gravy_RandomStandardPeptides_MatchOracleMean()
    {
        var rng = new Random(13572468);

        for (int iteration = 0; iteration < 2000; iteration++)
        {
            int len = rng.Next(1, 300);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = StandardResidues[rng.Next(StandardResidues.Length)];
            string seq = new string(chars);

            double gravy = SequenceStatistics.CalculateHydrophobicity(seq);

            gravy.Should().BeApproximately(GravyOracle(seq), Tolerance,
                "GRAVY = mean kd over recognized residues");
            AssertWellFormed(gravy);
        }
    }

    #endregion
}
