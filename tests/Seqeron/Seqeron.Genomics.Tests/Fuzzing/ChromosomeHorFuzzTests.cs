using System.Text;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Chromosome area — alpha-satellite higher-order-repeat (HOR)
/// structure detection (CHROM-HOR-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain input to a unit and
/// asserts that the code NEVER fails in an undisciplined way: no hang/infinite
/// loop, no state corruption, no nonsense output, and no *unhandled* runtime
/// exception (IndexOutOfRangeException from monomer/HOR windowing past the array
/// end, DivideByZeroException on a mean computed over zero monomers or with
/// monomerLength = 0, an HOR size / copy number reported negative, an identity
/// outside its documented [0,100] range, or a NaN where a real value is owed).
/// Every input must resolve to EITHER a well-defined, theory-correct result OR a
/// *documented, intentional* validation exception (ArgumentOutOfRangeException for
/// monomerLength &lt; 1). A raw runtime exception, a hang, or a spurious HOR call on
/// a monomeric / garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: CHROM-HOR-001 — higher-order-repeat detection
/// Checklist: docs/checklists/03_FUZZING.md, row 258 (the final row).
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: a MONOMERIC-ONLY array (the same monomer repeated, no
///          higher-order structure → period 1 / not a HOR), the EMPTY sequence
///          (→ defined no-structure result, no DivByZero), and a PERIOD that is
///          NOT a multiple of the monomer (sequence length not k·monomerLength,
///          i.e. a trailing partial monomer that must be ignored, not indexed
///          past the end), plus monomerLength = 0 / negative / > sequence length.
///   • MC = Malformed Content — non-ACGT characters (digits, IUPAC ambiguity
///          codes, whitespace, lowercase, unicode) fed through the same surface.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The HOR-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Centromeric alpha satellite is organised hierarchically: a block of N diverged
/// ~171 bp monomers forms a higher-order repeat (HOR) UNIT, and that unit is itself
/// tandemly repeated into a near-identical array (McNulty &amp; Sullivan 2018,
/// PMC6121732; docs/algorithms/Chromosome_Analysis/Higher_Order_Repeat_Detection.md
/// §1–§2). The detector is
///   ChromosomeAnalyzer.DetectHigherOrderRepeat(string sequence, int monomerLength = 171)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///    lines 821–905). It returns an EAGER value — the HorResult readonly record
///    struct (HasHigherOrderStructure, MonomersPerUnit, HorUnitLengthBp,
///    HorCopyNumber, MonomerCount, MeanInterHorIdentity, MeanIntraHorIdentity) —
///    so any exception or hang surfaces at the call itself; no `.ToList()` forcing
///    is needed.
///
/// How the detector works (ChromosomeAnalyzer.cs lines 821–905; Higher_Order_
/// Repeat_Detection.md §4.1, §3.3):
///   • monomerLength &lt; 1 → ArgumentOutOfRangeException (documented validation).
///   • null/empty → the no-structure result (false, 1, monomerLength, 0, 0, NaN, NaN).
///   • Split the uppercased array into M = ⌊len / monomerLength⌋ full monomers; the
///     trailing partial monomer (len not a multiple of monomerLength) is dropped.
///     M &lt; 2 → no-structure result with MonomerCount = M.
///   • HOR period = smallest k ∈ [1, ⌊M/2⌋] for which ≥ 90% of the k-periodic
///     monomer pairs (i, i+k) are ≥ 95% identical (gapped global alignment +
///     EMBOSS identity). k = 1 ⇒ homogeneous 1-mer array (NOT a multi-monomer HOR);
///     k ≥ 2 ⇒ genuine HOR ⇒ HasHigherOrderStructure = true.
///   • If no k qualifies → no structure (period 1, copy number = monomer count).
///   • HorUnitLengthBp = period × monomerLength; HorCopyNumber = ⌊M / period⌋;
///     MeanInterHorIdentity = mean identity(i, i+period); MeanIntraHorIdentity =
///     mean pairwise identity among the period distinct monomers of the first unit
///     (NaN when period = 1). Identities are percentages in [0,100].
///
/// Theory-derived invariants pinned below (derived from the doc rule + the test's
/// OWN construction, NOT read back off the code — INV-HOR-01..05, §2.4/§3.3):
///   • GENUINE HOR — concatenate K DISTINCT random 171-bp monomers into a unit
///     U = M1..MK, then tile U exactly N times. The k-periodic pairs at k = K are
///     EXACT copies (100% identity), so the smallest qualifying period is K:
///       HasHigherOrderStructure = true, MonomersPerUnit = K, HorCopyNumber = N,
///       HorUnitLengthBp = K × 171, MonomerCount = K·N, MeanInterHorIdentity = 100,
///       and (INV-HOR-04) MeanInterHorIdentity &gt; MeanIntraHorIdentity (distinct
///       random monomers are far below 95% mutually identical).
///   • MONOMERIC-ONLY — the SAME 171-bp monomer repeated is period 1:
///       HasHigherOrderStructure = false, MonomersPerUnit = 1 (NOT a multi-monomer
///       HOR), HorCopyNumber = MonomerCount, MeanInterHorIdentity = 100.
///   • DIVERGENT-MONOMERIC — K distinct random monomers, NO tiling (each appears
///     once): no period clears the 95% bar ⇒ no structure (period 1, IsHor false).
///   • EMPTY / &lt; 2 monomers / partial monomer — defined no-structure result, no
///     DivByZero and no IndexOutOfRange on the dropped tail.
///   • monomerLength &lt; 1 — ArgumentOutOfRangeException (documented), never DivByZero.
///   • Output bounds for EVERY input: MonomersPerUnit ≥ 1, HorUnitLengthBp ≥ 0,
///     HorCopyNumber ≥ 0, MonomerCount ≥ 0, identities NaN or ∈ [0,100],
///     HorUnitLengthBp = MonomersPerUnit × monomerLength (INV-HOR-02),
///     HorCopyNumber = ⌊MonomerCount / MonomersPerUnit⌋ (INV-HOR-03).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Citations
/// ───────────────────────────────────────────────────────────────────────────
/// • Algorithm doc: docs/algorithms/Chromosome_Analysis/Higher_Order_Repeat_Detection.md.
/// • Source: src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs
///   (DetectHigherOrderRepeat + HorResult + MeanPairwiseIdentity; lines 821–922).
/// • McNulty &amp; Sullivan 2018 (PMC6121732); Rosandić/Paar 2024 (PMC11050224);
///   Warburton &amp; Willard 1990 (JMB 216:3); Alkan 2007.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[Category("Fuzzing")]
public class ChromosomeHorFuzzTests
{
    #region Helpers

    /// <summary>Documented alphoid monomer length (ChromosomeAnalyzer.AlphaSatelliteMonomerLength = 171).</summary>
    private const int MonomerLength = ChromosomeAnalyzer.AlphaSatelliteMonomerLength;

    /// <summary>The documented inter-HOR identity bar (95%); k-periodic exact copies clear it at 100%.</summary>
    private const double InterHorBar = 95.0;

    /// <summary>A deterministic random ACGT monomer of <paramref name="length"/> bp (locally seeded; no shared static Rng).</summary>
    private static string RandomMonomer(int length, Random rng)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>Builds <paramref name="k"/> DISTINCT random 171-bp monomers (re-rolling any accidental duplicate).</summary>
    private static string[] DistinctMonomers(int k, int seed)
    {
        var rng = new Random(seed);
        var monomers = new string[k];
        for (int i = 0; i < k; i++)
        {
            string m;
            do { m = RandomMonomer(MonomerLength, rng); }
            while (monomers.Take(i).Contains(m));
            monomers[i] = m;
        }
        return monomers;
    }

    /// <summary>A genuine HOR array: a K-monomer unit (M1..MK) tiled EXACTLY <paramref name="copies"/> times.</summary>
    private static string GenuineHorArray(int k, int copies, int seed)
    {
        string unit = string.Concat(DistinctMonomers(k, seed));
        return string.Concat(Enumerable.Repeat(unit, copies));
    }

    /// <summary>Asserts the output-contract bounds that must hold for EVERY non-throwing input.</summary>
    private static void AssertInContract(ChromosomeAnalyzer.HorResult r, int monomerLength = MonomerLength)
    {
        r.MonomersPerUnit.Should().BeGreaterThanOrEqualTo(1, "the period search starts at k = 1 (INV-HOR-01)");
        r.MonomerCount.Should().BeGreaterThanOrEqualTo(0, "a monomer count is non-negative");
        r.HorCopyNumber.Should().BeGreaterThanOrEqualTo(0, "a copy number is non-negative");
        r.HorUnitLengthBp.Should().BeGreaterThanOrEqualTo(0, "a unit length in bp is non-negative");

        r.HorUnitLengthBp.Should().Be(r.MonomersPerUnit * monomerLength, "HorUnitLengthBp = period × monomerLength (INV-HOR-02)");
        r.HorCopyNumber.Should().Be(r.MonomerCount / r.MonomersPerUnit, "HorCopyNumber = ⌊MonomerCount / period⌋ (INV-HOR-03)");
        r.HasHigherOrderStructure.Should().Be(r.MonomersPerUnit >= 2, "HasHigherOrderStructure ⇔ period ≥ 2 (INV-HOR-01)");

        foreach (double identity in new[] { r.MeanInterHorIdentity, r.MeanIntraHorIdentity })
            if (!double.IsNaN(identity))
                identity.Should().BeInRange(0.0, 100.0, "a reported identity is a percentage ∈ [0,100], finite");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  CHROM-HOR-001 — higher-order-repeat detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region CHROM-HOR-001 — higher-order-repeat detection

    #region Positive control — a genuine K-monomer HOR is detected with the right size & copy number

    /// <summary>
    /// Theory anchor (INV-HOR-01..04, §2.2 worked example): a unit of K DISTINCT 171-bp
    /// monomers tiled exactly N times is a GENUINE HOR. The smallest period whose
    /// k-periodic monomer pairs are all exact copies is K, so:
    ///   HasHigherOrderStructure = true, MonomersPerUnit = K, HorCopyNumber = N,
    ///   HorUnitLengthBp = K×171, MonomerCount = K·N, MeanInterHorIdentity = 100,
    ///   and inter-HOR identity STRICTLY exceeds intra-HOR identity (the distinct
    ///   random monomers are far below 95% mutually identical). Expected values are
    ///   derived from the construction + the doc rule, NOT read off the code.
    /// </summary>
    [TestCase(2, 6)]
    [TestCase(3, 5)]
    [TestCase(4, 4)]
    [TestCase(5, 8)]
    public void DetectHigherOrderRepeat_GenuineHor_DetectedWithRightSizeAndCopyNumber(int k, int copies)
    {
        string array = GenuineHorArray(k, copies, seed: 7_000 + k * 31 + copies);

        var result = ChromosomeAnalyzer.DetectHigherOrderRepeat(array);

        AssertInContract(result);
        result.HasHigherOrderStructure.Should().BeTrue("a multi-monomer tiled unit is a genuine HOR");
        result.MonomersPerUnit.Should().Be(k, "the smallest exact-copy period equals the {0}-monomer unit size", k);
        result.HorCopyNumber.Should().Be(copies, "{0} monomers / {1}-monomer unit = {2} copies", k * copies, k, copies);
        result.HorUnitLengthBp.Should().Be(k * MonomerLength, "unit length = period × 171");
        result.MonomerCount.Should().Be(k * copies, "the array splits into exactly K·N full monomers");
        result.MeanInterHorIdentity.Should().Be(100.0, "the HOR copies are exact, so same-position monomers are identical");
        result.MeanIntraHorIdentity.Should().BeLessThan(result.MeanInterHorIdentity,
            "INV-HOR-04: distinct (random) intra-unit monomers are far less similar than the exact inter-HOR copies");
    }

    /// <summary>
    /// MC: the genuine-HOR verdict survives lowercasing the whole array — the detector
    /// uppercases internally before alignment (§3.3, INV-HOR-05 determinism).
    /// </summary>
    [Test]
    public void DetectHigherOrderRepeat_LowercaseGenuineHor_SameVerdict()
    {
        string array = GenuineHorArray(k: 3, copies: 5, seed: 4242);

        var upper = ChromosomeAnalyzer.DetectHigherOrderRepeat(array);
        var lower = ChromosomeAnalyzer.DetectHigherOrderRepeat(array.ToLowerInvariant());

        AssertInContract(lower);
        lower.Should().Be(upper, "the sequence is uppercased before alignment, so case cannot change the result");
        lower.MonomersPerUnit.Should().Be(3);
        lower.HasHigherOrderStructure.Should().BeTrue();
    }

    #endregion

    #region BE — Monomeric-only array is NOT a multi-monomer HOR (period 1)

    /// <summary>
    /// BE monomeric-only: the SAME 171-bp monomer repeated has period 1 (adjacent
    /// monomers are already exact copies), so it is a homogeneous 1-mer array, NOT a
    /// multi-monomer HOR: MonomersPerUnit = 1, HasHigherOrderStructure = false,
    /// HorCopyNumber = MonomerCount, MeanInterHorIdentity = 100 (§2.4, §6.1). This is
    /// the checklist's "monomeric-only seq" boundary — it must be classified, not crash.
    /// </summary>
    [TestCase(2)]
    [TestCase(5)]
    [TestCase(12)]
    public void DetectHigherOrderRepeat_MonomericOnly_IsNotAMultiMonomerHor(int copies)
    {
        string monomer = RandomMonomer(MonomerLength, new Random(909));
        string array = string.Concat(Enumerable.Repeat(monomer, copies));

        var result = ChromosomeAnalyzer.DetectHigherOrderRepeat(array);

        AssertInContract(result);
        result.HasHigherOrderStructure.Should().BeFalse("a single repeated monomer has no multi-monomer HOR structure");
        result.MonomersPerUnit.Should().Be(1, "adjacent monomers are exact copies ⇒ period 1");
        result.MonomerCount.Should().Be(copies);
        result.HorCopyNumber.Should().Be(copies, "with period 1 each monomer is its own unit");
        result.MeanInterHorIdentity.Should().Be(100.0, "period-1 neighbours are identical copies of the one monomer");
    }

    /// <summary>
    /// BE divergent-monomeric: K DISTINCT random monomers, each appearing ONCE (no
    /// tiling). No period clears the 95% inter-HOR bar, so the detector reports no
    /// structure (period 1, IsHor false) and MeanInterHorIdentity = NaN (§6.1
    /// "Monomeric, mutually divergent array"). Random monomers are ~25% identical.
    /// </summary>
    [TestCase(4, 31)]
    [TestCase(6, 32)]
    public void DetectHigherOrderRepeat_DivergentMonomers_NoPeriodNoStructure(int k, int seed)
    {
        string array = string.Concat(DistinctMonomers(k, seed));

        var result = ChromosomeAnalyzer.DetectHigherOrderRepeat(array);

        AssertInContract(result);
        result.HasHigherOrderStructure.Should().BeFalse("mutually divergent monomers form no high-identity period");
        result.MonomersPerUnit.Should().Be(1, "no k clears the 95% inter-HOR bar ⇒ reported as period 1");
        result.MonomerCount.Should().Be(k);
        double.IsNaN(result.MeanInterHorIdentity).Should().BeTrue("no accepted period ⇒ inter-HOR identity is undefined (NaN)");
    }

    #endregion

    #region BE — Empty / null / fewer than two monomers → defined no-structure result

    /// <summary>
    /// BE empty: the empty (and null) sequence short-circuits to the documented
    /// no-structure result (false, 1, monomerLength, 0, 0, NaN, NaN) — the split /
    /// mean arithmetic is never reached, so there is NO DivideByZero on zero monomers
    /// (§3.3, §6.1).
    /// </summary>
    [TestCase("")]
    [TestCase(null)]
    public void DetectHigherOrderRepeat_EmptyOrNull_ReturnsNoStructureResult(string? seq)
    {
        var act = () => ChromosomeAnalyzer.DetectHigherOrderRepeat(seq!);
        act.Should().NotThrow("empty/null is short-circuited before any indexing or division");

        var result = ChromosomeAnalyzer.DetectHigherOrderRepeat(seq!);

        AssertInContract(result);
        result.HasHigherOrderStructure.Should().BeFalse();
        result.MonomersPerUnit.Should().Be(1);
        result.HorUnitLengthBp.Should().Be(MonomerLength);
        result.HorCopyNumber.Should().Be(0);
        result.MonomerCount.Should().Be(0);
        double.IsNaN(result.MeanInterHorIdentity).Should().BeTrue();
        double.IsNaN(result.MeanIntraHorIdentity).Should().BeTrue();
    }

    /// <summary>
    /// BE fewer-than-two-monomers: 0 or 1 full monomers cannot show periodicity, so the
    /// detector returns a no-structure result with MonomerCount set accordingly and no
    /// DivByZero on the empty inter-HOR sum (§3.3, §6.1). Length 0 (one base) and
    /// length just under two monomers (2·171 − 1) are the boundaries.
    /// </summary>
    [TestCase(1)]                       // 0 full monomers
    [TestCase(MonomerLength)]           // exactly 1 monomer
    [TestCase(2 * MonomerLength - 1)]   // 1 full monomer + a partial (still < 2)
    public void DetectHigherOrderRepeat_FewerThanTwoMonomers_NoStructureNoDivByZero(int length)
    {
        string seq = RandomMonomer(length, new Random(5_000 + length));

        var act = () => ChromosomeAnalyzer.DetectHigherOrderRepeat(seq);
        act.Should().NotThrow("a sub-two-monomer array must be handled, not divided by zero");

        var result = ChromosomeAnalyzer.DetectHigherOrderRepeat(seq);

        AssertInContract(result);
        result.HasHigherOrderStructure.Should().BeFalse("fewer than two monomers cannot show periodicity");
        result.MonomersPerUnit.Should().Be(1);
        result.MonomerCount.Should().Be(length / MonomerLength, "MonomerCount = ⌊len / 171⌋");
        double.IsNaN(result.MeanInterHorIdentity).Should().BeTrue("no inter-HOR pairs ⇒ NaN, not 0/0");
    }

    #endregion

    #region BE — Period not a multiple of the monomer (trailing partial monomer ignored)

    /// <summary>
    /// BE partial-monomer: a genuine HOR array with arbitrary EXTRA bases appended (so
    /// the total length is NOT a multiple of monomerLength) must drop the trailing
    /// partial monomer and report the SAME HOR as the clean array — never IndexOutOfRange
    /// on the incomplete tail (§3.3, §6.1 "Trailing partial monomer"). The tail length
    /// is swept across the full [1, monomerLength-1] partial-window range.
    /// </summary>
    [TestCase(1)]
    [TestCase(57)]
    [TestCase(MonomerLength - 1)]
    public void DetectHigherOrderRepeat_TrailingPartialMonomer_IgnoredSameAsClean(int tailLen)
    {
        const int k = 3, copies = 5;
        string clean = GenuineHorArray(k, copies, seed: 8_080);
        string tail = RandomMonomer(tailLen, new Random(999));
        string withTail = clean + tail;
        (withTail.Length % MonomerLength).Should().NotBe(0, "the appended tail makes the length a non-multiple of 171");

        var act = () => ChromosomeAnalyzer.DetectHigherOrderRepeat(withTail);
        act.Should().NotThrow("a trailing partial monomer must be dropped, never indexed past the array end");

        var result = ChromosomeAnalyzer.DetectHigherOrderRepeat(withTail);
        var cleanResult = ChromosomeAnalyzer.DetectHigherOrderRepeat(clean);

        AssertInContract(result);
        result.Should().Be(cleanResult, "the dropped partial monomer leaves the full-monomer analysis unchanged");
        result.MonomersPerUnit.Should().Be(k);
        result.MonomerCount.Should().Be(k * copies, "only the K·N FULL monomers are counted");
    }

    #endregion

    #region BE — monomerLength boundaries (0 / negative / > sequence length)

    /// <summary>
    /// BE invalid monomerLength: monomerLength &lt; 1 (0 and negatives) is the documented
    /// validation boundary — it throws ArgumentOutOfRangeException, NOT a DivideByZero
    /// from a ⌊len / 0⌋ split (§3.1, §6.1). The exception is the *defined* contract, so
    /// catching it here is the passing outcome.
    /// </summary>
    [TestCase(0)]
    [TestCase(-1)]
    [TestCase(-171)]
    public void DetectHigherOrderRepeat_NonPositiveMonomerLength_ThrowsArgumentOutOfRange(int monomerLength)
    {
        string seq = GenuineHorArray(3, 4, seed: 12);

        var act = () => ChromosomeAnalyzer.DetectHigherOrderRepeat(seq, monomerLength);

        act.Should().Throw<ArgumentOutOfRangeException>("monomerLength < 1 is invalid, and must be a typed validation exception not a DivByZero")
            .Which.ParamName.Should().Be("monomerLength");
    }

    /// <summary>
    /// BE monomerLength &gt; sequence length: when one monomer is longer than the whole
    /// sequence, ⌊len / monomerLength⌋ = 0 full monomers, so the detector returns the
    /// no-structure result with MonomerCount 0 and no DivByZero — the degenerate split
    /// width is handled, not crashed (§3.3). The bounds use the supplied monomerLength.
    /// </summary>
    [TestCase(200)]
    [TestCase(5_000)]
    public void DetectHigherOrderRepeat_MonomerLengthExceedsSequence_ZeroMonomersNoCrash(int monomerLength)
    {
        string seq = RandomMonomer(MonomerLength, new Random(314)); // 171 bp < monomerLength

        var act = () => ChromosomeAnalyzer.DetectHigherOrderRepeat(seq, monomerLength);
        act.Should().NotThrow("monomerLength > length yields zero full monomers, not a crash");

        var result = ChromosomeAnalyzer.DetectHigherOrderRepeat(seq, monomerLength);

        AssertInContract(result, monomerLength);
        result.MonomerCount.Should().Be(0, "no full monomer of length {0} fits in a {1}-bp sequence", monomerLength, MonomerLength);
        result.HasHigherOrderStructure.Should().BeFalse();
        result.MonomersPerUnit.Should().Be(1);
        result.HorCopyNumber.Should().Be(0);
    }

    /// <summary>
    /// A SMALL custom monomerLength still detects a HOR built at that period: tile a
    /// 4-"monomer" unit of distinct 30-bp blocks, and detection must report period 4 at
    /// the 30-bp split, with unit length 4×30. Confirms monomerLength is honoured as a
    /// parameter, not hard-wired to 171, with no off-by-one in the windowing.
    /// </summary>
    [Test]
    public void DetectHigherOrderRepeat_CustomMonomerLength_DetectsHorAtThatPeriod()
    {
        const int len = 30, k = 4, copies = 5;
        var rng = new Random(271828);
        var blocks = new string[k];
        for (int i = 0; i < k; i++) blocks[i] = RandomMonomer(len, rng);
        string unit = string.Concat(blocks);
        string array = string.Concat(Enumerable.Repeat(unit, copies));

        var result = ChromosomeAnalyzer.DetectHigherOrderRepeat(array, monomerLength: len);

        AssertInContract(result, len);
        result.HasHigherOrderStructure.Should().BeTrue();
        result.MonomersPerUnit.Should().Be(k, "the smallest exact-copy period at the 30-bp split is the 4-block unit");
        result.HorUnitLengthBp.Should().Be(k * len);
        result.HorCopyNumber.Should().Be(copies);
        result.MonomerCount.Should().Be(k * copies);
    }

    #endregion

    #region MC — Malformed Content: non-ACGT characters never crash

    /// <summary>
    /// MC: a HOR array whose unit contains arbitrary non-ACGT characters (IUPAC ambiguity
    /// codes, digits, gaps, whitespace, unicode) must be processed without an unhandled
    /// exception and stay fully in-contract. The exact same non-ACGT bytes appear once per
    /// copy, so the period is still detected (the aligner treats the odd symbols verbatim),
    /// but the only hard guarantee fuzzing pins is: no crash, output in-contract.
    /// </summary>
    [TestCase("NRYWSKMBDHV")]
    [TestCase("-.* 0123\t")]
    [TestCase("ΩβΣ")]
    public void DetectHigherOrderRepeat_UnitWithNonAcgt_StaysInContract(string noise)
    {
        // Build a 3-block unit where each block is a 171-char mix of ACGT and the noise alphabet.
        var rng = new Random(noise.Length * 17 + 3);
        string alphabet = "ACGT" + noise;
        string Block()
        {
            var sb = new StringBuilder(MonomerLength);
            for (int i = 0; i < MonomerLength; i++) sb.Append(alphabet[rng.Next(alphabet.Length)]);
            return sb.ToString();
        }
        string unit = Block() + Block() + Block();
        string array = string.Concat(Enumerable.Repeat(unit, 5));

        var act = () => ChromosomeAnalyzer.DetectHigherOrderRepeat(array);
        act.Should().NotThrow("non-ACGT characters are aligned verbatim, never validated-then-indexed-out-of-range");

        AssertInContract(ChromosomeAnalyzer.DetectHigherOrderRepeat(array));
    }

    /// <summary>
    /// MC: a genuine HOR array sprinkled with random non-ACGT noise inserted at random
    /// positions (which also shifts monomer boundaries) must never crash and must stay
    /// in-contract. Determinism is per-test (locally seeded Random).
    /// </summary>
    [TestCase(11)]
    [TestCase(202)]
    public void DetectHigherOrderRepeat_GenuineHorWithInjectedNoise_StaysInContract(int seed)
    {
        const string noise = "NXRYWSKMBDHV-.* 0123\t";
        var rng = new Random(seed);
        var sb = new StringBuilder(GenuineHorArray(k: 4, copies: 6, seed: seed));
        for (int j = 0; j < 40; j++)
        {
            int pos = rng.Next(sb.Length + 1);
            sb.Insert(pos, noise[rng.Next(noise.Length)]);
        }

        var act = () => ChromosomeAnalyzer.DetectHigherOrderRepeat(sb.ToString());
        act.Should().NotThrow("malformed bases inside a real HOR array must not crash the detector");

        AssertInContract(ChromosomeAnalyzer.DetectHigherOrderRepeat(sb.ToString()));
    }

    #endregion

    #region BE — very long input terminates promptly

    /// <summary>
    /// BE scale: a large genuine HOR array (3-monomer unit × 800 copies = 2400 monomers,
    /// ≈ 410 kb) must terminate promptly — pairwise monomer identities are memoised — and
    /// remain in-contract with the correct period/copy number. Bounded by CancelAfter so a
    /// hang fails loudly rather than wedging the suite.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void DetectHigherOrderRepeat_VeryLongHorArray_TerminatesInContract()
    {
        const int k = 3, copies = 800;
        string array = GenuineHorArray(k, copies, seed: 1_234_567);

        var result = ChromosomeAnalyzer.DetectHigherOrderRepeat(array);

        AssertInContract(result);
        result.HasHigherOrderStructure.Should().BeTrue();
        result.MonomersPerUnit.Should().Be(k);
        result.HorCopyNumber.Should().Be(copies);
        result.MeanInterHorIdentity.Should().BeGreaterThanOrEqualTo(InterHorBar, "an exact HOR clears the 95% inter-HOR bar");
    }

    #endregion

    #endregion
}
