using FsCheck;
using FsCheck.Fluent;
using Seqeron.Genomics.Chromosome;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for chromosome-level analysis.
///
/// Test Units: CHROM-TELO-001 (this file is shared; CHROM-CENT/KARYO/ANEU/SYNT
/// extend it later in sibling regions).
///
/// CHROM-TELO-001 covers <see cref="ChromosomeAnalyzer.AnalyzeTelomeres"/> and
/// <see cref="ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio"/>:
///   R  positions/lengths valid (INV-01/02, exact repeat-unit multiples);
///   P  detected telomere is a repeat tract (purity, motif containment);
///   M  more repeats → strictly longer measured region (length = 6·k);
///   D  identical inputs → identical result.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Chromosome")]
public class ChromosomeProperties
{
    #region CHROM-TELO-001

    // ---- Domain constants (transcribed independently from Telomere_Analysis.md §2.1) ----

    /// <summary>Canonical vertebrate 3'-end telomere repeat (Meyne 1989).</summary>
    private const string Vertebrate3Prime = "TTAGGG";

    /// <summary>Reverse complement of <see cref="Vertebrate3Prime"/>, expected at the 5' end.</summary>
    private const string Vertebrate5Prime = "CCCTAA";

    // ===================== Generators =====================

    /// <summary>
    /// Generates a 3'-telomeric sequence: a non-telomeric prefix of A's whose length is a
    /// multiple of 6 (keeps repeat-sized windows aligned and never extends the tract — the
    /// last 6-base window "AAAAAA" matches TTAGGG at only 1/6 &lt; 0.7), followed by exactly
    /// <c>k</c> pure TTAGGG units placed at the very end. Yields (sequence, k).
    /// </summary>
    private static Arbitrary<(string seq, int k)> Telomeric3PrimeArbitrary() =>
        (from prefixUnits in Gen.Choose(0, 30)
         from k in Gen.Choose(1, 400)
         let prefix = new string('A', 6 * prefixUnits)
         let tract = string.Concat(Enumerable.Repeat(Vertebrate3Prime, k))
         select (prefix + tract, k)).ToArbitrary();

    /// <summary>
    /// Generates a 5'-telomeric sequence: exactly <c>k</c> pure CCCTAA units at the start,
    /// followed by a non-telomeric A-suffix. The first post-tract window "AAAAAA" matches
    /// CCCTAA at only 2/6 &lt; 0.7, so the tract is not extended. Yields (sequence, k).
    /// </summary>
    private static Arbitrary<(string seq, int k)> Telomeric5PrimeArbitrary() =>
        (from k in Gen.Choose(1, 400)
         from suffixLen in Gen.Choose(0, 120)
         let tract = string.Concat(Enumerable.Repeat(Vertebrate5Prime, k))
         select (tract + new string('A', suffixLen), k)).ToArbitrary();

    /// <summary>
    /// Generates non-telomeric DNA: only A and T. Neither motif can reach the 70% threshold —
    /// against TTAGGG the three G positions never match (max 3/6 = 0.5), and against CCCTAA the
    /// three C positions never match (max 3/6 = 0.5). Empty allowed.
    /// </summary>
    private static Arbitrary<string> NonTelomericArbitrary() =>
        Gen.Elements('A', 'T').ArrayOf().Select(a => new string(a)).ToArbitrary();

    /// <summary>Generates arbitrary DNA (A/C/G/T, possibly empty) for general invariants.</summary>
    private static Arbitrary<string> AnyDnaArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Select(a => new string(a)).ToArbitrary();

    /// <summary>Generates non-negative (tsRatio, refRatio &gt; 0, refLength) triples for the T/S oracle.</summary>
    private static Arbitrary<(double ts, double refRatio, double refLen)> TsRatioArbitrary() =>
        (from ts in Gen.Choose(0, 100_000).Select(i => i / 1000.0)
         from rr in Gen.Choose(1, 100_000).Select(i => i / 1000.0)
         from rl in Gen.Choose(0, 100_000).Select(i => (double)i)
         select (ts, rr, rl)).ToArbitrary();

    // ===================== Independent oracles =====================

    /// <summary>
    /// Independent T/S oracle, transcribed verbatim from Telomere_Analysis.md §2.2 /
    /// INV-04 — NOT routed through production.
    /// </summary>
    private static double ExpectedTsLength(double tsRatio, double referenceRatio, double referenceLength) =>
        referenceLength * tsRatio / referenceRatio;

    #endregion

    #region CHROM-TELO-001 — M: more repeats → longer region (central business invariant)

    /// <summary>
    /// INV-M (3' end): a tract of <c>k</c> pure TTAGGG units placed at the very end measures
    /// exactly <c>6·k</c> bases — the measured length is a linear function of repeat count,
    /// hence strictly increasing in k. The internal scan counts only complete 6-bp windows,
    /// and the A-prefix terminates the scan, so no over- or under-count occurs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_ThreePrimeLength_EqualsSixTimesRepeatCount()
    {
        return Prop.ForAll(Telomeric3PrimeArbitrary(), t =>
        {
            var (seq, k) = t;
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, minTelomereLength: 0);
            return (r.TelomereLength3Prime == 6 * k)
                .Label($"3' length {r.TelomereLength3Prime} != 6·{k} for tract of {k} units");
        });
    }

    /// <summary>
    /// INV-M (5' end): a tract of <c>k</c> pure CCCTAA units at the start measures exactly
    /// <c>6·k</c> bases (mirror of the 3' invariant against the reverse-complement motif).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_FivePrimeLength_EqualsSixTimesRepeatCount()
    {
        return Prop.ForAll(Telomeric5PrimeArbitrary(), t =>
        {
            var (seq, k) = t;
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, minTelomereLength: 0);
            return (r.TelomereLength5Prime == 6 * k)
                .Label($"5' length {r.TelomereLength5Prime} != 6·{k} for tract of {k} units");
        });
    }

    /// <summary>
    /// INV-M (strict monotonicity): adding one more repeat unit strictly lengthens the measured
    /// 3' tract. This is the literal "more repeats → longer region" claim, checked as a delta of
    /// exactly one repeat unit rather than relying on the 6·k identity alone.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_OneMoreRepeat_StrictlyLongerThreePrime()
    {
        return Prop.ForAll(Gen.Choose(1, 400).ToArbitrary(), k =>
        {
            string seqK = string.Concat(Enumerable.Repeat(Vertebrate3Prime, k));
            string seqK1 = string.Concat(Enumerable.Repeat(Vertebrate3Prime, k + 1));
            int lenK = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seqK, minTelomereLength: 0).TelomereLength3Prime;
            int lenK1 = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seqK1, minTelomereLength: 0).TelomereLength3Prime;
            return (lenK1 > lenK && lenK1 - lenK == 6)
                .Label($"k={k}: len {lenK} → {lenK1} (expected +6, strictly larger)");
        });
    }

    /// <summary>Exact anchors for the length=6·k law at fixed repeat counts.</summary>
    [TestCase(1, 6, TestName = "3' tract of 1 unit → length 6")]
    [TestCase(100, 600, TestName = "3' tract of 100 units → length 600")]
    [TestCase(200, 1200, TestName = "3' tract of 200 units → length 1200")]
    [Category("Property")]
    public void AnalyzeTelomeres_ThreePrimeLength_KnownRepeatCounts(int k, int expectedLength)
    {
        string seq = new string('A', 1000) + string.Concat(Enumerable.Repeat(Vertebrate3Prime, k));
        var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, minTelomereLength: 0);
        Assert.That(r.TelomereLength3Prime, Is.EqualTo(expectedLength),
            $"a 3' tract of {k} TTAGGG units must measure {expectedLength} bases");
    }

    #endregion

    #region CHROM-TELO-001 — P: detected telomere is a repeat tract

    /// <summary>
    /// INV-P (purity): a pure TTAGGG 3' tract yields perfect repeat purity (1.0) — every
    /// counted base matches the motif. Has3PrimeTelomere is true iff the measured length 6·k
    /// reaches the configured minTelomereLength.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_PureThreePrimeTract_PerfectPurityAndThresholdedFlag()
    {
        return Prop.ForAll(Telomeric3PrimeArbitrary(), t =>
        {
            var (seq, k) = t;
            const int minLen = 500;
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, minTelomereLength: minLen);
            bool purityPerfect = Math.Abs(r.RepeatPurity3Prime - 1.0) < 1e-12;
            bool flag = r.Has3PrimeTelomere == (6 * k >= minLen);
            return (purityPerfect && flag)
                .Label($"k={k}: purity={r.RepeatPurity3Prime}, has3'={r.Has3PrimeTelomere} (len {6 * k} vs {minLen})");
        });
    }

    /// <summary>
    /// INV-P (no false positives): non-telomeric DNA (only A/C, no window can reach 70%
    /// similarity to either motif) yields zero-length tracts and both presence flags false.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_NonTelomericDna_NoTelomereDetected()
    {
        return Prop.ForAll(NonTelomericArbitrary(), seq =>
        {
            if (seq.Length == 0) return true.ToProperty(); // empty is a separate edge case.
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq);
            return (r.TelomereLength3Prime == 0 && r.TelomereLength5Prime == 0
                    && !r.Has3PrimeTelomere && !r.Has5PrimeTelomere)
                .Label($"non-telomeric '{seq}' yielded len5={r.TelomereLength5Prime}, len3={r.TelomereLength3Prime}");
        });
    }

    /// <summary>
    /// INV-P (custom motif, ASM-02): a caller-supplied repeat (Arabidopsis TTTAGGG, 7 bp) is
    /// detected at the 3' end with length = 7·k when passed as <c>telomereRepeat</c>, proving
    /// the motif is parameterized rather than hard-coded.
    /// </summary>
    [Test]
    [Category("Property")]
    public void AnalyzeTelomeres_CustomArabidopsisMotif_DetectedWithCorrectLength()
    {
        const string motif = "TTTAGGG"; // 7 bp
        const int k = 100;
        // Prefix length a multiple of 7 keeps windows aligned; "AAAAAAA" vs TTTAGGG = 1/7 < 0.7.
        string seq = new string('A', 7 * 10) + string.Concat(Enumerable.Repeat(motif, k));
        var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, telomereRepeat: motif, minTelomereLength: 0);
        Assert.Multiple(() =>
        {
            Assert.That(r.TelomereLength3Prime, Is.EqualTo(7 * k), "7-bp motif tract must measure 7·k bases");
            Assert.That(r.RepeatPurity3Prime, Is.EqualTo(1.0).Within(1e-12), "pure custom tract → purity 1.0");
        });
    }

    #endregion

    #region CHROM-TELO-001 — R: positions/lengths valid (INV-01, INV-02)

    /// <summary>
    /// INV-01/INV-02 (R): for ANY DNA input, both measured lengths are ≥ 0, never exceed the
    /// scanned window <c>min(searchLength, |sequence|)</c>, are exact multiples of the repeat-unit
    /// length (only complete units counted), and both purities lie in [0, 1]. When a tract is
    /// non-empty its purity is ≥ 0.7 (every accepted window is ≥ 70% similar).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_LengthsAndPurities_AlwaysValid()
    {
        return Prop.ForAll(AnyDnaArbitrary(), seq =>
        {
            const int searchLength = 10000;
            const int unit = 6; // default TTAGGG
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, searchLength: searchLength);
            int cap = Math.Min(searchLength, seq.Length);

            bool nonNeg = r.TelomereLength5Prime >= 0 && r.TelomereLength3Prime >= 0;          // INV-01
            bool capped = r.TelomereLength5Prime <= cap && r.TelomereLength3Prime <= cap;
            bool multiples = r.TelomereLength5Prime % unit == 0 && r.TelomereLength3Prime % unit == 0;
            bool purityRange = r.RepeatPurity5Prime is >= 0.0 and <= 1.0                        // INV-02
                            && r.RepeatPurity3Prime is >= 0.0 and <= 1.0;
            bool accepted5 = r.TelomereLength5Prime == 0 || r.RepeatPurity5Prime >= 0.7 - 1e-12;
            bool accepted3 = r.TelomereLength3Prime == 0 || r.RepeatPurity3Prime >= 0.7 - 1e-12;

            return (nonNeg && capped && multiples && purityRange && accepted5 && accepted3)
                .Label($"invalid result for '{seq}': len5={r.TelomereLength5Prime}, len3={r.TelomereLength3Prime}, " +
                       $"pur5={r.RepeatPurity5Prime}, pur3={r.RepeatPurity3Prime}, cap={cap}");
        });
    }

    /// <summary>
    /// INV-03 (derived flags): Has5PrimeTelomere and Has3PrimeTelomere are exactly the
    /// comparison of the measured length to <c>minTelomereLength</c>, recomputed independently
    /// across random thresholds on a pure tract whose length is known.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_PresenceFlags_DerivedFromMinLength()
    {
        var arb = (from k in Gen.Choose(1, 400)
                   from minLen in Gen.Choose(0, 3000)
                   select (k, minLen)).ToArbitrary();
        return Prop.ForAll(arb, t =>
        {
            var (k, minLen) = t;
            string seq = string.Concat(Enumerable.Repeat(Vertebrate3Prime, k));
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, minTelomereLength: minLen);
            bool expected3 = r.TelomereLength3Prime >= minLen;
            bool expected5 = r.TelomereLength5Prime >= minLen;
            return (r.Has3PrimeTelomere == expected3 && r.Has5PrimeTelomere == expected5)
                .Label($"k={k}, min={minLen}: has3'={r.Has3PrimeTelomere} (len {r.TelomereLength3Prime}), " +
                       $"has5'={r.Has5PrimeTelomere} (len {r.TelomereLength5Prime})");
        });
    }

    #endregion

    #region CHROM-TELO-001 — IsCriticallyShort rule & short-sequence edge cases

    /// <summary>
    /// IsCriticallyShort rule: independently recomputed as
    /// <c>(Has5 &amp;&amp; len5 &lt; criticalLength) || (Has3 &amp;&amp; len3 &lt; criticalLength)</c>
    /// over random thresholds on a 3' tract. A detected-but-short telomere is critical; a long or
    /// absent one is not.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_IsCriticallyShort_MatchesRule()
    {
        var arb = (from k in Gen.Choose(1, 600)
                   from minLen in Gen.Choose(0, 1500)
                   from critical in Gen.Choose(1, 4000)
                   select (k, minLen, critical)).ToArbitrary();
        return Prop.ForAll(arb, t =>
        {
            var (k, minLen, critical) = t;
            string seq = string.Concat(Enumerable.Repeat(Vertebrate3Prime, k));
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq, minTelomereLength: minLen, criticalLength: critical);
            bool expected = (r.Has5PrimeTelomere && r.TelomereLength5Prime < critical)
                         || (r.Has3PrimeTelomere && r.TelomereLength3Prime < critical);
            return (r.IsCriticallyShort == expected)
                .Label($"k={k}, min={minLen}, crit={critical}: critical={r.IsCriticallyShort}, expected={expected}");
        });
    }

    /// <summary>
    /// Edge case (M4 / §6.1): empty or null sequence returns an all-zero result with
    /// IsCriticallyShort = true (special case).
    /// </summary>
    [Test]
    [Category("Property")]
    public void AnalyzeTelomeres_EmptyOrNull_AllZeroAndCriticallyShort()
    {
        Assert.Multiple(() =>
        {
            foreach (var seq in new[] { "", null })
            {
                var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq!);
                Assert.That(r.TelomereLength5Prime, Is.Zero);
                Assert.That(r.TelomereLength3Prime, Is.Zero);
                Assert.That(r.Has5PrimeTelomere, Is.False);
                Assert.That(r.Has3PrimeTelomere, Is.False);
                Assert.That(r.RepeatPurity5Prime, Is.Zero);
                Assert.That(r.RepeatPurity3Prime, Is.Zero);
                Assert.That(r.IsCriticallyShort, Is.True, "empty input is a special-case critically-short result");
            }
        });
    }

    /// <summary>
    /// Edge case (§6.1): a sequence shorter than the repeat unit cannot contain a complete
    /// window, so both measured lengths are 0 for any sub-6-base input.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_ShorterThanRepeatUnit_ZeroLengths()
    {
        var shortGen = (from len in Gen.Choose(1, 5)
                        from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
                        select new string(chars)).ToArbitrary();
        return Prop.ForAll(shortGen, seq =>
        {
            var r = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq);
            return (r.TelomereLength5Prime == 0 && r.TelomereLength3Prime == 0)
                .Label($"sub-unit '{seq}' gave len5={r.TelomereLength5Prime}, len3={r.TelomereLength3Prime}");
        });
    }

    #endregion

    #region CHROM-TELO-001 — INV-04: T/S ratio helper

    /// <summary>
    /// INV-04: <see cref="ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio"/> equals the
    /// independent proportional oracle <c>referenceLength · tsRatio / referenceRatio</c> within
    /// 1e-9, and is non-negative for non-negative inputs (Cawthon 2002).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EstimateTelomereLengthFromTSRatio_MatchesProportionalOracle()
    {
        return Prop.ForAll(TsRatioArbitrary(), t =>
        {
            var (ts, rr, rl) = t;
            double actual = ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio(ts, rr, rl);
            double expected = ExpectedTsLength(ts, rr, rl);
            return (Math.Abs(actual - expected) < 1e-9 && actual >= 0.0)
                .Label($"T/S({ts},{rr},{rl}) = {actual}, expected {expected}");
        });
    }

    /// <summary>Documented T/S anchors (TestSpec table): exact base-pair estimates.</summary>
    [TestCase(1.0, 1.0, 7000, 7000, TestName = "T/S=1 → reference length 7000")]
    [TestCase(1.5, 1.0, 7000, 10500, TestName = "T/S=1.5 → 10500")]
    [TestCase(0.5, 1.0, 7000, 3500, TestName = "T/S=0.5 → 3500")]
    [TestCase(2.0, 1.0, 7000, 14000, TestName = "T/S=2.0 → 14000")]
    [TestCase(0.0, 1.0, 7000, 0, TestName = "T/S=0 → 0")]
    [Category("Property")]
    public void EstimateTelomereLengthFromTSRatio_CanonicalValues(
        double ts, double refRatio, double refLen, double expected)
    {
        double actual = ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio(ts, refRatio, refLen);
        Assert.That(actual, Is.EqualTo(expected).Within(1e-9));
    }

    #endregion

    #region CHROM-TELO-001 — D: determinism

    /// <summary>
    /// INV-D: identical inputs produce an identical <see cref="ChromosomeAnalyzer.TelomereResult"/>
    /// (all eight fields), checked via record-struct value equality across arbitrary DNA.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeTelomeres_IsDeterministic()
    {
        return Prop.ForAll(AnyDnaArbitrary(), seq =>
        {
            var a = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq);
            var b = ChromosomeAnalyzer.AnalyzeTelomeres("chr", seq);
            return a.Equals(b).Label($"non-deterministic TelomereResult for '{seq}'");
        });
    }

    #endregion

    #region CHROM-CENT-001

    // ---- Domain constants (Centromere_Analysis.md §2.1 Levan 1964 arm-ratio thresholds) ----
    //
    // The checklist row's "P: AT-rich region" is BOGUS: AnalyzeCentromere does NOT use AT
    // content. Its centromere proxy is 15-mer REPEAT content × (1 − GC variability). The
    // genuine "P" tested here is detection on a strongly repetitive window vs. Unknown on a
    // non-repetitive / too-short input. The central rigorous test is the Levan classification
    // oracle recomputed from the RETURNED (Start, End, |sequence|).

    /// <summary>15-mer size used by <c>EstimateRepeatContent</c>; periodic units must be shorter.</summary>
    private const int CentKmerSize = 15;

    /// <summary>The six valid <c>CentromereType</c> values per INV-02 / §2.1 + Unknown.</summary>
    private static readonly string[] ValidCentromereTypes =
        { "Metacentric", "Submetacentric", "Subtelocentric", "Acrocentric", "Telocentric", "Unknown" };

    // ===================== Independent Levan oracle =====================

    /// <summary>
    /// Independent Levan (1964) classifier, transcribed verbatim from §2.1 / source
    /// <c>DetermineCentromereType</c> — NOT routed through production. Mirrors the exact
    /// integer-midpoint and threshold arithmetic: <c>centMid=(start+end)/2</c> (integer div),
    /// <c>pArm=min(centMid, len-centMid)</c>, <c>qArm=max(...)</c>; <c>pArm==0</c> ⇒ Telocentric;
    /// else ratio=qArm/pArm ⇒ ≤1.7 Metacentric, ≤3.0 Submetacentric, &lt;7.0 Subtelocentric,
    /// else Acrocentric.
    /// </summary>
    private static string ExpectedCentromereType(int start, int end, int length)
    {
        int centMid = (start + end) / 2;
        int pArm = Math.Min(centMid, length - centMid);
        int qArm = Math.Max(centMid, length - centMid);

        if (pArm == 0)
            return "Telocentric";

        double armRatio = (double)qArm / pArm;
        return armRatio switch
        {
            <= 1.7 => "Metacentric",
            <= 3.0 => "Submetacentric",
            < 7.0 => "Subtelocentric",
            _ => "Acrocentric"
        };
    }

    // ===================== Generators =====================

    /// <summary>
    /// Generates a strongly repetitive sequence: a short periodic unit (length 2–8, never a
    /// homopolymer, &lt; the 15-mer size) tiled to fill the whole sequence, plus a small
    /// windowSize. A perfectly periodic window of period <c>u</c> contains only <c>u</c> distinct
    /// 15-mers, each repeated, so <c>EstimateRepeatContent ≈ 1.0 &gt; minAlphaSatelliteContent</c> —
    /// the window is accepted and a region is detected. Yields (sequence, windowSize). The total
    /// length is always &gt; windowSize so the scan loop executes. Lengths are kept to a few
    /// thousand bp for speed.
    /// </summary>
    private static Arbitrary<(string seq, int windowSize)> RepetitiveArbitrary() =>
        (from windowSize in Gen.Choose(200, 1000)
         from unitLen in Gen.Choose(2, 8)
         from unitChars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(unitLen)
         from extraWindows in Gen.Choose(2, 6) // total length = (this+1)*windowSize → loop runs
         let unit = MakeNonHomopolymerUnit(unitChars)
         let totalLen = (extraWindows + 1) * windowSize
         let seq = Tile(unit, totalLen)
         select (seq, windowSize)).ToArbitrary();

    /// <summary>
    /// Generates a non-repetitive sequence shorter than its windowSize. Two guarantees combine to
    /// force <c>Unknown</c>: (1) length &lt; windowSize, so the scan loop body never runs; this is
    /// the §6.1 "shorter than the analysis window" edge case. Yields (sequence, windowSize).
    /// </summary>
    private static Arbitrary<(string seq, int windowSize)> ShortNonDetectableArbitrary() =>
        (from windowSize in Gen.Choose(200, 1000)
         from len in Gen.Choose(1, 199) // strictly &lt; minimum windowSize 200 ⇒ always &lt; windowSize
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
         select (new string(chars), windowSize)).ToArbitrary();

    /// <summary>Generates arbitrary DNA (A/C/G/T, possibly empty) for general invariants.</summary>
    private static Arbitrary<string> CentAnyDnaArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Select(a => new string(a)).ToArbitrary();

    // ---- Generator helpers (deterministic constructions) ----

    /// <summary>Ensures the periodic unit is not a homopolymer (a homopolymer yields 1 distinct
    /// 15-mer but identical behaviour; forcing variety keeps the unit a genuine period).</summary>
    private static string MakeNonHomopolymerUnit(char[] chars)
    {
        var s = new string(chars);
        if (s.Distinct().Count() == 1)
        {
            var c = s.ToCharArray();
            c[^1] = c[0] == 'A' ? 'C' : 'A';
            s = new string(c);
        }
        return s;
    }

    /// <summary>Tiles <paramref name="unit"/> to exactly <paramref name="totalLen"/> characters.</summary>
    private static string Tile(string unit, int totalLen)
    {
        var sb = new System.Text.StringBuilder(totalLen);
        while (sb.Length < totalLen)
            sb.Append(unit);
        return sb.ToString(0, totalLen);
    }

    // ===================== THE KEY TEST — Levan classification oracle =====================

    /// <summary>
    /// KEY (§2.1 Levan table, heuristic-independent): for any detected result the returned
    /// <c>CentromereType</c> EQUALS the class recomputed by the independent Levan oracle from the
    /// returned <c>(Start, End, |sequence|)</c>. This validates the classification table directly,
    /// regardless of where the heuristic placed the region. Driven over many randomized strongly
    /// repetitive sequences.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeCentromere_DetectedType_MatchesLevanOracle()
    {
        return Prop.ForAll(RepetitiveArbitrary(), t =>
        {
            var (seq, windowSize) = t;
            var r = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq, windowSize);
            if (!r.Start.HasValue || !r.End.HasValue)
                return true.ToProperty(); // detection asserted separately; oracle applies only when found.
            string expected = ExpectedCentromereType(r.Start.Value, r.End.Value, seq.Length);
            return (r.CentromereType == expected)
                .Label($"len={seq.Length} [{r.Start},{r.End}] → type '{r.CentromereType}', oracle '{expected}'");
        });
    }

    /// <summary>
    /// Exact Levan anchors: a region engineered to sit at a chosen midpoint, with the expected
    /// class recomputed BY HAND from §2.1. For total length L and a region [s,e] with
    /// midpoint m=(s+e)/2, pArm=min(m,L-m), qArm=max(m,L-m), ratio=qArm/pArm.
    ///   • centred (m≈L/2): ratio≈1 ⇒ Metacentric.
    ///   • m=L/4: pArm=L/4, qArm=3L/4, ratio=3 ⇒ Submetacentric.
    ///   • m near one end (ratio ≥7) ⇒ Acrocentric.
    ///   • m=0 (pArm=0) ⇒ Telocentric.
    /// </summary>
    [TestCase(0, 1000, 1000, "Metacentric", TestName = "centred region (mid 500/1000) → Metacentric")]
    [TestCase(400, 600, 1000, "Metacentric", TestName = "tight centred region (mid 500/1000) → Metacentric")]
    [TestCase(0, 500, 1000, "Submetacentric", TestName = "mid 250 of 1000 (ratio 3) → Submetacentric")]
    [TestCase(0, 280, 1000, "Subtelocentric", TestName = "mid 140 of 1000 (ratio ~6.1) → Subtelocentric")]
    [TestCase(0, 200, 1000, "Acrocentric", TestName = "mid 100 of 1000 (ratio 9) → Acrocentric")]
    [TestCase(0, 0, 1000, "Telocentric", TestName = "mid 0 (pArm 0) → Telocentric")]
    [Category("Property")]
    public void AnalyzeCentromere_LevanClassification_HandComputedAnchors(
        int start, int end, int length, string expected)
    {
        // Direct oracle anchor: confirms the hand-computed class equals the transcribed classifier.
        Assert.That(ExpectedCentromereType(start, end, length), Is.EqualTo(expected),
            $"Levan class for mid {(start + end) / 2} of {length}");
    }

    /// <summary>
    /// End-to-end Levan anchor: a real repetitive sequence whose detected region is centred yields a
    /// detected centromere classified consistently with the oracle on the returned boundaries.
    /// </summary>
    [Test]
    [Category("Property")]
    public void AnalyzeCentromere_CentredRepetitive_DetectedAndOracleConsistent()
    {
        const int windowSize = 400;
        string seq = Tile("ACGTGA", 4000); // fully periodic ⇒ repeat content ≈ 1 everywhere
        var r = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq, windowSize);
        Assert.That(r.Start, Is.Not.Null, "a fully repetitive sequence must yield a detected region");
        Assert.That(r.End, Is.Not.Null);
        string expected = ExpectedCentromereType(r.Start!.Value, r.End!.Value, seq.Length);
        Assert.That(r.CentromereType, Is.EqualTo(expected),
            "detected CentromereType must equal the Levan oracle on the returned boundaries");
    }

    #endregion

    #region CHROM-CENT-001 — R: index validity (INV-01) & INV-02/03/04

    /// <summary>
    /// R / INV-01: for a detected centromere, <c>0 ≤ Start ≤ End ≤ |sequence|</c> and
    /// <c>Length == End − Start</c>. Driven over randomized repetitive sequences.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeCentromere_DetectedBoundaries_Valid()
    {
        return Prop.ForAll(RepetitiveArbitrary(), t =>
        {
            var (seq, windowSize) = t;
            var r = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq, windowSize);
            if (!r.Start.HasValue || !r.End.HasValue)
                return (r.Length == 0).Label($"undetected but Length={r.Length} (expected 0)");
            bool ordered = 0 <= r.Start.Value && r.Start.Value <= r.End.Value && r.End.Value <= seq.Length;
            bool lengthOk = r.Length == r.End.Value - r.Start.Value;
            return (ordered && lengthOk)
                .Label($"len={seq.Length} [{r.Start},{r.End}] Length={r.Length}");
        });
    }

    /// <summary>
    /// INV-02 (type domain) + INV-03 (acrocentric flag) + INV-04 &amp; score range: for ANY DNA
    /// input the <c>CentromereType</c> is one of the six valid values, <c>IsAcrocentric</c> iff the
    /// type is "Acrocentric", and <c>AlphaSatelliteContent ∈ [0, 1]</c> (since repeatContent∈[0,1]
    /// and (1−gcVar)≤1).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeCentromere_TypeFlagAndScore_AlwaysValid()
    {
        return Prop.ForAll(CentAnyDnaArbitrary(), seq =>
        {
            var r = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq, windowSize: 300);
            bool typeOk = ValidCentromereTypes.Contains(r.CentromereType);                 // INV-02
            bool flagOk = r.IsAcrocentric == (r.CentromereType == "Acrocentric");          // INV-03
            bool scoreOk = r.AlphaSatelliteContent is >= 0.0 and <= 1.0 + 1e-12;           // INV-04
            return (typeOk && flagOk && scoreOk)
                .Label($"'{seq}': type='{r.CentromereType}', acro={r.IsAcrocentric}, score={r.AlphaSatelliteContent}");
        });
    }

    #endregion

    #region CHROM-CENT-001 — P: detection on repetitive vs non-repetitive

    /// <summary>
    /// P (detection, §6.1): a strongly repetitive sequence (a short unit tiled to fill windows,
    /// repeat content ≈ 1 &gt; threshold) is DETECTED — non-null Start/End, positive Length,
    /// non-Unknown type, score &gt; the default threshold proxy.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeCentromere_StronglyRepetitive_Detected()
    {
        return Prop.ForAll(RepetitiveArbitrary(), t =>
        {
            var (seq, windowSize) = t;
            var r = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq, windowSize);
            bool detected = r.Start.HasValue && r.End.HasValue && r.Length > 0
                            && r.CentromereType != "Unknown" && r.AlphaSatelliteContent > 0.0;
            return detected
                .Label($"repetitive (window {windowSize}, len {seq.Length}) not detected: " +
                       $"start={r.Start}, end={r.End}, type={r.CentromereType}, score={r.AlphaSatelliteContent}");
        });
    }

    /// <summary>
    /// P (non-detection, §6.1): a sequence shorter than the analysis window returns <c>Unknown</c>
    /// with null boundaries, <c>Length==0</c> — the scan loop body never executes.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeCentromere_ShorterThanWindow_Unknown()
    {
        return Prop.ForAll(ShortNonDetectableArbitrary(), t =>
        {
            var (seq, windowSize) = t;
            var r = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq, windowSize);
            return (r.Start is null && r.End is null && r.Length == 0 && r.CentromereType == "Unknown")
                .Label($"len {seq.Length} &lt; window {windowSize} gave start={r.Start}, type={r.CentromereType}");
        });
    }

    /// <summary>
    /// P (non-detection, §6.1): a full-length non-repetitive sequence returns <c>Unknown</c>. We use
    /// a single hand-constructed De Bruijn sequence of order 8 over {A,C,G,T} (every 8-mer occurs
    /// exactly once, hence — a fortiori — every 15-mer is unique), so <c>EstimateRepeatContent==0</c>
    /// in every window and no window reaches the threshold. This is an exact, repeat-free witness
    /// rather than a randomized generator (which cannot guarantee 15-mer uniqueness).
    /// </summary>
    [Test]
    [Category("Property")]
    public void AnalyzeCentromere_NonRepetitiveFullLength_Unknown()
    {
        const int windowSize = 400;
        string seq = DeBruijnOrder8(); // 4^8 = 65536 bases, every 8-mer unique ⇒ every 15-mer unique
        var r = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq, windowSize);
        Assert.Multiple(() =>
        {
            Assert.That(r.Start, Is.Null, "no repetitive window exists ⇒ no detection");
            Assert.That(r.End, Is.Null);
            Assert.That(r.Length, Is.Zero);
            Assert.That(r.CentromereType, Is.EqualTo("Unknown"));
        });
    }

    /// <summary>
    /// Builds a De Bruijn sequence B(4, 8): a cyclic string over {A,C,G,T} in which every length-8
    /// word appears exactly once. As a linear string of length 4^8 every 8-mer (and therefore every
    /// 15-mer) is distinct, giving a guaranteed repeat-content-zero, non-periodic test sequence.
    /// </summary>
    private static string DeBruijnOrder8()
    {
        const string alphabet = "ACGT";
        const int k = 8;
        const int n = 4;
        var a = new int[k * n];
        var seq = new System.Text.StringBuilder();

        void Db(int t, int p)
        {
            if (t > k)
            {
                if (k % p == 0)
                    for (int j = 1; j <= p; j++)
                        seq.Append(alphabet[a[j]]);
            }
            else
            {
                a[t] = a[t - p];
                Db(t + 1, p);
                for (int j = a[t - p] + 1; j < n; j++)
                {
                    a[t] = j;
                    Db(t + 1, t);
                }
            }
        }

        Db(1, 1);
        return seq.ToString();
    }

    #endregion

    #region CHROM-CENT-001 — edge cases & D: determinism

    /// <summary>
    /// Edge (§6.1): null or empty input ⇒ <c>Unknown</c>, null boundaries, <c>Length==0</c>,
    /// <c>AlphaSatelliteContent==0</c>, <c>IsAcrocentric==false</c>.
    /// </summary>
    [TestCase("", TestName = "empty sequence → Unknown")]
    [TestCase(null, TestName = "null sequence → Unknown")]
    [Category("Property")]
    public void AnalyzeCentromere_NullOrEmpty_UnknownAllZero(string? seq)
    {
        var r = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq!);
        Assert.Multiple(() =>
        {
            Assert.That(r.Start, Is.Null);
            Assert.That(r.End, Is.Null);
            Assert.That(r.Length, Is.Zero);
            Assert.That(r.CentromereType, Is.EqualTo("Unknown"));
            Assert.That(r.AlphaSatelliteContent, Is.Zero);
            Assert.That(r.IsAcrocentric, Is.False);
        });
    }

    /// <summary>
    /// P (case-insensitivity, §6.1): lowercase input produces an identical result to uppercase —
    /// the method uppercases before scoring. Checked via full record-struct value equality over
    /// repetitive sequences (where a detected region exercises every field).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeCentromere_LowercaseEqualsUppercase()
    {
        return Prop.ForAll(RepetitiveArbitrary(), t =>
        {
            var (seq, windowSize) = t;
            var upper = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq.ToUpperInvariant(), windowSize);
            var lower = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq.ToLowerInvariant(), windowSize);
            return upper.Equals(lower)
                .Label($"case mismatch (window {windowSize}, len {seq.Length}): {upper} vs {lower}");
        });
    }

    /// <summary>
    /// D (determinism): identical inputs produce an identical <c>CentromereResult</c> (all seven
    /// fields), checked via record-struct value equality across repetitive and arbitrary inputs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeCentromere_IsDeterministic()
    {
        return Prop.ForAll(RepetitiveArbitrary(), t =>
        {
            var (seq, windowSize) = t;
            var a = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq, windowSize);
            var b = ChromosomeAnalyzer.AnalyzeCentromere("chr", seq, windowSize);
            return a.Equals(b).Label($"non-deterministic CentromereResult (window {windowSize}, len {seq.Length})");
        });
    }

    #endregion

    #region CHROM-KARYO-001

    // Covers ChromosomeAnalyzer.AnalyzeKaryotype and ChromosomeAnalyzer.DetectPloidy
    // (Karyotype_Analysis.md §2.4 INV-01..INV-04, §4.2 ploidy formula + aneuploidy table, §6 edges).
    //
    //   R  chromosome count > 0 for non-empty input; depth-ploidy ∈ [1,8], confidence ∈ [0,1]
    //   D  identical input → identical result (both methods)
    //   P  every karyotype field recomputed by an independent oracle, including the
    //      per-autosome-group aneuploidy classification (P: each chrom group has a classification)
    //
    // All oracles below are transcribed INDEPENDENTLY from the doc/source semantics and are
    // NEVER routed through production.

    // ===================== Independent oracles =====================

    /// <summary>
    /// Independent base-name oracle, transcribed verbatim from the source helper
    /// <c>GetChromosomeBaseName</c> (Karyotype_Analysis.md §3.3 / §5.2): strip a trailing
    /// <c>_&lt;digits&gt;</c> suffix ONLY (e.g. <c>chr1_2 → chr1</c>); a letter suffix such as
    /// <c>chr1a</c> is left intact, and a leading-underscore name is not altered.
    /// </summary>
    private static string ExpectedBaseName(string name)
    {
        int underscoreIdx = name.LastIndexOf('_');
        if (underscoreIdx > 0 && underscoreIdx < name.Length - 1
            && int.TryParse(name[(underscoreIdx + 1)..], out _))
        {
            return name[..underscoreIdx];
        }
        return name;
    }

    /// <summary>
    /// Independent cytogenetic-term oracle, transcribed verbatim from the source helper
    /// <c>GetAneuploidyTerm</c> (§4.2 absolute-term table): 0→Nullisomy, 1→Monosomy, 2→Disomy,
    /// 3→Trisomy, 4→Tetrasomy, 5→Pentasomy, else <c>Polysomy ({n} copies)</c>.
    /// </summary>
    private static string ExpectedAneuploidyTerm(int copyCount) => copyCount switch
    {
        0 => "Nullisomy",
        1 => "Monosomy",
        2 => "Disomy",
        3 => "Trisomy",
        4 => "Tetrasomy",
        5 => "Pentasomy",
        _ => $"Polysomy ({copyCount} copies)"
    };

    /// <summary>
    /// Independent abnormality-set oracle (§2.4 INV-03 + §5.4 deviation #1): group ONLY the
    /// autosomes (sex chromosomes are excluded from grouping) by <see cref="ExpectedBaseName"/>;
    /// a group is abnormal iff its copy count differs from <paramref name="expectedPloidyLevel"/>;
    /// each abnormality string is <c>"{term} {baseName}"</c>. Returned as a set (production order
    /// is GroupBy-dependent and not part of the contract).
    /// </summary>
    private static HashSet<string> ExpectedAbnormalitySet(
        IReadOnlyList<(string Name, long Length, bool IsSexChromosome)> chroms,
        int expectedPloidyLevel)
    {
        return chroms
            .Where(c => !c.IsSexChromosome)
            .GroupBy(c => ExpectedBaseName(c.Name))
            .Where(g => g.Count() != expectedPloidyLevel)
            .Select(g => $"{ExpectedAneuploidyTerm(g.Count())} {g.Key}")
            .ToHashSet();
    }

    /// <summary>
    /// Independent ploidy oracle, transcribed verbatim from the source/§4.2 arithmetic.
    /// Median is the TRUE median of the sorted depths (average of the two middle values for an
    /// even count); <c>ratio = median / expectedDiploidDepth</c>;
    /// <c>ploidy = clamp(Math.Round(ratio*2), 1, 8)</c> using .NET banker's rounding
    /// (<see cref="MidpointRounding.ToEven"/>, the default of <see cref="Math.Round(double)"/>);
    /// <c>confidence = max(0, 1 − |ratio*2 − ploidy| · 2)</c> using the CLAMPED ploidy.
    /// Empty input ⇒ (2, 0).
    /// </summary>
    private static (int PloidyLevel, double Confidence) ExpectedPloidy(
        IReadOnlyList<double> depths, double expectedDiploidDepth)
    {
        if (depths.Count == 0)
            return (2, 0);

        var sorted = depths.OrderBy(d => d).ToList();
        double median = sorted.Count % 2 == 1
            ? sorted[sorted.Count / 2]
            : (sorted[sorted.Count / 2 - 1] + sorted[sorted.Count / 2]) / 2.0;
        double ratio = median / expectedDiploidDepth;

        int ploidy = (int)Math.Round(ratio * 2); // banker's rounding (round-half-to-even)
        ploidy = Math.Max(1, Math.Min(8, ploidy));

        double confidence = Math.Max(0.0, 1.0 - Math.Abs(ratio * 2 - ploidy) * 2);
        return (ploidy, confidence);
    }

    // ===================== Generators =====================

    /// <summary>
    /// Generates a chromosome tuple list with controllable base names, optional <c>_N</c> copy
    /// suffixes, sex flags, and positive lengths. Several distinct base names are emitted per
    /// group with varying copy counts (1..6) so aneuploidy grouping and classification are
    /// exercised across all term branches; a configurable subset is flagged as sex chromosomes
    /// (which must be excluded from abnormality grouping). Lengths are kept ≤ 250M and counts
    /// small so Σ never approaches <c>long</c> overflow. Empty list is NOT produced here (the
    /// empty edge case is a dedicated test).
    /// </summary>
    private static Arbitrary<List<(string Name, long Length, bool IsSexChromosome)>> KaryotypeArbitrary() =>
        (from groupCount in Gen.Choose(1, 5)
         from groups in GenGroups(groupCount)
         from sexCount in Gen.Choose(0, 3)
         let sex = Enumerable.Range(0, sexCount)
             .Select(i => ($"chrSex{i}", (long)(1_000_000 + i), true))
             .ToList()
         let all = groups.Concat(sex).ToList()
         select Shuffle(all)).ToArbitrary();

    /// <summary>Generates <paramref name="groupCount"/> autosome groups, each a base name with a
    /// random copy count (1..6) expanded into <c>base_1..base_n</c> tuples with positive lengths.</summary>
    private static Gen<List<(string Name, long Length, bool IsSexChromosome)>> GenGroups(int groupCount) =>
        from copies in Gen.Choose(1, 6).ListOf(groupCount)
        let tuples = copies
            .SelectMany((count, gi) =>
                Enumerable.Range(1, count).Select(ci =>
                    ($"chrA{gi}_{ci}", (long)(10_000_000 + gi * 1_000 + ci), false)))
            .ToList()
        select tuples;

    /// <summary>Deterministic Fisher-Yates-free shuffle by a fixed index permutation derived from
    /// the list — keeps the ordering nontrivial yet reproducible so generation stays deterministic.</summary>
    private static List<T> Shuffle<T>(List<T> items)
    {
        // Rotate by a length-derived offset: cheap, deterministic, and decouples sex/autosome order.
        if (items.Count <= 1) return items;
        int offset = items.Count / 2 + 1;
        return items.Skip(offset).Concat(items.Take(offset)).ToList();
    }

    /// <summary>
    /// Generates a non-empty normalized-depth list plus a positive <c>expectedDiploidDepth</c>.
    /// Depths and the divisor are finite positives in a sane range to avoid div-by-zero and
    /// floating noise.
    /// </summary>
    private static Arbitrary<(List<double> depths, double expectedDiploid)> DepthsArbitrary() =>
        (from n in Gen.Choose(1, 40)
         from raw in Gen.Choose(1, 8000).Select(i => i / 1000.0).ListOf(n)
         from divInt in Gen.Choose(500, 4000)
         select (raw, divInt / 1000.0)).ToArbitrary();

    // ===================== AnalyzeKaryotype — INV-01/02/03 field oracle =====================

    /// <summary>
    /// INV-01/02/03 (R + P): every field of the returned <see cref="ChromosomeAnalyzer.Karyotype"/>
    /// equals the independent oracle recomputed from the input list —
    /// <c>TotalChromosomes == input.Count == AutosomeCount + SexChromosomes.Count</c> (INV-01);
    /// <c>TotalGenomeSize == Σ lengths</c>, <c>MeanChromosomeLength == TotalGenomeSize/count</c>,
    /// autosome count and the ordered sex-chromosome names (INV-02);
    /// the abnormality SET and <c>HasAneuploidy ⟺ abnormalities non-empty</c> with sex
    /// chromosomes excluded from grouping (INV-03); and <c>PloidyLevel == expectedPloidyLevel</c>.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeKaryotype_AllFields_MatchIndependentOracle()
    {
        var arb = (from chroms in KaryotypeArbitrary().Generator
                   from ploidy in Gen.Choose(1, 4)
                   select (chroms, ploidy)).ToArbitrary();

        return Prop.ForAll(arb, t =>
        {
            var (chroms, expectedPloidy) = t;
            var k = ChromosomeAnalyzer.AnalyzeKaryotype(chroms, expectedPloidy);

            int count = chroms.Count;
            long totalSize = chroms.Sum(c => c.Length);
            double meanLen = totalSize / (double)count;
            var expectedSex = chroms.Where(c => c.IsSexChromosome).Select(c => c.Name).ToList();
            int expectedAutosomes = chroms.Count(c => !c.IsSexChromosome);
            var expectedAbn = ExpectedAbnormalitySet(chroms, expectedPloidy);

            bool inv01 = k.TotalChromosomes == count
                         && k.TotalChromosomes == k.AutosomeCount + k.SexChromosomes.Count;
            bool inv02 = k.TotalGenomeSize == totalSize
                         && Math.Abs(k.MeanChromosomeLength - meanLen) < 1e-9
                         && k.AutosomeCount == expectedAutosomes
                         && k.SexChromosomes.SequenceEqual(expectedSex);
            bool inv03 = k.Abnormalities.ToHashSet().SetEquals(expectedAbn)
                         && k.HasAneuploidy == (expectedAbn.Count > 0);
            bool ploidyEcho = k.PloidyLevel == expectedPloidy;

            return (inv01 && inv02 && inv03 && ploidyEcho)
                .Label($"count={count}, ploidy={expectedPloidy}: " +
                       $"INV01={inv01}, INV02={inv02}, INV03={inv03}, ploidyEcho={ploidyEcho}; " +
                       $"abn=[{string.Join("|", k.Abnormalities)}] vs oracle=[{string.Join("|", expectedAbn)}]");
        });
    }

    /// <summary>
    /// R (count &gt; 0): for any non-empty input the chromosome count is strictly positive and
    /// equals the input size — the checklist's literal "chromosome count &gt; 0" invariant.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeKaryotype_NonEmpty_PositiveCount()
    {
        return Prop.ForAll(KaryotypeArbitrary(), chroms =>
        {
            var k = ChromosomeAnalyzer.AnalyzeKaryotype(chroms);
            return (k.TotalChromosomes > 0 && k.TotalChromosomes == chroms.Count)
                .Label($"count={k.TotalChromosomes} for input of {chroms.Count}");
        });
    }

    /// <summary>
    /// INV-03 (sex exclusion, §5.4 deviation #1): a karyotype whose ONLY abnormal-count group is a
    /// sex chromosome produces NO abnormality. Here X is present three times (abnormal for diploid)
    /// but flagged as a sex chromosome, while every autosome group is disomic. Expect
    /// <c>HasAneuploidy == false</c> and an empty abnormality list.
    /// </summary>
    [Test]
    [Category("Property")]
    public void AnalyzeKaryotype_SexChromosomeAbnormalCount_NotReported()
    {
        var chroms = new List<(string, long, bool)>
        {
            ("chr1_1", 100, false), ("chr1_2", 100, false), // disomic autosome
            ("chrX_1", 50, true), ("chrX_2", 50, true), ("chrX_3", 50, true) // triple X, but sex → ignored
        };
        var k = ChromosomeAnalyzer.AnalyzeKaryotype(chroms, expectedPloidyLevel: 2);
        Assert.Multiple(() =>
        {
            Assert.That(k.HasAneuploidy, Is.False, "sex-chromosome copy abnormalities are excluded from grouping");
            Assert.That(k.Abnormalities, Is.Empty);
            Assert.That(k.SexChromosomes, Is.EquivalentTo(new[] { "chrX_1", "chrX_2", "chrX_3" }));
            Assert.That(k.AutosomeCount, Is.EqualTo(2));
            Assert.That(k.TotalChromosomes, Is.EqualTo(5));
        });
    }

    /// <summary>
    /// INV-03 labeling anchor (§4.2 table): a hand-built karyotype exercises every term branch —
    /// Monosomy (1 copy), Trisomy (3), Tetrasomy (4), Pentasomy (5), Polysomy (7), and a normal
    /// disomic group that produces no label. Asserts the exact abnormality set.
    /// </summary>
    [Test]
    [Category("Property")]
    public void AnalyzeKaryotype_AneuploidyTerms_ExactLabels()
    {
        List<(string, long, bool)> Copies(string b, int n) =>
            Enumerable.Range(1, n).Select(i => ($"{b}_{i}", (long)100, false)).ToList();

        var chroms = new List<(string, long, bool)>();
        chroms.AddRange(Copies("chrMono", 1));
        chroms.AddRange(Copies("chrDi", 2));   // normal → no label
        chroms.AddRange(Copies("chrTri", 3));
        chroms.AddRange(Copies("chrTetra", 4));
        chroms.AddRange(Copies("chrPenta", 5));
        chroms.AddRange(Copies("chrPoly", 7));

        var k = ChromosomeAnalyzer.AnalyzeKaryotype(chroms, expectedPloidyLevel: 2);
        var expected = new[]
        {
            "Monosomy chrMono",
            "Trisomy chrTri",
            "Tetrasomy chrTetra",
            "Pentasomy chrPenta",
            "Polysomy (7 copies) chrPoly"
        };
        Assert.Multiple(() =>
        {
            Assert.That(k.Abnormalities, Is.EquivalentTo(expected));
            Assert.That(k.HasAneuploidy, Is.True);
        });
    }

    /// <summary>
    /// Edge (§6.1): empty input ⇒ a fully zeroed <see cref="ChromosomeAnalyzer.Karyotype"/> —
    /// counts and sizes 0, empty sex/abnormality lists, PloidyLevel 0, no aneuploidy.
    /// </summary>
    [Test]
    [Category("Property")]
    public void AnalyzeKaryotype_Empty_ZeroedResult()
    {
        var k = ChromosomeAnalyzer.AnalyzeKaryotype(
            Enumerable.Empty<(string, long, bool)>(), expectedPloidyLevel: 2);
        Assert.Multiple(() =>
        {
            Assert.That(k.TotalChromosomes, Is.Zero);
            Assert.That(k.AutosomeCount, Is.Zero);
            Assert.That(k.SexChromosomes, Is.Empty);
            Assert.That(k.TotalGenomeSize, Is.Zero);
            Assert.That(k.MeanChromosomeLength, Is.Zero);
            Assert.That(k.PloidyLevel, Is.Zero);
            Assert.That(k.HasAneuploidy, Is.False);
            Assert.That(k.Abnormalities, Is.Empty);
        });
    }

    /// <summary>
    /// D (determinism): identical input ⇒ an identical <see cref="ChromosomeAnalyzer.Karyotype"/>.
    /// Record-struct equality covers value-typed fields; the list fields are compared element-wise.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AnalyzeKaryotype_IsDeterministic()
    {
        return Prop.ForAll(KaryotypeArbitrary(), chroms =>
        {
            var a = ChromosomeAnalyzer.AnalyzeKaryotype(chroms, 2);
            var b = ChromosomeAnalyzer.AnalyzeKaryotype(chroms, 2);
            bool same = a.TotalChromosomes == b.TotalChromosomes
                        && a.AutosomeCount == b.AutosomeCount
                        && a.SexChromosomes.SequenceEqual(b.SexChromosomes)
                        && a.TotalGenomeSize == b.TotalGenomeSize
                        && a.MeanChromosomeLength.Equals(b.MeanChromosomeLength)
                        && a.PloidyLevel == b.PloidyLevel
                        && a.HasAneuploidy == b.HasAneuploidy
                        && a.Abnormalities.SequenceEqual(b.Abnormalities);
            return same.Label("non-deterministic Karyotype for identical input");
        });
    }

    // ===================== DetectPloidy — INV-04 + formula oracle =====================

    /// <summary>
    /// INV-04 + §4.2 formula: production <c>DetectPloidy</c> equals the independent oracle exactly
    /// (PloidyLevel identical, Confidence within 1e-9), reproducing TRUE-median + banker's-rounding
    /// arithmetic. Also confirms the range guarantees <c>PloidyLevel ∈ [1,8]</c> and
    /// <c>Confidence ∈ [0,1]</c> for non-empty input.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DetectPloidy_MatchesFormulaOracle()
    {
        return Prop.ForAll(DepthsArbitrary(), t =>
        {
            var (depths, expectedDiploid) = t;
            var (ploidy, conf) = ChromosomeAnalyzer.DetectPloidy(depths, expectedDiploid);
            var (oPloidy, oConf) = ExpectedPloidy(depths, expectedDiploid);

            bool matches = ploidy == oPloidy && Math.Abs(conf - oConf) < 1e-9;
            bool ranges = ploidy is >= 1 and <= 8 && conf is >= 0.0 and <= 1.0 + 1e-12;
            return (matches && ranges)
                .Label($"depths(n={depths.Count}), diploid={expectedDiploid}: " +
                       $"got ({ploidy},{conf}), oracle ({oPloidy},{oConf})");
        });
    }

    /// <summary>
    /// INV-04 banker's-rounding witness: an EVEN number of depths whose median lands exactly on a
    /// rounding half-point. With median 1.25 and diploid 1.0, <c>ratio*2 = 2.5</c>; .NET's
    /// round-half-to-even gives ploidy 2 (NOT 3, which arithmetic floor(x+0.5) would give). Confirms
    /// the oracle and production both use <see cref="MidpointRounding.ToEven"/>.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DetectPloidy_HalfwayMedian_UsesBankersRounding()
    {
        var depths = new List<double> { 1.0, 1.5 }; // even count → median 1.25 → ratio*2 = 2.5
        var (ploidy, _) = ChromosomeAnalyzer.DetectPloidy(depths, expectedDiploidDepth: 1.0);
        Assert.That(ploidy, Is.EqualTo(2), "2.5 rounds to even (2), not up to 3");
        Assert.That(ExpectedPloidy(depths, 1.0).PloidyLevel, Is.EqualTo(2), "oracle agrees");
    }

    /// <summary>
    /// §4.2/§6.1 boundary anchors: ratio 0.5 ⇒ ploidy 1, conf 1; ratio 1.0 ⇒ ploidy 2, conf 1;
    /// ratio 2.0 ⇒ ploidy 4, conf 1; an extreme ratio (depth 10) ⇒ ploidy clamped to 8 with
    /// confidence 0. Each uses a single-element depth list so the median is the value itself.
    /// </summary>
    [TestCase(0.5, 1, 1.0, TestName = "ratio 0.5 → ploidy 1, conf 1")]
    [TestCase(1.0, 2, 1.0, TestName = "ratio 1.0 → ploidy 2, conf 1")]
    [TestCase(2.0, 4, 1.0, TestName = "ratio 2.0 → ploidy 4, conf 1")]
    [TestCase(10.0, 8, 0.0, TestName = "ratio 10 → ploidy clamped 8, conf 0")]
    [Category("Property")]
    public void DetectPloidy_BoundaryRatios(double depth, int expectedPloidy, double expectedConf)
    {
        var (ploidy, conf) = ChromosomeAnalyzer.DetectPloidy(new[] { depth }, expectedDiploidDepth: 1.0);
        Assert.Multiple(() =>
        {
            Assert.That(ploidy, Is.EqualTo(expectedPloidy));
            Assert.That(conf, Is.EqualTo(expectedConf).Within(1e-9));
        });
    }

    /// <summary>Edge (§6.1): empty depth list ⇒ exactly (2, 0).</summary>
    [Test]
    [Category("Property")]
    public void DetectPloidy_Empty_ReturnsDiploidZeroConfidence()
    {
        var (ploidy, conf) = ChromosomeAnalyzer.DetectPloidy(Enumerable.Empty<double>());
        Assert.Multiple(() =>
        {
            Assert.That(ploidy, Is.EqualTo(2));
            Assert.That(conf, Is.Zero);
        });
    }

    /// <summary>D (determinism): identical depth inputs ⇒ identical <c>(PloidyLevel, Confidence)</c>.</summary>
    [FsCheck.NUnit.Property]
    public Property DetectPloidy_IsDeterministic()
    {
        return Prop.ForAll(DepthsArbitrary(), t =>
        {
            var (depths, expectedDiploid) = t;
            var a = ChromosomeAnalyzer.DetectPloidy(depths, expectedDiploid);
            var b = ChromosomeAnalyzer.DetectPloidy(depths, expectedDiploid);
            return (a == b).Label($"non-deterministic DetectPloidy for n={depths.Count}");
        });
    }

    #endregion

    #region CHROM-ANEU-001

    // Covers ChromosomeAnalyzer.DetectAneuploidy and ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy
    // (Aneuploidy_Detection.md §2.4 INV-01..INV-03, §4.2 per-bin formula + ratio table +
    // whole-chromosome label map, §3 contract, §6 edge cases).
    //
    //   R  copy number ≥ 0 (INV-01: CopyNumber ∈ [0,10]); confidence ∈ [0,1] (INV-02)
    //   M  higher depth → higher CN (CopyNumber = clamp(round(2·mean/median),0,10) is monotone non-decreasing in mean)
    //   D  deterministic (identical inputs → identical output sequences for both methods)
    //
    // Key algebraic simplification used by the oracles: the source computes
    //   logRatio   = Math.Log2(meanDepth / medianDepth)
    //   copyNumber = round(2 ^ logRatio · 2)
    // and since 2 ^ log2(mean/median) = mean/median exactly, this is round(2 · mean/median).
    // observed = 2 ^ logRatio = mean/median, expected = copyNumber/2.0.
    // All oracles below are transcribed INDEPENDENTLY from the doc/source semantics and are
    // NEVER routed through production.

    // ===================== Independent oracles =====================

    /// <summary>
    /// Independent per-bin oracle for one bin, transcribed verbatim from Aneuploidy_Detection.md
    /// §4.2 / source: <c>observed = meanDepth / medianDepth</c>;
    /// <c>logRatio = Math.Log2(observed)</c>;
    /// <c>copyNumber = clamp((int)Math.Round(2 · observed), 0, 10)</c> using .NET banker's rounding
    /// (<see cref="MidpointRounding.ToEven"/>, the default of <see cref="Math.Round(double)"/>);
    /// <c>expected = copyNumber / 2.0</c>;
    /// <c>confidence = 1 - min(1, |expected - observed|)</c>;
    /// <c>Start = binIndex · binSize</c>, <c>End = (binIndex + 1) · binSize - 1</c>.
    /// </summary>
    private static ChromosomeAnalyzer.CopyNumberState ExpectedBinState(
        string chromosome, int binIndex, int binSize, double meanDepth, double medianDepth)
    {
        double observed = meanDepth / medianDepth;
        double logRatio = Math.Log2(observed);
        int copyNumber = (int)Math.Round(observed * 2.0); // round-half-to-even, matching source
        copyNumber = Math.Max(0, Math.Min(10, copyNumber));
        double expected = copyNumber / 2.0;
        double confidence = 1.0 - Math.Min(1.0, Math.Abs(expected - observed));
        return new ChromosomeAnalyzer.CopyNumberState(
            chromosome,
            binIndex * binSize,
            (binIndex + 1) * binSize - 1,
            copyNumber,
            logRatio,
            confidence);
    }

    /// <summary>
    /// Independent <c>DetectAneuploidy</c> oracle: replicate grouping by chromosome then by
    /// <c>Position / binSize</c>, average depths inside each bin, and emit one bin state per
    /// chromosome/bin in ascending bin order. Returns a flat list grouped per chromosome.
    /// </summary>
    private static List<ChromosomeAnalyzer.CopyNumberState> ExpectedAneuploidy(
        IReadOnlyList<(string Chromosome, int Position, double Depth)> data,
        double medianDepth,
        int binSize)
    {
        var result = new List<ChromosomeAnalyzer.CopyNumberState>();
        if (data.Count == 0 || medianDepth <= 0)
            return result;

        foreach (var chromGroup in data.GroupBy(d => d.Chromosome))
        {
            var bins = chromGroup
                .GroupBy(d => d.Position / binSize)
                .OrderBy(g => g.Key);
            foreach (var bin in bins)
            {
                double mean = bin.Average(d => d.Depth);
                result.Add(ExpectedBinState(chromGroup.Key, bin.Key, binSize, mean, medianDepth));
            }
        }
        return result;
    }

    /// <summary>
    /// Independent cytogenetic label oracle (§4.2 whole-chromosome map):
    /// 0→Nullisomy, 1→Monosomy, 3→Trisomy, 4→Tetrasomy, 5→Pentasomy, else <c>Copy number = {cn}</c>.
    /// </summary>
    private static string ExpectedWholeChromType(int cn) => cn switch
    {
        0 => "Nullisomy",
        1 => "Monosomy",
        3 => "Trisomy",
        4 => "Tetrasomy",
        5 => "Pentasomy",
        _ => $"Copy number = {cn}"
    };

    /// <summary>
    /// Independent <c>IdentifyWholeChromosomeAneuploidy</c> oracle (§2.4 INV-03): per chromosome,
    /// the dominant copy number is the one occupying the largest bin fraction; emit
    /// <c>(chrom, dominantCN, type)</c> only when <c>dominantFraction ≥ minFraction</c> AND
    /// <c>dominantCN != 2</c>. Returned as a set keyed by chromosome (production order is
    /// GroupBy-dependent and not part of the contract). Generators below avoid fraction ties so
    /// tie-break order is irrelevant.
    /// </summary>
    private static HashSet<(string Chromosome, int CopyNumber, string Type)> ExpectedWholeChrom(
        IReadOnlyList<ChromosomeAnalyzer.CopyNumberState> states, double minFraction)
    {
        var result = new HashSet<(string, int, string)>();
        foreach (var chromGroup in states.GroupBy(s => s.Chromosome))
        {
            int total = chromGroup.Count();
            var dominant = chromGroup
                .GroupBy(s => s.CopyNumber)
                .Select(g => (CopyNumber: g.Key, Fraction: g.Count() / (double)total))
                .OrderByDescending(g => g.Fraction)
                .First();
            if (dominant.Fraction >= minFraction && dominant.CopyNumber != 2)
                result.Add((chromGroup.Key, dominant.CopyNumber, ExpectedWholeChromType(dominant.CopyNumber)));
        }
        return result;
    }

    // ===================== Generators =====================

    /// <summary>
    /// Generates depth observations grouped into a few chromosomes and bins, plus a positive
    /// <c>medianDepth</c> and a modest positive <c>binSize</c>. Positions stay small (a handful of
    /// bins) so <c>binIndex · binSize</c> never overflows <c>int</c>; depths are finite positives.
    /// Yields (observations, medianDepth, binSize).
    /// </summary>
    private static Arbitrary<(List<(string Chromosome, int Position, double Depth)> data, double median, int binSize)>
        DepthDataArbitrary() =>
        (from binSize in Gen.Choose(10, 1000)
         from chromCount in Gen.Choose(1, 3)
         from n in Gen.Choose(1, 40)
         from median in Gen.Choose(1, 8000).Select(i => i / 100.0)
         from obs in
             (from chromIdx in Gen.Choose(0, chromCount - 1)
              from binIdx in Gen.Choose(0, 4)
              from offset in Gen.Choose(0, binSize - 1)
              from depthMilli in Gen.Choose(1, 2_000_000)
              select ($"chr{chromIdx}", binIdx * binSize + offset, depthMilli / 1000.0)).ListOf(n)
         select (obs, median, binSize)).ToArbitrary();

    /// <summary>
    /// Generates a <c>CopyNumberState</c> list for one or more chromosomes, each with a CLEAR
    /// dominant copy number: a controllable number of "dominant" bins of a chosen CN plus strictly
    /// fewer "noise" bins of a different CN, so no fraction tie can occur. The dominant CN and its
    /// fraction span disomy (⇒ nothing emitted), sub-threshold dominance (⇒ nothing emitted), and
    /// clearly-dominant non-disomic states (⇒ emitted). The non-state fields (Start/End/LogRatio/
    /// Confidence) are irrelevant to whole-chromosome classification and set to fixed dummies.
    /// Yields (states, minFraction).
    /// </summary>
    private static Arbitrary<(List<ChromosomeAnalyzer.CopyNumberState> states, double minFraction)>
        WholeChromArbitrary() =>
        (from chromCount in Gen.Choose(1, 3)
         from minFraction in Gen.Choose(50, 100).Select(i => i / 100.0)
         from perChrom in
             (from dominantCn in Gen.Choose(0, 6)
              from noiseDelta in Gen.Choose(1, 4)
              from dominantBins in Gen.Choose(2, 10)
              from noiseBins in Gen.Choose(0, 1) // strictly fewer ⇒ no tie; 0 allows a pure chromosome
              select (dominantCn, noiseCn: (dominantCn + noiseDelta) % 11, dominantBins, noiseBins))
             .ListOf(chromCount)
         let states = perChrom
             .SelectMany((p, ci) =>
                 Enumerable.Repeat(p.dominantCn, p.dominantBins)
                     .Concat(Enumerable.Repeat(p.noiseCn, p.noiseBins))
                     .Select((cn, bi) => new ChromosomeAnalyzer.CopyNumberState(
                         $"chr{ci}", bi * 1000, (bi + 1) * 1000 - 1, cn, 0.0, 0.5)))
             .ToList()
         select (states, minFraction)).ToArbitrary();

    // ===================== DetectAneuploidy — per-bin field oracle =====================

    /// <summary>
    /// §4.2 per-bin oracle (INV-01/INV-02 + coordinates): for arbitrary grouped depth data the
    /// production <c>DetectAneuploidy</c> output equals the independent oracle field-by-field —
    /// same chromosome ordering, same ascending bin order, integer fields exact, doubles
    /// (LogRatio, Confidence) within 1e-9. The oracle recomputes <c>CopyNumber =
    /// clamp(round(2·mean/median),0,10)</c>, <c>LogRatio = log2(mean/median)</c>, and
    /// <c>Confidence = 1 - min(1, |CopyNumber/2 - mean/median|)</c> independently of the code.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DetectAneuploidy_PerBin_MatchesOracle()
    {
        return Prop.ForAll(DepthDataArbitrary(), t =>
        {
            var (data, median, binSize) = t;
            var actual = ChromosomeAnalyzer.DetectAneuploidy(data, median, binSize).ToList();
            var expected = ExpectedAneuploidy(data, median, binSize);

            if (actual.Count != expected.Count)
                return false.Label($"count {actual.Count} != oracle {expected.Count}");

            for (int i = 0; i < actual.Count; i++)
            {
                var a = actual[i];
                var e = expected[i];
                bool same = a.Chromosome == e.Chromosome
                            && a.Start == e.Start
                            && a.End == e.End
                            && a.CopyNumber == e.CopyNumber
                            && Math.Abs(a.LogRatio - e.LogRatio) < 1e-9
                            && Math.Abs(a.Confidence - e.Confidence) < 1e-9;
                if (!same)
                    return false.Label($"bin {i}: got {a}, oracle {e}");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-01 (R: copy number ≥ 0) + INV-02 + valid coordinates: for arbitrary inputs every emitted
    /// state has <c>CopyNumber ∈ [0,10]</c>, <c>Confidence ∈ [0,1]</c>, <c>Start ≤ End</c>, and the
    /// bin coordinates satisfy <c>End = Start + binSize - 1</c> with <c>Start</c> a non-negative
    /// multiple of <c>binSize</c>.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DetectAneuploidy_Invariants_AlwaysHold()
    {
        return Prop.ForAll(DepthDataArbitrary(), t =>
        {
            var (data, median, binSize) = t;
            foreach (var s in ChromosomeAnalyzer.DetectAneuploidy(data, median, binSize))
            {
                bool cnOk = s.CopyNumber is >= 0 and <= 10;                       // INV-01
                bool confOk = s.Confidence is >= 0.0 and <= 1.0 + 1e-12;          // INV-02
                bool coordOk = s.Start <= s.End
                               && s.End == s.Start + binSize - 1
                               && s.Start >= 0
                               && s.Start % binSize == 0;
                if (!(cnOk && confOk && coordOk))
                    return false.Label($"invariant violated for {s} (binSize {binSize})");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// §4.2 ratio-table anchors (observed-ratio ⇒ copy number): a single bin whose mean depth is
    /// <c>ratio · median</c> yields <c>CopyNumber = round(2·ratio)</c>. Confirms 0.5⇒1, 1.0⇒2,
    /// 1.5⇒3, 2.0⇒4, 2.5⇒5 exactly. Each uses one observation so the bin mean is the value itself.
    /// </summary>
    [TestCase(0.5, 1, TestName = "ratio 0.5 → CN 1")]
    [TestCase(1.0, 2, TestName = "ratio 1.0 → CN 2")]
    [TestCase(1.5, 3, TestName = "ratio 1.5 → CN 3")]
    [TestCase(2.0, 4, TestName = "ratio 2.0 → CN 4")]
    [TestCase(2.5, 5, TestName = "ratio 2.5 → CN 5")]
    [Category("Property")]
    public void DetectAneuploidy_RatioTable_Anchors(double ratio, int expectedCn)
    {
        const double median = 30.0;
        var data = new[] { ("chr1", 0, ratio * median) };
        var states = ChromosomeAnalyzer.DetectAneuploidy(data, median, binSize: 1000).ToList();
        Assert.That(states, Has.Count.EqualTo(1));
        Assert.That(states[0].CopyNumber, Is.EqualTo(expectedCn),
            $"observed ratio {ratio} must map to copy number {expectedCn}");
    }

    // ===================== M: higher depth → higher CN (central business invariant) =====================

    /// <summary>
    /// M (monotonicity, central business invariant): with chromosome, position, median and binSize
    /// fixed, feeding a single bin with a lower mean depth never yields a HIGHER copy number than a
    /// higher mean depth — <c>CN(lo) ≤ CN(hi)</c> whenever <c>depthLo ≤ depthHi</c>. This holds
    /// because <c>CopyNumber = clamp(round(2·mean/median),0,10)</c> is monotone non-decreasing in
    /// the mean, and the clamp/round preserve order.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DetectAneuploidy_HigherDepth_NeverLowerCopyNumber()
    {
        var arb = (from median in Gen.Choose(1, 5000).Select(i => i / 100.0)
                   from dLo in Gen.Choose(1, 2_000_000).Select(i => i / 1000.0)
                   from dHi in Gen.Choose(1, 2_000_000).Select(i => i / 1000.0)
                   select (median, lo: Math.Min(dLo, dHi), hi: Math.Max(dLo, dHi))).ToArbitrary();
        return Prop.ForAll(arb, t =>
        {
            var (median, lo, hi) = t;
            int cnLo = ChromosomeAnalyzer.DetectAneuploidy(new[] { ("chr1", 0, lo) }, median, 1000)
                .Single().CopyNumber;
            int cnHi = ChromosomeAnalyzer.DetectAneuploidy(new[] { ("chr1", 0, hi) }, median, 1000)
                .Single().CopyNumber;
            return (cnLo <= cnHi)
                .Label($"median={median}: depth {lo}→CN {cnLo} but depth {hi}→CN {cnHi} (must be non-decreasing)");
        });
    }

    // ===================== DetectAneuploidy — edge cases =====================

    /// <summary>Edge (§6.1): empty depth data ⇒ no output, for any positive median/binSize.</summary>
    [Test]
    [Category("Property")]
    public void DetectAneuploidy_EmptyInput_NoOutput()
    {
        var states = ChromosomeAnalyzer.DetectAneuploidy(
            Enumerable.Empty<(string, int, double)>(), medianDepth: 30.0, binSize: 1000);
        Assert.That(states, Is.Empty);
    }

    /// <summary>Edge (§6.1): a non-positive <c>medianDepth</c> (0 or negative) ⇒ no output.</summary>
    [TestCase(0.0, TestName = "median 0 → no output")]
    [TestCase(-5.0, TestName = "median negative → no output")]
    [Category("Property")]
    public void DetectAneuploidy_NonPositiveMedian_NoOutput(double median)
    {
        var data = new[] { ("chr1", 0, 30.0), ("chr1", 100, 60.0) };
        var states = ChromosomeAnalyzer.DetectAneuploidy(data, median, binSize: 1000);
        Assert.That(states, Is.Empty);
    }

    // ===================== IdentifyWholeChromosomeAneuploidy — INV-03 oracle =====================

    /// <summary>
    /// INV-03 (§2.4 / §4.2): production <c>IdentifyWholeChromosomeAneuploidy</c> output equals the
    /// independent oracle as a set — a whole-chromosome call is emitted exactly for chromosomes whose
    /// dominant copy number occupies <c>≥ minFraction</c> of bins AND is non-disomic, with the label
    /// recomputed by the independent map. Generators guarantee a clear dominant CN (no fraction
    /// ties), and the random fraction/dominant-CN spans the suppressed cases (disomy dominant;
    /// dominant fraction below threshold) and the emitted cases.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IdentifyWholeChromosome_MatchesOracle()
    {
        return Prop.ForAll(WholeChromArbitrary(), t =>
        {
            var (states, minFraction) = t;
            var actual = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states, minFraction).ToHashSet();
            var expected = ExpectedWholeChrom(states, minFraction);
            return actual.SetEquals(expected)
                .Label($"minFraction={minFraction}: got [{string.Join("|", actual)}] vs oracle [{string.Join("|", expected)}]");
        });
    }

    /// <summary>
    /// INV-03 suppression anchor (dominant disomy): a chromosome whose dominant state is disomic
    /// (CopyNumber 2) emits NOTHING even when fully dominant, because the rule excludes
    /// <c>CopyNumber == 2</c>.
    /// </summary>
    [Test]
    [Category("Property")]
    public void IdentifyWholeChromosome_DominantDisomy_EmitsNothing()
    {
        var states = Enumerable.Range(0, 10)
            .Select(i => new ChromosomeAnalyzer.CopyNumberState("chr1", i * 1000, (i + 1) * 1000 - 1, 2, 0.0, 1.0))
            .ToList();
        var calls = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states, minFraction: 0.8);
        Assert.That(calls, Is.Empty, "a dominant disomic chromosome is never reported");
    }

    /// <summary>
    /// INV-03 suppression anchor (sub-threshold dominance): a chromosome whose dominant non-disomic
    /// state occupies only 6/10 = 0.6 of bins is NOT reported when <c>minFraction = 0.8</c>.
    /// </summary>
    [Test]
    [Category("Property")]
    public void IdentifyWholeChromosome_BelowMinFraction_EmitsNothing()
    {
        var states = new List<ChromosomeAnalyzer.CopyNumberState>();
        for (int i = 0; i < 6; i++)
            states.Add(new ChromosomeAnalyzer.CopyNumberState("chr1", i * 1000, (i + 1) * 1000 - 1, 3, 0.0, 1.0));
        for (int i = 6; i < 10; i++)
            states.Add(new ChromosomeAnalyzer.CopyNumberState("chr1", i * 1000, (i + 1) * 1000 - 1, 2, 0.0, 1.0));
        var calls = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states, minFraction: 0.8);
        Assert.That(calls, Is.Empty, "dominant fraction 0.6 < 0.8 ⇒ no whole-chromosome call");
    }

    /// <summary>
    /// INV-03 label anchors (§4.2 map): a fully dominant chromosome of each non-disomic copy number
    /// is reported with the exact cytogenetic label — 0 Nullisomy, 1 Monosomy, 3 Trisomy, 4
    /// Tetrasomy, 5 Pentasomy, and a fallback <c>Copy number = 6</c> for other non-disomic values.
    /// </summary>
    [TestCase(0, "Nullisomy", TestName = "CN 0 → Nullisomy")]
    [TestCase(1, "Monosomy", TestName = "CN 1 → Monosomy")]
    [TestCase(3, "Trisomy", TestName = "CN 3 → Trisomy")]
    [TestCase(4, "Tetrasomy", TestName = "CN 4 → Tetrasomy")]
    [TestCase(5, "Pentasomy", TestName = "CN 5 → Pentasomy")]
    [TestCase(6, "Copy number = 6", TestName = "CN 6 → fallback label")]
    [Category("Property")]
    public void IdentifyWholeChromosome_LabelAnchors(int cn, string expectedType)
    {
        var states = Enumerable.Range(0, 10)
            .Select(i => new ChromosomeAnalyzer.CopyNumberState("chr1", i * 1000, (i + 1) * 1000 - 1, cn, 0.0, 1.0))
            .ToList();
        var calls = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states, minFraction: 0.8).ToList();
        Assert.That(calls, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(calls[0].Chromosome, Is.EqualTo("chr1"));
            Assert.That(calls[0].CopyNumber, Is.EqualTo(cn));
            Assert.That(calls[0].Type, Is.EqualTo(expectedType));
        });
    }

    /// <summary>Edge: empty state input ⇒ no whole-chromosome output.</summary>
    [Test]
    [Category("Property")]
    public void IdentifyWholeChromosome_EmptyInput_NoOutput()
    {
        var calls = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(
            Enumerable.Empty<ChromosomeAnalyzer.CopyNumberState>());
        Assert.That(calls, Is.Empty);
    }

    // ===================== D: determinism =====================

    /// <summary>
    /// D (determinism): identical depth inputs ⇒ identical <c>DetectAneuploidy</c> output sequences
    /// (field-by-field, in order).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DetectAneuploidy_IsDeterministic()
    {
        return Prop.ForAll(DepthDataArbitrary(), t =>
        {
            var (data, median, binSize) = t;
            var a = ChromosomeAnalyzer.DetectAneuploidy(data, median, binSize).ToList();
            var b = ChromosomeAnalyzer.DetectAneuploidy(data, median, binSize).ToList();
            return a.SequenceEqual(b).Label($"non-deterministic DetectAneuploidy (n={data.Count}, binSize={binSize})");
        });
    }

    /// <summary>
    /// D (determinism): identical state inputs ⇒ identical <c>IdentifyWholeChromosomeAneuploidy</c>
    /// output sequences.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property IdentifyWholeChromosome_IsDeterministic()
    {
        return Prop.ForAll(WholeChromArbitrary(), t =>
        {
            var (states, minFraction) = t;
            var a = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states, minFraction).ToList();
            var b = ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(states, minFraction).ToList();
            return a.SequenceEqual(b).Label($"non-deterministic IdentifyWholeChromosomeAneuploidy (n={states.Count})");
        });
    }

    #endregion

    #region CHROM-SYNT-001

    // Covers ChromosomeAnalyzer.FindSyntenyBlocks and ChromosomeAnalyzer.DetectRearrangements
    // (Synteny_Analysis.md §2.4 INV-01..INV-04, §2.1/§4.2 rearrangement signature table,
    //  §3 contract — maxGap is in MEGABASES (×1,000,000), minGenes default 3 — §6 edges, §7.1 example).
    //
    //   R  block positions valid (INV-02): Species1Start ≤ Species1End, Species2Start ≤ Species2End,
    //      GeneCount ≥ minGenes; plus INV-01 (Strand ∈ {'+','-'}), INV-03 (SequenceIdentity is NaN).
    //   S  role-swap transposition symmetry on FORWARD-collinear input: swapping (genome1↔genome2)
    //      roles in the ortholog tuples transposes the single resulting block. NOT asserted for the
    //      general/reverse coordinate heuristic, which is documented as not globally symmetric.
    //   D  determinism: identical input → identical block / rearrangement sequences.
    //
    // The forward-block CONTENT oracle recomputes the expected block independently from the
    // generated run (first/last gene coordinates), proving the exact block contract — not echoing
    // the implementation. INV-04 is proven by the rearrangement Type ∈ the four-label set, with
    // hand-built inversion/translocation anchors. maxGap megabase semantics is proven by a split.

    // Ortholog tuple type used by FindSyntenyBlocks.
    private record struct Ortholog(
        string Chr1, int Start1, int End1, string Gene1,
        string Chr2, int Start2, int End2, string Gene2);

    /// <summary>One MB in base pairs — the factor <c>FindSyntenyBlocks</c> multiplies <c>maxGap</c> by.</summary>
    private const int OneMegabase = 1_000_000;

    /// <summary>The four valid rearrangement <c>Type</c> labels per INV-04 / §2.1.</summary>
    private static readonly string[] ValidRearrangementTypes =
        { "Inversion", "Translocation", "Deletion", "Duplication" };

    // ===================== Helpers (independent constructions) =====================

    /// <summary>
    /// Converts <see cref="Ortholog"/> records to the value-tuple shape expected by
    /// <c>FindSyntenyBlocks</c>.
    /// </summary>
    private static IEnumerable<(string Chr1, int Start1, int End1, string Gene1,
        string Chr2, int Start2, int End2, string Gene2)> ToTuples(IEnumerable<Ortholog> orthologs) =>
        orthologs.Select(o => (o.Chr1, o.Start1, o.End1, o.Gene1, o.Chr2, o.Start2, o.End2, o.Gene2));

    /// <summary>
    /// Swaps the genome roles of an ortholog list exactly as the symmetry contract specifies:
    /// <c>(Chr1,Start1,End1,Gene1,Chr2,Start2,End2,Gene2) → (Chr2,Start2,End2,Gene2,Chr1,Start1,End1,Gene1)</c>.
    /// </summary>
    private static List<Ortholog> SwapRoles(IEnumerable<Ortholog> orthologs) =>
        orthologs.Select(o => new Ortholog(
            o.Chr2, o.Start2, o.End2, o.Gene2,
            o.Chr1, o.Start1, o.End1, o.Gene1)).ToList();

    /// <summary>
    /// Treats two <see cref="ChromosomeAnalyzer.SyntenyBlock"/> values as equal with NaN==NaN for
    /// <c>SequenceIdentity</c> (record-struct equality already treats NaN==NaN via bit-equality, but
    /// this makes the determinism intent explicit and independent).
    /// </summary>
    private static bool BlocksEqual(
        ChromosomeAnalyzer.SyntenyBlock a, ChromosomeAnalyzer.SyntenyBlock b) =>
        a.Species1Chromosome == b.Species1Chromosome
        && a.Species1Start == b.Species1Start
        && a.Species1End == b.Species1End
        && a.Species2Chromosome == b.Species2Chromosome
        && a.Species2Start == b.Species2Start
        && a.Species2End == b.Species2End
        && a.Strand == b.Strand
        && a.GeneCount == b.GeneCount
        && (double.IsNaN(a.SequenceIdentity)
            ? double.IsNaN(b.SequenceIdentity)
            : a.SequenceIdentity == b.SequenceIdentity);

    // ===================== Generators =====================

    /// <summary>
    /// Generates one clean forward-collinear ortholog run on a single chromosome pair:
    /// <c>n ∈ [minGenes, 8]</c> genes strictly ascending in BOTH genomes, with every inter-gene gap
    /// well inside <c>maxGap·1e6</c> (spacing in the low thousands of bp, maxGap=10 ⇒ 10 Mb ceiling).
    /// Forward orientation is guaranteed because each <c>curr.Start2 &gt; prev.End2</c>. Coordinates
    /// stay modest (≤ ~100k) to avoid overflow when the implementation computes <c>maxGap·1e6</c>.
    /// Yields the ortholog list (the single expected block is recomputed independently by the oracle).
    /// </summary>
    private static Arbitrary<List<Ortholog>> ForwardCollinearRunArbitrary() =>
        (from n in Gen.Choose(3, 8)
         from geneLen in Gen.Choose(500, 2000)
         from gap in Gen.Choose(100, 5000) // « maxGap·1e6 (=10 Mb) so the run never splits
         from start1 in Gen.Choose(1_000, 50_000)
         from start2 in Gen.Choose(1_000, 50_000)
         let step = geneLen + gap
         let genes = Enumerable.Range(0, n).Select(i => new Ortholog(
             "chr1", start1 + i * step, start1 + i * step + geneLen, $"g1_{i}",
             "chrA", start2 + i * step, start2 + i * step + geneLen, $"gA_{i}")).ToList()
         select genes).ToArbitrary();

    // ===================== R / INV-01 / INV-02 / INV-03 — positions valid =====================

    /// <summary>
    /// INV-02 (R: positions valid) + INV-01 + INV-03 over ARBITRARY ortholog input: for EVERY emitted
    /// block <c>Species1Start ≤ Species1End</c>, <c>Species2Start ≤ Species2End</c>,
    /// <c>GeneCount ≥ minGenes</c> (INV-02); <c>Strand ∈ {'+','-'}</c> (INV-01); and
    /// <c>SequenceIdentity</c> is NaN (INV-03). The generator mixes chromosomes, forward and reverse
    /// orders, and gaps so the heuristic is exercised broadly — not just on clean runs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindSyntenyBlocks_EmittedBlocks_AlwaysValid()
    {
        var arb =
            (from n in Gen.Choose(0, 12)
             from rows in (from chrIdx in Gen.Choose(0, 1)
                           from s1 in Gen.Choose(1_000, 80_000)
                           from len1 in Gen.Choose(200, 3_000)
                           from s2 in Gen.Choose(1_000, 80_000)
                           from len2 in Gen.Choose(200, 3_000)
                           select new Ortholog(
                               $"chr{chrIdx}", s1, s1 + len1, "g",
                               "chrA", s2, s2 + len2, "gA")).ListOf(n)
             select rows).ToArbitrary();

        const int minGenes = 3;
        return Prop.ForAll(arb, rows =>
        {
            var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(ToTuples(rows), minGenes: minGenes).ToList();
            foreach (var b in blocks)
            {
                bool inv02 = b.Species1Start <= b.Species1End
                             && b.Species2Start <= b.Species2End
                             && b.GeneCount >= minGenes;
                bool inv01 = b.Strand is '+' or '-';
                bool inv03 = double.IsNaN(b.SequenceIdentity);
                if (!(inv02 && inv01 && inv03))
                    return false.Label(
                        $"invalid block: [{b.Species1Start},{b.Species1End}]→[{b.Species2Start},{b.Species2End}] " +
                        $"strand={b.Strand} genes={b.GeneCount} id={b.SequenceIdentity}");
            }
            return true.ToProperty();
        });
    }

    // ===================== Forward-block content oracle (exact contract) =====================

    /// <summary>
    /// Forward block CONTENT oracle: a single clean forward-collinear run of <c>n ≥ minGenes</c> genes
    /// (strictly ascending in both genomes, gaps ≤ maxGap·1e6) produces EXACTLY ONE block whose fields
    /// equal those recomputed independently from the run —
    /// <c>Species1Start = first.Start1</c>, <c>Species1End = last.End1</c>,
    /// <c>Species2Start = first.Start2</c>, <c>Species2End = last.End2</c> (ascending ⇒ first=min,
    /// last=max), <c>Strand = '+'</c>, <c>GeneCount = n</c>, <c>SequenceIdentity = NaN</c>.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindSyntenyBlocks_ForwardCollinearRun_MatchesContentOracle()
    {
        return Prop.ForAll(ForwardCollinearRunArbitrary(), genes =>
        {
            var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(ToTuples(genes), minGenes: 3, maxGap: 10).ToList();
            if (blocks.Count != 1)
                return false.Label($"expected exactly 1 block for a clean run of {genes.Count}, got {blocks.Count}");

            var b = blocks[0];
            var first = genes[0];
            var last = genes[^1];
            bool match = b.Species1Chromosome == "chr1"
                         && b.Species1Start == first.Start1
                         && b.Species1End == last.End1
                         && b.Species2Chromosome == "chrA"
                         && b.Species2Start == first.Start2
                         && b.Species2End == last.End2
                         && b.Strand == '+'
                         && b.GeneCount == genes.Count
                         && double.IsNaN(b.SequenceIdentity);
            return match.Label(
                $"block [{b.Species1Start},{b.Species1End}]→[{b.Species2Start},{b.Species2End}] " +
                $"strand={b.Strand} genes={b.GeneCount} vs run of {genes.Count} " +
                $"([{first.Start1},{last.End1}]→[{first.Start2},{last.End2}])");
        });
    }

    /// <summary>
    /// Exact forward anchor — §7.1 worked example: three collinear orthologs on chr1→chrA produce one
    /// forward block spanning the run with GeneCount 3 and NaN identity.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindSyntenyBlocks_WorkedExample_SingleForwardBlock()
    {
        var pairs = new List<(string, int, int, string, string, int, int, string)>
        {
            ("chr1", 1000, 2000, "gene1", "chrA", 1000, 2000, "geneA"),
            ("chr1", 3000, 4000, "gene2", "chrA", 3000, 4000, "geneB"),
            ("chr1", 5000, 6000, "gene3", "chrA", 5000, 6000, "geneC"),
        };
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(pairs, minGenes: 3, maxGap: 10).ToList();
        Assert.That(blocks, Has.Count.EqualTo(1));
        var b = blocks[0];
        Assert.Multiple(() =>
        {
            Assert.That(b.Species1Chromosome, Is.EqualTo("chr1"));
            Assert.That(b.Species2Chromosome, Is.EqualTo("chrA"));
            Assert.That(b.Species1Start, Is.EqualTo(1000));
            Assert.That(b.Species1End, Is.EqualTo(6000));
            Assert.That(b.Species2Start, Is.EqualTo(1000));
            Assert.That(b.Species2End, Is.EqualTo(6000));
            Assert.That(b.Strand, Is.EqualTo('+'));
            Assert.That(b.GeneCount, Is.EqualTo(3));
            Assert.That(double.IsNaN(b.SequenceIdentity), Is.True);
        });
    }

    /// <summary>
    /// Reverse block (INV-01 '-' branch): genes ascending in genome 1 but strictly DESCENDING in
    /// genome 2 form a reverse-collinear run ⇒ a single block with <c>Strand = '-'</c>. Boundaries
    /// use min/max so Species2Start/End remain ordered (INV-02). Hand-built anchor.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindSyntenyBlocks_ReverseCollinearRun_NegativeStrand()
    {
        // genome1 ascending; genome2 descending (each curr.Start2 < prev.End2 ⇒ not forward).
        var pairs = new List<(string, int, int, string, string, int, int, string)>
        {
            ("chr1", 1000, 2000, "g1", "chrA", 9000, 10000, "gA"),
            ("chr1", 3000, 4000, "g2", "chrA", 6000, 7000, "gB"),
            ("chr1", 5000, 6000, "g3", "chrA", 3000, 4000, "gC"),
        };
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(pairs, minGenes: 3, maxGap: 10).ToList();
        Assert.That(blocks, Has.Count.EqualTo(1));
        var b = blocks[0];
        Assert.Multiple(() =>
        {
            Assert.That(b.Strand, Is.EqualTo('-'), "descending genome-2 order ⇒ reverse strand");
            Assert.That(b.GeneCount, Is.EqualTo(3));
            Assert.That(b.Species1Start, Is.LessThanOrEqualTo(b.Species1End));
            Assert.That(b.Species2Start, Is.LessThanOrEqualTo(b.Species2End));
        });
    }

    // ===================== S — role-swap transposition symmetry =====================

    /// <summary>
    /// S (symmetry, provable case): for a clean forward-collinear single-group input, swapping the
    /// genome roles of every ortholog tuple TRANSPOSES the single resulting block —
    /// <c>swapped.Species1* == original.Species2*</c> and <c>swapped.Species2* == original.Species1*</c>
    /// (chromosome, start, end), with equal <c>GeneCount</c> and both blocks forward (<c>'+'</c>).
    /// Symmetry is asserted ONLY on this provable forward case, per the documented heuristic caveat.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindSyntenyBlocks_RoleSwap_TransposesForwardBlock()
    {
        return Prop.ForAll(ForwardCollinearRunArbitrary(), genes =>
        {
            var original = ChromosomeAnalyzer.FindSyntenyBlocks(ToTuples(genes), minGenes: 3, maxGap: 10).ToList();
            var swapped = ChromosomeAnalyzer.FindSyntenyBlocks(ToTuples(SwapRoles(genes)), minGenes: 3, maxGap: 10).ToList();

            if (original.Count != 1 || swapped.Count != 1)
                return false.Label($"expected one block each; original={original.Count}, swapped={swapped.Count}");

            var o = original[0];
            var s = swapped[0];
            bool transpose =
                s.Species1Chromosome == o.Species2Chromosome
                && s.Species1Start == o.Species2Start
                && s.Species1End == o.Species2End
                && s.Species2Chromosome == o.Species1Chromosome
                && s.Species2Start == o.Species1Start
                && s.Species2End == o.Species1End
                && s.GeneCount == o.GeneCount
                && o.Strand == '+' && s.Strand == '+';
            return transpose.Label(
                $"not a transpose: original [{o.Species1Start},{o.Species1End}]→[{o.Species2Start},{o.Species2End}] " +
                $"vs swapped [{s.Species1Start},{s.Species1End}]→[{s.Species2Start},{s.Species2End}]");
        });
    }

    // ===================== minGenes threshold & empty input (edges) =====================

    /// <summary>
    /// Edge (§6.1): fewer than <c>minGenes</c> total ortholog pairs ⇒ no blocks (the total count is
    /// checked before scanning). Driven over runs of length 0..minGenes-1.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindSyntenyBlocks_FewerThanMinGenes_NoBlocks()
    {
        var arb = (from n in Gen.Choose(0, 4)
                   select n).ToArbitrary();
        return Prop.ForAll(arb, n =>
        {
            const int minGenes = 5; // n ∈ [0,4] is always < minGenes
            var rows = Enumerable.Range(0, n).Select(i => new Ortholog(
                "chr1", 1000 + i * 4000, 2000 + i * 4000, $"g{i}",
                "chrA", 1000 + i * 4000, 2000 + i * 4000, $"gA{i}")).ToList();
            var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(ToTuples(rows), minGenes: minGenes).ToList();
            return (blocks.Count == 0).Label($"n={n} < minGenes={minGenes} produced {blocks.Count} blocks");
        });
    }

    /// <summary>
    /// Edge (§6.1): a collinear run STRICTLY shorter than <c>minGenes</c> (but with enough total pairs
    /// to pass the global count gate) is discarded. Two short runs of 2 genes each on different
    /// chromosome pairs total 4 ≥ minGenes=3, yet each run of 2 &lt; 3 ⇒ no blocks survive.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindSyntenyBlocks_RunShorterThanMinGenes_Discarded()
    {
        var pairs = new List<(string, int, int, string, string, int, int, string)>
        {
            ("chr1", 1000, 2000, "g1", "chrA", 1000, 2000, "gA"),
            ("chr1", 3000, 4000, "g2", "chrA", 3000, 4000, "gB"),
            ("chr2", 1000, 2000, "g3", "chrB", 1000, 2000, "gC"),
            ("chr2", 3000, 4000, "g4", "chrB", 3000, 4000, "gD"),
        };
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(pairs, minGenes: 3, maxGap: 10).ToList();
        Assert.That(blocks, Is.Empty, "each chromosome-pair run of 2 < minGenes=3 must be discarded");
    }

    /// <summary>Edge (§6.1): empty ortholog input ⇒ no blocks.</summary>
    [Test]
    [Category("Property")]
    public void FindSyntenyBlocks_EmptyInput_NoBlocks()
    {
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(
            Enumerable.Empty<(string, int, int, string, string, int, int, string)>());
        Assert.That(blocks, Is.Empty);
    }

    // ===================== maxGap megabase splitting =====================

    /// <summary>
    /// maxGap megabase semantics (§3.1 / §5.2 / §6.1): inserting a genome-1 gap GREATER than
    /// <c>maxGap·1e6</c> between two otherwise-collinear forward sub-runs SPLITS the run into two
    /// blocks. Built with <c>maxGap=1</c> (⇒ 1,000,000 bp ceiling): two 3-gene sub-runs are placed so
    /// the within-run gaps are tiny but the gap BETWEEN them exceeds 1 Mb. Confirms the value is read
    /// in megabases, not base pairs.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindSyntenyBlocks_GapOverMaxGap_SplitsRun()
    {
        const int maxGap = 1; // 1 Mb ceiling
        // First sub-run near origin; second sub-run shifted by > 1 Mb in BOTH genomes.
        int shift = 3 * OneMegabase; // 3 Mb >> 1 Mb ceiling
        var pairs = new List<(string, int, int, string, string, int, int, string)>
        {
            ("chr1", 1000, 2000, "g1", "chrA", 1000, 2000, "gA"),
            ("chr1", 3000, 4000, "g2", "chrA", 3000, 4000, "gB"),
            ("chr1", 5000, 6000, "g3", "chrA", 5000, 6000, "gC"),
            ("chr1", 1000 + shift, 2000 + shift, "g4", "chrA", 1000 + shift, 2000 + shift, "gD"),
            ("chr1", 3000 + shift, 4000 + shift, "g5", "chrA", 3000 + shift, 4000 + shift, "gE"),
            ("chr1", 5000 + shift, 6000 + shift, "g6", "chrA", 5000 + shift, 6000 + shift, "gF"),
        };
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(pairs, minGenes: 3, maxGap: maxGap).ToList();
        Assert.That(blocks, Has.Count.EqualTo(2),
            "a genome-1 gap > maxGap·1e6 must split the run into two blocks");
        Assert.Multiple(() =>
        {
            Assert.That(blocks.All(b => b.GeneCount == 3), Is.True, "each split block keeps 3 genes");
            Assert.That(blocks.All(b => b.Strand == '+'), Is.True);
        });
    }

    /// <summary>
    /// Counterpart to the split test: the SAME shifted layout, but with a generous <c>maxGap=10</c>
    /// (10 Mb &gt; 3 Mb inter-run gap), keeps everything in ONE block — confirming the split above is
    /// caused by the threshold, not the layout.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindSyntenyBlocks_GapWithinMaxGap_SingleBlock()
    {
        int shift = 3 * OneMegabase;
        var pairs = new List<(string, int, int, string, string, int, int, string)>
        {
            ("chr1", 1000, 2000, "g1", "chrA", 1000, 2000, "gA"),
            ("chr1", 3000, 4000, "g2", "chrA", 3000, 4000, "gB"),
            ("chr1", 5000, 6000, "g3", "chrA", 5000, 6000, "gC"),
            ("chr1", 1000 + shift, 2000 + shift, "g4", "chrA", 1000 + shift, 2000 + shift, "gD"),
            ("chr1", 3000 + shift, 4000 + shift, "g5", "chrA", 3000 + shift, 4000 + shift, "gE"),
            ("chr1", 5000 + shift, 6000 + shift, "g6", "chrA", 5000 + shift, 6000 + shift, "gF"),
        };
        var blocks = ChromosomeAnalyzer.FindSyntenyBlocks(pairs, minGenes: 3, maxGap: 10).ToList();
        Assert.That(blocks, Has.Count.EqualTo(1), "a 3 Mb gap is within the 10 Mb ceiling ⇒ one block");
        Assert.That(blocks[0].GeneCount, Is.EqualTo(6));
    }

    // ===================== DetectRearrangements — INV-04 type set & anchors =====================

    /// <summary>
    /// Builds a forward synteny block from explicit coordinates (SequenceIdentity = NaN, as the
    /// implementation emits) for assembling rearrangement anchors independently of block detection.
    /// </summary>
    private static ChromosomeAnalyzer.SyntenyBlock MakeBlock(
        string chr1, int s1, int e1, string chr2, int s2, int e2, char strand, int genes) =>
        new(chr1, s1, e1, chr2, s2, e2, strand, genes, double.NaN);

    /// <summary>
    /// INV-04 over ARBITRARY block lists: every emitted rearrangement <c>Type</c> is one of the four
    /// valid labels. Blocks are generated with mixed chromosomes, strands and overlapping intervals so
    /// all four detection branches can fire.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DetectRearrangements_TypesAlwaysValid()
    {
        var blockGen =
            from chr1Idx in Gen.Choose(0, 1)
            from chr2Idx in Gen.Choose(0, 1)
            from s1 in Gen.Choose(1_000, 50_000)
            from len1 in Gen.Choose(500, 20_000)
            from s2 in Gen.Choose(1_000, 50_000)
            from len2 in Gen.Choose(500, 20_000)
            from strandFwd in Gen.Elements(true, false)
            from genes in Gen.Choose(3, 8)
            select MakeBlock($"chr{chr1Idx}", s1, s1 + len1, $"chr{chr2Idx}A",
                s2, s2 + len2, strandFwd ? '+' : '-', genes);

        var arb = blockGen.ListOf().Select(l => l.ToList()).ToArbitrary();
        return Prop.ForAll(arb, blocks =>
        {
            var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();
            foreach (var e in events)
                if (!ValidRearrangementTypes.Contains(e.Type))
                    return false.Label($"invalid rearrangement Type '{e.Type}'");
            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-04 anchor — Inversion: two adjacent blocks on the SAME (Chr1, Chr2) with OPPOSITE strands
    /// (§2.1 signature: strand change within the same chromosome pair) ⇒ an "Inversion" event.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DetectRearrangements_OppositeStrandsSamePair_Inversion()
    {
        var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
        {
            MakeBlock("chr1", 1_000, 5_000, "chrA", 1_000, 5_000, '+', 3),
            MakeBlock("chr1", 10_000, 15_000, "chrA", 10_000, 15_000, '-', 3),
        };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();
        Assert.That(events.Any(e => e.Type == "Inversion"), Is.True,
            "opposite strands on the same chromosome pair must report an Inversion");
        Assert.That(events.All(e => ValidRearrangementTypes.Contains(e.Type)), Is.True);
    }

    /// <summary>
    /// INV-04 anchor — Translocation: two adjacent blocks on the SAME Chr1 but DIFFERENT Chr2 (§2.1
    /// signature: change in target chromosome) ⇒ a "Translocation" event.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DetectRearrangements_DifferentTargetChromosome_Translocation()
    {
        var blocks = new List<ChromosomeAnalyzer.SyntenyBlock>
        {
            MakeBlock("chr1", 1_000, 5_000, "chrA", 1_000, 5_000, '+', 3),
            MakeBlock("chr1", 10_000, 15_000, "chrB", 10_000, 15_000, '+', 3),
        };
        var events = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();
        Assert.That(events.Any(e => e.Type == "Translocation"), Is.True,
            "different target chromosomes on adjacent same-source blocks must report a Translocation");
        Assert.That(events.All(e => ValidRearrangementTypes.Contains(e.Type)), Is.True);
    }

    /// <summary>Edge (§6.1): empty block input ⇒ no rearrangements.</summary>
    [Test]
    [Category("Property")]
    public void DetectRearrangements_EmptyInput_NoEvents()
    {
        var events = ChromosomeAnalyzer.DetectRearrangements(
            Enumerable.Empty<ChromosomeAnalyzer.SyntenyBlock>());
        Assert.That(events, Is.Empty);
    }

    /// <summary>Edge (§6.1): a single block ⇒ no rearrangements (heuristics need ≥ 2 blocks).</summary>
    [FsCheck.NUnit.Property]
    public Property DetectRearrangements_SingleBlock_NoEvents()
    {
        var blockGen =
            from s1 in Gen.Choose(1_000, 50_000)
            from len1 in Gen.Choose(500, 5_000)
            from s2 in Gen.Choose(1_000, 50_000)
            from len2 in Gen.Choose(500, 5_000)
            from strandFwd in Gen.Elements(true, false)
            select MakeBlock("chr1", s1, s1 + len1, "chrA", s2, s2 + len2, strandFwd ? '+' : '-', 3);
        return Prop.ForAll(blockGen.ToArbitrary(), block =>
        {
            var events = ChromosomeAnalyzer.DetectRearrangements(new[] { block }).ToList();
            return (events.Count == 0).Label($"single block produced {events.Count} events");
        });
    }

    // ===================== D — determinism =====================

    /// <summary>
    /// D (determinism): identical ortholog input ⇒ identical <c>FindSyntenyBlocks</c> output sequences
    /// (all nine fields in order, NaN==NaN for SequenceIdentity).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindSyntenyBlocks_IsDeterministic()
    {
        return Prop.ForAll(ForwardCollinearRunArbitrary(), genes =>
        {
            var a = ChromosomeAnalyzer.FindSyntenyBlocks(ToTuples(genes), minGenes: 3, maxGap: 10).ToList();
            var b = ChromosomeAnalyzer.FindSyntenyBlocks(ToTuples(genes), minGenes: 3, maxGap: 10).ToList();
            bool same = a.Count == b.Count && a.Zip(b, BlocksEqual).All(x => x);
            return same.Label($"non-deterministic FindSyntenyBlocks (n={genes.Count})");
        });
    }

    /// <summary>
    /// D (determinism): identical block input ⇒ identical <c>DetectRearrangements</c> output sequences
    /// (record-struct value equality, field-by-field in order).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DetectRearrangements_IsDeterministic()
    {
        var blockGen =
            from chr1Idx in Gen.Choose(0, 1)
            from chr2Idx in Gen.Choose(0, 1)
            from s1 in Gen.Choose(1_000, 50_000)
            from len1 in Gen.Choose(500, 20_000)
            from s2 in Gen.Choose(1_000, 50_000)
            from len2 in Gen.Choose(500, 20_000)
            from strandFwd in Gen.Elements(true, false)
            select MakeBlock($"chr{chr1Idx}", s1, s1 + len1, $"chr{chr2Idx}A",
                s2, s2 + len2, strandFwd ? '+' : '-', 3);
        var arb = blockGen.ListOf().Select(l => l.ToList()).ToArbitrary();

        return Prop.ForAll(arb, blocks =>
        {
            var a = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();
            var b = ChromosomeAnalyzer.DetectRearrangements(blocks).ToList();
            return a.SequenceEqual(b).Label($"non-deterministic DetectRearrangements (n={blocks.Count})");
        });
    }

    #endregion

    #region CHROM-ALPHASAT-001: R: monomer period ≈ 171 bp; R: CENP-B boxes within monomers; D: deterministic

    // DetectAlphaSatellite / FindCenpBBoxes — alpha-satellite monomer detection (171-bp period, Willard 1985;
    // Waye & Willard 1987) and the CENP-B box consensus YTTCGTTGGAARCGGGA (Masumoto et al. 1989).

    private const int AlphaMonomer = ChromosomeAnalyzer.AlphaSatelliteMonomerLength; // 171
    private const string CenpBInstance = "CTTCGTTGGAAACGGGA"; // a valid CENP-B box (Y=C, R=A), 17 bp

    /// <summary>Builds an AT-rich 171-bp monomer carrying a CENP-B box at its 5' end, repeated m times.</summary>
    private static string BuildAlphaSatellite(int seed, int m)
    {
        var rng = new Random(seed);
        // 70%-AT filler so the overall AT content clears the >0.50 alpha-satellite gate.
        char[] atRich = { 'A', 'A', 'T', 'T', 'A', 'T', 'G', 'C', 'A', 'T' };
        var filler = new char[AlphaMonomer - CenpBInstance.Length];
        for (int i = 0; i < filler.Length; i++) filler[i] = atRich[rng.Next(atRich.Length)];
        string monomer = CenpBInstance + new string(filler);
        return string.Concat(Enumerable.Repeat(monomer, m));
    }

    private static Arbitrary<string> AlphaSatelliteArbitrary() =>
        (from seed in Gen.Choose(1, 1_000_000) from m in Gen.Choose(3, 5) select BuildAlphaSatellite(seed, m)).ToArbitrary();

    /// <summary>
    /// INV-1 (R): when a sequence is called alpha-satellite, its detected monomer period is ≈ 171 bp
    /// (within the ±5 bp tolerance window). Exercised on AT-rich 171-bp monomer tandems.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaSatellite_Period_NearMonomerLength()
    {
        return Prop.ForAll(AlphaSatelliteArbitrary(), seq =>
        {
            var r = ChromosomeAnalyzer.DetectAlphaSatellite(seq);
            if (!r.IsAlphaSatellite) return true.ToProperty();
            return (r.BestPeriod >= AlphaMonomer - 5 && r.BestPeriod <= AlphaMonomer + 5)
                .Label($"alpha-satellite period {r.BestPeriod} not ≈ {AlphaMonomer}");
        });
    }

    /// <summary>
    /// INV-2 (R): every CENP-B box found lies within the monomeric sequence — a valid 17-bp window
    /// (0 ≤ pos ≤ len − 17) — and the reported count matches the number of box positions.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AlphaSatellite_CenpBBoxes_WithinMonomers()
    {
        return Prop.ForAll(AlphaSatelliteArbitrary(), seq =>
        {
            var boxes = ChromosomeAnalyzer.FindCenpBBoxes(seq);
            var r = ChromosomeAnalyzer.DetectAlphaSatellite(seq);
            bool ok = boxes.All(p => p >= 0 && p + CenpBInstance.Length <= seq.Length)
                      && r.CenpBBoxCount == boxes.Count;
            return ok.Label($"CENP-B boxes invalid (count={boxes.Count}, reported={r.CenpBBoxCount})");
        });
    }

    /// <summary>INV-3 (D): alpha-satellite detection is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property AlphaSatellite_IsDeterministic()
    {
        return Prop.ForAll(AlphaSatelliteArbitrary(), seq =>
            (ChromosomeAnalyzer.DetectAlphaSatellite(seq) == ChromosomeAnalyzer.DetectAlphaSatellite(seq))
                .Label("DetectAlphaSatellite must be deterministic"));
    }

    #endregion

    #region CHROM-HOR-001: R: inter-HOR identity ≥ intra-monomer identity; R: HOR period = k×monomer; D: deterministic

    // DetectHigherOrderRepeat — alpha-satellite higher-order repeat (HOR) structure: a unit of k distinct
    // monomers tandemly repeated, so copies of the same unit are more identical (inter-HOR) than the distinct
    // monomers within a unit (intra-HOR). Source: Willard & Waye (1987); Miga et al. (2014, T2T).

    /// <summary>Builds a HOR: k distinct random 171-bp monomers forming a unit, repeated c times exactly.</summary>
    private static string BuildHor(int seed, int k, int c)
    {
        var rng = new Random(seed);
        char[] bases = { 'A', 'C', 'G', 'T' };
        var monomers = new string[k];
        for (int j = 0; j < k; j++)
        {
            var m = new char[AlphaMonomer];
            for (int i = 0; i < AlphaMonomer; i++) m[i] = bases[rng.Next(4)];
            monomers[j] = new string(m);
        }
        string unit = string.Concat(monomers);
        return string.Concat(Enumerable.Repeat(unit, c));
    }

    private static Arbitrary<string> HorArbitrary() =>
        (from seed in Gen.Choose(1, 1_000_000) from k in Gen.Choose(2, 3) from c in Gen.Choose(2, 3)
         select BuildHor(seed, k, c)).ToArbitrary();

    /// <summary>
    /// INV-1 (R): when a higher-order structure is detected, the HOR unit length equals MonomersPerUnit ×
    /// monomer length (the period is an integer number of monomers).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hor_UnitLength_IsMultipleOfMonomer()
    {
        return Prop.ForAll(HorArbitrary(), seq =>
        {
            var r = ChromosomeAnalyzer.DetectHigherOrderRepeat(seq, AlphaMonomer);
            if (!r.HasHigherOrderStructure) return true.ToProperty();
            return (r.HorUnitLengthBp == r.MonomersPerUnit * AlphaMonomer)
                .Label($"HOR unit {r.HorUnitLengthBp} ≠ {r.MonomersPerUnit}×{AlphaMonomer}");
        });
    }

    /// <summary>
    /// INV-2 (R): in a detected HOR, copies of the same unit are at least as identical as the distinct
    /// monomers within a unit — mean inter-HOR identity ≥ mean intra-HOR identity (the defining HOR property).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hor_InterIdentity_AtLeastIntra()
    {
        return Prop.ForAll(HorArbitrary(), seq =>
        {
            var r = ChromosomeAnalyzer.DetectHigherOrderRepeat(seq, AlphaMonomer);
            if (!r.HasHigherOrderStructure) return true.ToProperty();
            return (r.MeanInterHorIdentity >= r.MeanIntraHorIdentity - 1e-9)
                .Label($"inter-HOR {r.MeanInterHorIdentity} < intra-HOR {r.MeanIntraHorIdentity}");
        });
    }

    /// <summary>INV-3 (D): HOR detection is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property Hor_IsDeterministic()
    {
        return Prop.ForAll(HorArbitrary(), seq =>
            (ChromosomeAnalyzer.DetectHigherOrderRepeat(seq, AlphaMonomer) == ChromosomeAnalyzer.DetectHigherOrderRepeat(seq, AlphaMonomer))
                .Label("DetectHigherOrderRepeat must be deterministic"));
    }

    #endregion
}
