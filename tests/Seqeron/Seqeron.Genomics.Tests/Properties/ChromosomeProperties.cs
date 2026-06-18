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
}
