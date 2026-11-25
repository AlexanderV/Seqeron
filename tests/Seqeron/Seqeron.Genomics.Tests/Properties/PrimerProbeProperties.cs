using System.Text;
using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for primer/probe design: melting temperature, hairpin, dimer, probe validation.
///
/// Test Units: PRIMER-TM-001, PRIMER-DESIGN-001, PRIMER-STRUCT-001, PROBE-DESIGN-001, PROBE-VALID-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("MolTools")]
public class PrimerProbeProperties
{

    /// <summary>
    /// Generates random valid DNA primers (only A/C/G/T) within an explicit length window.
    /// </summary>
    private static Gen<string> ValidPrimerGen(int minLen, int maxLen) =>
        from len in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
        select new string(chars);

    /// <summary>Wraps <see cref="ValidPrimerGen"/> as an arbitrary.</summary>
    private static Arbitrary<string> ValidPrimerArbitrary(int minLen, int maxLen) =>
        ValidPrimerGen(minLen, maxLen).ToArbitrary();

    /// <summary>
    /// Generates arbitrary character soup: valid DNA bases mixed with case variants,
    /// IUPAC ambiguity codes, RNA uracil, and pure junk — including the empty string.
    /// Used to verify defined behaviors (Section 3 of the TestSpec).
    /// </summary>
    private static Arbitrary<string> MixedSequenceArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T', 'a', 'c', 'g', 't',
                     'N', 'n', 'U', 'u', 'R', 'Y', '-', ' ', 'x', '7')
            .ArrayOf()
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// Generates a length N and two G/C counts (lo ≤ hi) within [0, N], used to build
    /// equal-length sequences that differ only in GC composition.
    /// </summary>
    private static Arbitrary<(int n, int gcLo, int gcHi)> EqualLengthGcPairArbitrary() =>
        (from n in Gen.Choose(1, 60)
         from a in Gen.Choose(0, n)
         from b in Gen.Choose(0, n)
         select (n, Math.Min(a, b), Math.Max(a, b))).ToArbitrary();

    /// <summary>
    /// Generates a valid primer paired with two Na⁺ concentrations (mM) for salt-correction tests.
    /// </summary>
    private static Arbitrary<(string primer, int na1, int na2)> PrimerSaltPairArbitrary() =>
        (from p in ValidPrimerGen(8, 30)
         from na1 in Gen.Choose(1, 2000)
         from na2 in Gen.Choose(1, 2000)
         select (p, na1, na2)).ToArbitrary();

    /// <summary>
    /// Generates a clean ACGT primer together with a "dirty" copy that has non-ACGT
    /// characters sprinkled around every base. Tm must be identical (Section 3.1).
    /// </summary>
    private static Arbitrary<(string clean, string dirty)> CleanDirtyArbitrary() =>
        (from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= 1)
         // U/u is NOT junk: the basic Tm reads it as T (Biopython _check back-transcription), see
         // MeltingTemperature_UracilCountsAsThymine.
         from junk in Gen.Elements('N', 'n', 'R', 'Y', '-', ' ', 'x', '9', '.')
                         .ArrayOf().Where(a => a.Length >= 1)
         let clean = new string(chars)
         select (clean, Sprinkle(clean, new string(junk)))).ToArbitrary();

    /// <summary>Builds a string with a junk character inserted before and after every clean base.</summary>
    private static string Sprinkle(string clean, string junk)
    {
        var sb = new StringBuilder();
        sb.Append(junk[0]);
        for (int i = 0; i < clean.Length; i++)
        {
            sb.Append(clean[i]);
            sb.Append(junk[i % junk.Length]);
        }
        return sb.ToString();
    }

    /// <summary>
    /// Counts the bases the basic Tm uses (case-insensitive), partitioned into A/T and G/C; U counts
    /// as T (Biopython <c>MeltingTemp._check</c> back-transcription, ThermoConstants.CountBasicTmBases).
    /// </summary>
    private static (int at, int gc) CountAtGc(string s)
    {
        int at = 0, gc = 0;
        foreach (char ch in s.ToUpperInvariant())
        {
            if (ch is 'A' or 'T' or 'U') at++;
            else if (ch is 'G' or 'C') gc++;
        }
        return (at, gc);
    }

    /// <summary>
    /// Independent theory oracle for Tm, derived directly from the published equations
    /// (Thein &amp; Wallace 1986; Marmur &amp; Doty 1962) with literal constants — deliberately
    /// NOT routed through <c>ThermoConstants</c>, so a wrong production constant is caught.
    /// </summary>
    private static double ExpectedTm(string s)
    {
        var (at, gc) = CountAtGc(s);
        int n = at + gc;
        if (n == 0) return 0;
        // Wallace rule for short oligos (&lt; 14 valid bases): Tm = 2(A+T) + 4(G+C).
        if (n < 14) return 2 * at + 4 * gc;
        // Marmur-Doty for longer primers: Tm = 64.9 + 41(GC − 16.4)/N, floored at 0.
        return Math.Max(0, 64.9 + 41.0 * (gc - 16.4) / n);
    }

    /// <summary>Upper bound on Tm for any input: the Marmur-Doty asymptote 64.9 + 41 = 105.9°C.</summary>
    private static double TmUpperBound =>
        ThermoConstants.MarmurDotyBase + ThermoConstants.MarmurDotyGcCoefficient;



    /// <summary>
    /// INV-Wallace: For any primer with &lt; 14 valid bases, Tm equals the Wallace-rule value
    /// 2(A+T) + 4(G+C) computed independently from base counts.
    /// Evidence: Thein &amp; Wallace (1986). Validates regime selection, base counting, and
    /// the A/T = 2°C, G/C = 4°C contributions across hundreds of random short oligos.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_ShortPrimer_MatchesWallaceRule()
    {
        return Prop.ForAll(ValidPrimerArbitrary(1, 13), primer =>
        {
            double actual = PrimerDesigner.CalculateMeltingTemperature(primer);
            return (Math.Abs(actual - ExpectedTm(primer)) < 1e-9)
                .Label($"Wallace Tm mismatch for '{primer}': {actual} vs {ExpectedTm(primer)}");
        });
    }

    /// <summary>
    /// INV-MarmurDoty: For any primer with ≥ 14 valid bases, Tm equals
    /// max(0, 64.9 + 41(GC − 16.4)/N) computed independently from base counts.
    /// Evidence: Marmur &amp; Doty (1962). Validates the long-primer regime and GC/length scaling.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_LongPrimer_MatchesMarmurDoty()
    {
        return Prop.ForAll(ValidPrimerArbitrary(14, 40), primer =>
        {
            double actual = PrimerDesigner.CalculateMeltingTemperature(primer);
            return (Math.Abs(actual - ExpectedTm(primer)) < 1e-9)
                .Label($"Marmur-Doty Tm mismatch for '{primer}': {actual} vs {ExpectedTm(primer)}");
        });
    }

    /// <summary>
    /// INV-Threshold: The Wallace ↔ Marmur-Doty switch happens exactly at 14 valid bases.
    /// A 13-mer uses Wallace (integer multiple of 2/4); a 14-mer uses Marmur-Doty.
    /// Evidence: TestSpec Section 2 (threshold &lt; 14). Pins the boundary that the
    /// random-length properties exercise only stochastically.
    /// </summary>
    [Test]
    [Category("Property")]
    public void MeltingTemperature_RegimeThreshold_SwitchesAt14()
    {
        // 13 valid bases → Wallace: 7×A/T (2) + 6×G/C (4) = 14 + 24 = 38.
        double tm13 = PrimerDesigner.CalculateMeltingTemperature("ACGTACGTACGTA");
        // 14 valid bases → Marmur-Doty: 64.9 + 41×(7−16.4)/14 ≈ 37.36 (non-integer ⇒ not Wallace).
        double tm14 = PrimerDesigner.CalculateMeltingTemperature("ACGTACGTACGTAC");

        Assert.Multiple(() =>
        {
            Assert.That(tm13, Is.EqualTo(38.0).Within(1e-9), "13-mer must use the Wallace rule");
            Assert.That(tm14, Is.EqualTo(64.9 + 41.0 * (7 - 16.4) / 14).Within(1e-9),
                "14-mer must use the Marmur-Doty formula");
        });
    }



    /// <summary>
    /// Anchors Tm to absolute, literature-published values so that a self-consistent but
    /// wrong constant (shared between production and oracle) cannot pass unnoticed.
    /// Sources: Thein &amp; Wallace (1986) for short oligos; Marmur &amp; Doty (1962) for ≥14-mers.
    /// </summary>
    [TestCase("ATATATAT", 16.0, TestName = "Wallace: 8×A/T → 16°C")]
    [TestCase("GCGCGCGC", 32.0, TestName = "Wallace: 8×G/C → 32°C")]
    [TestCase("ACGTACGT", 24.0, TestName = "Wallace: 4 A/T + 4 G/C → 24°C")]
    [TestCase("ATATATATATATATATATAT", 31.28, TestName = "Marmur-Doty: 20bp 0% GC → 31.28°C")]
    [TestCase("ACGTACGTACGTACGTACGT", 51.78, TestName = "Marmur-Doty: 20bp 50% GC → 51.78°C")]
    [TestCase("GCGCGCGCGCGCGCGCGCGC", 72.28, TestName = "Marmur-Doty: 20bp 100% GC → 72.28°C")]
    [Category("Property")]
    public void MeltingTemperature_CanonicalLiteratureValues(string primer, double expected)
    {
        double actual = PrimerDesigner.CalculateMeltingTemperature(primer);
        Assert.That(actual, Is.EqualTo(expected).Within(0.01),
            $"Tm for '{primer}' must match the published value {expected}°C");
    }



    /// <summary>
    /// INV-1 (R): Tm is always within the physically meaningful band [0, 105.9) for ANY input,
    /// including junk, mixed case, and the empty string. The upper bound is the Marmur-Doty
    /// asymptote (64.9 + 41), which an all-G/C primer approaches but never reaches.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_AlwaysInPhysicalRange()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            double tm = PrimerDesigner.CalculateMeltingTemperature(seq);
            return (tm >= 0.0 && tm < TmUpperBound)
                .Label($"Tm {tm} outside [0, {TmUpperBound}) for '{seq}'");
        });
    }

    /// <summary>
    /// INV-2 (P): Tm is strictly positive iff the sequence contains at least one valid DNA base,
    /// and exactly 0 otherwise. This separates a correct DNA-aware implementation from one that
    /// would emit a spurious temperature for sequences like "NNNN" or "" (Section 3.1/3.3).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_PositiveIffValidBasePresent()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            var (at, gc) = CountAtGc(seq);
            bool hasValidBase = at + gc > 0;
            double tm = PrimerDesigner.CalculateMeltingTemperature(seq);
            return ((tm > 0.0) == hasValidBase)
                .Label($"Tm {tm} vs validBases={at + gc} for '{seq}'");
        });
    }



    /// <summary>
    /// INV-3 (M): For two sequences of equal valid length, the one with more G/C bases has a
    /// strictly higher Tm; equal G/C content yields an equal Tm. This is the central biological
    /// invariant — G·C pairs (3 hydrogen bonds) stabilize the duplex more than A·T pairs (2 H-bonds),
    /// so higher GC ⇒ higher melting temperature. Equal length keeps both sequences in the same
    /// formula regime, making the comparison rigorous rather than confounded by the 14-base threshold.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_HigherGc_StrictlyHigherTm_AtEqualLength()
    {
        return Prop.ForAll(EqualLengthGcPairArbitrary(), triple =>
        {
            var (n, gcLo, gcHi) = triple;
            string seqLo = new string('G', gcLo) + new string('A', n - gcLo);
            string seqHi = new string('G', gcHi) + new string('A', n - gcHi);
            double tmLo = PrimerDesigner.CalculateMeltingTemperature(seqLo);
            double tmHi = PrimerDesigner.CalculateMeltingTemperature(seqHi);

            bool ok = gcHi > gcLo ? tmHi > tmLo : Math.Abs(tmHi - tmLo) < 1e-9;
            return ok.Label($"n={n}, gc {gcLo}→{tmLo} vs {gcHi}→{tmHi}");
        });
    }

    /// <summary>
    /// INV-4 (M, length): a GC-rich elongation — appending a G or C base — strictly raises Tm.
    /// This is the "longer GC-rich → higher Tm" claim. In the Wallace regime each G/C adds +4 °C;
    /// in the Marmur-Doty regime Tm = 64.9 + 41·(gc−16.4)/N the increment from a G/C base is
    /// (at + 16.4)/(N·(N+1)) > 0 for all primers. (Restricted to the long-primer regime, length
    /// ≥ 14 before and after, so the comparison is not confounded by the Wallace↔Marmur-Doty switch
    /// at length 14. Note: appending an A/T base is NOT monotone — it lowers Marmur-Doty Tm for
    /// GC-rich primers — so only the GC-rich elongation is asserted.)
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_GcRichElongation_StrictlyHigherTm()
    {
        var gen = (from primer in ValidPrimerGen(14, 39)
                   from added in Gen.Elements('G', 'C')
                   select (primer, added)).ToArbitrary();

        return Prop.ForAll(gen, t =>
        {
            double before = PrimerDesigner.CalculateMeltingTemperature(t.primer);
            double after = PrimerDesigner.CalculateMeltingTemperature(t.primer + t.added);
            return (after > before)
                .Label($"Tm(len {t.primer.Length})={before:F3} → append '{t.added}' → {after:F3} must strictly increase");
        });
    }



    /// <summary>
    /// INV-4 (P): Tm is case-insensitive — uppercase, lowercase, and mixed-case spellings of the
    /// same primer all yield the identical temperature (Section 3.2).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_IsCaseInsensitive()
    {
        return Prop.ForAll(ValidPrimerArbitrary(1, 40), primer =>
        {
            char[] m = primer.ToCharArray();
            for (int i = 0; i < m.Length; i += 2) m[i] = char.ToLowerInvariant(m[i]);
            string mixed = new string(m);

            double up = PrimerDesigner.CalculateMeltingTemperature(primer.ToUpperInvariant());
            double lo = PrimerDesigner.CalculateMeltingTemperature(primer.ToLowerInvariant());
            double mi = PrimerDesigner.CalculateMeltingTemperature(mixed);

            return (Math.Abs(up - lo) < 1e-9 && Math.Abs(up - mi) < 1e-9)
                .Label($"Case sensitivity: up={up}, lo={lo}, mixed={mi} for '{primer}'");
        });
    }

    /// <summary>
    /// INV-5 (P): Non-ACGT characters are ignored — sprinkling N's, U's, gaps, and other junk
    /// around a primer leaves its Tm unchanged, because only valid DNA bases are counted
    /// (Section 3.1). This protects against an implementation that lets ambiguity codes leak
    /// into the temperature.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_IgnoresNonAcgtCharacters()
    {
        return Prop.ForAll(CleanDirtyArbitrary(), pair =>
        {
            var (clean, dirty) = pair;
            double tmClean = PrimerDesigner.CalculateMeltingTemperature(clean);
            double tmDirty = PrimerDesigner.CalculateMeltingTemperature(dirty);
            return (Math.Abs(tmClean - tmDirty) < 1e-9)
                .Label($"junk changed Tm: clean '{clean}'={tmClean}, dirty '{dirty}'={tmDirty}");
        });
    }

    /// <summary>
    /// INV-5b: U is read as T — the RNA spelling of a primer has the same basic Tm as its DNA spelling
    /// (ThermoConstants.CountBasicTmBases: Biopython <c>MeltingTemp._check</c> back-transcription).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_UracilCountsAsThymine()
    {
        return Prop.ForAll(CleanDirtyArbitrary(), pair =>
        {
            string dna = pair.clean;
            string rna = dna.Replace('T', 'U');
            double tmDna = PrimerDesigner.CalculateMeltingTemperature(dna);
            double tmRna = PrimerDesigner.CalculateMeltingTemperature(rna);
            return (Math.Abs(tmDna - tmRna) < 1e-9)
                .Label($"U not read as T: '{dna}'={tmDna}, '{rna}'={tmRna}");
        });
    }

    /// <summary>
    /// INV-6 (D): Tm is deterministic — repeated evaluation of the same input is bit-identical.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperature_IsDeterministic()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            double a = PrimerDesigner.CalculateMeltingTemperature(seq);
            double b = PrimerDesigner.CalculateMeltingTemperature(seq);
            return a.Equals(b).Label($"non-deterministic Tm for '{seq}': {a} vs {b}");
        });
    }



    /// <summary>
    /// INV-7 (P): Salt-adjusted Tm equals the OligoCalc "Salt Adjusted" formula (Kibbe 2007),
    /// written here independently of production code: N &lt; 14 → 2(A+T)+4(G+C)+16.6·log10([Na+]/0.050 M);
    /// N ≥ 14 → 100.5 + 41·GC/N − 820/N + 16.6·log10([Na+] M); rounded to one decimal.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperatureWithSalt_MatchesOligoCalcSaltAdjusted()
    {
        return Prop.ForAll(PrimerSaltPairArbitrary(), t =>
        {
            var (primer, na, _) = t;
            string u = primer.ToUpperInvariant();
            int at = u.Count(c => c is 'A' or 'T'), gc = u.Count(c => c is 'G' or 'C'), n = at + gc;
            double m = na / 1000.0;
            double oracle = n < 14
                ? 2 * at + 4 * gc + 16.6 * Math.Log10(m / 0.050)
                : 100.5 + 41.0 * gc / n - 820.0 / n + 16.6 * Math.Log10(m);
            double expected = Math.Round(oracle, 1);
            double actual = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, na);
            return (Math.Abs(actual - expected) < 1e-9)
                .Label($"salt Tm mismatch for '{primer}' @ {na}mM: {actual} vs {expected}");
        });
    }

    /// <summary>
    /// INV-8 (M): Higher salt stabilizes the duplex — increasing [Na⁺] never decreases the
    /// corrected Tm. Evidence: Owczarzy et al. (2004); higher cation concentration shields the
    /// phosphate backbone and raises melting temperature.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperatureWithSalt_IsMonotonicInSalt()
    {
        return Prop.ForAll(PrimerSaltPairArbitrary(), t =>
        {
            var (primer, na1, na2) = t;
            int lo = Math.Min(na1, na2), hi = Math.Max(na1, na2);
            double tmLo = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, lo);
            double tmHi = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, hi);
            return (tmHi >= tmLo)
                .Label($"salt monotonicity broken for '{primer}': {lo}mM→{tmLo}, {hi}mM→{tmHi}");
        });
    }

    /// <summary>
    /// INV-9 (P): The basic Wallace rule is defined at 50 mM Na+ (OligoCalc), so for short oligos
    /// (&lt; 14 valid bases) the salt-adjusted Tm at 50 mM equals the basic Tm — the salt term is not
    /// counted twice.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperatureWithSalt_ShortOligoAt50mM_EqualsBaseTm()
    {
        return Prop.ForAll(ValidPrimerArbitrary(1, 13), primer =>
        {
            double expected = Math.Round(PrimerDesigner.CalculateMeltingTemperature(primer), 1);
            double actual = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, 50);
            return (Math.Abs(actual - expected) < 1e-9)
                .Label($"@50mM '{primer}': {actual} vs base {expected}");
        });
    }

    /// <summary>
    /// INV-10 (M): Raising [Na+] from 50 mM to 1 M raises the salt-adjusted Tm by
    /// 16.6·log10(20) ≈ 21.6 °C for every primer (Schildkraut–Lifson slope), within rounding.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property MeltingTemperatureWithSalt_50mMTo1M_Adds16_6LogTwenty()
    {
        return Prop.ForAll(ValidPrimerArbitrary(8, 30), primer =>
        {
            double d = PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, 1000)
                       - PrimerDesigner.CalculateMeltingTemperatureWithSalt(primer, 50);
            return (Math.Abs(d - 16.6 * Math.Log10(20.0)) <= 0.1 + 1e-9)
                .Label($"'{primer}': ΔTm(50mM→1M) = {d}");
        });
    }

    /// <summary>
    /// INV-11 (Edge): Salt correction is undefined for an empty/null primer — both return 0,
    /// regardless of the requested concentration (Section 3.3).
    /// </summary>
    [Test]
    [Category("Property")]
    public void MeltingTemperatureWithSalt_EmptyOrNull_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateMeltingTemperatureWithSalt("", 50), Is.EqualTo(0.0));
            Assert.That(PrimerDesigner.CalculateMeltingTemperatureWithSalt(null!, 50), Is.EqualTo(0.0));
        });
    }



    /// <summary>
    /// Independent SantaLucia (1998) nearest-neighbor ΔG°37 oracle for the 3' terminal 5-mer.
    /// Constants are transcribed LITERALLY from Primer_Structure_Analysis.md §2.2/§4.2 and the
    /// SantaLucia (1998) unified table — deliberately NOT routed through production, so a wrong
    /// production constant is caught. ΔG = Σ(four NN dinucleotide steps) + two terminal initiation
    /// terms (+0.98 for a terminal G·C, +1.03 for a terminal A·T). Mirrors the doc anchors
    /// GCGCG = −6.86 (most stable) and TATAT = −0.86 (least stable).
    /// </summary>
    private static double ExpectedThreePrimeStability(string sequence)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 5)
            return 0;

        // Unified NN ΔG°37 (kcal/mol), SantaLucia (1998) — written here independently.
        var nn = new Dictionary<string, double>
        {
            ["AA"] = -1.00, ["TT"] = -1.00,
            ["AT"] = -0.88,
            ["TA"] = -0.58,
            ["CA"] = -1.45, ["TG"] = -1.45,
            ["GT"] = -1.44, ["AC"] = -1.44,
            ["CT"] = -1.28, ["AG"] = -1.28,
            ["GA"] = -1.30, ["TC"] = -1.30,
            ["CG"] = -2.17,
            ["GC"] = -2.24,
            ["GG"] = -1.84, ["CC"] = -1.84,
        };

        string last5 = sequence.ToUpperInvariant()[^5..];
        double dg = 0;
        for (int i = 0; i < last5.Length - 1; i++)
            if (nn.TryGetValue(last5.Substring(i, 2), out double step))
                dg += step;

        // Terminal initiation: +0.98 kcal/mol per terminal G·C, +1.03 per terminal A·T.
        dg += last5[0] is 'G' or 'C' ? 0.98 : 1.03;
        dg += last5[^1] is 'G' or 'C' ? 0.98 : 1.03;
        return dg;
    }

    /// <summary>Generates valid ACGT sequences with length ≥ 5 (the 3'-stability operating regime).</summary>
    private static Arbitrary<string> StabilitySequenceArbitrary() =>
        ValidPrimerGen(5, 40).ToArbitrary();

    /// <summary>
    /// Generates a valid ACGT sequence together with hairpin stem/loop knobs. The sequence length
    /// is drawn small enough to frequently fall below the structural minimum 2·stem + loop, so
    /// INV-01 is exercised on both sides of the threshold.
    /// </summary>
    private static Arbitrary<(string seq, int stem, int loop)> HairpinKnobArbitrary() =>
        (from stem in Gen.Choose(2, 6)
         from loop in Gen.Choose(1, 6)
         from len in Gen.Choose(1, 40)
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
         select (new string(chars), stem, loop)).ToArbitrary();

    /// <summary>
    /// Generates a pair of valid ACGT primers together with two complementarity thresholds, for
    /// HasPrimerDimer monotonicity. Primers are ≥ 8 bases so the 3' comparison window is full.
    /// </summary>
    private static Arbitrary<(string p1, string p2, int thrLo, int thrHi)> DimerThresholdArbitrary() =>
        (from p1 in ValidPrimerGen(8, 30)
         from p2 in ValidPrimerGen(8, 30)
         from a in Gen.Choose(1, 8)
         from b in Gen.Choose(1, 8)
         select (p1, p2, Math.Min(a, b), Math.Max(a, b))).ToArbitrary();

    /// <summary>
    /// Independent transcription of the production complement relation (A↔T, G↔C). Used to
    /// recompute the terminal complementary-pair count for HasPrimerDimer anchors.
    /// </summary>
    private static bool IsWatsonCrick(char a, char b) =>
        (a, b) is ('A', 'T') or ('T', 'A') or ('G', 'C') or ('C', 'G');

    /// <summary>Reverse complement built independently (no production helper) for dimer oracle work.</summary>
    private static string RevComp(string s)
    {
        var chars = s.ToUpperInvariant().ToCharArray();
        Array.Reverse(chars);
        for (int i = 0; i < chars.Length; i++)
            chars[i] = chars[i] switch { 'A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G', var c => c };
        return new string(chars);
    }



    /// <summary>
    /// INV-NN-oracle: for ANY ACGT sequence of length ≥ 5, production
    /// <see cref="PrimerDesigner.Calculate3PrimeStability"/> equals the independent SantaLucia
    /// nearest-neighbor oracle within 1e-9. This is the rigorous core: it pins every NN step and
    /// both initiation terms to literally transcribed published constants rather than the code's own.
    /// Source: SantaLucia (1998); Primer_Structure_Analysis.md §2.2/§4.2.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Calculate3PrimeStability_MatchesSantaLuciaOracle()
    {
        return Prop.ForAll(StabilitySequenceArbitrary(), seq =>
        {
            double actual = PrimerDesigner.Calculate3PrimeStability(seq);
            double expected = ExpectedThreePrimeStability(seq);
            return (Math.Abs(actual - expected) < 1e-9)
                .Label($"ΔG mismatch for '{seq}': {actual} vs oracle {expected}");
        });
    }

    /// <summary>
    /// INV-01 (R): ΔG ≤ 0 for every sequence of length ≥ 5. Every NN step is negative and the
    /// most positive initiation pair (+1.03 twice) cannot outweigh the four steps, so the duplex
    /// stability is never positive (the anchor TATAT = −0.86 is the documented least-stable case).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Calculate3PrimeStability_IsNonPositive_ForLengthAtLeast5()
    {
        return Prop.ForAll(StabilitySequenceArbitrary(), seq =>
        {
            double dg = PrimerDesigner.Calculate3PrimeStability(seq);
            return (dg <= 0.0 && double.IsFinite(dg))
                .Label($"ΔG {dg} not ≤ 0 (finite) for '{seq}'");
        });
    }

    /// <summary>
    /// INV-02 (Edge): a primer shorter than 5 nt is its own 3′ window — Primer3
    /// <c>end_oligodg(s, 5)</c> evaluates <c>oligodg</c> on the whole primer — and empty input gives 0.
    /// Checked against an independent oracle built from SantaLucia (1998) Table 1 unified ΔG°37 with
    /// the <c>oligodg</c> terms (initiation +1.96, +0.05 per terminal A·T, +0.43 if self-complementary).
    /// Source: Primer_Structure_Analysis.md §2.4 INV-02 and §6 ("Primer shorter than 5 nt").
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Calculate3PrimeStability_UsesWholePrimer_ForLengthBelow5()
    {
        var shortGen = (from len in Gen.Choose(0, 4)
                        from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
                        select new string(chars)).ToArbitrary();
        return Prop.ForAll(shortGen, seq =>
        {
            double dg = PrimerDesigner.Calculate3PrimeStability(seq);
            double expected = seq.Length == 0 ? 0.0 : OligoDgOracle(seq);
            return (Math.Abs(dg - expected) < 1e-9)
                .Label($"ΔG {dg} != whole-primer oligodg {expected} for short '{seq}'");
        });
    }

    /// <summary>
    /// Independent Primer3 <c>oligodg</c> (SantaLucia 1998) oracle, kcal/mol: Σ Table 1 ΔG°37 over the
    /// NN steps + 1.96 initiation + 0.05 per terminal A·T end (+0.43 for an even-length
    /// self-complementary sequence).
    /// </summary>
    private static double OligoDgOracle(string seq)
    {
        static double Step(char a, char b) => (a, b) switch
        {
            ('A', 'A') or ('T', 'T') => -1.00,
            ('A', 'T') => -0.88,
            ('T', 'A') => -0.58,
            ('C', 'A') or ('T', 'G') => -1.45,
            ('G', 'T') or ('A', 'C') => -1.44,
            ('C', 'T') or ('A', 'G') => -1.28,
            ('G', 'A') or ('T', 'C') => -1.30,
            ('C', 'G') => -2.17,
            ('G', 'C') => -2.24,
            ('G', 'G') or ('C', 'C') => -1.84,
            _ => throw new ArgumentException($"not ACGT: {a}{b}"),
        };
        double dg = 1.96;
        for (int i = 0; i + 1 < seq.Length; i++) dg += Step(seq[i], seq[i + 1]);
        if (seq[0] is 'A' or 'T') dg += 0.05;
        if (seq[^1] is 'A' or 'T') dg += 0.05;
        if (seq.Length % 2 == 0 && seq == RevComp(seq)) dg += 0.43;
        return dg;
    }

    /// <summary>
    /// INV-M (terminal G·C ⇒ not less stable): replacing the terminal A/T of an A/T-flanked 5-mer
    /// with G/C makes the duplex no LESS stable (ΔG not greater). G·C contributes a stronger NN
    /// step and a lower initiation penalty (+0.98 vs +1.03), so more terminal G/C ⇒ ΔG decreases.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Calculate3PrimeStability_MoreTerminalGc_NotLessStable()
    {
        // Inner trinucleotide is fixed; only the two terminal bases vary A/T → G/C.
        var arb = (from inner in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(3)
                   select new string(inner)).ToArbitrary();
        return Prop.ForAll(arb, inner =>
        {
            double dgAt = PrimerDesigner.Calculate3PrimeStability("A" + inner + "T");
            double dgGc = PrimerDesigner.Calculate3PrimeStability("G" + inner + "C");
            return (dgGc <= dgAt + 1e-9)
                .Label($"GC-terminal not ≤ AT-terminal for inner '{inner}': {dgGc} vs {dgAt}");
        });
    }

    /// <summary>
    /// INV-D (Determinism): repeated evaluation of the same sequence is bit-identical.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Calculate3PrimeStability_IsDeterministic()
    {
        return Prop.ForAll(StabilitySequenceArbitrary(), seq =>
        {
            double a = PrimerDesigner.Calculate3PrimeStability(seq);
            double b = PrimerDesigner.Calculate3PrimeStability(seq);
            return a.Equals(b).Label($"non-deterministic ΔG for '{seq}': {a} vs {b}");
        });
    }

    /// <summary>
    /// Evidence anchors: the two canonical 5-mers documented in §2.2 with literal published values —
    /// GCGCG = −6.86 kcal/mol (most stable) and TATAT = −0.86 (least stable). A shared-but-wrong
    /// constant in production+oracle could pass the oracle property; these absolute values cannot.
    /// </summary>
    [TestCase("GCGCG", -6.86, TestName = "SantaLucia: GCGCG → −6.86 (most stable)")]
    [TestCase("TATAT", -0.86, TestName = "SantaLucia: TATAT → −0.86 (least stable)")]
    [Category("Property")]
    public void Calculate3PrimeStability_CanonicalAnchors(string seq, double expected)
    {
        double actual = PrimerDesigner.Calculate3PrimeStability(seq);
        Assert.That(actual, Is.EqualTo(expected).Within(0.01),
            $"ΔG for '{seq}' must match the documented SantaLucia value {expected} kcal/mol");
    }



    /// <summary>
    /// INV-01 (R): HasHairpinPotential returns false whenever the sequence is shorter than the
    /// structural minimum 2·minStemLength + minLoopLength — a valid stem-loop cannot fit. Exercised
    /// over randomized sequences AND randomized stem/loop knobs. Source: §2.4 INV-01.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HasHairpinPotential_BelowMinStructureSize_IsFalse()
    {
        return Prop.ForAll(HairpinKnobArbitrary(), t =>
        {
            var (seq, stem, loop) = t;
            if (seq.Length >= 2 * stem + loop)
                return true.ToProperty(); // property only constrains the sub-minimum regime.
            return (!PrimerDesigner.HasHairpinPotential(seq, stem, loop))
                .Label($"hairpin reported below min size: len {seq.Length} < 2·{stem}+{loop} for '{seq}'");
        });
    }

    /// <summary>
    /// INV-M (monotonicity in permissiveness): a shorter required stem is strictly more permissive,
    /// so if a hairpin is detected at stem = k it must still be detected at stem = k−1 (loop fixed).
    /// Detection sets are nested: lowering the stem threshold never turns true into false.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HasHairpinPotential_LowerStem_NeverLosesDetection()
    {
        var arb = (from stem in Gen.Choose(3, 6)
                   from loop in Gen.Choose(1, 4)
                   from len in Gen.Choose(2 * stem + loop, 60)
                   from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
                   select (new string(chars), stem, loop)).ToArbitrary();

        return Prop.ForAll(arb, t =>
        {
            var (seq, stem, loop) = t;
            bool atK = PrimerDesigner.HasHairpinPotential(seq, stem, loop);
            bool atKMinus1 = PrimerDesigner.HasHairpinPotential(seq, stem - 1, loop);
            // detected at stricter stem ⇒ must be detected at more permissive (shorter) stem.
            return (!atK || atKMinus1)
                .Label($"monotonicity broken for '{seq}' loop={loop}: stem {stem}={atK}, stem {stem - 1}={atKMinus1}");
        });
    }

    /// <summary>
    /// INV-P (case-insensitivity): upper-, lower-, and mixed-case spellings of the same sequence
    /// give the identical hairpin verdict (the source upper-cases before comparison). Source: §3.3.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HasHairpinPotential_IsCaseInsensitive()
    {
        return Prop.ForAll(ValidPrimerArbitrary(1, 40), seq =>
        {
            char[] m = seq.ToCharArray();
            for (int i = 0; i < m.Length; i += 2) m[i] = char.ToLowerInvariant(m[i]);
            bool up = PrimerDesigner.HasHairpinPotential(seq.ToUpperInvariant());
            bool lo = PrimerDesigner.HasHairpinPotential(seq.ToLowerInvariant());
            bool mi = PrimerDesigner.HasHairpinPotential(new string(m));
            return (up == lo && up == mi)
                .Label($"case sensitivity for '{seq}': up={up}, lo={lo}, mixed={mi}");
        });
    }

    /// <summary>
    /// INV-D (Determinism): the hairpin verdict is identical across repeated calls.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HasHairpinPotential_IsDeterministic()
    {
        return Prop.ForAll(HairpinKnobArbitrary(), t =>
        {
            var (seq, stem, loop) = t;
            bool a = PrimerDesigner.HasHairpinPotential(seq, stem, loop);
            bool b = PrimerDesigner.HasHairpinPotential(seq, stem, loop);
            return (a == b).Label($"non-deterministic hairpin for '{seq}'");
        });
    }

    /// <summary>
    /// Constructed positive anchor: a 5' stem, a loop, then the reverse complement of the stem as
    /// the 3' stem forms a genuine stem-loop. Stem "GGGG" (positions 0-3), loop "AAA" (3 nt ≥
    /// minLoop), 3' stem "CCCC" = revcomp("GGGG"). With minStem=4/minLoop=3 this must be detected.
    /// Source: §2.2 stem-loop model; §6.1.
    /// </summary>
    [Test]
    [Category("Property")]
    public void HasHairpinPotential_EngineeredStemLoop_ReturnsTrue()
    {
        // 5' stem GGGG | loop AAA | 3' stem CCCC (= reverse complement of GGGG read 3'→5').
        const string hairpin = "GGGGAAACCCC";
        Assert.That(PrimerDesigner.HasHairpinPotential(hairpin, minStemLength: 4, minLoopLength: 3),
            Is.True, "engineered stem(GGGG)+loop(AAA)+revcomp-stem(CCCC) must form a hairpin");
    }

    /// <summary>
    /// Constructed negative anchor: a homopolymer has no self-complementary stem (A never pairs with
    /// A), so no hairpin is possible regardless of length. Source: §6.1 "No self-complementary regions".
    /// </summary>
    [Test]
    [Category("Property")]
    public void HasHairpinPotential_Homopolymer_ReturnsFalse()
    {
        Assert.That(PrimerDesigner.HasHairpinPotential(new string('A', 40)),
            Is.False, "a poly-A run has no complementary stem and cannot form a hairpin");
    }



    /// <summary>
    /// INV-Edge: an empty or null primer on either side yields false — no terminal complementarity
    /// can be evaluated. Source: §3.3 / §6.1.
    /// </summary>
    [Test]
    [Category("Property")]
    public void HasPrimerDimer_EmptyOrNull_ReturnsFalse()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.HasPrimerDimer("", "ACGTACGT"), Is.False);
            Assert.That(PrimerDesigner.HasPrimerDimer("ACGTACGT", ""), Is.False);
            Assert.That(PrimerDesigner.HasPrimerDimer(null!, "ACGTACGT"), Is.False);
            Assert.That(PrimerDesigner.HasPrimerDimer("ACGTACGT", null!), Is.False);
        });
    }

    /// <summary>
    /// INV-M (monotonicity in threshold): a higher minComplementarity is stricter, so if a dimer is
    /// flagged at the high threshold it must also be flagged at the lower one (same primers). Raising
    /// the threshold can never INCREASE detections. Source: §2.2 (threshold on terminal window).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HasPrimerDimer_HigherThreshold_NeverIncreasesDetection()
    {
        return Prop.ForAll(DimerThresholdArbitrary(), t =>
        {
            var (p1, p2, lo, hi) = t;
            bool atHi = PrimerDesigner.HasPrimerDimer(p1, p2, hi);
            bool atLo = PrimerDesigner.HasPrimerDimer(p1, p2, lo);
            // detected at the stricter (higher) threshold ⇒ must be detected at the looser (lower) one.
            return (!atHi || atLo)
                .Label($"threshold monotonicity broken for '{p1}'/'{p2}': thr {hi}={atHi}, thr {lo}={atLo}");
        });
    }

    /// <summary>
    /// INV-D (Determinism): the dimer verdict is identical across repeated calls.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property HasPrimerDimer_IsDeterministic()
    {
        return Prop.ForAll(DimerThresholdArbitrary(), t =>
        {
            var (p1, p2, _, thr) = t;
            bool a = PrimerDesigner.HasPrimerDimer(p1, p2, thr);
            bool b = PrimerDesigner.HasPrimerDimer(p1, p2, thr);
            return (a == b).Label($"non-deterministic dimer for '{p1}'/'{p2}'");
        });
    }

    /// <summary>
    /// Constructed positive anchor with independent recomputation: a 3′ primer-dimer forms when the
    /// 3′ terminus of primer1 anneals ANTIPARALLEL to the 3′ terminus of primer2, i.e. primer1's
    /// 3′ k-mer is the reverse complement of primer2's 3′ k-mer (equivalently: it equals the FIRST k
    /// bases of revcomp(primer2)). Primer3's alignment-mode <c>compl_end</c> (end-anchored
    /// <c>align(p1, revcomp(p2))</c>, +1 per pair) then scores k. We build such a pair with k = 8 and
    /// verify (a) the independently recomputed antiparallel pair count is 8, (b) production's
    /// <c>compl_end</c> is 8 (primer3-py 2.3.1 <c>check_primers</c>,
    /// <c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0</c>: <c>PRIMER_PAIR_0_COMPL_END = 8.0</c>) and (c) the
    /// dimer is flagged at minComplementarity = 4. Source: Primer_Structure_Analysis.md §2.2 / §5.2.
    /// </summary>
    [Test]
    [Category("Property")]
    public void HasPrimerDimer_EngineeredComplementaryEnds_ReturnsTrue()
    {
        // primer2 fixed; choose primer1 so its 3' 8-mer == first 8 of revcomp(primer2), i.e. the
        // reverse complement of primer2's 3' 8-mer ⇒ the two 3' ends pair antiparallel over 8 bases.
        const string primer2 = "TTTTTTTTGGCCAATT";
        string rc2 = RevComp(primer2);                 // independent reverse complement
        string window2 = rc2.Substring(0, 8);          // first 8 bases of revcomp(primer2)
        string primer1 = "ACGTACGT" + window2;          // length 16, last 8 = revcomp(primer2's 3' 8-mer)

        // Independent recomputation: primer1[n1-8+i] pairs antiparallel with primer2[n2-1-i].
        string end1 = primer1.Substring(primer1.Length - 8);
        int comp = 0;
        for (int i = 0; i < 8; i++)
            if (IsWatsonCrick(end1[i], primer2[primer2.Length - 8 + (7 - i)])) comp++;

        Assert.Multiple(() =>
        {
            Assert.That(comp, Is.EqualTo(8), "engineered 3' ends must pair antiparallel over all 8 bases");
            Assert.That(PrimerDesigner.CalculatePrimerDimerEndComplementarity(primer1, primer2), Is.EqualTo(8.0),
                "Primer3 compl_end of the engineered pair (primer3-py 2.3.1: 8.0)");
            Assert.That(PrimerDesigner.HasPrimerDimer(primer1, primer2, minComplementarity: 4),
                Is.True, "fully complementary 3' ends must be flagged as a primer-dimer");
        });
    }

    /// <summary>Independent GC% oracle: 100 × (#G + #C) / length, computed from raw base counts.</summary>
    private static double ExpectedGcPercent(string seq)
    {
        if (seq.Length == 0) return 0;
        var (_, gc) = CountAtGc(seq);
        return 100.0 * gc / seq.Length;
    }



    /// <summary>
    /// Independent GC fraction oracle: (#G + #C) / length on the uppercased sequence, computed
    /// from raw base counts — NOT routed through any production helper. Returns 0 for empty input.
    /// (The production code stores GcContent as a fraction in [0,1], cf. INV-02.)
    /// </summary>
    private static double ExpectedGcFraction(string seq)
    {
        if (seq.Length == 0) return 0;
        var (_, gc) = CountAtGc(seq);
        return (double)gc / seq.Length;
    }



    /// <summary>
    /// Independent reverse-complement oracle over an uppercased ACGT string, built char-by-char
    /// from the complement map A↔T, C↔G and then reversed — deliberately NOT routed through any
    /// production helper, so a wrong complement table is caught. Non-ACGT chars are mapped to 'N'.
    /// </summary>
    private static string OracleReverseComplement(string seq)
    {
        var chars = new char[seq.Length];
        for (int i = 0; i < seq.Length; i++)
        {
            char c = char.ToUpperInvariant(seq[i]);
            chars[seq.Length - 1 - i] = c switch
            {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                _ => 'N'
            };
        }
        return new string(chars);
    }

    /// <summary>
    /// Independent self-complementarity oracle (INV-02): fraction of positions i where
    /// <c>seq[i] == reverseComplement(seq)[i]</c>, over the uppercased probe. Mirrors the
    /// definition in ProbeDesigner §CalculateSelfComplementarity but with its OWN revcomp.
    /// </summary>
    private static double OracleSelfComplementarity(string seq)
    {
        string up = seq.ToUpperInvariant();
        string rc = OracleReverseComplement(up);
        int matches = 0;
        for (int i = 0; i < up.Length; i++)
            if (up[i] == rc[i]) matches++;
        return up.Length == 0 ? 0.0 : matches / (double)up.Length;
    }

    /// <summary>
    /// Independent off-target counter (INV-03): substitution-only, fixed-length, ungapped sliding
    /// window over each reference (both probe and reference uppercased), counting every start
    /// position whose Hamming distance to the probe over the window is ≤ maxMismatches, summed
    /// across all references. Re-implemented from the doc §4.1 / §5.2 — NOT from production
    /// internals — so it is a true oracle for <c>FindApproximateMatches</c>.
    /// </summary>
    private static int OracleOffTargetHits(string probe, IEnumerable<string> references, int maxMismatches)
    {
        string p = probe.ToUpperInvariant();
        int total = 0;
        foreach (var reference in references)
        {
            string text = reference.ToUpperInvariant();
            for (int i = 0; i + p.Length <= text.Length; i++)
            {
                int mismatches = 0;
                for (int j = 0; j < p.Length; j++)
                {
                    if (text[i + j] != p[j]) mismatches++;
                    if (mismatches > maxMismatches) break;
                }
                if (mismatches <= maxMismatches) total++;
            }
        }
        return total;
    }

    /// <summary>
    /// Independent specificity map (INV-01 / INV-04 / §2.2): 0 hits ⇒ 0.0, 1 hit ⇒ 1.0,
    /// N&gt;1 hits ⇒ 1.0/N. Re-derived purely from a hit count, not from production.
    /// </summary>
    private static double OracleSpecificity(int offTargetHits) =>
        offTargetHits == 0 ? 0.0 : offTargetHits == 1 ? 1.0 : 1.0 / offTargetHits;

    /// <summary>
    /// Generates a probe-validation scenario: a random ACGT probe (8-20 bp), a set of 0-4 reference
    /// sequences, and a mismatch tolerance. References are deliberately MIXED in kind so the
    /// off-target counter is exercised across the interesting regimes:
    /// (a) references that EMBED the probe verbatim 0/1/several times (exact hits),
    /// (b) references that embed a NEAR-match (the probe with 1-3 substitutions — within tolerance),
    /// (c) unrelated filler (no hits). Lengths are kept short so hundreds of cases run quickly.
    /// </summary>
    private static Arbitrary<(string probe, string[] refs, int maxMismatches)>
        ProbeValidationScenarioArbitrary() =>
        (from plen in Gen.Choose(8, 20)
         from pchars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(plen)
         let probe = new string(pchars)
         from maxMismatches in Gen.Choose(0, 4)
         from refCount in Gen.Choose(0, 4)
         from refs in BuildReferenceGen(probe, maxMismatches).ArrayOf(refCount)
         select (probe, refs, maxMismatches)).ToArbitrary();

    /// <summary>
    /// Builds a single reference sequence around the probe: random ACGT flanks plus, with varying
    /// probability, an embedded exact copy of the probe or a perturbed copy carrying a controlled
    /// number of substitutions (0..probeLen). The perturbed copy may or may not be within tolerance,
    /// so both "counted" and "not counted" near-matches are produced.
    /// </summary>
    private static Gen<string> BuildReferenceGen(string probe, int maxMismatches)
    {
        return
            from leftLen in Gen.Choose(0, 8)
            from rightLen in Gen.Choose(0, 8)
            from left in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(leftLen)
            from right in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(rightLen)
            from kind in Gen.Choose(0, 3) // 0=no insert, 1=exact, 2=near-match, 3=double exact
            from subs in Gen.Choose(0, probe.Length)
            let middle = kind switch
            {
                1 => probe,
                2 => PerturbProbe(probe, subs),
                3 => probe + probe,
                _ => string.Empty
            }
            select new string(left) + middle + new string(right);
    }

    /// <summary>
    /// Returns a copy of the probe with the first <paramref name="count"/> positions substituted to
    /// a DIFFERENT base (deterministic 'A'→'C', else→'A'), producing exactly <paramref name="count"/>
    /// guaranteed mismatches at known positions.
    /// </summary>
    private static string PerturbProbe(string probe, int count)
    {
        var chars = probe.ToCharArray();
        for (int i = 0; i < count && i < chars.Length; i++)
            chars[i] = chars[i] == 'A' ? 'C' : 'A';
        return new string(chars);
    }

    /// <summary>
    /// CheckSpecificity suffix-tree anchor (0 / 1 / 1/hits exact-hit map): a genome containing the
    /// probe EXACTLY TWICE ⇒ 0.5; a unique probe ⇒ 1.0; an absent probe ⇒ 0.0. Uses a real
    /// <c>global::SuffixTree.SuffixTree.Build(...)</c> index, with the occurrence count cross-checked
    /// via <c>FindAllOccurrences</c>. Source: ProbeDesigner.cs §CheckSpecificity.
    /// </summary>
    [Test]
    [Category("Property")]
    public void CheckSpecificity_SuffixTree_ExactHitMap()
    {
        // "GATTACA" appears twice; "CCCCCCC" is absent; "ACAGTC" unique.
        const string genome = "GATTACAGTCAAGATTACATTTT";
        var index = global::SuffixTree.SuffixTree.Build(genome);

        Assert.Multiple(() =>
        {
            int twiceCount = index.FindAllOccurrences("GATTACA").Count;
            Assert.That(twiceCount, Is.EqualTo(2), "GATTACA occurs exactly twice");
            Assert.That(ProbeDesigner.CheckSpecificity("GATTACA", index),
                Is.EqualTo(0.5).Within(1e-9), "2 hits ⇒ 1/2");

            Assert.That(index.FindAllOccurrences("CCCCCCC").Count, Is.EqualTo(0), "absent probe");
            Assert.That(ProbeDesigner.CheckSpecificity("CCCCCCC", index),
                Is.EqualTo(0.0).Within(1e-9), "0 hits ⇒ 0.0");

            Assert.That(index.FindAllOccurrences("ACAGTC").Count, Is.EqualTo(1), "unique window");
            Assert.That(ProbeDesigner.CheckSpecificity("ACAGTC", index),
                Is.EqualTo(1.0).Within(1e-9), "1 hit ⇒ 1.0");
        });
    }



    // CalculateMeltingTemperatureNN — SantaLucia & Hicks (2004) unified nearest-neighbour duplex Tm with
    // an Owczarzy (2004) monovalent salt correction. CalculateMeltingTemperatureNNMismatch scores an
    // imperfect duplex with the internal-mismatch NN parameters.

    /// <summary>
    /// INV-1 (R): the nearest-neighbour Tm is a finite number for every valid ACGT primer of length ≥ 2
    /// (a duplex needs at least one NN stack).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property NnTm_IsFinite_ForValidPrimer()
    {
        return Prop.ForAll(ValidPrimerArbitrary(2, 40), primer =>
        {
            double tm = PrimerDesigner.CalculateMeltingTemperatureNN(primer);
            return double.IsFinite(tm).Label($"NN Tm={tm} must be finite for '{primer}'");
        });
    }

    /// <summary>
    /// INV-2 (M): raising the monovalent salt concentration never lowers the duplex Tm — Na⁺ shields the
    /// backbone charge and stabilises the duplex (Owczarzy 2004). Compared over the calibrated [50 mM,1 M] range.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property NnTm_HigherSodium_NotLowerTm()
    {
        var gen = (from primer in ValidPrimerGen(6, 40)
                   from naLoMm in Gen.Choose(50, 1000)
                   from naHiMm in Gen.Choose(50, 1000)
                   select (primer, lo: Math.Min(naLoMm, naHiMm) / 1000.0, hi: Math.Max(naLoMm, naHiMm) / 1000.0)).ToArbitrary();

        return Prop.ForAll(gen, t =>
        {
            double tmLo = PrimerDesigner.CalculateMeltingTemperatureNN(t.primer, sodiumMolar: t.lo);
            double tmHi = PrimerDesigner.CalculateMeltingTemperatureNN(t.primer, sodiumMolar: t.hi);
            return (tmHi >= tmLo - 1e-9).Label($"[Na+] {t.lo}→{tmLo:F2}, {t.hi}→{tmHi:F2} must not lower Tm");
        });
    }

    /// <summary>
    /// INV-3 (M): a perfectly complementary duplex melts no lower than the same duplex carrying an internal
    /// <b>destabilising</b> mismatch — breaking a Watson-Crick pair generally replaces two favourable stacks
    /// with a destabilising mismatch.
    /// <para>
    /// DOMAIN: the introduced mismatch is restricted to types that are thermodynamically destabilising. The
    /// internal <b>G·A</b> sheared/imino mismatch is a documented EXCEPTION — Allawi &amp; SantaLucia (1998)
    /// Biochemistry 37:9435 report anomalously stable G·A NN stacks (e.g. the GG/CG term is ΔH°=−6.0,
    /// ΔS°=−15.8, more stabilising than the matched stack it replaces), so a single internal G·A can raise the
    /// duplex Tm above the perfect match in GC-rich contexts. That is the NN model behaving correctly, not a
    /// defect, so the monotonicity theorem is asserted only on the destabilising-mismatch domain (G·A excluded).
    /// Verified empirically: excluding G·A removes every counterexample across 150k random duplexes; including
    /// it produces ~0.08% violations, all G·A. The G·A stack parameters themselves are validated against the
    /// published equation by <c>PrimerDesigner_NearestNeighborTm_Tests</c>.
    /// </para>
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property NnTmMismatch_PerfectMatch_NotBelowSingleDestabilisingMismatch()
    {
        var gen = (from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= 8 && a.Length <= 24)
                   let top = new string(chars)
                   from pos in Gen.Choose(2, top.Length - 3) // internal position, away from both ends
                   select (top, pos)).ToArbitrary();

        return Prop.ForAll(gen, t =>
        {
            // The mismatch model aligns the two strands index-to-index (top[i] pairs with bottom[i]),
            // so a perfect duplex's bottom strand is the PER-POSITION complement of top (no reversal).
            string perfect = OracleComplement(t.top);
            var mm = perfect.ToCharArray();
            mm[t.pos] = mm[t.pos] switch { 'A' => 'C', 'C' => 'A', 'G' => 'T', _ => 'G' }; // break the WC pair at top[pos]

            // Skip the anomalously-stable internal G·A mismatch (top[pos] paired with mm[pos] = a G·A pair):
            // it is documented to stabilise, so it lies outside this monotonicity invariant's domain.
            char topBase = char.ToUpperInvariant(t.top[t.pos]);
            char botBase = mm[t.pos];
            bool isGaShearedMismatch = (topBase == 'A' && botBase == 'G') || (topBase == 'G' && botBase == 'A');
            if (isGaShearedMismatch) return true.ToProperty();

            double tmPerfect = PrimerDesigner.CalculateMeltingTemperatureNNMismatch(t.top, perfect);
            double tmMismatch = PrimerDesigner.CalculateMeltingTemperatureNNMismatch(t.top, new string(mm));
            if (!double.IsFinite(tmPerfect) || !double.IsFinite(tmMismatch)) return true.ToProperty();
            return (tmPerfect >= tmMismatch - 1e-9)
                .Label($"perfect Tm={tmPerfect:F2} < mismatch Tm={tmMismatch:F2}");
        });
    }

    /// <summary>INV-4 (D): the nearest-neighbour Tm is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property NnTm_IsDeterministic()
    {
        return Prop.ForAll(ValidPrimerArbitrary(2, 40), primer =>
            (PrimerDesigner.CalculateMeltingTemperatureNN(primer) == PrimerDesigner.CalculateMeltingTemperatureNN(primer))
                .Label("CalculateMeltingTemperatureNN must be deterministic"));
    }



    // FindMostStableHairpin returns the minimum-ΔG°37 intramolecular hairpin (SantaLucia & Hicks 2004) or
    // null when no stem ≥2 bp can close a loop ≥3 nt. NOTE: it returns the *most stable* hairpin even when
    // that minimum ΔG is slightly positive (weak stem + large loop), so ΔG ≤ 0 is asserted on the
    // meaningful domain of a genuine GC stem rather than universally.

    /// <summary>Per-position DNA complement (no reversal) — the bottom strand of an index-aligned duplex.</summary>
    private static string OracleComplement(string seq)
    {
        var chars = new char[seq.Length];
        for (int i = 0; i < seq.Length; i++)
            chars[i] = char.ToUpperInvariant(seq[i]) switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N' };
        return new string(chars);
    }

    /// <summary>A GC-stem hairpin: a G/C stem of length L, an A-loop (3..6 nt), then the stem's reverse complement.</summary>
    private static Gen<(string seq, int stemLen)> GcStemHairpinGen(int minStem) =>
        from stemLen in Gen.Choose(minStem, 9)
        from stemChars in Gen.Elements('G', 'C').ArrayOf(stemLen)
        from loopLen in Gen.Choose(3, 6)
        let stem = new string(stemChars)
        select (stem + new string('A', loopLen) + OracleReverseComplement(stem), stemLen);

    /// <summary>
    /// INV-1 (R, well-formed): a sequence built as GC-stem + loop + reverse-complement has a hairpin whose
    /// stem ≥ minStemLength, loop ≥ 3 nt, and positions are in bounds.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hairpin_PlantedGcStem_IsWellFormed()
    {
        return Prop.ForAll(GcStemHairpinGen(2).ToArbitrary(), t =>
        {
            var hp = PrimerDesigner.FindMostStableHairpin(t.seq, minStemLength: 2);
            if (hp is null) return false.Label($"no hairpin found in planted construct '{t.seq}'");
            var h = hp.Value;
            bool ok = h.StemLength >= 2 && h.LoopSize >= 3
                      && h.StemStart >= 0 && h.StemEnd < t.seq.Length && h.StemStart < h.StemEnd;
            return ok.Label($"stem={h.StemLength}, loop={h.LoopSize}");
        });
    }

    /// <summary>
    /// INV-1b (R, ΔG ≤ 0): a genuinely stable hairpin (GC stem ≥ 5 bp, enough nearest-neighbour stacks to
    /// overcome the loop-initiation penalty) has ΔG°37 ≤ 0. A weaker 2-bp stem need not — the method returns
    /// the most stable hairpin even when its minimum ΔG is slightly positive, so this is the "real stem" domain.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hairpin_StableGcStem_HasNegativeDeltaG()
    {
        return Prop.ForAll(GcStemHairpinGen(5).ToArbitrary(), t =>
        {
            var hp = PrimerDesigner.FindMostStableHairpin(t.seq, minStemLength: 2);
            if (hp is null) return false.Label($"no hairpin found in planted construct '{t.seq}'");
            return (hp.Value.DeltaG37 <= 1e-9).Label($"stem≥5 ΔG37={hp.Value.DeltaG37:F3} must be ≤ 0");
        });
    }

    /// <summary>
    /// INV-2 (M): extending the stem by one GC pair makes the best hairpin more stable — ΔG°37 decreases
    /// (one extra stabilising nearest-neighbour stack), at a fixed loop.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Hairpin_LongerStem_LowerDeltaG()
    {
        var gen = (from stemLen in Gen.Choose(2, 7)
                   from stemChars in Gen.Elements('G', 'C').ArrayOf(stemLen)
                   from extra in Gen.Elements('G', 'C')
                   let stem = new string(stemChars)
                   let loop = new string('A', 4)
                   select (shortStem: stem, longStem: stem + extra, loop)).ToArbitrary();

        return Prop.ForAll(gen, t =>
        {
            string seqShort = t.shortStem + t.loop + OracleReverseComplement(t.shortStem);
            string seqLong = t.longStem + t.loop + OracleReverseComplement(t.longStem);
            var hpShort = PrimerDesigner.FindMostStableHairpin(seqShort, minStemLength: 2);
            var hpLong = PrimerDesigner.FindMostStableHairpin(seqLong, minStemLength: 2);
            if (hpShort is null || hpLong is null) return false.Label("planted hairpin missing");
            return (hpLong.Value.DeltaG37 <= hpShort.Value.DeltaG37 + 1e-9)
                .Label($"longer stem not more stable: ΔG(L+1)={hpLong.Value.DeltaG37:F3} > ΔG(L)={hpShort.Value.DeltaG37:F3}");
        });
    }

    /// <summary>INV-3 (D): the most-stable-hairpin search is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property Hairpin_IsDeterministic()
    {
        return Prop.ForAll(ValidPrimerArbitrary(4, 40), seq =>
            (PrimerDesigner.FindMostStableHairpin(seq) == PrimerDesigner.FindMostStableHairpin(seq))
                .Label("FindMostStableHairpin must be deterministic"));
    }



    // FindMostStableDimer returns the most stable intermolecular duplex (Primer3/ntthal NN model) or null
    // when no contiguous run of ≥2 bp exists. As with hairpins, the most stable dimer can have a slightly
    // positive ΔG for a weak 2-bp run, so ΔG ≤ 0 is asserted on a genuine GC complementary run (length ≥ 3).

    /// <summary>
    /// INV-1 (R): two strands that are reverse complements over a GC run (length 3..8) form a stable dimer —
    /// non-null, ≥ 2 base pairs, in-bounds offsets, and ΔG°37 ≤ 0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Dimer_PlantedGcRun_IsStableAndWellFormed()
    {
        var gen = (from runLen in Gen.Choose(3, 8)
                   from runChars in Gen.Elements('G', 'C').ArrayOf(runLen)
                   let s1 = new string(runChars)
                   select (s1, s2: OracleReverseComplement(s1))).ToArbitrary();

        return Prop.ForAll(gen, t =>
        {
            var d = PrimerDesigner.FindMostStableDimer(t.s1, t.s2);
            if (d is null) return false.Label($"no dimer found for complementary GC run '{t.s1}'/'{t.s2}'");
            var dim = d.Value;
            bool ok = dim.BasePairs >= 2 && dim.Strand1Start >= 0 && dim.Strand2Start >= 0
                      && dim.DeltaG37 <= 1e-9;
            return ok.Label($"bp={dim.BasePairs}, ΔG37={dim.DeltaG37:F3}");
        });
    }

    /// <summary>
    /// INV-2 (M): a longer complementary run yields a more stable dimer — extending the GC run by one base
    /// (and its complement) lowers ΔG°37 (one extra stabilising stack).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Dimer_LongerRun_LowerDeltaG()
    {
        var gen = (from runLen in Gen.Choose(3, 7)
                   from runChars in Gen.Elements('G', 'C').ArrayOf(runLen)
                   from extra in Gen.Elements('G', 'C')
                   let shortRun = new string(runChars)
                   select (shortRun, longRun: shortRun + extra)).ToArbitrary();

        return Prop.ForAll(gen, t =>
        {
            var dShort = PrimerDesigner.FindMostStableDimer(t.shortRun, OracleReverseComplement(t.shortRun));
            var dLong = PrimerDesigner.FindMostStableDimer(t.longRun, OracleReverseComplement(t.longRun));
            if (dShort is null || dLong is null) return false.Label("planted dimer missing");
            return (dLong.Value.DeltaG37 <= dShort.Value.DeltaG37 + 1e-9)
                .Label($"longer run not more stable: ΔG(L+1)={dLong.Value.DeltaG37:F3} > ΔG(L)={dShort.Value.DeltaG37:F3}");
        });
    }

    /// <summary>INV-3 (D): the most-stable-dimer search is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property Dimer_IsDeterministic()
    {
        var gen = (from a in ValidPrimerGen(4, 30) from b in ValidPrimerGen(4, 30) select (a, b)).ToArbitrary();
        return Prop.ForAll(gen, t =>
            (PrimerDesigner.FindMostStableDimer(t.a, t.b) == PrimerDesigner.FindMostStableDimer(t.a, t.b))
                .Label("FindMostStableDimer must be deterministic"));
    }



    // CalculateMeltingTemperatureNNLna applies the McTigue (2004) LNA nearest-neighbour increments, which
    // are all stabilising, so locking a base as LNA never lowers the duplex Tm. EvaluateMgbProbeDesign
    // returns the qualitative Kutyavin (2000) MGB design rules (a boolean + a reasons list).

    /// <summary>
    /// INV-1 (R): adding an <b>isolated</b> LNA-<b>C</b> substitution never lowers the duplex Tm: Tm(S ∪ {p}) ≥ Tm(S) when the
    /// locked base at p is C and neither neighbour (p−1, p+1) is already locked.
    /// <para>
    /// Two facts bound the claim. (1) The published LNA nearest-neighbour increments are mixed-sign by design (McTigue,
    /// Peterson &amp; Kahn 2004; e.g. the G_L-A step ΔΔH°=+3.16, ΔΔS°=+10.5), so the blanket "any added LNA never lowers Tm"
    /// is false: a single LNA-A, and LNA-G in a few steps, can shave a few tenths of a °C off in some sequence contexts.
    /// LNA-C is the reliably stabilising case. (2) The model used here by default, Owczarzy et al. (2011) (the MELTING 5
    /// default), publishes separate complete step parameters for a step with one locked base ("single") and for a step with
    /// two locked bases ("consecutive"). Adding an LNA next to an existing LNA therefore does not add one stabilising
    /// increment; it swaps two "single" steps for "consecutive" steps, and that can lower Tm. Observed with this
    /// implementation (found by the CI property run of 2026-11-04; not independently re-checked against MELTING 5):
    /// CCCCGCTCGCCCCGCC with LNA at {1,2,3,4,5,7,8,10,11,12,13,14}, adding LNA-C at 9, moves Tm 100.30 → 100.21 °C.
    /// The generator therefore only adds an LNA whose neighbours are not LNA; for those no violation was found in
    /// 1,000,000 random duplexes (the earlier ~11k-sample check used the unrestricted generator and missed the
    /// consecutive-LNA case).
    /// </para>
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property LnaTm_AddingLnaCPosition_NeverLowersTm()
    {
        var gen = (from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= 8 && a.Length <= 30)
                   let seq = new string(chars)
                   where seq.Contains('C')
                   from p in Gen.Elements(Enumerable.Range(0, seq.Length).Where(i => seq[i] == 'C').ToArray())
                   from baseFlags in Gen.Elements(true, false).ArrayOf(seq.Length)
                   where (p == 0 || !baseFlags[p - 1]) && (p == seq.Length - 1 || !baseFlags[p + 1])
                   select (seq, p, baseSet: Enumerable.Range(0, seq.Length).Where(i => baseFlags[i] && i != p).ToList()))
                   .ToArbitrary();

        return Prop.ForAll(gen, t =>
        {
            double without = PrimerDesigner.CalculateMeltingTemperatureNNLna(t.seq, t.baseSet);
            double with = PrimerDesigner.CalculateMeltingTemperatureNNLna(t.seq, t.baseSet.Append(t.p).ToList());
            if (!double.IsFinite(without) || !double.IsFinite(with)) return true.ToProperty();
            return (with >= without - 1e-9).Label($"adding LNA-C@{t.p} lowered Tm: {without:F2} → {with:F2}");
        });
    }

    /// <summary>INV-2 (D): the LNA Tm is deterministic for a fixed sequence and LNA position set.</summary>
    [FsCheck.NUnit.Property]
    public Property LnaTm_IsDeterministic()
    {
        return Prop.ForAll(ValidPrimerArbitrary(8, 30), seq =>
        {
            var lna = new[] { 0, seq.Length / 2 };
            // double.Equals treats NaN == NaN as true: NaN (no LNA thermodynamics) is still deterministic.
            return PrimerDesigner.CalculateMeltingTemperatureNNLna(seq, lna)
                    .Equals(PrimerDesigner.CalculateMeltingTemperatureNNLna(seq, lna))
                .Label("CalculateMeltingTemperatureNNLna must be deterministic");
        });
    }

    /// <summary>
    /// INV-3 (MGB rules → boolean + reasons): EvaluateMgbProbeDesign always returns a non-empty guidance
    /// list with the 3'-attachment rule, reports the input length, uppercases the sequence, and its
    /// LengthInMgbRange boolean is consistent with whether an "outside the MGB window" warning was emitted.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Mgb_ReturnsBooleanAndConsistentReasons()
    {
        return Prop.ForAll(ValidPrimerArbitrary(1, 40), seq =>
        {
            var d = ProbeDesigner.EvaluateMgbProbeDesign(seq);
            bool warned = d.Guidance.Any(g => g.Contains("outside", StringComparison.OrdinalIgnoreCase));
            bool ok = d.Guidance.Count >= 1
                      && d.Guidance.Any(g => g.Contains("3'"))
                      && d.MgbAttachmentEnd == "3'"
                      && d.Length == seq.Length
                      && d.Sequence == seq.ToUpperInvariant()
                      && d.LengthInMgbRange == !warned;
            return ok.Label($"len={d.Length}, inRange={d.LengthInMgbRange}, warned={warned}, reasons={d.Guidance.Count}");
        });
    }



    // ComputeKarlinAltschul: E = K·m·n·e^{−λS}, bit score S' = (λS − ln K)/ln 2 (Karlin & Altschul 1990;
    // Altschul et al. 1990). All of K, m, n, e^{−λS} are positive, so E ≥ 0.

    private static Arbitrary<(double score, int m, long n)> KaProblemArbitrary() =>
        (from scoreCenti in Gen.Choose(0, 10000)
         from m in Gen.Choose(1, 5000)
         from nKb in Gen.Choose(1, 1_000_000)
         select (scoreCenti / 100.0, m, (long)nKb * 1000)).ToArbitrary();

    /// <summary>INV-1 (R): the Karlin–Altschul E-value is non-negative (E = K·m·n·e^{−λS}, all factors ≥ 0).</summary>
    [FsCheck.NUnit.Property]
    public Property KarlinAltschul_EValue_IsNonNegative()
    {
        return Prop.ForAll(KaProblemArbitrary(), t =>
        {
            double e = ProbeDesigner.ComputeKarlinAltschul(t.score, t.m, t.n).EValue;
            return (e >= 0.0).Label($"E-value {e} must be ≥ 0");
        });
    }

    /// <summary>
    /// INV-2 (M): a higher raw alignment score yields a higher bit score but a LOWER E-value — E decays as
    /// e^{−λS} while the bit score grows linearly in S (Altschul et al. 1990).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property KarlinAltschul_HigherScore_LowerEValue_HigherBitScore()
    {
        var gen = (from m in Gen.Choose(1, 5000)
                   from nKb in Gen.Choose(1, 1_000_000)
                   from s1 in Gen.Choose(0, 9000)
                   from delta in Gen.Choose(1, 1000)
                   select (m, n: (long)nKb * 1000, lo: s1 / 100.0, hi: (s1 + delta) / 100.0)).ToArbitrary();

        return Prop.ForAll(gen, t =>
        {
            var loStat = ProbeDesigner.ComputeKarlinAltschul(t.lo, t.m, t.n);
            var hiStat = ProbeDesigner.ComputeKarlinAltschul(t.hi, t.m, t.n);
            return (hiStat.EValue <= loStat.EValue + 1e-12 && hiStat.BitScore >= loStat.BitScore - 1e-12)
                .Label($"score {t.lo}→{t.hi}: E {loStat.EValue}→{hiStat.EValue}, bits {loStat.BitScore}→{hiStat.BitScore}");
        });
    }

    /// <summary>
    /// INV-3 (M): a larger search space (longer database) yields a higher E-value — E is linear in the
    /// search space m·n (Karlin & Altschul 1990).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property KarlinAltschul_LargerSearchSpace_HigherEValue()
    {
        var gen = (from scoreCenti in Gen.Choose(0, 10000)
                   from m in Gen.Choose(1, 5000)
                   from n1Kb in Gen.Choose(1, 500_000)
                   from extraKb in Gen.Choose(1, 500_000)
                   select (score: scoreCenti / 100.0, m, n1: (long)n1Kb * 1000, n2: (long)(n1Kb + extraKb) * 1000)).ToArbitrary();

        return Prop.ForAll(gen, t =>
        {
            double eSmall = ProbeDesigner.ComputeKarlinAltschul(t.score, t.m, t.n1).EValue;
            double eLarge = ProbeDesigner.ComputeKarlinAltschul(t.score, t.m, t.n2).EValue;
            return (eLarge >= eSmall - 1e-12).Label($"db {t.n1}→{t.n2}: E {eSmall}→{eLarge} must not decrease");
        });
    }

}
