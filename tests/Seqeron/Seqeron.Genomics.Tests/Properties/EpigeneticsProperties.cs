using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for epigenetics algorithms.
///
/// Test Units: EPIGEN-CPG-001 (CpG site detection, observed/expected ratio, CpG islands), EPIGEN-AGE-001, EPIGEN-BISULF-001, EPIGEN-CHROM-001, EPIGEN-DMR-001, EPIGEN-METHYL-001.
/// Future siblings (EPIGEN-AGE/BISULF/CHROM/DMR/METHYL-001) extend this fixture in their own regions.
///
/// Theory: CpG dinucleotide scanning (Gardiner-Garden &amp; Frommer 1987; Wikipedia "CpG site");
/// observed/expected ratio O/E = CpG / ((C·G)/L); island criteria length≥minLength, GC≥minGc, O/E≥minCpGRatio.
/// Doc: docs/algorithms/Epigenetics/CpG_Site_Detection.md.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Epigenetics")]
public class EpigeneticsProperties
{
    #region EPIGEN-CPG-001 — Generators &amp; Independent Oracles

    /// <summary>Generates random valid DNA (A/C/G/T only) within an explicit length window.</summary>
    private static Gen<string> DnaGen(int minLen, int maxLen) =>
        from len in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
        select new string(chars);

    /// <summary>Wraps <see cref="DnaGen"/> as an arbitrary.</summary>
    private static Arbitrary<string> DnaArbitrary(int minLen, int maxLen) =>
        DnaGen(minLen, maxLen).ToArbitrary();

    /// <summary>
    /// Generates arbitrary character soup: valid bases, case variants, IUPAC codes, RNA, and junk,
    /// including the empty string. Used to stress range/edge invariants on undefined alphabets.
    /// </summary>
    private static Arbitrary<string> MixedSequenceArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T', 'a', 'c', 'g', 't',
                     'N', 'n', 'U', 'u', 'R', 'Y', '-', ' ', 'x', '7')
            .ArrayOf()
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// Generates a fixed multiset of bases (cCount C's, gCount G's, plus filler A's) and TWO
    /// arrangements of it: a "max-adjacency" string ("CGCG…" then leftovers) and a
    /// "min-adjacency" string ("CC…GG…AA…"). Both share identical C-count, G-count and length,
    /// so O/E differs only through the CpG count — the substrate for the more-CG monotonicity test.
    /// </summary>
    private static Arbitrary<(string maxAdj, string minAdj)> SameMultisetPairArbitrary() =>
        (from k in Gen.Choose(2, 30)            // ≥2 CG pairs so max- and min-adjacency strings differ
         from fillA in Gen.Choose(0, 20)        // extra A filler (changes L but kept equal in both)
         let cCount = k
         let gCount = k
         // max adjacency: k copies of "CG" gives exactly k CpG dinucleotides, then A filler.
         let maxAdj = string.Concat(Enumerable.Repeat("CG", k)) + new string('A', fillA)
         // min adjacency: all C's first, then all G's, then filler → exactly ONE "CG" junction.
         let minAdj = new string('C', cCount) + new string('G', gCount) + new string('A', fillA)
         select (maxAdj, minAdj)).ToArbitrary();

    /// <summary>
    /// Independent CpG-site oracle: 0-based indices i where upper(seq)[i]=='C' &amp;&amp; upper(seq)[i+1]=='G'.
    /// Recomputed from scratch — never calls production. INV-01 reference.
    /// </summary>
    private static List<int> ExpectedCpGSites(string seq)
    {
        var result = new List<int>();
        if (string.IsNullOrEmpty(seq))
            return result;
        string s = seq.ToUpperInvariant();
        for (int i = 0; i < s.Length - 1; i++)
            if (s[i] == 'C' && s[i + 1] == 'G')
                result.Add(i);
        return result;
    }

    /// <summary>
    /// Independent O/E oracle: cpg / ((C·G)/L) on the uppercased sequence, exactly 0.0 when
    /// null/empty/len&lt;2 or expected==0 (C==0 or G==0). Mirrors Gardiner-Garden &amp; Frommer (1987)
    /// with literal counting — never routed through production. INV-02 reference.
    /// </summary>
    private static double ExpectedObservedExpected(string seq)
    {
        if (string.IsNullOrEmpty(seq) || seq.Length < 2)
            return 0.0;
        string s = seq.ToUpperInvariant();
        int c = s.Count(ch => ch == 'C');
        int g = s.Count(ch => ch == 'G');
        int cpg = 0;
        for (int i = 0; i < s.Length - 1; i++)
            if (s[i] == 'C' && s[i + 1] == 'G')
                cpg++;
        double expected = (c * (double)g) / s.Length;
        return expected > 0 ? cpg / expected : 0.0;
    }

    /// <summary>
    /// Independent GC fraction oracle: (#G+#C) / (#valid ACGT bases) on the uppercased window.
    /// Matches the production GC-fraction definition (valid-base denominator); for pure-ACGT
    /// island windows this equals GC/length. Used to recompute island GcContent. INV-04 reference.
    /// </summary>
    private static double ExpectedGcFraction(string seq)
    {
        if (string.IsNullOrEmpty(seq))
            return 0.0;
        string s = seq.ToUpperInvariant();
        int gc = 0, valid = 0;
        foreach (char ch in s)
        {
            if (ch is 'A' or 'C' or 'G' or 'T') valid++;
            if (ch is 'G' or 'C') gc++;
        }
        return valid == 0 ? 0.0 : (double)gc / valid;
    }

    #endregion

    #region EPIGEN-CPG-001 — FindCpGSites (INV-01)

    /// <summary>
    /// INV-01 (exact): the reported CpG positions EXACTLY equal the independently recomputed set
    /// { i : upper(seq)[i]=='C' &amp;&amp; upper(seq)[i+1]=='G' } — same elements, same order. This pins the
    /// defining predicate of a CpG dinucleotide (Wikipedia "CpG site"; doc §2.4 INV-01) across hundreds
    /// of random ACGT sequences, catching off-by-one, GpC confusion, or missed adjacent CpGs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindCpGSites_MatchesIndependentOracle()
    {
        return Prop.ForAll(DnaArbitrary(0, 60), seq =>
        {
            var actual = EpigeneticsAnalyzer.FindCpGSites(seq).ToList();
            var expected = ExpectedCpGSites(seq);
            return actual.SequenceEqual(expected)
                .Label($"sites [{string.Join(",", actual)}] != oracle [{string.Join(",", expected)}] for '{seq}'");
        });
    }

    /// <summary>
    /// INV-01 (predicate): every reported position p has 'C' at p and 'G' at p+1 in the uppercased
    /// sequence; the position is in range [0, len-2]. Operates on mixed-case/junk soup so the
    /// case-insensitive normalisation and the C-then-G order are both exercised.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindCpGSites_EveryReportedSiteIsCThenG()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            seq ??= "";
            string s = seq.ToUpperInvariant();
            foreach (int p in EpigeneticsAnalyzer.FindCpGSites(seq))
            {
                if (p < 0 || p + 1 >= s.Length || s[p] != 'C' || s[p + 1] != 'G')
                    return false.Label($"bad site {p} in '{seq}'");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-D (Determinism): identical input yields an identical site list across repeated calls.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindCpGSites_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(0, 60), seq =>
        {
            var a = EpigeneticsAnalyzer.FindCpGSites(seq).ToList();
            var b = EpigeneticsAnalyzer.FindCpGSites(seq).ToList();
            return a.SequenceEqual(b).Label($"non-deterministic sites for '{seq}'");
        });
    }

    /// <summary>
    /// Anchors for the C-then-G predicate (doc §6.1): adjacent CGCG reports BOTH cytosine positions
    /// (0 and 2); GpC ("GC") is NEVER reported; case is ignored; null/empty yield no sites.
    /// </summary>
    [TestCase("CGCG", new[] { 0, 2 }, TestName = "FindCpGSites: adjacent CGCG → {0,2}")]
    [TestCase("GC", new int[0], TestName = "FindCpGSites: GpC is not a CpG → {}")]
    [TestCase("ACGT", new[] { 1 }, TestName = "FindCpGSites: single CpG at index 1")]
    [TestCase("cg", new[] { 0 }, TestName = "FindCpGSites: lowercase cg → {0}")]
    [TestCase("AAAA", new int[0], TestName = "FindCpGSites: no C/G → {}")]
    [TestCase("", new int[0], TestName = "FindCpGSites: empty → {}")]
    [Category("Property")]
    public void FindCpGSites_Anchors(string seq, int[] expected)
    {
        Assert.That(EpigeneticsAnalyzer.FindCpGSites(seq).ToArray(), Is.EqualTo(expected));
    }

    /// <summary>Edge: null input yields no sites (doc §6.1 guard).</summary>
    [Test]
    [Category("Property")]
    public void FindCpGSites_Null_YieldsNoSites()
    {
        Assert.That(EpigeneticsAnalyzer.FindCpGSites(null!).ToArray(), Is.Empty);
    }

    #endregion

    #region EPIGEN-CPG-001 — CalculateCpGObservedExpected (R: ratio ≥ 0, INV-02)

    /// <summary>
    /// INV-02 (oracle): for ANY sequence, production O/E equals the independent
    /// cpg / ((C·G)/L) oracle within 1e-9 (Gardiner-Garden &amp; Frommer 1987). This is the rigorous core:
    /// it pins the formula and the counting to a literally re-implemented reference, not the code's own.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CalculateCpGObservedExpected_MatchesOracle()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            double actual = EpigeneticsAnalyzer.CalculateCpGObservedExpected(seq);
            double expected = ExpectedObservedExpected(seq);
            return (Math.Abs(actual - expected) < 1e-9)
                .Label($"O/E mismatch for '{seq}': {actual} vs oracle {expected}");
        });
    }

    /// <summary>
    /// INV-R (ratio ≥ 0): the O/E ratio is always non-negative and finite for ANY input — the
    /// checklist range invariant. CpG count, C, G, and L are all non-negative, so the quotient cannot
    /// be negative; the guards keep it finite.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CalculateCpGObservedExpected_IsNonNegativeAndFinite()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            double oe = EpigeneticsAnalyzer.CalculateCpGObservedExpected(seq);
            return (oe >= 0.0 && double.IsFinite(oe))
                .Label($"O/E {oe} not ≥ 0 / finite for '{seq}'");
        });
    }

    /// <summary>
    /// INV-Edge (R: returns exactly 0.0): O/E is exactly 0.0 whenever the expected count is zero,
    /// i.e. the sequence has no C or no G (doc §6.1). Built over ACGT sequences that deliberately
    /// omit C (or omit G) so the expected denominator collapses.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CalculateCpGObservedExpected_NoCOrNoG_IsZero()
    {
        var arb = (from len in Gen.Choose(2, 40)
                   from dropC in Gen.Elements(true, false)
                   // Omit either C or G; remaining alphabet still has the other so it is a real test.
                   let alphabet = dropC ? new[] { 'A', 'G', 'T' } : new[] { 'A', 'C', 'T' }
                   from chars in Gen.Elements(alphabet).ArrayOf(len)
                   select new string(chars)).ToArbitrary();

        return Prop.ForAll(arb, seq =>
            (EpigeneticsAnalyzer.CalculateCpGObservedExpected(seq) == 0.0)
                .Label($"O/E not 0.0 for no-C-or-no-G '{seq}'"));
    }

    /// <summary>
    /// INV-D (Determinism): identical input yields a bit-identical ratio across repeated calls.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CalculateCpGObservedExpected_IsDeterministic()
    {
        return Prop.ForAll(MixedSequenceArbitrary(), seq =>
        {
            double a = EpigeneticsAnalyzer.CalculateCpGObservedExpected(seq);
            double b = EpigeneticsAnalyzer.CalculateCpGObservedExpected(seq);
            return a.Equals(b).Label($"non-deterministic O/E for '{seq}': {a} vs {b}");
        });
    }

    /// <summary>
    /// Hand-computed anchors (doc §2.2 formula): "CGCG" → C=2,G=2,L=4,cpg=2 ⇒ 2/((2·2)/4)=2.0;
    /// "CCGG" → cpg=1 ⇒ 1/((2·2)/4)=1.0; "CGCGCG" → C=3,G=3,L=6,cpg=3 ⇒ 3/((3·3)/6)=2.0.
    /// Absolute values cannot pass a self-consistent-but-wrong constant.
    /// </summary>
    [TestCase("CGCG", 2.0, TestName = "O/E: CGCG = 2/((2·2)/4) = 2.0")]
    [TestCase("CCGG", 1.0, TestName = "O/E: CCGG = 1/((2·2)/4) = 1.0")]
    [TestCase("CGCGCG", 2.0, TestName = "O/E: CGCGCG = 3/((3·3)/6) = 2.0")]
    [Category("Property")]
    public void CalculateCpGObservedExpected_Anchors(string seq, double expected)
    {
        Assert.That(EpigeneticsAnalyzer.CalculateCpGObservedExpected(seq),
            Is.EqualTo(expected).Within(1e-9));
    }

    /// <summary>Edge: null/empty/length-1 inputs all return exactly 0.0 (doc §6.1).</summary>
    [TestCase(null, TestName = "O/E: null → 0.0")]
    [TestCase("", TestName = "O/E: empty → 0.0")]
    [TestCase("C", TestName = "O/E: length-1 → 0.0")]
    [TestCase("G", TestName = "O/E: length-1 → 0.0")]
    [Category("Property")]
    public void CalculateCpGObservedExpected_EdgeInputs_AreZero(string? seq)
    {
        Assert.That(EpigeneticsAnalyzer.CalculateCpGObservedExpected(seq!), Is.EqualTo(0.0));
    }

    #endregion

    #region EPIGEN-CPG-001 — More CG ⇒ Higher Ratio (M, central business invariant)

    /// <summary>
    /// INV-M (more CG → higher ratio): holding C-count, G-count and length FIXED, more CpG
    /// dinucleotides yield a STRICTLY higher O/E, because O/E = cpgCount / constant when the
    /// composition is fixed. The generator produces two arrangements of the SAME base multiset —
    /// "CGCG…" (k CpG junctions, maximal) versus "CC…GG…" (a single junction, minimal) — and asserts
    /// the max-adjacency O/E strictly exceeds the min-adjacency O/E. This is the checklist's core
    /// monotonicity claim, verified on identical compositions so the comparison is confounder-free.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CalculateCpGObservedExpected_MoreCG_StrictlyHigherRatio()
    {
        return Prop.ForAll(SameMultisetPairArbitrary(), pair =>
        {
            var (maxAdj, minAdj) = pair;
            double oeMax = EpigeneticsAnalyzer.CalculateCpGObservedExpected(maxAdj);
            double oeMin = EpigeneticsAnalyzer.CalculateCpGObservedExpected(minAdj);

            // Sanity: both strings share C-count, G-count, and length (fixed composition).
            int Count(string s, char ch) => s.Count(c => c == ch);
            bool sameComposition =
                maxAdj.Length == minAdj.Length &&
                Count(maxAdj, 'C') == Count(minAdj, 'C') &&
                Count(maxAdj, 'G') == Count(minAdj, 'G');

            return (sameComposition && oeMax > oeMin)
                .Label($"more-CG monotonicity broken: max '{maxAdj}'={oeMax} vs min '{minAdj}'={oeMin}");
        });
    }

    /// <summary>
    /// Anchor for the fixed-composition monotonicity: "CGCG" (max adjacency, cpg=2 ⇒ O/E=2.0) versus
    /// "CCGG" (same C=2/G=2/L=4 multiset, cpg=1 ⇒ O/E=1.0). More CG ⇒ strictly higher ratio.
    /// </summary>
    [Test]
    [Category("Property")]
    public void CalculateCpGObservedExpected_MoreCG_Anchor_CGCG_GreaterThan_CCGG()
    {
        double cgcg = EpigeneticsAnalyzer.CalculateCpGObservedExpected("CGCG");
        double ccgg = EpigeneticsAnalyzer.CalculateCpGObservedExpected("CCGG");
        Assert.Multiple(() =>
        {
            Assert.That(cgcg, Is.EqualTo(2.0).Within(1e-9));
            Assert.That(ccgg, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(cgcg, Is.GreaterThan(ccgg),
                "more CG adjacencies at fixed composition must give a strictly higher O/E");
        });
    }

    #endregion

    #region EPIGEN-CPG-001 — FindCpGIslands (R: GC% ∈ [0,1], INV-04)

    /// <summary>
    /// INV-04 (R + recompute): for EVERY returned island, recomputing on the island substring
    /// independently must reproduce the reported metrics and satisfy all three criteria:
    /// GcContent ∈ [0,1] and == ExpectedGcFraction(window); CpGRatio == CalculateCpGObservedExpected(window)
    /// (the doc defines the island ratio via that same O/E); length End−Start ≥ minLength;
    /// GcContent ≥ minGc; CpGRatio ≥ minCpGRatio. The 'End' coordinate is treated as exclusive
    /// (doc §3.3). Exercised over GC/CpG-rich random templates so islands actually appear.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindCpGIslands_EveryIsland_SatisfiesCriteriaAndRecomputes()
    {
        // GC/CpG-rich substrate: bias toward C and G (and many CG pairs) so islands are produced.
        var richGen = (from len in Gen.Choose(200, 400)
                       from chars in Gen.Elements('C', 'G', 'C', 'G', 'A', 'T').ArrayOf(len)
                       select new string(chars)).ToArbitrary();

        return Prop.ForAll(richGen, seq =>
        {
            const int minLength = 200;
            const double minGc = 0.5;
            const double minCpGRatio = 0.6;

            foreach (var (start, end, gc, ratio) in
                     EpigeneticsAnalyzer.FindCpGIslands(seq, minLength, minGc, minCpGRatio))
            {
                if (start < 0 || end > seq.Length || end - start < minLength)
                    return false.Label($"bad span [{start},{end}) in len {seq.Length}");

                string window = seq.Substring(start, end - start);
                double oracleGc = ExpectedGcFraction(window);
                double oracleRatio = EpigeneticsAnalyzer.CalculateCpGObservedExpected(window);

                if (gc < 0.0 || gc > 1.0)
                    return false.Label($"GcContent {gc} outside [0,1]");
                if (Math.Abs(gc - oracleGc) > 1e-9)
                    return false.Label($"GcContent {gc} != oracle {oracleGc}");
                if (Math.Abs(ratio - oracleRatio) > 1e-9)
                    return false.Label($"CpGRatio {ratio} != recomputed {oracleRatio}");
                if (gc < minGc || ratio < minCpGRatio)
                    return false.Label($"criteria violated: gc={gc}, ratio={ratio}");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// INV-D (Determinism): identical input yields an identical island list across repeated calls.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FindCpGIslands_IsDeterministic()
    {
        var richGen = (from len in Gen.Choose(200, 360)
                       from chars in Gen.Elements('C', 'G', 'A', 'T').ArrayOf(len)
                       select new string(chars)).ToArbitrary();

        return Prop.ForAll(richGen, seq =>
        {
            var a = EpigeneticsAnalyzer.FindCpGIslands(seq).ToList();
            var b = EpigeneticsAnalyzer.FindCpGIslands(seq).ToList();
            return a.SequenceEqual(b).Label($"non-deterministic islands for length {seq.Length}");
        });
    }

    /// <summary>
    /// Positive anchor: "CG" repeated 120× = 240 bp has GC=1.0 and a very high O/E, so it MUST yield
    /// at least one island, and the reported island must cover ≥200 bp with GcContent==1.0 and
    /// CpGRatio recomputed on the window. Source: doc §2.2 criteria; §7 examples.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindCpGIslands_GcCpGRichSequence_YieldsIsland()
    {
        string seq = string.Concat(Enumerable.Repeat("CG", 120)); // 240 bp, GC=1.0
        var islands = EpigeneticsAnalyzer.FindCpGIslands(seq).ToList();

        Assert.That(islands, Is.Not.Empty, "a 240 bp GC+CpG-rich sequence must contain an island");
        foreach (var (start, end, gc, ratio) in islands)
        {
            string window = seq.Substring(start, end - start);
            Assert.Multiple(() =>
            {
                Assert.That(end - start, Is.GreaterThanOrEqualTo(200), "island ≥ minLength");
                Assert.That(gc, Is.EqualTo(1.0).Within(1e-9), "all-CG window has GC=1.0");
                Assert.That(ratio, Is.EqualTo(EpigeneticsAnalyzer.CalculateCpGObservedExpected(window)).Within(1e-9),
                    "reported CpGRatio must equal the O/E recomputed on the window");
            });
        }
    }

    /// <summary>
    /// Negative anchor: a 240 bp AT-only sequence has GC=0.0 and no CpG dinucleotides, so it fails
    /// every island criterion and MUST yield no islands. Source: doc §2.2 criteria.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindCpGIslands_AtOnlySequence_YieldsNoIslands()
    {
        string seq = new string('A', 120) + new string('T', 120); // 240 bp, GC=0.0
        Assert.That(EpigeneticsAnalyzer.FindCpGIslands(seq).ToList(), Is.Empty);
    }

    /// <summary>
    /// Edge: a sequence shorter than minLength yields no islands; null/empty yield no islands
    /// (doc §6.1). 199 bp of all-CG (otherwise island-qualifying) is rejected purely on length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindCpGIslands_ShorterThanMinLength_OrNullEmpty_YieldNoIslands()
    {
        string shortRich = string.Concat(Enumerable.Repeat("CG", 99)) + "C"; // 199 bp
        Assert.Multiple(() =>
        {
            Assert.That(EpigeneticsAnalyzer.FindCpGIslands(shortRich).ToList(), Is.Empty, "len 199 < 200");
            Assert.That(EpigeneticsAnalyzer.FindCpGIslands(null!).ToList(), Is.Empty, "null");
            Assert.That(EpigeneticsAnalyzer.FindCpGIslands("").ToList(), Is.Empty, "empty");
        });
    }

    #endregion

    #region EPIGEN-AGE-001: M: more methylation at (positively-weighted) clock sites → higher age; D: deterministic

    // CalculateEpigeneticAge = Horvath (2013) linear predictor (intercept + Σ coef·β) mapped through
    // the monotone anti-transform F⁻¹. NOTE: the anti-transform infimum is −1 (as x→−∞), so the
    // theoretical lower bound is age > −1, not the checklist's ≥ 0 — we test the true bound.

    /// <summary>Clock with positive coefficients and two profiles where methHigh ≥ methLow everywhere.</summary>
    private static Arbitrary<(Dictionary<string, double> coeffs, Dictionary<string, double> low, Dictionary<string, double> high)>
        ClockArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int n = 1 + rng.Next(5);
            var coeffs = new Dictionary<string, double>();
            var low = new Dictionary<string, double>();
            var high = new Dictionary<string, double>();
            for (int i = 0; i < n; i++)
            {
                string cg = $"cg{i}";
                coeffs[cg] = 0.1 + rng.NextDouble();           // strictly positive weight
                double a = rng.NextDouble();                    // β in [0,1]
                low[cg] = a;
                high[cg] = a + (1.0 - a) * rng.NextDouble();    // ≥ a, ≤ 1
            }
            return (coeffs, low, high);
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (M): with positive clock coefficients, raising methylation at the clock CpGs does not
    /// lower the predicted age.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EpigeneticAge_MoreMethylation_RaisesAge()
    {
        return Prop.ForAll(ClockArbitrary(), c =>
        {
            double ageLow = EpigeneticsAnalyzer.CalculateEpigeneticAge(c.low, c.coeffs);
            double ageHigh = EpigeneticsAnalyzer.CalculateEpigeneticAge(c.high, c.coeffs);
            return (ageHigh >= ageLow - 1e-9).Label($"age dropped with more methylation: {ageHigh} < {ageLow}");
        });
    }

    /// <summary>
    /// INV-2 (R, true bound): the predicted age is always greater than −1 (the anti-transform infimum).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EpigeneticAge_IsAboveAntiTransformInfimum()
    {
        return Prop.ForAll(ClockArbitrary(), c =>
        {
            double age = EpigeneticsAnalyzer.CalculateEpigeneticAge(c.low, c.coeffs);
            return (age > -1.0).Label($"age {age} ≤ −1");
        });
    }

    /// <summary>
    /// INV-3 (monotone calibration): the Horvath anti-transform is non-decreasing and equals 20 at 0.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AntiTransform_IsMonotone()
    {
        var pairs = Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            double x1 = (rng.NextDouble() - 0.5) * 20;
            double x2 = x1 + rng.NextDouble() * 10;
            return (x1, x2);
        }).ToArbitrary();

        return Prop.ForAll(pairs, p =>
            (EpigeneticsAnalyzer.HorvathAntiTransform(p.x1) <= EpigeneticsAnalyzer.HorvathAntiTransform(p.x2) + 1e-9)
                .Label("anti-transform must be non-decreasing"));
    }

    /// <summary>
    /// INV-4 (D + boundary): age is deterministic; an empty coefficient table is rejected.
    /// </summary>
    [Test]
    [Category("Property")]
    public void EpigeneticAge_DeterministicAndBoundary()
    {
        var beta = new Dictionary<string, double> { ["cg0"] = 0.5 };
        var coeffs = new Dictionary<string, double> { ["cg0"] = 1.0 };
        double age1 = EpigeneticsAnalyzer.CalculateEpigeneticAge(beta, coeffs);
        double age2 = EpigeneticsAnalyzer.CalculateEpigeneticAge(beta, coeffs);
        Assert.Multiple(() =>
        {
            Assert.That(age2, Is.EqualTo(age1), "deterministic");
            Assert.That(EpigeneticsAnalyzer.HorvathAntiTransform(0.0), Is.EqualTo(20.0).Within(1e-9));
            Assert.Throws<ArgumentException>(
                () => EpigeneticsAnalyzer.CalculateEpigeneticAge(beta, new Dictionary<string, double>()));
        });
    }

    #endregion

    #region EPIGEN-BISULF-001: P: unmethylated C→T, methylated C preserved; P: length preserved; D: deterministic

    // SimulateBisulfiteConversion (Frommer et al. 1992): unmethylated cytosine → thymine; methylated
    // (protected) cytosine and all non-cytosine bases are unchanged; the strand length is preserved.

    /// <summary>A DNA sequence and a random subset of indices marked as protected (methylated).</summary>
    private static Arbitrary<(string seq, HashSet<int> methylated)> BisulfiteInputArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            const string bases = "ACGT";
            int len = 10 + rng.Next(20);
            var c = new char[len];
            for (int i = 0; i < len; i++) c[i] = bases[rng.Next(4)];
            var methylated = new HashSet<int>();
            for (int i = 0; i < len; i++) if (rng.Next(2) == 0) methylated.Add(i);
            return (new string(c), methylated);
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (P): each position follows the conversion rule — unmethylated C→T, methylated C kept,
    /// non-cytosine unchanged — and the converted strand has the same length.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Bisulfite_AppliesConversionRule()
    {
        return Prop.ForAll(BisulfiteInputArbitrary(), input =>
        {
            var (seq, methylated) = input;
            string conv = EpigeneticsAnalyzer.SimulateBisulfiteConversion(seq, methylated);
            if (conv.Length != seq.Length) return false.Label("length changed");
            for (int i = 0; i < seq.Length; i++)
            {
                char s = seq[i], r = conv[i];
                bool isC = s is 'C' or 'c';
                char expected = isC ? (methylated.Contains(i) ? s : (s == 'C' ? 'T' : 't')) : s;
                if (r != expected) return false.Label($"pos {i}: '{s}'→'{r}', expected '{expected}'");
            }
            return true.Label("conversion rule holds");
        });
    }

    /// <summary>
    /// INV-2 (P): with no protected positions, every cytosine is converted (none remain) and all
    /// other bases are unchanged.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Bisulfite_NoMethylation_ConvertsAllCytosines()
    {
        return Prop.ForAll(BisulfiteInputArbitrary(), input =>
        {
            var (seq, _) = input;
            string conv = EpigeneticsAnalyzer.SimulateBisulfiteConversion(seq);
            return (!conv.Contains('C') && !conv.Contains('c'))
                .Label("an unmethylated cytosine survived conversion");
        });
    }

    /// <summary>
    /// INV-3 (involution-like): re-converting the output with the same protected set is a no-op (the
    /// only remaining cytosines are protected).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Bisulfite_Reconversion_IsStable()
    {
        return Prop.ForAll(BisulfiteInputArbitrary(), input =>
        {
            var (seq, methylated) = input;
            string once = EpigeneticsAnalyzer.SimulateBisulfiteConversion(seq, methylated);
            string twice = EpigeneticsAnalyzer.SimulateBisulfiteConversion(once, methylated);
            return (once == twice).Label("re-conversion changed the strand");
        });
    }

    /// <summary>
    /// INV-4 (D + boundary): conversion is deterministic; empty input yields empty output.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Bisulfite_DeterministicAndBoundary()
    {
        var methylated = new HashSet<int> { 2 };
        string a = EpigeneticsAnalyzer.SimulateBisulfiteConversion("ACCGT", methylated);
        string b = EpigeneticsAnalyzer.SimulateBisulfiteConversion("ACCGT", methylated);
        Assert.Multiple(() =>
        {
            Assert.That(b, Is.EqualTo(a), "deterministic");
            Assert.That(a, Is.EqualTo("ATCGT"), "pos1 C→T, pos2 C protected, others unchanged");
            Assert.That(EpigeneticsAnalyzer.SimulateBisulfiteConversion(""), Is.EqualTo(""));
        });
    }

    #endregion

    #region EPIGEN-CHROM-001: P: each region assigned a state; R: positions preserved; D: deterministic

    // AnnotateHistoneModifications labels each region by its single histone mark's canonical Roadmap
    // chromatin state, or LowSignal when the mark's signal is below the presence threshold.

    private const double ChromThreshold = 1.0;

    private static EpigeneticsAnalyzer.ChromatinState ExpectedState(string mark, double signal) =>
        signal < ChromThreshold ? EpigeneticsAnalyzer.ChromatinState.LowSignal : mark.ToUpperInvariant() switch
        {
            "H3K4ME3" => EpigeneticsAnalyzer.ChromatinState.ActivePromoter,
            "H3K4ME1" => EpigeneticsAnalyzer.ChromatinState.WeakEnhancer,
            "H3K27AC" => EpigeneticsAnalyzer.ChromatinState.ActiveEnhancer,
            "H3K36ME3" => EpigeneticsAnalyzer.ChromatinState.Transcribed,
            "H3K27ME3" => EpigeneticsAnalyzer.ChromatinState.Repressed,
            "H3K9ME3" => EpigeneticsAnalyzer.ChromatinState.Heterochromatin,
            "H3K9AC" => EpigeneticsAnalyzer.ChromatinState.ActivePromoter,
            _ => EpigeneticsAnalyzer.ChromatinState.LowSignal,
        };

    private static readonly string[] ChromMarks =
        { "H3K4me3", "H3K4me1", "H3K27ac", "H3K36me3", "H3K27me3", "H3K9me3", "H3K9ac", "H3Kunknown" };

    private static Arbitrary<(int Start, int End, string Mark, double Signal)[]> HistoneRegionsArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int n = 1 + rng.Next(6);
            var arr = new (int, int, string, double)[n];
            for (int i = 0; i < n; i++)
                arr[i] = (i * 100, i * 100 + 99, ChromMarks[rng.Next(ChromMarks.Length)], rng.NextDouble() * 2.0);
            return arr;
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (P + R): one annotation per region, preserving Start/End/Mark/Signal, with the predicted
    /// state matching the present-mark rule (independently computed).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Chromatin_EachRegion_AssignedCorrectState()
    {
        return Prop.ForAll(HistoneRegionsArbitrary(), regions =>
        {
            var result = EpigeneticsAnalyzer.AnnotateHistoneModifications(regions, ChromThreshold).ToList();
            if (result.Count != regions.Length) return false.Label("region count changed");
            for (int i = 0; i < regions.Length; i++)
            {
                var r = result[i];
                var (s, e, m, sig) = regions[i];
                if (r.Start != s || r.End != e || r.Mark != m || r.Signal != sig)
                    return false.Label($"region {i} fields not preserved");
                if (r.PredictedState != ExpectedState(m, sig))
                    return false.Label($"region {i} state {r.PredictedState} ≠ expected {ExpectedState(m, sig)}");
            }
            return true.Label("all regions correctly annotated");
        });
    }

    /// <summary>
    /// INV-2 (P, golden): canonical marks above threshold map to their Roadmap states; below-threshold
    /// or unknown marks map to LowSignal.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Chromatin_KnownMarks_MapToStates()
    {
        var regions = new (int, int, string, double)[]
        {
            (0, 99, "H3K4me3", 5.0),
            (100, 199, "H3K27ac", 5.0),
            (200, 299, "H3K9me3", 5.0),
            (300, 399, "H3K4me3", 0.1),   // below threshold
            (400, 499, "H3Kbogus", 5.0),  // unknown mark
        };
        var result = EpigeneticsAnalyzer.AnnotateHistoneModifications(regions, ChromThreshold).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(result[0].PredictedState, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.ActivePromoter));
            Assert.That(result[1].PredictedState, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.ActiveEnhancer));
            Assert.That(result[2].PredictedState, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.Heterochromatin));
            Assert.That(result[3].PredictedState, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.LowSignal));
            Assert.That(result[4].PredictedState, Is.EqualTo(EpigeneticsAnalyzer.ChromatinState.LowSignal));
        });
    }

    /// <summary>
    /// INV-3 (D): Annotation is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Chromatin_IsDeterministic()
    {
        return Prop.ForAll(HistoneRegionsArbitrary(), regions =>
        {
            var a = EpigeneticsAnalyzer.AnnotateHistoneModifications(regions, ChromThreshold).Select(r => r.PredictedState).ToList();
            var b = EpigeneticsAnalyzer.AnnotateHistoneModifications(regions, ChromThreshold).Select(r => r.PredictedState).ToList();
            return a.SequenceEqual(b).Label("AnnotateHistoneModifications must be deterministic");
        });
    }

    #endregion

    #region EPIGEN-DMR-001: R: start < end; M: lower threshold → ≥ DMRs; P: |Δmethylation| > threshold; D: deterministic

    // FindDMRs (methylKit tiling model): positions are grouped into fixed windows; a window is a DMR
    // when |mean(sample2 − sample1)| exceeds the cutoff and it has ≥ minCpGCount covered cytosines.

    private static EpigeneticsAnalyzer.MethylationSite Site(int pos, double level) =>
        new(pos, EpigeneticsAnalyzer.MethylationType.CpG, "CG", level, 10);

    /// <summary>Two single-window samples of 3..8 CpGs with random per-site methylation levels.</summary>
    private static Arbitrary<(EpigeneticsAnalyzer.MethylationSite[] s1, EpigeneticsAnalyzer.MethylationSite[] s2)>
        DmrSamplesArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int n = 3 + rng.Next(6);
            var s1 = new EpigeneticsAnalyzer.MethylationSite[n];
            var s2 = new EpigeneticsAnalyzer.MethylationSite[n];
            for (int i = 0; i < n; i++)
            {
                s1[i] = Site(i, rng.NextDouble());
                s2[i] = Site(i, rng.NextDouble());
            }
            return (s1, s2);
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R + P): every DMR has Start &lt; End, at least minCpGCount cytosines, an absolute mean
    /// difference exceeding the cutoff, and a hyper/hypo annotation consistent with its sign.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Dmr_RegionsAreValid()
    {
        return Prop.ForAll(DmrSamplesArbitrary(), input =>
        {
            var (s1, s2) = input;
            const double minDiff = 0.25;
            var dmrs = EpigeneticsAnalyzer.FindDMRs(s1, s2, windowSize: 1000, minDifference: minDiff, minCpGCount: 3).ToList();
            bool ok = dmrs.All(d =>
                d.Start < d.End && d.CpGCount >= 3 &&
                Math.Abs(d.MeanDifference) > minDiff &&
                d.Annotation == (d.MeanDifference > 0 ? "Hypermethylated" : "Hypomethylated"));
            return ok.Label("a DMR was invalid");
        });
    }

    /// <summary>
    /// INV-2 (M): a lower difference cutoff reports at least as many DMRs — the strict-cutoff region
    /// starts are a subset of the loose-cutoff region starts.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Dmr_LowerThreshold_IsSuperset()
    {
        return Prop.ForAll(DmrSamplesArbitrary(), input =>
        {
            var (s1, s2) = input;
            var loose = EpigeneticsAnalyzer.FindDMRs(s1, s2, 1000, 0.1, 3).Select(d => d.Start).ToHashSet();
            var strict = EpigeneticsAnalyzer.FindDMRs(s1, s2, 1000, 0.5, 3).Select(d => d.Start).ToHashSet();
            return (strict.IsSubsetOf(loose) && strict.Count <= loose.Count)
                .Label($"strict DMRs ({strict.Count}) not ⊆ loose DMRs ({loose.Count})");
        });
    }

    /// <summary>
    /// INV-3 (D + positive control): DMR calling is deterministic; a clear hypo→hyper shift is called
    /// as one Hypermethylated DMR.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Dmr_DeterministicAndGolden()
    {
        var control = new[] { Site(0, 0.1), Site(1, 0.1), Site(2, 0.1), Site(3, 0.1) };
        var treatment = new[] { Site(0, 0.9), Site(1, 0.9), Site(2, 0.9), Site(3, 0.9) };
        var a = EpigeneticsAnalyzer.FindDMRs(control, treatment, 1000, 0.25, 3).ToList();
        var b = EpigeneticsAnalyzer.FindDMRs(control, treatment, 1000, 0.25, 3).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(a.Select(d => d.Start), Is.EqualTo(b.Select(d => d.Start)), "deterministic");
            Assert.That(a, Has.Count.EqualTo(1));
            Assert.That(a[0].Annotation, Is.EqualTo("Hypermethylated"));
            Assert.That(a[0].MeanDifference, Is.EqualTo(0.8).Within(1e-9));
        });
    }

    #endregion

    #region EPIGEN-METHYL-001: R: level ∈ [0,1]; P: level = methylated/total; D: deterministic

    // CalculateMethylationFromBisulfite (Bismark call rule): at each reference CpG, a read base C is a
    // methylated call and T an unmethylated call; level = C / (C+T), Coverage = valid C/T calls.

    /// <summary>A reference and a few reads aligned at random start positions.</summary>
    private static Arbitrary<(string reference, (string, int)[] reads)> MethylInputArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            const string bases = "ACGT";
            int refLen = 20;
            var rc = new char[refLen];
            for (int i = 0; i < refLen; i++) rc[i] = bases[rng.Next(4)];
            string reference = new string(rc);

            int n = 1 + rng.Next(4);
            var reads = new (string, int)[n];
            for (int r = 0; r < n; r++)
            {
                int rl = 5 + rng.Next(6);
                var c = new char[rl];
                for (int i = 0; i < rl; i++) c[i] = bases[rng.Next(4)];
                reads[r] = (new string(c), rng.Next(refLen));
            }
            return (reference, reads);
        }).ToArbitrary();

    /// <summary>Independent (Bismark-rule) recomputation of (methylated, total) at a CpG position.</summary>
    private static (int meth, int total) RecomputeCall(string reference, (string seq, int start)[] reads, int site)
    {
        int meth = 0, total = 0;
        foreach (var (seq, start) in reads)
        {
            string read = seq.ToUpperInvariant();
            for (int i = 0; i < read.Length && start + i < reference.Length - 1; i++)
            {
                if (start + i != site) continue;
                if (read[i] == 'C') { meth++; total++; }
                else if (read[i] == 'T') { total++; }
            }
        }
        return (meth, total);
    }

    /// <summary>
    /// INV-1 (R + P): every returned site is a reference CpG with level ∈ [0,1] equal to its
    /// methylated/total Bismark call fraction and Coverage equal to the valid C/T call count.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Methylation_LevelEqualsMethylatedOverTotal()
    {
        return Prop.ForAll(MethylInputArbitrary(), input =>
        {
            var (reference, reads) = input;
            var sites = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(reference, reads).ToList();
            bool ok = sites.All(s =>
            {
                bool isCpg = s.Position + 1 < reference.Length && reference[s.Position] == 'C' && reference[s.Position + 1] == 'G';
                var (meth, total) = RecomputeCall(reference, reads, s.Position);
                return isCpg && s.Coverage > 0 && total == s.Coverage &&
                       s.MethylationLevel is >= 0.0 and <= 1.0 &&
                       Math.Abs(s.MethylationLevel - (double)meth / total) < 1e-9 &&
                       s.Type == EpigeneticsAnalyzer.MethylationType.CpG;
            });
            return ok.Label("a methylation call did not equal methylated/total at a CpG");
        });
    }

    /// <summary>
    /// INV-2 (D + golden): all-C reads give level 1, all-T reads give 0, a half/half mix gives 0.5;
    /// calling is deterministic.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Methylation_DeterministicAndGolden()
    {
        const string reference = "TCGAA"; // CpG at index 1
        var allC = new[] { ("TCGAA", 0), ("TCGAA", 0) };  // 'C' at index 1
        var allT = new[] { ("TTGAA", 0), ("TTGAA", 0) };  // 'T' at index 1
        var mix = new[] { ("TCGAA", 0), ("TTGAA", 0) };

        Assert.Multiple(() =>
        {
            Assert.That(EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(reference, allC).Single().MethylationLevel,
                Is.EqualTo(1.0).Within(1e-9), "all methylated");
            Assert.That(EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(reference, allT).Single().MethylationLevel,
                Is.EqualTo(0.0).Within(1e-9), "all unmethylated");
            Assert.That(EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(reference, mix).Single().MethylationLevel,
                Is.EqualTo(0.5).Within(1e-9), "half methylated");
            var once = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(reference, mix).Select(s => s.MethylationLevel).ToList();
            var twice = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(reference, mix).Select(s => s.MethylationLevel).ToList();
            Assert.That(twice, Is.EqualTo(once), "deterministic");
        });
    }

    #endregion
}
