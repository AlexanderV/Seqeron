using FsCheck;
using FsCheck.Fluent;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for epigenetics algorithms.
///
/// Test Units: EPIGEN-CPG-001 (CpG site detection, observed/expected ratio, CpG islands).
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
}
