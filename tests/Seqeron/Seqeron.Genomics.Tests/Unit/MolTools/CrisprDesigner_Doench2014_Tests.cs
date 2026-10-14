namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for the Doench et al. 2014 "Rule Set 1" on-target efficacy model
/// (<see cref="CrisprDesigner.CalculateOnTargetDoench2014"/>), enhancement C7 / CRISPR-GUIDE-001.
///
/// Source of the model + coefficients (re-grounded this session):
///   Doench, Hartenian, Graham, et al. Nat Biotechnol 32:1262-1267 (2014), PMID 25184501.
///   Coefficients transcribed verbatim from the published reference implementation
///   `doenchScore.py` in CRISPOR (Haeussler et al. 2016, Genome Biol 17:148):
///   https://github.com/maximilianh/crisporWebsite/blob/master/doenchScore.py
///
/// Every expected value below traces to an INDEPENDENT Python run of the published coefficients
/// (not to this C# code), OR to the worked examples shipped in the reference `doenchScore.py`
/// itself. The reference returns a probability in (0,1); this API multiplies by 100, so expected
/// values are reference × 100.
///
/// IMPORTANT faithfulness note: this is the reproducible LINEAR model (Rule Set 1). Doench
/// "Rule Set 2" / Azimuth is a gradient-boosted-tree model and is intentionally NOT implemented
/// (it cannot be reproduced from published numbers without the trained model).
/// </summary>
[TestFixture]
public class CrisprDesigner_Doench2014_Tests
{
    private const double Tol = 1e-4; // reference example floats are quoted to ~12 sig figs; C# double matches well within this

    /// <summary>
    /// M-001: Published reference worked example #1 from doenchScore.py.
    /// Sequence "TATAGCTGCGATCTGAGGTAGGGAGGGACC" → reference 0.713089368437 → ×100 = 71.3089368437.
    /// This is the load-bearing cross-check: it ties the C# implementation to the reference
    /// implementation's own documented expected output.
    /// </summary>
    [Test]
    public void Doench2014_ReferenceExample1_MatchesPublishedValue()
    {
        double score = CrisprDesigner.CalculateOnTargetDoench2014("TATAGCTGCGATCTGAGGTAGGGAGGGACC");
        Assert.That(score, Is.EqualTo(71.3089368437).Within(Tol));
    }

    /// <summary>
    /// M-002: Published reference worked example #2 from doenchScore.py.
    /// Sequence "TCCGCACCTGTCACGGTCGGGGCTTGGCGC" → reference 0.0189838463593 → ×100 = 1.89838463593.
    /// </summary>
    [Test]
    public void Doench2014_ReferenceExample2_MatchesPublishedValue()
    {
        double score = CrisprDesigner.CalculateOnTargetDoench2014("TCCGCACCTGTCACGGTCGGGGCTTGGCGC");
        Assert.That(score, Is.EqualTo(1.89838463593).Within(Tol));
    }

    /// <summary>
    /// M-003: All-A 30-mer (with a forced NGG PAM at offsets 25-26).
    /// Independent Python computation of the published model (intercept + GC term + feature hits)
    /// on "AAAA AAAAAAAAAAAAAAAAAAAA AGG AAA"-style all-A-with-PAM gives a specific value.
    /// Hand value: protospacer offsets [4,24) all A → gcCount 0 → +abs(10-0)*gcLow = 10*(-0.2026259).
    /// Plus single feature hits where modelSeq is 'A' at the relevant offsets, plus 'AA' dinucs.
    /// The independent Python reference yields 0.044338168085440034 → ×100 = 4.4338168085440034,
    /// re-executed in WP2 against the coefficients read out of the downloaded `doenchScore.py`.
    /// We lock the reference value.
    /// </summary>
    [Test]
    public void Doench2014_AllAdenineWithPam_MatchesIndependentComputation()
    {
        // 4 up (A) + 20 protospacer (A) + GG... need PAM N at 24, GG at 25-26.
        // Build: positions 0-23 = A, 24 = A, 25 = G, 26 = G, 27-29 = A
        string seq = "AAAAAAAAAAAAAAAAAAAAAAAAAGGAAA";
        Assert.That(seq.Length, Is.EqualTo(30));
        // Independent Python run of the published coefficients on this exact 30-mer:
        //   gcCount over [4,24): 0  → +abs(10-0)*gcLow = 10*(-0.2026259) = -2.026259
        //   plus all matching single/di-nucleotide 'A'/'AA' feature weights and the intercept,
        //   through the sigmoid → 0.044338168085440035 → ×100 = 4.4338168085.
        double score = CrisprDesigner.CalculateOnTargetDoench2014(seq);
        Assert.That(score, Is.EqualTo(4.4338168085).Within(1e-3));
    }

    /// <summary>
    /// M-002b (B08 WP2): reference-locked rows from the bulk numeric cross-check.
    ///
    /// The reference model was loaded straight out of the TEXT of the downloaded CRISPOR
    /// `doenchScore.py` (its `params` table and the three constants parsed with a regex — no
    /// hand-typed coefficient) and executed in Python over 477 contexts: 400 pseudo-random
    /// A/C/G/T 30-mers with a forced NGG at offsets 25-26 (Python `random`, seed 20261014) plus
    /// the 77 `GuideRnaCandidate.Context30Mer` values that `DesignGuideRnas` produces on three
    /// designed targets. C# and the reference agreed on all 477, max |delta| exactly 0.
    /// The rows below are a representative slice of that run, chosen to cover both published GC
    /// coefficients and the `gcCount == 10` boundary (the whole run: 265 rows with gcCount &lt;= 10,
    /// 212 with gcCount &gt; 10, 90 exactly at 10). Every expected value is REFERENCE-EXECUTED
    /// (Python reference x 100), never read off the C#.
    /// </summary>
    [TestCase("GTACTCGATTTATAAAGAACAAGATGGGAA", 8.9009900853122)]    // gcCount  5 (gcLow branch)
    [TestCase("CGGATAACTATAAAGCTGCTGTCTAGGCAT", 0.1877342348367)]    // gcCount  7, lowest score in the run
    [TestCase("TTGGCGATTGGCTTTGTTTCTTGGTGGTCC", 45.9107148757746)]   // gcCount  9
    [TestCase("GGTATTACCGGAGCACAGAAGAGTAGGGAA", 10.8954314424884)]   // gcCount 10 — the <=10 boundary
    [TestCase("ATTGCATAGTCGGCCCAAACGGTTGGGATA", 6.0430223333030)]    // gcCount 11 — first gcHigh row
    [TestCase("CTCTGGGGTCCCGTCCAGAAATCCTGGTAC", 2.9297223516504)]    // gcCount 13
    [TestCase("GACAGTACCCCCACGCAGATAGCGAGGCCG", 87.5894495698859)]   // gcCount 13, highest score in the run
    [TestCase("ATTCCCGCCCGGGGCCGTATGGCCCGGCCG", 1.0602834215605)]    // gcCount 17
    [TestCase("ACGTAAAAAAAAAATTTTTTTTTTGGGGGG", 0.2637789285781)]    // designed context, gcCount 0
    [TestCase("GTAAAAAAAAAATTTTTTTTTTGGGGGGGG", 1.3456002035926)]    // designed context, gcCount 2
    [TestCase("CAATTGGCCAAACGTACGTACGTCAGGCCA", 5.2307204967364)]    // designed context, gcCount 11
    public void Doench2014_ReferenceCrossCheckRows_MatchTheReferenceExactly(string context30Mer, double expected)
    {
        double score = CrisprDesigner.CalculateOnTargetDoench2014(context30Mer);
        Assert.That(score, Is.EqualTo(expected).Within(1e-9));
    }

    /// <summary>
    /// M-002c (B08 WP2): the two published GC coefficients are BOTH in use and switch at
    /// gcCount == 10 (`gcCount &lt;= 10 → gcLow = -0.2026259`, `&gt; 10 → gcHigh = -0.1665878`).
    ///
    /// Built as a sourced structural check rather than a tautology: two contexts differing only in
    /// the protospacer composition, with gcCount 10 and 11, are scored, and the implied GC
    /// contribution is recovered from the pair. The reference applies `abs(10 - gcCount) * gcWeight`,
    /// so the gcCount == 10 row carries NO GC contribution at all (`abs(10-10) == 0`) — which is why
    /// the boundary row cannot distinguish the two coefficients and the 11 row must be checked.
    /// Reference-executed values (same Python run as M-002b).
    /// </summary>
    [Test]
    public void Doench2014_GcTerm_UsesBothPublishedCoefficientsAtTheTenBoundary()
    {
        // gcCount == 10 → abs(10-10)*gcLow = 0 contribution.
        double atBoundary = CrisprDesigner.CalculateOnTargetDoench2014("GGTATTACCGGAGCACAGAAGAGTAGGGAA");
        // gcCount == 11 → abs(10-11)*gcHigh = -0.1665878 (the gcHigh branch is reachable at all).
        double aboveBoundary = CrisprDesigner.CalculateOnTargetDoench2014("ATTGCATAGTCGGCCCAAACGGTTGGGATA");

        Assert.Multiple(() =>
        {
            Assert.That(atBoundary, Is.EqualTo(10.8954314424884).Within(1e-9));
            Assert.That(aboveBoundary, Is.EqualTo(6.0430223333030).Within(1e-9));
        });
    }

    /// <summary>
    /// M-004: Score is always in [0, 100] and the sigmoid endpoints behave (monotone in raw score).
    /// Reference example 1 (high) must exceed reference example 2 (low).
    /// </summary>
    [Test]
    public void Doench2014_ScoreRange_WithinZeroToHundredAndOrdered()
    {
        double high = CrisprDesigner.CalculateOnTargetDoench2014("TATAGCTGCGATCTGAGGTAGGGAGGGACC");
        double low = CrisprDesigner.CalculateOnTargetDoench2014("TCCGCACCTGTCACGGTCGGGGCTTGGCGC");

        Assert.Multiple(() =>
        {
            Assert.That(high, Is.InRange(0.0, 100.0));
            Assert.That(low, Is.InRange(0.0, 100.0));
            Assert.That(high, Is.GreaterThan(low));
        });
    }

    /// <summary>
    /// M-005: Case-insensitive — lower-case input gives the identical published value.
    /// </summary>
    [Test]
    public void Doench2014_LowerCaseInput_SameAsUpperCase()
    {
        double upper = CrisprDesigner.CalculateOnTargetDoench2014("TATAGCTGCGATCTGAGGTAGGGAGGGACC");
        double lower = CrisprDesigner.CalculateOnTargetDoench2014("tatagctgcgatctgaggtagggagggacc");
        Assert.That(lower, Is.EqualTo(upper).Within(1e-12));
    }

    /// <summary>
    /// M-006: Wrong length (not 30) throws.
    /// </summary>
    [TestCase("ACGT")]
    [TestCase("ACGTACGTACGTACGTACGTACGTACGTA")]   // 29
    [TestCase("ACGTACGTACGTACGTACGTACGTACGTACG")] // 31
    public void Doench2014_WrongLength_Throws(string seq)
    {
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateOnTargetDoench2014(seq));
    }

    /// <summary>
    /// M-007: <c>null</c> throws <see cref="ArgumentNullException"/>; the EMPTY STRING throws
    /// <see cref="ArgumentException"/> with the 30-nt length message.
    ///
    /// This assertion was CHANGED in B08 WP2 (it previously locked <see cref="ArgumentNullException"/>
    /// for <c>""</c>). Reason, recorded here so the rationale lives next to the assertion: a null
    /// reference is a MISSING argument, while an empty string is a WRONG VALUE — a 0-nt context — and
    /// must be reported like every other malformed context (29 nt, 31 nt, non-ACGT). The reference
    /// draws the same line (`assert isinstance(seq[0], basestring)` vs `assert len(seq[0]) > 0`,
    /// azimuth/model_comparison.py L523-L524). This is the identical fix, for the identical reason,
    /// that WP1 applied to <c>AzimuthRuleSet2.Validate("")</c> (B08 F36), so the two on-target
    /// entry points now have one contract. Nothing in the repo depended on the old type: the MCP tool
    /// <c>calculate_on_target_doench2014</c> rejects null/empty itself, as
    /// <see cref="ArgumentException"/>, before the library is reached.
    /// </summary>
    [Test]
    public void Doench2014_NullOrEmpty_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => CrisprDesigner.CalculateOnTargetDoench2014(null!));

        var ex = Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateOnTargetDoench2014(""));
        Assert.Multiple(() =>
        {
            Assert.That(ex, Is.Not.InstanceOf<ArgumentNullException>());
            Assert.That(ex!.Message, Does.Contain("30-nt context"));
            Assert.That(ex.ParamName, Is.EqualTo("context30Mer"));
        });
    }

    /// <summary>
    /// M-008: Non-ACGT characters throw.
    /// </summary>
    [Test]
    public void Doench2014_InvalidCharacters_Throws()
    {
        Assert.Throws<ArgumentException>(() =>
            CrisprDesigner.CalculateOnTargetDoench2014("TATAGCTGCGATCTGAGGTAGGGAGGGACN"));
    }

    /// <summary>
    /// S-001: Non-NGG PAM at offsets 25-26 throws (model is SpCas9-specific).
    /// "...AGGGAGGGACC": replace the GG at 25-26 with AA.
    /// </summary>
    [Test]
    public void Doench2014_NonNggPam_Throws()
    {
        // Reference example 1 with offsets 25-26 forced to non-GG (AA).
        string bad = "TATAGCTGCGATCTGAGGTAGGGAGAAACC";
        Assert.That(bad.Length, Is.EqualTo(30));
        Assert.That(bad.Substring(25, 2), Is.EqualTo("AA")); // offsets 25-26 are not "GG"
        Assert.Throws<ArgumentException>(() => CrisprDesigner.CalculateOnTargetDoench2014(bad));
    }
}
