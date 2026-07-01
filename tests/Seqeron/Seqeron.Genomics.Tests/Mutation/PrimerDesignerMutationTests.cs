using NUnit.Framework;
using FluentAssertions;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// Targeted mutation-killing tests for PrimerDesigner.cs (checklist 04 rows 21/22/23:
/// PRIMER-TM-001, PRIMER-DESIGN-001, PRIMER-STRUCT-001).
///
/// These pin the documented business logic that the canonical suite left as Stryker
/// "survivors": the primer-quality SCORING formula, the 3' GC-clamp rule, the
/// dinucleotide-repeat acceptance boundary, and the candidate-generation window
/// boundaries. Each test encodes the published rule (not the production result) so an
/// arithmetic/relational/logical mutant of the source diverges and is killed.
///
/// Scoring model (PrimerDesigner.CalculatePrimerScore, transcribed from source docs):
///   score = 100
///         − |len − OptimalLength| · 2
///         − |Tm  − OptimalTm|     · 2
///         − |GC  − 50|            · 0.5
///         − longestHomopolymer    · 5
///         + (3'-base ∈ {G,C} ? 5 : 0)      , clamped to ≥ 0, rounded to 2 dp.
/// </summary>
[TestFixture]
public class PrimerDesignerMutationTests
{
    // Independent ground-truth scoring formula (the theory), built only from the
    // separately-tested primitive metrics. Kills the arithmetic/clamp mutants in
    // CalculatePrimerScore because the test arithmetic is never mutated.
    private static double ExpectedScore(string seq, PrimerParameters p)
    {
        double score = 100;
        score -= Math.Abs(seq.Length - p.OptimalLength) * 2;

        double tm = PrimerDesigner.CalculateMeltingTemperature(seq);
        score -= Math.Abs(tm - p.OptimalTm) * 2;

        double gc = PrimerDesigner.CalculateGcContent(seq);
        score -= Math.Abs(gc - 50) * 0.5;

        int homopolymer = PrimerDesigner.FindLongestHomopolymer(seq);
        score -= homopolymer * 5;

        if (seq.Length >= 2 && (seq[^1] == 'G' || seq[^1] == 'C'))
            score += 5;

        return Math.Round(Math.Max(0, score), 2);
    }

    // ── Scoring: exact value pins every penalty term and the GC-clamp bonus ──────────
    // len ≠ OptimalLength, GC ≠ 50, homopolymer ≥ 2 in each construct so the *2 / *0.5 / *5
    // operators have non-zero contributions (×→÷ etc. all change the result).

    [Test]
    [TestCase("ACGGTACGGTACGGTACGGTAG")] // 22 bp, GC=60%, ends in G  → clamp applied
    [TestCase("ACGGTACGGTACGGTACGGTAC")] // 22 bp, GC=60%, ends in C  → clamp applied
    [TestCase("ACGGTACGGTACGGTACGGTAA")] // 22 bp, GC=55%, ends in A  → clamp NOT applied
    public void EvaluatePrimer_Score_MatchesPublishedScoringFormula(string seq)
    {
        var candidate = PrimerDesigner.EvaluatePrimer(seq, position: 0, isForward: true);

        candidate.Score.Should().Be(ExpectedScore(seq, PrimerDesigner.DefaultParameters));
    }

    [Test]
    public void EvaluatePrimer_GcClampBonus_OnlyAppliesToGorCThreePrimeEnd()
    {
        // Two primers identical except for the 3' base. The only score difference must be
        // the +5 GC-clamp bonus (kills the `last=='G' || last=='C'` logical/equality mutants).
        const string baseBody = "ACGGTACGGTACGGTACGGTA"; // 21 bp body
        var endsG = PrimerDesigner.EvaluatePrimer(baseBody + "G", 0, true);
        var endsT = PrimerDesigner.EvaluatePrimer(baseBody + "T", 0, true);

        // G/C ending and A/T ending share length(22)+homopolymer; GC% differs by one base
        // so compare against the formula rather than a raw constant.
        endsG.Score.Should().Be(ExpectedScore(baseBody + "G", PrimerDesigner.DefaultParameters));
        endsT.Score.Should().Be(ExpectedScore(baseBody + "T", PrimerDesigner.DefaultParameters));

        // The G-ending primer must score strictly higher once GC%/Tm differences are removed:
        // construct a controlled pair where only the clamp differs.
        const string body2 = "ATGCATGCATGCATGCATGCA"; // 21 bp, balanced
        var pG = PrimerDesigner.EvaluatePrimer(body2 + "G", 0, true);
        var pA = PrimerDesigner.EvaluatePrimer(body2 + "A", 0, true);
        (pG.Score - pA.Score).Should().BeGreaterThan(0,
            "a G/C 3' terminus earns the +5 clamp bonus that an A/T terminus does not");
    }

    // ── 3' GC clamp rule (EvaluatePrimer issue list, Avoid3PrimeGC = true) ────────────

    private static PrimerParameters WithGcClampRequired() =>
        PrimerDesigner.DefaultParameters with { Avoid3PrimeGC = true };

    [Test]
    public void EvaluatePrimer_Avoid3PrimeGc_FlagsPrimerWithNoGcInLastTwoBases()
    {
        // last2 = "AT" → gcCount == 0 → "No GC clamp" issue MUST be present.
        // Kills the `gcCount == 0` → `!= 0` mutant and the `c=='G'||c=='C'` → `!='G'`/`!='C'` mutants.
        var primer = "GCGCGCGCGCGCGCGCGCAT"; // 20 bp, ends "AT"
        var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true, WithGcClampRequired());

        candidate.Issues.Should().Contain(i => i.Contains("GC clamp"),
            "a 3' end without G or C must be flagged when Avoid3PrimeGC is requested");
    }

    [Test]
    public void EvaluatePrimer_Avoid3PrimeGc_AcceptsPrimerWithGcInLastTwoBases()
    {
        // last2 = "GC" → gcCount == 2 → NO "No GC clamp" issue.
        // Kills the `c=='G' || c=='C'` → `&&` mutant (which would force gcCount to 0).
        var primer = "ATATATATATATATATATGC"; // 20 bp, ends "GC"
        var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true, WithGcClampRequired());

        candidate.Issues.Should().NotContain(i => i.Contains("GC clamp"),
            "a 3' end containing G or C satisfies the GC clamp");
    }

    // ── Dinucleotide-repeat acceptance boundary (EvaluatePrimer, strict '>') ──────────

    [Test]
    public void EvaluatePrimer_DinucleotideRepeat_ExactlyMaxIsAccepted()
    {
        // MaxDinucleotideRepeats = 4. A primer whose longest dinucleotide repeat is EXACTLY 4
        // must NOT be flagged (rule is strict '>'). Kills the `> Max` → `>= Max` mutant.
        const string primer = "GCTATATATATCGGCTAGCAT"; // contains AT×4 ("ATATATAT"), nothing longer
        PrimerDesigner.FindLongestDinucleotideRepeat(primer).Should().Be(4,
            "guard: the construct must have a longest dinucleotide repeat of exactly 4");

        var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true);

        candidate.Issues.Should().NotContain(i => i.Contains("Dinucleotide repeat"),
            "exactly MaxDinucleotideRepeats (4) is within tolerance, not exceeding it");
    }

    [Test]
    public void EvaluatePrimer_DinucleotideRepeat_AboveMaxIsRejected()
    {
        // 5 repeats > 4 → must be flagged (anchors the other side of the boundary).
        const string primer = "GCTATATATATATCGGCTAGCA"; // AT×5
        PrimerDesigner.FindLongestDinucleotideRepeat(primer).Should().Be(5);

        var candidate = PrimerDesigner.EvaluatePrimer(primer, 0, true);

        candidate.Issues.Should().Contain(i => i.Contains("Dinucleotide repeat"));
    }

    // ── GeneratePrimerCandidates window boundaries (inclusive '<=') ───────────────────

    [Test]
    public void GeneratePrimerCandidates_EmitsLengthsUpToMaxLengthInclusive()
    {
        // Region [0,30) with Min=18, Max=25. A candidate of length EXACTLY MaxLength(25) must
        // be produced. Kills `len <= MaxLength` → `len < MaxLength`.
        var template = new DnaSequence(new string('A', 10) + "CGTACGTACGTACGTACGTACGTACGTACGTA"); // ≥ 42 bp
        var candidates = PrimerDesigner.GeneratePrimerCandidates(template, 0, 30).ToList();

        candidates.Should().Contain(c => c.Length == 25, "MaxLength is an inclusive bound");
        candidates.Should().Contain(c => c.Length == 18, "MinLength candidates are emitted too");
    }

    [Test]
    public void GeneratePrimerCandidates_LastStartReachingRegionEndIsInclusive()
    {
        // The final window must be allowed to end exactly at regionEnd, and the final start at
        // regionEnd-MinLength must still emit. Kills `start + len <= regionEnd` → `<` and
        // `start + MinLength <= regionEnd` → `<`.
        var template = new DnaSequence(new string('A', 10) + "CGTACGTACGTACGTACGTACGTACGTACGTA");
        int regionEnd = 30;
        var candidates = PrimerDesigner.GeneratePrimerCandidates(template, 0, regionEnd).ToList();

        candidates.Should().Contain(c => c.Position + c.Length == regionEnd,
            "a window may end exactly at regionEnd (inclusive upper bound)");
        candidates.Should().Contain(c => c.Position == regionEnd - PrimerDesigner.DefaultParameters.MinLength,
            "the last legal start (regionEnd − MinLength) must still produce a candidate");
    }

    // ── Hairpin detection: exercise BOTH the O(n²) and suffix-tree paths ──────────────

    [Test]
    public void HasHairpinPotential_ShortSelfComplementaryStem_DetectsHairpin()
    {
        // stem GGGGG / loop AAA / stem CCCCC (revcomp(CCCCC)=GGGGG) → hairpin, O(n²) path (<100 bp).
        PrimerDesigner.HasHairpinPotential("GGGGGAAACCCCC").Should().BeTrue();
    }

    [Test]
    public void HasHairpinPotential_ShortNonComplementary_NoHairpin()
    {
        PrimerDesigner.HasHairpinPotential("AAAAAAAAAAAAA").Should().BeFalse();
    }

    [Test]
    public void HasHairpinPotential_LongSequenceWithPlantedHairpin_UsesSuffixTreePath()
    {
        // ≥100 bp routes through the suffix-tree implementation. Poly-A padding cannot form
        // spurious stems (A pairs only with T, absent), so detection is driven solely by the
        // planted GGGGG…CCCCC stem-loop. Kills the NoCoverage block on the suffix-tree branch.
        string seq = new string('A', 60) + "GGGGGAAACCCCC" + new string('A', 60); // 133 bp
        PrimerDesigner.HasHairpinPotential(seq).Should().BeTrue();
    }

    [Test]
    public void HasHairpinPotential_LongSequenceNoComplementarity_NoHairpin()
    {
        PrimerDesigner.HasHairpinPotential(new string('A', 120)).Should().BeFalse();
    }
}
