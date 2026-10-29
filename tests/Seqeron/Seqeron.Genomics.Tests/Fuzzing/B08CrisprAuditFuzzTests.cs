using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Heavy-tier (fuzz) coverage for the CRISPR surface that B08 audit rounds 4-5 added or tightened and the older
/// fuzz suites do not reach: the MIT / CFD genome aggregates, the SaCas9 hit score and the validation of the
/// aggregate arguments.
///
/// Every expectation is a published definition or a documented argument contract, not an observed output:
///   - the MIT/Hsu 2013 and CFD (Doench 2016) off-target models are SpCas9-specific, so a PAM-before-target
///     (Cas12a-class / CasX) system has no score and must be refused, never scored with a wrong geometry
///     (CRISPOR reports no MIT or CFD off-target score for these systems);
///   - the genome aggregates take the same 0-5 mismatch range as <c>FindOffTargets</c>;
///   - the Tycko 2018 SaCas9 model gives 100 x the product of per-position coefficients, with no contribution from
///     the 5'-most position ("Mismatches at the 5' position are not truly determined by the model").
///
/// Test Units: CRISPR-OFF-001 (+ CRISPR-PAM-001 for the cut-site range).
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("MolTools")]
[Category("CRISPR-OFF-001")]
public class B08CrisprAuditFuzzTests
{
    private const string Guide20 = "GACGCATAAAGATGAGACGC";
    private const string Guide21 = "GACGCATAAAGATGAGACGCT";

    private static readonly DnaSequence Genome = new(
        "TTGACGCATAAAGATGAGACGCTGGAATTCCTAGGACGCATAAAGATGAGACGCAGGTTACCGGATCCAAGGACGCATAAAGATGAGACGCTAGCCAAA");

    private static string GuideFor(CrisprSystemType type) =>
        CrisprDesigner.GetSystem(type).GuideLength == 21 ? Guide21 : Guide20;

    public static IEnumerable<CrisprSystemType> AllSystems =>
        Enum.GetValues<CrisprSystemType>();

    public static IEnumerable<CrisprSystemType> PamBeforeSystems =>
        AllSystems.Where(t => !CrisprDesigner.GetSystem(t).PamAfterTarget);

    public static IEnumerable<CrisprSystemType> PamAfterSystems =>
        AllSystems.Where(t => CrisprDesigner.GetSystem(t).PamAfterTarget);

    public static IEnumerable<CrisprSystemType> PamAfter20NtSystems =>
        PamAfterSystems.Where(t => CrisprDesigner.GetSystem(t).GuideLength == 20);

    public static IEnumerable<CrisprSystemType> PamAfter21NtSystems =>
        PamAfterSystems.Where(t => CrisprDesigner.GetSystem(t).GuideLength == 21);

    #region MIT / CFD are refused for every PAM-before-target system (split out, not relaxed)

    /// <summary>
    /// The MIT/Hsu 2013 model is indexed PAM-proximal-last, so for a PAM-before-target system every position weight would
    /// land on the wrong base. All five genome aggregates (plain, alternative-PAM, CFD plain, CFD alternative-PAM) and the
    /// PAM-annotated hit-list overload refuse such a system with ArgumentException instead of returning a number.
    /// </summary>
    [TestCaseSource(nameof(PamBeforeSystems))]
    public void Aggregates_PamBeforeTargetSystem_AreRefused_NotScored(CrisprSystemType type)
    {
        string guide = GuideFor(type);

        var calls = new (string Name, Action Call)[]
        {
            ("MIT genome", () => CrisprDesigner.CalculateMitSpecificityScore(guide, Genome, 4, type)),
            ("MIT alt-PAM", () => CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(guide, Genome, 4, type)),
            ("CFD genome", () => CrisprDesigner.CalculateCfdSpecificityScore(guide, Genome, 4, type)),
            ("CFD alt-PAM", () => CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(guide, Genome, 4, type)),
            ("MIT hit list", () => CrisprDesigner.CalculateMitSpecificityScore(guide, Array.Empty<PamAnnotatedOffTarget>(), type)),
        };

        foreach (var (name, call) in calls)
            call.Should().Throw<ArgumentException>($"{name} must refuse {type}: the model is SpCas9-specific")
                .Which.ParamName.Should().Be("systemType", $"{name} blames the system argument");
    }

    /// <summary>
    /// The refusal is a split, not a blanket: every PAM-after-target system (the SpCas9 family and the SaCas9 family)
    /// still scores through the MIT genome aggregate with a guide of its own length, and the result is a finite value in
    /// (0, 100] - 100 / (100 + sum) x 100 with a finite non-negative sum.
    /// </summary>
    [TestCaseSource(nameof(PamAfterSystems))]
    public void MitAggregate_PamAfterTargetSystem_IsScored_InTheDocumentedRange(CrisprSystemType type)
    {
        double score = CrisprDesigner.CalculateMitSpecificityScore(GuideFor(type), Genome, 4, type);

        double.IsFinite(score).Should().BeTrue();
        score.Should().BeGreaterThan(0).And.BeLessThanOrEqualTo(100);
    }

    /// <summary>The CFD model is defined for a 20-nt protospacer only: the 20-nt PAM-after systems score.</summary>
    [TestCaseSource(nameof(PamAfter20NtSystems))]
    public void CfdAggregate_TwentyNtPamAfterSystem_IsScored(CrisprSystemType type)
    {
        double score = CrisprDesigner.CalculateCfdSpecificityScore(Guide20, Genome, 4, type);

        score.Should().BeGreaterThan(0).And.BeLessThanOrEqualTo(100);
    }

    /// <summary>... and the 21-nt SaCas9 family is refused by the CFD aggregates, blaming the system.</summary>
    [TestCaseSource(nameof(PamAfter21NtSystems))]
    public void CfdAggregate_TwentyOneNtSystem_IsRefused(CrisprSystemType type)
    {
        var plain = () => CrisprDesigner.CalculateCfdSpecificityScore(Guide21, Genome, 4, type);
        var altPam = () => CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(Guide21, Genome, 4, type);

        plain.Should().Throw<ArgumentException>().Which.ParamName.Should().Be("systemType");
        altPam.Should().Throw<ArgumentException>().Which.ParamName.Should().Be("systemType");
    }

    #endregion

    #region maxMismatches outside 0-5

    private static Arbitrary<int> OutOfRangeMismatches() =>
        Gen.OneOf(Gen.Choose(int.MinValue, -1), Gen.Choose(6, int.MaxValue)).ToArbitrary();

    /// <summary>
    /// All four genome aggregates reject any <c>maxMismatches</c> outside 0-5 eagerly with ArgumentOutOfRangeException
    /// (the same range as <c>FindOffTargets</c>), whatever the value - including the int extremes.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GenomeAggregates_MaxMismatchesOutsideZeroToFive_ThrowArgumentOutOfRange()
    {
        return Prop.ForAll(OutOfRangeMismatches(), mm =>
        {
            Action[] calls =
            [
                () => CrisprDesigner.CalculateMitSpecificityScore(Guide20, Genome, mm),
                () => CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(Guide20, Genome, mm),
                () => CrisprDesigner.CalculateCfdSpecificityScore(Guide20, Genome, mm),
                () => CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(Guide20, Genome, mm),
            ];

            return calls.All(call =>
            {
                try { call(); return false; }
                catch (ArgumentOutOfRangeException) { return true; }
            }).Label($"maxMismatches = {mm}");
        });
    }

    /// <summary>The boundary values 0 and 5 are inside the range: they score (no throw), each in (0, 100].</summary>
    [TestCase(0)]
    [TestCase(5)]
    public void GenomeAggregates_MaxMismatchesAtTheRangeEnds_AreAccepted(int maxMismatches)
    {
        double mit = CrisprDesigner.CalculateMitSpecificityScore(Guide20, Genome, maxMismatches);
        double cfd = CrisprDesigner.CalculateCfdSpecificityScore(Guide20, Genome, maxMismatches);

        mit.Should().BeGreaterThan(0).And.BeLessThanOrEqualTo(100);
        cfd.Should().BeGreaterThan(0).And.BeLessThanOrEqualTo(100);
    }

    #endregion

    #region SaCas9 hit score

    private static readonly char[] Bases = ['A', 'C', 'G', 'T'];

    private static Arbitrary<(string Guide, string OffTarget)> SaPairs() =>
        Gen.Elements(Bases).ArrayOf(21)
            .SelectMany(g => Gen.Elements(Bases).ArrayOf(21).Select(o => (new string(g), new string(o))))
            .ToArbitrary();

    /// <summary>
    /// Range, case invariance and the 5' rule of the Tycko 2018 hit score, over random 21-mer pairs:
    /// the score is in [0, 100]; it does not depend on the case of either argument; and it is exactly 100 when the pair
    /// differs at index 0 only (that position contributes nothing) - including identical pairs.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property SaCas9HitScore_RangeCaseInvariance_AndFivePrimeRule()
    {
        return Prop.ForAll(SaPairs(), pair =>
        {
            var (guide, off) = pair;
            double score = CrisprDesigner.CalculateSaCas9HitScore(guide, off);
            double lower = CrisprDesigner.CalculateSaCas9HitScore(guide.ToLowerInvariant(), off.ToLowerInvariant());
            double mixed = CrisprDesigner.CalculateSaCas9HitScore(guide.ToLowerInvariant(), off);

            // Force the pair to differ at index 0 only.
            char other = guide[0] == 'A' ? 'C' : 'A';
            string onlyFirstDiffers = other + guide[1..];
            double first = CrisprDesigner.CalculateSaCas9HitScore(guide, onlyFirstDiffers);

            return (score >= 0 && score <= 100).Label($"range: {score}")
                .And((score == lower && score == mixed).Label("case invariance"))
                .And((first == 100.0).Label($"5' position contributes nothing: {first}"));
        });
    }

    /// <summary>
    /// "Equals 100 only when the sole mismatches are at index 0": a mismatch at any other position strictly lowers the
    /// activity (every trained coefficient is below 1), for every position 1..20 and every substitution.
    /// </summary>
    [Test]
    public void SaCas9HitScore_AnyMismatchBeyondIndexZero_StrictlyLowersTheScore()
    {
        for (int i = 1; i < 21; i++)
            foreach (char sub in Bases.Where(b => b != Guide21[i]))
            {
                string off = Guide21[..i] + sub + Guide21[(i + 1)..];

                double score = CrisprDesigner.CalculateSaCas9HitScore(Guide21, off);

                score.Should().BeLessThan(100, $"substitution {Guide21[i]}>{sub} at index {i}");
                score.Should().BeGreaterThan(0);
            }
    }

    /// <summary>A gap, an ambiguity code or a non-nucleotide character anywhere is refused (the reference refuses the gap).</summary>
    [TestCase('-')]
    [TestCase('N')]
    [TestCase('X')]
    [TestCase('1')]
    public void SaCas9HitScore_InvalidBase_IsRefusedAtAnyPosition(char junk)
    {
        for (int i = 0; i < 21; i++)
        {
            string bad = Guide21[..i] + junk + Guide21[(i + 1)..];

            var inGuide = () => CrisprDesigner.CalculateSaCas9HitScore(bad, Guide21);
            var inOff = () => CrisprDesigner.CalculateSaCas9HitScore(Guide21, bad);

            inGuide.Should().Throw<ArgumentException>($"index {i} of the guide");
            inOff.Should().Throw<ArgumentException>($"index {i} of the off-target");
        }
    }

    #endregion

    #region GetCutSite stays inside the target sequence

    private static Arbitrary<string> Dna(int min, int max) =>
        Gen.Elements(Bases).ArrayOf()
            .Where(a => a.Length >= min && a.Length <= max)
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// For every system and every PAM site FindPamSites reports in a random sequence, the predicted cut lies strictly
    /// inside the sequence, in [1, Length - 1] - a cut at 0 or at Length would be a cut outside the molecule. This is the
    /// narrowing of the older [0, Length] bound.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GetCutSite_ForEveryReportedPamSite_IsInside_OneToLengthMinusOne()
    {
        return Prop.ForAll(Dna(60, 160), seq =>
        {
            var dna = new DnaSequence(seq);
            foreach (var type in Enum.GetValues<CrisprSystemType>())
                foreach (var site in CrisprDesigner.FindPamSites(dna, type))
                {
                    int cut = CrisprDesigner.GetCutSite(site);
                    if (cut < 1 || cut > seq.Length - 1)
                        return false.Label($"{type} {(site.IsForwardStrand ? "+" : "-")} PAM at {site.Position}: cut {cut} outside [1, {seq.Length - 1}]");
                }

            return true.ToProperty();
        });
    }

    /// <summary>
    /// The cut follows the documented geometry for every system: the data-driven <c>CutOffset</c> c gives
    /// <c>Position - c</c> / <c>Position + pamLength + c - 1</c> for PAM-after systems and <c>Position + pamLength + c</c> /
    /// <c>Position - c - 1</c> for PAM-before systems (forward / reverse). Independent re-derivation of GetCutSite.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property GetCutSite_MatchesTheDocumentedGeometry_ForEverySystem()
    {
        return Prop.ForAll(Dna(60, 160), seq =>
        {
            var dna = new DnaSequence(seq);
            foreach (var type in Enum.GetValues<CrisprSystemType>())
            {
                var system = CrisprDesigner.GetSystem(type);
                foreach (var site in CrisprDesigner.FindPamSites(dna, type))
                {
                    int pam = site.PamSequence.Length;
                    int c = system.CutOffset;
                    int expected = system.PamAfterTarget
                        ? (site.IsForwardStrand ? site.Position - c : site.Position + pam + c - 1)
                        : (site.IsForwardStrand ? site.Position + pam + c : site.Position - c - 1);
                    if (CrisprDesigner.GetCutSite(site) != expected)
                        return false.Label($"{type}: expected {expected}, got {CrisprDesigner.GetCutSite(site)}");
                }
            }

            return true.ToProperty();
        });
    }

    #endregion
}
