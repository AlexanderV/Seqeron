using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// MIRNA-SEED-001: Seed Sequence Analysis
/// Canonical test file for GetSeedSequence, CreateMiRna, CompareSeedRegions.
/// Evidence: miRBase, TargetScan, Bartel (2009), Lewis (2005).
/// </summary>
[TestFixture]
public class MiRnaAnalyzer_SeedAnalysis_Tests
{
    #region Reference Data (miRBase)

    // hsa-let-7a-5p: UGAGGUAGUAGGUUGUAUAGUU → seed GAGGUAG
    private const string Let7a_Sequence = "UGAGGUAGUAGGUUGUAUAGUU";
    private const string Let7a_Seed = "GAGGUAG";

    // hsa-let-7b-5p: UGAGGUAGUAGGUUGUGUGGUU → seed GAGGUAG
    private const string Let7b_Sequence = "UGAGGUAGUAGGUUGUGUGGUU";
    private const string Let7b_Seed = "GAGGUAG";

    // hsa-let-7c-5p: UGAGGUAGUAGGUUGUAUGGUU → seed GAGGUAG
    private const string Let7c_Sequence = "UGAGGUAGUAGGUUGUAUGGUU";
    private const string Let7c_Seed = "GAGGUAG";

    // hsa-miR-21-5p: UAGCUUAUCAGACUGAUGUUGA → seed AGCUUAU
    private const string MiR21_Sequence = "UAGCUUAUCAGACUGAUGUUGA";
    private const string MiR21_Seed = "AGCUUAU";

    // hsa-miR-155-5p: UUAAUGCUAAUUGUGAUAGGGGU → seed UAAUGCU
    private const string MiR155_Sequence = "UUAAUGCUAAUUGUGAUAGGGGU";
    private const string MiR155_Seed = "UAAUGCU";

    // hsa-miR-1-3p: UGGAAUGUAAAGAAGUAUGUAU → seed GGAAUGU
    private const string MiR1_Sequence = "UGGAAUGUAAAGAAGUAUGUAU";
    private const string MiR1_Seed = "GGAAUGU";

    #endregion

    #region M-001: GetSeedSequence — Known miRNA returns correct seed

    [Test]
    public void GetSeedSequence_Let7a_ReturnsGAGGUAG()
    {
        // Evidence: miRBase hsa-let-7a-5p, seed = positions 2-8
        string seed = GetSeedSequence(Let7a_Sequence);

        Assert.That(seed, Is.EqualTo(Let7a_Seed));
    }

    #endregion

    #region M-002: GetSeedSequence — Multiple known miRNAs

    [TestCase("UAGCUUAUCAGACUGAUGUUGA", "AGCUUAU", TestName = "miR-21")]
    [TestCase("UUAAUGCUAAUUGUGAUAGGGGU", "UAAUGCU", TestName = "miR-155")]
    [TestCase("UGGAAUGUAAAGAAGUAUGUAU", "GGAAUGU", TestName = "miR-1")]
    public void GetSeedSequence_KnownMiRNAs_ReturnsExpectedSeed(string sequence, string expectedSeed)
    {
        // Evidence: miRBase reference sequences
        string seed = GetSeedSequence(sequence);

        Assert.That(seed, Is.EqualTo(expectedSeed));
    }

    #endregion

    #region M-003: GetSeedSequence — Seed is always 7 nucleotides

    [TestCase("UGAGGUAGUAGGUUGUAUAGUU", TestName = "let-7a (22 nt)")]
    [TestCase("UAGCUUAUCAGACUGAUGUUGA", TestName = "miR-21 (22 nt)")]
    [TestCase("UUAAUGCUAAUUGUGAUAGGGGU", TestName = "miR-155 (23 nt)")]
    [TestCase("ABCDEFGH", TestName = "minimum 8 nt")]
    [TestCase("ABCDEFGHIJKLMNOP", TestName = "16 nt")]
    public void GetSeedSequence_ValidInput_ReturnsExactly7Characters(string sequence)
    {
        // Evidence: TargetScan: seed region is always 7 nt (positions 2-8)
        string seed = GetSeedSequence(sequence);

        Assert.That(seed, Has.Length.EqualTo(7));
    }

    #endregion

    #region M-004: GetSeedSequence — Empty/null input returns empty

    [Test]
    public void GetSeedSequence_NullInput_ReturnsEmpty()
    {
        Assert.That(GetSeedSequence(null!), Is.Empty);
    }

    [Test]
    public void GetSeedSequence_EmptyString_ReturnsEmpty()
    {
        Assert.That(GetSeedSequence(""), Is.Empty);
    }

    #endregion

    #region M-005: GetSeedSequence — Short sequence (< 8 nt) returns empty

    [TestCase("U", TestName = "1 nt")]
    [TestCase("UA", TestName = "2 nt")]
    [TestCase("UAGCA", TestName = "5 nt")]
    [TestCase("UAGCAUU", TestName = "7 nt")]
    public void GetSeedSequence_ShortSequence_ReturnsEmpty(string sequence)
    {
        // Evidence: Cannot extract 7-nt seed from sequence shorter than 8 nt
        string seed = GetSeedSequence(sequence);

        Assert.That(seed, Is.Empty);
    }

    #endregion

    #region M-006: GetSeedSequence — Case normalization

    [Test]
    public void GetSeedSequence_LowercaseInput_ReturnsUppercaseSeed()
    {
        // Evidence: Implementation normalizes to uppercase
        string seed = GetSeedSequence("ugagguaguagguuguauaguu");

        Assert.That(seed, Is.EqualTo("GAGGUAG"));
    }

    [Test]
    public void GetSeedSequence_MixedCaseInput_ReturnsUppercaseSeed()
    {
        string seed = GetSeedSequence("uGaGgUaGuAgGuUgUaUaGuu");

        Assert.That(seed, Is.EqualTo("GAGGUAG"));
    }

    #endregion

    #region M-007: CreateMiRna — Factory produces correct record

    [Test]
    public void CreateMiRna_Let7a_ProducesCorrectRecord()
    {
        // Evidence: miRBase hsa-let-7a-5p
        var mirna = CreateMiRna("let-7a", Let7a_Sequence);

        Assert.Multiple(() =>
        {
            Assert.That(mirna.Name, Is.EqualTo("let-7a"));
            Assert.That(mirna.Sequence, Is.EqualTo(Let7a_Sequence));
            Assert.That(mirna.SeedSequence, Is.EqualTo(Let7a_Seed));
            Assert.That(mirna.SeedStart, Is.EqualTo(1));
            Assert.That(mirna.SeedEnd, Is.EqualTo(7));
        });
    }

    #endregion

    #region M-008: CreateMiRna — DNA input converted to RNA

    [Test]
    public void CreateMiRna_DnaInput_ConvertsToRna()
    {
        // Evidence: Standard T→U conversion
        string dnaSequence = "TGAGGTAGTAGGTTGTATAGTT"; // DNA version of let-7a
        var mirna = CreateMiRna("let-7a-dna", dnaSequence);

        Assert.Multiple(() =>
        {
            Assert.That(mirna.Sequence, Does.Not.Contain("T"));
            Assert.That(mirna.Sequence, Does.Contain("U"));
            Assert.That(mirna.Sequence, Is.EqualTo(Let7a_Sequence));
            Assert.That(mirna.SeedSequence, Is.EqualTo(Let7a_Seed));
        });
    }

    #endregion

    #region M-009: CompareSeedRegions — Identical seeds

    [Test]
    public void CompareSeedRegions_IdenticalSeeds_ZeroMismatchesSameFamily()
    {
        // Evidence: TargetScan: same seed = same family
        var mirna1 = CreateMiRna("let-7a", Let7a_Sequence);
        var mirna2 = CreateMiRna("let-7b", Let7b_Sequence);

        var result = CompareSeedRegions(mirna1, mirna2);

        Assert.Multiple(() =>
        {
            Assert.That(result.Matches, Is.EqualTo(7));
            Assert.That(result.Mismatches, Is.EqualTo(0));
            Assert.That(result.IsSameFamily, Is.True);
        });
    }

    #endregion

    #region M-010: CompareSeedRegions — Different seeds

    [Test]
    public void CompareSeedRegions_Let7a_Vs_MiR21_CorrectHammingDistance()
    {
        // Evidence: Hamming distance between GAGGUAG and AGCUUAU
        var let7a = CreateMiRna("let-7a", Let7a_Sequence);
        var mir21 = CreateMiRna("miR-21", MiR21_Sequence);

        var result = CompareSeedRegions(let7a, mir21);

        // Manually compute: G≠A, A=A, G≠C, G≠U, U≠U→match? G vs U: G≠U, A≠A→match? A vs A: match, G≠U
        // GAGGUAG vs AGCUUAU:
        // G≠A, A≠G, G≠C, G≠U, U≠U→U=U match, A=A match, G≠U
        // Pos 0: G vs A → mismatch
        // Pos 1: A vs G → mismatch
        // Pos 2: G vs C → mismatch
        // Pos 3: G vs U → mismatch
        // Pos 4: U vs U → match
        // Pos 5: A vs A → match
        // Pos 6: G vs U → mismatch
        // Matches=2, Mismatches=5

        Assert.Multiple(() =>
        {
            Assert.That(result.Matches, Is.EqualTo(2));
            Assert.That(result.Mismatches, Is.EqualTo(5));
            Assert.That(result.IsSameFamily, Is.False);
        });
    }

    #endregion

    #region M-011: CompareSeedRegions — Single mismatch

    [Test]
    public void CompareSeedRegions_SingleMismatch_CorrectlyReported()
    {
        // Create two miRNAs that differ in exactly one seed position
        // Seed of "XGAGGUAGUYYY..." = GAGGUAG (let-7a seed)
        // Modify one position: change position 2-8 to GAGUUAG (pos 4: G→U)
        // Build sequences: U + seed + padding
        var mirna1 = CreateMiRna("ref", "UGAGGUAGXXXXXXXXXXXXXXX");
        var mirna2 = CreateMiRna("mut", "UGAGUUAGXXXXXXXXXXXXXXX");

        // mirna1 seed: GAGGUAG
        // mirna2 seed: GAGUUAG
        // Diff at pos 2: G vs U → 1 mismatch, 6 matches

        var result = CompareSeedRegions(mirna1, mirna2);

        Assert.Multiple(() =>
        {
            Assert.That(result.Matches, Is.EqualTo(6));
            Assert.That(result.Mismatches, Is.EqualTo(1));
            Assert.That(result.IsSameFamily, Is.False);
        });
    }

    #endregion

    #region M-012: CompareSeedRegions — Empty seed handling

    [Test]
    public void CompareSeedRegions_EmptySeed_ReturnsGracefulResult()
    {
        // Evidence: Defensive programming (ASSUMPTION)
        var validMirna = CreateMiRna("let-7a", Let7a_Sequence);
        var shortMirna = CreateMiRna("short", "UAGCA"); // Too short → empty seed

        var result = CompareSeedRegions(validMirna, shortMirna);

        Assert.Multiple(() =>
        {
            Assert.That(result.Matches, Is.EqualTo(0));
            Assert.That(result.Mismatches, Is.EqualTo(0));
            Assert.That(result.IsSameFamily, Is.False);
        });
    }

    [Test]
    public void CompareSeedRegions_BothEmptySeeds_ReturnsGracefulResult()
    {
        var short1 = CreateMiRna("s1", "UGA");
        var short2 = CreateMiRna("s2", "UAG");

        var result = CompareSeedRegions(short1, short2);

        Assert.Multiple(() =>
        {
            Assert.That(result.Matches, Is.EqualTo(0));
            Assert.That(result.Mismatches, Is.EqualTo(0));
            Assert.That(result.IsSameFamily, Is.False);
        });
    }

    #endregion

    #region S-001: GetSeedSequence — Boundary (exactly 8 nt)

    [Test]
    public void GetSeedSequence_Exactly8Nucleotides_ReturnsSeed()
    {
        // Evidence: Minimum viable input for seed extraction
        string seed = GetSeedSequence("ABCDEFGH");

        Assert.That(seed, Is.EqualTo("BCDEFGH"));
    }

    #endregion

    #region S-002: CreateMiRna — Mixed case DNA input

    [Test]
    public void CreateMiRna_MixedCaseDna_ConvertsCorrectly()
    {
        // Evidence: Robustness requirement
        string mixedDna = "tGaGgTaGtAgGtTgTaTaGtT";
        var mirna = CreateMiRna("mixed", mixedDna);

        Assert.Multiple(() =>
        {
            Assert.That(mirna.Sequence, Is.EqualTo(Let7a_Sequence));
            Assert.That(mirna.SeedSequence, Is.EqualTo(Let7a_Seed));
        });
    }

    #endregion

    #region S-003: let-7 family members share identical seed

    [Test]
    public void GetSeedSequence_Let7Family_AllShareSameSeed()
    {
        // Evidence: miRBase let-7a/b/c all have seed GAGGUAG
        string seedA = GetSeedSequence(Let7a_Sequence);
        string seedB = GetSeedSequence(Let7b_Sequence);
        string seedC = GetSeedSequence(Let7c_Sequence);

        Assert.Multiple(() =>
        {
            Assert.That(seedA, Is.EqualTo(Let7a_Seed));
            Assert.That(seedB, Is.EqualTo(Let7b_Seed));
            Assert.That(seedC, Is.EqualTo(Let7c_Seed));
            Assert.That(seedA, Is.EqualTo(seedB));
            Assert.That(seedB, Is.EqualTo(seedC));
        });
    }

    #endregion

    #region S-004: CompareSeedRegions — Completely different seeds

    [Test]
    public void CompareSeedRegions_CompletelyDifferentSeeds_MaxMismatches()
    {
        // Evidence: Maximal Hamming distance = seed length
        // Construct two miRNAs with maximally different seeds
        // Seed1: AAAAAAA, Seed2: UUUUUUU
        var mirna1 = CreateMiRna("all-A", "UAAAAAAAXXXXXXXXXXX");
        var mirna2 = CreateMiRna("all-U", "UUUUUUUUXXXXXXXXXXX");

        var result = CompareSeedRegions(mirna1, mirna2);

        Assert.Multiple(() =>
        {
            Assert.That(result.Matches, Is.EqualTo(0));
            Assert.That(result.Mismatches, Is.EqualTo(7));
            Assert.That(result.IsSameFamily, Is.False);
        });
    }

    #endregion

    #region C-001: Seed extraction is a pure function (deterministic)

    [Test]
    public void GetSeedSequence_CalledMultipleTimes_ReturnsSameResult()
    {
        // Evidence: Implementation requirement — pure function
        string result1 = GetSeedSequence(Let7a_Sequence);
        string result2 = GetSeedSequence(Let7a_Sequence);
        string result3 = GetSeedSequence(Let7a_Sequence);

        Assert.Multiple(() =>
        {
            Assert.That(result1, Is.EqualTo(Let7a_Seed));
            Assert.That(result2, Is.EqualTo(result1));
            Assert.That(result3, Is.EqualTo(result1));
        });
    }

    #endregion
}
