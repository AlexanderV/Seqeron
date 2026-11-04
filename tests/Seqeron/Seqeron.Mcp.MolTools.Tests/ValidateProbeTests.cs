using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class ValidateProbeTests
{
    // From Seqeron.Genomics.Tests PROBE-VALID-001.
    private const string UniqueProbe = "ATCGATCGATCGATCGATCG";
    private static readonly string[] SingleMatchReference = { "NNNNNATCGATCGATCGATCGATCGNNNN" };

    [Test]
    public void ValidateProbe_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.validate_probe(UniqueProbe, SingleMatchReference));
        Assert.Throws<ArgumentException>(() => MolToolsTools.validate_probe(null!, SingleMatchReference));
        Assert.Throws<ArgumentException>(() => MolToolsTools.validate_probe(UniqueProbe, null!));
        Assert.Throws<ArgumentException>(() => MolToolsTools.validate_probe(UniqueProbe, SingleMatchReference, max_mismatches: -1));
    }

    [Test]
    public void ValidateProbe_UniqueProbe_SpecificityOne()
    {
        // 1 exact hit -> specificity 1.0.
        var v = MolToolsTools.validate_probe(UniqueProbe, SingleMatchReference, max_mismatches: 0);
        Assert.Multiple(() =>
        {
            Assert.That(v.OffTargetHits, Is.EqualTo(1));
            Assert.That(v.SpecificityScore, Is.EqualTo(1.0).Within(1e-9));
        });
    }

    [Test]
    public void ValidateProbe_MultipleHits_SpecificityIsOneOverCount()
    {
        // 10-mer poly-A in a 34-mer poly-A: 34-10+1 = 25 exact-match positions -> specificity 1/25.
        var v = MolToolsTools.validate_probe("AAAAAAAAAA",
            new[] { "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" }, max_mismatches: 0);
        Assert.Multiple(() =>
        {
            Assert.That(v.OffTargetHits, Is.EqualTo(25));
            Assert.That(v.SpecificityScore, Is.EqualTo(1.0 / 25).Within(1e-9));
        });
    }

    [Test]
    public void ValidateProbe_EmptyProbe_InvalidZeroSpecificity()
    {
        // Empty probe is a degenerate input: invalid, zero specificity, no hits.
        var v = MolToolsTools.validate_probe("", SingleMatchReference);
        Assert.Multiple(() =>
        {
            Assert.That(v.IsValid, Is.False);
            Assert.That(v.SpecificityScore, Is.EqualTo(0.0).Within(1e-9));
            Assert.That(v.OffTargetHits, Is.EqualTo(0));
        });
    }

    [Test]
    public void ValidateProbe_ThermodynamicScreen_ReportsPrimer3NtthalTm()
    {
        // primer3-py 2.3.1 calc_homodimer / calc_hairpin at mv 50, dv 0, dntp 0, dna 50 (Primer3 probe conditions).
        var v = MolToolsTools.validate_probe("GCGCGCGCGCGCGCGCGCGC", Array.Empty<string>());
        Assert.Multiple(() =>
        {
            Assert.That(v.ThermodynamicScreen, Is.True);
            Assert.That(v.SelfDimerTm!.Value, Is.EqualTo(78.85652531616256).Within(1e-9));
            Assert.That(v.HairpinTm!.Value, Is.EqualTo(87.30265612393043).Within(1e-9));
            Assert.That(v.IsValid, Is.False);
        });
    }

    [Test]
    public void ValidateProbe_NonTargets_KaneCriteria()
    {
        // Kane et al. (2000): identity > 75 % or > 15 contiguous identical bases. Non-target = probe with a
        // substitution at every 5th position (Biopython local alignment: 40/50 identical → 80 %).
        const string probe = "TATGCCTCCGGTACATCAACTACAGTTAGCCTTAAGAGAAAAATCCCAAA";
        const string nonTarget = "CCGCACCATGAGACTGTTTCTATGGCTCCTGTACCTCAAGTACATTTAGGCTTACGAGACAAATGCCAACCACATCGGCTTCGCACGTCT";
        var v = MolToolsTools.validate_probe(probe, new[] { probe }, non_target_sequences: new[] { nonTarget });
        Assert.Multiple(() =>
        {
            Assert.That(v.CrossHybridization, Has.Count.EqualTo(2));
            Assert.That(v.CrossHybridization[0].Identity, Is.EqualTo(0.8).Within(1e-12));
            Assert.That(v.CrossHybridization[0].CrossHybridizes, Is.True);
            Assert.That(v.IsValid, Is.False);
        });
        Assert.Throws<ArgumentException>(() => MolToolsTools.validate_probe(probe, new[] { probe }, max_non_target_identity: 1.5));
        Assert.Throws<ArgumentException>(() => MolToolsTools.validate_probe(probe, new[] { probe }, max_contiguous_match: -1));
    }

    [Test]
    public void ValidateProbe_LongProbe_ReportsPrimer3AlignmentSelfAny()
    {
        // > 60 nt: fallback self-dimer criterion = Primer3 alignment-mode self_any/self_end (dpal.c: (ACGT)16 → 64.00
        // > PRIMER_INTERNAL_MAX_SELF_ANY 12.00).
        var v = MolToolsTools.validate_probe(string.Concat(Enumerable.Repeat("ACGT", 16)), Array.Empty<string>());
        Assert.Multiple(() =>
        {
            Assert.That(v.ThermodynamicScreen, Is.False);
            Assert.That(v.SelfAny, Is.EqualTo(64.0));
            Assert.That(v.SelfEnd, Is.EqualTo(64.0));
            Assert.That(v.Issues, Has.Some.EqualTo("Self-complementarity: Primer3 self_any 64.00 exceeds 12.00"));
        });
    }

    [Test]
    public void ValidateProbe_LongProbe_FallbackStemFlagIsWarning()
    {
        // Audit round 5, A5-5: random 61-mer (thal.c -DTHAL_MAX_ALIGN=10000: hairpin 33.98 °C, ANY / END1 below 0 °C)
        // carries exact 4-bp stems → the fallback stem-loop flag is reported as a warning but isValid stays true;
        // with thermodynamic_screen_max_length = 61 the sourced ntthal hairpin Tm decides (33.98 °C ≤ 47 °C).
        const string random61 = "AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGT";
        var v = MolToolsTools.validate_probe(random61, new[] { random61 });
        var thermo = MolToolsTools.validate_probe(random61, new[] { random61 }, thermodynamic_screen_max_length: 61);
        Assert.Multiple(() =>
        {
            Assert.That(v.ThermodynamicScreen, Is.False);
            Assert.That(v.HasSecondaryStructure, Is.True);
            Assert.That(v.Warnings, Is.EqualTo(new[] { "Potential secondary structure formation" }));
            Assert.That(v.Issues, Is.Empty);
            Assert.That(v.IsValid, Is.True);
            Assert.That(thermo.HairpinTm!.Value, Is.EqualTo(33.980529935122263).Within(1e-6));
            Assert.That(thermo.IsValid, Is.True);
        });
    }

    [Test]
    public void ValidateProbe_ReactionConditions_DelegateToNtthalScreen()
    {
        // primer3-py 2.3.1 calc_homodimer / calc_end_stability / calc_hairpin of GCGC…(20 nt) at mv 100, dv 2, dntp 0.2,
        // dna 250 (default conditions mv 50 / dv 0 / dntp 0 / dna 50: 78.85652531616256 / 87.30265612393043).
        var v = MolToolsTools.validate_probe("GCGCGCGCGCGCGCGCGCGC", Array.Empty<string>(),
            monovalent_mm: 100, divalent_mm: 2, dntp_mm: 0.2, dna_conc_nm: 250);
        Assert.Multiple(() =>
        {
            Assert.That(v.ThermodynamicScreen, Is.True);
            Assert.That(v.SelfDimerTm!.Value, Is.EqualTo(88.57933031369095).Within(1e-9));
            Assert.That(v.SelfEndDimerTm!.Value, Is.EqualTo(88.57933031369095).Within(1e-9));
            Assert.That(v.HairpinTm!.Value, Is.EqualTo(93.4845818217363).Within(1e-9));
            Assert.That(v.IsValid, Is.False);
        });
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.validate_probe(UniqueProbe, SingleMatchReference, monovalent_mm: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.validate_probe(UniqueProbe, SingleMatchReference, dna_conc_nm: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.validate_probe(UniqueProbe, SingleMatchReference, divalent_mm: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.validate_probe(UniqueProbe, SingleMatchReference, dntp_mm: -1));
    }

    [Test]
    public void ValidateProbe_ScreenLimits_DelegateToLibrary()
    {
        // max_structure_tm above the ntthal values (78.86 / 78.86 / 87.30 °C) → no structure issue.
        var lenient = MolToolsTools.validate_probe("GCGCGCGCGCGCGCGCGCGC", Array.Empty<string>(), max_structure_tm: 90);
        // thermodynamic_screen = false → the fallback screens (Primer3 alignment-mode self_any / self_end).
        var fallback = MolToolsTools.validate_probe("GCGCGCGCGCGCGCGCGCGC", Array.Empty<string>(), thermodynamic_screen: false);
        // (ACGT)16: dpal self_any = self_end = 64.00 (dpal.c) → issue at the default 12, none at 64.
        string acgt16 = string.Concat(Enumerable.Repeat("ACGT", 16));
        var limit64 = MolToolsTools.validate_probe(acgt16, Array.Empty<string>(), max_self_any: 64, max_self_end: 64);
        Assert.Multiple(() =>
        {
            Assert.That(lenient.HairpinTm!.Value, Is.EqualTo(87.30265612393043).Within(1e-9));
            Assert.That(lenient.HasSecondaryStructure, Is.False);
            Assert.That(lenient.IsValid, Is.True);
            Assert.That(fallback.ThermodynamicScreen, Is.False);
            Assert.That(fallback.SelfDimerTm, Is.Null);
            Assert.That(limit64.SelfAny, Is.EqualTo(64.0));
            Assert.That(limit64.Issues, Has.None.StartWith("Self-complementarity"));
        });
    }

    [Test]
    public void ValidateProbe_ReferenceSites_JudgedByKaneCriteria()
    {
        // 12-mer: exact site + 3-mismatch site (9/12 = 0.75, not > 0.75) → no off-target issue; + 2-mismatch site
        // (10/12 > 0.75) → 2 Kane sites → issue (Python brute-force Hamming oracle).
        const string probe = "GATCCGACGCTA";
        var one = MolToolsTools.validate_probe(probe, new[] { "TTTTTGATCCGACGCTATTTTTGCTCCTACGGTATTTTT" });
        var two = MolToolsTools.validate_probe(probe, new[] { "TTTTTGATCCGACGCTATTTTTGCTCCTACGGTATTTTT", "GGGGGGAACCGAGGCTAGGGGG" });
        Assert.Multiple(() =>
        {
            Assert.That((one.OffTargetHits, one.CrossHybridizingHits), Is.EqualTo((2, 1)));
            Assert.That(one.Issues, Has.None.Contain("off-target"));
            Assert.That((two.OffTargetHits, two.CrossHybridizingHits), Is.EqualTo((3, 2)));
            Assert.That(two.Issues, Has.Some.StartWith("2 potential off-target sites"));
        });
    }

    [Test]
    public void ValidateProbe_BothStrands_DelegatesReverseComplementScan()
    {
        // Reference = exact site + rc(probe with 2 substitutions) + rc(probe with 5 substitutions) (random.seed(20261008));
        // independent Python brute force over both strands: radius 3 → (1, 1) given strand, (2, 2) both strands;
        // radius 5 → (3, 2) both strands (the 5-substitution site is 15/20 = 0.75, longest run 3).
        const string probe = "GATCCGACGCTATATGCCGT";
        string[] refs = { "ACAGTTTTAAGATAGGATCCGACGCTATATGCCGTAGCGAAAGCGCAGACACGGCATAGAGCGTCGCATCAATAAATAATCCGTAACCGCAGATACCGTAGGAGCGGGAGACCTGGCACA" };
        var single = MolToolsTools.validate_probe(probe, refs);
        var both = MolToolsTools.validate_probe(probe, refs, both_strands: true);
        var both5 = MolToolsTools.validate_probe(probe, refs, max_mismatches: 5, both_strands: true);
        Assert.Multiple(() =>
        {
            Assert.That((single.OffTargetHits, single.CrossHybridizingHits, single.IsValid), Is.EqualTo((1, 1, true)));
            Assert.That((both.OffTargetHits, both.CrossHybridizingHits, both.IsValid), Is.EqualTo((2, 2, false)));
            Assert.That(both.SpecificityScore, Is.EqualTo(0.5).Within(1e-12));
            Assert.That((both5.OffTargetHits, both5.CrossHybridizingHits), Is.EqualTo((3, 2)));
            // Reverse-palindromic probe: one site, not two.
            var pal = MolToolsTools.validate_probe("GAATTCCGGAATTC", new[] { "TTTTTGAATTCCGGAATTCAAAAA" }, both_strands: true);
            Assert.That((pal.OffTargetHits, pal.CrossHybridizingHits), Is.EqualTo((1, 1)));
        });
    }
}
