using System.Text;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class DesignPrimersTests
{
    // Standard template from Seqeron.Genomics.Tests PrimerDesigner_PrimerDesign_Tests:
    // 100 bp GAACTCGT-unit forward region + 50 bp poly-T target + 100 bp TCCGAAGT-unit
    // reverse region = 258 bp. Yields a valid primer pair for target [100,150) (half-open).
    private static string StandardTemplate()
    {
        var sb = new StringBuilder();
        while (sb.Length < 100) sb.Append("GAACTCGT");
        sb.Append(new string('T', 50));
        int revStart = sb.Length;
        while (sb.Length - revStart < 100) sb.Append("TCCGAAGT");
        return sb.ToString();
    }

    [Test]
    public void DesignPrimers_Schema_ValidatesCorrectly()
    {
        var t = StandardTemplate();
        Assert.DoesNotThrow(() => MolToolsTools.design_primers(t, 100, 150));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers("", 100, 150));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(null!, 100, 150));
        // Negative start.
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, -1, 150));
        // End at/after template length.
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 100, t.Length));
        // start >= end.
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 150, 100));
    }

    [Test]
    public void DesignPrimers_Binding_InvokesSuccessfully()
    {
        var t = StandardTemplate();
        Assert.That(t.Length, Is.EqualTo(258));

        var r = MolToolsTools.design_primers(t, 100, 150);

        Assert.Multiple(() =>
        {
            // Documented behavior of PrimerDesigner.DesignPrimers on this fixture.
            Assert.That(r.IsValid, Is.True);
            Assert.That(r.Forward, Is.Not.Null);
            Assert.That(r.Reverse, Is.Not.Null);

            // Forward is upstream of the target and ends at/before target start.
            Assert.That(r.Forward!.Position, Is.LessThan(100));
            Assert.That(r.Forward.Position + r.Forward.Length, Is.LessThanOrEqualTo(100));
            // Reverse is downstream of target end.
            Assert.That(r.Reverse!.Position, Is.GreaterThanOrEqualTo(150));

            // Exact selected pair = primer3-py 2.3.1 design_primers on this template (sizes 18/20/25,
            // Tm 57/60/63, GC 40-60, poly-X 4, pair ΔTm ≤ 5, PRIMER_PRODUCT_SIZE_RANGE 100-300 (default),
            // SEQUENCE_TARGET 100,50): PRIMER_LEFT_0 = [76,20] TCGTGAACTCGTGAACTCGT, PRIMER_RIGHT_0 = [181,20]
            // CGGAACTTCGGAACTTCGGA, PRIMER_PAIR_0_PENALTY = 0.6983187292252637, PRIMER_PAIR_0_PRODUCT_SIZE 106,
            // PRIMER_PAIR_0_PRODUCT_TM 74.01606258306629.
            Assert.That(r.Forward.Position, Is.EqualTo(76));
            Assert.That(r.Forward.Length, Is.EqualTo(20));
            Assert.That(r.Forward.Sequence, Is.EqualTo("TCGTGAACTCGTGAACTCGT"));
            Assert.That(r.Reverse.Position, Is.EqualTo(162));
            Assert.That(r.Reverse.Length, Is.EqualTo(20));
            Assert.That(r.Reverse.Sequence, Is.EqualTo("CGGAACTTCGGAACTTCGGA"));
            Assert.That(r.PairPenalty!.Value, Is.EqualTo(0.6983187292252637).Within(1e-9));
            Assert.That(r.ProductTm!.Value, Is.EqualTo(74.01606258306629).Within(1e-9));

            // Product size = reverse.Position + reverse.Length - forward.Position = 162 + 20 - 76
            // (= Primer3 PRIMER_PAIR_0_PRODUCT_SIZE 181 - 76 + 1).
            Assert.That(r.ProductSize, Is.EqualTo(106));
            Assert.That(r.ProductSize,
                Is.EqualTo(r.Reverse.Position + r.Reverse.Length - r.Forward.Position));
            Assert.That(r.Pairs, Has.Count.EqualTo(1));

            // Both primers in the documented GC (40-60%) and Tm (57-63°C) windows.
            Assert.That(r.Forward.GcContent, Is.InRange(40.0, 60.0));
            Assert.That(r.Reverse.GcContent, Is.InRange(40.0, 60.0));
            Assert.That(r.Forward.MeltingTemperature, Is.InRange(57.0, 63.0));
            Assert.That(r.Reverse.MeltingTemperature, Is.InRange(57.0, 63.0));
            // Pair Tm difference within 5°C (compatibility rule).
            Assert.That(Math.Abs(r.Forward.MeltingTemperature - r.Reverse.MeltingTemperature),
                Is.LessThanOrEqualTo(5.0));
        });
    }

    [Test]
    public void DesignPrimers_NumReturnAndProductSizeRange_MatchPrimer3()
    {
        // primer3-py 2.3.1 (same settings, PRIMER_NUM_RETURN 3, PRIMER_PRODUCT_SIZE_RANGE "150-250 100-149"):
        // PRIMER_LEFT_k = [76,20]; PRIMER_RIGHT_k = [229,20], [230,20], [237,20]; product sizes 154/155/162.
        var r = MolToolsTools.design_primers(StandardTemplate(), 100, 150,
            product_size_range: "150-250 100-149", num_return: 3);

        Assert.Multiple(() =>
        {
            Assert.That(r.IsValid, Is.True);
            Assert.That(r.Pairs, Has.Count.EqualTo(3));
            Assert.That(r.Pairs.Select(p => p.Reverse!.Position + p.Reverse.Length - 1), Is.EqualTo(new[] { 229, 230, 237 }));
            Assert.That(r.Pairs.Select(p => p.ProductSize), Is.EqualTo(new[] { 154, 155, 162 }));
            Assert.That(r.Pairs.All(p => p.Forward!.Position == 76), Is.True);
            Assert.That(r.ProductSize, Is.EqualTo(154));
            Assert.That(r.ProductTm!.Value, Is.EqualTo(78.6737445237671).Within(1e-9));
        });
    }

    [Test]
    public void DesignPrimers_PickInternalOligo_MatchesPrimer3()
    {
        // primer3-py 2.3.1 (same settings, PRIMER_PICK_INTERNAL_OLIGO 1): PRIMER_LEFT_0 = [76,20],
        // PRIMER_RIGHT_0 = [197,20], PRIMER_INTERNAL_0 = [154,24] TCCGAAGTTCCGAAGTTCCGAAGT,
        // PRIMER_INTERNAL_0_PENALTY 5.875422099275511 (pairs without an internal oligo between the primers fail:
        // "no internal oligo 10").
        var r = MolToolsTools.design_primers(StandardTemplate(), 100, 150, pick_internal_oligo: true);

        Assert.Multiple(() =>
        {
            Assert.That(r.IsValid, Is.True);
            Assert.That(r.Reverse!.Position + r.Reverse.Length - 1, Is.EqualTo(197));
            Assert.That(r.InternalOligo, Is.Not.Null);
            Assert.That(r.InternalOligo!.Value.Start, Is.EqualTo(154));
            Assert.That(r.InternalOligo.Value.Length, Is.EqualTo(24));
            Assert.That(r.InternalOligo.Value.Sequence, Is.EqualTo("TCCGAAGTTCCGAAGTTCCGAAGT"));
            Assert.That(r.InternalOligo.Value.Penalty, Is.EqualTo(5.875422099275511).Within(1e-9));
        });
    }

    [Test]
    public void DesignPrimers_NewOptions_Validated()
    {
        var t = StandardTemplate();
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 100, 150, product_size_range: "300-100"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 100, 150, product_size_range: "abc"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 100, 150, product_size_range: " "));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 100, 150, num_return: 0));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 100, 150, max_tm_difference: -1));

        // A 0 °C Tm-difference limit leaves only equal-Tm pairs: none here.
        var strict = MolToolsTools.design_primers(t, 100, 150, max_tm_difference: 0);
        Assert.That(strict.IsValid, Is.False);
        Assert.That(strict.Pairs, Is.Empty);
    }

    [Test]
    public void DesignPrimers_Primer3AlignmentScreen_MatchesPrimer3()
    {
        // primer3-py 2.3.1 design_primers(SEQUENCE_TARGET [47, 17], PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT 0):
        // PRIMER_LEFT_0 (26, 20), PRIMER_RIGHT_0 (172, 20), PENALTY 0.353112963359024, COMPL_ANY 4, COMPL_END 2,
        // PRIMER_LEFT_0_SELF_ANY 5, _SELF_END 3.
        const string t = "GATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTCGATGTCAAACCCCGGGGGGAGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATTGGTGACGAAAGGTTGTAAGTAGCT";
        var p = Seqeron.Genomics.MolTools.PrimerDesigner.Primer3DefaultParameters with
        {
            StructureScreen = Seqeron.Genomics.MolTools.PrimerStructureScreen.Primer3Alignment,
        };
        var r = MolToolsTools.design_primers(t, 47, 64, p, max_tm_difference: 100);

        Assert.Multiple(() =>
        {
            Assert.That(r.Forward!.Position, Is.EqualTo(26));
            Assert.That(r.Reverse!.Position + r.Reverse.Length - 1, Is.EqualTo(172));
            Assert.That(r.PairPenalty!.Value, Is.EqualTo(0.353112963359024).Within(1e-9));
            Assert.That(r.ComplAny, Is.EqualTo(4.0));
            Assert.That(r.ComplEnd, Is.EqualTo(2.0));
            Assert.That(r.ComplAnyTh, Is.Null);
            Assert.That(r.Forward.SelfAny, Is.EqualTo(5.0));
            Assert.That(r.Forward.SelfEnd, Is.EqualTo(3.0));
        });
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 47, 64, p, pair_max_compl_any: -1));
    }
    [Test]
    public void DesignPrimers_ReactionConditions_MatchPrimer3()
    {
        // primer3-py 2.3.1 design_primers(SEQUENCE_TARGET [100, 30], PRIMER_PAIR_MAX_DIFF_TM 100, PRIMER_SALT_MONOVALENT 20,
        // PRIMER_SALT_DIVALENT 3.0, PRIMER_DNTP_CONC 0.8, PRIMER_DNA_CONC 250, PRIMER_OPT_GC_PERCENT 45,
        // PRIMER_WT_GC_PERCENT_GT 0.5, _LT 1.0, PRIMER_PICK_INTERNAL_OLIGO 1, PRIMER_INTERNAL_SALT_MONOVALENT 120,
        // _SALT_DIVALENT 2, _DNTP_CONC 0.2, _DNA_CONC 100, PRIMER_INTERNAL_OPT_GC_PERCENT 40, _WT_GC_PERCENT_GT/_LT 0.5):
        // PRIMER_LEFT_0 [74,20], PRIMER_RIGHT_0 [278,20], PRIMER_PAIR_0_PENALTY 0.7052466341623926,
        // PRODUCT_TM 87.49740375825182, COMPL_END_TH 17.360051093275388, PRIMER_INTERNAL_0 [122,20] PENALTY 2.1790501883629076.
        const string t = "GATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTCGATGTCAAACCCCGGGGGGAGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATTGGTGACGAAAGGTTGTAAGTAGCT";
        var r = MolToolsTools.design_primers(t, 100, 130, Seqeron.Genomics.MolTools.PrimerDesigner.Primer3DefaultParameters,
            max_tm_difference: 100, num_return: 5, pick_internal_oligo: true,
            salt_monovalent: 20, salt_divalent: 3.0, dntp_conc: 0.8, dna_conc: 250,
            opt_gc_percent: 45, wt_gc_percent_gt: 0.5, wt_gc_percent_lt: 1.0,
            internal_salt_monovalent: 120, internal_salt_divalent: 2.0, internal_dntp_conc: 0.2, internal_dna_conc: 100,
            internal_opt_gc_percent: 40, internal_wt_gc_percent_gt: 0.5, internal_wt_gc_percent_lt: 0.5);

        Assert.Multiple(() =>
        {
            Assert.That(r.IsValid, Is.True);
            Assert.That(r.Pairs, Has.Count.EqualTo(5));
            Assert.That(r.Forward!.Position, Is.EqualTo(74));
            Assert.That(r.Reverse!.Position + r.Reverse.Length - 1, Is.EqualTo(278));
            Assert.That(r.PairPenalty!.Value, Is.EqualTo(0.7052466341623926).Within(1e-9));
            Assert.That(r.ProductTm!.Value, Is.EqualTo(87.49740375825182).Within(1e-9));
            Assert.That(r.ComplEndTh!.Value, Is.EqualTo(17.360051093275388).Within(1e-9));
            Assert.That(r.InternalOligo!.Value.Start, Is.EqualTo(122));
            Assert.That(r.InternalOligo.Value.Penalty, Is.EqualTo(2.1790501883629076).Within(1e-9));
            Assert.That(r.Pairs[4].PairPenalty!.Value, Is.EqualTo(1.1878536407214142).Within(1e-9));
        });

        // Primer3 _pr_data_control: salt / DNA concentration must be > 0, divalent / dNTP ≥ 0.
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.design_primers(t, 100, 130, salt_monovalent: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.design_primers(t, 100, 130, dna_conc: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.design_primers(t, 100, 130, salt_divalent: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            MolToolsTools.design_primers(t, 100, 130, pick_internal_oligo: true, internal_dna_conc: 0));
    }

    [Test]
    public void DesignPrimers_ThreePrimeEndChecksAndDistance_MatchPrimer3()
    {
        // primer3-py 2.3.1 design_primers(SEQUENCE_TARGET [54, 21], PRIMER_PAIR_MAX_DIFF_TM 100):
        // PRIMER_GC_CLAMP 2 → PRIMER_LEFT_0 [34,18], PRIMER_RIGHT_0 [143,21], PRIMER_PAIR_0_PENALTY 3.8176034305392363;
        // PRIMER_MAX_END_GC 2 → [19,19] / [136,20], 2.2130093738348933; PRIMER_MAX_END_STABILITY 3.0 → no pair;
        // PRIMER_MIN_THREE_PRIME_DISTANCE 3 → PRIMER_PAIR_1 [35,19] / [140,20], 2.570253724637439;
        // PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE 5 + _RIGHT_ 0 → PRIMER_PAIR_2 [16,18] / [140,20], 4.073043969173341.
        const string t = "GGTCCTAATTGGAGCGCCCAGTTACCGGCCGAGTGCTACGGGCACTCGTTGGTAGTGGGCTCCCTAAGTCGGCGCATCCGTTCCTAGCTTTAAAATATCCGTTGAAAGAATGTTCTGAGTCTCGCCTAGTGAAAGCCAACTCCTTTGGATTTTGTCATA";
        var p3 = Seqeron.Genomics.MolTools.PrimerDesigner.Primer3DefaultParameters;
        var clamp = MolToolsTools.design_primers(t, 54, 75, p3, max_tm_difference: 100, num_return: 5, gc_clamp: 2);
        var endGc = MolToolsTools.design_primers(t, 54, 75, p3, max_tm_difference: 100, num_return: 5, max_end_gc: 2);
        var stab = MolToolsTools.design_primers(t, 54, 75, p3, max_tm_difference: 100, num_return: 5, max_end_stability: 3.0);
        var alias = MolToolsTools.design_primers(t, 54, 75, p3, max_tm_difference: 100, num_return: 5, min_three_prime_distance: 3);
        var lr = MolToolsTools.design_primers(t, 54, 75, p3, max_tm_difference: 100, num_return: 5,
            min_left_three_prime_distance: 5, min_right_three_prime_distance: 0);

        Assert.Multiple(() =>
        {
            Assert.That(clamp.Forward!.Position, Is.EqualTo(34));
            Assert.That(clamp.Reverse!.Position + clamp.Reverse.Length - 1, Is.EqualTo(143));
            Assert.That(clamp.PairPenalty!.Value, Is.EqualTo(3.8176034305392363).Within(1e-9));
            Assert.That(endGc.Forward!.Position, Is.EqualTo(19));
            Assert.That(endGc.PairPenalty!.Value, Is.EqualTo(2.2130093738348933).Within(1e-9));
            Assert.That(stab.IsValid, Is.False);
            Assert.That(stab.Pairs, Is.Empty);
            Assert.That(alias.Pairs, Has.Count.EqualTo(5));
            Assert.That(alias.Pairs[1].Forward!.Position, Is.EqualTo(35));
            Assert.That(alias.Pairs[1].Reverse!.Position + alias.Pairs[1].Reverse!.Length - 1, Is.EqualTo(140));
            Assert.That(alias.Pairs[1].PairPenalty!.Value, Is.EqualTo(2.570253724637439).Within(1e-9));
            Assert.That(lr.Pairs[2].Forward!.Position, Is.EqualTo(16));
            Assert.That(lr.Pairs[2].PairPenalty!.Value, Is.EqualTo(4.073043969173341).Within(1e-9));
        });

        // read_boulder: "Both PRIMER_MIN_THREE_PRIME_DISTANCE and PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE specified";
        // _pr_data_control: PRIMER_MAX_END_GC ∉ [0, 5], PRIMER_GC_CLAMP > PRIMER_MIN_SIZE, distance < −1.
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 54, 75, min_three_prime_distance: 3, min_left_three_prime_distance: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.design_primers(t, 54, 75, max_end_gc: 6));
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.design_primers(t, 54, 75, p3, gc_clamp: 19));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 54, 75, min_three_prime_distance: -2));
    }

    [Test]
    public void DesignPrimers_TemplateMispriming_MatchesPrimer3()
    {
        // primer3-py 2.3.1 design_primers(SEQUENCE_TARGET [100, 20]) on a repetitive template:
        // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT 0, PRIMER_MAX_TEMPLATE_MISPRIMING 9, PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING 16,
        // PRIMER_WT_TEMPLATE_MISPRIMING 0.5, PRIMER_PAIR_WT_TEMPLATE_MISPRIMING 0.2 → PRIMER_LEFT_0 [59,21], PRIMER_RIGHT_0
        // [237,20], PRIMER_PAIR_0_PENALTY 10.14045393659626, LEFT/RIGHT/PAIR_0_TEMPLATE_MISPRIMING 6 / 4 / 10;
        // PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT 1, PRIMER_MAX_TEMPLATE_MISPRIMING_TH 45, PRIMER_WT_TEMPLATE_MISPRIMING_TH 0.5,
        // PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH 0, PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH 0.2 → [58,20] / [208,20],
        // PRIMER_PAIR_0_PENALTY 1.1607262752978076, PRIMER_PAIR_0_TEMPLATE_MISPRIMING_TH 29.881724966610534.
        const string t = "AGACTTTCAATGTCAAACCAATCTACCCCCTTCCGTATTATTTGTTACCAATTCTCATTGTGTTTCGGAACTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCGATCCGTTCCTAATAAGGAATGGTGATTCCCTGTCATACCAATCTACCCCCTTCCGAAACACAATGAGAATTGGTAACAAACAGCGCAGCGGCAGATCAAGCAGGAGGCGGAATGTAAACA";
        var p3 = Seqeron.Genomics.MolTools.PrimerDesigner.Primer3DefaultParameters;
        var align = MolToolsTools.design_primers(t, 100, 120,
            p3 with { StructureScreen = Seqeron.Genomics.MolTools.PrimerStructureScreen.Primer3Alignment },
            max_tm_difference: 100, num_return: 5, max_template_mispriming: 9, pair_max_template_mispriming: 16,
            wt_template_mispriming: 0.5, pair_wt_template_mispriming: 0.2);
        var thermo = MolToolsTools.design_primers(t, 100, 120, p3, max_tm_difference: 100, num_return: 5,
            thermodynamic_template_alignment: true, max_template_mispriming_th: 45, wt_template_mispriming_th: 0.5,
            pair_max_template_mispriming_th: 0, pair_wt_template_mispriming_th: 0.2);
        var none = MolToolsTools.design_primers(t, 100, 120, p3, max_tm_difference: 100);
        Assert.Multiple(() =>
        {
            Assert.That(none.Forward!.Position, Is.EqualTo(16));
            Assert.That(none.TemplateMispriming, Is.Null);
            Assert.That(align.Forward!.Position, Is.EqualTo(59));
            Assert.That(align.Forward.Length, Is.EqualTo(21));
            Assert.That(align.Reverse!.Position + align.Reverse.Length - 1, Is.EqualTo(237));
            Assert.That(align.PairPenalty!.Value, Is.EqualTo(10.14045393659626).Within(1e-9));
            Assert.That(align.Forward.TemplateMispriming, Is.EqualTo(6.0));
            Assert.That(align.Reverse.TemplateMispriming, Is.EqualTo(4.0));
            Assert.That(align.TemplateMispriming, Is.EqualTo(10.0));
            Assert.That(align.Pairs, Has.Count.EqualTo(5));
            Assert.That(thermo.Forward!.Position, Is.EqualTo(58));
            Assert.That(thermo.Reverse!.Position + thermo.Reverse.Length - 1, Is.EqualTo(208));
            Assert.That(thermo.PairPenalty!.Value, Is.EqualTo(1.1607262752978076).Within(1e-9));
            Assert.That(thermo.TemplateMispriming!.Value, Is.EqualTo(29.881724966610534).Within(1e-9));
        });
    }

    [Test]
    public void DesignPrimers_BoundAndPositionPenalty_MatchPrimer3()
    {
        // primer3-py 2.3.1 design_primers(SEQUENCE_TARGET [140, 20]) on random.Random(5) 300 nt:
        // PRIMER_ANNEALING_TEMP 60, PRIMER_MIN_BOUND 50, PRIMER_OPT_BOUND 90, PRIMER_WT_BOUND_LT/GT 0.1 → PRIMER_LEFT_0
        // [116,20], RIGHT_0 [233,20], PAIR_0_PENALTY 8.061085955165812, LEFT/RIGHT_0_BOUND 53.45996470116724 / 50.93138362747912;
        // PRIMER_INSIDE_PENALTY 2, PRIMER_OUTSIDE_PENALTY 0.1 → [108,20] / [211,20], 5.6992044449417225, POSITION_PENALTY
        // 1.2000000000000002 / 3.2; PRIMER_PICK_INTERNAL_OLIGO + PRIMER_ANNEALING_TEMP 60 + PRIMER_INTERNAL_MIN_BOUND 60 →
        // PRIMER_INTERNAL_0 [87,20], PENALTY 1.1533159279721872, BOUND 60.29348152377712.
        const string t = "GGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGCTTGTCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCTAATCGACTATGCTACTGCGGTTAACG";
        var p3 = Seqeron.Genomics.MolTools.PrimerDesigner.Primer3DefaultParameters;
        var bound = MolToolsTools.design_primers(t, 140, 160, p3, max_tm_difference: 100, num_return: 5,
            annealing_temp: 60, min_bound: 50, opt_bound: 90, wt_bound_lt: 0.1, wt_bound_gt: 0.1);
        var pos = MolToolsTools.design_primers(t, 140, 160, p3, max_tm_difference: 100, inside_penalty: 2, outside_penalty: 0.1);
        var intl = MolToolsTools.design_primers(t, 140, 160, p3, max_tm_difference: 100, pick_internal_oligo: true,
            annealing_temp: 60, internal_min_bound: 60);
        Assert.Multiple(() =>
        {
            Assert.That(bound.Forward!.Position, Is.EqualTo(116));
            Assert.That(bound.Reverse!.Position + bound.Reverse.Length - 1, Is.EqualTo(233));
            Assert.That(bound.PairPenalty!.Value, Is.EqualTo(8.061085955165812).Within(1e-9));
            Assert.That(bound.Forward.Bound!.Value, Is.EqualTo(53.45996470116724).Within(1e-9));
            Assert.That(bound.Reverse.Bound!.Value, Is.EqualTo(50.93138362747912).Within(1e-9));
            Assert.That(bound.Pairs, Has.Count.EqualTo(5));
            Assert.That(pos.Forward!.Position, Is.EqualTo(108));
            Assert.That(pos.Reverse!.Position + pos.Reverse.Length - 1, Is.EqualTo(211));
            Assert.That(pos.PairPenalty!.Value, Is.EqualTo(5.6992044449417225).Within(1e-9));
            Assert.That(pos.Forward.PositionPenalty!.Value, Is.EqualTo(1.2).Within(1e-9));
            Assert.That(pos.Reverse.PositionPenalty!.Value, Is.EqualTo(3.2).Within(1e-9));
            Assert.That(intl.InternalOligo!.Value.Start, Is.EqualTo(87));
            Assert.That(intl.InternalOligo.Value.Penalty, Is.EqualTo(1.1533159279721872).Within(1e-9));
            Assert.That(intl.InternalOligo.Value.Bound!.Value, Is.EqualTo(60.29348152377712).Within(1e-9));
            // _pr_data_control errors.
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.design_primers(t, 140, 160, p3, opt_bound: 120));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.design_primers(t, 140, 160, p3, pick_internal_oligo: true, internal_opt_bound: -20));
            // A negative pair penalty (default inside −1 with outside 0.05): primer3-py aborts (obj_fn PR_ASSERT).
            Assert.Throws<InvalidOperationException>(() => MolToolsTools.design_primers(t, 140, 160, p3, max_tm_difference: 100, outside_penalty: 0.05));
        });
    }

    [Test]
    public void DesignPrimers_SequenceQuality_MatchesPrimer3()
    {
        // primer3-py 2.3.1 design_primers(SEQUENCE_TARGET [140, 20], SEQUENCE_QUALITY q[i] = 12 if i % 61 == 5, 28 if
        // i % 23 == 7, else 35 + i % 6) on random.Random(5) 300 nt: PRIMER_MIN_QUALITY 20, PRIMER_MIN_END_QUALITY 30,
        // PRIMER_WT_SEQ_QUAL 0.1 → PRIMER_LEFT_0 [46,20], RIGHT_0 [236,20], PAIR_0_PENALTY 14.273993465730893,
        // LEFT/RIGHT_0_MIN_SEQ_QUALITY 28 / 35; PRIMER_PICK_INTERNAL_OLIGO + PRIMER_INTERNAL_MIN_QUALITY 25 +
        // PRIMER_INTERNAL_WT_SEQ_QUAL 0.2 + PRIMER_PAIR_WT_IO_PENALTY 1 + PRIMER_QUALITY_RANGE_MAX 60 + PRIMER_WT_SEQ_QUAL 0.05 →
        // PAIR_0_PENALTY 9.380853713532826, PRIMER_INTERNAL_0 [31,20], PENALTY 5.3889380835424845, MIN_SEQ_QUALITY 35.
        const string t = "GGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGCTTGTCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCTAATCGACTATGCTACTGCGGTTAACG";
        int[] q = Enumerable.Range(0, t.Length).Select(i => i % 61 == 5 ? 12 : i % 23 == 7 ? 28 : 35 + i % 6).ToArray();
        var p3 = Seqeron.Genomics.MolTools.PrimerDesigner.Primer3DefaultParameters;
        var r = MolToolsTools.design_primers(t, 140, 160, p3, max_tm_difference: 100, num_return: 5,
            sequence_quality: q, min_quality: 20, min_end_quality: 30, wt_seq_qual: 0.1, wt_end_qual: 7);
        var intl = MolToolsTools.design_primers(t, 140, 160, p3, max_tm_difference: 100, pick_internal_oligo: true,
            sequence_quality: q, quality_range_max: 60, wt_seq_qual: 0.05, internal_min_quality: 25, internal_wt_seq_qual: 0.2,
            pair_wt_io_penalty: 1);
        Assert.Multiple(() =>
        {
            Assert.That(r.Forward!.Position, Is.EqualTo(46));
            Assert.That(r.Reverse!.Position + r.Reverse.Length - 1, Is.EqualTo(236));
            Assert.That(r.PairPenalty!.Value, Is.EqualTo(14.273993465730893).Within(1e-9));
            Assert.That(r.Forward.MinSequenceQuality, Is.EqualTo(28));
            Assert.That(r.Reverse.MinSequenceQuality, Is.EqualTo(35));
            Assert.That(r.Pairs, Has.Count.EqualTo(5));
            Assert.That(intl.PairPenalty!.Value, Is.EqualTo(9.380853713532826).Within(1e-9));
            Assert.That(intl.InternalOligo!.Value.Start, Is.EqualTo(31));
            Assert.That(intl.InternalOligo.Value.Penalty, Is.EqualTo(5.3889380835424845).Within(1e-9));
            Assert.That(intl.InternalOligo.Value.MinSequenceQuality, Is.EqualTo(35));
            // _pr_data_control errors.
            Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 140, 160, p3, min_quality: 20));
            Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 140, 160, p3, wt_seq_qual: 1));
            Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 140, 160, p3, sequence_quality: q[..^1]));
            Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 140, 160, p3, pair_wt_io_penalty: 0.5));
        });
    }

    [Test]
    public void DesignPrimers_MaskTemplate_MatchesPrimer3()
    {
        // primer3-py 2.3.1 design_primers(SEQUENCE_TARGET [140, 20], PRIMER_MASK_TEMPLATE 1, PRIMER_MASK_KMERLIST_PATH with
        // GenomeTester4 lists of the counts below, PRIMER_WT_MASK_FAILURE_RATE 1) on random.Random(5) 300 nt: rank 3 =
        // [117,20] / [233,20], PAIR_PENALTY 0.5401089974063126, LEFT_PENALTY 0.43329038433534206; with
        // PRIMER_MASK_FAILURE_RATE 0.05, PRIMER_MASK_5P_DIRECTION 2, PRIMER_MASK_3P_DIRECTION 1 (no weight): rank 0 = [117,20] /
        // [233,20], PAIR_PENALTY 0.5003559657240544.
        const string t = "GGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGCTTGTCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCTAATCGACTATGCTACTGCGGTTAACG";
        var k11 = new Dictionary<string, int>();
        for (int i = 0; i <= t.Length - 11; i++)
            if (i % 3 == 0)
            {
                string k = t.Substring(i, 11);
                k11[i % 2 == 1 ? Seqeron.Genomics.Core.DnaSequence.GetReverseComplementString(k) : k] = 2 + i * 37 % 1000;
            }
        var k16 = new Dictionary<string, int>();
        for (int i = 0; i <= t.Length - 16; i++)
            if (i % 7 == 0)
                k16[t.Substring(i, 16)] = (int)Math.Pow(10, 1 + i % 5);
        var p3 = Seqeron.Genomics.MolTools.PrimerDesigner.Primer3DefaultParameters;
        var r = MolToolsTools.design_primers(t, 140, 160, p3, max_tm_difference: 100, num_return: 5,
            mask_template: true, mask_kmers_11: k11, mask_kmers_16: k16, wt_mask_failure_rate: 1);
        var r2 = MolToolsTools.design_primers(t, 140, 160, p3, max_tm_difference: 100,
            mask_template: true, mask_kmers_11: k11, mask_kmers_16: k16, mask_failure_rate: 0.05, mask_5p_direction: 2, mask_3p_direction: 1);
        Assert.Multiple(() =>
        {
            Assert.That(r.Pairs, Has.Count.EqualTo(5));
            Assert.That(r.Pairs[3].Forward!.Position, Is.EqualTo(117));
            Assert.That(r.Pairs[3].PairPenalty!.Value, Is.EqualTo(0.5401089974063126).Within(1e-9));
            Assert.That(r.Pairs[3].Forward!.Penalty, Is.EqualTo(0.43329038433534206).Within(1e-9));
            Assert.That(r.Pairs[3].Forward!.MaskFailureRate, Is.Not.Null);
            Assert.That(r2.Forward!.Position, Is.EqualTo(117));
            Assert.That(r2.PairPenalty!.Value, Is.EqualTo(0.5003559657240544).Within(1e-9));
            Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 140, 160, p3, mask_template: true));
        });
    }

    [Test]
    public void DesignPrimers_MisprimingLibrary_MatchesPrimer3()
    {
        // primer3-py 2.3.1 design_primers(SEQUENCE_TARGET [100, 30], PRIMER_PAIR_MAX_DIFF_TM 100, misprime_lib = lib):
        // default PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS 0 → PRIMER_LEFT_0 [26,20] (the library-free [68,20] scores 20 > 12),
        // PRIMER_RIGHT_0 [180,20], PRIMER_PAIR_0_PENALTY 0.42552744112833807, PRIMER_LEFT_0_LIBRARY_MISPRIMING
        // (4.0, reverse site1), PRIMER_RIGHT_0_ (10.0, site2*0.5), PRIMER_PAIR_0_ (11.0, site2*0.5);
        // PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS 1 → PRIMER_LEFT_0_LIBRARY_MISPRIMING (9.0, iupac), PRIMER_PAIR_0_ (17.0, iupac).
        const string t = "GATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTCGATGTCAAACCCCGGGGGGAGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATTGGTGACGAAAGGTTGTAAGTAGCT";
        var lib = new Dictionary<string, string>
        {
            ["site1"] = "TGTATAGTCCCACCTGGTGATCCTATGCTTGTGAG",
            ["site2*0.5"] = "gccagaaggctgcaactcatcgactctatg",
            ["iupac"] = "NNRYCCAGAAAATAGCGWSKMBDHV",
        };
        var p3 = Seqeron.Genomics.MolTools.PrimerDesigner.Primer3DefaultParameters;
        var r = MolToolsTools.design_primers(t, 100, 130, p3, max_tm_difference: 100, num_return: 5, mispriming_library: lib);
        var cons = MolToolsTools.design_primers(t, 100, 130, p3, max_tm_difference: 100, num_return: 5, mispriming_library: lib,
            lib_ambiguity_codes_consensus: true);
        var noLib = MolToolsTools.design_primers(t, 100, 130, p3, max_tm_difference: 100);

        Assert.Multiple(() =>
        {
            Assert.That(noLib.Forward!.Position, Is.EqualTo(68));
            Assert.That(noLib.LibraryMispriming, Is.Null);
            Assert.That(r.Forward!.Position, Is.EqualTo(26));
            Assert.That(r.Reverse!.Position + r.Reverse.Length - 1, Is.EqualTo(180));
            Assert.That(r.PairPenalty!.Value, Is.EqualTo(0.42552744112833807).Within(1e-9));
            Assert.That(r.Forward.LibraryMispriming, Is.EqualTo(4.0));
            Assert.That(r.Forward.LibraryMisprimingName, Is.EqualTo("reverse site1"));
            Assert.That(r.Reverse.LibraryMispriming, Is.EqualTo(10.0));
            Assert.That(r.Reverse.LibraryMisprimingName, Is.EqualTo("site2*0.5"));
            Assert.That(r.LibraryMispriming, Is.EqualTo(11.0));
            Assert.That(r.LibraryMisprimingName, Is.EqualTo("site2*0.5"));
            Assert.That(r.Pairs[4].LibraryMispriming, Is.EqualTo(7.0));
            Assert.That(cons.Forward!.LibraryMispriming, Is.EqualTo(9.0));
            Assert.That(cons.Forward.LibraryMisprimingName, Is.EqualTo("iupac"));
            Assert.That(cons.LibraryMispriming, Is.EqualTo(17.0));
            Assert.That(System.Text.Json.JsonSerializer.Serialize(r), Does.Contain("\"LibraryMisprimingName\":\"site2*0.5\""));
        });

        // _pr_data_control: a library weight without a library; add_seq_to_seq_lib: illegal weight / empty sequence.
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 100, 130, p3, wt_library_mispriming: 1));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 100, 130, p3, pair_wt_library_mispriming: 1));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 100, 130, p3,
            mispriming_library: new Dictionary<string, string> { ["x*101"] = "ACGT" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 100, 130, p3,
            mispriming_library: new Dictionary<string, string> { ["x"] = "" }));
    }

    [Test]
    public void DesignPrimers_InternalOligoMishybLibrary_MatchesPrimer3()
    {
        // primer3-py 2.3.1 design_primers(SEQUENCE_TARGET [100, 30], PRIMER_PAIR_MAX_DIFF_TM 100, PRIMER_PICK_INTERNAL_OLIGO 1,
        // PRIMER_INTERNAL_MAX_LIBRARY_MISHYB 9, mishyb_lib = lib): PRIMER_INTERNAL_0 [108,19] (the library-free [108,20]
        // scores 10 > 9), PRIMER_INTERNAL_0_PENALTY 1.2417931164713423, PRIMER_INTERNAL_0_LIBRARY_MISPRIMING (9.0, s1);
        // PRIMER_INTERNAL_4 [194,21] (5.0, s1).
        const string t = "GATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTCGATGTCAAACCCCGGGGGGAGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATTGGTGACGAAAGGTTGTAAGTAGCT";
        var lib = new Dictionary<string, string>
        {
            ["s1"] = "GCGGTGTTAAGTGTCGAGCTACATCACTTCTC",
            ["s2*0.5"] = "atgtagccagaaggctgcaactcatcgactctatg",
            ["iu"] = "GGTGTTAAGTGTCRAGCTACAYC",
        };
        var p3 = Seqeron.Genomics.MolTools.PrimerDesigner.Primer3DefaultParameters;
        var r = MolToolsTools.design_primers(t, 100, 130, p3, max_tm_difference: 100, num_return: 5, pick_internal_oligo: true,
            internal_mishyb_library: lib, internal_max_library_mishyb: 9);
        var noLib = MolToolsTools.design_primers(t, 100, 130, p3, max_tm_difference: 100, pick_internal_oligo: true);

        Assert.Multiple(() =>
        {
            Assert.That((noLib.InternalOligo!.Value.Start, noLib.InternalOligo.Value.Length), Is.EqualTo((108, 20)));
            Assert.That(noLib.InternalOligo.Value.LibraryMishyb, Is.Null);
            var io = r.InternalOligo!.Value;
            Assert.That((io.Start, io.Length), Is.EqualTo((108, 19)));
            Assert.That(io.Penalty, Is.EqualTo(1.2417931164713423).Within(1e-9));
            Assert.That(io.LibraryMishyb, Is.EqualTo(9.0));
            Assert.That(io.LibraryMishybName, Is.EqualTo("s1"));
            var io4 = r.Pairs[4].InternalOligo!.Value;
            Assert.That((io4.Start, io4.Length, io4.LibraryMishyb), Is.EqualTo((194, 21, (double?)5.0)));
        });

        // _pr_data_control: an internal-oligo mishyb weight without a mishyb library.
        Assert.Throws<ArgumentException>(() => MolToolsTools.design_primers(t, 100, 130, p3, pick_internal_oligo: true,
            internal_wt_library_mishyb: 1));
    }

    [Test]
    public void DesignPrimers_LowercaseMaskingAndGcOptimum_MatchPrimer3()
    {
        // primer3-py 2.3.1 design_primers(SEQUENCE_TARGET [100,20], PRIMER_LOWERCASE_MASKING 1, PRIMER_PICK_INTERNAL_OLIGO 1)
        // on a mixed-case template (B07 F46): PRIMER_PAIR_1 [14,20] / [232,20] 0.936221987820943 (without masking the
        // second pair is [14,20] / [188,20] 0.9201203739298762 — its right primer ends on a lower-case base).
        const string t = "AGACTTTCAAagatatgctgggtagaggtcGAGGTTATTAtTTGTTAcCAATtCTCATTGTGTTTCGGAActtgCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCGatCCGTTCcTAaTAAGGAATGGTGATTCCCtgtcataccaatctaccccctgttaTGCGCGTTTGTCGTTaGACCAaTGtCAGCGcAGCGGCAGATCAAGCAGGAGGCGGAATGTAAACA";
        var p3 = Seqeron.Genomics.MolTools.PrimerDesigner.Primer3DefaultParameters;
        var on = MolToolsTools.design_primers(t, 100, 120, p3, max_tm_difference: 100, num_return: 5, pick_internal_oligo: true,
            lowercase_masking: true);
        var off = MolToolsTools.design_primers(t, 100, 120, p3, max_tm_difference: 100, num_return: 5, pick_internal_oligo: true);
        Assert.Multiple(() =>
        {
            Assert.That(on.Pairs.Select(x => (x.Forward!.Position, x.Reverse!.Position + x.Reverse.Length - 1)),
                Is.EqualTo(new[] { (14, 235), (14, 232), (16, 235), (14, 205), (16, 232) }));
            Assert.That(on.Pairs[1].PairPenalty!.Value, Is.EqualTo(0.936221987820943).Within(1e-9));
            Assert.That(on.Pairs[3].InternalOligo!.Value.Start, Is.EqualTo(104));
            Assert.That(off.Pairs[1].Reverse!.Position + off.Pairs[1].Reverse!.Length - 1, Is.EqualTo(188));
            Assert.That(off.Pairs[1].PairPenalty!.Value, Is.EqualTo(0.9201203739298762).Within(1e-9));
        });

        // Primer3 _pr_data_control: a GC weight without PRIMER_[INTERNAL_]OPT_GC_PERCENT (undefined by default).
        Assert.That(() => MolToolsTools.design_primers(t, 100, 120, p3, wt_gc_percent_gt: 0.5),
            Throws.ArgumentException.With.Message.StartsWith("Primer GC content is part of objective function while optimum gc_content is not defined"));
        Assert.That(() => MolToolsTools.design_primers(t, 100, 120, p3, internal_wt_gc_percent_lt: 0.5),
            Throws.ArgumentException.With.Message.StartsWith("Hyb probe GC content is part of objective function while optimum gc_content is not defined"));
    }
}
