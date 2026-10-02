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
}
