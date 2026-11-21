namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A5-9b / A5-9c, B09 F188): the multi-segment (paired-read) mode of minimap2 2.26: <see cref="Minimap2Program.Run"/> with two query files or an interleaved file (collect_minimizers with segment ids, the multi-segment step score of
/// <c>comput_sc</c>, <c>mm_select_sub_multi</c>, <c>mm_seg_gen</c>, per-segment <c>mm_set_parent</c> / mapq, <c>pe_ori</c> flipping, <c>mm_pair</c> under <c>-c</c>) against the installed program, the whole standard output of 2 100 command lines
/// (TestData/Alignment/minimap2_pair_oracle.jsonl.gz, scripts/alignment/generate_minimap2_pair_oracle.py): FR pairs with inserts of 100 - 1500 bp and errors, unmapped / wrongly oriented / distant mates, overlapping mates, adapters, N bases, repeats (secondary
/// hits), one to three fragments per file, FASTA and FASTQ, <c>/1</c> <c>/2</c> suffixes, three reads with one name (<c>FI</c> tags).
/// </summary>
[TestFixture]
public class Minimap2Pair_Tests
{
    private static readonly Lazy<List<Minimap2Oracle.Case>> Oracle = new(() => Minimap2Oracle.Read("minimap2_pair_oracle.jsonl.gz"));

    private static void Check(string set, int cases)
    {
        var all = Oracle.Value.Where(c => c.Set == set).ToList();
        Assert.That(all, Has.Count.EqualTo(cases), set);
        var failures = Minimap2Oracle.Failures(all);
        Assert.That(failures, Is.Empty, $"{failures.Count} of {all.Count}:\n" + string.Join("\n", failures.Take(3)));
    }

    [Test] public void TwoFiles_Paf_EqualsMinimap2() => Check("pe_files_paf", 300);
    [Test] public void TwoFiles_PafWithCigar_Pairing_EqualsMinimap2() => Check("pe_files_pafc", 300);
    [Test] public void TwoFiles_Sam_EqualsMinimap2() => Check("pe_files_sam", 300);
    [Test] public void Interleaved_PafWithCigar_EqualsMinimap2() => Check("pe_inter_paf", 300);
    [Test] public void Interleaved_Sam_EqualsMinimap2() => Check("pe_inter_sam", 300);
    [Test] public void FragmentMode_WithoutPreset_ForwardForwardOrientation_EqualsMinimap2() => Check("pe_ff", 300);
    [Test] public void Options_EqualMinimap2() => Check("pe_options", 300);

    [Test]
    public void OracleContainsProperPairsMateFieldsAndSupplementaryHits()
    {
        string sam = string.Concat(Oracle.Value.Where(c => c.Set is "pe_files_sam" or "pe_inter_sam").Select(c => c.StdOut));
        var flags = sam.Split('\n').Where(l => l.Length > 0 && l[0] != '@').Select(l => int.Parse(l.Split('\t')[1])).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(flags.Count(f => (f & 0x2) != 0), Is.GreaterThan(100), "proper pairs");
            Assert.That(flags.Count(f => (f & 0x8) != 0), Is.GreaterThan(50), "mate unmapped");
            Assert.That(flags.Count(f => (f & 0x20) != 0), Is.GreaterThan(50), "mate reversed");
            Assert.That(flags.Count(f => (f & 0x40) != 0), Is.GreaterThan(100), "first in pair");
            Assert.That(sam, Does.Contain("FI:i:"), "three reads with one name");
        });
    }
}
