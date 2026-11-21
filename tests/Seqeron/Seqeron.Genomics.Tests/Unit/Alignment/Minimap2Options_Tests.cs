namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A3-9c2 / A5-9b, B09 F189): the options of the mapping and alignment stages of minimap2 2.26 beyond the alignment itself and the paired-read mode (<see cref="Minimap2Program"/>: <c>main.c</c> options with getopt semantics, <c>-x</c>
/// applied first and cumulatively, <c>mm_check_opt</c>, the reader, the index in parts, <c>--split-prefix</c> merge, <c>--junc-bed</c>, <c>--alt</c>, <c>--qstrand</c>, the output flags), the whole standard output (and for the error group the exit status and the
/// standard error at <c>-v 2</c>) of 2 403 command lines against the installed program (TestData/Alignment/minimap2_options_oracle.jsonl.gz, scripts/alignment/generate_minimap2_options_oracle.py). <c>-t -K -2 --bucket-bits --no-kalloc --cap-kalloc
/// --print-*</c> are in the sets and have no effect on the standard output; <c>-I</c> has one (the index is split into parts of about the given size and every part is mapped on its own).
/// </summary>
[TestFixture]
public class Minimap2Options_Tests
{
    private static readonly Lazy<List<Minimap2Oracle.Case>> Oracle = new(() => Minimap2Oracle.Read("minimap2_options_oracle.jsonl.gz"));

    private static void Check(string set, int cases)
    {
        var all = Oracle.Value.Where(c => c.Set == set).ToList();
        Assert.That(all, Has.Count.EqualTo(cases), set);
        var failures = Minimap2Oracle.Failures(all);
        Assert.That(failures, Is.Empty, $"{failures.Count} of {all.Count}:\n" + string.Join("\n", failures.Take(3)));
    }

    [Test] public void IndexParts_BatchSizeAndSplitPrefixMerge_EqualMinimap2() => Check("opt_index_parts", 300);
    [Test] public void ChainingAndSeedingOptions_WithCigar_EqualMinimap2() => Check("opt_chain", 300);
    [Test] public void OutputOptions_EqualMinimap2() => Check("opt_output", 300);
    [Test] public void AltContigs_EqualMinimap2() => Check("opt_alt", 300);
    [Test] public void QueryStrandMode_EqualsMinimap2() => Check("opt_qstrand", 300);
    [Test] public void AllVersusAll_WithCigar_EqualsMinimap2() => Check("opt_ava", 300);
    [Test] public void UnusualInput_EqualsMinimap2() => Check("opt_input", 300);
    [Test] public void LongCigar_CgTag_EqualsMinimap2() => Check("opt_long_cigar", 3);
    [Test] public void OptionChecksHelpVersionAndErrors_EqualMinimap2() => Check("opt_errors", 300);
}
