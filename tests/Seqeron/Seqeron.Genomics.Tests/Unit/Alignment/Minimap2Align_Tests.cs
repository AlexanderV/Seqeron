namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A3-9c2, B09 F187): the base-level alignment of minimap2 2.26 (<c>minimap2 -c</c> / <c>-a</c>, <c>align.c</c>): <see cref="Minimap2Program.Run"/> against the installed program, the WHOLE standard output of
/// 2 700 command lines compared as text (TestData/Alignment/minimap2_align_oracle.jsonl.gz, scripts/alignment/generate_minimap2_align_oracle.py): PAF with <c>cg</c> / <c>cs</c> / <c>MD</c> / <c>--eqx</c> and the tags <c>NM ms AS nn ts de zd</c>, SAM with
/// flags, MAPQ, <c>SA</c>, clipping and the header. Groups: the six presets (150 cases each: related queries with substitutions, short and long indels, N bases, both strands, repeats, chimeras), random alignment options (300), SAM (300), the short-read mode (300),
/// spliced alignment (300), repeat families with secondary hits (300) and Z-drop / inversion rescue (300).
/// </summary>
[TestFixture]
public class Minimap2Align_Tests
{
    private static readonly Lazy<List<Minimap2Oracle.Case>> Oracle = new(() => Minimap2Oracle.Read("minimap2_align_oracle.jsonl.gz"));

    private static void Check(string set, int cases)
    {
        var all = Oracle.Value.Where(c => c.Set == set).ToList();
        Assert.That(all, Has.Count.EqualTo(cases), set);
        var failures = Minimap2Oracle.Failures(all);
        Assert.That(failures, Is.Empty, $"{failures.Count} of {all.Count}:\n" + string.Join("\n", failures.Take(3)));
    }

    [Test] public void PafMapOnt_EqualsMinimap2() => Check("paf_map_ont", 150);
    [Test] public void PafMapHifi_EqualsMinimap2() => Check("paf_map_hifi", 150);
    [Test] public void PafMapPb_HomopolymerCompressed_EqualsMinimap2() => Check("paf_map_pb", 150);
    [Test] public void PafAsm5_EqualsMinimap2() => Check("paf_asm5", 150);
    [Test] public void PafAsm10_EqualsMinimap2() => Check("paf_asm10", 150);
    [Test] public void PafAsm20_EqualsMinimap2() => Check("paf_asm20", 150);
    [Test] public void PafAlignmentOptions_EqualsMinimap2() => Check("paf_options", 300);
    [Test] public void Sam_EqualsMinimap2() => Check("sam", 300);
    [Test] public void ShortReads_UngappedAlignment_EqualsMinimap2() => Check("sr", 300);
    [Test] public void Splice_EqualsMinimap2() => Check("splice", 300);
    [Test] public void SecondaryAndSupplementaryHits_EqualMinimap2() => Check("secondary", 300);
    [Test] public void ZDropAndInversionRescue_EqualMinimap2() => Check("inversion", 300);

    [Test]
    public void OracleContainsTheFeaturesItClaims()
    {
        string all = string.Concat(Oracle.Value.Select(c => c.StdOut));
        Assert.Multiple(() =>
        {
            Assert.That(all, Does.Contain("\ttp:A:I\t"), "inversion hit");
            Assert.That(all, Does.Contain("\tzd:i:"), "Z-drop split");
            Assert.That(all, Does.Contain("\tSA:Z:"), "supplementary alignments");
            Assert.That(all, Does.Contain("\tts:A:+"), "transcript strand");
            Assert.That(all, Does.Contain("\tcs:Z:"), "cs tag");
            Assert.That(all, Does.Contain("\tMD:Z:"), "MD tag");
            Assert.That(all, Does.Contain("~"), "intron in cs");
            Assert.That(Oracle.Value.Count(c => c.Args.Contains("-a")), Is.GreaterThan(500), "SAM cases");
        });
    }
}
