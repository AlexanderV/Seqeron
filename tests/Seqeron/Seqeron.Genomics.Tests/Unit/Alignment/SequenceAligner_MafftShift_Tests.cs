using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-s5, B09 WP74, F343): <c>mafft --globalpair --shiftpenalty x</c> against the installed Debian mafft 7.505-1 (<c>scripts/alignment/generate_mafft_wp74_oracle.py shift</c>): return code, standard output, the files and the whole
/// standard error of every run. With x below 10 the programs switch to the "trywarp" mode (<c>constants.c</c>: <c>penalty_shift_factor &lt; 10</c>): <c>A__align</c> (the progressive step, with and without the blocks of the pair library), <c>G__align11</c> (the pairwise
/// alignments of the pair library) and <c>partA__align</c> (the iterative refinement; <c>dvtditr</c> fixes every x below 100 to 2) may jump from the best cell of the part of the matrix above and to the left of the cell to the cell, at the price of the shift penalty
/// (<c>penalty_shift = x * penalty</c> plus the extension penalty of the skipped columns); the trace-back reads the jump targets (<c>warpis</c>, <c>warpjs</c>) from the cells that hold <c>warpbase + n</c>. The script allows the option only with <c>--globalpair</c>
/// (and a pair library weight of 1 or more); every other method stops with its message.
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftShift_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Cases = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_wp74_shift_oracle.tsv"));

    private static double Value(SequenceAligner_MafftProgram_Tests.Case c)
    {
        int i = Array.IndexOf(c.Args, "--shiftpenalty");
        return i >= 0 && i + 1 < c.Args.Length && double.TryParse(c.Args[i + 1], System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out double v) ? v : double.NaN;
    }

    [Test]
    public void Oracle_Shape()
    {
        var all = Cases.Value;
        Assert.Multiple(() =>
        {
            Assert.That(all.Count, Is.GreaterThanOrEqualTo(200));
            var warp = all.Where(c => c.Rc == 0 && Value(c) < 10).ToList();
            Assert.That(warp.Count, Is.GreaterThanOrEqualTo(120), "x below 10: the warp of A__align / G__align11 / partA__align");
            Assert.That(warp.Count(c => c.Args.Contains("--maxiterate")), Is.GreaterThanOrEqualTo(20), "partA__align (dvtditr fixes x to 2)");
            Assert.That(warp.Count(c => c.Args.Contains("--retree")), Is.GreaterThanOrEqualTo(10), "--retree");
            Assert.That(warp.Count(c => c.Args.Contains("--large") || c.Args.Contains("--memsavetree") || c.Args.Contains("--stepadd")), Is.GreaterThanOrEqualTo(5), "nodepair");
            Assert.That(warp.Count(c => c.Args.Contains("--weighti")), Is.GreaterThanOrEqualTo(3), "--weighti");
            Assert.That(all.Count(c => c.Rc == 0 && Value(c) >= 10), Is.GreaterThanOrEqualTo(50), "x of 10 or more: no warp");
            Assert.That(all.Count(c => c.Rc != 0 && c.Err.Contains("supported only with the --globalpair option", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(10), "the other methods");
            Assert.That(all.Count(c => c.Rc != 0 && c.Err.Contains("shiftpenalty must be >1.", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(3), "x below 1");
        });
    }

    [Test]
    public void Oracle_EqualsMafft()
    {
        var failures = new List<string>();
        foreach (var c in Cases.Value)
        {
            if (SequenceAligner_MafftProgram_Tests.Verify(c) is { } failure)
                failures.Add(failure);
        }
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    [Test]
    public void ShiftPenaltyOfOne_ChangesTheAlignment()
    {
        // installed mafft 7.505-1: the same five DNA sequences, --globalpair with and without --shiftpenalty 1.0 (the jump to an earlier cell of the matrix)
        const string input = ">s0\ntaaacaggcttgggtctggtcgaatcctagagccctatgag\n>s1\ncaacagccttgggtctggtctaatcgctacacacatgtgag\n>s2\nTAACGTGACGCTGTTACGTGCCAGCTCATACATTCATCGTATGTATAGTGACCCATTGAGAAGAATCATA\n"
            + ">s3\ncaacaccttggggtctggttaaatactagatcccggtgat\n>s4\nCTTTTTAATCGGGAACCAGTACAATGTCCATCTACATCGCATGCGAGCCCACAAAGAACAATATTTCCAA\n";
        var files = new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = input };
        var shifted = MafftProgram.Run(["--quiet", "--globalpair", "--shiftpenalty", "1.0", SequenceAligner_MafftProgram_Tests.InputName], files);
        var plain = MafftProgram.Run(["--quiet", "--globalpair", "--shiftpenalty", "100", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(shifted.ExitCode, Is.Zero);
        Assert.That(plain.ExitCode, Is.Zero);
        Assert.That(shifted.Stdout, Is.EqualTo(">s0\ntaaacaggcttgggtctggtcgaatc-ctagagccctatga-------------------\n----------------g\n>s1\nc-aacagccttgggtctggtctaatcgctacacacatgtga-------------------\n----------------g\n"
            + ">s2\ntaacgtgacgctgttacgtgccagctcatacattcatcgtatgtatagtgacccattgag\naagaat-------cata\n>s3\nc-aacaccttggggtctggttaaata-ctagatcccggtga-------------------\n----------------t\n"
            + ">s4\nctttttaatcgggaaccagtacaatgtccatctacatcgcatgcg----agcccacaaag\naacaatatttcca---a\n"));
        Assert.That(plain.Stdout, Is.EqualTo(">s0\ntaaacaggcttgggtctggtcgaatc-ctagagccct-----------------------\n---------atgag\n>s1\nc-aacagccttgggtctggtctaatcgctacacacat-----------------------\n---------gtgag\n"
            + ">s2\ntaacgtgacgctgttacgtgccagctcatacattcatcgtatgtatagtgacccattgag\naagaa----tcata\n>s3\nc-aacaccttggggtctggttaaata-ctagatcccg-----------------------\n---------gtgat\n"
            + ">s4\nctttttaatcgggaaccagtacaatgtccatctacatcgcatgcg----agcccacaaag\naacaatatttccaa\n"));
    }

    [Test]
    public void ShiftPenaltyBelowOne_StopsTheScript()
    {
        var files = new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = ">a\nMKVLAAGIVG\n>b\nMKVLAGIVG\n" };
        var r = MafftProgram.Run(["--globalpair", "--shiftpenalty", "0.5", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.EqualTo(1));
        Assert.That(r.Stderr, Is.EqualTo("\nshiftpenalty must be >1.\n"));
    }

    [Test]
    public void ShiftPenaltyWithLocalPair_StopsTheScript()
    {
        var files = new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = ">a\nMKVLAAGIVG\n>b\nMKVLAGIVG\n" };
        var r = MafftProgram.Run(["--localpair", "--shiftpenalty", "2", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.EqualTo(1));
        Assert.That(r.Stderr, Is.EqualTo("\nAt present, --unalignlevel # or --allowshift is supported only with the --globalpair option.\n\n"));
    }
}
