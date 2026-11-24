using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-s5, B09 WP74, F342): <c>mafft --seed file</c> / <c>--seedtable file</c> against the installed Debian mafft 7.505-1 (<c>scripts/alignment/generate_mafft_wp74_oracle.py</c> seed, seed-multi, seedtable, seed-err):
/// return code, standard output, the files and the whole standard error of every run.
/// <list type="bullet">
/// <item><c>multi2hat3s</c> writes the seed sequences (<c>&gt;_seed_name</c>, with their gaps) in front of the input and the ungapped blocks of every pair of the seed alignment as lines of kind 'k' (weight <c>nhomologs^2 * 100</c>, numbers shifted by the
/// seed offset); <c>disttbfast -E 1</c> aligns everything without the seeds' columns, <c>tbfast -J</c> (the later passes of the script: distances from the alignment, the structure weights <c>fastconjuction_noname_kozo</c> / <c>effarr_kozo</c>,
/// no reuse of profiles) and <c>dvtditr -l</c> use the blocks (<c>readlocalhomtable2_half</c>: <c>kozoarivec</c>); with a pair library the blocks are added to the chains of the pair stage.</item>
/// <item><c>seed-multi</c> - two or three seed files (the offsets, the numbers of the files count down); <c>--anysymbol</c> - <c>replaceu -o seedoffset</c> of every seed and of the input, <c>restoreu</c> reads the seeds and the input from <c>orig</c>;
/// <c>seed-err</c> - a single sequence, unaligned seeds, both options, a seed table with a sequence number outside the input, a sequence in the table twice, an unknown letter, <c>--anysymbol</c>, a pair library, CRLF.</item>
/// </list>
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftSeed_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Cases = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_wp74_seed_oracle.tsv"));

    public static IEnumerable<TestCaseData> Groups() => Cases.Value.Select(c => c.Group).Distinct().Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    [Test]
    public void Oracle_HasEveryGroup()
    {
        Assert.Multiple(() =>
        {
            foreach (var g in Cases.Value.GroupBy(c => c.Group))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            Assert.That(Cases.Value.Select(c => c.Group).Distinct().Order(), Is.EqualTo(new[] { "seed", "seed-err", "seed-large", "seed-multi", "seed-opt", "seedtable" }));
            var seed = Cases.Value.Where(c => c.Group == "seed").ToList();
            Assert.That(seed.Count(c => c.Rc == 0 && c.Out.Contains(">_seed_", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(150), "the seeds are part of the alignment");
            Assert.That(seed.Count(c => c.Args.Contains("--anysymbol")), Is.GreaterThanOrEqualTo(20), "replaceu of the seed file and the input, restoreu");
            Assert.That(seed.Count(c => c.Args.Any(a => a is "--localpair" or "--globalpair" or "--genafpair")), Is.GreaterThanOrEqualTo(15), "the seed blocks are added to the pair library");
            Assert.That(seed.Count(c => c.Args.Contains("--maxiterate")), Is.GreaterThanOrEqualTo(20), "dvtditr -l with the seed blocks");
            Assert.That(seed.Count(c => c.OutFiles.Keys.Any(k => k.EndsWith(".tree", StringComparison.Ordinal))), Is.GreaterThanOrEqualTo(5), "--treeout");
            Assert.That(Cases.Value.Count(c => c.Err.Contains("Progressive alignment method is incompatible with the --seed option.", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(10), "the warning of the script");
            Assert.That(Cases.Value.Count(c => c.Err.Contains("multi2hat3s (aa) Version 7.505", StringComparison.Ordinal) || c.Err.Contains("multi2hat3s (nuc) Version 7.505", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(5), "the standard error of multi2hat3s");
            var multi = Cases.Value.Where(c => c.Group == "seed-multi").ToList();
            Assert.That(multi.Count(c => c.Args.Count(a => a == "--seed") >= 2 && c.Rc == 0), Is.GreaterThanOrEqualTo(180), "two or three seed files");
            var opt = Cases.Value.Where(c => c.Group == "seed-opt").ToList();
            Assert.That(opt.Count(c => c.Rc == 0 && c.Out.Contains(">_seed_", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(100), "the seed with the other options of the script");
            Assert.That(opt.Count(c => c.Args.Contains("--treein")), Is.GreaterThanOrEqualTo(30), "a guide tree file with the seeds");
            Assert.That(opt.Count(c => c.Args.Contains("--merge")), Is.GreaterThanOrEqualTo(30), "subalignments with the seeds (the stopped disttbfast, then tbfast -J on the empty file)");
            Assert.That(opt.Count(c => c.Args.Contains("--add") || c.Args.Contains("--addprofile")), Is.GreaterThanOrEqualTo(10), "--add with --seed: the script says 'Use either ONE'");
            var large = Cases.Value.Where(c => c.Group == "seed-large").ToList();
            Assert.That(large.Count(c => c.Rc == 0 && c.Out.Contains(">_seed_", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(170), "nodepair and tbfast -U with the seeds");
            Assert.That(large.Count(c => c.Args.Count(a => a == "--seed") >= 2), Is.GreaterThanOrEqualTo(60), "two or three seed files");
            Assert.That(large.Count(c => c.Args.Contains("--large")), Is.GreaterThanOrEqualTo(50), "--large");
            Assert.That(large.Count(c => c.Args.Contains("--memsavetree")), Is.GreaterThanOrEqualTo(50), "--memsavetree");
            Assert.That(large.Count(c => c.Args.Contains("--stepadd")), Is.GreaterThanOrEqualTo(50), "--stepadd");
            Assert.That(large.Count(c => c.Args.Contains("--localpair")), Is.GreaterThanOrEqualTo(40), "L-INS");
            Assert.That(large.Count(c => c.Args.Contains("--globalpair")), Is.GreaterThanOrEqualTo(40), "G-INS");
            Assert.That(large.Count(c => c.Args.Contains("--genafpair") || c.Args.Contains("--oldgenafpair")), Is.GreaterThanOrEqualTo(80), "E-INS");
            Assert.That(large.Count(c => c.Args.Contains("--anysymbol")), Is.GreaterThanOrEqualTo(20), "replaceu / restoreu of the seeds");
            var table = Cases.Value.Where(c => c.Group == "seedtable").ToList();
            Assert.That(table.Count(c => c.Rc == 0), Is.GreaterThanOrEqualTo(200), "--seedtable");
            Assert.That(table.Count(c => c.Files.Values.Any(v => v.Contains(" h\n", StringComparison.Ordinal))), Is.GreaterThanOrEqualTo(20), "blocks of kind 'h' in a seed table");
            var err = Cases.Value.Where(c => c.Group == "seed-err").ToList();
            Assert.That(err.Count(c => c.Rc != 0), Is.GreaterThanOrEqualTo(60), "the stops");
            Assert.That(err.Count(c => c.Err.Contains("At least 2 sequences should be input!", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(2), "a seed file with one sequence");
            Assert.That(err.Count(c => c.Err.Contains("Use either one of seedtable and seed.  Not both.", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(5), "both options");
            Assert.That(err.Count(c => c.Rc != 0 && c.Args.Contains("--seedtable")), Is.GreaterThanOrEqualTo(10), "a seed table line with j <= i or a sequence number outside the input: tbfast stops with 'Check hat3.'");
        });
    }

    [Test]
    public void Multi2Hat3s_WritesTheSeedsAndTheBlocksOfEveryPair()
    {
        // installed multi2hat3s 7.505-1: multi2hat3s -t 5 -o 2 -i seed1
        const string seed = ">s1\nMKV-LAAGIVG\n>s2\nMKVALAAGIVG\n>s3\nMRVLL-AGIVG\n";
        var r = SequenceAligner.MafftRun.Multi2Hat3s(seed, 5, 2, "seed1");
        Assert.That(r.Failure, Is.Null);
        Assert.That(r.Rows, Is.EqualTo(">_seed_s1\nMKV-LAAGIVG\n>_seed_s2\nMKVALAAGIVG\n>_seed_s3\nMRVLL-AGIVG\n"));
        Assert.That(r.Hat3, Is.EqualTo("2 3 10 14632.917 0 2 0 2 k\n2 3 10 14632.917 3 9 4 10 k\n2 4 9 14132.130 0 2 0 2 k\n2 4 9 14132.130 3 3 4 4 k\n2 4 9 14132.130 5 9 5 9 k\n3 4 10 12786.583 0 4 0 4 k\n3 4 10 12786.583 6 10 5 9 k\n"));
        Assert.That(r.Log, Is.EqualTo("nhomologs = 5\nseedoffset = 2\nseed = seed1\nrescale = 1\nGap Penalty = -1.53, +0.00, -0.12\ntsuyosa = 2500.000000\n\nadding 0-1\radding 0-2\radding 1-2\r\nmulti2hat3s (aa) Version 7.505\nalg=A, model=BLOSUM62, 1.53, +0.12, -0.00, noshift, amax=0.0\n1 thread(s)\n\n"));
        Assert.That(SequenceAligner.MafftRun.Multi2Hat3s(">only\nACGT\n", 5, 0, "seed1").Failure!.Text, Is.EqualTo("At least 2 sequences should be input!\nOnly 1 sequence found.\n"));
    }

    [Test]
    public void Seed_KeepsTheSeedsAndAlignsTheInputToThem()
    {
        // installed mafft 7.505-1: mafft --quiet --seed s1.fa in.fa
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">x1\nMKVLAGIVGLLAX\n>x2\nMRVLAAGIVGLLL\n",
            ["s1.fa"] = ">s1\nMKV-LAAGIVG\n>s2\nMKVALAAGIVG\n>s3\nMRVLL-AGIVG\n",
        };
        var r = MafftProgram.Run(["--quiet", "--seed", "s1.fa", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.Zero);
        Assert.That(r.Stdout, Is.EqualTo(">_seed_s1\nMKV-LAAGIVG----\n>_seed_s2\nMKVALAAGIVG----\n>_seed_s3\nMRVLL-AGIVG----\n>x1\nMKV-L-AGIVGLLAX\n>x2\nMRV-LAAGIVGLLL-\n"));
        Assert.That(r.Stderr, Does.StartWith("############################################################################\n# Warning:\n#   Progressive alignment method is incompatible with the --seed option.\n"));
    }

    [Test]
    public void Clustal_ColumnOfCAndU_GetsNoMark()
    {
        // f2cl.c: weaker[2] = "CU" but nweaker = 2 - the group is never tested (found by the oracle case seed-0146: a seed with 'U'); installed mafft 7.505-1: mafft --quiet --clustalout
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">a\nACGCACGCAC\n>b\nACGUACGUAC\n>c\nACGCACGCAC\n",
        };
        var r = MafftProgram.Run(["--quiet", "--clustalout", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.Zero);
        Assert.That(r.Stdout, Is.EqualTo("CLUSTAL format alignment by MAFFT FFT-NS-2 (v7.505)\n\n\na               acgcacgcac\nb               acguacguac\nc               acgcacgcac\n                *** *** **\n"));
    }

    [TestCase("--parttree")]
    [TestCase("--dpparttree")]
    public void SeedWithPartTree_IsImpossible(string option)
    {
        // installed mafft 7.505-1: the script stops before any program ("Impossible", return code 1)
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">x1\nMKVLAGIVGLLAX\n>x2\nMRVLAAGIVGLLL\n>x3\nMRVLAAGIVGLLL\n",
            ["s1.fa"] = ">s1\nMKV-LAAGIVG\n>s2\nMKVALAAGIVG\n>s3\nMRVLL-AGIVG\n",
        };
        var r = MafftProgram.Run(["--quiet", "--seed", "s1.fa", option, SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.EqualTo(1));
        Assert.That(r.Stdout, Is.Empty);
        Assert.That(r.Stderr, Is.EqualTo("Impossible\n"));
    }

    [Test]
    public void SeedWithMemSave_StopsInTbfast()
    {
        // installed mafft 7.505-1: tbfast -J -M with the blocks of the seeds (constraint 2): "Memory saving mode is not supported." and exit status 1
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">x1\nMKVLAGIVGLLAX\n>x2\nMRVLAAGIVGLLL\n>x3\nMRVLAAGIVGLLI\n",
            ["s1.fa"] = ">s1\nMKV-LAAGIVG\n>s2\nMKVALAAGIVG\n>s3\nMRVLL-AGIVG\n",
        };
        var r = MafftProgram.Run(["--seed", "s1.fa", "--memsave", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.EqualTo(1));
        Assert.That(r.Stdout, Is.Empty);
        Assert.That(r.Stderr, Does.EndWith("Progressive alignment ... \n\rSTEP     1 /5 \n\nMemory saving mode is not supported.\n\n"));
    }

    [Test]
    public void TenMerPairOnProtein_StopsWithTheMessageOfDisttbfast()
    {
        // installed mafft 7.505-1: disttbfast "tuplesize must be 6 for aa sequence", exit status 1 (an ordinary stop of the program, not an unported refusal)
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">x1\nMKVLAGIVGLLAX\n>x2\nMRVLAAGIVGLLL\n>x3\nMRVLAAGIVGLLI\n",
        };
        var r = MafftProgram.Run(["--10merpair", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.EqualTo(1));
        Assert.That(r.Stdout, Is.Empty);
        Assert.That(r.Stderr, Does.EndWith("Gap Penalty = -1.53, +0.00, +0.00\n\ntuplesize must be 6 for aa sequence\n"));
    }

    [Test]
    public void SeedTable_WithAShortLine_IsRefused()
    {
        // valgrind on the installed tbfast with hat3 = "x y z": "Conditional jump or move depends on uninitialised value(s)" - readlocalhomtable2_half reads i, j and the kind of a line that sscanf did not fill
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">a\nMKVLAAGIVG\n>b\nMKVLAGIVG\n>c\nMKVIAAGLVG\n",
            ["t.txt"] = "x y z\n",
        };
        Assert.Throws<NotSupportedException>(() => MafftProgram.Run(["--quiet", "--seedtable", "t.txt", SequenceAligner_MafftProgram_Tests.InputName], files));
    }

    [TestCaseSource(nameof(Groups))]
    public void Group_EqualsMafft(string group)
    {
        var failures = new List<string>();
        int n = 0;
        foreach (var c in Cases.Value.Where(x => x.Group == group))
        {
            n++;
            if (SequenceAligner_MafftProgram_Tests.Verify(c) is { } failure)
                failures.Add(failure);
        }
        Assert.That(n, Is.GreaterThan(0));
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }
}
