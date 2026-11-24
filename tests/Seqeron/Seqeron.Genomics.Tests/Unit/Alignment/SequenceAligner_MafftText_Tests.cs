using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c-s5, B09 WP74, F341): <c>mafft --text</c> / <c>--textmatrix file</c> (an alphabet of 256 characters: <c>constants.c</c> <c>extendedmtx</c>, <c>overridematrix</c> of <c>blosum.c</c>, <c>calcfreq_extended</c>,
/// <c>charfilter</c>, <c>f2cl -E</c>) against the installed Debian mafft 7.505-1 (<c>scripts/alignment/generate_mafft_wp74_oracle.py</c> text, text-format, text-pair, text-iter, text-matrix, text-err): return code, standard output, the
/// files and the whole standard error of every run.
/// <list type="bullet">
/// <item><c>text</c> / <c>text-format</c> / <c>text-pair</c> / <c>text-iter</c> - texts over alphabets of letters, digits, punctuation, bytes above 127 and control characters; the progressive method, the options of the scoring and trees, the output
/// formats (the marks of CLUSTAL are '*' for a column of identical characters), the pair libraries and the iterative refinement (<c>setcore</c> reads the alignment again in the plain mode: <c>onlyAlpha_upper</c>).</item>
/// <item><c>text-matrix</c> - <c>--textmatrix</c>: the lines "0x41 0x42 3.5" of <c>_aamtx</c> replace entries of the matrix (+1 on the diagonal, -1 elsewhere), the "Score(...)" lines of the program, the stops of <c>overridematrix</c>.</item>
/// <item><c>text-err</c> - '=', '&lt;' and '&gt;' in a sequence, <c>--text</c> with <c>--anysymbol</c>, a missing or malformed matrix file, a letter that cannot be used (NUL, blank, '-', '&lt;', '=', '&gt;'), a negative average of the diagonal.</item>
/// </list>
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftText_Tests
{
    internal static readonly Lazy<List<SequenceAligner_MafftProgram_Tests.Case>> Cases = new(() => SequenceAligner_MafftProgram_Tests.Load("mafft_wp74_text_oracle.tsv"));

    public static IEnumerable<TestCaseData> Groups() => Cases.Value.Select(c => c.Group).Distinct().Select(g => new TestCaseData(g).SetArgDisplayNames(g));

    [Test]
    public void Oracle_HasEveryGroup()
    {
        Assert.Multiple(() =>
        {
            foreach (var g in Cases.Value.GroupBy(c => c.Group))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            Assert.That(Cases.Value.Select(c => c.Group).Distinct().Order(), Is.EqualTo(new[] { "text", "text-err", "text-format", "text-iter", "text-matrix", "text-pair" }));
            // the standard error is stored in full for every 25th case and for the short ones only
            var text = Cases.Value.Where(c => c.Group == "text").ToList();
            Assert.That(Cases.Value.Count(c => c.Rc == 0 && c.Out.Length > 0 && c.Err.Contains("(text) Version", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(30), "the version line of the programs says text");
            Assert.That(text.Count(c => c.Input.Any(ch => ch > 127)), Is.GreaterThanOrEqualTo(30), "bytes above 127");
            Assert.That(text.Count(c => c.Input.Contains('\t', StringComparison.Ordinal) || c.Input.Contains('\u0001', StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(10), "control characters");
            Assert.That(text.Count(c => c.OutFiles.Keys.Any(k => k.EndsWith(".tree", StringComparison.Ordinal))), Is.GreaterThanOrEqualTo(2), "--treeout");
            var format = Cases.Value.Where(c => c.Group == "text-format").ToList();
            Assert.That(format.Count(c => c.Out.StartsWith("CLUSTAL", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(80), "CLUSTAL");
            Assert.That(format.Count(c => c.Args.Contains("--phylipout")), Is.GreaterThanOrEqualTo(30), "PHYLIP");
            var pair = Cases.Value.Where(c => c.Group == "text-pair").ToList();
            Assert.That(pair.Count(c => c.Rc == 0), Is.GreaterThanOrEqualTo(200), "a pair library of the extended alphabet");
            var iter = Cases.Value.Where(c => c.Group == "text-iter").ToList();
            Assert.That(iter.Count(c => c.Args.Contains("--core")), Is.GreaterThanOrEqualTo(5), "setcore reads the alignment as plain protein: the bytes that are no letters are dropped");
            Assert.That(iter.Count(c => c.Args.Contains("--maxiterate") && c.Rc == 0), Is.GreaterThanOrEqualTo(150), "dvtditr with the extended alphabet");
            var matrix = Cases.Value.Where(c => c.Group == "text-matrix").ToList();
            Assert.That(matrix.Count(c => c.Rc == 0), Is.GreaterThanOrEqualTo(200), "the entries that the file sets");
            Assert.That(Cases.Value.Count(c => c.Err.Contains("Score(", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(8), "the entries that the file sets");
            Assert.That(matrix.Count(c => c.Rc != 0), Is.GreaterThanOrEqualTo(3), "a matrix with a negative average of the diagonal");
            var err = Cases.Value.Where(c => c.Group == "text-err").ToList();
            Assert.That(err.Count(c => c.Err.Contains("Characters '= < >' can be used only in the title lines", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(15), "'=' '<' '>' in a sequence");
            Assert.That(err.Count(c => c.Err.Contains("The combination of --text and ( --anysymbol or --preservecase ) is impossible.", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(15), "--text with --anysymbol");
            Assert.That(err.Count(c => c.Err.Contains("cannot be used (", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(30), "a letter that cannot be a letter");
            Assert.That(err.Count(c => c.Err.Contains("Format error in this line?", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(30), "a malformed matrix line");
            Assert.That(err.Count(c => c.Err.Contains("Cannot open missing.txt", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(10), "no matrix file");
            Assert.That(err.Count(c => c.Err.Contains("Unrealistic scoring matrix", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(5), "a negative average of the diagonal");
        });
    }

    [Test]
    public void Text_AlignsEveryCharacterAsALetter()
    {
        // installed mafft 7.505-1: mafft --quiet --text in.fa   (case counts, digits and punctuation are letters, blanks are removed)
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">a\nhello world\n>b\nhelo wrld\n>c\nHallo world\n",
        };
        var r = MafftProgram.Run(["--quiet", "--text", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.Zero);
        Assert.That(r.Stdout, Is.EqualTo(">a\nhelloworld\n>b\n-helow-rld\n>c\nHalloworld\n"));
        Assert.That(r.Stderr, Is.Empty);
        var clustal = MafftProgram.Run(["--quiet", "--text", "--clustalout", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(clustal.Stdout, Is.EqualTo("CLUSTAL format alignment by MAFFT NW-NS-2 (v7.505)\n\n\na               helloworld\nb               -helow-rld\nc               Halloworld\n                   *** ***\n"), "only '*' marks the columns of identical characters (the case counts)");
    }

    [Test]
    public void TextMatrix_ScoresTheFileSets()
    {
        // installed mafft 7.505-1: the scores of the file replace the entries (and the matrix is rescaled to an average diagonal of 600): the stderr of the program names every entry
        var files = new Dictionary<string, string>
        {
            [SequenceAligner_MafftProgram_Tests.InputName] = ">a\nabcab\n>b\nabab\n>c\nbcab\n",
            ["m.txt"] = "0x61 0x61 5 # a\n0x62 0x63 2.5\n",
        };
        var r = MafftProgram.Run(["--textmatrix", "m.txt", SequenceAligner_MafftProgram_Tests.InputName], files);
        Assert.That(r.ExitCode, Is.Zero);
        Assert.That(r.Stderr, Does.Contain("Score(a=0x61,a=0x61)=5.000000\nScore(b=0x62,c=0x63)=2.500000\n"));
    }

    [Test]
    public void OneIterationWithText_IsRefused()
    {
        // installed mafft 7.505-1: 24 of 24 inputs make dooneiteration crash (valgrind: "Invalid write of size 8", "stack overflow in your program's main thread"), the message of the script names a process id
        var files = new Dictionary<string, string> { [SequenceAligner_MafftProgram_Tests.InputName] = ">a\nabcab\n>b\nabab\n>c\nbcab\n" };
        Assert.Throws<NotSupportedException>(() => MafftProgram.Run(["--text", "--oneiteration", SequenceAligner_MafftProgram_Tests.InputName], files));
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
