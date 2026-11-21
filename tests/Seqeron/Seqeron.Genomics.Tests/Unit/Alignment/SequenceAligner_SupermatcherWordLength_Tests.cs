using System.Diagnostics;
using System.Globalization;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-LOCAL-001 (audit A6-5, B09 F180): <see cref="SequenceAligner.Supermatcher"/> with <c>-wordlen</c> around and above the targets against the installed supermatcher 6.6.0
/// (TestData/Alignment/supermatcher_wordlen_oracle.tsv, scripts/alignment/generate_supermatcher_wordlen_oracle.py): the program prints <c>Error: no kmers found</c> on stderr whenever its word table is empty (every word length above the longest
/// target, or targets that hold only N / X) and lists every pair in its error file; the port returns the same text without running the Rabin-Karp loops (<c>-wordlen 2147483647</c> took 30 s in the port, 9 s in the program, which then crashes).
/// </summary>
[TestFixture]
public class SequenceAligner_SupermatcherWordLength_Tests
{
    private sealed record Row(int Line, string Type, string FastaA, string FastaB, int WordLength, string StdErr, string ErrFile, string Alignments);

    private static string Unescape(string s)
    {
        var sb = new System.Text.StringBuilder();
        bool escaped = false;
        foreach (char c in s)
        {
            if (escaped)
            {
                sb.Append(c switch { 'n' => '\n', 't' => '\t', _ => c });
                escaped = false;
            }
            else if (c == '\\')
            {
                escaped = true;
            }
            else
            {
                sb.Append(c);
            }
        }

        return sb.ToString();
    }

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.supermatcher_wordlen_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var rows = new List<Row>();
        string? line;
        int n = 1;
        while ((line = reader.ReadLine()) is not null)
        {
            n++;
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            rows.Add(new Row(n, f[1], Unescape(f[2]), Unescape(f[3]), int.Parse(f[4], CultureInfo.InvariantCulture), Unescape(f[5]), Unescape(f[6]), f[7]));
        }

        return rows;
    });

    private static List<EmbossSequence> ParseFasta(string fasta)
    {
        var list = new List<EmbossSequence>();
        foreach (string block in fasta.Split('>', StringSplitOptions.RemoveEmptyEntries))
        {
            var lines = block.TrimEnd('\n').Split('\n');
            string h = lines[0];
            int sp = h.IndexOf(' ', StringComparison.Ordinal);
            list.Add(new EmbossSequence(sp < 0 ? h : h[..sp], string.Concat(lines.Skip(1)), sp < 0 ? "" : h[(sp + 1)..]));
        }

        return list;
    }

    [Test]
    public void Supermatcher_WordLengthAroundAndAboveTheTargets_EqualsTheProgram_StderrErrorFileAndAlignments()
    {
        var failures = new List<string>();
        int noKmers = 0, aligned = 0;
        foreach (var r in Oracle.Value)
        {
            var result = SequenceAligner.Supermatcher(ParseFasta(r.FastaA), ParseFasta(r.FastaB), wordLength: r.WordLength);
            string alignments = string.Join(';', result.Alignments.Select(a => $"{a.First.Name}|{a.Second.Name}|{a.Score.ToString("0.0", CultureInfo.InvariantCulture)}"));
            string expected = string.Join(';', r.Alignments.Split(';', StringSplitOptions.RemoveEmptyEntries).Select(x => x.Split('|')).Select(p => $"{p[0]}|{p[1]}|{double.Parse(p[2], CultureInfo.InvariantCulture).ToString("0.0", CultureInfo.InvariantCulture)}"));
            if (result.StandardError != r.StdErr || result.ErrorText != r.ErrFile || alignments != expected)
                failures.Add($"line {r.Line} -wordlen {r.WordLength}\n   stderr   expected [{r.StdErr.Replace("\n", "\\n")}] actual [{result.StandardError.Replace("\n", "\\n")}]\n   errfile  expected [{r.ErrFile.Replace("\n", "\\n")}] actual [{result.ErrorText.Replace("\n", "\\n")}]\n   aligned  expected {expected} actual {alignments}");
            if (r.StdErr.Length > 0)
                noKmers++;
            else
                aligned++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {Oracle.Value.Count}:\n" + string.Join("\n", failures.Take(3)));
        Assert.That(noKmers, Is.GreaterThan(80), "rows where the word table is empty");
        Assert.That(aligned, Is.GreaterThan(60), "rows with a word table");
        Assert.That(Oracle.Value.Count(r => r.WordLength >= 100_000), Is.GreaterThan(5));
    }

    [Test]
    public void Supermatcher_AbsurdWordLength_ReturnsTheProgramsOutcome_WithoutTheWordLoops()
    {
        // two 23-base sequences, -wordlen 2147483647: the program loops ~9 s, prints "Error: no kmers found" and crashes; the port spent 30 s in the Rabin-Karp loops (no poll) and printed only the error file text
        var q = new[] { new EmbossSequence("q", "ACGTACGTTAGCCGATTAGGCAT") };
        var t = new[] { new EmbossSequence("t", "TTGCATGCATTAGCGATCGATCG") };
        var sw = Stopwatch.StartNew();
        foreach (int w in new[] { 24, 100, 100_000_000, int.MaxValue })
        {
            var r = SequenceAligner.Supermatcher(q, t, wordLength: w);
            Assert.That(r.Alignments, Is.Empty, $"-wordlen {w}");
            Assert.That(r.StandardError, Is.EqualTo("Error: no kmers found\n"), $"-wordlen {w}");
            Assert.That(r.ErrorText, Is.EqualTo("No wordmatch start points for q vs t. No alignment\n"), $"-wordlen {w}");
        }

        // -wordlen 23 (= the target) is an ordinary run: no stderr line
        Assert.That(SequenceAligner.Supermatcher(q, t, wordLength: 23).StandardError, Is.Empty);
        // generous guard against the loops (the old code needed 30 s for int.MaxValue alone)
        Assert.That(sw.Elapsed, Is.LessThan(TimeSpan.FromSeconds(5)));
    }

    [Test]
    public void Supermatcher_TargetsOfNOnly_HaveNoWords_LikeTheProgram()
    {
        // supermatcher q.fa tN.fa -wordlen 3: "Error: no kmers found" and the pair listed in the error file (installed 6.6.0)
        var r = SequenceAligner.Supermatcher([new EmbossSequence("q", "ACGTACGTTAGCCGATTAGGCAT")], [new EmbossSequence("t", new string('N', 30))], wordLength: 3);
        Assert.That((r.StandardError, r.ErrorText, r.Alignments.Count), Is.EqualTo(("Error: no kmers found\n", "No wordmatch start points for q vs t. No alignment\n", 0)));
    }

    [Test]
    public void SeqMatchAll_MegaMerge_DotMatcher_AbsurdWordSizes_DoNotSpin()
    {
        // the same parameter of the sibling word programs (audit A6-5): a word size above every sequence finishes at once (no match / no plot / the program's refusal), never a loop over the word length
        const string a = "ACGTACGTTAGCCGATTAGGCATCGATCGG", b = "TTGCATGCATTAGCGATCGATCGGATTACAGG";
        var sw = Stopwatch.StartNew();
        Assert.That(GenomicAnalyzer.SeqMatchAll([a, b], int.MaxValue), Is.Empty);
        Assert.That(GenomicAnalyzer.SeqMatchAll([a, b], EmbossSequenceQualifiers.None, int.MaxValue), Is.Empty);
        Assert.That(DotMatcher.Compute(a, b, windowSize: int.MaxValue, threshold: 10), Is.Empty);
        Assert.That(Assert.Throws<ArgumentException>(() => GenomicAnalyzer.MegaMerge(a, b, wordSize: int.MaxValue))!.Message, Does.Contain("No match found"));       // the program's own refusal
        Assert.That(sw.Elapsed, Is.LessThan(TimeSpan.FromSeconds(5)));
    }
}
