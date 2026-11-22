using System.Globalization;
using System.IO.Compression;
using System.Security.Cryptography;
using System.Text;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19b, B09 F283): <see cref="SequenceAligner.MultipleAlignMafft"/> / <see cref="SequenceAligner.RunMafft"/> for three and more sequences (the progressive strategies FFT-NS-1 / -2, NW-NS-1 / -2
/// of MAFFT 7.505: 6-mer distance, UPGMA guide tree, <c>treebase</c> with <c>A__align</c> and <c>Falign</c>) against the installed MAFFT 7.505 (Debian 7.505-1; scripts/alignment/generate_mafft_progressive_oracle.py,
/// mafft_progressive_oracle.tsv.gz). Every committed run must give the program's rows and its standard error (SHA-256 of the whole text, the full text for every 25th run).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftProgressive_Tests
{
    internal sealed record Row(string Id, string Strategy, string Options, string[] Seqs, int Rc, string[] Rows, string Sha, string Stderr);

    internal static readonly Lazy<List<Row>> Oracle = new(() => LoadOracle("mafft_progressive_oracle.tsv.gz"));

    internal static List<Row> LoadOracle(string file)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment." + file;
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException(name);
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], f[2], f[3].Split('|'), int.Parse(f[4], CultureInfo.InvariantCulture), f[5].Length == 0 ? [] : f[5].Split('|'), f[6], Unescape(f[7])));
        }
        return list;
    }

    private static string Unescape(string s)
    {
        var sb = new StringBuilder(s.Length);
        int i = 0;
        while (i < s.Length)
        {
            if (s[i] == '\\' && i + 1 < s.Length)
            {
                sb.Append(s[i + 1] switch { 'n' => '\n', 'r' => '\r', 'b' => '\b', 't' => '\t', _ => s[i + 1] });
                i += 2;
            }
            else
            {
                sb.Append(s[i]);
                i++;
            }
        }
        return sb.ToString();
    }

    internal static MafftOptions ParseOptions(string options)
    {
        var o = new MafftOptions { AnySymbol = false };
        if (options.Length == 0)
            return o;
        var t = options.Split(' ');
        int i = 0;
        string Next() => t[i++];
        double D() => double.Parse(Next(), CultureInfo.InvariantCulture);
        int I() => int.Parse(Next(), CultureInfo.InvariantCulture);
        while (i < t.Length)
        {
            switch (Next())
            {
                case "--quiet": o = o with { Quiet = true }; break;
                case "--anysymbol": o = o with { AnySymbol = true }; break;
                case "--op": o = o with { GapOpen = D() }; break;
                case "--ep": o = o with { Offset = D() }; break;
                case "--exp": o = o with { GapExtension = D() }; break;
                case "--opdist": o = o with { GapOpenDistance = D() }; break;
                case "--bl": o = o with { Blosum = I() }; break;
                case "--jtt": case "--kimura": o = o with { ScoringModel = MafftScoringModel.Jtt, Pam = I() }; break;
                case "--tm": o = o with { ScoringModel = MafftScoringModel.Tm, Pam = I() }; break;
                case "--kappa": o = o with { Kappa = I() }; break;
                case "--fmodel": o = o with { FrequencyModel = true }; break;
                case "--nwildcard": o = o with { WildcardN = true }; break;
                case "--nuc": o = o with { SequenceType = MafftSequenceType.Nucleotide }; break;
                case "--amino": o = o with { SequenceType = MafftSequenceType.Amino }; break;
                case "--retree": o = o with { Retree = I() }; break;
                case "--nofft": o = o with { Fft = false }; break;
                case "--fft": o = o with { Fft = true, ForceFft = true }; break;
                case "--6merpair": o = o with { TupleSize = 6 }; break;
                case "--10merpair": o = o with { TupleSize = 10 }; break;
                case "--leavegappyregion": o = o with { LegacyGapPenalty = true }; break;
                case "--termgappenalty": o = o with { TerminalGapPenalty = true }; break;
                case "--unweight": o = o with { Unweight = true }; break;
                case "--averagelinkage": o = o with { TreeMixing = 1.0 }; break;
                case "--minimumlinkage": o = o with { TreeMixing = 0.0 }; break;
                case "--mixedlinkage": o = o with { TreeMixing = D() }; break;
                case "--reorder": o = o with { Reorder = true }; break;
                case "--thread": o = o with { Threads = I() }; break;
                case "--threadit": o = o with { ThreadsIterate = I() }; break;
                case "--maxiterate": o = o with { MaxIterate = I() }; break;
                case "--globalpair": o = o with { PairLibrary = MafftPairLibrary.Global }; break;
                case "--localpair": o = o with { PairLibrary = MafftPairLibrary.Local }; break;
                case "--genafpair": o = o with { PairLibrary = MafftPairLibrary.Genaf }; break;
                case "--weighti": o = o with { WeightI = D() }; break;
                case "--gop": o = o with { GlobalPairOpen = D() }; break;
                case "--gep": o = o with { GlobalPairOffset = D() }; break;
                case "--gexp": o = o with { GlobalPairExtension = D() }; break;
                case "--lop": o = o with { LocalPairOpen = D() }; break;
                case "--lep": o = o with { LocalPairOffset = D() }; break;
                case "--lexp": o = o with { LocalPairExtension = D() }; break;
                case "--LOP": o = o with { GenafOpen = D() }; break;
                case "--LEXP": o = o with { GenafExtension = D() }; break;
                default: throw new InvalidOperationException(t[i - 1]);
            }
        }
        return o;
    }

    public static IEnumerable<TestCaseData> Strategies() => Oracle.Value.Select(r => r.Strategy).Distinct().Select(s => new TestCaseData(s).SetArgDisplayNames(s));

    [Test]
    public void Oracle_HasEveryStrategy()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(5000));
            foreach (var g in rows.GroupBy(r => r.Strategy))
                Assert.That(g.Count(), Is.GreaterThanOrEqualTo(200), g.Key);
            Assert.That(rows.Select(r => r.Strategy).Distinct().Count(), Is.GreaterThanOrEqualTo(25));
            Assert.That(rows.Count(r => r.Rc != 0), Is.GreaterThanOrEqualTo(10), "the plain command stops on letters outside the alphabet");
            Assert.That(rows.Count(r => r.Stderr.Length > 0), Is.GreaterThanOrEqualTo(150), "full texts of the standard error");
            Assert.That(rows.Max(r => r.Seqs.Length), Is.GreaterThanOrEqualTo(40));
        });
    }

    [TestCaseSource(nameof(Strategies))]
    public void RowsAndStandardError_EqualMafft(string strategy) => CompareStrategy(Oracle.Value, strategy);

    internal static void CompareStrategy(List<Row> oracle, string strategy)
    {
        var failures = new List<string>();
        int n = 0;
        foreach (var r in oracle.Where(x => x.Strategy == strategy))
        {
            n++;
            var o = ParseOptions(r.Options);
            if (r.Rc != 0)
            {
                if (!Throws<ArgumentException>(() => SequenceAligner.RunMafft(r.Seqs, o)))
                    failures.Add(r.Id + ": the program refuses (rc " + r.Rc + "), the port does not");
                continue;
            }
            MafftRunResult res;
            try
            {
                res = SequenceAligner.RunMafft(r.Seqs, o);
            }
            catch (Exception e)
            {
                failures.Add(r.Id + ": " + e.GetType().Name + " " + e.Message);
                continue;
            }
            if (!res.Alignment.AlignedSequences.SequenceEqual(r.Rows))
            {
                failures.Add(r.Id + " (" + r.Options + "): rows differ");
                continue;
            }
            if (r.Sha != "-")
            {
                string sha = Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(res.StandardError))).ToLowerInvariant();
                if (sha != r.Sha)
                    failures.Add(r.Id + " (" + r.Options + "): standard error differs");
                else if (r.Stderr.Length > 0 && r.Stderr != res.StandardError)
                    failures.Add(r.Id + ": standard error text differs although the hash is equal");
            }
        }
        Assert.That(n, Is.GreaterThan(0));
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    internal static bool Throws<T>(Action a) where T : Exception
    {
        try
        {
            a();
            return false;
        }
        catch (T)
        {
            return true;
        }
    }

    // Hand run of the installed mafft 7.505 (mafft --quiet in.fa) on three short proteins and three short DNA sequences.
    [Test]
    public void HandRun_ThreeProteins_FftNs2()
    {
        var rows = SequenceAligner.MultipleAlignMafft(["MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ"], new MafftOptions { AnySymbol = false }).AlignedSequences;
        Assert.That(rows.Select(r => r.Replace("-", "")), Is.EqualTo(new[] { "MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ" }));
        Assert.That(rows.Select(r => r.Length).Distinct().Count(), Is.EqualTo(1));
    }

    [Test]
    public void Threads_DoNotChangeTheAlignment()
    {
        // the program is deterministic with several threads (verified on 220 committed sets: strategy 'thread'); the port is a single-thread run
        foreach (var r in Oracle.Value.Where(x => x.Strategy == "thread" && x.Rc == 0).Take(40))
        {
            var o = ParseOptions(r.Options);
            Assert.That(SequenceAligner.MultipleAlignMafft(r.Seqs, o).AlignedSequences, Is.EqualTo(SequenceAligner.MultipleAlignMafft(r.Seqs, o with { Threads = 0 }).AlignedSequences), r.Id);
        }
    }

    [Test]
    public void Reorder_PutsTheSequencesInTheOrderOfTheGuideTree()
    {
        foreach (var r in Oracle.Value.Where(x => x.Strategy == "reorder" && x.Rc == 0).Take(20))
        {
            var plain = SequenceAligner.MultipleAlignMafft(r.Seqs, ParseOptions("")).AlignedSequences;
            var reordered = SequenceAligner.MultipleAlignMafft(r.Seqs, ParseOptions("--reorder")).AlignedSequences;
            Assert.That(reordered.OrderBy(x => x, StringComparer.Ordinal), Is.EqualTo(plain.OrderBy(x => x, StringComparer.Ordinal)), "the same rows in another order: " + r.Id);
        }
    }

    [Test]
    public void Contracts()
    {
        var three = new[] { "ACGTACGTAGCT", "ACGTACGAGCT", "ACGTCGTAGCT" };
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.MultipleAlignMafft(three).AlignedSequences, Has.Length.EqualTo(3));
            Assert.Throws<NotSupportedException>(() => SequenceAligner.MultipleAlignMafft(three, new MafftOptions { Retree = 0 }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignMafft(three, new MafftOptions { TupleSize = 10, SequenceType = MafftSequenceType.Amino }), "tuplesize must be 6 for aa sequence");
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignMafft([.. three, "ACG-"], new MafftOptions { Blosum = 63 }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignMafft([.. three, "---"]), "no residues");
            Assert.Throws<OperationCanceledException>(() => SequenceAligner.MultipleAlignMafft(three, new MafftOptions { CancellationToken = new CancellationToken(true) }));
            Assert.That(SequenceAligner.RunMafft(three, new MafftOptions { Quiet = true }).StandardError, Is.Empty);
            Assert.That(SequenceAligner.RunMafft(three).StandardError, Does.Contain("Strategy:\n FFT-NS-2 (Fast but rough)\n Progressive method (guide trees were built 2 times.)\n"));
        });
    }
}
