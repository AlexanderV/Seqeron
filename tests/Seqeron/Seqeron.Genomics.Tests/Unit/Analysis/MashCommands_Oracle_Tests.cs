using System.IO.Compression;
using System.Text.Json;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// A5-10 (B09 F174): the text of <c>mash bounds</c>, <c>mash dist</c> (<c>-t -v -d -i -n -k -s -S</c>) and <c>mash triangle</c> (<c>-E -v -d -i</c>) from <see cref="MashCommands"/> against the installed Mash 2.3, byte for byte,
/// on 220 random cases per command (TestData/Alignment/mash_commands_oracle.jsonl.gz, scripts/alignment/generate_mash_commands_oracle.py).
/// </summary>
[TestFixture]
public class MashCommands_Oracle_Tests
{
    private static readonly Lazy<List<JsonElement>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.mash_commands_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var rows = new List<JsonElement>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
            rows.Add(JsonDocument.Parse(line).RootElement.Clone());
        return rows;
    });

    private static MashInput[] Files(JsonElement r, string property) => [.. r.GetProperty(property).EnumerateArray().Select(f => new MashInput(f.GetProperty("name").GetString()!,
        [.. f.GetProperty("records").EnumerateArray().Select(x => new MashRecord(x.GetProperty("id").GetString()!, x.GetProperty("comment").GetString()!, x.GetProperty("seq").GetString()!))]))];

    private static double Fraction(JsonElement e) => e.ValueKind == JsonValueKind.Null ? 1.0 : e.GetDouble();

    private static string Produce(JsonElement r)
    {
        string cmd = r.GetProperty("cmd").GetString()!;
        if (cmd == "bounds")
            return MashCommands.Bounds(r.GetProperty("k").GetInt32(), r.GetProperty("p").GetDouble());
        int k = r.GetProperty("k").GetInt32(), s = r.GetProperty("s").GetInt32();
        uint seed = r.GetProperty("seed").GetUInt32();
        bool canonical = !r.GetProperty("noncanonical").GetBoolean(), individual = r.GetProperty("individual").GetBoolean();
        double v = Fraction(r.GetProperty("v")), d = Fraction(r.GetProperty("d"));
        return cmd == "dist"
            ? MashCommands.Dist(Files(r, "refs"), Files(r, "queries"), k, s, seed, canonical, individual, r.GetProperty("table").GetBoolean(), v, d)
            : MashCommands.Triangle(Files(r, "refs"), k, s, seed, canonical, individual, r.GetProperty("edge").GetBoolean(), v, d);
    }

    private static void Check(string command)
    {
        var rows = Oracle.Value.Where(r => r.GetProperty("cmd").GetString() == command).ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(200), command);
        var failures = new List<string>();
        int filtered = 0;
        foreach (var r in rows)
        {
            string expected = r.GetProperty("expected").GetString()!;
            string got = Produce(r);
            if (got != expected)
                failures.Add($"{command} {r.GetRawText()[..Math.Min(120, r.GetRawText().Length)]}\n   expected {expected[..Math.Min(200, expected.Length)]}\n   actual   {got[..Math.Min(200, got.Length)]}");
            if (command != "bounds" && (r.GetProperty("v").ValueKind != JsonValueKind.Null || r.GetProperty("d").ValueKind != JsonValueKind.Null))
                filtered++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(2)));
        if (command != "bounds")
            Assert.That(filtered, Is.GreaterThan(rows.Count / 10), "rows with a -v / -d filter");
    }

    [Test]
    public void FilesWithoutARecordOfKSymbols_AreRefused_AsTheProgramRefusesThem_OnAllErrorRows()
    {
        // A6-6: `mash dist` / `mash triangle` exit with status 1 and "WARNING: All fasta records in input files were shorter than the k-mer size (k)." (or "ERROR: Did not find fasta records in "input files"." for a file without records);
        // with -i the short record is simply not sketched, but a dist -i with no reference sequence at all crashes the program (SIGFPE, kept by the generator as a negative status)
        var rows = Oracle.Value.Where(r => r.GetProperty("cmd").GetString() == "error").ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(80));
        var statuses = rows.Select(r => r.GetProperty("status").GetInt32()).Distinct().ToList();
        Assert.That(statuses, Does.Contain(1).And.Contain(-8), "refusals with status 1 and the crash of dist -i");
        int warnings = 0, errors = 0;
        foreach (var r in rows)
        {
            string sub = r.GetProperty("sub").GetString()!;
            int k = r.GetProperty("k").GetInt32(), s = r.GetProperty("s").GetInt32();
            uint seed = r.GetProperty("seed").GetUInt32();
            bool canonical = !r.GetProperty("noncanonical").GetBoolean(), individual = r.GetProperty("individual").GetBoolean();
            var ex = Assert.Throws<ArgumentException>(() => _ = sub == "dist"
                ? MashCommands.Dist(Files(r, "refs"), Files(r, "queries"), k, s, seed, canonical, individual, r.GetProperty("table").GetBoolean())
                : MashCommands.Triangle(Files(r, "refs"), k, s, seed, canonical, individual, r.GetProperty("edge").GetBoolean()), r.GetRawText()[..Math.Min(160, r.GetRawText().Length)]);
            if (r.GetProperty("status").GetInt32() == 1)
            {
                Assert.That(ex!.Message, Does.StartWith(r.GetProperty("message").GetString()!), r.GetRawText()[..Math.Min(160, r.GetRawText().Length)]);
                if (ex.Message.StartsWith("WARNING", StringComparison.Ordinal))
                    warnings++;
                else
                    errors++;
            }
        }

        Assert.That(warnings, Is.GreaterThan(30));
        Assert.That(errors, Is.GreaterThan(10));
    }

    [Test] public void Bounds_EqualsMash_OnAllRows() => Check("bounds");
    [Test] public void Dist_TableAndFilters_EqualsMash_OnAllRows() => Check("dist");
    [Test] public void Triangle_MatrixAndEdgeList_EqualsMash_OnAllRows() => Check("triangle");

    [Test]
    public void Bounds_DefaultParameters_FirstRowsOfTheProgram()
    {
        // mash bounds -k 21 -p 0.99 (2.3)
        string text = MashCommands.Bounds();
        Assert.That(text, Does.StartWith("\nParameters (run with -h for details):\n   k:   21\n   p:   0.99\n\n\tMash distance\nSketch\t0.05\t0.1\t0.15\t0.2\t0.25\t0.3\t0.35\t0.4\n100\t0.0270708\t0.0867606\tinf\tinf"));
        Assert.That(text, Does.Contain("\n1000\t0.00675661\t0.0157591\t0.0322692\t0.0630219\tinf"));
        Assert.Throws<ArgumentOutOfRangeException>(() => MashCommands.Bounds(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => MashCommands.Bounds(21, 1.0));
    }

    [Test]
    public void Dist_ShortSequenceInIndividualMode_IsNotSketched()
    {
        var files = new[] { new MashInput("a.fa", [new MashRecord("long", "", new string('A', 40) + "CGTACGTAGCTAGCTAGCATCGATCGAT"), new MashRecord("short", "", "ACG")]) };
        string text = MashCommands.Dist(files, files, k: 9, sketchSize: 20, individual: true, table: true);
        Assert.That(text, Does.StartWith("#query\tlong\n"));
        Assert.That(text, Does.Not.Contain("short"));
    }

    private static MashInput Good(string name) => new(name, [new MashRecord("g", "", "ACGTTGCAAGCTTGACCGTAGCTAGGCTAACGTAGGCATCGATCAGCTAGCATCGACTAGC")]);

    [Test]
    public void Dist_And_Triangle_RefuseNullInputs_WithTheParameterName()
    {
        // A6-4: a null element of the lists, a MashInput without Records / Name, a default MashRecord (null Sequence) - typed exceptions that name the parameter, never a NullReferenceException
        var ok = new[] { Good("a.fa") };
        Assert.That(Assert.Throws<ArgumentNullException>(() => MashCommands.Dist([null!], ok, 11, 20))!.ParamName, Is.EqualTo("references"));
        Assert.That(Assert.Throws<ArgumentNullException>(() => MashCommands.Dist(ok, [Good("q"), null!], 11, 20))!.ParamName, Is.EqualTo("queries"));
        Assert.That(Assert.Throws<ArgumentNullException>(() => MashCommands.Triangle([Good("a"), null!], 11, 20))!.ParamName, Is.EqualTo("inputs"));
        Assert.That(Assert.Throws<ArgumentException>(() => MashCommands.Triangle([Good("a"), new MashInput("b", null!)], 11, 20))!.Message, Does.StartWith("Input 1 ('b') has no Records list."));
        Assert.That(Assert.Throws<ArgumentException>(() => MashCommands.Triangle([Good("a"), new MashInput(null!, [])], 11, 20))!.ParamName, Is.EqualTo("inputs"));
        var ex = Assert.Throws<ArgumentException>(() => MashCommands.Dist(ok, [new MashInput("q", [default])], 11, 20))!;
        Assert.That(ex.ParamName, Is.EqualTo("queries"));
        Assert.That(ex.Message, Does.StartWith("Record 0 of input 0 ('q') has no Sequence."));
        Assert.That(Assert.Throws<ArgumentException>(() => MashCommands.Dist(ok, [new MashInput("q", [new MashRecord(null!, "", "ACGTACGTACGTACGT")])], 11, 20, individual: true))!.Message, Does.Contain("has no Id"));
        Assert.Throws<ArgumentNullException>(() => MashCommands.Dist(null!, ok));
        Assert.Throws<ArgumentNullException>(() => MashCommands.Triangle(null!));
    }

    [Test]
    public void Dist_And_Triangle_RefuseOutOfRangeSketchParameters_BeforeReadingTheInputs()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => MashCommands.Dist([], [], 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => MashCommands.Dist([], [], 33));
        Assert.Throws<ArgumentOutOfRangeException>(() => MashCommands.Triangle([], 21, 0));
    }

    [Test]
    public void Dist_IndividualWithoutAnyReferenceSequence_And_TriangleMatrixWithoutSequences_AreRefused()
    {
        // mash dist -i short.fa g.fa -> SIGFPE; mash triangle -i short.fa short.fa -> std::out_of_range abort; the edge list of the latter is empty (status 0)
        var shortOnly = new MashInput("s.fa", [new MashRecord("s", "", "ACG")]);
        Assert.Throws<ArgumentException>(() => MashCommands.Dist([shortOnly], [Good("q")], 11, 20, individual: true));
        Assert.Throws<ArgumentException>(() => MashCommands.Triangle([shortOnly, shortOnly], 11, 20, individual: true));
        Assert.That(MashCommands.Triangle([shortOnly, shortOnly], 11, 20, individual: true, edgeList: true), Is.Empty);
        Assert.That(MashCommands.Dist([Good("r")], [shortOnly], 11, 20, individual: true), Is.Empty);                    // status 0, empty output
    }

    [Test]
    public void Dist_And_Triangle_ObserveTheCancellationToken_PerInputAndPerPair()
    {
        // cancelled while the third input's records are being enumerated (the token is polled at the start of the next input): nothing is returned
        using var cts = new CancellationTokenSource();
        var inputs = new List<MashInput> { Good("a"), Good("b"), new("c", new CancellingRecords(cts)), Good("d") };
        Assert.Throws<OperationCanceledException>(() => MashCommands.Triangle(inputs, 11, 20, cancellationToken: cts.Token));

        using var pre = new CancellationTokenSource();
        pre.Cancel();
        Assert.Throws<OperationCanceledException>(() => MashCommands.Dist([Good("r")], [Good("q")], 11, 20, cancellationToken: pre.Token));
        Assert.Throws<OperationCanceledException>(() => MashCommands.Triangle([Good("r"), Good("q")], 11, 20, cancellationToken: pre.Token));
    }

    private sealed class CancellingRecords(CancellationTokenSource cts) : IReadOnlyList<MashRecord>
    {
        public int Count => 1;
        public MashRecord this[int index] { get { cts.Cancel(); return new MashRecord("c", "", "ACGTTGCAAGCTTGACCGTAGCTAGGCTAACGTAGG"); } }
        public IEnumerator<MashRecord> GetEnumerator() { yield return this[0]; }
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() => GetEnumerator();
    }
}
