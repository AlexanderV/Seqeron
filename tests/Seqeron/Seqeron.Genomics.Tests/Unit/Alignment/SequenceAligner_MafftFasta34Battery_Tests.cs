using System.Globalization;
using System.IO.Compression;
using System.Text.RegularExpressions;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19c, B09 F323): the port of FASTA 3.4 that <c>splittbfast -S</c> uses (<c>SequenceAligner.MafftFasta34</c>) against the RAW output of the real <c>fasta34 -z3 -m10 -p|-n -Q -H -b n -E n -d0 query library 2|6</c>
/// (FASTA 34.26.5, sha256 7491909b..., the user's Docker run <c>handoff/A2-19c-fasta34-battery</c>, scripts/alignment/convert_fasta34_battery_returned.py, fasta34_battery_oracle.tsv.gz): the 20 probes of <c>getfastascoremtx</c> and 320 random calls as
/// <c>splittbfast.c</c> makes them (protein and DNA, 1 - 60 library sequences of 6 - 400 letters, unknown letters, low complexity, duplicates, reverse-complement frames of DNA). Every header value (query length, library residues / sequences, <c>join</c>,
/// <c>opt</c> cut-off) and every printed line of the list of the best scores (library sequence, length, frame, <c>opt</c>, bit score <c>%.1f</c>, E-value <c>%7.2g</c>, the ORDER of the lines, the E-value cut and the "No library sequences" case) must be equal.
/// </summary>
[TestFixture]
[Parallelizable(ParallelScope.Children)]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MafftFasta34Battery_Tests
{
    internal sealed record RawHit(int Sequence, int Length, string Frame, int Opt, string Bits, string E);

    internal sealed record Call(string Id, bool Dna, string Query, string[] Library, int Rc, string Raw)
    {
        public int? QueryLength => Regex.Match(Raw, @"^\s+1>>>\S+\s+(\d+) (?:aa|nt) - (\d+) (?:aa|nt)", RegexOptions.Multiline) is { Success: true } m ? int.Parse(m.Groups[2].Value, CultureInfo.InvariantCulture) : null;

        public (long Residues, int Sequences)? Library_ => Regex.Match(Raw, @"^\s+(\d+) residues in\s+(\d+) sequences", RegexOptions.Multiline) is { Success: true } m
            ? (long.Parse(m.Groups[1].Value, CultureInfo.InvariantCulture), int.Parse(m.Groups[2].Value, CultureInfo.InvariantCulture)) : null;

        public (int Join, int OptCut)? Cuts => Regex.Match(Raw, @"^ join: (\d+), opt: (\d+),", RegexOptions.Multiline) is { Success: true } m
            ? (int.Parse(m.Groups[1].Value, CultureInfo.InvariantCulture), int.Parse(m.Groups[2].Value, CultureInfo.InvariantCulture)) : null;

        public bool NoSignificant => Raw.Contains("!! No library sequences with E() <", StringComparison.Ordinal);

        public List<RawHit> Hits
        {
            get
            {
                var list = new List<RawHit>();
                var lines = Raw.Split('\n');
                int start = Array.FindIndex(lines, l => l.StartsWith("The best scores are:", StringComparison.Ordinal));
                if (start < 0)
                    return list;
                for (int i = start + 1; i < lines.Length && lines[i].Length > 0; i++)
                {
                    var m = Regex.Match(lines[i], @"^\+===========\+(\d+)\s+\(\s*(\d+)\)(?: \[(\w)\])?\s+(\d+) (\S+)\s+(\S+)$");
                    Assert.That(m.Success, Is.True, $"{Id}: unreadable line '{lines[i]}'");
                    list.Add(new RawHit(int.Parse(m.Groups[1].Value, CultureInfo.InvariantCulture), int.Parse(m.Groups[2].Value, CultureInfo.InvariantCulture), m.Groups[3].Value,
                        int.Parse(m.Groups[4].Value, CultureInfo.InvariantCulture), m.Groups[5].Value, m.Groups[6].Value));
                }

                return list;
            }
        }
    }

    internal static readonly Lazy<List<Call>> Oracle = new(() =>
    {
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream("Seqeron.Genomics.Tests.TestData.Alignment.fasta34_battery_oracle.tsv.gz") ?? throw new InvalidOperationException("oracle");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var list = new List<Call>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Call(f[0], f[1] == "n", f[2], f[3].Split('|'), int.Parse(f[4], CultureInfo.InvariantCulture), Unescape(f[5])));
        }

        return list;
    });

    private static string Unescape(string s) => Regex.Replace(s, @"\\(.)", m => m.Groups[1].Value switch { "n" => "\n", "t" => "\t", "r" => "\r", var other => other }, RegexOptions.Singleline);

    /// <summary>C <c>printf("%.2g")</c> of a non-negative double (the E-values of the list are printed with <c>%7.2g</c>).</summary>
    internal static string FormatG2(double v)
    {
        if (v == 0.0)
            return "0";
        string e = v.ToString("E1", CultureInfo.InvariantCulture);        // d.dE+ddd, correctly rounded
        int x = int.Parse(e[(e.IndexOf('E') + 1)..], CultureInfo.InvariantCulture);
        if (x < -4 || x >= 2)
        {
            string mant = e[..e.IndexOf('E')];
            if (mant.Contains('.'))
                mant = mant.TrimEnd('0').TrimEnd('.');
            return mant + "e" + (x < 0 ? "-" : "+") + Math.Abs(x).ToString("D2", CultureInfo.InvariantCulture);
        }

        string fixedText = v.ToString("F" + (1 - x).ToString(CultureInfo.InvariantCulture), CultureInfo.InvariantCulture);
        if (fixedText.Contains('.'))
            fixedText = fixedText.TrimEnd('0').TrimEnd('.');
        return fixedText;
    }

    [Test]
    public void Oracle_Shape()
    {
        var calls = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(calls, Has.Count.EqualTo(340));
            Assert.That(calls.Count(c => c.Id.StartsWith("probe-", StringComparison.Ordinal)), Is.EqualTo(20));
            Assert.That(calls.All(c => c.Rc == 0));
            Assert.That(calls.Count(c => c.Dna), Is.EqualTo(160));
            Assert.That(calls.Count(c => c.NoSignificant), Is.GreaterThan(0), "calls with no significant library sequence");
            Assert.That(calls.Sum(c => c.Hits.Count), Is.GreaterThan(1000));
            Assert.That(calls.Count(c => c.Hits.Any(h => h.Frame == "r")), Is.GreaterThan(20), "DNA calls with a reverse-complement line");
            Assert.That(FormatG2(0.034), Is.EqualTo("0.034"));
            Assert.That(FormatG2(0.41), Is.EqualTo("0.41"));
            Assert.That(FormatG2(2.7), Is.EqualTo("2.7"));
            Assert.That(FormatG2(1.2e-5), Is.EqualTo("1.2e-05"));
            Assert.That(FormatG2(30.0), Is.EqualTo("30"));
            Assert.That(FormatG2(123.0), Is.EqualTo("1.2e+02"));
        });
    }

    [Test]
    public void EveryRawOutput_EqualsThePort_HeaderListAndOrder()
    {
        var failures = new List<string>();
        int hits = 0, none = 0;
        foreach (var c in Oracle.Value)
        {
            var port = SequenceAligner.MafftFasta34.SearchListing(c.Query, c.Library, c.Dna);
            var raw = c.Hits;
            if (c.QueryLength is { } n0 && n0 != port.QueryLength)
                failures.Add($"{c.Id}: query length {port.QueryLength}, fasta34 {n0}");
            if (c.Library_ is { } lib && (lib.Residues != port.LibraryResidues || lib.Sequences != port.LibrarySequences))
                failures.Add($"{c.Id}: library {port.LibraryResidues} residues in {port.LibrarySequences} sequences, fasta34 {lib.Residues} in {lib.Sequences}");
            if (c.Cuts is { } cuts && (cuts.Join != port.Join || cuts.OptCut != port.OptCut))
                failures.Add($"{c.Id}: join/opt {port.Join}/{port.OptCut}, fasta34 {cuts.Join}/{cuts.OptCut}");
            if (c.NoSignificant != port.NoSignificant)
                failures.Add($"{c.Id}: no-significant {port.NoSignificant}, fasta34 {c.NoSignificant}");
            if (c.NoSignificant)
                none++;
            if (raw.Count != port.Hits.Count)
            {
                failures.Add($"{c.Id}: {port.Hits.Count} lines, fasta34 {raw.Count}");
                continue;
            }

            for (int i = 0; i < raw.Count; i++)
            {
                hits++;
                var r = raw[i];
                var p = port.Hits[i];
                string frame = c.Dna ? (p.Frame == 0 ? "f" : "r") : "";
                string bits = p.Bits.ToString("F1", CultureInfo.InvariantCulture);
                string e = FormatG2(p.EValue);
                if (r.Sequence != p.Sequence || r.Length != p.Length || r.Frame != frame || r.Opt != p.Opt || r.Bits != bits || r.E != e)
                    failures.Add($"{c.Id} line {i}: port seq {p.Sequence} len {p.Length} [{frame}] opt {p.Opt} bits {bits} E {e}; fasta34 seq {r.Sequence} len {r.Length} [{r.Frame}] opt {r.Opt} bits {r.Bits} E {r.E}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} differences; first: {failures.FirstOrDefault()}");
        TestContext.Out.WriteLine($"fasta34 battery: {Oracle.Value.Count} calls, {hits} printed lines equal, {none} calls with no significant sequence");
    }

    [Test]
    public void ScoresOfTheListAreWhatSplittbfastReads()
    {
        // ReadFasta34m10_scoreonly: the opt of a library sequence = the first line of that sequence in the list (the forward frame for DNA), 0 for the sequences that are not on the list
        var failures = new List<string>();
        foreach (var c in Oracle.Value)
        {
            var expected = new double[c.Library.Length];
            var seen = new HashSet<int>();
            foreach (var h in c.Hits)
                if ((!c.Dna || h.Frame == "f") && seen.Add(h.Sequence))
                    expected[h.Sequence] += h.Opt;
            var got = SequenceAligner.MafftFasta34.Search(c.Query, c.Library, c.Dna);
            if (!got.SequenceEqual(expected))
                failures.Add($"{c.Id}: scores differ");
        }

        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    [Test]
    public void ProbesOfGetfastascoremtx_GiveTheDiagonalOfTheMatrix()
    {
        // getfastascoremtx: the score of six equal residues against the probe library, divided by six: the 20 diagonal values printed on standard error ("A: 5.000000", ...)
        var table = SequenceAligner.MafftFasta34.ScoreMatrix(false, _ => { });
        foreach (var c in Oracle.Value.Where(c => c.Id.StartsWith("probe-", StringComparison.Ordinal)))
        {
            char aa = c.Id[^1];
            var hit = c.Hits.Single();
            Assert.That(table[aa, aa] * 6, Is.EqualTo(hit.Opt), c.Id);
        }
    }
}
