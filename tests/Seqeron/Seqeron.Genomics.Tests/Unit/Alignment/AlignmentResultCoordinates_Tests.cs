using System.Text.Json.Nodes;
using System.Text.RegularExpressions;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 / ALIGN-SEMI-001 (audit R14-9, B09 F302): the bridge from <see cref="AlignmentResult"/> to the coordinate world
/// (<see cref="AlignmentResultCoordinates"/>: block coordinates, BED / PSL / chain / MAF / SAM / Exonerate text, CIGAR (extended), MD, cs). Oracle 1: 210 alignments Biopython 1.88 <c>PairwiseAligner</c> finds (global, affine, local, fitting,
/// overlap; <c>alignment.coordinates</c> and <c>alignment.format(...)</c> of every format, and of <c>alignment[:, first:last]</c> for the trimmed end gaps), described by their rows; oracle 2: the extended CIGAR of parasail 1.3.4
/// (<c>nw_trace</c>, <c>sg_trace</c>, <c>sg_qx_trace</c>, <c>sg_dx_trace</c>); oracle 3: <c>cs</c> / <c>MD</c> / <c>--eqx</c> CIGAR of minimap2 2.x (see <c>AlignmentResultCigarOracle</c>).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class AlignmentResultCoordinates_Tests
{
    public sealed record Case(string Id, string Mode, string TName, string QName, string T, string Q, string Row1, string Row2, int Start1, int Start2, double Score, long[][] Coordinates, JsonObject Formats, JsonObject? Trimmed);

    private static string Unesc(string s) => Regex.Replace(s, @"\\(.)", m => m.Groups[1].Value switch { "t" => "\t", "n" => "\n", "r" => "\r", var c => c });

    private static readonly Lazy<List<Case>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.alignment_result_bridge_oracle.tsv";
        using var stream = TestDataResource.Open(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Case>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            var coords = ((JsonArray)JsonNode.Parse(Unesc(f[11]))!).Select(r => ((JsonArray)r!).Select(v => v!.GetValue<long>()).ToArray()).ToArray();
            list.Add(new Case(f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], int.Parse(f[8]), int.Parse(f[9]), double.Parse(f[10], System.Globalization.CultureInfo.InvariantCulture), coords,
                (JsonObject)JsonNode.Parse(Unesc(f[12]))!, JsonNode.Parse(Unesc(f[13])) as JsonObject));
        }

        return list;
    });

    private static IEnumerable<TestCaseData> Cases() => Oracle.Value.Select(c => new TestCaseData(c).SetName($"Bridge_{c.Id}_{c.Mode}"));

    private static AlignmentResult ResultOf(Case c, bool swap = false)
    {
        string r1 = swap ? c.Row2 : c.Row1, r2 = swap ? c.Row1 : c.Row2;
        int s1 = swap ? c.Start2 : c.Start1, s2 = swap ? c.Start1 : c.Start2;
        var type = c.Mode.StartsWith("local", StringComparison.Ordinal) ? AlignmentType.Local : c.Mode.StartsWith("global", StringComparison.Ordinal) ? AlignmentType.Global : AlignmentType.SemiGlobal;
        int n1 = r1.Count(ch => ch != '-'), n2 = r2.Count(ch => ch != '-');
        return new AlignmentResult(r1, r2, (int)Math.Round(c.Score, MidpointRounding.AwayFromZero), type, s1, s2, s1 + n1 - 1, s2 + n2 - 1) { ExactScore = c.Score };
    }

    private static readonly string[] FormatNames = ["bed", "bed3", "bed6", "psl", "chain", "maf", "sam", "sam:md", "exonerate", "exonerate:cigar"];

    [TestCaseSource(nameof(Cases))]
    public void Coordinates_EqualBiopython(Case c)
    {
        var a = ResultOf(c).ToBioCoordAlignment(c.TName, c.T, c.QName, c.Q);
        Assert.That(a.Coordinates, Is.EqualTo(c.Coordinates));
        Assert.That(a.Sequences[0].Length, Is.EqualTo(c.T.Length));
        Assert.That(a.Sequences[1].Length, Is.EqualTo(c.Q.Length));
    }

    [TestCaseSource(nameof(Cases))]
    public void Formats_EqualBiopython(Case c)
    {
        var result = ResultOf(c);
        foreach (string f in FormatNames)
            Assert.That(result.Format(f, c.TName, c.T, c.QName, c.Q), Is.EqualTo(c.Formats[f]!.GetValue<string>()), $"{c.Id} {c.Mode} {f}");
    }

    [TestCaseSource(nameof(Cases))]
    public void Formats_FirstIsQuery_EqualBiopython(Case c)
    {
        // SemiGlobalAlign(query, reference): the same alignment with the rows and the names swapped and the roles declared
        var result = ResultOf(c, swap: true);
        var options = new AlignmentBridgeOptions { FirstIsQuery = true };
        foreach (string f in FormatNames)
            Assert.That(result.Format(f, c.QName, c.Q, c.TName, c.T, options), Is.EqualTo(c.Formats[f]!.GetValue<string>()), $"{c.Id} {c.Mode} {f}");
    }

    [TestCaseSource(nameof(Cases))]
    public void Formats_TrimEndGaps_EqualBiopythonColumnSlice(Case c)
    {
        var result = ResultOf(c);
        var options = new AlignmentBridgeOptions { TrimEndGaps = true };
        foreach (string f in FormatNames)
            Assert.That(result.Format(f, c.TName, c.T, c.QName, c.Q, options), Is.EqualTo(c.Trimmed![f]!.GetValue<string>()), $"{c.Id} {c.Mode} {f}");
    }

    [TestCaseSource(nameof(Cases))]
    public void CigarAndMd_EqualTheSamLineOfBiopython(Case c)
    {
        var result = ResultOf(c);
        var sam = c.Formats["sam:md"]!.GetValue<string>().TrimEnd('\n').Split('\t');
        Assert.That(result.ToCigar(c.T.Length, c.Q.Length), Is.EqualTo(sam[5]));
        string md = sam.First(x => x.StartsWith("MD:Z:", StringComparison.Ordinal))[5..];
        // the alignment as Biopython's SAM writer sees it: the same MD, from the rows alone (no sequence needed)
        Assert.That(result.ToMd(), Is.EqualTo(md));
        var trimmed = c.Trimmed!["sam:md"]!.GetValue<string>().TrimEnd('\n').Split('\t');
        var trim = new AlignmentBridgeOptions { TrimEndGaps = true };
        Assert.That(result.ToCigar(c.T.Length, c.Q.Length, options: trim), Is.EqualTo(trimmed[5]));
        Assert.That(result.ToMd(trim), Is.EqualTo(trimmed.First(x => x.StartsWith("MD:Z:", StringComparison.Ordinal))[5..]));
    }

    [Test]
    public void Oracle_HasTheAnnouncedInventory()
    {
        Assert.That(Oracle.Value, Has.Count.EqualTo(210));
        Assert.That(Oracle.Value.Select(c => c.Mode).Distinct().OrderBy(x => x), Is.EqualTo(new[] { "fitting", "fitting_affine", "global", "global_affine", "local", "local_affine", "overlap" }));
        Assert.That(Oracle.Value.Count(c => c.Trimmed is not null), Is.EqualTo(210));
    }

    [Test]
    public void FreeEndGaps_StayAsGapBlocks_AndTrimToTheSamRecord()
    {
        // R14-1 example: the fitting result keeps the flanks in the rows; Biopython's coordinates keep them as D blocks, a SAM record trims them
        var r = SequenceAligner.SemiGlobalAlign("ACGT", "TTACGTTT");
        var options = new AlignmentBridgeOptions { FirstIsQuery = true };
        Assert.That(r.ToBioCoordAlignment("q", "ACGT", "r", "TTACGTTT", options).Coordinates, Is.EqualTo(new[] { new long[] { 0, 2, 6, 8 }, [0, 0, 4, 4] }));
        Assert.That(r.ToCigar(4, 8, options: options), Is.EqualTo("2D4M2D"));
        var trim = new AlignmentBridgeOptions { FirstIsQuery = true, TrimEndGaps = true };
        Assert.That(r.ToCigar(4, 8, options: trim), Is.EqualTo("4M"));
        Assert.That(r.Format("sam", "q", "ACGT", "r", "TTACGTTT", trim), Is.EqualTo("q\t0\tr\t3\t255\t4M\t*\t0\t0\tACGT\t*\tAS:i:4\n"));
    }

    [Test]
    public void LocalResult_SoftClipsTheQueryFlanks_AndKnowsOnlyTheAlignedResidues()
    {
        var r = SequenceAligner.LocalAlign("TTACGTTT", "GGACGTCC");
        Assert.That(r.ToCigar(8, 8), Is.EqualTo("2S4M2S"));
        var lengthsOnly = r.ToBioCoordAlignment("a", 8, "b", 8);
        Assert.That(lengthsOnly.Coordinates, Is.EqualTo(new[] { new long[] { 2, 6 }, [2, 6] }));
        Assert.That(lengthsOnly.Sequences[0].Slice(2, 6), Is.EqualTo("ACGT"));
        Assert.That(lengthsOnly.Sequences[0].Slice(0, 8), Is.Null);
        Assert.That(r.Format("sam", "a", "TTACGTTT", "b", "GGACGTCC"), Does.Contain("\t3\t255\t2S4M2S\t"));
        Assert.That(r.ToCigar(), Is.EqualTo("4M"));
        Assert.That(r.ToCigar(extended: true), Is.EqualTo("4="));
    }

    [Test]
    public void EmptyAlignment_HasNoCoordinates_AndStarCigar()
    {
        var r = SequenceAligner.LocalAlign("AAAA", "CCCC");
        Assert.That(r.ToCigar(), Is.EqualTo("*"));
        Assert.That(r.ToCigar(4, 4), Is.EqualTo("*"));
        Assert.Throws<InvalidOperationException>(() => r.ToBioCoordAlignment("a", "AAAA", "b", "CCCC"));
        Assert.Throws<InvalidOperationException>(() => r.ToMd());
        Assert.That(r.ToCs(), Is.EqualTo(""));
    }

    [Test]
    public void BadInput_IsRefused()
    {
        var r = SequenceAligner.GlobalAlign("ACGT", "ACGT");
        Assert.Throws<ArgumentException>(() => r.Format("fasta", "a", "ACGT", "b", "ACGT"));
        Assert.Throws<ArgumentException>(() => r.Format("bed13", "a", "ACGT", "b", "ACGT"));
        Assert.Throws<ArgumentException>(() => r.ToBioCoordAlignment("a", "ACG", "b", "ACGT"));
        Assert.Throws<ArgumentException>(() => new AlignmentResult("A-", "A", 0, AlignmentType.Global, 0, 0, 0, 0).ToCigar());
        Assert.Throws<ArgumentException>(() => new AlignmentResult("A-", "A-", 0, AlignmentType.Global, 0, 0, 0, 0).ToCigar());
        Assert.That(AlignmentResultCoordinates.IsFormatName("sam:md"), Is.True);
        Assert.That(AlignmentResultCoordinates.IsFormatName("bed2"), Is.False);
        Assert.That(AlignmentResultCoordinates.IsFormatName(null), Is.False);
    }

    [Test]
    public void ScoreIsTheExactScore()
    {
        var r = SequenceAligner.GlobalAlign("ACGT", "ACGT", new MatrixScoring(SubstitutionMatrix.Get("GONNET1992"), -2.5, -0.5));
        Assert.That(r.ToBioCoordAlignment("a", "ACGT", "b", "ACGT").Attributes["score"]!.GetValue<double>(), Is.EqualTo(r.ExactScore));
    }
}
