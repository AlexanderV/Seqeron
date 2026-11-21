using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-10b2 (B09 F215): <c>mash taxscreen</c> of Mash 2.3 (<see cref="MashTaxScreen.Run"/>: the Kraken-style report of the clade and taxon hash counts, the lowest common ancestor of the references that share a hash) against the installed program on
/// 300 random runs (TestData/Alignment/mash_taxscreen_oracle.jsonl.gz, scripts/alignment/generate_mash_taxscreen_oracle.py): synthetic taxonomies in the NCBI dump format (4 - 300 taxa), nucleotide and amino-acid (6-frame) queries, taxonomy IDs in the comment or the
/// mapping file, whole stdout compared as TEXT. Covers the program's own containers whose order decides the output (robin_hood hash counts, the libstdc++ node list that is walked while it grows, std::sort of the children).
/// </summary>
[TestFixture]
public class MashTaxScreen_Oracle_Tests
{
    [Test]
    public void TaxScreen_EqualsMash()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.mash_taxscreen_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        int rows = 0, withTaxa = 0, mapped = 0, large = 0;
        var failures = new List<string>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            var pools = j.GetProperty("files").EnumerateArray().Select(f => MashCommands.ParseFasta(f.GetProperty("name").GetString()!, f.GetProperty("fasta").GetString()!)).ToList();
            var query = MashSketchFile.Read(Convert.FromBase64String(j.GetProperty("qmsh").GetString()!));
            string? map = j.GetProperty("map").ValueKind == JsonValueKind.Null ? null : j.GetProperty("map").GetString();
            string expected = j.GetProperty("out").GetString()!;
            string got;
            try
            {
                got = MashTaxScreen.Run([new MashInput("q.msh", [], query)], pools, new MashTaxonomy(j.GetProperty("names").GetString()!, j.GetProperty("nodes").GetString()!, map));
            }
            catch (Exception e)
            {
                got = "EXCEPTION " + e.GetType().Name + ": " + e.Message;
            }

            if (got != expected)
                failures.Add($"row {j.GetProperty("id")} {string.Join(' ', j.GetProperty("args").EnumerateArray().Select(a => a.GetString()))}\n--- expected\n{expected}--- actual\n{got}");
            rows++;
            if (expected.Count(c => c == '\n') > 1)
                withTaxa++;
            if (map is not null)
                mapped++;
            if (expected.Count(c => c == '\n') > 25)
                large++;
        }

        Assert.That(rows, Is.GreaterThanOrEqualTo(250));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(2)));
        Assert.That(withTaxa, Is.GreaterThanOrEqualTo(200), "reports with at least one taxon");
        Assert.That(mapped, Is.GreaterThanOrEqualTo(40), "runs with a mapping file");
        Assert.That(large, Is.GreaterThanOrEqualTo(8), "reports of more than 25 taxa");
    }
}
