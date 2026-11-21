using System.Text.Json.Nodes;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>Helpers of the Bio.Align coordinate-format oracle (scripts/alignment/generate_biopython_coord_oracle.py, biopython_coord_oracle.tsv): the JSON description of an alignment and its comparison.</summary>
internal static class BiopythonCoordOracle
{
    internal sealed record Row(string Id, string Format, string Kind, JsonObject Options, string Data1, string Data2);

    internal static List<Row> Load()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.biopython_coord_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], f[2], (JsonObject)JsonNode.Parse(Unesc(f[3]))!, Unesc(f[4]), Unesc(f[5])));
        }

        return list;
    }

    internal static string Unesc(string s) => System.Text.RegularExpressions.Regex.Replace(s, @"\\(.)", m => m.Groups[1].Value switch { "t" => "\t", "n" => "\n", "r" => "\r", var c => c });

    internal static BioCoordAlignment FromSpec(JsonNode spec) => Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.FromJson(spec);

    internal static JsonNode ToSpec(BioCoordAlignment a) => Seqeron.Genomics.Alignment.AlignmentCoordinateFormats.ToJson(a);

    /// <summary>Structural equality in which numbers compare by value (Python prints 1 and 1.0 differently, the description does not distinguish them).</summary>
    internal static bool Same(JsonNode? a, JsonNode? b)
    {
        if (a is null || b is null)
            return a is null && b is null;
        switch (a)
        {
            case JsonObject oa when b is JsonObject ob:
                return oa.Count == ob.Count && oa.All(kv => ob.TryGetPropertyValue(kv.Key, out var v) && Same(kv.Value, v));
            case JsonArray aa when b is JsonArray ab:
                return aa.Count == ab.Count && aa.Zip(ab).All(p => Same(p.First, p.Second));
            case JsonValue va when b is JsonValue vb:
                if (va.TryGetValue<string>(out var sa))
                    return vb.TryGetValue<string>(out var sb) && sa == sb;
                if (va.TryGetValue<bool>(out var ba))
                    return vb.TryGetValue<bool>(out var bb) && ba == bb;
                return !vb.TryGetValue<string>(out _) && ToDouble(va) == ToDouble(vb);
            default:
                return false;
        }
    }

    private static double ToDouble(JsonValue v) => v.TryGetValue<long>(out var l) ? l : v.GetValue<double>();
}
