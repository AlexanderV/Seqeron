using System.IO.Compression;
using System.Text.Json;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>Reads an embedded gzip-compressed JSON-lines oracle of the MUMmer programs (TestData/Alignment/*.jsonl.gz).</summary>
public static class MummerOracleData
{
    public static IEnumerable<JsonElement> Rows(string resource)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment." + resource;
        string? dir = Environment.GetEnvironmentVariable("SEQERON_ORACLE_DIR");
        using Stream stream = dir != null
            ? File.OpenRead(Path.Combine(dir, resource))
            : System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        string? line;
        while ((line = reader.ReadLine()) is not null)
            yield return JsonDocument.Parse(line).RootElement;
    }

    public static string Str(this JsonElement e, string property) => e.GetProperty(property).GetString()!;

    public static string? StrOrNull(this JsonElement e, string property) => e.TryGetProperty(property, out var v) && v.ValueKind == JsonValueKind.String ? v.GetString() : null;

    public static List<string> StrList(this JsonElement e, string property) => e.GetProperty(property).EnumerateArray().Select(x => x.GetString()!).ToList();

    /// <summary>The position and a short context of the first difference of two texts.</summary>
    public static string Diff(string ours, string expected)
    {
        int i = 0;
        while (i < Math.Min(ours.Length, expected.Length) && ours[i] == expected[i])
            i++;
        static string Ctx(string s, int at) => at >= s.Length ? "<end>" : s.Substring(Math.Max(0, at - 40), Math.Min(s.Length - Math.Max(0, at - 40), 100)).Replace("\n", "\\n", StringComparison.Ordinal);
        return $"@{i}: ours '{Ctx(ours, i)}' expected '{Ctx(expected, i)}'";
    }
}
