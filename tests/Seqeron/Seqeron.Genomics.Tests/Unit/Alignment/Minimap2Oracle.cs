using System.IO.Compression;
using System.Text.Json;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// The oracle files of the minimap2 program tests (B09 WP45): JSON lines <c>{"set", "args", "files", "stdout"}</c> (+ <c>"stderr"</c>, <c>"exit"</c>, <c>"outfiles"</c> for the option / error sets) produced by the installed minimap2 2.26
/// (scripts/alignment/generate_minimap2_{align,pair,options}_oracle.py); a reference text that occurs in several cases is stored once as a <c>{"pool": N, "text"}</c> line and referred to as <c>{"$": N}</c>. A case is checked by running
/// <see cref="Minimap2Program.Run"/> with the same arguments and files and comparing the whole standard output (and, where recorded, the exit status, the standard error and the output files) as text.
/// </summary>
internal static class Minimap2Oracle
{
    internal sealed record Case(string Set, string[] Args, Dictionary<string, string> Files, string StdOut, string? StdErr, int? Exit, Dictionary<string, string> OutFiles, bool StdErrRecorded);

    internal static List<Case> Read(string file)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment." + file;
        using var stream = typeof(Minimap2Oracle).Assembly.GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var pool = new Dictionary<int, string>();
        var cases = new List<Case>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            using var doc = JsonDocument.Parse(line);
            var r = doc.RootElement;
            if (r.TryGetProperty("pool", out var id))
            {
                pool[id.GetInt32()] = r.GetProperty("text").GetString()!;
                continue;
            }

            var files = r.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.ValueKind == JsonValueKind.Object ? pool[p.Value.GetProperty("$").GetInt32()] : p.Value.GetString()!);
            bool hasExit = r.TryGetProperty("exit", out var exit);
            string? stderr = hasExit && r.GetProperty("stderr").ValueKind == JsonValueKind.String ? r.GetProperty("stderr").GetString() : null;
            var outFiles = hasExit ? r.GetProperty("outfiles").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!) : [];
            cases.Add(new Case(r.GetProperty("set").GetString()!, [.. r.GetProperty("args").EnumerateArray().Select(e => e.GetString()!)], files, r.GetProperty("stdout").GetString()!, stderr, hasExit ? exit.GetInt32() : null, outFiles, hasExit));
        }

        return cases;
    }

    /// <summary>Run every case of the set and return the descriptions of those whose result differs from the program's.</summary>
    internal static List<string> Failures(IEnumerable<Case> cases)
    {
        var failures = new List<string>();
        foreach (var c in cases)
        {
            string? problem = null;
            try
            {
                var result = Minimap2Program.Run(c.Args, c.Files);
                if (result.StdOut != c.StdOut)
                    problem = $"stdout\n   expected {Show(c.StdOut)}\n   actual   {Show(result.StdOut)}";
                else if (c.StdErrRecorded && result.ExitCode != c.Exit)
                    problem = $"exit {result.ExitCode}, expected {c.Exit}";
                else if (c.StdErr is not null && result.StdErr != c.StdErr)
                    problem = $"stderr\n   expected {Show(c.StdErr)}\n   actual   {Show(result.StdErr)}";
                else if (c.StdErrRecorded && !c.OutFiles.OrderBy(k => k.Key).SequenceEqual(result.Files.OrderBy(k => k.Key)))
                    problem = "output files differ";
            }
            catch (Exception ex) when (ex is not OperationCanceledException)
            {
                problem = $"{ex.GetType().Name}: {ex.Message}";
            }

            if (problem is not null)
                failures.Add($"{string.Join(' ', c.Args)}: {problem}");
        }

        return failures;
    }

    private static string Show(string s)
    {
        string t = s.Replace("\n", " || ").Replace("\r", "<CR>");
        return t.Length > 600 ? t[..600] + "..." : t;
    }
}
