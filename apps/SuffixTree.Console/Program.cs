using System.Diagnostics;
using System.Text;
using SuffixTree;
using SuffixTree.Persistent;

var text = BuildText(length: 120_000);
var tempFile = Path.Combine(Path.GetTempPath(), $"suffix-tree-persistent-perf-{Guid.NewGuid():N}.tree");

try
{
    var buildSw = Stopwatch.StartNew();
    using (PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text), tempFile))
    {
        // build + dispose only — the point is to measure persistent-build wall time
    }
    buildSw.Stop();

    using var tree = PersistentSuffixTreeFactory.LoadPersistent(tempFile);

    string existingContains = "TTGACCAT";
    string missingContains = "ZZZZZZZZ";
    string multiHitPattern = "ACGTAC";

    // Warm-up
    _ = tree.Contains(existingContains);
    _ = tree.Contains(missingContains);
    _ = tree.CountOccurrences(multiHitPattern);
    _ = tree.FindAllOccurrences(multiHitPattern);
    _ = tree.LongestRepeatedSubstring();
    _ = tree.LongestRepeatedSubstring();

    var metrics = new Dictionary<string, double>
    {
        ["build_ms"] = buildSw.Elapsed.TotalMilliseconds,
        ["contains_hit_us"] = Measure(40_000, () => tree.Contains(existingContains)),
        ["contains_miss_us"] = Measure(40_000, () => tree.Contains(missingContains)),
        ["count_us"] = Measure(30_000, () => tree.CountOccurrences(multiHitPattern)),
        ["find_all_us"] = Measure(2_000, () => _ = tree.FindAllOccurrences(multiHitPattern).Count),
        ["lrs_first_us"] = MeasureColdLrs(tempFile, text),
        ["lrs_cached_us"] = Measure(50_000, () => _ = tree.LongestRepeatedSubstring().Length),
    };

    foreach (var metric in metrics.OrderBy(pair => pair.Key, StringComparer.Ordinal))
    {
        Console.WriteLine($"{metric.Key}={metric.Value:F4}");
    }
}
finally
{
    Cleanup(tempFile);
}

static double Measure(int iterations, Action action)
{
    var sw = Stopwatch.StartNew();
    for (int i = 0; i < iterations; i++)
    {
        action();
    }

    sw.Stop();
    return sw.Elapsed.TotalMilliseconds * 1000.0 / iterations;
}

static double MeasureColdLrs(string filePath, string expectedText)
{
    const int Iterations = 25;
    var sw = Stopwatch.StartNew();
    for (int i = 0; i < Iterations; i++)
    {
        using var tree = PersistentSuffixTreeFactory.LoadPersistent(filePath);
        var lrs = tree.LongestRepeatedSubstring();
        if (lrs.Length == 0 || !expectedText.Contains(lrs, StringComparison.Ordinal))
        {
            throw new InvalidOperationException("Unexpected LRS result during benchmark.");
        }
    }

    sw.Stop();
    return sw.Elapsed.TotalMilliseconds * 1000.0 / Iterations;
}

static string BuildText(int length)
{
    const string motifA = "ACGTACGTTTGACCATGGAACCTA";
    const string motifB = "GATTACAGGCTTAGGACCTTACGA";
    const string motifC = "TTGACCATACGTACGGATTACAAC";

    var builder = new StringBuilder(length + 64);
    int seed = 17;

    while (builder.Length < length)
    {
        seed = unchecked(seed * 1103515245 + 12345);
        builder.Append((seed & 1) == 0 ? motifA : motifB);
        builder.Append(motifC);
        builder.Append((char)('A' + ((seed >> 8) & 0x0F) % 20));
        builder.Append((char)('A' + ((seed >> 13) & 0x0F) % 20));
    }

    return builder.ToString(0, length);
}

static void Cleanup(string treePath)
{
    TryDelete(treePath);
    TryDelete(treePath + ".children.tmp");
    TryDelete(treePath + ".depth.tmp");
}

static void TryDelete(string path)
{
    try
    {
        if (File.Exists(path))
        {
            File.Delete(path);
        }
    }
    catch (IOException)
    {
        // best-effort temp cleanup; ignore if the file is locked/gone
    }
    catch (UnauthorizedAccessException)
    {
        // best-effort temp cleanup; ignore if access is denied
    }
}
