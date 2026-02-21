using System.Diagnostics;
using System.IO.Compression;
using System.Text;
using SuffixTree;
using SuffixTree.Persistent;

/// <summary>
/// Focused profiler for persistent suffix tree build.
/// Reads a subset of chr1.fa.gz and measures build with GC/memory diagnostics.
/// Usage: dotnet run --project apps\SuffixTree.GenomeDemo -- --profile [chars]
/// </summary>
static class ProfileBuild
{
    internal static int Run(string[] args)
    {
        int charLimit = 1_000_000; // default 1M
        if (args.Length > 1 && int.TryParse(args[1], out int n))
            charLimit = n;

        var dataDir = Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..", "data"));
        var fastaGz = Path.Combine(dataDir, "chr1.fa.gz");
        var treePath = Path.Combine(Path.GetTempPath(), $"ST_profile_{charLimit}.dat");

        if (!File.Exists(fastaGz))
        {
            Console.Error.WriteLine($"FASTA not found: {fastaGz}");
            return 1;
        }

        Console.WriteLine($"Profile: building persistent suffix tree from first {charLimit:N0} chars of chr1");
        Console.WriteLine($"Tree file: {treePath}");
        Console.WriteLine();

        // ── Parse FASTA (subset) ──
        var sw = Stopwatch.StartNew();
        var sb = new StringBuilder(charLimit + 1024);
        using (var fs = File.OpenRead(fastaGz))
        using (var gz = new GZipStream(fs, CompressionMode.Decompress))
        using (var reader = new StreamReader(gz, Encoding.ASCII))
        {
            string? line;
            while ((line = reader.ReadLine()) != null && sb.Length < charLimit)
            {
                if (line.Length == 0 || line[0] == '>') continue;
                int take = Math.Min(line.Length, charLimit - sb.Length);
                sb.Append(line, 0, take);
            }
        }
        var sequence = sb.ToString();
        sb = null;
        var parseTime = sw.Elapsed;
        Console.WriteLine($"[Parse]   {parseTime.TotalMilliseconds,10:F1} ms  ({sequence.Length:N0} chars)");

        // ── Force GC baseline ──
        GC.Collect(2, GCCollectionMode.Aggressive, true, true);
        GC.WaitForPendingFinalizers();
        GC.Collect(2, GCCollectionMode.Aggressive, true, true);
        long memBefore = GC.GetTotalMemory(true);
        long gen0Before = GC.CollectionCount(0);
        long gen1Before = GC.CollectionCount(1);
        long gen2Before = GC.CollectionCount(2);

        // ── Build ──
        if (File.Exists(treePath)) File.Delete(treePath);
        var textSource = new StringTextSource(sequence);

        sw.Restart();
        ISuffixTree tree;
        try
        {
            tree = PersistentSuffixTreeFactory.Create(textSource, treePath);
        }
        catch (Exception ex)
        {
            Console.Error.WriteLine($"Build FAILED: {ex.Message}");
            Console.Error.WriteLine(ex.StackTrace);
            return 2;
        }
        var buildTime = sw.Elapsed;

        long memAfter = GC.GetTotalMemory(false);
        long gen0After = GC.CollectionCount(0);
        long gen1After = GC.CollectionCount(1);
        long gen2After = GC.CollectionCount(2);

        long treeFileSize = File.Exists(treePath) ? new FileInfo(treePath).Length : 0;

        Console.WriteLine($"[Build]   {buildTime.TotalMilliseconds,10:F1} ms");
        Console.WriteLine();
        Console.WriteLine("── Results ──────────────────────────────────────");
        Console.WriteLine($"  Text length:     {sequence.Length,15:N0} chars");
        Console.WriteLine($"  Tree file:       {treeFileSize,15:N0} bytes ({treeFileSize / 1_048_576.0:F1} MB)");
        Console.WriteLine($"  Nodes:           {tree.NodeCount,15:N0}");
        Console.WriteLine($"  Managed heap Δ:  {(memAfter - memBefore),15:N0} bytes ({(memAfter - memBefore) / 1_048_576.0:F1} MB)");
        Console.WriteLine($"  GC gen0:         {gen0After - gen0Before,15:N0} collections");
        Console.WriteLine($"  GC gen1:         {gen1After - gen1Before,15:N0} collections");
        Console.WriteLine($"  GC gen2:         {gen2After - gen2Before,15:N0} collections");
        Console.WriteLine($"  Build rate:      {sequence.Length / buildTime.TotalSeconds,15:N0} chars/sec");
        Console.WriteLine($"  Bytes/char:      {(treeFileSize > 0 ? (double)treeFileSize / sequence.Length : 0),15:F1}");
        Console.WriteLine();

        // ── Quick query sanity check ──
        sw.Restart();
        bool found = tree.Contains("GATTACA");
        Console.WriteLine($"[Query]   Contains(\"GATTACA\") = {found}  ({sw.Elapsed.TotalMilliseconds:F3} ms)");

        sw.Restart();
        int count = tree.CountOccurrences("ATG");
        Console.WriteLine($"[Query]   CountOccurrences(\"ATG\") = {count:N0}  ({sw.Elapsed.TotalMilliseconds:F3} ms)");

        sw.Restart();
        string lrs = tree.LongestRepeatedSubstring();
        Console.WriteLine($"[Query]   LRS length = {lrs.Length:N0}  ({sw.Elapsed.TotalMilliseconds:F3} ms)");

        // ── Cleanup ──
        (tree as IDisposable)?.Dispose();
        try { File.Delete(treePath); } catch { }

        Console.WriteLine();
        Console.WriteLine($"Total: {parseTime + buildTime}");
        return 0;
    }
}
