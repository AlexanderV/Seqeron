#if NET9_0_OR_GREATER
using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Jobs;
using SuffixTree.Persistent;

namespace SuffixTree.Benchmarks;

/// <summary>
/// Benchmarks for PersistentSuffixTree operations (memory-mapped file storage).
/// Measures build, load, and query performance with unsafe pointer access,
/// capacity pre-estimation, and PrefetchVirtualMemory optimizations.
///
/// Run:  dotnet run -c Release -f net9.0 -- --filter "*Persistent*"
/// </summary>
[SimpleJob(RuntimeMoniker.Net90)]
[MemoryDiagnoser]
[RankColumn]
public class PersistentSuffixTreeBenchmarks
{
    private string _mediumText = null!;
    private string _longText = null!;
    private string _dnaText = null!;

    private string _mediumPattern = null!;
    private string _longPattern = null!;
    private string _dnaPattern = null!;

    // Pre-built trees for query benchmarks
    private ISuffixTree _mmfMediumTree = null!;
    private ISuffixTree _mmfLongTree = null!;
    private ISuffixTree _mmfDnaTree = null!;

    private ISuffixTree _heapMediumTree = null!;
    private ISuffixTree _heapDnaTree = null!;

    // Pre-built files for Load benchmarks
    private string _loadFile = null!;

    // Fixed paths for Build benchmarks — reused across invocations so that
    // Windows Defender scans them once during GlobalSetup, not on every iteration.
    private string _buildMediumFile = null!;
    private string _buildDnaFile = null!;
    private string _buildLongFile = null!;

    private readonly List<string> _tempFiles = [];

    private string AllocTempFile()
    {
        var path = Path.Combine(Path.GetTempPath(), $"ST_bench_{Guid.NewGuid():N}.dat");
        _tempFiles.Add(path);
        return path;
    }

    private string AllocFixedFile(string name)
    {
        var path = Path.Combine(Path.GetTempPath(), $"ST_bench_{name}.dat");
        _tempFiles.Add(path);
        return path;
    }

    [GlobalSetup]
    public void Setup()
    {
        var loremBase = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. " +
                        "Sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. " +
                        "Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris. ";
        _mediumText = string.Concat(Enumerable.Repeat(loremBase, 50));   // ~10K chars
        _longText = string.Concat(Enumerable.Repeat(loremBase, 500));    // ~100K chars

        var random = new Random(42);
        var dnaChars = new char[50_000];
        var bases = "ACGT";
        for (int i = 0; i < dnaChars.Length; i++)
            dnaChars[i] = bases[random.Next(4)];
        _dnaText = new string(dnaChars);

        _mediumPattern = "consectetur adipiscing";
        _longPattern = "tempor incididunt ut labore";
        _dnaPattern = _dnaText.Substring(25000, 20);

        // Pre-build MMF trees for query benchmarks
        _mmfMediumTree = PersistentSuffixTreeFactory.Create(new StringTextSource(_mediumText), AllocTempFile());
        _mmfLongTree = PersistentSuffixTreeFactory.Create(new StringTextSource(_longText), AllocTempFile());
        _mmfDnaTree = PersistentSuffixTreeFactory.Create(new StringTextSource(_dnaText), AllocTempFile());

        // Pre-build Heap trees for comparison
        _heapMediumTree = PersistentSuffixTreeFactory.Create(new StringTextSource(_mediumText));
        _heapDnaTree = PersistentSuffixTreeFactory.Create(new StringTextSource(_dnaText));

        // Pre-build file for Load benchmarks
        _loadFile = AllocTempFile();
        using (PersistentSuffixTreeFactory.Create(new StringTextSource(_dnaText), _loadFile) as IDisposable ?? throw new InvalidOperationException())
        {
        }

        // Fixed paths for Build benchmarks — create once so Windows Defender
        // scans them during setup, not during measured iterations.
        _buildMediumFile = AllocFixedFile("build_medium");
        _buildDnaFile = AllocFixedFile("build_dna");
        _buildLongFile = AllocFixedFile("build_long");

        // Warm the paths: create and dispose a small tree to trigger any AV scan now.
        foreach (var path in new[] { _buildMediumFile, _buildDnaFile, _buildLongFile })
        {
            using var warmup = PersistentSuffixTreeFactory.Create(new StringTextSource("warmup"), path) as IDisposable;
        }
    }

    [GlobalCleanup]
    public void Cleanup()
    {
        (_mmfMediumTree as IDisposable)?.Dispose();
        (_mmfLongTree as IDisposable)?.Dispose();
        (_mmfDnaTree as IDisposable)?.Dispose();
        (_heapMediumTree as IDisposable)?.Dispose();
        (_heapDnaTree as IDisposable)?.Dispose();

        foreach (var f in _tempFiles)
        {
            try { File.Delete(f); } catch { /* cleanup */ }
        }
    }

    // ──────────────── Build: MMF vs Heap ────────────────

    [Benchmark]
    [BenchmarkCategory("PersistentBuild")]
    public ISuffixTree Build_Mmf_Medium()
    {
        var tree = PersistentSuffixTreeFactory.Create(new StringTextSource(_mediumText), _buildMediumFile);
        (tree as IDisposable)?.Dispose();
        return tree;
    }

    [Benchmark]
    [BenchmarkCategory("PersistentBuild")]
    public ISuffixTree Build_Heap_Medium()
    {
        var tree = PersistentSuffixTreeFactory.Create(new StringTextSource(_mediumText));
        (tree as IDisposable)?.Dispose();
        return tree;
    }

    [Benchmark]
    [BenchmarkCategory("PersistentBuild")]
    public ISuffixTree Build_Mmf_DNA()
    {
        var tree = PersistentSuffixTreeFactory.Create(new StringTextSource(_dnaText), _buildDnaFile);
        (tree as IDisposable)?.Dispose();
        return tree;
    }

    [Benchmark]
    [BenchmarkCategory("PersistentBuild")]
    public ISuffixTree Build_Heap_DNA()
    {
        var tree = PersistentSuffixTreeFactory.Create(new StringTextSource(_dnaText));
        (tree as IDisposable)?.Dispose();
        return tree;
    }

    [Benchmark]
    [BenchmarkCategory("PersistentBuild")]
    public ISuffixTree Build_Mmf_Long()
    {
        var tree = PersistentSuffixTreeFactory.Create(new StringTextSource(_longText), _buildLongFile);
        (tree as IDisposable)?.Dispose();
        return tree;
    }

    // ──────────────── Load from file (with Prefetch) ────────────────

    [Benchmark]
    [BenchmarkCategory("PersistentLoad")]
    public ISuffixTree Load_Mmf_DNA()
    {
        var tree = PersistentSuffixTreeFactory.Load(_loadFile);
        (tree as IDisposable)?.Dispose();
        return tree;
    }

    // ──────────────── Contains: MMF vs Heap ────────────────

    [Benchmark(Baseline = true)]
    [BenchmarkCategory("PersistentContains")]
    public bool Contains_Mmf_Medium() => _mmfMediumTree.Contains(_mediumPattern);

    [Benchmark]
    [BenchmarkCategory("PersistentContains")]
    public bool Contains_Heap_Medium() => _heapMediumTree.Contains(_mediumPattern);

    [Benchmark]
    [BenchmarkCategory("PersistentContains")]
    public bool Contains_Mmf_DNA() => _mmfDnaTree.Contains(_dnaPattern);

    [Benchmark]
    [BenchmarkCategory("PersistentContains")]
    public bool Contains_Heap_DNA() => _heapDnaTree.Contains(_dnaPattern);

    [Benchmark]
    [BenchmarkCategory("PersistentContains")]
    public bool Contains_Mmf_Long() => _mmfLongTree.Contains(_longPattern);

    // ──────────────── Count: MMF ────────────────

    [Benchmark]
    [BenchmarkCategory("PersistentCount")]
    public int Count_Mmf_DNA() => _mmfDnaTree.CountOccurrences("ACGT");

    [Benchmark]
    [BenchmarkCategory("PersistentCount")]
    public int Count_Mmf_Long() => _mmfLongTree.CountOccurrences("dolor");

    // ──────────────── FindAll: MMF ────────────────

    [Benchmark]
    [BenchmarkCategory("PersistentFindAll")]
    public IReadOnlyList<int> FindAll_Mmf_DNA() => _mmfDnaTree.FindAllOccurrences("ACGT");

    [Benchmark]
    [BenchmarkCategory("PersistentFindAll")]
    public IReadOnlyList<int> FindAll_Mmf_Long() => _mmfLongTree.FindAllOccurrences("dolor");

    // ──────────────── LRS: MMF ────────────────

    [Benchmark]
    [BenchmarkCategory("PersistentLRS")]
    public string LRS_Mmf_Medium() => _mmfMediumTree.LongestRepeatedSubstring();

    [Benchmark]
    [BenchmarkCategory("PersistentLRS")]
    public string LRS_Mmf_DNA() => _mmfDnaTree.LongestRepeatedSubstring();
}
#endif
