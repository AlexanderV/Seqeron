using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Jobs;
using Seqeron.Genomics;

namespace SuffixTree.Benchmarks;

/// <summary>
/// Benchmarks for SuffixTree operations.
/// Run with: dotnet run -c Release
/// </summary>
[SimpleJob(RuntimeMoniker.Net80)]
[MemoryDiagnoser]
[RankColumn]
public class SuffixTreeBenchmarks
{
    private string _shortText = null!;
    private string _mediumText = null!;
    private string _longText = null!;
    private string _dnaText = null!;

    private SuffixTree _shortTree = null!;
    private SuffixTree _mediumTree = null!;
    private SuffixTree _longTree = null!;
    private SuffixTree _dnaTree = null!;

    private string _shortPattern = null!;
    private string _mediumPattern = null!;
    private string _longPattern = null!;
    private string _dnaPattern = null!;

    [GlobalSetup]
    public void Setup()
    {
        // Short text: ~100 chars
        _shortText = "The quick brown fox jumps over the lazy dog. " +
                     "Pack my box with five dozen liquor jugs.";

        // Medium text: ~10K chars (Lorem ipsum repeated)
        var loremBase = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. " +
                        "Sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. " +
                        "Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris. ";
        _mediumText = string.Concat(Enumerable.Repeat(loremBase, 50));

        // Long text: ~100K chars
        _longText = string.Concat(Enumerable.Repeat(loremBase, 500));

        // DNA sequence: ~50K chars (only ACGT)
        var random = new Random(42);
        var dnaChars = new char[50_000];
        var bases = "ACGT";
        for (int i = 0; i < dnaChars.Length; i++)
        {
            dnaChars[i] = bases[random.Next(4)];
        }
        _dnaText = new string(dnaChars);

        // Build trees
        _shortTree = SuffixTree.Build(_shortText);
        _mediumTree = SuffixTree.Build(_mediumText);
        _longTree = SuffixTree.Build(_longText);
        _dnaTree = SuffixTree.Build(_dnaText);

        // Patterns that exist in texts
        _shortPattern = "fox jumps";
        _mediumPattern = "consectetur adipiscing";
        _longPattern = "tempor incididunt ut labore";
        _dnaPattern = _dnaText.Substring(25000, 20); // 20-char pattern from middle
    }

    #region Build Benchmarks

    [Benchmark]
    [BenchmarkCategory("Build")]
    public SuffixTree Build_Short() => SuffixTree.Build(_shortText);

    [Benchmark]
    [BenchmarkCategory("Build")]
    public SuffixTree Build_Medium() => SuffixTree.Build(_mediumText);

    [Benchmark]
    [BenchmarkCategory("Build")]
    public SuffixTree Build_Long() => SuffixTree.Build(_longText);

    [Benchmark]
    [BenchmarkCategory("Build")]
    public SuffixTree Build_DNA() => SuffixTree.Build(_dnaText);

    #endregion

    #region Contains Benchmarks

    [Benchmark]
    [BenchmarkCategory("Contains")]
    public bool Contains_Short() => _shortTree.Contains(_shortPattern);

    [Benchmark]
    [BenchmarkCategory("Contains")]
    public bool Contains_Medium() => _mediumTree.Contains(_mediumPattern);

    [Benchmark]
    [BenchmarkCategory("Contains")]
    public bool Contains_Long() => _longTree.Contains(_longPattern);

    [Benchmark]
    [BenchmarkCategory("Contains")]
    public bool Contains_DNA() => _dnaTree.Contains(_dnaPattern);

    [Benchmark]
    [BenchmarkCategory("Contains")]
    public bool Contains_NotFound() => _longTree.Contains("xyz123notfound");

    #endregion

    #region FindAllOccurrences Benchmarks

    [Benchmark]
    [BenchmarkCategory("FindAll")]
    public IReadOnlyList<int> FindAll_Short() => _shortTree.FindAllOccurrences(_shortPattern);

    [Benchmark]
    [BenchmarkCategory("FindAll")]
    public IReadOnlyList<int> FindAll_Medium() => _mediumTree.FindAllOccurrences("Lorem");

    [Benchmark]
    [BenchmarkCategory("FindAll")]
    public IReadOnlyList<int> FindAll_Long() => _longTree.FindAllOccurrences("dolor");

    [Benchmark]
    [BenchmarkCategory("FindAll")]
    public IReadOnlyList<int> FindAll_DNA() => _dnaTree.FindAllOccurrences("ACGT");

    #endregion

    #region CountOccurrences Benchmarks

    [Benchmark]
    [BenchmarkCategory("Count")]
    public int Count_Medium() => _mediumTree.CountOccurrences("Lorem");

    [Benchmark]
    [BenchmarkCategory("Count")]
    public int Count_Long() => _longTree.CountOccurrences("dolor");

    [Benchmark]
    [BenchmarkCategory("Count")]
    public int Count_DNA() => _dnaTree.CountOccurrences("ACGT");

    #endregion

    #region LongestRepeatedSubstring Benchmarks

    [Benchmark]
    [BenchmarkCategory("LRS")]
    public string LRS_Short() => _shortTree.LongestRepeatedSubstring();

    [Benchmark]
    [BenchmarkCategory("LRS")]
    public string LRS_Medium() => _mediumTree.LongestRepeatedSubstring();

    [Benchmark]
    [BenchmarkCategory("LRS")]
    public string LRS_Long() => _longTree.LongestRepeatedSubstring();

    [Benchmark]
    [BenchmarkCategory("LRS")]
    public string LRS_DNA() => _dnaTree.LongestRepeatedSubstring();

    #endregion

    #region LongestCommonSubstring Benchmarks

    [Benchmark]
    [BenchmarkCategory("LCS")]
    public string LCS_Short() => _shortTree.LongestCommonSubstring("quick brown fox");

    [Benchmark]
    [BenchmarkCategory("LCS")]
    public string LCS_Medium() => _mediumTree.LongestCommonSubstring("dolor sit amet consectetur");

    [Benchmark]
    [BenchmarkCategory("LCS")]
    public string LCS_Long() => _longTree.LongestCommonSubstring("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed do eiusmod");

    [Benchmark]
    [BenchmarkCategory("LCS")]
    public string LCS_DNA() => _dnaTree.LongestCommonSubstring(_dnaText.Substring(10000, 100));

    #endregion

    #region HasHairpinPotential Benchmarks (Genomics)

    private string _shortPrimer = null!;
    private string _mediumSequence = null!;
    private string _longSequence = null!;

    [GlobalSetup(Target = nameof(Hairpin_ShortPrimer))]
    public void SetupHairpin()
    {
        // Typical primer length (20bp)
        _shortPrimer = "ACGTACGTACGTACGTACGT";

        // Medium sequence (50bp) - below threshold
        _mediumSequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";

        // Long sequence (200bp) - above suffix tree threshold
        var sb = new System.Text.StringBuilder();
        sb.Append("ACGTACGTACGT"); // stem region
        sb.Append(new string('A', 176));
        sb.Append("ACGTACGTACGT"); // complementary stem
        _longSequence = sb.ToString();
    }

    [Benchmark]
    [BenchmarkCategory("Hairpin")]
    public bool Hairpin_ShortPrimer() => PrimerDesigner.HasHairpinPotential(_shortPrimer);

    [Benchmark]
    [BenchmarkCategory("Hairpin")]
    public bool Hairpin_MediumSequence() => PrimerDesigner.HasHairpinPotential(_mediumSequence);

    [Benchmark]
    [BenchmarkCategory("Hairpin")]
    public bool Hairpin_LongSequence() => PrimerDesigner.HasHairpinPotential(_longSequence);

    #endregion
}
