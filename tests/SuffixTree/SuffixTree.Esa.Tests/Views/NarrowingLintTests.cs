namespace SuffixTree.Esa.Tests.Views;

/// <summary>
/// The rule of work package S2.4 as a test (see <see cref="NarrowingLint"/>): no cast to <c>int</c> or <c>uint</c> in the files of the search side that is not behind a proof, and the 32-bit loop of S1 is the text that was gated. Written
/// before the implementation: against the S1 code and an empty allow-list it is red, and the failure lists the casts that the package has to remove or to justify.
/// </summary>
[TestFixture]
public class NarrowingLintTests
{
    /// <summary>The casts that stay and why each cannot wrap (key: <see cref="NarrowingLint.Cast.Key"/>). A new entry needs a reason in the log of the package.</summary>
    private static readonly Dictionary<string, string> AllowList = new(StringComparer.Ordinal)
    {
        ["MappedSearch.cs|(int)file.TextReader.Length"] = "the constructor of the 32-bit loop: reached only for SesaFile.IsNarrow, a text of fewer than 2^31 - 1 symbols (the long descent serves every other container)",
        ["MappedSearch.cs|(int)_lcp.Length"] = "the same branch: the LCP array has n + 1 entries, which fit an int for a narrow container",
        ["MappedSearch.cs|(int)((((long)pattern.Length * _bits) + 63) / 64)"] = "the line of S1: the 64-bit words of the packed pattern, at most 8 bits per symbol of an int-length span",
        ["MappedSearch.cs|(int)(((w << (int)(bit & 7)) & 0xFFFF) >> (16 - bits))"] = "Symbol (S1): a 16-bit window shifted into a symbol of at most 8 bits",
        ["MappedSearch.cs|(int)(bit & 7)"] = "Symbol (S1): a shift count, masked to 0..7",
        ["EsaIndex.Search.cs|(int)available"] = "PageOf: taken only when available is below the capacity, which is an int",
        ["EsaIndex.Search.cs|(uint)p"] = "DescendMapped, the 32-bit readable reference of the 32-bit loop (S1): ranks and positions are ints there, tested as unsigned numbers; the long reference is DescendMappedWide",
        ["EsaIndex.Search.cs|(uint)n"] = "DescendMapped, the 32-bit reference: the unsigned bound of an int position",
        ["EsaIndex.Search.cs|(uint)start"] = "DescendMapped, the 32-bit reference: an int position tested as an unsigned number",
        ["EsaIndex.Search.cs|(uint)from"] = "DescendMapped, the 32-bit reference: an int position tested as an unsigned number",
        ["EsaIndex.Search.cs|(int)text.Symbol(p)"] = "a dense symbol of at most 8 bits read from the packed text (in both references)",
        ["EsaIndex.Analysis.cs|(int)n"] = "CountDistinctSubstringsByLength: taken only when n is below maxLength, which is an int",
        ["EsaIndex.Analysis.cs|(int)value"] = "CountDistinctSubstringsByLength: taken only when the value is below m, which is an int",
        ["EsaIndex.Traverse.cs|(int)(position & 63)"] = "MarkLeaf: the shift count of a bit of the position, masked to 0..63",
        ["ChildTable.cs|(int)rank"] = "ChildTable.Rank: taken only when (ulong)rank <= int.MaxValue on the same line (a managed table has 32-bit ranks)",
        ["External/Sesa/SesaFile.cs|(int)(n + 1)"] = "the narrow branch of the constructor (S2.3): n < int.MaxValue is the condition of the branch",
        ["External/Sesa/SesaFile.cs|(int)parsed.Lcp.Param0"] = "the narrow branch: the number of overflow pairs is at most n + 1 (the header reader checks it), which fits an int",
        ["External/Sesa/SesaFile.cs|(int)section.Length"] = "Payload: after the test that refuses a section of more than int.MaxValue bytes (a span cannot hold one)",
        ["IntBridge.cs|(int)length"] = "IntBridge.Length: the checked narrowing itself, after the test of the line above",
        ["IntBridge.cs|(int)value"] = "IntBridge.Value: the checked narrowing itself, after the test of the line above",
    };

    /// <summary>The number of entries that the log of the package says (the list does not grow without a reason in the log).</summary>
    private const int AllowListSize = 20;

    private static List<NarrowingLint.Cast> Scan()
    {
        var casts = new List<NarrowingLint.Cast>();
        foreach (string file in NarrowingLint.Files)
        {
            string source = File.ReadAllText(Path.Combine(NarrowingLint.EsaDirectory(), file));
            if (file == "MappedSearch.cs") source = NarrowingLint.SplitMappedSearch(source).Other;
            casts.AddRange(NarrowingLint.Find(file, source));
        }
        return casts;
    }

    [Test]
    public void TheScanner_FindsACastToIntOrUint_AndIgnoresCommentsStringsTypeofAndCheckedGroups()
    {
        const string source = """
            class C
            {
                int A(long x) => (int)x;                                  // found
                uint B(long x, long y) => (uint)(x + y);                  // found
                int Dd(long x) => checked((int)x);                        // checked: not found
                int E(long[] v, long i) => (int)v[i] + (int) Math.Min(i, 3); // two found
                System.Type F => typeof(int);                             // not a cast
                string G => "(int)x" + '(';                               // a string
                // (int)y in a comment
                /* (uint)z in a block comment */
                int H(long x) => (int)-x;                                 // found, with its unary minus
                int I(long x) => unchecked((int)x);                       // found: unchecked is no proof
            }
            """;
        List<NarrowingLint.Cast> found = [.. NarrowingLint.Find("C.cs", source)];
        Assert.That(found.Select(c => c.Key), Is.EqualTo(new[]
        {
            "C.cs|(int)x", "C.cs|(uint)(x + y)", "C.cs|(int)v[i]", "C.cs|(int)Math.Min(i, 3)", "C.cs|(int)-x", "C.cs|(int)x",
        }));
        Assert.That(found.Select(c => c.Line), Is.EqualTo(new[] { 3, 4, 6, 6, 11, 12 }));
    }

    [Test]
    public void TheThirtyTwoBitLoopOfMappedSearch_IsTheTextThatS1Gated_ByteForByte()
    {
        string source = File.ReadAllText(Path.Combine(NarrowingLint.EsaDirectory(), "MappedSearch.cs"));
        (string pinned, _) = NarrowingLint.SplitMappedSearch(source);
        Assert.That(NarrowingLint.Hash(pinned), Is.EqualTo(NarrowingLint.PinnedLoop),
            "the readers and cells of the codecs, Run, Window, Equal, PatternWindow and PackPattern are the code that the A/B gates measured; a change is a change of the S1 hot path: say so in the log and pin the new text");
    }

    [Test]
    public void EveryCastToIntOrUint_InTheFilesOfTheSearchSide_IsOnTheAllowList()
    {
        List<string> unlisted = [.. Scan().Where(c => !AllowList.ContainsKey(c.Key)).Select(c => $"{c.File}:{c.Line} ({c.Type}){c.Operand}")];
        Assert.That(unlisted, Is.Empty, "a cast that can wrap a 64-bit number: remove it (use the long, or IntBridge, or checked(...)) or put it on the allow-list with the reason why it cannot wrap:\n" + string.Join("\n", unlisted));
    }

    [Test]
    public void EveryEntryOfTheAllowList_MatchesACast_AndGivesAReason_AndTheListHasTheSizeOfTheLog()
    {
        HashSet<string> keys = [.. Scan().Select(c => c.Key)];
        foreach ((string key, string reason) in AllowList)
        {
            Assert.That(keys, Does.Contain(key), $"the allow-list entry '{key}' matches no cast: remove it");
            Assert.That(reason.Length, Is.GreaterThan(20), $"the entry '{key}' needs a reason");
        }
        Assert.That(AllowList, Has.Count.EqualTo(AllowListSize), "the allow-list changed: the log of the package says the number and the reasons");
    }
}
