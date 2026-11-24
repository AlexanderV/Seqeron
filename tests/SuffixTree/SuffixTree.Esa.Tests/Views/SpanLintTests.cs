namespace SuffixTree.Esa.Tests.Views;

/// <summary>
/// Risk R8 of ADR-005 (work package S2.3): no span is made from a pointer or a section with a narrowed length. The scanner (<see cref="SpanLint"/>) is pinned on synthetic sources first (what it finds, what it does not), then run over
/// <c>src/SuffixTree/Algorithms/SuffixTree.Esa</c>: every creation whose length argument narrows to an integer or names a section must be on <see cref="AllowList"/> with its reason, and every entry of the list must still match.
/// </summary>
[TestFixture]
public class SpanLintTests
{
    /// <summary>
    /// The accounted-for sites: file, call, length expression, and why the length cannot exceed <c>int.MaxValue</c> (or why a longer one is refused). A site is identified by its text, not by its line, so an edit elsewhere in the file does
    /// not disturb the list; a second site with the same text is the same entry (the count says how many).
    /// </summary>
    internal static readonly (string File, string Call, string Length, int Count, string Reason)[] AllowList = [];

    private static string Src(string code) => "class C { void M() { " + code + " } }";

    private static IEnumerable<SpanLint.Site> Hazards(string code) => SpanLint.Hazards("t.cs", Src(code));

    // ---------------------------------------------------------------------------------------------- the scanner

    [Test]
    public void TheScanner_FindsACastToInt_InEveryFormOfSpanCreation()
    {
        Assert.That(Hazards("var a = new Span<byte>(p, (int)length);").Select(s => s.Length), Is.EqualTo(new[] { "(int)length" }));
        Assert.That(Hazards("var a = new ReadOnlySpan<byte>(p + offset, (int)(end - start));").Count(), Is.EqualTo(1));
        Assert.That(Hazards("var a = new ReadOnlySpan<int>(array, from, (int)count);").Select(s => s.Length), Is.EqualTo(new[] { "(int)count" }));
        Assert.That(Hazards("var a = buffer.AsSpan(0, (int)take);").Count(), Is.EqualTo(1));
        Assert.That(Hazards("var a = MemoryMarshal.CreateSpan(ref r, (int)n);").Count(), Is.EqualTo(1));
        Assert.That(Hazards("var a = MemoryMarshal.CreateReadOnlySpan(ref r, (uint)n);").Count(), Is.EqualTo(1));
        Assert.That(Hazards("var a = new Span<T>(_p, checked((int)length));").Count(), Is.EqualTo(1));
        Assert.That(Hazards("var a = new ReadOnlySpan<byte>(data + at, checked(1 + n));").Count(), Is.EqualTo(1), "checked(...) is a narrowing in disguise");
    }

    [Test]
    public void TheScanner_FindsALengthThatComesFromASection()
    {
        Assert.That(Hazards("var a = new ReadOnlySpan<byte>(_data + s.Offset, section.Length);").Count(), Is.EqualTo(1));
        Assert.That(Hazards("var a = new ReadOnlySpan<byte>(p, Section.Length);").Count(), Is.EqualTo(1));
        Assert.That(Hazards("var a = x.AsSpan(0, section.Count);").Count(), Is.EqualTo(1));
    }

    [Test]
    public void TheScanner_IgnoresABoundedLength_AndTheCallsThatHaveNone()
    {
        Assert.That(Hazards("var a = new ReadOnlySpan<byte>(p, 8);"), Is.Empty);
        Assert.That(Hazards("var a = new Span<byte>(buffer, _blockBytes);"), Is.Empty);
        Assert.That(Hazards("var a = new ReadOnlySpan<ushort>(p + alphabet.Offset, sigma);"), Is.Empty);
        Assert.That(Hazards("var a = array.AsSpan(start);"), Is.Empty);
        Assert.That(Hazards("var a = array.AsSpan();"), Is.Empty);
        Assert.That(Hazards("var a = new Span<byte>(array);"), Is.Empty);
        Assert.That(Hazards("var a = text.AsSpan(1, 5);"), Is.Empty);
        Assert.That(Hazards("var a = (int)length; var b = new Span<byte>(p, 4);"), Is.Empty, "a cast elsewhere on the line is not the length");
        Assert.That(Hazards("// new Span<byte>(p, (int)n)").Count(), Is.EqualTo(1), "a comment is scanned like code: say so in the allow-list or remove it");
    }

    [Test]
    public void TheScanner_ReadsNestedAndMultiLineArguments_AndStringsWithParentheses()
    {
        var sites = Hazards("var a = new ReadOnlySpan<byte>(\n    file.Data + end,\n    checked((int)(s.Offset - end)))\n    .IndexOfAnyExcept((byte)0);").ToList();
        Assert.That(sites, Has.Count.EqualTo(1));
        Assert.That(sites[0].Length, Is.EqualTo("checked((int)(s.Offset - end))"));
        Assert.That(Hazards("var a = Math.Max(1, new Span<byte>(p, (int)n).Length);").Count(), Is.EqualTo(1));
        Assert.That(Hazards("var a = new Span<byte>(F(\"a,(b\", 1), (int)n);").Select(s => s.Length), Is.EqualTo(new[] { "(int)n" }), "a comma and a parenthesis inside a string literal");
        Assert.That(Hazards("var a = new ReadOnlySpan<byte>(p, ").Count(), Is.EqualTo(0), "an unterminated call is not a site");
    }

    [Test]
    public void TheScanner_NamesTheCall_TheLine_AndTheKey()
    {
        SpanLint.Site site = SpanLint.Hazards("src/X.cs", "\n\nclass C { object M() => new ReadOnlySpan<int>(cells + done, (int)piece); }").Single();
        Assert.That(site.File, Is.EqualTo("src/X.cs"));
        Assert.That(site.Line, Is.EqualTo(3));
        Assert.That(site.Call, Is.EqualTo("new ReadOnlySpan<int>"));
        Assert.That(site.Length, Is.EqualTo("(int)piece"));
        Assert.That(site.Key, Is.EqualTo("src/X.cs|new ReadOnlySpan<int>|(int)piece"));
        Assert.That(SpanLint.Hazards("a.cs", "x.AsSpan(0, (int)n)").Single().Call, Is.EqualTo(".AsSpan"));
    }

    // ---------------------------------------------------------------------------------------------- the project

    private static string EsaDirectory => Path.Combine(SpanLint.RepositoryRoot(), "src", "SuffixTree", "Algorithms", "SuffixTree.Esa");

    [Test]
    public void TheEsaProject_HasNoSpanWithANarrowedOrSectionLength_ThatIsNotOnTheAllowList()
    {
        List<SpanLint.Site> sites = SpanLint.ScanDirectory(EsaDirectory);
        var allowed = AllowList.ToDictionary(e => $"{e.File}|{e.Call}|{e.Length}", e => e.Count);
        string[] unlisted = [.. sites.GroupBy(s => s.Key).Where(g => !allowed.TryGetValue(g.Key, out int c) || c != g.Count()).Select(g => $"{g.Key} x{g.Count()} (line {g.First().Line})")];
        Assert.That(unlisted, Is.Empty, "a span is made with a length that is a cast to int or comes from a section: bound it, refuse the length that does not fit, or add the site to SpanLintTests.AllowList with the reason why it is bounded");
    }

    [Test]
    public void TheAllowList_HasNoEntryThatNoLongerMatches_AndEveryEntryHasAReason()
    {
        List<SpanLint.Site> sites = SpanLint.ScanDirectory(EsaDirectory);
        var found = sites.GroupBy(s => s.Key).ToDictionary(g => g.Key, g => g.Count());
        foreach ((string file, string call, string length, int count, string reason) in AllowList)
        {
            string key = $"{file}|{call}|{length}";
            Assert.That(found.TryGetValue(key, out int actual) ? actual : 0, Is.EqualTo(count), $"the allow-list entry {key} (stale, or the number of sites changed)");
            Assert.That(reason.Length, Is.GreaterThan(30), $"{key}: the reason says why the length is bounded");
            Assert.That(count, Is.GreaterThan(0));
        }
        Assert.That(AllowList.Select(e => (e.File, e.Call, e.Length)).Distinct().Count(), Is.EqualTo(AllowList.Length), "one entry per site text");
    }

    [Test]
    public void TheProjectScan_FindsTheSitesThatItMustFind()
    {
        // the scan is not vacuous: the container reader and the verifier contain spans over sections
        List<SpanLint.Site> all = Directory.EnumerateFiles(EsaDirectory, "*.cs", SearchOption.AllDirectories)
            .Where(p => !p.Contains("/obj/") && !p.Contains("/bin/"))
            .SelectMany(p => SpanLint.Find(Path.GetRelativePath(EsaDirectory, p), File.ReadAllText(p))).ToList();
        Assert.That(all.Count, Is.GreaterThan(50), "the project has more than fifty span creations");
        Assert.That(all.Any(s => s.File.EndsWith("SesaVerifier.cs", StringComparison.Ordinal) && s.Call == "new ReadOnlySpan<byte>"), Is.True);
        Assert.That(all.Any(s => s.File.EndsWith("MappedTextReader.cs", StringComparison.Ordinal)), Is.True);
    }
}
