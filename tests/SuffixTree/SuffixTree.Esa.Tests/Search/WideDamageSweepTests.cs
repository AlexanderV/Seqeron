using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Search;

/// <summary>
/// The sweeps of the audit of work package S2.4 (round 1): the long descent over the sections that a container hands to it (<c>SesaFile.WideSections()</c>, the access mode and the stored width included) answers a damaged table exactly as the
/// 32-bit descent of S1 (the reference of the tests, <see cref="EsaIndex.DescendMapped"/>) does, for EVERY value of a cell, not for a sample: a boundary that is one above the right bound, a child as deep as the text plus one, an LCP value
/// of 32 bits that is the length of the text plus one or two. The mutants that moved those bounds by one survived the random damage and the fixed values of <see cref="WideDescentTests"/>.
/// </summary>
[TestFixture]
public unsafe class WideDamageSweepTests
{
    private const string Text = "mississippimississippi";

    private static readonly EsaOpenOptions Checked = new() { Access = EsaAccessMode.Checked };

    private static byte[] Dense(char[] alphabet, string pattern)
    {
        var dense = new byte[pattern.Length];
        for (int i = 0; i < dense.Length; i++)
        {
            int s = Array.IndexOf(alphabet, pattern[i]);
            if (s < 0) return [];
            dense[i] = (byte)s;
        }
        return dense;
    }

    private static string Outcome(Func<(bool Found, long Lo, long Hi)> run)
    {
        try
        {
            (bool found, long lo, long hi) = run();
            return found ? $"found [{lo}..{hi}]" : "absent";
        }
        catch (InvalidDataException)
        {
            return "InvalidDataException";
        }
    }

    private static string[] Patterns(string text, int maxLength)
    {
        var all = new HashSet<string>(StringComparer.Ordinal);
        for (int i = 0; i < text.Length; i++)
        {
            for (int len = 1; len <= Math.Min(maxLength, text.Length - i); len++) all.Add(text.Substring(i, len));
        }
        return [.. all.Order(StringComparer.Ordinal)];
    }

    /// <summary>Opens the damaged image in Checked mode and compares the long descent (over the sections of the container) with the 32-bit reference for every pattern; returns how many patterns ended as <c>InvalidDataException</c> (-1: refused at the open).</summary>
    private static int Compare(byte[] image, bool flush, string[] patterns, bool narrowerBranch, string where)
    {
        using var guarded = new SesaReaderHostileTests.Guarded(image, flush);
        SesaFile file;
        try
        {
            file = SesaFile.FromMemory(guarded.Data, guarded.Length, Checked);
        }
        catch (InvalidDataException)
        {
            return -1;
        }
        using (file)
        {
            WideSearcher wide = WideSearcher.Create(file.WideSections());
            int refused = 0;
            foreach (string pattern in patterns)
            {
                byte[] dense = Dense(file.TextView.Alphabet, pattern);
                string expected = Outcome(() => (EsaIndex.DescendMapped(file.Child.Mapped, file.Sa, file.Lcp, file.TextReader, dense, true, out int lo, out int hi), lo, hi));
                string actual = Outcome(() => (wide.TryMatch(dense, out long lo, out long hi), lo, hi));
                if (narrowerBranch)
                {
                    // a value above 2^31 is a negative number to the 32-bit descent; to the long descent it is a rank or a length that no text of this size has: damage, refused (stricter, never a different answer)
                    Assert.That(actual, Is.EqualTo(expected).Or.EqualTo("InvalidDataException"), $"{where}, pattern '{pattern}'");
                }
                else
                {
                    Assert.That(actual, Is.EqualTo(expected), $"{where}, pattern '{pattern}'");
                }
                if (actual == "InvalidDataException") refused++;
            }
            return refused;
        }
    }

    [TestCase(false)]
    [TestCase(true)]
    public void ACellOfTheChildTable_AtEveryValueAndEveryKind_IsAnsweredAsTheThirtyTwoBitDescentAnswersIt(bool tagged)
    {
        ChildCodec codec = tagged ? ChildCodec.Tagged : ChildCodec.Packed;
        byte[] clean = SesaTestFile.Image(LcpChildCase.Of(Text), codec: codec);
        OracleSection child = SesaOracle.FromBytes(clean).Get(SesaFormat.KindChild);
        int cells = (int)child.Count;
        string[] patterns = Patterns(Text, 10);
        Assert.That(Compare(clean, true, patterns, false, "the clean file"), Is.EqualTo(0));
        int refused = 0;
        for (int cell = 0; cell < cells; cell++)
        {
            for (uint kind = 0; kind < 4; kind++)
            {
                for (uint value = 0; value <= (uint)cells + 1; value++)
                {
                    byte[] image = (byte[])clean.Clone();
                    BitConverter.TryWriteBytes(image.AsSpan((int)child.Offset + (4 * cell)), value | (kind << 30));
                    int r = Compare(image, (cell + value) % 2 == 0, patterns, !tagged && kind >= 2, $"{codec} cell {cell} = {value} | kind {kind}");
                    if (r > 0) refused++;
                }
            }
        }
        Assert.That(refused, Is.GreaterThan(0), "some cell values are refused");
    }

    [Test]
    public void AnLcpValueOf8Bits_FromTheLengthOfTheTextDown_IsAnsweredAsTheThirtyTwoBitDescentAnswersIt()
    {
        LcpChildCase c = LcpChildCase.Of(Text);
        Assert.That(c.Lcp.Kind, Is.EqualTo(LcpKind.U8));
        byte[] clean = SesaTestFile.Image(c);
        OracleSection lcp = SesaOracle.FromBytes(clean).Get(SesaFormat.KindLcp);
        long n = c.N;
        string[] patterns = Patterns(Text, 12);
        int refused = 0;
        for (int rank = 0; rank <= n; rank++)
        {
            foreach (long value in new[] { n - 3, n - 2, n - 1, n, n + 1, n + 2, n + 3, 100, 254 })
            {
                byte[] image = (byte[])clean.Clone();
                image[lcp.Offset + rank] = (byte)value;
                int r = Compare(image, rank % 2 == 0, patterns, false, $"U8 lcp[{rank}] = {value}");
                if (r > 0) refused++;
            }
        }
        Assert.That(refused, Is.GreaterThan(0), "some LCP values are refused");
    }

    [Test]
    public void AnLcpValueOf32Bits_AroundTheLengthOfTheText_IsAnsweredAsTheThirtyTwoBitDescentAnswersIt()
    {
        string text = new('a', 200_000);
        LcpChildCase c = LcpChildCase.Of(text);
        Assert.That(c.Lcp.Kind, Is.EqualTo(LcpKind.U32), "the text is long enough for 32-bit LCP entries");
        byte[] clean = SesaTestFile.Image(c);
        OracleSection lcp = SesaOracle.FromBytes(clean).Get(SesaFormat.KindLcp);
        long n = c.N;
        string[] patterns = [.. new[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100, 1000 }.Select(k => new string('a', k))];
        int refused = 0;
        foreach (int rank in new[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 101, 1001, (int)(n / 2), (int)n })
        {
            foreach (uint value in new[] { (uint)n - 2, (uint)n - 1, (uint)n, (uint)n + 1, (uint)n + 2, (uint)n + 5, 0x7FFFFFFFu, 0x80000000u, 0xFFFFFFFFu })
            {
                byte[] image = (byte[])clean.Clone();
                BitConverter.TryWriteBytes(image.AsSpan((int)lcp.Offset + (4 * rank)), value);
                int r = Compare(image, rank % 2 == 0, patterns, value >= 0x80000000u, $"U32 lcp[{rank}] = {value}");
                if (r > 0) refused++;
            }
        }
        Assert.That(refused, Is.GreaterThan(0), "some LCP values are refused");
    }
}
