using System.Diagnostics;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// Work package S2.1, deliverable F: the red suites of the packages S2.3 - S2.8, only those that compile against the S1 head and fail because of a limit. Each test opens a sparse container of 2^31 + 10 symbols (POS_WIDTH 4) or
/// 2^32 + 10 symbols (POS_WIDTH 5, ADR-005 section 6.1) through an entry point that exists today and asserts the TARGET behaviour of stage S2; the property <c>TurnsGreenAt</c> names the package that makes it pass. A test that
/// needs an API that does not exist yet (<c>OpenIndex</c>, <c>IPosReader</c>) is not here: each later package writes its own tests first. The tests marked "guard" pass today and must keep passing (the behaviour that S2 must not change).
/// <para>The containers are shells (zeros and the ends of <c>a^n</c>, see <see cref="SparseContainer"/>): valid at the <c>Header</c> and <c>Sections</c> levels, damaged at the <c>Structure</c> level, which is what the tests rely on.</para>
/// Category <c>Large</c>, explicit: <c>dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter FullyQualifiedName~LargeRedTests</c>.
/// </summary>
[TestFixture]
[Category("Large")]
[Explicit(ExplicitReason)]
public class LargeRedTests : LargeFixtureBase
{
    private const long TwoTo31 = 1L << 31;
    private const long TwoTo32 = 1L << 32;

    /// <summary>The time in which a container opens at the Header level (ADR-005 row S2.1: 5 s).</summary>
    private const double HeaderBudgetSeconds = 5;

    /// <summary>
    /// The time in which the Sections level of the 78 GB shell of 2^32 + 10 symbols must finish: the checksum pass reads 78 GB of holes at 3.6 - 4.6 GB/s (measured on the S1 reader below 2^31: 7.2 s for 26 GB), 17 - 22 s;
    /// the 5 s of the ADR row cannot be met by a serial pass over that many bytes (see the record of S2.1), so the budget is 30 s.
    /// </summary>
    private const double SectionsBudgetSeconds = 30;

    /// <summary>The anonymous memory that opening may add to the process (ADR-005 row S2.1: RssAnon below 100 MB).</summary>
    private const double AnonBudgetMb = 100;

    private string Shell(long n)
    {
        ClosedFormText f = new PowerOfA(n);
        SparseSpec spec = SparseSpec.For(f);
        return WriteContainer("red", spec, SparseContainer.EndPatches(f, SparseContainer.Plan(spec), 4096));
    }

    private static EsaOpenOptions At(EsaValidationLevel level, EsaAccessMode access = EsaAccessMode.Checked) => new() { Validation = level, Access = access };

    // ---------------------------------------------------------------------------------- n = 2^31 + 10, POS_WIDTH 4

    [Test]
    [Property("Kind", "guard")]
    public void Width4_TheHeaderLevelOpens_InTimeAndMemory()
    {
        string path = Shell(TwoTo31 + 10);
        var clock = Stopwatch.StartNew();
        using SesaFile file = SesaFile.Open(path, At(EsaValidationLevel.Header));
        double seconds = clock.Elapsed.TotalSeconds;
        (double start, double peak) = RssAnonSoFar();
        TestContext.Out.WriteLine($"OPEN level=Header n={TwoTo31 + 10} seconds={seconds:F4} rss_anon_added_mb={peak - start:F1}");
        Assert.That(file.Header.TextLength, Is.EqualTo(TwoTo31 + 10));
        Assert.That(seconds, Is.LessThan(HeaderBudgetSeconds));
        Assert.That(peak - start, Is.LessThan(AnonBudgetMb));
    }

    [Test]
    [Property("Kind", "guard")]
    public void Width4_TheSectionsLevelOpens_AndReadsEveryByteOfTheSections()
    {
        string path = Shell(TwoTo31 + 10);
        var clock = Stopwatch.StartNew();
        using SesaFile file = SesaFile.Open(path, At(EsaValidationLevel.Sections));
        double seconds = clock.Elapsed.TotalSeconds;
        (double start, double peak) = RssAnonSoFar();
        TestContext.Out.WriteLine($"OPEN level=Sections n={TwoTo31 + 10} seconds={seconds:F3} section_bytes={file.Report.SectionBytesChecksummed} rss_anon_added_mb={peak - start:F1}");
        Assert.That(file.Report.SectionBytesChecksummed, Is.GreaterThan(3 * 4 * TwoTo31));
        Assert.That(peak - start, Is.LessThan(AnonBudgetMb));
    }

    [Test]
    [Property("TurnsGreenAt", "S2.3")]
    public void Width4_TheViewsOfAnOpenedContainerExist_NoNullAtAnySize()
    {
        string path = Shell(TwoTo31 + 10);
        using SesaFile file = SesaFile.Open(path, At(EsaValidationLevel.Header));
        Assert.That(file.Lcp, Is.Not.Null, "SesaFile.cs:62: _lcp = null! for a text of 2^31 - 1 symbols or more");
        Assert.That(file.Child, Is.Not.Null, "SesaFile.cs:63: _child = null! for a text of 2^31 - 1 symbols or more");
    }

    [Test]
    [Property("TurnsGreenAt", "S2.4")]
    public void Width4_TheStatisticsOfTheHeaderAreNotNarrowedToInt()
    {
        string path = Shell(TwoTo31 + 10);
        using SesaFile file = SesaFile.Open(path, At(EsaValidationLevel.Header));
        Assert.That((long)file.Statistics.Length, Is.EqualTo(TwoTo31 + 10), "SesaFile.Statistics: (int)h.TextLength wraps to a negative number");
        Assert.That((long)file.Statistics.LrsDepth, Is.EqualTo(TwoTo31 + 9), "LRS_DEPTH of a^n is n - 1");
    }

    [Test]
    [Property("TurnsGreenAt", "S2.4")]
    public void Width4_TheDescentOverTheContainerAnswersContains_AtTheEndsOfTheText()
    {
        // the ends of a^n are real in the shell (ranks 0 .. 4095 and the last 4096): the descent of a^k with k < 4000 reads only those cells
        string path = Shell(TwoTo31 + 10);
        using SesaFile file = SesaFile.Open(path, At(EsaValidationLevel.Header));
        EsaIndex index = EsaIndex.FromFile(file);
        Assert.That(index.Contains("aaaa"), Is.True, "a^4 occurs n - 3 times (S1: no search structure above 2^31 - 2 symbols)");
        Assert.That(index.Contains("b"), Is.False);
        Assert.That(index.Contains("ab"), Is.False);
    }

    [Test]
    [Property("TurnsGreenAt", "S2.8")]
    public void Width4_TheStructureLevelIsNotRefusedForTheLengthOfTheText_ItFindsTheDamageOfTheShell()
    {
        string path = Shell(TwoTo31 + 10);
        Assert.Throws<InvalidDataException>(() => SesaFile.Open(path, At(EsaValidationLevel.Structure)).Dispose(), "the Structure level must run above 2^31 symbols and say that the zero suffix array is damaged (S1: NotSupportedException, SesaVerifier.RequireSupportedLength)");
    }

    [Test]
    [Property("Kind", "guard")]
    public void Width4_TheIntSurfaceRefusesTheContainer_NamingTheLongApi()
    {
        string path = Shell(TwoTo31 + 10);
        NotSupportedException e = Assert.Throws<NotSupportedException>(() => SuffixIndexFactory.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Header }))!;
        Assert.That(e.Message, Does.Contain("long API"), "ISuffixTree stays 32-bit: Open refuses with the exception that names the long API (ADR-005 section 3.2)");
        Assert.Throws<NotSupportedException>(() => EsaSuffixTreeFactory.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Header }));
    }

    // ---------------------------------------------------------------------------------- n = 2^32 + 10, POS_WIDTH 5

    [Test]
    [Property("TurnsGreenAt", "S2.7")]
    public void Width5_TheHeaderLevelOpens_InUnderFiveSecondsAndUnder100MbOfAnonymousMemory()
    {
        string path = Shell(TwoTo32 + 10);
        var clock = Stopwatch.StartNew();
        using SesaFile file = SesaFile.Open(path, At(EsaValidationLevel.Header));
        double seconds = clock.Elapsed.TotalSeconds;
        (double start, double peak) = RssAnonSoFar();
        TestContext.Out.WriteLine($"OPEN level=Header n={TwoTo32 + 10} seconds={seconds:F4} rss_anon_added_mb={peak - start:F1}");
        Assert.That(file.Header.TextLength, Is.EqualTo(TwoTo32 + 10));
        Assert.That(file.Header.PosWidth, Is.EqualTo(5));
        Assert.That(seconds, Is.LessThan(HeaderBudgetSeconds));
        Assert.That(peak - start, Is.LessThan(AnonBudgetMb));
    }

    [Test]
    [Property("TurnsGreenAt", "S2.7")]
    public void Width5_TheSectionsLevelOpens_AndTheChecksumPassOver78GbOfHolesIsFast()
    {
        string path = Shell(TwoTo32 + 10);
        var clock = Stopwatch.StartNew();
        using SesaFile file = SesaFile.Open(path, At(EsaValidationLevel.Sections));
        double seconds = clock.Elapsed.TotalSeconds;
        (double start, double peak) = RssAnonSoFar();
        TestContext.Out.WriteLine($"OPEN level=Sections n={TwoTo32 + 10} seconds={seconds:F3} section_bytes={file.Report.SectionBytesChecksummed} rss_anon_added_mb={peak - start:F1}");
        Assert.That(file.Report.SectionBytesChecksummed, Is.GreaterThan(18 * TwoTo32));
        Assert.That(seconds, Is.LessThan(SectionsBudgetSeconds));
        Assert.That(peak - start, Is.LessThan(AnonBudgetMb));
    }

    [Test]
    [Property("TurnsGreenAt", "S2.7")]
    public void Width5_TheIntSurfaceRefusesTheContainerAsTooLong_NotAsDamaged()
    {
        string path = Shell(TwoTo32 + 10);
        NotSupportedException e = Assert.Throws<NotSupportedException>(() => SuffixIndexFactory.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Header }))!;
        Assert.That(e.Message, Does.Contain("long API"));
    }

    [Test]
    [Property("TurnsGreenAt", "S2.7")]
    public void Width5_TheDescentOverTheContainerAnswersContains_AtTheEndsOfTheText()
    {
        string path = Shell(TwoTo32 + 10);
        using SesaFile file = SesaFile.Open(path, At(EsaValidationLevel.Header));
        EsaIndex index = EsaIndex.FromFile(file);
        Assert.That(index.Contains("aaaa"), Is.True);
        Assert.That(index.Contains("b"), Is.False);
    }

    [Test]
    [Property("TurnsGreenAt", "S2.8")]
    public void Width5_TheStructureLevelFindsTheDamageOfTheShell_NotALimit()
    {
        string path = Shell(TwoTo32 + 10);
        InvalidDataException e = Assert.Throws<InvalidDataException>(() => SesaFile.Open(path, At(EsaValidationLevel.Structure)).Dispose())!;
        Assert.That(e.Message, Does.Not.Contain("POS_WIDTH"), "the damage found must be in the suffix array of the shell, not the refusal of the width (S1: POS_WIDTH is 5: 4 or 8)");
    }
}
