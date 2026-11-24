using System.Diagnostics;
using System.Globalization;

namespace SuffixTree.Esa.Tests.Large;

/// <summary>
/// The base of the tests of the <c>Large</c> category (work package S2.1, deliverable D). A derived fixture carries <c>[Category("Large")]</c> and <c>[Explicit(...)]</c> itself (never in the default run); the base
/// <list type="bullet">
/// <item>checks the disk once per fixture: the directory (<c>$SEQERON_LARGE_DIR</c> or a subdirectory of the system temporary directory) must make sparse files (<see cref="LargeDisk.RequireSparse"/>), else every test of the fixture is
/// inconclusive with the message, and nothing is written;</item>
/// <item>tracks every file that a test asks for with <see cref="NewPath"/> and deletes them all in the tear-down, whatever the outcome;</item>
/// <item>prints one time/space report line per test (<c>LARGE-REPORT</c>: seconds, the peak of RssAnon sampled every 20 ms, the disk that the tracked files occupy and their nominal size).</item>
/// </list>
/// </summary>
public abstract class LargeFixtureBase
{
    private readonly List<string> _files = [];
    private Stopwatch _clock = new();
    private RssSampler? _sampler;
    private long _startKb;

    /// <summary>The directory of the sparse files.</summary>
    protected static string Directory => LargeDisk.Directory;

    /// <summary>Makes the whole fixture inconclusive when the filesystem does not make sparse files.</summary>
    [OneTimeSetUp]
    public void RequireASparseFilesystem()
    {
        try
        {
            LargeDisk.RequireSparse(Directory);
        }
        catch (SparseUnavailableException e)
        {
            Assert.Inconclusive(e.Message);
        }
    }

    [SetUp]
    public void StartTheClock()
    {
        _files.Clear();
        _startKb = ReadRssAnonKb();
        _sampler = new RssSampler();
        _clock = Stopwatch.StartNew();
    }

    [TearDown]
    public void ReportAndCleanUp()
    {
        double seconds = _clock.Elapsed.TotalSeconds;
        double peak = _sampler?.Stop() ?? 0;
        long allocated = 0, nominal = 0;
        foreach (string path in _files)
        {
            if (!File.Exists(path)) continue;
            nominal += new FileInfo(path).Length;
            allocated += Math.Max(0, LargeDisk.AllocatedBytes(path));
        }
        foreach (string path in _files)
        {
            try
            {
                File.Delete(path);
            }
            catch (IOException)
            {
                TestContext.Out.WriteLine($"LARGE-CLEANUP could not delete {path}");
            }
        }
        string line = string.Create(CultureInfo.InvariantCulture, $"LARGE-REPORT test={TestContext.CurrentContext.Test.Name} outcome={TestContext.CurrentContext.Result.Outcome.Status} seconds={seconds:F2} rss_anon_peak_mb={peak:F1} disk_mb={allocated / 1048576.0:F2} nominal_gb={nominal / 1073741824.0:F2}");
        TestContext.Out.WriteLine(line);
        TestContext.Progress.WriteLine(line);
    }

    /// <summary>RssAnon at the start of the test and its peak so far (MB), sampled every 20 ms: what the test added to the anonymous memory of the process.</summary>
    internal (double StartMb, double PeakMb) RssAnonSoFar() => (_startKb / 1024.0, Math.Max(_sampler?.PeakKb ?? 0, ReadRssAnonKb()) / 1024.0);

    /// <summary>A path for a sparse file; the file is deleted in the tear-down.</summary>
    protected string NewPath(string stem)
    {
        string path = LargeDisk.NewPath(Directory, stem);
        _files.Add(path);
        return path;
    }

    /// <summary>Writes a sparse container and returns its path (deleted in the tear-down).</summary>
    internal string WriteContainer(string stem, SparseSpec spec, IReadOnlyList<SparsePatch> patches)
    {
        string path = NewPath(stem);
        SparseContainer.Write(path, spec, patches);
        return path;
    }

    /// <summary>The attributes of a <c>Large</c> fixture, in words (the attributes themselves are on each class).</summary>
    protected const string ExplicitReason = "Heavy: sparse files of 8 - 80 GB nominal (a few MB on disk); run with --filter Category=Large (see tests/SuffixTree/SuffixTree.Esa.Tests/Large/README.md)";

    /// <summary>Samples <c>RssAnon</c> of <c>/proc/self/status</c> every 20 ms and keeps the largest value.</summary>
    private sealed class RssSampler
    {
        private readonly Thread _thread;
        private volatile bool _stop;
        private long _peakKb;

        internal RssSampler()
        {
            _peakKb = ReadRssAnonKb();
            _thread = new Thread(Run) { IsBackground = true, Name = "rss-sampler" };
            _thread.Start();
        }

        internal long PeakKb => Interlocked.Read(ref _peakKb);

        internal double Stop()
        {
            _stop = true;
            _thread.Join();
            return Math.Max(Interlocked.Read(ref _peakKb), ReadRssAnonKb()) / 1024.0;
        }

        private void Run()
        {
            while (!_stop)
            {
                long kb = ReadRssAnonKb();
                if (kb > Interlocked.Read(ref _peakKb)) Interlocked.Exchange(ref _peakKb, kb);
                Thread.Sleep(20);
            }
        }
    }

    /// <summary><c>RssAnon</c> of this process in KiB (0 where <c>/proc</c> has none).</summary>
    protected static long ReadRssAnonKb()
    {
        try
        {
            foreach (string line in File.ReadLines("/proc/self/status"))
            {
                if (!line.StartsWith("RssAnon:", StringComparison.Ordinal)) continue;
                string[] parts = line.Split(' ', StringSplitOptions.RemoveEmptyEntries);
                return long.Parse(parts[1], CultureInfo.InvariantCulture);
            }
        }
        catch (Exception e) when (e is IOException or FormatException or UnauthorizedAccessException)
        {
            // no /proc: the report says 0
        }
        return 0;
    }
}
