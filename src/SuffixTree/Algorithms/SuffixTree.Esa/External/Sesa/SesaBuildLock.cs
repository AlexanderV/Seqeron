namespace SuffixTree.Esa.External;

/// <summary>
/// The lock that keeps two resumable builds of one target from running at once (work package 1b.9e): the file <c>&lt;target&gt;.lock</c> opened with <see cref="FileShare.None"/> (managed code: the runtime takes an advisory
/// <c>flock</c> on Unix and a sharing mode on Windows, and the operating system gives it up when the process dies, so a killed build never leaves the target locked). The holder removes the file when it is done; a file left
/// by a killed process is only a name, and the next build locks it again. Two builds that open the file in the instant between a finishing holder's last use and its removal can both believe that they hold it (one holds
/// a removed file): a window of microseconds that needs a third process to matter, documented in ADR-004 section 5a "1b.9e" and not closed. Builds with <see cref="SesaBuildOptions.Resume"/> off do not take the lock.
/// </summary>
internal sealed class SesaBuildLock : IDisposable
{
    private FileStream? _stream;

    private SesaBuildLock(FileStream stream, string path)
    {
        _stream = stream;
        Path = path;
    }

    /// <summary>The path of the lock file.</summary>
    internal string Path { get; }

    /// <summary>Takes the lock.</summary>
    /// <exception cref="IOException">Another build of the target holds it (or the file cannot be created).</exception>
    internal static SesaBuildLock Acquire(string path)
    {
        ArgumentNullException.ThrowIfNull(path);
        try
        {
            var stream = new FileStream(path, FileMode.OpenOrCreate, FileAccess.ReadWrite, FileShare.None, 1);
            return new SesaBuildLock(stream, path);
        }
        catch (IOException e) when (File.Exists(path))
        {
            throw new IOException($"another build of this target holds the lock {path}: two resumable builds of one container cannot run at once", e);
        }
    }

    public void Dispose()
    {
        FileStream? stream = Interlocked.Exchange(ref _stream, null);
        if (stream is null) return;
        try
        {
            File.Delete(Path); // while it is still held: a file that is deleted after the lock is given up could be one that the next holder has just created
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // a platform that does not delete an open file: the name stays, the next build locks it again
        }
        stream.Dispose();
    }
}
