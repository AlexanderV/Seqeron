namespace SuffixTree.Persistent;

public sealed partial class PersistentSuffixTree
{
    /// <summary>Disposes the tree, releasing storage and text source resources.</summary>
    public void Dispose()
    {
        if (Interlocked.CompareExchange(ref _disposed, 1, 0) != 0) return;
        GC.SuppressFinalize(this);

        // Dispose textSource BEFORE storage: MemoryMappedTextSource must
        // ReleasePointer() while the underlying accessor is still alive.
        // Use try/finally to guarantee storage is disposed even if textSource throws.
        try
        {
            if (_ownsTextSource && _textSource is IDisposable disposableText)
                disposableText.Dispose();
        }
        finally
        {
            _storage.Dispose();
        }
    }

    private void ThrowIfDisposed()
    {
        ObjectDisposedException.ThrowIf(Volatile.Read(ref _disposed) != 0, this);
    }
}
