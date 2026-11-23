namespace SuffixTree.Esa.External;

/// <summary>
/// The checked build's list of the files that this process has mapped for reading (the packed text and the rank file of PH1 and PH3): a bounds check cannot see the one way in which a read through a mapping that is right in every
/// index goes wrong, a file that is made shorter under it. Creating a file again (<c>FileMode.Create</c>, the reserved <see cref="MappedOutputFile"/>) truncates it, and the next read of the mapping is a bus error that ends the process
/// with an <see cref="AccessViolationException"/>: two levels of the sample recursion that were given the same file name did exactly that (work package 1b.9h, the mutant <c>L10</c> of 1b.9c). In the checked build a creation of a
/// file that is mapped is an <see cref="EsaMemoryViolationException"/> that names the file; the production build has no such list (the file does not exist in it).
/// </summary>
internal static class MappedFileRegistry
{
#if ESA_BOUNDS
    private static readonly Dictionary<string, int> s_mapped = new(StringComparer.Ordinal);

    /// <summary>The file is mapped now (one more mapping of it); returns the name to give <see cref="Unmap"/>.</summary>
    internal static string Map(string path)
    {
        string full = Path.GetFullPath(path);
        lock (s_mapped) s_mapped[full] = s_mapped.GetValueOrDefault(full) + 1;
        return full;
    }

    /// <summary>One mapping of the file ended.</summary>
    internal static void Unmap(string full)
    {
        lock (s_mapped)
        {
            if (!s_mapped.TryGetValue(full, out int count)) return;
            if (count <= 1) s_mapped.Remove(full);
            else s_mapped[full] = count - 1;
        }
    }

    /// <summary>The file is about to be created (made empty): no mapping of it may be alive.</summary>
    internal static void RequireNotMapped(string path)
    {
        string full = Path.GetFullPath(path);
        bool mapped;
        lock (s_mapped) mapped = s_mapped.ContainsKey(full);
        if (mapped)
        {
            throw new EsaMemoryViolationException(MemoryViolationKind.OutOfBounds, $"the file {full} is mapped for reading by this process and is about to be created again: that truncates it under the mapping, and the next read of the mapping is a bus error");
        }
    }

    /// <summary>How many files are mapped now (the tests look at it: nothing may be left once the readers are disposed).</summary>
    internal static int Count
    {
        get
        {
            lock (s_mapped) return s_mapped.Count;
        }
    }
#else
    /// <summary>The production build keeps no list.</summary>
    internal static int Count => 0;
#endif
}
