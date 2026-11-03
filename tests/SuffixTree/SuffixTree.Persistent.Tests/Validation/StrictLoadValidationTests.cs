using System.Buffers.Binary;
using System.Text;

namespace SuffixTree.Persistent.Tests.Validation;

/// <summary>
/// ST-S0 WP 0.5: <see cref="PersistentSuffixTree.Load(IStorageProvider)"/> validates only the header, so a body
/// that is corrupt behind a valid header is accepted and later queried through raw pointers. The opt-in strict
/// validation (<see cref="PersistentSuffixTree.ValidateStructure"/>, <c>Load(storage, validateStructure: true)</c>,
/// <c>LoadPersistent(path, validateStructure: true)</c>) must accept every tree the builder produces and reject
/// every structurally inconsistent file with <see cref="InvalidOperationException"/> — never an access violation,
/// another exception type, or a hang.
/// </summary>
[TestFixture]
[Category("Validation")]
public class StrictLoadValidationTests
{
    // Compact (v6) layout facts used to corrupt images by hand.
    private const int RootOffset = PersistentConstants.HEADER_SIZE_V6;
    private const int NodeStart = 0, NodeEnd = 4, NodeLink = 8, NodeLeafCount = 12, NodeHead = 16, NodeChildCount = 20;
    private const int EntrySize = 8;

    // ───────────────────────── image helpers ─────────────────────────

    private static byte[] BuildImage(string text, long? compactOffsetLimit = null)
    {
        using var storage = new HeapStorageProvider();
        var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact);
        if (compactOffsetLimit is { } limit)
            builder.CompactOffsetLimit = limit;
        builder.Build(new StringTextSource(text));
        var bytes = new byte[storage.Size];
        storage.ReadBytes(0, bytes, 0, bytes.Length);
        return bytes;
    }

    private static HeapStorageProvider Open(byte[] image)
    {
        var storage = new HeapStorageProvider(Math.Max(image.Length, 1));
        storage.Allocate(image.Length);
        storage.WriteBytes(0, image, 0, image.Length);
        return storage;
    }

    private static uint U32(byte[] b, long o) => BinaryPrimitives.ReadUInt32LittleEndian(b.AsSpan((int)o));
    private static int I32(byte[] b, long o) => BinaryPrimitives.ReadInt32LittleEndian(b.AsSpan((int)o));
    private static long I64(byte[] b, long o) => BinaryPrimitives.ReadInt64LittleEndian(b.AsSpan((int)o));
    private static void SetU32(byte[] b, long o, uint v) => BinaryPrimitives.WriteUInt32LittleEndian(b.AsSpan((int)o), v);
    private static void SetI32(byte[] b, long o, int v) => BinaryPrimitives.WriteInt32LittleEndian(b.AsSpan((int)o), v);

    private static bool IsLeaf(byte[] b, long node) => U32(b, node + NodeEnd) == uint.MaxValue;
    private static int ChildCount(byte[] b, long node) => I32(b, node + NodeChildCount) & int.MaxValue;
    private static long EntryOffset(byte[] b, long node, int i) => U32(b, node + NodeHead) + (long)i * EntrySize;
    private static long ChildAt(byte[] b, long node, int i) => U32(b, EntryOffset(b, node, i) + 4);

    /// <summary>All nodes below <paramref name="node"/> in depth-first order (compact, non-hybrid images only).</summary>
    private static List<long> Descendants(byte[] b, long node)
    {
        var result = new List<long>();
        var stack = new Stack<long>();
        stack.Push(node);
        while (stack.Count > 0)
        {
            long current = stack.Pop();
            for (int i = ChildCount(b, current) - 1; i >= 0; i--)
            {
                long child = ChildAt(b, current, i);
                result.Add(child);
                if (!IsLeaf(b, child))
                    stack.Push(child);
            }
        }
        return result;
    }

    private static long FirstLeaf(byte[] b) => Descendants(b, RootOffset).First(n => IsLeaf(b, n));
    private static long SecondLeaf(byte[] b) => Descendants(b, RootOffset).Where(n => IsLeaf(b, n)).Skip(1).First();
    private static long FirstInternal(byte[] b) => Descendants(b, RootOffset).First(n => !IsLeaf(b, n));
    private static long LastInternal(byte[] b) => Descendants(b, RootOffset).Last(n => !IsLeaf(b, n));

    // ───────────────────────── valid trees are accepted ─────────────────────────

    private static string RandomText(int length, int seed, string alphabet)
    {
        var rng = new Random(seed);
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(alphabet[rng.Next(alphabet.Length)]);
        return sb.ToString();
    }

    private static IEnumerable<string> ValidTexts()
    {
        yield return "";
        yield return "a";
        yield return "ab";
        yield return "aa";
        yield return "banana";
        yield return "mississippi";
        yield return "abracadabra";
        yield return new string('a', 200);
        yield return "abcdefghijklmnopqrstuvwxyz";
        yield return string.Concat(Enumerable.Repeat("ACGT", 50));
        yield return RandomText(3000, 1, "ACGT");
        yield return RandomText(2000, 2, "ACDEFGHIKLMNPQRSTVWY");
        yield return RandomText(500, 3, "ab");
        yield return "αβγαβγ日本語日本語";
        yield return "a\0b\0a\0b￿￿";
    }

    [TestCaseSource(nameof(ValidTexts))]
    public void ValidateStructure_AcceptsEveryBuiltTree_Heap(string text)
    {
        using var storage = Open(BuildImage(text));
        using var tree = PersistentSuffixTree.Load(storage, validateStructure: true);

        Assert.DoesNotThrow(() => tree.ValidateStructure());
        Assert.That(tree.Text.Length, Is.EqualTo(text.Length));
    }

    [TestCaseSource(nameof(ValidTexts))]
    public void LoadPersistent_StrictValidation_AcceptsEveryBuiltTree_MappedFile(string text)
    {
        string path = Path.Combine(Path.GetTempPath(), $"strict_{Guid.NewGuid():N}.tree");
        try
        {
            PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text), path).Dispose();

            using var tree = PersistentSuffixTreeFactory.LoadPersistent(path, validateStructure: true);

            Assert.That(tree.Text.Length, Is.EqualTo(text.Length));
            Assert.That(tree.Contains(text.Length > 3 ? text.Substring(1, 2) : string.Empty), Is.True);
        }
        finally
        {
            try { File.Delete(path); } catch (IOException) { }
        }
    }

    [TestCase("mississippi")]
    [TestCase("abracadabra")]
    [TestCase("aababcabcdabcde")]
    [TestCase("aaaaaaaaaa")]
    public void ValidateStructure_AcceptsHybridTreesAtEveryTransitionBoundary(string text)
    {
        long maxLimit = BuildImage(text).Length;
        int hybridTrees = 0;
        for (long limit = RootOffset + NodeLayout.Compact.NodeSize; limit < maxLimit; limit += NodeLayout.Compact.NodeSize)
        {
            using var storage = new HeapStorageProvider();
            var builder = new PersistentSuffixTreeBuilder(storage, NodeLayout.Compact) { CompactOffsetLimit = limit };
            builder.Build(new StringTextSource(text));
            if (!builder.IsHybrid)
                continue;
            hybridTrees++;

            using var tree = PersistentSuffixTree.Load(storage, validateStructure: false);
            Assert.DoesNotThrow(() => tree.ValidateStructure(), $"hybrid tree for '{text}' with compact limit {limit}");
        }
        Assert.That(hybridTrees, Is.GreaterThan(0), "the sweep must produce hybrid trees");
    }

    // ───────────────────────── default Load is unchanged: header only ─────────────────────────

    public static IEnumerable<TestCaseData> BodyCorruptions()
    {
        // Every case damages the node area only; the header stays valid, so the default Load still succeeds.
        static TestCaseData Case(string name, Action<byte[]> corrupt) => new TestCaseData(corrupt).SetName(name);

        yield return Case("RootChildCountHuge", b => SetI32(b, RootOffset + NodeChildCount, 0x7FFFFFFF));
        yield return Case("RootChildCountZero", b => SetI32(b, RootOffset + NodeChildCount, 0));
        yield return Case("RootMarkedAsLeaf", b => SetU32(b, RootOffset + NodeEnd, uint.MaxValue));
        yield return Case("RootJumpFlagInNonHybridTree", b => SetI32(b, RootOffset + NodeChildCount, ChildCount(b, RootOffset) | int.MinValue));
        yield return Case("ChildOffsetBeyondStorage", b => SetU32(b, EntryOffset(b, RootOffset, 1) + 4, (uint)b.Length + 4096));
        yield return Case("ChildOffsetInsideHeader", b => SetU32(b, EntryOffset(b, RootOffset, 1) + 4, 40));
        yield return Case("ChildOffsetMisaligned", b => SetU32(b, EntryOffset(b, RootOffset, 1) + 4, (uint)ChildAt(b, RootOffset, 1) + 1));
        yield return Case("ChildOffsetNull", b => SetU32(b, EntryOffset(b, RootOffset, 1) + 4, uint.MaxValue));
        yield return Case("ChildPointsToRoot_Cycle", b => SetU32(b, EntryOffset(b, FirstInternal(b), 0) + 4, RootOffset));
        yield return Case("ChildPointsToAncestor_Cycle", b =>
        {
            long inner = LastInternal(b);
            SetU32(b, EntryOffset(b, inner, 0) + 4, (uint)FirstInternal(b));
        });
        yield return Case("TwoEntriesShareOneChild", b => SetU32(b, EntryOffset(b, RootOffset, 2) + 4, (uint)ChildAt(b, RootOffset, 1)));
        yield return Case("ChildKeysUnsorted", b =>
        {
            long e1 = EntryOffset(b, RootOffset, 1), e2 = EntryOffset(b, RootOffset, 2);
            uint k1 = U32(b, e1), k2 = U32(b, e2);
            SetU32(b, e1, k2);
            SetU32(b, e2, k1);
        });
        yield return Case("ChildKeyDoesNotMatchEdge", b => SetU32(b, EntryOffset(b, RootOffset, 2), U32(b, EntryOffset(b, RootOffset, 2)) + 1));
        yield return Case("LeafCountOfLeafIsTwo", b => SetU32(b, FirstLeaf(b) + NodeLeafCount, 2));
        yield return Case("LeafCountOfInternalWrong", b => SetU32(b, FirstInternal(b) + NodeLeafCount, U32(b, FirstInternal(b) + NodeLeafCount) + 1));
        yield return Case("RootLeafCountWrong", b => SetU32(b, RootOffset + NodeLeafCount, U32(b, RootOffset + NodeLeafCount) + 1));
        yield return Case("HeaderNodeCountWrong", b => SetI32(b, PersistentConstants.HEADER_OFFSET_NODE_COUNT, I32(b, PersistentConstants.HEADER_OFFSET_NODE_COUNT) + 1));
        yield return Case("SuffixLinkBeyondStorage", b => SetU32(b, FirstInternal(b) + NodeLink, 0xFFFFFFF0u));
        yield return Case("SuffixLinkToLeaf", b => SetU32(b, FirstInternal(b) + NodeLink, (uint)FirstLeaf(b)));
        yield return Case("SuffixLinkMisaligned", b => SetU32(b, FirstInternal(b) + NodeLink, (uint)FirstInternal(b) + 2));
        yield return Case("InternalEdgeEmpty", b => SetU32(b, FirstInternal(b) + NodeStart, U32(b, FirstInternal(b) + NodeEnd)));
        yield return Case("InternalEdgeEndsBeyondText", b => SetU32(b, FirstInternal(b) + NodeEnd, 1000));
        yield return Case("InternalEdgeStartsAfterEnd", b => SetU32(b, FirstInternal(b) + NodeStart, U32(b, FirstInternal(b) + NodeEnd) + 5));
        yield return Case("LeafStartBeyondText", b => SetU32(b, FirstLeaf(b) + NodeStart, 1000));
        yield return Case("TwoLeavesForOneSuffix", b => SetU32(b, SecondLeaf(b) + NodeStart, U32(b, FirstLeaf(b) + NodeStart)));
        yield return Case("ChildArrayHeadBeyondStorage", b => SetU32(b, FirstInternal(b) + NodeHead, 0xFFFFFFFEu));
        yield return Case("ChildArrayHeadMisaligned", b => SetU32(b, FirstInternal(b) + NodeHead, U32(b, FirstInternal(b) + NodeHead) + 2));
        yield return Case("ChildArrayRunsPastStorage", b => SetI32(b, FirstInternal(b) + NodeChildCount, 100_000));
        yield return Case("InternalNodeWithoutChildren", b => SetI32(b, FirstInternal(b) + NodeChildCount, 0));
        yield return Case("LeafWithChildren", b => SetI32(b, FirstLeaf(b) + NodeChildCount, 3));
        yield return Case("DeepestNodeIsALeaf", b => BinaryPrimitives.WriteInt64LittleEndian(b.AsSpan(PersistentConstants.HEADER_OFFSET_DEEPEST_NODE), FirstLeaf(b)));
        yield return Case("DeepestNodeIsNotTheDeepest", b =>
        {
            long shallow = Descendants(b, RootOffset).First(n => !IsLeaf(b, n));
            BinaryPrimitives.WriteInt64LittleEndian(b.AsSpan(PersistentConstants.HEADER_OFFSET_DEEPEST_NODE), shallow);
        });
        yield return Case("LrsDepthWrong", b => SetI32(b, PersistentConstants.HEADER_OFFSET_LRS_DEPTH, I32(b, PersistentConstants.HEADER_OFFSET_LRS_DEPTH) + 1));
    }

    [TestCaseSource(nameof(BodyCorruptions))]
    public void BodyCorruption_IsAcceptedByHeaderOnlyLoad_ButRejectedByStrictValidation(Action<byte[]> corrupt)
    {
        byte[] image = BuildImage("mississippi");
        corrupt(image);

        // Header-only validation (unchanged default): the damaged body is not noticed.
        using (var storage = Open(image))
        using (var tree = PersistentSuffixTree.Load(storage))
            Assert.That(tree, Is.Not.Null);

        // Strict validation names the problem.
        using var strictStorage = Open(image);
        var ex = Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(strictStorage, validateStructure: true));
        Assert.That(ex!.Message, Does.StartWith("Invalid storage format:"));
    }

    [Test]
    public void BodyCorruption_JumpTableSlotBeyondStorage_IsRejected()
    {
        // A hybrid tree stores some child arrays and suffix links behind 8-byte jump-table slots.
        const string text = "mississippi";
        long maxLimit = BuildImage(text).Length;
        byte[]? image = null;
        for (long limit = RootOffset + 2 * NodeLayout.Compact.NodeSize; limit < maxLimit && image is null; limit += NodeLayout.Compact.NodeSize)
        {
            byte[] candidate = BuildImage(text, limit);
            long transition = I64(candidate, PersistentConstants.HEADER_OFFSET_TRANSITION);
            long jumpStart = I64(candidate, PersistentConstants.HEADER_OFFSET_JUMP_START);
            if (transition >= 0 && jumpStart >= 0)
                image = candidate;
        }
        Assert.That(image, Is.Not.Null, "need a hybrid image with a jump table");

        long start = I64(image!, PersistentConstants.HEADER_OFFSET_JUMP_START);
        long end = I64(image!, PersistentConstants.HEADER_OFFSET_JUMP_END);
        Assert.That(end - start, Is.GreaterThanOrEqualTo(8));

        using (var ok = Open(image!))
            Assert.DoesNotThrow(() => PersistentSuffixTree.Load(ok, validateStructure: true).Dispose(), "the untouched hybrid image is valid");

        // Damage every slot in turn: a slot is either a child-array pointer or a suffix-link target, and in both
        // cases pointing beyond the storage must be reported (not dereferenced).
        int rejected = 0;
        for (long slot = start; slot + 8 <= end; slot += 8)
        {
            byte[] damaged = (byte[])image!.Clone();
            BinaryPrimitives.WriteInt64LittleEndian(damaged.AsSpan((int)slot), damaged.LongLength + 10_000);
            using var storage = Open(damaged);
            try
            {
                PersistentSuffixTree.Load(storage, validateStructure: true).Dispose();
            }
            catch (InvalidOperationException)
            {
                rejected++;
            }
        }
        Assert.That(rejected, Is.GreaterThan(0), "at least the child-array slots must be rejected");
    }

    [Test]
    public void StrictLoad_RejectedStructure_DisposesTheStorage()
    {
        byte[] image = BuildImage("banana");
        SetI32(image, RootOffset + NodeChildCount, 0x7FFFFFFF);
        var storage = Open(image);

        Assert.Throws<InvalidOperationException>(() => PersistentSuffixTree.Load(storage, validateStructure: true));

        Assert.Throws<ObjectDisposedException>(() => _ = storage.ReadInt32(0), "the storage is released when the structure is rejected");
    }

    [Test]
    public void StrictLoad_WithFalse_IsTheHeaderOnlyLoad()
    {
        byte[] image = BuildImage("banana");
        SetI32(image, RootOffset + NodeChildCount, 0x7FFFFFFF);
        using var storage = Open(image);

        using var tree = PersistentSuffixTree.Load(storage, validateStructure: false);

        Assert.Throws<InvalidOperationException>(() => tree.ValidateStructure());
    }

    [Test]
    public void ValidateStructure_AfterDispose_ThrowsObjectDisposed()
    {
        using var storage = Open(BuildImage("banana"));
        var tree = PersistentSuffixTree.Load(storage);
        tree.Dispose();

        Assert.Throws<ObjectDisposedException>(() => tree.ValidateStructure());
    }

    [Test]
    public void LoadPersistent_StrictValidation_RejectsCorruptFile_AndReleasesIt()
    {
        string path = Path.Combine(Path.GetTempPath(), $"strict_{Guid.NewGuid():N}.tree");
        try
        {
            PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource("mississippi"), path).Dispose();
            byte[] image = File.ReadAllBytes(path);
            SetU32(image, EntryOffset(image, RootOffset, 1) + 4, (uint)image.Length + 4096);
            File.WriteAllBytes(path, image);

            using (PersistentSuffixTreeFactory.LoadPersistent(path)) { /* header-only: loads */ }
            Assert.Throws<InvalidOperationException>(() => PersistentSuffixTreeFactory.LoadPersistent(path, validateStructure: true));

            // The failed strict load must not keep the file open or mapped.
            Assert.DoesNotThrow(() => File.Delete(path));
        }
        finally
        {
            try { File.Delete(path); } catch (IOException) { }
        }
    }

    // ───────────────────────── fuzz: flipped bytes ─────────────────────────

    /// <summary>
    /// Runs reads that use raw-pointer traversal on a (possibly corrupt but validated) tree. Results may be wrong
    /// for a corrupt tree; the contract being tested is that nothing throws and nothing hangs or crashes.
    /// </summary>
    private static void QueryBattery(PersistentSuffixTree tree, string text)
    {
        for (int i = 0; i < text.Length; i++)
        {
            for (int len = 1; len <= Math.Min(6, text.Length - i); len++)
            {
                string p = text.Substring(i, len);
                _ = tree.Contains(p);
                _ = tree.CountOccurrences(p);
                _ = tree.FindAllOccurrences(p);
            }
        }
        _ = tree.Contains("zzz");
        _ = tree.CountOccurrences("ss");
        _ = tree.LongestRepeatedSubstring();
        _ = tree.NodeCount;
        _ = tree.LeafCount;
        _ = tree.MaxDepth;
        tree.Traverse(new CountingVisitor());
        _ = tree.LongestCommonSubstring(text);
        _ = tree.FindExactMatchAnchors(text, 2);
        _ = tree.PrintTree();
    }

    /// <summary>Delegating storage that records every byte the validator reads.</summary>
    private sealed class RecordingStorage(IStorageProvider inner, HashSet<long> touched) : IStorageProvider
    {
        private T Note<T>(long offset, int size, T value)
        {
            for (int i = 0; i < size; i++)
            {
                touched.Add(offset + i);
            }
            return value;
        }

        public long Size => inner.Size;
        public void EnsureCapacity(long capacity) => inner.EnsureCapacity(capacity);
        public int ReadInt32(long offset) => Note(offset, 4, inner.ReadInt32(offset));
        public void WriteInt32(long offset, int value) => inner.WriteInt32(offset, value);
        public uint ReadUInt32(long offset) => Note(offset, 4, inner.ReadUInt32(offset));
        public void WriteUInt32(long offset, uint value) => inner.WriteUInt32(offset, value);
        public long ReadInt64(long offset) => Note(offset, 8, inner.ReadInt64(offset));
        public void WriteInt64(long offset, long value) => inner.WriteInt64(offset, value);
        public char ReadChar(long offset) => Note(offset, 2, inner.ReadChar(offset));
        public void WriteChar(long offset, char value) => inner.WriteChar(offset, value);
        public void ReadBytes(long offset, byte[] buffer, int start, int count)
        {
            Note(offset, count, 0);
            inner.ReadBytes(offset, buffer, start, count);
        }
        public void WriteBytes(long offset, byte[] buffer, int start, int count) => inner.WriteBytes(offset, buffer, start, count);
        public long Allocate(int size) => inner.Allocate(size);
        public void Dispose() => inner.Dispose();
    }

    [Test]
    [CancelAfter(300_000)]
    public void FlippingAByteTheValidatorReads_IsMostlyDetected([Values(0xFF, 0x01)] int mask)
    {
        const string text = "mississippi";
        byte[] original = BuildImage(text);

        // Bytes the strict validation reads on the pristine image = the live part of the file
        // (the area between the child arrays and the text contains unused space).
        var live = new HashSet<long>();
        using (var recording = new RecordingStorage(Open(original), live))
        using (var tree = PersistentSuffixTree.Load(recording, validateStructure: true))
            Assert.That(tree, Is.Not.Null);
        long textOffset = I64(original, PersistentConstants.HEADER_OFFSET_TEXT_OFF);
        var liveBytes = live.Where(o => o >= RootOffset && o < textOffset).OrderBy(o => o).ToList();
        Assert.That(liveBytes.Count, Is.GreaterThan(300), "the node and child-array bytes must be read");

        int detected = 0;
        var undetected = new List<long>();
        foreach (long offset in liveBytes)
        {
            byte[] image = (byte[])original.Clone();
            image[offset] ^= (byte)mask;
            using var storage = Open(image);
            try
            {
                PersistentSuffixTree.Load(storage, validateStructure: true).Dispose();
                undetected.Add(offset);
            }
            catch (InvalidOperationException)
            {
                detected++;
            }
        }

        TestContext.Out.WriteLine($"mask {mask:X2}: {liveBytes.Count} live bytes, {detected} detected, {undetected.Count} undetected: {string.Join(",", undetected.Take(60))}");
        // The root's Start/End fields (the first 8 bytes of the root node) are the only inspected bytes whose value
        // is irrelevant: the root has no incoming edge. Every other inspected byte must be protected.
        Assert.That(undetected.Where(o => o < RootOffset || o >= RootOffset + 8), Is.Empty,
            "flipping a byte that the validator inspects must be detected (except the root's unused Start/End)");
        Assert.That(detected, Is.GreaterThan(liveBytes.Count * 3 / 4));
    }

    private sealed class CountingVisitor : ISuffixTreeVisitor
    {
        public int Nodes;
        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) => Nodes++;
        public void EnterBranch(int key) { }
        public void ExitBranch() { }
    }

    [Test]
    [CancelAfter(300_000)]
    public void FlippingAnyByte_NeverCrashesOrHangs_AndStrictValidationOnlyThrowsInvalidOperation([Values(0xFF, 0x01, 0x80)] int mask)
    {
        const string text = "mississippi";
        byte[] original = BuildImage(text);
        int rejected = 0, accepted = 0;

        for (int i = 0; i < original.Length; i++)
        {
            byte[] image = (byte[])original.Clone();
            image[i] ^= (byte)mask;

            using var storage = Open(image);
            PersistentSuffixTree? tree = null;
            try
            {
                tree = PersistentSuffixTree.Load(storage, validateStructure: true);
            }
            catch (InvalidOperationException)
            {
                rejected++;
                continue;
            }

            accepted++;
            using (tree)
            {
                // Validation passed: the structure is memory-safe, so queries must complete without exceptions
                // (answers may differ from the original tree for a flip that kept the structure consistent).
                Assert.DoesNotThrow(() => QueryBattery(tree, text), $"queries on the tree with byte {i} flipped (mask {mask:X2})");
            }
        }

        TestContext.Out.WriteLine($"mask {mask:X2}: {rejected} of {original.Length} flips rejected, {accepted} accepted");
        // Most of the image is unused space between the child arrays and the text, so the rate over all bytes is low;
        // FlippingAByteTheValidatorReads_IsMostlyDetected checks the detection rate over the bytes that matter.
        Assert.That(rejected, Is.GreaterThan(300), "byte flips in the node area must be detected");
    }

    [Test]
    [CancelAfter(300_000)]
    public void RandomCorruptionOfMappedFile_NeverCrashesAfterStrictValidation()
    {
        // Same property on a memory-mapped file, where queries use the unchecked raw-pointer fast path.
        foreach (var (text, seed) in new[] { ("mississippi", 1), (RandomText(300, 5, "ACGT"), 2), ("abracadabraabracadabra", 3) })
        {
            string path = Path.Combine(Path.GetTempPath(), $"strict_{Guid.NewGuid():N}.tree");
            try
            {
                PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text), path).Dispose();
                byte[] original = File.ReadAllBytes(path);
                var rng = new Random(seed);

                for (int round = 0; round < 150; round++)
                {
                    byte[] image = (byte[])original.Clone();
                    int flips = 1 + rng.Next(3);
                    for (int f = 0; f < flips; f++)
                        image[rng.Next(image.Length)] ^= (byte)(1 << rng.Next(8));
                    File.WriteAllBytes(path, image);

                    PersistentSuffixTree tree;
                    try
                    {
                        tree = PersistentSuffixTreeFactory.LoadPersistent(path, validateStructure: true);
                    }
                    catch (InvalidOperationException)
                    {
                        continue;
                    }

                    using (tree)
                        Assert.DoesNotThrow(() => QueryBattery(tree, text.Length > 30 ? text[..30] : text), $"round {round}, seed {seed}");
                }
            }
            finally
            {
                try { File.Delete(path); } catch (IOException) { }
            }
        }
    }
}
