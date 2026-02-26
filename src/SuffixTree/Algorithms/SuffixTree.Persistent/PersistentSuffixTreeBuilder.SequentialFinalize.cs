namespace SuffixTree.Persistent;

#pragma warning disable CA1001 // _childStore and _depthStore ownership finalized in FinalizeTree/SequentialFinalize
public partial class PersistentSuffixTreeBuilder
#pragma warning restore CA1001
{
    // ──────────────── Sequential multi-pass finalize ────────────────

    /// <summary>
    /// Replaces the old DFS-based finalization with sequential disk passes to
    /// eliminate random I/O that caused 13+ minute DFS phases on genome-scale builds.
    ///
    /// <b>Pass 1</b> — Sequential scan of main storage (19 GB):
    ///   For each node in allocation order: read child linked list from child store,
    ///   sort children, write sorted child array (batch-allocated), update
    ///   ChildrenHead/ChildCount in node header. Build parent-graph in temp file.
    ///
    /// <b>Pass 2</b> — Bottom-up leaf count propagation (temp files only, ~7 GB):
    ///   BFS from leaves → parents using temp file + sequential queue file.
    ///   Zero main storage I/O.
    ///
    /// <b>Pass 3</b> — Sequential write-back (19 GB):
    ///   Write LeafCount from temp file. Track deepest internal node via depth store.
    ///
    /// Disk: ~7 GB temp files (auto-deleted). RAM: near-zero (just OS page cache).
    /// </summary>
    private unsafe void SequentialFinalize(long rootOffset)
    {
        var mmfMain = _mmfStorage;
        byte* childBase = _mmfChildStore != null ? _mmfChildStore.RawPointer : null;
        int nodeCount = _nodeCount;
        int hdr = HeaderSize;

        string tempNodePath = Path.Combine(Path.GetTempPath(), $"sfx_{Guid.NewGuid():N}_node.tmp");
        string tempQueuePath = Path.Combine(Path.GetTempPath(), $"sfx_{Guid.NewGuid():N}_queue.tmp");

        try
        {
            // ── Create temp file for node graph metadata ──
            long tempNodeCap = (long)nodeCount * TEMP_ENTRY;
            using var tempNodeMmf = new MappedFileStorageProvider(tempNodePath, Math.Max(tempNodeCap, 65536));
            tempNodeMmf.SetSize(tempNodeCap);
            byte* tp = tempNodeMmf.RawPointer;

            // Initialize parent = -1 (leafCount=0 and remaining=0 are zero-filled by OS)
            for (long i = 0; i < nodeCount; i++)
                *(int*)(tp + i * TEMP_ENTRY + TEMP_OFF_PARENT) = -1;

            // ═══════════ PASS 1: Sequential scan → child arrays + parent graph ═══════════
            _progress?.Report(("Pass 1: child arrays", 92.0));

            // Prefetch child store into OS page cache before linked list walks.
            // This converts random page faults into sequential read-ahead.
            if (_mmfChildStore != null)
                _mmfChildStore.PrefetchForBuild();

            // Batch allocation for child arrays
            long estimatedChildArrayBytes = (long)nodeCount * 10;
            int batchSize = (int)Math.Min(estimatedChildArrayBytes, int.MaxValue - 65536);
            batchSize = Math.Max(batchSize, 4096);
            long batchStart = _storage.Allocate(batchSize);
            long batchEnd = batchStart + batchSize;
            long batchPos = batchStart;

            // Reusable child buffer (max 256 children per node)
            const int CHILD_BUF_MAX = 256;
            uint* cKeys = stackalloc uint[CHILD_BUF_MAX];
            long* cOffs = stackalloc long[CHILD_BUF_MAX];

            long compactEnd = IsHybrid
                ? (_compactNodesEnd >= 0 ? _compactNodesEnd : _transitionOffset)
                : _maxNodeEndOffset;
            int compactNodeSize = _initialLayout.NodeSize;
            int totalNodes = nodeCount;
            int reportInterval = Math.Max(totalNodes / 200, 1);
            int reportCountdown = reportInterval;

            int idx = 0;

            // Compact zone — sequential scan
            for (long off = hdr; off < compactEnd; off += compactNodeSize, idx++)
            {
                uint nodeStart, nodeEnd;
                if (mmfMain != null)
                {
                    nodeStart = mmfMain.ReadUInt32Unchecked(off);
                    nodeEnd = mmfMain.ReadUInt32Unchecked(off + 4);
                }
                else
                {
                    nodeStart = _storage.ReadUInt32(off);
                    nodeEnd = _storage.ReadUInt32(off + 4);
                }

                if (nodeEnd == PersistentConstants.BOUNDLESS)
                {
                    // Leaf: leafCount = 1, remaining = 0 (already zero)
                    *(uint*)(tp + (long)idx * TEMP_ENTRY + TEMP_OFF_LEAFCOUNT) = 1;
                }
                else
                {
                    // Internal node: process children
                    Pass1Internal(off, idx, _initialLayout, /* useMmfUnchecked */ mmfMain != null,
                        mmfMain, childBase, tp, cKeys, cOffs, CHILD_BUF_MAX,
                        ref batchPos, ref batchStart, ref batchEnd, batchSize);
                }

                if (--reportCountdown <= 0)
                {
                    reportCountdown = reportInterval;
                    _progress?.Report(("Pass 1: child arrays", 92.0 + (double)idx / totalNodes * 2.0));
                }
            }

            // Large zone — sequential scan (hybrid only)
            if (IsHybrid)
            {
                int largeNodeSize = NodeLayout.Large.NodeSize;
                for (long off = _transitionOffset; off < _maxNodeEndOffset; off += largeNodeSize, idx++)
                {
                    uint nodeEnd = _storage.ReadUInt32(off + 4);

                    if (nodeEnd == PersistentConstants.BOUNDLESS)
                    {
                        *(uint*)(tp + (long)idx * TEMP_ENTRY + TEMP_OFF_LEAFCOUNT) = 1;
                    }
                    else
                    {
                        Pass1Internal(off, idx, NodeLayout.Large, /* useMmfUnchecked */ false,
                            mmfMain, childBase, tp, cKeys, cOffs, CHILD_BUF_MAX,
                            ref batchPos, ref batchStart, ref batchEnd, batchSize);
                    }

                    if (--reportCountdown <= 0)
                    {
                        reportCountdown = reportInterval;
                        _progress?.Report(("Pass 1: child arrays", 92.0 + (double)idx / totalNodes * 2.0));
                    }
                }
            }

            // Dispose child store — free its page cache for subsequent passes
            childBase = null;
            if (_ownsChildStore) _childStore.Dispose();

            RunPass2LeafCountPropagation(tp, nodeCount, tempQueuePath);
            RunPass3WriteBackAndTrackLrs(tp, nodeCount, hdr, rootOffset, compactEnd, compactNodeSize, reportInterval, mmfMain);

            // Dispose depth store after Pass 3 (needed for deepest node tracking)
            if (_ownsDepthStore) _depthStore?.Dispose();
        }
        finally
        {
            TryDeleteTempFile(tempNodePath);
            TryDeleteTempFile(tempQueuePath);
        }
    }

    private unsafe void RunPass2LeafCountPropagation(byte* tp, int nodeCount, string tempQueuePath)
    {
        _progress?.Report(("Pass 2: leaf counts", 94.0));

        long queueCap = (long)nodeCount * 4;
        using var tempQueueMmf = new MappedFileStorageProvider(tempQueuePath, Math.Max(queueCap, 65536));
        tempQueueMmf.SetSize(queueCap);
        byte* qp = tempQueueMmf.RawPointer;

        // Enqueue all leaves (remaining == 0 AND leafCount == 1)
        long qWrite = 0;
        for (int i = 0; i < nodeCount; i++)
        {
            byte* e = tp + (long)i * TEMP_ENTRY;
            if (*(int*)(e + TEMP_OFF_REMAINING) == 0 && *(uint*)(e + TEMP_OFF_LEAFCOUNT) > 0)
            {
                *(int*)(qp + qWrite) = i;
                qWrite += 4;
            }
        }

        // BFS: propagate leaf counts bottom-up
        long qRead = 0;
        while (qRead < qWrite)
        {
            int ni = *(int*)(qp + qRead);
            qRead += 4;

            byte* e = tp + (long)ni * TEMP_ENTRY;
            int parentIdx = *(int*)(e + TEMP_OFF_PARENT);
            if (parentIdx < 0) continue; // root — no parent

            uint myLeaves = *(uint*)(e + TEMP_OFF_LEAFCOUNT);
            byte* pe = tp + (long)parentIdx * TEMP_ENTRY;
            *(uint*)(pe + TEMP_OFF_LEAFCOUNT) += myLeaves;
            int rem = *(int*)(pe + TEMP_OFF_REMAINING) - 1;
            *(int*)(pe + TEMP_OFF_REMAINING) = rem;

            if (rem == 0)
            {
                *(int*)(qp + qWrite) = parentIdx;
                qWrite += 4;
            }
        }
    }

    private unsafe void RunPass3WriteBackAndTrackLrs(
        byte* tp,
        int nodeCount,
        int headerSize,
        long rootOffset,
        long compactEnd,
        int compactNodeSize,
        int reportInterval,
        MappedFileStorageProvider? mmfMain)
    {
        _progress?.Report(("Pass 3: write-back", 96.0));

        long deepestOffset = rootOffset;
        int maxLrsDepth = 0;
        int idx = 0;
        int reportCountdown = reportInterval;

        // Compact zone
        for (long off = headerSize; off < compactEnd; off += compactNodeSize, idx++)
        {
            uint lc = *(uint*)(tp + (long)idx * TEMP_ENTRY + TEMP_OFF_LEAFCOUNT);
            if (mmfMain != null)
                mmfMain.WriteUInt32Unchecked(off + (uint)_initialLayout.OffsetLeafCount, lc);
            else
                _storage.WriteUInt32(off + _initialLayout.OffsetLeafCount, lc);

            // Track deepest internal node (read Start/End sequentially)
            uint nodeEnd = mmfMain != null
                ? mmfMain.ReadUInt32Unchecked(off + 4)
                : _storage.ReadUInt32(off + 4);
            if (nodeEnd != PersistentConstants.BOUNDLESS)
            {
                uint nodeStart = mmfMain != null
                    ? mmfMain.ReadUInt32Unchecked(off)
                    : _storage.ReadUInt32(off);
                int edgeLen = (int)(nodeEnd - nodeStart);
                uint dfr = GetBuildDepth(off);
                int sd = (int)dfr + edgeLen;
                if (sd > maxLrsDepth) { maxLrsDepth = sd; deepestOffset = off; }
            }

            if (--reportCountdown <= 0)
            {
                reportCountdown = reportInterval;
                _progress?.Report(("Pass 3: write-back", 96.0 + (double)idx / nodeCount * 1.5));
            }
        }

        // Large zone
        if (IsHybrid)
        {
            int largeNodeSize = NodeLayout.Large.NodeSize;
            for (long off = _transitionOffset; off < _maxNodeEndOffset; off += largeNodeSize, idx++)
            {
                uint lc = *(uint*)(tp + (long)idx * TEMP_ENTRY + TEMP_OFF_LEAFCOUNT);
                _storage.WriteUInt32(off + NodeLayout.Large.OffsetLeafCount, lc);

                uint nodeEnd = _storage.ReadUInt32(off + 4);
                if (nodeEnd != PersistentConstants.BOUNDLESS)
                {
                    uint nodeStart = _storage.ReadUInt32(off);
                    int edgeLen = (int)(nodeEnd - nodeStart);
                    uint dfr = GetBuildDepth(off);
                    int sd = (int)dfr + edgeLen;
                    if (sd > maxLrsDepth) { maxLrsDepth = sd; deepestOffset = off; }
                }

                if (--reportCountdown <= 0)
                {
                    reportCountdown = reportInterval;
                    _progress?.Report(("Pass 3: write-back", 96.0 + (double)idx / nodeCount * 1.5));
                }
            }
        }

        _deepestInternalNodeOffset = deepestOffset;
        _lrsDepth = maxLrsDepth;
    }

    private static void TryDeleteTempFile(string path)
    {
        try { File.Delete(path); } catch (IOException) { }
    }

}
