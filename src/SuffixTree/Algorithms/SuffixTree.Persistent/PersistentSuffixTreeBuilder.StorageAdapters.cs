using System.Runtime.CompilerServices;

namespace SuffixTree.Persistent;

public partial class PersistentSuffixTreeBuilder
{
    private sealed class ChildStoreAdapter
    {
        private readonly IStorageProvider _store;
        private readonly MappedFileStorageProvider? _mmf;

        public ChildStoreAdapter(IStorageProvider store, MappedFileStorageProvider? mmf)
        {
            _store = store;
            _mmf = mmf;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static long EntryOffset(int index) => (long)index * CHILD_ENTRY_SIZE;

        public void ReadEntry(int index, out uint key, out int nextIndex, out long childOffset)
        {
            long entryOff = EntryOffset(index);
            if (_mmf != null)
            {
                key = _mmf.ReadUInt32Unchecked(entryOff + CE_OFF_KEY);
                nextIndex = _mmf.ReadInt32Unchecked(entryOff + CE_OFF_NEXT);
                childOffset = _mmf.ReadInt64Unchecked(entryOff + CE_OFF_CHILD);
                return;
            }

            key = _store.ReadUInt32(entryOff + CE_OFF_KEY);
            nextIndex = _store.ReadInt32(entryOff + CE_OFF_NEXT);
            childOffset = _store.ReadInt64(entryOff + CE_OFF_CHILD);
        }

        public void WriteChildOffset(int index, long childOffset)
        {
            long entryOff = EntryOffset(index);
            if (_mmf != null)
            {
                _mmf.WriteInt64Unchecked(entryOff + CE_OFF_CHILD, childOffset);
                return;
            }

            _store.WriteInt64(entryOff + CE_OFF_CHILD, childOffset);
        }

        public int AddEntry(int nextIndex, uint key, long childOffset)
        {
            long newOff = _store.Allocate(CHILD_ENTRY_SIZE);
            int newIndex = (int)(newOff / CHILD_ENTRY_SIZE);

            if (_mmf != null)
            {
                _mmf.WriteUInt32Unchecked(newOff + CE_OFF_KEY, key);
                _mmf.WriteInt32Unchecked(newOff + CE_OFF_NEXT, nextIndex);
                _mmf.WriteInt64Unchecked(newOff + CE_OFF_CHILD, childOffset);
                return newIndex;
            }

            _store.WriteUInt32(newOff + CE_OFF_KEY, key);
            _store.WriteInt32(newOff + CE_OFF_NEXT, nextIndex);
            _store.WriteInt64(newOff + CE_OFF_CHILD, childOffset);
            return newIndex;
        }
    }

    private sealed class DepthStoreAdapter
    {
        private readonly IStorageProvider _store;
        private readonly MappedFileStorageProvider? _mmf;

        public DepthStoreAdapter(IStorageProvider store, MappedFileStorageProvider? mmf)
        {
            _store = store;
            _mmf = mmf;
        }

        public long Size => _store.Size;

        public void EnsureSizeAtLeast(long size)
        {
            while (_store.Size < size)
                _store.Allocate(4);
        }

        public uint ReadUInt32(long offset)
        {
            if (_mmf != null)
                return _mmf.ReadUInt32Unchecked(offset);
            return _store.ReadUInt32(offset);
        }

        public void WriteUInt32(long offset, uint value)
        {
            if (_mmf != null)
            {
                _mmf.WriteUInt32Unchecked(offset, value);
                return;
            }

            _store.WriteUInt32(offset, value);
        }
    }
}
