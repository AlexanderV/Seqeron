namespace SuffixTree.Persistent;

public partial class PersistentSuffixTreeBuilder
{
    private sealed class HybridTransitionState
    {
        public long TransitionOffset = -1;
        public long JumpTableStart = -1;
        public long JumpTableEnd = -1;
        public long CompactNodesEnd = -1;
        public long JumpTableReserveStart = -1;
        public long JumpTableReserveEnd = -1;
    }

    private sealed class UkkonenRuntimeState
    {
        public long ActiveNodeOffset;
        public int ActiveEdgeIndex;
        public int ActiveLength;
        public int Remainder;
        public int Position = -1;
        public long LastCreatedInternalNodeOffset = PersistentConstants.NULL_OFFSET;
    }
}
