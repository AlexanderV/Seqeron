namespace SuffixTree.Persistent;

public static class PersistentConstants
{
    public const int NODE_SIZE = 40;
    public const int CHILD_ENTRY_SIZE = 12;
    public const int HEADER_SIZE = 64;
    public const uint BOUNDLESS = uint.MaxValue;
    public const uint TERMINATOR_KEY = uint.MaxValue;  // Same bit pattern as -1 when cast to int
    public const long NULL_OFFSET = -1L;
    public const long MAGIC_NUMBER = 0x5452454558494646L; // "SUFFIXTR"

    // Header Layout Offsets
    public const int HEADER_OFFSET_MAGIC = 0;
    public const int HEADER_OFFSET_VERSION = 8;
    public const int HEADER_OFFSET_ROOT = 16;
    public const int HEADER_OFFSET_TEXT_OFF = 24;
    public const int HEADER_OFFSET_TEXT_LEN = 32;
    public const int HEADER_OFFSET_NODE_COUNT = 36;
    public const int HEADER_OFFSET_SIZE = 44;

    // Node Layout Offsets
    public const int OFFSET_START = 0;
    public const int OFFSET_END = 4;
    public const int OFFSET_SUFFIX_LINK = 8;
    public const int OFFSET_DEPTH = 16;
    public const int OFFSET_LEAF_COUNT = 20;
    public const int OFFSET_CHILDREN_HEAD = 24;
    public const int OFFSET_CHILD_COUNT = 32;

    // Child Entry Layout Offsets (sorted contiguous array: Key(4) + ChildNodeOffset(8) = 12 bytes)
    public const int CHILD_OFFSET_KEY = 0;
    public const int CHILD_OFFSET_NODE = 4;
}
