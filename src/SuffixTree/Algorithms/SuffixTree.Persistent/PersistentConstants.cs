namespace SuffixTree.Persistent;

public static class PersistentConstants
{
    public const int NODE_SIZE = 40;
    public const int CHILD_ENTRY_SIZE = 12;
    public const int HEADER_SIZE = 48;
    public const uint BOUNDLESS = uint.MaxValue;
    public const uint TERMINATOR_KEY = uint.MaxValue;  // Same bit pattern as -1 when cast to int
    public const long NULL_OFFSET = -1L;
    public const long MAGIC_NUMBER = 0x5452454558494646L; // "SUFFIXTR"
    public const int CURRENT_VERSION = 3;

    // Header Layout (48 bytes, packed, naturally aligned)
    //  0-7  : MAGIC      (int64)
    //  8-11 : VERSION    (int32)
    // 12-15 : TEXT_LEN   (int32)
    // 16-23 : ROOT       (int64)
    // 24-31 : TEXT_OFF   (int64)
    // 32-35 : NODE_COUNT (int32)
    // 36-39 : reserved
    // 40-47 : SIZE       (int64)
    public const int HEADER_OFFSET_MAGIC = 0;
    public const int HEADER_OFFSET_VERSION = 8;
    public const int HEADER_OFFSET_TEXT_LEN = 12;
    public const int HEADER_OFFSET_ROOT = 16;
    public const int HEADER_OFFSET_TEXT_OFF = 24;
    public const int HEADER_OFFSET_NODE_COUNT = 32;
    public const int HEADER_OFFSET_SIZE = 40;

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
