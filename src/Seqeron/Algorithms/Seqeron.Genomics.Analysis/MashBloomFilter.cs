namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The Bloom filter of <c>mash sketch -b</c> (<c>bloom_filter.hpp</c> of Mash 2.3 as <c>MinHashHeap</c> builds it; B09 A5-10b2, F197): <c>MinHashHeap</c> asks for the optimal parameters with a projected count of 10^9 and a false-positive
/// probability of 0, for which the search of <c>compute_optimal_parameters</c> gives -infinity bits for one hash function; the conversion of -infinity to <c>unsigned long long</c> is an out-of-range value that is cut to <c>maximum_size</c> = 8 x the
/// memory bound, so the filter has ONE hash function (the salt 0xAAAAAAAA squared plus the low 32 bits of the seeded random seed) over exactly <c>bound x 8</c> bits (measured on the x86-64 build: files byte-equal). The hash is the AP hash of the 4 (32-bit hash) or 8 (64-bit hash) little-endian
/// bytes of the k-mer hash.
/// </summary>
internal sealed class MashBloomFilter
{
    private readonly byte[] _table;
    private readonly ulong _tableBits;
    private readonly uint _salt;

    /// <summary>Creates the filter for a memory bound of <paramref name="memoryBoundBytes"/> bytes.</summary>
    /// <exception cref="ArgumentOutOfRangeException">The bound is not positive.</exception>
    /// <exception cref="NotSupportedException">The table would exceed the 2 GiB of a .NET array (the program allocates it).</exception>
    public MashBloomFilter(long memoryBoundBytes)
    {
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(memoryBoundBytes, 0);
        if (memoryBoundBytes > Array.MaxLength)
            throw new NotSupportedException($"A Bloom filter of {memoryBoundBytes} bytes does not fit into a .NET array (limit {Array.MaxLength}).");
        _table = new byte[memoryBoundBytes];
        _tableBits = (ulong)memoryBoundBytes * 8UL;
        unchecked
        {
            // random_seed_ = (0xA5A5A5A55A5A5A5A * 0xA5A5A5A5) + 1 (64 bit); only its low 32 bits enter the salt
            uint seed = 0x5A5A5A5Au * 0xA5A5A5A5u + 1u;
            _salt = 0xAAAAAAAAu * 0xAAAAAAAAu + seed;
        }
    }

    /// <summary><c>bloom_filter::contains</c> of a hash (8 bytes for a 64-bit sketch, else 4).</summary>
    public bool Contains(ulong hash, bool use64)
    {
        ulong index = Hash(hash, use64) % _tableBits;
        return (_table[index / 8] & (1 << (int)(index % 8))) != 0;
    }

    /// <summary><c>bloom_filter::insert</c>.</summary>
    public void Insert(ulong hash, bool use64)
    {
        ulong index = Hash(hash, use64) % _tableBits;
        _table[index / 8] |= (byte)(1 << (int)(index % 8));
    }

    // hash_ap over the little-endian bytes of the k-mer hash
    private uint Hash(ulong key, bool use64)
    {
        unchecked
        {
            uint hash = _salt;
            if (use64)
            {
                uint i1 = (uint)key, i2 = (uint)(key >> 32);
                hash ^= (hash << 7) ^ i1 * (hash >> 3) ^ ~((hash << 11) + (i2 ^ (hash >> 5)));
                return hash;
            }

            uint i = (uint)key;
            hash ^= ~((hash << 11) + (i ^ (hash >> 5)));
            return hash;
        }
    }
}
