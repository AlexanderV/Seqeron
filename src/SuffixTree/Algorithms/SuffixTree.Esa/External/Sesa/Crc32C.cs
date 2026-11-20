using System.Buffers.Binary;
using System.Runtime.Intrinsics.Arm;
using System.Runtime.Intrinsics.X86;

namespace SuffixTree.Esa.External;

/// <summary>
/// CRC-32C (Castagnoli, polynomial 0x1EDC6F41, reflected 0x82F63B78, initial value and final xor 0xFFFFFFFF) as ADR-003 section 6 wants it for the section and header checksums: the hardware instruction where the CPU
/// has it (SSE 4.2 on x86-64, the CRC extension on arm64), a slice-by-8 managed table path everywhere else; both give the same value. The functions chain: <c>Append(Append(0, a), b)</c> is the CRC of <c>a</c> followed
/// by <c>b</c>, and <c>Append(0, empty)</c> is 0.
/// </summary>
internal static unsafe class Crc32C
{
    private const uint Reflected = 0x82F63B78;

    /// <summary>Eight tables of 256 entries: table k holds the CRC of a byte followed by k zero bytes (slice-by-8).</summary>
    private static readonly uint[] Tables = BuildTables();

    /// <summary>True where <see cref="AppendHardware(uint, ReadOnlySpan{byte})"/> can run.</summary>
    internal static bool HardwareAccelerated { get; } = Sse42.X64.IsSupported || Crc32.Arm64.IsSupported;

    /// <summary>The CRC of <paramref name="data"/> continued from <paramref name="crc"/> (0 to start), on the hardware path where there is one.</summary>
    internal static uint Append(uint crc, ReadOnlySpan<byte> data)
    {
        fixed (byte* p = data) return Append(crc, p, data.Length);
    }

    /// <summary>The same for memory at <paramref name="data"/> (a mapped view); <paramref name="length"/> may exceed 2^31.</summary>
    internal static uint Append(uint crc, byte* data, long length) => HardwareAccelerated ? Hardware(crc, data, length) : Table(crc, data, length);

    /// <summary>The managed slice-by-8 path.</summary>
    internal static uint AppendTable(uint crc, ReadOnlySpan<byte> data)
    {
        fixed (byte* p = data) return Table(crc, p, data.Length);
    }

    /// <summary>The hardware path.</summary>
    /// <exception cref="PlatformNotSupportedException">The CPU has no CRC-32C instruction.</exception>
    internal static uint AppendHardware(uint crc, ReadOnlySpan<byte> data)
    {
        if (!HardwareAccelerated) throw new PlatformNotSupportedException("this CPU has no CRC-32C instruction");
        fixed (byte* p = data) return Hardware(crc, p, data.Length);
    }

    /// <summary>The CRC of <paramref name="data"/>.</summary>
    internal static uint Compute(ReadOnlySpan<byte> data) => Append(0, data);

    private static uint[] BuildTables()
    {
        var t = new uint[8 * 256];
        for (uint i = 0; i < 256; i++)
        {
            uint c = i;
            for (int k = 0; k < 8; k++) c = (c & 1) != 0 ? (c >> 1) ^ Reflected : c >> 1;
            t[i] = c;
        }
        for (int k = 1; k < 8; k++)
        {
            for (int i = 0; i < 256; i++)
            {
                uint previous = t[((k - 1) * 256) + i];
                t[(k * 256) + i] = (previous >> 8) ^ t[previous & 0xFF];
            }
        }
        return t;
    }

    private static uint Table(uint crc, byte* data, long length)
    {
        uint c = ~crc;
        uint[] t = Tables;
        long at = 0;
        for (; at + 8 <= length; at += 8)
        {
            uint lo = BinaryPrimitives.ReadUInt32LittleEndian(new ReadOnlySpan<byte>(data + at, 4)) ^ c;
            uint hi = BinaryPrimitives.ReadUInt32LittleEndian(new ReadOnlySpan<byte>(data + at + 4, 4));
            c = t[(7 * 256) + (lo & 0xFF)] ^ t[(6 * 256) + ((lo >> 8) & 0xFF)] ^ t[(5 * 256) + ((lo >> 16) & 0xFF)] ^ t[(4 * 256) + (lo >> 24)]
                ^ t[(3 * 256) + (hi & 0xFF)] ^ t[(2 * 256) + ((hi >> 8) & 0xFF)] ^ t[256 + ((hi >> 16) & 0xFF)] ^ t[hi >> 24];
        }
        for (; at < length; at++) c = t[(c ^ data[at]) & 0xFF] ^ (c >> 8);
        return ~c;
    }

    private static uint Hardware(uint crc, byte* data, long length)
    {
        uint c = ~crc;
        long at = 0;
        if (Sse42.X64.IsSupported)
        {
            ulong c64 = c;
            for (; at + 32 <= length; at += 32)
            {
                c64 = Sse42.X64.Crc32(c64, BinaryPrimitives.ReadUInt64LittleEndian(new ReadOnlySpan<byte>(data + at, 8)));
                c64 = Sse42.X64.Crc32(c64, BinaryPrimitives.ReadUInt64LittleEndian(new ReadOnlySpan<byte>(data + at + 8, 8)));
                c64 = Sse42.X64.Crc32(c64, BinaryPrimitives.ReadUInt64LittleEndian(new ReadOnlySpan<byte>(data + at + 16, 8)));
                c64 = Sse42.X64.Crc32(c64, BinaryPrimitives.ReadUInt64LittleEndian(new ReadOnlySpan<byte>(data + at + 24, 8)));
            }
            for (; at + 8 <= length; at += 8) c64 = Sse42.X64.Crc32(c64, BinaryPrimitives.ReadUInt64LittleEndian(new ReadOnlySpan<byte>(data + at, 8)));
            c = (uint)c64;
            for (; at < length; at++) c = Sse42.Crc32(c, data[at]);
        }
        else
        {
            c = ArmLoops(new ArmCrcOps(), c, data, length);
        }
        return ~c;
    }

    /// <summary>The two CRC-32C steps of the arm64 path (raw updates, without the initial and the final inversion).</summary>
    internal interface IArmCrcOps
    {
        uint Step64(uint crc, ulong value);

        uint Step8(uint crc, byte value);
    }

    private readonly struct ArmCrcOps : IArmCrcOps
    {
        public uint Step64(uint crc, ulong value) => Crc32.Arm64.ComputeCrc32C(crc, value);

        public uint Step8(uint crc, byte value) => Crc32.ComputeCrc32C(crc, value);
    }

    /// <summary>
    /// The loops of the arm64 path, over the steps of <paramref name="ops"/> (a struct, so that the production call is specialised and inlined): eight bytes at a time, then the tail byte by byte. A test supplies the steps in
    /// software, so that the loops run, and are checked against the table path, on a CPU that has no such instruction.
    /// </summary>
    internal static uint ArmLoops<TOps>(TOps ops, uint c, byte* data, long length)
        where TOps : struct, IArmCrcOps
    {
        long at = 0;
        for (; at + 8 <= length; at += 8) c = ops.Step64(c, BinaryPrimitives.ReadUInt64LittleEndian(new ReadOnlySpan<byte>(data + at, 8)));
        for (; at < length; at++) c = ops.Step8(c, data[at]);
        return c;
    }
}
