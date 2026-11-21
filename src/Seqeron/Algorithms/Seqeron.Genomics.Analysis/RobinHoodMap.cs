namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Port of the flat <c>robin_hood::unordered_map&lt;K, uint32_t&gt;</c> of robin_hood.h 3.9.1 (martinus/robin-hood-hashing, the header bundled with Mash 2.3: <c>HashSet.h</c> keeps
/// the sketch's hashes and counts in it, <c>HashSet::toCounts</c> lists them in the ITERATION ORDER of the table): <c>robin_hood::hash</c> (<c>hash_int</c>: two 64-bit multiplications and a 32-bit rotation, the
/// integer key widened to 64 bits), a power-of-two table of at least 8 buckets with the 80 % load factor and the per-bucket info byte (the low 5 bits carry hash bits, the rest the distance from the home bucket;
/// <c>try_increase_info</c> halves the info increment when the distance would overflow), insertion by <c>doCreateByKey</c> (skip forward while the info is less, <c>shiftUp</c> over the occupied run),
/// the rehash by <c>insert_move</c> in the old bucket order, and the backward-shift erase (<c>shiftDown</c>). Only the operations Mash performs are ported (insert-or-add, count, erase, clear, enumeration in bucket order).
/// B09 A5-10b3 (F194).
/// </summary>
internal sealed class RobinHoodMap
{
    private const int InitialNumElements = 8;
    private const int InitialInfoNumBits = 5;
    private const uint InitialInfoInc = 1U << InitialInfoNumBits;
    private const uint InfoMask = InitialInfoInc - 1U;
    private const int MaxLoadFactor100 = 80;

    private ulong[] _keys = [];
    private uint[] _values = [];
    private byte[] _info = new byte[16];        // the never-allocated table reads as zeros (it points at its own mask)
    private ulong _numElements;
    private ulong _mask;
    private ulong _maxNumElementsAllowed;
    private uint _infoInc = InitialInfoInc;
    private int _infoHashShift;

    /// <summary>The number of elements.</summary>
    public int Count => (int)_numElements;

    /// <summary><c>robin_hood::hash_int</c> (also the hash of an unsigned 32-bit or 64-bit key).</summary>
    internal static ulong HashInt(ulong x)
    {
        ulong h1 = x * 0xA24BAED4963EE407UL;
        ulong h2 = ulong.RotateRight(x, 32) * 0x9FB21C651E98DF25UL;
        return ulong.RotateRight(h1 + h2, 32);
    }

    /// <summary>The value of the key, 0 when it is absent (<c>HashSet::count</c>).</summary>
    public uint Get(ulong key)
    {
        int idx = FindIdx(key);
        return idx < 0 ? 0 : _values[idx];
    }

    /// <summary>Whether the key is present.</summary>
    public bool Contains(ulong key) => FindIdx(key) >= 0;

    /// <summary><c>HashSet::insert</c>: the key's value grows by <paramref name="count"/>, a new key is created (<c>operator[]</c>) with the value <paramref name="count"/>.</summary>
    public void Insert(ulong key, uint count)
    {
        int idx = FindIdx(key);
        if (idx >= 0)
        {
            _values[idx] += count;
            return;
        }

        int created = CreateByKey(key);      // (the arrays are replaced inside: index after the call)
        _values[created] = count;
    }

    /// <summary>Removes the key (the backward-shift erase); true when it was present.</summary>
    public bool Erase(ulong key)
    {
        KeyToIdx(key, out ulong idx, out uint info);
        do
        {
            if (info == _info[idx] && _keys[idx] == key)
            {
                ShiftDown(idx);
                _numElements--;
                return true;
            }

            idx++;
            info += _infoInc;
        } while (info <= _info[idx]);

        return false;
    }

    /// <summary>Removes every element and keeps the allocation (<c>Table::clear</c>).</summary>
    public void Clear()
    {
        if (_numElements == 0)
            return;
        Array.Clear(_keys);
        Array.Clear(_values);
        ulong withBuffer = NumElementsWithBuffer(_mask + 1);
        Array.Clear(_info, 0, _info.Length);
        _info[withBuffer] = 1;
        _numElements = 0;
        _infoInc = InitialInfoInc;
        _infoHashShift = 0;
    }

    /// <summary>The values in the order of the table's iterator (bucket order), as <c>HashSet::toCounts</c> lists them.</summary>
    public List<uint> ValuesInIterationOrder()
    {
        var list = new List<uint>((int)_numElements);
        Enumerate((_, v) => list.Add(v));
        return list;
    }

    /// <summary>The keys in the order of the table's iterator (<c>HashSet::toHashList</c>).</summary>
    public List<ulong> KeysInIterationOrder()
    {
        var list = new List<ulong>((int)_numElements);
        Enumerate((k, _) => list.Add(k));
        return list;
    }

    private void Enumerate(Action<ulong, uint> visit)
    {
        if (_numElements == 0)
            return;
        ulong withBuffer = NumElementsWithBuffer(_mask + 1);
        for (ulong i = 0; i < withBuffer; i++)
            if (_info[i] != 0)
                visit(_keys[i], _values[i]);
    }

    private static ulong CalcMaxNumElementsAllowed(ulong maxElements) => maxElements * MaxLoadFactor100 / 100;

    private static ulong NumElementsWithBuffer(ulong numElements) => numElements + Math.Min(CalcMaxNumElementsAllowed(numElements), 0xFFUL);

    private void KeyToIdx(ulong key, out ulong idx, out uint info)
    {
        ulong h = HashInt(key);
        info = _infoInc + (uint)((h & InfoMask) >> _infoHashShift);
        idx = (h >> InitialInfoNumBits) & _mask;
    }

    private int FindIdx(ulong key)
    {
        if (_numElements == 0)
            return -1;
        KeyToIdx(key, out ulong idx, out uint info);
        do
        {
            if (info == _info[idx] && _keys[idx] == key)
                return (int)idx;
            idx++;
            info += _infoInc;
            if (info == _info[idx] && _keys[idx] == key)
                return (int)idx;
            idx++;
            info += _infoInc;
        } while (info <= _info[idx]);

        return -1;
    }

    private void ShiftUp(ulong startIdx, ulong insertionIdx)
    {
        ulong idx = startIdx;
        _keys[idx] = _keys[idx - 1];
        _values[idx] = _values[idx - 1];
        while (--idx != insertionIdx)
        {
            _keys[idx] = _keys[idx - 1];
            _values[idx] = _values[idx - 1];
        }

        idx = startIdx;
        while (idx != insertionIdx)
        {
            _info[idx] = (byte)(_info[idx - 1] + _infoInc);
            if (_info[idx] + _infoInc > 0xFF)
                _maxNumElementsAllowed = 0;
            idx--;
        }
    }

    private void ShiftDown(ulong idx)
    {
        while (_info[idx + 1] >= 2 * _infoInc)
        {
            _info[idx] = (byte)(_info[idx + 1] - _infoInc);
            _keys[idx] = _keys[idx + 1];
            _values[idx] = _values[idx + 1];
            idx++;
        }

        _info[idx] = 0;
        _keys[idx] = 0;
        _values[idx] = 0;
    }

    // doCreateByKey for a key that is absent: returns the bucket of the new element
    private int CreateByKey(ulong key)
    {
        while (true)
        {
            KeyToIdx(key, out ulong idx, out uint info);
            while (info < _info[idx])
            {
                idx++;
                info += _infoInc;
            }

            while (info == _info[idx])
            {
                idx++;
                info += _infoInc;
            }

            if (_numElements >= _maxNumElementsAllowed)
            {
                IncreaseSize();
                continue;
            }

            ulong insertionIdx = idx;
            uint insertionInfo = info;
            if (insertionInfo + _infoInc > 0xFF)
                _maxNumElementsAllowed = 0;
            while (_info[idx] != 0)
            {
                idx++;
                info += _infoInc;
            }

            if (idx != insertionIdx)
                ShiftUp(idx, insertionIdx);
            _keys[insertionIdx] = key;
            _values[insertionIdx] = 0;
            _info[insertionIdx] = (byte)insertionInfo;
            _numElements++;
            return (int)insertionIdx;
        }
    }

    private bool TryIncreaseInfo()
    {
        if (_infoInc <= 2)
            return false;
        _infoInc >>= 1;
        _infoHashShift++;
        ulong withBuffer = NumElementsWithBuffer(_mask + 1);
        for (ulong i = 0; i < withBuffer; i += 8)
        {
            ulong val = BitConverter.ToUInt64(_info, (int)i);
            val = (val >> 1) & 0x7f7f7f7f7f7f7f7fUL;
            BitConverter.TryWriteBytes(_info.AsSpan((int)i, 8), val);
        }

        _info[withBuffer] = 1;
        _maxNumElementsAllowed = CalcMaxNumElementsAllowed(_mask + 1);
        return true;
    }

    private void IncreaseSize()
    {
        if (_mask == 0)
        {
            InitData(InitialNumElements);
            return;
        }

        ulong maxAllowed = CalcMaxNumElementsAllowed(_mask + 1);
        if (_numElements < maxAllowed && TryIncreaseInfo())
            return;
        if (_numElements * 2 < CalcMaxNumElementsAllowed(_mask + 1))
            throw new InvalidOperationException("robin_hood::map overflow");
        RehashPowerOfTwo((_mask + 1) * 2);
    }

    private void InitData(ulong maxElements)
    {
        _numElements = 0;
        _mask = maxElements - 1;
        _maxNumElementsAllowed = CalcMaxNumElementsAllowed(maxElements);
        ulong withBuffer = NumElementsWithBuffer(maxElements);
        _keys = new ulong[withBuffer];
        _values = new uint[withBuffer];
        _info = new byte[withBuffer + 8];
        _info[withBuffer] = 1;
        _infoInc = InitialInfoInc;
        _infoHashShift = 0;
    }

    private void RehashPowerOfTwo(ulong numBuckets)
    {
        var oldKeys = _keys;
        var oldValues = _values;
        var oldInfo = _info;
        ulong oldWithBuffer = _mask == 0 ? 1 : NumElementsWithBuffer(_mask + 1);
        InitData(numBuckets);
        if (oldWithBuffer > 1)
            for (ulong i = 0; i < oldWithBuffer; i++)
                if (oldInfo[i] != 0)
                    InsertMove(oldKeys[i], oldValues[i]);
    }

    private void InsertMove(ulong key, uint value)
    {
        if (_maxNumElementsAllowed == 0 && !TryIncreaseInfo())
            throw new InvalidOperationException("robin_hood::map overflow");
        KeyToIdx(key, out ulong idx, out uint info);
        while (info <= _info[idx])
        {
            idx++;
            info += _infoInc;
        }

        ulong insertionIdx = idx;
        uint insertionInfo = (byte)info;
        if (insertionInfo + _infoInc > 0xFF)
            _maxNumElementsAllowed = 0;
        while (_info[idx] != 0)
        {
            idx++;
            info += _infoInc;
        }

        if (idx != insertionIdx)
            ShiftUp(idx, insertionIdx);
        _keys[insertionIdx] = key;
        _values[insertionIdx] = value;
        _info[insertionIdx] = (byte)insertionInfo;
        _numElements++;
    }
}
