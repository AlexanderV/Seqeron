namespace SuffixTree.Esa.External;

/// <summary>
/// SplitMix64 (Steele, Lea, Flood 2014; Vigna's public-domain constants) with the rule <c>Next() % bound</c> for a bound: the generator of the prototype's splitter sampling and of the pruned difference covers, so that the
/// bucket trees and the covers are those that E3, E3b and E6 verified. Not thread-safe (every task has its own).
/// </summary>
internal sealed class SeededRandom(ulong seed)
{
    private ulong _state = seed;

    internal ulong Next()
    {
        unchecked
        {
            _state += 0x9E3779B97F4A7C15UL;
            ulong z = _state;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
            return z ^ (z >> 31);
        }
    }

    /// <summary>A value in <c>[0, bound)</c>: <c>Next() % bound</c>.</summary>
    internal int NextInt(int bound)
    {
        ArgumentOutOfRangeException.ThrowIfLessThanOrEqual(bound, 0);
        return (int)(Next() % (ulong)bound);
    }
}

/// <summary>
/// LIFO work list with termination detection (depth-first order keeps the number of pending tasks, hence of scratch files, small): a task that splits a bucket pushes its children before it calls <see cref="Done"/>, so the
/// list is finished only when it is empty and no task is running. A failure ends every <see cref="Take"/> at once and keeps the first exception for the caller.
/// </summary>
internal sealed class WorkPool<T>
    where T : class
{
    private readonly Stack<T> _stack = new();
    private int _active;
    private Exception? _error;

    /// <summary>The first failure that was reported, or null.</summary>
    internal Exception? Error
    {
        get
        {
            lock (_stack) return _error;
        }
    }

    internal void Push(T task)
    {
        lock (_stack)
        {
            _stack.Push(task);
            Monitor.Pulse(_stack);
        }
    }

    /// <summary>Blocks until a task is available (true; the caller must call <see cref="Done"/>) or all work is finished or failed (false).</summary>
    internal bool Take(out T task)
    {
        lock (_stack)
        {
            while (true)
            {
                if (_error is null && _stack.Count > 0)
                {
                    task = _stack.Pop();
                    _active++;
                    return true;
                }
                if (_error is not null || _active == 0)
                {
                    task = null!;
                    Monitor.PulseAll(_stack);
                    return false;
                }
                Monitor.Wait(_stack);
            }
        }
    }

    internal void Done()
    {
        lock (_stack)
        {
            _active--;
            Monitor.PulseAll(_stack);
        }
    }

    internal void Fail(Exception e)
    {
        lock (_stack)
        {
            _error ??= e;
            Monitor.PulseAll(_stack);
        }
    }
}
