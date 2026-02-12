using System.Collections.Generic;
using System.Text;

namespace Seqeron.Genomics.Tests.Helpers;

/// <summary>
/// Shared test helpers for reverse complement and other sequence operations.
/// </summary>
internal static class TestSequenceHelpers
{
    /// <summary>
    /// Gets reverse complement of a DNA sequence (supports N/n).
    /// </summary>
    internal static string GetReverseComplement(string sequence)
    {
        var complement = new Dictionary<char, char>
        {
            ['A'] = 'T',
            ['T'] = 'A',
            ['C'] = 'G',
            ['G'] = 'C',
            ['a'] = 't',
            ['t'] = 'a',
            ['c'] = 'g',
            ['g'] = 'c',
            ['N'] = 'N',
            ['n'] = 'n'
        };
        var sb = new StringBuilder(sequence.Length);
        for (int i = sequence.Length - 1; i >= 0; i--)
        {
            char c = sequence[i];
            sb.Append(complement.GetValueOrDefault(c, c));
        }
        return sb.ToString();
    }
}
