// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.Repeats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Emboss.Repeats;

using Seqeron.Genomics.Analysis;

/// <summary>
/// Port of the Tarhio-Ukkonen-Bleasby k-mismatch search of EMBOSS <c>nucleus/embpat.c</c> (<c>embPatTUBInit</c> / <c>embPatTUBSearch</c>), the engine <c>embPatGetTypeII</c> picks (type 6) for a pattern without a range
/// that has a class <c>[..]</c> or an exclusion <c>{..}</c> when mismatches are allowed (<c>fuzzpro</c> / <c>fuzztran</c> / <c>fuzznuc -pmismatch</c>, B09 A5-8c2, F195). The skip table knows the letters A-Z only: a text symbol that is no
/// letter (the stop <c>*</c> of a translation) gets the default skip although an exclusion or <c>x</c> matches it, so the search may jump over windows that contain it (the output differs from a plain Hamming scan there; measured on fuzztran).
/// </summary>
public static class EmbossTubSearch
{
    private const int Alpha = 256;

    /// <summary><c>embPatTUBInit</c>: the skip table <c>skipm[m][256]</c> for the pattern text (letters, <c>?</c>, <c>[..]</c>, <c>{..}</c>), <paramref name="m"/> real elements, <paramref name="k"/> mismatches.</summary>
    public static uint[][] Init(string pat, int m, int k)
    {
        var skipm = new uint[m][];
        for (int i = 0; i < m; i++)
            skipm[i] = new uint[Alpha];
        var ready = new uint[Alpha];
        for (int i = 0; i < Alpha; i++)
        {
            ready[i] = (uint)m;
            for (int j = m - k - 1; j < m; j++)
                skipm[j][i] = (uint)Math.Max(m - k - 1, 1);
        }

        int p = pat.Length - 1;
        if (pat[p] is '}' or ']')
        {
            while (pat[p] != '{' && pat[p] != '[')
                p--;
            p--;
        }
        else
        {
            p--;
        }

        for (int j = m - 2; j > -1; --j)
        {
            int jj = m - k - 1;
            int x = Math.Max(j + 1, jj);
            if (pat[p] == '?')
            {
                for (int z = 'A'; z <= 'Z'; ++z)
                {
                    for (int i = (int)ready[z] - 1; i >= x; --i)
                        skipm[i][z] = (uint)(i - j);
                    ready[z] = (uint)x;
                }

                p--;
                continue;
            }

            if (pat[p] == ']')
            {
                p--;
                while (pat[p] != '[')
                {
                    int c = pat[p];
                    for (int i = (int)ready[c] - 1; i >= x; --i)
                        skipm[i][c] = (uint)(i - j);
                    ready[c] = (uint)x;
                    p--;
                }

                p--;
                continue;
            }

            if (pat[p] == '}')
            {
                int s = --p;
                for (int z = 'A'; z <= 'Z'; ++z)
                {
                    int q = s;
                    bool flag = false;
                    while (pat[q] != '{')
                    {
                        if (pat[q] == z)
                        {
                            flag = true;
                            break;
                        }

                        q--;
                    }

                    if (!flag)
                    {
                        for (int i = (int)ready[z] - 1; i >= x; --i)
                            skipm[i][z] = (uint)(i - j);
                        ready[z] = (uint)x;
                    }
                }

                while (pat[p] != '{')
                    p--;
                p--;
                continue;
            }

            int ch = pat[p];
            for (int i = (int)ready[ch] - 1; i >= x; --i)
                skipm[i][ch] = (uint)(i - j);
            ready[ch] = (uint)x;
            p--;
        }

        return skipm;
    }

    /// <summary><c>embPatTUBSearch</c>: the windows (start, length m, mismatches) in the order the program pushes them.</summary>
    public static List<(int Start, int Length, int Mismatches)> Search(string pat, string text, uint[][] skipm, int m, int k, bool amino, bool carboxyl)
    {
        var hits = new List<(int, int, int)>();
        int slen = text.Length;
        int plen = pat.Length;
        long i = m - 1;
        while (i < slen)
        {
            long h = i;
            int j = m - 1;
            int p = plen - 1;
            uint skip = (uint)(m - k);
            int mm = 0;
            while (j > -1 && mm <= k)
            {
                int jj = m - k - 1;
                int a = text[(int)h];
                if (j >= jj)
                    skip = Math.Min(skip, skipm[j][a & 0xFF]);
                if (pat[p] != '?')
                {
                    if (pat[p] == ']')
                    {
                        bool flag = false;
                        p--;
                        while (pat[p] != '[')
                        {
                            if (a == pat[p])
                                flag = true;
                            p--;
                        }

                        if (!flag)
                            mm++;
                    }
                    else if (pat[p] == '}')
                    {
                        bool flag = false;
                        p--;
                        while (pat[p] != '{')
                        {
                            if (a == pat[p])
                                flag = true;
                            p--;
                        }

                        if (flag)
                            mm++;
                    }
                    else if (a != pat[p])
                    {
                        mm++;
                    }
                }

                p--;
                h--;
                j--;
            }

            if (mm <= k)
            {
                if (amino && h + 1 != 0)
                    return hits;
                if (!carboxyl || h + 1 == slen - m)
                    hits.Add(((int)(h + 1), m, mm));
            }

            i += skip;
        }

        return hits;
    }
}
