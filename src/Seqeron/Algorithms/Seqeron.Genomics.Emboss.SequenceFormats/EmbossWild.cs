// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

/// <summary>
/// The wildcard comparisons of the AJAX string library (<c>ajstr.c</c>: <c>ajCharCmpWild</c>, <c>ajCharCmpWildCase</c>, <c>ajCharMatchWildNextC</c>, <c>ajCharMatchWildNextCaseC</c>, <c>ajCharMatchWildWordC</c>, <c>ajCharMatchWildWordCaseC</c>)
/// with the C strings of the library (the text up to a NUL character): <c>*</c> matches any run, <c>?</c> one character; the case-less forms compare with <c>toupper</c> of the C locale (ASCII letters only).
/// A5-8c3h, B09 F255.
/// </summary>
public static class EmbossWild
{
    // C strings: the text from the first NUL character on is not seen
    private static string CStr(string s)
    {
        int nul = s.IndexOf('\0');
        return nul >= 0 ? s[..nul] : s;
    }

    private static char Up(char c) => c is >= 'a' and <= 'z' ? (char)(c - 32) : c;

    private static bool IsSpace(char c) => c is ' ' or '\t' or '\n' or '\v' or '\f' or '\r';

    private static bool IsAlnum(char c) => c is (>= 'a' and <= 'z') or (>= 'A' and <= 'Z') or (>= '0' and <= '9');

    /// <summary><c>ajStrMatchWildS</c>: <paramref name="text"/> matches the wildcard pattern.</summary>
    public static bool Match(string text, string pattern) => Cmp(CStr(text), 0, CStr(pattern), 0, false) == 0;

    /// <summary><c>ajStrMatchWildCaseS</c>.</summary>
    public static bool MatchCase(string text, string pattern) => Cmp(CStr(text), 0, CStr(pattern), 0, true) == 0;

    // ajCharCmpWild / ajCharCmpWildCase: txt at q, pattern at p (both are C strings, the end is the end of the string)
    private static int Cmp(string txt, int q, string pat, int p, bool nocase)
    {
        if (p >= pat.Length && q >= txt.Length)
            return 0;
        if (p >= pat.Length)
            return -1;
        while (p < pat.Length)
        {
            if (q >= txt.Length && pat[p] != '*')
                return 1;
            switch (pat[p])
            {
                case '?':
                    p++;
                    q++;
                    break;
                case '*':
                    p++;
                    if (p >= pat.Length)
                        return 0;
                    if (q >= txt.Length)
                        return Cmp(txt, q, pat, p, nocase);
                    while (q < txt.Length)
                    {
                        if (Cmp(txt, q, pat, p, nocase) == 0)
                            return 0;
                        q++;
                    }

                    return 1;
                default:
                    {
                        char a = nocase ? Up(pat[p]) : pat[p], b = nocase ? Up(txt[q]) : txt[q];
                        if (a != b)
                            return a > b ? -1 : 1;
                        p++;
                        if (q < txt.Length)
                            q++;
                        break;
                    }
            }
        }

        if (p < pat.Length)
            return -1;
        if (q < txt.Length)
            return 1;
        return 0;
    }

    /// <summary><c>ajStrMatchWildWordCaseS</c>: the pattern matches a word of the text (a run of letters and digits), any case.</summary>
    public static bool MatchWordCase(string text, string pattern)
    {
        string txt = CStr(text), pat = CStr(pattern);
        if (pat.Length == 0 && txt.Length == 0)
            return true;
        if (pat.Length == 0)
            return false;
        bool word = false;
        for (int i = 0; i < txt.Length; i++)
        {
            if (!word)
            {
                if (IsAlnum(txt[i]))
                {
                    word = true;
                    if (Next(txt, i, pat, 0, true))
                        return true;
                }
            }
            else if (!IsAlnum(txt[i]))
            {
                word = false;
            }
        }

        return false;
    }

    // ajCharMatchWildNextC / ajCharMatchWildNextCaseC: the pattern ends at its first white space; the match must end at a word boundary
    private static bool Next(string txt, int q, string pat, int p, bool nocase)
    {
        char lastch = '\0';
        if (p >= pat.Length && q >= txt.Length)
            return true;
        if (p >= pat.Length)
            return false;
        while (p < pat.Length && !IsSpace(pat[p]))
        {
            if (q >= txt.Length && pat[p] != '*')
                return false;
            switch (pat[p])
            {
                case '?':
                    lastch = q < txt.Length ? txt[q] : '\0';
                    p++;
                    q++;
                    break;
                case '*':
                    p++;
                    if (p >= pat.Length)
                        return true;
                    if (q >= txt.Length)
                    {
                        int savecp = p;
                        while (p < pat.Length && pat[p] == '*')
                            savecp = p++;
                        if (p >= pat.Length)
                            return true;
                        return Next(txt, q, pat, savecp, nocase);
                    }

                    while (q < txt.Length)
                    {
                        if (Next(txt, q, pat, p, nocase))
                            return true;
                        q++;
                    }

                    return false;
                default:
                    {
                        char a = nocase ? Up(pat[p]) : pat[p], b = nocase ? Up(txt[q]) : txt[q];
                        if (a != b)
                            return false;
                        p++;
                        if (q < txt.Length)
                        {
                            lastch = txt[q];
                            q++;
                        }

                        break;
                    }
            }
        }

        if (!IsAlnum(lastch))
            return false;
        if (p < pat.Length)
            return false;
        if (q < txt.Length)
            return !IsAlnum(txt[q]);
        return true;
    }
}
