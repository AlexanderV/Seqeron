// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

/// <summary>
/// The command-line qualifier resolution of the EMBOSS 6.6.0 ACD library (<c>ajacd.c</c>: <c>acdArgsParse</c>, <c>acdIsQual</c>, <c>acdQualParse</c>, <c>acdFindQual</c>, <c>acdFindQualDetail</c>, <c>acdFindQualAssoc</c>,
/// <c>acdFindQualMaster</c>, <c>acdTestQualC</c>; EMBOSS-6.6.0 at emboss.open-bio.org via raw.githubusercontent.com mirrors, read for audit A6-2, B09 F177) for the qualifiers of the sequence inputs of a program, run on the qualifier text
/// that follows the sequence arguments of the command line: <c>needle a.fa b.fa -sbegin 3 ...</c> (the text is what comes after <c>a.fa b.fa</c>). The rules the program applies, all reproduced:
/// <list type="bullet">
/// <item>A token starting with <c>-</c> or <c>/</c> (optionally doubled, <c>--sbegin</c>) is a qualifier when its name is letters / digits / underscore, or when it contains <c>=</c> (<c>-sbegin=5</c>, also without the dash); names are case-insensitive.</item>
/// <item>The name is split into <c>name[_master][number]</c>; <c>no</c> before a name is the negation of a boolean.</item>
/// <item>A name matches an item exactly, or as a prefix when that is unique: the prefix is tested against EVERY item of the program (parameters, qualifiers, the associated qualifiers of each parameter, the general qualifiers), so <c>-se</c> is
/// ambiguous in a one-sequence program (<c>-send</c> / <c>-sequence</c>) and <c>-s</c> always; the same name of two numbered inputs (<c>-sbegin1</c> / <c>-sbegin2</c>) is one name.</item>
/// <item>The qualifiers of a sequence are tied to the parameter that PRECEDES them on the command line (the "master"): an unnumbered <c>-sbegin 3</c> after the second sequence is <c>-sbegin2 3</c>; a numbered one for another input, any qualifier that is
/// not an associated qualifier of the master, and an abbreviation that the general search would resolve to another item, end that tie; an unnumbered qualifier found without a master (or by a unique prefix that the tie did not accept) applies to ALL inputs.</item>
/// <item>A non-boolean qualifier takes the next token as its value unless that token names a known qualifier ("Value required"); a boolean takes the next token only when it is <c>y yes n no t true f false 0 1</c>; <c>=value</c> and <c>-no</c> + value are refused.</item>
/// <item>Integer values are <c>strtol</c> decimal (sign, digits; saturated to the 32-bit range, a 64-bit overflow is refused), boolean values <c>ajStrToBool</c> (<c>y yes t true 1</c>, <c>n no f false</c>, any number: non-zero = true).</item>
/// </list>
/// </summary>
public sealed class AcdCommandLine
{
    private static readonly Regex QualName = new("^([a-z]+)(_([a-z]+))?([0-9]+)?$", RegexOptions.CultureInvariant);

    private readonly AcdProfile _profile;
    private AcdItem? _masterQual;
    private readonly Dictionary<AcdItem, string> _values = new(ReferenceEqualityComparer.Instance);
    private readonly List<string> _warnings = [];

    /// <summary>Starts a resolution after the sequence parameters of <paramref name="profile"/> were given positionally (the last one is the master of the qualifiers that follow).</summary>
    public AcdCommandLine(AcdProfile profile)
    {
        _profile = profile;
        _masterQual = profile.Inputs.Count > 0 ? profile.Inputs[^1] : null;
    }

    /// <summary>The value the command line gave to an item (the last assignment), or null.</summary>
    public string? Value(AcdItem item) => _values.TryGetValue(item, out var v) ? v : null;

    /// <summary>Every item the command line set, with its value (the last assignment of each).</summary>
    public IEnumerable<KeyValuePair<AcdItem, string>> Assignments => _values;

    private ArgumentException Die(string message)
        => new(_warnings.Count == 0 ? message : string.Join(" ", _warnings.Select(w => "Warning: " + w)) + " " + message);

    private void Warn(string message) => _warnings.Add(message);

    private readonly System.Text.StringBuilder _saved = new("-auto");
    private readonly List<string> _entries = ["-auto"];

    /// <summary>The same arguments one entry per argument (<c>-auto</c>, <c>-sequence a.fa</c> ...): the lines of the command line a report header prints (<c>ajUtilGetCmdline</c>, entries separated by a new line; B09 WP63, F282).</summary>
    public IReadOnlyList<string> SavedEntries => _entries;

    /// <summary>
    /// The arguments of the run as the library keeps them (<c>acdArgSave</c>, <c>ajUtilGetCmdline</c>, tokens separated by a blank here): every qualifier under its full resolved name with its number or its master
    /// (<c>-sbegin1</c>, <c>-osformat_outseq</c> ...), a boolean without its value when the value is the one the name says, a value that is no word (or holds one of <c>* ? [ ] { } | ! &amp; ^</c>) in double quotes.
    /// </summary>
    public string SavedCommandLine => _saved.ToString();

    // acdArgsParse: the part of the loop that appends the token to acdArgSave
    public void Save(AcdItem acd, string qual, string noqual, string value, int number, string master)
    {
        int before = _saved.Length;
        SaveText(acd, qual, noqual, value, number, master);
        _entries.Add(_saved.ToString(before + 1, _saved.Length - before - 1));
    }

    private void SaveText(AcdItem acd, string qual, string noqual, string value, int number, string master)
    {
        _saved.Append(" -").Append(noqual);
        if (number != 0)
            _saved.Append(number.ToString(CultureInfo.InvariantCulture));
        else if (master.Length > 0)
            _saved.Append('_').Append(master);
        bool isWord = value.Length > 0 && !value.Any(char.IsWhiteSpace);
        if (acd.IsBoolean)
        {
            string argval = value;
            if (qual == noqual)
            {
                if (value == "Y")
                    argval = "";
            }
            else if (value == "N")
            {
                argval = "";
            }

            if (argval.Length > 0)
                _saved.Append(' ').Append(argval);
        }
        else if (isWord && value.IndexOfAny(['*', '?', '[', ']', '{', '}', '|', '!', '&', '^']) < 0)
        {
            _saved.Append(' ').Append(value);
        }
        else if (qual == noqual)
        {
            _saved.Append(" \"").Append(value).Append('"');
        }
    }

    /// <summary>Resolves the tokens one after the other (<c>acdArgsParse</c>).</summary>
    /// <exception cref="ArgumentException">The program would die: unknown / ambiguous qualifier, missing value, "no" prefix with a value; or a token that is not a qualifier (the text holds qualifiers only).</exception>
    public void Run(IReadOnlyList<string> tokens)
    {
        int i = 0;
        while (i < tokens.Count)
        {
            string cp = tokens[i];
            string? cq = i + 1 < tokens.Count ? tokens[i + 1] : null;
            int j = IsQual(cp, cq, out int jparam, out string qual, out string noqual, out string value, out int number, out string master, out AcdItem? acd);
            if (j == 0)
                throw Die($"Unexpected argument '{cp}': the text takes qualifiers of the sequence inputs only (for example -sbegin1 3 -sreverse2).");
            // acdFindQualDetail dies (Unknown qualifier) if acd is not set, so acd is set here
            _values[acd!] = value;
            Save(acd!, qual, noqual, value, number, master);
            if (master.Length == 0 && number == 0 && jparam == 0 && acd!.PNum != 0)
            {
                for (int itest = acd.PNum + 1; itest <= _profile.ParamCount; itest++)
                {
                    int dummy = 0;
                    var acd2 = FindQualDetail(qual, null, null, itest, ref dummy);
                    if (acd2 is not null)
                        _values[acd2] = value;
                }
            }

            i += j;
        }
    }

    // ---- name parsing ---------------------------------------------------------------------------------------------------------

    private static bool Prefix(string name, string? prefix) => !string.IsNullOrEmpty(prefix) && name.StartsWith(prefix, StringComparison.Ordinal);

    private static string LowerAscii(string s)
    {
        var sb = new StringBuilder(s.Length);
        foreach (char c in s)
            sb.Append(c is >= 'A' and <= 'Z' ? (char)(c + 32) : c);
        return sb.ToString();
    }

    private static string AmbigApp(string list, string name) => list.Length > 0 ? list + "," + name : name;

    /// <summary><c>acdQualParse</c>: lower-cases and splits <c>name[_master][number]</c>; a name that does not fit clears everything.</summary>
    private void QualParse(ref string qual, out string noqual, out string master, out int number)
    {
        qual = LowerAscii(qual);
        var m = QualName.Match(qual);
        if (!m.Success)
        {
            qual = "";
            noqual = "";
            master = "";
            number = 0;
            return;
        }

        qual = m.Groups[1].Value;
        master = m.Groups[3].Value;
        string num = m.Groups[4].Value;
        noqual = qual.StartsWith("no", StringComparison.Ordinal) ? qual[2..] : "";
        number = 0;
        if (num.Length > 0 && long.TryParse(num, NumberStyles.None, CultureInfo.InvariantCulture, out long n))
            number = (int)Math.Min(n, int.MaxValue);
        if (master.Length > 0)
            FindQual(ref master);
    }

    /// <summary><c>acdFindQual</c>: the full name of a qualifier / parameter from an exact name or a prefix that matches exactly one item.</summary>
    private void FindQual(ref string qual)
    {
        if (qual.Length == 0)
            return;
        AcdItem? ret = null;
        int ifound = 0;
        foreach (var pa in _profile.List)
        {
            if (!Prefix(pa.Name, qual))
                continue;
            if (pa.Name == qual)
                return;
            ifound++;
            ret = pa;
        }

        if (ifound == 1)
            qual = ret!.Name;
    }

    // ---- resolution ------------------------------------------------------------------------------------------------------------

    /// <summary><c>acdFindQualDetail</c>: the item for a name and an optional number (0 = unnumbered: the first match, the caller sets the others).</summary>
    private AcdItem? FindQualDetail(string qual, string? noqual, string? master, int pnum, ref int iparam)
    {
        if (!string.IsNullOrEmpty(master))
        {
            iparam = 0;
            return FindQualMaster(qual, noqual, master, pnum);
        }

        AcdItem? ret = null;
        string ambig = "";
        int ifound = 0;
        bool isparam = false;
        foreach (var pa in _profile.List)
        {
            if (!pa.IsParam)
            {
                if (!Prefix(pa.Name, qual) && !Prefix(pa.Name, noqual))
                    continue;
                bool found = false;
                if (pnum != 0)
                    found = pnum == pa.PNum;
                else if (pa.PNum != 0)
                    found = ifound == 0 || pa.Name != ambig;
                else
                    found = true;
                if (!found)
                    continue;
                if (pa.Name == qual || (!string.IsNullOrEmpty(noqual) && pa.Name == noqual))
                    return pa;
                ambig = AmbigApp(ambig, pa.Name);
                ifound++;
                ret = pa;
            }
            else if (Prefix(pa.Name, qual) || Prefix(pa.Name, noqual))
            {
                if (pa.Name == qual || (!string.IsNullOrEmpty(noqual) && pa.Name == noqual))
                {
                    iparam = pa.PNum;
                    return pa;
                }

                ambig = AmbigApp(ambig, pa.Name);
                ifound++;
                isparam = true;
                ret = pa;
            }
        }

        if (ifound == 1)
        {
            if (isparam)
                iparam = ret!.PNum;
            return ret;
        }

        if (ifound > 1)
            Warn($"ambiguous qualifier '{qual}' ({ambig})");
        return null;
    }

    /// <summary><c>acdFindQualAssoc</c>: the associated qualifier of the master the previous parameter set (a qualifier written for another input number does not match).</summary>
    private AcdItem? FindQualAssoc(AcdItem thys, string qual, string? noqual, string? master, int pnum)
    {
        if (thys.Assoc.Count == 0 || (pnum != 0 && thys.Assoc[0].PNum != pnum))
            return null;
        AcdItem? ret = null;
        string ambig = "";
        int ifound = 0;
        foreach (var pa in thys.Assoc)
        {
            if (!Prefix(pa.Name, qual) && !Prefix(pa.Name, noqual))
                continue;
            if (pa.Name == qual || (!string.IsNullOrEmpty(noqual) && pa.Name == noqual))
                return pa;
            ifound++;
            ret = pa;
            ambig = AmbigApp(ambig, pa.Name);
        }

        if (ifound == 1)
        {
            int iparam = 0;
            var pa = FindQualDetail(qual, noqual, master, pnum, ref iparam);
            if (ReferenceEquals(pa, ret))
                return ret;
        }

        if (ifound > 1)
            Warn($"ambiguous qualifier {qual} ({ambig})");
        return null;
    }

    /// <summary><c>acdFindQualMaster</c>: <c>-qualifier_master</c>, the associated qualifier of a named master (parameter or qualifier).</summary>
    private AcdItem? FindQualMaster(string qual, string? noqual, string master, int pnum)
    {
        AcdItem? ret = null;
        string ambig = "";
        int ifound = 0;
        foreach (var pa in _profile.List)
        {
            if (!pa.IsParam)
            {
                if (!Prefix(pa.Name, master))
                    continue;
                bool found = false;
                if (pnum != 0)
                    found = pnum == pa.PNum;
                else if (pa.PNum != 0)
                    found = ifound == 0 || pa.Name != ambig;
                else
                    found = true;
                if (!found)
                    continue;
                if (pa.Name == master)
                {
                    ret = pa;
                    ifound = 1;
                    break;
                }

                ambig = AmbigApp(ambig, pa.Name);
                ifound++;
                ret = pa;
            }
            else if (Prefix(pa.Name, master))
            {
                if (pa.Name == master)
                {
                    ret = pa;
                    ifound = 1;
                    break;
                }

                ambig = AmbigApp(ambig, pa.Name);
                ifound++;
                ret = pa;
            }
        }

        if (ifound > 1)
        {
            Warn($"Ambiguous master qualifier '{master}' in {qual}_{master} ({ambig})");
            return null;
        }

        if (ifound == 0)
            return null;

        ifound = 0;
        foreach (var pa in ret!.Assoc)
        {
            if (!Prefix(pa.Name, qual) && !Prefix(pa.Name, noqual))
                continue;
            bool found = false;
            if (pnum != 0)
                found = pnum == pa.PNum;
            else if (pa.PNum != 0)
                found = ifound == 0 || pa.Name != ambig;
            else
                found = true;
            if (!found)
                continue;
            if (pa.Name == qual || (!string.IsNullOrEmpty(noqual) && pa.Name == noqual))
                return pa;
            ambig = AmbigApp(ambig, pa.Name);
            ifound++;
            ret = pa;
        }

        if (ifound == 1)
            return ret;
        if (ifound > 1)
            Warn($"ambiguous associated qualifier {qual}_{master} ({ambig})");
        return null;
    }

    /// <summary><c>acdTestAssoc</c>.</summary>
    private static AcdItem? TestAssoc(AcdItem thys, string name, string noname)
    {
        AcdItem? ret = null;
        int ifound = 0;
        foreach (var pa in thys.Assoc)
        {
            if (!Prefix(pa.Name, name) && !Prefix(pa.Name, noname))
                continue;
            if (pa.Name == name || (noname.Length > 0 && pa.Name == noname))
                return pa;
            ifound++;
            ret = pa;
        }

        return ifound == 1 ? ret : null;
    }

    /// <summary><c>acdTestQualC</c>: does the next token name a qualifier (then it is not the value of the previous one)? The count of an abbreviation is over ITEMS (two numbered inputs are two).</summary>
    private bool TestQualC(string name)
    {
        if (name.Length == 0 || (name[0] != '-' && name[0] != '/' && !name.Contains('=')))
            return false;
        string qstr = name[1..];
        int i = qstr.IndexOf('=');
        if (i > 0)
            qstr = qstr[..i];
        QualParse(ref qstr, out string qnostr, out string qmaster, out _);
        if (qmaster.Length > 0)
        {
            foreach (var pa in _profile.List)
            {
                if (pa.Name == qmaster)
                    return TestAssoc(pa, qstr, qnostr) is not null;
            }

            return false;
        }

        int ifound = 0;
        foreach (var pa in _profile.List)
        {
            if (pa.Name == qstr)
                return true;
            if (Prefix(pa.Name, qstr))
                ifound++;
        }

        return ifound == 1;
    }

    /// <summary><c>acdValIsBool</c>: the next token is a boolean value (<c>y yes n no t true f false 0 1</c>, any case).</summary>
    private static bool ValIsBool(string? arg)
    {
        if (arg is null || arg.Length == 0)
            return false;
        return LowerAscii(arg) is "y" or "yes" or "n" or "no" or "t" or "true" or "f" or "false" or "0" or "1";
    }

    /// <summary><c>acdIsQual</c>: is the argument a qualifier, which item, which value; the number of arguments consumed (0 = a parameter, 1, 2).</summary>
    private int IsQual(string arg, string? arg2, out int iparam, out string pqual, out string pnoqual, out string pvalue, out int number, out string pmaster, out AcdItem? acd)
    {
        iparam = 0;
        pqual = "";
        pnoqual = "";
        pvalue = "";
        pmaster = "";
        number = 0;
        acd = null;
        if (arg == "-" || arg.Length == 0 || arg == "--")
            return 0;

        int pos = 0;
        bool qstart = false;
        if (arg[0] is '-' or '/')
        {
            pos++;
            qstart = true;
            if (pos < arg.Length && arg[pos] == '-')
                pos++;
        }

        if (pos >= arg.Length)
            return 0;

        string qual = arg[pos..];
        int ret = 1;
        bool gotvalue = false;
        int eq = qual.IndexOf('=');
        if (eq >= 0)
        {
            pvalue = eq + 1 == qual.Length ? "" : qual[(eq + 1)..];
            qual = qual[..eq];
            gotvalue = true;
        }
        else
        {
            if (!qstart)
                return 0;
            if (qual.Any(c => c != '_' && !(c is >= '0' and <= '9' || c is >= 'a' and <= 'z' || c is >= 'A' and <= 'Z')))
                return 0;
        }

        QualParse(ref qual, out string noqual, out string master, out number);
        pqual = qual;
        pmaster = master;
        if (master.Length > 0)
            _masterQual = null;

        if (_masterQual is not null)
        {
            acd = FindQualAssoc(_masterQual, qual, noqual, master, number);
            if (acd is null)
                _masterQual = null;
            else
                number = _masterQual.PNum;
        }

        if (_masterQual is null)
            acd = FindQualDetail(qual, noqual, master, number, ref iparam);

        if (acd is null)
            throw Die($"Unknown qualifier {arg}");

        pqual = acd.Name;
        pnoqual = acd.Name;
        if (acd.Assoc.Count > 0)
            _masterQual = acd;

        if (gotvalue)
        {
            if (Prefix(acd.Name, noqual))
                throw Die($"'no' prefix used with value for '{arg}'");
            return ret;
        }

        if (Prefix(acd.Name, noqual))
        {
            pnoqual = "no" + acd.Name;
            if (acd.IsBoolean)
            {
                pvalue = "N";
                return 1;
            }

            throw Die($"'no' prefix invalid for '{arg}'");
        }

        if (acd.IsBoolean)
        {
            if (ValIsBool(arg2))
            {
                ret = 2;
                pvalue = arg2!;
            }
            else
            {
                pvalue = "Y";
            }

            return ret;
        }

        if (arg2 is null)
            throw Die($"Value required for '{arg}'");
        if (arg2.Length > 0 && arg2[0] == '-' && TestQualC(arg2))
            throw Die($"Value required for '{arg}' before '{arg2}'");

        pvalue = arg2;
        return 2;
    }

    // ---- values ----------------------------------------------------------------------------------------------------------------

    /// <summary><c>ajStrToInt</c>: <c>strtol(.., 10)</c> over the whole string, saturated to the 32-bit range; a 64-bit overflow or trailing text is invalid.</summary>
    public static bool TryToInt(string text, out int value)
    {
        value = 0;
        if (text.Length == 0)
            return false;
        int i = 0;
        while (i < text.Length && (text[i] is ' ' or '\t' or '\n' or '\v' or '\f' or '\r'))
            i++;
        bool negative = false;
        if (i < text.Length && text[i] is '+' or '-')
        {
            negative = text[i] == '-';
            i++;
        }

        int digits = 0;
        System.Numerics.BigInteger n = 0;
        for (; i < text.Length && text[i] is >= '0' and <= '9'; i++, digits++)
            n = n * 10 + (text[i] - '0');
        if (digits == 0 || i != text.Length)
            return false;
        if (negative)
            n = -n;
        if (n > long.MaxValue || n < long.MinValue)
            return false;
        value = (int)Math.Clamp((long)n, int.MinValue, int.MaxValue);
        return true;
    }

    /// <summary><c>ajStrToBool</c> (the value of a boolean qualifier): <c>y yes t true 1</c> true, <c>n no f false</c> false, a number (<c>strtod</c>, at most FLT_MAX) true unless it is zero; anything else is invalid.</summary>
    public static bool TryToBool(string text, out bool value)
    {
        value = false;
        if (text.Length == 0)
            return false;
        char c = text[0];
        if ("YyTt1".Contains(c))
        {
            value = true;
            if (text.Length == 1)
                return true;
            return string.Equals(text, "yes", StringComparison.OrdinalIgnoreCase) || string.Equals(text, "true", StringComparison.OrdinalIgnoreCase);
        }

        if ("NnFf".Contains(c))
        {
            if (text.Length == 1)
                return true;
            return string.Equals(text, "no", StringComparison.OrdinalIgnoreCase) || string.Equals(text, "false", StringComparison.OrdinalIgnoreCase);
        }

        if ("123456789".Contains(c))
        {
            value = true;
            return IsFloat(text);
        }

        if ("0+-".Contains(c))
        {
            value = text.Any(ch => ch is >= '1' and <= '9');
            return IsFloat(text);
        }

        return false;
    }

    /// <summary><c>ajStrIsFloat</c>: the whole string is a <c>strtod</c> number (decimal, hexadecimal, inf / nan) within the single-precision range, without ERANGE.</summary>
    private static bool IsFloat(string text)
    {
        if (text.Length == 0)
            return false;
        int i = 0;
        while (i < text.Length && (text[i] is ' ' or '\t' or '\n' or '\v' or '\f' or '\r'))
            i++;
        bool negative = false;
        if (i < text.Length && text[i] is '+' or '-')
        {
            negative = text[i] == '-';
            i++;
        }

        string rest = text[i..];
        string lower = LowerAscii(rest);
        if (lower is "inf" or "infinity")
            return false;                                      // |d| > FLT_MAX
        if (lower == "nan" || (lower.StartsWith("nan(", StringComparison.Ordinal) && lower.EndsWith(')') && lower[4..^1].All(ch => ch == '_' || char.IsAsciiLetterOrDigit(ch))))
            return true;                                       // NaN fails neither comparison
        double d;
        if (lower.StartsWith("0x", StringComparison.Ordinal))
        {
            if (!TryHexFloat(lower[2..], out d))
                return false;
        }
        else
        {
            var m = Regex.Match(rest, @"^(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$", RegexOptions.CultureInvariant);
            if (!m.Success || !double.TryParse(rest, NumberStyles.Float, CultureInfo.InvariantCulture, out d))
                return false;
            if (double.IsInfinity(d))
                return false;                                  // ERANGE
            if (d == 0 && rest.Any(ch => ch is >= '1' and <= '9') && m.Groups[1].Value.Any(ch => ch is >= '1' and <= '9'))
                return false;                                  // underflow to zero: ERANGE
        }

        if (negative)
            d = -d;
        return d <= float.MaxValue && d >= -float.MaxValue;
    }

    private static bool TryHexFloat(string text, out double value)
    {
        value = 0;
        var m = Regex.Match(text, @"^([0-9a-f]*)(?:\.([0-9a-f]*))?(?:p([+-]?\d+))?$", RegexOptions.CultureInvariant);
        if (!m.Success || (m.Groups[1].Length == 0 && m.Groups[2].Length == 0))
            return false;
        double mantissa = 0;
        foreach (char ch in m.Groups[1].Value)
            mantissa = mantissa * 16 + Convert.ToInt32(ch.ToString(), 16);
        double scale = 1.0 / 16;
        foreach (char ch in m.Groups[2].Value)
        {
            mantissa += Convert.ToInt32(ch.ToString(), 16) * scale;
            scale /= 16;
        }

        int exp = m.Groups[3].Success && int.TryParse(m.Groups[3].Value, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out int e) ? e : 0;
        value = mantissa * Math.Pow(2, exp);
        return !double.IsInfinity(value);
    }
}
