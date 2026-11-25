// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.Motifs, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127, S1905, S2681, S1121
using System.Globalization;

namespace Seqeron.Genomics.Emboss.Motifs;

using Seqeron.Genomics.Analysis;

/// <summary>The regular expression text of the library of EMBOSS 6.6.0 (<c>ajRegComp</c> = <c>pcre_compile</c> of the bundled PCRE 7, no options) is not accepted by PCRE.</summary>
public sealed class EmbossPcreException(string message, int position) : Exception(message)
{
    /// <summary>The offset in the expression where PCRE found the error (<c>errptr</c> / <c>erroffset</c>).</summary>
    public int Position { get; } = position;
}

/// <summary>
/// The regular expression engine behind <c>ajRegComp</c> / <c>ajRegExecC</c> of EMBOSS 6.6.0 as the programs that search PROSITE patterns use it (<c>embPatMatchFind</c>: <c>patmatdb</c>, <c>patmatmotifs</c>): PCRE 7 (bundled
/// <c>libepcre</c>) with no options, Perl-compatible matching (leftmost match, greedy / lazy / possessive quantifiers, backtracking) on 8-bit characters. Supported syntax: literals, <c>.</c>, character classes (ranges, POSIX classes,
/// escapes), escapes (<c>\d \D \s \S \w \W \h \H \v \V \R \b \B \A \z \Z \G \K \C \n \r \t \f \e \a \0nn \xhh \x{h} \cX \Q..\E</c>, back references), groups <c>( )</c>, <c>(?: )</c>, <c>(?= )</c>, <c>(?! )</c>, <c>(?&gt; )</c>,
/// <c>(?# )</c>, alternation, the quantifiers <c>* + ? {n} {n,} {n,m}</c> with lazy and possessive suffixes, <c>^</c> and <c>$</c>. The errors of the compiler are the messages and positions of PCRE 7 (checked against the library on the
/// expressions that the PROSITE conversion can produce); the constructs that PROSITE text cannot produce (look-behind, option settings, named groups, recursion, conditionals, callouts) are not ported.
/// </summary>
internal sealed class EmbossPcre
{
    private interface Node
    {
        bool M(Matcher m, int i, Func<int, bool> k);
    }

    private sealed class Matcher(string s)
    {
        public readonly string S = s;
        public readonly int[] Caps = new int[2 * 100];
        public int StartOfMatch;
    }

    private sealed class CharNode(char c) : Node
    {
        public bool M(Matcher m, int i, Func<int, bool> k) => i < m.S.Length && m.S[i] == c && k(i + 1);
    }

    private sealed class SetNode(bool[] set) : Node
    {
        public bool M(Matcher m, int i, Func<int, bool> k) => i < m.S.Length && m.S[i] < 256 && set[m.S[i]] && k(i + 1);
    }

    private sealed class AnyNode(bool all) : Node
    {
        public bool M(Matcher m, int i, Func<int, bool> k) => i < m.S.Length && (all || m.S[i] != '\n') && k(i + 1);
    }

    private sealed class SeqNode(List<Node> items) : Node
    {
        public List<Node> Items { get; } = items;

        public bool M(Matcher m, int i, Func<int, bool> k) => Step(m, 0, i, k);

        private bool Step(Matcher m, int idx, int i, Func<int, bool> k) => idx == Items.Count ? k(i) : Items[idx].M(m, i, j => Step(m, idx + 1, j, k));
    }

    private sealed class AltNode(List<Node> alts) : Node
    {
        public bool M(Matcher m, int i, Func<int, bool> k)
        {
            foreach (var a in alts)
                if (a.M(m, i, k))
                    return true;
            return false;
        }
    }

    private sealed class GroupNode(Node inner, int index) : Node
    {
        public bool M(Matcher m, int i, Func<int, bool> k)
        {
            if (index < 0)
                return inner.M(m, i, k);
            int os = m.Caps[2 * index], oe = m.Caps[2 * index + 1];
            bool ok = inner.M(m, i, j =>
            {
                int ps = m.Caps[2 * index], pe = m.Caps[2 * index + 1];
                m.Caps[2 * index] = i;
                m.Caps[2 * index + 1] = j;
                if (k(j))
                    return true;
                m.Caps[2 * index] = ps;
                m.Caps[2 * index + 1] = pe;
                return false;
            });
            if (!ok)
            {
                m.Caps[2 * index] = os;
                m.Caps[2 * index + 1] = oe;
            }

            return ok;
        }
    }

    private sealed class AtomicNode(Node inner) : Node
    {
        public bool M(Matcher m, int i, Func<int, bool> k)
        {
            int end = -1;
            if (!inner.M(m, i, j => { end = j; return true; }))
                return false;
            return k(end);
        }
    }

    private sealed class LookNode(Node inner, bool negative) : Node
    {
        public bool M(Matcher m, int i, Func<int, bool> k)
        {
            bool found = inner.M(m, i, _ => true);
            return found != negative && k(i);
        }
    }

    private sealed class RepeatNode(Node inner, int min, int max, int mode) : Node
    {
        // mode 0 greedy, 1 lazy, 2 possessive
        public bool M(Matcher m, int i, Func<int, bool> k)
        {
            if (mode == 2)
            {
                var atom = new RepeatNode(inner, min, max, 0);
                return new AtomicNode(atom).M(m, i, k);
            }

            return mode == 0 ? Greedy(m, i, 0, k) : Lazy(m, i, 0, k);
        }

        private bool Greedy(Matcher m, int i, int count, Func<int, bool> k)
        {
            if (count < max && (inner.M(m, i, j => (j != i || count < min) && Greedy(m, j, count + 1, k))))
                return true;
            return count >= min && k(i);
        }

        private bool Lazy(Matcher m, int i, int count, Func<int, bool> k)
        {
            if (count >= min && k(i))
                return true;
            return count < max && inner.M(m, i, j => (j != i || count < min) && Lazy(m, j, count + 1, k));
        }
    }

    private sealed class AssertNode(char kind) : Node
    {
        public bool M(Matcher m, int i, Func<int, bool> k)
        {
            var s = m.S;
            bool ok;
            switch (kind)
            {
                case '^':
                case 'A':
                    ok = i == 0;
                    break;
                case '$':
                case 'Z':
                    ok = i == s.Length || (i == s.Length - 1 && s[i] == '\n');
                    break;
                case 'z':
                    ok = i == s.Length;
                    break;
                case 'G':
                    ok = i == 0;
                    break;
                case 'b':
                case 'B':
                    {
                        bool a = i > 0 && IsWord(s[i - 1]), b = i < s.Length && IsWord(s[i]);
                        ok = (a != b) == (kind == 'b');
                        break;
                    }

                default:
                    ok = true;
                    break;
            }

            return ok && k(i);
        }
    }

    private sealed class KeepNode : Node
    {
        public bool M(Matcher m, int i, Func<int, bool> k)
        {
            int old = m.StartOfMatch;
            m.StartOfMatch = i;
            if (k(i))
                return true;
            m.StartOfMatch = old;
            return false;
        }
    }

    private sealed class BackrefNode(int n) : Node
    {
        public bool M(Matcher m, int i, Func<int, bool> k)
        {
            int s = m.Caps[2 * n], e = m.Caps[2 * n + 1];
            if (e < 0 || s < 0)
                return false;                                                     // an unset group makes the reference fail
            int len = e - s;
            if (i + len > m.S.Length || string.CompareOrdinal(m.S, s, m.S, i, len) != 0)
                return false;
            return k(i + len);
        }
    }

    private sealed class StubNode : Node
    {
        public bool M(Matcher m, int i, Func<int, bool> k) => false;
    }

    private sealed class CrLfNode : Node
    {
        public bool M(Matcher m, int i, Func<int, bool> k)
        {
            var s = m.S;
            if (i >= s.Length)
                return false;
            if (s[i] == '\r' && i + 1 < s.Length && s[i + 1] == '\n' && k(i + 2))
                return true;
            char c = s[i];
            return (c is '\n' or '\u000b' or '\u000c' or '\r' or '\u0085') && k(i + 1);
        }
    }

    private static bool IsWord(char c) => c < 128 && (char.IsAsciiLetterOrDigit(c) || c == '_');

    private readonly Node _root;

    private EmbossPcre(Node root)
    {
        _root = root;
    }

    /// <summary>The result of a search: the offset of the match in the subject and its length (<c>ajRegOffset</c>, <c>ajRegLenI(0)</c>).</summary>
    internal readonly record struct Match(int Offset, int Length);

    /// <summary><c>pcre_exec</c> of the subject <paramref name="s"/> from its first character: the first match found in the order of Perl (a start position tried after another), or null.</summary>
    internal Match? Exec(string s)
    {
        var m = new Matcher(s);
        for (int start = 0; start <= s.Length; start++)
        {
            Array.Fill(m.Caps, -1);
            m.StartOfMatch = start;
            int end = -1;
            if (_root.M(m, start, j => { end = j; return true; }))
                return new Match(m.StartOfMatch, end - m.StartOfMatch);
        }

        return null;
    }

    // ---------------------------------------------------------------------------------------------------------------- the compiler

    /// <summary><c>pcre_compile</c> with no options.</summary>
    /// <exception cref="EmbossPcreException">The expression is refused (the message and the offset of PCRE).</exception>
    internal static EmbossPcre Compile(string pattern)
    {
        var p = new PatternReader(pattern);
        var root = p.ParseAlternation(0);
        p.Finish();                                                              // the references that name nothing (found while the second pass of PCRE compiles) come before the check of the end
        if (p.Pos < pattern.Length)
            throw new EmbossPcreException("unmatched parentheses", p.Pos);
        p.FinishNumbered();
        return new EmbossPcre(root);
    }

    private sealed class PatternReader(string p)
    {
        public int Pos;
        public int Groups;
        private int _maxBackref;
        private readonly string _p = p;

        public void Finish()
        {
            // the first reference that names nothing, by position
            int firstPos = int.MaxValue;
            bool unsupported = false;
            if (_duplicateNamePos is { } dup)
                throw new EmbossPcreException("two named subpatterns have the same name", dup);
            foreach (var (name, pos, subroutine) in _namedRefs)
            {
                if (!_names.ContainsKey(name))
                    firstPos = Math.Min(firstPos, pos);
                else if (subroutine)
                    unsupported = true;
            }

            foreach (var (number, pos) in _subroutineRefs)
            {
                if (number > Groups)
                    firstPos = Math.Min(firstPos, pos);
                else
                    unsupported = true;
            }

            if (firstPos != int.MaxValue)
                throw new EmbossPcreException("reference to non-existent subpattern", firstPos);
            if (unsupported)
                throw new NotSupportedException("A subroutine call of PCRE is not ported.");
        }

        public void FinishNumbered()
        {
            if (_maxBackref > Groups)
                throw new EmbossPcreException("reference to non-existent subpattern", _p.Length);
        }

        private bool More => Pos < _p.Length;

        // alternation: items | items ... ; stops at ')' (depth > 0) or at the end
        public Node ParseAlternation(int depth)
        {
            var alts = new List<Node>();
            alts.Add(ParseSequence(depth));
            while (More && _p[Pos] == '|')
            {
                Pos++;
                alts.Add(ParseSequence(depth));
            }

            return alts.Count == 1 ? alts[0] : new AltNode(alts);
        }

        private Node ParseSequence(int depth)
        {
            var items = new List<Node>();
            Node? previous = null;                                              // the item a quantifier applies to (null: nothing to repeat)
            while (More)
            {
                char c = _p[Pos];
                if (c == '|')
                    break;
                if (c == ')')
                    break;                                                      // closes the group; at the top level it is reported after the references are checked

                Pos++;
                switch (c)
                {
                    case '^':
                        items.Add(new AssertNode('^'));
                        previous = null;
                        break;
                    case '$':
                        items.Add(new AssertNode('$'));
                        previous = null;
                        break;
                    case '.':
                        items.Add(previous = new AnyNode(false));
                        break;
                    case '[':
                        if (Pos < _p.Length && _p[Pos] is ':' or '.' or '=' && CheckPosixSyntax(Pos - 1, out _))
                            throw new EmbossPcreException(_p[Pos] == ':' ? "POSIX named classes are supported only within a class" : "POSIX collating elements are not supported", Pos - 1);
                        items.Add(previous = ParseClass());
                        break;
                    case '*':
                    case '+':
                    case '?':
                        Quantify(items, ref previous, c == '*' ? 0 : c == '+' ? 1 : 0, c == '?' ? 1 : -1, c);
                        break;
                    case '{':
                        {
                            if (IsCountedRepeat(Pos))
                            {
                                var (min, max) = ReadRepeat();
                                Quantify(items, ref previous, min, max, '{');
                            }
                            else
                            {
                                items.Add(previous = new CharNode('{'));
                            }

                            break;
                        }

                    case '(':
                        {
                            var g = ParseGroup(depth, out bool skip);
                            if (!skip)
                                items.Add(previous = g);
                            break;
                        }

                    case '\\':
                        {
                            var node = ParseEscape();
                            if (node is null)
                                break;                                          // \E: nothing
                            if (node is AssertNode || node is KeepNode)
                            {
                                items.Add(node);
                                previous = null;
                            }
                            else if (node is SeqNode quoted)
                            {
                                items.AddRange(quoted.Items);                    // \Q...\E: the literal characters; a quantifier applies to the last
                                previous = quoted.Items.Count > 0 ? quoted.Items[^1] : previous;
                            }
                            else
                            {
                                items.Add(previous = node);
                            }

                            break;
                        }

                    default:
                        items.Add(previous = new CharNode(c));
                        break;
                }
            }

            return items.Count == 1 ? items[0] : new SeqNode(items);
        }

        // the quantifier applies to the last item
        private void Quantify(List<Node> items, ref Node? previous, int min, int max, char kind)
        {
            if (previous is null || items.Count == 0 || !ReferenceEquals(items[^1], previous))
                throw new EmbossPcreException("nothing to repeat", Pos - 1);
            if (kind == '?' || kind == '*' || kind == '+')
            {
                // * -> {0,inf}, + -> {1,inf}, ? -> {0,1}
                if (kind == '*') { min = 0; max = int.MaxValue; }
                else if (kind == '+') { min = 1; max = int.MaxValue; }
                else { min = 0; max = 1; }
            }
            else if (max < 0)
            {
                max = int.MaxValue;
            }

            int mode = 0;
            if (More && _p[Pos] == '+')
            {
                mode = 2;
                Pos++;
            }
            else if (More && _p[Pos] == '?')
            {
                mode = 1;
                Pos++;
            }

            items[^1] = new RepeatNode(previous, min, max, mode);
            previous = null;                                                     // END_REPEAT: another quantifier has nothing to repeat
        }

        private bool IsCountedRepeat(int p)
        {
            int q = p;
            if (q >= _p.Length || !char.IsAsciiDigit(_p[q]))
                return false;
            while (q < _p.Length && char.IsAsciiDigit(_p[q]))
                q++;
            if (q < _p.Length && _p[q] == '}')
                return true;
            if (q >= _p.Length || _p[q] != ',')
                return false;
            q++;
            while (q < _p.Length && char.IsAsciiDigit(_p[q]))
                q++;
            return q < _p.Length && _p[q] == '}';
        }

        private (int Min, int Max) ReadRepeat()
        {
            int min = ReadNumber();
            int max = min;
            if (_p[Pos] == ',')
            {
                Pos++;
                max = _p[Pos] == '}' ? -1 : ReadNumber();
            }

            if (_p[Pos] != '}')
                throw new EmbossPcreException("internal error", Pos);
            Pos++;
            if (min > 65535 || max > 65535)
                throw new EmbossPcreException("number too big in {} quantifier", Pos - 1);
            if (max >= 0 && max < min)
                throw new EmbossPcreException("numbers out of order in {} quantifier", Pos - 1);
            return (min, max);
        }

        private int ReadNumber()
        {
            long v = 0;
            while (Pos < _p.Length && char.IsAsciiDigit(_p[Pos]))
            {
                v = Math.Min(v * 10 + (_p[Pos] - '0'), 1_000_000);
                Pos++;
            }

            return (int)v;
        }

        private readonly Dictionary<string, int> _names = new(StringComparer.Ordinal);
        private int? _duplicateNamePos;
        private readonly List<(string Name, int Pos, bool Subroutine)> _namedRefs = [];
        private readonly List<(int Number, int Pos)> _subroutineRefs = [];

        // the name of a group: word characters up to the terminator; Pos at the first character of the name
        private string ReadName(char terminator, string emptyMessage)
        {
            int start = Pos;
            while (Pos < _p.Length && (char.IsAsciiLetterOrDigit(_p[Pos]) || _p[Pos] == '_'))
                Pos++;
            if (Pos == start)
                throw new EmbossPcreException(emptyMessage, Pos);
            if (Pos >= _p.Length || _p[Pos] != terminator)
                throw new EmbossPcreException("syntax error in subpattern name (missing terminator)", Pos);
            if (Pos - start > 32)
                throw new EmbossPcreException("subpattern name is too long (maximum 32 characters)", Pos);
            string name = _p[start..Pos];
            Pos++;
            return name;
        }

        private int NewGroup(string? name)
        {
            int index = ++Groups;
            if (Groups >= 99)
                throw new NotSupportedException("More than 98 groups are not ported.");
            if (name is not null && !_names.TryAdd(name, index))
                _duplicateNamePos ??= Pos;                                        // found by the second pass of PCRE: an error of the first pass comes first

            return index;
        }

        private Node ParseGroup(int depth, out bool skip)
        {
            skip = false;
            int index = -1;
            int kind = 0;                                                       // 0 capture, 1 non-capture, 2 atomic, 3 lookahead, 4 negative lookahead
            bool unsupported = false;
            if (More && _p[Pos] == '?')
            {
                Pos++;
                if (!More)
                    throw new EmbossPcreException("unrecognized character after (? or (?-", Pos);
                char t = _p[Pos];
                switch (t)
                {
                    case '#':
                        {
                            int close = _p.IndexOf(')', Pos);
                            if (close < 0)
                                throw new EmbossPcreException("missing ) after comment", _p.Length);
                            Pos = close + 1;
                            skip = true;
                            return new SeqNode([]);
                        }

                    case ':':
                        Pos++;
                        kind = 1;
                        break;
                    case '|':
                        Pos++;
                        kind = 1;
                        unsupported = true;                                     // a branch reset group
                        break;
                    case '>':
                        Pos++;
                        kind = 2;
                        break;
                    case '=':
                        Pos++;
                        kind = 3;
                        break;
                    case '!':
                        Pos++;
                        kind = 4;
                        break;
                    case '<':
                        Pos++;
                        if (!More)
                            throw new EmbossPcreException("unrecognized character after (?<", Pos);
                        if (_p[Pos] is '=' or '!')
                        {
                            Pos++;
                            kind = 1;
                            unsupported = true;                                 // a look-behind assertion
                            break;
                        }

                        index = NewGroupNamed(ReadName('>', "unrecognized character after (?<"));
                        break;
                    case '\'':
                        Pos++;
                        index = NewGroupNamed(ReadName('\'', "unrecognized character after (?'"));
                        break;
                    case 'P':
                        Pos++;
                        if (!More || _p[Pos] is not ('<' or '=' or '>'))
                            throw new EmbossPcreException("unrecognized character after (?P", Pos);
                        if (_p[Pos] == '<')
                        {
                            Pos++;
                            index = NewGroupNamed(ReadName('>', "syntax error in subpattern name (missing terminator)"));
                            break;
                        }

                        {
                            bool subroutine = _p[Pos] == '>';
                            Pos++;
                            string refName = ReadName(')', "syntax error in subpattern name (missing terminator)");
                            _namedRefs.Add((refName, Pos - 1, subroutine));
                            skip = false;
                            return subroutine ? new StubNode() : new NamedBackrefNode(this, refName);
                        }

                    case 'R':
                    case 'C':
                    case '(':
                    case '+':
                    case >= '0' and <= '9':
                        throw new NotSupportedException($"The group form '(?{t}' of PCRE (recursion, callout, condition) is not ported.");
                    case '-' when Pos + 1 < _p.Length && char.IsAsciiDigit(_p[Pos + 1]):
                        throw new NotSupportedException("The relative recursion '(?-n)' of PCRE is not ported.");
                    default:
                        {
                            int q = Pos;
                            while (q < _p.Length && "imsxUXJ-".Contains(_p[q], StringComparison.Ordinal))
                                q++;
                            if (q >= _p.Length || (_p[q] != ')' && _p[q] != ':'))
                                throw new EmbossPcreException("unrecognized character after (? or (?-", q);
                            if (_p[q] == ')')
                            {
                                if (q == Pos)
                                {
                                    Pos = q + 1;                                // (?) sets nothing
                                    skip = true;
                                    return new SeqNode([]);
                                }

                                throw new NotSupportedException("The option setting '(?" + _p[Pos..q] + ")' of PCRE is not ported.");
                            }

                            Pos = q + 1;
                            kind = 1;
                            unsupported = true;                                 // a group with options
                            break;
                        }
                }
            }
            else
            {
                index = NewGroup(null);
            }

            var inner = ParseAlternation(depth + 1);
            if (!More || _p[Pos] != ')')
                throw new EmbossPcreException("missing )", _p.Length);
            Pos++;
            if (unsupported)
                throw new NotSupportedException("A look-behind assertion, branch reset or group with options of PCRE is not ported.");
            return kind switch
            {
                2 => new AtomicNode(inner),
                3 => new LookNode(inner, false),
                4 => new LookNode(inner, true),
                _ => new GroupNode(inner, index),
            };
        }

        private int NewGroupNamed(string name) => NewGroup(name);

        private sealed class NamedBackrefNode(PatternReader parser, string name) : Node
        {
            public bool M(Matcher m, int i, Func<int, bool> k) => new BackrefNode(parser._names[name]).M(m, i, k);
        }

        // ------------------------------------------------------------------------------------------------------------ escapes

        // returns the node of the escape at Pos (after the backslash); null for \E
        private Node? ParseEscape()
        {
            int start = Pos;
            if (!More)
                throw new EmbossPcreException("\\ at end of pattern", Pos);
            char c = _p[Pos++];
            switch (c)
            {
                case 'A':
                case 'B':
                case 'b':
                case 'G':
                case 'Z':
                case 'z':
                    return new AssertNode(c);
                case 'K':
                    return new KeepNode();
                case 'E':
                    return null;
                case 'Q':
                    {
                        var lits = new List<Node>();
                        while (More)
                        {
                            if (_p[Pos] == '\\' && Pos + 1 < _p.Length && _p[Pos + 1] == 'E')
                            {
                                Pos += 2;
                                break;
                            }

                            lits.Add(new CharNode(_p[Pos++]));
                        }

                        return new SeqNode(lits);
                    }

                case 'C':
                    return new AnyNode(true);
                case 'R':
                    return new CrLfNode();
                case 'X':
                case 'p':
                case 'P':
                    throw new EmbossPcreException("support for \\P, \\p, and \\X has not been compiled", Pos - 1);
                case 'L':
                case 'l':
                case 'N':
                case 'U':
                case 'u':
                    throw new EmbossPcreException("PCRE does not support \\L, \\l, \\N, \\U, or \\u", Pos - 1);
                case 'g':
                    return ParseGEscape();
                case 'k':
                    return ParseKEscape();
                case >= '1' and <= '9':
                    {
                        // a back reference \n when the number has one or two digits and is below 10 or a group that exists... decided at the end (PCRE 7: \1 - \9 always; larger numbers when that many groups exist before)
                        int save = Pos - 1;
                        int n = c - '0';
                        int q = Pos;
                        while (q < _p.Length && char.IsAsciiDigit(_p[q]) && n < 100000)
                        {
                            n = n * 10 + (_p[q] - '0');
                            q++;
                        }

                        if (n < 10 || n <= Groups)
                        {
                            Pos = q;
                            if (n > _maxBackref)
                                _maxBackref = n;

                            return new BackrefNode(n);
                        }

                        Pos = save;
                        return Literal(ReadOctalOrLiteralDigit());
                    }

                default:
                    break;
            }

            Pos = start;
            var cls = new bool[256];
            if (SimpleEscapeClass(_p[Pos], cls))
            {
                Pos++;
                return new SetNode(cls);
            }

            return Literal(ReadCharEscape(false));
        }

        private static Node Literal(int c) => new CharNode((char)c);

        // \g: \gn, \g{n}, \g{-n} back references; \g<n>, \g'n' subroutine calls (Pos after the g)
        private Node ParseGEscape()
        {
            int gpos = Pos - 1;
            if (!More)
                throw new EmbossPcreException("a numbered reference must not be zero", gpos);
            char c = _p[Pos];
            bool subroutine = false;
            string? name = null;
            int number = 0;
            bool relative = false;
            if (c is '{' or '<' or '\'')
            {
                char close = c == '{' ? '}' : c == '<' ? '>' : '\'';
                subroutine = c != '{';
                Pos++;
                if (More && _p[Pos] == close)
                    throw new EmbossPcreException("a numbered reference must not be zero", Pos);
                if (More && (char.IsAsciiDigit(_p[Pos]) || _p[Pos] is '-' or '+'))
                {
                    relative = _p[Pos] is '-' or '+';
                    bool plus = _p[Pos] == '+';
                    if (relative)
                        Pos++;
                    int q = Pos;
                    while (q < _p.Length && char.IsAsciiDigit(_p[q]))
                        q++;
                    if (q == Pos || q >= _p.Length || _p[q] != close)
                        throw new EmbossPcreException("syntax error in subpattern name (missing terminator)", Math.Min(q, _p.Length));
                    number = int.Parse(_p[Pos..q], System.Globalization.CultureInfo.InvariantCulture);
                    Pos = q + 1;
                    if (relative)
                        number = plus ? Groups + number : Groups - number + 1;
                }
                else
                {
                    name = ReadName(close, "subpattern name expected");
                }
            }
            else if (c == '-' || char.IsAsciiDigit(c))
            {
                relative = c == '-';
                int q = relative ? Pos + 1 : Pos;
                int start = q;
                while (q < _p.Length && char.IsAsciiDigit(_p[q]))
                    q++;
                if (q == start)
                    throw new EmbossPcreException("a numbered reference must not be zero", q - 1);
                number = int.Parse(_p[start..q], System.Globalization.CultureInfo.InvariantCulture);
                Pos = q;
                if (relative)
                    number = Groups - number + 1;
            }
            else
            {
                throw new EmbossPcreException("a numbered reference must not be zero", gpos);
            }

            if (name is null)
            {
                if (relative && number < 1)
                    throw new EmbossPcreException("reference to non-existent subpattern", Pos - 1);
                if (number == 0)
                    throw new EmbossPcreException("a numbered reference must not be zero", Pos - 1);
                if (subroutine)
                {
                    _subroutineRefs.Add((number, Pos - 1));
                    return new StubNode();
                }

                if (number > _maxBackref)
                    _maxBackref = number;
                return new BackrefNode(number);
            }

            _namedRefs.Add((name, Pos - 1, subroutine));
            return subroutine ? new StubNode() : new NamedBackrefNode(this, name);
        }

        // \g in a class (Pos after the g): the braced / angle / quoted forms are consumed and add nothing, a plain number is left to be read as characters; only the errors of the syntax remain
        private void ClassGEscape()
        {
            int gpos = Pos - 1;
            if (More && _p[Pos] is '{' or '<' or '\'')
            {
                char close = _p[Pos] == '{' ? '}' : _p[Pos] == '<' ? '>' : '\'';
                int end = _p.IndexOf(close, Pos + 1);
                if (end < 0)
                    throw new EmbossPcreException("\\g is not followed by a braced, angle-bracketed, or quoted name/number or by a plain number", _p.Length);
                string text = _p[(Pos + 1)..end];
                Pos = end + 1;
                if (text.Length > 0 && text.All(ch => char.IsAsciiDigit(ch) || ch is '-' or '+'))
                {
                    int n = int.Parse(text.TrimStart('-', '+'), System.Globalization.CultureInfo.InvariantCulture);
                    if (text[0] == '-' && Groups - n + 1 < 1)
                        throw new EmbossPcreException("reference to non-existent subpattern", Pos - 1);
                    if (n == 0)
                        throw new EmbossPcreException("a numbered reference must not be zero", Pos - 1);
                }

                throw new NotSupportedException("A back reference or subroutine call inside a class is not ported.");
            }

            if (!More || !(char.IsAsciiDigit(_p[Pos]) && _p[Pos] != '0'))
                throw new EmbossPcreException("a numbered reference must not be zero", gpos);
            throw new NotSupportedException("A back reference inside a class is not ported (PCRE 7 reads \\g1 there as the digit).");
        }

        // \k<name> \k{name} \k'name': a back reference by name; any other \k is the letter of the library's table (PCRE 7 reads it as U)
        private Node ParseKEscape()
        {
            if (More && _p[Pos] is '<' or '{' or '\'')
            {
                char open = _p[Pos];
                char close = open == '{' ? '}' : open == '<' ? '>' : '\'';
                Pos++;
                string name = ReadName(close, "subpattern name expected");
                _namedRefs.Add((name, Pos - 1, false));
                return new NamedBackrefNode(this, name);
            }

            throw new NotSupportedException("A bare '\\k' of PCRE 7 is not ported.");
        }

        // \0nn, \1-\7 octal (up to 3 digits), \8 \9 literal digits: Pos at the first digit
        private int ReadOctalOrLiteralDigit()
        {
            char d = _p[Pos];
            if (d >= '8')
                return 0;                                                        // \8 and \9: a binary zero, the digit is the next literal (PCRE 7 follows the Perl of its time)

            int v = 0, i = 0;
            while (i < 3 && Pos < _p.Length && _p[Pos] >= '0' && _p[Pos] <= '7')
            {
                v = v * 8 + (_p[Pos] - '0');
                Pos++;
                i++;
            }

            return v & 0xFF;
        }

        // the sets of \d \D \s \S \w \W \h \H \v \V; Pos at the letter
        private static bool SimpleEscapeClass(char c, bool[] set)
        {
            Func<int, bool>? f = c switch
            {
                'd' => x => x is >= '0' and <= '9',
                'D' => x => !(x is >= '0' and <= '9'),
                's' => x => x is ' ' or '\t' or '\n' or '\u000b' or '\u000c' or '\r',
                'S' => x => !(x is ' ' or '\t' or '\n' or '\u000b' or '\u000c' or '\r'),
                'w' => x => IsWord((char)x),
                'W' => x => !IsWord((char)x),
                'h' => x => x is '\t' or ' ' or 0xA0,
                'H' => x => !(x is '\t' or ' ' or 0xA0),
                'v' => x => x is '\n' or '\u000b' or '\u000c' or '\r' or 0x85,
                'V' => x => !(x is '\n' or '\u000b' or '\u000c' or '\r' or 0x85),
                _ => null,
            };
            if (f is null)
                return false;
            for (int i = 0; i < 256; i++)
                set[i] = f(i);
            return true;
        }

        // a single-character escape: Pos at the character after the backslash
        private int ReadCharEscape(bool inClass)
        {
            if (!More)
                throw new EmbossPcreException("\\ at end of pattern", Pos);
            char c = _p[Pos++];
            switch (c)
            {
                case 'a': return 7;
                case 'e': return 27;
                case 'f': return 12;
                case 'n': return 10;
                case 'r': return 13;
                case 't': return 9;
                case '0':
                    {
                        int v = 0, i = 0;
                        while (i < 2 && Pos < _p.Length && _p[Pos] >= '0' && _p[Pos] <= '7')
                        {
                            v = v * 8 + (_p[Pos] - '0');
                            Pos++;
                            i++;
                        }

                        return v & 0xFF;
                    }

                case 'x':
                    {
                        if (More && _p[Pos] == '{')
                        {
                            int q = Pos + 1, v = 0, digits = 0;
                            while (q < _p.Length && Uri.IsHexDigit(_p[q]))
                            {
                                v = v * 16 + Convert.ToInt32(_p[q].ToString(), 16);
                                digits++;
                                q++;
                                if (v > 0xFF)
                                    throw new EmbossPcreException("character value in \\x{...} sequence is too large", q);
                            }

                            if (q < _p.Length && _p[q] == '}' && digits > 0)
                            {
                                Pos = q + 1;
                                return v;
                            }

                            return 0;                                           // "\x{" not followed by hex digits and }: \x is NUL followed by the text
                        }

                        int val = 0, n = 0;
                        while (n < 2 && Pos < _p.Length && Uri.IsHexDigit(_p[Pos]))
                        {
                            val = val * 16 + Convert.ToInt32(_p[Pos].ToString(), 16);
                            Pos++;
                            n++;
                        }

                        return val;
                    }

                case 'c':
                    {
                        if (!More)
                            throw new EmbossPcreException("\\c at end of pattern", Pos);
                        char x = _p[Pos++];
                        if (x >= 128)
                            throw new EmbossPcreException("\\c must be followed by an ASCII character", Pos);
                        if (x is >= 'a' and <= 'z')
                            x = (char)(x - 32);
                        return x ^ 0x40;
                    }

                default:
                    if (c is >= '1' and <= '9')
                    {
                        Pos--;
                        return ReadOctalOrLiteralDigit();
                    }

                    if (char.IsAsciiLetter(c))
                        return c;                                                // handled for the letters PCRE reads as literals; the refused letters are decided by the callers
                    return c;
            }
        }

        // ------------------------------------------------------------------------------------------------------------ classes

        private static readonly string[] PosixNames = ["alpha", "lower", "upper", "alnum", "ascii", "blank", "cntrl", "digit", "graph", "print", "punct", "space", "word", "xdigit"];

        private static bool PosixMember(string name, int c) => name switch
        {
            "alpha" => char.IsAsciiLetter((char)c),
            "lower" => c is >= 'a' and <= 'z',
            "upper" => c is >= 'A' and <= 'Z',
            "alnum" => char.IsAsciiLetterOrDigit((char)c),
            "ascii" => c < 128,
            "blank" => c is ' ' or '\t',
            "cntrl" => c < 32 || c == 127,
            "digit" => c is >= '0' and <= '9',
            "graph" => c is > 32 and < 127,
            "print" => c is >= 32 and < 127,
            "punct" => c is > 32 and < 127 && !char.IsAsciiLetterOrDigit((char)c),
            "space" => c is ' ' or '\t' or '\n' or '\u000b' or '\u000c' or '\r',
            "word" => IsWord((char)c),
            "xdigit" => Uri.IsHexDigit((char)c),
            _ => false,
        };

        // Pos is after '['
        private Node ParseClass()
        {
            var set = new bool[256];
            bool negate = false;
            if (More && _p[Pos] == '^')
            {
                negate = true;
                Pos++;
            }

            bool first = true;
            bool inQuote = false;
            int lastChar = -1;                                                   // the last single character (a range start); -1 when none
            bool pendingRange = false;
            while (true)
            {
                if (!More)
                    throw new EmbossPcreException("missing terminating ] for character class", _p.Length);
                char c = _p[Pos];
                if (c == ']' && !first && !inQuote)
                {
                    Pos++;
                    break;
                }

                if (inQuote)
                {
                    if (c == '\\' && Pos + 1 < _p.Length && _p[Pos + 1] == 'E')
                    {
                        Pos += 2;
                        inQuote = false;
                        continue;
                    }

                    Pos++;
                    first = false;
                    AddChar(set, c, ref lastChar, ref pendingRange);
                    continue;
                }

                if (c == '[' && !pendingRange && Pos + 1 < _p.Length && _p[Pos + 1] is ':' or '.' or '=' && CheckPosixSyntax(Pos, out int term))
                {
                    first = false;
                    if (_p[Pos + 1] != ':')
                        throw new EmbossPcreException("POSIX collating elements are not supported", Pos);
                    int ns = Pos + 2;
                    bool neg = false;
                    if (_p[ns] == '^')
                    {
                        neg = true;
                        ns++;
                    }

                    string name = _p.Substring(ns, term - ns);
                    if (Array.IndexOf(PosixNames, name) < 0)
                        throw new EmbossPcreException("unknown POSIX class name", ns);
                    for (int i = 0; i < 256; i++)
                        if (PosixMember(name, i) != neg)
                            set[i] = true;
                    Pos = term + 2;
                    lastChar = -1;
                    continue;
                }

                if (c == '\\')
                {
                    Pos++;
                    if (!More)
                        throw new EmbossPcreException("\\ at end of pattern", Pos);
                    char e = _p[Pos];
                    if (e == 'E')
                    {
                        Pos++;
                        continue;
                    }

                    if (e == 'Q')
                    {
                        Pos++;
                        inQuote = true;
                        continue;
                    }

                    first = false;
                    if (pendingRange && e is 'A' or 'B' or 'C' or 'G' or 'K' or 'Z' or 'z')
                    {
                        set['-'] = true;                                           // PCRE 7: a hyphen before an escape that is no character is a hyphen; the escape is read again as a single character
                        pendingRange = false;
                    }

                    var sub = new bool[256];
                    if (SimpleEscapeClass(e, sub))
                    {
                        Pos++;
                        if (pendingRange)
                        {
                            set['-'] = true;                                       // PCRE 7: a hyphen before a class escape is a hyphen
                            pendingRange = false;
                        }

                        for (int i = 0; i < 256; i++)
                            if (sub[i])
                                set[i] = true;
                        lastChar = -1;
                        continue;
                    }

                    int ch;
                    switch (e)
                    {
                        case 'b':
                            Pos++;
                            ch = 8;
                            break;
                        case 'B':
                        case 'R':
                        case 'X':
                        case 'p':
                        case 'P':
                            Pos++;
                            ch = e;                                              // PCRE 7 reads \B, \R, \X, \p and \P in a class as the letters
                            break;
                        case 'L':
                        case 'l':
                        case 'N':
                        case 'U':
                        case 'u':
                            Pos++;
                            throw new EmbossPcreException("PCRE does not support \\L, \\l, \\N, \\U, or \\u", Pos - 1);
                        case 'g':
                            Pos++;
                            ClassGEscape();
                            continue;
                        default:
                            ch = ReadCharEscape(true);
                            break;
                    }

                    AddChar(set, ch, ref lastChar, ref pendingRange);
                    continue;
                }

                if (c == '-' && lastChar >= 0 && !pendingRange && Pos + 1 < _p.Length && _p[Pos + 1] != ']')
                {
                    Pos++;
                    pendingRange = true;
                    continue;
                }

                Pos++;
                first = false;
                AddChar(set, c, ref lastChar, ref pendingRange);
            }

            if (pendingRange)
                set['-'] = true;
            if (negate)
                for (int i = 0; i < 256; i++)
                    set[i] = !set[i];
            return new SetNode(set);
        }

        // check_posix_syntax: the text from the '[' at <from> is [:name:] / [.x.] / [=x=] (a terminator followed by ']' before any ']')
        private bool CheckPosixSyntax(int from, out int end)
        {
            end = -1;
            char terminator = _p[from + 1];
            for (int q = from + 2; q < _p.Length; q++)
            {
                if (_p[q] == '\\' && q + 1 < _p.Length && _p[q + 1] is ']' or '\\')
                {
                    q++;
                    continue;
                }

                if (_p[q] == ']')
                    return false;
                if (_p[q] == terminator && q + 1 < _p.Length && _p[q + 1] == ']')
                {
                    end = q;
                    return true;
                }
            }

            return false;
        }

        private void AddChar(bool[] set, int ch, ref int lastChar, ref bool pendingRange)
        {
            ch &= 0xFF;
            if (pendingRange)
            {
                if (ch < lastChar)
                    throw new EmbossPcreException("range out of order in character class", Pos - 1);
                for (int i = lastChar; i <= ch; i++)
                    set[i] = true;
                pendingRange = false;
                lastChar = -1;
                return;
            }

            set[ch] = true;
            lastChar = ch;
        }
    }
}
