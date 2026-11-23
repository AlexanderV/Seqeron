using System.Collections.Concurrent;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Contract;

/// <summary>
/// The definitional oracle of a text with the answers that the random walks ask for again and again cached (the same text is walked by every engine): occurrences, loci, parents and children of a string.
/// </summary>
internal sealed class WalkOracle
{
    private static readonly ConcurrentDictionary<string, WalkOracle> Known = new(StringComparer.Ordinal);

    private readonly ConcurrentDictionary<string, List<int>> _occurrences = new(StringComparer.Ordinal);
    private readonly ConcurrentDictionary<string, IndexLocus> _loci = new(StringComparer.Ordinal);
    private readonly ConcurrentDictionary<string, int[]?> _parents = new(StringComparer.Ordinal);
    private readonly ConcurrentDictionary<string, List<(int Key, int[] Child)>> _children = new(StringComparer.Ordinal);
    private readonly Lazy<List<int[]>> _all;

    private WalkOracle(string text)
    {
        Oracle = new IndexContractOracle(text);
        _all = new Lazy<List<int[]>>(() => AliasCatalog.Loci(Oracle));
    }

    internal static WalkOracle Of(string text) => Known.GetOrAdd(text, t => new WalkOracle(t));

    internal IndexContractOracle Oracle { get; }

    /// <summary>Every locus string of the text (the root first).</summary>
    internal List<int[]> AllLoci => _all.Value;

    internal List<int> Occurrences(int[] w) => _occurrences.GetOrAdd(IndexContractOracle.Key(w), _ => Oracle.Occ(w));

    internal bool Exists(int[] w) => Occurrences(w).Count > 0;

    internal bool IsExplicit(int[] w) => w.Length == 0 || Oracle.IsExplicit(w);

    internal IndexLocus LocusOf(int[] w) => _loci.GetOrAdd(IndexContractOracle.Key(w), _ => Oracle.LocusOf(w));

    internal int[]? Parent(int[] w) => _parents.GetOrAdd(IndexContractOracle.Key(w), _ => Oracle.Parent(w));

    internal List<(int Key, int[] Child)> Children(int[] w) => _children.GetOrAdd(IndexContractOracle.Key(w), _ => Oracle.Children(w));

    /// <summary>The string with <paramref name="symbol"/> appended when that string occurs, else null (a value that is no code unit never occurs).</summary>
    internal int[]? ExtendRight(int[] w, int symbol)
    {
        if (symbol is < 0 or > 0xFFFF) return null;
        int[] extended = [.. w, symbol];
        return Exists(extended) ? extended : null;
    }

    /// <summary>The string with <paramref name="symbol"/> prepended when that string occurs, else null.</summary>
    internal int[]? ExtendLeft(int symbol, int[] w)
    {
        if (symbol is < 0 or > 0xFFFF) return null;
        int[] extended = [symbol, .. w];
        return Exists(extended) ? extended : null;
    }
}

/// <summary>
/// Random walks of loci for the contract of the aliased result (work package 1b.12, finding R01 of the review of 1b.10): on every text of a family (the texts of the differential tests: random DNA, the Fibonacci word, a periodic
/// text, a^n, a protein-like text, Thue-Morse, runs, abracadabra, mississippi, banana, the empty text, one letter, surrogates and the characters 0 and FFFF) and for every engine, a few hundred walks start at the root, take about a dozen
/// steps and jump to a random locus now and then. Each step picks a member (<c>TryExtendRight</c> with a symbol or a pattern, <c>TryExtendLeft</c>, <c>TryParent</c>, <c>TryChild</c>, <c>TrySuffixLink</c>, <c>GetChildren</c>) and
/// an argument that continues the current string (most of the time), a symbol that does not occur or a value that is no code unit, calls it with two variables and with ONE variable for the locus and the result, and checks that the answers
/// (the bool and the locus, the number and the content of the children) are equal, and equal to the answer of the definitional oracle (<see cref="IndexContractOracle"/>); the engine that has no navigation (the Ukkonen adapter of the oracle)
/// is checked against the oracle on its count. The walk is driven by the aliased results. "The end of the text" is reached by the walks that follow a suffix to its leaf (a pattern that ends at the terminator, <c>TryChild(node, -1)</c>, the
/// leaf of the empty suffix) and by the symbols that are not code units; a locus that is not an explicit node is given to the members that need one, and both forms must refuse it. The seeds are fixed.
/// </summary>
[TestFixtureSource(typeof(IndexEngineFactories), nameof(IndexEngineFactories.All))]
public class IndexAliasingPropertyTests(IndexEngineFactory factory)
{
    private const int WalksPerText = 40;
    private const int StepsPerWalk = 12;

    private static readonly (string Name, string Text)[] Texts =
    [
        ("abracadabra", "abracadabra"),
        ("mississippi", "mississippi"),
        ("banana", "banana"),
        ("a^30", new string('a', 30)),
        ("dna-110", StringFamilies.RandomDna(110, 17)),
        ("fibonacci-90", StringFamilies.Fibonacci(90)),
        ("periodic-5", StringFamilies.Periodic(100, 5, 3)),
        ("protein-100", StringFamilies.RandomProtein(100, 9)),
        ("thue-morse-80", StringFamilies.ThueMorse(80)),
        ("runs-90", StringFamilies.Runs(90, 5, 4)),
        ("one-letter", "a"),
        ("empty", string.Empty),
        ("surrogates", "\uD800a\uD800"),
        ("nul-and-ffff", "\0￿\0a"),
    ];

    private enum Member { ExtendRightSymbol, ExtendRightPattern, ExtendLeft, Parent, Child, SuffixLink, GetChildren }

    [Test]
    public void RandomWalks_TheAliasedCallAndTheSeparateCallGiveTheSameAnswer_AndTheAnswerOfTheOracle()
    {
        var ok = new Dictionary<Member, int>();
        var failed = new Dictionary<Member, int>();
        int refused = 0, walks = 0;
        for (int t = 0; t < Texts.Length; t++)
        {
            (string name, string text) = Texts[t];
            if (!factory.AcceptsText(text)) continue;
            ISuffixIndex index = factory.Create(text);
            var walker = new Walker(index, WalkOracle.Of(text), $"[{factory.Name}] {name}", ok, failed);
            for (int walk = 0; walk < WalksPerText; walk++)
            {
                walker.Walk(new SplitMix64((ulong)((t * 1000) + walk + 1)), walk);
                walks++;
            }
            refused += walker.Refused;
        }
        TestContext.Out.WriteLine($"[{factory.Name}] {walks} walks, {StepsPerWalk} steps each; answered / not answered per member: {string.Join(", ", Enum.GetValues<Member>().Select(m => $"{m} {ok.GetValueOrDefault(m)}/{failed.GetValueOrDefault(m)}"))}; refused by both forms: {refused}");
        Assert.That(walks, Is.GreaterThan(300), "a few hundred walks");
        bool navigation = factory.Name != "ukkonen-adapter";
        foreach (Member member in Enum.GetValues<Member>())
        {
            if (!navigation && member > Member.ExtendRightPattern) continue;
            Assert.That(ok.GetValueOrDefault(member), Is.GreaterThanOrEqualTo(20), $"[{factory.Name}] walks that took {member} and got an answer");
            Assert.That(failed.GetValueOrDefault(member), Is.GreaterThanOrEqualTo(5), $"[{factory.Name}] walks that took {member} and got none");
        }
        if (navigation) Assert.That(refused, Is.GreaterThanOrEqualTo(20), $"[{factory.Name}] calls on a locus that is no explicit node, refused by both forms");
    }

    /// <summary>One walker per engine and text: the state of a walk (the string and the locus of the engine) and the checks of a step.</summary>
    private sealed class Walker
    {
        private static readonly int[] Absent = ['z', 1, 0xFFFE];
        private static readonly int[] NotCodeUnits = [-1, -2, 65536, int.MaxValue, int.MinValue];

        private readonly ISuffixIndex _index;
        private readonly ISuffixTreeNavigation? _nav;
        private readonly WalkOracle _oracle;
        private readonly string _text;
        private readonly string _context;
        private readonly bool _esaLoci;
        private readonly int[] _alphabet;
        private readonly Dictionary<Member, int> _ok;
        private readonly Dictionary<Member, int> _failed;
        private readonly List<string> _trail = [];
        private SplitMix64 _rng = new(1);
        private int[] _w = [];
        private IndexLocus _locus;
        private int _walk;

        internal Walker(ISuffixIndex index, WalkOracle oracle, string context, Dictionary<Member, int> ok, Dictionary<Member, int> failed)
        {
            _index = index;
            _nav = index as ISuffixTreeNavigation;
            _oracle = oracle;
            _text = oracle.Oracle.Text;
            _context = context;
            _esaLoci = _nav is not null; // the other engine of the suite numbers its loci by text position
            int[] alphabet = oracle.Oracle.Alphabet();
            _alphabet = alphabet.Length > 0 ? alphabet : ['a'];
            _ok = ok;
            _failed = failed;
        }

        /// <summary>The number of calls that both forms refused because the locus was no explicit node.</summary>
        internal int Refused { get; private set; }

        internal void Walk(SplitMix64 rng, int walk)
        {
            _rng = rng;
            _walk = walk;
            _trail.Clear();
            _w = [];
            _locus = _index.Root;
            for (int step = 0; step < StepsPerWalk; step++) Step();
        }

        private void Step()
        {
            int roll = _rng.NextInt(_nav is null ? 6 : 14);
            switch (roll)
            {
                case < 3: ExtendRightSymbol(); break;
                case < 5: ExtendRightPattern(); break;
                case < 7:
                    if (_nav is null) Jump();
                    else ExtendLeft();
                    break;
                case 7: Parent(); break;
                case < 10: Child(); break;
                case 10: SuffixLink(); break;
                case 11: GetChildren(); break;
                default: Jump(); break;
            }
        }

        private void Jump()
        {
            _trail.Add("jump");
            if (!_esaLoci || _rng.NextInt(3) == 0)
            {
                _w = [];
                _locus = _index.Root;
                return;
            }
            _w = _oracle.AllLoci[_rng.NextInt(_oracle.AllLoci.Count)];
            _locus = _oracle.LocusOf(_w);
        }

        private void Fail(string message)
            => Assert.Fail($"{_context} walk {_walk}, step {_trail.Count}: at {IndexContractOracle.Render(_w)} {_locus}: {message}. Trail: {string.Join(" ", _trail)}");

        private void Tally(Member member, bool answered)
        {
            Dictionary<Member, int> counts = answered ? _ok : _failed;
            counts[member] = counts.GetValueOrDefault(member) + 1;
        }

        /// <summary>The answer of one call: equal for the two forms, and equal to the oracle (<paramref name="expected"/>, null for none).</summary>
        private void Verify(Member member, string call, Outcome apart, Outcome alias, int[]? expected)
        {
            _trail.Add(call);
            if (!alias.Equals(apart)) Fail($"{call}: one variable for the locus and the result gave {alias}, two variables gave {apart}");
            if (_esaLoci)
            {
                Outcome wanted = expected is null ? new Outcome(false, default) : new Outcome(true, _oracle.LocusOf(expected));
                if (!apart.Equals(wanted)) Fail($"{call}: the oracle says {wanted}, the engine says {apart}");
            }
            else
            {
                if (apart.Ok != (expected is not null)) Fail($"{call}: the oracle says {(expected is not null)}, the engine says {apart.Ok}");
                if (!apart.Ok && !apart.Locus.Equals(default(IndexLocus))) Fail($"{call}: a failed call must leave the default locus, got {apart.Locus}");
                if (apart.Ok && (_index.Count(apart.Locus) != _oracle.Oracle.Count(expected!) || apart.Locus.Length != expected!.Length)) Fail($"{call}: the locus {apart.Locus} has {_index.Count(apart.Locus)} occurrences and depth {apart.Locus.Length}, the oracle has {_oracle.Oracle.Count(expected!)} and {expected!.Length}");
            }
            Tally(member, apart.Ok);
            if (alias.Ok)
            {
                _w = expected!;
                _locus = alias.Locus; // the walk goes on with the aliased result
            }
        }

        private void RefusedByBoth(string call, Action apart, Action alias)
        {
            _trail.Add(call + "!");
            if (!Refuses(apart) || !Refuses(alias)) Fail($"{call}: a locus that is no explicit node must be refused with an ArgumentException by both forms");
            Refused++;
        }

        private static bool Refuses(Action call)
        {
            try
            {
                call();
                return false;
            }
            catch (ArgumentException)
            {
                return true;
            }
        }

        private int PickSymbol(bool right)
        {
            int roll = _rng.NextInt(20);
            if (roll < 11)
            {
                List<int> next = [];
                foreach (int i in _oracle.Occurrences(_w))
                {
                    int at = right ? i + _w.Length : i - 1;
                    if (at >= 0 && at < _text.Length) next.Add(_text[at]);
                }
                if (next.Count > 0) return next[_rng.NextInt(next.Count)];
            }
            if (roll < 16) return _alphabet[_rng.NextInt(_alphabet.Length)];
            if (roll < 18) return Absent[_rng.NextInt(Absent.Length)];
            return NotCodeUnits[_rng.NextInt(NotCodeUnits.Length)];
        }

        private void ExtendRightSymbol()
        {
            int symbol = PickSymbol(right: true);
            Outcome apart = AliasProbe.ExtendRightApart(_index, _locus, symbol);
            Outcome alias = AliasProbe.ExtendRightAlias(_index, _locus, symbol);
            Verify(Member.ExtendRightSymbol, $"R{symbol}", apart, alias, _oracle.ExtendRight(_w, symbol));
        }

        private void ExtendRightPattern()
        {
            string pattern;
            List<int> starts = [.. _oracle.Occurrences(_w).Where(i => i + _w.Length < _text.Length)];
            if (starts.Count > 0 && _rng.NextInt(4) > 0)
            {
                int at = starts[_rng.NextInt(starts.Count)] + _w.Length;
                pattern = _text.Substring(at, Math.Min(_rng.NextInt(6), _text.Length - at));
                if (pattern.Length > 0 && _rng.NextInt(5) == 0) pattern = pattern[..^1] + 'z';
            }
            else
            {
                char[] garbage = new char[_rng.NextInt(4)];
                for (int k = 0; k < garbage.Length; k++) garbage[k] = (char)_alphabet[_rng.NextInt(_alphabet.Length)];
                pattern = new string(garbage);
            }
            int[] extended = [.. _w, .. IndexContractOracle.Symbols(pattern)];
            Outcome apart = AliasProbe.ExtendRightApart(_index, _locus, pattern);
            Outcome alias = AliasProbe.ExtendRightAlias(_index, _locus, pattern);
            Verify(Member.ExtendRightPattern, $"R\"{pattern}\"", apart, alias, _oracle.Exists(extended) ? extended : null);
        }

        private void ExtendLeft()
        {
            int symbol = PickSymbol(right: false);
            Outcome apart = AliasProbe.ExtendLeftApart(_index, _locus, symbol);
            Outcome alias = AliasProbe.ExtendLeftAlias(_index, _locus, symbol);
            Verify(Member.ExtendLeft, $"L{symbol}", apart, alias, _oracle.ExtendLeft(symbol, _w));
        }

        private void Parent()
        {
            ISuffixTreeNavigation nav = _nav!;
            Outcome apart = AliasProbe.ParentApart(nav, _locus);
            Outcome alias = AliasProbe.ParentAlias(nav, _locus);
            Verify(Member.Parent, "P", apart, alias, _oracle.Parent(_w));
        }

        private int PickKey()
        {
            List<(int Key, int[] Child)> children = _oracle.Children(_w);
            int roll = _rng.NextInt(10);
            if (roll < 6 && children.Count > 0) return children[_rng.NextInt(children.Count)].Key;
            if (roll < 7) return -1;
            if (roll < 9) return _alphabet[_rng.NextInt(_alphabet.Length)];
            return roll == 9 && _rng.NextInt(2) == 0 ? Absent[_rng.NextInt(Absent.Length)] : NotCodeUnits[_rng.NextInt(NotCodeUnits.Length)];
        }

        private void Child()
        {
            ISuffixTreeNavigation nav = _nav!;
            int key = PickKey();
            IndexLocus node = _locus;
            if (!_oracle.IsExplicit(_w))
            {
                RefusedByBoth($"C{key}", () => AliasProbe.ChildApart(nav, node, key), () => AliasProbe.ChildAlias(nav, node, key));
                return;
            }
            int[]? child = _oracle.Children(_w).Where(k => k.Key == key).Select(k => k.Child).FirstOrDefault();
            Verify(Member.Child, $"C{key}", AliasProbe.ChildApart(nav, node, key), AliasProbe.ChildAlias(nav, node, key), child);
        }

        private void SuffixLink()
        {
            ISuffixTreeNavigation nav = _nav!;
            IndexLocus node = _locus;
            if (!_oracle.IsExplicit(_w))
            {
                RefusedByBoth("S", () => AliasProbe.SuffixLinkApart(nav, node), () => AliasProbe.SuffixLinkAlias(nav, node));
                return;
            }
            int[]? link = IndexContractOracle.SuffixLink(_w);
            bool isInternal = link is not null && _oracle.Children(_w).Count > 0;
            Verify(Member.SuffixLink, "S", AliasProbe.SuffixLinkApart(nav, node), AliasProbe.SuffixLinkAlias(nav, node), isInternal ? link : null);
        }

        private void GetChildren()
        {
            ISuffixTreeNavigation nav = _nav!;
            IndexLocus node = _locus;
            if (!_oracle.IsExplicit(_w))
            {
                RefusedByBoth("K", () => AliasProbe.ChildrenApart(nav, node, 3), () => AliasProbe.ChildrenAlias(nav, node, 3, 1));
                return;
            }
            IndexLocus[] children = [.. _oracle.Children(_w).Select(k => _oracle.LocusOf(k.Child))];
            int[] sizes = [children.Length, children.Length, children.Length + 1, children.Length + 3, Math.Max(0, children.Length - 1)];
            int size = sizes[_rng.NextInt(sizes.Length)];
            _trail.Add($"K{size}/{children.Length}");
            ChildrenOutcome apart = AliasProbe.ChildrenApart(nav, node, size);
            CheckChildren("two variables", apart, children, size, default, -1);
            if (size > 0)
            {
                int slot = _rng.NextInt(size);
                CheckChildren($"the node in slot {slot} of the destination", AliasProbe.ChildrenAlias(nav, node, size, slot), children, size, node, slot);
            }
            Tally(Member.GetChildren, children.Length > 0);
        }

        private void CheckChildren(string form, ChildrenOutcome got, IndexLocus[] children, int size, IndexLocus node, int slot)
        {
            IndexLocus[] expected = AliasProbe.ExpectedBuffer(children, size, node, slot);
            if (got.Count != children.Length) Fail($"GetChildren({form}, size {size}) returned {got.Count}, the oracle has {children.Length} children");
            if (!got.Buffer.SequenceEqual(expected)) Fail($"GetChildren({form}, size {size}): expected {new ChildrenOutcome(children.Length, expected)}, got {got}");
        }
    }
}
