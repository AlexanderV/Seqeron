// RNA-DOTBRACKET-001 — Dot-Bracket (extended WUSS) Notation
// Evidence: docs/Evidence/RNA-DOTBRACKET-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-DOTBRACKET-001.md
// Source: ViennaRNA Package — RNA Structure Notations / Dot-Bracket / WUSS (Lorenz et al. 2011);
//         Nawrocki & Eddy (2013) Infernal 1.1, Bioinformatics 29(22):2933-2935; Rfam glossary (WUSS).

using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_ParseDotBracket_Tests
{
    private static HashSet<(int, int)> PairsOf(string s) =>
        ParseDotBracket(s).ToHashSet();

    #region ParseDotBracket

    // M1 — Simple nested hairpin "((((....))))": 4 bp, outermost-with-outermost.
    // Evidence: ViennaRNA basic notation, example "((((....))))"; INV-01 (i<j), INV-03 (#pairs == #openers).
    [Test]
    public void ParseDotBracket_SimpleHairpin_ReturnsExactNestedPairs()
    {
        var pairs = PairsOf("((((....))))");

        Assert.That(pairs, Is.EquivalentTo(new[] { (0, 11), (1, 10), (2, 9), (3, 8) }),
            "((((....)))) has 4 nested pairs; position i pairs with 11-i — exact full set required");
    }

    // M2 — Crossing families "([)]" must be matched on independent stacks: (0,2) and (1,3).
    // A single shared stack would wrongly yield (1,2). Evidence: ViennaRNA extended notation /
    // WUSS — different bracket families need not nest; INV-02 (same family per pair).
    [Test]
    public void ParseDotBracket_CrossingFamilies_PairsEachFamilyIndependently()
    {
        var pairs = PairsOf("([)]");

        Assert.That(pairs, Is.EquivalentTo(new[] { (0, 2), (1, 3) }),
            "'(' pairs with ')' at (0,2) and '[' with ']' at (1,3); a shared stack would give the wrong (1,2)");
    }

    // M3 — Equivalence of bracket families and letter pairs: the two retrieved equivalent encodings
    // of the same crossing structure must produce identical pair sets.
    // Evidence: ViennaRNA verbatim examples "<<<<[[[[....>>>>]]]]" and "((((AAAA....))))aaaa".
    [Test]
    public void ParseDotBracket_BracketAndLetterEncodings_AreEquivalent()
    {
        // Within each family brackets nest (LIFO): the outermost '<' (idx 0) pairs with the
        // outermost '>' (idx 15), and the outermost '[' (idx 4) with the outermost ']' (idx 19).
        var expected = new[]
        {
            (0, 15), (1, 14), (2, 13), (3, 12), // <<<< ... >>>> helix (nested)
            (4, 19), (5, 18), (6, 17), (7, 16), // [[[[ ... ]]]] crossing helix (nested)
        };

        var byBrackets = PairsOf("<<<<[[[[....>>>>]]]]");
        var byLetters = PairsOf("((((AAAA....))))aaaa");

        Assert.Multiple(() =>
        {
            Assert.That(byBrackets, Is.EquivalentTo(expected),
                "<>/[] crossing helices: <<<< pairs with >>>> and [[[[ with ]]]] per ViennaRNA example");
            Assert.That(byLetters, Is.EquivalentTo(expected),
                "letter-pair encoding ((((AAAA....))))aaaa must yield the same pairs as the bracket encoding");
            Assert.That(byLetters, Is.EquivalentTo(byBrackets),
                "ViennaRNA states the two encodings are equivalent representations of one structure");
        });
    }

    // M6 — Pair count equals number of opening symbols for a well-formed string.
    // Evidence: ViennaRNA balanced requirement; INV-03. "(([[]]))" = 4 openers → 4 pairs with positions.
    [Test]
    public void ParseDotBracket_NestedMixedFamilies_ReturnsFourPairsWithPositions()
    {
        var pairs = PairsOf("(([[]]))");

        Assert.That(pairs, Is.EquivalentTo(new[] { (0, 7), (1, 6), (2, 5), (3, 4) }),
            "(([[]])) nests ( ( [ [ ] ] ) ): pairs (0,7),(1,6) round and (2,5),(3,4) square");
    }

    // S1 — Letter direction: uppercase opens (5'), lowercase closes (3'); "AAAA....aaaa" nests like brackets.
    // Evidence: ViennaRNA equivalent example uses AAAA on the 5' side and aaaa on the 3' side.
    [Test]
    public void ParseDotBracket_UppercaseLowercaseLetters_UppercaseIsOpener()
    {
        var pairs = PairsOf("AAAA....aaaa");

        Assert.That(pairs, Is.EquivalentTo(new[] { (0, 11), (1, 10), (2, 9), (3, 8) }),
            "uppercase A is the 5' opener and lowercase a the 3' closer, nested like ((((....))))");
    }

    // S2 — Non-bracket WUSS symbols (-, ,) are single-stranded and must not break pairing.
    // Evidence: Rfam glossary — '-', ',', ':' are single-stranded residues; only brackets pair.
    [Test]
    public void ParseDotBracket_NonBracketWussSymbols_TreatedAsUnpaired()
    {
        var dashes = PairsOf("<<<-->>>");
        var commas = PairsOf("((,,))");

        Assert.Multiple(() =>
        {
            Assert.That(dashes, Is.EquivalentTo(new[] { (0, 7), (1, 6), (2, 5) }),
                "'-' is single-stranded; <<<-->>> pairs the angle brackets across the dashes");
            Assert.That(commas, Is.EquivalentTo(new[] { (0, 5), (1, 4) }),
                "',' is single-stranded; ((,,)) pairs only the parentheses");
        });
    }

    // S3 — Best-effort parse: a stray closing bracket is dropped, not thrown.
    // Evidence: documented contract (Evidence Assumption 1).
    [Test]
    public void ParseDotBracket_StrayCloser_DroppedNotThrown()
    {
        var pairs = PairsOf("())");

        Assert.That(pairs, Is.EquivalentTo(new[] { (0, 1) }),
            "the matched pair (0,1) is returned and the trailing unmatched ')' is dropped");
    }

    // C1 — Empty / all-dots / null inputs produce no pairs.
    // Evidence: dots are unpaired (ViennaRNA); empty/null contract (Evidence Assumption 2).
    [Test]
    public void ParseDotBracket_EmptyDotsOrNull_ReturnsNoPairs()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ParseDotBracket("").Any(), Is.False, "empty string encodes no base pairs");
            Assert.That(ParseDotBracket(".....").Any(), Is.False, "all-dots string is fully single-stranded");
            Assert.That(ParseDotBracket(null!).Any(), Is.False, "null is treated as an empty structure");
        });
    }

    // S4 — Letter pairs mix freely with brackets and are independent systems.
    // Reference: ViennaRNA 2.7.2 RNA.ptable(s, RNA.BRACKETS_ANY) → "(A)a" pairs (1,3),(2,4) 1-based
    // = (0,2),(1,3) 0-based; "AB.ba" → (0,4),(1,3); "((..#..))" (unknown '#' is unpaired) → (0,8),(1,7).
    [Test]
    public void ParseDotBracket_LettersMixedWithBrackets_MatchViennaRnaPtable()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PairsOf("(A)a"), Is.EquivalentTo(new[] { (0, 2), (1, 3) }),
                "'(' and 'A' are independent families: crossing pairs (0,2) and (1,3) — ViennaRNA ptable");
            Assert.That(PairsOf("AB.ba"), Is.EquivalentTo(new[] { (0, 4), (1, 3) }),
                "A/a and B/b are distinct letter families — ViennaRNA ptable");
            Assert.That(PairsOf("((..#..))"), Is.EquivalentTo(new[] { (0, 8), (1, 7) }),
                "non-bracket symbol '#' is unpaired — ViennaRNA ptable ignores it");
        });
    }

    // S5 — Only ASCII A–Z / a–z are pseudoknot letters. ViennaRNA 2.7.2 vrna_ptable_from_string
    // (structures/structure_pairtable.c, VRNA_BRACKETS_ALPHA) loops i = 65..90 and pairs (char)i
    // with (char)(i+32); any other character is unpaired. Reference (non-ASCII masked to '.'):
    // RNA.ptable("ÉÉ..éé") → no pairs; "Σ((..))σ" → (1,6),(2,5) only.
    // Regression: char.IsLetter previously paired 'É'/'é' and 'Σ'/'σ'.
    [Test]
    public void ParseDotBracket_NonAsciiLetters_AreUnpaired()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ParseDotBracket("ÉÉ..éé").Any(), Is.False,
                "'É'/'é' are not ViennaRNA letter families (A–Z only) → unpaired");
            Assert.That(PairsOf("Σ((..))σ"), Is.EquivalentTo(new[] { (1, 6), (2, 5) }),
                "Greek Σ/σ are unpaired; only the parentheses pair");
        });
    }

    // S6 — Round trip: the dot-bracket rendered by the MFE predictor (canonical BuildDotBracket
    // renderer) parses back to exactly the predictor's own base-pair list (vrna_db_from_ptable ∘
    // vrna_ptable is the identity for nested structures).
    [TestCase("GGGAAACCC")]
    [TestCase("GGGGAAAACCCCAUAUGGGGAAAACCCC")]
    [TestCase("GCGCUUCGGCGCAUCGAUCGAAAACGAUCGAU")]
    public void ParseDotBracket_MfeDotBracket_RoundTripsToBasePairs(string seq)
    {
        var mfe = CalculateMfeStructure(seq);

        Assert.Multiple(() =>
        {
            Assert.That(ValidateDotBracket(mfe.DotBracket), Is.True, "rendered MFE structure is well-formed");
            Assert.That(PairsOf(mfe.DotBracket), Is.EquivalentTo(mfe.BasePairs),
                "parse(render(pairs)) == pairs");
        });
    }

    #endregion

    #region ValidateDotBracket

    // M4 — Well-formed strings (balanced/nested AND crossing families) validate true.
    // Evidence: ViennaRNA balanced requirement; WUSS crossing families "([)]"; INV-04.
    [Test]
    public void ValidateDotBracket_WellFormed_ReturnsTrue()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ValidateDotBracket("(((...)))"), Is.True, "balanced, nested parentheses");
            Assert.That(ValidateDotBracket("(([[]]))"), Is.True, "nested mixed families, each balanced");
            Assert.That(ValidateDotBracket("([)]"), Is.True,
                "crossing families: '(' and '[' each have a matching closer — valid pseudoknot per WUSS");
            Assert.That(ValidateDotBracket("...."), Is.True, "all unpaired");
            Assert.That(ValidateDotBracket(""), Is.True, "empty structure is balanced");
            Assert.That(ValidateDotBracket(null!), Is.True, "null treated as empty structure");
        });
    }

    // M5 — Malformed strings validate false, including the mismatched-family case "(]" that a
    // single-counter validator would wrongly accept.
    // Evidence: ViennaRNA balanced requirement; WUSS "partners must match up"; INV-04.
    [Test]
    public void ValidateDotBracket_Malformed_ReturnsFalse()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ValidateDotBracket("(((...)"), Is.False, "two '(' left unclosed");
            Assert.That(ValidateDotBracket("...)"), Is.False, "')' has no opening partner");
            Assert.That(ValidateDotBracket(")("), Is.False, "')' before any '(' — closer precedes opener");
            Assert.That(ValidateDotBracket("(]"), Is.False,
                "mismatched families: '(' is unclosed and ']' unopened — partners must match up");
        });
    }

    // M7 — Validation agrees with ViennaRNA 2.7.2 RNA.ptable(s, RNA.BRACKETS_ANY) != None
    // (vrna_ptable_from_string returns NULL on unbalanced brackets of any family).
    // Reference: "aaaa....AAAA" → NULL (closer before opener); "A.a" → valid; "(A)a" → valid;
    // "É..." → valid (non-ASCII is unpaired, not an unclosed opener).
    [Test]
    public void ValidateDotBracket_LetterFamilies_MatchViennaRnaPtable()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ValidateDotBracket("aaaa....AAAA"), Is.False,
                "lowercase closers precede their uppercase openers — ViennaRNA ptable returns NULL");
            Assert.That(ValidateDotBracket("A.a"), Is.True, "A/a balanced — ViennaRNA ptable valid");
            Assert.That(ValidateDotBracket("(A)a"), Is.True, "crossing ( and A families — ViennaRNA ptable valid");
            Assert.That(ValidateDotBracket("É..."), Is.True,
                "non-ASCII 'É' is unpaired (ViennaRNA letters are A–Z only), not an unclosed opener");
        });
    }

    #endregion
}
