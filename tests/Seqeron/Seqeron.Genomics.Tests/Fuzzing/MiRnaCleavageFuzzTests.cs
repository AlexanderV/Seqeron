using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for MIRNA-CLEAVAGE-001 — Drosha/Dicer cleavage-site prediction on a folded pri-miRNA
/// hairpin (<see cref="MiRnaAnalyzer.PredictDroshaDicerCleavage(string, int)"/> and the structure /
/// counting-rule overload).
///
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8): degenerate, boundary and out-of-domain
/// inputs never hang or throw an unhandled runtime exception, and every result is either a well-defined
/// record satisfying the published contract, the documented null ("cannot predict"), or the documented
/// ArgumentException for a malformed structure. Checklist: docs/checklists/03_FUZZING.md, row 255
/// ("precursor too short, no stem, non-ACGU, empty seq"). Strategies: BE = boundary exploitation,
/// MC = malformed content.
///
/// The contract (re-derived from the sources, not read off the code):
///   • Han 2006: count 11 base PAIRS up the unbranched folded stem from the junction pair; the 5' end of the
///     pre-miRNA is the 5' base of the 12th pair (bulged nucleotides do not count).
///   • RNase III 2-nt 3' overhang: Drosha 3' end = partner(5' end) + 2 (miRDeep2 find_star).
///   • Park 2011 5' counting: 5p = [5' end, 5' end + 21]; the 3p lies 3' of the 5p (opposite arm).
///   INV-CV-1: null/empty, junction ∉ [0,len), junction not the 5' base of a pair, fewer than 12 stem
///             pairs, or products off the sequence / not on opposite arms ⇒ null, never an exception.
///   INV-CV-2: every coordinate ∈ [0,len); spans ordered; MatureEnd &lt; StarStart; sub-sequences equal the
///             indexed spans of the normalised (upper, T→U) sequence; ThreePrimeOverhang = 2.
///   INV-CV-3: DroshaCut5' is the 12th pair's 5' base on the reported structure; DroshaCut3' = partner + 2;
///             5p length 22.
///   INV-CV-4: a malformed structure (unbalanced or wrong length) ⇒ ArgumentException.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MiRnaCleavageFuzzTests
{
    #region Helpers

    private static string RandomRna(int length, int seed)
    {
        const string bases = "ACGU";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    private const int LowerStemBp = 11; // Han 2006
    private const int MatureLength = 22; // Park 2011
    private const int Overhang = 2;      // Lee 2003 / Han 2006

    // Independent pair table (simple '(' ')' stack — the structures here use only round brackets).
    private static int[] Partners(string db)
    {
        var p = Enumerable.Repeat(-1, db.Length).ToArray();
        var st = new Stack<int>();
        for (int i = 0; i < db.Length; i++)
        {
            if (db[i] == '(') st.Push(i);
            else if (db[i] == ')') { int j = st.Pop(); p[i] = j; p[j] = i; }
        }
        return p;
    }

    /// <summary>INV-CV-2/3: re-derives the 12th stem pair and the overhang from the reported structure.</summary>
    private static void AssertCleavageContract(DroshaDicerCleavage c, int junction)
    {
        int len = c.Sequence.Length;
        c.Structure.Length.Should().Be(len);
        var p = Partners(c.Structure);

        // Walk 11 pairs up the stem: next pair = first paired base 3' of i whose partner is the last paired
        // base 5' of j (no branch).
        int i = junction, j = p[junction];
        j.Should().BeGreaterThan(i, "the junction is the 5' base of a pair (INV-CV-3)");
        for (int k = 0; k < LowerStemBp; k++)
        {
            int a = i + 1; while (p[a] < 0) a++;
            int b = j - 1; while (p[b] < 0) b--;
            p[a].Should().Be(b, "the lower stem is unbranched");
            i = a; j = b;
        }
        c.BasalJunction.Should().Be(junction);
        c.DroshaCut5Prime.Should().Be(i, "the 5' base of the 12th pair (Han 2006, pairs not nucleotides)");
        c.DroshaCut3Prime.Should().Be(j + Overhang, "partner + 2-nt 3' overhang");
        c.MatureStart.Should().Be(c.DroshaCut5Prime);
        (c.MatureEnd - c.MatureStart + 1).Should().Be(MatureLength, "Park 2011 5' counting");
        c.StarEnd.Should().Be(c.DroshaCut3Prime);
        c.MatureEnd.Should().BeLessThan(c.StarStart, "products on opposite arms");
        foreach (int coord in new[] { c.BasalJunction, c.DroshaCut5Prime, c.DroshaCut3Prime, c.MatureStart, c.MatureEnd, c.StarStart, c.StarEnd })
            coord.Should().BeInRange(0, len - 1);
        c.MatureSequence.Should().Be(c.Sequence.Substring(c.MatureStart, MatureLength));
        c.StarSequence.Should().Be(c.Sequence.Substring(c.StarStart, c.StarEnd - c.StarStart + 1));
        c.ThreePrimeOverhang.Should().Be(Overhang);
        c.Sequence.Should().NotContain("T");
    }

    // Real pri-miR-21: miRBase MI0000077 ± 30 nt hg38 flank; RNAfold 2.8 structure; junction 23 (11 bp below
    // the annotated 5p start). Expected: 5p [37,58] = hsa-miR-21-5p 22-nt form, 3p [76,96].
    private const string Mir21Pri =
        "CAUCGUGACAUCUCCAUGGCUGUACCACCUUGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACAUUUUGGUAUCUUUCAUCUGACCAUCCAUAU";
    private const string Mir21Fold =
        "..(((((.......)))))..((((((...((((((((((((((((.(((((.(((((.((((.((...))))))))))).)))))))))))))))))))))...)))))).....................";
    private const string HsaMir21_5p = "UAGCUUAUCAGACUGAUGUUGA";

    // Perfect hairpin: 5-nt flank, 35 G·C pairs, 8-nt loop, 25-nt flank (108 nt).
    private static readonly string Hairpin = new string('A', 5) + new string('G', 35) + new string('A', 8) + new string('C', 35) + new string('A', 25);
    private static readonly string HairpinFold = new string('.', 5) + new string('(', 35) + new string('.', 8) + new string(')', 35) + new string('.', 25);

    #endregion

    #region BE: empty / null sequence

    [Test]
    public void Predict_EmptyOrNull_ReturnsNullNeverThrows()
    {
        var actEmpty = () => PredictDroshaDicerCleavage("", 0);
        actEmpty.Should().NotThrow().Subject.Should().BeNull("empty ⇒ null (INV-CV-1)");
        var actNull = () => PredictDroshaDicerCleavage(null!, 0);
        actNull.Should().NotThrow().Subject.Should().BeNull("null ⇒ null (INV-CV-1)");
        PredictDroshaDicerCleavage("", 5).Should().BeNull();
        PredictDroshaDicerCleavage(null!, 0, null).Should().BeNull();
    }

    #endregion

    #region BE: precursor too short / stem too short

    /// <summary>
    /// A stem needs ≥ 12 pairs above the junction (≥ 24 paired nt + a 3-nt loop), so every homopolymer or
    /// random sequence of length 1..30 is null without throwing; on the perfect hairpin, junctions 5..28
    /// keep ≥ 12 pairs and opposite-arm products, junctions with ≤ 11 pairs above are null.
    /// </summary>
    [Test]
    public void Predict_PrecursorOrStemTooShort_ReturnsNullNeverIndexOverflow()
    {
        for (int len = 1; len <= 30; len++)
        {
            foreach (string s in new[] { new string('G', len), RandomRna(len, len) })
            {
                var act = () => PredictDroshaDicerCleavage(s, 0);
                act.Should().NotThrow($"len {len}").Subject.Should().BeNull($"len {len} cannot hold a 12-bp stem + 22-nt products");
            }
        }

        for (int j = 5; j < 40; j++)
        {
            var act = () => PredictDroshaDicerCleavage(Hairpin, j, HairpinFold);
            var r = act.Should().NotThrow($"junction {j}").Subject;
            if (r is not null) AssertCleavageContract(r.Value, j);
        }
        PredictDroshaDicerCleavage(Hairpin, 5, HairpinFold).Should().NotBeNull("35 pairs above junction 5");
        PredictDroshaDicerCleavage(Hairpin, 30, HairpinFold).Should().BeNull("only 10 pairs from junction 30");
    }

    #endregion

    #region BE: junction out of range / not the 5' base of a pair

    [Test]
    public void Predict_JunctionOutOfRangeOrUnpaired_ReturnsNullNeverThrows()
    {
        foreach (int j in new[] { -1, -11, int.MinValue, Hairpin.Length, Hairpin.Length + 1, int.MaxValue })
        {
            var act = () => PredictDroshaDicerCleavage(Hairpin, j, HairpinFold);
            act.Should().NotThrow($"junction {j}").Subject.Should().BeNull("out of range ⇒ null (INV-CV-1)");
        }
        foreach (int j in new[] { 0, 4, 42, 50, 60, 82, 107 }) // flank, loop, 3' arm, 3' flank
            PredictDroshaDicerCleavage(Hairpin, j, HairpinFold).Should().BeNull($"index {j} is not the 5' base of a pair");
    }

    #endregion

    #region MC: malformed structure

    [Test]
    public void Predict_MalformedStructure_ThrowsArgumentException()
    {
        foreach (string bad in new[] { "", "(", HairpinFold + ".", HairpinFold[1..], ")" + HairpinFold[1..], HairpinFold[..^1] + "(" })
        {
            var act = () => PredictDroshaDicerCleavage(Hairpin, 5, bad);
            act.Should().Throw<ArgumentException>($"structure of length {bad.Length} (INV-CV-4)");
        }
    }

    #endregion

    #region MC: non-ACGU / no-stem content

    [Test]
    public void Predict_NonAcguAndNoStem_NormalizedNeverCrash()
    {
        // DNA / lower-case spelling ⇒ identical record (T→U, upper) on the same structure.
        var rna = PredictDroshaDicerCleavage(Mir21Pri, 23, Mir21Fold)!.Value;
        var dna = PredictDroshaDicerCleavage(Mir21Pri.Replace('U', 'T').ToLowerInvariant(), 23, Mir21Fold)!.Value;
        dna.Should().Be(rna);

        foreach (string junk in new[]
                 {
                     new string('A', 80), new string('N', 80), new string('U', 80),
                     new string('N', 20) + "12345!@#$ \t%^&*" + new string('N', 20),
                     "GCGC" + new string('N', 40) + "GCGC",
                 })
        {
            foreach (int j in new[] { 0, 3, 20 })
            {
                var act = () => PredictDroshaDicerCleavage(junk, j);
                var r = act.Should().NotThrow($"junk len {junk.Length}, j {j}").Subject;
                if (r is not null) AssertCleavageContract(r.Value, j);
            }
        }
    }

    #endregion

    #region Discriminative pins

    [Test]
    public void Predict_Mir21PriMiRna_ReproducesMiRBaseMature()
    {
        var c = PredictDroshaDicerCleavage(Mir21Pri, 23, Mir21Fold)!.Value;
        c.DroshaCut5Prime.Should().Be(37);
        c.MatureSequence.Should().Be(HsaMir21_5p, "11 bp up the RNAfold stem + 22 nt = hsa-miR-21-5p");
        (c.StarStart, c.StarEnd).Should().Be((76, 96), "3p on the 3' arm (miRBase 3p = pri 75–97)");
        AssertCleavageContract(c, 23);
    }

    /// <summary>
    /// Random sequences folded by the canonical MFE with every 5'-paired base as junction: never throws, and
    /// every non-null record satisfies the contract re-derived on its own reported structure.
    /// </summary>
    [Test]
    public void Predict_RandomPrecursors_ContractHoldsAndCrashFree()
    {
        foreach (int seed in new[] { 1, 7, 42, 2026 })
        {
            foreach (int len in new[] { 60, 90, 120 })
            {
                string raw = RandomRna(len, seed);
                string fold = Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMfeStructure(raw).DotBracket;
                var p = Partners(fold);
                for (int j = 0; j < len; j++)
                {
                    var act = () => PredictDroshaDicerCleavage(raw, j, fold);
                    var r = act.Should().NotThrow($"len {len}, j {j}, seed {seed}").Subject;
                    if (p[j] <= j) r.Should().BeNull("junction not the 5' base of a pair");
                    if (r is null) continue;
                    AssertCleavageContract(r.Value, j);
                    PredictDroshaDicerCleavage(raw, j).Should().Be(r, "the 2-arg overload uses the same canonical fold");
                }
            }
        }
    }

    #endregion
}
