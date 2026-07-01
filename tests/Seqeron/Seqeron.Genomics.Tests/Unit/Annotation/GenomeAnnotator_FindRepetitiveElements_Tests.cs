// ANNOT-REPEAT-001 — Repetitive Element Detection and Classification
// Evidence: docs/Evidence/ANNOT-REPEAT-001-Evidence.md
// TestSpec: tests/TestSpecs/ANNOT-REPEAT-001.md
// Source: Wikipedia "Tandem repeat" (cites Duitama et al. 2014); Wikipedia "Inverted repeat";
//         Hampson et al. (2021) IUPACpal, BMC Bioinformatics 22:51, doi:10.1186/s12859-021-03983-2;
//         RepeatMasker documentation (Smit/Hubley/Green; Repbase).

namespace Seqeron.Genomics.Tests.Unit.Annotation;

using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Annotation;

[TestFixture]
public class GenomeAnnotator_FindRepetitiveElements_Tests
{
    #region FindRepetitiveElements

    // M1 — Tandem repeat: "ATTCGATTCGATTCG" = unit "ATTCG" x3, head-to-tail (Wikipedia worked example).
    [Test]
    public void FindRepetitiveElements_ThreeAdjacentCopies_ReturnsSingleTandemRepeatSpanningWhole()
    {
        var elements = GenomeAnnotator.FindRepetitiveElements("ATTCGATTCGATTCG", minRepeatLength: 5, minCopies: 2)
            .Where(e => e.type == "tandem_repeat")
            .ToList();

        var tr = elements.SingleOrDefault(e => e.start == 0 && e.end == 15);
        Assert.Multiple(() =>
        {
            Assert.That(tr.type, Is.EqualTo("tandem_repeat"),
                "A head-to-tail repeat of ATTCG x3 must be reported as a tandem_repeat.");
            Assert.That(tr.start, Is.EqualTo(0), "Tandem array starts at index 0.");
            Assert.That(tr.end, Is.EqualTo(15), "End is exclusive; 3 copies x 5 bp span [0,15).");
            Assert.That(tr.sequence, Is.EqualTo("ATTCGATTCGATTCG"),
                "Reported sequence must equal the array slice (INV-01).");
        });
    }

    // M2 — Single occurrence of a motif is NOT a tandem repeat (Wikipedia: "two or more").
    [Test]
    public void FindRepetitiveElements_MotifAppearsOnce_NoTandemRepeatForThatMotif()
    {
        // "ATTCG" appears once; the only run is AAAAA (unit A x5).
        var tandem = GenomeAnnotator.FindRepetitiveElements("ATTCGAAAAA", minRepeatLength: 5, minCopies: 2)
            .Where(e => e.type == "tandem_repeat")
            .ToList();

        Assert.That(tandem.Any(e => e.sequence.Contains("ATTCGATTCG")), Is.False,
            "A motif occurring only once must not be reported as a tandem repeat.");
        Assert.That(tandem.Any(e => e.start == 0 && e.sequence.StartsWith("ATTCG", StringComparison.Ordinal)
                                    && e.end - e.start >= 10), Is.False,
            "No tandem array of ATTCG exists (it appears a single time).");
    }

    // M3 — Inverted repeat, gap 0: "GAATTC" -> left "GAA", right "TTC" = revcomp(GAA) (Wikipedia + IUPACpal WW^R).
    [Test]
    public void FindRepetitiveElements_ReverseComplementPalindrome_ReturnsInvertedRepeat()
    {
        var inverted = GenomeAnnotator.FindRepetitiveElements("GAATTC", minRepeatLength: 3, minCopies: 2)
            .Where(e => e.type == "inverted_repeat")
            .ToList();

        var ir = inverted.SingleOrDefault(e => e.start == 0 && e.end == 6);
        Assert.Multiple(() =>
        {
            Assert.That(ir.type, Is.EqualTo("inverted_repeat"),
                "GAATTC is a reverse-complement palindrome (arms GAA / TTC, gap 0).");
            Assert.That(ir.start, Is.EqualTo(0), "Left arm starts at index 0.");
            Assert.That(ir.end, Is.EqualTo(6), "Right arm ends at exclusive index 6.");
            Assert.That(ir.sequence, Is.EqualTo("GAATTC"),
                "Gap-0 inverted repeat composite is arm1+arm2 with no bracketed gap.");
        });
    }

    // M4 — Gapped inverted repeat: "TTACGAAAAAACGTAA": left "TTACG", gap "AAAAAA" (6), right "CGTAA"=revcomp (WGW^R).
    [Test]
    public void FindRepetitiveElements_GappedInvertedRepeat_ReturnsArmsWithBracketedGap()
    {
        var inverted = GenomeAnnotator.FindRepetitiveElements("TTACGAAAAAACGTAA", minRepeatLength: 5, minCopies: 2)
            .Where(e => e.type == "inverted_repeat" && e.start == 0 && e.end == 16)
            .ToList();

        var ir = inverted.SingleOrDefault();
        Assert.Multiple(() =>
        {
            Assert.That(ir.type, Is.EqualTo("inverted_repeat"),
                "Left arm TTACG and right arm CGTAA=revcomp(TTACG) form a gapped inverted repeat.");
            Assert.That(ir.start, Is.EqualTo(0), "Left arm starts at index 0.");
            Assert.That(ir.end, Is.EqualTo(16), "Right arm ends at exclusive index 16.");
            Assert.That(ir.sequence, Is.EqualTo("TTACG[AAAAAA]CGTAA"),
                "Composite must embed the 6-nt gap in brackets: arm1[gap]arm2.");
        });
    }

    // INV-3 — For every reported inverted repeat, the right arm equals the reverse complement of the
    //         left arm (Wikipedia "Inverted repeat"; IUPACpal WGW^R, Hampson et al. 2021).
    //         Hand-checked: revcomp(GAA)=TTC, revcomp(TTACG)=CGTAA, revcomp(GGATCC)=GGATCC (palindrome).
    [Test]
    public void FindRepetitiveElements_InvertedRepeats_RightArmIsReverseComplementOfLeftArm()
    {
        // (sequence, minArm, expectedLeftArm, expectedRightArm)
        var cases = new[]
        {
            ("GAATTC", 3, "GAA", "TTC"),                 // EcoRI palindrome, gap 0
            ("TTACGAAAAAACGTAA", 5, "TTACG", "CGTAA"),   // gapped, gap 6 (Wikipedia example)
            ("GGATCC", 3, "GGA", "TCC"),                 // BamHI palindrome, gap 0
        };

        foreach (var (seq, minArm, leftArm, rightArm) in cases)
        {
            int armLen = leftArm.Length;
            var ir = GenomeAnnotator.FindRepetitiveElements(seq, minRepeatLength: minArm, minCopies: 2)
                .Where(e => e.type == "inverted_repeat" && e.start == 0)
                .OrderBy(e => e.end - e.start)
                .First(e => e.end == seq.Length);

            string actualLeft = seq.Substring(0, armLen);
            string actualRight = seq.Substring(seq.Length - armLen, armLen);
            string revComp = ReverseComplement(actualLeft);

            Assert.Multiple(() =>
            {
                Assert.That(actualLeft, Is.EqualTo(leftArm), $"Left arm of '{seq}'.");
                Assert.That(actualRight, Is.EqualTo(rightArm), $"Right arm of '{seq}'.");
                Assert.That(actualRight, Is.EqualTo(revComp),
                    $"INV-3: right arm must equal reverse complement of left arm for '{seq}'.");
            });
        }
    }

    // Test helper: A<->T, C<->G complement then reverse (independent of library code under test).
    private static string ReverseComplement(string s)
    {
        var chars = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
        {
            char c = s[s.Length - 1 - i];
            chars[i] = c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c };
        }
        return new string(chars);
    }

    // S1 — Primitive unit preference: "AAAAAA" reported as mononucleotide A, not AA/AAA (Wikipedia non-primitive).
    [Test]
    public void FindRepetitiveElements_HomopolymerRun_ReportedAsPrimitiveMononucleotide()
    {
        var tandem = GenomeAnnotator.FindRepetitiveElements("AAAAAA", minRepeatLength: 4, minCopies: 2)
            .Where(e => e.type == "tandem_repeat")
            .ToList();

        var array = tandem.SingleOrDefault(e => e.start == 0 && e.end == 6);
        Assert.Multiple(() =>
        {
            Assert.That(array.sequence, Is.EqualTo("AAAAAA"),
                "The maximal homopolymer array spans the whole input.");
            // Only the primitive period (A) array is reported once; no AA/AAA duplicate spans.
            Assert.That(tandem.Count(e => e.start == 0 && e.end == 6), Is.EqualTo(1),
                "Non-primitive units (AA, AAA) must collapse to the primitive A; reported once (de-dup).");
        });
    }

    // S2 — null sequence throws (input contract, mirrors sibling analyzer methods).
    [Test]
    public void FindRepetitiveElements_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() => GenomeAnnotator.FindRepetitiveElements(null!).ToList(),
            "A null sequence must raise ArgumentNullException.");
    }

    // S2b — minCopies < 2 throws: a tandem repeat requires two or more copies (Wikipedia).
    [Test]
    public void FindRepetitiveElements_MinCopiesBelowTwo_ThrowsArgumentOutOfRange()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => GenomeAnnotator.FindRepetitiveElements("ACGT", minRepeatLength: 1, minCopies: 1).ToList(),
            "minCopies < 2 is invalid: a tandem repeat needs at least two adjacent copies.");
    }

    // S3 — empty sequence yields empty result (no throw).
    [Test]
    public void FindRepetitiveElements_EmptySequence_ReturnsEmpty()
    {
        var elements = GenomeAnnotator.FindRepetitiveElements("").ToList();

        Assert.That(elements, Is.Empty, "An empty sequence has no repetitive elements.");
    }

    // C1 — INV-01/INV-02 property: every tandem repeat span equals an integer number of >= minCopies copies
    //      and its reported sequence equals the slice. Deterministic inputs (no randomness).
    [Test]
    public void FindRepetitiveElements_TandemRepeats_SatisfyStructuralInvariants()
    {
        const int minCopies = 2;
        string[] inputs = { "CACACACAGT", "ATATATATAT", "GGGGGGTTAA", "ACGTACGTACGTNN", "AAAAAACCCCCC" };

        foreach (var seq in inputs)
        {
            var tandem = GenomeAnnotator.FindRepetitiveElements(seq, minRepeatLength: 4, minCopies: minCopies)
                .Where(e => e.type == "tandem_repeat");

            foreach (var (start, end, _, sequence) in tandem)
            {
                Assert.Multiple(() =>
                {
                    Assert.That(sequence, Is.EqualTo(seq.ToUpperInvariant().Substring(start, end - start)),
                        $"INV-01: reported sequence must equal the slice [{start},{end}) of '{seq}'.");
                    // INV-02: span is an integer number of >= minCopies copies of a primitive unit.
                    int span = end - start;
                    int unitLen = FindPrimitiveUnitLength(sequence);
                    Assert.That(span % unitLen, Is.EqualTo(0),
                        "INV-01: span must be an integer multiple of the unit length.");
                    Assert.That(span / unitLen, Is.GreaterThanOrEqualTo(minCopies),
                        "INV-02: array must contain at least minCopies copies.");
                });
            }
        }
    }

    // Test helper: smallest period of a perfect tandem array string.
    private static int FindPrimitiveUnitLength(string array)
    {
        int n = array.Length;
        for (int period = 1; period <= n; period++)
        {
            if (n % period != 0) continue;
            bool ok = true;
            for (int i = period; i < n && ok; i++)
                if (array[i] != array[i % period]) ok = false;
            if (ok) return period;
        }
        return n;
    }

    #endregion

    #region ClassifyRepeat

    // M5 — Classify by library best-match: query contains an Alu SINE library entry -> "SINE/Alu" (RepeatMasker).
    [Test]
    public void ClassifyRepeat_QueryContainsLibraryElement_ReturnsLibraryClass()
    {
        var db = new Dictionary<string, string>
        {
            ["GGCCGGGCGCGGTGGCTCAC"] = "SINE/Alu",
            ["TTAGGG"] = "Satellite",
        };

        string cls = GenomeAnnotator.ClassifyRepeat("CCCCGGCCGGGCGCGGTGGCTCACGGGG", db);

        Assert.That(cls, Is.EqualTo("SINE/Alu"),
            "The class of the matching library element (best/longest match) must be assigned.");
    }

    // M6 — Classify fallback to simple repeat: dinucleotide CA repeat absent from library -> "Simple_repeat".
    [Test]
    public void ClassifyRepeat_NoLibraryMatchButSimpleMotif_ReturnsSimpleRepeat()
    {
        var db = new Dictionary<string, string> { ["GGGGGGGGGGGGGGGGGGGG"] = "Satellite" };

        string cls = GenomeAnnotator.ClassifyRepeat("CACACACA", db);

        Assert.That(cls, Is.EqualTo("Simple_repeat"),
            "An unmatched 1-6 bp tandem motif (CA dinucleotide) falls back to Simple_repeat (STR).");
    }

    // M7 — Classify with no match and non-simple sequence -> "Unknown" (RepeatMasker Unclassified).
    [Test]
    public void ClassifyRepeat_NoLibraryMatchAndNotSimple_ReturnsUnknown()
    {
        var db = new Dictionary<string, string> { ["GGCCGGGCGCGGTGGCTCAC"] = "SINE/Alu" };

        string cls = GenomeAnnotator.ClassifyRepeat("ACGTACTGGATCCAGTTGCAC", db);

        Assert.That(cls, Is.EqualTo("Unknown"),
            "A non-matching, non-simple sequence is unclassified (Unknown).");
    }

    // M7b — A trivially short query must NOT be classified as a library class merely because
    //       a longer consensus happens to contain those letters. RepeatMasker screens the
    //       query for occurrences of known elements (element within query), not the reverse;
    //       a single base "A" is not a SINE just because an Alu consensus contains an A.
    //       (Defect found in validation: bidirectional containment misclassified short queries.)
    [Test]
    public void ClassifyRepeat_SingleBaseSubstringOfConsensus_DoesNotInheritLibraryClass()
    {
        // Both library consensi contain the letter 'A'; neither equals nor is contained in "A".
        var db = new Dictionary<string, string>
        {
            ["GGCCGGGCGCGGTGGCTCAC"] = "SINE/Alu", // contains an 'A'
            ["TTAGGG"] = "Satellite",              // contains an 'A'
        };

        // "A" is a 1-bp homopolymer: not a library element occurrence, but a simple mononucleotide repeat
        // is only recognised with >=2 copies, so a lone base is unclassified (Unknown), never SINE/Alu.
        string cls = GenomeAnnotator.ClassifyRepeat("A", db);

        Assert.That(cls, Is.EqualTo("Unknown"),
            "A 1-bp query is not an occurrence of any library element; it must not inherit a library class.");
    }

    // M7c — A short query that is a *substring of* a longer library consensus (but does not itself
    //       contain that consensus) is not a library-element occurrence. RepeatMasker reports a
    //       known element found *within* the query, so a 7-bp fragment of a 20-bp Alu consensus,
    //       containing no full library element, is not classified as that family by containment.
    [Test]
    public void ClassifyRepeat_QueryIsSubstringOfConsensus_NotClassifiedByReverseContainment()
    {
        var db = new Dictionary<string, string> { ["GGCCGGGCGCGGTGGCTCAC"] = "SINE/Alu" };

        // "GGCCGGG" is a substring of the consensus but contains no full library element and is not simple.
        string cls = GenomeAnnotator.ClassifyRepeat("GGCCGGG", db);

        Assert.That(cls, Is.EqualTo("Unknown"),
            "A fragment that is contained in a consensus (but does not contain a full library element) is Unknown, not SINE/Alu.");
    }

    // M5b — Best (longest) match wins when several library elements occur in the query (RepeatMasker
    //       reports the best match). A query harbouring both a 6-bp Satellite motif and a 20-bp Alu
    //       consensus is classified by the longer Alu element.
    [Test]
    public void ClassifyRepeat_MultipleElementsOccur_ReturnsLongestMatch()
    {
        var db = new Dictionary<string, string>
        {
            ["GGCCGGGCGCGGTGGCTCAC"] = "SINE/Alu", // 20 bp
            ["TTAGGG"] = "Satellite",              // 6 bp
        };

        // Query contains both elements; the 20-bp Alu is the better (longer) match.
        string cls = GenomeAnnotator.ClassifyRepeat("TTAGGGAAGGCCGGGCGCGGTGGCTCACTT", db);

        Assert.That(cls, Is.EqualTo("SINE/Alu"),
            "When several library elements occur, the longest (best) match's class is assigned.");
    }

    // M6b — Empty query is a degenerate simple-repeat fallback (no library element can occur in it).
    [Test]
    public void ClassifyRepeat_EmptySequence_ReturnsSimpleRepeat()
    {
        var db = new Dictionary<string, string> { ["GGCCGGGCGCGGTGGCTCAC"] = "SINE/Alu" };

        string cls = GenomeAnnotator.ClassifyRepeat("", db);

        Assert.That(cls, Is.EqualTo("Simple_repeat"),
            "An empty query has no library occurrence; the documented fallback is Simple_repeat.");
    }

    // S4 — null arguments throw (input contract).
    [Test]
    public void ClassifyRepeat_NullArguments_ThrowArgumentNullException()
    {
        var db = new Dictionary<string, string>();
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => GenomeAnnotator.ClassifyRepeat(null!, db),
                "null sequence must throw.");
            Assert.Throws<ArgumentNullException>(() => GenomeAnnotator.ClassifyRepeat("ACGT", null!),
                "null repeatDb must throw.");
        });
    }

    #endregion
}
