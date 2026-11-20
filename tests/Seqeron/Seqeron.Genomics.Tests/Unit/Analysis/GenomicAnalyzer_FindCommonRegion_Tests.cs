// GENOMIC-COMMON-001 — Common Region Detection (Longest Common Substring + all common regions)
// Evidence: docs/Evidence/GENOMIC-COMMON-001-Evidence.md
// TestSpec: tests/TestSpecs/GENOMIC-COMMON-001.md
// Source: Gusfield D. (1997). Algorithms on Strings, Trees and Sequences. Cambridge Univ. Press, ISBN 0-521-58519-8;
//         Wikipedia "Longest common substring"; GeeksforGeeks "Suffix Tree Application 5" — all accessed 2026-06-13.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
// Audit A3-15: runs alone - 0.5 GB peak working set (measured, LOCAL_RUN.md section 5).
[NonParallelizable]
public class GenomicAnalyzer_FindCommonRegion_Tests
{
    #region FindLongestCommonRegion

    // M1 — Definition (Wikipedia): longest CONTIGUOUS substring of both. ACGTACGT vs TTACGTGG -> TACGT (len 5).
    // Brute-force-verified (enumerate s[i:j] in t): maximal is TACGT at seq1 pos 3, seq2 pos 1.
    [Test]
    public void FindLongestCommonRegion_UniqueLcs_ReturnsTacgt()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("TTACGTGG");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(region.Sequence, Is.EqualTo("TACGT"),
                "Longest contiguous substring common to ACGTACGT and TTACGTGG is TACGT.");
            Assert.That(region.Length, Is.EqualTo(5), "TACGT has length 5; no length-6 common substring exists.");
            Assert.That(region.PositionInFirst, Is.EqualTo(3), "TACGT starts at 0-based index 3 in ACGTACGT.");
            Assert.That(region.PositionInSecond, Is.EqualTo(1), "TACGT starts at 0-based index 1 in TTACGTGG.");
        });
    }

    // M2 — Tie (Wikipedia: BADANAT/CANADAS share ADA & ANA). DNA analogue CACAGAG vs TACATAGAT shares
    // two distinct length-3 substrings ACA and AGA; documented tie-break returns the one first found in
    // sequence2 (TACATAGAT): ACA ends earlier than AGA -> ACA.
    [Test]
    public void FindLongestCommonRegion_Tie_ReturnsFirstInOther()
    {
        var seq1 = new DnaSequence("CACAGAG");
        var seq2 = new DnaSequence("TACATAGAT");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(region.Sequence, Is.EqualTo("ACA"),
                "Both ACA and AGA are maximal (len 3); tie-break returns the one first found in sequence2 -> ACA.");
            Assert.That(region.Length, Is.EqualTo(3), "Maximal common substring length is 3.");
            Assert.That(region.PositionInFirst, Is.EqualTo(1), "ACA starts at index 1 in CACAGAG.");
            Assert.That(region.PositionInSecond, Is.EqualTo(1), "ACA starts at index 1 in TACATAGAT.");
        });
    }

    // M3 — Identical sequences: a string is a substring of itself, so the whole sequence is the LCS at 0/0.
    [Test]
    public void FindLongestCommonRegion_Identical_ReturnsWholeSequence()
    {
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("ACGT");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(region.Sequence, Is.EqualTo("ACGT"), "Identical sequences share the whole sequence.");
            Assert.That(region.Length, Is.EqualTo(4), "Whole 4-base sequence is the LCS.");
            Assert.That(region.PositionInFirst, Is.EqualTo(0), "Starts at 0 in the first sequence.");
            Assert.That(region.PositionInSecond, Is.EqualTo(0), "Starts at 0 in the second sequence.");
        });
    }

    // M4 — No shared character: only the empty string qualifies -> length-0 LCS -> CommonRegion.None.
    [Test]
    public void FindLongestCommonRegion_NoCommon_ReturnsNone()
    {
        var seq1 = new DnaSequence("AAAA");
        var seq2 = new DnaSequence("GGGG");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(region.IsEmpty, Is.True, "No shared character -> empty LCS (length 0).");
            Assert.That(region.Sequence, Is.EqualTo(string.Empty), "Empty substring is returned.");
            Assert.That(region.Length, Is.EqualTo(0), "Length is 0.");
            Assert.That(region.PositionInFirst, Is.EqualTo(-1), "Sentinel position -1 when no region.");
            Assert.That(region.PositionInSecond, Is.EqualTo(-1), "Sentinel position -1 when no region.");
        });
    }

    // S1 — Empty first sequence -> CommonRegion.None (length-0 LCS).
    [Test]
    public void FindLongestCommonRegion_EmptyFirst_ReturnsNone()
    {
        var seq1 = new DnaSequence("");
        var seq2 = new DnaSequence("ACGT");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.That(region.IsEmpty, Is.True, "Empty first sequence shares nothing -> CommonRegion.None.");
    }

    // S2 — Empty second sequence -> CommonRegion.None.
    [Test]
    public void FindLongestCommonRegion_EmptySecond_ReturnsNone()
    {
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.That(region.IsEmpty, Is.True, "Empty second sequence shares nothing -> CommonRegion.None.");
    }

    // S4 — Single-character overlap: minimal non-empty LCS. "A" vs "TTTAT" -> A at seq1 pos 0, seq2 pos 3.
    [Test]
    public void FindLongestCommonRegion_SingleCharOverlap_ReturnsA()
    {
        var seq1 = new DnaSequence("A");
        var seq2 = new DnaSequence("TTTAT");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(region.Sequence, Is.EqualTo("A"), "Only the single base A is shared.");
            Assert.That(region.Length, Is.EqualTo(1), "Length-1 LCS.");
            Assert.That(region.PositionInFirst, Is.EqualTo(0), "A is at index 0 in the first sequence.");
            Assert.That(region.PositionInSecond, Is.EqualTo(3), "First A in TTTAT is at index 3.");
        });
    }

    // C1 — INV-01 property: the returned substring occurs in both sequences at the reported positions.
    [Test]
    public void FindLongestCommonRegion_ReturnedSubstring_OccursInBothAtPositions()
    {
        var seq1 = new DnaSequence("GATTACACGT");
        var seq2 = new DnaSequence("CCTTACAGG");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(region.IsEmpty, Is.False, "These sequences share TTACA.");
            Assert.That(seq1.Sequence.Substring(region.PositionInFirst, region.Length), Is.EqualTo(region.Sequence),
                "INV-01: substring must appear at PositionInFirst in sequence1.");
            Assert.That(seq2.Sequence.Substring(region.PositionInSecond, region.Length), Is.EqualTo(region.Sequence),
                "INV-01: substring must appear at PositionInSecond in sequence2.");
        });
    }

    #endregion

    #region FindCommonRegions

    // M6 — Right-maximal common regions of length >= 3 in ACGTACGT vs TTACGTGG.
    // FindCommonRegions returns, per start position in seq2, the SINGLE LONGEST substring also in seq1
    // (deduplicated), NOT every common substring. Hand-enumerating the longest match at each start of
    // TTACGTGG (positions 0..5) and deduplicating gives exactly {TACGT@(3,1), ACGT@(0,2), CGT@(1,3)}.
    // The FULL set of all common substrings >=3 would also include the prefixes TAC, TACG, ACG — those
    // are intentionally OMITTED (see FindCommonRegions_DoesNotReturnPrefixesOfLongerMatches below).
    [Test]
    public void FindCommonRegions_MinLengthThree_ReturnsThreeRegions()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("TTACGTGG");

        var regions = GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength: 3).ToList();
        var tuples = regions.Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();

        Assert.That(tuples, Is.EquivalentTo(new[]
        {
            ("TACGT", 3, 1),
            ("ACGT", 0, 2),
            ("CGT", 1, 3),
        }), "minLength=3 yields exactly the right-maximal set TACGT@(3,1), ACGT@(0,2), CGT@(1,3).");
    }

    // M5 — Same pair, minLength=4 drops CGT (length 3); leaves the right-maximal TACGT and ACGT.
    // Hand-verified: at each start of TTACGTGG, the longest substring also in ACGTACGT with length >=4 is
    // TACGT (start 1) and ACGT (start 2); all other starts have no match >=4.
    [Test]
    public void FindCommonRegions_MinLengthFour_ReturnsTwoRegions()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("TTACGTGG");

        var regions = GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength: 4).ToList();
        var tuples = regions.Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();

        Assert.That(tuples, Is.EquivalentTo(new[]
        {
            ("TACGT", 3, 1),
            ("ACGT", 0, 2),
        }), "minLength=4 yields exactly TACGT@(3,1) and ACGT@(0,2).");
    }

    // M7 — Contract lock: FindCommonRegions reports ONLY the longest match per start position in seq2
    // (right-maximal), NOT every common substring. For ACGTACGT vs TTACGTGG (minLength 3) the prefixes
    // TAC, TACG (prefixes of TACGT) and ACG (prefix of ACGT) are common substrings of length >= 3 but are
    // intentionally NOT returned. Hand-verified by full enumeration of all common substrings >= 3:
    // {TAC, TACG, TACGT, ACG, ACGT, CGT}; the returned set is the right-maximal subset {TACGT, ACGT, CGT}.
    [Test]
    public void FindCommonRegions_DoesNotReturnPrefixesOfLongerMatches()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("TTACGTGG");

        var returned = GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength: 3)
            .Select(r => r.Sequence)
            .ToList();

        Assert.Multiple(() =>
        {
            // These ARE common substrings of length >= 3, yet are deliberately omitted (not right-maximal).
            Assert.That(seq1.Sequence.Contains("TAC") && seq2.Sequence.Contains("TAC"), Is.True,
                "Sanity: TAC is a genuine common substring.");
            Assert.That(returned, Does.Not.Contain("TAC"), "Prefix TAC of TACGT must be omitted.");
            Assert.That(returned, Does.Not.Contain("TACG"), "Prefix TACG of TACGT must be omitted.");
            Assert.That(returned, Does.Not.Contain("ACG"), "Prefix ACG of ACGT must be omitted.");
            // The returned set is exactly the right-maximal one.
            Assert.That(returned, Is.EquivalentTo(new[] { "TACGT", "ACGT", "CGT" }),
                "Only the longest match per start position in seq2 is reported.");
        });
    }

    // S5 — minLength < 1 is treated as 1 (a region is non-empty). For ACGT vs ACGT the right-maximal
    // match at each start of seq2 is ACGT@(0,0), CGT@(1,1), GT@(2,2), T@(3,3). Hand-verified; minLength=0
    // and minLength=-5 must give the same result as minLength=1 (no empty/length-0 regions).
    [Test]
    public void FindCommonRegions_MinLengthBelowOne_TreatedAsOne()
    {
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("ACGT");

        var expected = new[] { ("ACGT", 0, 0), ("CGT", 1, 1), ("GT", 2, 2), ("T", 3, 3) };

        var atZero = GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength: 0)
            .Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();
        var atNeg = GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength: -5)
            .Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(atZero, Is.EquivalentTo(expected), "minLength=0 is treated as 1; no empty regions.");
            Assert.That(atNeg, Is.EquivalentTo(expected), "minLength<0 is treated as 1; no empty regions.");
            Assert.That(atZero.All(t => t.Sequence.Length >= 1), Is.True,
                "No length-0 / empty region is ever returned.");
        });
    }

    // S3 — No shared substring -> empty enumeration.
    [Test]
    public void FindCommonRegions_NoMatch_ReturnsEmpty()
    {
        var seq1 = new DnaSequence("AAAA");
        var seq2 = new DnaSequence("GGGG");

        var regions = GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength: 2).ToList();

        Assert.That(regions, Is.Empty, "No shared base -> no common regions.");
    }

    // C2 — INV-04 property (O(n*m) algorithm): every region has length >= minLength and is contained in sequence1.
    [Test]
    public void FindCommonRegions_AllRegions_SatisfyMinLengthAndContainment()
    {
        var seq1 = new DnaSequence("GATTACAGATTACA");
        var seq2 = new DnaSequence("TTACAGGGATTAC");
        const int MinLength = 3;

        var regions = GenomicAnalyzer.FindCommonRegions(seq1, seq2, MinLength).ToList();

        Assert.That(regions, Is.Not.Empty, "Sequences share several substrings of length >= 3.");
        Assert.Multiple(() =>
        {
            foreach (var r in regions)
            {
                Assert.That(r.Length, Is.GreaterThanOrEqualTo(MinLength),
                    $"INV-04: region '{r.Sequence}' must be at least minLength={MinLength}.");
                Assert.That(seq1.Sequence.Contains(r.Sequence), Is.True,
                    $"INV-04: region '{r.Sequence}' must be a substring of sequence1.");
                Assert.That(seq2.Sequence.Substring(r.PositionInSecond, r.Length), Is.EqualTo(r.Sequence),
                    $"INV-01: region '{r.Sequence}' must appear at PositionInSecond in sequence2.");
            }
        });
    }

    #endregion

    #region B09 review: first-occurrence positions, null contract, oracle, scalability

    // F15 — PositionInFirst must be the FIRST (lowest) occurrence in sequence1, independent of the
    // suffix-tree's internal (DFS) occurrence order. Reference values: Python str.find (independent oracle,
    // scripts/analysis/generate_common_region_oracle.py).
    [Test]
    public void FindCommonRegions_PositionInFirst_IsLowestOccurrence()
    {
        // The final region GA occurs in seq1 at 2 and 15; the suffix tree lists 15 first (DFS order), the first occurrence is 2.
        var seq1 = new DnaSequence("AGGAAGAAGGAGGCCGATTACT");
        var seq2 = new DnaSequence("CGTGGAAGGAGGCCGA");

        var regions = GenomicAnalyzer.FindCommonRegions(seq1, seq2, minLength: 2).ToList();

        Assert.That(regions.Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)), Is.EqualTo(new[]
        {
            ("CG", 14, 0), ("GGAAG", 1, 3), ("GAAGGAGGCCGA", 5, 4), ("AAGGAGGCCGA", 6, 5), ("AGGAGGCCGA", 7, 6),
            ("GGAGGCCGA", 8, 7), ("GAGGCCGA", 9, 8), ("AGGCCGA", 10, 9), ("GGCCGA", 11, 10), ("GCCGA", 12, 11),
            ("CCGA", 13, 12), ("CGA", 14, 13), ("GA", 2, 14),
        }));
    }

    [Test]
    public void FindLongestCommonRegion_PositionInFirst_IsLowestOccurrence()
    {
        // LCS "GA" (length 2) occurs at 4, 15 (and others) in seq1; str.find = 4.
        var seq1 = new DnaSequence("CACAGATAAAGGACTGAGGCTCTACG");
        var seq2 = new DnaSequence("GA");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.That((region.Sequence, region.PositionInFirst, region.PositionInSecond), Is.EqualTo(("GA", 4, 0)));
    }

    // Tie: first-found-in-sequence2 representative; both positions are first occurrences (poly-repeat seq1).
    [Test]
    public void FindLongestCommonRegion_RepeatedOccurrences_ReportsFirstInBoth()
    {
        // ACAC occurs in seq1 at 0, 2, 4, 6 and in seq2 only at 2 (Python str.find: 0 and 2).
        var seq1 = new DnaSequence("ACACACACAC");
        var seq2 = new DnaSequence("GGACACGG");

        CommonRegion region = GenomicAnalyzer.FindLongestCommonRegion(seq1, seq2);

        Assert.That((region.Sequence, region.PositionInFirst, region.PositionInSecond), Is.EqualTo(("ACAC", 0, 2)));
    }

    [Test]
    public void CommonRegionMethods_NullSequences_ThrowArgumentNullException()
    {
        var seq = new DnaSequence("ACGT");

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindLongestCommonRegion(null!, seq));
            Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindLongestCommonRegion(seq, null!));
            // Eager validation: the exception is thrown by the call itself, not on first enumeration.
            Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindCommonRegions(null!, seq, 2));
            Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindCommonRegions(seq, null!, 2));
        });
    }

    // Lower-case input is folded by DnaSequence, so matching is case-insensitive (ties/positions unchanged).
    [Test]
    public void CommonRegionMethods_LowerCaseInput_SameAsUpperCase()
    {
        var upper1 = new DnaSequence("ACGTACGT"); var upper2 = new DnaSequence("TTACGTGG");
        var lower1 = new DnaSequence("acgtacgt"); var lower2 = new DnaSequence("ttacgtgg");

        Assert.Multiple(() =>
        {
            Assert.That(GenomicAnalyzer.FindLongestCommonRegion(lower1, lower2),
                Is.EqualTo(GenomicAnalyzer.FindLongestCommonRegion(upper1, upper2)));
            Assert.That(GenomicAnalyzer.FindCommonRegions(lower1, lower2, 3),
                Is.EqualTo(GenomicAnalyzer.FindCommonRegions(upper1, upper2, 3)));
        });
    }

    // Non-ACGT symbols (IUPAC N, RNA U) never reach the algorithm: DnaSequence rejects them.
    [TestCase("ACGN")]
    [TestCase("ACGU")]
    public void CommonRegionInput_NonAcgt_RejectedByDnaSequence(string text)
    {
        Assert.Throws<ArgumentException>(() => new DnaSequence(text));
    }

    // minLength = int.MinValue / 0 / 1 all mean "at least one base"; huge minLength yields nothing.
    [Test]
    public void FindCommonRegions_ExtremeMinLength_Handled()
    {
        var a = new DnaSequence("ACGTACGT"); var b = new DnaSequence("TTACGTGG");
        var atOne = GenomicAnalyzer.FindCommonRegions(a, b, 1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(GenomicAnalyzer.FindCommonRegions(a, b, int.MinValue), Is.EqualTo(atOne));
            Assert.That(GenomicAnalyzer.FindCommonRegions(a, b, 0), Is.EqualTo(atOne));
            Assert.That(GenomicAnalyzer.FindCommonRegions(a, b, int.MaxValue), Is.Empty);
            Assert.That(GenomicAnalyzer.FindCommonRegions(a, b, 6), Is.Empty, "longest common substring is 5");
            Assert.That(GenomicAnalyzer.FindCommonRegions(a, b, 5).Select(r => r.Sequence), Is.EqualTo(new[] { "TACGT" }));
        });
    }

    // Differential test against the independent Python oracle (DP + suffix automaton for the LCS,
    // brute-force matching statistics for the regions): 284 pairs incl. empty inputs, binary alphabet, tandem
    // repeats, planted (point-mutated) shared blocks up to 250 bp that exercise the maximal-exact-match path.
    [Test]
    public void FindLongestCommonRegionAndFindCommonRegions_MatchIndependentOracle()
    {
        var rows = ReadOracle();
        var failures = new List<string>();
        foreach (var (id, s1, s2, minLength, lcs, regions) in rows)
        {
            var a = new DnaSequence(s1); var b = new DnaSequence(s2);

            var l = GenomicAnalyzer.FindLongestCommonRegion(a, b);
            string gotLcs = l.IsEmpty ? "-" : $"{l.Sequence}@{l.PositionInFirst}@{l.PositionInSecond}";
            if (gotLcs != lcs) failures.Add($"#{id} LCS {gotLcs} != {lcs}");

            var r = GenomicAnalyzer.FindCommonRegions(a, b, minLength).ToList();
            string gotRegions = r.Count == 0 ? "-" : string.Join(",", r.Select(x => $"{x.Length}:{x.PositionInFirst}:{x.PositionInSecond}"));
            bool sequencesOk = r.All(x => s2.Substring(x.PositionInSecond, x.Length) == x.Sequence
                                          && s1.Substring(x.PositionInFirst, x.Length) == x.Sequence);
            if (gotRegions != regions || !sequencesOk) failures.Add($"#{id} regions(min={minLength}) {gotRegions} != {regions} (seq ok: {sequencesOk})");
        }

        Assert.That(rows, Has.Count.EqualTo(283));
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    // Matching statistics from the oracle's own definition, checked on fresh random pairs through a brute-force
    // C# reference (independent of the suffix tree): longest-common-substring length, region set and positions.
    [Test]
    public void FindCommonRegions_RandomPairs_MatchBruteForce()
    {
        var rnd = new Random(20261118);
        for (int k = 0; k < 150; k++)
        {
            string alphabet = k % 3 == 0 ? "AC" : "ACGT";
            string s = RandomText(rnd, rnd.Next(0, 60), alphabet);
            string t = RandomText(rnd, rnd.Next(0, 60), alphabet);
            if (k % 2 == 0 && s.Length > 10) t = t.Insert(rnd.Next(t.Length + 1), s.Substring(rnd.Next(s.Length - 8), 8));
            int min = rnd.Next(0, 9);

            var expected = BruteForceRegions(s, t, min);
            var actual = GenomicAnalyzer.FindCommonRegions(new DnaSequence(s), new DnaSequence(t), min)
                .Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();

            Assert.That(actual, Is.EqualTo(expected), $"s={s} t={t} min={min}");

            var lcs = GenomicAnalyzer.FindLongestCommonRegion(new DnaSequence(s), new DnaSequence(t));
            int bestLength = BruteForceRegions(s, t, 1).Select(x => x.Item1.Length).DefaultIfEmpty(0).Max();
            Assert.That(lcs.Length, Is.EqualTo(bestLength), $"LCS length s={s} t={t}");
        }
    }

    // Complexity guard (INV: no quadratic time/memory). A 1 Mbp random pair; the former per-start binary
    // search over Contains copied up to m characters per probe (~O(m^2): minutes and tens of GB for 1 Mbp).
    [Test]
    public void FindCommonRegions_OneMegabasePair_IsLinearAndCorrect()
    {
        var rnd = new Random(7);
        const int N = 1_000_000;
        string s = RandomText(rnd, N, "ACGT");
        string t = RandomText(rnd, N, "ACGT");
        // Plant a 5 kb shared block (point-mutated in the middle) so a long match must be found.
        string block = s.Substring(250_000, 5_000);
        t = t.Remove(600_000, 5_000).Insert(600_000, block);
        var a = new DnaSequence(s); var b = new DnaSequence(t);

        var sw = System.Diagnostics.Stopwatch.StartNew();
        var lcs = GenomicAnalyzer.FindLongestCommonRegion(a, b);
        var big = GenomicAnalyzer.FindCommonRegions(a, b, 1000).ToList();
        var small = GenomicAnalyzer.FindCommonRegions(a, b, 1).Take(200).ToList();
        sw.Stop();

        Assert.Multiple(() =>
        {
            Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(60), "linear-time behaviour expected");
            Assert.That(lcs.Length, Is.GreaterThanOrEqualTo(5_000));
            Assert.That(t.Substring(lcs.PositionInSecond, lcs.Length), Is.EqualTo(lcs.Sequence));
            Assert.That(string.CompareOrdinal(s, lcs.PositionInFirst, lcs.Sequence, 0, lcs.Length), Is.EqualTo(0));
            // Regions >= 1000: the planted block (t[600000..605000) = s[250000..255000)) and its right-maximal suffix starts.
            Assert.That(big[0].PositionInSecond, Is.LessThanOrEqualTo(600_000));
            Assert.That(big.All(r => r.Length >= 1000), Is.True);
            Assert.That(big.First(r => r.PositionInSecond == 600_000).Length, Is.GreaterThanOrEqualTo(5_000));
            foreach (var r in small)
            {
                Assert.That(string.CompareOrdinal(t, r.PositionInSecond, r.Sequence, 0, r.Length), Is.EqualTo(0));
                Assert.That(string.CompareOrdinal(s, r.PositionInFirst, r.Sequence, 0, r.Length), Is.EqualTo(0));
                Assert.That(s.IndexOf(r.Sequence, StringComparison.Ordinal), Is.EqualTo(r.PositionInFirst), "first occurrence");
            }
        });
    }

    private static string RandomText(Random rnd, int length, string alphabet)
    {
        var c = new char[length];
        for (int i = 0; i < length; i++) c[i] = alphabet[rnd.Next(alphabet.Length)];
        return new string(c);
    }

    // Definition-level reference: ms(i) by brute-force substring search, distinct strings by first start.
    private static List<(string, int, int)> BruteForceRegions(string s, string t, int minLength)
    {
        int min = Math.Max(1, minLength);
        var seen = new HashSet<string>();
        var result = new List<(string, int, int)>();
        for (int i = 0; i < t.Length; i++)
        {
            int l = 0;
            while (i + l < t.Length && s.Contains(t.AsSpan(i, l + 1), StringComparison.Ordinal)) l++;
            if (l >= min && seen.Add(t.Substring(i, l)))
                result.Add((t.Substring(i, l), s.IndexOf(t.Substring(i, l), StringComparison.Ordinal), i));
        }
        return result;
    }

    private static List<(int Id, string S1, string S2, int Min, string Lcs, string Regions)> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.common_region_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<(int, string, string, int, string, string)>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#') continue;
            var f = line.Split('\t');
            rows.Add((int.Parse(f[0]), f[1] == "." ? "" : f[1], f[2] == "." ? "" : f[2], int.Parse(f[3]), f[4], f[5]));
        }
        return rows;
    }

    #endregion
}
