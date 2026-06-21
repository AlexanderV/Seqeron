namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Analysis area (GenomicAnalyzer).
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of combinatorial testing.
/// Each grid cell carries a real business assertion; small grids use the exhaustive
/// <c>[Combinatorial]</c> product (a strict superset of pairwise).
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Analysis")]
public class AnalysisCombinatorialTests
{
    /// <summary>Deterministic well-mixed ACGT sequence (LCG).</summary>
    private static string DiverseDna(int n, uint seed)
    {
        const string bases = "ACGT";
        var chars = new char[n];
        uint state = seed;
        for (int i = 0; i < n; i++)
        {
            state = state * 1664525u + 1013904223u;
            chars[i] = bases[(int)((state >> 16) & 3u)];
        }
        return new string(chars);
    }

    private static int CountOverlapping(string text, string sub)
    {
        int count = 0, idx = 0;
        while ((idx = text.IndexOf(sub, idx, StringComparison.Ordinal)) >= 0) { count++; idx++; }
        return count;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: GENOMIC-COMMON-001 — Common regions between two sequences (Analysis)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 175.
    // Spec: tests/TestSpecs/GENOMIC-COMMON-001.md (canonical FindLongestCommonRegion / FindCommonRegions).
    // ADVANCED §10.
    // Dimensions: nSeqs(3) × minLen(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Gusfield 1997; Wikipedia "Longest common substring"): the longest common substring (LCS)
    // of two strings is found via a generalized suffix tree; FindCommonRegions reports, per start
    // position in seq2, the single longest substring of length ≥ minLength that also occurs in seq1.
    //
    // Axis mapping (documented — these are 2-sequence methods): nSeqs → the planted common-block
    // LENGTH class (6/9/12); minLen → the minLength threshold. Engineered construct: a unique block W
    // shared by seq1 (poly-C flanks) and seq2 (poly-A flanks), so W is the UNIQUE LCS (no incidental
    // common run exceeds 1). The combinatorial point: the LCS equals W at every cell, while the block
    // appears in FindCommonRegions exactly when minLen ≤ block length.
    // ═══════════════════════════════════════════════════════════════════════

    private const string CommonBlock = "GATTACAGGCAT"; // 12 nt; prefixes give the 6/9/12 length classes

    private static int BruteLcsLength(string a, string b)
    {
        int best = 0;
        for (int i = 0; i < a.Length; i++)
            for (int j = 0; j < b.Length; j++)
            {
                int len = 0;
                while (i + len < a.Length && j + len < b.Length && a[i + len] == b[j + len]) len++;
                if (len > best) best = len;
            }
        return best;
    }

    [Test, Combinatorial]
    public void GenomicCommon_LcsAndRegions_AcrossBlockLengthAndMinLen(
        [Values(6, 9, 12)] int blockLen,
        [Values(3, 6, 9)] int minLen)
    {
        string w = CommonBlock[..blockLen];
        string s1 = new string('C', 6) + w + new string('C', 6);
        string s2 = new string('A', 6) + w + new string('A', 6);
        int wPosInS2 = 6;

        // Self-check: W is the unique longest common substring (poly-C/poly-A share no run > 1).
        BruteLcsLength(s1, s2).Should().Be(blockLen, "the planted block is the LCS");

        var dna1 = new DnaSequence(s1);
        var dna2 = new DnaSequence(s2);

        var longest = GenomicAnalyzer.FindLongestCommonRegion(dna1, dna2);
        longest.Length.Should().Be(blockLen, "LCS length equals the block length");
        longest.Sequence.Should().Be(w);
        s1.Should().Contain(longest.Sequence);
        s2.Should().Contain(longest.Sequence);

        var regions = GenomicAnalyzer.FindCommonRegions(dna1, dna2, minLen).ToList();
        regions.Should().OnlyContain(r => r.Length >= minLen, "every region meets minLength");
        regions.Should().OnlyContain(r => s1.Contains(r.Sequence) && s2.Contains(r.Sequence), "regions occur in both sequences");
        regions.Should().OnlyContain(r => s2.Substring(r.PositionInSecond, r.Length) == r.Sequence, "PositionInSecond is correct");

        bool expectBlock = minLen <= blockLen;
        regions.Any(r => r.PositionInSecond == wPosInS2 && r.Sequence == w)
            .Should().Be(expectBlock, "the planted block is a region iff minLen ≤ its length");
    }

    /// <summary>
    /// Interaction witness — the minLength floor gates FindCommonRegions while FindLongestCommonRegion
    /// (which has no threshold) always returns the LCS.
    /// </summary>
    [Test]
    public void GenomicCommon_MinLengthGatesRegionsNotLcs()
    {
        var s1 = new DnaSequence("CCCCCC" + CommonBlock + "CCCCCC");
        var s2 = new DnaSequence("AAAAAA" + CommonBlock + "AAAAAA");

        GenomicAnalyzer.FindLongestCommonRegion(s1, s2).Sequence.Should().Be(CommonBlock, "LCS ignores minLength");

        GenomicAnalyzer.FindCommonRegions(s1, s2, 6).Should().Contain(r => r.Sequence == CommonBlock);
        GenomicAnalyzer.FindCommonRegions(s1, s2, CommonBlock.Length + 1)
            .Should().BeEmpty("no common substring is longer than the block");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: GENOMIC-ORF-001 — Six-frame ORF detection (Analysis)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 177.
    // Spec: tests/TestSpecs/GENOMIC-ORF-001.md (canonical FindOpenReadingFrames). ADVANCED §10.
    // Dimensions: minLen(3) × frame(3) × startCodon(2). Grid 3×3×2 = 18 (full, exhaustive).
    //
    // Model (Rosalind ORF; NCBI ORFfinder): scan all six frames; at every ATG find the first in-frame
    // stop (TAA/TAG/TGA) downstream and report that ATG→stop span (length divisible by 3). ONLY ATG is
    // a start codon — a GTG start is never reported (documented: StartCodon = "ATG").
    //
    // The combinatorial point: minimum length, the frame/strand the ORF sits in, and the start codon
    // interact — a planted 36-nt ORF is detected exactly when the start is ATG AND minLen ≤ 36, on
    // whichever strand/frame it was placed; every reported ORF starts with ATG, ends with a stop, and
    // has length ≥ minLen divisible by 3.
    // ═══════════════════════════════════════════════════════════════════════

    public enum OrfPlacement { Forward0, Forward1, Reverse }

    private static string RevComp(string s)
    {
        var chars = new char[s.Length];
        for (int i = 0; i < s.Length; i++)
            chars[s.Length - 1 - i] = s[i] switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N' };
        return new string(chars);
    }

    [Test, Combinatorial]
    public void GenomicOrf_DetectsPlantedOrf_AcrossMinLenFramePlacementStart(
        [Values(12, 30, 60)] int minLen,
        [Values(OrfPlacement.Forward0, OrfPlacement.Forward1, OrfPlacement.Reverse)] OrfPlacement placement,
        [Values("ATG", "GTG")] string startCodon)
    {
        string orf = startCodon + string.Concat(Enumerable.Repeat("CAA", 10)) + "TAA"; // 36 nt, no internal stop
        const int orfLen = 36;

        string text = placement switch
        {
            OrfPlacement.Forward0 => new string('T', 6) + orf + new string('T', 6),
            OrfPlacement.Forward1 => new string('T', 7) + orf + new string('T', 6), // shift reading frame
            OrfPlacement.Reverse => new string('T', 6) + RevComp(orf) + new string('T', 6),
            _ => orf,
        };
        var dna = new DnaSequence(text);

        var orfs = GenomicAnalyzer.FindOpenReadingFrames(dna, minLen).ToList();

        // Every reported ORF respects the ATG…stop contract.
        orfs.Should().OnlyContain(o => o.Sequence.StartsWith("ATG", StringComparison.Ordinal), "ORFs start at ATG");
        orfs.Should().OnlyContain(o => o.Length % 3 == 0 && o.Length >= minLen, "length divisible by 3 and ≥ minLen");
        orfs.Should().OnlyContain(o =>
            o.Sequence.EndsWith("TAA", StringComparison.Ordinal) ||
            o.Sequence.EndsWith("TAG", StringComparison.Ordinal) ||
            o.Sequence.EndsWith("TGA", StringComparison.Ordinal), "ORFs end at a stop codon");

        bool expectDetected = startCodon == "ATG" && orfLen >= minLen;
        bool onReverse = placement == OrfPlacement.Reverse;
        orfs.Any(o => o.Sequence == "ATG" + string.Concat(Enumerable.Repeat("CAA", 10)) + "TAA" && o.IsReverseComplement == onReverse)
            .Should().Be(expectDetected, "the planted ORF is found iff start is ATG and it clears minLen");
    }

    /// <summary>
    /// Interaction witness — only ATG starts an ORF (a GTG start yields nothing), and the minLength
    /// floor excludes shorter ORFs.
    /// </summary>
    [Test]
    public void GenomicOrf_OnlyAtgStarts_AndMinLengthFilters()
    {
        string body = string.Concat(Enumerable.Repeat("CAA", 10)) + "TAA";
        var atg = new DnaSequence("TTTTTT" + "ATG" + body + "TTTTTT");
        var gtg = new DnaSequence("TTTTTT" + "GTG" + body + "TTTTTT");

        GenomicAnalyzer.FindOpenReadingFrames(atg, 12).Should().Contain(o => o.Sequence == "ATG" + body);
        GenomicAnalyzer.FindOpenReadingFrames(gtg, 12).Should().NotContain(o => o.Sequence == "GTG" + body, "GTG is not a start codon");

        GenomicAnalyzer.FindOpenReadingFrames(atg, 100).Should().NotContain(o => o.Sequence == "ATG" + body,
            "a 36-nt ORF is excluded at minLength 100");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: GENOMIC-REPEAT-001 — Repeated substrings (Analysis)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 178.
    // Spec: tests/TestSpecs/GENOMIC-REPEAT-001.md (canonical FindRepeats / FindLongestRepeat). ADVANCED §10.
    // Dimensions: minLen(3) × seqLen(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Gusfield 1997; CMU 15-451): a repeat is a substring occurring ≥ 2 times; FindRepeats
    // returns every distinct such substring of length ≥ minLength, and FindLongestRepeat returns the
    // longest repeated substring (LRS).
    //
    // The combinatorial point: minLength and sequence length interact — FindRepeats equals the
    // independent brute-force set of substrings occurring ≥ 2 times with length ≥ minLength (with
    // matching overlapping occurrence positions), and the LRS length equals the maximum repeat length.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void GenomicRepeat_AllRepeatsAndLrs_AcrossMinLenAndLength(
        [Values(3, 5, 8)] int minLen,
        [Values(32, 48, 64)] int seqLen)
    {
        // A periodic text (16-nt block tiled) guarantees rich, non-vacuous repeats.
        string block = "ACGTTGCAACGGTACC";
        string text = string.Concat(Enumerable.Repeat(block, seqLen / block.Length + 1))[..seqLen];
        var dna = new DnaSequence(text);

        // Independent brute force: every distinct substring of length ≥ minLen occurring ≥ 2 times.
        var brute = new Dictionary<string, int>();
        for (int len = minLen; len <= seqLen; len++)
            for (int i = 0; i + len <= seqLen; i++)
            {
                string sub = text.Substring(i, len);
                if (!brute.ContainsKey(sub)) brute[sub] = CountOverlapping(text, sub);
            }
        var bruteRepeats = brute.Where(kv => kv.Value >= 2).Select(kv => kv.Key).ToHashSet();

        var repeats = GenomicAnalyzer.FindRepeats(dna, minLen).ToList();

        repeats.Select(r => r.Sequence).Should().BeEquivalentTo(bruteRepeats,
            "FindRepeats = every distinct substring (len ≥ minLen) occurring ≥ 2 times");
        foreach (var r in repeats)
        {
            r.Length.Should().BeGreaterThanOrEqualTo(minLen);
            r.Count.Should().BeGreaterThanOrEqualTo(2, "a repeat occurs at least twice");
            r.Count.Should().Be(CountOverlapping(text, r.Sequence), "Count = overlapping occurrence count");
        }

        var lrs = GenomicAnalyzer.FindLongestRepeat(dna);
        lrs.Count.Should().BeGreaterThanOrEqualTo(2);
        lrs.Length.Should().Be(bruteRepeats.Max(s => s.Length), "the LRS is the longest repeated substring");
    }

    /// <summary>
    /// Interaction witness — the minLength floor is monotone (raising it can only drop repeats), and a
    /// sequence with no repeat (all-distinct characters region) yields none above that length.
    /// </summary>
    [Test]
    public void GenomicRepeat_MinLengthMonotone_AndNoRepeatCase()
    {
        var dna = new DnaSequence(string.Concat(Enumerable.Repeat("ACGTTGCAACGGTACC", 4)));
        var at3 = GenomicAnalyzer.FindRepeats(dna, 3).Select(r => r.Sequence).ToHashSet();
        var at6 = GenomicAnalyzer.FindRepeats(dna, 6).Select(r => r.Sequence).ToHashSet();
        at6.Should().BeSubsetOf(at3, "a higher minLength retains no more repeats");

        // No substring of length ≥ 4 repeats in a De Bruijn-like all-distinct 4-mer run.
        GenomicAnalyzer.FindLongestRepeat(new DnaSequence("ACGT")).IsEmpty.Should().BeTrue("no substring repeats in ACGT");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: GENOMIC-TANDEM-001 — Tandem repeats (Analysis)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 180.
    // Spec: tests/TestSpecs/GENOMIC-TANDEM-001.md (canonical FindTandemRepeats; same method as REP-TANDEM-001).
    // ADVANCED §10.
    // Dimensions: unitLen(3) × minReps(3) × seqLen(3). Grid 3×3×3 = 27 (full, exhaustive).
    //
    // Model (Benson 1999 Tandem Repeats Finder): a tandem repeat is ≥ minRepetitions consecutive copies
    // of a unit of length ≥ minUnitLength. FindTandemRepeats reports (Unit, Position, Repetitions).
    //
    // The combinatorial point: minUnitLength, minRepetitions and sequence length interact — every
    // reported tandem is a genuine run (Unit repeated Repetitions times exactly at Position in the text,
    // honouring both floors), and a planted "ACGG"×5 run is detected whenever both floors admit it.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void GenomicTandem_ValidRunsAndPlantedDetection_AcrossUnitRepsLength(
        [Values(2, 3, 4)] int minUnitLength,
        [Values(2, 3, 4)] int minRepetitions,
        [Values(40, 60, 80)] int seqLen)
    {
        int fillLen = (seqLen - 24) / 2; // planted run is 2+20+2 = 24 nt
        string text = DiverseDna(fillLen, 0x111u) + "TT" + string.Concat(Enumerable.Repeat("ACGG", 5)) + "TT"
            + DiverseDna(fillLen, 0x222u);
        int plantedPos = fillLen + 2;
        var dna = new DnaSequence(text);

        var tandems = GenomicAnalyzer.FindTandemRepeats(dna, minUnitLength, minRepetitions).ToList();

        // Every reported tandem is a real run honouring both floors.
        foreach (var t in tandems)
        {
            t.Unit.Length.Should().BeGreaterThanOrEqualTo(minUnitLength, "unit meets minUnitLength");
            t.Repetitions.Should().BeGreaterThanOrEqualTo(minRepetitions, "run meets minRepetitions");
            t.FullSequence.Should().Be(string.Concat(Enumerable.Repeat(t.Unit, t.Repetitions)), "FullSequence is the unit tiled");
            text.Substring(t.Position, t.TotalLength).Should().Be(t.FullSequence, "the run actually occurs at Position");
        }

        // The planted ACGG×5 run (unit 4, 5 copies) is detected when both floors admit it.
        bool expectPlanted = minUnitLength <= 4 && minRepetitions <= 5;
        tandems.Any(t => t.Position == plantedPos && t.Unit == "ACGG" && t.Repetitions == 5)
            .Should().Be(expectPlanted, "the planted ACGG×5 tandem is found iff both floors admit it");
    }

    /// <summary>
    /// Interaction witness — each floor independently gates the planted tandem: too-long minUnitLength
    /// or too-many minRepetitions removes it.
    /// </summary>
    [Test]
    public void GenomicTandem_EachFloorGatesDetection()
    {
        var dna = new DnaSequence("GTGT" + "TT" + string.Concat(Enumerable.Repeat("ACGG", 5)) + "TTGTGT");

        bool Found(int minUnit, int minReps) =>
            GenomicAnalyzer.FindTandemRepeats(dna, minUnit, minReps).Any(t => t.Unit == "ACGG" && t.Repetitions == 5);

        Found(2, 2).Should().BeTrue("permissive floors find the ACGG×5 run");
        Found(5, 2).Should().BeFalse("minUnitLength 5 > unit length 4");
        Found(4, 6).Should().BeFalse("minRepetitions 6 > 5 copies");
    }
}
