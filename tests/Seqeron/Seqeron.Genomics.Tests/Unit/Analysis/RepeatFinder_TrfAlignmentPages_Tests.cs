// REP-APPROX-001 — TRF 4.10.0 alignment pages (<prefix>.<parameters>.N.txt.html) (B04 audit WP17, F64).
// Evidence: docs/Evidence/REP-APPROX-001-Evidence.md (§WP17); TestSpec: tests/TestSpecs/REP-APPROX-001.md (O13..O19).
// Sources: TRF 4.10.0 README ("Alignment explanation", "-f"); TRF 4.10.0 source read for the layout only
//          (trfrun.h alignment-file heading, tr30dat.c print_alignment_headings / alt3_print_alignment /
//          get_statistics / printECtoAlignments / print_flanking_sequence, trfclean.h CleanAlignments / BreakAlignments).
//
// Every expected page below was written by the compiled TRF 4.10.0 binary (github.com/Benson-Genomics-Lab/TRF, commit
// 355c1f9) run without -h; the multi-page case is locked by the SHA-256 of TRF's two pages. TRF's apparent-size table is
// not embedded: every case here gives the same repeats with the library's exact table.

using System.Security.Cryptography;
using System.Text;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RepeatFinder_TrfAlignmentPages_Tests
{
    // ">U1" (183 bp, two N).
    private const string U1 =
        "CTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTCATTGGTCATTGGTCANTGGTCATTGGTCATTGGTCATTNGTCATTGGTCATTGGTCATTGGTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGA";

    // set700 sequence 274 (">s273", 642 bp): a period-6 repeat with a deletion (copies separated by blanks).
    private const string S273 =
        "TTATTACCAGCCGGAGAGGTTTCGCTTATACGTTAGCAGGTAGCAAGAGGTTTAGTCAGCCTCGTCTTTTAAAGACACACGTATCATGGGATGTACACTT"
        + "TTTACAACACACTCATAGCAGGCTACGGCTCAAATTGATGGTCCCAAGTCAGCCGCGATGTCCGAGGACGGTCTATCAACGGGTGGTTTATCCCAAGAGA"
        + "AAAGTAGGTTTCATACATCCGTCGGCCGCTGGGCCGCAGTAAGAGCACTGGTGAGGGTGAGGTTGAGGTTGAGGTTGAGGTTGAGGTTGAGGTTGAGGTT"
        + "GAGGTTGAGGTTGAGGTTGAGGTTGCGGTTGAGGTTGAGGTTGAGGAGCCCGTCTTGCGCGTCGGACAGCGATGAGCGGATGGCATGACCGCTGTACATT"
        + "TACTCAGCCAAGTACTTTAACTGCTCTCGAATTAGGCTAGCTGTGCGCGGCTTCCCAACTCATTAATCCGAAAGTGAGCTTGTCTCGAGATCCCAGCACT"
        + "GAACAGTACTCCAGATACCTTTGCCATGACTCCTCTTCCATAATACTCGATAACTGTAAGGCGGTTCACAGGCTGGGGTCTGTAAATTGCTTAGCTTTAG"
        + "CGTACACCTGGGTTATTTCCGCCCTGAATCTATGAAATTCGT";

    // ">s12" (332 bp, lower case): one period-34 repeat with insertions and deletions; TRF reports two alignments and
    // drops the first as redundant.
    private const string S12 =
        "ggtgtcgtgagatcactcgcggttctcccggggggctcttgtaggcgtgagctcgaactgttcgccaaacagaccacgtgacacagtagctaggagcgtt"
        + "gatagtctatctctgtgtggaagggaagacttgcttagtccatcctgtgtggagaggaagacgttgcagtagtctatcctgtgtggaagggatgacgttg"
        + "catccgctcgaaagggccacatatttgtcggtaaatcccggtagtctctatccagattaatagccggcagcagagtttactcacgtcagcagtggtgttt"
        + "tggacaccagtactactattcgatcggacggt";

    private static TandemRepeatsFinderParameters Set(int match, int mismatch, int delta, int pm, int pi, int minScore, int maxPeriod) =>
        new()
        {
            MatchWeight = match, MismatchPenalty = mismatch, IndelPenalty = delta,
            MatchProbability = pm, IndelProbability = pi, MinScore = minScore, MaxPeriod = maxPeriod,
        };

    // O13: `trf U1.fa 2 7 7 80 10 50 500` → U1.fa.2.7.7.80.10.50.500.1.txt.html (period 7 → one row per copy, '*' under
    // the N mismatches, 10-base context, distance table; a later alignment was dropped → blank line before "Done.").
    [Test]
    public void U1_MatchesTrfByteForByte()
    {
        const string trf = "<HTML><HEAD><TITLE>U1.fa.2.7.7.80.10.50.500.txt.html</TITLE></HEAD><BODY bgcolor=\"#FBF8BC\"><PRE>\n"
        + "Tandem Repeats Finder Program written by:\n"
        + "\n"
        + "                 Gary Benson\n"
        + "      Program in Bioinformatics\n"
        + "          Boston University\n"
        + "\n"
        + "Version 4.10.0\n"
        + "\n"
        + "Sequence: U1\n"
        + "\n"
        + "Parameters: 2 7 7 80 10 50 500\n"
        + "\n"
        + "Pmatch=0.80,Pindel=0.10\n"
        + "tuple sizes 0,4,5,7\n"
        + "tuple distances 0, 29, 159, 200\n"
        + "\n"
        + "Length: 183\n"
        + "ACGTcount: A:0.20, C:0.25, G:0.23, T:0.31\n"
        + "\n"
        + "Warning! 2 characters in sequence are not A, C, G, or T\n"
        + "\n"
        + "\n"
        + "Found at i:71 original size:7 final size:7\n"
        + "\n"
        + "<A NAME=\"61--124,7,9.1,7,1\"></A><A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#alignment\" target =\"explanation\">Alignment explanation</A><BR><BR>\n"
        + "    Indices: 61--124  Score: 110\n"
        + "    Period size: 7  Copynumber: 9.1  Consensus size: 7\n"
        + "\n"
        + "         51 TCATTTCCGC\n"
        + "\n"
        + "                   \n"
        + "         61 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                   \n"
        + "         68 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "               *   \n"
        + "         75 TCANTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                   \n"
        + "         82 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                   \n"
        + "         89 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                 * \n"
        + "         96 TCATTNG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                   \n"
        + "        103 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                   \n"
        + "        110 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                   \n"
        + "        117 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "             \n"
        + "        124 T\n"
        + "          1 T\n"
        + "\n"
        + "        125 AGACATAATC\n"
        + "\n"
        + "\n"
        + "Statistics\n"
        + "Matches: 53,  Mismatches: 4, Indels: 0\n"
        + "        0.93            0.07        0.00\n"
        + "\n"
        + "Matches are distributed among these distances:\n"
        + "   7   53  1.00\n"
        + "\n"
        + "ACGTcount: A:0.14, C:0.14, G:0.27, T:0.42\n"
        + "\n"
        + "\n"
        + "Consensus pattern (7 bp):   \n"
        + "TCATTGG\n"
        + "\n"
        + "Done.\n"
        + "</PRE></BODY></HTML>\n";
        var p = TandemRepeatsFinderParameters.Recommended;
        var found = RepeatFinder.FindApproximateTandemRepeats(U1, p, out int outputCount);
        Assert.Multiple(() =>
        {
            Assert.That(found, Has.Count.EqualTo(1));
            Assert.That(found[0].DetectionPosition + 1, Is.EqualTo(71), "Found at i:71");
            Assert.That(found[0].DetectionDistance, Is.EqualTo(7), "original size:7");
            Assert.That(outputCount, Is.EqualTo(3), "periods 7, 14, 21 reported; 14 and 21 dropped as redundant");
            Assert.That(found[0].OutputCount, Is.EqualTo(3));
        });
        var pages = RepeatFinder.FormatTrfAlignmentPages(U1, found, "U1", p, "U1.fa");
        Assert.That(pages, Has.Count.EqualTo(1));
        Assert.That(pages[0].FileName, Is.EqualTo("U1.fa.2.7.7.80.10.50.500.1.txt.html"));
        Assert.That(pages[0].Html, Is.EqualTo(trf));
        Assert.That(RepeatFinder.FormatTrfAlignmentPages(U1, found, "U1", p, "U1.fa", outputCount)[0].Html, Is.EqualTo(trf));
    }

    // O14: `trf U1.fa 2 7 7 80 10 50 500 -r -f` → all three alignments (periods 7, 14, 21), flanking sequences (TRF -f = 500).
    [Test]
    public void U1_NoRedundancyElimination_WithFlanks_MatchesTrf()
    {
        const string trf = "<HTML><HEAD><TITLE>U1.fa.2.7.7.80.10.50.500.txt.html</TITLE></HEAD><BODY bgcolor=\"#FBF8BC\"><PRE>\n"
        + "Tandem Repeats Finder Program written by:\n"
        + "\n"
        + "                 Gary Benson\n"
        + "      Program in Bioinformatics\n"
        + "          Boston University\n"
        + "\n"
        + "Version 4.10.0\n"
        + "\n"
        + "Sequence: U1\n"
        + "\n"
        + "Parameters: 2 7 7 80 10 50 500\n"
        + "\n"
        + "Pmatch=0.80,Pindel=0.10\n"
        + "tuple sizes 0,4,5,7\n"
        + "tuple distances 0, 29, 159, 200\n"
        + "\n"
        + "Length: 183\n"
        + "ACGTcount: A:0.20, C:0.25, G:0.23, T:0.31\n"
        + "\n"
        + "Warning! 2 characters in sequence are not A, C, G, or T\n"
        + "\n"
        + "\n"
        + "Found at i:71 original size:7 final size:7\n"
        + "\n"
        + "<A NAME=\"61--124,7,9.1,7,1\"></A><A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#alignment\" target =\"explanation\">Alignment explanation</A><BR><BR>\n"
        + "    Indices: 61--124  Score: 110\n"
        + "    Period size: 7  Copynumber: 9.1  Consensus size: 7\n"
        + "\n"
        + "         51 TCATTTCCGC\n"
        + "\n"
        + "                   \n"
        + "         61 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                   \n"
        + "         68 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "               *   \n"
        + "         75 TCANTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                   \n"
        + "         82 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                   \n"
        + "         89 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                 * \n"
        + "         96 TCATTNG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                   \n"
        + "        103 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                   \n"
        + "        110 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "                   \n"
        + "        117 TCATTGG\n"
        + "          1 TCATTGG\n"
        + "\n"
        + "             \n"
        + "        124 T\n"
        + "          1 T\n"
        + "\n"
        + "        125 AGACATAATC\n"
        + "\n"
        + "\n"
        + "Statistics\n"
        + "Matches: 53,  Mismatches: 4, Indels: 0\n"
        + "        0.93            0.07        0.00\n"
        + "\n"
        + "Matches are distributed among these distances:\n"
        + "   7   53  1.00\n"
        + "\n"
        + "ACGTcount: A:0.14, C:0.14, G:0.27, T:0.42\n"
        + "\n"
        + "\n"
        + "Consensus pattern (7 bp):   \n"
        + "TCATTGG\n"
        + "\n"
        + "Left flanking sequence: Indices 1 -- 60\n"
        + "CTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGC\n"
        + "\n"
        + "\n"
        + "Right flanking sequence: Indices 125 -- 183\n"
        + "AGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGA\n"
        + "\n"
        + "\n"
        + "\n"
        + "Found at i:84 original size:14 final size:14\n"
        + "\n"
        + "<A NAME=\"61--124,14,4.6,14,2\"></A><A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#alignment\" target =\"explanation\">Alignment explanation</A><BR><BR>\n"
        + "    Indices: 61--124  Score: 110\n"
        + "    Period size: 14  Copynumber: 4.6  Consensus size: 14\n"
        + "\n"
        + "         51 TCATTTCCGC\n"
        + "\n"
        + "                          \n"
        + "         61 TCATTGGTCATTGG\n"
        + "          1 TCATTGGTCATTGG\n"
        + "\n"
        + "               *          \n"
        + "         75 TCANTGGTCATTGG\n"
        + "          1 TCATTGGTCATTGG\n"
        + "\n"
        + "                        * \n"
        + "         89 TCATTGGTCATTNG\n"
        + "          1 TCATTGGTCATTGG\n"
        + "\n"
        + "                          \n"
        + "        103 TCATTGGTCATTGG\n"
        + "          1 TCATTGGTCATTGG\n"
        + "\n"
        + "                    \n"
        + "        117 TCATTGGT\n"
        + "          1 TCATTGGT\n"
        + "\n"
        + "        125 AGACATAATC\n"
        + "\n"
        + "\n"
        + "Statistics\n"
        + "Matches: 46,  Mismatches: 4, Indels: 0\n"
        + "        0.92            0.08        0.00\n"
        + "\n"
        + "Matches are distributed among these distances:\n"
        + "  14   46  1.00\n"
        + "\n"
        + "ACGTcount: A:0.14, C:0.14, G:0.27, T:0.42\n"
        + "\n"
        + "\n"
        + "Consensus pattern (14 bp):   \n"
        + "TCATTGGTCATTGG\n"
        + "\n"
        + "Left flanking sequence: Indices 1 -- 60\n"
        + "CTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGC\n"
        + "\n"
        + "\n"
        + "Right flanking sequence: Indices 125 -- 183\n"
        + "AGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGA\n"
        + "\n"
        + "\n"
        + "\n"
        + "Found at i:87 original size:21 final size:21\n"
        + "\n"
        + "<A NAME=\"61--124,21,3.0,21,3\"></A><A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#alignment\" target =\"explanation\">Alignment explanation</A><BR><BR>\n"
        + "    Indices: 61--124  Score: 110\n"
        + "    Period size: 21  Copynumber: 3.0  Consensus size: 21\n"
        + "\n"
        + "         51 TCATTTCCGC\n"
        + "\n"
        + "                             *   \n"
        + "         61 TCATTGGTCATTGGTCANTGG\n"
        + "          1 TCATTGGTCATTGGTCATTGG\n"
        + "\n"
        + "                               * \n"
        + "         82 TCATTGGTCATTGGTCATTNG\n"
        + "          1 TCATTGGTCATTGGTCATTGG\n"
        + "\n"
        + "                                 \n"
        + "        103 TCATTGGTCATTGGTCATTGG\n"
        + "          1 TCATTGGTCATTGGTCATTGG\n"
        + "\n"
        + "             \n"
        + "        124 T\n"
        + "          1 T\n"
        + "\n"
        + "        125 AGACATAATC\n"
        + "\n"
        + "\n"
        + "Statistics\n"
        + "Matches: 40,  Mismatches: 3, Indels: 0\n"
        + "        0.93            0.07        0.00\n"
        + "\n"
        + "Matches are distributed among these distances:\n"
        + "  21   40  1.00\n"
        + "\n"
        + "ACGTcount: A:0.14, C:0.14, G:0.27, T:0.42\n"
        + "\n"
        + "\n"
        + "Consensus pattern (21 bp):   \n"
        + "TCATTGGTCATTGGTCATTGG\n"
        + "\n"
        + "Left flanking sequence: Indices 1 -- 60\n"
        + "CTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGC\n"
        + "\n"
        + "\n"
        + "Right flanking sequence: Indices 125 -- 183\n"
        + "AGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGA\n"
        + "\n"
        + "\n"
        + "Done.\n"
        + "</PRE></BODY></HTML>\n";
        var p = TandemRepeatsFinderParameters.Recommended with { EliminateRedundancy = false, FlankLength = 500 };
        var found = RepeatFinder.FindApproximateTandemRepeats(U1, p).ToList();
        Assert.That(RepeatFinder.FormatTrfAlignmentPages(U1, found, "U1", p, "U1.fa")[0].Html, Is.EqualTo(trf));
    }

    // O15: `trf set700.fa 2 3 3 80 20 50 500` (sequence 274) → set700.fa.s274…1.txt.html: period ≤ 6, so copies share a
    // row separated by blanks while two more fit; the sequence row shows the deletion ('-') and its index skips nothing.
    [Test]
    public void SmallPeriod_CopiesShareRows_MatchesTrf()
    {
        const string trf = "<HTML><HEAD><TITLE>set700.fa.s274.2.3.3.80.20.50.500.txt.html</TITLE></HEAD><BODY bgcolor=\"#FBF8BC\"><PRE>\n"
        + "Tandem Repeats Finder Program written by:\n"
        + "\n"
        + "                 Gary Benson\n"
        + "      Program in Bioinformatics\n"
        + "          Boston University\n"
        + "\n"
        + "Version 4.10.0\n"
        + "\n"
        + "Sequence: s273\n"
        + "\n"
        + "Parameters: 2 3 3 80 20 50 500\n"
        + "\n"
        + "Pmatch=0.80,Pindel=0.20\n"
        + "tuple sizes 0,4,5,7\n"
        + "tuple distances 0, 29, 159, 385\n"
        + "\n"
        + "Length: 642\n"
        + "ACGTcount: A:0.23, C:0.22, G:0.28, T:0.27\n"
        + "\n"
        + "\n"
        + "Found at i:261 original size:6 final size:6\n"
        + "\n"
        + "<A NAME=\"250--346,6,16.3,6,1\"></A><A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#alignment\" target =\"explanation\">Alignment explanation</A><BR><BR>\n"
        + "    Indices: 250--346  Score: 181\n"
        + "    Period size: 6  Copynumber: 16.3  Consensus size: 6\n"
        + "\n"
        + "        240 TAAGAGCACT\n"
        + "\n"
        + "                     *                                             \n"
        + "        250 GG-TGA GGGTGA GGTTGA GGTTGA GGTTGA GGTTGA GGTTGA GGTTGA\n"
        + "          1 GGTTGA GGTTGA GGTTGA GGTTGA GGTTGA GGTTGA GGTTGA GGTTGA\n"
        + "\n"
        + "                                             *                     \n"
        + "        297 GGTTGA GGTTGA GGTTGA GGTTGA GGTTGC GGTTGA GGTTGA GGTTGA\n"
        + "          1 GGTTGA GGTTGA GGTTGA GGTTGA GGTTGA GGTTGA GGTTGA GGTTGA\n"
        + "\n"
        + "              \n"
        + "        345 GG\n"
        + "          1 GG\n"
        + "\n"
        + "        347 AGCCCGTCTT\n"
        + "\n"
        + "\n"
        + "Statistics\n"
        + "Matches: 88,  Mismatches: 3, Indels: 1\n"
        + "        0.96            0.03        0.01\n"
        + "\n"
        + "Matches are distributed among these distances:\n"
        + "   5    2  0.02\n"
        + "   6   86  0.98\n"
        + "\n"
        + "ACGTcount: A:0.15, C:0.01, G:0.53, T:0.31\n"
        + "\n"
        + "\n"
        + "Consensus pattern (6 bp):   \n"
        + "GGTTGA\n"
        + "\n"
        + "Done.\n"
        + "</PRE></BODY></HTML>\n";
        var p = Set(2, 3, 3, 80, 20, 50, 500);
        var found = RepeatFinder.FindApproximateTandemRepeats(S273, p).ToList();
        var pages = RepeatFinder.FormatTrfAlignmentPages(S273, found, "s273", p, "set700.fa.s274");
        Assert.That(pages.Single().FileName, Is.EqualTo("set700.fa.s274.2.3.3.80.20.50.500.1.txt.html"));
        Assert.That(pages.Single().Html, Is.EqualTo(trf));
    }

    // O16: `trf s12.fa 2 7 7 80 10 50 500` (lower-case input) → period 34: insertions in both rows, the matching-distance
    // table 31..35, OUTPUTcount label 2 (alignment 1 removed as redundant), no blank line before "Done.".
    [Test]
    public void LargePeriodWithIndels_MatchesTrf()
    {
        const string trf = "<HTML><HEAD><TITLE>q12.fa.2.7.7.80.10.50.500.txt.html</TITLE></HEAD><BODY bgcolor=\"#FBF8BC\"><PRE>\n"
        + "Tandem Repeats Finder Program written by:\n"
        + "\n"
        + "                 Gary Benson\n"
        + "      Program in Bioinformatics\n"
        + "          Boston University\n"
        + "\n"
        + "Version 4.10.0\n"
        + "\n"
        + "Sequence: s12\n"
        + "\n"
        + "Parameters: 2 7 7 80 10 50 500\n"
        + "\n"
        + "Pmatch=0.80,Pindel=0.10\n"
        + "tuple sizes 0,4,5,7\n"
        + "tuple distances 0, 29, 159, 200\n"
        + "\n"
        + "Length: 332\n"
        + "ACGTcount: A:0.22, C:0.22, G:0.30, T:0.26\n"
        + "\n"
        + "\n"
        + "Found at i:180 original size:34 final size:34\n"
        + "\n"
        + "<A NAME=\"103--202,34,3.0,34,2\"></A><A HREF=\"http://tandem.bu.edu/trf/trf.definitions.html#alignment\" target =\"explanation\">Alignment explanation</A><BR><BR>\n"
        + "    Indices: 103--202  Score: 134\n"
        + "    Period size: 34  Copynumber: 3.0  Consensus size: 34\n"
        + "\n"
        + "         93 GGAGCGTTGA\n"
        + "\n"
        + "                                              *\n"
        + "        103 TAGTCTATCTCTGTGTGGAAGGGAAGAC-TTGC-T\n"
        + "          1 TAGTCTATC-CTGTGTGGAAGGGAAGACGTTGCAG\n"
        + "\n"
        + "                 *                             \n"
        + "        136 TAGTCCATCCTGTGTGG-AGAGGAAGACGTTGCAG\n"
        + "          1 TAGTCTATCCTGTGTGGAAG-GGAAGACGTTGCAG\n"
        + "\n"
        + "                                   *         \n"
        + "        170 TAGTCTATCCTGTGTGGAAGGGATGACGTTGCA\n"
        + "          1 TAGTCTATCCTGTGTGGAAGGGAAGACGTTGCA\n"
        + "\n"
        + "        203 TCCGCTCGAA\n"
        + "\n"
        + "\n"
        + "Statistics\n"
        + "Matches: 59,  Mismatches: 4, Indels: 7\n"
        + "        0.84            0.06        0.10\n"
        + "\n"
        + "Matches are distributed among these distances:\n"
        + "  31    2  0.03\n"
        + "  32   15  0.25\n"
        + "  33   12  0.20\n"
        + "  34   28  0.47\n"
        + "  35    2  0.03\n"
        + "\n"
        + "ACGTcount: A:0.22, C:0.16, G:0.33, T:0.29\n"
        + "\n"
        + "\n"
        + "Consensus pattern (34 bp):   \n"
        + "TAGTCTATCCTGTGTGGAAGGGAAGACGTTGCAG\n"
        + "Done.\n"
        + "</PRE></BODY></HTML>\n";
        var p = TandemRepeatsFinderParameters.Recommended;
        var found = RepeatFinder.FindApproximateTandemRepeats(S12, p, out int outputCount);
        Assert.That(found.Single().OutputIndex, Is.EqualTo(2));
        Assert.That(outputCount, Is.EqualTo(2));
        Assert.That(RepeatFinder.FormatTrfAlignmentPages(S12, found, "s12", p, "q12.fa").Single().Html, Is.EqualTo(trf));
    }

    // O17: `trf s12.fa 2 7 7 80 10 50 3` → no repeat kept (period 34 > MaxPeriod) but two were reported: the heading, one
    // extra blank line, "Done.". Without the output count the blank line is missing.
    [Test]
    public void NoRepeatKept_UsesOutputCount()
    {
        const string trf = "<HTML><HEAD><TITLE>q12.fa.2.7.7.80.10.50.3.txt.html</TITLE></HEAD><BODY bgcolor=\"#FBF8BC\"><PRE>\n"
        + "Tandem Repeats Finder Program written by:\n"
        + "\n"
        + "                 Gary Benson\n"
        + "      Program in Bioinformatics\n"
        + "          Boston University\n"
        + "\n"
        + "Version 4.10.0\n"
        + "\n"
        + "Sequence: s12\n"
        + "\n"
        + "Parameters: 2 7 7 80 10 50 3\n"
        + "\n"
        + "Pmatch=0.80,Pindel=0.10\n"
        + "tuple sizes 0,4,5,7\n"
        + "tuple distances 0, 29, 159, 200\n"
        + "\n"
        + "Length: 332\n"
        + "ACGTcount: A:0.22, C:0.22, G:0.30, T:0.26\n"
        + "\n"
        + "\n"
        + "Done.\n"
        + "</PRE></BODY></HTML>\n";
        var p = TandemRepeatsFinderParameters.Recommended with { MaxPeriod = 3 };
        var found = RepeatFinder.FindApproximateTandemRepeats(S12, p, out int outputCount);
        Assert.That(found, Is.Empty);
        Assert.That(outputCount, Is.EqualTo(2));
        Assert.That(RepeatFinder.FormatTrfAlignmentPages(S12, found, "s12", p, "q12.fa", outputCount).Single().Html, Is.EqualTo(trf));
        string withoutCount = RepeatFinder.FormatTrfAlignmentPages(S12, found, "s12", p, "q12.fa").Single().Html;
        Assert.That(withoutCount, Is.EqualTo(trf.Replace("\n\n\nDone.", "\n\nDone.")));
    }

    // O18: 150 repeats → two pages (TRF BreakAlignments: heading repeated, "File k of 2", 120 + 30 alignments). Sequence:
    // 150 blocks of 40 pseudo-random bases + a random trinucleotide × 12 (LCG x ← 1103515245·x + 12345 mod 2^31, base =
    // (x >> 16) & 3, seed 12345). `trf big.fa 2 7 7 80 10 50 500` page hashes below.
    [Test]
    public void TwoPages_MatchTrfHashes()
    {
        var sb = new StringBuilder();
        long x = 12345;
        char Next()
        {
            x = (x * 1103515245 + 12345) % 2147483648;
            return "ACGT"[(int)((x >> 16) & 3)];
        }

        for (int b = 0; b < 150; b++)
        {
            for (int k = 0; k < 40; k++)
                sb.Append(Next());
            string unit = new([Next(), Next(), Next()]);
            for (int k = 0; k < 12; k++)
                sb.Append(unit);
        }

        string seq = sb.ToString();
        var p = TandemRepeatsFinderParameters.Recommended;
        var found = RepeatFinder.FindApproximateTandemRepeats(seq, p).ToList();
        Assert.That(found, Has.Count.EqualTo(150));
        var pages = RepeatFinder.FormatTrfAlignmentPages(seq, found, "big synthetic", p, "big.fa");
        static string Sha(string s) => Convert.ToHexString(SHA256.HashData(Encoding.ASCII.GetBytes(s))).ToLowerInvariant();
        Assert.Multiple(() =>
        {
            Assert.That(pages.Select(q => q.FileName), Is.EqualTo(new[] { "big.fa.2.7.7.80.10.50.500.1.txt.html", "big.fa.2.7.7.80.10.50.500.2.txt.html" }));
            Assert.That(pages[0].Html, Does.Contain("\nFile 1 of 2\n\nFound at i:"));
            Assert.That(pages[1].Html, Does.Contain("\nFile 2 of 2\n\nFound at i:9169 "));
            Assert.That(Sha(pages[0].Html), Is.EqualTo("1cc6f5cc235c5f41e000f1397f3aa890d596eeb4cfb28f5bdf58ee3371dd0b2e"));
            Assert.That(Sha(pages[1].Html), Is.EqualTo("e5a66df271ca931dbf826234d779ffd6e83697156f9d056e273776988b556685"));
        });

        // Every table link resolves to an anchor of the alignment page with the same number.
        var tables = RepeatFinder.FormatTrfHtmlTables(seq, found, "big synthetic", p, "big.fa");
        for (int k = 0; k < tables.Count; k++)
        {
            foreach (System.Text.RegularExpressions.Match m in System.Text.RegularExpressions.Regex.Matches(tables[k].Html, "HREF=\"(big[^\"#]+)#([^\"]+)\""))
            {
                Assert.That(m.Groups[1].Value, Is.EqualTo(pages[k].FileName));
                Assert.That(pages[k].Html, Does.Contain("<A NAME=\"" + m.Groups[2].Value + "\">"));
            }
        }
    }

    // O19: TRF's page split (trfclean.h BreakAlignments) reads 199-character line chunks; a chunk starting with 'F'
    // starts a new alignment and a chunk starting with 'D' ends the page. The reader keeps every letter (trfrun.h
    // LoadSequenceFromFileBenson: A..Z / a..z, upper-cased) and -f flank lines start at column 0, so an IUPAC 'D' or a
    // letter 'F' at the start of a flank line cuts the pages. Sequence: 150 blocks of 40 background letters from a
    // 32-letter alphabet (7 × ACGT + A, C, D, F, so D and F are 1/32 each) + a random trinucleotide × 12 (same LCG as
    // O18, background letter = alphabet[(x >> 16) mod 32]). `trf df.fa 2 7 7 80 10 50 500 -f` writes two pages. Page 1
    // holds only 4 "Found at": two 'F' flank lines each used up one of the 120 slots, then a 'D' flank line ended the
    // page. Page 2 holds just the two lines from that 'D' chunk to the next 'D' chunk. The other 146 alignments are lost,
    // and their table links (e.g. #345--380,3,12.0,3,13) dangle in TRF too. TRF's own apparent-size table and the exact
    // table give the same pages.
    [Test]
    public void FlankLinesStartingWithDOrF_CutPagesLikeTrf()
    {
        const string alphabet = "ACGTACGTACGTACGTACGTACGTACGTACDF";
        var sb = new StringBuilder();
        long x = 12345;
        char Next(string letters)
        {
            x = (x * 1103515245 + 12345) % 2147483648;
            return letters[(int)((x >> 16) % letters.Length)];
        }

        for (int b = 0; b < 150; b++)
        {
            for (int k = 0; k < 40; k++)
                sb.Append(Next(alphabet));
            string unit = new([Next("ACGT"), Next("ACGT"), Next("ACGT")]);
            for (int k = 0; k < 12; k++)
                sb.Append(unit);
        }

        string seq = sb.ToString();
        var p = TandemRepeatsFinderParameters.Recommended with { FlankLength = 500 };
        var found = RepeatFinder.FindApproximateTandemRepeats(seq, p).ToList();
        var pages = RepeatFinder.FormatTrfAlignmentPages(seq, found, "df", p, "df.fa");
        static string Sha(string s) => Convert.ToHexString(SHA256.HashData(Encoding.ASCII.GetBytes(s))).ToLowerInvariant();
        static int Count(string s, string what) => (s.Length - s.Replace(what, "").Length) / what.Length;

        const string page2 = "<HTML><HEAD><TITLE>df.fa.2.7.7.80.10.50.500.txt.html</TITLE></HEAD><BODY bgcolor=\"#FBF8BC\"><PRE>\n"
            + "Tandem Repeats Finder Program written by:\n"
            + "\n"
            + "                 Gary Benson\n"
            + "      Program in Bioinformatics\n"
            + "          Boston University\n"
            + "\n"
            + "Version 4.10.0\n"
            + "\n"
            + "Sequence: df\n"
            + "\n"
            + "Parameters: 2 7 7 80 10 50 500\n"
            + "\n"
            + "Pmatch=0.80,Pindel=0.10\n"
            + "tuple sizes 0,4,5,7\n"
            + "tuple distances 0, 29, 159, 500\n"
            + "\n"
            + "Length: 11400\n"
            + "ACGTcount: A:0.26, C:0.23, G:0.23, T:0.24\n"
            + "\n"
            + "Warning! 404 characters in sequence are not A, C, G, or T\n"
            + "\n"
            + "\n"
            + "File 2 of 2\n"
            + "\n"
            + "DGGGCGADAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGTAACGDGCCAGACTCTACG\n"
            + "TCGAGGAATTFGADGTGTGTTATTATTATTATTATTATTATTATTATTATTATTATATDCCAGCA\n"
            + "\n"
            + "Done.\n"
            + "</PRE></BODY></HTML>\n";

        Assert.Multiple(() =>
        {
            Assert.That(found, Has.Count.EqualTo(150));
            Assert.That(pages.Select(q => q.FileName), Is.EqualTo(new[] { "df.fa.2.7.7.80.10.50.500.1.txt.html", "df.fa.2.7.7.80.10.50.500.2.txt.html" }));
            Assert.That(Count(pages[0].Html, "Found at i:"), Is.EqualTo(4));
            Assert.That(pages[0].Html, Does.Contain("\nFound at i:276 ").And.Not.Contain("<A NAME=\"345--380,3,12.0,3,13\">"));
            Assert.That(pages[0].Html, Does.Contain("\nFGAACAGAAACFCGADCGGATAATAATAATAATAATAATAATAATAATAATAATAGGGCGAFCTG\n"));
            Assert.That(pages[1].Html, Is.EqualTo(page2));
            Assert.That(Sha(pages[0].Html), Is.EqualTo("c9a03351303444b7043c5a7b10306e3dfdf13592027abd20525a6f3b2eddc691"));
            Assert.That(Sha(pages[1].Html), Is.EqualTo("9110867d714b0c42020e0a85e249cd0c97ec62dee193f4e6e160788a2fe061f0"));
        });
    }

    [Test]
    public void AlignmentPages_RejectInvalidInput()
    {
        var p = TandemRepeatsFinderParameters.Recommended;
        var noRows = new[] { new ApproximateTandemRepeatResult(0, 8, 2, 2, "AC", 4, 100, 0, 16) };
        var badRows = new[] { noRows[0] with { AlignedSequence = "ACAC", AlignedConsensus = "ACAC" } };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => RepeatFinder.FormatTrfAlignmentPages("ACACACAC", noRows, "x", p, "x.fa"));
            Assert.Throws<ArgumentException>(() => RepeatFinder.FormatTrfAlignmentPages("ACACACAC", badRows, "x", p, "x.fa"));
            Assert.Throws<ArgumentException>(() => RepeatFinder.FormatTrfAlignmentPages("", [], "x", p, "x.fa"));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FormatTrfAlignmentPages("AC", [], "x", p, null!));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FormatTrfAlignmentPages("AC", [], "x", null!, "x.fa"));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FormatTrfAlignmentPages("AC", [], "x", p, "x.fa", -1));
        });
    }
}
