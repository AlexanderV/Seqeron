// TRANS-PROT-001 / TRANS-SIXFRAME-001 — OrfResult.ToGetOrfRange (B02 FIN WP19, audit A2-11)
// Evidence: docs/Evidence/TRANS-SIXFRAME-001-Evidence.md
// TestSpec: tests/TestSpecs/TRANS-PROT-001.md
// Source: EMBOSS getorf.c getorf_WriteORF (1-based header [s - e], STOP excluded, reverse strand
//         s = L - s + 1, e = L - e + 1); reference values are the literal output of the EMBOSS 6.6.0
//         binary: getorf -sequence sN.fa -find 1 -table 1 -minsize 9 -outseq stdout -auto
//         (-minsize 9 nt = 3 residues = FindOrfs minLength 3; EMBOSS -table 1 = NCBI table 1 with the
//         alternative TTG/CTG starts, as GeneticCode.Standard; -reverse Y, -methionine Y).

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class Translator_OrfResultGetOrfRange_Tests
{
    // Random A/C/G/T sequences (Python random.Random(19)), L = 90, 120, 150.
    private const string S0 =
        "AACTGGCGAGTGGAGGACACATTAATAATTTGCTCACTCCCTATATTATTGTCACAATTTAGCTGTTGCTTGGAACGTATATTACTGAAC";
    private const string S1 =
        "CTGTACTATCATGAGGCGGAAGCACGTCGGGGTTACCATAAGTTGGGCTGATTAGAGTGTGAGGGACTACCTAAATACGCAGCCCGATCATGTCGACGCGTCTGATACATCTGTGTATAC";
    private const string S2 =
        "GCATTCCCTCTAGTTAACAGCTCTGTTACACTAGTGCGCAGAAAAATAAATTGGGGCAAAGCAGTGGTACACGTGACTAATCCCCTCACGACTTCATGCGCCAAGGTTCATATGGTTGATCTCATCCAAGCTTGATACATCGTTCGCGGG";

    private static IEnumerable<TestCaseData> GetorfOutputs()
    {
        // getorf 6.6.0 headers + proteins, in getorf output order.
        yield return new TestCaseData(S0, new[]
        {
            (3, 59, "MASGGHINNLLTPYIIVTI"),               // >s0_1
            (63, 89, "MLLGTYITE"),                        // >s0_2
            (49, 90, "MSQFSCCLERILLN"),                   // >s0_3, open ORF
            (22, 2, "MCPPLAS"),                           // >s0_4 (REVERSE SENSE), open ORF
        }).SetName("OrfResult_ToGetOrfRange_S0_MatchesGetorf");
        yield return new TestCaseData(S1, new[]
        {
            (1, 48, "MYYHEAEARRGYHKLG"),                  // >s1_1 (CTG start read as M)
            (11, 52, "MRRKHVGVTISWAD"),                   // >s1_2
            (48, 59, "MIRV"),                             // >s1_3
            (90, 119, "MSTRLIHLCI"),                      // >s1_4, open ORF
            (91, 74, "MIGLRI"),                           // >s1_5 (REVERSE SENSE)
            (110, 36, "MYQTRRHDRAAYLGSPSHSNQPNLW"),       // >s1_6 (REVERSE SENSE)
            (39, 1, "MVTPTCFRLMIVQ"),                     // >s1_7 (REVERSE SENSE), open ORF
        }).SetName("OrfResult_ToGetOrfRange_S1_MatchesGetorf");
        yield return new TestCaseData(S2, new[]
        {
            (23, 31, "MLH"),                              // >s2_1
            (51, 77, "MGQSSGTRD"),                        // >s2_2
            (96, 116, "MRQGSYG"),                         // >s2_3
            (112, 132, "MVDLIQA"),                        // >s2_4
            (116, 148, "MISSKLDTSFA"),                    // >s2_5
            (132, 149, "MIHRSR"),                         // >s2_6, open ORF
            (140, 111, "MYQAWMRSTI"),                     // >s2_7 (REVERSE SENSE)
            (129, 97, "MDEINHMNLGA"),                     // >s2_8 (REVERSE SENSE, TTG start)
            (97, 89, "MKS"),                              // >s2_9 (REVERSE SENSE)
            (64, 29, "MLCPNLFFCALV"),                     // >s2_10 (REVERSE SENSE)
            (104, 3, "MAHEVVRGLVTCTTALPQFIFLRTSVTELLTRGN"), // >s2_11 (REVERSE SENSE)
        }).SetName("OrfResult_ToGetOrfRange_S2_MatchesGetorf");
    }

    [TestCaseSource(nameof(GetorfOutputs))]
    public void ToGetOrfRange_LegacyFindOrfs_MatchesGetorfHeaders(string sequence, (int Start, int End, string Protein)[] expected)
    {
        var orfs = Translator.FindOrfs(new DnaSequence(sequence), minLength: 3).ToList();

        var actual = orfs.Select(o =>
        {
            var (s, e) = o.ToGetOrfRange(sequence.Length);
            return (s, e, o.Protein.Sequence);
        });

        // The legacy overload yields per frame (forward then reverse), getorf in closing order.
        Assert.That(actual, Is.EquivalentTo(expected));
    }

    [Test]
    public void ToGetOrfRange_EqualsOrfFeatureMapping_ForTheSameOrfs()
    {
        foreach (var sequence in new[] { S0, S1, S2 })
        {
            var legacy = Translator.FindOrfs(new DnaSequence(sequence), minLength: 3)
                .Select(o => (o.ToGetOrfRange(sequence.Length), o.IsTerminated, o.Protein.Sequence));
            var modes = Translator.FindOrfs(new DnaSequence(sequence), OrfFindMode.TranslationStartToStop, minSize: 9)
                .Select(o => (o.ToGetOrfRange(sequence.Length), o.IsTerminated, o.Sequence));

            Assert.That(legacy, Is.EquivalentTo(modes), sequence);
        }
    }

    [Test]
    public void IsTerminated_DistinguishesStopTerminatedAndOpenOrfs()
    {
        // s1_2 [11 - 52] ends with a STOP (stop-inclusive End = 54, 0-based); s1_4 [90 - 119] runs off the end.
        var orfs = Translator.FindOrfs(new DnaSequence(S1), minLength: 3, searchBothStrands: false).ToList();

        var terminated = orfs.Single(o => o.Protein.Sequence == "MRRKHVGVTISWAD");
        var open = orfs.Single(o => o.Protein.Sequence == "MSTRLIHLCI");

        Assert.Multiple(() =>
        {
            Assert.That((terminated.StartPosition, terminated.EndPosition, terminated.IsTerminated), Is.EqualTo((10, 54, true)));
            Assert.That((open.StartPosition, open.EndPosition, open.IsTerminated), Is.EqualTo((89, 118, false)));
        });
    }

    [TestCase(0)]
    [TestCase(-5)]
    public void ToGetOrfRange_NonPositiveLength_Throws(int length)
    {
        var orf = Translator.FindOrfs(new DnaSequence(S1), minLength: 3).First();

        Assert.Throws<ArgumentOutOfRangeException>(() => orf.ToGetOrfRange(length));
    }
}
