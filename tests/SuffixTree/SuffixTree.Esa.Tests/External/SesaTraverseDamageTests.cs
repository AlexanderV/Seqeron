using System.Buffers.Binary;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// A defect that the container corruption fuzz of work package 1b.8 found in the file-backed arm of <c>Traverse</c>: a suffix array whose entries are in range but not a permutation (damaged and sealed, so that neither the header
/// nor the checksums see it; only the <c>Structure</c> level does) makes <c>Traverse</c> report the same leaf twice, and <c>FindMaximalRepeatedPairs</c>, which builds linked lists of leaves by text position from that stream, closes a
/// cycle and appends pairs to a list without end (an <see cref="OutOfMemoryException"/> after many gigabytes, found as <c>SealedWords</c> of the suffix array in the fuzz). In the <c>Checked</c> mode (ADR-003 section 6: never a hang)
/// the walk must refuse a repeated leaf with <see cref="InvalidDataException"/>.
/// </summary>
[TestFixture]
public class SesaTraverseDamageTests
{
    private sealed class Counting : ISuffixTreeVisitor
    {
        internal long Nodes;

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth) => Nodes++;

        public void EnterBranch(int key)
        {
        }

        public void ExitBranch()
        {
        }
    }

    /// <summary>The image of the text with the suffix-array entry of rank <paramref name="rank"/> replaced by the entry of rank <paramref name="from"/> (two ranks now name the same text position), sealed.</summary>
    private static byte[] RepeatedEntry(string text, int rank, int from)
    {
        byte[] image = SesaTestFile.Image(LcpChildCase.Of(text));
        OracleSection sa = SesaTestFile.Section(image, SesaFormat.KindSa);
        int value = BinaryPrimitives.ReadInt32LittleEndian(image.AsSpan((int)sa.Offset + (4 * from)));
        BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan((int)sa.Offset + (4 * rank)), value);
        SesaTestFile.Reseal(image);
        return image;
    }

    private static void Open(TempDir dir, byte[] image, EsaAccessMode access, out EsaSuffixTree tree)
    {
        string path = dir.File("d.sesa");
        File.WriteAllBytes(path, image);
        tree = SesaTestFile.Open(path, EsaValidationLevel.Header, access);
    }

    private static T WithTimeout<T>(Func<T> work, int seconds = 20)
    {
        Task<T> task = Task.Run(work);
        Assert.That(Task.WaitAny([task], TimeSpan.FromSeconds(seconds)), Is.Zero, "the call does not return (the leaf lists are cyclic)");
        return task.GetAwaiter().GetResult();
    }

    [TestCase("abracadabra", 3, 4)]
    [TestCase("aaaaaaaaaaaa", 5, 6)]
    [TestCase("mississippimississippi", 9, 12)]
    [TestCase("abababababababababababab", 10, 11)]
    public void Traverse_OfASuffixArrayWithARepeatedEntry_IsInvalidData_InCheckedMode(string text, int rank, int from)
    {
        using var dir = new TempDir();
        Open(dir, RepeatedEntry(text, rank, from), EsaAccessMode.Checked, out EsaSuffixTree tree);
        using (tree)
        {
            Assert.That(() => tree.Traverse(new Counting()), Throws.TypeOf<InvalidDataException>().With.Message.Contain("suffix array"));
        }
    }

    [TestCase(3, 100)]
    [TestCase(5, int.MaxValue)]
    [TestCase(2, -1)]
    [TestCase(7, int.MinValue)]
    public void Traverse_OfASuffixArrayWithAnEntryOutsideTheText_IsInvalidData_NeverAnIndexError(int rank, int value)
    {
        // found by the corruption fuzz in the first version of the repeated-entry check (an index error of the bitmap): a position that is not in [0, n] is a damaged file as well
        byte[] image = SesaTestFile.Image(LcpChildCase.Of("abracadabra"));
        OracleSection sa = SesaTestFile.Section(image, SesaFormat.KindSa);
        BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan((int)sa.Offset + (4 * rank)), value);
        SesaTestFile.Reseal(image);
        using var dir = new TempDir();
        Open(dir, image, EsaAccessMode.Checked, out EsaSuffixTree tree);
        using (tree)
        {
            Assert.That(() => tree.Traverse(new Counting()), Throws.TypeOf<InvalidDataException>());
        }
    }

    [TestCase("a")]
    [TestCase("aaa")]
    [TestCase("aaaaaa")]
    [TestCase("abracadabra")]
    public void Traverse_OfALcpArrayWithAnyByteAtAnyRank_IsAnAnswerOrInvalidData_NeverAnotherException(string text)
    {
        // found by the corruption fuzz (SealedByteFlip and SealedWords of the LCP section): a value that makes the interval stack of LeftmostOccurrence pop its base frame ended in an IndexOutOfRangeException of a managed array
        using var dir = new TempDir();
        byte[] clean = SesaTestFile.Image(LcpChildCase.Of(text));
        OracleSection lcp = SesaTestFile.Section(clean, SesaFormat.KindLcp);
        var failures = new List<string>();
        for (int rank = 0; rank < Math.Min(lcp.Length, text.Length + 1); rank++)
        {
            for (int value = 0; value < 256; value++)
            {
                byte[] image = (byte[])clean.Clone();
                image[(int)lcp.Offset + rank] = (byte)value;
                SesaTestFile.Reseal(image);
                try
                {
                    Open(dir, image, EsaAccessMode.Checked, out EsaSuffixTree tree);
                    using (tree) tree.Traverse(new Counting());
                }
                catch (InvalidDataException)
                {
                    // allowed
                }
                catch (Exception e)
                {
                    failures.Add($"byte {value} at rank {rank}: {e.GetType().Name}");
                }
            }
        }
        Assert.That(failures, Is.Empty);
    }

    [Test]
    public void FindMaximalRepeatedPairs_OfASuffixArrayWithARepeatedEntry_IsInvalidData_NotAnEndlessList()
    {
        using var dir = new TempDir();
        Open(dir, RepeatedEntry("aaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaab", 4, 5), EsaAccessMode.Checked, out EsaSuffixTree tree);
        using (tree)
        {
            Assert.That(() => WithTimeout(() => tree.FindMaximalRepeatedPairs(3)), Throws.TypeOf<InvalidDataException>());
        }
    }

    [Test]
    public void Traverse_OfACleanFile_ReportsEveryLeafOnce_InBothModes()
    {
        string text = StringFamilies.Random(300, 3, 9) + "abababababababab";
        using EsaSuffixTree memory = EsaSuffixTree.Build(text);
        var expected = new Counting();
        memory.Traverse(expected);
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, text);
        foreach (EsaAccessMode access in new[] { EsaAccessMode.Checked, EsaAccessMode.Trusted })
        {
            using EsaSuffixTree tree = SesaTestFile.Open(path, EsaValidationLevel.Header, access);
            var counting = new Counting();
            tree.Traverse(counting);
            Assert.That(counting.Nodes, Is.EqualTo(expected.Nodes), access.ToString());
            Assert.That(tree.FindMaximalRepeatedPairs(5), Is.EqualTo(memory.FindMaximalRepeatedPairs(5)), access.ToString());
        }
    }

    [Test]
    public void ARepeatedEntry_IsRefusedAtTheStructureLevel_AsBefore()
    {
        using var dir = new TempDir();
        string path = dir.File("s.sesa");
        File.WriteAllBytes(path, RepeatedEntry("abracadabra", 3, 4));
        Assert.That(() => SesaTestFile.Open(path, EsaValidationLevel.Structure), Throws.TypeOf<InvalidDataException>());
    }
}
