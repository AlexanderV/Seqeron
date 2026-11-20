using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// The reader count of a tree over a container (work package 1b.8) once more, for the survivor of the audit of 1b.9b: the count is striped by the managed thread id modulo 8, and a caller on any of the eight stripes must leave
/// it again, or the disposal would never find the tree idle and the mapping would stay open.
/// </summary>
[TestFixture]
public class FacadeReaderCountTests
{
    [Test]
    public void ADisposedFileBackedTree_ReleasesItsMapping_WhicheverStripeItsCallersCounted()
    {
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, StringFamilies.RandomDna(2_000, 3));
        EsaSuffixTree tree = SesaTestFile.Open(path);
        SesaFile file = tree.File!;
        var covered = new HashSet<int>();
        // managed thread ids grow with every thread that is made, so eight threads in a row cover the eight residues
        for (int attempts = 0; attempts < 400 && covered.Count < 8; attempts++)
        {
            int residue = -1;
            var thread = new Thread(() =>
            {
                residue = Environment.CurrentManagedThreadId & 7;
                Assert.That(tree.CountOccurrences("ACG"), Is.GreaterThanOrEqualTo(0));
            });
            thread.Start();
            thread.Join();
            covered.Add(residue);
        }
        Assert.That(covered, Has.Count.EqualTo(8), "a caller on every stripe");
        Assert.That(() => file.Payload(file.Parsed.Text), Throws.Nothing, "the mapping is open before the disposal");
        tree.Dispose();
        Assert.That(() => file.Payload(file.Parsed.Text), Throws.TypeOf<ObjectDisposedException>(), "the disposal found the tree idle and released the mapping");
    }

    [Test]
    public void ADisposedFileBackedTree_WithACallerInside_ReleasesItsMapping_WhenTheLastCallerLeaves()
    {
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, StringFamilies.RandomDna(2_000, 4));
        EsaSuffixTree tree = SesaTestFile.Open(path);
        SesaFile file = tree.File!;
        using var inside = new ManualResetEventSlim();
        using var go = new ManualResetEventSlim();
        Task reader = Task.Run(() => tree.Traverse(new Blocking(inside, go)));
        Assert.That(inside.Wait(TimeSpan.FromSeconds(30)), Is.True);
        tree.Dispose();
        Assert.That(() => file.Payload(file.Parsed.Text), Throws.Nothing, "a caller is inside: the mapping stays");
        go.Set();
        Assert.That(reader.Wait(TimeSpan.FromSeconds(30)), Is.True);
        Assert.That(() => file.Payload(file.Parsed.Text), Throws.TypeOf<ObjectDisposedException>(), "the last caller released the mapping when it left");
    }

    private sealed class Blocking(ManualResetEventSlim inside, ManualResetEventSlim go) : ISuffixTreeVisitor
    {
        private int _nodes;

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
        {
            if (_nodes++ == 3)
            {
                inside.Set();
                go.Wait(TimeSpan.FromSeconds(30));
            }
        }

        public void EnterBranch(int key)
        {
        }

        public void ExitBranch()
        {
        }
    }
}
