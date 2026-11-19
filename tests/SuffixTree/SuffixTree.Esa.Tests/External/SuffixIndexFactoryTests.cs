using SuffixTree.Esa.External;
using SuffixTree.Persistent;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The dispatcher <see cref="SuffixIndexFactory.Open"/> (work package 1b.7f, ADR-003 section 8): the first eight bytes decide between a SESA container and a v6 persistent tree; the v6 loader keeps its
/// exception for a container; everything else is "unknown suffix index kind".
/// </summary>
[TestFixture]
public class SuffixIndexFactoryTests
{
    private static string V6File(TempDir dir, string text, string name = "v6.suffix")
    {
        string path = dir.File(name);
        using (PersistentSuffixTree tree = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text), path)) Assert.That(tree.IsEmpty, Is.False);
        return path;
    }

    [Test]
    public void ASesaContainer_IsOpened_AsTheEsaTree()
    {
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, "mississippi");
        ISuffixTree tree = SuffixIndexFactory.Open(path);
        using var treeHandle = (IDisposable)tree;
        Assert.That(tree, Is.InstanceOf<EsaSuffixTree>());
        Assert.That(((EsaSuffixTree)tree).File!.Report.Level, Is.EqualTo(EsaValidationLevel.Header));
        Assert.That(tree.CountOccurrences("issi"), Is.EqualTo(2));
        Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo("issi"));
    }

    [Test]
    public void AV6File_IsOpened_AsThePersistentTree()
    {
        using var dir = new TempDir();
        string path = V6File(dir, "mississippi");
        ISuffixTree tree = SuffixIndexFactory.Open(path);
        using var treeHandle = (IDisposable)tree;
        Assert.That(tree, Is.InstanceOf<PersistentSuffixTree>());
        Assert.That(tree.CountOccurrences("issi"), Is.EqualTo(2));
    }

    [Test]
    public void TheOptions_AreIgnoredForAV6File_AndApplyToAContainer()
    {
        using var dir = new TempDir();
        string v6 = V6File(dir, "banana");
        ISuffixTree v6Tree = SuffixIndexFactory.Open(v6, new EsaOpenOptions { Validation = EsaValidationLevel.Structure, Access = EsaAccessMode.Trusted });
        using (var v6Handle = (IDisposable)v6Tree) Assert.That(v6Tree.CountOccurrences("an"), Is.EqualTo(2));
        string sesa = SesaTestFile.Write(dir, "banana");
        ISuffixTree esa = SuffixIndexFactory.Open(sesa, new EsaOpenOptions { Validation = EsaValidationLevel.Structure, Access = EsaAccessMode.Trusted });
        using var esaHandle = (IDisposable)esa;
        Assert.That(((EsaSuffixTree)esa).File!.Report, Has.Property("Level").EqualTo(EsaValidationLevel.Structure).And.Property("Access").EqualTo(EsaAccessMode.Trusted));
    }

    [Test]
    public void AFileOfNeitherKind_IsUnknown()
    {
        using var dir = new TempDir();
        var cases = new Dictionary<string, byte[]>
        {
            ["text"] = "this is not a suffix index at all, but it is longer than the header page".ToCharArray().Select(c => (byte)c).Concat(new byte[5000]).ToArray(),
            ["empty"] = [],
            ["seven bytes"] = [0x89, 0x53, 0x45, 0x53, 0x41, 0x0D, 0x0A],
            ["zeros"] = new byte[8192],
            ["sesa magic with a flipped byte"] = [0x89, 0x53, 0x45, 0x53, 0x41, 0x0D, 0x0A, 0x1B, 0, 0, 0, 0],
        };
        foreach ((string name, byte[] content) in cases)
        {
            string path = dir.File(name.Replace(' ', '_'));
            File.WriteAllBytes(path, content);
            Assert.That(() => SuffixIndexFactory.Open(path), Throws.TypeOf<InvalidDataException>().With.Message.Contains("unknown suffix index kind"), name);
        }
    }

    [Test]
    public void ATmpFile_WithoutMagic_IsRejected_AsUnknown()
    {
        // the writer's state between the header and the magic: the whole file is there, the first eight bytes are zero (ADR-003 section 6)
        using var dir = new TempDir();
        byte[] image = SesaTestFile.Image(LcpChildCase.Of("banana"));
        image.AsSpan(0, 8).Clear();
        string path = dir.File("x.sesa.tmp");
        File.WriteAllBytes(path, image);
        Assert.That(() => SuffixIndexFactory.Open(path), Throws.TypeOf<InvalidDataException>().With.Message.Contains("unknown suffix index kind"));
        Assert.That(() => EsaSuffixTreeFactory.Open(path), Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void TheV6Loader_StillRejectsAContainer_WithItsExistingException()
    {
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, "mississippi");
        Assert.That(() => PersistentSuffixTreeFactory.Load(path), Throws.TypeOf<InvalidOperationException>().With.Message.Contains("Magic number mismatch"));
        Assert.That(() => PersistentSuffixTreeFactory.LoadPersistent(path), Throws.TypeOf<InvalidOperationException>().With.Message.Contains("Magic number mismatch"));
    }

    [Test]
    public void EsaSuffixTreeFactoryOpen_RefusesAV6File_AsNotAContainer()
    {
        using var dir = new TempDir();
        string path = V6File(dir, "mississippi");
        Assert.That(() => EsaSuffixTreeFactory.Open(path), Throws.TypeOf<InvalidDataException>());
    }

    [Test]
    public void TheArguments_AreChecked()
    {
        Assert.That(() => SuffixIndexFactory.Open(null!), Throws.InstanceOf<ArgumentException>());
        Assert.That(() => SuffixIndexFactory.Open(string.Empty), Throws.InstanceOf<ArgumentException>());
        Assert.That(() => EsaSuffixTreeFactory.Open(null!), Throws.InstanceOf<ArgumentException>());
        using var dir = new TempDir();
        Assert.That(() => SuffixIndexFactory.Open(dir.File("missing.sesa")), Throws.TypeOf<FileNotFoundException>());
        Assert.That(() => SuffixIndexFactory.Open(dir.Path), Throws.InstanceOf<IOException>().Or.InstanceOf<UnauthorizedAccessException>());
    }

    [Test]
    public void ADeviceThatIsNotARegularFile_IsRefused_NeverHangs()
    {
        if (!OperatingSystem.IsLinux()) Assert.Ignore("Linux devices");
        foreach (string device in new[] { "/dev/null", "/dev/zero" })
        {
            if (!File.Exists(device)) continue;
            string d = device;
            var task = Task.Run(() => Assert.Catch(() => SuffixIndexFactory.Open(d)));
            Assert.That(task.Wait(TimeSpan.FromSeconds(20)), Is.True, $"{d} hangs");
            Assert.That(task.Result, Is.InstanceOf<InvalidDataException>().Or.InstanceOf<NotSupportedException>().Or.InstanceOf<IOException>(), d);
        }
    }

    [Test]
    public void TheMagicIsTheOnlyDecision_AV6FileWithADamagedHeaderIsTheV6LoadersBusiness()
    {
        using var dir = new TempDir();
        string path = V6File(dir, "mississippi");
        byte[] bytes = File.ReadAllBytes(path);
        bytes[40] ^= 0xFF; // far from the magic
        File.WriteAllBytes(path, bytes);
        try
        {
            ((IDisposable)SuffixIndexFactory.Open(path)).Dispose();
        }
        catch (InvalidDataException e)
        {
            Assert.That(e.Message, Does.Not.Contain("unknown suffix index kind"), "the v6 magic decides: a v6 file is never 'unknown'");
        }
        catch (InvalidOperationException)
        {
            // the v6 loader's own refusal of a damaged header
        }
    }

    [Test]
    public void TheDispatcher_DoesNotMapTheFileForTheDecision()
    {
        using var dir = new TempDir();
        string path = dir.File("junk.bin");
        File.WriteAllBytes(path, new byte[100_000]);
        Assert.That(() => SuffixIndexFactory.Open(path), Throws.TypeOf<InvalidDataException>());
        if (ProcMaps.Contains(path) is { } mapped) Assert.That(mapped, Is.False, "a file of an unknown kind is never mapped");
    }

    [Test]
    public void ASesaFileOpenedTwice_ThroughTheDispatcher_GivesIndependentTrees()
    {
        using var dir = new TempDir();
        string path = SesaTestFile.Write(dir, "abracadabra");
        ISuffixTree a = SuffixIndexFactory.Open(path);
        ISuffixTree b = SuffixIndexFactory.Open(path, new EsaOpenOptions { Validation = EsaValidationLevel.Sections });
        using var bHandle = (IDisposable)b;
        ((IDisposable)a).Dispose();
        Assert.That(b.CountOccurrences("abra"), Is.EqualTo(2));
        Assert.That(() => a.CountOccurrences("abra"), Throws.TypeOf<ObjectDisposedException>());
    }

    [Test]
    public void TheSesaMagic_IsTheOneOfAdr003()
    {
        Assert.That(SesaFormat.Magic.ToArray(), Is.EqualTo(new byte[] { 0x89, 0x53, 0x45, 0x53, 0x41, 0x0D, 0x0A, 0x1A }));
        Assert.That(BitConverter.ToInt64(SesaFormat.Magic), Is.EqualTo(0x1A0A0D4153455389L));
    }
}
