using System;
using System.IO;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests.Algorithms
{
    [TestFixture]
    [Category("Algorithms")]
    public class SuffixTreeFuzzTests
    {
        private string _tempDir;

        [SetUp]
        public void SetUp()
        {
            _tempDir = Path.Combine(Path.GetTempPath(), "SuffixTreeFuzzTests_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_tempDir);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_tempDir))
            {
                try { Directory.Delete(_tempDir, true); } catch { /* best effort */ }
            }
        }

        private string CreateValidTree()
        {
            string path = Path.Combine(_tempDir, "valid.st");
            using var tree = (IDisposable)PersistentSuffixTreeFactory.Create(new StringTextSource("ACGTACGT$"), path);
            return path;
        }

        [Test]
        public void Fuzz_Header_MagicNumber_Throws()
        {
            string path = CreateValidTree();
            CorruptFile(path, 0, new byte[] { 0x00, 0x00, 0x00, 0x00 });

            Assert.Throws<InvalidOperationException>(() =>
            {
                PersistentSuffixTreeFactory.Load(path);
            }, "Invalid storage format: Magic number mismatch.");
        }

        [Test]
        public void Fuzz_Header_VersionMismatch_Throws()
        {
            string path = CreateValidTree();
            // Corrupt version at offset 8 (magic is 8 bytes)
            CorruptFile(path, 8, BitConverter.GetBytes(9999));

            Assert.Throws<InvalidOperationException>(() =>
            {
                PersistentSuffixTreeFactory.Load(path);
            });
        }

        [Test]
        public void Fuzz_Header_TruncatedSize_Throws()
        {
            string path = CreateValidTree();
            long actualSize = new FileInfo(path).Length;
            CorruptFile(path, 40, BitConverter.GetBytes(actualSize + 5000));

            Assert.Throws<InvalidOperationException>(() =>
            {
                PersistentSuffixTreeFactory.Load(path);
            });
        }
        
        [Test]
        public void Fuzz_Header_RootOutOfBounds_Throws()
        {
            string path = CreateValidTree();
            long actualSize = new FileInfo(path).Length;
            CorruptFile(path, 16, BitConverter.GetBytes(actualSize + 1000));

            Assert.Throws<InvalidOperationException>(() =>
            {
                PersistentSuffixTreeFactory.Load(path);
            });
        }

        [Test]
        public void Fuzz_Header_TextLengthNegative_Throws()
        {
            string path = CreateValidTree();
            CorruptFile(path, 12, BitConverter.GetBytes(-50));

            Assert.Throws<InvalidOperationException>(() =>
            {
                PersistentSuffixTreeFactory.Load(path);
            });
        }

        [Test]
        public void Fuzz_Header_TextRegionOutOfBounds_Throws()
        {
            string path = CreateValidTree();
            long actualSize = new FileInfo(path).Length;
            CorruptFile(path, 24, BitConverter.GetBytes(actualSize + 1));

            Assert.Throws<InvalidOperationException>(() =>
            {
                PersistentSuffixTreeFactory.Load(path);
            });
        }

        [Test]
        public void Fuzz_Header_JumpTableInverted_Throws()
        {
            string path = CreateValidTree();
            long actualSize = new FileInfo(path).Length;
            
            // 48: TransitionOffset, 56: JumpTableStart, 64: JumpTableEnd
            CorruptFile(path, 48, BitConverter.GetBytes(100L));
            
            // Invert the jump table limits
            CorruptFile(path, 56, BitConverter.GetBytes(1000L));
            CorruptFile(path, 64, BitConverter.GetBytes(500L));

            Assert.Throws<InvalidOperationException>(() =>
            {
                PersistentSuffixTreeFactory.Load(path);
            });
        }

        [Test]
        public void Fuzz_Header_JumpTableOutOfBounds_Throws()
        {
            string path = CreateValidTree();
            long actualSize = new FileInfo(path).Length;
            
            CorruptFile(path, 48, BitConverter.GetBytes(100L)); // Transition >= 0
            CorruptFile(path, 56, BitConverter.GetBytes(actualSize + 100)); // Start OOB
            CorruptFile(path, 64, BitConverter.GetBytes(actualSize + 500)); // End OOB

            Assert.Throws<InvalidOperationException>(() =>
            {
                PersistentSuffixTreeFactory.Load(path);
            });
        }

        private void CorruptFile(string path, long offset, byte[] payload)
        {
            using var fs = new FileStream(path, FileMode.Open, FileAccess.ReadWrite);
            fs.Seek(offset, SeekOrigin.Begin);
            fs.Write(payload, 0, payload.Length);
        }
    }
}
