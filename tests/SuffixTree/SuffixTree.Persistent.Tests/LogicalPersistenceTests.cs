using System;
using System.IO;
using System.Linq;
using NUnit.Framework;
using SuffixTree;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests
{
    [TestFixture]
    public class LogicalPersistenceTests
    {
        private string _tempExportFile = "test_export.bin";
        private string _tempMmfFile = "test_mmf_import.tree";
        private string _tempSaveFile = "test_save.tree";

        [TearDown]
        public void Cleanup()
        {
            if (File.Exists(_tempExportFile)) File.Delete(_tempExportFile);
            if (File.Exists(_tempMmfFile)) File.Delete(_tempMmfFile);
            if (File.Exists(_tempSaveFile)) File.Delete(_tempSaveFile);
        }

        [Test]
        public void Checksum_IsLayoutIndependent()
        {
            string text = "abracadabra";

            // Reference tree (Heap)
            var reference = global::SuffixTree.SuffixTree.Build(text);
            var refHash = SuffixTreeSerializer.CalculateLogicalHash(reference);

            // Persistent tree (Heap Storage)
            using (var heapTree = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable)
            {
                var heapHash = SuffixTreeSerializer.CalculateLogicalHash((ISuffixTree)heapTree!);
                Assert.That(heapHash, Is.EqualTo(refHash), "Heap Persistent and Reference hashes should match");
            }

            // Persistent tree (MMF Storage)
            using (var mmfTree = PersistentSuffixTreeFactory.Create(new StringTextSource(text), _tempMmfFile) as IDisposable)
            {
                var mmfHash = SuffixTreeSerializer.CalculateLogicalHash((ISuffixTree)mmfTree!);
                Assert.That(mmfHash, Is.EqualTo(refHash), "MMF and Reference hashes should match");
            }
        }

        [Test]
        public void DeterministicExport_IsByteIdentical()
        {
            string text = "banana";
            var st = PersistentSuffixTreeFactory.Create(new StringTextSource(text));

            byte[] export1;
            byte[] export2;

            using (var ms1 = new MemoryStream())
            {
                SuffixTreeSerializer.Export(st, ms1);
                export1 = ms1.ToArray();
            }

            using (var ms2 = new MemoryStream())
            {
                SuffixTreeSerializer.Export(st, ms2);
                export2 = ms2.ToArray();
            }

            Assert.That(export2, Is.EqualTo(export1), "Two exports of the same tree must be byte-identical");
        }

        [Test]
        public void ExportImport_Parity()
        {
            string text = "mississippi";
            var original = PersistentSuffixTreeFactory.Create(new StringTextSource(text));
            var originalHash = SuffixTreeSerializer.CalculateLogicalHash(original);

            using (var ms = new MemoryStream())
            {
                SuffixTreeSerializer.Export(original, ms);
                ms.Position = 0;

                // Import into a new storage provider — rebuilds via Ukkonen
                var importedStorage = new HeapStorageProvider();
                var imported = SuffixTreeSerializer.Import(ms, importedStorage);

                var importedHash = SuffixTreeSerializer.CalculateLogicalHash(imported);
                Assert.That(importedHash, Is.EqualTo(originalHash), "Imported tree must have the same logical hash");

                // Functional parity
                Assert.That(imported.Contains("ssi"), Is.True);
                Assert.That(imported.CountOccurrences("i"), Is.EqualTo(4));
                Assert.That(imported.LongestRepeatedSubstring(), Is.EqualTo("issi"));

                // FindExactMatchAnchors works on imported tree (suffix links intact)
                var anchors = imported.FindExactMatchAnchors("mississippi", 3);
                Assert.That(anchors.Count, Is.GreaterThan(0),
                    "FindExactMatchAnchors must work on imported tree (suffix links rebuilt by Ukkonen)");

                var refAnchors = original.FindExactMatchAnchors("mississippi", 3);
                Assert.That(anchors.Count, Is.EqualTo(refAnchors.Count), "Anchor count must match");
            }
        }

        [Test]
        public void ExportImport_FromInMemoryTree()
        {
            string text = "abracadabra";
            var inMemory = global::SuffixTree.SuffixTree.Build(text);

            using (var ms = new MemoryStream())
            {
                SuffixTreeSerializer.Export(inMemory, ms);
                ms.Position = 0;

                var storage = new HeapStorageProvider();
                var imported = SuffixTreeSerializer.Import(ms, storage);

                Assert.That(imported.Contains("abra"), Is.True);
                Assert.That(imported.CountOccurrences("a"), Is.EqualTo(5));
                Assert.That(imported.LongestRepeatedSubstring(), Is.EqualTo("abra"));

                // Full anchor parity with in-memory tree
                string query = "cadabracadabra";
                var refAnchors = inMemory.FindExactMatchAnchors(query, 3);
                var impAnchors = imported.FindExactMatchAnchors(query, 3);
                Assert.That(impAnchors.Count, Is.EqualTo(refAnchors.Count));
                for (int i = 0; i < refAnchors.Count; i++)
                {
                    Assert.That(impAnchors[i].PositionInQuery, Is.EqualTo(refAnchors[i].PositionInQuery));
                    Assert.That(impAnchors[i].Length, Is.EqualTo(refAnchors[i].Length));
                }
            }
        }

        [Test]
        public void SaveToFile_LoadFromFile_FullFunctionality()
        {
            string text = "the quick brown fox jumps over the lazy dog";
            var inMemory = global::SuffixTree.SuffixTree.Build(text);

            // Save in-memory tree to MMF file
            using (var saved = SuffixTreeSerializer.SaveToFile(inMemory, _tempSaveFile) as IDisposable)
            {
                var savedTree = (ISuffixTree)saved!;
                Assert.That(savedTree.Contains("brown fox"), Is.True);
                Assert.That(savedTree.LongestRepeatedSubstring(), Is.EqualTo("the "));
            }

            // Load from file — must have full functionality
            using (var loaded = SuffixTreeSerializer.LoadFromFile(_tempSaveFile) as IDisposable)
            {
                var tree = (ISuffixTree)loaded!;

                Assert.That(tree.Text.ToString(), Is.EqualTo(text));
                Assert.That(tree.Contains("lazy dog"), Is.True);
                Assert.That(tree.CountOccurrences("the"), Is.EqualTo(2));
                Assert.That(tree.LongestRepeatedSubstring(), Is.EqualTo("the "));

                // FindExactMatchAnchors — the key test
                string query = "the brown lazy fox";
                var anchors = tree.FindExactMatchAnchors(query, 3);
                Assert.That(anchors.Count, Is.GreaterThan(0),
                    "FindExactMatchAnchors must work on tree loaded from MMF file");

                // Verify anchors are valid substrings
                foreach (var anchor in anchors)
                {
                    string inText = text.Substring(anchor.PositionInText, anchor.Length);
                    string inQuery = query.Substring(anchor.PositionInQuery, anchor.Length);
                    Assert.That(inText, Is.EqualTo(inQuery),
                        $"Anchor mismatch at text[{anchor.PositionInText}] vs query[{anchor.PositionInQuery}] len={anchor.Length}");
                }
            }
        }

        [Test]
        public void SaveToFile_LoadFromFile_HashParity()
        {
            string text = "repetitive-repetitive-repetitive";
            var reference = global::SuffixTree.SuffixTree.Build(text);
            var refHash = SuffixTreeSerializer.CalculateLogicalHash(reference);

            using (var saved = SuffixTreeSerializer.SaveToFile(reference, _tempSaveFile) as IDisposable)
            {
                var savedHash = SuffixTreeSerializer.CalculateLogicalHash((ISuffixTree)saved!);
                Assert.That(savedHash, Is.EqualTo(refHash), "Saved tree hash must match reference");
            }

            using (var loaded = SuffixTreeSerializer.LoadFromFile(_tempSaveFile) as IDisposable)
            {
                var loadedHash = SuffixTreeSerializer.CalculateLogicalHash((ISuffixTree)loaded!);
                Assert.That(loadedHash, Is.EqualTo(refHash), "Loaded tree hash must match reference");
            }
        }
    }
}
