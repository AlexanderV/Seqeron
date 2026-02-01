using System;
using System.IO;
using NUnit.Framework;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests
{
    [TestFixture]
    public class LogicalPersistenceTests
    {
        private string _tempExportFile = "test_export.bin";
        private string _tempMmfFile = "test_mmf_import.tree";

        [TearDown]
        public void Cleanup()
        {
            if (File.Exists(_tempExportFile)) File.Delete(_tempExportFile);
            if (File.Exists(_tempMmfFile)) File.Delete(_tempMmfFile);
        }

        [Test]
        public void Checksum_IsLayoutIndependent()
        {
            string text = "abracadabra";
            
            // Reference tree (Heap)
            var reference = global::SuffixTree.SuffixTree.Build(text);
            var refHash = SuffixTreeSerializer.CalculateLogicalHash(reference);
            
            // Persistent tree (Heap Storage)
            using (var heapTree = PersistentSuffixTreeFactory.Create(text) as IDisposable)
            {
                var heapHash = SuffixTreeSerializer.CalculateLogicalHash((ISuffixTree)heapTree!);
                Assert.That(heapHash, Is.EqualTo(refHash), "Heap Persistent and Reference hashes should match");
            }
            
            // Persistent tree (MMF Storage)
            using (var mmfTree = PersistentSuffixTreeFactory.Create(text, _tempMmfFile) as IDisposable)
            {
                var mmfHash = SuffixTreeSerializer.CalculateLogicalHash((ISuffixTree)mmfTree!);
                Assert.That(mmfHash, Is.EqualTo(refHash), "MMF and Reference hashes should match");
            }
        }

        [Test]
        public void DeterministicExport_IsByteIdentical()
        {
            string text = "banana";
            var st = PersistentSuffixTreeFactory.Create(text);
            
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
            var original = PersistentSuffixTreeFactory.Create(text);
            var originalHash = SuffixTreeSerializer.CalculateLogicalHash(original);
            
            using (var ms = new MemoryStream())
            {
                SuffixTreeSerializer.Export(original, ms);
                ms.Position = 0;
                
                // Import into a new storage provider
                var importedStorage = new HeapStorageProvider();
                var imported = SuffixTreeSerializer.Import(ms, importedStorage);
                
                var importedHash = SuffixTreeSerializer.CalculateLogicalHash(imported);
                Assert.That(importedHash, Is.EqualTo(originalHash), "Imported tree must have the same logical hash");
                
                // Functional parity
                Assert.That(imported.Contains("ssi"), Is.True);
                Assert.That(imported.CountOccurrences("i"), Is.EqualTo(4));
                Assert.That(imported.LongestRepeatedSubstring(), Is.EqualTo("issi"));
            }
        }
    }
}
