using System.Diagnostics;
using System.IO.Compression;
using System.Text;
using SuffixTree;
using SuffixTree.Persistent;

// ────────────────────────────────────────────────────────────
// Human Chromosome 1 (GRCh38) — Suffix Tree Demo
// ────────────────────────────────────────────────────────────
// Downloads chr1 FASTA (~249 Mbp), builds a persistent suffix
// tree on a memory-mapped file, and runs biologically meaningful
// queries with timing for each step.
// ────────────────────────────────────────────────────────────

var dataDir = Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, "..", "..", "..", "..", "..", "data"));
var fastaGz = Path.Combine(dataDir, "chr1.fa.gz");
var treePath = Path.Combine(dataDir, "chr1.suffixtree.dat");

if (!File.Exists(fastaGz))
{
    Console.Error.WriteLine($"FASTA not found: {fastaGz}");
    Console.Error.WriteLine("Download it first:");
    Console.Error.WriteLine("  curl -L -o data/chr1.fa.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr1.fa.gz");
    return 1;
}

// ── Step 1: Parse FASTA ─────────────────────────────────────
Console.WriteLine("═══════════════════════════════════════════════════════════════");
Console.WriteLine("  Human Chromosome 1 (GRCh38) — Suffix Tree on MMF");
Console.WriteLine("═══════════════════════════════════════════════════════════════");
Console.WriteLine();

Console.Write("[1/3] Parsing FASTA (gzipped)...");
var sw = Stopwatch.StartNew();

var sb = new StringBuilder(250_000_000);
using (var fs = File.OpenRead(fastaGz))
using (var gz = new GZipStream(fs, CompressionMode.Decompress))
using (var reader = new StreamReader(gz, Encoding.ASCII))
{
    string? line;
    while ((line = reader.ReadLine()) != null)
    {
        if (line.Length == 0 || line[0] == '>') continue;
        sb.Append(line);
    }
}

var sequence = sb.ToString();
sb.Clear();
sb = null; // free ~500 MB

var parseDuration = sw.Elapsed;
Console.WriteLine($" done in {parseDuration.TotalSeconds:F1}s");
Console.WriteLine($"    Sequence length: {sequence.Length:N0} bp ({sequence.Length / 1_000_000.0:F1} Mbp)");
Console.WriteLine($"    Estimated tree:  ~{sequence.Length * 80L / 1_073_741_824.0:F1} GB");
Console.WriteLine();

// ── Step 2: Build Suffix Tree ────────────────────────────────
Console.Write($"[2/3] Building suffix tree → {Path.GetFileName(treePath)}...");
Console.Out.Flush();

// Delete old tree if exists
if (File.Exists(treePath)) File.Delete(treePath);

sw.Restart();
var textSource = new StringTextSource(sequence);

// Force GC before build to get clean memory baseline
GC.Collect(2, GCCollectionMode.Aggressive, true, true);
var memBefore = GC.GetTotalMemory(false);

var tree = PersistentSuffixTreeFactory.Create(textSource, treePath);
var buildDuration = sw.Elapsed;

var memAfter = GC.GetTotalMemory(false);
var treeFileSize = new FileInfo(treePath).Length;

Console.WriteLine($" done in {buildDuration}");
Console.WriteLine($"    Tree file size:  {treeFileSize / 1_073_741_824.0:F2} GB ({treeFileSize:N0} bytes)");
Console.WriteLine($"    Node count:      {tree.NodeCount:N0}");
Console.WriteLine($"    Leaf count:      {tree.LeafCount:N0}");
Console.WriteLine($"    Max depth:       {tree.MaxDepth:N0}");
Console.WriteLine($"    Managed heap Δ:  {(memAfter - memBefore) / 1_048_576.0:F1} MB");
Console.WriteLine();

// ── Step 3: Queries ──────────────────────────────────────────
Console.WriteLine("[3/3] Running queries on chr1 suffix tree:");
Console.WriteLine("─────────────────────────────────────────────────────────────");
Console.WriteLine();

// Helper to time a query
void TimeQuery(string label, string description, Action action)
{
    Console.Write($"  {label,-42}");
    Console.Out.Flush();
    sw.Restart();
    action();
    var elapsed = sw.Elapsed;
    Console.WriteLine($"  {elapsed.TotalMilliseconds,10:F3} ms");
    if (description.Length > 0)
        Console.WriteLine($"    → {description}");
}

// ── 3a: Pattern existence ────────────────────────────────────
Console.WriteLine("  Pattern Search (Contains):");

TimeQuery("\"GATTACA\" (movie motif)", "", () =>
{
    bool found = tree.Contains("GATTACA");
    Console.Write(found ? "FOUND" : "not found");
});

TimeQuery("\"TTAGGGTTAGGGTTAGGG\" (telomere x3)", "", () =>
{
    bool found = tree.Contains("TTAGGGTTAGGGTTAGGG");
    Console.Write(found ? "FOUND" : "not found");
});

// TP53 exon 4 partial (a real tumor suppressor gene sequence)
const string tp53Exon = "TACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAG";
TimeQuery("TP53 exon 4 (40 bp)", "", () =>
{
    bool found = tree.Contains(tp53Exon);
    Console.Write(found ? "FOUND" : "not found");
});

Console.WriteLine();

// ── 3b: Count occurrences ────────────────────────────────────
Console.WriteLine("  Count Occurrences:");

TimeQuery("\"ATG\" (start codon)", "", () =>
{
    int count = tree.CountOccurrences("ATG");
    Console.Write($"{count,10:N0}");
});

TimeQuery("\"TATAAA\" (TATA box)", "", () =>
{
    int count = tree.CountOccurrences("TATAAA");
    Console.Write($"{count,10:N0}");
});

TimeQuery("\"GAATTC\" (EcoRI restriction site)", "", () =>
{
    int count = tree.CountOccurrences("GAATTC");
    Console.Write($"{count,10:N0}");
});

TimeQuery("\"TTAGGG\" (telomeric repeat)", "", () =>
{
    int count = tree.CountOccurrences("TTAGGG");
    Console.Write($"{count,10:N0}");
});

TimeQuery("\"CAGCAG\" (Huntington CAG repeat)", "", () =>
{
    int count = tree.CountOccurrences("CAGCAG");
    Console.Write($"{count,10:N0}");
});

TimeQuery("\"AAAAAAAAAAAA\" (poly-A 12-mer)", "", () =>
{
    int count = tree.CountOccurrences("AAAAAAAAAAAA");
    Console.Write($"{count,10:N0}");
});

Console.WriteLine();

// ── 3c: Find all occurrences ─────────────────────────────────
Console.WriteLine("  Find All Occurrences:");

TimeQuery("\"GAATTC\" (EcoRI, all positions)", "", () =>
{
    var positions = tree.FindAllOccurrences("GAATTC");
    Console.Write($"{positions.Count,10:N0}");
});

TimeQuery("\"GGATCC\" (BamHI restriction site)", "", () =>
{
    var positions = tree.FindAllOccurrences("GGATCC");
    Console.Write($"{positions.Count,10:N0}");
});

Console.WriteLine();

// ── 3d: Longest repeated substring ───────────────────────────
Console.WriteLine("  Longest Repeated Substring:");

TimeQuery("LRS computation", "", () =>
{
    string lrs = tree.LongestRepeatedSubstring();
    Console.Write($"len={lrs.Length,7:N0}");
});

Console.WriteLine();

// ── 3e: Longest common substring ─────────────────────────────
Console.WriteLine("  Longest Common Substring:");

// A 200-bp probe from BRCA1 gene (chr17, so NOT on chr1 — interesting to see partial match)
const string brca1Probe =
    "ATGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCA" +
    "TCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAAC";

TimeQuery("BRCA1 probe (150 bp, from chr17)", "", () =>
{
    var lcs = tree.LongestCommonSubstring(brca1Probe);
    Console.Write($"len={lcs.Length,4}");
});

// A probe FROM chr1 — should match fully
var chr1Probe = sequence.Substring(100_000_000, 200);
TimeQuery("chr1 self-probe (200 bp @ 100M)", "", () =>
{
    var lcs = tree.LongestCommonSubstring(chr1Probe);
    Console.Write($"len={lcs.Length,4}");
});

Console.WriteLine();
Console.WriteLine("═══════════════════════════════════════════════════════════════");
Console.WriteLine($"  Total time: parse {parseDuration.TotalSeconds:F1}s + build {buildDuration} + queries");
Console.WriteLine("═══════════════════════════════════════════════════════════════");

// Cleanup
(tree as IDisposable)?.Dispose();

return 0;
