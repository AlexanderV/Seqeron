using System.Diagnostics;
using SuffixTree.Esa;

// The time to answer the first queries in a fresh process (work package 1b.9g, item 2: what leaving the mapped descent to tiered compilation costs): the container is opened, then <blocks> blocks of <block size> queries
// (count of the P patterns, or locate of the L patterns) are timed one by one; the line shows the total, the mean block time of the first block, of blocks 2-5, 6-10 and of every group of ten after, and the last block per query.
// usage: <container> <pattern file> <P|L> <blocks> <block size>
var tree = EsaSuffixTreeFactory.Open(args[0]);
string op = args[2];
string[] pats = File.ReadAllLines(args[1]).Where(l => l[0] == op[0]).Select(l => l[2..]).ToArray();
int blocks = int.Parse(args[3]), size = int.Parse(args[4]);
long sink = 0; double total = 0;
var times = new List<double>();
for (int b = 0; b < blocks; b++)
{
    long t0 = Stopwatch.GetTimestamp();
    for (int i = 0; i < size; i++) { string p = pats[(b * size + i) % pats.Length]; sink += op == "L" ? tree.FindAllOccurrences(p).Count : tree.CountOccurrences(p); }
    double ms = Stopwatch.GetElapsedTime(t0).TotalMilliseconds; times.Add(ms); total += ms;
}
// the line: the total, then the mean block time (ms) of the first block alone, of blocks 2-5, 6-10 and then of every group of ten, and the last block per query
var groups = new List<string> { times[0].ToString("F1") };
int at = 1;
foreach (int end in new[] { 5, 10 }) { if (at < times.Count) groups.Add(times.Skip(at).Take(end - at).Average().ToString("F2")); at = end; }
for (; at < times.Count; at += 10) groups.Add(times.Skip(at).Take(10).Average().ToString("F2"));
Console.WriteLine($"{op} total {total:F1} ms for {blocks * size} queries; mean block ms (block 1; 2-5; 6-10; then groups of ten): " + string.Join(" ", groups) + $" ... last block {times[^1] / size * 1e6:F0} ns/query  sink {sink}");
