using System.Diagnostics;
using SuffixTree.Esa;

// time to answer the first queries in a fresh process: usage warm <container> <pattern file> <op: P|L> <blocks> <block size>
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
Console.WriteLine($"{op} total {total:F1} ms for {blocks * size} queries; per block of {size} (ms): " + string.Join(" ", times.Take(12).Select(x => x.ToString("F2"))) + " ... last " + times[^1].ToString("F2") + $" ({times[^1] / size * 1e6:F0} ns/query)  sink {sink}");
