using Seqeron.Genomics;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

Console.WriteLine("=== DIAGNOSTIC: Actual values for weak tests ===\n");

// Test 1: PredictSignalPeptide_ClassicSignal
{
    string protein = "MKRLLLLLLLLLLLLLLLLLLASAGDDDEEEFFF";
    var s = PredictSignalPeptide(protein);
    Console.WriteLine($"SignalPeptide1: IsNull={s == null}");
    if (s != null)
    {
        Console.WriteLine($"  CleavagePosition={s.Value.CleavagePosition}");
        Console.WriteLine($"  NRegion=\"{s.Value.NRegion}\" len={s.Value.NRegion.Length}");
        Console.WriteLine($"  HRegion=\"{s.Value.HRegion}\" len={s.Value.HRegion.Length}");
        Console.WriteLine($"  CRegion=\"{s.Value.CRegion}\" len={s.Value.CRegion.Length}");
        Console.WriteLine($"  Score={s.Value.Score}");
        Console.WriteLine($"  Probability={s.Value.Probability}");
    }
}

// Test 2: PredictSignalPeptide_ReturnsRegions (different protein!)
{
    string protein = "MKRLLLLLLLLLLLLLLLLLLLASAGDDDEEEFFF";
    var s = PredictSignalPeptide(protein);
    Console.WriteLine($"\nSignalPeptide2: IsNull={s == null}");
    if (s != null)
    {
        Console.WriteLine($"  CleavagePosition={s.Value.CleavagePosition}");
        Console.WriteLine($"  NRegion=\"{s.Value.NRegion}\" len={s.Value.NRegion.Length}");
        Console.WriteLine($"  HRegion=\"{s.Value.HRegion}\" len={s.Value.HRegion.Length}");
        Console.WriteLine($"  CRegion=\"{s.Value.CRegion}\" len={s.Value.CRegion.Length}");
        Console.WriteLine($"  Score={s.Value.Score}");
    }
}

// Test 3: Transmembrane
{
    string protein = "AAAA" + new string('L', 22) + "EEEE";
    var h = PredictTransmembraneHelices(protein, windowSize: 19, threshold: 1.0).ToList();
    Console.WriteLine($"\nTM_hydrophobic: count={h.Count}");
    foreach (var (start, end, score) in h)
        Console.WriteLine($"  Start={start} End={end} Score={score:F4}");
}

// Test 4: TM multiple
{
    string loop = "EEEEEEEEEEE";
    string tm = new string('L', 22);
    string protein = tm + loop + tm + loop + tm;
    var h = PredictTransmembraneHelices(protein, threshold: 1.0).ToList();
    Console.WriteLine($"\nTM_multiple: count={h.Count}");
    foreach (var (start, end, score) in h)
        Console.WriteLine($"  Start={start} End={end} Score={score:F4}");
}

// Test 5: Disorder
{
    string protein = "LLLLVVVV" + "PEKSPEKSPPEKSPEKS" + "LLLLVVVV";
    var r = PredictDisorderedRegions(protein, threshold: 0.4).ToList();
    Console.WriteLine($"\nDisorder: count={r.Count}");
    foreach (var (start, end, score) in r)
        Console.WriteLine($"  Start={start} End={end} Score={score:F4}");
}

// Test 6: CoiledCoil
{
    string coiledCoil = "";
    for (int i = 0; i < 6; i++) coiledCoil += "LAEALEK";
    var c = PredictCoiledCoils(coiledCoil, threshold: 0.3).ToList();
    Console.WriteLine($"\nCoiledCoil: count={c.Count} seqLen={coiledCoil.Length}");
    foreach (var (start, end, score) in c)
        Console.WriteLine($"  Start={start} End={end} Score={score:F4}");
}

// Test 7: LowComplexity poly-A
{
    string protein = "MKKK" + new string('A', 15) + "VVVV";
    var r = FindLowComplexityRegions(protein, threshold: 0.5).ToList();
    Console.WriteLine($"\nLowComplex_polyA: count={r.Count}");
    foreach (var (start, end, dominant, freq) in r)
        Console.WriteLine($"  Start={start} End={end} DominantAa={dominant} Freq={freq:F4}");
}

// Test 8: LowComplexity multiple
{
    string protein = new string('G', 12) + "MKLVFP" + new string('S', 12);
    var r = FindLowComplexityRegions(protein, threshold: 0.5).ToList();
    Console.WriteLine($"\nLowComplex_multi: count={r.Count}");
    foreach (var (start, end, dominant, freq) in r)
        Console.WriteLine($"  Start={start} End={end} DominantAa={dominant} Freq={freq:F4}");
}

// Test 9: FindDomains zinc finger
{
    string protein = "AAAACXXCXXXLXXXXXXXXHXXXHAAA";
    var d = FindDomains(protein).ToList();
    Console.WriteLine($"\nDomains_ZF: count={d.Count}");
    foreach (var dm in d)
        Console.WriteLine($"  Name={dm.Name} Start={dm.Start} End={dm.End}");
}

// Test 10: FindDomains kinase
{
    string protein = "AAAAGXXXXGKSAAAA";
    var d = FindDomains(protein).ToList();
    Console.WriteLine($"\nDomains_kinase: count={d.Count}");
    foreach (var dm in d)
        Console.WriteLine($"  Name={dm.Name} Start={dm.Start} End={dm.End}");
}

// Test 11: Integration - full workflow motifs
{
    string protein = "MKRLLLLLLLLLLLLLLLLLLASAG" + "NFTAAAA" + "SARK" + "RGDAAA" + new string('L', 22) + "EEEEE";
    var motifs = FindCommonMotifs(protein).ToList();
    Console.WriteLine($"\nFullWorkflow: motifCount={motifs.Count}");
    foreach (var m in motifs)
        Console.WriteLine($"  {m.MotifName} Start={m.Start} End={m.End} Seq={m.Sequence}");
    var signal = PredictSignalPeptide(protein);
    Console.WriteLine($"  signalIsNull={signal == null}");
    if (signal != null) Console.WriteLine($"  cleavage={signal.Value.CleavagePosition}");
}

// Test 12: Large protein workflow
{
    var random = new Random(42);
    var aas = "ACDEFGHIKLMNPQRSTVWY";
    var protein = new string(Enumerable.Range(0, 500).Select(_ => aas[random.Next(aas.Length)]).ToArray());
    var motifs = FindCommonMotifs(protein).ToList();
    Console.WriteLine($"\nLargeProtein: motifCount={motifs.Count}");
    var signal = PredictSignalPeptide(protein);
    Console.WriteLine($"  signalIsNull={signal == null}");
    var tm = PredictTransmembraneHelices(protein).ToList();
    Console.WriteLine($"  tmCount={tm.Count}");
    var disorder = PredictDisorderedRegions(protein).ToList();
    Console.WriteLine($"  disorderCount={disorder.Count}");
    var domains = FindDomains(protein).ToList();
    Console.WriteLine($"  domainCount={domains.Count}");
    foreach (var dm in domains)
        Console.WriteLine($"    {dm.Name} Start={dm.Start} End={dm.End}");
}

Console.WriteLine("\n=== DONE ===");

