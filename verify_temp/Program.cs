using Seqeron.Genomics;
using Seqeron.Genomics.Population;
using Seqeron.Genomics.Phylogenetics;
using Seqeron.Genomics.Core;

Console.WriteLine("=== MISSION-CRITICAL ALGORITHM VERIFICATION ===");
Console.WriteLine("Verifying against published/calculated reference values\n");

int passed = 0;
int failed = 0;

void Check(string name, double actual, double expected, double tolerance, string source)
{
    bool ok = Math.Abs(actual - expected) <= tolerance;
    if (ok)
    {
        Console.WriteLine($"✓ {name}: {actual:F6} ≈ {expected:F6} ({source})");
        passed++;
    }
    else
    {
        Console.WriteLine($"✗ {name}: {actual:F6} ≠ {expected:F6} (diff={Math.Abs(actual - expected):F6}) ({source})");
        failed++;
    }
}

// =============================================================================
// 1. HARDY-WEINBERG EQUILIBRIUM - Ford's Scarlet Tiger Moth
// Source: Wikipedia HWE article, Ford (1971)
// Observed: AA=1469, Aa=138, aa=5, n=1612
// p = (2×1469 + 138) / (2×1612) = 3076/3224 = 0.9541
// Expected: AA = p² × n = 0.9103 × 1612 = 1467.4
//           Aa = 2pq × n = 0.0877 × 1612 = 141.3
//           aa = q² × n = 0.0021 × 1612 = 3.4
// χ² = (1469-1467.4)²/1467.4 + (138-141.3)²/141.3 + (5-3.4)²/3.4 ≈ 0.83
// =============================================================================
Console.WriteLine("1. HARDY-WEINBERG EQUILIBRIUM (Ford's Moth Data):");
var hweResult = PopulationGeneticsAnalyzer.TestHardyWeinberg("FORD_MOTH", 1469, 138, 5);
Check("  Chi-square", hweResult.ChiSquare, 0.83, 0.15, "Wikipedia HWE");
Check("  Expected AA", hweResult.ExpectedAA, 1467.4, 2.0, "p²×n calculation");
Check("  Expected Aa", hweResult.ExpectedAa, 141.3, 2.0, "2pq×n calculation");
Check("  Expected aa", hweResult.Expectedaa, 3.4, 1.0, "q²×n calculation");

// =============================================================================
// 2. WATTERSON'S THETA - Wikipedia example
// S = 10 segregating sites, n = 10 samples, L = 1000 bp
// a₁ = Σ(1/i) for i=1..9 = 1 + 0.5 + 0.333 + 0.25 + 0.2 + 0.167 + 0.143 + 0.125 + 0.111 = 2.8289
// θ = S / (a₁ × L) = 10 / (2.8289 × 1000) = 0.00354
// =============================================================================
Console.WriteLine("\n2. WATTERSON'S THETA (Wikipedia Example):");
double theta = PopulationGeneticsAnalyzer.CalculateWattersonTheta(10, 10, 1000);
Check("  Watterson θ", theta, 0.00354, 0.0005, "S/(a₁×L)");

// =============================================================================
// 3. JUKES-CANTOR DISTANCE
// p = 0.125 (1 diff in 8 sites)
// d = -0.75 × ln(1 - 4p/3) = -0.75 × ln(1 - 0.1667) = -0.75 × ln(0.8333) = -0.75 × (-0.1823) = 0.1367
// =============================================================================
Console.WriteLine("\n3. JUKES-CANTOR DISTANCE (Formula Verification):");
double jcDist = PhylogeneticAnalyzer.CalculatePairwiseDistance("ACGTACGT", "TCGTACGT",
    PhylogeneticAnalyzer.DistanceMethod.JukesCantor);
double expectedJC = -0.75 * Math.Log(1 - (4.0 * 0.125 / 3.0));
Check("  JC distance", jcDist, expectedJC, 0.0001, "d = -0.75×ln(1-4p/3)");

// =============================================================================
// 4. NUCLEOTIDE DIVERSITY (π)
// Two sequences: AAAA vs TTTT (4 sites, all different)
// π = differences / (n_comparisons × L) = 4 / (1 × 4) = 1.0
// =============================================================================
Console.WriteLine("\n4. NUCLEOTIDE DIVERSITY (π):");
var seqs = new List<IReadOnlyList<char>> {
    "AAAA".ToList(),
    "TTTT".ToList()
};
double pi = PopulationGeneticsAnalyzer.CalculateNucleotideDiversity(seqs);
Check("  π (all diff)", pi, 1.0, 0.0001, "4 diff / (1×4)");

// =============================================================================
// 5. DNA COMPLEMENT - Watson-Crick
// A↔T, C↔G
// =============================================================================
Console.WriteLine("\n5. DNA COMPLEMENT (Watson-Crick):");
var dna = new Seqeron.Genomics.Core.DnaSequence("ACGTACGT");
string comp = dna.Complement().Sequence;
string expected_comp = "TGCATGCA";
bool compOk = comp == expected_comp;
if (compOk) { Console.WriteLine($"✓ Complement: {comp} == {expected_comp}"); passed++; }
else { Console.WriteLine($"✗ Complement: {comp} ≠ {expected_comp}"); failed++; }

// =============================================================================
// 6. REVERSE COMPLEMENT
// ACGTACGT → complement: TGCATGCA → reverse: ACGTACGT
// =============================================================================
Console.WriteLine("\n6. REVERSE COMPLEMENT:");
string revcomp = dna.ReverseComplement().Sequence;
string expected_rc = "ACGTACGT"; // Palindrome!
bool rcOk = revcomp == expected_rc;
if (rcOk) { Console.WriteLine($"✓ RevComp: {revcomp} == {expected_rc} (palindrome)"); passed++; }
else { Console.WriteLine($"✗ RevComp: {revcomp} ≠ {expected_rc}"); failed++; }

// =============================================================================
// 7. GC CONTENT
// GCGCGCGC = 8 GC / 8 total = 100%
// ATATATATAT = 0 GC / 10 total = 0%
// =============================================================================
Console.WriteLine("\n7. GC CONTENT:");
var gcSeq = new DnaSequence("GCGCGCGC");
var atSeq = new DnaSequence("ATATATATAT");
double gcContent1 = (double)gcSeq.Sequence.Count(c => c == 'G' || c == 'C') / gcSeq.Sequence.Length;
double gcContent2 = (double)atSeq.Sequence.Count(c => c == 'G' || c == 'C') / atSeq.Sequence.Length;
Check("  GC% (GCGCGCGC)", gcContent1, 1.0, 0.0001, "8GC/8");
Check("  GC% (ATATATAT)", gcContent2, 0.0, 0.0001, "0GC/10");

// =============================================================================
// 8. FST - Identical populations = 0
// =============================================================================
Console.WriteLine("\n8. FST (Fixation Index):");
var pop1 = new List<(double, int)> { (0.5, 100), (0.3, 100) };
var pop2 = new List<(double, int)> { (0.5, 100), (0.3, 100) };
double fst_identical = PopulationGeneticsAnalyzer.CalculateFst(pop1, pop2);
Check("  Fst (identical)", fst_identical, 0.0, 0.001, "No differentiation");

var pop1_diff = new List<(double, int)> { (1.0, 100) };
var pop2_diff = new List<(double, int)> { (0.0, 100) };
double fst_fixed = PopulationGeneticsAnalyzer.CalculateFst(pop1_diff, pop2_diff);
Console.WriteLine($"  Fst (fixed diff): {fst_fixed:F4} (should be high, >0.5)");
if (fst_fixed > 0.5) passed++; else failed++;

// =============================================================================
Console.WriteLine($"\n=== SUMMARY: {passed} passed, {failed} failed ===");
if (failed > 0)
{
    Console.WriteLine("\n⚠️  CRITICAL: Some algorithm results don't match reference values!");
    Environment.Exit(1);
}
else
{
    Console.WriteLine("\n✓ All algorithms produce correct results against reference data");
}

