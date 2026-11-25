// RNA-PKPREDICT-001 / RNA-PKRECURSIVE-001 — oracle-locked rows (ViennaRNA-primitive independent oracle).
//
// Each row: sequence, ViennaRNA 2.7.2 MFE (NNDB-precision Turner 2004 parameter file, md.dangles = 2),
//           pknotsRG-over-Turner-2004 optimum of the single-knot space, of the full recursive grammar,
//           and the dot-bracket of both optima.
// The two optima are the energies of an independent Python implementation of the documented grammar
// that takes EVERY nested energy (hairpin / interior / multiloop / exterior stem / stack / helix-end mismatch)
// from ViennaRNA `E_Hairpin`, `E_IntLoop`, `E_MLstem`, `E_ExtLoop` (docs/Evidence/RNA-PKPREDICT-001-Evidence.md §6).
// Sequences containing NNDB-only terms ViennaRNA lacks (special C bulge, 1-nt bulge degeneracy, GGUC stack,
// special GU closure, all-C loops) are excluded, so the rows agree with the oracle to 0.00.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

internal static class PknotsRgOracleData
{
    // { sequence, viennaMfe, singleKnotOptimum, recursiveOptimum, singleDotBracket, recursiveDotBracket }
    internal static readonly object[][] Rows =
    {
        new object[] { "GGGGAACCCCAACCCCAAGGGG", -7.84, -13.56, -13.56, "((((..[[[[..))))..]]]]", "((((..[[[[..))))..]]]]" },
        new object[] { "GGGCGAAAGCCCAAGGGGAACCCCAACCCCAAGGGG", -15.98, -21.70, -21.70, "((((....))))..((((..[[[[..))))..]]]]", "((((....))))..((((..[[[[..))))..]]]]" },
        new object[] { "GGCGCGGCACCGUCCGCGGAACAAACGG", -7.90, -10.51, -10.51, "..(((((..[[[[)))))......]]]]", "..(((((..[[[[)))))......]]]]" },
        new object[] { "UCAAGUAUUCCGAAGCUCAACGGGAAAAUGAGCUA", -6.54, -12.31, -12.31, ".......(((((.[[[[[[.)))))...]]]]]].", ".......(((((.[[[[[[.)))))...]]]]]]." },
        new object[] { "UCAUGUUAGCUCAUGAAUGGCU", -4.22, -5.93, -5.93, "(((((..[[[[)))))..]]]]", "(((((..[[[[)))))..]]]]" },
        new object[] { "UGCCGUCGUAAGCAAAACGACUCG", -5.14, -7.52, -7.52, "(((.[[[[[..)))..]]]]]...", "(((.[[[[[..)))..]]]]]..." },
        new object[] { "ACCGCAAUGGGUGAAAGCGAUACACCCA", -7.53, -11.31, -11.31, "..(((..[[[[[[...)))...]]]]]]", "..(((..[[[[[[...)))...]]]]]]" },
        new object[] { "UCCCACCCAACGUGUGGGUGUUAUACACGGACC", -8.63, -12.94, -12.94, "...((((((.[[[[)))))).....]]]]....", "...((((((.[[[[)))))).....]]]]...." },
        new object[] { "AGACUAGCUCCGCGAGUCUCGGAGC", -8.27, -10.33, -10.33, ".((((.[[[[[...))))..]]]]]", ".((((.[[[[[...))))..]]]]]" },
        new object[] { "CGGAGGCCCAAUGCAAGUGGUUUGCCCGAAGCC", -5.69, -6.81, -6.81, "(((.[[[.....((((.....)))))))..]]]", "(((.[[[.....((((.....)))))))..]]]" },
        new object[] { "CGCACCGCGGUGGCGACCGAACGCCGCGAACGG", -12.98, -15.30, -15.30, "(((..[[[[[[[)))......]]]]]]].....", "(((..[[[[[[[)))......]]]]]]]....." },
        new object[] { "ACCCAGGGAUCCGGAAGCCCUGAACUCCGGCGU", -7.46, -14.69, -14.69, "...(((((..[[[[[..)))))...]]]]]...", "...(((((..[[[[[..)))))...]]]]]..." },
        new object[] { "UGCUGCGGCGCGUGCUAUUCGCCGUC", -8.27, -8.88, -8.88, ".((.[[[[[[...))....]]]]]].", ".((.[[[[[[...))....]]]]]]." },
        new object[] { "CAGUCCCGGGGCAUCGGGAUUCUGCAU", -10.27, -11.90, -11.90, "..(((((((.[[[)))))))..]]]..", "..(((((((.[[[)))))))..]]].." },
        new object[] { "CUGGAUGUGUCCCUAAAGCGAG", -0.96, -0.96, -0.96, "..(((....)))..........", "..(((....))).........." },
        new object[] { "CAAGCGCAGAGCUUAAUCUU", -3.38, -3.38, -3.38, ".((((.....))))......", ".((((.....))))......" },
        new object[] { "CUGCACUGCACUCGGAGUGGUU", -5.13, -5.13, -5.13, "...((((.(....).))))...", "...((((.(....).))))..." },
        new object[] { "UUAAAGCCUAACGGAUAGCCG", -1.82, -1.82, -1.82, "...........(((....)))", "...........(((....)))" },
        new object[] { "GCCCACAGACAGUUGGACGUGCAGUGCACC", -5.09, -5.09, -5.09, "..(((........)))..((((...)))).", "..(((........)))..((((...))))." },
        new object[] { "CAACAUUUUCAAUAAUGAUCGAUGG", -2.29, -2.29, -2.29, "...((((.(((....)))..)))).", "...((((.(((....)))..))))." },
        new object[] { "CCCAUCUGGUAGGGAACCCAGUCCGGAACGCUCCGUUGCG", -9.82, -11.36, -14.27, "(((..[[[[..)))...]]]]..((.((((...)))).))", "(((..[[[[..)))...]]]]..((((.[[[))))..]]]" },
        new object[] { "GGCAUCCGAGCCUCCGGAUCGGCCGCAAGCCAAGCG", -10.78, -12.28, -13.90, "(((.[[[[.)))..]]]]..(((.....))).....", "(((.[[[[.)))..]]]]..(((.[[..)))..]]." },
        new object[] { "CUCGUAUGCGACGAGUAACGCAUGCCGGAGCCCACCGGUUGGC", -14.39, -16.36, -17.48, "(((((.[[[[)))))...]]]].(((((......)))))....", "(((((.[[[[)))))...]]]]..((((.[[[..))))..]]]" },
        new object[] { "GCCAGCGGCCUUACGGGAUAGGCAUACCCGGGCAACGC", -10.96, -11.56, -12.15, "(((.[[[......((((.((....)))))))))..]]]", "(((.[[[((((..[[[[..))))...]]]])))..]]]" },
        new object[] { "CGGAGCCCCUACCGAAGGGGGUACCCGAUAAGCUCGGGAAAGCUU", -17.08, -20.10, -22.78, "(((.[[[[[[.)))..]]]]]].(((((.....))))).......", "(((.[[[[[[.)))..]]]]]].(((((.[[[[)))))...]]]]" },
        new object[] { "CGGAAGGGGUACCGAUACCCCAUUCGAGAGAGACUCGAUUCUC", -14.99, -16.12, -17.91, ".....((((((....))))))...((((.[[[[))))..]]]]", "(((..[[[[[.)))..]]]]]...((((.[[[[))))..]]]]" },
        new object[] { "CCGCUUGGAGUGCGGAAACUCCUUUCGGAAGCCGUCCGUACGGC", -17.43, -18.85, -22.38, "((((..[[[[[))))..]]]]]........(((((....)))))", "((((..[[[[[))))..]]]]]...((((.[[[[))))..]]]]" },
        new object[] { "CGCCCGGGCAAGGCGGGAAACCGCCAAGCGGGAGCCCCCCCCCGCCUGGGGGCAAGCCCGGGCG", -43.46, -43.46, -46.66, "(((((((((..(((((....)))))........((((((........))))))..)))))))))", "(((((((((..(((((....)))))..(((((.[[[[[[.)))))..]]]]]]..)))))))))" },
        new object[] { "CGCCGGGGAACCCGGGAAACCGGGAACGGGGCUGCGCCAUCCCCGUAAGGCGCAACCCCGGCG", -36.26, -36.26, -40.47, "((((((((..(((((....))))).((((((.((...)).)))))).........))))))))", "((((((((..(((((....))))).((((((.[[[[[[..))))))..]]]]]].))))))))" },
    };

    internal static IEnumerable<object[]> RowSource() => Rows;
}
