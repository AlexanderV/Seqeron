using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Phylogenetics;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Phylogenetic area — phylogenetic distance calculation
/// (PHYLO-DIST-001): the pairwise evolutionary-distance routine and the symmetric
/// distance matrix built from it — and phylogenetic tree construction (PHYLO-TREE-001):
/// the distance-matrix-to-tree builders (UPGMA / Neighbor-Joining).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds malformed, out-of-domain and boundary inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang or infinite loop, no
/// nonsense output, and no *unhandled* runtime exception (DivideByZeroException, a NaN
/// or ±Infinity leaking from log-of-non-positive, IndexOutOfRangeException). Every input
/// must resolve to EITHER a well-defined, theory-correct result OR a *documented,
/// intentional* validation exception (ArgumentException / ArgumentNullException). For a
/// distance function the central theory contract is that it behaves like a (pseudo)metric:
/// d(a,a) = 0 (identity of indiscernibles), d(a,b) = d(b,a) (symmetry), and d ≥ 0
/// (non-negativity). A raw runtime exception, a hang, a NaN, or a negative distance on
/// garbage input is a bug, not a passing test. — docs/ADVANCED_TESTING_CHECKLIST.md §8
/// "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PHYLO-DIST-001 — phylogenetic distance
/// Checklist: docs/checklists/03_FUZZING.md, row 39.
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the empty sequence (no comparable sites: a
///          division-by-zero hazard on the p = differences/comparableSites step), the
///          single-sequence matrix (a degenerate 1×1 matrix with only the 0 diagonal),
///          and the saturation boundary p ≥ 0.75 where the Jukes-Cantor / Kimura
///          logarithm argument turns non-positive (log of ≤ 0 → NaN/−Inf hazard).
///   • MC = Malformed Content — non-DNA / gap / ambiguous characters in otherwise
///          aligned sequences, which the implementation skips by pairwise deletion.
/// — docs/checklists/03_FUZZING.md §Description (BE = boundary values: 0, empty;
///   MC = malformed content); row 39 targets:
///   "Identical seqs, empty seqs, single seq, non-DNA chars".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The distance contract under test (Distance_Matrix.md)
/// ───────────────────────────────────────────────────────────────────────────
/// PhylogeneticAnalyzer exposes ONE shared pairwise scan over aligned sequences behind
/// four distance models selected by the DistanceMethod enum (Distance_Matrix.md §1, §2):
///   • Hamming  — raw mismatch COUNT over comparable sites:  d_H = Σ 1[s1[i] ≠ s2[i]].
///   • PDistance — proportion of differing comparable sites:  p = differences / comparable.
///   • JukesCantor (JC69) — d = −3/4 · ln(1 − 4p/3)            (Distance_Matrix.md §2.C).
///   • Kimura2Parameter (K80) — d = −1/2 · ln((1 − 2S − V)·√(1 − 2V))  (§2.D),
///       S = transitions/comparable, V = transversions/comparable.
///
/// API entry points (Distance_Matrix.md §5.1; PhylogeneticAnalyzer.cs):
///   • CalculatePairwiseDistance(string seq1, string seq2, DistanceMethod = JukesCantor)
///     (lines 223–270) — one scalar distance. null seq → ArgumentNullException;
///     unequal lengths → ArgumentException (lines 226–229).
///   • CalculateDistanceMatrix(IReadOnlyList<string>, DistanceMethod = JukesCantor)
///     (lines 199–218) — a symmetric n×n double[,] with a 0 diagonal, filled only for
///     i &lt; j and mirrored (j,i); the diagonal is left at its default 0.0.
///
/// THE FOUR ROW-39 FUZZ TARGETS, mapped to the theory-correct contract:
///   • Identical seqs (KEY — identity of indiscernibles): two equal sequences have zero
///     differences, so p = 0 and EVERY method returns exactly 0 — Hamming 0, PDistance 0,
///     JC −3/4·ln(1) = 0, K2P −1/2·ln(1·√1) = 0 (Distance_Matrix.md §6.1 "Identical
///     comparable sequences → distance 0 for all methods"). Pinned as the metric's
///     defining property, and as the zero diagonal of the matrix.
///   • Empty seqs (BE — div-by-zero hazard): two empty (or all-gap / all-ambiguous)
///     sequences have comparableSites = 0. The implementation GUARDS this with an explicit
///     `if (comparableSites == 0) return 0;` (PhylogeneticAnalyzer.cs line 256) BEFORE the
///     p = differences/comparableSites division (line 258), so the denominator is never
///     zero — the theory-correct boundary is a finite 0, NOT a DivideByZeroException, NaN,
///     or Infinity (Distance_Matrix.md §3.3, §6.1 "no comparable sites → 0"). This is the
///     central BE probe.
///   • Single seq (BE — degenerate matrix): a 1-element sequence list yields a 1×1 matrix.
///     The fill loop `for j = i+1` never runs, so only the default-0.0 diagonal remains —
///     a well-formed 1×1 matrix whose single entry [0,0] is 0, no throw, no out-of-range
///     indexing (CalculateDistanceMatrix lines 206–215). The pairwise scalar surface
///     instead requires TWO sequences and is exercised via self-distance d(a,a) = 0.
///   • Non-DNA chars (MC — pairwise deletion): gaps ('-') and any non-A/C/G/T symbol are
///     SKIPPED at a site rather than crashing or being miscounted (lines 242–243,
///     IsStandardBase). Junk therefore lowers the comparable-site count but never throws,
///     never invents a difference, and never produces a NaN; an all-junk pair collapses to
///     the comparableSites = 0 → 0 boundary above (Distance_Matrix.md §3.3, §5.2). Case is
///     irrelevant: each site is upper-cased before inspection (line 238–239).
///
/// THE SATURATION BOUNDARY (BE — log of non-positive): for a CORRECTED model the hazard is
/// p → 0.75, where 1 − 4p/3 → 0 and the JC logarithm argument turns non-positive. The
/// source GUARDS this with `if (arg <= 0) return double.PositiveInfinity;` (line 287; the
/// K2P helper guards both of its arguments, lines 296). So a fully-saturated pair returns a
/// DEFINED +Infinity sentinel — NEVER a NaN and never a −Infinity from ln(negative)
/// (Distance_Matrix.md §6.1 "p ≥ 0.75 in JC69 → positive infinity"; INV-JC-02, INV-K2P-02).
/// +Infinity is a finite-arithmetic-safe, non-negative, intentional saturation marker, so it
/// satisfies non-negativity; we pin that it is +Infinity and specifically NOT NaN.
///
/// Documented invariants pinned (Distance_Matrix.md §2): INV-HAMMING-02 / INV-PDIST /
/// identical → 0; INV-PDIST-01 0 ≤ p ≤ 1; INV-JC-01 d_JC ≥ p; INV-JC-02 / INV-K2P-02
/// saturation → +Infinity (not NaN). The pinned exact value: for p = 1/4 (one mismatch in
/// four comparable A/C/G/T sites) JC = −3/4·ln(1 − 1/3) = −3/4·ln(2/3) ≈ 0.3040988.
/// CalculatePairwiseDistance and CalculateDistanceMatrix are pure (no iterators), so every
/// probe calls them directly; deterministic fuzz inputs use a locally fixed-seed Random
/// (no shared static RNG).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PHYLO-TREE-001 — phylogenetic tree construction
/// Checklist: docs/checklists/03_FUZZING.md, row 40.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the count boundaries of the taxon set: 0 taxa,
///          1 taxon (both below the documented 2-taxon minimum), the smallest valid set
///          (2 taxa), a large set (100+ taxa, an O(n³) completion probe), and the
///          degenerate all-identical-sequence input where every pairwise distance is 0.
/// — docs/checklists/03_FUZZING.md §Description (BE = boundary values: 0, empty);
///   row 40 targets: "0 seqs, 1 seq, 2 seqs, 100+ seqs, all identical seqs".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The tree-construction contract under test (Tree_Construction.md)
/// ───────────────────────────────────────────────────────────────────────────
/// PhylogeneticAnalyzer turns a distance matrix into a tree over the taxa via two classical
/// agglomerative builders selected by the TreeMethod enum (Tree_Construction.md §2):
///   • UPGMA — rooted ultrametric clustering; child branch lengths are clamped
///     non-negative via `Math.Max(0, newHeight − childHeight)` (INV-UPGMA-01/03, §5.2).
///   • NeighborJoining — non-clock topology; branch lengths from the standard NJ formula,
///     which MAY be negative (INV-NJ-02); the in-memory result is rooted by convention so
///     it fits the binary `PhyloNode` API, with the final unrooted centre kept as an
///     N-ary trifurcation of the last three taxa (Saitou & Nei 1987; §4.B, §5.2).
///
/// API entry points (Tree_Construction.md §5.1; PhylogeneticAnalyzer.cs):
///   • BuildTree(IReadOnlyDictionary&lt;string,string&gt; sequences, DistanceMethod, TreeMethod)
///     (lines 136–164) — computes the distance matrix then dispatches to the builder.
///     &lt; 2 sequences OR a null dictionary → ArgumentException; unaligned (unequal-length)
///     sequences → ArgumentException (lines 141–150; §3.3).
///   • BuildTreeFromMatrix(IReadOnlyList&lt;string&gt; taxa, double[,], TreeMethod)
///     (lines 174–194) — builds directly from a precomputed matrix; &lt; 2 taxa, a null
///     matrix, or a dimension mismatch → ArgumentException (lines 179–184; §3.3).
/// The returned PhylogeneticTree carries the Root `PhyloNode` (Children/IsLeaf N-ary model),
/// the Taxa list, the DistanceMatrix used, and the Method name (§3.2).
///
/// THE FIVE ROW-40 FUZZ TARGETS, mapped to the theory-correct contract:
///   • 0 seqs (BE — below the minimum): an empty dictionary / empty taxa list is &lt; 2 and is
///     REJECTED with a documented ArgumentException — never an IndexOutOfRange or a malformed
///     empty tree (BuildTree line 141; BuildTreeFromMatrix line 179; §3.3). The defined
///     contract is "throw", which we pin as the theory-correct boundary, not a crash.
///   • 1 seq (BE — below the minimum): a single-taxon input is also &lt; 2 and likewise throws
///     ArgumentException. A "trivial single-leaf tree" is NOT this implementation's contract;
///     the doc fixes the minimum at two taxa, so we VERIFY and pin the documented rejection
///     rather than inventing a one-leaf tree (§3.1, §3.3).
///   • 2 seqs (BE — smallest valid set): the minimal tree is a single internal root joining
///     two leaves. We pin: exactly 2 leaves, the root has 2 children that are both leaves
///     carrying the two taxon names, and the branch lengths are finite and non-negative
///     (UPGMA join height d/2 ≥ 0; INV-UPGMA-01/03).
///   • 100+ seqs (BE — large set, O(n³) completion): 120 SHORT deterministic sequences must
///     BUILD without hang or crash and yield a tree with exactly 120 leaves and finite, ≥0
///     UPGMA branch lengths. Guarded by `[CancelAfter]`; sequences kept short so the O(n³)
///     builder completes quickly (§4.A complexity).
///   • all-identical seqs (BE — KEY div-by-zero/NaN hazard): identical sequences ⇒ every
///     pairwise distance is 0 ⇒ the join height d/2 = 0 and every UPGMA branch length is the
///     clamped 0. This must NOT divide by zero (the n−2 NJ divisor is only reached for ≥4
///     active taxa) nor leak a NaN/Inf branch length: a valid tree with the right leaf count
///     and all-zero, finite, non-negative branch lengths (INV-UPGMA-03; §5.2).
///
/// Structural invariants pinned for every built tree (TreeBranchLengthsValid helper):
/// the leaf count equals the taxon count, and every branch length is finite (no NaN/±Inf).
/// For UPGMA we additionally pin non-negativity (INV-UPGMA-03); NJ branch lengths may be
/// negative by design (INV-NJ-02) so only finiteness is asserted there. The builders are
/// pure static methods, so every probe calls them directly with deterministic, locally
/// fixed-seed Random inputs (no shared static RNG).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PHYLO-NEWICK-001 — Newick parsing / serialization
/// Checklist: docs/checklists/03_FUZZING.md, row 41.
/// Fuzz strategies exercised for THIS unit:
///   • TF  = Truncated Fields — a Newick string cut off mid-tree: an opened descendant
///           list that is never closed ("(A,B"), a trailing colon with no number,
///           a tree truncated before its terminator.
///   • MC  = Malformed Content — structurally invalid Newick: unbalanced parentheses
///           (more "(" than ")" and vice-versa), trailing garbage after a complete tree,
///           and the empty / whitespace string.
///   • INJ = Injection — pathological / special-character payloads: a DEEPLY NESTED paren
///           stack (a recursive-descent-parser StackOverflow / DoS vector), embedded
///           NUL bytes, and non-ASCII unicode label characters.
/// — docs/checklists/03_FUZZING.md §Description (TF = truncated fields; MC = malformed
///   content; INJ = injection special chars / null bytes / unicode); row 41 targets:
///   "Malformed Newick, unbalanced parens, missing semicolon, empty string".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The Newick contract under test (Newick_Format.md; PHYLO-NEWICK-001.md)
/// ───────────────────────────────────────────────────────────────────────────
/// PhylogeneticAnalyzer is a recursive-descent Newick reader/writer over the PhyloNode tree
/// (Newick_Format.md §2.2 grammar: Tree → Subtree ";"; Internal → "(" BranchSet ")" Name;
/// BranchSet → Branch | Branch "," BranchSet; Length → empty | ":" number):
///   • ToNewick(PhyloNode, bool includeBranchLengths = true) (lines 566–608) — serializes a
///     tree, always terminating with ";" (INV-01); branch lengths rendered with "." via
///     InvariantCulture (INV-02); internal names emitted only when valid unquoted labels
///     (INV-03). A null node returns "" (Newick_Format.md §6.1; N7).
///   • ParseNewick(string) (lines 643–671) — recursive descent over "(", ",", ")", labels and
///     ":"-prefixed numbers; trims input and strips ONE trailing ";" if present; accepts an
///     optional root branch length (INV-04); any unconsumed trailing input → FormatException.
///
/// THE FOUR ROW-41 FUZZ TARGETS, mapped to the theory-correct contract (empirically verified):
///   • Malformed Newick (MC — garbage / trailing junk): a structurally broken string is
///     REJECTED with a documented parse exception, never an unhandled crash. Trailing input
///     after a complete tree ("(A,B);extra") and a "(...)" followed by an extra ")"
///     ("(A,B));") both throw FormatException (ParseNewick line 668). A bare token with no
///     parentheses ("@#$%") is — by the permissive raw-label rule (Newick_Format.md §3.3,
///     §5.2) — read as a single LEAF named verbatim, not a crash: we pin that documented
///     leniency rather than inventing a rejection the parser does not implement.
///   • Unbalanced parens (MC/TF — KEY): an opening "(" with no matching ")" ("(A,B",
///     "(((A)") throws FormatException ("unbalanced parentheses …", ParseNewickRecursive),
///     and a stray extra ")" likewise throws — DETERMINISTICALLY, never a hang and never a
///     silently-accepted truncated tree (Newick_Format.md §3.3).
///   • Missing semicolon (VERIFIED tolerance): the parser strips a trailing ";" ONLY when
///     present, so a terminator-less but otherwise valid tree ("(A,B)") parses to the SAME
///     tree as "(A,B);" — the missing ";" is TOLERATED, not rejected (ParseNewick lines
///     648–650). We pin the verified behavior, not a guessed rejection.
///   • Empty string (MC/BE): null, "", and whitespace-only input throw ArgumentException
///     ("Newick string is empty.", ParseNewick line 645) — a documented, intentional
///     rejection, never a NullReference (Newick_Format.md §3.1, §6.1; N6).
///
/// THE INJECTION / DEEP-NESTING PROBE (INJ — DoS / StackOverflow, KEY): ParseNewick is
/// recursive descent, so each level of nested parentheses consumes one call-stack frame.
/// An UNGUARDED parser overflows the stack on a deeply nested payload — an UNCATCHABLE
/// StackOverflowException that terminates the process (a denial-of-service bug). The source
/// is depth-guarded: PhylogeneticAnalyzer.MaxParseDepth caps the nesting and a payload
/// deeper than the cap is rejected with a catchable FormatException rather than overflowing
/// (ParseNewickRecursive depth guard). The fuzz test asserts that a payload FAR beyond the
/// cap (here 50 000 nested "(") is GRACEFULLY rejected — no StackOverflow, no hang under
/// [CancelAfter] — and that the deepest still-accepted depth (MaxParseDepth) parses. NOTE: a
/// StackOverflowException cannot be asserted with Throws — it kills the runner — so the test
/// stays strictly on the guarded side of the cap. FINDING: the unguarded parser was measured
/// to overflow at ≈3000 frames on a default 1 MB thread stack; the cap was set conservatively
/// below that (and the guard added to the source) precisely so this DoS vector is closed.
///
/// THE ROUND-TRIP POSITIVE SANITY (N2/N3/N4): a hand-written valid Newick string with branch
/// lengths parses to the EXACT expected leaf set, and ToNewick(ParseNewick(s)) re-parses to a
/// structurally equivalent tree with the same leaves — the serializer/parser are inverse on a
/// well-formed binary tree (Newick_Format.md §4; PHYLO-NEWICK-001.md M05–M10). ParseNewick /
/// ToNewick are pure, so every probe calls them directly; all fuzz inputs are deterministic
/// (fixed strings or locally fixed-seed Random — no shared static RNG).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PHYLO-COMP-001 — tree comparison (Robinson–Foulds distance)
/// Checklist: docs/checklists/03_FUZZING.md, row 42 (the LAST Phylogenetic-area unit).
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate tree-pair boundaries: the SAME tree
///          compared with itself (the identity-of-a-metric KEY), trees over DIFFERENT
///          leaf sets (the mismatched-taxon boundary — must not KeyNotFound/NullReference),
///          and the empty / single-leaf / null tree (no comparable clade — never crash).
/// — docs/checklists/03_FUZZING.md §Description (BE = boundary values: 0, empty);
///   row 42 targets: "Same tree vs same tree, different leaf sets, empty tree".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The tree-comparison contract under test (Tree_Comparison.md)
/// ───────────────────────────────────────────────────────────────────────────
/// PhylogeneticAnalyzer compares two trees by ROOTED Robinson–Foulds distance — the size of
/// the symmetric difference of the two trees' sets of non-trivial CLADES (Tree_Comparison.md
/// §2.A; the rooted-clade interpretation, NOT the unrooted-bipartition variant):
///   RF(T1,T2) = |S1 △ S2| = |S1 \ S2| + |S2 \ S1|,  S = clades(T) (Tree_Comparison.md §2.A).
/// A clade is the set of taxa beneath an internal node; the implementation keeps only
/// NON-TRIVIAL clades — strictly more than one taxon (not a leaf) and strictly fewer than all
/// (not the full-tree root) — each materialised as a sorted, "|"-joined taxon-name string and
/// compared as a string set (CollectClades, lines 1034–1058; §5.2).
///
/// API entry point (Tree_Comparison.md §5.1; PhylogeneticAnalyzer.cs):
///   • RobinsonFouldsDistance(PhyloNode tree1, PhyloNode tree2) (lines 868–875) — returns the
///     raw rooted RF distance as an int. Branch lengths are IGNORED (topology only; §3.3).
///
/// THE THREE ROW-42 FUZZ TARGETS, mapped to the theory-correct contract (empirically verified):
///   • Same tree vs same tree (BE — INV-RF-01, identity, KEY): clades(T) = clades(T) so the
///     symmetric difference is empty and RF(T,T) = 0 — the defining metric property
///     (Tree_Comparison.md §2.A INV-RF-01). Pinned for hand-built and Newick-parsed trees and
///     for an equal-content-but-distinct-object pair (the comparison is by clade STRING, not
///     reference identity), and alongside it INV-RF-02 symmetry and INV-RF-03 non-negativity.
///   • Different leaf sets (BE — mismatched-taxon boundary, KEY): the RAW rooted RF surface is
///     defined as a string-set symmetric difference, so trees over DIFFERENT (partially- or
///     fully-disjoint) leaf sets do NOT throw and do NOT KeyNotFound/NullReference — they yield
///     a finite, non-negative, symmetric clade-difference count (verified: a one-taxon-swapped
///     four-taxon pair → 2; fully-disjoint four-taxon trees → 4). NOTE: this differs from the
///     SEPARATE CalculateUnrootedRobinsonFoulds surface, which DOES reject mismatched leaf sets
///     with ArgumentException; the row-42 unit under test is the raw rooted RobinsonFouldsDistance,
///     whose documented behaviour is the tolerant symmetric-difference count (§3.3, §5.2,
///     §5.3.A "raw rooted RF"). We pin the VERIFIED tolerant behaviour, not a guessed rejection.
///   • Empty tree (BE — no comparable clade): a null root, a single-leaf tree, and a two-leaf
///     tree each carry NO non-trivial clade (GetLeaves(null) yields nothing; CollectClades guards
///     null; a 1-/2-leaf tree produces an empty clade set). RF therefore stays a finite 0 (or a
///     plain clade count vs a richer tree) — NEVER a NullReferenceException, KeyNotFound, or
///     DivideByZero (GetClades/CollectClades null-guard; §6.1 "no non-trivial clade → 0").
///
/// THE POSITIVE SANITY (INV-RF-01/02/05): two four-taxon trees differing in exactly ONE clade —
/// the balanced cherry ((A,B),(C,D)) vs the caterpillar (((A,B),C),D) — have RF = 2 (each tree
/// contributes one clade absent from the other), the comparison is symmetric (RF(a,b)=RF(b,a)),
/// and RF stays within the rooted-binary bound 2(n−2) = 4 for n = 4 (Tree_Comparison.md §2.A
/// INV-RF-05). RobinsonFouldsDistance is a pure static method, so every probe calls it directly
/// with deterministic, hand-built or fixed-string Newick trees (no shared static RNG).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PHYLO-BOOT-001 — phylogenetic bootstrap (Felsenstein 1985)
/// Checklist: docs/checklists/03_FUZZING.md, row 221.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries of the bootstrap: a SINGLE
///          minimal tree (the smallest 2-taxon and 3-taxon alignments — a 2-taxon input has
///          NO non-trivial clade to score, a 3-taxon input has exactly one), ZERO (and
///          negative) replicates (the denominator B; B &lt; 1 must be a documented throw, never
///          a divide-by-zero), and ALL-IDENTICAL sequences (no phylogenetic signal — a
///          zero-distance matrix that must not NaN/divide-by-zero, yet still yields defined
///          support in [0,1]).
/// — docs/checklists/03_FUZZING.md §Description (BE = boundary values: 0, -1, empty);
///   row 221 targets: "single tree, 0 replicates, identical sequences".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The bootstrap contract under test (Bootstrap_Analysis.md)
/// ───────────────────────────────────────────────────────────────────────────
/// PhylogeneticAnalyzer.Bootstrap is Felsenstein's (1985) non-parametric phylogenetic bootstrap
/// [Bootstrap_Analysis.md §2.2, refs 1–3]: build the reference tree T₀ from the original
/// alignment and record its set of non-trivial clades C(T₀); then for each of B replicates draw
/// L column indices uniformly WITH REPLACEMENT from {0..L−1}, assemble a pseudo-alignment of the
/// SAME length L keeping ALL taxa, rebuild a replicate tree Tᵣ, and for every reference clade
/// c ∈ C(T₀) count whether some clade of Tᵣ has the identical leaf-name set. The returned support
/// is the PROPORTION
///   support(c) = #{ r : c ∈ C(Tᵣ) } / B ∈ [0,1]   (Bootstrap_Analysis.md §2.2; ×100 = percentage).
/// A clade is keyed by its sorted, '|'-joined leaf names; only NON-TRIVIAL clades (more than one
/// taxon, fewer than all) are scored — exactly C(T₀) (INV-03; CollectClades).
///
/// API entry point (Bootstrap_Analysis.md §3.1, §5.1; PhylogeneticAnalyzer.cs lines 1173–1235):
///   • Bootstrap(IReadOnlyDictionary&lt;string,string&gt; sequences, int replicates = 100,
///       DistanceMethod = JukesCantor, TreeMethod = UPGMA, int seed = 42)
///     → IReadOnlyDictionary&lt;string,double&gt;  (clade key → support in [0,1]).
///     null sequences → ArgumentNullException; &lt; 2 sequences → ArgumentException;
///     replicates &lt; 1 → ArgumentException (§3.3, §6.1).
/// The `seed` parameter makes the otherwise-stochastic resampling REPRODUCIBLE (INV-04); every
/// probe passes a FIXED seed (no shared static RNG).
///
/// THE ROW-221 FUZZ TARGETS, mapped to the theory-correct contract:
///   • single tree (BE — minimal taxon set): a 2-taxon alignment has only the full-tree root
///     as an internal node, which is the TRIVIAL clade {both taxa} and is NOT scored — so the
///     support map is EMPTY, never a crash and never a spurious entry (INV-03; CollectThe
///     non-trivial filter). The smallest input that yields one scored clade is 3 taxa, where
///     exactly ONE non-trivial clade (the closest pair) is reported. Both pinned.
///   • 0 replicates (BE — KEY divide-by-zero hazard): the support denominator is B; the source
///     VALIDATES `replicates &lt; 1` and throws ArgumentException BEFORE the count/B division
///     (line 1184) — so B = 0 (and B = −1, the BE −1 probe) is a documented rejection, NEVER a
///     DivideByZeroException, a 0/0 = NaN, or an Infinity (Bootstrap_Analysis.md §3.3, §6.1).
///     The smallest VALID B = 1 is pinned to return quantized 0-or-1 support (count/1 ∈ {0,1}).
///   • identical sequences (BE — no signal, NaN / div-by-zero hazard): all-identical sequences
///     give a zero-distance matrix; every replicate's resampled pseudo-alignment is ALSO
///     all-identical (resampling identical columns yields identical columns), so each replicate
///     rebuilds the same degenerate topology. The build must not NaN or divide by zero, and the
///     returned support is DEFINED and ∈ [0,1] for every reported clade (Bootstrap_Analysis.md
///     §6.1 "all-identical sequences → every reported clade has support 1.0", INV-05).
///
/// THE POSITIVE-SANITY ANCHOR (worked example, Bootstrap_Analysis.md §7.1): the doc's hand-checkable
/// four-taxon alignment A=B="AAAAAAAAAA", C=D="GGGGGGGGGG" has two well-separated groups {A,B} and
/// {C,D}; because d(A,B)=d(C,D)=0 and A/B-vs-C/D is saturated under EVERY column resample, every
/// replicate recovers the SAME {A,B},{C,D} topology, so support["A|B"] = support["C|D"] = 1.0
/// (count = B ⇒ B/B = 1; INV-05, ref 1). Pinned verbatim against the documented expected values.
///
/// THE INVARIANTS PINNED (Bootstrap_Analysis.md §2.4): INV-01 0 ≤ support(c) ≤ 1 for every reported
/// clade; INV-02 support(c) = k/B for an integer k (quantization); INV-03 reported clades = the
/// reference tree's non-trivial clades; INV-04 fixed (alignment, B, methods, seed) ⇒ identical
/// result across runs; INV-05 a clade in every replicate has support 1.0. A randomized boundary
/// sweep over locally-seeded random alignments asserts the no-crash / [0,1] / quantization
/// contract across many shapes. The fuzz bar (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing"):
/// no crash / hang / NaN / Infinity / corruption, and the real algorithmic contract (support ∈
/// [0,1], one per non-trivial reference clade, = fraction of replicates recovering it) is pinned.
/// Bootstrap is a pure static method; every probe passes a FIXED seed (no shared static RNG) and
/// hang-sensitive / randomized-sweep probes carry [CancelAfter].
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PHYLO-STATS-001 — tree statistics (leaves, total tree length, tree height/depth)
/// Checklist: docs/checklists/03_FUZZING.md, row 222.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate TOPOLOGY boundaries of a rooted tree:
///          a SINGLE LEAF (the smallest non-empty tree — height 0, one leaf, length = its
///          own branch), a STAR tree (one internal root over all leaves as direct children
///          — height 1, the shallowest non-trivial shape), a DEEP LADDER / caterpillar (the
///          maximally IMBALANCED binary tree — height = N−1, the deepest shape on N leaves),
///          and the EMPTY (null) tree (height −1, no leaves, length 0, the empty-tree
///          convention) — plus the all-zero-branch-length boundary (length 0).
/// — docs/checklists/03_FUZZING.md §Description (BE = boundary values: 0, -1, empty);
///   row 222 targets: "single leaf, star tree, deep ladder".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The tree-statistics contract under test (Tree_Statistics.md)
/// ───────────────────────────────────────────────────────────────────────────
/// PhylogeneticAnalyzer exposes THREE exact, deterministic descriptive measures over a rooted
/// tree, each a single O(n) traversal (Tree_Statistics.md §2.2, §4.1):
///   • GetLeaves(root)            → the terminal nodes (a leaf is "a vertex with no children",
///                                   IsLeaf); pre-order order. Null root → ∅ (INV-01/02/06).
///   • CalculateTreeLength(root)  → Σ over ALL nodes of node.BranchLength — the sum-of-edge-
///                                   lengths total tree length (DendroPy Tree.length(),
///                                   Biopython total_branch_length). Null → 0; all-zero → 0
///                                   (INV-03/04/06).
///   • GetTreeDepth(root)         → the tree HEIGHT in EDGES: the number of edges on the longest
///                                   root→leaf path. A single-node (leaf) tree = 0; a null/empty
///                                   tree = −1 by the cited graph-theory convention (INV-05/06).
///
/// API entry points (Tree_Statistics.md §3, §5.1; PhylogeneticAnalyzer.cs lines 801–872):
///   • GetLeaves(PhyloNode root)            → IEnumerable&lt;PhyloNode&gt;   (yield/pre-order).
///   • CalculateTreeLength(PhyloNode root)  → double                    (recursive Σ BranchLength).
///   • GetTreeDepth(PhyloNode root)         → int                       (recursive 1 + max child).
/// All three are pure static methods; every probe calls them directly with hand-built or
/// locally-seeded-random trees (no shared static RNG).
///
/// THE ROW-222 FUZZ TARGETS, mapped to the theory-correct contract (independently hand-derived
/// from Tree_Statistics.md §2.2/§6.1 and the graph-theory definitions, NOT echoed from the code):
///   • single leaf (BE — degenerate smallest tree): a lone leaf node with no children has
///     EXACTLY one leaf (itself), height 0 ("a tree with only a single node has height zero",
///     §6.1), and total length = its OWN BranchLength (sum over the one node). No crash, no −1
///     (that is reserved for the null tree). Both the default-0-branch and a non-zero-branch
///     single leaf are pinned.
///   • star tree (BE — shallowest non-trivial shape): one internal root whose children are ALL
///     leaves directly. GetLeaves returns exactly the N children; height = 1 (root is not a leaf,
///     so 1 + max(child height = 0) = 1); total length = root.BranchLength + Σ child branch
///     lengths. Pinned for several N, including N=2 (the binary minimum) and a wide N (polytomy,
///     exercising the N-ary Children traversal of every method, not a first-two-children
///     shortcut).
///   • deep ladder / caterpillar (BE — maximally imbalanced, deepest shape): a binary caterpillar
///     on N leaves — each internal node has one leaf child and one internal child, nested N−1
///     deep. GetLeaves returns N; height = N−1 (one edge per internal level on the longest path);
///     total length = Σ of all N+(N−1) edges' branch lengths. The N−1 height is the documented
///     maximum for a rooted binary tree on N leaves and is the central deep-recursion probe
///     ([CancelAfter] guards against any traversal hang).
///   • empty / null tree (BE −1 boundary): GetLeaves(null) yields NOTHING, CalculateTreeLength(null)
///     returns 0, and GetTreeDepth(null) returns −1 — the empty-tree convention (§2.4 INV-06,
///     §6.1). NEVER a NullReferenceException. This is the −1 BE probe distinct from the leaf's 0.
///
/// THE POSITIVE-SANITY ANCHOR (worked example, Tree_Statistics.md §7.1, hand-checkable): the
/// balanced four-taxon tree ((A:1,B:1):1,(C:1,D:1):1) has leaves {A,B,C,D} (count 4), total
/// length = the six unit edges = 6.0, and height 2 (root→internal→leaf = 2 edges). Pinned verbatim
/// against the documented expected values, on a hand-built tree (independent of ParseNewick), and
/// cross-checked against a ParseNewick build of the SAME doc string.
///
/// THE INVARIANTS PINNED (Tree_Statistics.md §2.4): INV-01 every GetLeaves node IsLeaf; INV-02 an
/// N-leaf tree returns exactly N leaves; INV-03 CalculateTreeLength = Σ node.BranchLength;
/// INV-04 length ≥ 0 when all branches ≥ 0; INV-05 GetTreeDepth = edges on the longest root→leaf
/// path, leaf-only = 0; INV-06 the null triple (∅, 0, −1). A randomized boundary sweep over
/// locally-seeded random caterpillars and stars asserts the no-crash / leaf-count / height /
/// length contract across many shapes and sizes. The fuzz bar (docs/ADVANCED_TESTING_CHECKLIST.md
/// §8 "Fuzzing"): no crash / hang / NaN / Infinity / corruption, AND the real algorithmic contract
/// (the documented leaf set, total length and topological height for known topologies) is pinned to
/// EXACT values. Hang-sensitive / randomized-sweep probes carry [CancelAfter]; all random inputs use
/// a locally fixed-seed Random (no shared static RNG).
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class PhylogeneticFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    private static readonly PhylogeneticAnalyzer.DistanceMethod[] AllMethods =
    {
        PhylogeneticAnalyzer.DistanceMethod.Hamming,
        PhylogeneticAnalyzer.DistanceMethod.PDistance,
        PhylogeneticAnalyzer.DistanceMethod.JukesCantor,
        PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter,
    };

    /// <summary>Collects all leaf nodes (Children.Count == 0) of a tree by pre-order traversal.</summary>
    private static List<PhylogeneticAnalyzer.PhyloNode> CollectLeaves(PhylogeneticAnalyzer.PhyloNode root)
    {
        var leaves = new List<PhylogeneticAnalyzer.PhyloNode>();
        var stack = new Stack<PhylogeneticAnalyzer.PhyloNode>();
        stack.Push(root);
        while (stack.Count > 0)
        {
            var node = stack.Pop();
            if (node.IsLeaf)
                leaves.Add(node);
            else
                foreach (var child in node.Children)
                    stack.Push(child);
        }
        return leaves;
    }

    /// <summary>Enumerates every node of a tree (pre-order) for branch-length inspection.</summary>
    private static IEnumerable<PhylogeneticAnalyzer.PhyloNode> AllNodes(PhylogeneticAnalyzer.PhyloNode root)
    {
        var stack = new Stack<PhylogeneticAnalyzer.PhyloNode>();
        stack.Push(root);
        while (stack.Count > 0)
        {
            var node = stack.Pop();
            yield return node;
            foreach (var child in node.Children)
                stack.Push(child);
        }
    }

    /// <summary>
    /// Structural invariant assertion shared by every PHYLO-TREE-001 probe: the tree carries
    /// exactly <paramref name="expectedLeafCount"/> leaves and every branch length is finite
    /// (no NaN, no ±Infinity). When <paramref name="requireNonNegative"/> is set (UPGMA), branch
    /// lengths are additionally pinned ≥ 0 (INV-UPGMA-03); NJ is allowed negatives (INV-NJ-02).
    /// </summary>
    private static void TreeBranchLengthsValid(
        PhylogeneticAnalyzer.PhyloNode root, int expectedLeafCount, bool requireNonNegative)
    {
        CollectLeaves(root).Count.Should().Be(
            expectedLeafCount, "a tree on n taxa has exactly n leaves");

        foreach (var node in AllNodes(root))
        {
            double bl = node.BranchLength;
            double.IsNaN(bl).Should().BeFalse("a branch length must never be NaN");
            double.IsInfinity(bl).Should().BeFalse("a branch length must never be ±Infinity");
            if (requireNonNegative)
                bl.Should().BeGreaterThanOrEqualTo(0.0, "UPGMA clamps branch lengths non-negative (INV-UPGMA-03)");
        }
    }

    /// <summary>Builds an ordered, named sequence dictionary "T0..Tn" over the given sequences.</summary>
    private static Dictionary<string, string> NamedSequences(IReadOnlyList<string> sequences)
    {
        var dict = new Dictionary<string, string>();
        for (int i = 0; i < sequences.Count; i++)
            dict[$"T{i}"] = sequences[i];
        return dict;
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PHYLO-DIST-001 — phylogenetic distance : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PHYLO-DIST-001 — phylogenetic distance

    #region Positive sanity — known distance + metric properties

    /// <summary>
    /// Positive-sanity anchor: a hand-checkable pair pins the EXACT documented distance for
    /// every model, and the three metric properties (identical → 0, symmetry, non-negativity).
    /// seq1 = "AAAA", seq2 = "AAAG": one mismatch (A→G, a transition) in four comparable A/C/G/T
    /// sites ⇒ differences = 1, comparable = 4, p = 1/4, transitions = 1, transversions = 0.
    ///   • Hamming  = 1 (raw count).
    ///   • PDistance = 0.25.
    ///   • JC69     = −3/4·ln(1 − 4·0.25/3) = −3/4·ln(2/3) ≈ 0.3040988.
    ///   • K2P      = −1/2·ln((1 − 2·0.25 − 0)·√(1 − 0)) = −1/2·ln(0.5) ≈ 0.3465736.
    /// </summary>
    [Test]
    public void PairwiseDistance_KnownPair_MatchesDocumentedFormula_AndIsAMetric()
    {
        const string a = "AAAA";
        const string b = "AAAG";

        double hamming = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.Hamming);
        double pDist = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.PDistance);
        double jc = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.JukesCantor);
        double k2p = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter);

        hamming.Should().Be(1.0, "one differing site over the comparable A/C/G/T positions (INV-HAMMING-01)");
        pDist.Should().BeApproximately(0.25, 1e-12, "p = differences/comparable = 1/4 (Distance_Matrix.md §2.B)");
        jc.Should().BeApproximately(0.3040988, 1e-6, "JC69 = −3/4·ln(1 − 4p/3) at p = 1/4 (Distance_Matrix.md §2.C)");
        k2p.Should().BeApproximately(0.3465736, 1e-6, "K2P = −1/2·ln((1 − 2S − V)·√(1 − 2V)) at S = 1/4, V = 0 (§2.D)");

        // INV-JC-01: the corrected distance never under-estimates the raw proportion.
        jc.Should().BeGreaterThanOrEqualTo(pDist, "JC69 corrects upward for hidden substitutions (INV-JC-01)");

        foreach (var method in AllMethods)
        {
            double d = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, method);

            // Symmetry: d(a,b) = d(b,a).
            double reversed = PhylogeneticAnalyzer.CalculatePairwiseDistance(b, a, method);
            reversed.Should().Be(d, "phylogenetic distance is symmetric: d(a,b) = d(b,a)");

            // Non-negativity and finiteness for this non-saturated pair.
            d.Should().BeGreaterThanOrEqualTo(0.0, "an evolutionary distance is non-negative");
            double.IsNaN(d).Should().BeFalse("a well-defined pair must not produce NaN");

            // Identity of indiscernibles: d(a,a) = 0.
            double self = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, a, method);
            self.Should().Be(0.0, "identical sequences have distance 0 (identity of indiscernibles)");
        }
    }

    #endregion

    #region BE — Identical sequences (identity of indiscernibles)

    /// <summary>
    /// Identical aligned sequences ⇒ zero differences ⇒ p = 0 ⇒ distance 0 for EVERY method,
    /// both as the scalar self-distance and as the zero diagonal of the matrix
    /// (Distance_Matrix.md §6.1; INV-HAMMING-02). Uses a deterministic, locally-seeded random
    /// sequence so the identity holds for arbitrary content, not a hand-picked string.
    /// </summary>
    [Test]
    public void IdenticalSequences_AreDistanceZero_ForAllMethods()
    {
        string s = RandomDna(40, seed: 39_001);

        foreach (var method in AllMethods)
        {
            double d = PhylogeneticAnalyzer.CalculatePairwiseDistance(s, s, method);
            d.Should().Be(0.0, "identical comparable sequences differ at no site → distance 0 ({0})", method);
        }

        // The matrix diagonal is the self-distance: an identical/equal-content matrix has a 0 diagonal.
        var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(new[] { s, s, s }, PhylogeneticAnalyzer.DistanceMethod.JukesCantor);
        for (int i = 0; i < 3; i++)
            matrix[i, i].Should().Be(0.0, "the distance-matrix diagonal is the self-distance, which is 0");
        // All off-diagonal pairs are identical content ⇒ also 0.
        matrix[0, 1].Should().Be(0.0, "two equal sequences have distance 0 off the diagonal too");
        matrix[1, 0].Should().Be(0.0, "the matrix is symmetric");
    }

    #endregion

    #region BE — Empty sequences (no comparable sites → div-by-zero hazard)

    /// <summary>
    /// Two empty sequences (and the all-gap / all-ambiguous degenerate that also leaves zero
    /// comparable sites) hit the `comparableSites == 0` guard and return a finite 0 for every
    /// method — NEVER a DivideByZeroException, NaN, or Infinity from the p = differences/0 step
    /// (Distance_Matrix.md §3.3, §6.1). This is the central boundary probe.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void EmptySequences_NoComparableSites_ReturnZeroNeverNaN_ForAllMethods()
    {
        foreach (var method in AllMethods)
        {
            // Both empty: length 0 ⇒ no sites ⇒ comparableSites == 0 ⇒ guarded 0.
            Action act = () => PhylogeneticAnalyzer.CalculatePairwiseDistance(string.Empty, string.Empty, method);
            act.Should().NotThrow("zero comparable sites are guarded before the division ({0})", method);

            double empty = PhylogeneticAnalyzer.CalculatePairwiseDistance(string.Empty, string.Empty, method);
            empty.Should().Be(0.0, "no comparable sites → distance 0 (no division by zero)");
            double.IsNaN(empty).Should().BeFalse("the zero-site boundary must not produce NaN");
            double.IsInfinity(empty).Should().BeFalse("the zero-site boundary must not produce Infinity");

            // All-gap pair: every site is skipped by pairwise deletion ⇒ also zero comparable sites.
            double allGap = PhylogeneticAnalyzer.CalculatePairwiseDistance("----", "----", method);
            allGap.Should().Be(0.0, "an all-gap pair leaves no comparable site → guarded 0");
        }
    }

    #endregion

    #region BE — Single sequence (degenerate 1×1 matrix)

    /// <summary>
    /// A single-element sequence list yields a well-formed 1×1 matrix whose only entry is the
    /// default-0.0 diagonal: the i&lt;j fill loop never runs, so no division and no out-of-range
    /// indexing occur (CalculateDistanceMatrix). The scalar surface instead requires two
    /// sequences, so a "single sequence" is exercised there as the self-distance d(a,a) = 0.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void SingleSequence_MatrixIsOneByOneZero_AndSelfDistanceIsZero()
    {
        string s = RandomDna(25, seed: 39_002);

        foreach (var method in AllMethods)
        {
            Action act = () => PhylogeneticAnalyzer.CalculateDistanceMatrix(new[] { s }, method);
            act.Should().NotThrow("a 1×1 matrix has no pair to compare; the fill loop never runs ({0})", method);

            var matrix = PhylogeneticAnalyzer.CalculateDistanceMatrix(new[] { s }, method);
            matrix.GetLength(0).Should().Be(1, "one taxon → a 1×1 matrix");
            matrix.GetLength(1).Should().Be(1);
            matrix[0, 0].Should().Be(0.0, "the lone diagonal entry is the self-distance, 0");

            // The single-sequence notion on the scalar surface: distance to itself is 0.
            double self = PhylogeneticAnalyzer.CalculatePairwiseDistance(s, s, method);
            self.Should().Be(0.0, "a single sequence compared to itself has distance 0");
        }
    }

    #endregion

    #region MC — Non-DNA / gap / ambiguous characters (pairwise deletion)

    /// <summary>
    /// Gaps and non-A/C/G/T symbols are SKIPPED at a site (pairwise deletion), never crash and
    /// never invent a difference: the distance equals that of the cleaned, comparable-only
    /// sub-alignment. Here only positions 0,1,2 are comparable A/C/G/T (one mismatch at index 2),
    /// while positions 3+ carry junk in at least one row and are excluded ⇒ p = 1/3 regardless of
    /// the junk (Distance_Matrix.md §3.3, §5.2). Case-insensitivity is pinned alongside.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void NonDnaCharacters_AreSkipped_NoCrashNoNaN_NoSpuriousDifference()
    {
        //               idx: 0123456
        const string a = "ACG-N*x";
        const string b = "ACT9 Zq";
        // Comparable A/C/G/T-vs-A/C/G/T positions: 0 (A=A), 1 (C=C), 2 (G≠T). Index 3+ excluded.
        // ⇒ differences = 1, comparable = 3, p = 1/3.

        foreach (var method in AllMethods)
        {
            Action act = () => PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, method);
            act.Should().NotThrow("non-DNA/gap symbols are skipped, never crash ({0})", method);

            double d = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, method);
            double.IsNaN(d).Should().BeFalse("pairwise-deleted junk must not produce NaN ({0})", method);
            d.Should().BeGreaterThanOrEqualTo(0.0, "distance stays non-negative under pairwise deletion");
        }

        // Exact pins on the uncorrected models prove junk is excluded, not miscounted.
        PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.Hamming)
            .Should().Be(1.0, "exactly one comparable mismatch (index 2); all junk positions are deleted");
        PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.PDistance)
            .Should().BeApproximately(1.0 / 3.0, 1e-12, "p = 1 difference / 3 comparable sites");

        // An all-junk pair collapses to the zero-comparable-sites boundary → 0, not a crash.
        PhylogeneticAnalyzer.CalculatePairwiseDistance("****", "----", PhylogeneticAnalyzer.DistanceMethod.JukesCantor)
            .Should().Be(0.0, "no comparable site survives pairwise deletion → guarded 0");

        // Case-insensitive: lower-case bases compare equal to their upper-case counterparts.
        PhylogeneticAnalyzer.CalculatePairwiseDistance("acgt", "ACGT", PhylogeneticAnalyzer.DistanceMethod.Hamming)
            .Should().Be(0.0, "each site is upper-cased before inspection (case-insensitive)");
    }

    #endregion

    #region BE — Saturation boundary (log of non-positive → +Infinity, not NaN)

    /// <summary>
    /// The corrected-model saturation boundary: when p ≥ 0.75 the JC argument 1 − 4p/3 ≤ 0 and
    /// the K2P argument 1 − 2S − V ≤ 0. The source guards both and returns +Infinity rather than
    /// taking ln of a non-positive value (Distance_Matrix.md §6.1; INV-JC-02, INV-K2P-02). The
    /// theory-correct sentinel is +Infinity — specifically NOT NaN and NOT −Infinity — while the
    /// uncorrected models stay finite (PDistance saturates at 1, Hamming counts).
    /// </summary>
    [Test]
    public void SaturatedPair_CorrectedModelsReturnPositiveInfinity_NeverNaN()
    {
        // All four sites differ as transitions ⇒ p = 1.0 (≥ 0.75), S = 1.0, V = 0.0.
        const string a = "AAAA";
        const string b = "GGGG";

        double jc = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.JukesCantor);
        double k2p = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.Kimura2Parameter);

        double.IsPositiveInfinity(jc).Should().BeTrue("JC69 with p ≥ 0.75 has a non-positive log argument → +Infinity (INV-JC-02)");
        double.IsNaN(jc).Should().BeFalse("the saturation guard must return +Infinity, never NaN from ln(≤0)");
        double.IsPositiveInfinity(k2p).Should().BeTrue("K2P with a non-positive log argument → +Infinity (INV-K2P-02)");
        double.IsNaN(k2p).Should().BeFalse("the K2P saturation guard must return +Infinity, never NaN");

        // The uncorrected models stay finite at saturation.
        double p = PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.PDistance);
        p.Should().Be(1.0, "every site differs → p saturates at 1, finite (INV-PDIST-01 upper bound)");
        PhylogeneticAnalyzer.CalculatePairwiseDistance(a, b, PhylogeneticAnalyzer.DistanceMethod.Hamming)
            .Should().Be(4.0, "all four comparable sites differ");
    }

    #endregion

    #region Validation — null and unequal-length inputs (documented exceptions)

    /// <summary>
    /// The two documented validation throws on the scalar surface: a null sequence →
    /// ArgumentNullException; unequal lengths (an unaligned pair) → ArgumentException
    /// (PhylogeneticAnalyzer.cs lines 226–229; Distance_Matrix.md §6.1). These are
    /// INTENTIONAL, contract-defined rejections — not raw runtime crashes.
    /// </summary>
    [Test]
    public void NullOrUnequalLength_ThrowDocumentedValidationExceptions()
    {
        Action nullFirst = () => PhylogeneticAnalyzer.CalculatePairwiseDistance(null!, "ACGT");
        Action nullSecond = () => PhylogeneticAnalyzer.CalculatePairwiseDistance("ACGT", null!);
        nullFirst.Should().Throw<ArgumentNullException>("a null sequence is an explicit, documented rejection");
        nullSecond.Should().Throw<ArgumentNullException>("a null sequence is an explicit, documented rejection");

        Action unequal = () => PhylogeneticAnalyzer.CalculatePairwiseDistance("ACG", "ACGT");
        unequal.Should().Throw<ArgumentException>("pairwise distance requires aligned (equal-length) sequences");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PHYLO-TREE-001 — tree construction : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PHYLO-TREE-001 — tree construction

    #region Positive sanity — known topology (closest pair grouped)

    /// <summary>
    /// Positive-sanity anchor: four taxa with a hand-controlled distance structure build a tree
    /// with the EXACTLY-expected leaf count and the closest pair grouped together. Taxa A and B
    /// are nearly identical (one mismatch) while C and D diverge sharply, so UPGMA must first
    /// merge {A,B}: the deepest internal node grouping A and B forms before any other join, and
    /// A and B end up as siblings under a common ancestor that does NOT contain C or D.
    /// Branch lengths are pinned finite and non-negative (INV-UPGMA-03).
    /// </summary>
    [Test]
    public void BuildTree_KnownTopology_GroupsClosestPair_AndHasExpectedLeafCount()
    {
        // A vs B: 1 mismatch (tail). A/B vs C/D: many mismatches. C vs D: a few mismatches.
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "AAAAAAAAAA",
            ["B"] = "AAAAAAAAAG", // one mismatch from A → A,B are the closest pair
            ["C"] = "TTTTTTTTTT",
            ["D"] = "TTTTTTTGCC", // diverges from C, but far from A,B
        };

        var tree = PhylogeneticAnalyzer.BuildTree(
            sequences,
            PhylogeneticAnalyzer.DistanceMethod.PDistance,
            PhylogeneticAnalyzer.TreeMethod.UPGMA);

        tree.Root.Should().NotBeNull("a valid four-taxon input must build a tree");
        TreeBranchLengthsValid(tree.Root, expectedLeafCount: 4, requireNonNegative: true);

        // The closest pair {A,B} must appear together under one internal node that excludes C,D.
        var abClade = AllNodes(tree.Root)
            .FirstOrDefault(n => !n.IsLeaf
                && n.Taxa.Contains("A") && n.Taxa.Contains("B")
                && !n.Taxa.Contains("C") && !n.Taxa.Contains("D"));
        abClade.Should().NotBeNull("UPGMA merges the closest pair {A,B} first, forming an {A,B}-only clade");
    }

    #endregion

    #region BE — 0 sequences (below the 2-taxon minimum → documented throw)

    /// <summary>
    /// An empty taxon set is below the documented two-taxon minimum and is REJECTED with an
    /// ArgumentException on BOTH surfaces — never an IndexOutOfRange/NullReference crash and never
    /// a malformed empty tree (BuildTree line 141; BuildTreeFromMatrix line 179; §3.3). The
    /// null-dictionary degenerate is folded into the same guard.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void ZeroSequences_ThrowDocumentedArgumentException_NeverCrashNeverTree()
    {
        var empty = new Dictionary<string, string>();
        Action buildFromSeqs = () => PhylogeneticAnalyzer.BuildTree(empty);
        buildFromSeqs.Should().Throw<ArgumentException>("zero sequences is below the 2-taxon minimum (§3.3)");

        Action buildFromNull = () => PhylogeneticAnalyzer.BuildTree(null!);
        buildFromNull.Should().Throw<ArgumentException>("a null sequence dictionary is rejected by the same guard");

        var noTaxa = Array.Empty<string>();
        var emptyMatrix = new double[0, 0];
        Action buildFromMatrix = () => PhylogeneticAnalyzer.BuildTreeFromMatrix(noTaxa, emptyMatrix);
        buildFromMatrix.Should().Throw<ArgumentException>("zero taxa is below the 2-taxon minimum on the matrix surface (§3.3)");
    }

    #endregion

    #region BE — 1 sequence (below the 2-taxon minimum → documented throw)

    /// <summary>
    /// A single-taxon input is also below the two-taxon minimum and throws ArgumentException —
    /// the implementation does NOT manufacture a trivial single-leaf tree; the doc fixes the
    /// minimum at two taxa (§3.1, §3.3). VERIFIED behavior: rejection, never a crash, never a
    /// one-leaf tree. Pinned on both the sequence and the precomputed-matrix surfaces.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void SingleSequence_ThrowsDocumentedArgumentException_NotASingleLeafTree()
    {
        var one = new Dictionary<string, string> { ["Only"] = RandomDna(20, seed: 40_001) };
        Action buildFromSeqs = () => PhylogeneticAnalyzer.BuildTree(one);
        buildFromSeqs.Should().Throw<ArgumentException>("one sequence is below the documented 2-taxon minimum (§3.3)");

        var oneTaxon = new[] { "Only" };
        var oneByOne = new double[1, 1]; // a well-formed 1×1 zero matrix is still < 2 taxa
        Action buildFromMatrix = () => PhylogeneticAnalyzer.BuildTreeFromMatrix(oneTaxon, oneByOne);
        buildFromMatrix.Should().Throw<ArgumentException>("one taxon is below the minimum on the matrix surface (§3.3)");
    }

    #endregion

    #region BE — 2 sequences (smallest valid tree: one root, two leaves)

    /// <summary>
    /// The smallest valid input — two taxa — builds the minimal tree: a single internal root
    /// joining two leaves carrying the two taxon names. Pinned: exactly 2 leaves, the root is
    /// internal with two leaf children, the leaf names are the two taxa, and branch lengths are
    /// finite and non-negative (UPGMA join height d/2 ≥ 0; INV-UPGMA-01/03).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void TwoSequences_BuildMinimalTwoLeafTree_JoinedAtOneRoot()
    {
        var sequences = new Dictionary<string, string>
        {
            ["X"] = RandomDna(30, seed: 40_002),
            ["Y"] = RandomDna(30, seed: 40_003),
        };

        var tree = PhylogeneticAnalyzer.BuildTree(
            sequences,
            PhylogeneticAnalyzer.DistanceMethod.JukesCantor,
            PhylogeneticAnalyzer.TreeMethod.UPGMA);

        tree.Root.Should().NotBeNull();
        tree.Root.IsLeaf.Should().BeFalse("two taxa join under a single internal root");
        tree.Root.Children.Count.Should().Be(2, "the root bifurcates into the two taxa");

        var leaves = CollectLeaves(tree.Root);
        leaves.Should().HaveCount(2, "a two-taxon tree has exactly two leaves");
        leaves.All(l => l.IsLeaf).Should().BeTrue();
        leaves.Select(l => l.Name).Should().BeEquivalentTo(new[] { "X", "Y" }, "both taxa appear as leaves");

        TreeBranchLengthsValid(tree.Root, expectedLeafCount: 2, requireNonNegative: true);
    }

    #endregion

    #region BE — 100+ sequences (O(n³) completion under [CancelAfter])

    /// <summary>
    /// A large taxon set (120 SHORT deterministic sequences) must BUILD to completion without
    /// hanging or crashing and yield a structurally valid tree: exactly 120 leaves with finite,
    /// non-negative UPGMA branch lengths (§4.A O(n³) complexity). Sequences are kept short so the
    /// cubic builder finishes quickly; `[CancelAfter]` fails the test if it ever hangs.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void ManySequences_BuildCompletes_WithCorrectLeafCount_AndFiniteBranches()
    {
        const int n = 120;
        var sequences = new Dictionary<string, string>();
        for (int i = 0; i < n; i++)
            sequences[$"T{i}"] = RandomDna(40, seed: 40_100 + i);

        PhylogeneticAnalyzer.PhylogeneticTree tree = default;
        Action act = () => tree = PhylogeneticAnalyzer.BuildTree(
            sequences,
            PhylogeneticAnalyzer.DistanceMethod.JukesCantor,
            PhylogeneticAnalyzer.TreeMethod.UPGMA);
        act.Should().NotThrow("a large but valid input must build without crashing");

        tree.Root.Should().NotBeNull("the O(n³) builder completes for 120 taxa");
        TreeBranchLengthsValid(tree.Root, expectedLeafCount: n, requireNonNegative: true);
    }

    #endregion

    #region BE — All-identical sequences (every distance 0 → div-by-zero / NaN hazard)

    /// <summary>
    /// KEY div-by-zero / NaN probe: when every sequence is identical, every pairwise distance is 0,
    /// so each UPGMA join height d/2 = 0 and every branch length is the clamped 0. The build must
    /// NOT divide by zero (the n−2 NJ divisor is only reached for ≥4 active taxa) and must NOT leak
    /// a NaN/±Inf branch length — it produces a valid tree with the right leaf count and all-zero,
    /// finite, non-negative branch lengths (INV-UPGMA-03; §5.2). Exercised for BOTH builders;
    /// UPGMA additionally pins every branch length is exactly 0.
    /// </summary>
    [Test]
    [CancelAfter(10_000)]
    public void AllIdenticalSequences_BuildValidZeroLengthTree_NoNaNNoDivByZero()
    {
        const int n = 8;
        string identical = RandomDna(30, seed: 40_500);
        var sequences = new Dictionary<string, string>();
        for (int i = 0; i < n; i++)
            sequences[$"T{i}"] = identical;

        // UPGMA: all branch lengths clamp to exactly 0, no NaN, no division by zero.
        PhylogeneticAnalyzer.PhylogeneticTree upgma = default;
        Action buildUpgma = () => upgma = PhylogeneticAnalyzer.BuildTree(
            sequences, PhylogeneticAnalyzer.DistanceMethod.JukesCantor, PhylogeneticAnalyzer.TreeMethod.UPGMA);
        buildUpgma.Should().NotThrow("zero distances must not divide by zero or crash the builder");
        TreeBranchLengthsValid(upgma.Root, expectedLeafCount: n, requireNonNegative: true);
        AllNodes(upgma.Root).Select(node => node.BranchLength)
            .Should().OnlyContain(bl => bl == 0.0, "identical taxa → every UPGMA branch length is exactly 0");

        // Neighbor-Joining: also completes with finite branch lengths on the all-zero matrix
        // (the n−2 divisor is ≥ 2 while joining, and the 3-/2-taxon closures are pure arithmetic).
        PhylogeneticAnalyzer.PhylogeneticTree nj = default;
        Action buildNj = () => nj = PhylogeneticAnalyzer.BuildTree(
            sequences, PhylogeneticAnalyzer.DistanceMethod.JukesCantor, PhylogeneticAnalyzer.TreeMethod.NeighborJoining);
        buildNj.Should().NotThrow("an all-zero distance matrix must not crash Neighbor-Joining");
        // NJ branch lengths may be negative by design (INV-NJ-02) — only finiteness is required.
        TreeBranchLengthsValid(nj.Root, expectedLeafCount: n, requireNonNegative: false);
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PHYLO-NEWICK-001 — Newick parsing / serialization : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PHYLO-NEWICK-001 — Newick parsing

    /// <summary>Sorted leaf names of a parsed tree, for order-independent leaf-set comparison.</summary>
    private static List<string> LeafNames(PhylogeneticAnalyzer.PhyloNode root) =>
        CollectLeaves(root).Select(l => l.Name).OrderBy(n => n, StringComparer.Ordinal).ToList();

    #region Positive sanity — valid Newick parses to the right leaf set and round-trips

    /// <summary>
    /// Positive-sanity anchor: a hand-written, well-formed Newick string with branch lengths
    /// parses to the EXACT expected leaf set, and ToNewick(ParseNewick(s)) re-parses to a
    /// structurally equivalent tree (same leaves) — the parser and serializer are inverse on a
    /// valid binary tree (Newick_Format.md §4; PHYLO-NEWICK-001.md M05–M10, N2/N3/N4).
    /// </summary>
    [Test]
    public void Newick_ValidString_ParsesToLeafSet_AndRoundTrips()
    {
        const string s = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);";

        var tree = PhylogeneticAnalyzer.ParseNewick(s);
        tree.Should().NotBeNull("a well-formed Newick string must parse");

        LeafNames(tree).Should().Equal(new[] { "A", "B", "C", "D" },
            "the four leaf Name productions are A,B,C,D (Newick_Format.md §2.2; N2)");

        // The structure is two cherries: {A,B} and {C,D} each under one internal node.
        var ab = AllNodes(tree).FirstOrDefault(n => !n.IsLeaf
            && n.Taxa.Contains("A") && n.Taxa.Contains("B")
            && !n.Taxa.Contains("C") && !n.Taxa.Contains("D"));
        ab.Should().NotBeNull("(A,B) forms an internal node that excludes C,D");

        // Round-trip: serialize then re-parse → same leaf set, structurally equivalent (N3/N4).
        string serialized = PhylogeneticAnalyzer.ToNewick(tree);
        serialized.Should().EndWith(";", "ToNewick always terminates with ';' (INV-01; N1)");

        var reparsed = PhylogeneticAnalyzer.ParseNewick(serialized);
        LeafNames(reparsed).Should().Equal(LeafNames(tree),
            "ToNewick→ParseNewick preserves the leaf set (round-trip N4)");
    }

    #endregion

    #region MC/BE — Empty / whitespace / null string (documented ArgumentException)

    /// <summary>
    /// The empty-input boundary: null, the empty string, and whitespace-only input each hit the
    /// `string.IsNullOrWhiteSpace` guard and throw ArgumentException ("Newick string is empty.")
    /// — a documented, intentional rejection, NEVER a NullReferenceException or an empty/corrupt
    /// tree (Newick_Format.md §3.1, §6.1; PHYLO-NEWICK-001.md N6).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void EmptyOrNull_ThrowsDocumentedArgumentException_NeverCrash()
    {
        foreach (string bad in new[] { null!, "", "   ", "\t", "\n", "  \t \n " })
        {
            Action act = () => PhylogeneticAnalyzer.ParseNewick(bad);
            act.Should().Throw<ArgumentException>(
                "null/empty/whitespace Newick is an explicit, documented rejection (N6) — got [{0}]",
                bad is null ? "null" : "whitespace");
        }
    }

    #endregion

    #region MC/TF — Unbalanced parentheses (deterministic FormatException, no hang)

    /// <summary>
    /// Unbalanced parentheses — an opening "(" with no matching ")" (truncated descendant list)
    /// OR a stray extra ")" — are REJECTED deterministically with a FormatException, never a hang,
    /// an infinite loop, or a silently-accepted truncated tree (ParseNewickRecursive "unbalanced
    /// parentheses" guard; Newick_Format.md §3.3). [CancelAfter] pins termination. The positive
    /// control "(A,B);" parses, proving the rejection is specific to the imbalance.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void UnbalancedParentheses_RejectedWithFormatException_NoHang()
    {
        // Too many "(" (truncated / never closed) and too many ")" (excess close).
        string[] unbalanced =
        {
            "(A,B",        // opened list never closed (truncated field)
            "(((A)",       // three opens, one close
            "((A,B)",      // missing outer close
            "(A,B));",     // one extra close after a complete subtree
            "A)",          // a leaf then a stray close
            "(A,(B,C);",   // inner group never closed before ';'
        };

        foreach (string s in unbalanced)
        {
            Action act = () => PhylogeneticAnalyzer.ParseNewick(s);
            act.Should().Throw<FormatException>(
                "unbalanced parentheses are a malformed tree and must be rejected, not hang [{0}]", s);
        }

        // Positive control: the balanced counterpart parses cleanly.
        Action balanced = () => PhylogeneticAnalyzer.ParseNewick("(A,B);");
        balanced.Should().NotThrow("a balanced tree parses — the rejection is specific to imbalance");
    }

    #endregion

    #region MC — Malformed content / trailing garbage (documented parse exception)

    /// <summary>
    /// Malformed content that is NOT a well-formed tree is rejected with a documented parse
    /// exception rather than an unhandled crash: trailing input after a complete tree
    /// ("(A,B);extra") and a ":"-branch-length with no number ("(A:,B);"→ trailing junk) throw
    /// FormatException (ParseNewick trailing-input guard). A bare token with no parentheses
    /// ("@#$%") is — per the parser's permissive raw-label rule (Newick_Format.md §3.3, §5.2) —
    /// read as a single LEAF named verbatim; we pin that documented leniency (no crash, no NaN),
    /// not a rejection the parser does not implement.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void MalformedContent_RejectedOrLeniencyDocumented_NeverUnhandledCrash()
    {
        // Trailing garbage after a complete, terminated tree → FormatException.
        foreach (string s in new[] { "(A,B);extra", "(A,B);(C,D);", "(A,B):0.1 junk" })
        {
            Action act = () => PhylogeneticAnalyzer.ParseNewick(s);
            act.Should().Throw<FormatException>(
                "unconsumed trailing input after a complete tree is malformed [{0}]", s);
        }

        // A bare token with no structure: documented leniency → a single leaf named verbatim.
        var bare = PhylogeneticAnalyzer.ParseNewick("@#$%");
        bare.IsLeaf.Should().BeTrue("a bare token with no parentheses is parsed as a single leaf");
        bare.Name.Should().Be("@#$%", "leaf labels are read as a raw run (Newick_Format.md §5.2)");

        // Every probe above either threw a documented FormatException or returned a leaf — no
        // unhandled exception type leaks; assert that no malformed input causes a non-parse crash.
        foreach (string s in new[] { "()", "(,)", "(A,)", ",", ":::" })
        {
            Action act = () => PhylogeneticAnalyzer.ParseNewick(s);
            act.Should().NotThrow<NullReferenceException>("malformed Newick must not NullReference [{0}]", s);
            act.Should().NotThrow<IndexOutOfRangeException>("malformed Newick must not IndexOutOfRange [{0}]", s);
            act.Should().NotThrow<ArgumentOutOfRangeException>("malformed Newick must not throw AOORE [{0}]", s);
        }
    }

    #endregion

    #region VERIFIED — Missing semicolon is tolerated (not rejected)

    /// <summary>
    /// VERIFIED tolerance: the parser strips a trailing ";" only when present, so a valid tree
    /// WITHOUT a terminator ("(A,B)") parses to the SAME tree as the terminated form ("(A,B);")
    /// — the missing semicolon is TOLERATED, not rejected (ParseNewick trims and conditionally
    /// strips ";"; Newick_Format.md §3.3). We pin the verified behavior rather than guessing a
    /// rejection; the leaf set is identical with and without the terminator.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void MissingSemicolon_IsTolerated_SameTreeAsTerminated()
    {
        var withSemi = PhylogeneticAnalyzer.ParseNewick("((A,B),(C,D));");
        var noSemi = PhylogeneticAnalyzer.ParseNewick("((A,B),(C,D))");

        LeafNames(noSemi).Should().Equal(new[] { "A", "B", "C", "D" },
            "a terminator-less tree still yields the full leaf set");
        LeafNames(noSemi).Should().Equal(LeafNames(withSemi),
            "the missing ';' is tolerated — same tree as the terminated form");

        // A bare leaf without a terminator is likewise tolerated as a single leaf.
        var leaf = PhylogeneticAnalyzer.ParseNewick("Taxon1");
        leaf.IsLeaf.Should().BeTrue("a bare leaf token without ';' parses as one leaf");
        leaf.Name.Should().Be("Taxon1");
    }

    #endregion

    #region INJ — Deep nesting (StackOverflow / DoS) + null bytes + unicode

    /// <summary>
    /// KEY injection probe — recursive-descent StackOverflow / DoS: ParseNewick recurses once per
    /// nested "(", so an unbounded depth would overflow the call stack (an UNCATCHABLE
    /// StackOverflowException that kills the process). The source is depth-guarded
    /// (PhylogeneticAnalyzer.MaxParseDepth): a payload deeper than the cap is rejected with a
    /// catchable FormatException instead of overflowing. We assert a payload FAR beyond the cap
    /// (50 000 nested "(") is GRACEFULLY rejected — no StackOverflow, no hang under [CancelAfter]
    /// — and that the deepest still-accepted depth (MaxParseDepth) parses to a single leaf.
    /// (A StackOverflowException cannot be asserted with Throws, so the test stays on the guarded
    /// side of the cap; the guard was added precisely because the unguarded parser overflowed at
    /// ≈3000 frames on a default 1 MB thread stack.)
    /// </summary>
    [Test]
    [CancelAfter(15_000)]
    public void DeepNesting_RejectedGracefully_NoStackOverflowNoHang()
    {
        int cap = PhylogeneticAnalyzer.MaxParseDepth;

        // FAR beyond the cap: must be rejected with a catchable FormatException, never overflow.
        int deep = Math.Max(cap * 5, 50_000);
        string overDeep = new string('(', deep) + "A" + new string(')', deep) + ";";
        Action act = () => PhylogeneticAnalyzer.ParseNewick(overDeep);
        act.Should().Throw<FormatException>(
            "nesting beyond MaxParseDepth ({0}) must be rejected gracefully, never StackOverflow", cap);

        // Just over the cap → also rejected (boundary of the guard).
        string justOver = new string('(', cap + 1) + "A" + new string(')', cap + 1) + ";";
        Action overByOne = () => PhylogeneticAnalyzer.ParseNewick(justOver);
        overByOne.Should().Throw<FormatException>("depth cap + 1 is over the limit");

        // At the cap → still accepted (the deepest legal depth parses to its single leaf).
        string atCap = new string('(', cap) + "A" + new string(')', cap) + ";";
        var atCapTree = PhylogeneticAnalyzer.ParseNewick(atCap);
        CollectLeaves(atCapTree).Select(l => l.Name).Should().Equal(new[] { "A" },
            "nesting exactly at MaxParseDepth is the deepest accepted tree (single leaf A)");
    }

    /// <summary>
    /// Injection of special characters into LABEL positions: embedded NUL bytes and non-ASCII
    /// unicode are read as raw label characters (the permissive raw-label rule, Newick_Format.md
    /// §5.2) — they must not crash, hang, or corrupt the structural parse. The leaf COUNT and
    /// tree shape stay correct; only the label text carries the injected bytes.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void InjectedSpecialCharacters_InLabels_NoCrashNoCorruption()
    {
        // NUL byte inside a label, plus a unicode-named pair — both must parse to 2 leaves.
        var nul = PhylogeneticAnalyzer.ParseNewick("(A\0X,B);");
        CollectLeaves(nul).Should().HaveCount(2, "a NUL byte in a label is a label char, not a delimiter");

        var unicode = PhylogeneticAnalyzer.ParseNewick("(αβγ,δεζ);");
        var names = CollectLeaves(unicode).Select(l => l.Name).ToList();
        names.Should().HaveCount(2, "unicode-labelled tree parses to two leaves");
        names.Should().Contain("αβγ").And.Contain("δεζ", "non-ASCII label text is preserved verbatim");

        // A large injected-character label run must not hang or crash.
        string longLabel = new string('Z', 50_000);
        Action act = () => PhylogeneticAnalyzer.ParseNewick($"({longLabel},B);");
        act.Should().NotThrow("a very long label run is read linearly, no crash or hang");
    }

    #endregion

    #region INV — ToNewick null node returns empty; terminator invariant

    /// <summary>
    /// Serializer-side defensive contract: ToNewick(null) returns the empty string (Newick_Format.md
    /// §6.1; N7), and every non-null serialization ends with ";" (INV-01; N1). Pinned alongside the
    /// fuzz probes so the writer's null/edge behavior is locked next to the reader's.
    /// </summary>
    [Test]
    public void ToNewick_NullNode_ReturnsEmpty_AndOutputAlwaysTerminated()
    {
        PhylogeneticAnalyzer.ToNewick(null!).Should().Be("",
            "a null root serializes to the empty string (N7)");

        var tree = PhylogeneticAnalyzer.ParseNewick("(A:0.1,B:0.2);");
        PhylogeneticAnalyzer.ToNewick(tree).Should().EndWith(";",
            "every serialized tree ends with ';' (INV-01; N1)");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PHYLO-COMP-001 — tree comparison (Robinson–Foulds) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PHYLO-COMP-001 — tree comparison

    #region Positive sanity — one-clade-difference RF + symmetry + bound

    /// <summary>
    /// Positive-sanity anchor: two four-taxon trees differing in EXACTLY one clade — the balanced
    /// cherry ((A,B),(C,D)) versus the caterpillar (((A,B),C),D) — have rooted RF = 2 (each tree
    /// contributes one non-trivial clade absent from the other). The comparison is symmetric
    /// (INV-RF-02) and stays within the rooted-binary upper bound 2(n−2) = 4 for n = 4 (INV-RF-05;
    /// Tree_Comparison.md §2.A). The shared clade {A,B} cancels in the symmetric difference, so only
    /// the differing clades — {A,B,C} in the caterpillar and {C,D} in the cherry — are counted.
    /// </summary>
    [Test]
    public void RobinsonFoulds_OneCladeDifference_IsTwo_AndSymmetric_WithinBound()
    {
        var cherry = PhylogeneticAnalyzer.ParseNewick("((A,B),(C,D));");
        var caterpillar = PhylogeneticAnalyzer.ParseNewick("(((A,B),C),D);");

        int rf = PhylogeneticAnalyzer.RobinsonFouldsDistance(cherry, caterpillar);
        rf.Should().Be(2, "the two trees differ in exactly one clade each → symmetric difference 2");

        // INV-RF-02: symmetry.
        PhylogeneticAnalyzer.RobinsonFouldsDistance(caterpillar, cherry)
            .Should().Be(rf, "rooted RF is symmetric: RF(T1,T2) = RF(T2,T1) (INV-RF-02)");

        // INV-RF-03 / INV-RF-05: non-negative and within the rooted-binary maximum 2(n−2).
        rf.Should().BeGreaterThanOrEqualTo(0, "RF is a symmetric-difference count → non-negative (INV-RF-03)");
        const int n = 4;
        rf.Should().BeLessThanOrEqualTo(2 * (n - 2), "RF ≤ 2(n−2) for rooted binary trees on n taxa (INV-RF-05)");

        // INV-RF-04: even for binary trees with the same leaf set.
        (rf % 2).Should().Be(0, "RF is even for binary trees over the same leaf set (INV-RF-04)");
    }

    #endregion

    #region BE — Same tree vs same tree (identity, INV-RF-01, KEY)

    /// <summary>
    /// KEY identity probe (INV-RF-01): a tree compared with itself has identical clade sets, so the
    /// symmetric difference is empty and RF(T,T) = 0 (Tree_Comparison.md §2.A). Pinned for a
    /// hand-built tree, the SAME reference re-used, and — crucially — two DISTINCT objects with
    /// equal content (RF compares clade STRINGS, not reference identity), so a structurally equal
    /// re-parse also yields 0. Symmetry and non-negativity are pinned alongside.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void SameTreeVsSameTree_DistanceIsZero_AndSymmetric()
    {
        var tree = PhylogeneticAnalyzer.ParseNewick("(((A,B),C),(D,E));");

        // Same reference: trivially identical clade sets.
        PhylogeneticAnalyzer.RobinsonFouldsDistance(tree, tree)
            .Should().Be(0, "a tree compared with itself has an empty clade symmetric difference (INV-RF-01)");

        // Two DISTINCT objects of equal content: comparison is by clade string, not by reference.
        var twin = PhylogeneticAnalyzer.ParseNewick("(((A,B),C),(D,E));");
        int rf = PhylogeneticAnalyzer.RobinsonFouldsDistance(tree, twin);
        rf.Should().Be(0, "equal-content but distinct tree objects still have identical clade sets → 0");

        // INV-RF-02 symmetry, INV-RF-03 non-negativity on the identity pair.
        PhylogeneticAnalyzer.RobinsonFouldsDistance(twin, tree)
            .Should().Be(rf, "RF is symmetric even on the identity pair (INV-RF-02)");
        rf.Should().BeGreaterThanOrEqualTo(0, "RF is non-negative (INV-RF-03)");
    }

    #endregion

    #region BE — Different leaf sets (mismatched-taxon boundary, KEY)

    /// <summary>
    /// KEY mismatched-taxon boundary: the RAW rooted <see cref="PhylogeneticAnalyzer.RobinsonFouldsDistance"/>
    /// is defined as a clade-STRING symmetric difference, so trees over DIFFERENT leaf sets do NOT
    /// throw, NOT KeyNotFound, and NOT NullReference — they yield a finite, non-negative, symmetric
    /// clade-difference count (Tree_Comparison.md §5.2, §5.3.A "raw rooted RF"). VERIFIED values:
    /// a one-taxon-swapped four-taxon pair ((A,B),(C,D)) vs ((A,B),(C,E)) → 2 (the {C,D}/{C,E} clades
    /// differ); fully-disjoint trees ((A,B),(C,D)) vs ((W,X),(Y,Z)) → 4 (no clade in common). This
    /// is the documented RAW rooted contract — distinct from the SEPARATE
    /// CalculateUnrootedRobinsonFoulds surface, which intentionally REJECTS mismatched leaf sets.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void DifferentLeafSets_NoCrash_FiniteSymmetricNonNegativeCount()
    {
        var baseTree = PhylogeneticAnalyzer.ParseNewick("((A,B),(C,D));");
        var oneSwapped = PhylogeneticAnalyzer.ParseNewick("((A,B),(C,E));"); // D → E: partial overlap
        var disjoint = PhylogeneticAnalyzer.ParseNewick("((W,X),(Y,Z));");   // no taxon in common

        foreach (var other in new[] { oneSwapped, disjoint })
        {
            Action act = () => PhylogeneticAnalyzer.RobinsonFouldsDistance(baseTree, other);
            act.Should().NotThrow("the raw rooted RF tolerates mismatched leaf sets (no KeyNotFound/NullReference)");

            int rf = PhylogeneticAnalyzer.RobinsonFouldsDistance(baseTree, other);
            rf.Should().BeGreaterThanOrEqualTo(0, "a clade symmetric-difference count is non-negative (INV-RF-03)");
            PhylogeneticAnalyzer.RobinsonFouldsDistance(other, baseTree)
                .Should().Be(rf, "RF stays symmetric even across mismatched leaf sets (INV-RF-02)");
        }

        // Exact pins prove the differing clades are counted, not miscounted or skipped.
        PhylogeneticAnalyzer.RobinsonFouldsDistance(baseTree, oneSwapped)
            .Should().Be(2, "the partially-overlapping {C,D} vs {C,E} clades each differ → 2");
        PhylogeneticAnalyzer.RobinsonFouldsDistance(baseTree, disjoint)
            .Should().Be(4, "fully-disjoint leaf sets share no clade → all clades differ");
    }

    #endregion

    #region BE — Empty / single-leaf / null tree (no comparable clade → never crash)

    /// <summary>
    /// The empty / degenerate-tree boundary: a null root, a single leaf, and a two-leaf tree each
    /// carry NO non-trivial clade (GetLeaves(null) yields nothing; CollectClades null-guards; a 1-/
    /// 2-leaf tree produces an empty clade set). RF therefore stays a finite, non-negative count and
    /// NEVER NullReferences, KeyNotFounds, or divides by zero (Tree_Comparison.md §6.1 "no non-trivial
    /// clade → 0"). The empty-vs-empty and empty-vs-populated pairs are both pinned.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void EmptyOrSingleLeafOrNullTree_NoCrash_DistanceIsCladeCount()
    {
        var populated = PhylogeneticAnalyzer.ParseNewick("((A,B),(C,D));");
        var singleLeaf = new PhylogeneticAnalyzer.PhyloNode { Name = "A", Taxa = new List<string> { "A" } };
        var twoLeaf = PhylogeneticAnalyzer.ParseNewick("(A,B);");

        // null vs null: both clade sets empty → 0, no NullReference.
        Action nullPair = () => PhylogeneticAnalyzer.RobinsonFouldsDistance(null!, null!);
        nullPair.Should().NotThrow("two null trees have empty clade sets → guarded, never NullReference");
        PhylogeneticAnalyzer.RobinsonFouldsDistance(null!, null!)
            .Should().Be(0, "no clade on either side → symmetric difference 0");

        // null vs populated: only the populated side contributes clades; no crash, symmetric.
        int nullVsPop = PhylogeneticAnalyzer.RobinsonFouldsDistance(null!, populated);
        nullVsPop.Should().BeGreaterThanOrEqualTo(0, "a one-sided clade set yields a non-negative count");
        PhylogeneticAnalyzer.RobinsonFouldsDistance(populated, null!)
            .Should().Be(nullVsPop, "RF stays symmetric with one empty operand (INV-RF-02)");

        // single leaf: no non-trivial clade (a lone taxon is trivial) → self-RF 0, no crash.
        PhylogeneticAnalyzer.RobinsonFouldsDistance(singleLeaf, singleLeaf)
            .Should().Be(0, "a single-leaf tree has no non-trivial clade → RF(T,T)=0 (INV-RF-01)");
        Action leafVsPop = () => PhylogeneticAnalyzer.RobinsonFouldsDistance(singleLeaf, populated);
        leafVsPop.Should().NotThrow("a single-leaf vs a populated tree must not crash");

        // two-leaf tree: still no non-trivial clade (the only internal node is the full-tree root).
        PhylogeneticAnalyzer.RobinsonFouldsDistance(twoLeaf, twoLeaf)
            .Should().Be(0, "a two-leaf tree carries no non-trivial clade → RF 0");
        PhylogeneticAnalyzer.RobinsonFouldsDistance(twoLeaf, populated)
            .Should().BeGreaterThanOrEqualTo(0, "two-leaf vs richer tree is a non-negative clade count, no crash");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PHYLO-BOOT-001 — phylogenetic bootstrap (Felsenstein 1985) : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PHYLO-BOOT-001 — phylogenetic bootstrap

    #region Helpers — bootstrap

    /// <summary>
    /// Independent re-derivation of a clade key as the implementation and the doc define it
    /// (Bootstrap_Analysis.md §4.2): the subtree's leaf names, sorted ascending and '|'-joined.
    /// Computed here from the TAXON NAMES alone — NOT read back off the algorithm's own output —
    /// so an expected key is derived from the spec, not echoed from the code.
    /// </summary>
    private static string CladeKey(params string[] leafNames) =>
        string.Join("|", leafNames.OrderBy(n => n, StringComparer.Ordinal));

    /// <summary>
    /// Pins the universal documented support contract on a returned map (Bootstrap_Analysis.md
    /// §2.4): every value is finite (no NaN/±Inf), lies in [0,1] (INV-01), and is quantized to
    /// k/replicates for an integer k (INV-02). Returns nothing; asserts.
    /// </summary>
    private static void SupportMapObeysContract(
        IReadOnlyDictionary<string, double> support, int replicates)
    {
        foreach (var kvp in support)
        {
            double s = kvp.Value;
            double.IsNaN(s).Should().BeFalse("a bootstrap support value must never be NaN (clade {0})", kvp.Key);
            double.IsInfinity(s).Should().BeFalse("a bootstrap support value must never be ±Infinity (clade {0})", kvp.Key);
            s.Should().BeGreaterThanOrEqualTo(0.0, "support = count/B ≥ 0 (INV-01, clade {0})", kvp.Key);
            s.Should().BeLessThanOrEqualTo(1.0, "support = count/B ≤ 1 (INV-01, clade {0})", kvp.Key);

            // INV-02: support = k/B for an integer k ⇒ s*B is an integer.
            double scaled = s * replicates;
            scaled.Should().BeApproximately(Math.Round(scaled), 1e-9,
                "support is quantized to (integer count)/replicates (INV-02, clade {0})", kvp.Key);
        }
    }

    #endregion

    #region Positive sanity — documented worked example (two well-separated groups → 1.0)

    /// <summary>
    /// Positive-sanity anchor reproducing the doc's hand-checkable worked example verbatim
    /// (Bootstrap_Analysis.md §7.1): the four-taxon alignment A=B="AAAAAAAAAA",
    /// C=D="GGGGGGGGGG". Because d(A,B)=d(C,D)=0 and the A/B-vs-C/D split is fully saturated,
    /// EVERY column-resample leaves the same {A,B},{C,D} structure, so every replicate recovers
    /// both clades ⇒ support["A|B"] = support["C|D"] = 1.0 (count = B ⇒ B/B = 1; INV-05, ref 1).
    /// The expected keys and values are derived from the SPEC, not read off the algorithm's map.
    /// </summary>
    [Test]
    [CancelAfter(15_000)]
    public void Bootstrap_TwoWellSeparatedGroups_DocumentedExample_BothCladesFullySupported()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "AAAAAAAAAA",
            ["B"] = "AAAAAAAAAA",
            ["C"] = "GGGGGGGGGG",
            ["D"] = "GGGGGGGGGG",
        };

        var support = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 100, seed: 42);

        SupportMapObeysContract(support, replicates: 100);

        // The two non-trivial clades of the reference tree, keyed independently from the taxon names.
        string ab = CladeKey("A", "B");
        string cd = CladeKey("C", "D");

        support.Should().ContainKey(ab, "{{A,B}} is a non-trivial clade of the reference tree (INV-03)");
        support.Should().ContainKey(cd, "{{C,D}} is a non-trivial clade of the reference tree (INV-03)");
        support[ab].Should().Be(1.0, "every resampled replicate recovers {{A,B}} → support 1.0 (Bootstrap_Analysis.md §7.1, INV-05)");
        support[cd].Should().Be(1.0, "every resampled replicate recovers {{C,D}} → support 1.0 (Bootstrap_Analysis.md §7.1, INV-05)");

        // INV-03: the reported keys are EXACTLY the reference tree's non-trivial clades — for a
        // 4-taxon rooted tree split into two cherries, that is precisely {A,B} and {C,D}.
        support.Keys.Should().BeEquivalentTo(new[] { ab, cd },
            "only the reference tree's non-trivial clades are scored (INV-03)");
    }

    /// <summary>
    /// INV-04 (reproducibility): for a FIXED (alignment, replicates, methods, seed) the support
    /// map is byte-for-byte identical across independent calls — resampling is the only randomness
    /// and it is seeded (Bootstrap_Analysis.md §2.4 INV-04). Two different seeds are allowed to
    /// differ but must each independently obey the [0,1]/quantization contract.
    /// </summary>
    [Test]
    [CancelAfter(15_000)]
    public void Bootstrap_FixedSeed_IsReproducible_AcrossRuns()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGTACGTAC",
            ["B"] = "ACGTACGTACGTAG", // close to A
            ["C"] = "TGCATGCATGCATG",
            ["D"] = "TGCATGCATGCATC", // close to C
        };

        var run1 = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 50, seed: 12345);
        var run2 = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 50, seed: 12345);

        run2.Should().BeEquivalentTo(run1, "a fixed seed makes the resampling reproducible (INV-04)");
        SupportMapObeysContract(run1, replicates: 50);
    }

    #endregion

    #region BE — single tree (minimal taxon set: 2 taxa → empty, 3 taxa → one clade)

    /// <summary>
    /// BE "single tree" — the minimal taxon set. A 2-taxon alignment builds a single tree whose
    /// only internal node is the full-tree root: that is the TRIVIAL clade {both taxa} (all taxa),
    /// which the non-trivial filter excludes (CollectClades: subtreeTaxa.Count &lt; totalLeaves) —
    /// so the support map is EMPTY (no clade is scorable), never a crash and never a spurious
    /// entry (INV-03). The smallest input yielding ONE scored clade is 3 taxa (the closest pair
    /// forms the single non-trivial clade). Both pinned; the value still obeys [0,1]/quantization.
    /// </summary>
    [Test]
    [CancelAfter(10_000)]
    public void Bootstrap_SingleMinimalTree_TwoTaxaEmpty_ThreeTaxaOneClade()
    {
        // 2 taxa: only the root (trivial, = all taxa) exists ⇒ no non-trivial clade ⇒ empty map.
        var twoTaxa = new Dictionary<string, string>
        {
            ["X"] = "ACGTACGTAC",
            ["Y"] = "ACGTACGTTT",
        };
        IReadOnlyDictionary<string, double> twoSupport = null!;
        Action twoAct = () => twoSupport = PhylogeneticAnalyzer.Bootstrap(twoTaxa, replicates: 20, seed: 221_001);
        twoAct.Should().NotThrow("a valid 2-taxon alignment builds a tree; it simply has no non-trivial clade");
        twoSupport.Should().BeEmpty("a 2-taxon tree's only internal node is the trivial full-set root → no scored clade (INV-03)");

        // 3 taxa with a clear closest pair (P,Q) ⇒ exactly ONE non-trivial clade {P,Q}.
        var threeTaxa = new Dictionary<string, string>
        {
            ["P"] = "AAAAAAAAAA",
            ["Q"] = "AAAAAAAAAG", // closest to P
            ["R"] = "GGGGGGGGGG", // far outgroup
        };
        var threeSupport = PhylogeneticAnalyzer.Bootstrap(threeTaxa, replicates: 30, seed: 221_002);
        SupportMapObeysContract(threeSupport, replicates: 30);
        threeSupport.Should().HaveCount(1, "a 3-taxon rooted tree has exactly one non-trivial clade (INV-03)");
        string pq = CladeKey("P", "Q");
        threeSupport.Should().ContainKey(pq, "the closest pair {{P,Q}} is the single non-trivial clade");
        threeSupport[pq].Should().Be(1.0,
            "{{P,Q}} are identical-but-one and saturated against R, so every replicate recovers it → 1.0 (INV-05)");
    }

    #endregion

    #region BE — 0 (and -1) replicates (denominator B → documented throw, never div-by-zero)

    /// <summary>
    /// BE "0 replicates" — the KEY divide-by-zero hazard. Support = count/B, so B = 0 would make
    /// every value 0/0 = NaN (or a DivideByZeroException on an integer path). The source GUARDS
    /// this: `replicates &lt; 1` throws ArgumentException BEFORE any resampling or division
    /// (PhylogeneticAnalyzer.cs line 1184; Bootstrap_Analysis.md §3.3, §6.1). Both the 0 and the
    /// BE −1 boundary are pinned as documented rejections — never a crash, NaN, or Infinity. The
    /// smallest VALID denominator B = 1 is pinned to return support quantized to {0,1} = count/1.
    /// </summary>
    [Test]
    [CancelAfter(10_000)]
    public void Bootstrap_ZeroOrNegativeReplicates_ThrowDocumentedException_NeverDivideByZero()
    {
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "AAAAAAAAAA",
            ["B"] = "AAAAAAAAAG",
            ["C"] = "GGGGGGGGGG",
            ["D"] = "GGGGGGGGGC",
        };

        Action zero = () => PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 0, seed: 221_010);
        zero.Should().Throw<ArgumentException>("0 replicates is an invalid denominator → documented rejection, not 0/0 NaN (§3.3, §6.1)");

        Action negative = () => PhylogeneticAnalyzer.Bootstrap(sequences, replicates: -1, seed: 221_010);
        negative.Should().Throw<ArgumentException>("the BE −1 replicate boundary is below the ≥1 minimum → documented rejection (§3.3)");

        // Smallest valid denominator: B = 1 ⇒ every support is count/1 ∈ {0,1}, finite, no crash.
        var one = PhylogeneticAnalyzer.Bootstrap(sequences, replicates: 1, seed: 221_011);
        SupportMapObeysContract(one, replicates: 1);
        one.Values.Should().OnlyContain(s => s == 0.0 || s == 1.0,
            "with B = 1 each support is count/1, an integer 0 or 1 (INV-02)");
    }

    /// <summary>
    /// The remaining documented validation throws on the bootstrap surface (Bootstrap_Analysis.md
    /// §3.3, §6.1): null sequences → ArgumentNullException; fewer than 2 sequences (a tree needs
    /// ≥2 taxa) → ArgumentException. Intentional, contract-defined rejections — not raw crashes.
    /// </summary>
    [Test]
    public void Bootstrap_NullOrTooFewSequences_ThrowDocumentedValidationExceptions()
    {
        Action nullSeqs = () => PhylogeneticAnalyzer.Bootstrap(null!, replicates: 10);
        nullSeqs.Should().Throw<ArgumentNullException>("a null sequence dictionary is an explicit, documented rejection (§3.3)");

        var single = new Dictionary<string, string> { ["Only"] = "ACGTACGT" };
        Action tooFew = () => PhylogeneticAnalyzer.Bootstrap(single, replicates: 10);
        tooFew.Should().Throw<ArgumentException>("fewer than 2 sequences cannot form a tree (§3.3)");

        var empty = new Dictionary<string, string>();
        Action none = () => PhylogeneticAnalyzer.Bootstrap(empty, replicates: 10);
        none.Should().Throw<ArgumentException>("zero sequences is below the 2-taxon minimum (§3.3)");
    }

    #endregion

    #region BE — identical sequences (no phylogenetic signal → defined, no NaN/div-by-zero)

    /// <summary>
    /// BE "identical sequences" — no phylogenetic signal. All-identical sequences give a
    /// zero-distance matrix; resampling identical columns yields identical pseudo-alignments, so
    /// every replicate rebuilds the same degenerate topology. This must NOT NaN or divide by zero,
    /// and the returned support is DEFINED and ∈ [0,1] for every reported clade
    /// (Bootstrap_Analysis.md §6.1, INV-05). Whatever non-trivial clades the deterministic
    /// reference tie-break happens to produce are recovered identically each replicate ⇒ their
    /// support is 1.0; the contract (finite, [0,1], quantized) is pinned regardless of the tie-break.
    /// </summary>
    [Test]
    [CancelAfter(10_000)]
    public void Bootstrap_AllIdenticalSequences_NoSignal_DefinedSupportNoNaN()
    {
        const int n = 6;
        string identical = RandomDna(20, seed: 221_020);
        var sequences = new Dictionary<string, string>();
        for (int i = 0; i < n; i++)
            sequences[$"T{i}"] = identical;

        IReadOnlyDictionary<string, double> support = null!;
        Action act = () => support = PhylogeneticAnalyzer.Bootstrap(
            sequences, replicates: 40, seed: 221_021);
        act.Should().NotThrow("a zero-distance (no-signal) matrix must not divide by zero or crash the bootstrap");

        SupportMapObeysContract(support, replicates: 40);

        // INV-05: every reference clade is recovered in every replicate (each replicate is the same
        // degenerate tree), so every reported support is exactly 1.0 — no signal, but fully defined.
        support.Values.Should().OnlyContain(s => s == 1.0,
            "identical sequences reproduce one topology each replicate → every reported clade has support 1.0 (§6.1, INV-05)");
    }

    #endregion

    #region BE — randomized boundary sweep (no crash / [0,1] / quantization across many shapes)

    /// <summary>
    /// Randomized boundary sweep: over many locally-seeded random alignments of varying taxon
    /// counts, sequence lengths, distance/tree methods and small replicate counts, the bootstrap
    /// must NEVER crash/hang/NaN and must ALWAYS satisfy the documented support contract
    /// (Bootstrap_Analysis.md §2.4): finite values in [0,1] (INV-01), quantized to count/B
    /// (INV-02), and keyed exactly by the reference tree's non-trivial clades (INV-03). A locally
    /// fixed master seed makes the whole sweep reproducible (no shared static RNG); [CancelAfter]
    /// fails the test on any hang.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void Bootstrap_RandomizedSweep_NeverCrashes_AlwaysObeysSupportContract()
    {
        var master = new Random(221_777);
        var methods = AllMethods;
        var treeMethods = new[]
        {
            PhylogeneticAnalyzer.TreeMethod.UPGMA,
            PhylogeneticAnalyzer.TreeMethod.NeighborJoining,
        };

        for (int iter = 0; iter < 40; iter++)
        {
            int taxaCount = 2 + master.Next(6);   // 2..7 taxa (includes the 2-taxon empty-clade boundary)
            int length = 1 + master.Next(30);     // 1..30 columns (includes the single-column boundary)
            int replicates = 1 + master.Next(15); // 1..15 replicates (small B near the boundary)
            var distanceMethod = methods[master.Next(methods.Length)];
            var treeMethod = treeMethods[master.Next(treeMethods.Length)];
            int seed = master.Next();

            var sequences = new Dictionary<string, string>();
            for (int t = 0; t < taxaCount; t++)
                sequences[$"S{t}"] = RandomDna(length, seed: master.Next());

            IReadOnlyDictionary<string, double> support = null!;
            Action act = () => support = PhylogeneticAnalyzer.Bootstrap(
                sequences, replicates, distanceMethod, treeMethod, seed);
            act.Should().NotThrow(
                "bootstrap must not crash on a random {0}-taxon, {1}-col, B={2}, {3}/{4} input",
                taxaCount, length, replicates, distanceMethod, treeMethod);

            SupportMapObeysContract(support, replicates);

            // INV-03: every reported key is a non-trivial clade — its leaf set is a strict, &gt;1-taxon
            // subset of the full taxon set. Re-derived from the key string itself, not from the code.
            var allTaxa = new HashSet<string>(sequences.Keys);
            foreach (var key in support.Keys)
            {
                var members = key.Split('|');
                members.Length.Should().BeGreaterThan(1, "a non-trivial clade has more than one taxon (INV-03)");
                members.Length.Should().BeLessThan(allTaxa.Count, "a non-trivial clade is a strict subset of all taxa (INV-03)");
                members.All(allTaxa.Contains).Should().BeTrue("every clade member is one of the input taxa (INV-03)");
            }
        }
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PHYLO-STATS-001 — tree statistics : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region PHYLO-STATS-001 — tree statistics (leaves, total length, height)

    #region Helpers (PHYLO-STATS-001)

    /// <summary>A leaf node carrying the given name and branch length (no children).</summary>
    private static PhylogeneticAnalyzer.PhyloNode Leaf(string name, double branchLength = 0.0)
        => new(name) { BranchLength = branchLength };

    /// <summary>An internal node with the given branch length and children.</summary>
    private static PhylogeneticAnalyzer.PhyloNode Internal(double branchLength, params PhylogeneticAnalyzer.PhyloNode[] children)
        => new() { BranchLength = branchLength, Children = children.ToList() };

    /// <summary>
    /// Builds a STAR tree: one internal root (branch length <paramref name="rootBranch"/>) whose
    /// children are <paramref name="leafCount"/> leaves L0..L(n-1), each with the same
    /// <paramref name="leafBranch"/>. The shallowest non-trivial shape: height 1, all leaves direct.
    /// </summary>
    private static PhylogeneticAnalyzer.PhyloNode BuildStar(int leafCount, double rootBranch, double leafBranch)
    {
        var children = new PhylogeneticAnalyzer.PhyloNode[leafCount];
        for (int i = 0; i < leafCount; i++)
            children[i] = Leaf($"L{i}", leafBranch);
        return Internal(rootBranch, children);
    }

    /// <summary>
    /// Builds a binary CATERPILLAR (deep ladder) on <paramref name="leafCount"/> leaves: each
    /// internal node has one leaf child and one internal child, nested (leafCount-1) levels deep.
    /// Every edge (internal and leaf) carries <paramref name="edgeBranch"/>. This is the maximally
    /// imbalanced rooted binary tree — height = leafCount - 1.
    /// </summary>
    private static PhylogeneticAnalyzer.PhyloNode BuildCaterpillar(int leafCount, double edgeBranch)
    {
        if (leafCount < 2) throw new ArgumentOutOfRangeException(nameof(leafCount));

        // Deepest internal node joins the last two leaves.
        var node = Internal(edgeBranch, Leaf($"L{leafCount - 2}", edgeBranch), Leaf($"L{leafCount - 1}", edgeBranch));
        // Each step up adds one leaf and one internal level.
        for (int i = leafCount - 3; i >= 0; i--)
            node = Internal(edgeBranch, Leaf($"L{i}", edgeBranch), node);
        return node;
    }

    #endregion

    #region Positive sanity — documented worked example (§7.1)

    /// <summary>
    /// Positive-sanity anchor (Tree_Statistics.md §7.1, hand-checkable): the balanced four-taxon
    /// tree ((A:1,B:1):1,(C:1,D:1):1) has leaf set {A,B,C,D} (count 4), total length = the SIX unit
    /// edges = 6.0, and topological height 2 (root→internal→leaf = 2 edges). Expected values are
    /// derived independently from the doc, not from the code. Pinned on a HAND-BUILT tree and
    /// cross-checked against ParseNewick of the SAME documented string.
    /// </summary>
    [Test]
    public void TreeStatistics_DocumentedBalancedExample_LeavesLengthHeight_MatchSpec()
    {
        // ((A:1,B:1):1,(C:1,D:1):1)  — six unit-length edges.
        var ab = Internal(1.0, Leaf("A", 1.0), Leaf("B", 1.0));
        var cd = Internal(1.0, Leaf("C", 1.0), Leaf("D", 1.0));
        var root = Internal(0.0, ab, cd); // root edge has no length in the doc string

        var leafNames = PhylogeneticAnalyzer.GetLeaves(root).Select(l => l.Name).ToList();
        leafNames.Should().BeEquivalentTo(new[] { "A", "B", "C", "D" }, "the four taxa are the leaves (§7.1)");
        PhylogeneticAnalyzer.GetLeaves(root).Count().Should().Be(4, "an N-leaf tree returns exactly N leaves (INV-02)");
        PhylogeneticAnalyzer.GetLeaves(root).All(l => l.IsLeaf).Should().BeTrue("every returned node IsLeaf (INV-01)");

        PhylogeneticAnalyzer.CalculateTreeLength(root)
            .Should().BeApproximately(6.0, 1e-10, "Σ of the six unit edges = 6.0 (§7.1; INV-03)");

        PhylogeneticAnalyzer.GetTreeDepth(root)
            .Should().Be(2, "root→internal→leaf is 2 edges (§7.1; INV-05)");

        // Cross-check the SAME tree parsed from the documented Newick string.
        var parsed = PhylogeneticAnalyzer.ParseNewick("((A:1,B:1):1,(C:1,D:1):1);");
        PhylogeneticAnalyzer.GetLeaves(parsed).Select(l => l.Name).Should()
            .BeEquivalentTo(new[] { "A", "B", "C", "D" }, "ParseNewick recovers the same four leaves");
        PhylogeneticAnalyzer.CalculateTreeLength(parsed).Should().BeApproximately(6.0, 1e-10, "parsed total length = 6.0");
        PhylogeneticAnalyzer.GetTreeDepth(parsed).Should().Be(2, "parsed height = 2");
    }

    #endregion

    #region BE — Empty / null tree (the −1 boundary)

    /// <summary>
    /// The empty-tree convention (Tree_Statistics.md §2.4 INV-06, §6.1): a null root yields NO
    /// leaves, total length 0, and height −1 — NEVER a NullReferenceException. This is the BE −1
    /// probe, distinct from a single leaf's height 0.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void NullTree_YieldsEmptyLeaves_ZeroLength_AndMinusOneHeight()
    {
        Action act = () =>
        {
            _ = PhylogeneticAnalyzer.GetLeaves(null!).ToList();
            _ = PhylogeneticAnalyzer.CalculateTreeLength(null!);
            _ = PhylogeneticAnalyzer.GetTreeDepth(null!);
        };
        act.Should().NotThrow("the null/empty tree is a defined boundary, never a NullReference");

        PhylogeneticAnalyzer.GetLeaves(null!).Should().BeEmpty("GetLeaves(null) yields nothing (INV-06)");
        PhylogeneticAnalyzer.CalculateTreeLength(null!).Should().Be(0.0, "an empty tree has no edges → length 0 (INV-06)");
        PhylogeneticAnalyzer.GetTreeDepth(null!).Should().Be(-1, "an empty tree has height −1 by convention (INV-05/06)");
    }

    #endregion

    #region BE — Single leaf (degenerate smallest tree: height 0)

    /// <summary>
    /// A single leaf (Tree_Statistics.md §6.1): exactly one leaf (itself), height 0 ("a tree with
    /// only a single node has height zero"), and total length = its OWN BranchLength (sum over the
    /// one node) — specifically NOT −1 (that is the null/empty case). Pinned for both a default-0
    /// branch and a non-zero branch.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void SingleLeaf_OneLeaf_HeightZero_LengthEqualsOwnBranch()
    {
        // Default branch (0): length 0, height 0, one leaf.
        var bare = Leaf("Only");
        var bareLeaves = PhylogeneticAnalyzer.GetLeaves(bare).ToList();
        bareLeaves.Should().ContainSingle("a single leaf is one terminal node (INV-02)");
        bareLeaves[0].Name.Should().Be("Only");
        bareLeaves[0].IsLeaf.Should().BeTrue("the lone node IsLeaf (INV-01)");
        PhylogeneticAnalyzer.CalculateTreeLength(bare).Should().Be(0.0, "a default-branch leaf has length 0 (§6.1)");
        PhylogeneticAnalyzer.GetTreeDepth(bare).Should().Be(0, "a single node has height 0, NOT −1 (§6.1; INV-05)");

        // Non-zero branch: length = its own branch length; still one leaf, height 0.
        var withBranch = Leaf("Solo", 2.5);
        PhylogeneticAnalyzer.GetLeaves(withBranch).Should().ContainSingle();
        PhylogeneticAnalyzer.CalculateTreeLength(withBranch)
            .Should().BeApproximately(2.5, 1e-10, "total length of a single leaf = its own BranchLength (§6.1; INV-03)");
        PhylogeneticAnalyzer.GetTreeDepth(withBranch).Should().Be(0, "a single node has height 0 regardless of its branch length");
    }

    #endregion

    #region BE — Star tree (one internal root, all leaves direct: height 1)

    /// <summary>
    /// A STAR tree — one internal root whose children are ALL leaves directly (the shallowest
    /// non-trivial shape). For N leaves: GetLeaves returns exactly N; height = 1 (root not a leaf →
    /// 1 + max(child height 0)); total length = root.BranchLength + Σ child branch lengths. Pinned
    /// for N=2 (binary minimum) and a WIDE N (polytomy — exercises the N-ary Children traversal of
    /// every method, not a first-two-children shortcut).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void StarTree_AllLeavesDirectChildren_HeightOne_LengthIsSumOfEdges()
    {
        // N=2 binary star: root branch 0.5, two leaf branches 1.0 each → length 2.5, height 1.
        var star2 = BuildStar(leafCount: 2, rootBranch: 0.5, leafBranch: 1.0);
        PhylogeneticAnalyzer.GetLeaves(star2).Count().Should().Be(2, "a 2-leaf star has 2 leaves (INV-02)");
        PhylogeneticAnalyzer.GetTreeDepth(star2).Should().Be(1, "root→leaf is exactly 1 edge in a star (INV-05)");
        PhylogeneticAnalyzer.CalculateTreeLength(star2)
            .Should().BeApproximately(0.5 + 2 * 1.0, 1e-10, "length = rootBranch + Σ leaf branches (INV-03)");

        // Wide polytomy star: 7 leaves. root branch 0 → length = 7 × 0.3; height still 1.
        const int n = 7;
        var star7 = BuildStar(leafCount: n, rootBranch: 0.0, leafBranch: 0.3);
        var leaves7 = PhylogeneticAnalyzer.GetLeaves(star7).ToList();
        leaves7.Should().HaveCount(n, "every direct child of the star is returned — all {0} leaves (INV-02)", n);
        leaves7.All(l => l.IsLeaf).Should().BeTrue("each returned node IsLeaf (INV-01)");
        leaves7.Select(l => l.Name).Should().OnlyHaveUniqueItems("the leaf names are distinct L0..L6");
        PhylogeneticAnalyzer.GetTreeDepth(star7).Should().Be(1, "a star of any width has height 1 (INV-05)");
        PhylogeneticAnalyzer.CalculateTreeLength(star7)
            .Should().BeApproximately(n * 0.3, 1e-10, "Σ of the {0} leaf edges (root edge is 0); traverses ALL children (INV-03)", n);
    }

    #endregion

    #region BE — Deep ladder / caterpillar (maximally imbalanced: height N−1)

    /// <summary>
    /// A binary CATERPILLAR (deep ladder) on N leaves — the maximally imbalanced rooted binary tree
    /// (Tree_Statistics.md §2.1 "how many levels deep"; INV-05). Independently hand-derived: N leaves,
    /// height = N−1 (one edge per internal level on the longest root→leaf path), total length = Σ over
    /// all N+(N−1) unit edges. Pinned for several N, including a DEEP one to exercise the recursive
    /// traversals to depth (guarded by [CancelAfter] against any hang).
    /// </summary>
    [Test]
    [CancelAfter(10_000)]
    public void Caterpillar_DeepLadder_HeightIsNMinusOne_AndLeafCountIsN()
    {
        foreach (int n in new[] { 2, 3, 5, 10, 64 })
        {
            var tree = BuildCaterpillar(n, edgeBranch: 1.0);

            PhylogeneticAnalyzer.GetLeaves(tree).Count().Should().Be(n, "a caterpillar on {0} leaves returns {0} leaves (INV-02)", n);
            PhylogeneticAnalyzer.GetLeaves(tree).All(l => l.IsLeaf).Should().BeTrue("every returned node IsLeaf (INV-01)");

            PhylogeneticAnalyzer.GetTreeDepth(tree).Should().Be(
                n - 1, "a caterpillar on {0} leaves is maximally deep: height = N−1 (INV-05)", n);

            // Edges: N leaf edges + (N−1) internal edges, all unit length → total = 2N−1.
            PhylogeneticAnalyzer.CalculateTreeLength(tree)
                .Should().BeApproximately(2.0 * n - 1.0, 1e-10, "Σ of the N + (N−1) unit edges (INV-03)");
        }
    }

    #endregion

    #region BE — Randomized boundary sweep (no crash, exact contract over random shapes)

    /// <summary>
    /// A randomized boundary sweep over locally-seeded random STARS and CATERPILLARS asserts the
    /// no-crash / leaf-count / height / length contract across many sizes and branch lengths. Expected
    /// values are computed independently from the topology (star: height 1, leaf count = N, length =
    /// root + Σ leaves; caterpillar: height N−1, leaf count = N, length = Σ all edges) — NOT read off
    /// the code. No NaN, no Infinity, no negative length for non-negative branches (INV-04).
    /// </summary>
    [Test]
    [CancelAfter(15_000)]
    public void RandomizedSweep_StarsAndCaterpillars_ObeyExactStatisticsContract()
    {
        var master = new Random(222_001);

        for (int iter = 0; iter < 60; iter++)
        {
            int n = 2 + master.Next(20);                 // 2..21 leaves
            double rootBranch = master.NextDouble() * 3;
            double leafBranch = master.NextDouble() * 3;
            bool star = master.Next(2) == 0;

            PhylogeneticAnalyzer.PhyloNode tree;
            int expectedHeight;
            double expectedLength;

            if (star)
            {
                tree = BuildStar(n, rootBranch, leafBranch);
                expectedHeight = 1;
                expectedLength = rootBranch + n * leafBranch;
            }
            else
            {
                double edge = leafBranch; // uniform edge length for the caterpillar
                tree = BuildCaterpillar(n, edge);
                expectedHeight = n - 1;
                expectedLength = (2.0 * n - 1.0) * edge; // N leaf + (N−1) internal edges
            }

            int leafCount = 0;
            double length = 0;
            int height = 0;
            Action act = () =>
            {
                leafCount = PhylogeneticAnalyzer.GetLeaves(tree).Count();
                length = PhylogeneticAnalyzer.CalculateTreeLength(tree);
                height = PhylogeneticAnalyzer.GetTreeDepth(tree);
            };
            act.Should().NotThrow("statistics must not crash on a random {0} of {1} leaves", star ? "star" : "caterpillar", n);

            leafCount.Should().Be(n, "leaf count equals the leaf count of the constructed topology (INV-02)");
            height.Should().Be(expectedHeight, "height matches the independently-derived topology height (INV-05)");
            double.IsNaN(length).Should().BeFalse("total length must never be NaN");
            double.IsInfinity(length).Should().BeFalse("total length must never be ±Infinity");
            length.Should().BeApproximately(expectedLength, 1e-9, "length = Σ edges, independently derived (INV-03)");
            length.Should().BeGreaterThanOrEqualTo(0.0, "Σ of non-negative branches is non-negative (INV-04)");
        }
    }

    #endregion

    #endregion
}
