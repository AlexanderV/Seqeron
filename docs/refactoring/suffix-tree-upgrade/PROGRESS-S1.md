# S1 progress (ESA backend behind ISuffixTree, managed, no native code)

Rules: [PROCESS.md](PROCESS.md). Gates: [PLAN.md](PLAN.md) section S1 plus **peak RSS during build <= 24 B/char** (about 240 MB at 10 Mbp; Ukkonen peaks at 1.76 GB).
Pre-registered S1a gate (10 Mbp random DNA, 1 thread, this machine class): SA+LCP build <= 3.3 s (<= 10x libsais 0.328 s), index <= 16 B/char incl. text.
Targets: ceiling goal <= 1.0 s (3x), aspiration <= 0.5 s (1.5x); libsais reference 0.328 s, libdivsufsort 0.662 s, Go 0.73 s.

Decisions already taken by the owner: managed only; positions are `long` in the new contract from the start; contract-first (stable contract, replaceable engines);
Ukkonen stays untouched (oracle, online); ESA persistent container is a separate format kind (name decision pending); native accelerator not allowed.

| step | content | status |
|---|---|---|
| S1-0 | Contract + ADRs (docs only), SA construction approach and optimization roadmap, S1a/S1b work-package list with oracles | [ ] |
| S1a | SA-IS (int32, then generic), LCP (Phi/PLCP), oracles (brute force n<=2000, Ukkonen, libsais vectors, adversarial strings), differential fuzz, early speed/memory measurement vs the gate | [ ] |
| S1a-opt | Iterative optimization per the stopping rule in PROCESS.md | [ ] |
| go/no-go | Owner-visible decision after S1a numbers | [ ] |
| S1b | Child table, LCP-interval navigator reusing `SuffixTreeAlgorithms.*<TNode,TNav>`, public navigation/handle/lazy-locate API, factory, mmap container | [ ] |

## Log
(append newest last)
