| input | kills (signal landed) | boundary | random | chained | finished before the signal | target after the kill: none / complete / damaged | hash equal | nothing else left |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dna_100m (n = 100000000) | 28 | 20 | 8 | 1 | 0 | 26 / 2 / 0 | 28 of 28 | 28 of 28 |
| dna_10m (n = 10000000) | 73 | 20 | 53 | 5 | 3 | 62 / 14 / 0 | 76 of 76 | 76 of 76 |
| prot_10m (n = 10000000) | 74 | 20 | 54 | 7 | 2 | 72 / 4 / 0 | 76 of 76 | 76 of 76 |
| rep_100m (n = 100000000) | 8 | 0 | 8 | 1 | 0 | 6 / 2 / 0 | 8 of 8 | 8 of 8 |
| rep_10m (n = 10000000) | 76 | 20 | 56 | 5 | 0 | 72 / 4 / 0 | 76 of 76 | 76 of 76 |

Total: 259 kills landed in 264 cases; hash equal in 264, nothing else left in 264. RssAnon of the finishing runs: max 252.1 MB over 264 runs (the gate is 320).

Kills by the phase that the process was in (the last marker that it printed before the signal) and what the resume skipped (the most common set of stages). The seconds of the resumes are in the raw rows but are not a measurement of the saving: the resumes of these campaigns run with 1, 2 or 4 threads at random and right after the kill; the saving is measured by run_1b9e_timing.py.

| input | ph0 | ph1 | ph23 | ph4 | bwt | psi | ph7 | seal |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dna_100m | 4 ((none)) | 4 (packed) | 2 (packed,ranks) | 6 (packed,sa) | 2 (packed,sa,lcp) | 1 (packed,sa,lcp,bwt) | 6 (packed,sa,lcp,psi) | 3 ((none)) |
| dna_10m | 10 ((none)) | 10 (packed) | 8 (packed,ranks) | 11 (packed,sa) | 4 (packed,sa,lcp) | 10 (packed,sa,lcp,bwt) | 8 (packed,sa,lcp,psi) | 12 ((none)) |
| prot_10m | 10 ((none)) | 8 (packed) | 10 (packed,ranks) | 21 (packed,sa) | 12 (packed,sa,lcp) | 4 (packed,sa,lcp,bwt) | 6 (packed,sa,lcp,psi) | 3 (packed,sa,lcp,psi) |
| rep_100m | 1 (packed,ranks) | 1 (packed) | 1 (packed,ranks) | 1 (packed,sa) | 1 (packed,sa,lcp) | 0 | 1 (packed,sa,lcp,psi) | 2 ((none)) |
| rep_10m | 10 ((none)) | 9 (packed) | 10 (packed,ranks) | 19 (packed,sa) | 13 (packed,sa,lcp) | 3 (packed,sa,lcp,bwt) | 7 (packed,sa,lcp,psi) | 5 ((none)) |
