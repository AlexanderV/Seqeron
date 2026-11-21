| input | kills (signal landed) | boundary | random | chained | finished before the signal | target after the kill: none / complete / damaged | hash equal | nothing else left |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dna_100m (n = 100000000) | 6 | 0 | 6 | 4 | 2 | 6 / 2 / 0 | 8 of 8 | 8 of 8 |
| dna_10m (n = 10000000) | 14 | 0 | 14 | 2 | 2 | 11 / 5 / 0 | 16 of 16 | 16 of 16 |
| prot_10m (n = 10000000) | 10 | 0 | 10 | 2 | 6 | 9 / 7 / 0 | 16 of 16 | 16 of 16 |
| rep_100m (n = 100000000) | 8 | 0 | 8 | 4 | 0 | 8 / 0 / 0 | 8 of 8 | 8 of 8 |
| rep_10m (n = 10000000) | 15 | 0 | 15 | 2 | 1 | 13 / 3 / 0 | 16 of 16 | 16 of 16 |

Total: 53 kills landed in 64 cases; hash equal in 64, nothing else left in 64. RssAnon of the finishing runs: max 255.8 MB over 64 runs (the gate is 320).

Kills by the phase that the process was in (the last marker that it printed before the signal) and what the resume skipped (the most common set of stages). The seconds of the resumes are in the raw rows but are not a measurement of the saving: the resumes of these campaigns run with 1, 2 or 4 threads at random and right after the kill; the saving is measured by run_1b9e_timing.py.

| input | ph0 | ph1 | ph23 | ph4 | bwt | psi | ph7 | seal |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| dna_100m | 1 (packed,ranks) | 1 (packed,sa) | 1 (packed,ranks) | 1 (packed,sa,lcp) | 1 ((none)) | 1 (packed,sa,lcp,bwt) | 0 | 0 |
| dna_10m | 2 ((none)) | 2 (packed) | 2 (packed,ranks) | 1 (packed,sa) | 1 (packed,sa,lcp) | 3 (packed,sa,lcp,bwt) | 0 | 3 ((none)) |
| prot_10m | 2 ((none)) | 2 (packed) | 2 (packed,ranks) | 3 (packed,sa) | 0 | 0 | 0 | 1 ((none)) |
| rep_100m | 1 (packed,ranks) | 1 (packed,sa) | 1 (packed,ranks) | 3 ((none)) | 1 (packed,sa,lcp) | 1 (packed,sa,lcp,bwt) | 0 | 0 |
| rep_10m | 2 ((none)) | 2 (packed) | 2 (packed,ranks) | 3 (packed,sa) | 1 (packed,sa,lcp) | 1 (packed,sa,lcp,bwt) | 2 (packed,sa,lcp,psi) | 2 (packed,ranks) |
