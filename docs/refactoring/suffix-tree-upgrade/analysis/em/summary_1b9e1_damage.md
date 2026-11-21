| marker | damage | skipped (before / after) | bytes that differ from the uninterrupted build, before the fix | after the fix |
| --- | --- | --- | --- | --- |
| lcp-child-built | child_last | packed,sa / packed,sa | 9 | 0 |
| lcp-child-built | tail_a5 | packed,sa,lcp / packed,sa,lcp | 0 | 0 |
| lcp-child-built | fill_a5 | packed,sa,lcp / packed,sa,lcp | 11047 | 0 |
| packed | tail_a5 | packed / packed | 10971 | 0 |
| sa-built | tail_a5 | packed,sa / packed,sa | 8507 | 0 |
| psi-built | tail_a5 | packed,sa,lcp,psi / packed,sa,lcp,psi | 0 | 0 |
| psi-built | fill_a5 | packed,sa,lcp,psi / packed,sa,lcp,psi | 11047 | 0 |
| info-written | fill_a5 | packed,sa,lcp,psi / packed,sa,lcp,psi | 11047 | 0 |

| file | resumes (each in a new process) | markers | kinds of damage | exit 0 | SHA-256 equal | nothing but the target left |
| --- | --- | --- | --- | --- | --- | --- |
| matrix_1m.tsv | 495 | 18 | 11 | 495 | 495 | 495 |
| matrix_10m.tsv | 198 | 18 | 11 | 198 | 198 | 198 |
| matrix_100m.tsv | 24 | 4 | 6 | 24 | 24 | 24 |
