| experiment | label | n pairs | seconds A (median [min-max]) | seconds B | B / A per pair (median [min-max]) | device MB written A | device MB written B | fsync calls B | recording s (B) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Resume off (A) vs on (B) | dna_100m.fsync | 5 | 17.171 [16.544-17.891] | 16.715 [15.820-18.101] | 1.006 [0.921-1.015] | 2330 | 2330 | 14 | 0.320 [0.306-0.336] |
| Resume off (A) vs on (B) | dna_100m.nofsync | 5 | 17.041 [16.313-17.704] | 17.533 [16.905-17.921] | 1.046 [0.958-1.099] | 2329 | 2330 | 0 | 0.323 [0.305-0.338] |
| Resume off (A) vs on (B) | dna_10m.fsync | 8 | 2.073 [1.791-2.171] | 2.125 [2.010-2.178] | 1.033 [0.972-1.122] | 194 | 194 | 14 | 0.043 [0.042-0.046] |
| Resume off (A) vs on (B) | dna_10m.nofsync | 8 | 1.822 [1.715-2.148] | 1.911 [1.777-2.003] | 1.011 [0.902-1.146] | 194 | 194 | 0 | 0.038 [0.036-0.040] |
| Resume off (A) vs on (B) | prot_10m.fsync | 8 | 2.032 [1.875-2.120] | 2.163 [1.920-2.193] | 1.057 [0.914-1.154] | 198 | 198 | 14 | 0.041 [0.036-0.042] |
| Resume off (A) vs on (B) | prot_10m.nofsync | 8 | 2.151 [2.041-2.188] | 2.024 [1.882-2.212] | 0.952 [0.892-1.019] | 198 | 198 | 0 | 0.039 [0.034-0.042] |
| Resume off (A) vs on (B) | rep_10m.fsync | 8 | 2.439 [2.210-2.599] | 2.529 [2.229-2.600] | 1.008 [0.973-1.139] | 209 | 209 | 14 | 0.043 [0.040-0.045] |
| Resume off (A) vs on (B) | rep_10m.nofsync | 8 | 2.469 [2.240-2.584] | 2.435 [2.261-2.560] | 1.007 [0.903-1.034] | 209 | 209 | 0 | 0.041 [0.037-0.044] |
| the build before the package (A) vs this commit, Resume off (B) | dna_100m.head.fsync | 5 | 17.426 [16.187-17.804] | 16.639 [15.136-17.130] | 0.962 [0.879-0.989] | | | | |
| the build before the package (A) vs this commit, Resume off (B) | dna_100m.head.nofsync | 5 | 15.017 [14.864-15.683] | 15.277 [15.068-15.793] | 1.017 [0.961-1.056] | | | | |
| the build before the package (A) vs this commit, Resume off (B) | dna_10m.head.fsync | 8 | 1.811 [1.717-2.004] | 1.826 [1.736-2.051] | 1.020 [0.966-1.088] | | | | |
| the build before the package (A) vs this commit, Resume off (B) | dna_10m.head.nofsync | 8 | 1.927 [1.844-2.083] | 1.896 [1.804-2.129] | 0.985 [0.866-1.110] | | | | |
| the build before the package (A) vs this commit, Resume off (B) | prot_10m.head.fsync | 8 | 1.949 [1.871-2.249] | 1.878 [1.843-2.235] | 0.981 [0.858-1.130] | | | | |
| the build before the package (A) vs this commit, Resume off (B) | prot_10m.head.nofsync | 8 | 2.143 [1.914-2.267] | 2.123 [1.954-2.190] | 0.999 [0.862-1.043] | | | | |
| the build before the package (A) vs this commit, Resume off (B) | rep_10m.head.fsync | 8 | 2.262 [2.135-2.501] | 2.214 [2.133-2.514] | 1.001 [0.908-1.126] | | | | |
| the build before the package (A) vs this commit, Resume off (B) | rep_10m.head.nofsync | 8 | 2.450 [2.329-2.532] | 2.332 [2.220-2.517] | 0.965 [0.897-1.020] | | | | |
