**STATUS: DOABLE-AWAITING-USER (created 2026-11-24 by B09 WP74; result to be returned as `returned/A2-19c-dash-out.zip`). `mafft --dash` calls the Go program `dash_client` (`core/dash_client.go`, 700 lines; `core/Makefile`: `DASH_CLIENT = dash_client`, "Go language is required"), which the Debian binary package `mafft` 7.505-1 does not contain, and the DASH server `https://sysimm.org/dash/REST1.0/` (HTTP 200 from the audit machine, 2026-11-23) answers a live database: the client's output cannot be produced or checked on the audit machine (never compile or run downloaded code there). The package builds the real `dash_client` and the real MAFFT 7.505 source tree in a container, records EVERY HTTP exchange of the client with the server (a logging proxy: `exchanges.jsonl`, replayable offline) and runs `mafft --dash` with its options.**

# TEMPORARY hand-off package A2-19c-dash (for the user's Mac)

**Why.** B09 item A2-19c-s5 / `--dash`: `mafft --dash` (the structural information of DASH: `newdash=1`, `distance=global`, `iterate=3` when none is given; script lines 2413 - 2480) writes the input without gaps to `dashin`, calls
`dash_client -url "$dashserver" -i dashin -sequences dashsequences -hat3 hat3.seed`, prepends the template sequences (`>DASH_...` renamed `>DASH|...`) to the input, uses `hat3.seed` (the blocks of the structural alignments of the templates, kind `k`) as the pair library of the
iterative refinement (the `--seed` machinery: `tbfast -J`, `dvtditr -l`) and removes the `DASH|` sequences from the output (`f2cl -d`, `restoreu`). To port `dash_client` (template selection `chain_search_sequence`, the `domains` / `domain_alignments` queries, slicing, the filter 22.5 %, `equivalence`,
the `hat3` writer) and to verify the ported script, the raw behaviour is needed, and because the server's answers change over time, EVERY exchange is logged so the next session replays the same answers offline.

**What it does.** (1) `dash_client` itself on each of the 37 sets of `inputs/sets.tsv` (14 well-known proteins with PDB homologs - ubiquitin, lysozymes, myoglobin, cytochrome c, BPTI, GFP, calmodulin, thioredoxin, barnase, RNase A, haemoglobin chains, staphylococcal nuclease -, mixtures, mutated families, a nucleotide set, unknown letters,
random sequences without homologs; `make_inputs.py`, deterministic): `dashin`, standard output / error, return code, `dashsequences`, `hat3.seed`, the exchanges. (2) `mafft <options> in.fa` with the Debian-equivalent serial build (`ENABLE_MULTITHREAD=`: "nthread = 0" as the Debian binary prints) for the 18 option lines of `inputs/options.tsv`
(`--dash`, with `--quiet`, `--retree 1`, `--localpair`, `--globalpair`, `--maxiterate 0 / 10`, `--anysymbol`, `--originalseqonly`, `--excludehomologs`, `--reorder --clustalout`, `--treeout --distout`, `--op --ep`, `--dashserver URL` with and without the final slash, `--merge table`, `--amino`, `--nuc`): the first line on every set, the others on 5 sets;
return code, rows, standard error, `in.fa.tree` / `in.fa.hat2` / `GuideTree`, the exchanges. The client talks to a logging proxy on 127.0.0.1:8765 (`dashproxy.py`, standard library, forwards to the real server with the neutral User-Agent `Seqeron-review/1.0`); the proxy's log has the request and response bodies (base64).
Source: the upstream tarball `https://mafft.cbrc.jp/alignment/software/mafft-7.505-without-extensions-src.tgz` (its SHA-256 goes to `out/versions.txt`; the audit machine's copy has `9b9ee0dea1d493f0fad6f9c227bcc291ddedd045b5558b7fea30dc5932db33b8`), Debian 12 `golang-go`. No account data, no secrets, no e-mail / name anywhere.

## Run it (Apple Silicon: the image is x86-64 and is emulated; 30 - 90 minutes, mostly the server)

```bash
git pull
cd docs/Validation/review-2026-09/handoff/A2-19c-dash
./run.sh
ls -l returned/          # A2-19c-dash-out.zip : commit it under docs/Validation/review-2026-09/handoff/A2-19c-dash/returned/
```
The machine needs internet access to `sysimm.org` and `mafft.cbrc.jp`. If the server is unreachable the client's message is the result (the package records it too).

## How the next session integrates it

1. `scripts/alignment/convert_dash_returned.py` (to be written then) embeds the runs (`dash_oracle.tsv.gz`: input, options, return code, rows, standard error, files) and the exchanges (`dash_exchanges.jsonl.gz`).
2. Port `dash_client.go` (the HTTP part is replaced by a replay of the logged exchanges: the same request -> the logged response; the client's own logic - the sorting of templates, the slicing, the filter, `hat3` writer, `dashsequences` writer, the messages "Querying DASH for domains - [n%]" ... - is verified against `dashsequences` / `hat3.seed` / standard output of the direct runs), then the script's `--dash` branch
   on top of the `--seed` machinery (B09.md F342 ff.): `MafftProgram.Run` with a `dashserver` replay (a files entry that names the exchanges) must equal the `mafft --dash` runs in return code, rows, standard error and files.
3. The refusal of `--dash` in `MafftScriptRun` (today: the script's "Install DASH client" message, which an installation without `dash_client` prints) stays for an installation without the client.
