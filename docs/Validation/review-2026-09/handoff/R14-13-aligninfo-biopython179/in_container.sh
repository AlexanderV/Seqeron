#!/usr/bin/env bash
# Runs inside the image: makes the cases, runs them on Biopython 1.79, writes /work/out/{cases.jsonl,results.jsonl,environment.txt}.
set -u
cd /work
rm -rf out; mkdir -p out
python -I make_cases.py out/cases.jsonl
python -I run_cases.py out/cases.jsonl out/results.jsonl
{ echo "uname -m: $(uname -m)"; python --version; python -c "import Bio; print('Biopython', Bio.__version__)"; sha256sum out/cases.jsonl out/results.jsonl; date -u; } > out/environment.txt
cat out/environment.txt
if [ -n "${HOST_UID:-}" ]; then chown -R "$HOST_UID:$HOST_GID" out 2>/dev/null || true; fi
