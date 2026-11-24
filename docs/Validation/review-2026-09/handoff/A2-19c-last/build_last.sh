#!/usr/bin/env bash
# Builds / installs candidate LAST versions and selects one: /opt/last/bin -> the first candidate for which `mafft --lastpair /opt/smoke.fa` (with the logging-free real programs) ends within 60 s with return code 0.
# Never fails: every step is logged (the image build must not stop; the run on the Mac records which candidate was selected).
set +e
PREFIX_MAFFT=/usr/lib/mafft/lib/mafft
mkdir -p /opt/last-cand /opt/last
: > /opt/last-candidates.txt
UA="Seqeron-review/1.0"

smoke() {   # $1 = directory with lastdb / lastal
  d=$(mktemp -d)
  for t in lastdb lastal; do printf '#!/bin/sh\nexec %s/%s "$@"\n' "$1" "$t" > $PREFIX_MAFFT/$t; chmod +x $PREFIX_MAFFT/$t; done
  ( cd $d && timeout 60 mafft --lastpair /opt/smoke.fa > out.txt 2> err.txt ); rc=$?
  echo "candidate $1: rc=$rc, output lines $(wc -l < $d/out.txt), $(head -c 120 $d/err.txt | tr '\n' ' ')" | tee -a /opt/last-candidates.txt
  [ $rc -eq 0 ] && [ "$(grep -c '^>' $d/out.txt)" = 4 ]
}

# 1. the Debian package
echo "== Debian package last-align" | tee -a /opt/last-candidates.txt
apt-get update > /dev/null 2>&1; DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends last-align > /opt/last-apt.log 2>&1
if command -v lastal > /dev/null; then
  mkdir -p /opt/last-cand/debian/bin; for t in lastal lastdb; do printf '#!/bin/sh\nexec %s "$@"\n' "$(command -v $t)" > /opt/last-cand/debian/bin/$t; chmod +x /opt/last-cand/debian/bin/$t; done
  lastal --version 2>&1 | head -1 | tee -a /opt/last-candidates.txt
fi
# 2. release sources
for v in $LAST_VERSIONS; do
  echo "== LAST $v" | tee -a /opt/last-candidates.txt
  mkdir -p /opt/last-cand/$v && cd /opt/last-cand/$v
  for url in "https://gitlab.com/mcfrith/last/-/archive/$v/last-$v.tar.gz" "https://gitlab.com/mcfrith/last/-/archive/$v.tar.gz"; do
    curl -fsSL -A "$UA" -o src.tar.gz "$url" && break
  done
  [ -s src.tar.gz ] || { echo "no source for $v" | tee -a /opt/last-candidates.txt; continue; }
  tar xzf src.tar.gz && cd last-* && make -j2 CXXFLAGS="-O3" > ../make.log 2>&1 && make install prefix=/opt/last-cand/$v > ../install.log 2>&1 || echo "build of $v failed (see make.log)" | tee -a /opt/last-candidates.txt
done
# 3. select
for c in /opt/last-cand/debian /opt/last-cand/*; do
  [ -x $c/bin/lastal ] || continue
  if smoke $c/bin; then rm -rf /opt/last; ln -s $c /opt/last; echo "selected $c" | tee -a /opt/last-candidates.txt; break; fi
done
ls -l /opt/last 2>&1 | tee -a /opt/last-candidates.txt
exit 0
