#!/usr/bin/env bash
# Builds / installs candidate LAST versions and selects one: /opt/last -> the first candidate (LAST 1256, 1250, 1282, 1411, 1542 built from source, then the Debian package) whose real lastal call of mafft (-m 10 -e 5000 -f 0 -s 1 -p matrix -a 3599 -b 179)
# on smoke_db.fa / smoke_q.fa answers within SMOKE_FAST seconds (default 90) with at least one alignment; none fast: the fastest one. (Run 1 selected by "mafft exits 0 within 60 s", which no candidate did - see README.)
# Never fails: every step is logged (the image build must not stop; the run on the Mac records which candidate was selected).
set +e
PREFIX_MAFFT=/usr/lib/mafft/lib/mafft
mkdir -p /opt/last-cand
: > /opt/last-candidates.txt
UA="Seqeron-review/1.0"

smoke() {   # $1 = directory with lastdb / lastal: the real call of mafft 7.505 (core/pairlocalalign.c lastcallthread) for the last query of smoke.fa, timed. Selection is by the TIME of this call (and a non-empty answer), NOT by
            # `mafft --lastpair`: run 1 of this package selected by "mafft exits 0 within 60 s" and found no candidate - see README (the program mafft itself crashes for >= 3 related sequences, and LAST >= 1447 needs ~26 s per call).
  d=$(mktemp -d); cp /opt/smoke_db.fa $d/_db3d.fa; cp /opt/smoke_q.fa $d/_q3; cp /opt/smoke.mat $d/_scoringmatrixforlast
  ( cd $d && $1/lastdb _db3d _db3d.fa > lastdb.out 2>&1; s=$(date +%s); timeout ${SMOKE_MAX:-900} $1/lastal -m 10 -e 5000 -f 0 -s 1 -p _scoringmatrixforlast -a 3599 -b 179 _db3d _q3 > lastres 2> lastal.err; rc=$?; e=$(date +%s); echo "rc=$rc seconds=$((e - s))" > timing )
  hits=$(grep -vc '^#' $d/lastres)
  echo "candidate $1: $(cat $d/timing), hit lines $hits, $(head -c 100 $d/lastal.err | tr '\n' ' ')" | tee -a /opt/last-candidates.txt
  secs=$(sed 's/.*seconds=//' $d/timing | cut -d. -f1)
  echo "$secs $1" >> /opt/last-times.txt
  [ "$hits" -ge 1 ] && [ "${secs:-9999}" -le "${SMOKE_FAST:-90}" ]
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
# 3. select: the contemporary sources first (mafft 7.505 is of April 2022; LAST 1256 of January 2022), then the other sources, the Debian package last; a candidate is taken if the real lastal call of mafft answers within SMOKE_FAST seconds
: > /opt/last-times.txt
for c in $(for v in $LAST_VERSIONS; do echo /opt/last-cand/$v; done) /opt/last-cand/debian; do
  [ -x $c/bin/lastal ] || continue
  if smoke $c/bin; then rm -rf /opt/last; ln -s $c /opt/last; echo "selected $c" | tee -a /opt/last-candidates.txt; break; fi
done
if [ ! -e /opt/last ]; then   # none is fast: take the fastest that answered
  best=$(sort -n /opt/last-times.txt | head -1 | cut -d' ' -f2); [ -n "$best" ] && ln -s "$(dirname $best)" /opt/last && echo "no fast candidate; selected the fastest: $best" | tee -a /opt/last-candidates.txt
fi
ls -l /opt/last 2>&1 | tee -a /opt/last-candidates.txt
exit 0
