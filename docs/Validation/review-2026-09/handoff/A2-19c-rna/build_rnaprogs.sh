#!/usr/bin/env bash
# Run 3 of the package A2-19c-rna: the two programs that the MAFFT source archive does NOT contain and that run 2 therefore stopped at ("Install DAFS" for --dafs / --dafspair, "Install CONTRAfold" for --contrafold):
#   CONTRAfold 2.02  (Do et al. 2006; BSD licence)   -> /opt/mafft-ext/libexec/mafft/contrafold   (mafft: `contrafoldwrap -d $prefix` runs `env PATH=$prefix contrafold predict _contrafoldin --posteriors 0.01 _contrafoldout`)
#   DAFS             (Sato et al. 2012)              -> /opt/mafft-ext/libexec/mafft/dafs          (mafft: `mccaskillwrap -G` / `pairlocalalign -G` run `$prefix/dafs --mafft-out ...` / `--mafft-in ...`)
# Every step is logged to stdout (the Dockerfile writes /opt/build_rnaprogs.log); nothing here ever fails the image build (exit 0): in_container.sh lists the programs that exist (out/versions.txt) and the sets that need a missing program
# simply end with the script's own "Install ..." message again. Third-party SOURCE is built here, inside the image, never on the audit machine. Downloads use the neutral User-Agent "Seqeron-review/1.0".
set -u
UA="Seqeron-review/1.0"
LIBEXEC=/opt/mafft-ext/libexec/mafft
VIENNA=/opt/vienna
mkdir -p /opt/rsrc "$LIBEXEC" && cd /opt/rsrc

echo "######## CONTRAfold 2.02"
CF_SHA=9e2657262797f9fcba5c325deaca565401698d682026439476774abee08163cd
if curl -fsSL -A "$UA" --retry 3 -o cf.tgz https://contra.stanford.edu/contrafold/contrafold_v2_02.tar.gz; then
  sha256sum cf.tgz; echo "expected  $CF_SHA (the value seen when the package was written; a different value is only reported)"
  rm -rf contrafold && tar xzf cf.tgz && cd contrafold/src || true
  # the 2008 sources need a few headers that newer g++ no longer includes implicitly; three flag sets, the first that links wins
  for FL in "" "-fpermissive -include cstring -include cstdlib -include cstdio -include climits" "-fpermissive -std=gnu++98 -include cstring -include cstdlib -include cstdio -include climits"; do
    [ -x "$LIBEXEC/contrafold" ] && break
    echo "-- make contrafold CXXFLAGS+= $FL"
    make clean > /dev/null 2>&1
    make contrafold "CXXFLAGS=-O3 -DNDEBUG -pipe -w $FL" > make.log 2>&1 && [ -x contrafold ] && cp contrafold "$LIBEXEC/contrafold" || { echo "make failed:"; tail -15 make.log; }
  done
  cd /opt/rsrc
fi
if [ ! -x "$LIBEXEC/contrafold" ]; then
  # fall back to the bioconda binary of the same version (linux-64, sha256 of the conda repodata when the package was written; a different value is only reported)
  echo "-- fallback: bioconda contrafold 2.02 binary"
  if curl -fsSL -A "$UA" --retry 3 -o cfb.tar.bz2 https://conda.anaconda.org/bioconda/linux-64/contrafold-2.02-h9f5acd7_2.tar.bz2; then
    sha256sum cfb.tar.bz2; echo "expected  40e99eae83b346fe07baeb75237a3181ef3c2ed18b7b110959f6ae53717d469a"
    mkdir -p cfb && tar xjf cfb.tar.bz2 -C cfb && find cfb -name 'contrafold*' -type f | head
    B=$(find cfb -type f -name contrafold -perm -u+x | head -1)
    [ -n "$B" ] && cp "$B" "$LIBEXEC/contrafold"
  fi
fi
[ -x "$LIBEXEC/contrafold" ] && { chmod 755 "$LIBEXEC/contrafold"; ls -l "$LIBEXEC/contrafold"; "$LIBEXEC/contrafold" 2>&1 | head -3; } || echo "CONTRAfold: NOT AVAILABLE"

echo "######## DAFS"
# the repository of the authors (GitHub, satoken/dafs); the first source that downloads and builds wins. DAFS needs the ViennaRNA library (built by build_vienna.sh into $VIENNA).
dafs_build() {
  d="$1"; cd "$d" || return 1
  ls | head -30
  if [ -f CMakeLists.txt ]; then
    echo "-- cmake"; mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$VIENNA" -DCMAKE_CXX_FLAGS="-O2 -fpermissive -include cstdint" > cmake.log 2>&1 && make -j2 > make.log 2>&1 || { tail -20 cmake.log make.log 2>/dev/null; cd ..; }
    [ -d build ] && cd build; B=$(find . -type f -name dafs -perm -u+x | head -1); [ -n "$B" ] && cp "$B" "$LIBEXEC/dafs"; cd "$d"
  fi
  if [ ! -x "$LIBEXEC/dafs" ] && { [ -f configure ] || [ -f configure.ac ]; }; then
    echo "-- autotools"; [ -f configure ] || autoreconf -i > autoreconf.log 2>&1
    CXXFLAGS="-O2 -fpermissive -include cstdint" CPPFLAGS="-I$VIENNA/include" LDFLAGS="-L$VIENNA/lib" ./configure --with-vienna-rna="$VIENNA" > configure.log 2>&1 && make -j2 > make.log 2>&1 || { tail -20 configure.log make.log 2>/dev/null; }
    B=$(find . -type f -name dafs -perm -u+x | head -1); [ -n "$B" ] && cp "$B" "$LIBEXEC/dafs"
  fi
  cd /opt/rsrc
}
for P in \
  "https://github.com/satoken/dafs/archive/refs/tags/v0.0.3.tar.gz dafs-0.0.3" \
  "https://github.com/satoken/dafs/archive/refs/heads/master.tar.gz dafs-master" \
  "http://www.ncrna.org/software/dafs/dafs-0.0.3.tar.gz dafs-0.0.3-ncrna"; do
  set -- $P
  [ -x "$LIBEXEC/dafs" ] && break
  echo "== $1"
  rm -rf "$2" && mkdir -p "$2" && curl -fsSL -A "$UA" --retry 3 -o "$2.tgz" "$1" && sha256sum "$2.tgz" && tar xzf "$2.tgz" -C "$2" --strip-components=1 && dafs_build "/opt/rsrc/$2" || echo "failed: $1"
done
if [ -x "$LIBEXEC/dafs" ]; then chmod 755 "$LIBEXEC/dafs"; ls -l "$LIBEXEC/dafs"; "$LIBEXEC/dafs" --help 2>&1 | head -5; else echo "DAFS: NOT AVAILABLE"; fi

echo "######## DAFS PROBE (A2-19c-rna run 3 analysis, WP86): does this DAFS know the options that mccaskillwrap / pairlocalalign give it?"
# mccaskillwrap.c runs   dafs --mafft-out _mccaskillout _mccaskillin     and pairlocalalign.c runs   dafs --mafft-in _bpp _dafsin ;  the sets of run 3 stop with "ERROR IN mccaskill_mea" (the call fails)
if [ -x "$LIBEXEC/dafs" ]; then
  mkdir -p /tmp/dafsprobe && cd /tmp/dafsprobe && printf '>in\nGGGAAAGCCCUUUGGGAAACCC\n' > one.fa && printf '>a\nGGGAAAGCCCUUUGGGAAACCC\n>b\nGGGAAAGCCCUUUGGGAAAGCC\n' > two.fa
  echo "-- dafs --help (whole)"; "$LIBEXEC/dafs" --help 2>&1
  echo "-- dafs --mafft-out o one.fa"; "$LIBEXEC/dafs" --mafft-out o one.fa > d1 2> d2; echo "rc=$?"; echo "stdout:"; cat d1; echo "stderr:"; cat d2; ls
  echo "-- dafs --mafft-in b two.fa"; "$LIBEXEC/dafs" --mafft-in b two.fa > d1 2> d2; echo "rc=$?"; echo "stdout:"; cat d1; echo "stderr:"; cat d2
  echo "-- dafs one.fa (one sequence, no mafft option)"; "$LIBEXEC/dafs" one.fa > d1 2> d2; echo "rc=$?"; echo "stdout:"; cat d1; echo "stderr:"; cat d2
  echo "-- dafs two.fa"; "$LIBEXEC/dafs" two.fa > d1 2> d2; echo "rc=$?"; echo "stdout:"; cat d1; echo "stderr:"; cat d2
  echo "-- strings dafs | grep -i mafft"; strings "$LIBEXEC/dafs" | grep -i mafft | head
  echo "-- source tree: grep -i mafft"; grep -rIn -i mafft /opt/rsrc/dafs-*/src 2>/dev/null | head
  cd /opt/rsrc
fi
exit 0
