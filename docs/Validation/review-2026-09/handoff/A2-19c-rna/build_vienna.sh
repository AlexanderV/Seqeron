#!/usr/bin/env bash
# Obtains ViennaRNA (RNAalifold) for the image. Debian 12 has NO installable `vienna-rna` in `main` (run 1 of this package: "E: Package 'vienna-rna' has no installation candidate"; the package is in `non-free`), so
#   step 1: build the upstream release 2.5.1 (the version of the Debian package) from the GitHub release tarball with CPU-friendly flags (no Perl / Python / SWIG / OpenMP / LTO, -O2, no -march=native);
#   step 2: if that fails, the same for 2.4.18 and 2.6.4;
#   step 3: if that fails, `apt-get install vienna-rna` with `non-free` enabled.
# Every step is logged to stdout (the Dockerfile writes it to /opt/build_vienna.log); nothing here ever fails the image build (exit 0); the result is judged by in_container.sh (RNAalifold --version in out/versions.txt).
set -u
UA="Seqeron-review/1.0"
PREFIX=/opt/vienna
mkdir -p /opt/vsrc && cd /opt/vsrc
have() { [ -x "$PREFIX/bin/RNAalifold" ] || command -v RNAalifold > /dev/null 2>&1; }
for V in 2.5.1 2.4.18 2.6.4; do
  have && break
  echo "== ViennaRNA $V"
  U="https://github.com/ViennaRNA/ViennaRNA/releases/download/v$V/ViennaRNA-$V.tar.gz"
  rm -rf "ViennaRNA-$V" v.tgz
  curl -fsSL -A "$UA" --retry 3 -o v.tgz "$U" || { echo "download failed: $U"; continue; }
  sha256sum v.tgz
  tar xzf v.tgz || { echo "unpack failed"; continue; }
  cd "ViennaRNA-$V" || continue
  for FLAGS in "--without-perl --without-python --without-python2 --without-swig --without-kinfold --without-forester --without-rnalocmin --without-svm --without-gsl --without-json --disable-openmp --disable-lto --disable-doc" \
               "--without-perl --without-python --without-swig --disable-lto" \
               "--without-perl --without-python"; do
    have && break
    echo "-- configure $FLAGS"
    make distclean > /dev/null 2>&1
    CFLAGS="-O2 -g0" CXXFLAGS="-O2 -g0" ./configure --prefix=$PREFIX $FLAGS > configure.log 2>&1 || { echo "configure failed:"; tail -15 configure.log; continue; }
    grep -i "warning: unrecognized" configure.log | head -5
    ( make -j2 > make.log 2>&1 && make install > install.log 2>&1 ) || { echo "make failed:"; tail -25 make.log; continue; }
  done
  cd /opt/vsrc
  have && echo "built ViennaRNA $V"
done
if [ -x "$PREFIX/bin/RNAalifold" ]; then ln -sf "$PREFIX/bin/RNAalifold" /usr/local/bin/RNAalifold; ln -sf "$PREFIX/bin/RNAfold" /usr/local/bin/RNAfold; fi
if ! have; then
  echo "== apt non-free vienna-rna"
  echo "deb http://deb.debian.org/debian bookworm main non-free non-free-firmware" > /etc/apt/sources.list.d/nonfree.list
  apt-get update 2>&1 | tail -3
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends vienna-rna 2>&1 | tail -10
fi
echo "== result"; command -v RNAalifold && RNAalifold --version 2>&1 | head -1
exit 0
