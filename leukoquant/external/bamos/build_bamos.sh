#!/bin/bash
# Builds source/ (a private CMIC fork of NiftySeg 0.9.4 -- see
# leukoquant/external/bamos/source/README.txt) into the two BaMoS-specific
# applications, Seg_BiASM and Seg_Analysis. These do not exist in mainline
# NiftySeg (confirmed via diff against niftyseg/source_latest/seg-apps/), so
# this build must come from this fork specifically -- nothing else in the
# repo can produce them. Everything else this fork's seg-apps/ could build
# (seg_maths, seg_EM, seg_LabFusion, seg_LoAd, seg_stats) is deliberately
# NOT built here (see source/seg-apps/CMakeLists.txt) since BaMoS_WMH_*.sh
# sources those from the separately built leukoquant/external/niftyseg/bin/
# instead; see niftyseg/build_niftyseg.sh.
#
# Usage: build_bamos.sh <install_prefix>
#   install_prefix: a scratch install location (this script's `make install`
#   populates install_prefix/bin/{Seg_BiASM,Seg_Analysis}). After a
#   successful build, copy those two binaries into
#   leukoquant/external/bamos/bin/, replacing what's there.
set -euo pipefail

INSTALL_PREFIX="${1:?Usage: build_bamos.sh <install_prefix>}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$SCRIPT_DIR/source"

export PATH=/share/apps/gcc-9.2.0/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/gcc-9.2.0/lib64:${LD_LIBRARY_PATH:-}

BUILD_DIR="$(mktemp -d)"
mkdir -p "$BUILD_DIR/build"
cd "$BUILD_DIR/build"
cmake "$SOURCE_DIR" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=/share/apps/gcc-9.2.0/bin/gcc \
  -DCMAKE_CXX_COMPILER=/share/apps/gcc-9.2.0/bin/g++ \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX"
make -j4
make install

rm -rf "$BUILD_DIR"
echo BUILD_AND_INSTALL_DONE
