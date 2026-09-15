#!/bin/bash
# Builds the latest vanilla NiftySeg (source_latest/) into a bin/ of standard
# NiftySeg tools (seg_maths, seg_EM, seg_LabFusion, seg_LoAd, seg_stats, etc.).
# These are the tools BaMoS_WMH_*.sh sources for anything that isn't
# BaMoS-specific -- see leukoquant/external/bamos/build_bamos.sh for the two
# BaMoS-specific applications (Seg_BiASM, Seg_Analysis) that this build does
# NOT produce (they don't exist in mainline NiftySeg at all).
#
# Usage: build_niftyseg.sh <install_prefix> <eigen_include_dir> <niftyreg_install_dir>
#   install_prefix: a scratch install location (this script's `make install`
#   populates install_prefix/{bin,lib,include,priors}, matching a normal
#   NiftySeg install tree). Only install_prefix/bin/* is meant to be
#   committed -- after a successful build, copy those binaries into
#   leukoquant/external/niftyseg/bin/, replacing what's there.
#   eigen_include_dir: a directory containing an Eigen/ header subdirectory
#   (a plain Eigen checkout/extraction works; only the Eigen/ headers are
#   needed, not unsupported/ or test/). Passed via -DUSE_SYSTEM_EIGEN=ON to
#   avoid NiftySeg's own ExternalProject_Add(Eigen) download, which on a
#   quota-limited filesystem can exhaust the inode quota (Eigen ships
#   thousands of files under unsupported/ and test/ that are never used here).
#   niftyreg_install_dir: a NiftyReg install tree (bin/lib/include), same as
#   what build_gif.sh needs. [leukoquant patch] source_latest/nifti/{nifti1_io,
#   nifti1,znzlib}.h now forward to this install's copies instead of defining
#   their own, so NiftySeg's compiled library and NiftyReg-linked code (e.g.
#   GIF) agree on one struct identity -- see those files for why. This
#   argument's include/ dir needs to be on the compiler's search path for
#   those forwards to resolve.
set -euo pipefail

INSTALL_PREFIX="${1:?Usage: build_niftyseg.sh <install_prefix> <eigen_include_dir> <niftyreg_install_dir>}"
EIGEN_INCLUDE_DIR="${2:?Usage: build_niftyseg.sh <install_prefix> <eigen_include_dir> <niftyreg_install_dir>}"
NIFTYREG_INSTALL_DIR="${3:?Usage: build_niftyseg.sh <install_prefix> <eigen_include_dir> <niftyreg_install_dir>}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$SCRIPT_DIR/source_latest"

export PATH=/share/apps/gcc-9.2.0/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/gcc-9.2.0/lib64:${LD_LIBRARY_PATH:-}

BUILD_DIR="$(mktemp -d)"
mkdir -p "$BUILD_DIR/build"
cd "$BUILD_DIR/build"
cmake "$SOURCE_DIR" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=/share/apps/gcc-9.2.0/bin/gcc \
  -DCMAKE_CXX_COMPILER=/share/apps/gcc-9.2.0/bin/g++ \
  -DUSE_SYSTEM_EIGEN=ON \
  -DEIGEN_INCLUDE_DIR="$EIGEN_INCLUDE_DIR" \
  -DCMAKE_CXX_STANDARD=17 \
  -DCMAKE_CXX_STANDARD_REQUIRED=ON \
  -DCMAKE_C_FLAGS="-I${NIFTYREG_INSTALL_DIR}/include" \
  -DCMAKE_CXX_FLAGS="-I${NIFTYREG_INSTALL_DIR}/include" \
  -DCMAKE_EXE_LINKER_FLAGS="-L${NIFTYREG_INSTALL_DIR}/lib -lz" \
  -DCMAKE_SHARED_LINKER_FLAGS="-L${NIFTYREG_INSTALL_DIR}/lib -lz" \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX"
make -j4
make install

rm -rf "$BUILD_DIR"
echo BUILD_AND_INSTALL_DONE
