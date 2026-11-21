#!/bin/bash
# Builds gif_source/ into bin_gpu/seg_GIF. Depends on two other builds having
# already been done first:
#   1. leukoquant/external/niftyreg/build_niftyreg_cuda.sh -> a NiftyReg
#      install tree (bin/lib/include). Must be the CUDA-enabled build, even
#      if you only intend to run seg_GIF with -platf 0 (CPU) -- there is only
#      ever one seg_GIF binary, and it always links against this same
#      CUDA-enabled NiftyReg.
#   2. leukoquant/external/niftyseg/build_niftyseg.sh -> a NiftySeg install
#      tree (bin/lib/include).
#
# Run on a GPU node (the NiftyReg dependency's own checkCudaCard.cpp sanity
# check only matters at NiftyReg's own build time, but linking against CUDA
# libraries generally wants the toolkit's runtime libs on LD_LIBRARY_PATH,
# which are only guaranteed present via the cuda module/PATH setup below):
#   qsub -l gpu=true -pe gpu 1 -l tmem=8G,h_vmem=8G -l h_rt=1:0:0 build_gif.sh ...
#
# Usage: build_gif.sh <install_prefix> <niftyreg_install_dir> <niftyseg_install_dir> <eigen_include_dir>
#   install_prefix: a scratch install location (this script's `make install`
#   populates install_prefix/bin/seg_GIF). After a successful build, copy
#   that one binary into leukoquant/external/gif/bin_gpu/seg_GIF, replacing
#   what's there (libcudart.so.11.0 and libcuda.so.1 alongside it in bin_gpu/
#   don't need rebuilding -- they're just copied CUDA runtime/driver-stub
#   libraries, not part of this build).
set -euo pipefail

INSTALL_PREFIX="${1:?Usage: build_gif.sh <install_prefix> <niftyreg_install_dir> <niftyseg_install_dir> <eigen_include_dir>}"
NIFTYREG_INSTALL_DIR="${2:?Usage: build_gif.sh <install_prefix> <niftyreg_install_dir> <niftyseg_install_dir> <eigen_include_dir>}"
NIFTYSEG_INSTALL_DIR="${3:?Usage: build_gif.sh <install_prefix> <niftyreg_install_dir> <niftyseg_install_dir> <eigen_include_dir>}"
EIGEN_INCLUDE_DIR="${4:?Usage: build_gif.sh <install_prefix> <niftyreg_install_dir> <niftyseg_install_dir> <eigen_include_dir>}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$SCRIPT_DIR/gif_source"

export PATH=/share/apps/gcc-9.2.0/bin:/share/apps/cuda-11.8/bin:/share/apps/python-3.13.0a6-shared/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/gcc-9.2.0/lib64:/share/apps/cuda-11.8/lib64:/share/apps/python-3.13.0a6-shared/lib:${LD_LIBRARY_PATH:-}
# LIBRARY_PATH (not LD_LIBRARY_PATH -- that's a runtime-only dynamic loader
# variable) is what gcc/ld consult to resolve -lcudart/-lcuda at link time,
# needed by FindNIFTYREG.cmake's explicit `cudart`/`cuda` links against the
# CUDA-enabled NiftyReg .a libs. libcuda.so itself only exists as a link-time
# stub (no real GPU driver on the build node) under cuda-11.8's own
# lib64/stubs/ -- the exact same stub bundled at runtime as
# leukoquant/external/gif/bin_gpu/libcuda.so.1.
export LIBRARY_PATH=/share/apps/cuda-11.8/lib64:/share/apps/cuda-11.8/lib64/stubs:${LIBRARY_PATH:-}

BUILD_DIR="$(mktemp -d)"
mkdir -p "$BUILD_DIR/build"
cd "$BUILD_DIR/build"
cmake "$SOURCE_DIR" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=/share/apps/gcc-9.2.0/bin/gcc \
  -DCMAKE_CXX_COMPILER=/share/apps/gcc-9.2.0/bin/g++ \
  -DNIFTYREG_DIR="$NIFTYREG_INSTALL_DIR" \
  -DNIFTYSEG_DIR="$NIFTYSEG_INSTALL_DIR" \
  -DUSE_SYSTEM_EIGEN=ON \
  -DEIGEN_INCLUDE_DIR="$EIGEN_INCLUDE_DIR" \
  -DUSE_OPENMP=ON \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX"
make -j4
make install

rm -rf "$BUILD_DIR"
echo BUILD_AND_INSTALL_DONE
