#!/bin/bash
# Builds the CUDA-enabled NiftyReg used by GIF's bin_gpu/seg_GIF and, separately, the
# downloaded niftyreg/gpu tarball for metrics/z-score/BaMoS (see container_utils.py's
# ensure_niftyreg_gpu()). Source lives in the sibling directory source/ (kept
# outside it deliberately -- that directory carries upstream's own inherited
# .gitignore, whose 'build*' pattern would otherwise swallow this script), vendored
# from https://github.com/KCL-BMEIS/niftyreg, master branch (currently at commit
# cd5440fd). This vendored copy now uses upstream's RNifti-based NIfTI I/O
# throughout (a switch from the older classic nifticlib API), which required
# reintroducing classic-style include guards in nifti1_io.h/nifti1.h/znzlib.h so
# GIF can still link against both NiftyReg and NiftySeg in the same binary
# without struct redefinition collisions -- NiftySeg's own build now forwards
# its copies of those 3 files to these ones for the same reason (see
# leukoquant/external/niftyseg/source_latest/nifti/). All patches applied to
# this vendored copy are documented inline at each patch site (search for
# "[leukoquant patch]"); this script only performs the build/install, no
# further source modification.
#
# Run on a GPU node so the checkCudaCard.cpp sanity check (upstream CMake) succeeds.
# GPU jobs on this cluster request tmem only, not h_vmem -- h_vmem on a GPU job has
# caused CUDA OOM elsewhere in this project:
#   qsub -l gpu=true -pe gpu 1 -l tmem=8G -l h_rt=1:0:0 build_niftyreg_cuda.sh
#
# Usage: build_niftyreg_cuda.sh <install_prefix>
set -euo pipefail

INSTALL_PREFIX="${1:?Usage: build_niftyreg_cuda.sh <install_prefix>}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="$SCRIPT_DIR/source"

export PATH=/share/apps/gcc-9.2.0/bin:/share/apps/cuda-11.8/bin:/share/apps/python-3.13.0a6-shared/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/gcc-9.2.0/lib64:/share/apps/cuda-11.8/lib64:/share/apps/python-3.13.0a6-shared/lib:${LD_LIBRARY_PATH:-}

BUILD_DIR="$(mktemp -d)"
mkdir -p "$BUILD_DIR/build"
cd "$BUILD_DIR/build"
# [leukoquant patch] reg-lib/cuda/CMakeLists.txt's CHECK_GPU option (ON by
# default) runs a small program at configure time that detects the *build
# node's own* GPU and overwrites CMAKE_CUDA_ARCHITECTURES with just that one
# card's capability, silently discarding any value passed here on the
# command line -- confirmed directly: passing our own multi-arch list still
# produced a binary with only sm_75 embedded (whatever GPU the build
# happened to land on), which then failed with "no kernel image is
# available for execution on the device" on every other architecture in
# the fleet (confirmed on a real Ada Lovelace/RTX 4070 Ti Super node).
# CHECK_GPU=OFF routes to that same CMakeLists.txt's own *other* branch,
# which sets a sensible broad default itself
# ("60-real;61-real;70-real;75-real;80-real;86-real;89") -- letting
# upstream's own intended fallback do this is more robust than us
# hand-maintaining a separate list that can silently stop applying if this
# override logic changes again.
cmake "$SOURCE_DIR" -DCMAKE_BUILD_TYPE=Release -DUSE_CUDA=ON \
  -DCHECK_GPU=OFF \
  -DCMAKE_C_COMPILER=/share/apps/gcc-9.2.0/bin/gcc \
  -DCMAKE_CXX_COMPILER=/share/apps/gcc-9.2.0/bin/g++ \
  -DCMAKE_CUDA_HOST_COMPILER=/share/apps/gcc-9.2.0/bin/g++ \
  -DCMAKE_CUDA_FLAGS="--fmad=false" \
  -DCMAKE_C_FLAGS="-ffp-contract=off" \
  -DCMAKE_CXX_FLAGS="-ffp-contract=off" \
  -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX"
make -j4
make install

# reg-lib/CMakeLists.txt's install() rules only cover a subset of the headers
# actually needed by downstream consumers of Platform.h (GIF's _seg_GIF.cpp
# transitively pulls in ContentCreatorFactory.h, MeasureFactory.h, etc., none of
# which have their own install() entry -- confirmed via trial builds, not something
# worth hand-maintaining a matching install(FILES) list for). Bulk-copy the full
# top-level reg-lib header set so any consumer of NIFTYREG_DIR/include has
# everything it transitively needs.
cp "$SOURCE_DIR"/reg-lib/*.h "$INSTALL_PREFIX/include/"

rm -rf "$BUILD_DIR"
echo BUILD_AND_INSTALL_DONE
