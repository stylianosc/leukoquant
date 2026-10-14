#!/bin/bash
# Builds dcm2niix with -DmyInstanceNumberOrderIsNotSpatial, the variant that
# dcm2niix itself suggests when it stops with "Recompiling with
# '-DmyInstanceNumberOrderIsNotSpatial' might help". It orders slices by
# position instead of by DICOM InstanceNumber, which rescues series whose
# InstanceNumbers do not follow slice order (seen on ADNI-3 DWI).
#
# dwi_utils._convert_dicom_with_fallback() looks for the result as
# leukoquant/external/noddi/dcm2niix_nonspatial and tries it after the stock
# dcm2niix (leukoquant/external/noddi/dcm2niix) fails.
#
# Usage: build_dcm2niix_nonspatial.sh <output_binary_path>
# Run it as an SGE job, not on the login node, inside the job container
# (needs only git, g++ and make). After a successful build,
# commit the binary to leukoquant/external/noddi/dcm2niix_nonspatial like the
# other external binaries (see build_niftyseg.sh for the convention).
set -euo pipefail

OUTPUT_BINARY="${1:?Usage: build_dcm2niix_nonspatial.sh <output_binary_path>}"
# Pinned so the build is reproducible.
DCM2NIIX_TAG="v1.0.20260724"

BUILD_DIR="$(mktemp -d)"
trap 'rm -rf "$BUILD_DIR"' EXIT

git clone --quiet --depth 1 --branch "$DCM2NIIX_TAG" https://github.com/rordenlab/dcm2niix.git "$BUILD_DIR/src"

# dcm2niix's own Makefile needs only g++ and make (the job container has no
# cmake). CFLAGS is overridden wholesale, so its default flags are repeated.
make -C "$BUILD_DIR/src/console" CFLAGS="-O3 -s -flto -DmyInstanceNumberOrderIsNotSpatial"

cp "$BUILD_DIR/src/console/dcm2niix" "$OUTPUT_BINARY"
chmod +x "$OUTPUT_BINARY"
"$OUTPUT_BINARY" -v | head -1
