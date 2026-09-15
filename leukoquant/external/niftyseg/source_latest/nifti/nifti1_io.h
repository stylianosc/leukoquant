// [leukoquant patch] Forwards to the vendored NiftyReg's RNifti-aware copy
// instead of NiftySeg's own classic one, so NiftySeg's compiled library and
// GIF's code (which links against NiftyReg headers) agree on the same
// underlying struct identity (NiftyReg's version aliases nifti_image to a
// distinct nifti1_image struct tag; NiftySeg's own copy defines the struct
// directly as nifti_image -- two different type identities at the ABI
// level despite being structurally identical, which breaks linking
// anything that mixes NiftySeg's precompiled library with NiftyReg-based
// code). Needs the NiftyReg install's include/ directory on the compiler's
// include search path (added by build_niftyseg.sh's NIFTYREG_INSTALL_DIR
// argument) for this to resolve.
#include "niftilib/nifti1_io.h"
