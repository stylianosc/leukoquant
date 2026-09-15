#/*============================================================================
#
#  NifTK: A software platform for medical image computing.
#
#  Copyright (c) University College London (UCL). All rights reserved.
#
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.
#
#  See LICENSE.txt in the top level directory for details.
#
#============================================================================*/

# Construct consistent error messages for use below.
set(NIFTYREG_DIR_DESCRIPTION "directory containing the file 'include/_reg_tools.h'. This is either the root of the build tree, or PREFIX for an installation.")
set(NIFTYREG_DIR_MESSAGE "NIFTYREG not found.  Set the NIFTYREG_DIR cmake cache entry to the ${NIFTYREG_DIR_DESCRIPTION}")


if(NOT NIFTYREG_FOUND)

  # Look for reg_tools.h in build trees or under <prefix>/include.
  find_path(NIFTYREG_DIR
    NAMES include/_reg_f3d.h
    HINTS ENV NIFTYREG_DIR
    PATHS

    # Help the user find it if we cannot.
    DOC "The ${NIFTYREG_DIR_DESCRIPTION}"
    )

  if(NIFTYREG_DIR)
    if(EXISTS ${NIFTYREG_DIR}/include/_reg_f3d.h)
      set(NIFTYREG_FOUND 1)
    else()
      set(NIFTYREG_DIR "NIFTYREG_DIR-NOTFOUND" CACHE PATH "The ${NIFTYREG_DIR_DESCRIPTION}" FORCE)
    endif()
  endif()

endif()


if (NIFTYREG_FOUND)

  # Library list updated for the current (2026, post-#130 "reg_resample cuda
  # enabling" refactor) NiftyReg source tree, which has a different internal
  # module layout than what this Find module was originally written against.
  # Confirmed by listing the actual .a files an install of this tree
  # produces: _reg_femTrans and _reg_common_cuda no longer exist as separate
  # libraries (their functionality was absorbed into other consolidated
  # targets during that refactor); _reg_content, _reg_platform, and znz are
  # new libraries this tree didn't have before. Wrapped in --start-group/
  # --end-group so GNU ld resolves inter-archive symbol references
  # regardless of order, since this Find module (unlike NiftyReg's own CMake
  # targets) links these as plain -l flags with no CMake-tracked transitive
  # dependency graph.
  set(NIFTYREG_LIBRARIES
    -Wl,--start-group
    _reg_aladin
    _reg_f3d
    _reg_measure
    _reg_blockMatching
    _reg_globalTrans
    _reg_localTrans
    _reg_resampling
    _reg_ReadWriteImage
    _reg_tools
    _reg_maths
    _reg_compute
    _reg_content
    _reg_platform
    _reg_kernels
    _reg_cuda_kernels
    _reg_cudainfo
    reg_png
    reg_nifti
    znz
    z
    -Wl,--end-group
    # [leukoquant patch] The CUDA-enabled NiftyReg .a libs above (built with
    # -DUSE_CUDA=ON) contain real calls into the CUDA runtime (cudaLaunchKernel,
    # cudaPeekAtLastError, cudaGetErrorString, etc.) that nothing in this list
    # otherwise resolves, since these are plain -l links (not CMake-target-aware,
    # so they don't pull in NiftyReg's own PUBLIC_LINK_LIBRARIES transitively).
    # Without this, linking seg_GIF against a CUDA NiftyReg build fails with
    # "undefined reference to `cudaLaunchKernel'" and similar.
    cudart
    # Separately, CudaContextSingleton (lib_reg_cuda_kernels.a) calls the CUDA
    # DRIVER API directly (cuInit, cuCtxCreate_v2, etc. -- the `cu`-prefixed
    # symbols, distinct from the `cuda`-prefixed runtime API above), which
    # only libcuda.so (not libcudart.so) provides.
    cuda
    # CudaLts.cu (least-squares/SVD point-set solver, new in this refactor)
    # calls cusolverDn* functions directly -- only libcusolver.so provides
    # these, separate from both the CUDA runtime and driver libs above.
    cusolver
    )

  set(NIFTYREG_INCLUDE_DIR
    ${NIFTYREG_DIR}/include
    )
  set(NIFTYREG_LIBRARY_DIR
    ${NIFTYREG_DIR}/lib
    )

else()
  # Eigen not found, explain to the user how to specify its location.
  if(NIFTYREG_FIND_REQUIRED)
    message(FATAL_ERROR ${NIFTYREG_DIR_MESSAGE})
  else()
    if(NOT NIFTYREG_FIND_QUIETLY)
      message(STATUS ${NIFTYREG_DIR_MESSAGE})
    endif()
  endif()

endif()
