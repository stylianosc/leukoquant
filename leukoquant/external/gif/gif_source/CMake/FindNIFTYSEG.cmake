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
set(NIFTYSEG_DIR_DESCRIPTION "directory containing the file 'include/_seg_common.h'. This is either the root of the build tree, or PREFIX for an installation.")
set(NIFTYSEG_DIR_MESSAGE "NIFTYSEG not found.  Set the NIFTYSEG_DIR cmake cache entry to the ${NIFTYSEG_DIR_DESCRIPTION}")


if(NOT NIFTYSEG_FOUND)
  
  # Look for reg_tools.h in build trees or under <prefix>/include.
  find_path(NIFTYSEG_DIR
    NAMES include/_seg_common.h
    HINTS ENV NIFTYSEG_DIR
    PATHS

    # Help the user find it if we cannot.
    DOC "The ${NIFTYSEG_DIR_DESCRIPTION}"
    )

  if(NIFTYSEG_DIR)
    if(EXISTS ${NIFTYSEG_DIR}/include/_seg_common.h)
      set(NIFTYSEG_FOUND 1)
    else()
      set(NIFTYSEG_DIR "NIFTYSEG_DIR-NOTFOUND" CACHE PATH "The ${NIFTYSEG_DIR_DESCRIPTION}" FORCE)
    endif()
  endif()

endif()


if (NIFTYSEG_FOUND)

  set(NIFTYSEG_LIBRARIES
    _seg_tools
    _seg_EM
    _seg_tools
    # [leukoquant patch] _seg_nifti is now an INTERFACE library in the
    # vendored NiftySeg build (see its nifti/CMakeLists.txt) -- it produces
    # no .a archive of its own to link against with a plain -l flag from an
    # external consumer like this one. Link its real underlying libraries
    # (NiftyReg's own reg_nifti/znz, which it forwards to) directly instead.
    reg_nifti
    znz
    _seg_FMM
    z
    )

  set(NIFTYSEG_INCLUDE_DIR 
    ${NIFTYSEG_DIR}/include
    )
  set(NIFTYSEG_LIBRARY_DIR 
    ${NIFTYSEG_DIR}/lib
    )
  
else()

  # Eigen not found, explain to the user how to specify its location.
  if(NIFTYSEG_FIND_REQUIRED)
    message(FATAL_ERROR ${NIFTYSEG_DIR_MESSAGE})
  else()
    if(NOT NIFTYSEG_FIND_QUIETLY)
      message(STATUS ${NIFTYSEG_DIR_MESSAGE})
    endif()
  endif()

endif()


