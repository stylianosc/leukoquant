#-----------------------------------------------------------------------------
option(USE_NII_NAN "To enable NaN support with nifti images" ON)
mark_as_advanced(FORCE USE_NII_NAN)
if(USE_NII_NAN)
    add_definitions(-DUSE_NII_NAN)
endif(USE_NII_NAN)
add_definitions(-DHAVE_ZLIB)

set(NAME _seg_nifti)

# On Windows none of this is exported, so leaving it on STATIC
# otherwise there will be no .lib file generated
if(WIN32)
  SET(LIBRARY_TYPE STATIC)
endif()
 
# [leukoquant patch] nifti1_io.c/znzlib.c here are NiftySeg's own reference
# nifticlib implementation, matched to NiftySeg's own classic nifti1_io.h.
# Since nifti1_io.h/nifti1.h/znzlib.h in this directory now forward to the
# vendored NiftyReg's RNifti-aware headers instead (see those files), this
# .c file would need the classic API surface (nifti_type_list, g_opts, etc)
# those newer headers don't provide -- it doesn't compile against them.
# Rather than porting nifti1_io.c to a different API, just link against
# NiftyReg's own already-compiled implementation of the exact same
# functions (built against the exact headers we're now using), instead of
# building a redundant, now-incompatible copy here.
add_library(${NAME} INTERFACE)
target_link_libraries(${NAME} INTERFACE reg_nifti znz z)
install(FILES nifti1_io.h znzlib.h nifti1.h DESTINATION include COMPONENT Development)
#add_dependencies(${NAME} Eigen)
#-----------------------------------------------------------------------------
