
set(SEG_GIF_LIBRARIES 
  _seg_GIF
)


foreach(library ${SEG_GIF_LIBRARIES} )

  add_library(${library} ${${PROJECT_NAME}_BUILD_TYPE}
    ${library}.cpp
    # [leukoquant patch] Provides reg_mat44_expm() and the rest of the old
    # _reg_maths_eigen.h API, which the current vendored NiftyReg no longer
    # has (removed in its post-#130 refactor) -- see
    # _reg_maths_eigen_compat.h for details.
    _reg_maths_eigen_compat.cpp
    )

  target_link_libraries(${library} 
    ${NIFTYREG_LIBRARIES}
    ${NIFTYSEG_LIBRARIES}
    )  

  install(TARGETS ${library}
    BUNDLE DESTINATION  lib
    RUNTIME DESTINATION lib
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

endforeach()


if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  file(GLOB __files1 "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
  file(GLOB __files2 "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
  file(GLOB __files3 "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
  install(FILES ${__files1} ${__files2} ${__files3}
    DESTINATION include
    COMPONENT   Development
    )
endif(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
