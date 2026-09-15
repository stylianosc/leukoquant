set(SEG_DEV_APPLICATIONS
  seg_CalcTopNCC
  seg_DetectSulci 
  seg_PBSI  
  seg_Denoise
  seg_GetLayers 
  seg_DTI_seed 
  )

foreach( executable ${SEG_DEV_APPLICATIONS} )
  add_executable(${executable}
    ${executable}.cpp
    )
  target_link_libraries(${executable}
    ${NIFTYREG_LIBRARIES}
    ${NIFTYSEG_LIBRARIES}
    )
  install(TARGETS ${executable}
    BUNDLE DESTINATION bin
    RUNTIME DESTINATION bin
    )
endforeach()


add_executable(seg_GIF
  seg_GIF.cpp
  )
target_link_libraries(seg_GIF
  _seg_GIF
  )
install(TARGETS seg_GIF
  BUNDLE DESTINATION bin
  RUNTIME DESTINATION bin
  )

