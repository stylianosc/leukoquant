#-----------------------------------------------------------------------------

    ADD_EXECUTABLE(seg_LoAd seg_LoAd.cpp)
        TARGET_LINK_LIBRARIES(seg_LoAd _seg_LoAd _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_LoAd)

    ADD_EXECUTABLE(seg_stats seg_Stats.cpp)
        TARGET_LINK_LIBRARIES(seg_stats _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_stats)

    ADD_EXECUTABLE(seg_EM seg_EM.cpp)
        TARGET_LINK_LIBRARIES(seg_EM _seg_EM _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_EM)

    ADD_EXECUTABLE(seg_LabFusion seg_LabFusion.cpp)
        TARGET_LINK_LIBRARIES(seg_LabFusion _seg_LabFusion _seg_tools _seg_Topo _seg_MRF _seg_FMM _seg_BiasCorrection _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_LabFusion)

    ADD_EXECUTABLE(seg_maths seg_maths.cpp)
        TARGET_LINK_LIBRARIES(seg_maths _seg_tools _seg_Topo _seg_MRF _seg_FMM _seg_BiasCorrection _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_maths)

    ADD_EXECUTABLE(Seg_BiASM Seg_BiASM.cpp)
        TARGET_LINK_LIBRARIES(Seg_BiASM _Seg_InitAndRun _seg_tools _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin Seg_BiASM)

    ADD_EXECUTABLE(Seg_Analysis Seg_Analysis.cpp)
        TARGET_LINK_LIBRARIES(Seg_Analysis _TreeEM_new _Seg_InitAndRun _seg_tools _seg_nifti _Seg_Analysis ${ZLIB})
    INSTALL_TARGETS(/bin Seg_Analysis)

#-----------------------------------------------------------------------------
