#-----------------------------------------------------------------------------
    # [leukoquant patch] This fork only ever needs to build the two BaMoS-specific
    # applications below (Seg_BiASM, Seg_Analysis) -- they don't exist in mainline
    # NiftySeg (confirmed via diff against niftyseg/source_latest/seg-apps/).
    # Everything else this fork's seg-apps directory could build (seg_LoAd,
    # seg_stats, seg_EM, seg_LabFusion, seg_maths) is a standard NiftySeg tool with
    # no BaMoS-specific changes; BaMoS_WMH_*.sh sources those from the separately
    # built leukoquant/external/niftyseg/bin/ instead, so building duplicate copies
    # of them from this old fork (with an old NiftySeg codebase) served no purpose.

    ADD_EXECUTABLE(Seg_BiASM Seg_BiASM.cpp)
        TARGET_LINK_LIBRARIES(Seg_BiASM _Seg_InitAndRun _seg_tools _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin Seg_BiASM)

    ADD_EXECUTABLE(Seg_Analysis Seg_Analysis.cpp)
        TARGET_LINK_LIBRARIES(Seg_Analysis _TreeEM_new _Seg_InitAndRun _seg_tools _seg_nifti _Seg_Analysis ${ZLIB})
    INSTALL_TARGETS(/bin Seg_Analysis)

#-----------------------------------------------------------------------------
