art_make(
MODULE_LIBRARIES
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_ROOT_IO_TFILE_SUPPORT}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  nusimdata_SimulationBase
  cetlib cetlib_except
  larsim_MCCheater_ParticleInventoryService_service
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${ROOT_BASIC_LIB_LIST}

  larcore_Geometry_Geometry_service
  larcorealg_Geometry
  lardataobj_RawData
  lardataobj_RecoBase
  lardataobj_MCBase
  lardataobj_Simulation
  nusimdata_SimulationBase
  nug4_ParticleNavigation

  sbndcode_OpDetSim
  sbndcode_OpDetReco_OpFlash_FlashFinder

  larsim_Utils
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_PhotonBackTrackerService_service
)
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)

add_subdirectory(fcl)
