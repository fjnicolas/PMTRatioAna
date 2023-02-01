set(      MODULE_LIBRARIES
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    nusimdata::SimulationBase
    larsim::MCCheater_ParticleInventoryService_service
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${ROOT_BASIC_LIB_LIST}

    art_root_io::tfile_support
    art_root_io::TFileService_service

    ROOT::Tree
    ROOT::Core

    larcore::Geometry_Geometry_service
    larcorealg::Geometry
    lardataobj::RawData
    lardataobj::RecoBase
    lardataobj::MCBase
    lardataobj::Simulation
    nusimdata::SimulationBase
    nug4::ParticleNavigation

    sbndcode::OpDetSim

    larsim::Utils
    larsim::MCCheater_BackTrackerService_service
    larsim::MCCheater_PhotonBackTrackerService_service
)


cet_build_plugin(PMTRatioAna art::Module SOURCE PMTRatioAna_module.cc LIBRARIES ${MODULE_LIBRARIES} )

add_subdirectory(fcl)

install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)
