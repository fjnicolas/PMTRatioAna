////////////////////////////////////////////////////////////////////////
// Class:       PMTRatioAna
// File:        PMTRatioAna_module.cc
//
// Input: single events
////////////////////////////////////////////////////////////////////////


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "fhiclcpp/ParameterSet.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "larcore/Geometry/Geometry.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "TVector3.h"

//#include "larsim/MCCheater/ParticleInventory.h"
//#include "larsim/MCCheater/ParticleInventoryService.h"


#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"


#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"



#include "TTree.h"
#include "TFile.h"
#include "TInterpreter.h"
#include "TTimeStamp.h"

#include <vector>
#include <limits>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>

#define fXFidCut 200
#define fYFidCut 200
#define fZFidCut1 5
#define fZFidCut2 495

namespace pmtratio {
  class PMTRatioAna
;
}


class pmtratio::PMTRatioAna : public art::EDAnalyzer {
public:
  explicit PMTRatioAna
(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PMTRatioAna(PMTRatioAna const&) = delete;
  PMTRatioAna(PMTRatioAna&&) = delete;
  PMTRatioAna & operator=(PMTRatioAna const&) = delete;
  PMTRatioAna & operator=(PMTRatioAna &&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  void resetVars();
  std::string GetTPCLabel(int Cryo, int Tpc);
  void FillEnergyDepositions(std::vector<art::Ptr< sim::SimEnergyDeposit> > SimED);

  std::string fMCTruthLabel;
  std::string fSimEnergyDepositLabel;
  std::string fSimEnergyDepositInstanceLabel;
  std::vector<std::string> fSimPhotonsModuleLabel;
  std::vector<std::string> fOpFlashesModuleLabel;
  bool fSaveTruth;
  bool fSaveSimED;
  bool fApplyFiducialCut;

  TTree* fTree;
  int fEventID, fRunID, fSubRunID;

  //True variables
  std::vector<int> fTruePrimariesPDG;
  std::vector<double> fTruePrimariesE;
  double fTrueVx;
  double fTrueVy;
  double fTrueVz;
  double fTrueVt;
  int fTrueVU;
  int fTrueVV;
  int fTrueVC;
  int fTrueVTimeTick;
  double fTrueVEnergy;

  //True SimEnergyDeposits
  std::vector<double> fEnDepE;
  std::vector<double> fEnDepX;
  std::vector<double> fEnDepY;
  std::vector<double> fEnDepZ;
  std::vector<double> fEnDepT;

  //SimPhotons
  double fSimPhotonsCoated;
  double fSimPhotonsUncoated;
  double fPMTRatioSimPhotons;

  //OpFlash
  double fPECoated;
  double fPEUncoated;
  double fPMTRatioPE;

  int fNAnalyzedEvents;

  double fTotalEnergyDep;
  double fdEPromX, fdEPromY, fdEPromZ;
  double fdESpreadX, fdESpreadY, fdESpreadZ;
  double fdETPCBalance;

  /*std::vector<double> fdEtpc, fdEpromx, fdEpromy, fdEpromz;
  std::vector<double> fdEspreadx, fdEspready, fdEspreadz;
  std::vector<std::vector<double>> fdElowedges, fdEmaxedges;*/

  //Geometry
  const geo::GeometryCore* fGeom = art::ServiceHandle<geo::Geometry>()->provider();
  //PDS map
  opdet::sbndPDMapAlg fPDSMap; //map for photon detector types
};


pmtratio::PMTRatioAna::PMTRatioAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fMCTruthLabel( p.get<std::string>("MCTruthLabel", "generator") ),
  fSimEnergyDepositLabel( p.get<std::string>("SimEnergyDepositLabel", "ionandscint") ),
  fSimEnergyDepositInstanceLabel( p.get<std::string>("SimEnergyDepositInstanceLabel", "priorSCE") ),
  fSimPhotonsModuleLabel ( p.get<std::vector<std::string>>("SimPhotonsModuleLabel",   {"pdfastsim", "pdfastsimout"}) ),
  fOpFlashesModuleLabel ( p.get<std::vector<std::string>>("OpFlashesModuleLabel",   {"opflashtpc0", "opflashtpc1"}) ),
  fSaveTruth( p.get<bool>("SaveTruth", "true") ),
  fSaveSimED( p.get<bool>("SaveSimED", "true") ),
  fApplyFiducialCut( p.get<bool>("ApplyFiducialCut", "true") )
  // More initializers here.
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();
  //auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();



}


void pmtratio::PMTRatioAna::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  std::cout<<"Running PMTRatioAna---run="<< e.id().run()<<" --subrun="<< e.id().subRun()<<" --event="<<e.id().event()<<"\n";

  //............................Event General Info
  fNAnalyzedEvents++;
  fRunID = e.id().run();
  fSubRunID = e.id().subRun();
  fEventID = e.id().event();
  //Reset tree variables
  resetVars();

  //............................Read Truth Objects
  if(fSaveTruth){
    art::Handle<std::vector<simb::MCTruth>> mctruths;
    e.getByLabel(fMCTruthLabel, mctruths);
    std::cout<<" --- Saving MCTruth\n";

    for (auto const& truth : *mctruths) {
      for (int p = 0; p < truth.NParticles(); p++){
        simb::MCParticle const& TruePart = truth.GetParticle(p);
        if( TruePart.StatusCode()==1 ){
          fTruePrimariesPDG.push_back( TruePart.PdgCode() );
          fTruePrimariesE.push_back( TruePart.E() );
        }

        if( TruePart.Mother()==-1 && abs(TruePart.PdgCode())==13 ){
          fTrueVx=TruePart.EndX();
          fTrueVy=TruePart.EndY();
          fTrueVz=TruePart.EndZ();
          fTrueVt=TruePart.T();
          fTrueVEnergy=TruePart.E();
          std::cout<<"  -- Vertex: "<<fTrueVx<<" "<<fTrueVy<<" "<<fTrueVz<<std::endl;
        }
      }
    }
  }


  //............................Read SimEnergyDeposits
  art::Handle<std::vector<sim::SimEnergyDeposit> > SimEDHandle;
  std::vector<art::Ptr< sim::SimEnergyDeposit> > SimEDList;
  //if (e.getByLabel(fSimEnergyDepositLabel, fSimEnergyDepositInstanceLabel, SimEDHandle)) {
  if (e.getByLabel(fSimEnergyDepositLabel, fSimEnergyDepositInstanceLabel, SimEDHandle)) {
    art::fill_ptr_vector(SimEDList, SimEDHandle);
  }
  std::cout<<"  ---- Reading SimEnergyDeposition from handle: "<<SimEDHandle.provenance()->moduleLabel();
  std::cout<<":"<<SimEDHandle.provenance()->productInstanceName()<<" ----\n";

  FillEnergyDepositions(SimEDList);

  if(fSaveSimED){
    for (auto const& SimED : SimEDList){
      fEnDepE.push_back(SimED->Energy());
      fEnDepX.push_back(SimED->MidPointX());
      fEnDepY.push_back(SimED->MidPointY());
      fEnDepZ.push_back(SimED->MidPointZ());
      fEnDepT.push_back( (SimED->StartT()+SimED->EndT())/2. );
    }
  }




  std::vector<art::Handle<std::vector<sim::SimPhotonsLite> >> fLitePhotonHandle_list;
  fLitePhotonHandle_list = e.getMany<std::vector<sim::SimPhotonsLite>>();
  for ( const art::Handle<std::vector<sim::SimPhotonsLite>>& fLitePhotonHandle: (fLitePhotonHandle_list) ){
    if( std::find(fSimPhotonsModuleLabel.begin(), fSimPhotonsModuleLabel.end(), fLitePhotonHandle.provenance()->moduleLabel() )==fSimPhotonsModuleLabel.end()) continue;
    std::vector<art::Ptr<sim::SimPhotonsLite> > fLitePhotonList;
    art::fill_ptr_vector(fLitePhotonList, fLitePhotonHandle);
    bool reflected = (fLitePhotonHandle.provenance()->productInstanceName()== "Reflected");
    std::cout<<"*********SimPhotonLite List: "<<fLitePhotonHandle.provenance()->moduleLabel()<<" "<<fLitePhotonHandle.provenance()->productInstanceName()<<" "<<reflected<<std::endl;

    for ( auto const& fLitePhotons : (*fLitePhotonHandle) ){
      int fOpCh=fLitePhotons.OpChannel;
      std::string pd_type=fPDSMap.pdType(fOpCh);

      if(pd_type=="xarapuca_vuv" || pd_type=="xarapuca_vis") continue;
      if(!reflected && pd_type=="pmt_uncoated") continue;

      std::map<int, int> fLitePhotons_map = fLitePhotons.DetectedPhotons;

      for(auto& fphoton : fLitePhotons_map){
        if(pd_type=="pmt_coated") fSimPhotonsCoated+=fphoton.second;
        else if(pd_type=="pmt_uncoated") fSimPhotonsUncoated+=fphoton.second;
      }
    }
  }
  if(fSimPhotonsCoated!=0)
    fPMTRatioSimPhotons=4*fSimPhotonsUncoated/fSimPhotonsCoated;
  else
    fPMTRatioSimPhotons=-1;

  //Saving OpFlashes
  std::cout<<"Reading OpFlashes..."<<std::endl;
  art::Handle< std::vector<recob::OpFlash> > opflashListHandle;

  for (size_t s = 0; s < fOpFlashesModuleLabel.size(); s++) {
    std::cout<<"  --OpFlash Module Label:"<<fOpFlashesModuleLabel[s]<<std::endl;
    e.getByLabel(fOpFlashesModuleLabel[s], opflashListHandle);

    for (unsigned int i = 0; i < opflashListHandle->size(); ++i) {
      // Get OpFlash
      art::Ptr<recob::OpFlash> FlashPtr(opflashListHandle, i);
      recob::OpFlash Flash = *FlashPtr;

      double flash_t0 = Flash.AbsTime();
      std::vector<double> FlashPE_v = Flash.PEs();

      std::cout<<"   +++ OpFlash Time[ns]="<<1000*flash_t0<<" PE="<<Flash.TotalPE()<<std::endl;

      for(size_t OpCh=0; OpCh<FlashPE_v.size(); OpCh++){
        std::string pd_type=fPDSMap.pdType(OpCh);
        if(pd_type=="xarapuca_vuv" || pd_type=="xarapuca_vis") continue;

        if(pd_type=="pmt_coated") fPECoated+=FlashPE_v[OpCh];
        else if(pd_type=="pmt_uncoated") fPEUncoated+=FlashPE_v[OpCh];

      }
    }
  }
  if(fPECoated!=0)
    fPMTRatioPE=4*fPEUncoated/fPECoated;
  else
    fPMTRatioPE=-1;

  fTree->Fill();
}

void pmtratio::PMTRatioAna::resetVars()
{
  if(fSaveTruth){
    fTruePrimariesPDG.clear();
    fTruePrimariesE.clear();
    fTrueVx=-1e3;
    fTrueVy=-1e3;
    fTrueVz=-1e3;
    fTrueVt=-1e3;
    fTrueVEnergy=-1e3;
  }

  if(fSaveSimED){
    fEnDepE.clear();
    fEnDepX.clear();
    fEnDepY.clear();
    fEnDepZ.clear();
    fEnDepT.clear();
  }

  fSimPhotonsCoated=0;
  fSimPhotonsUncoated=0;
  fPMTRatioSimPhotons=-1;

  fPECoated=0;
  fPEUncoated=0;
  fPMTRatioPE=-1;

  fdEPromX=fdEPromY=fdEPromZ=0;
  fdESpreadX=fdESpreadY=fdESpreadZ=0;
  fTotalEnergyDep=0;

}



void pmtratio::PMTRatioAna::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  fTree=tfs->make<TTree>("PMTRatioTree", "Tree for PMTRatio calibration");
  fTree->Branch("RunID", &fRunID, "RunID/I");
  fTree->Branch("SubRunID", &fSubRunID, "SubRunID/I");
  fTree->Branch("EventID", &fEventID, "EventID/I");

  if(fSaveTruth){
    fTree->Branch("TruePrimariesPDG", &fTruePrimariesPDG);
    fTree->Branch("TruePrimariesE", &fTruePrimariesE);
    fTree->Branch("TrueVx", &fTrueVx, "TrueVx/D");
    fTree->Branch("TrueVy", &fTrueVy, "TrueVy/D");
    fTree->Branch("TrueVz", &fTrueVz, "TrueVz/D");
    fTree->Branch("TrueVt", &fTrueVt, "TrueVt/D");
    fTree->Branch("TrueVEnergy", &fTrueVEnergy, "TrueVEnergy/D");
  }

  if(fSaveSimED){
    fTree->Branch("EnDepE", &fEnDepE);
    fTree->Branch("EnDepX", &fEnDepX);
    fTree->Branch("EnDepY", &fEnDepY);
    fTree->Branch("EnDepZ", &fEnDepZ);
    fTree->Branch("EnDepT", &fEnDepT);
  }

  fTree->Branch("SimPhotonsCoated", &fSimPhotonsCoated);
  fTree->Branch("SimPhotonsUncoated", &fSimPhotonsUncoated);
  fTree->Branch("PMTRatioSimPhotons", &fPMTRatioSimPhotons);

  fTree->Branch("PECoated", &fPECoated);
  fTree->Branch("PEUncoated", &fPEUncoated);
  fTree->Branch("PMTRatioPE", &fPMTRatioPE);

  fTree->Branch("TotalEnergyDep", &fTotalEnergyDep);
  fTree->Branch("dETPCBalance", &fdETPCBalance);
  fTree->Branch("dEPromX", &fdEPromX);
  fTree->Branch("dEPromY", &fdEPromY);
  fTree->Branch("dEPromZ", &fdEPromZ);
  fTree->Branch("dESpreadX", &fdESpreadX);
  fTree->Branch("dESpreadY", &fdESpreadY);
  fTree->Branch("dESpreadZ", &fdESpreadZ);

  fNAnalyzedEvents=0;
}

void pmtratio::PMTRatioAna::endJob(){
}


std::string pmtratio::PMTRatioAna::GetTPCLabel(int Cryo, int Tpc){
  return "C"+std::to_string(Cryo)+"_TPC"+std::to_string(Tpc);
}

void pmtratio::PMTRatioAna::FillEnergyDepositions(std::vector<art::Ptr< sim::SimEnergyDeposit> > SimED){
  //TODO: change this harcoded
  std::unordered_map<std::string, double> EnergyDep;
  std::unordered_map<std::string, double> EnergyDepX;
  std::unordered_map<std::string, double> EnergyDepY;
  std::unordered_map<std::string, double> EnergyDepZ;
  double EnergyDepSpX=0, EnergyDepSpY=0, EnergyDepSpZ=0;
  EnergyDep["C0_TPC0"]=0; EnergyDep["C0_TPC1"]=0;
  EnergyDepX["C0_TPC0"]=0; EnergyDepX["C0_TPC1"]=0;
  EnergyDepY["C0_TPC0"]=0; EnergyDepY["C0_TPC1"]=0;
  EnergyDepZ["C0_TPC0"]=0; EnergyDepZ["C0_TPC1"]=0;

  for(auto & ed:SimED){
    //double const pos[]={ed->MidPointX(),ed->MidPointY(),ed->MidPointZ()};
    geo::Point_t const pos{ed->MidPointX(),ed->MidPointY(),ed->MidPointZ()};
    std::string tpclabel=GetTPCLabel(fGeom->FindTPCAtPosition(pos).Cryostat, fGeom->FindTPCAtPosition(pos).TPC);
    EnergyDep[tpclabel]+=ed->Energy();
    EnergyDepX[tpclabel]+=pos.X()*ed->Energy();
    EnergyDepY[tpclabel]+=pos.Y()*ed->Energy();
    EnergyDepZ[tpclabel]+=pos.Z()*ed->Energy();
    EnergyDepSpX+=pos.X()*pos.X()*ed->Energy();
    EnergyDepSpY+=pos.Y()*pos.Y()*ed->Energy();
    EnergyDepSpZ+=pos.Z()*pos.Z()*ed->Energy();
  }



  fTotalEnergyDep=0;
  for(auto & tpc_endep:EnergyDep){
    std::cout<<"JU "<<EnergyDepX[tpc_endep.first]<<std::endl;
    fTotalEnergyDep+=tpc_endep.second;
    fdEPromX+=EnergyDepX[tpc_endep.first];
    fdEPromY+=EnergyDepY[tpc_endep.first];
    fdEPromZ+=EnergyDepZ[tpc_endep.first];
    EnergyDepX[tpc_endep.first]/=tpc_endep.second;
    EnergyDepY[tpc_endep.first]/=tpc_endep.second;
    EnergyDepZ[tpc_endep.first]/=tpc_endep.second;
  }

  fdEPromX/=fTotalEnergyDep;
  fdEPromY/=fTotalEnergyDep;
  fdEPromZ/=fTotalEnergyDep;
  fdETPCBalance=std::abs(EnergyDep["C0_TPC0"]-EnergyDep["C0_TPC1"])/fTotalEnergyDep;

  fdESpreadX = std::sqrt(EnergyDepSpX/fTotalEnergyDep-fdEPromX*fdEPromX);
  fdESpreadY = std::sqrt(EnergyDepSpY/fTotalEnergyDep-fdEPromY*fdEPromY);
  fdESpreadZ = std::sqrt(EnergyDepSpZ/fTotalEnergyDep-fdEPromZ*fdEPromZ);
}

DEFINE_ART_MODULE(pmtratio::PMTRatioAna)
