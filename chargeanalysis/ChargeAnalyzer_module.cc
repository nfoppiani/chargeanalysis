////////////////////////////////////////////////////////////////////////
// Class:       ChargeAnalyzer
// Plugin Type: analyzer (art v2_05_01)
// File:        ChargeAnalyzer_module.cc
//
// Generated at Wed Feb 14 10:17:20 2018 by Nicolo Foppiani using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "uboone/Database/TPCEnergyCalib/TPCEnergyCalibService.h"

#include "EnergyHelper.h"
#include "GeometryHelper.h"
#include "McPfpMatch.h"

#include "TTree.h"

class ChargeAnalyzer;

class ChargeAnalyzer : public art::EDAnalyzer
{
public:
  explicit ChargeAnalyzer(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ChargeAnalyzer(ChargeAnalyzer const &) = delete;
  ChargeAnalyzer(ChargeAnalyzer &&) = delete;
  ChargeAnalyzer &operator=(ChargeAnalyzer const &) = delete;
  ChargeAnalyzer &operator=(ChargeAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;
  void reconfigure(fhicl::ParameterSet const &p) override;
  void clear();

  // double GetChargeCorrection(int plane, double x, double y, double z);

private:
  std::string _pfp_producer;

  // reco-true matching tools
  ubana::McPfpMatch _mcpfpMatcher;
  std::string _spacepointLabel;
  std::string _hitfinderLabel;
  std::string _geantModuleLabel;
  std::string _mcpHitAssLabel;
  bool _use_premade_ass;

  // dqdx tools
  lee::EnergyHelper energyHelper;
  lee::GeometryHelper geoHelper;

  TTree *fChargeTree;
  int fRun, fSubrun, fEvent;
  int fPdgCode, fPdgCodeParent;
  double fStartx, fStarty, fStartz;
  double fEndx, fEndy, fEndz;
  double fDirectionx, fDirectiony, fDirectionz;
  int fNhits, fNclusters;
  int fNhitsU, fNhitsV, fNhitsY;
  double fLength, fSpLength;

  std::vector<double> fSpace_points_path;
  std::vector<double> fX;
  std::vector<double> fY;
  std::vector<double> fZ;
  std::vector<double> fChargeU;
  std::vector<double> fChargeV;
  std::vector<double> fChargeY;
  std::vector<double> fSigmaChargeU;
  std::vector<double> fSigmaChargeV;
  std::vector<double> fSigmaChargeY;

  // const lariov::TPCEnergyCalibProvider &energyCalibProvider = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();
  // std::vector<double> fChargeU_corrected;
  // std::vector<double> fChargeV_corrected;
  // std::vector<double> fChargeY_corrected;
  // std::vector<double> fSigmaChargeU_corrected;
  // std::vector<double> fSigmaChargeV_corrected;
  // std::vector<double> fSigmaChargeY_corrected;

  // reco-true matching information
  int fMatchedPdgCode;
  double fMatchedE;
  double fMatchedPx, fMatchedPy, fMatchedPz;
  double fMatchedVx, fMatchedVy, fMatchedVz;
  double fMatchedEndx, fMatchedEndy, fMatchedEndz;

  // dqdx information
  double _dQdxRectangleLength;
  double _dQdxRectangleWidth;
  std::vector<double> fDQdx_hits;
  std::vector<int> fDQdx_wires;
  std::vector<double> fDQdx;
  double fDQdx_U, fDQdx_V, fDQdx_Y;
  int fn_hits_dQdx;
  double fBox_start_z, fBox_start_x, fBox_direction_z, fBox_direction_x;
  double fReco_energy_U, fReco_energy_V, fReco_energy_Y;
};

ChargeAnalyzer::ChargeAnalyzer(fhicl::ParameterSet const &p)
    : EDAnalyzer(p) // ,
                    // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;

  this->reconfigure(p);

  fChargeTree = tfs->make<TTree>("Charge", "Charge PF Particles Tree");
  fChargeTree->Branch("event", &fEvent, "event/i");
  fChargeTree->Branch("run", &fRun, "run/i");
  fChargeTree->Branch("subrun", &fSubrun, "subrun/i");
  fChargeTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
  fChargeTree->Branch("pdg_code_parent", &fPdgCodeParent, "pdg_code_parent/i");
  fChargeTree->Branch("start_x", &fStartx, "start_x/d");
  fChargeTree->Branch("start_y", &fStarty, "start_y/d");
  fChargeTree->Branch("start_z", &fStartz, "start_z/d");
  fChargeTree->Branch("end_x", &fEndx, "end_x/d");
  fChargeTree->Branch("end_y", &fEndy, "end_y/d");
  fChargeTree->Branch("end_z", &fEndz, "end_z/d");
  fChargeTree->Branch("dir_x", &fDirectionx, "dir_x/d");
  fChargeTree->Branch("dir_y", &fDirectiony, "dir_y/d");
  fChargeTree->Branch("dir_z", &fDirectionz, "dir_z/d");
  fChargeTree->Branch("n_hits", &fNhits, "n_hits/i");
  fChargeTree->Branch("n_hitsU", &fNhitsU, "n_hitsU/i");
  fChargeTree->Branch("n_hitsV", &fNhitsV, "n_hitsV/i");
  fChargeTree->Branch("n_hitsY", &fNhitsY, "n_hitsY/i");
  fChargeTree->Branch("n_clusters", &fNclusters, "n_clusters/i");
  
  fChargeTree->Branch("length", &fLength, "length/d");
  fChargeTree->Branch("sp_length", &fSpLength, "sp_length/d");
  
  fChargeTree->Branch("space_points_path", "std::vector<double>", &fSpace_points_path);
  fChargeTree->Branch("X", "std::vector<double>", &fX);
  fChargeTree->Branch("Y", "std::vector<double>", &fY);
  fChargeTree->Branch("Z", "std::vector<double>", &fZ);
  fChargeTree->Branch("chargeU", "std::vector<double>", &fChargeU);
  fChargeTree->Branch("chargeV", "std::vector<double>", &fChargeV);
  fChargeTree->Branch("chargeY", "std::vector<double>", &fChargeY);
  fChargeTree->Branch("sigma_chargeU", "std::vector<double>", &fSigmaChargeU);
  fChargeTree->Branch("sigma_chargeV", "std::vector<double>", &fSigmaChargeV);
  fChargeTree->Branch("sigma_chargeY", "std::vector<double>", &fSigmaChargeY);

  // // fChargeTree->Branch("chargeU_corrected", "std::vector<double>", &fChargeU_corrected);
  // // fChargeTree->Branch("chargeV_corrected", "std::vector<double>", &fChargeV_corrected);
  // // fChargeTree->Branch("chargeY_corrected", "std::vector<double>", &fChargeY_corrected);
  // // fChargeTree->Branch("sigma_chargeU_corrected", "std::vector<double>", &fSigmaChargeU_corrected);
  // // fChargeTree->Branch("sigma_chargeV_corrected", "std::vector<double>", &fSigmaChargeV_corrected);
  // // fChargeTree->Branch("sigma_chargeY_corrected", "std::vector<double>", &fSigmaChargeY_corrected);

  // reco-true matching information
  fChargeTree->Branch("matched_pdg_code", &fMatchedPdgCode, "matched_pdg_code/i");
  fChargeTree->Branch("matched_E", &fMatchedE, "matched_E/d");
  fChargeTree->Branch("matched_px", &fMatchedPx, "matched_px/d");
  fChargeTree->Branch("matched_py", &fMatchedPy, "matched_py/d");
  fChargeTree->Branch("matched_pz", &fMatchedPz, "matched_pz/d");
  fChargeTree->Branch("matched_vx", &fMatchedVx, "matched_vx/d");
  fChargeTree->Branch("matched_vy", &fMatchedVy, "matched_vy/d");
  fChargeTree->Branch("matched_vz", &fMatchedVz, "matched_vz/d");
  fChargeTree->Branch("matched_endx", &fMatchedEndx, "matched_endx/d");
  fChargeTree->Branch("matched_endy", &fMatchedEndy, "matched_endy/d");
  fChargeTree->Branch("matched_endz", &fMatchedEndz, "matched_endz/d");

  // dqdx information
  fChargeTree->Branch("dQdx_hits", "std::vector<double>", &fDQdx_hits);
  fChargeTree->Branch("dQdx_wires", "std::vector<int>", &fDQdx_wires);
  fChargeTree->Branch("dQdx", "std::vector<double>", &fDQdx);
  fChargeTree->Branch("dQdx_U", &fDQdx_U, "dQdx_U/d");
  fChargeTree->Branch("dQdx_V", &fDQdx_V, "dQdx_V/d");
  fChargeTree->Branch("dQdx_Y", &fDQdx_Y, "dQdx_Y/d");
  fChargeTree->Branch("n_hits_dQdx", &fn_hits_dQdx, "n_hits_dQdx/i");
  fChargeTree->Branch("box_start_z", &fBox_start_z, "box_start_z/d");
  fChargeTree->Branch("box_start_x", &fBox_start_x, "box_start_x/d");
  fChargeTree->Branch("box_direction_z", &fBox_direction_z, "box_direction_z/d");
  fChargeTree->Branch("box_direction_x", &fBox_direction_x, "box_direction_x/d");
  fChargeTree->Branch("reco_energy_U", &fReco_energy_U, "reco_energy_U/d");
  fChargeTree->Branch("reco_energy_V", &fReco_energy_V, "reco_energy_V/d");
  fChargeTree->Branch("reco_energy_Y", &fReco_energy_Y, "reco_energy_Y/d");
  
}

// double ChargeAnalyzer::GetChargeCorrection(int plane, double x, double y, double z)
// {
//   double x_correction, yz_correction, correction;
//   yz_correction = energyCalibProvider.YZdqdxCorrection(plane, y, z);
//   x_correction = energyCalibProvider.XdqdxCorrection(plane, x);
//   if (!yz_correction) yz_correction = 1.0;
//   if (!x_correction) x_correction = 1.0;
//   correction = yz_correction * x_correction;
//   return correction;
// }

void ChargeAnalyzer::clear()
{
  fSpace_points_path.clear();
  fX.clear();
  fY.clear();
  fZ.clear();
  fChargeU.clear();
  fChargeV.clear();
  fChargeY.clear();
  fSigmaChargeU.clear();
  fSigmaChargeV.clear();
  fSigmaChargeY.clear();

  // fChargeU_corrected.clear();
  // fChargeV_corrected.clear();
  // fChargeY_corrected.clear();
  // fSigmaChargeU_corrected.clear();
  // fSigmaChargeV_corrected.clear();
  // fSigmaChargeY_corrected.clear();

  fDQdx_hits.clear();
  fDQdx_wires.clear();
  fDQdx.clear();
}

void ChargeAnalyzer::reconfigure(fhicl::ParameterSet const &p)
{
  _pfp_producer                   = p.get<std::string>("PFParticleProducer", "pandoraNu");
  _hitfinderLabel                 = p.get<std::string>("HitProducer", "pandoraCosmicHitRemoval");
  _geantModuleLabel               = p.get<std::string>("GeantModule", "largeant");
  _spacepointLabel                = p.get<std::string>("SpacePointProducer", "pandoraNu");
  _mcpHitAssLabel                 = p.get<std::string>("MCPHitAssProducer", "crHitRemovalTruthMatch");
  _use_premade_ass                = p.get<bool>("UsePremadeMCPHitAss", true);

  _dQdxRectangleWidth             = p.get<double>("dQdxRectangleWidth", 1);
  _dQdxRectangleLength            = p.get<double>("dQdxRectangleLength", 4);
}

void ChargeAnalyzer::analyze(art::Event const &evt)
{
  clear();
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();
  // std::cout << "Event: " << fEvent << ", Run: " << fRun << ", subRun: " << fSubrun << std::endl;

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  
  std::vector<art::Ptr<recob::PFParticle>> pfp_v;
  art::fill_ptr_vector(pfp_v, pfparticle_handle);
  
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);

  art::FindManyP<recob::Shower> showers_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Track> tracks_per_pfpart(pfparticle_handle, evt, _pfp_producer);

  art::FindManyP<recob::Vertex> vertices_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_cluster(cluster_handle, evt, _pfp_producer);
  auto const &spacepoint_handle =
        evt.getValidHandle<std::vector<recob::SpacePoint>>(_pfp_producer);

  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, _pfp_producer);

  // for (size_t i_pfp = 0; i_pfp < pfp_v.size(); i_pfp++)
  // {
  //   art::Ptr<recob::PFParticle> const pfparticle = pfp_v.at(i_pfp);
  //   if (pfparticle->Self() != i_pfp)
  //   {
  //     std::cout << "Self != pfparticle, " << pfparticle->Self() << " " <<  i_pfp << std::endl;
  //     std::cout << "Event: " << fEvent << ", Run: " << fRun << ", subRun: " << fSubrun << std::endl;
  //   }
  // }
  // std::cout << "event ok " << std::endl;

  for (size_t i_pfp = 0; i_pfp < pfp_v.size(); i_pfp++)
  {
    art::Ptr<recob::PFParticle> const pfparticle = pfp_v.at(i_pfp);
    fPdgCode = pfparticle->PdgCode();

    if (fPdgCode != 11 && fPdgCode != 13)
    {
      continue;
    }
    
    try
    {
      double i_parent = pfparticle->Parent();
      fPdgCodeParent = pfp_v.at(i_parent)->PdgCode();
    }
    catch (...)
    {
      fPdgCodeParent = 1000000.;
    }

    // std::cout << "\n\nPfparticle:  " << i_pfp << ", pdgcode: " << fPdgCode << std::endl;

    // starting point and end point
    if (fPdgCode == 11)
    {
      std::vector<art::Ptr<recob::Shower>> pf_objs = showers_per_pfpart.at(i_pfp);
      // std::cout << "Found shower or tracks per pf part with lenght " << fNshowers_or_tracks << ", for pdg = " << fPdgCode << std::endl;
      if (pf_objs.size() != 0)
      {
        fStartx = pf_objs[0]->ShowerStart().X();
        fStarty = pf_objs[0]->ShowerStart().Y();
        fStartz = pf_objs[0]->ShowerStart().Z();
        fEndx = 1000000.;
        fEndy = 1000000.;
        fEndz = 1000000.;
        fLength = pf_objs[0]->Length();
        fDirectionx = pf_objs[0]->Direction().X();
        fDirectiony = pf_objs[0]->Direction().Y();
        fDirectionz = pf_objs[0]->Direction().Z();
        // std::cout << "shower dir: " << pf_objs[0]->Direction().Z() << " , " << pf_objs[0]->Direction().X() << " , " << pf_objs[0]->Direction().X()/pf_objs[0]->Direction().Z() << std::endl;
      }
      else
      {
        fStartx = 1000000.;
        fStarty = 1000000.;
        fStartz = 1000000.;
        fEndx = 1000000.;
        fEndy = 1000000.;
        fEndz = 1000000.;
        fLength = 1000000.;
        fDirectionx = 1000000.;
        fDirectiony = 1000000.;
        fDirectionz = 1000000.;
        // std::cout << "No start point found for shower " << fPdgCode << std::endl;
      }
    }
    else if (fPdgCode == 13)
    {
      std::vector<art::Ptr<recob::Track>> pf_objs = tracks_per_pfpart.at(i_pfp);
      // std::cout << "Found shower or tracks per pf part with lenght " << fNshowers_or_tracks << ", for pdg = " << fPdgCode << std::endl;
      if (pf_objs.size() != 0)
      {
        fStartx = pf_objs[0]->Start().X();
        fStarty = pf_objs[0]->Start().Y();
        fStartz = pf_objs[0]->Start().Z();
        fEndx = pf_objs[0]->End().X();
        fEndy = pf_objs[0]->End().Y();
        fEndz = pf_objs[0]->End().Z();
        fLength = pf_objs[0]->Length();
        fDirectionx = pf_objs[0]->StartDirection().X();
        fDirectiony = pf_objs[0]->StartDirection().Y();
        fDirectionz = pf_objs[0]->StartDirection().Z();
        // std::cout << "track dir: " << pf_objs[0]->StartDirection().Z() << " , " << pf_objs[0]->StartDirection().X() << " , " << pf_objs[0]->StartDirection().X()/pf_objs[0]->StartDirection().Z() << std::endl;
      }
      else
      {
        fStartx = 1000000.;
        fStarty = 1000000.;
        fStartz = 1000000.;
        fEndx = 1000000.;
        fEndy = 1000000.;
        fEndz = 1000000.;
        fLength = 1000000.;
        fDirectionx = 1000000.;
        fDirectiony = 1000000.;
        fDirectionz = 1000000.;
        // std::cout << "No start point found for track " << fPdgCode << std::endl;
      }
    }
    else
    {
      fStartx = 1000000.;
      fStarty = 1000000.;
      fStartz = 1000000.;
      fEndx = 1000000.;
      fEndy = 1000000.;
      fEndz = 1000000.;
      fLength = 1000000.;
      fDirectionx = 1000000.;
      fDirectiony = 1000000.;
      fDirectionz = 1000000.;
    }

    // std::cout << "Startz: " << fStartz <<  ", Startx:  " << fStartx << std::endl;
    // clusters
    fNhits = 0;
    fNhitsU = 0;
    fNhitsV = 0;
    fNhitsY = 0;
    fNclusters = 0;

    std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(i_pfp);
    fNclusters = clusters.size();
    int j = 0;
    for (art::Ptr<recob::Cluster> &cluster : clusters)
    {
      fNhits += cluster->NHits();
      // if (cluster->Plane().Plane !=2)
      //   continue;
      // std::cout << "Cluster:  " << j << ", plane: " << cluster->Plane().Plane << ", nhits: " << fNhits << std::endl;
      // std::cout << "Cluster_start:  " << 0.3*cluster->StartWire() << ", cluster end: " << 0.3*cluster->EndWire() << std::endl;
      // if (cluster->EndWire() <= cluster->StartWire())
      //   std::cout << "endwire smaller than startwire" << std::endl;

      // std::cout << "Cluster_start time:  " << cluster->StartTick() << ", cluster end time: " << cluster->EndTick() << std::endl;
      // if (cluster->EndTick() <= cluster->StartTick())
      //   std::cout << "endtick smaller than starttick" << std::endl;

      // std::cout << "Cluster_startaxis:  " << cos(cluster->StartAngle()) << " , " << sin(cluster->StartAngle()) << " , " << tan(cluster->StartAngle()) << std::endl;
      // if (cos(cluster->StartAngle()) <= 0)
      //   std::cout << "cluster Startangle < 0" << std::endl;
      
      // std::cout << "Cluster_endaxis:  " << cos(cluster->EndAngle()) << " , " << sin(cluster->EndAngle()) << " , " << tan(cluster->EndAngle()) << std::endl;
      // if (cos(cluster->EndAngle()) <= 0)
      //   std::cout << "cluster Endangle < 0" << std::endl;

      std::vector<art::Ptr<recob::Hit>> hits = hits_per_cluster.at(cluster.key());
      for (art::Ptr<recob::Hit> &hit : hits)
      {
        // std::cout << "hit on wire " << hit->WireID().Wire << " , " << 0.3*hit->WireID().Wire << " , " << hit->PeakTime() << std::endl;
        int fPlane = hit->WireID().Plane;
        if (fPlane == 0)
        {
          fNhitsU += 1;
        }
        else if (fPlane == 1)
        {
          fNhitsV += 1;
        }
        else if (fPlane == 2)
        {
          fNhitsY += 1;
        }
      }
      j++;
    }

    // dqdx
    fDQdx = {-1., -1., -1.};
    std::vector<double> box_start = {1000000., 1000000.};
    std::vector<double> box_direction = {1000000., 1000000.};
    energyHelper.dQdx(i_pfp, evt, fDQdx, fDQdx_hits, fDQdx_wires, box_start, box_direction, _dQdxRectangleLength, _dQdxRectangleWidth, _pfp_producer);
    fDQdx_U = fDQdx[0];
    fDQdx_V = fDQdx[1];
    fDQdx_Y = fDQdx[2];
    fn_hits_dQdx = fDQdx_hits.size();
    fBox_start_z = box_start[0];
    fBox_start_x = box_start[1];
    fBox_direction_z = box_direction[0];
    fBox_direction_x = box_direction[1];

    std::vector<double> this_energy;
    std::vector<int> this_nhits;
    energyHelper.energyFromHits(*pfparticle, this_nhits, this_energy, evt, _pfp_producer);
    fReco_energy_U = this_energy[0];
    fReco_energy_V = this_energy[1];
    fReco_energy_Y = this_energy[2];

    // reco-true matching configuration 
    bool _is_data = evt.isRealData();
    if (_is_data) 
    {
      std::cout << "[RecoTrueMatcher] Running on a real data file. No MC-PFP matching will be attempted." << std::endl;
    } 
  
    if (_use_premade_ass)   
      _mcpfpMatcher.Configure(evt, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel, _mcpHitAssLabel, lar_pandora::LArPandoraHelper::kAddDaughters);
    else 
      _mcpfpMatcher.Configure(evt, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel);

    // This is a map: PFParticle to matched MCParticle: std::map<art::Ptr<recob::PFParticle>, art::Ptr<simb::MCParticle> >
    lar_pandora::PFParticlesToMCParticles matched_pfp_to_mcp_map;
    _mcpfpMatcher.GetRecoToTrueMatches(matched_pfp_to_mcp_map);

    // spacepoints loop
    double x, y, z, x_prev, y_prev, z_prev;
    int i = 0;

    std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(i_pfp);
    for (art::Ptr<recob::SpacePoint> &sps : spcpnts)
    {
      x_prev = x;
      y_prev = y;
      z_prev = z;
      auto xyz = sps->XYZ();
      x = xyz[0];
      y = xyz[1];
      z = xyz[2];

      fX.push_back(x);
      fY.push_back(y);
      fZ.push_back(z);

      if (i == 0)
      {
        fSpace_points_path.push_back(0);
      }
      else
      {
        double aux_d = std::sqrt((x-x_prev)*(x-x_prev) +
                                 (y-y_prev)*(y-y_prev) +
                                 (z-z_prev)*(z-z_prev));
        aux_d += fSpace_points_path.back();
        fSpace_points_path.push_back(aux_d);
      }

      std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(sps.key());
      double aux_ChargeU = 0;
      double aux_ChargeV = 0;
      double aux_ChargeY = 0;
      double aux_SigmaChargeU = 0;
      double aux_SigmaChargeV = 0;
      double aux_SigmaChargeY = 0;
      // double aux_ChargeU_corrected = 0;
      // double aux_ChargeV_corrected = 0;
      // double aux_ChargeY_corrected = 0;
      // double aux_SigmaChargeU_corrected = 0;
      // double aux_SigmaChargeV_corrected = 0;
      // double aux_SigmaChargeY_corrected = 0;
      for (art::Ptr<recob::Hit> &hit : hits)
      {
        double hit_plane = hit->WireID().Plane;
        double hit_integral = hit->Integral();
        double hit_sigma_integral = hit->SigmaIntegral();
        // double correction = GetChargeCorrection(hit_plane, x, y, z);
        if (hit_plane == 0)
        {
          aux_ChargeU += hit_integral;          
          // aux_ChargeU_corrected += aux_ChargeU*correction;
          aux_SigmaChargeU += hit_sigma_integral;
          // aux_SigmaChargeU_corrected += aux_SigmaChargeU*correction;
        }
        else if (hit_plane == 1)
        {
          aux_ChargeV += hit_integral;          
          // aux_ChargeV_corrected += aux_ChargeV*correction;
          aux_SigmaChargeV += hit_sigma_integral;
          // aux_SigmaChargeV_corrected += aux_SigmaChargeV*correction;
        }
        else if (hit_plane == 2)
        {
          aux_ChargeY += hit_integral;          
          // aux_ChargeY_corrected += aux_ChargeY*correction;
          aux_SigmaChargeY += hit_sigma_integral;
          // aux_SigmaChargeY_corrected += aux_SigmaChargeY*correction;
        }
        else
        {
          std::cout << "hit plane != 0, 1, 2, but " << hit_plane << std::endl;
        }
      }
      fChargeU.push_back(aux_ChargeU);
      fChargeV.push_back(aux_ChargeV);
      fChargeY.push_back(aux_ChargeY);
      fSigmaChargeU.push_back(aux_SigmaChargeU);
      fSigmaChargeV.push_back(aux_SigmaChargeV);
      fSigmaChargeY.push_back(aux_SigmaChargeY);

      // // fChargeU_corrected.push_back(aux_ChargeU_corrected);
      // // fChargeV_corrected.push_back(aux_ChargeV_corrected);
      // // fChargeY_corrected.push_back(aux_ChargeY_corrected);
      // // fSigmaChargeU_corrected.push_back(aux_SigmaChargeU_corrected);
      // // fSigmaChargeV_corrected.push_back(aux_SigmaChargeV_corrected);
      // // fSigmaChargeY_corrected.push_back(aux_SigmaChargeY_corrected);

      i++;
    }

    fSpLength = fSpace_points_path.back();

    // reco-true matching
    if (!_is_data) 
    {
      auto iter = matched_pfp_to_mcp_map.find(pfparticle);
      if (iter == matched_pfp_to_mcp_map.end()) 
      {
        fMatchedPdgCode = 1000000.;
        fMatchedE = 1000000.;
        fMatchedPx = 1000000.;
        fMatchedPy = 1000000.;
        fMatchedPz = 1000000.;
        fMatchedVx = 1000000.;
        fMatchedVy = 1000000.;
        fMatchedVz = 1000000.;
        fMatchedEndx = 1000000.;
        fMatchedEndy = 1000000.;
        fMatchedEndz = 1000000.;
      } 
      else 
      {
        art::Ptr<simb::MCParticle> mc_part = iter->second;
        fMatchedPdgCode = mc_part->PdgCode();
        fMatchedE = mc_part->E();
        fMatchedPx = mc_part->Px();
        fMatchedPy = mc_part->Py();
        fMatchedPz = mc_part->Pz();
        fMatchedVx = mc_part->Vx();
        fMatchedVy = mc_part->Vy();
        fMatchedVz = mc_part->Vz();
        fMatchedEndx = mc_part->EndX();
        fMatchedEndy = mc_part->EndY();
        fMatchedEndz = mc_part->EndZ();
      }
    }  
    else
    {
      fMatchedPdgCode = 1000000.;
      fMatchedE = 1000000.;
      fMatchedPx = 1000000.;
      fMatchedPy = 1000000.;
      fMatchedPz = 1000000.;
      fMatchedVx = 1000000.;
      fMatchedVy = 1000000.;
      fMatchedVz = 1000000.;
      fMatchedEndx = 1000000.;
      fMatchedEndy = 1000000.;
      fMatchedEndz = 1000000.;
    }

    fChargeTree->Fill();
    clear();
  }
}

DEFINE_ART_MODULE(ChargeAnalyzer)
