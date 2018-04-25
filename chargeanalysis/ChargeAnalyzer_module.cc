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

#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "SortingFunctions.cxx"
#include "MathFunctions.cxx"

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
  void recoTrueMatching(art::Event const &evt, art::Ptr<recob::PFParticle> const &pfparticle);
  void trueNeutrinoInformation(art::Event const &evt);
  void spacePointTrajectory(art::Event const &evt, size_t const i_pfp, double const dir[3]);


private:
  std::string _pfp_producer;

  // reco-true matching tools
  ubana::McPfpMatch _mcpfpMatcher;
  std::string _spacepointLabel;
  std::string _hitfinderLabel;
  std::string _geantModuleLabel;
  std::string _mcpHitAssLabel;
  bool _use_premade_ass;

  std::string _mctruthLabel;

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
  double fDistanceFromTrue;

  double fTrue_vx, fTrue_vy, fTrue_vz;
  std::vector<double> fTrue_v;
  int fTrue_nu_pdg, fTrue_ccnc;
  double fTrue_nu_energy, fTrue_theta, fTrue_vx_sce, fTrue_vy_sce, fTrue_vz_sce;
  std::vector<double> fTrue_v_sce;

  // space points
  std::vector<double> fDistance_projection;
  std::vector<double> fX;
  std::vector<double> fY;
  std::vector<double> fZ;
  std::vector<double> fChargeU;
  std::vector<double> fChargeV;
  std::vector<double> fChargeY;
  std::vector<double> fDqdx_U;
  std::vector<double> fDqdx_V;
  std::vector<double> fDqdx_Y;
  std::vector<int> fN_hits_U;
  std::vector<int> fN_hits_V;
  std::vector<int> fN_hits_Y;
  std::vector<double> fTrajectory_range;
  std::vector<double> fRadial_distance;
  std::vector<double> fScattering_angle;
  double fTrajectory_length;
  double fMean_angle;
  double fRms_angle;
  double fMean_radial_distance;
  double fRms_radial_distance;

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
  double fDistanceFromMatched;

  // dqdx information
  double _dQdxRectangleLength;
  double _dQdxRectangleWidth;
  std::vector<double> fDQdx_hits;
  std::vector<int> fDQdx_wires;
  std::vector<double> fDQdx;
  double fDQdx_U, fDQdx_V, fDQdx_Y;
  int fn_hits_dQdx;
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


  fChargeTree->Branch("true_nu_pdg", &fTrue_nu_pdg, "true_nu_pdg/d");
  fChargeTree->Branch("true_ccnc", &fTrue_ccnc, "true_ccnc/d");
  fChargeTree->Branch("true_nu_energy", &fTrue_nu_energy, "true_nu_energy/d");
  fChargeTree->Branch("true_theta", &fTrue_theta, "true_theta/d");

  fChargeTree->Branch("true_vx", &fTrue_vx, "true_vx/d");
  fChargeTree->Branch("true_vy", &fTrue_vy, "true_vy/d");
  fChargeTree->Branch("true_vz", &fTrue_vz, "true_vz/d");

  fChargeTree->Branch("true_vx_sce", &fTrue_vx_sce, "true_vx_sce/d");
  fChargeTree->Branch("true_vy_sce", &fTrue_vy_sce, "true_vy_sce/d");
  fChargeTree->Branch("true_vz_sce", &fTrue_vz_sce, "true_vz_sce/d");
  fChargeTree->Branch("distance_from_true", &fDistanceFromTrue, "distance_from_true/d");

  // space points
  fChargeTree->Branch("distance_projection", "std::vector<double>", &fDistance_projection);
  fChargeTree->Branch("x", "std::vector<double>", &fX);
  fChargeTree->Branch("y", "std::vector<double>", &fY);
  fChargeTree->Branch("z", "std::vector<double>", &fZ);
  fChargeTree->Branch("chargeU", "std::vector<double>", &fChargeU);
  fChargeTree->Branch("chargeV", "std::vector<double>", &fChargeV);
  fChargeTree->Branch("chargeY", "std::vector<double>", &fChargeY);
  fChargeTree->Branch("dqdx_U", "std::vector<double>", &fDqdx_U);
  fChargeTree->Branch("dqdx_V", "std::vector<double>", &fDqdx_V);
  fChargeTree->Branch("dqdx_Y", "std::vector<double>", &fDqdx_Y);
  fChargeTree->Branch("n_hits_U", "std::vector<int>", &fN_hits_U);
  fChargeTree->Branch("n_hits_V", "std::vector<int>", &fN_hits_V);
  fChargeTree->Branch("n_hits_Y", "std::vector<int>", &fN_hits_Y);
  fChargeTree->Branch("trajectory_range", "std::vector<double>", &fTrajectory_range);
  fChargeTree->Branch("radial_distance", "std::vector<double>", &fRadial_distance);
  fChargeTree->Branch("scattering_angle", "std::vector<double>", &fScattering_angle);
  
  fChargeTree->Branch("trajectory_length", &fTrajectory_length, "trajectory_length/d");
  fChargeTree->Branch("mean_angle", &fMean_angle, "mean_angle/d");
  fChargeTree->Branch("rms_angle", &fRms_angle, "rms_angle/d");
  fChargeTree->Branch("mean_radial_distance", &fMean_radial_distance, "mean_radial_distance/d");
  fChargeTree->Branch("rms_radial_distance", &fRms_radial_distance, "rms_radial_distance/d");

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
  fChargeTree->Branch("distance_from_matched", &fDistanceFromMatched, "distance_from_matched/d");

  // dqdx information
  fChargeTree->Branch("dQdx_hits", "std::vector<double>", &fDQdx_hits);
  fChargeTree->Branch("dQdx_wires", "std::vector<int>", &fDQdx_wires);
  fChargeTree->Branch("dQdx", "std::vector<double>", &fDQdx);
  fChargeTree->Branch("dQdx_U", &fDQdx_U, "dQdx_U/d");
  fChargeTree->Branch("dQdx_V", &fDQdx_V, "dQdx_V/d");
  fChargeTree->Branch("dQdx_Y", &fDQdx_Y, "dQdx_Y/d");
  fChargeTree->Branch("n_hits_dQdx", &fn_hits_dQdx, "n_hits_dQdx/i");
  fChargeTree->Branch("reco_energy_U", &fReco_energy_U, "reco_energy_U/d");
  fChargeTree->Branch("reco_energy_V", &fReco_energy_V, "reco_energy_V/d");
  fChargeTree->Branch("reco_energy_Y", &fReco_energy_Y, "reco_energy_Y/d");
  
}

void ChargeAnalyzer::clear()
{
  fDistance_projection.clear();
  fX.clear();
  fY.clear();
  fZ.clear();
  fChargeU.clear();
  fChargeV.clear();
  fChargeY.clear();
  fDqdx_U.clear();
  fDqdx_V.clear();
  fDqdx_Y.clear();
  fN_hits_U.clear();
  fN_hits_V.clear();
  fN_hits_Y.clear();
  fTrajectory_range.clear();
  fRadial_distance.clear();
  fScattering_angle.clear();


  fDQdx_hits.clear();
  fDQdx_wires.clear();
  fDQdx.clear();
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

void ChargeAnalyzer::recoTrueMatching(art::Event const &evt, art::Ptr<recob::PFParticle> const &pfparticle)
{
  bool _is_data = evt.isRealData();
  if (_is_data) 
  {
    std::cout << "[RecoTrueMatcher] Running on a real data file. No MC-PFP matching will be attempted." << std::endl;
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
    if (_use_premade_ass)   
      _mcpfpMatcher.Configure(evt, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel, _mcpHitAssLabel, lar_pandora::LArPandoraHelper::kAddDaughters);
    else 
      _mcpfpMatcher.Configure(evt, _pfp_producer, _spacepointLabel, _hitfinderLabel, _geantModuleLabel);
    
    // This is a map: PFParticle to matched MCParticle: std::map<art::Ptr<recob::PFParticle>, art::Ptr<simb::MCParticle> >
    lar_pandora::PFParticlesToMCParticles matched_pfp_to_mcp_map;
    _mcpfpMatcher.GetRecoToTrueMatches(matched_pfp_to_mcp_map);

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
      std::vector<double> fStart_true = {fMatchedVx, fMatchedVy, fMatchedVz};
      fDistanceFromMatched = geoHelper.distance(fStart_true, fTrue_v);
    }
  }  
}

void ChargeAnalyzer::trueNeutrinoInformation(art::Event const &evt)
{
  auto const &generator_handle =
    evt.getValidHandle<std::vector<simb::MCTruth>>(_mctruthLabel);
  auto const &generator(*generator_handle);

  // bool there_is_a_neutrino = false;
  //std::cout << "[PandoraLEEAnalyzer] Generator size " << generator.size() << std::endl;
  if (generator.size() > 0)
  {
    for (auto &gen : generator)
    {
    //std::cout << "[PandoraLEEAnalyzer] Generator origin " << gen.Origin() << std::endl;

      if (gen.Origin() == simb::kBeamNeutrino)
      {
        // there_is_a_neutrino = true;
        fTrue_nu_pdg = gen.GetNeutrino().Nu().PdgCode();
        fTrue_nu_energy = gen.GetNeutrino().Nu().E();
        fTrue_ccnc = gen.GetNeutrino().CCNC();
        // _qsqr = gen.GetNeutrino().QSqr();
        fTrue_theta = gen.GetNeutrino().Theta();

        // if (_ccnc == simb::kNC)
        // {
        //   fcategory = k_nc;
        // }

        fTrue_vx = gen.GetNeutrino().Nu().Vx();
        fTrue_vy = gen.GetNeutrino().Nu().Vy();
        fTrue_vz = gen.GetNeutrino().Nu().Vz();

        fTrue_v = {fTrue_vx, fTrue_vy, fTrue_vz};
        // _interaction_type = gen.GetNeutrino().InteractionType();

        auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
        if (SCE->GetPosOffsets(fTrue_vx, fTrue_vy, fTrue_vz).size() == 3)
        {
          fTrue_vx_sce =
          fTrue_vx - SCE->GetPosOffsets(fTrue_vx, fTrue_vy, fTrue_vz)[0] + 0.7;
          fTrue_vy_sce =
          fTrue_vy + SCE->GetPosOffsets(fTrue_vx, fTrue_vy, fTrue_vz)[1];
          fTrue_vz_sce =
          fTrue_vz + SCE->GetPosOffsets(fTrue_vx, fTrue_vy, fTrue_vz)[2];

          fTrue_v_sce = {fTrue_vx_sce, fTrue_vy_sce, fTrue_vz_sce};
        }
        else
        {
          std::cout << "[PandoraLEEAnalyzer] "
          << "Space Charge service offset size < 3" << std::endl;
          continue;
        }

        // if (!geoHelper.isActive(true_neutrino_vertex))
        // {
        //   _category = k_dirt;
        // }
      }
    }
  }

  // if (!there_is_a_neutrino)
  // _category = k_cosmic;

  // auto const &mcparticles_handle = evt.getValidHandle<std::vector<simb::MCParticle>>("largeant");
  // auto const &mcparticles(*mcparticles_handle);

  // for (auto &mcparticle : mcparticles)
  // {
  // if (!(mcparticle.Process() == "primary" &&
  //   mcparticle.T() != 0 &&
  //   mcparticle.StatusCode() == 1))
  // continue;

  // const auto mc_truth = pandoraHelper.TrackIDToMCTruth(evt, "largeant", mcparticle.TrackId());
  // if (mc_truth->Origin() == simb::kBeamNeutrino)
  // {
  // _nu_daughters_E.push_back(mcparticle.E());
  // _nu_daughters_pdg.push_back(mcparticle.PdgCode());

  // _nu_daughters_px.push_back(mcparticle.Px());
  // _nu_daughters_py.push_back(mcparticle.Py());
  // _nu_daughters_pz.push_back(mcparticle.Pz());

  // _nu_daughters_vx.push_back(mcparticle.Vx());
  // _nu_daughters_vy.push_back(mcparticle.Vy());
  // _nu_daughters_vz.push_back(mcparticle.Vz());

  // _nu_daughters_endx.push_back(mcparticle.EndX());
  // _nu_daughters_endy.push_back(mcparticle.EndY());
  // _nu_daughters_endz.push_back(mcparticle.EndZ());
  // }
  // }

  // //Insert block to save the start point of the MCshower object for all showers that have a neutrino as mother and a kbeamneutrino as origin
  // auto const &mcshower_handle = evt.getValidHandle<std::vector<sim::MCShower>>("mcreco");
  // for (size_t _i_mcs = 0; _i_mcs < mcshower_handle->size(); _i_mcs++)
  // {
  // int pdg_mother = mcshower_handle->at(_i_mcs).MotherPdgCode();
  // int origin = mcshower_handle->at(_i_mcs).Origin();

  // if ((pdg_mother == 22 || pdg_mother == 11) && origin == 1)
  // {
  // _true_shower_pdg.push_back(mcshower_handle->at(_i_mcs).AncestorPdgCode());
  // _true_shower_depE.push_back(mcshower_handle->at(_i_mcs).DetProfile().E());

  // double x_det = mcshower_handle->at(_i_mcs).Start().X();
  // double y_det = mcshower_handle->at(_i_mcs).Start().Y();
  // double z_det = mcshower_handle->at(_i_mcs).Start().Z();

  // if (pdg_mother == 22)
  // { //For photons take the end of the shower
  // x_det = mcshower_handle->at(_i_mcs).End().X();
  // y_det = mcshower_handle->at(_i_mcs).End().Y();
  // z_det = mcshower_handle->at(_i_mcs).End().Z();
  // }

  // std::vector<double> dqdx = mcshower_handle->at(_i_mcs).dQdx();
  // //std::vector< double > chrg = mcshower_handle->at(_i_mcs).Charge();

  // //unsigned int maxindex= (dqdx.size() > chrg.size())? chrg.size() : dqdx.size();
  // //std::cout << "[PandoraLEE] " << "dqdx.size(): " << dqdx.size() << "\t chrg.size(): " << chrg.size() << std::endl;
  // //for(unsigned int j=0; j<maxindex; j++){
  // //  std::cout << "[PandoraLEE] " << j << " dqdx: " << dqdx.at(j) << "\t chrg: " << chrg.at(j) << std::endl;
  // //}

  // auto const *SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
  // _true_shower_x_sce.push_back(x_det - SCE->GetPosOffsets(x_det, y_det, z_det)[0] + 0.7);
  // _true_shower_y_sce.push_back(y_det + SCE->GetPosOffsets(x_det, y_det, z_det)[1]);
  // _true_shower_z_sce.push_back(z_det + SCE->GetPosOffsets(x_det, y_det, z_det)[2]);

  // //std::cout << "[PandoraLEE] "
  // //    << "MCShower End: (" << x_det - SCE->GetPosOffsets(x_det, y_det, z_det)[0] + 0.7
  // //    << "," << y_det + SCE->GetPosOffsets(x_det, y_det, z_det)[1]
  // //    << "," << z_det + SCE->GetPosOffsets(x_det, y_det, z_det)[2] << ")" << std::endl;

  // //std::cout << "[PandoraLEE] "
  // //    << "TrueVTX: (" << _true_vx_sce << "," << _true_vy_sce << "," << _true_vz_sce << ")" << std::endl;
  // }
  // }

  // if (_category != k_cosmic && _category != k_dirt && _category != k_nc)
  // {
  // if (abs(_nu_pdg) == 12)
  // {
  // _category = k_nu_e;
  // }
  // if (abs(_nu_pdg) == 14)
  // {
  // _category = k_nu_mu;
  // }
  // }
  // }
  // else
  // {
  //   _gain = 240;
  //   _category = k_data;
  // }
}

void ChargeAnalyzer::spacePointTrajectory(art::Event const &evt, size_t const i_pfp, double const dir[3])
{
  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  auto const &spacepoint_handle =
        evt.getValidHandle<std::vector<recob::SpacePoint>>(_pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, _pfp_producer);
  art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt, _pfp_producer);

  auto const &vertex_obj = vertex_per_pfpart.at(i_pfp);
  double vertex[3];
  vertex_obj->XYZ(vertex);
  
  double norm_dir[3];
  geoHelper.normalize(dir, norm_dir);
  TVector3 t_dir = TVector3(norm_dir);
  double pitch[3];
  for(size_t ii=0; ii<3; ii++)
    pitch[ii] = geoHelper.getPitch(t_dir, ii);

  std::vector<double> distance_projection, x, y, z;
  std::vector<double> charge[3], dqdx[3];
  std::vector<int> n_hits[3];
  std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(i_pfp);

  for(art::Ptr<recob::SpacePoint> &sps : spcpnts)
  {
    double xyz[3] = {sps->XYZ()[0] - vertex[0],
                     sps->XYZ()[1] - vertex[1],
                     sps->XYZ()[2] - vertex[2]};
    double aux_distance = geoHelper.dotProduct(xyz, norm_dir);
    distance_projection.push_back(aux_distance);
    x.push_back(sps->XYZ()[0]);
    y.push_back(sps->XYZ()[1]);
    z.push_back(sps->XYZ()[2]);

    int aux_n_hits[3] = {0, 0, 0};
    double aux_charge[3] = {0, 0, 0};
    std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(sps.key());
    for(art::Ptr<recob::Hit> &hit : hits)
    {
      int hit_plane = hit->WireID().Plane;
      double hit_integral = hit->Integral();

      aux_n_hits[hit_plane] ++;
      aux_charge[hit_plane] += hit_integral; 
    }
    for(size_t ii=0; ii<3; ii++)
    {
      double aux_dqdx = aux_charge[ii]/(aux_n_hits[ii]*pitch[ii]);
      charge[ii].push_back(aux_charge[ii]);
      dqdx[ii].push_back(aux_dqdx);
      n_hits[ii].push_back(aux_n_hits[ii]);
    }
  }

  // order spacepoints
  sortVectors(distance_projection, std::less<double>(), 
              distance_projection, x, y, z, 
              charge[0], charge[1], charge[2], 
              dqdx[0], dqdx[1], dqdx[2], 
              n_hits[0], n_hits[1], n_hits[2]);
  
  std::vector<double> trajectory_range, radial_distance, scattering_angle;
  double aux_previous_point[3];
  for(size_t i=0; i<distance_projection.size(); i++)
  {
    double this_point[3] = {x[i], y[i], z[i]};
    if (i == 0)
      trajectory_range.push_back(0);
    else
      aux_previous_point[0] = x[i-1];
      aux_previous_point[1] = y[i-1];
      aux_previous_point[2] = z[i-1];
      double aux_trajectory_distance = geoHelper.distance(this_point, aux_previous_point);
      trajectory_range.push_back(aux_trajectory_distance);

    double xyz[3] = {x[i] - vertex[0],
                     y[i] - vertex[1],
                     z[i] - vertex[2]};

    double aux_distance = geoHelper.distance(xyz, vertex);    
    double aux_radial_distance = sqrt(pow(aux_distance, 2) - pow(distance_projection[i], 2));
    radial_distance.push_back(aux_radial_distance);

    if (i != 0 && i != distance_projection.size())
    {
      double aux_next_point[3] = {x[i+1], y[i+1], z[i+1]};
      double aux_angle = geoHelper.scatteringAngle(aux_previous_point, this_point, aux_next_point);
      scattering_angle.push_back(aux_angle);
    }
  }
  double trajectory_length = trajectory_range.back();

  double mean_angle =  mean(scattering_angle);
  double rms_angle =  std_dev(scattering_angle);
  double mean_radial_distance =  mean(radial_distance);
  double rms_radial_distance =  std_dev(radial_distance);

  fDistance_projection = distance_projection;
  fX = x;
  fY = y;
  fZ = z;
  fChargeU = charge[0];
  fChargeV = charge[1];
  fChargeY = charge[2];
  fDqdx_U = dqdx[0];
  fDqdx_V = dqdx[1];
  fDqdx_Y = dqdx[2];
  fN_hits_U = n_hits[0];
  fN_hits_V = n_hits[1];
  fN_hits_Y = n_hits[2];
  fTrajectory_range = trajectory_range;
  fRadial_distance = radial_distance;
  fScattering_angle = scattering_angle;
  fTrajectory_length = trajectory_length;
  fMean_angle = mean_angle;
  fRms_angle = rms_angle;
  fMean_radial_distance = mean_radial_distance;
  fRms_radial_distance = rms_radial_distance;
}

void ChargeAnalyzer::reconfigure(fhicl::ParameterSet const &p)
{
  _pfp_producer                   = p.get<std::string>("PFParticleProducer", "pandoraNu");
  _hitfinderLabel                 = p.get<std::string>("HitProducer", "pandoraCosmicHitRemoval");
  _geantModuleLabel               = p.get<std::string>("GeantModule", "largeant");
  _spacepointLabel                = p.get<std::string>("SpacePointProducer", "pandoraNu");
  _mcpHitAssLabel                 = p.get<std::string>("MCPHitAssProducer", "crHitRemovalTruthMatch");
  _use_premade_ass                = p.get<bool>("UsePremadeMCPHitAss", true);
  _mctruthLabel                   = p.get<std::string>("MCTruthLabel", "generator");

  _dQdxRectangleWidth             = p.get<double>("dQdxRectangleWidth", 1);
  _dQdxRectangleLength            = p.get<double>("dQdxRectangleLength", 4);
}

void ChargeAnalyzer::analyze(art::Event const &evt)
{
  clear();
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();
  std::cout << "Event: " << fEvent << ", Run: " << fRun << ", subRun: " << fSubrun << std::endl;

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

    std::cout << "\n\nPfparticle:  " << i_pfp << ", pdgcode: " << fPdgCode << std::endl;

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
        continue;
      }
    }
    else
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
        continue;
      }
    }

    std::vector<double> fStart = {fStartx, fStarty, fStartz};
    fDistanceFromTrue = geoHelper.distance(fStart, fTrue_v_sce);
    // std::cout << "Startz: " << fStartz <<  ", Startx:  " << fStartx << std::endl;
    // clusters
    fNhits = 0;
    fNhitsU = 0;
    fNhitsV = 0;
    fNhitsY = 0;
    fNclusters = 0;

    std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(i_pfp);
    fNclusters = clusters.size();
    // int j = 0;
    for (art::Ptr<recob::Cluster> &cluster : clusters)
    {
      int aux_n_hits = cluster->NHits();
      fNhits += aux_n_hits;
      int fPlane = cluster->Plane().Plane;
      if (fPlane == 0)
      {
        fNhitsU += aux_n_hits;
      }
      else if (fPlane == 1)
      {
        fNhitsV += aux_n_hits;
      }
      else if (fPlane == 2)
      {
        fNhitsY += aux_n_hits;
      }
      if (cluster->Plane().Plane !=2)
        continue;
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

      // std::vector<art::Ptr<recob::Hit>> hits = hits_per_cluster.at(cluster.key());
      // for (art::Ptr<recob::Hit> &hit : hits)
      // {
      //   // std::cout << "hit on wire " << hit->WireID().Wire << " , " << 0.3*hit->WireID().Wire << " , " << hit->PeakTime() << std::endl;
      //   int fPlane = hit->WireID().Plane;
      //   if (fPlane == 0)
      //   {
      //     fNhitsU += 1;
      //   }
      //   else if (fPlane == 1)
      //   {
      //     fNhitsV += 1;
      //   }
      //   else if (fPlane == 2)
      //   {
      //     fNhitsY += 1;
      //   }
      // }
      // j++;
    }

    // dqdx
    fDQdx = {-1., -1., -1.};
    double collection_plane_cluster_length = fLength * 
      sqrt((pow(fDirectionx, 2) + pow(fDirectionz, 2))/(pow(fDirectionx, 2) + pow(fDirectiony, 2) + pow(fDirectionz, 2)));
    energyHelper.dQdx(i_pfp, 
                          evt, 
                          fDQdx, 
                          fDQdx_hits, 
                          fDQdx_wires, 
                          collection_plane_cluster_length, 
                          _dQdxRectangleWidth, 
                          _pfp_producer);
    fDQdx_U = fDQdx[0];
    fDQdx_V = fDQdx[1];
    fDQdx_Y = fDQdx[2];
    fn_hits_dQdx = fDQdx_hits.size();

    // storing energy
    std::vector<double> this_energy;
    std::vector<int> this_nhits;
    energyHelper.energyFromHits(*pfparticle, this_nhits, this_energy, evt, _pfp_producer);
    fReco_energy_U = this_energy[0];
    fReco_energy_V = this_energy[1];
    fReco_energy_Y = this_energy[2];

    recoTrueMatching(evt, pfparticle);

    double direction[3] = {fDirectionx, fDirectiony, fDirectionz};
    spacePointTrajectory(evt, i_pfp, direction);
    // // spacepoints loop
    // double x, y, z, x_prev, y_prev, z_prev;
    // int i = 0;

    // // double pitch[3];
    // // TVector3 direction = TVector3(fDirectionx, fDirectiony, fDirectionz);
    // // for(int j=0; j<3; j++)
    // // {
    // //   pitch[j] = geoHelper.getPitch(direction, j);
    // // }

    // std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(i_pfp);
    // for (art::Ptr<recob::SpacePoint> &sps : spcpnts)
    // {
    //   x_prev = x;
    //   y_prev = y;
    //   z_prev = z;
    //   auto xyz = sps->XYZ();
    //   x = xyz[0];
    //   y = xyz[1];
    //   z = xyz[2];

    //   fX.push_back(x);
    //   fY.push_back(y);
    //   fZ.push_back(z);

    //   if (i == 0)
    //   {
    //     fSpace_points_path.push_back(0);
    //   }
    //   else
    //   {
    //     double aux_d = std::sqrt((x-x_prev)*(x-x_prev) +
    //                              (y-y_prev)*(y-y_prev) +
    //                              (z-z_prev)*(z-z_prev));
    //     aux_d += fSpace_points_path.back();
    //     fSpace_points_path.push_back(aux_d);
    //   }

    //   std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(sps.key());
    //   double aux_ChargeU = 0;
    //   double aux_ChargeV = 0;
    //   double aux_ChargeY = 0;
    //   double aux_SigmaChargeU = 0;
    //   double aux_SigmaChargeV = 0;
    //   double aux_SigmaChargeY = 0;
    //   // double aux_ChargeU_corrected = 0;
    //   // double aux_ChargeV_corrected = 0;
    //   // double aux_ChargeY_corrected = 0;
    //   // double aux_SigmaChargeU_corrected = 0;
    //   // double aux_SigmaChargeV_corrected = 0;
    //   // double aux_SigmaChargeY_corrected = 0;
    //   for (art::Ptr<recob::Hit> &hit : hits)
    //   {
    //     double hit_plane = hit->WireID().Plane;
    //     double hit_integral = hit->Integral();
    //     double hit_sigma_integral = hit->SigmaIntegral();
    //     // double correction = GetChargeCorrection(hit_plane, x, y, z);
    //     if (hit_plane == 0)
    //     {
    //       aux_ChargeU += hit_integral;          
    //       // aux_ChargeU_corrected += aux_ChargeU*correction;
    //       aux_SigmaChargeU += hit_sigma_integral;
    //       // aux_SigmaChargeU_corrected += aux_SigmaChargeU*correction;
    //     }
    //     else if (hit_plane == 1)
    //     {
    //       aux_ChargeV += hit_integral;          
    //       // aux_ChargeV_corrected += aux_ChargeV*correction;
    //       aux_SigmaChargeV += hit_sigma_integral;
    //       // aux_SigmaChargeV_corrected += aux_SigmaChargeV*correction;
    //     }
    //     else if (hit_plane == 2)
    //     {
    //       aux_ChargeY += hit_integral;          
    //       // aux_ChargeY_corrected += aux_ChargeY*correction;
    //       aux_SigmaChargeY += hit_sigma_integral;
    //       // aux_SigmaChargeY_corrected += aux_SigmaChargeY*correction;
    //     }
    //     else
    //     {
    //       std::cout << "hit plane != 0, 1, 2, but " << hit_plane << std::endl;
    //     }
    //   }
    //   fChargeU.push_back(aux_ChargeU);
    //   fChargeV.push_back(aux_ChargeV);
    //   fChargeY.push_back(aux_ChargeY);
    //   fSigmaChargeU.push_back(aux_SigmaChargeU);
    //   fSigmaChargeV.push_back(aux_SigmaChargeV);
    //   fSigmaChargeY.push_back(aux_SigmaChargeY);

    //   // // fChargeU_corrected.push_back(aux_ChargeU_corrected);
    //   // // fChargeV_corrected.push_back(aux_ChargeV_corrected);
    //   // // fChargeY_corrected.push_back(aux_ChargeY_corrected);
    //   // // fSigmaChargeU_corrected.push_back(aux_SigmaChargeU_corrected);
    //   // // fSigmaChargeV_corrected.push_back(aux_SigmaChargeV_corrected);
    //   // // fSigmaChargeY_corrected.push_back(aux_SigmaChargeY_corrected);

    //   i++;
    // }

    // fSpLength = fSpace_points_path.back();



    fChargeTree->Fill();
    clear();
  }
}

DEFINE_ART_MODULE(ChargeAnalyzer)
