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
  std::vector<double> fResidual_range;
  std::vector<double> fRadial_distance;
  std::vector<double> fScattering_angle;
  double fRange;
  double fMean_angle;
  double fRms_angle;
  double fMean_radial_distance;
  double fRms_radial_distance;
  int fN_spacepoints;
  double fRange_length_ratio;
  double fSpacepoints_length_ratio;

  double fProjected_length_range_ratio;
  double fProjected_length_length_ratio;
  double fProjected_length;

  double fPitch_U, fPitch_V, fPitch_Y;
  double fDqdx_median_U, fDqdx_median_V, fDqdx_median_Y;
  double fDqdx_rms_U, fDqdx_rms_V, fDqdx_rms_Y;
  double fDqdx_median_first6_U, fDqdx_median_first6_V, fDqdx_median_first6_Y;
  double fDqdx_median_last8_U, fDqdx_median_last8_V, fDqdx_median_last8_Y;
  double fDqdx_median_last8_without_first5_U, fDqdx_median_last8_without_first5_V, fDqdx_median_last8_without_first5_Y;

  // reco-true matching information
  int fMatchedPdgCode;
  double fMatchedE;
  double fMatchedPx, fMatchedPy, fMatchedPz;
  double fMatchedVx, fMatchedVy, fMatchedVz;
  double fMatchedEndx, fMatchedEndy, fMatchedEndz;
  double fDistanceFromMatched, fMC_reco_costheta;

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
  fChargeTree->Branch("residual_range", "std::vector<double>", &fResidual_range);
  fChargeTree->Branch("radial_distance", "std::vector<double>", &fRadial_distance);
  fChargeTree->Branch("scattering_angle", "std::vector<double>", &fScattering_angle);

  fChargeTree->Branch("range", &fRange, "range/d");
  fChargeTree->Branch("mean_angle", &fMean_angle, "mean_angle/d");
  fChargeTree->Branch("rms_angle", &fRms_angle, "rms_angle/d");
  fChargeTree->Branch("mean_radial_distance", &fMean_radial_distance, "mean_radial_distance/d");
  fChargeTree->Branch("rms_radial_distance", &fRms_radial_distance, "rms_radial_distance/d");
  fChargeTree->Branch("projected_length", &fProjected_length, "projected_length/d");

  fChargeTree->Branch("n_space_points", &fN_spacepoints, "n_space_points/i");
  fChargeTree->Branch("range_length_ratio", &fRange_length_ratio, "range_length_ratio/d");
  fChargeTree->Branch("space_points_length_ratio", &fSpacepoints_length_ratio, "space_points_length_ratio/d");
  fChargeTree->Branch("projected_length_range_ratio", &fProjected_length_range_ratio, "projected_length_range_ratio/d");
  fChargeTree->Branch("projected_length_length_ratio", &fProjected_length_length_ratio, "projected_length_length_ratio/d");

  fChargeTree->Branch("pitch_U", &fPitch_U, "pitch_U/d");
  fChargeTree->Branch("pitch_V", &fPitch_V, "pitch_V/d");
  fChargeTree->Branch("pitch_Y", &fPitch_Y, "pitch_Y/d");

  fChargeTree->Branch("dqdx_median_U", &fDqdx_median_U, "dqdx_median_U/d");
  fChargeTree->Branch("dqdx_median_V", &fDqdx_median_V, "dqdx_median_V/d");
  fChargeTree->Branch("dqdx_median_Y", &fDqdx_median_Y, "dqdx_median_Y/d");

  fChargeTree->Branch("dqdx_rms_U", &fDqdx_rms_U, "dqdx_rms_U/d");
  fChargeTree->Branch("dqdx_rms_V", &fDqdx_rms_V, "dqdx_rms_V/d");
  fChargeTree->Branch("dqdx_rms_Y", &fDqdx_rms_Y, "dqdx_rms_Y/d");

  fChargeTree->Branch("dqdx_median_first6_U", &fDqdx_median_first6_U, "dqdx_median_first6_U/d");
  fChargeTree->Branch("dqdx_median_first6_V", &fDqdx_median_first6_V, "dqdx_median_first6_V/d");
  fChargeTree->Branch("dqdx_median_first6_Y", &fDqdx_median_first6_Y, "dqdx_median_first6_Y/d");

  fChargeTree->Branch("dqdx_median_last8_U", &fDqdx_median_last8_U, "dqdx_median_last8_U/d");
  fChargeTree->Branch("dqdx_median_last8_V", &fDqdx_median_last8_V, "dqdx_median_last8_V/d");
  fChargeTree->Branch("dqdx_median_last8_Y", &fDqdx_median_last8_Y, "dqdx_median_last8_Y/d");

  fChargeTree->Branch("dqdx_median_last8_without_first5_U", &fDqdx_median_last8_without_first5_U, "dqdx_median_last8_without_first5_U/d");
  fChargeTree->Branch("dqdx_median_last8_without_first5_V", &fDqdx_median_last8_without_first5_V, "dqdx_median_last8_without_first5_V/d");
  fChargeTree->Branch("dqdx_median_last8_without_first5_Y", &fDqdx_median_last8_without_first5_Y, "dqdx_median_last8_without_first5_Y/d");

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
  fChargeTree->Branch("mc_reco_costheta", &fMC_reco_costheta, "mc_reco_costheta/d");

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
  fResidual_range.clear();
  fRadial_distance.clear();
  fScattering_angle.clear();

  fDQdx_hits.clear();
  fDQdx_wires.clear();
  fDQdx.clear();
}

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

  if (generator.size() > 0)
  {
    for (auto &gen : generator)
    {
      if (gen.Origin() == simb::kBeamNeutrino)
      {
        // there_is_a_neutrino = true;
        fTrue_nu_pdg = gen.GetNeutrino().Nu().PdgCode();
        fTrue_nu_energy = gen.GetNeutrino().Nu().E();
        fTrue_ccnc = gen.GetNeutrino().CCNC();
        fTrue_theta = gen.GetNeutrino().Theta();

        fTrue_vx = gen.GetNeutrino().Nu().Vx();
        fTrue_vy = gen.GetNeutrino().Nu().Vy();
        fTrue_vz = gen.GetNeutrino().Nu().Vz();

        fTrue_v = {fTrue_vx, fTrue_vy, fTrue_vz};

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
      }
    }
  }
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
  for (size_t ii = 0; ii < 3; ii++)
  {
    pitch[ii] = geoHelper.getPitch(t_dir, ii);
  }

  std::vector<double> distance_projection, x, y, z;
  std::vector<double> charge[3], dqdx[3];
  std::vector<int> n_hits[3];
  std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(i_pfp);

  for (art::Ptr<recob::SpacePoint> &sps : spcpnts)
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
    for (art::Ptr<recob::Hit> &hit : hits)
    {
      int hit_plane = hit->WireID().Plane;
      double hit_integral = hit->Integral();

      aux_n_hits[hit_plane]++;
      aux_charge[hit_plane] += energyHelper.integral2charge(hit_integral, hit_plane);
    }
    for (size_t ii = 0; ii < 3; ii++)
    {
      double aux_dqdx = aux_charge[ii] / (aux_n_hits[ii] * pitch[ii]);
      aux_dqdx = energyHelper.dQdx2dEdx(aux_dqdx);
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
  for (size_t i = 0; i < distance_projection.size(); i++)
  {
    double this_point[3] = {x[i], y[i], z[i]};
    if (i == 0)
    {
      trajectory_range.push_back(0);
    }
    else
    {
      aux_previous_point[0] = x[i - 1];
      aux_previous_point[1] = y[i - 1];
      aux_previous_point[2] = z[i - 1];
      double aux_trajectory_distance = geoHelper.distance(this_point, aux_previous_point);
      aux_trajectory_distance += trajectory_range.back();
      trajectory_range.push_back(aux_trajectory_distance);
    }

    double xyz[3] = {x[i], y[i], z[i]};
    double aux_distance = geoHelper.distance(xyz, vertex);
    double aux_radial_distance = sqrt(pow(aux_distance, 2) - pow(distance_projection[i], 2));
    radial_distance.push_back(aux_radial_distance);

    if (i != 0 && i != (distance_projection.size() - 1))
    {
      double aux_next_point[3] = {x[i + 1], y[i + 1], z[i + 1]};
      double aux_angle = geoHelper.scatteringAngle(aux_previous_point, this_point, aux_next_point);
      scattering_angle.push_back(aux_angle);
    }
  }

  double range = trajectory_range.back();
  std::vector<double> residual_range;
  std::vector<double> dqdx_for_median_first6[3], dqdx_for_median_last8[3], dqdx_for_median_last8_without_first5[3];
  for (size_t i = 0; i < trajectory_range.size(); i++)
  {
    residual_range.push_back(range - trajectory_range[i]);
    if (trajectory_range[i] < 6)
    {
      dqdx_for_median_first6[0].push_back(dqdx[0][i]);
      dqdx_for_median_first6[1].push_back(dqdx[1][i]);
      dqdx_for_median_first6[2].push_back(dqdx[2][i]);
    }
    if (residual_range[i] < 8)
    {
      dqdx_for_median_last8[0].push_back(dqdx[0][i]);
      dqdx_for_median_last8[1].push_back(dqdx[1][i]);
      dqdx_for_median_last8[2].push_back(dqdx[2][i]);
      if (trajectory_range[i] > 5)
      {
        dqdx_for_median_last8_without_first5[0].push_back(dqdx[0][i]);
        dqdx_for_median_last8_without_first5[1].push_back(dqdx[1][i]);
        dqdx_for_median_last8_without_first5[2].push_back(dqdx[2][i]);
      }
    }
  }

  double mean_angle = mean(scattering_angle);
  double rms_angle = std_dev(scattering_angle);
  double mean_radial_distance = mean(radial_distance);
  double rms_radial_distance = std_dev(radial_distance);

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
  fResidual_range = residual_range;
  fRadial_distance = radial_distance;
  fScattering_angle = scattering_angle;
  fRange = range;
  fMean_angle = mean_angle;
  fRms_angle = rms_angle;
  fMean_radial_distance = mean_radial_distance;
  fRms_radial_distance = rms_radial_distance;
  fN_spacepoints = x.size();

  fPitch_U = pitch[0],
  fPitch_V = pitch[1],
  fPitch_Y = pitch[2];
  fDqdx_median_U = median(fDqdx_U);
  fDqdx_median_V = median(fDqdx_V);
  fDqdx_median_Y = median(fDqdx_Y);
  fDqdx_rms_U = std_dev(fDqdx_U);
  fDqdx_rms_V = std_dev(fDqdx_V);
  fDqdx_rms_Y = std_dev(fDqdx_Y);
  fDqdx_median_first6_U = median(dqdx_for_median_first6[0]);
  fDqdx_median_first6_V = median(dqdx_for_median_first6[1]);
  fDqdx_median_first6_Y = median(dqdx_for_median_first6[2]);
  fDqdx_median_last8_U = median(dqdx_for_median_last8[0]);
  fDqdx_median_last8_V = median(dqdx_for_median_last8[1]);
  fDqdx_median_last8_Y = median(dqdx_for_median_last8[2]);
  fDqdx_median_last8_without_first5_U = median(dqdx_for_median_last8_without_first5[0]);
  fDqdx_median_last8_without_first5_V = median(dqdx_for_median_last8_without_first5[1]);
  fDqdx_median_last8_without_first5_Y = median(dqdx_for_median_last8_without_first5[2]);

  fProjected_length = distance_projection.back() - distance_projection.front();
}

void ChargeAnalyzer::reconfigure(fhicl::ParameterSet const &p)
{
  _pfp_producer = p.get<std::string>("PFParticleProducer", "pandoraNu");
  _hitfinderLabel = p.get<std::string>("HitProducer", "pandoraCosmicHitRemoval");
  _geantModuleLabel = p.get<std::string>("GeantModule", "largeant");
  _spacepointLabel = p.get<std::string>("SpacePointProducer", "pandoraNu");
  _mcpHitAssLabel = p.get<std::string>("MCPHitAssProducer", "crHitRemovalTruthMatch");
  _use_premade_ass = p.get<bool>("UsePremadeMCPHitAss", true);
  _mctruthLabel = p.get<std::string>("MCTruthLabel", "generator");

  _dQdxRectangleWidth = p.get<double>("dQdxRectangleWidth", 1);
  _dQdxRectangleLength = p.get<double>("dQdxRectangleLength", 4);
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

  trueNeutrinoInformation(evt);
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
        fLength = pf_objs[0]->Length();
        fDirectionx = pf_objs[0]->Direction().X();
        fDirectiony = pf_objs[0]->Direction().Y();
        fDirectionz = pf_objs[0]->Direction().Z();
        fEndx = fStartx + fDirectionx * fLength;
        fEndy = fStarty + fDirectiony * fLength;
        fEndz = fStartz + fDirectionz * fLength;
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

    recoTrueMatching(evt, pfparticle);
    double reco_dir[3] = {fDirectionx, fDirectiony, fDirectionz};
    double mc_dir[3] = {fMatchedPx, fMatchedPy, fMatchedPz};
    fMC_reco_costheta = geoHelper.costheta(reco_dir, mc_dir);

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
      if (cluster->Plane().Plane != 2)
        continue;
    }

    // dqdx
    fDQdx = {-1., -1., -1.};
    double collection_plane_cluster_length = fLength *
                                             sqrt((pow(fDirectionx, 2) + pow(fDirectionz, 2)) / (pow(fDirectionx, 2) + pow(fDirectiony, 2) + pow(fDirectionz, 2)));
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

    spacePointTrajectory(evt, i_pfp, reco_dir);
    fRange_length_ratio = fRange / fLength;
    fSpacepoints_length_ratio = fN_spacepoints / fLength;
    fProjected_length_range_ratio = fProjected_length / fRange;
    fProjected_length_length_ratio = fProjected_length / fLength;

    fChargeTree->Fill();
    clear();
  }
}

DEFINE_ART_MODULE(ChargeAnalyzer)
