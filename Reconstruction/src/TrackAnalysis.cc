/*
TrackAnalysis.cc - Implementation of TrackAnalysis class for muon track reconstruction.
Units in Input Data:
Positions are in mm (as you mentioned)
Momenta are in MeV/c (as you mentioned)
Units in GENFIT:
GENFIT generally expects positions in cm (not mm)
GENFIT expects momenta in GeV/c (not MeV/c)
Magnetic field strength should be in Tesla
*/

#include "TrackAnalysis.hh"
#include <RKTrackRep.h>
#include <FieldManager.h>
#include <DAF.h>
#include <PlanarMeasurement.h>
#include <SpacepointMeasurement.h>

#include <StateOnPlane.h>
#include <FitStatus.h> 
#include <TMath.h>
#include <iostream>

#include "TDecompChol.h" 
#include "TDecompLU.h"
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TFile.h>
#include <TTree.h>
#include <map>
#include <cmath>
#include <tuple>

#include <TCanvas.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TPolyLine3D.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLine.h>  
#include <TMultiGraph.h>

#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
#include <TGeoShape.h>
#include <algorithm>
#include <functional>
#include <memory> // for unique_ptr
#include <iomanip>
#include <numeric>

using namespace genfit;
///////////////////////////////////////////////////////// 
///// TrackAnalysis class implementation.         ///////
///////////////////////////////////////////////////////// 
TrackAnalysis::TrackAnalysis(TGeoManager* geo) : trackRepMinus_(nullptr), trackRepPlus_(nullptr), magneticField_(nullptr), geo_(geo) {
    // Initialize GENFIT track representation for muons
    trackRepMinus_ = new genfit::RKTrackRep(13);   // mu-
    trackRepPlus_  = new genfit::RKTrackRep(-13);  // mu+
    
    // Set up magnetic field
    if (geo_) {
        geo_ = geo;
        BuildGeometryMaps();
        // Use geometry-based field
        magneticField_ = new GeoField(this);
        genfit::FieldManager::getInstance()->init(magneticField_);
    } else {
        std::cerr << "Warning: No geometry provided, using default uniform magnetic field." << std::endl;
        // Default to uniform field as fallback
        magneticField_ = new MagneticField(this);
        genfit::FieldManager::getInstance()->init(magneticField_);
    }
}
///////////////////////////////////////////////////////// 
//// Destructor for TrackAnalysis class.          ///////
/////////////////////////////////////////////////////////
TrackAnalysis::~TrackAnalysis() {
    delete trackRepMinus_;
    delete trackRepPlus_;
    delete magneticField_;
}
///////////////////////////////////////////////////////// 
//// BuildGeometryMaps method implementation.     ///////
/////////////////////////////////////////////////////////
void TrackAnalysis::BuildGeometryMaps() {
    magnet_z_ranges_.clear();
    scifi_layer_z_.clear();

    auto nodes = geo_->GetTopVolume()->GetNodes();
    for (int i = 0; i < nodes->GetEntries(); ++i) {
        TGeoNode* node = (TGeoNode*)nodes->At(i);
        std::string name = node->GetName();
        TGeoBBox* shape = dynamic_cast<TGeoBBox*>(node->GetVolume()->GetShape());
        if (!shape) continue;

        double zc = node->GetMatrix()->GetTranslation()[2];
        double dz = shape->GetDZ();

        if (name.find("Magnet") != std::string::npos) {
            magnet_z_ranges_.push_back({zc - dz, zc + dz});
        }
        if (name.find("SciFiLayer") != std::string::npos || name.find("InitialSciFiLayer") != std::string::npos) {
            scifi_layer_z_.push_back(zc);
        }
    }
    // Sort the ranges and layer positions
    std::sort(magnet_z_ranges_.begin(), magnet_z_ranges_.end());
    std::sort(scifi_layer_z_.begin(), scifi_layer_z_.end());
    // For debugging:
    if(verbose >=2){
        std::cout << "[Geom] Magnet Z-ranges:" << std::endl;
        for (auto& r : magnet_z_ranges_) 
            std::cout << "   " << r.first << " to " << r.second << " mm" << std::endl;
        std::cout << "[Geom] SciFiLayer Z positions: ";
        for (auto& z : scifi_layer_z_) 
            std::cout << z << " ";
        std::cout << std::endl;
    }
}
/////////////////////////////////////////////////////////
//// IsInMagnet method implementation.            ///////
/////////////////////////////////////////////////////////
bool TrackAnalysis::IsInMagnet(double z) const {
    z = z / 10.0; // Convert to cm
    for (const auto& rng : magnet_z_ranges_)
        if (z > rng.first && z < rng.second)
            return true;
    return false;
}
/////////////////////////////////////////////////////////
//// isMagnetExist method implementation.         ///////
/////////////////////////////////////////////////////////
bool TrackAnalysis::IsMagnetExist(double z1,double z2) const { 
    double z_low = std::min(z1, z2)/10.; // convert to cm
    double z_high = std::max(z1, z2)/10.; // convert to cm  
    bool segment_crosses_magnet = false;
    
    for (const auto& rng : magnet_z_ranges_)
        {
            if(verbose>= 3)
                std::cout << "[DEBUG] Checking magnet range: " << rng.first << " to " << rng.second 
                          << " against segment [" << z_low << ", " << z_high << "]" << std::endl;
            if (rng.first<z_high && rng.second >z_low)
                return true;
        }
    return false;
}
/////////////////////////////////////////////////////////
//// GetFieldX method implementation.            ///////
/////////////////////////////////////////////////////////
double TrackAnalysis::GetFieldX(double x, double y, double z) const {
   if (verbose >= 3) {
        std::cout << "[GetFieldX] Input: x=" << x << ", y=" << y << ", z=" << z << std::endl;
    }    
    // x, y and z are in mm
    // magnet_z_ranges_ are in cm
    z = z / 10.0; // Convert to cm
    y = y / 10.0; // Convert to cm
    // Check if the position is within any magnet range
    for (const auto& rng : magnet_z_ranges_) {
        if (z > rng.first && z < rng.second) {
            if (std::abs(y) >= 25. && std::abs(y) <= 50.) return MagnetFieldStrength; // Tesla
            else if (std::abs(y) < 25.) return -1*MagnetFieldStrength; // Tesla
            else return 0.0; // Outside the magnet field    
        }
    }
    return 0.0;
}
double TrackAnalysis::GetFieldBX(double y) const {
    y = y / 10.0; // Convert to cm
    // Check if the position is within any magnet range
    if (std::abs(y) >= 25. && std::abs(y) <= 50.) return MagnetFieldStrength; // Tesla
    else if (std::abs(y) < 25.) return -1*MagnetFieldStrength; // Tesla
    else return 0.0; // Outside the magnet field        
    return 0.0;
}
/////////////////////////////////////////////////////////
//// GeoField class implementation.               ///////
/////////////////////////////////////////////////////////
/*
TVector3 TrackAnalysis::GeoField::get(const TVector3& position) const {
    // GENFIT gives position in cm! 
    double x_mm = position.X() * 10.0;
    double y_mm = position.Y() * 10.0;
    double z_mm = position.Z() * 10.0;
    double bx = ta_->GetFieldX(x_mm, y_mm, z_mm); // in Tesla
    if (ta_->verbose >= 3)
        std::cout << "[GeoField] Query at (x=" << x_mm << ", y=" << y_mm << ", z=" << z_mm
                  << ") mm --> Bx=" << bx << " T" << std::endl;
    // Only Bx, By = 0, Bz = 0
    return TVector3(bx, 0, 0);
}*/
/////////////////////////////////////////////////////////
//// Get IronThickness                            ///////
/////////////////////////////////////////////////////////
double TrackAnalysis::GetIronThickness() const {
    // the iron volume is named "Magnet" in the geometry
    TGeoVolume* ironVol = geo_ ? geo_->FindVolumeFast("Magnet") : nullptr;
    if (ironVol) {
        TGeoShape* shape = ironVol->GetShape();
        // box: TGeoBBox* box = dynamic_cast<TGeoBBox*>(shape);
        if (auto box = dynamic_cast<TGeoBBox*>(shape)) {
            // Iron thickness is half-length in Z times 2 (full length)
            return box->GetDZ() * 2.0 / 10.0;  // [cm] to [m]
        }
    }
    return IronThickness; // 15 cm default
}
///////////////////////////////////////////////////////// 
/////   ProcessFile method implementation.        ///////
/////////////////////////////////////////////////////////
void TrackAnalysis::ProcessFile(const char* filename) {
    TFile file(filename);
    TTree* tree = (TTree*)file.Get("Hits");
    if (!tree) {
        std::cerr << "Error: Could not find tree 'Hits'" << std::endl;
        return;
    }

    // Branch variables
    int eventID;
    std::vector<double> *x=0, *y=0, *z=0, *px=0, *py=0, *pz=0;
    std::vector<int> *pdg=0, *trackID=0, *layerID=0, *stationID=0;

    tree->SetBranchAddress("eventID",    &eventID);
    tree->SetBranchAddress("x",          &x);
    tree->SetBranchAddress("y",          &y);
    tree->SetBranchAddress("z",          &z);
    tree->SetBranchAddress("px",         &px);
    tree->SetBranchAddress("py",         &py);
    tree->SetBranchAddress("pz",         &pz);
    tree->SetBranchAddress("pdg",        &pdg);
    tree->SetBranchAddress("trackID",    &trackID);
    tree->SetBranchAddress("layerID",    &layerID);
    tree->SetBranchAddress("stationID",  &stationID);
    int nEvents = tree->GetEntries();
    for (int i = 0; i < nEvents; ++i) {
        tree->GetEntry(i);
        ProcessEvent(eventID, *trackID, *x, *y, *z, *px, *py, *pz, *layerID, *stationID, *pdg);
    }
    return;
}
///////////////////////////////////////////////////////// 
///// ProcessEvent method implementation.         ///////
/////////////////////////////////////////////////////////

void TrackAnalysis::ProcessEvent(
    int eventID, 
    const std::vector<int>& trackID,
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& z,
    const std::vector<double>& px,
    const std::vector<double>& py,
    const std::vector<double>& pz,
    const std::vector<int>& layerID,
    const std::vector<int>& stationID,
    const std::vector<int>& pdg)
{
    //if (eventID > 1 )
    //    return; // Skip events that are not the first one

    // 1. Convert raw input to hits
    std::vector<Hit> hits;
    for (size_t i = 0; i < x.size(); ++i) {
        hits.push_back({
            .stationID = stationID[i],
            .layerID   = layerID[i],
            .trackID   = trackID[i],
            .pdg       = pdg[i],
            .position  = TVector3(x[i], y[i], z[i]),
            .momentum  = TVector3(px[i], py[i], pz[i])
        });
    }
    std::cout << "############################################" << std::endl;
    std::cout << "[DEBUG] Processing event " << eventID 
              << " with " << hits.size() << " hits." << std::endl; 
    std::cout << "############################################" << std::endl;

    // 2. Group hits by (trackID, pdg)
    std::map<std::pair<int,int>, std::vector<Hit>> hits_by_track;
    for (const auto& hit : hits) {
        hits_by_track[{hit.trackID, hit.pdg}].push_back(hit);
    }

    // 3. Loop over tracks
    double iron_thickness = GetIronThickness(); // in m
    //double B_field = -1.8; // Tesla

    for (const auto& entry : hits_by_track) {
        // -- Extract trackID, pdg and hits
        // entry.first is a pair of (trackID, pdg)
        // entry.second is a vector of Hits for that track
        int this_trackID = entry.first.first;
        int this_pdg     = entry.first.second;
        const std::vector<Hit>& track_hits = entry.second;
        if (track_hits.size() < 4) continue;
        results_.emplace_back();
        results_.back().eventID = eventID;
        results_.back().trackID = this_trackID;
        results_.back().pdg = this_pdg;
        //Initialize all charge fields to 0
        results_.back().sagittaCharge = 0;
        results_.back().kasaCharge = 0;
        results_.back().taubinCharge = 0;
        results_.back().trackletDeflectionCharge = 0;
        results_.back().kalmanCharge = 0;
        results_.back().genfitCharge = 0;
        //Initialize resolution fields
        results_.back().sagittaResolution = -1;
        results_.back().kasaResolution = -1;
        results_.back().taubinResolution = -1;
        results_.back().trackletDeflectionMomentumErr = -1;       
        results_.back().sagittaResolution_trk = -1;
        results_.back().kasaResolution_trk = -1;
        results_.back().taubinResolution_trk = -1;
        results_.back().trackletDeflectionMomentumErr_trk = -1; // Initialize tracklet deflection momentum error
        
        // -- Sort hits by Z
        std::vector<Hit> sorted_hits = track_hits;
        std::sort(sorted_hits.begin(), sorted_hits.end(),
                  [](const Hit& a, const Hit& b) { return a.position.Z() < b.position.Z(); });

        // -- Truth momentum from MC using the first hit
        results_.back().truthMomentum = sorted_hits.front().momentum.Mag() * MeV_to_GeV;
        // -- Print debug information
        std::cout << "[DEBUG] TrackID: " << this_trackID 
                  << ", PDG: " << this_pdg 
                  << ", Truth Momentum: " << results_.back().truthMomentum << " GeV/c" << std::endl;

        if(this_pdg != 13 && this_pdg != -13) {
            std::cerr << "[SKIP] It is not muon track: " << this_pdg << " for trackID=" << this_trackID << std::endl;
            continue; // Skip unsupported PDG codes
        }
        ///////////////////////////////////////////////////////
        // using only hits
        // -- Calculate sagitta momentum using hits   
        double sagittaMomentum = CalculateSagittaMomentum(sorted_hits);
        results_.back().sagittaMomentum = sagittaMomentum;
        //if (verbose == 1) 
            std::cout << "[SAGITTA] TrackID: " << this_trackID
                      << ", PDG: " << this_pdg
                      << ", Sagitta Momentum: " << sagittaMomentum << " GeV/c" << std::endl;
        std::cout << "#####%%%%%%#######" << std::endl;
        // -- 2D Circle fit using Kasa method
        auto kasaResult = CircleFitYZ_Kasa(sorted_hits, MagnetFieldStrength);
        results_.back().kasaMomentum = kasaResult.p;
        results_.back().kasaCharge = kasaResult.charge;
        results_.back().kasaResolution = kasaResult.sigma; // Store sigma as resolution
        //if (verbose == 1) 
            std::cout << "[CIRCLE_Kasa_Hits] TrackID: " << this_trackID 
                      << ", PDG: " << this_pdg 
                      << ", Circle Fit Center: (" << kasaResult.yc << ", " << kasaResult.zc<< ") mm" 
                      << ", Circle Fit Radius: " << kasaResult.R << ", Sigma: " << kasaResult.sigma
                      << ", Circle Fit Momentum: " << kasaResult.p << " GeV/c ± " << kasaResult.p_err  << " GeV/c" 
                      << ", circleCharge: " << kasaResult.charge
                      << std::endl;
        std::cout << "#####%%%%%%#######" << std::endl;
        // -- 2D Circle fit using Taubin method
        auto taubinResult = CircleFitYZ_Taubin(sorted_hits, MagnetFieldStrength);
        results_.back().taubinMomentum = taubinResult.p;
        results_.back().taubinCharge = taubinResult.charge;
        results_.back().taubinResolution = taubinResult.sigma; // Store sigma as resolution
        //if (verbose == 1) 
            std::cout << "[CIRCLE_Taubin_Hits] TrackID: " << this_trackID 
                      << ", PDG: " << this_pdg 
                      << ", Circle Fit Center: (" << taubinResult.yc << ", " << taubinResult.zc << ") mm" 
                      << ", Circle Fit Radius: " << taubinResult.R << ", Sigma: " << taubinResult.sigma 
                      << ", Circle Fit Momentum: " << taubinResult.p << " GeV/c ± " << taubinResult.p_err  << " GeV/c" 
                      << ", circleCharge: " << taubinResult.charge
                      << std::endl;
        std::cout << "#####%%%%%%#######" << std::endl;
        ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////
        // -- constructing tracklets from hits
        ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////
        // Fit tracklets from hits
        auto tracklets = FitStationTracklets(sorted_hits);
        if (tracklets.size() < 2) {
            std::cerr << "[ERROR] Not enough tracklets for Kalman/Genfit on trackID=" << this_trackID << std::endl;
            continue;
        }

        // -- Calculate sagitta momentum using hits   
        double sagittaMomentum_trk = CalculateSagittaMomentum(tracklets);
        results_.back().sagittaMomentum_trk = sagittaMomentum_trk;
        //if (verbose == 1) 
            std::cout << "[SAGITTA] TrackID: " << this_trackID
                      << ", PDG: " << this_pdg
                      << ", Sagitta Momentum: " << sagittaMomentum_trk  << " GeV/c" << std::endl;
        std::cout << "#####%%%%%%#######" << std::endl;
        

        auto kasaResult_trk = CircleFitYZ_Kasa(tracklets, MagnetFieldStrength);
        results_.back().kasaMomentum_trk = kasaResult_trk.p;
        results_.back().kasaCharge_trk = kasaResult_trk.charge;
        results_.back().kasaResolution_trk = kasaResult_trk.p_err;
        //if (verbose == 1) 
            std::cout << "[CIRCLE_Kasa_Tracklets] TrackID: " << this_trackID 
                      << ", PDG: " << this_pdg 
                      << ", Circle Fit Center: (" << kasaResult_trk.yc << ", " << kasaResult_trk.zc << ") mm" 
                      << ", Circle Fit Radius: " << kasaResult_trk.R << ", Sigma: " << kasaResult_trk.sigma
                      << ", Circle Fit Momentum: " << kasaResult_trk.p << " GeV/c ± " << kasaResult_trk.p_err  << " GeV/c" 
                      << std::endl;
        std::cout << "#####%%%%%%#######" << std::endl;
        // -- 2D Circle fit using Taubin method
        auto taubinResult_trk = CircleFitYZ_Taubin(tracklets, MagnetFieldStrength);
        results_.back().taubinMomentum_trk = taubinResult_trk.p;
        results_.back().taubinCharge_trk = taubinResult_trk.charge;
        results_.back().taubinResolution_trk = taubinResult_trk.p_err;
        //if (verbose == 1) 
            std::cout << "[CIRCLE_Taubin_Tracklets] TrackID: " << this_trackID 
                      << ", PDG: " << this_pdg 
                      << ", Circle Fit Center: (" << taubinResult_trk.yc << ", " << taubinResult_trk.zc << ") mm" 
                      << ", Circle Fit Radius: " << taubinResult_trk.R << ", Sigma: " << taubinResult_trk.sigma
                      << ", Circle Fit Momentum: " << taubinResult_trk.p << " GeV/c ± " << taubinResult_trk.p_err  << " GeV/c" 
                      << std::endl;
        std::cout << "#####%%%%%%#######" << std::endl;
        //use the tracklet deflection method in 3D to get momentum
        auto reco = MomentumFromTrackletDeflection(tracklets, MagnetFieldStrength, iron_thickness, true);
        //double trackletMomentum = reco.first;
        //double trackletMomentumErr = reco.second;
        double trackletMomentum = std::get<0>(reco);
        double trackletMomentumErr = std::get<1>(reco);
        int trackletCharge = std::get<2>(reco);
        results_.back().trackletDeflectionMomentum = trackletMomentum;
        results_.back().trackletDeflectionMomentumErr = trackletMomentumErr;
        results_.back().trackletDeflectionCharge = trackletCharge;
        // Print the tracklet deflection momentum
        std::cout << "[DEFLECTION] Tracklet Deflection Momentum: " 
                  << trackletMomentum << " GeV/c, error = " << trackletMomentumErr 
                  << ", charge = " << trackletCharge
                  << std::endl;        
        std::cout << "#####%%%%%%#######" << std::endl;

        double effective_iron_thickness = iron_thickness * (tracklets.size() - 1);
        results_.back().trackletDeflectionMomentum_trk = GlobalDirectionDeflectionMomentum(tracklets, MagnetFieldStrength, effective_iron_thickness);
        results_.back().trackletDeflectionMomentumErr_trk = trackletMomentumErr; // Use the same error as for the deflection method
        results_.back().trackletDeflectionCharge_trk = trackletCharge;
        std::cout << "[GLOBAL] Global Direction Deflection Momentum: " 
                  << results_.back().trackletDeflectionMomentum_trk << " GeV/c" 
                  << ", error = " << results_.back().trackletDeflectionMomentumErr_trk 
                  << ", charge = " << results_.back().trackletDeflectionCharge_trk
                  << std::endl;
        std::cout << "#####%%%%%%#######" << std::endl;
        
        //auto smoothed_positions = SmoothedTrackletPositions(tracklets);
        //std::vector<Tracklet> smoothed_tracklets = tracklets; // Make a copy
        //for (size_t i = 0; i < tracklets.size(); ++i)
        //    smoothed_tracklets[i].position = smoothed_positions[i];
        //auto res = CircleFitYZ_Taubin(smoothed_tracklets);
        //double p = 0.3 * 1.8 * res.R; // p in GeV/c  
        //std::cout << "[CIRCLE_Taubin Tracklets] TrackID: " << this_trackID 
        //          << ", PDG: " << this_pdg 
        //          << ", Circle Fit Momentum: " << p << " GeV/c" << std::endl;   
        //bool collinear = true;
        //if (smoothed_positions.size() >= 3) {
        //    TVector3 v1 = smoothed_positions[1] - smoothed_positions[0];
        //    TVector3 v2 = smoothed_positions[2] - smoothed_positions[1];
        //    double angle = v1.Angle(v2);
        //    if (angle > 1e-3) collinear = false;
       // }
        //if (!collinear) {
        //    std::vector<Tracklet> smoothed_tracklets = tracklets;
        //    for (size_t i = 0; i < tracklets.size(); ++i)
        //        smoothed_tracklets[i].position = smoothed_positions[i];
        //    auto res = CircleFitYZ_Taubin(smoothed_tracklets);
        //    double p = 0.3 * MagnetFieldStrength * res.R;
        //    std::cout << "[CIRCLE_Taubin Tracklets] TrackID: " << this_trackID 
        //            << ", PDG: " << this_pdg 
        //            << ", Circle Fit Momentum: " << p << " GeV/c" << std::endl;
        //} else {
        //    std::cout << "[CIRCLE_Taubin Tracklets] Track is collinear after smoothing, circle fit not reliable." << std::endl;
       // }
        

       
        // -- Try both seeds for Kalman/Genfit fits
        struct SeedOption { double value; std::string name; };
        std::vector<SeedOption> seeds = {
            {kasaResult.p, "kasaFit"},
            {taubinResult.p, "tabuinFit"},
        };

        for (const auto& seed : seeds) {
            std::cout << "[DEBUG] Kalman Seed (" << seed.name << "): " << seed.value << std::endl;
            if (std::isnan(seed.value) || seed.value < 0.1) {
                std::cerr << "[ERROR] Invalid seed value (" << seed.value << "), skipping Kalman for this seed!" << std::endl;
                continue;
            }

            std::cout << "[DEBUG] Running Kalman on " << tracklets.size() << " tracklets for trackID " << this_trackID << std::endl;
            if(verbose>=3)
                for (const auto& t : tracklets) {
                    std::cout << "   Tracklet: z=" << t.position.Z() << ", station=" << t.station << std::endl;
                }

            // ----- Run Kalman filter
            RunKalmanFilterTracklets(tracklets, seed.value);
            Result& res = results_.back();
            std::cout << "[KALMANFILTER] results for trackID " << this_trackID 
                      << ": Momentum = " << res.kalmanMomentum 
                      << " GeV/c, Chi2 = " << res.kalmanChi2 
                      << " GeV/c, Error = " << res.kalmanMomentumErr << std::endl;
            std::cout << "#####%%%%%%#######" << std::endl;

            // ----- Run Genfit
            RunGenfitTracklets(tracklets, seed.value);
            std::cout << "[GENFIT] results for trackID " << this_trackID 
                      << ": Momentum = " << res.genfitMomentum 
                      << " GeV/c, Chi2 = " << res.genfitChi2 
                      << " GeV/c, Error = " << res.genfitMomentumErr << std::endl;
            std::cout << "#####%%%%%%#######" << std::endl;
        
        }
    }


}

/////////////////////////////////////////////////////////
// SmoothedTrackletPositions method implementation.  //// 
/////////////////////////////////////////////////////////
std::vector<TVector3> TrackAnalysis::SmoothedTrackletPositions(const std::vector<Tracklet>& tracklets) {
    // Fit a line in Z-X and Z-Y
    std::vector<TVector3> positions;
    for (const auto& t : tracklets) positions.push_back(t.position);
    double slope_xz = 0, intcpt_xz = 0, slope_yz = 0, intcpt_yz = 0;
    LinearFit(positions, slope_xz, intcpt_xz, 'x');
    LinearFit(positions, slope_yz, intcpt_yz, 'y');

    std::vector<TVector3> smoothed;
    for (const auto& t : tracklets) {
        double z = t.position.Z();
        double x = slope_xz * z + intcpt_xz;
        double y = slope_yz * z + intcpt_yz;
        smoothed.emplace_back(x, y, z);
    }
    return smoothed;
}
/////////////////////////////////////////////////////////////
// Plot CircleFit YZ together with hits implementation.  ////
/////////////////////////////////////////////////////////////
void TrackAnalysis::PlotCircleFitYZ(const std::vector<Hit>& hits, double yc, double zc, double R, const std::string& title) {
    TCanvas* c = new TCanvas(("c_"+title).c_str(), title.c_str(), 800, 700);
    TGraph* gHits = new TGraph(hits.size());
    for (size_t i = 0; i < hits.size(); ++i)
        gHits->SetPoint(i, hits[i].position.Z(), hits[i].position.Y());
    gHits->SetMarkerStyle(20);
    gHits->SetMarkerColor(kRed);
    gHits->SetTitle((title+";Z [mm];Y [mm]").c_str());
    gHits->Draw("AP");

    // Draw fitted circle
    const int N = 200;
    double theta0 = 0, theta1 = 2*TMath::Pi();
    TGraph* gCircle = new TGraph(N);
    for (int i = 0; i < N; ++i) {
        double theta = theta0 + (theta1-theta0)*i/(N-1);
        double y = yc + R * std::sin(theta);
        double z = zc + R * std::cos(theta);
        gCircle->SetPoint(i, z, y);
    }
    gCircle->SetLineColor(kBlue);
    gCircle->SetLineWidth(2);
    gCircle->Draw("L SAME");

    TLegend* leg = new TLegend(0.15,0.8,0.45,0.93);
    leg->AddEntry(gHits, "Hits (Y,Z)", "p");
    leg->AddEntry(gCircle, "Fitted circle", "l");
    leg->Draw();

    c->SetGrid();
    c->Modified(); c->Update();
    c->SaveAs((title+"_CircleFitYZ.pdf").c_str());
}
/////////////////////////////////////////////////////////
// Linear interpolation between SciFi hits for a given Z/
/////////////////////////////////////////////////////////
TVector3 InterpolateAtZ(const std::vector<TrackAnalysis::Hit>& hits, double z_query) {
    // hits must be sorted in Z
    for (size_t i = 1; i < hits.size(); ++i) {
        double z1 = hits[i-1].position.Z();
        double z2 = hits[i].position.Z();
        if ((z_query >= z1 && z_query <= z2) || (z_query >= z2 && z_query <= z1)) {
            double t = (z_query - z1) / (z2 - z1);
            TVector3 pos = hits[i-1].position + t * (hits[i].position - hits[i-1].position);
            return pos;
        }
    }
    // If out of range, return first or last
    if (z_query < hits.front().position.Z()) return hits.front().position;
    if (z_query > hits.back().position.Z()) return hits.back().position;
    return hits.front().position;
}
///////////////////////////////////////////////////////// 
///// RunKalmanFilter method using hits.          ///////  
/////////////////////////////////////////////////////////
void TrackAnalysis::RunKalmanFilter(const std::vector<Hit>& hits, double seed_p) {

     if (hits.size() < 4) {
        results_.back().kalmanMomentum = 0;
        results_.back().kalmanChi2 = 9999;
        return;
    }
    
    const int nStates = 5; // x, y, z, tx, ty
    const int nMeas = 3;   // x, y, z measurements
    
    TMatrixD state(nStates, 1);

    //TMatrixD cov(nStates, nStates);
    TMatrixDSym cov(nStates);
    state(0,0) = hits[0].position.X();
    state(1,0) = hits[0].position.Y();
    state(2,0) = hits[0].position.Z();
    double dz = hits[1].position.Z() - hits[0].position.Z();
    state(3,0) = (fabs(dz) > 1e-6) ? (hits[1].position.X() - hits[0].position.X()) / dz : 0.0;
    state(4,0) = (fabs(dz) > 1e-6) ? (hits[1].position.Y() - hits[0].position.Y()) / dz : 0.0;


    // Initialize covariance matrix
    cov.Zero();
    cov(0,0) = cov(1,1) = cov(2,2) = 0.01; // position variance
    cov(3,3) = cov(4,4) = 1e-4;            // slopes variance

    double chi2 = 0;
    int nUpdates = 0;

    // For momentum estimation
    double total_Bdl = 0.0;
    double ty_entry = 0.0, ty_exit = 0.0;
    bool found_entry = false, found_exit = false;

    for (size_t i = 1; i < hits.size(); ++i) {
        double z_prev = state(2,0);
        double z_next = hits[i].position.Z();
        double dz = z_next - z_prev;

        // Check if segment crosses a magnet
        bool crosses_magnet = IsMagnetExist(z_prev, z_next);
        if (crosses_magnet) {
            double x_avg = 0.5 * (state(0,0) + hits[i].position.X());
            double y_avg = 0.5 * (state(1,0) + hits[i].position.Y());
            double z_avg = 0.5 * (z_prev + z_next);
            double bx = GetFieldX(x_avg, y_avg, z_avg);

            double dz_m = std::abs(dz) * 1e-3; // mm to m
            total_Bdl += std::abs(bx) * dz_m;

            // Use kasaMomentum as the momentum for propagation
            double p_est = seed_p;
            double dtheta = 0.3 * bx * dz_m / p_est; // radians
            state(4,0) += dtheta; // update ty (bending in y-z plane)

            // Record entry/exit slopes
            if (!found_entry) {
                ty_entry = state(4,0);
                found_entry = true;
            }
            ty_exit = state(4,0);
        }

        // Prediction
       //double dz = hits[i].position.Z() - state(2,0);
        state(0,0) += dz * state(3,0); // x += dz*tx
        state(1,0) += dz * state(4,0); // y += dz*ty
        state(2,0) = z_next;

        // Measurement matrix H (3x5)
        TMatrixD H(nMeas, nStates);
        H.Zero();
        H(0,0) = H(1,1) = H(2,2) = 1.0; // maps states to measurement x, y, z
        
        TVectorD meas(nMeas);
        meas(0) = hits[i].position.X();
        meas(1) = hits[i].position.Y();
        meas(2) = hits[i].position.Z();
        
        // Calculate H * state first
        TVectorD Hstate(nMeas);
        //Hstate = H * state;
        for (int m = 0; m < nMeas; m++) {
            Hstate(m) = 0;
            for (int s = 0; s < nStates; s++) {
                Hstate(m) += H(m,s) * state(s,0);
            }
        }
        // Calculate residual
        TVectorD residual = meas - Hstate; // residual = meas - H * state

        // Measurement covariance (3x3)
        TMatrixDSym R(nMeas);
        R.Zero();
        R(0,0) = R(1,1) = R(2,2) = 0.01; // 1 cm^2 resolution
        // R is diagonal, so we can use TMatrixDSym for efficiency

        // Innovation covariance S = H * cov * H.T() + R
        TMatrixD H_T(TMatrixD::kTransposed, H);
        TMatrixD temp = H * cov;
        TMatrixD S = temp * H_T;
        S += R; // Add measurement noise covariance

        if (S.Determinant() == 0) {
            std::cerr << "Error: S matrix is singular!" << std::endl;
            continue;
        }
        // Invert S using TMatrixD
        TMatrixD S_inv = S;
        S_inv.Invert();

         // Then compute residual properly
        //TVectorD residual(nMeas);
        //for (int j = 0; j < nMeas; ++j) {
        //    residual(j) = meas(j) - Hstate(j);
       // }

        // Kalman gain  K = cov * H.T() * S_inv
        TMatrixD K = cov * H_T; // (5x5)*(5x3) = (5x3)
         K *= S_inv; // (5x3)*(3x3) = (5x3)
        
        // Update state
        TVectorD state_update = K * residual;  // K (5x3) * residual (3x1) → (5x1) vector
        for (int j = 0; j < nStates; ++j) {
            state(j, 0) += state_update(j);    // Update state using vector indexing
        }

        // Update covariance - Joseph form for numerical stability
        TMatrixD I(TMatrixD::kUnit, cov); // Identity matrix (5x5)
        TMatrixD KH = K * H;              // (5x3)*(3x5) → (5x5)
        TMatrixD IminusKH = I - KH;       // (5x5) - (5x5) → (5x5)

        // First term: IminusKH * cov * IminusKH^T
        TMatrixD tempCov = IminusKH * cov;
        TMatrixD IminusKH_T(TMatrixD::kTransposed, IminusKH);
        TMatrixD cov1 = tempCov * IminusKH_T;  // (5x5)*(5x5) → (5x5)
        
        // Second term: K * R * K^T
        TMatrixD tempKR = K * R;               // (5x3)*(3x3) → (5x3)
        TMatrixD K_T(TMatrixD::kTransposed, K); // (3x5)
        TMatrixD cov2 = tempKR * K_T;          // (5x3)*(3x5) → (5x5)

        // Final covariance update
       for (int ii = 0; ii < nStates; ++ii) {
            for (int jj = 0; jj < nStates; ++jj) {
                cov(ii, jj) = cov1(ii, jj) + cov2(ii, jj);
            }
        }
        int DEBUG = 0; // Set to 1 to enable debug output
        // Debug output
        if (DEBUG) {
            std::cout << "[DEBUG] Kalman update for hit " << i << std::endl;
            std::cout << "H: " << H.GetNrows() << "x" << H.GetNcols() << std::endl;
            std::cout << "R: " << R.GetNrows() << "x" << R.GetNcols() << std::endl;
            std::cout << "S: " << S.GetNrows() << "x" << S.GetNcols() << std::endl;
            std::cout << "S_inv: " << S_inv.GetNrows() << "x" << S_inv.GetNcols() << std::endl;
            //std::cout << "residual: " << residual.GetNrows() << std::endl;
            std::cout << "K: " << K.GetNrows() << "x" << K.GetNcols() << std::endl;
            std::cout << "KH:  " << KH.GetNrows() << "x" << KH.GetNcols() << std::endl;
            std::cout << "cov: " << cov.GetNrows() << "x" << cov.GetNcols() << std::endl;
        }
        // Chi2 calculation (equivalent to residual^T * S_inv * residual)
        double chi2Cont = 0.0;
        for (int i = 0; i < nMeas; ++i) {
            for (int j = 0; j < nMeas; ++j) {
                chi2Cont += residual(i) * S_inv(i, j) * residual(j);
            }
        }
        chi2 += chi2Cont;
        nUpdates++;
    }

    // Final momentum estimation
    double kalman_p = 0;
    double dtheta = ty_exit - ty_entry;
    if (found_entry && fabs(dtheta) > 1e-8 && fabs(total_Bdl) > 1e-8) {
        kalman_p = 0.3 * fabs(total_Bdl) / fabs(dtheta);
        if(verbose >=2)
            std::cout << "[KALMAN] ty_entry = " << ty_entry
                  << ", ty_exit = " << ty_exit
                  << ", dtheta = " << dtheta
                  << ", total_Bdl = " << total_Bdl
                  << ", p = " << kalman_p << " GeV/c" << std::endl;
    } else {
        kalman_p = 0;
    }
    results_.back().kalmanMomentum = kalman_p;
    results_.back().kalmanChi2 = (nUpdates > 0) ? chi2/nUpdates : 9999;

}
///////////////////////////////////////////////////////// 
///// RunGenfit method using hits                 ///////
/////////////////////////////////////////////////////////
void TrackAnalysis::RunGenfit(const std::vector<Hit>& hits) 
{
    if (hits.size() < 4) {
        results_.back().genfitMomentum = 0;
        results_.back().genfitChi2 = 9999;
        results_.back().genfitCharge = 0;
        return;
    }

    // Sort hits by Z position (critical for GENFIT)
    std::vector<Hit> sortedHits = hits;
    std::sort(sortedHits.begin(), sortedHits.end(), 
        [](const Hit& a, const Hit& b) { return a.position.Z() < b.position.Z(); });

    // Print hit information for debugging
    // Debug print hit positions
    std::cout << "\n==== GENFIT Setup ====\n";
    std::cout << "Track ID: " << sortedHits[0].trackID << "\n";
    for (const auto& hit : sortedHits) {
        std::cout << "[HIT] (x,y,z in mm): " 
              << hit.position.X() << ", " 
              << hit.position.Y() << ", " 
              << hit.position.Z() << std::endl;
    }


    std::vector<Hit> augmentedHits = sortedHits; // Start with measured hits
    // For each magnet region, insert a virtual hit at entry and exit
    for (const auto& rng : magnet_z_ranges_) {
        double z_in = rng.first;
        double z_out = rng.second;
        // Add entry point
        TVector3 pos_in = InterpolateAtZ(sortedHits, z_in);
        Hit vh_in = {
            .stationID = -99, // or another flag value
            .layerID = -99,
            .trackID = sortedHits.front().trackID,
            .pdg = sortedHits.front().pdg,
            .position = pos_in,
            .momentum = TVector3() // not needed for fitting
        };
        augmentedHits.push_back(vh_in);
        // Add exit point
        TVector3 pos_out = InterpolateAtZ(sortedHits, z_out);
        Hit vh_out = vh_in;
        vh_out.position = pos_out;
        augmentedHits.push_back(vh_out);
    }
    // Now sort augmentedHits by Z again
    std::sort(augmentedHits.begin(), augmentedHits.end(), 
        [](const Hit& a, const Hit& b) { return a.position.Z() < b.position.Z(); });

    // Print augmented hits for debugging
    std::cout << "\n==== Augmented Hits ====\n";
    for (const auto& hit : augmentedHits) {
        std::cout << "[AUGMENTED HIT] (x,y,z in mm): " 
              << hit.position.X() << ", " 
              << hit.position.Y() << ", " 
              << hit.position.Z() << std::endl;
    }


    //Replace every usage of sortedHits with augmentedHits

    // ---- Seed state: position from first hit, momentum (direction) from first/last ----
    TVector3 dir = (augmentedHits.back().position - augmentedHits.front().position).Unit();
    double truthP = augmentedHits.front().momentum.Mag() * MeV_to_GeV;
    //double seed_p = 10.0; // 5 GeV/c; // Use a fixed seed momentum for now, can be improved later
    double seed_p = truthP; // 5 GeV/c; // Use a fixed seed momentum for now, can be improved later

    for (const auto& hit : augmentedHits) {
    double bx = GetFieldX(hit.position.X(), hit.position.Y(), hit.position.Z());
    std::cout << "Hit at z=" << hit.position.Z() << "mm: Bx=" << bx << " T" << std::endl;
    }

    // Positions in cm, momenta in GeV/c
    TVectorD stateSeed(6);
    stateSeed[0] = augmentedHits[0].position.X()* mm_to_cm;;
    stateSeed[1] = augmentedHits[0].position.Y()* mm_to_cm;;
    stateSeed[2] = augmentedHits[0].position.Z()* mm_to_cm;;
    stateSeed[3] = seed_p * dir.X();
    stateSeed[4] = seed_p * dir.Y();
    stateSeed[5] = seed_p * dir.Z();
    //stateSeed[3] = dir.X() * truthP;
    //stateSeed[4] = dir.Y() * truthP;
    //stateSeed[5] = dir.Z() * truthP;

    // ---- Covariance: generous errors ----
    TMatrixDSym covSeed(6);
    covSeed.Zero();
    // Position errors in cm^2, momentum errors in (GeV/c)^2
    covSeed(0,0) = covSeed(1,1) = covSeed(2,2) = 0.01;  // 1 mm2 --> 0.01 cm^2  resolution
    covSeed(3,3) = covSeed(4,4) = covSeed(5,5) = 1.;  // 1 (GeV/c)^2

    // Before fitting, print initial state:
    std::cout << "Initial seed state (cm, GeV/c):\n";
    for (int i=0; i<6; ++i) 
        std::cout << "  " << stateSeed[i] << "\n";
    std::cout << "Initial cov:\n";
    covSeed.Print();

    // Fit function for both charges
    auto fitCharge = [&](int pdg, double& momentum, double& chi2, bool& success) {
        success = false;
        momentum = 0;
        chi2 = 1e9;
        
        try {
           // Create track with representation
            genfit::Track* track = new genfit::Track(new genfit::RKTrackRep(pdg), stateSeed, covSeed);

            // Add measurements (convert mm→cm)
            for (const auto& hit : augmentedHits) {
                TVectorD pos(3);
                pos(0) = hit.position.X() * mm_to_cm;
                pos(1) = hit.position.Y() * mm_to_cm;
                pos(2) = hit.position.Z() * mm_to_cm;
                
                // ~1 mm resolution
                TMatrixDSym cov(3); 
                cov.UnitMatrix();
                //cov *= 0.01; // 0.01 cm² = (1 mm)²
                if (hit.stationID < 0) {
                   cov *= 1.0;   // 1 cm^2 uncertainty for virtual hits
                } else {
                    cov *= 0.01;  // 1 mm^2 for measured SciFi hits
                }
                track->insertMeasurement(new genfit::SpacepointMeasurement(pos, cov, hit.stationID, hit.layerID, nullptr));
            }

            // Create and configure fitter
            genfit::DAF daf;
            daf.setMaxIterations(20);
            daf.processTrack(track);

            const genfit::FitStatus* status = track->getFitStatus();
            if (status && status->isFitConverged()) {
                genfit::MeasuredStateOnPlane mop = track->getFittedState();
                // Extract momentum and chi2                
                momentum = mop.getMom().Mag();
                chi2 = status->getChi2() / status->getNdf();
                success = true;

                // Debug output
                std::cout << "\nSuccessful fit for PDG " << pdg << ":\n";
                std::cout << "  Momentum: " << momentum << " GeV/c\n";
                std::cout << "  Chi2/NDF: " << chi2 << "\n";
                std::cout << "  Final State:\n";
                mop.Print();
            } else {
                std::cout << "\nFit failed for PDG " << pdg << ": ";
                if (status) {
                    std::cout << "  Reason: " << status->isFitConverged() << "\n";
                    std::cout << "  Chi2: " << status->getChi2() << "\n";
                }
            }
            // Clean up
            delete track;
        } catch (const std::exception& e) {
            std::cerr << "Exception in fit: " << e.what() << "\n";
        }
    };
    // Fit both charges
    double p_minus, chi2_minus;
    double p_plus, chi2_plus;
    bool ok_minus, ok_plus;
    
    fitCharge(13, p_minus, chi2_minus, ok_minus);    // mu-
    fitCharge(-13, p_plus, chi2_plus, ok_plus);      // mu+

    // Store results
    if (ok_minus || ok_plus) {
        if (ok_minus && (!ok_plus || chi2_minus < chi2_plus)) {
            results_.back().genfitMomentum = p_minus;
            results_.back().genfitChi2 = chi2_minus;
            results_.back().genfitCharge = -1;
        } else {
            results_.back().genfitMomentum = p_plus;
            results_.back().genfitChi2 = chi2_plus;
            results_.back().genfitCharge = +1;
        }
    } else {
        results_.back().genfitMomentum = 0;
        results_.back().genfitChi2 = 9999;
        results_.back().genfitCharge = 0;
    }

}
///////////////////////////////////////////////////////////// 
///// CalculateSagittaMomentum method using hits.     ///////
/////////////////////////////////////////////////////////////
double TrackAnalysis::CalculateSagittaMomentum(const std::vector<Hit>& hits) {
    if (hits.size() < 3) return 0;

    const TVector3& p0 = hits[0].position;
    const TVector3& p1 = hits[hits.size()/2].position;
    const TVector3& p2 = hits.back().position;

    // Convert mm to meters!
    //double L = (p2.Z() - p0.Z()) * 0.001;  // mm to m
    // in the following i need to use the total magnetized length, not just the distance between first and last hit
    // Compute number of stations
    std::set<int> stations;
    for (const auto& hit : hits) stations.insert(hit.stationID);
    int nStations = stations.size();

    // Use iron thickness from class constant or geometry
    double iron_thickness = GetIronThickness(); // in meters

    // Effective L: total magnetized length
    double L = (nStations - 1) * iron_thickness;

    double sagitta = (p1.Y() - 0.5*(p0.Y() + p2.Y())) * 0.001; // mm to m
    //double B = std::abs(GetFieldX(p1.X(), p1.Y(), p1.Z())); // Tesla
    double B = std::abs(MagnetFieldStrength); // Tesla

    // Formula: p = 0.3 * B * L^2 / (8 * sagitta)
    double p = (std::abs(sagitta) > 1e-8) ? (0.3 * B * L * L / (8. * std::abs(sagitta))) : 0.0;
    if(verbose >= 2)
        std::cout << "[SagittaMomentum] p = " << p 
                  << " GeV/c, B = " << B 
                  << " T, L = " << L 
                  << " m, sagitta = " << sagitta 
                  << " m" << std::endl;

     // Determine charge
    double local_bx = GetFieldBX(p1.Y()); // Get field with sign
    int charge = 0;
    if (std::abs(sagitta) > 1e-8) {
        // For Bx field: positive sagitta with positive Bx means negative charge
        charge = (sagitta * local_bx > 0) ? -1 : +1;
    }
    results_.back().sagittaCharge = charge;
    // Calculate resolution (error propagation)
    double sigma_sagitta = 0.05e-3; // 50 μm resolution in meters
    double resolution = (std::abs(sagitta) > 1e-8) ? 
        p * sigma_sagitta / std::abs(sagitta) : -1;
    results_.back().sagittaResolution = resolution;
    if(verbose >= 2)
        std::cout << "[SagittaMomentum] Charge = " << charge 
              << " local Bx = " << local_bx
              << ", Resolution = " << resolution 
              << " GeV/c" << std::endl;

    return p; // (std::abs(sagitta) > 1e-8) ? (0.3 * B * L * L / (8. * std::abs(sagitta))) : 0.0; // GeV/c
}
///////////////////////////////////////////////////////////// 
///// CalculateSagittaMomentum method using tracklet. ///////
/////////////////////////////////////////////////////////////
double TrackAnalysis::CalculateSagittaMomentum(const std::vector<Tracklet>& hits) {
    if (hits.size() < 3) return 0;

    const TVector3& p0 = hits[0].position;
    const TVector3& p1 = hits[hits.size()/2].position;
    const TVector3& p2 = hits.back().position;

    // Convert mm to meters!
    //double L = (p2.Z() - p0.Z()) * 0.001;  // mm to m
    // in the following i need to use the total magnetized length, not just the distance between first and last hit
    // Compute number of stations
    std::set<int> stations;
    for (const auto& hit : hits) stations.insert(hit.station);
    int nStations = stations.size();

    // Use iron thickness from class constant or geometry
    double iron_thickness = GetIronThickness(); // in meters

    // Effective L: total magnetized length
    double L = (nStations - 1) * iron_thickness;
    // i have to use Y coordinate for sagitta calculation
    double sagitta = (p1.Y() - 0.5*(p0.Y() + p2.Y())) * 0.001; // mm to m
    double B = std::abs(MagnetFieldStrength); // Use absolute value for momentum
    //double B = std::abs(GetFieldX(p1.X(), p1.Y(), p1.Z())); // Tesla
    // Calculate momentum
    // Formula: p = 0.3 * B * L^2 / (8 * sagitta)
    double p = (std::abs(sagitta) > 1e-8) ? (0.3 * B * L * L / (8. * std::abs(sagitta))) : 0.0;
    //Determine charge using local signed field
    double local_bx = GetFieldBX(p1.Y()); // Get field with sign
    int charge = 0;
    if (std::abs(sagitta) > 1e-8 && std::abs(local_bx) > 1e-6) {
        charge = (sagitta * local_bx > 0) ? -1 : +1;
    }
    results_.back().sagittaCharge_trk = charge; 

    //Calculate resolution (error propagation)
    double sigma_sagitta = 0.05e-3; // 50 μm resolution in meters
    double resolution = (std::abs(sagitta) > 1e-8) ? 
        p * sigma_sagitta / std::abs(sagitta) : -1;
    results_.back().sagittaResolution_trk = resolution;
    if(verbose >= 2)
        std::cout << "[SagittaMomentum] p = " << p 
                  << " GeV/c, B = " << B 
                  << " T, L = " << L 
                  << " m, sagitta = " << sagitta 
                  << " m" << std::endl;
    return p;
}
///////////////////////////////////////////////////////// 
/////  FitStationTracklets method implementation. ///////
/////////////////////////////////////////////////////////
std::vector<TrackAnalysis::Tracklet> TrackAnalysis::FitStationTracklets(const std::vector<Hit>& hits) {
    //For every unique (stationID, trackID) pair in the event:
    //Collect all hits belonging to that station and track
    //Fit a straight line to the hits in x(z) and y(z), to find the local direction the particle took through that station.
    //Calculate the average position of the tracklet.
    //Store a Tracklet struct containing: station, average position, direction, trackID, pdg, and truth momentum.
    //Returns all these fitted tracklets as a vector.
    // Each tracklet gives you the local track direction before/after a magnetized block.
    // Group hits by stationID 
    //std::map<int, std::vector<const Hit*>> stationHits;
    //for (const auto& hit : hits) {
    //    stationHits[hit.stationID].push_back(&hit);
    //}
    // group hits by stationID and trackID
    std::map<std::pair<int, int>, std::vector<const Hit*>> tracklets;
    for (const auto& hit : hits) {
        tracklets[{hit.stationID, hit.trackID}].push_back(&hit);
    }

    // Prepare a vector to hold the tracklets
    std::vector<Tracklet> stationTracklets;

    if(verbose >= 2)
        std::cout << "---- Station-by-station Tracklet Fits ----" << std::endl;
    
    for (const auto& entry : tracklets) {
        int station = entry.first.first;
        int trackID = entry.first.second; 
        const auto& hitvec = entry.second;

        if (hitvec.size() < 2) continue; // Need at least two hits for a fit
        if( verbose >= 3)
            std::cout << "Processing tracklet for station " << station 
                      << " with trackID " << trackID 
                      << " with PDG" << hitvec.front()->pdg 
                      << std::endl;
        
        // sometimes i see more than 4 hits in a station for muon ... 
        // --- Median selection per layer ---
        std::map<int, std::vector<const Hit*>> hitsByLayer;
        for (const auto* hit : hitvec) {
            hitsByLayer[hit->layerID].push_back(hit);
        }
        std::vector<const Hit*> selectedHits;
        for (auto& layerHits : hitsByLayer) {
            auto& layer_vec = layerHits.second;
            std::sort(layer_vec.begin(), layer_vec.end(),
                      [](const Hit* a, const Hit* b) { return a->position.X() < b->position.X(); });
            // Pick median hit (if odd, it's clear; if even, pick lower middle)
            selectedHits.push_back(layer_vec[layer_vec.size()/2]);
        }
        if (selectedHits.size() < 2) continue; // Need at least two hits for a fit

        // Prepare sums for linear fit: x(z) and y(z)
        double sumZ = 0, sumZZ = 0, sumX = 0, sumY = 0, sumXZ = 0, sumYZ = 0;
        for (const auto* hit : selectedHits) {
            double x = hit->position.X(); // X coordinate in mm
            double y = hit->position.Y(); // Y coordinate in mm
            double z = hit->position.Z(); // Z coordinate in mm
            sumZ  += z;
            sumZZ += z*z;
            sumX  += x;
            sumY  += y;
            sumXZ += x*z;
            sumYZ += y*z;
        }
        double N = selectedHits.size();
        double denom = N * sumZZ - sumZ * sumZ;
        if (std::abs(denom) < 1e-12) continue; // avoid division by zero

        double slopeX = (N * sumXZ - sumZ * sumX) / denom;
        double interceptX = (sumX - slopeX * sumZ) / N;
        double slopeY = (N * sumYZ - sumZ * sumY) / denom;
        double interceptY = (sumY - slopeY * sumZ) / N;

        // --- Remove hits with large residuals and refit ---
        std::vector<const Hit*> inliers;
        double max_residual = 2.0; // mm
        for (const auto* hit : selectedHits) {
            double pred_x = slopeX * hit->position.Z() + interceptX;
            double pred_y = slopeY * hit->position.Z() + interceptY;
            double dx = hit->position.X() - pred_x;
            double dy = hit->position.Y() - pred_y;
            double res = std::sqrt(dx*dx + dy*dy);
            if (res < max_residual) inliers.push_back(hit);
        }
        // If enough inliers, refit
        if (inliers.size() >= 2 && inliers.size() < selectedHits.size()) {
            sumZ = sumZZ = sumX = sumY = sumXZ = sumYZ = 0;
            for (const auto* hit : inliers) {
                double x = hit->position.X();
                double y = hit->position.Y();
                double z = hit->position.Z();
                sumZ  += z;
                sumZZ += z*z;
                sumX  += x;
                sumY  += y;
                sumXZ += x*z;
                sumYZ += y*z;
            }
            N = inliers.size();
            denom = N * sumZZ - sumZ * sumZ;
            if (std::abs(denom) < 1e-12) continue;
            slopeX = (N * sumXZ - sumZ * sumX) / denom;
            interceptX = (sumX - slopeX * sumZ) / N;
            slopeY = (N * sumYZ - sumZ * sumY) / denom;
            interceptY = (sumY - slopeY * sumZ) / N;
            selectedHits = inliers;
        }
        // --- Error estimation for slopeX and slopeY ---
        double sum_res2_x = 0, sum_res2_y = 0;
        for (const auto* hit : selectedHits) {
            double pred_x = slopeX * hit->position.Z() + interceptX;
            double pred_y = slopeY * hit->position.Z() + interceptY;
            sum_res2_x += std::pow(hit->position.X() - pred_x, 2);
            sum_res2_y += std::pow(hit->position.Y() - pred_y, 2);
        }
        double sigma2_x = (N > 2) ? sum_res2_x / (N - 2) : 0;
        double sigma2_y = (N > 2) ? sum_res2_y / (N - 2) : 0;
        double slopeX_err = (denom > 0 && N > 2) ? std::sqrt(sigma2_x / (sumZZ - sumZ*sumZ/N)) : 0;
        double slopeY_err = (denom > 0 && N > 2) ? std::sqrt(sigma2_y / (sumZZ - sumZ*sumZ/N)) : 0;
        // Print the fit results
        if (verbose >= 3) {
            std::cout << "Fit results for station " << station 
                      << " with trackID " << trackID
                      << ": N=" << N
                      << ", slopeX=" << slopeX
                      << " ± " << slopeX_err
                      << ", interceptX=" << interceptX
                      << ", slopeY=" << slopeY
                      << " ± " << slopeY_err
                      << ", interceptY=" << interceptY
                      << std::endl;     
        }   

        // The direction vector 
        TVector3 dir(slopeX, slopeY, 1.0);
        dir = dir.Unit();

        // Average position (unit in mm)
        TVector3 position(0,0,0);
        for (const auto* hit : selectedHits)
            position += hit->position;
        position *= 1.0 / selectedHits.size();
        
        int pdg = selectedHits.front()->pdg;
        double truthP = selectedHits.front()->momentum.Mag() * MeV_to_GeV; // If in MeV/c, convert to GeV/c

        // Store the tracklet
        stationTracklets.push_back({station, position, dir, pdg, trackID, truthP});

        // Print the results
         if (verbose >= 3)
            std::cout << "Tracklet fit for station " << station
                      << " with trackID " << trackID
                      << ": N=" << N
                      << ", slopeX=" << slopeX << " ± " << slopeX_err << " mm/mm"
                      << ", slopeY=" << slopeY << " ± " << slopeY_err << " mm/mm"
                      << ", direction=(" << dir.X() << ", " << dir.Y() << ", " << dir.Z() << ")"
                      << ", avg pos (mm): (" << position.X() << ", " << position.Y() << ", " << position.Z() << ")"
                      << std::endl;
       
        // print all hits in this station
        if (verbose >=3)  
        for (size_t i = 0; i < hitvec.size(); ++i) {
            const Hit* hit = hitvec[i];
            std::cout << "    all hits in given station " << i
                      << " trackID=" << hit->trackID
                      << ", pdg=" << hit->pdg
                      << ": layer=" << hit->layerID 
                      << ", (" << hit->position.X() 
                      << ", " << hit->position.Y() 
                      << ", " << hit->position.Z() << ")" 
                      << std::endl;
        }
    }
        return stationTracklets;
}
///////////////////////////////////////////////////////// 
///// use tracklets to estimate momentum.         ///////
/////////////////////////////////////////////////////////
//[tracklet1]--[iron block]--[tracklet2]--[iron block]--[tracklet3]--...
//Tracklet N: SciFi hits, no field, fit direction before the iron.
//Tracklet N+1: SciFi hits, after the iron, fit direction after the iron.
//Iron block: Bx field, known thickness, known Bx value.
//The deflection angle deltaTheta is just the angle between the direction vectors 
//before and after the iron block in the bending plane (here, y-z, since the field is in x).
// Keep in mind that the magnetic field is changing sign in the iron block, 
// check the positions of the tarcklets to define the magnetic field sign.
//std::pair<double, double> 
std::tuple<double, double, int> TrackAnalysis::MomentumFromTrackletDeflection(
    const std::vector<Tracklet>& tracklets,
    double bx_field,         // e.g. -1.5 T
    double iron_thickness_m,  // e.g. 0.10 (10 cm)
    bool use_3D_angle        // true = 3D angle, false = 2D bending plane
) {
    if (tracklets.size() < 2) return {-999., -1., 0};

    std::vector<double> momenta, theta_errs;
    std::vector<double> weights;
    double sum_w = 0, sum_pw = 0;
    double ang_res = 0.0058; // radians, for 4 layers at 2.5mm spacing, 50um resolution    //double Bld; // Bld = Bx * iron_thickness_m
    
    for (size_t i = 0; i + 1 < tracklets.size(); ++i) {
        double z1 = tracklets[i].position.Z();
        double z2 = tracklets[i+1].position.Z();
        
        if(verbose >= 4)    
        std::cout << "[Tracklet-1] " << i
                  << " trackID=" << tracklets[i].trackID  
                  << ", pdg=" << tracklets[i].pdg 
                  << " at z1=" << z1 
                  << ", position: (" 
                  << tracklets[i].position.X() << ", " 
                  << tracklets[i].position.Y() << ", " 
                  << tracklets[i].position.Z() << ")" 
                  << std::endl;
        if(verbose >= 4)    
        std::cout << "[Tracklet-2] " << i+1 
                  << " trackID=" << tracklets[i+1].trackID
                  << ", pdg=" << tracklets[i+1].pdg
                  << " at z2=" << z2 
                  << ", position: (" 
                  << tracklets[i+1].position.X() << ", " 
                  << tracklets[i+1].position.Y() << ", " 
                  << tracklets[i+1].position.Z() << ")"
                  << std::endl;
        // Only consider segments that cross a magnet
        bool has_magnet = IsMagnetExist(z1, z2);
        if(verbose >= 2)
        std::cout << "[TRACKLET DEFLECTION] Processing tracklet " << i 
                  << " to " << i+1 << ", z1=" << z1 << ", z2=" << z2 
                  << ", has_magnet_in_between=" << has_magnet << std::endl;
        if (!has_magnet) continue; // No magnet between these tracklets

        // Compute average position for field query
        double y_avg = 0.5 * (tracklets[i].position.Y() + tracklets[i+1].position.Y()); // in mm
        double x_avg = 0.5 * (tracklets[i].position.X() + tracklets[i+1].position.X()); // in mm
        double z_avg = 0.5 * (tracklets[i].position.Z() + tracklets[i+1].position.Z()); // in mm      
       
       // Get local Bx at this segment (in Tesla)
        double bx_local = GetFieldX(x_avg, y_avg, z_avg);
        if (verbose >= 4)
            std::cout << "[TRACKLET DEFLECTION] Bx at segment: " << bx_local 
                      << " T, positions of " << x_avg << ", " << y_avg << ", " << z_avg << std::endl;        
        if (verbose >= 4)
        std::cout << "[DEBUG] z1 = " << z1 << ", z2 = " << z2 << ", dz = " << (z2-z1) << std::endl;
        
        // Segment length inside magnetic field is the same as the iron thickness
        //double dz_m = std::abs(z2 - z1) * 1e-3;
        //double dz_m = iron_thickness_m; // in meters, e.g. 0.10 m for 10 cm iron block
        double dz_m = 0.0;
        for (const auto& rng : magnet_z_ranges_) {
            double zin = std::max(rng.first*10.0, std::min(z1, z2));   // cm to mm
            double zout = std::min(rng.second*10.0, std::max(z1, z2)); // cm to mm
            if (zout > zin) dz_m += (zout - zin) * 1e-3; // mm to m
        }
        if (dz_m < 1e-6) continue; // skip if no overlap

        // Local Bdl for this segment
        double Bdl = bx_local * dz_m; // Tesla*m
        if (verbose >= 4)
            std::cout << "[TRACKLET DEFLECTION] Bdl for segment: " << Bdl 
                      << " T*m, dz_m = " << dz_m << " m" << std::endl;

        // Calculate the deflection angle in the bending plane (y-z)
        double delta_theta = 0;
        if (use_3D_angle) {
            // Using full 3D angle between direction vectors
            const TVector3& dir1 = tracklets[i].dir;
            const TVector3& dir2 = tracklets[i+1].dir;
            double cos_angle = dir1.Dot(dir2) / (dir1.Mag() * dir2.Mag());
            cos_angle = std::max(-1.0, std::min(1.0, cos_angle)); 
            delta_theta = std::acos(cos_angle); 
        } else {
            // Using change in slope in bending (y-z) plane
            double tan_in  = tracklets[i].dir.Y() / tracklets[i].dir.Z();
            double tan_out = tracklets[i+1].dir.Y() / tracklets[i+1].dir.Z();
            delta_theta = tan_out - tan_in; 
        }

        if (std::abs(delta_theta) > 1e-7 && std::abs(Bdl) > 1e-8) {
            double p = 0.3 * std::abs(Bdl) / std::abs(delta_theta); // Use abs(Bdl)
            // Calculate error on angle
            double theta_err = ang_res;
            double p_err = p * theta_err / std::abs(delta_theta);
            // Outlier rejection: skip segments with tiny angle or huge momentum
            if (std::abs(delta_theta) < 0.005 || p > 200.0) {
                if (verbose >= 1)
                    std::cout << "[TRACKLET DEFLECTION] Skipping outlier: Δθ = " << delta_theta << ", p = " << p << std::endl;
                continue;
            }
            
            double w = 1. / (p_err*p_err);
            momenta.push_back(p);
            theta_errs.push_back(p_err);
            weights.push_back(w);
            sum_w += w;
            sum_pw += w * p;
            if(verbose >=1)
                std::cout << "[TRACKLET DEFLECTION] Δθ = " << delta_theta 
                          << " rad, Bdl = " << Bdl 
                          << ", p = " << p << " GeV/c, p_err = " << p_err 
                          << ", weight = " << w 
                          << std::endl;
        }
    }

    // Return average, median, 
    if (momenta.empty()) return {-999., -1., 0};
    //double avg = sum_pw / sum_w;
    //double var = 0;
    //for (size_t i = 0; i < momenta.size(); ++i)
    //    var += weights[i] * (momenta[i] - avg) * (momenta[i] - avg);
    //double rms = sqrt(var / sum_w);
    std::sort(momenta.begin(), momenta.end());
    double median = momenta[momenta.size()/2];
    double rms = 0.0;
    for (auto p : momenta) rms += (p - median) * (p - median);
    rms = std::sqrt(rms / momenta.size());

    // Add charge determination from overall bending
    int overall_charge = 0;
    if (tracklets.size() >= 2) {
        TVector3 dir_first = tracklets.front().dir;
        TVector3 dir_last = tracklets.back().dir;
        // Change in Y slope indicates charge for Bx field
        double delta_ty = dir_last.Y()/dir_last.Z() - dir_first.Y()/dir_first.Z();
        // Get local field with sign at middle of track
        double z_mid = 0.5 * (tracklets.front().position.Z() + tracklets.back().position.Z());
        double y_mid = 0.5 * (tracklets.front().position.Y() + tracklets.back().position.Y());
        double x_mid = 0.5 * (tracklets.front().position.X() + tracklets.back().position.X());
        //double local_bx = GetFieldX(x_mid, y_mid, z_mid);
        double local_bx = GetFieldBX(y_mid); // Get field with sign

        overall_charge = (delta_ty * local_bx > 0) ? -1 : +1;
    }
    //std::cout << "[TRACKLET DEFLECTION] Weighted average momentum: " << avg << " GeV/c, RMS: " << rms << std::endl;
    std::cout << "[TRACKLET DEFLECTION] Median momentum: " << median << " GeV/c, RMS: " << rms << std::endl;
    std::cout << "[TRACKLET DEFLECTION] Overall charge: " << overall_charge << std::endl;
    //return {avg, rms};
    //return {median, rms, overall_charge}; // Return median, RMS, and overall charge as a tuple
    return std::make_tuple(median, rms, overall_charge); // Return median and overall charge
}

///////////////////////////////////////////////////////// 
///// estimate momentum using global direction method ///
/////////////////////////////////////////////////////////
double TrackAnalysis::GlobalDirectionDeflectionMomentum(
    const std::vector<Tracklet>& tracklets,
    double bx_field,           // e.g. -1.5 T
    double total_iron_thickness_m // e.g. 1.0 (10 blocks × 0.10m)
) {
    if (tracklets.size() < 2) return -999.;

   // Fit a global line to all tracklet positions (Z-X and Z-Y)
    std::vector<TVector3> positions;
    for (const auto& t : tracklets) positions.push_back(t.position);

    double slope_xz = 0, intcpt_xz = 0, slope_yz = 0, intcpt_yz = 0;
    LinearFit(positions, slope_xz, intcpt_xz, 'x');
    LinearFit(positions, slope_yz, intcpt_yz, 'y');

    // Build direction vectors at the first and last tracklet positions
    double z_first = tracklets.front().position.Z();
    double z_last  = tracklets.back().position.Z();
    TVector3 dir_first(slope_xz, slope_yz, 1.0); dir_first = dir_first.Unit();
    TVector3 dir_last(slope_xz, slope_yz, 1.0);  dir_last = dir_last.Unit();

    //fitted directions from the first and last tracklet
    // TVector3 dir_first = tracklets.front().dir;
    // TVector3 dir_last  = tracklets.back().dir;

    // Compute the total deflection angle in 3D
    double cos_angle = dir_first.Dot(dir_last) / (dir_first.Mag() * dir_last.Mag());
    cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
    double delta_theta = std::acos(cos_angle);

    // Total Bdl (field × thickness)
    double Bdl = bx_field * total_iron_thickness_m; // all blocks together

    // Momentum formula
    double p = (std::abs(delta_theta) > 1e-7 && std::abs(Bdl) > 1e-8)
        ? (0.3 * std::abs(Bdl) / std::abs(delta_theta))
        : -999.;

    if (verbose >= 1)
        std::cout << "[GLOBAL DIR DEFLECTION] Δθ = " << delta_theta
                  << " rad, Bdl = " << Bdl
                  << ", p = " << p << " GeV/c" << std::endl;

    return p;
}
////////////////////////////////////////////////////////////// 
/////CircleFitYZ_Taubin method implementation on Tracklet/////
//////////////////////////////////////////////////////////////
// Circle fit in (y, z) plane for Bx field
// Returns the momentum in GeV/c
// fitting the (Y, Z) plane for Bx field.
// The circle fit reconstructs the trajectory's curvature, which is directly related to the momentum via the above formula.
TrackAnalysis::CircleFitResult TrackAnalysis::CircleFitYZ_Taubin(const std::vector<Tracklet>& tracklets, double bx_field) {
    CircleFitResult result = {0, 0, 0, 0, 0, 0, 0, 0, false};
    if (tracklets.size() < 3) return result;
    // Fit a circle in the (Y, Z) plane using Taubin's method
    // Returns the circle center (yc, zc), radius R, RMS error sigma,
    std::vector<double> y, z;
    for (const auto& t : tracklets) {
        y.push_back(t.position.Y() * 1e-3); // Convert mm to m
        z.push_back(t.position.Z() * 1e-3);
    }
    int N = y.size();

    // -- Compute means
    double meanY = 0, meanZ = 0;
    for (int i=0; i<N; ++i) { meanY += y[i]; meanZ += z[i]; }
    meanY /= N; meanZ /= N;

    // -- Build moment sums
    double Mzz=0, Myy=0, Myz=0, Myz2=0, Mz2=0, My2=0, Myz_y=0, Myz_z=0;
    for (int i=0; i<N; ++i) {
        double Yi = y[i] - meanY, Zi = z[i] - meanZ, Ri2 = Yi*Yi + Zi*Zi;
        Myy += Yi*Yi;      // variance Y
        Mzz += Zi*Zi;      // variance Z
        Myz += Yi*Zi;
        Myz2 += Yi*Ri2;
        Myz_z += Zi*Ri2;
        Myz_y += Yi*Ri2;
        Mz2  += Ri2;
        My2  += Ri2*Ri2;
    }
    Myy /= N; 
    Mzz /= N; 
    Myz /= N; 
    Myz2 /= N; 
    Myz_z /= N; 
    Myz_y /= N; 
    Mz2 /= N; 
    My2 /= N;

    // -- Taubin algebraic fit 
    double Mz = Myy + Mzz;
    double Cov_yz = Myy * Mzz - Myz * Myz;
    double Var_r = My2 - Mz*Mz;
    double A3 = 4*Mz;
    double A2 = -3*Mz*Mz - My2;
    double A1 = Var_r*Mz + 4*Cov_yz*Mz - Myz2*Myz2 - Myz_z*Myz_z;
    double A0 = Myz2*(Myz2*Mzz - Myz_z*Myz) + Myz_z*(Myz_z*Myy - Myz2*Myz) - Var_r*Cov_yz;
    double x = 0, ynew = A0, xnew, Dy;
    int iter, IterMAX=99;
    for (iter=0; iter<IterMAX; ++iter) {
        Dy = A1 + x*(2*A2 + 3*A3*x);
        if (Dy == 0) break;
        xnew = x - ynew/Dy;
        if (!std::isfinite(xnew) || fabs(xnew-x)<1e-12) break;
        ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
        if (fabs(ynew)>=fabs(A0)) break;
        x = xnew;
    }
    double DET = x*x - x*Mz + Cov_yz;
    if (fabs(DET) < 1e-16) return {0,0,0,0,0,0,false};
    double yc = (Myz2*(Mzz-x) - Myz_z*Myz) / DET / 2 + meanY;
    double zc = (Myz_z*(Myy-x) - Myz2*Myz) / DET / 2 + meanZ;
    double R = 0;
    std::vector<double> residuals;
    for (int i=0; i<N; ++i) {
        double dist = std::sqrt((y[i]-yc)*(y[i]-yc) + (z[i]-zc)*(z[i]-zc));
        R += dist;
        residuals.push_back(dist);
    }
    R /= N;

    // -- RMS error
    double sigma = 0;
    for (int i=0; i<N; ++i) {
        sigma += (residuals[i] - R)*(residuals[i] - R);
    }
    sigma = sqrt(sigma/N);
    int ndf = N - 3;

    // Calculate momentum and its error
    double p = 0.3 * std::abs(bx_field) * R; // GeV/c
    double p_err = 0.3 * std::abs(bx_field) * sigma;
    
    // Determine charge from bending direction
    int charge = 0;
    if (tracklets.size() >= 3) {
        // Use first, middle, last points to determine bending direction
        TVector3 p1 = tracklets.front().position;
        TVector3 p2 = tracklets[tracklets.size()/2].position;
        TVector3 p3 = tracklets.back().position;

        // Calculate sagitta in Y direction (for Bx field)
        double sagitta_y = p2.Y() - 0.5 * (p1.Y() + p3.Y());
        
        // Get local field with sign at middle point
        double local_bx = GetFieldBX(p2.Y()); // Get field with sign

        if (std::abs(local_bx) > 1e-6) {
            charge = (sagitta_y * local_bx > 0) ? -1 : +1;
        }
    }
    
    results_.back().taubinCharge_trk = charge;  // ADD this line  
    results_.back().taubinResolution_trk = p_err;
    
    result = {yc, zc, R, sigma, ndf, p, p_err, charge, true};
    return result;
}
////////////////////////////////////////////////////////////// 
/////CircleFitYZ_Taubin method implementation on Hits/////
//////////////////////////////////////////////////////////////
// Circle fit in (y, z) plane for Bx field
// Returns the momentum in GeV/c
// fitting the (Y, Z) plane for Bx field.
TrackAnalysis::CircleFitResult TrackAnalysis::CircleFitYZ_Taubin(const std::vector<Hit>& hits, double bx_field) {
    CircleFitResult result = {0, 0, 0, 0, 0, 0, 0, false};
    if (hits.size() < 3) return result;
    // Fit a circle in the (Y, Z) plane using Taubin's method
    // Returns the circle center (yc, zc), radius R, RMS error sigma,
    std::vector<double> y, z;
    for (const auto& t : hits) {
        y.push_back(t.position.Y() * 1e-3); // Convert mm to m
        z.push_back(t.position.Z() * 1e-3);
    }
    int N = y.size();

    // -- Compute means
    double meanY = 0, meanZ = 0;
    for (int i=0; i<N; ++i) { meanY += y[i]; meanZ += z[i]; }
    meanY /= N; meanZ /= N;

    // -- Build moment sums
    double Mzz=0, Myy=0, Myz=0, Myz2=0, Mz2=0, My2=0, Myz_y=0, Myz_z=0;
    for (int i=0; i<N; ++i) {
        double Yi = y[i] - meanY, Zi = z[i] - meanZ, Ri2 = Yi*Yi + Zi*Zi;
        Myy += Yi*Yi;      // variance Y
        Mzz += Zi*Zi;      // variance Z
        Myz += Yi*Zi;
        Myz2 += Yi*Ri2;
        Myz_z += Zi*Ri2;
        Myz_y += Yi*Ri2;
        Mz2  += Ri2;
        My2  += Ri2*Ri2;
    }
    Myy /= N; 
    Mzz /= N; 
    Myz /= N; 
    Myz2 /= N; 
    Myz_z /= N; 
    Myz_y /= N; 
    Mz2 /= N; 
    My2 /= N;

    // -- Taubin algebraic fit 
    double Mz = Myy + Mzz;
    double Cov_yz = Myy * Mzz - Myz * Myz;
    double Var_r = My2 - Mz*Mz;
    double A3 = 4*Mz;
    double A2 = -3*Mz*Mz - My2;
    double A1 = Var_r*Mz + 4*Cov_yz*Mz - Myz2*Myz2 - Myz_z*Myz_z;
    double A0 = Myz2*(Myz2*Mzz - Myz_z*Myz) + Myz_z*(Myz_z*Myy - Myz2*Myz) - Var_r*Cov_yz;
    double x = 0, ynew = A0, xnew, Dy;
    int iter, IterMAX=99;
    for (iter=0; iter<IterMAX; ++iter) {
        Dy = A1 + x*(2*A2 + 3*A3*x);
        if (Dy == 0) break;
        xnew = x - ynew/Dy;
        if (!std::isfinite(xnew) || fabs(xnew-x)<1e-12) break;
        ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
        if (fabs(ynew)>=fabs(A0)) break;
        x = xnew;
    }
    double DET = x*x - x*Mz + Cov_yz;
    if (fabs(DET) < 1e-16) return {0,0,0,0,0,0,false};
    double yc = (Myz2*(Mzz-x) - Myz_z*Myz) / DET / 2 + meanY;
    double zc = (Myz_z*(Myy-x) - Myz2*Myz) / DET / 2 + meanZ;
    double R = 0;
    std::vector<double> residuals;
    for (int i=0; i<N; ++i) {
        double dist = std::sqrt((y[i]-yc)*(y[i]-yc) + (z[i]-zc)*(z[i]-zc));
        R += dist;
        residuals.push_back(dist);
    }
    R /= N;

    // -- RMS error
    double sigma = 0;
    for (int i=0; i<N; ++i) {
        sigma += (residuals[i] - R)*(residuals[i] - R);
    }
    sigma = sqrt(sigma/N);
    
    int ndf = N - 3;

    // Calculate momentum and its error
    double p = 0.3 * std::abs(bx_field) * R; // GeV/c
    double p_err = 0.3 * std::abs(bx_field) * sigma;

    // Determine charge from bending direction
    int charge = 0;
    if (hits.size() >= 3) {
        TVector3 p1 = hits.front().position;
        TVector3 p2 = hits[hits.size()/2].position;
        TVector3 p3 = hits.back().position;
        
        double sagitta_y = p2.Y() - 0.5 * (p1.Y() + p3.Y());
        
        // Get local field with sign at middle point
        //double local_bx = GetFieldX(p2.X(), p2.Y(), p2.Z());
        double local_bx = GetFieldBX(p2.Y()); // Get field with sign

        if (std::abs(local_bx) > 1e-6) {
            charge = (sagitta_y * local_bx > 0) ? -1 : +1;
        }
    }
    //std::cout << "[CircleFitYZ_Taubin] yc=" << yc << ", zc=" << zc
    //          << ", R=" << R << ", sigma=" << sigma
    //          << ", ndf=" << ndf << ", p=" << p << ", p_err=" << p_err
    //          << ", charge=" << charge << std::endl;

    result = {yc, zc, R, sigma, ndf, p, p_err, charge, true};
    return result;
}
///////////////////////////////////////////////////////// 
/////CircleFitYZ_Kasa method implementation on Hits//////
/////////////////////////////////////////////////////////
// Circle fit in (y, z) plane for Bx field
// Returns the momentum in GeV/c
// fitting the (Y, Z) plane for Bx field.
TrackAnalysis::CircleFitResult TrackAnalysis::CircleFitYZ_Kasa(const std::vector<Hit>& hits, double bx_field) {
    CircleFitResult result = {0, 0, 0, 0, 0, 0, 0., false};
    if (hits.size() < 3) return result;

    // Kåsa method sums
    double sum_y = 0, sum_z = 0, sum_yy = 0, sum_zz = 0, sum_yz = 0, sum_r2 = 0;
    double sum_y_r2 = 0, sum_z_r2 = 0;
    int n = 0;

    for (const auto& hit : hits) {
        double y = hit.position.Y() * 1e-3; // mm to m
        double z = hit.position.Z() * 1e-3; // mm to m
        double r2 = y*y + z*z;
        sum_y += y;
        sum_z += z;
        sum_yy += y*y;
        sum_zz += z*z;
        sum_yz += y*z;
        sum_r2 += r2;
        sum_y_r2 += y*r2;
        sum_z_r2 += z*r2;
        n++;
    }

    double A = n * sum_yy - sum_y * sum_y;
    double B = n * sum_yz - sum_y * sum_z;
    double C = n * sum_zz - sum_z * sum_z;
    double D = 0.5 * (n * sum_y_r2 - sum_y * sum_r2);
    double E = 0.5 * (n * sum_z_r2 - sum_z * sum_r2);

    double det = A*C - B*B;
    if (std::abs(det) < 1e-12) {
        std::cout << "[CircleFit] Singular matrix, cannot fit circle." << std::endl;
        return result;
    }

    double yc = (D*C - B*E) / det;
    double zc = (A*E - B*D) / det;

    // Compute radius as average distance to center
    double sum_r = 0;
    std::vector<double> residuals;
    for (const auto& hit : hits) {
        double dy = hit.position.Y()* 1e-3 - yc;
        double dz = hit.position.Z()* 1e-3 - zc;
        double r = std::sqrt(dy*dy + dz*dz);
        sum_r += r;
        residuals.push_back(r);
    }
    double R = sum_r / n; // in m

    // RMS error of radius
    double s = 0;
    for (auto r : residuals) s += (r - R) * (r - R);
    double sigma = std::sqrt(s / n); // in meters
    int ndf = n - 3;

    // Calculate momentum and its error
    double p = 0.3 * std::abs(bx_field) * R; // GeV/c
    double p_err = 0.3 * std::abs(bx_field) * sigma;

    // Determine charge from bending direction
    int charge = 0;
    if (hits.size() >= 3) {
        TVector3 p1 = hits.front().position;
        TVector3 p2 = hits[hits.size()/2].position;
        TVector3 p3 = hits.back().position;
        
        double sagitta_y = p2.Y() - 0.5 * (p1.Y() + p3.Y());
        
        // Get local field with sign at middle point
        //double local_bx = GetFieldX(p2.X(), p2.Y(), p2.Z());
        double local_bx = GetFieldBX(p2.Y()); // Get field with sign
        if (std::abs(local_bx) > 1e-6) {
            charge = (sagitta_y * local_bx > 0) ? -1 : +1;
        }
    }
    //std::cout << "[CircleFitYZ_Kasa] yc=" << yc << ", zc=" << zc
    //          << ", R=" << R << ", sigma=" << sigma
    //          << ", ndf=" << ndf << ", p=" << p << ", p_err=" << p_err
    //          << ", charge=" << charge << std::endl;
    result = {yc, zc, R, sigma, ndf, p, p_err, charge, true};

    return result;
}
///////////////////////////////////////////////////////////// 
/////CircleFitYZ_Kasa method implementation on Tracklet//////
/////////////////////////////////////////////////////////////
// Circle fit in (y, z) plane for Bx field
// Returns the momentum in GeV/c
// fitting the (Y, Z) plane for Bx field.
TrackAnalysis::CircleFitResult TrackAnalysis::CircleFitYZ_Kasa(const std::vector<Tracklet>& tracklets, double bx_field) {
    CircleFitResult result = {0, 0, 0, 0, 0, 0, 0, 0, false};
    if (tracklets.size() < 3) return result;

    // Kåsa method sums
    double sum_y = 0, sum_z = 0, sum_yy = 0, sum_zz = 0, sum_yz = 0, sum_r2 = 0;
    double sum_y_r2 = 0, sum_z_r2 = 0;
    int n = 0;

    for (const auto& hit : tracklets) {
        double y = hit.position.Y() * 1e-3; // mm to m
        double z = hit.position.Z() * 1e-3; // mm to m
        double r2 = y*y + z*z;
        sum_y += y;
        sum_z += z;
        sum_yy += y*y;
        sum_zz += z*z;
        sum_yz += y*z;
        sum_r2 += r2;
        sum_y_r2 += y*r2;
        sum_z_r2 += z*r2;
        n++;
    }

    double A = n * sum_yy - sum_y * sum_y;
    double B = n * sum_yz - sum_y * sum_z;
    double C = n * sum_zz - sum_z * sum_z;
    double D = 0.5 * (n * sum_y_r2 - sum_y * sum_r2);
    double E = 0.5 * (n * sum_z_r2 - sum_z * sum_r2);

    double det = A*C - B*B;
    if (std::abs(det) < 1e-12) {
        std::cout << "[CircleFit] Singular matrix, cannot fit circle." << std::endl;
        return result;
    }

    double yc = (D*C - B*E) / det;
    double zc = (A*E - B*D) / det;

    // Compute radius as average distance to center
    double sum_r = 0;
    std::vector<double> residuals;
    for (const auto& hit : tracklets) {
        double dy = hit.position.Y()* 1e-3 - yc;
        double dz = hit.position.Z()* 1e-3 - zc;
        double r = std::sqrt(dy*dy + dz*dz);
        sum_r += r;
        residuals.push_back(r);
    }
    double R = sum_r / n; // in mm

    // RMS error of radius
    double s = 0;
    for (auto r : residuals) s += (r - R) * (r - R);
    double sigma = std::sqrt(s / n); // in meters
    int ndf = n - 3;

    // Calculate momentum and its error
    double p = 0.3 * std::abs(bx_field) * R; // GeV/c
    double p_err = 0.3 * std::abs(bx_field) * sigma;

    int charge = 0;
    if (tracklets.size() >= 3) {
        // Use first, middle, last points to determine bending direction
        TVector3 p1 = tracklets.front().position;
        TVector3 p2 = tracklets[tracklets.size()/2].position;
        TVector3 p3 = tracklets.back().position;
        
        // Calculate sagitta in Y direction (for Bx field)
        double sagitta_y = p2.Y() - 0.5 * (p1.Y() + p3.Y());
        
        // Get local field with sign at middle point
        //double local_bx = GetFieldX(p2.X(), p2.Y(), p2.Z());
        double local_bx = GetFieldBX(p2.Y()); // Get field with sign

        if (std::abs(local_bx) > 1e-6) {
            charge = (sagitta_y * local_bx > 0) ? -1 : +1;
        }
    }
    
    //results_.back().kasaCharge_trk = charge;
    //results_.back().kasaResolution_trk = p_err;
    
    result = {yc, zc, R, sigma, ndf, p, p_err, charge, true};
    return result;
}
///////////////////////////////////////////////////////// 
///// Linear Fit function                       /////////
/////////////////////////////////////////////////////////
// -- Fit a line to N points in (z,x) or (z,y) for direction and intercept
void TrackAnalysis::LinearFit(const std::vector<TVector3>& pts, double& slope, double& intercept, char coord) {
    // coord: 'x' or 'y'
    double sumz = 0, sumz2 = 0, sumx = 0, sumzx = 0;
    int N = pts.size();
    for (const auto& p : pts) {
        double z = p.Z();
        double x = (coord == 'x') ? p.X() : p.Y();
        sumz  += z;
        sumz2 += z*z;
        sumx  += x;
        sumzx += z*x;
    }
    double denom = N * sumz2 - sumz * sumz;
    if (denom == 0) { slope = 0; intercept = 0; return; }
    slope = (N * sumzx - sumz * sumx) / denom;
    intercept = (sumx - slope * sumz) / N;
}
///////////////////////////////////////////////////////// 
/////MomentumFromSlopeChange method implementation///////
/////////////////////////////////////////////////////////
double TrackAnalysis::MomentumFromSlopeChange(
    const std::vector<TVector3>& avgs,
    std::function<double(double, double, double)> GetField,
    const std::vector<std::pair<double, double>>& magnet_z_ranges){
    // Find indices before and after the field region(s)
    int i_first = -1, i_last = -1;
    double z_first = 0, z_last = 0;
    for (size_t i = 0; i < avgs.size(); ++i) {
        double z = avgs[i].Z();
        for (const auto& rng : magnet_z_ranges) {
            if (z > rng.first && z < rng.second) {
                if (i_first < 0 && i > 0) { i_first = i-1; z_first = avgs[i-1].Z(); }
                i_last = i; z_last = z;
            }
        }
    }
    if (i_first < 0 || i_last < 0 || i_first == i_last) return 0.0; // Not enough points

    // Collect points before and after field region
    std::vector<TVector3> before, after;
    for (size_t i = 0; i <= i_first; ++i) before.push_back(avgs[i]);
    for (size_t i = i_last; i < avgs.size(); ++i) after.push_back(avgs[i]);

    // Fit before and after
    double slope_xz1, intcpt_xz1, slope_yz1, intcpt_yz1;
    double slope_xz2, intcpt_xz2, slope_yz2, intcpt_yz2;
    LinearFit(before, slope_xz1, intcpt_xz1, 'x');
    LinearFit(before, slope_yz1, intcpt_yz1, 'y');
    LinearFit(after,  slope_xz2, intcpt_xz2, 'x');
    LinearFit(after,  slope_yz2, intcpt_yz2, 'y');
    TVector3 dir1(slope_xz1, slope_yz1, 1.0); dir1 = dir1.Unit();
    TVector3 dir2(slope_xz2, slope_yz2, 1.0); dir2 = dir2.Unit();
    double dtheta = dir2.Angle(dir1); // radians

    // Estimate B at center of magnet (or average if needed)
    double z_center = 0.5*(z_first + z_last);
    double B = GetField(0, 0, z_center); // assumes field along x

    double L = (z_last - z_first) * 0.001; // mm -> m

    if (std::abs(dtheta) < 1e-7) return 0.0;
    return 0.3 * std::abs(B) * L / dtheta; // GeV/c
}
///////////////////////////////////////////////////////// 
/////ReconstructGlobalTrack method implementation.///////
/////////////////////////////////////////////////////////
void TrackAnalysis::ReconstructGlobalTrack(const std::vector<Tracklet>& tracklets)
{
    // Group tracklets by (trackID, pdg)
    std::map<std::pair<int,int>, std::vector<Tracklet>> tkmap;
    for (const auto& t : tracklets)
        tkmap[{t.trackID, t.pdg}].push_back(t);

    for (const auto& entry : tkmap) {
        int trackID = entry.first.first;
        int pdg     = entry.first.second;
        const auto& vec = entry.second;

        // Skip tracks with too few stations/tracklets
        if (vec.size() < 3) continue;

        // Collect avg positions of all tracklets
        std::vector<TVector3> avgs;
        for (const auto& t : vec) avgs.push_back(t.position);

        // Fit lines in Z-X and Z-Y to get the overall direction
        double slope_xz, intcpt_xz, slope_yz, intcpt_yz;
        LinearFit(avgs, slope_xz, intcpt_xz, 'x');
        LinearFit(avgs, slope_yz, intcpt_yz, 'y');

        // Build a direction vector (not normalized)
        TVector3 dir(slope_xz, slope_yz, 1.0); 
        dir = dir.Unit();

        // Calculate the average momentum using the circle fit in (Z,Y) plane
        // This is a simplified approach, assuming the track is in the Bx field
        // and using the average of the tracklet positions as points.
        // Sagitta/circle fit for momentum if in field (use avgs as points!)
        // For Bx field, fit (Z,Y), for By, fit (Z,X) -- here for Bx:
        double momentum = 0.0;
        if (avgs.size() >= 3) {
            // Simple circle fit in (Z,Y)
            // 1. Take first, middle, last
            TVector3 p0 = avgs.front();
            TVector3 pm = avgs[avgs.size()/2];
            TVector3 p2 = avgs.back();
            
            // Lever arm, L_total, total path length inside the magnetic field(s) traversed by the track.
            // Find the path length INSIDE all magnet regions
            double L_total = 0.0;
            for (const auto& rng : magnet_z_ranges_) {
                // Only consider the part of the track between p0.Z() and p2.Z()
                double z_in  = std::max(rng.first, std::min(p0.Z(), p2.Z()));
                double z_out = std::min(rng.second, std::max(p0.Z(), p2.Z()));
                if (z_out > z_in) {
                    L_total += (z_out - z_in);
                }
            }
            L_total *= 0.01; // cm → m
            // Calculate sagitta
            // Sagitta is the distance from the midpoint of the chord to the circle center
            // Sagitta = (Y_mid - Y_center) * 0.001; //
            // where Y_mid is the midpoint of the chord (p0.Y() + p2.Y()) / 2
            // and Y_center is the Y coordinate of the circle center, which we can estimate as
            // Y_center = Y_mid + sagitta
            // Here we use the midpoint of the chord (p0, p2) and the middle point pm
            // to calculate the sagitta and then the momentum.
            double sagitta = (pm.Y() - 0.5 * (p0.Y() + p2.Y())) * 0.001; // mm -> m
            double B = 1.5; // Tesla
            momentum = (std::abs(sagitta) > 1e-8) ? (0.3 * B * L_total * L_total / (8. * std::abs(sagitta))) : 0.0;
            
            std::cout << "[RECO] TrackID: " << trackID
                  << ", PDG: " << pdg
                  << ", Nstation=" << vec.size()
                  << ", Direction: (" << dir.X() << "," << dir.Y() << "," << dir.Z() << ")"
                  << ", Lever arm L: " << L_total << " m"
                  << ", Sagitta: " << sagitta << " m"
                  << ", Momentum: " << momentum << " GeV/c"
                  << std::endl;
       
            // Build a list of "pseudo-hits" from tracklet positions
            std::vector<Hit> pseudoHits;
            for (const auto& trk : vec) {
                pseudoHits.push_back({
                    .stationID = trk.station,
                    .layerID   = -1, // no specific layer
                    .trackID   = trk.trackID,
                    .pdg       = trk.pdg,
                    .position  = trk.position,
                    .momentum  = TVector3() // leave blank; can fill with dummy if needed
                });
            }
            /////// Let's use Kalman filter here to refine the momentum estimate ///////
            //results_.emplace_back(); // Add a result entry
            //RunKalmanFilter(pseudoHits);
            //RunGenfit(pseudoHits);
            // Print results
            //const auto& r = results_.back();
           // std::cout << "[RECO][KF/GF] TrackID: " << trackID
            //      << ", Stations: " << vec.size()
            //      << ", Kalman P: " << r.kalmanMomentum << " GeV/c"
            //      << ", Kalman chi2: " << r.kalmanChi2
            //      << ", GENFIT P: " << r.genfitMomentum << " GeV/c"
            //      << ", GENFIT chi2: " << r.genfitChi2
            //      << ", GENFIT charge: " << r.genfitCharge
            //      << std::endl;
       
        }
        //double p_slope = MomentumFromSlopeChange(
        //avgs, [this](double x, double y, double z){ return this->GetFieldX(x, y, z); }, magnet_z_ranges_);
        //std::cout << "[SLOPE] Angle method: momentum = " << p_slope << " GeV/c\n";  
    }
}
///////////////////////////////////////////////////////// 
/////   WriteResults method implementation.       /////// 
/////////////////////////////////////////////////////////
void TrackAnalysis::WriteResults(const char* filename) {
    TFile outfile(filename, "RECREATE");
    TTree tree("results", "Reconstruction Results");

    double truth, kalman, genfit;
    double sagittaMomentum, sagittaMomentum_trk;
    double kasaMomentum, kasaMomentumTrk;
    double taubinMomentum, taubinMomentumTrk;
    double trackletMomentum, trackletMomentumGlobal;


    double kalmanMomentumErr, kalmanChi2;
    double genfitChi2, genfitMomentumErr, genfitCharge;
    int event_, trackid_, pdg_;

    // ADD: New charge variables
    int sagittaCharge, kasaCharge, taubinCharge, trackletCharge, trackletChargeGlobal, kalmanCharge;
    int sagittaCharge_trk, kasaCharge_trk, taubinCharge_trk;

    // ADD: New resolution variables
    double sagittaRes, kasaRes, taubinRes, trackletRes, trackletResGlobal;
    double sagittaRes_trk, kasaRes_trk, taubinRes_trk;

    // Branches for each result
    tree.Branch("event", &event_, "event/I");
    tree.Branch("trackid", &trackid_, "trackid/I");
    tree.Branch("pdg", &pdg_, "pdg/I");
    tree.Branch("truth", &truth, "truth/D");
    // Sagitta momentum and charge branches
    tree.Branch("sagitta", &sagittaMomentum, "sagitta/D");
    tree.Branch("sagitta_trk", &sagittaMomentum_trk, "sagitta_trk/D");
    tree.Branch("sagitta_charge", &sagittaCharge, "sagitta_charge/I");
    tree.Branch("sagitta_charge_trk", &sagittaCharge_trk, "sagitta_charge_trk/I");
    tree.Branch("sagitta_res", &sagittaRes, "sagitta_res/D");
    tree.Branch("sagitta_res_trk", &sagittaRes_trk, "sagitta_res_trk/D");

    // Kasa and Taubin momentum branches
    tree.Branch("kasa", &kasaMomentum, "kasa/D");
    tree.Branch("kasa_trk", &kasaMomentumTrk, "kasa_trk/D");
    tree.Branch("kasa_charge", &kasaCharge, "kasa_charge/I");
    tree.Branch("kasa_charge_trk", &kasaCharge_trk, "kasa_charge_trk/I");
    tree.Branch("kasa_res", &kasaRes, "kasa_res/D");
    tree.Branch("kasa_res_trk", &kasaRes_trk, "kasa_res_trk/D");


    tree.Branch("taubin", &taubinMomentum, "taubin/D");
    tree.Branch("taubin_trk", &taubinMomentumTrk, "taubin_trk/D");
    tree.Branch("taubin_charge", &taubinCharge, "taubin_charge/I");
    tree.Branch("taubin_charge_trk", &taubinCharge_trk, "taubin_charge_trk/I");
    tree.Branch("taubin_res", &taubinRes, "taubin_res/D");
    tree.Branch("taubin_res_trk", &taubinRes_trk, "taubin_res_trk/D");

    tree.Branch("tracklet", &trackletMomentum, "tracklet/D");
    tree.Branch("trackletGlobal", &trackletMomentumGlobal, "trackletGlobal/D");
    tree.Branch("tracklet_charge", &trackletCharge, "tracklet_charge/I");
    tree.Branch("tracklet_charge_global", &trackletChargeGlobal, "tracklet_charge_global/I");
    tree.Branch("tracklet_res", &trackletRes, "tracklet_res/D");
    tree.Branch("tracklet_res_global", &trackletResGlobal, "tracklet_res_global/D");

    tree.Branch("kalman", &kalman, "kalman/D");
    tree.Branch("kalman_err", &kalmanMomentumErr, "kalman_err/D");
    tree.Branch("kalman_chi2", &kalmanChi2, "kalman_chi2/D");
    tree.Branch("kalman_charge", &kalmanCharge, "kalman_charge/I");

    // Branch for genfit momentum and chi2
    tree.Branch("genfit", &genfit, "genfit/D");
    tree.Branch("genfit_chi2", &genfitChi2, "genfit_chi2/D");
    tree.Branch("genfit_err", &genfitMomentumErr, "genfit_err/D");
    tree.Branch("genfit_charge", &genfitCharge, "genfit_charge/D");

    
    for (const auto& r : results_) {
        event_ = r.eventID;
        trackid_ = r.trackID;
        pdg_ = r.pdg;
        truth = r.truthMomentum;

        sagittaMomentum = r.sagittaMomentum;
        sagittaMomentum_trk = r.sagittaMomentum_trk;
        sagittaCharge = r.sagittaCharge;
        sagittaCharge_trk = r.sagittaCharge_trk;
        sagittaRes = r.sagittaResolution;
        sagittaRes_trk = r.sagittaResolution_trk;

        kalman = r.kalmanMomentum;
        kalmanMomentumErr = r.kalmanMomentumErr;
        kalmanChi2 = r.kalmanChi2;
        kalmanCharge = r.kalmanCharge;

        trackletMomentum = r.trackletDeflectionMomentum;
        trackletMomentumGlobal = r.trackletDeflectionMomentum_trk;
        trackletCharge = r.trackletDeflectionCharge;
        trackletChargeGlobal = r.trackletDeflectionCharge_trk;
        trackletRes = r.trackletDeflectionMomentumErr;
        trackletResGlobal = r.trackletDeflectionMomentumErr_trk;

        kasaMomentum = r.kasaMomentum;
        kasaMomentumTrk = r.kasaMomentum_trk;
        kasaCharge = r.kasaCharge;
        kasaCharge_trk = r.kasaCharge_trk;
        kasaRes = r.kasaResolution;
        kasaRes_trk = r.kasaResolution_trk;
        
        taubinMomentum = r.taubinMomentum;
        taubinMomentumTrk = r.taubinMomentum_trk;
        taubinCharge = r.taubinCharge;
        taubinCharge_trk = r.taubinCharge_trk;
        taubinRes = r.taubinResolution;
        taubinRes_trk = r.taubinResolution_trk;


        genfit = r.genfitMomentum;
        genfitChi2 = r.genfitChi2;
        genfitCharge = r.genfitCharge;
        genfitMomentumErr = r.genfitMomentumErr;

        tree.Fill();
    }

    tree.Write();
    outfile.Close();
}
///////////////////////////////////////////////////////// 
/////   PrintOut method implementation.           /////// 
/////////////////////////////////////////////////////////
void TrackAnalysis::PrintOut(const std::vector<Hit> &hits)
{
   std::cout << "Processing event with hits: " << std::endl;
    for (size_t i = 0; i < hits.size(); ++i) 
    {
        std::cout << i << " "
          << hits[i].pdg << " "
          << hits[i].stationID << " "
          << hits[i].layerID << " " 
          << hits[i].position.X() << " " 
          << hits[i].position.Y() << " " 
          << hits[i].position.Z() << " " 
          << hits[i].momentum.X() << " " 
          << hits[i].momentum.Y() << " " 
          << hits[i].momentum.Z() << " " 
          << std::endl;
    } 

    std::set<std::pair<int,int>> layers_seen;
    for (const auto& hit : hits) {
        std::cout << "station=" << hit.stationID << " layer=" << hit.layerID
                << " trackID=" << hit.trackID << " pdg=" << hit.pdg << std::endl;
        layers_seen.insert({hit.stationID, hit.layerID});
    }
    std::cout << "Unique (station,layer) pairs: " << layers_seen.size() << std::endl;

    std::map<std::pair<int,int>, int> layer_counts;
    for (const auto& hit : hits) {
        if (hit.trackID == 1 && std::abs(hit.pdg)==13) {
            auto key = std::make_pair(hit.stationID, hit.layerID);
            layer_counts[key]++;
        }
    }
    for (const auto& entry : layer_counts) {
        if (entry.second > 1)
            std::cout << "Duplicate hit: station=" << entry.first.first
                    << ", layer=" << entry.first.second
                    << ", count=" << entry.second << std::endl;
    }

}
///////////////////////////////////////////////////////// 
/////   PrintStationFits method implementation.   /////// 
/////////////////////////////////////////////////////////
void TrackAnalysis::PrintStationFits(const std::vector<Hit>& hits) {
    // Group hits by stationID
    std::map<int, std::vector<TVector3>> stationHits;
    for (const auto& hit : hits) {
        stationHits[hit.stationID].push_back(hit.position);
    }

    std::cout << "---- Station Line Fits ----" << std::endl;
    for (const auto& entry : stationHits) {
        int station = entry.first;
        const auto& vecs = entry.second;
        if (vecs.size() < 2) continue; // Need at least 2 points

        // Linear fit: x = a + b*z, y = c + d*z
        double sumZ = 0, sumX = 0, sumY = 0, sumZZ = 0, sumXZ = 0, sumYZ = 0;
        for (const auto& pos : vecs) {
            double x = pos.X(), y = pos.Y(), z = pos.Z();
            sumZ  += z;
            sumX  += x;
            sumY  += y;
            sumZZ += z * z;
            sumXZ += x * z;
            sumYZ += y * z;
        }
        double N = vecs.size();
        double denom = N * sumZZ - sumZ * sumZ;
        double slopeX = (N * sumXZ - sumZ * sumX) / denom;
        double slopeY = (N * sumYZ - sumZ * sumY) / denom;
        double interceptX = (sumX - slopeX * sumZ) / N;
        double interceptY = (sumY - slopeY * sumZ) / N;

        // Direction vector (dx/dz, dy/dz, 1)
        TVector3 dir(slopeX, slopeY, 1.0);
        dir = dir.Unit(); // Normalize

        std::cout << "Station " << station
                  << ": direction = (" << dir.X() << ", " << dir.Y() << ", " << dir.Z() << ")" << std::endl;

        // (Optional) Print points used in fit
        for (size_t i = 0; i < vecs.size(); ++i) {
            std::cout << "    hit " << i << ": ("
                      << vecs[i].X() << ", " << vecs[i].Y() << ", " << vecs[i].Z() << ")" << std::endl;
        }
    }
}

///////////////////////////////////////////////////////// 
/////   PlotHitsAndTrack method implementation.   /////// 
/////////////////////////////////////////////////////////
void TrackAnalysis::PlotHitsAndTrack(const std::vector<Hit>& hits,
                                     const std::vector<TVector3>& fittedTrackPoints,
                                     const std::string& title)
{
    std::string canvasName = "c_" + title;

    TCanvas* c1 = new TCanvas(canvasName.c_str(), title.c_str(), 800, 600);

    // Plot hits as red points
    TGraph2D* gHits = new TGraph2D(hits.size());
    gHits->SetName(("Graph2D_" + std::to_string(hits.front().trackID)).c_str());
    gHits->SetTitle(title.c_str());

    for (size_t i = 0; i < hits.size(); ++i)
        gHits->SetPoint(i, hits[i].position.X(), hits[i].position.Y(), hits[i].position.Z());
    gHits->SetMarkerStyle(20);
    gHits->SetMarkerSize(1.2);
    gHits->SetMarkerColor(kRed);
    gHits->Draw("P");

    // --- Axis labels ---
    gHits->GetXaxis()->SetTitle("X [cm]");
    gHits->GetYaxis()->SetTitle("Y [cm]");
    gHits->GetZaxis()->SetTitle("Z [cm]");
    gHits->GetXaxis()->SetTitleOffset(1.2);
    gHits->GetYaxis()->SetTitleOffset(1.2);
    gHits->GetZaxis()->SetTitleOffset(1.2);

    c1->SetGrid();

    // Plot fitted track if available
    if (!fittedTrackPoints.empty()) {
        TPolyLine3D* pl = new TPolyLine3D(fittedTrackPoints.size());
        for (size_t i = 0; i < fittedTrackPoints.size(); ++i) {
            pl->SetPoint(i,
                fittedTrackPoints[i].X(),
                fittedTrackPoints[i].Y(),
                fittedTrackPoints[i].Z());
        }
        pl->SetLineColor(kBlue);
        pl->SetLineWidth(3);
        pl->Draw("same");
    }
    // Legend
    TLegend* leg = new TLegend(0.15,0.8,0.5,0.93);
    leg->AddEntry(gHits, "Measured hits", "p");
    leg->AddEntry((TObject*)0, "Fitted track", "l");
    leg->Draw();

    c1->Modified();
    c1->Update();

    // Save canvas as PDF, file name based on title/track ID
    std::string pdfName = title + ".pdf";
    c1->SaveAs(pdfName.c_str());
}
///////////////////////////////////////////////////////// 
/////   PlotMCViews method implementation.        /////// 
/////////////////////////////////////////////////////////
void TrackAnalysis::PlotMCViews(const std::vector<Hit>& hits, const std::string& basename)
{
    // Group by PDG
    std::vector<double> x_mu, y_mu, z_mu;
    std::vector<double> x_e, y_e, z_e;
    std::vector<double> x_oth, y_oth, z_oth;

    for (const auto& hit : hits) {
        if (std::abs(hit.pdg) == 13) { // muon
            x_mu.push_back(hit.position.X());
            y_mu.push_back(hit.position.Y());
            z_mu.push_back(hit.position.Z());
        } else if (std::abs(hit.pdg) == 11) { // electron/positron
            x_e.push_back(hit.position.X());
            y_e.push_back(hit.position.Y());
            z_e.push_back(hit.position.Z());
        } else { // other
            x_oth.push_back(hit.position.X());
            y_oth.push_back(hit.position.Y());
            z_oth.push_back(hit.position.Z());
        }
    }

    auto plot_proj = [&](const std::string& title, 
                         const std::string& xlab, const std::string& ylab,
                         const std::vector<double>& x1, const std::vector<double>& y1, Color_t c1,
                         const std::vector<double>& x2, const std::vector<double>& y2, Color_t c2,
                         const std::vector<double>& x3, const std::vector<double>& y3, Color_t c3,
                         const std::string& outname)
    {
        TCanvas* c = new TCanvas((basename+"_"+title).c_str(), (basename+"_"+title).c_str(), 900, 700);

        TGraph* gr_mu = nullptr;
        TGraph* gr_e = nullptr;
        TGraph* gr_oth = nullptr;

        if (!x1.empty()) {
            gr_mu = new TGraph(x1.size(), x1.data(), y1.data());
            gr_mu->SetMarkerStyle(20); gr_mu->SetMarkerSize(1.4); gr_mu->SetMarkerColor(c1);
            gr_mu->SetTitle((title+";"+xlab+";"+ylab).c_str());
            gr_mu->Draw("AP");
        }
        if (!x2.empty()) {
            gr_e = new TGraph(x2.size(), x2.data(), y2.data());
            gr_e->SetMarkerStyle(21); gr_e->SetMarkerSize(1.3); gr_e->SetMarkerColor(c2);
            if (!gr_mu) gr_e->SetTitle((title+";"+xlab+";"+ylab).c_str());
            gr_e->Draw(gr_mu ? "P SAME" : "AP");
        }
        if (!x3.empty()) {
            gr_oth = new TGraph(x3.size(), x3.data(), y3.data());
            gr_oth->SetMarkerStyle(22); gr_oth->SetMarkerSize(1.1); gr_oth->SetMarkerColor(c3);
            if (!gr_mu && !gr_e) gr_oth->SetTitle((title+";"+xlab+";"+ylab).c_str());
            gr_oth->Draw((gr_mu || gr_e) ? "P SAME" : "AP");
        }

        // Legend creation (only for non-empty graphs)
        TLegend* leg = new TLegend(0.13, 0.81, 0.43, 0.93);
        if (gr_mu) leg->AddEntry(gr_mu, "#mu^{#pm}", "p");
        if (gr_e)  leg->AddEntry(gr_e, "e^{#pm}", "p");
        if (gr_oth) leg->AddEntry(gr_oth, "other", "p");
        leg->Draw();
        c->SetGrid();

        // Save PDF
        c->SaveAs(outname.c_str());

        // Interactive pause 
        c->WaitPrimitive();
        if (gr_mu) delete gr_mu;
        if (gr_e) delete gr_e;
        if (gr_oth) delete gr_oth;
        delete leg;
        delete c;

    };

    // XY (standard, X horizontal, Y vertical)
    plot_proj("XY", "X [mm]", "Y [mm]",
        x_mu, y_mu, kGreen+2, x_e, y_e, kRed, x_oth, y_oth, kGray+2,
        basename+"_XY.pdf"
    );

    // XZ (Z on X axis, X on Y axis)
    plot_proj("XZ", "Z [mm]", "X [mm]",
        z_mu, x_mu, kGreen+2, z_e, x_e, kRed, z_oth, x_oth, kGray+2,
        basename+"_XZ.pdf"
    );

    // YZ (Z on X axis, Y on Y axis)
    plot_proj("YZ", "Z [mm]", "Y [mm]",
        z_mu, y_mu, kGreen+2, z_e, y_e, kRed, z_oth, y_oth, kGray+2,
        basename+"_YZ.pdf"
    );
}
///////////////////////////////////////////////////////// 
/////PlotHitsAndTrackLets method implementation.  /////// 
/////////////////////////////////////////////////////////
void TrackAnalysis::PlotHitsAndTracklets(const std::vector<Hit>& hits, const std::vector<Tracklet>& tracklets, const std::string& title)
{
    TCanvas* c3d = new TCanvas(("c3d_"+title).c_str(), ("3D "+title).c_str(), 1000, 800);
    TGraph2D* gHits = new TGraph2D(hits.size());
    for (size_t i = 0; i < hits.size(); ++i)
        gHits->SetPoint(i, hits[i].position.X(), hits[i].position.Y(), hits[i].position.Z());
    gHits->SetMarkerStyle(20);
    gHits->SetMarkerSize(1.2);
    gHits->SetMarkerColor(kRed+1);
    gHits->SetTitle((title+";X [mm];Y [mm];Z [mm]").c_str());
    gHits->Draw("P");
    c3d->SetGrid();

    // Draw station tracklets
    double L = 100.0; // mm, length of the drawn tracklet
    for (const auto& trk : tracklets) {
        TVector3 p1 = trk.position - 0.5*L*trk.dir;
        TVector3 p2 = trk.position + 0.5*L*trk.dir;
        TPolyLine3D* line = new TPolyLine3D(2);
        line->SetPoint(0, p1.X(), p1.Y(), p1.Z());
        line->SetPoint(1, p2.X(), p2.Y(), p2.Z());
        line->SetLineColor(kViolet+1);
        line->SetLineWidth(4);
        line->Draw("same");
    }
    c3d->Modified(); c3d->Update();
    c3d->SaveAs((title+"_3DTracklets.pdf").c_str());
}
///////////////////////////////////////////////////////// 
/////PlotTrackletsViews method implementation.  /////// 
/////////////////////////////////////////////////////////
void TrackAnalysis::PlotTrackletsViews(const std::vector<Hit>& hits, const std::vector<Tracklet>& tracklets, const std::string& basename)
{
    double L = 100.0; // mm: visual length of tracklet line
    auto colorByPDG = [](int pdg) -> int {
        if (std::abs(pdg) == 13) return kGreen+2;     // Muon
        if (std::abs(pdg) == 11) return kRed;         // Electron/positron
        return kGray+2;                               // Other
    };

    // ---- XZ projection (Z horizontal, X vertical) ----
    TCanvas* cXZ = new TCanvas((basename+"_XZ").c_str(), (basename+" XZ").c_str(), 900, 700);
    // Hits
    TGraph* gr = new TGraph(hits.size());
    for (size_t i = 0; i < hits.size(); ++i)
        gr->SetPoint(i, hits[i].position.Z(), hits[i].position.X());
    gr->SetMarkerStyle(20); gr->SetMarkerSize(1.1); gr->SetMarkerColor(kBlue);
    gr->SetTitle((basename+";Z [mm];X [mm]").c_str());
    gr->Draw("AP");
    cXZ->SetGrid();

    // Tracklets
    for (const auto& t : tracklets) {
        TVector3 p1 = t.position - 0.5*L*t.dir;
        TVector3 p2 = t.position + 0.5*L*t.dir;
        TLine* line = new TLine(p1.Z(), p1.X(), p2.Z(), p2.X());
        line->SetLineColor(colorByPDG(t.pdg));
        line->SetLineWidth(3);
        line->Draw("same");
    }
    
    // Legend
    TLegend* legXZ = new TLegend(0.13, 0.81, 0.43, 0.93);
    TLine* l_mu = new TLine(0,0,1,0); l_mu->SetLineColor(kGreen+2); l_mu->SetLineWidth(3);
    TLine* l_e = new TLine(0,0,1,0);  l_e->SetLineColor(kRed);      l_e->SetLineWidth(3);
    TLine* l_o = new TLine(0,0,1,0);  l_o->SetLineColor(kGray+2);   l_o->SetLineWidth(3);
    legXZ->AddEntry(l_mu, "#mu^{#pm} tracklet", "l");
    legXZ->AddEntry(l_e, "e^{#pm} tracklet", "l");
    legXZ->AddEntry(l_o, "other", "l");
    legXZ->Draw();

    cXZ->Modified(); cXZ->Update();
    cXZ->SaveAs((basename+"_XZ_tracklets.pdf").c_str());

    // ---- YZ projection (Z horizontal, Y vertical) ----
    TCanvas* cYZ = new TCanvas((basename+"_YZ").c_str(), (basename+" YZ").c_str(), 900, 700);
    TGraph* gr2 = new TGraph(hits.size());
    for (size_t i = 0; i < hits.size(); ++i)
        gr2->SetPoint(i, hits[i].position.Z(), hits[i].position.Y());
    gr2->SetMarkerStyle(20); gr2->SetMarkerSize(1.1); gr2->SetMarkerColor(kBlue);
    gr2->SetTitle((basename+";Z [mm];Y [mm]").c_str());
    gr2->Draw("AP");
    cYZ->SetGrid();

    for (const auto& t : tracklets) {
        TVector3 p1 = t.position - 0.5*L*t.dir;
        TVector3 p2 = t.position + 0.5*L*t.dir;
        TLine* line = new TLine(p1.Z(), p1.Y(), p2.Z(), p2.Y());
        line->SetLineColor(colorByPDG(t.pdg));
        line->SetLineWidth(3);
        line->Draw("same");
    }

    // Legend for YZ
    TLegend* legYZ = new TLegend(0.13, 0.81, 0.43, 0.93);
    TLine* l_mu2 = new TLine(0,0,1,0); l_mu2->SetLineColor(kGreen+2); l_mu2->SetLineWidth(3);
    TLine* l_e2 = new TLine(0,0,1,0);  l_e2->SetLineColor(kRed);      l_e2->SetLineWidth(3);
    TLine* l_o2 = new TLine(0,0,1,0);  l_o2->SetLineColor(kGray+2);   l_o2->SetLineWidth(3);
    legYZ->AddEntry(l_mu2, "#mu tracklet", "l");
    legYZ->AddEntry(l_e2, "e^{#pm} tracklet", "l");
    legYZ->AddEntry(l_o2, "other", "l");
    legYZ->Draw();

    cYZ->Modified(); cYZ->Update();
    cYZ->SaveAs((basename+"_YZ_tracklets.pdf").c_str());
    
    // Clean up 
    //delete gr; delete gr2;
    //delete l_mu; delete l_e; delete l_o; delete legXZ;
    //delete l_mu2; delete l_e2; delete l_o2; delete legYZ;

}

static int colorByPDG(int pdg) {
    if (std::abs(pdg) == 13) return kGreen+2;     // Muon
    if (std::abs(pdg) == 11) return kRed;         // Electron/positron
    return kGray+2;                               // Other
}
static std::string pdgLabel(int pdg) {
    if (std::abs(pdg) == 13) return "#mu^{#pm} track";
    if (std::abs(pdg) == 11) return "e^{#pm} track";
    return "other";
}
/////////////////////////////////////////////////////////////
/////   PlotGlobalTrackletFit method implementation. ///////
/////////////////////////////////////////////////////////////
void TrackAnalysis::PlotGlobalTrackletFit(const std::vector<Tracklet>& tracklets, const std::string& basename)
{
    // Group by (trackID, pdg)
    std::map<std::pair<int,int>, std::vector<Tracklet>> tkmap;
    for (const auto& t : tracklets)
        tkmap[{t.trackID, t.pdg}].push_back(t);

    // Store for legend entries to only add once per particle type
    std::set<int> legendAdded;

    // --- XZ Projection ---
    TCanvas* cXZ = new TCanvas((basename+"_XZ").c_str(), (basename+" XZ").c_str(), 950, 700);
    TMultiGraph* mgXZ = new TMultiGraph();

    // --- YZ Projection ---
    TCanvas* cYZ = new TCanvas((basename+"_YZ").c_str(), (basename+" YZ").c_str(), 950, 700);
    TMultiGraph* mgYZ = new TMultiGraph();

    TLegend* legXZ = new TLegend(0.15, 0.77, 0.45, 0.93);
    TLegend* legYZ = new TLegend(0.15, 0.77, 0.45, 0.93);

    for (const auto& entry : tkmap) {
        int trackID = entry.first.first;
        int pdg     = entry.first.second;
        const auto& vec = entry.second;

        if (vec.size() < 3) continue;

        // -- Plot tracklet average positions
        std::vector<double> z, x, y;
        for (const auto& t : vec) {
            z.push_back(t.position.Z());
            x.push_back(t.position.X());
            y.push_back(t.position.Y());
        }
        TGraph* grXZ = new TGraph(z.size(), z.data(), x.data());
        grXZ->SetMarkerStyle(20);
        grXZ->SetMarkerColor(colorByPDG(pdg));
        grXZ->SetMarkerSize(1.4);
        mgXZ->Add(grXZ);

        TGraph* grYZ = new TGraph(z.size(), z.data(), y.data());
        grYZ->SetMarkerStyle(20);
        grYZ->SetMarkerColor(colorByPDG(pdg));
        grYZ->SetMarkerSize(1.4);
        mgYZ->Add(grYZ);

        if (legendAdded.count(pdg)==0) {
            legXZ->AddEntry(grXZ, pdgLabel(pdg).c_str(), "p");
            legYZ->AddEntry(grYZ, pdgLabel(pdg).c_str(), "p");
            legendAdded.insert(pdg);
        }

        // ---- Fit and draw global line ----
        double slope_xz, intcpt_xz, slope_yz, intcpt_yz;
        std::vector<TVector3> avgs;
        for (const auto& t : vec) avgs.push_back(t.position);
        LinearFit(avgs, slope_xz, intcpt_xz, 'x');
        LinearFit(avgs, slope_yz, intcpt_yz, 'y');

        // Draw line from min to max Z
        double minZ = *std::min_element(z.begin(), z.end());
        double maxZ = *std::max_element(z.begin(), z.end());
        double minX = slope_xz * minZ + intcpt_xz, maxX = slope_xz * maxZ + intcpt_xz;
        double minY = slope_yz * minZ + intcpt_yz, maxY = slope_yz * maxZ + intcpt_yz;
        TLine* lineXZ = new TLine(minZ, minX, maxZ, maxX);
        TLine* lineYZ = new TLine(minZ, minY, maxZ, maxY);
        lineXZ->SetLineColor(colorByPDG(pdg));
        lineXZ->SetLineWidth(2);
        lineXZ->SetLineStyle(2);
        lineYZ->SetLineColor(colorByPDG(pdg));
        lineYZ->SetLineWidth(2);
        lineYZ->SetLineStyle(2);

        // Draw after multigraph
        cXZ->cd(); lineXZ->Draw();
        cYZ->cd(); lineYZ->Draw();

        // ---- Momentum estimate (sagitta/circle) using avg positions ----
        double momentum = 0;
        if (avgs.size() >= 3) {
            const TVector3& p0 = avgs.front();
            const TVector3& pm = avgs[avgs.size()/2];
            const TVector3& p2 = avgs.back();
            double L = 1;//(p2.Z() - p0.Z()) * 0.001;  // mm -> m
            double sagitta = (pm.Y() - 0.5 * (p0.Y() + p2.Y())) * 0.001; // mm -> m
            double B = std::abs(MagnetFieldStrength); //1.5; // Tesla
            momentum = (std::abs(sagitta) > 1e-8) ? (0.3 * B * L * L / (8. * std::abs(sagitta))) : 0.0;
        }
        // Print track info
        std::cout << "[RECO] TrackID: " << trackID
                  << ", PDG: " << pdg
                  << ", Nstations: " << vec.size()
                  << ", Dir: (" << slope_xz << "," << slope_yz << ",1)"
                  << ", Momentum: " << momentum << " GeV/c"
                  << std::endl;
    }

    cXZ->cd();
    mgXZ->SetTitle((basename+";Z [mm];X [mm]").c_str());
    mgXZ->Draw("AP");
    legXZ->Draw();
    cXZ->SetGrid();
    cXZ->Modified(); cXZ->Update();
    cXZ->SaveAs((basename+"_XZ_globaltrack.pdf").c_str());

    cYZ->cd();
    mgYZ->SetTitle((basename+";Z [mm];Y [mm]").c_str());
    mgYZ->Draw("AP");
    legYZ->Draw();
    cYZ->SetGrid();
    cYZ->Modified(); cYZ->Update();
    cYZ->SaveAs((basename+"_YZ_globaltrack.pdf").c_str());
}
///////////////////////////////////////////////////////// 
///RunKalmanFilterTracklets method implementation /////// 
/////////////////////////////////////////////////////////
void TrackAnalysis::RunKalmanFilterTracklets(const std::vector<Tracklet>& hits, double seed_p) {

    if (hits.size() < 2) {
        results_.back().kalmanMomentum = 0;
        results_.back().kalmanChi2 = 9999;
        return;
    }

    // Debug print
    if(verbose>=4){
        std::cout << "[DEBUG] Tracklets for Kalman: " << hits.size() << std::endl;
        for (size_t i = 0; i < hits.size(); ++i)
            std::cout << "  Tracklet " << i << ": Z = " << hits[i].position.Z() << std::endl;
    }

    const int nStates = 5; // x, y, z, tx, ty
    const int nMeas = 3;   // x, y, z measurements
    
    TMatrixD state(nStates, 1);
    TMatrixDSym cov(nStates);

    // seed state from first tracklet
    state(0,0) = hits[0].position.X();
    state(1,0) = hits[0].position.Y();
    state(2,0) = hits[0].position.Z();
    double dz = hits[1].position.Z() - hits[0].position.Z();
    state(3,0) = (fabs(dz) > 1e-6) ? (hits[1].position.X() - hits[0].position.X()) / dz : 0.0;
    state(4,0) = (fabs(dz) > 1e-6) ? (hits[1].position.Y() - hits[0].position.Y()) / dz : 0.0;

    // Initialize covariance - larger errors for initial state
    cov.Zero();
    cov(0,0) = cov(1,1) = cov(2,2) = 1.0; // 1 mm position uncertainty
    cov(3,3) = cov(4,4) = 0.1;             // slope uncertainty

    double chi2 = 0;
    int nUpdates = 0;

     // For momentum estimation
    double total_Bdl = 0.0;
    double ty_entry = 0.0, ty_exit = 0.0;
    bool found_entry = false, found_exit = false;

    for (size_t i = 1; i < hits.size(); ++i) {
        // Get current and next positions
        double z_prev = state(2,0);
        double z_next = hits[i].position.Z();
        double dz = z_next - z_prev;
        
        // Check if segment crosses a magnet
        bool crosses_magnet = IsMagnetExist(z_prev, z_next);
        if (crosses_magnet) {
            double x_avg = 0.5 * (state(0,0) + hits[i].position.X());
            double y_avg = 0.5 * (state(1,0) + hits[i].position.Y());
            double z_avg = 0.5 * (z_prev + z_next);
            double bx = GetFieldX(x_avg, y_avg, z_avg);

            double dz_m = std::abs(dz) * 1e-3; // mm to m
            total_Bdl += std::abs(bx) * dz_m;
            // Estimate momentum for propagation (use truth if available, or a fixed guess)
            //double p_est = hits[0].truthMomentum > 0 ? hits[0].truthMomentum : 10.0; // GeV/c
            //double p_est = 30.0; // might worh for momentum from 10 to 100GeV
            double p_est = seed_p;
            double dtheta = 0.3 * bx * dz_m / p_est; // radians
            state(4,0) += dtheta; // update ty (bending in y-z plane)

            // Record entry/exit slopes
            if (!found_entry) {
                ty_entry = state(4,0);
                found_entry = true;
            }
            ty_exit = state(4,0);
        }

        // Propagate state to next hit
        state(0,0) += dz * state(3,0);
        state(1,0) += dz * state(4,0);
        state(2,0) = z_next;
        
        // Measurement matrix H (3x5)
        TMatrixD H(nMeas, nStates);
        H.Zero();
        H(0,0) = H(1,1) = H(2,2) = 1.0; // maps states to measurement x, y, z
        TVectorD meas(nMeas);
        meas(0) = hits[i].position.X();
        meas(1) = hits[i].position.Y();
        meas(2) = hits[i].position.Z();
        
        // Calculate H * state first
        TVectorD Hstate(nMeas);
        //Hstate = H * state;
        for (int m = 0; m < nMeas; m++) {
            Hstate(m) = 0;
            for (int s = 0; s < nStates; s++) {
                Hstate(m) += H(m,s) * state(s,0);
            }
        }
        // Calculate residual
        TVectorD residual = meas - Hstate; // residual = meas - H * state

        // Measurement covariance (3x3)
        TMatrixDSym R(nMeas);
        R.Zero();

        R(0,0) = R(1,1) = R(2,2) = 0.0025; // scifi resolutions 50um
        // R is diagonal, so we can use TMatrixDSym for efficiency

        // Innovation covariance S = H * cov * H.T() + R
        TMatrixD H_T(TMatrixD::kTransposed, H);
        TMatrixD temp = H * cov;
        TMatrixD S = temp * H_T;
        S += R; // Add measurement noise covariance

        //if (S.Determinant() == 0) {
        if (fabs(S.Determinant()) < 1e-12) {
            std::cerr << "Error: S matrix is singular!" << std::endl;
            continue;
        }
        // Invert S using TMatrixD
        TMatrixD S_inv = S;
        S_inv.Invert();
        // Kalman gain  K = cov * H.T() * S_inv
        TMatrixD K = cov * H_T; // (5x5)*(5x3) = (5x3)
         K *= S_inv; // (5x3)*(3x3) = (5x3)
        // Update state
        TVectorD state_update = K * residual;  // K (5x3) * residual (3x1) → (5x1) vector
        if (state_update.GetNrows() != nStates) {
            std::cerr << "[ERROR] Kalman state update vector size mismatch: "
                      << state_update.GetNrows() << " vs. " << nStates << std::endl;
            continue;
        }
        for (int j = 0; j < nStates; ++j) {
            state(j, 0) += state_update(j);    // Update state using vector indexing
        }
        // Update covariance - Joseph form for numerical stability
        TMatrixD I(TMatrixD::kUnit, cov); // Identity matrix (5x5)
        TMatrixD KH = K * H;              // (5x3)*(3x5) → (5x5)
        TMatrixD IminusKH = I - KH;       // (5x5) - (5x5) → (5x5)
        // First term: IminusKH * cov * IminusKH^T
        TMatrixD tempCov = IminusKH * cov;
        TMatrixD IminusKH_T(TMatrixD::kTransposed, IminusKH);
        TMatrixD cov1 = tempCov * IminusKH_T;  // (5x5)*(5x5) → (5x5)
        // Second term: K * R * K^T
        TMatrixD tempKR = K * R;               // (5x3)*(3x3) → (5x3)
        TMatrixD K_T(TMatrixD::kTransposed, K); // (3x5)
        TMatrixD cov2 = tempKR * K_T;          // (5x3)*(3x5) → (5x5)
        // Final covariance update
        for (int ii = 0; ii < nStates; ++ii) {
            for (int jj = 0; jj < nStates; ++jj) {
                cov(ii, jj) = cov1(ii, jj) + cov2(ii, jj);
            }
        }
        // Chi2 calculation (equivalent to residual^T * S_inv * residual)
        double chi2Cont = 0.0;
        for (int i = 0; i < nMeas; ++i) {
            for (int j = 0; j < nMeas; ++j) {
                chi2Cont += residual(i) * S_inv(i, j) * residual(j);
            }
        }
        chi2 += chi2Cont;
        nUpdates++;
    }
    
    // Final momentum estimation
    double kalman_p = 0;
    double dtheta = ty_exit - ty_entry;
    double dtheta_err = 0.0058; // typical angular resolution in radians for 4 layers of SciFi with 50um resolution
    double kalman_p_err = 0.0;

    if(verbose >=4)
        std::cout << "[KALMAN DEBUG] found_entry: " << found_entry
              << ", ty_entry: " << ty_entry
              << ", ty_exit: " << ty_exit
              << ", dtheta: " << dtheta
              << ", total_Bdl: " << total_Bdl << std::endl;
    
    if (found_entry && fabs(dtheta) > 1e-8 && fabs(total_Bdl) > 1e-8) {
        kalman_p = 0.3 * fabs(total_Bdl) / fabs(dtheta);
        // Error propagation: deltap/p = deltatheta/theta
        kalman_p_err = kalman_p * dtheta_err / fabs(dtheta);
         if(verbose >=4)
            std::cout << "[KALMAN] ty_entry = " << ty_entry
                  << ", ty_exit = " << ty_exit
                  << ", dtheta = " << dtheta
                  << ", total_Bdl = " << total_Bdl
                  << ", p = " << kalman_p << " GeV/c" 
                  << ", dp = " << kalman_p_err << " GeV/c" 
                  << std::endl;
    } else {
        kalman_p = 0;
        kalman_p_err = -1;
    }

    int kalman_charge = 0;
    if (found_entry && fabs(dtheta) > 1e-8) {
        // Get local field with sign at middle of magnetic region
        double z_mid = 0.5 * (hits.front().position.Z() + hits.back().position.Z());
        double y_mid = 0.5 * (hits.front().position.Y() + hits.back().position.Y());
        double x_mid = 0.5 * (hits.front().position.X() + hits.back().position.X());
        //double local_bx = GetFieldX(x_mid, y_mid, z_mid);
        double local_bx = GetFieldBX(y_mid); // Get field with sign


        kalman_charge = (dtheta * local_bx > 0) ? -1 : +1;
    }


    // Store results
    results_.back().kalmanMomentum = kalman_p;
    results_.back().kalmanMomentumErr = kalman_p_err;
    results_.back().kalmanChi2 = (nUpdates > 0) ? chi2/nUpdates : 9999;
    results_.back().kalmanCharge = kalman_charge;


}
///////////////////////////////////////////////////////////
///////   RunGenfitTracklets method implementation /////// 
///////////////////////////////////////////////////////////
void TrackAnalysis::RunGenfitTracklets(const std::vector<Tracklet>& tracklets, double seed_p) {
    if (tracklets.size() < 3) {
        results_.back().genfitMomentum = 0;
        results_.back().genfitChi2 = 9999;
        results_.back().genfitCharge = 0;
        return;
    }

    // Sort by Z
    std::vector<Tracklet> sorted = tracklets;
    std::sort(sorted.begin(), sorted.end(), [](const Tracklet& a, const Tracklet& b) {
        return a.position.Z() < b.position.Z();
    });

    // Debug: print sorted tracklets
    if(verbose >=4)
        std::cout << "[GENFIT] Sorted tracklets (Z order):" << std::endl;
    if(verbose >=4)
    for (const auto& t : sorted) {
        std::cout << "  Station: " << t.station
                  << ", Z: " << t.position.Z()
                  << ", X: " << t.position.X()
                  << ", Y: " << t.position.Y()
                  << ", Dir: (" << t.dir.X() << ", " << t.dir.Y() << ", " << t.dir.Z() << ")"
                  << std::endl;
    }

    // Seed state: position from first, direction from first/last
    //TVector3 dir = (sorted.back().position - sorted.front().position).Unit();
    TVector3 dir(0,0,0);
    for (const auto& t : sorted) dir += t.dir;
    dir = dir.Unit();

    TVectorD stateSeed(6);
    stateSeed[0] = sorted.front().position.X() * 0.1; // mm to cm
    stateSeed[1] = sorted.front().position.Y() * 0.1;
    stateSeed[2] = sorted.front().position.Z() * 0.1;
    stateSeed[3] = seed_p * dir.X();
    stateSeed[4] = seed_p * dir.Y();
    stateSeed[5] = seed_p * dir.Z();

    // Debug: print seed state
    if(verbose >=4)
    {   
        std::cout << "[GENFIT] Seed state (cm, GeV/c): ";
        for (int i = 0; i < 6; ++i) std::cout << stateSeed[i] << " ";
        std::cout << std::endl;
    }

    TMatrixDSym covSeed(6);
    covSeed.Zero();
    covSeed(0,0) = covSeed(1,1) = covSeed(2,2) = 0.01;  // 1 mm^2 in cm^2
    covSeed(3,3) = covSeed(4,4) = covSeed(5,5) = 1.0; // (GeV/c)^2


    auto fitCharge = [&](int pdg, double& momentum, double& chi2, double& p_err, bool& success) {
        success = false;
        momentum = 0;
        chi2 = 1e9;
        p_err = -1;
        try {
            genfit::Track* track = new genfit::Track(new genfit::RKTrackRep(pdg), stateSeed, covSeed);
            for (const auto& t : sorted) {
                TVectorD pos(3);
                pos(0) = t.position.X() * 0.1; // mm to cm
                pos(1) = t.position.Y() * 0.1;
                pos(2) = t.position.Z() * 0.1;
                TMatrixDSym cov(3); 
                //cov.UnitMatrix(); 
                //cov.Zero();
                // CORRECTED: Larger uncertainty for tracklets (they're averaged positions)
                double sciFiSigma = 0.01; // 100um resolutions considered for tracklets instead of 50um
                cov(0,0) = sciFiSigma * sciFiSigma;
                cov(1,1) = sciFiSigma * sciFiSigma;
                cov(2,2) = sciFiSigma * sciFiSigma;
                //cov *= 0.000025; // (50um)^2 = 0.000025 cm^2
                track->insertMeasurement(new genfit::SpacepointMeasurement(pos, cov, t.station, -1, nullptr));
            }
            genfit::DAF daf;
            // Set DAF parameters: more conservative for tracklets
            daf.setMaxIterations(20);
            daf.setAnnealingScheme(0.8, 0.2, 5);
            daf.processTrack(track);

            const genfit::FitStatus* status = track->getFitStatus();
            if (status && status->isFitConverged()) {
                genfit::MeasuredStateOnPlane mop = track->getFittedState();
                TVector3 pvec = mop.getMom();
                //momentum = pvec.Mag();
                // following added to check curvature
                TVector3 pos = mop.getPos();
    
                // Get the track representation to access state parameters
                genfit::AbsTrackRep* rep = track->getTrackRep(0);
                
                // Get the 5D state vector at the fitted plane
                TVectorD state5D = mop.getState();
                
                // For RKTrackRep, the state vector is typically:
                // state[0] = x/y (depending on plane orientation)
                // state[1] = y/x 
                // state[2] = tx = dx/dz (slope in x-z)
                // state[3] = ty = dy/dz (slope in y-z)  
                // state[4] = q/p (charge over momentum)
                
                double qOverP = state5D(4);  // q/p in (e*c)/(GeV)
                double charge = (qOverP > 0) ? 1.0 : -1.0;
                double momentum_from_state = std::abs(1.0 / qOverP);  // GeV/c
                
                // Calculate curvature in the bending plane (y-z for Bx field)
                double tx = state5D(2);  // dx/dz
                double ty = state5D(3);  // dy/dz
                
                // Momentum components
                double px = pvec.X();
                double py = pvec.Y(); 
                double pz = pvec.Z();
                double p_total = pvec.Mag();
                
                // Transverse momentum in bending plane (y-z for Bx field)
                double pt_yz = std::sqrt(py*py + pz*pz);
                
                // Curvature = 1/R, where R is radius of curvature
                // For a charged particle in magnetic field: R = p_t / (0.3 * |B| * |q|)
                // So curvature = 0.3 * |B| * |q| / p_t
                
                // Get magnetic field at current position
                TVector3 bField = magneticField_->get(pos);  // pos is in cm for GENFIT
                double bx = bField.X();  // Tesla
                
                // Calculate curvature (1/R) in m^-1
                double curvature = 0.0;
                if (pt_yz > 0 && std::abs(bx) > 1e-6) {
                    double R_meters = pt_yz / (0.3 * std::abs(bx));  // radius in meters
                    curvature = 1.0 / R_meters;  // curvature in m^-1
                }
                
                // Alternative: calculate curvature directly from q/p and field
                double curvature_alt = 0.0;
                if (std::abs(qOverP) > 1e-8 && std::abs(bx) > 1e-6) {
                    curvature_alt = 0.3 * std::abs(bx) * std::abs(qOverP);  // m^-1
                }
                
                // Debug output
                std::cout << "[DEBUG] GENFIT results:" << std::endl;
                std::cout << "  Position (cm): (" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << ")" << std::endl;
                std::cout << "  Momentum (GeV/c): (" << px << ", " << py << ", " << pz << ")" << std::endl;
                std::cout << "  |p| = " << p_total << " GeV/c" << std::endl;
                std::cout << "  |p| from q/p = " << momentum_from_state << " GeV/c" << std::endl;
                std::cout << "  q/p = " << qOverP << " (e*c)/GeV" << std::endl;
                std::cout << "  charge = " << charge << std::endl;
                std::cout << "  slopes: tx = " << tx << ", ty = " << ty << std::endl;
                std::cout << "  p_t(yz) = " << pt_yz << " GeV/c" << std::endl;
                std::cout << "  B_x = " << bx << " T" << std::endl;
                std::cout << "  Curvature (from p_t/B): " << curvature << " m^-1" << std::endl;
                std::cout << "  Curvature (from q/p): " << curvature_alt << " m^-1" << std::endl;
                if (curvature > 0) {
                    std::cout << "  Radius of curvature: " << 1.0/curvature << " m" << std::endl;
                    std::cout << "  Momentum from R: " << (0.3 * std::abs(bx) * (1.0/curvature)) << " GeV/c" << std::endl;
                }
                
                momentum = p_total;
                

                chi2 = status->getChi2() / status->getNdf();
                // error propagation for momentum
                const TMatrixDSym& covMom = mop.get6DCov();
                //double px = pvec.X(), py = pvec.Y(), pz = pvec.Z();
                double p = momentum;
                if (p > 0 && covMom.GetNrows() >= 6) {
                    double dpx2 = covMom(3,3);
                    double dpy2 = covMom(4,4);
                    double dpz2 = covMom(5,5);
                    double dpxpy = covMom(3,4);
                    double dpxpz = covMom(3,5);
                    double dpypz = covMom(4,5);
                    // Error propagation for |p| = sqrt(px^2 + py^2 + pz^2)
                    p_err = (1.0/p) * std::sqrt(
                        px*px*dpx2 + py*py*dpy2 + pz*pz*dpz2 +
                        2*px*py*dpxpy + 2*px*pz*dpxpz + 2*py*pz*dpypz
                    );
                } else {
                    p_err = -1;
                }
                success = true;
            }
            delete track;
        } catch (...) {}
    };

    double p_minus, chi2_minus, p_plus, chi2_plus, p_err_minus, p_err_plus;
    bool ok_minus, ok_plus;
    fitCharge(13, p_minus, chi2_minus, p_err_minus, ok_minus);    // mu-
    fitCharge(-13, p_plus, chi2_plus, p_err_plus, ok_plus);      // mu+

    // Store results
    if (ok_minus || ok_plus) {
        if (ok_minus && (!ok_plus || chi2_minus < chi2_plus)) {
            results_.back().genfitMomentum = p_minus;
            results_.back().genfitChi2 = chi2_minus;
            results_.back().genfitCharge = -1;
            results_.back().genfitMomentumErr = p_err_minus;
        } else {
            results_.back().genfitMomentum = p_plus;
            results_.back().genfitChi2 = chi2_plus;
            results_.back().genfitCharge = +1;
            results_.back().genfitMomentumErr = p_err_plus;
        }
    } else {
        results_.back().genfitMomentum = 0;
        results_.back().genfitChi2 = 9999;
        results_.back().genfitCharge = 0;
        results_.back().genfitMomentumErr = -1;
    }
    
}