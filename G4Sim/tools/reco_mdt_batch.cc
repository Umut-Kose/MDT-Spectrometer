// =============================================================================
// reco_mdt.cc
// Global Chi2 track reconstruction with Local Left-Right Ambiguity Solving
// Inspired by ATLAS TrkGlobalChi2Fitter (Thijs Cornelissen Thesis methodology)
// =============================================================================

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TDecompLU.h>
#include <TVector3.h>

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TGeoManager.h>
#include <TString.h>
#include "TCanvas.h"
#include "TH2F.h"
#include "TBox.h"
#include "TLatex.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "TLine.h"
#include "TLegend.h"
#include "TString.h"
#include "TObject.h"
#include "TPad.h"
#include "TROOT.h"

// -----------------------------------------------------------------------------
// Core Experimental Configuration
// -----------------------------------------------------------------------------
const double zm_fit[3] = {-531.962, 0.0, +531.962}; // Magnet centers (mm)
const double magnetHalfZ_mm = 200.0;                 // 40 cm total thickness per magnet
const double B_T = 1.5;                              // Field strength (Tesla)
const double sigma_MDT_mm = 0.080;                   // 80 micron intrinsic resolution

// Data structures from your existing environment
struct MDTMeas {
    double wireY_mm;
    double wireZ_mm;
    double r_meas_mm;
    double y_mm; // Set after L/R resolution
    int sta;     // Station ID (0 to 3)
    int pln;     // Plane Layer ID
    int side;    // Left-Right side assignment (-1 or +1)
};

struct GlobalTrackState {
    double y0;       // Position at reference Z = 0 (mm)
    double slope0;   // Slope dY/dZ at reference Z = 0
    double qOverP;   // Charge over momentum (1/GeV)
};

struct RecoOut {
    int eventID;
    int entry;
    int trackID;
    int muonIndex;
    int nMuons;

    int nHits;
    int nStations;
    int fitOK;
    int lrOK;
    int pdg;
    int charge;

    double pTrue;
    double pReco;
    double pErr;
    double pAnalytic;
    double qAnalytic;

    double invPTrue;
    double invPReco;
    double invPResidual;

    double relRes;
    double absRes;
    double relErr;

    double chi2;
    double ndf;
    double pval;
    double manualChi2;
    double manualChi2Ndf;

    double bendAnalytic;
    double bestLRChi2;
    double bestLRNdf;

    double qOverP;
    double qOverPErr;
    double qOverPSignificance;

    double px;
    double py;
    double pz;
};

void resetOut(RecoOut& out)
{
    out.eventID = -1;
    out.entry = -1;
    out.trackID = -1;
    out.muonIndex = -1;
    out.nMuons = 0;

    out.nHits = 0;
    out.nStations = 0;
    out.fitOK = 0;
    out.lrOK = 0;
    out.pdg = 0;
    out.charge = 0;

    out.pTrue = -999.;
    out.pReco = -999.;
    out.pErr = -999.;
    out.pAnalytic = -999.;
    out.qAnalytic = 0.0;

    out.invPTrue = -999.;
    out.invPReco = -999.;
    out.invPResidual = -999.;

    out.relRes = -999.;
    out.absRes = -999.;
    out.relErr = -999.;

    out.chi2 = -999.;
    out.ndf = -999.;
    out.pval = -999.;
    out.manualChi2 = -999.;
    out.manualChi2Ndf = -999.;

    out.bendAnalytic = -999.;
    out.bestLRChi2 = -999.;
    out.bestLRNdf = -999.;

    out.qOverP = -999.;
    out.qOverPErr = -999.;
    out.qOverPSignificance = -999.;

    out.px = -999.;
    out.py = -999.;
    out.pz = -999.;
}

// -----------------------------------------------------------------------------
// Magnetic Field Configuration
// -----------------------------------------------------------------------------
double getBxTesla(double x_mm, double y_mm, double z_mm)
{
    bool inMagnet = false;
    for (double zc : zm_fit) {
        if (std::fabs(z_mm - zc) <= magnetHalfZ_mm) {
            inMagnet = true;
            break;
        }
    }
    if (!inMagnet) return 0.0;
    if (std::fabs(x_mm) > 500.0) return 0.0;

    const double ay = std::fabs(y_mm);
    if (ay < 250.0)  return -1.5;
    if (ay <= 500.0) return +1.5;
    return 0.0;
}

// -----------------------------------------------------------------------------
// Step 1: Local Left-Right Ambiguity Combinatorial Solver
// -----------------------------------------------------------------------------
struct LocalTrackSegment {
    double slope;
    double intercept;
    double chi2;
    std::vector<int> bestSides;
};

LocalTrackSegment solveLocalStationLR(const std::vector<MDTMeas>& stationMeas) 
{
    LocalTrackSegment result;
    result.chi2 = 1e12; // High initialization penalty
    
    int nHits = stationMeas.size();
    if (nHits < 3) {
        result.slope = 0.0; result.intercept = 0.0;
        result.bestSides.assign(nHits, 1);
        return result;
    }

    int nCombos = 1 << nHits; // 2^3 = 8 combinations for 3 layers
    
    for (int combo = 0; combo < nCombos; ++combo) {
        std::vector<int> testSides(nHits);
        for (int i = 0; i < nHits; ++i) {
            testSides[i] = ((combo >> i) & 1) ? 1 : -1;
        }

        double sumZ = 0, sumZ2 = 0, sumY = 0, sumZY = 0;
        std::vector<double> testY(nHits);

        for (int i = 0; i < nHits; ++i) {
            testY[i] = stationMeas[i].wireY_mm + testSides[i] * stationMeas[i].r_meas_mm;
            double z = stationMeas[i].wireZ_mm;
            sumZ  += z;
            sumZ2 += z * z;
            sumY  += testY[i];
            sumZY += z * testY[i];
        }

        double denom = nHits * sumZ2 - sumZ * sumZ;
        if (std::fabs(denom) < 1e-6) continue;

        double slope = (nHits * sumZY - sumZ * sumY) / denom;
        double intercept = (sumY - slope * sumZ) / nHits;

        double currentChi2 = 0.0;
        const double sigma2 = sigma_MDT_mm * sigma_MDT_mm; 

        for (int i = 0; i < nHits; ++i) {
            double expectedY = intercept + slope * stationMeas[i].wireZ_mm;
            double residual = testY[i] - expectedY;
            currentChi2 += (residual * residual) / sigma2;
        }

        if (currentChi2 < result.chi2) {
            result.chi2 = currentChi2;
            result.slope = slope;
            result.intercept = intercept;
            result.bestSides = testSides;
        }
    }
    return result;
}

// -----------------------------------------------------------------------------
// Step 2: Track Seeder (Generates initial GlobalTrackState estimation)
// -----------------------------------------------------------------------------
GlobalTrackState generateInitialSeed(const std::map<int, LocalTrackSegment>& resolvedSegments)
{
    GlobalTrackState seed;
    seed.y0 = 0.0;
    seed.slope0 = 0.0;
    seed.qOverP = 0.01; // Initial default guess of 100 GeV

    // Locate outer bounding tracking layers (Stations 0 and 3)
    auto it0 = resolvedSegments.find(0);
    auto it3 = resolvedSegments.find(3);

    if (it0 != resolvedSegments.end() && it3 != resolvedSegments.end()) {
        // Calculate average linear path across the entire system footprint
        double y_sta0 = it0->second.intercept; // evaluated around regional center
        double y_sta3 = it3->second.intercept;
        
        // Approximate geometric slope anchoring positions
        seed.slope0 = (y_sta3 - y_sta0) / 1500.0; // Assume global baseline length scale
        seed.y0 = (y_sta0 + y_sta3) / 2.0;

        // Extract bend difference between station entry (0) and station exit (3)
        double deltaSlope = it3->second.slope - it0->second.slope;
        double totalBdl = 3.0 * B_T * (2.0 * magnetHalfZ_mm * 1e-3); // 1.8 T*m
        
        if (std::fabs(deltaSlope) > 1e-5) {
            seed.qOverP = deltaSlope / (0.299792458 * totalBdl);
        }
    }
    return seed;
}

// -----------------------------------------------------------------------------
// Step 3: Runge-Kutta Numerical Track Propagator & Jacobian Matrix Generator
// -----------------------------------------------------------------------------
/*
double propagateTrack(const GlobalTrackState& state, double targetZ, 
                       double& dY_dy0, double& dY_dslope0, double& dY_dqOverP) 
{
    double currentZ = 0.0; // Reference origin plane
    double currentY = state.y0;
    double currentSlope = state.slope0;
    
    // Core Variational Matrix Initializations at reference boundary Z = 0
    dY_dy0 = 1.0;     dY_dslope0 = 0.0;       dY_dqOverP = 0.0;
    double dSlope_dy0 = 0.0; double dSlope_dslope0 = 1.0; double dSlope_dqOverP = 0.0;
    
    double stepSize = (targetZ > currentZ) ? 10.0 : -10.0; // mm iteration frames
    
    while (std::fabs(targetZ - currentZ) > 1e-3) {
        if (std::fabs(targetZ - currentZ) < std::fabs(stepSize)) {
            stepSize = targetZ - currentZ;
        }
        
        double Bx = getBxTesla(0.0, currentY, currentZ); 
        double secTheta = std::sqrt(1.0 + currentSlope * currentSlope);
        double curvatureFactor = 0.299792458 * 1e-3 * Bx * secTheta; 
        
        // Advance physics boundaries
        currentY += currentSlope * stepSize;
        currentSlope += curvatureFactor * state.qOverP * stepSize;
        
        // Accumulate directional matrix derivatives (The Tracking Jacobian)
        dY_dy0 += dSlope_dy0 * stepSize;
        dY_dslope0 += dSlope_dslope0 * stepSize;
        dY_dqOverP += dSlope_dqOverP * stepSize;
        
        dSlope_dqOverP += curvatureFactor * stepSize; 
        
        currentZ += stepSize;
    }
    
    return currentY;
}
*/
// 
// -----------------------------------------------------------------------------
// Step 3: Runge-Kutta Numerical Track Propagator with dE/dx slowdown
// -----------------------------------------------------------------------------
double propagateTrack(const GlobalTrackState& state, double targetZ, 
                       double& dY_dy0, double& dY_dslope0, double& dY_dqOverP) 
{
    double currentZ = 0.0; // Reference origin plane
    double currentY = state.y0;
    double currentSlope = state.slope0;
    
    // Track momentum evolves as it loses energy in iron
    // Start with the initial qOverP at the Z=0 reference plane
    double currentQOverP = state.qOverP; 

    // Core Variational Matrix Initializations at reference boundary Z = 0
    dY_dy0 = 1.0;     dY_dslope0 = 0.0;       dY_dqOverP = 0.0;
    double dSlope_dy0 = 0.0; double dSlope_dslope0 = 1.0; double dSlope_dqOverP = 0.0;
    
    double stepSize = (targetZ > currentZ) ? 10.0 : -10.0; // mm iteration frames
    
    // Constant for Bethe-Bloch ionization loss in Iron: ~1.14 MeV/mm
    const double dEdx_GeV_per_mm = 0.00114; 

    while (std::fabs(targetZ - currentZ) > 1e-3) {
        if (std::fabs(targetZ - currentZ) < std::fabs(stepSize)) {
            stepSize = targetZ - currentZ;
        }
        
        double Bx = getBxTesla(0.0, currentY, currentZ); 
        double secTheta = std::sqrt(1.0 + currentSlope * currentSlope);
        
        // Calculate the actual path length element: ds = dZ * sec(theta)
        double ds = std::fabs(stepSize) * secTheta;

        // If inside the iron core (indicated by active magnetic field presence)
        if (std::fabs(Bx) > 0.1) {
            double currentP = 1.0 / std::fabs(currentQOverP);
            
            // Deduct ionization energy loss from the absolute momentum scale
            // Ensure the particle doesn't range out/stop entirely
            if (currentP > 0.5) { 
                double newP = currentP - dEdx_GeV_per_mm * ds;
                int chargeSign = (currentQOverP > 0) ? 1 : -1;
                currentQOverP = chargeSign / newP;
            }
        }

        double curvatureFactor = 0.299792458 * 1e-3 * Bx * secTheta; 
        
        // Advance physics boundaries
        currentY += currentSlope * stepSize;
        currentSlope += curvatureFactor * currentQOverP * stepSize;
        
        // Accumulate directional matrix derivatives (The Tracking Jacobian)
        dY_dy0 += dSlope_dy0 * stepSize;
        dY_dslope0 += dSlope_dslope0 * stepSize;
        dY_dqOverP += dSlope_dqOverP * stepSize;
        
        // The derivative step takes into account the locally degraded currentQOverP scaling
        double dQOverP_dqOverP0 = currentQOverP / (state.qOverP + 1e-15); 
        dSlope_dqOverP += curvatureFactor * dQOverP_dqOverP0 * stepSize; 
        
        currentZ += stepSize;
    }
    
    return currentY;
}

// -----------------------------------------------------------------------------
// Step 4: Global Chi2 Fitter Core Engine
// -----------------------------------------------------------------------------
/*
bool fitGlobalChi2(const std::vector<MDTMeas>& measurements, GlobalTrackState& trackTrack, TMatrixDSym& covMatrix, double& finalChi2) 
{
    const int nTrackParams = 3; // Params: y0, slope0, qOverP
    int nHits = measurements.size();
    
    // Gauss-Newton iterating convergence configuration
    for (int iter = 0; iter < 8; ++iter) {
        TMatrixDSym A(nTrackParams); 
        TVectorD b(nTrackParams);    
        A.Zero(); b.Zero();
        finalChi2 = 0.0;
        
        for (int i = 0; i < nHits; ++i) {
            double dY_dy0, dY_dslope0, dY_dqOverP;
            double predY = propagateTrack(trackTrack, measurements[i].wireZ_mm, dY_dy0, dY_dslope0, dY_dqOverP);
            
            double residual = measurements[i].y_mm - predY;
            
            // Multiple scattering error scale inflation formulation inside iron absorbers
            double sigma_MS = 0.0;   
            if (std::fabs(measurements[i].wireZ_mm) > 200.0) {
                // Radiation length index projection scaling for iron
                sigma_MS = (0.0136 / (0.3 * 1.5)) * std::sqrt(400.0 / 17.6) * std::fabs(measurements[i].wireZ_mm - 200.0) * 1e-3;
            }
            double totalWeight = 1.0 / (sigma_MDT_mm * sigma_MDT_mm + sigma_MS * sigma_MS);
            
            finalChi2 += residual * residual * totalWeight;
            
            double D[3] = {dY_dy0, dY_dslope0, dY_dqOverP};
            for (int r = 0; r < nTrackParams; ++r) {
                for (int c = 0; c < nTrackParams; ++c) {
                    A(r, c) += D[r] * totalWeight * D[c];
                }
                b(r) += D[r] * totalWeight * residual;
            }
        }
        
        TDecompLU lu(A);
        if (!lu.Decompose()) return false; 

        // lu.Solve modifies 'b' in-place so that 'b' now holds the delta_p values
        if (!lu.Solve(b)) return false; 

        trackTrack.y0     += b(0);
        trackTrack.slope0 += b(1);
        trackTrack.qOverP += b(2);

        // Check convergence using the updated b vector
        if (std::fabs(b(2)) < 1e-7) {
            covMatrix = A.Invert(); 
            return true; 
        }
    }
    return false;
}
*/
// -----------------------------------------------------------------------------
// Step 4: Upgraded Global Chi2 Fitter with Full Covariance Scattering Matrix
// -----------------------------------------------------------------------------
bool fitGlobalChi2(const std::vector<MDTMeas>& measurements, GlobalTrackState& trackTrack, TMatrixDSym& covMatrix, double& finalChi2) 
{
    const int nTrackParams = 3; // Params: y0, slope0, qOverP
    int nHits = measurements.size();
    
    // Gauss-Newton iterating convergence configuration
    for (int iter = 0; iter < 8; ++iter) {
        
        // 1. Build the full Hit Covariance Matrix (V) of dimension (nHits x nHits)
        TMatrixDSym V(nHits);
        V.Zero();
        
        // Add independent intrinsic MDT resolution to the diagonal
        for (int i = 0; i < nHits; ++i) {
            V(i, i) = sigma_MDT_mm * sigma_MDT_mm;
        }
        
        // Inject non-diagonal forward-correlated multiple scattering terms
        double currentP = 1.0 / std::max(1e-3, std::fabs(trackTrack.qOverP));
        
        // Compute standard Highland scattering angle per 400mm iron block magnet
        // Radiation length of Iron (X0) = 17.6 mm. Thickness (L) = 400 mm.
        double theta_ms = (0.0136 / currentP) * std::sqrt(400.0 / 17.6) * (1.0 + 0.038 * std::log(400.0 / 17.6));
        double theta_ms2 = theta_ms * theta_ms;
        // Calibration factor K to reconcile the analytical model with the 4-station geometry constraints
        // Since full active matrix yielded 0.07 and zero matrix yielded 1.39, the ideal scale factor is ~0.15
        const double K_calibration = 0.15; 

        for (int i = 0; i < nHits; ++i) {
            for (int j = 0; j < nHits; ++j) {
                double scatteringCov = 0.0;
                
                // Accumulate correlated errors from each of the 3 magnets passed
                for (double z_magnet : zm_fit) {
                    // Elements only correlate if BOTH hits are downstream of the scattering material
                    if (measurements[i].wireZ_mm > z_magnet && measurements[j].wireZ_mm > z_magnet) {
                        double leverArmI = measurements[i].wireZ_mm - z_magnet;
                        double leverArmJ = measurements[j].wireZ_mm - z_magnet;
                        // Apply the calibration scaling factor here
                        scatteringCov += K_calibration *theta_ms2 * leverArmI * leverArmJ;
                    }
                }
                V(i, j) += scatteringCov;
            }
        }
        
        // 2. Invert the Hit Covariance Matrix to obtain the full Weight Matrix (W = V^-1)
        TMatrixDSym W(V);
        W.Invert();
        
        // 3. Accumulate information matrices using the complete tracking Jacobian
        TMatrixDSym A(nTrackParams); 
        TVectorD b(nTrackParams);    
        A.Zero(); b.Zero();
        
        // Temporary storage for track predictions and Jacobians per hit
        std::vector<double> residuals(nHits);
        std::vector<std::vector<double>> Jacobians(nHits, std::vector<double>(nTrackParams));
        
        for (int i = 0; i < nHits; ++i) {
            double dY_dy0, dY_dslope0, dY_dqOverP;
            double predY = propagateTrack(trackTrack, measurements[i].wireZ_mm, dY_dy0, dY_dslope0, dY_dqOverP);
            
            residuals[i] = measurements[i].y_mm - predY;
            Jacobians[i][0] = dY_dy0;
            Jacobians[i][1] = dY_dslope0;
            Jacobians[i][2] = dY_dqOverP;
        }
        
        // Perform non-diagonal matrix multiplications: A = D^T * W * D and b = D^T * W * r
        for (int i = 0; i < nHits; ++i) {
            for (int j = 0; j < nHits; ++j) {
                double weightElement = W(i, j);
                
                for (int r = 0; r < nTrackParams; ++r) {
                    for (int c = 0; c < nTrackParams; ++c) {
                        A(r, c) += Jacobians[i][r] * weightElement * Jacobians[j][c];
                    }
                    b(r) += Jacobians[i][r] * weightElement * residuals[j];
                }
            }
        }
        
        // 4. Compute the true global multi-point Chi2 math block
        finalChi2 = 0.0;
        for (int i = 0; i < nHits; ++i) {
            for (int j = 0; j < nHits; ++j) {
                finalChi2 += residuals[i] * W(i, j) * residuals[j];
            }
        }
        
        // 5. Solve linear step updates using ROOT LU decomposition
        TDecompLU lu(A);
        if (!lu.Decompose()) return false; 
        if (!lu.Solve(b)) return false;
        
        trackTrack.y0     += b(0);
        trackTrack.slope0 += b(1);
        trackTrack.qOverP += b(2);
        
        // Convergence criteria check
        if (std::fabs(b(2)) < 1e-7) {
            covMatrix = A.Invert(); 
            return true; 
        }
    }
    return false;
}

// Global pointers for the output tree (place near the top of your script or inside main)
TFile* fOut = nullptr;
TTree* tOut = nullptr;
RecoOut outData; // Using the RecoOut struct from your original code snippet

void initializeOutputTree(const char* outputFilename) {
    fOut = TFile::Open(outputFilename, "RECREATE");
    tOut = new TTree("Reco", "Reconstructed Muon Tracks");
    
    tOut->Branch("eventID",    &outData.eventID,    "eventID/I");
    tOut->Branch("trackID",    &outData.trackID,    "trackID/I");
    tOut->Branch("nHits",      &outData.nHits,      "nHits/I");
    tOut->Branch("fitOK",      &outData.fitOK,      "fitOK/I");
    tOut->Branch("pReco",      &outData.pReco,      "pReco/D");
    tOut->Branch("pErr",       &outData.pErr,       "pErr/D");
    tOut->Branch("charge",     &outData.charge,     "charge/I");
    tOut->Branch("chi2",       &outData.chi2,       "chi2/D");
    tOut->Branch("ndf",        &outData.ndf,        "ndf/D");
    tOut->Branch("qOverP",     &outData.qOverP,     "qOverP/D");
}

// -----------------------------------------------------------------------------
// Master Execution Coordinator Workflow
// -----------------------------------------------------------------------------
void processEvent(const std::vector<MDTMeas>& rawEventHits, int evID, int trkID) {
    resetOut(outData); // Clear stale fields
    outData.eventID = evID;
    outData.trackID = trkID;
    outData.nHits   = rawEventHits.size();

    // 1. Sort raw hits by station mapping indices
    std::map<int, std::vector<MDTMeas>> stationGroups;
    for (const auto& hit : rawEventHits) {
        stationGroups[hit.sta].push_back(hit);
    }

    // 2. Perform Local L/R Resolving updates
    std::map<int, LocalTrackSegment> resolvedSegments;
    std::vector<MDTMeas> resolvedEventHits;
    for (auto& pair : stationGroups) {
        int stationID = pair.first;
        LocalTrackSegment seg = solveLocalStationLR(pair.second);
        resolvedSegments[stationID] = seg;
        for (size_t i = 0; i < pair.second.size(); ++i) {
            MDTMeas hit = pair.second[i];
            hit.side = seg.bestSides[i];
            hit.y_mm = hit.wireY_mm + hit.side * hit.r_meas_mm; // Resolved 2D spacepoint
            resolvedEventHits.push_back(hit);
        }
    }

    // 3. Extract parameter trajectories seeding anchors
    GlobalTrackState trackState = generateInitialSeed(resolvedSegments);
    TMatrixDSym covariance(3);
    double trackChi2 = 0.0;
    
    bool fitSuccess = fitGlobalChi2(resolvedEventHits, trackState, covariance, trackChi2);
    
    if (fitSuccess) {
        outData.fitOK     = 1;
        outData.pReco     = 1.0 / std::fabs(trackState.qOverP);
        outData.pErr      = std::sqrt(covariance(2,2)) * (outData.pReco * outData.pReco);
        outData.charge    = (trackState.qOverP > 0) ? 1 : -1;
        outData.chi2      = trackChi2;
        outData.ndf       = resolvedEventHits.size() - 3;
        outData.qOverP    = trackState.qOverP;
    } else {
        outData.fitOK     = 0;
    }
    
    // Fill the current row into the output buffer
    if (tOut) tOut->Fill();
}


int main(int argc, char** argv)
{

    const char* inputFile  = (argc > 1) ? argv[1] : "mdt_hits.root";
    const char* outputFile = (argc > 2) ? argv[2] : "reco_results.root";
    const char* gdmlFile   = (argc > 3) ? argv[3] : "detector_geometry.gdml";
    double smear_mm        = (argc > 4) ? std::atof(argv[4]) : 0.080;
    int useMaterial        = (argc > 5) ? std::atoi(argv[5]) : 1;


    std::cout << "\n============================================================\n";
    std::cout << "  Batch MDT GenFit reconstruction, per trackID\n";
    std::cout << "  Input:       " << inputFile << "\n";
    std::cout << "  Output:      " << outputFile << "\n";
    std::cout << "  GDML:        " << gdmlFile << "\n";
    std::cout << "  sigma_r:     " << smear_mm * 1000.0 << " um\n";
    std::cout << "  Material:    " << (useMaterial ? "ON" : "OFF") << "\n";
    std::cout << "============================================================\n\n";


    initializeOutputTree(outputFile);

    // -------------------------------------------------------------------------
    // Open input ROOT file
    // -------------------------------------------------------------------------

    TFile* fin = TFile::Open(inputFile);
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open input file " << inputFile << "\n";
        return 1;
    }

    TTree* tree = (TTree*)fin->Get("Hits");
    if (!tree) {
        std::cerr << "ERROR: cannot find TTree 'Hits'\n";
        fin->Close();
        return 1;
    }

    int eventID = -1;
    int nDumped = 1;
    std::vector<double>* trueDriftRadius = nullptr;
    std::vector<double>* tubeCenterY = nullptr;
    std::vector<double>* tubeCenterZ = nullptr;
    std::vector<double>* tubeCenterX = nullptr;
    std::vector<double>* trueY = nullptr;
    std::vector<double>* trueZ = nullptr;
    std::vector<double>* pz_br = nullptr;
    std::vector<double>* px_br = nullptr;
    std::vector<double>* py_br = nullptr;
    std::vector<double>* hitX_br = nullptr;   // local X along tube (after MDTSD fix)

    std::vector<int>* stationID = nullptr;
    std::vector<int>* planeID = nullptr;
    std::vector<int>* pdg_br = nullptr;
    std::vector<int>* trackID_br = nullptr;

    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("trueDriftRadius", &trueDriftRadius);
    tree->SetBranchAddress("tubeCenterY", &tubeCenterY);
    tree->SetBranchAddress("tubeCenterZ", &tubeCenterZ);
    if (tree->GetBranch("tubeCenterX")) {
        tree->SetBranchAddress("tubeCenterX", &tubeCenterX);
    }
    tree->SetBranchAddress("trueY", &trueY);
    tree->SetBranchAddress("trueZ", &trueZ);
    tree->SetBranchAddress("stationID", &stationID);
    tree->SetBranchAddress("planeID", &planeID);
    tree->SetBranchAddress("pdg", &pdg_br);

    if (!tree->GetBranch("trackID")) {
        std::cerr << "ERROR: branch 'trackID' not found. Needed for multi-muon separation.\n";
        fin->Close();
        return 1;
    }

    tree->SetBranchAddress("trackID", &trackID_br);

    if (tree->GetBranch("pz")) {
        tree->SetBranchAddress("pz", &pz_br);
    }
    if (tree->GetBranch("px")) {
        tree->SetBranchAddress("px", &px_br);
    }
    if (tree->GetBranch("py")) {
        tree->SetBranchAddress("py", &py_br);
    }
    if (tree->GetBranch("hitX")) {
        tree->SetBranchAddress("hitX", &hitX_br);
    }

    // -------------------------------------------------------------------------
    // Event processing loop
    // -------------------------------------------------------------------------
    Long64_t nEntries = tree->GetEntries();
    std::cout << "Processing " << nEntries << " entries from TTree...\n\n";

    // Set up a random engine if you want to smear the true drift radius into a measured one
    std::mt19937 gen(42);
    
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        int nHitsInTree = stationID->size();
        if (nHitsInTree == 0) continue;

        std::cout << "--- Processing EventID: " << eventID << " (Tree Entry: " << entry << ") with " << nHitsInTree << " hits ---\n";

        // Group hit indices by trackID to separate individual particles/muons
        std::map<int, std::vector<int>> trackHitIndices;
        for (int i = 0; i < nHitsInTree; ++i) {
            // Optional: Filter for muons/anti-muons only (PDG 13 or -13)
            if (pdg_br && std::abs((*pdg_br)[i]) != 13) continue;
            
            int tid = (*trackID_br)[i];
            trackHitIndices[tid].push_back(i);
        }

        // Process each track independently
        for (auto const& [trackID, indices] : trackHitIndices) {
            std::vector<MDTMeas> trackHits;

            for (int idx : indices) {
                MDTMeas hit;
                hit.wireY_mm   = (*tubeCenterY)[idx];
                hit.wireZ_mm   = (*tubeCenterZ)[idx];
                
                // Smear the true drift radius to simulate realistic experimental conditions
                double r_true  = (*trueDriftRadius)[idx];
                std::normal_distribution<double> dist(r_true, smear_mm);
                hit.r_meas_mm  = std::max(0.0, dist(gen)); // Ensure radius stays positive
                
                hit.sta        = (*stationID)[idx];
                hit.pln        = (*planeID)[idx];
                hit.side       = 0;   // To be resolved by solveLocalStationLR
                hit.y_mm       = 0.0; // To be populated by solveLocalStationLR

                trackHits.push_back(hit);
            }

            // Skip tracks with insufficient hits to execute the 5-parameter track fit
            if (trackHits.size() < 6) {
                std::cout << "  [TrackID " << trackID << "] Skipped: insufficient hits (" << trackHits.size() << ").\n";
                continue;
            }

            std::cout << "  [TrackID " << trackID << "] Processing reconstruction with " << trackHits.size() << " precision hits:\n";
            
            // Pass the vector and its tracking identities into your master coordination loop
            // This function processes the fit, populates outData, and calls tOut->Fill() internally
            processEvent(trackHits, eventID, trackID);
        }
    }

    // Clean up pointers and close out files cleanly
    fin->Close();
    delete fin;
     // Flush the tracking TTree branches to disk and close out the output file safely
    if (tOut && fOut) {
        fOut->cd();
        tOut->Write();
        fOut->Close();
        std::cout << "\nReconstruction completed cleanly. Results written to: " << outputFile << "\n";
    }
    return 0;
}
