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

TrackAnalysis::TrackAnalysis() : trackRep_(nullptr), magneticField_(nullptr) {
    // Initialize GENFIT track representation (muon)
    trackRep_ = new genfit::RKTrackRep(13); // PDG code for muon
    
    // Set up magnetic field
    magneticField_ = new MagneticField();
    genfit::FieldManager::getInstance()->init(magneticField_);
}

TrackAnalysis::~TrackAnalysis() {
    delete trackRep_;
    delete magneticField_;
}


void TrackAnalysis::ProcessFile(const char* filename) {
    TFile file(filename);
    TTree* tree = (TTree*)file.Get("Hits");
    if (!tree) {
        std::cerr << "Error: Could not find tree 'Hits'" << std::endl;
        return;
    }

    // Branch variables
    int eventID, trackID;
    std::vector<double> *x=0, *y=0, *z=0, *px=0, *py=0, *pz=0;
    std::vector<int> *layerID=0;

    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("eventID", &trackID);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("z", &z);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("layerID", &layerID);

    int nEvents = tree->GetEntries();
    for (int i = 0; i < nEvents; ++i) {
        tree->GetEntry(i);
        ProcessEvent(eventID, trackID, *x, *y, *z, *px, *py, *pz, *layerID);
    }
}

void TrackAnalysis::ProcessEvent(int eventID, 
                               int trackID,
                               const std::vector<double>& x,
                               const std::vector<double>& y,
                               const std::vector<double>& z,
                               const std::vector<double>& px,
                               const std::vector<double>& py,
                               const std::vector<double>& pz,
                               const std::vector<int>& layerID) {
    // Convert to Hit structures
    std::vector<Hit> hits;
    for (size_t i = 0; i < x.size(); ++i) {
        // Determine station ID from z-position 
        int station = static_cast<int>(z[i] / 100.0);
        
        hits.push_back({
            .stationID = station,
            .layerID = layerID[i],
            .position = TVector3(x[i], y[i], z[i]),
            .momentum = TVector3(px[i], py[i], pz[i])
        });
    }

    // Skip events with too few hits
    if (hits.size() < 4) return;

    // Calculate truth momentum from average of first and last hit
    double truthP = 0.5 * (hits.front().momentum.Mag() + hits.back().momentum.Mag()) / 1000.;

    results_.emplace_back();
    results_.back().truthMomentum = truthP;

    RunKalmanFilter(hits);
    RunGenfit(hits);

    // Print event summary
    const auto& r = results_.back();
    printf("Event %d: Truth=%.3f GeV | Kalman=%.3f GeV (χ²=%.1f) | GENFIT=%.3f GeV (χ²=%.1f)\n",
           eventID, r.truthMomentum, r.kalmanMomentum, r.kalmanChi2, 
           r.genfitMomentum, r.genfitChi2);
}

void TrackAnalysis::RunKalmanFilter(const std::vector<Hit>& hits) {
    const int nStates = 5; // x, y, z, tx, ty
    TMatrixD state(nStates, 1);
    TMatrixD cov(nStates, nStates);

    // Initialize with first hit
    state(0,0) = hits[0].position.X();
    state(1,0) = hits[0].position.Y();
    state(2,0) = hits[0].position.Z();
    
    double chi2 = 0;
    int nUpdates = 0;

    for (size_t i = 1; i < hits.size(); ++i) {
        // Prediction
        double dz = hits[i].position.Z() - state(2,0);
        state(0,0) += dz * state(3,0); // x += dz*tx
        state(1,0) += dz * state(4,0); // y += dz*ty
        state(2,0) = hits[i].position.Z();

        // Measurement update
        TVector3 meas = hits[i].position;
        TMatrixD H(3, nStates);
        H(0,0) = H(1,1) = H(2,2) = 1.0;

        TMatrixD residual(3,1);
        residual(0,0) = meas.X() - state(0,0);
        residual(1,0) = meas.Y() - state(1,0);
        residual(2,0) = meas.Z() - state(2,0);

        TMatrixD R(3,3);
        R(0,0) = R(1,1) = R(2,2) = 0.01; // 100 μm resolution

        TMatrixD S = H * cov * H.T() + R;
        TMatrixD K = cov * H.T() * S.Invert();

        state += K * residual;
        cov -= K * H * cov;

        // Chi2 calculation
        TMatrixD chi2Cont = residual.T() * S.Invert() * residual;
        chi2 += chi2Cont(0,0);
        nUpdates++;
    }

    results_.back().kalmanMomentum = CalculateSagittaMomentum(hits);
    results_.back().kalmanChi2 = (nUpdates > 0) ? chi2/nUpdates : 9999;
}

void TrackAnalysis::RunGenfit(const std::vector<Hit>& hits) {
    genfit::Track* track = new genfit::Track(trackRep_, 
                                           hits[0].position, 
                                           TVector3(0,0,10)); // Seed with 10 GeV

    // Create measurements
    for (const auto& hit : hits) {
        TVectorD hitCoords(3); // considering we have x,y,z positions
        hitCoords(0) = hit.position.X();
        hitCoords(1) = hit.position.Y();
        hitCoords(2) = hit.position.Z();

        TMatrixDSym hitCov(3);
        hitCov(0,0) = hitCov(1,1) = hitCov(2,2) = 0.01; // 100 μm resolution
        hitCov(0,1) = hitCov(0,2) = hitCov(1,0) = hitCov(1,2) = hitCov(2,0) = hitCov(2,1) = 0.0;

        std::cout << "[DEBUG] hitCoords size: " << hitCoords.GetNrows() 
          << ", hitCov size: " << hitCov.GetNrows() << "x" << hitCov.GetNcols() << std::endl;
        
          for (int i = 0; i < hitCoords.GetNrows(); ++i)
            std::cout << "  hitCoords[" << i << "] = " << hitCoords[i] << std::endl;


        genfit::SpacepointMeasurement* meas = new genfit::SpacepointMeasurement(
            hitCoords, hitCov, hit.stationID, hit.layerID,nullptr);
        track->insertMeasurement(meas);
    }

    // Fit with DAF
    genfit::DAF daf;
    try {
        daf.processTrack(track);
        
        //const genfit::MeasuredStateOnPlane& mop = track->getFittedState(0, trackRep_);
        //track->getFittedState(mop, 0, trackRep_);
       // For the first error:
        genfit::MeasuredStateOnPlane mop;
        // Correct way to get the fitted state is to assign the return value:
        mop = track->getFittedState(); // Use default arguments (id=0)

        // For the fitStatus errors:
    const genfit::FitStatus* fitStatus = track->getFitStatus();
    if (fitStatus) {
        results_.back().genfitChi2 = fitStatus->getChi2() / fitStatus->getNdf();
    }
    else {
            results_.back().genfitChi2 = 9999;
        }
    } catch (...) {
        results_.back().genfitMomentum = 0;
        results_.back().genfitChi2 = 9999;
    }

    delete track;
}

double TrackAnalysis::CalculateSagittaMomentum(const std::vector<Hit>& hits) {
    if (hits.size() < 3) return 0;

    const TVector3& p0 = hits[0].position;
    const TVector3& p1 = hits[hits.size()/2].position;
    const TVector3& p2 = hits.back().position;

    double L = p2.Z() - p0.Z();
    double sagitta = p1.X() - 0.5*(p0.X() + p2.X());
    double B = 1.5; // Tesla

    return 0.3 * B * L*L / (8 * std::abs(sagitta)); // GeV/c
}

void TrackAnalysis::WriteResults(const char* filename) {
    TFile outfile(filename, "RECREATE");
    TTree tree("results", "Reconstruction Results");

    double truth, kalman, genfit, chi2_kf, chi2_gf;
    tree.Branch("truth", &truth, "truth/D");
    tree.Branch("kalman", &kalman, "kalman/D");
    tree.Branch("genfit", &genfit, "genfit/D");
    tree.Branch("chi2_kf", &chi2_kf, "chi2_kf/D");
    tree.Branch("chi2_gf", &chi2_gf, "chi2_gf/D");

    for (const auto& r : results_) {
        truth = r.truthMomentum;
        kalman = r.kalmanMomentum;
        genfit = r.genfitMomentum;
        chi2_kf = r.kalmanChi2;
        chi2_gf = r.genfitChi2;
        tree.Fill();
    }

    tree.Write();
    outfile.Close();
}
