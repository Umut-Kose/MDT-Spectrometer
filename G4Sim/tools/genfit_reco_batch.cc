// =============================================================================
// genfit_reco_batch.cc
// Batch ATLAS-MDT-like GenFit reconstruction, one fit per muon trackID
// =============================================================================
//
// Build:
//   make -f Makefile.genfit_batch
//
// Run:
//   ./genfit_reco_batch mdt_hits.root reco_results.root detector_geometry.gdml 0.080 1
//
// Arguments:
//   argv[1] input ROOT file      default: mdt_hits.root
//   argv[2] output ROOT file     default: reco_results.root
//   argv[3] GDML geometry file   default: detector_geometry.gdml
//   argv[4] smear sigma [mm]     default: 0.080
//   argv[5] use material         default: 1
//
// Output tree:
//   reco_results.root / Reco
//
// One output row = one reconstructed muon trackID
// =============================================================================

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TGeoManager.h>
#include <TString.h>

#include <RKTrackRep.h>
#include <Track.h>
#include <TrackPoint.h>
#include <PlanarMeasurement.h>
#include <KalmanFitterRefTrack.h>
#include <FitStatus.h>
#include <FieldManager.h>
#include <MaterialEffects.h>
#include <TGeoMaterialInterface.h>
#include <AbsBField.h>
#include <Exception.h>
#include <MeasuredStateOnPlane.h>
#include <DetPlane.h>

const double zm_fit[3] = {-531.962, 0.0, +531.962}; // mm
const double magnetHalfZ_mm = 200.0;
const double B_T = 1.5;
const int nMagnets = 3;

const double totalMagLength_m = nMagnets * 2.0 * magnetHalfZ_mm * 1e-3; // 1.2 m

double trackY(double z, const double* p)
{
    double y = p[0] + p[1] * z;
    for (int m = 0; m < 3; m++) {
        if (z > zm_fit[m]) y += p[m + 2] * (z - zm_fit[m]);
    }
    return y;
}

void basisRow(double z, double* A)
{
    A[0] = 1.0;
    A[1] = z;
    for (int m = 0; m < 3; m++) {
        A[m + 2] = (z > zm_fit[m]) ? (z - zm_fit[m]) : 0.0;
    }
}

bool solve5(double ATA[5][5], double ATy[5], double p[5])
{
    double M[5][6];

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) M[i][j] = ATA[i][j];
        M[i][5] = ATy[i];
    }

    for (int col = 0; col < 5; col++) {
        int mx = col;
        for (int r = col + 1; r < 5; r++) {
            if (std::fabs(M[r][col]) > std::fabs(M[mx][col])) mx = r;
        }

        if (std::fabs(M[mx][col]) < 1e-15) return false;

        for (int j = 0; j <= 5; j++) std::swap(M[col][j], M[mx][j]);

        for (int r = 0; r < 5; r++) {
            if (r == col) continue;
            double f = M[r][col] / M[col][col];
            for (int j = col; j <= 5; j++) M[r][j] -= f * M[col][j];
        }
    }

    for (int i = 0; i < 5; i++) p[i] = M[i][5] / M[i][i];
    return true;
}

class MDTMagneticField : public genfit::AbsBField {
public:
    TVector3 get(const TVector3& pos_cm) const override
    {
        const double x_mm = pos_cm.X() * 10.0;
        const double z_mm = pos_cm.Z() * 10.0;
        const double y_mm = pos_cm.Y() * 10.0;

        bool inMagnet = false;
        for (int m = 0; m < 3; m++) {
            if (std::fabs(z_mm - zm_fit[m]) < magnetHalfZ_mm) {
                inMagnet = true;
                break;
            }
        }

        if (!inMagnet) return TVector3(0.0, 0.0, 0.0);
        if (std::fabs(x_mm) > 500.0) return TVector3(0.0, 0.0, 0.0);

        const double ay = std::fabs(y_mm);

        if      (ay < 250.0) return TVector3(-15.0, 0.0, 0.0);
        else if (ay < 500.0) return TVector3(+15.0, 0.0, 0.0);

        return TVector3(0.0, 0.0, 0.0);
    }
};

double getBxTesla(double x_mm, double y_mm, double z_mm)
{
    // Same geometry as MDTMagneticField
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

struct Hit {
    double wx;   // tube center global X
    double wy;
    double wz;
    double r;
    double rtrue;
    double ty;
    double tz;
    double hitX;   // local X along tube (from MDTSD: closestPos.x - tubeCenter.x)
    int sta;
    int pln;
};

struct MDTMeas {
    double x_mm;      // global X of hit = tubeCenterX + hitX (local)
    double y_mm;
    double z_mm;
    double wireY_mm;
    double wireZ_mm;
    double r_meas_mm;
    double r_true_mm;
    double trueY_mm;
    double trueZ_mm;
    int sta;
    int pln;
    int side;
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

void dumpEventDebug(const RecoOut& out, const std::vector<Hit>& hits, const std::vector<MDTMeas>& meas, const double* bestP, double signedBdl_Tm
)
{
    std::cout << "\n============================================================\n";
    std::cout << "DEBUG EVENT DUMP\n";
    std::cout << "eventID=" << out.eventID
              << " entry=" << out.entry
              << " trackID=" << out.trackID
              << " pdg=" << out.pdg
              << " chargeTruth=" << out.charge
              << " nHits=" << out.nHits
              << " nStations=" << out.nStations
              << "\n";

    std::cout << "pTrue      = " << out.pTrue << " GeV\n";
    std::cout << "pAnalytic  = " << out.pAnalytic << " GeV\n";
    std::cout << "pReco      = " << out.pReco << " GeV\n";
    std::cout << "pErr       = " << out.pErr << " GeV\n";
    std::cout << "relRes     = " << out.relRes << "\n";
    std::cout << "invPRes    = " << out.invPResidual << " 1/GeV\n";
    std::cout << "bend       = " << out.bendAnalytic << " rad\n";
    std::cout << "signedBdl  = " << signedBdl_Tm << " Tm\n";
    std::cout << "qAnalytic  = " << out.qAnalytic << "\n";
    std::cout << "qOverP     = " << out.qOverP << "\n";
    std::cout << "qOverPErr  = " << out.qOverPErr << "\n";
    std::cout << "qSignif    = " << out.qOverPSignificance << "\n";

    std::cout << "\nAnalytic fit parameters:\n";
    std::cout << "p0 y0      = " << bestP[0] << " mm\n";
    std::cout << "p1 slope0  = " << bestP[1] << "\n";
    std::cout << "bend1      = " << bestP[2] << "\n";
    std::cout << "bend2      = " << bestP[3] << "\n";
    std::cout << "bend3      = " << bestP[4] << "\n";

    std::cout << "\nHits:\n";
    std::cout << " i sta pln"
              << " wireY wireZ"
              << " rTrue rReco"
              << " trueY trueZ"
              << " recoY side"
              << " yFit residual_um Bx_T\n";

    for (size_t i = 0; i < meas.size(); ++i) {
        double yFit = trackY(meas[i].z_mm, bestP);
        double residual_um = (yFit - meas[i].y_mm) * 1000.0;
        double Bx_T = getBxTesla(0.0, yFit, meas[i].z_mm);

        std::cout << i << " "
                  << meas[i].sta << " "
                  << meas[i].pln << " "
                  << meas[i].wireY_mm << " "
                  << meas[i].wireZ_mm << " "
                  << meas[i].r_true_mm << " "
                  << meas[i].r_meas_mm << " "
                  << meas[i].trueY_mm << " "
                  << meas[i].trueZ_mm << " "
                  << meas[i].y_mm << " "
                  << meas[i].side << " "
                  << yFit << " "
                  << residual_um << " "
                  << Bx_T << "\n";
    }

    std::cout << "============================================================\n";
}


int main(int argc, char** argv)
{
    if (argc < 6) {
        std::cerr
            << "Usage:\n"
            << "  " << argv[0]
            << " mdt_hits.root"
            << " reco_results.root"
            << " detector_geometry.gdml"
            << " driftSigma_mm"
            << " useMaterial(0|1)"
            << " [nDebugDump]\n\n"
            << "Example:\n"
            << "  ./genfit_reco_batch "
            << "mdt_hits.root "
            << "reco_results.root "
            << "detector_geometry.gdml "
            << "0.080 "
            << "1 "
            << "10\n\n"
            << "Arguments:\n"
            << "  driftSigma_mm : MDT resolution in mm\n"
            << "  useMaterial   : 0=no material effects, 1=Highland MSC\n"
            << "  nDebugDump    : optional number of events to dump (default=5)\n";
        return 1;
    }
    const char* inputFile  = (argc > 1) ? argv[1] : "mdt_hits.root";
    const char* outputFile = (argc > 2) ? argv[2] : "reco_results.root";
    const char* gdmlFile   = (argc > 3) ? argv[3] : "detector_geometry.gdml";
    double smear_mm        = (argc > 4) ? std::atof(argv[4]) : 0.080;
    int useMaterial        = (argc > 5) ? std::atoi(argv[5]) : 1;
    int nDebugDump = (argc > 6) ? std::atoi(argv[6]) : 5;

    std::cout << "\n============================================================\n";
    std::cout << "  Batch MDT GenFit reconstruction, per trackID\n";
    std::cout << "  Input:       " << inputFile << "\n";
    std::cout << "  Output:      " << outputFile << "\n";
    std::cout << "  GDML:        " << gdmlFile << "\n";
    std::cout << "  sigma_r:     " << smear_mm * 1000.0 << " um\n";
    std::cout << "  Material:    " << (useMaterial ? "ON" : "OFF") << "\n";
    std::cout << "  Debug dumps: " << nDebugDump << "\n";
    std::cout << "============================================================\n\n";

    TGeoManager::Import(gdmlFile);
    if (!gGeoManager) {
        std::cerr << "ERROR: could not import GDML file: " << gdmlFile << "\n";
        return 1;
    }

    gGeoManager->CloseGeometry();
    gGeoManager->SetVerboseLevel(0);

    genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
    genfit::MaterialEffects::getInstance()->setDebugLvl(0);

    if (useMaterial) {
        genfit::MaterialEffects::getInstance()->setNoEffects(false);
        genfit::MaterialEffects::getInstance()->setMscModel("Highland");
    } else {
        genfit::MaterialEffects::getInstance()->setNoEffects(true);
    }

    genfit::FieldManager::getInstance()->init(new MDTMagneticField());

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
    int nDumped = 0;

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

    TFile* fout = new TFile(outputFile, "RECREATE");
    TTree* tout = new TTree("Reco", "Per-muon MDT GenFit reconstruction results");

    RecoOut out;
    resetOut(out);

    tout->Branch("eventID", &out.eventID, "eventID/I");
    tout->Branch("entry", &out.entry, "entry/I");
    tout->Branch("trackID", &out.trackID, "trackID/I");
    tout->Branch("muonIndex", &out.muonIndex, "muonIndex/I");
    tout->Branch("nMuons", &out.nMuons, "nMuons/I");

    tout->Branch("nHits", &out.nHits, "nHits/I");
    tout->Branch("nStations", &out.nStations, "nStations/I");
    tout->Branch("fitOK", &out.fitOK, "fitOK/I");
    tout->Branch("lrOK", &out.lrOK, "lrOK/I");
    tout->Branch("pdg", &out.pdg, "pdg/I");
    tout->Branch("charge", &out.charge, "charge/I");

    tout->Branch("pTrue", &out.pTrue, "pTrue/D");
    tout->Branch("pReco", &out.pReco, "pReco/D");
    tout->Branch("pErr", &out.pErr, "pErr/D");
    tout->Branch("pAnalytic", &out.pAnalytic, "pAnalytic/D");
    tout->Branch("qAnalytic", &out.qAnalytic, "qAnalytic/D");

    tout->Branch("invPTrue",&out.invPTrue,"invPTrue/D");
    tout->Branch("invPReco",&out.invPReco,"invPReco/D");
    tout->Branch("invPResidual",&out.invPResidual,"invPResidual/D");

    tout->Branch("relRes", &out.relRes, "relRes/D");
    tout->Branch("absRes", &out.absRes, "absRes/D");
    tout->Branch("relErr", &out.relErr, "relErr/D");

    tout->Branch("chi2", &out.chi2, "chi2/D");
    tout->Branch("ndf", &out.ndf, "ndf/D");
    tout->Branch("pval", &out.pval, "pval/D");
    tout->Branch("manualChi2", &out.manualChi2, "manualChi2/D");
    tout->Branch("manualChi2Ndf", &out.manualChi2Ndf, "manualChi2Ndf/D");

    tout->Branch("bendAnalytic", &out.bendAnalytic, "bendAnalytic/D");
    tout->Branch("bestLRChi2", &out.bestLRChi2, "bestLRChi2/D");
    tout->Branch("bestLRNdf", &out.bestLRNdf, "bestLRNdf/D");

    tout->Branch("qOverP", &out.qOverP, "qOverP/D");
    tout->Branch("qOverPErr", &out.qOverPErr, "qOverPErr/D");
    tout->Branch("qOverPSignificance", &out.qOverPSignificance, "qOverPSignificance/D");

    tout->Branch("px", &out.px, "px/D");
    tout->Branch("py", &out.py, "py/D");
    tout->Branch("pz", &out.pz, "pz/D");

    const double tubeInnerRadius = 14.6;
    const double detRes_cm = (smear_mm > 0.0) ? smear_mm / 10.0 : 0.001;
    const double sigmaX_cm = 1.0;

    // smearRadius is defined per-muon inside the loop using a deterministic seed

    Long64_t nEntries = tree->GetEntries();

    Long64_t nMuonRows = 0;
    Long64_t nFitOK = 0;

    for (Long64_t iev = 0; iev < nEntries; iev++) {
        tree->GetEntry(iev);


        if (iev < 10) {
            std::cout << "\nEntry " << iev
                    << " eventID=" << eventID
                    << " nHits=" << (pdg_br ? pdg_br->size() : 0)
                    << std::endl;

            if (pdg_br) {
                std::map<int,int> pdgCount;
                for (int pdg : *pdg_br) pdgCount[pdg]++;

                for (const auto& kv : pdgCount) {
                    std::cout << "  pdg=" << kv.first
                            << " count=" << kv.second << std::endl;
                }
            }
        }


        if (!trueDriftRadius || !tubeCenterY || !tubeCenterZ ||
            !trueY || !trueZ || !stationID || !planeID || !pdg_br || !trackID_br) {
            continue;
        }

        std::map<int, std::vector<int>> muonGroups;

        for (int i = 0; i < (int)pdg_br->size(); i++) {
            if (std::abs((*pdg_br)[i]) != 13) continue;
            muonGroups[(*trackID_br)[i]].push_back(i);
        }

        int nMuons = (int)muonGroups.size();
        int muonIndex = 0;

        for (auto& kv : muonGroups) {
            int tid = kv.first;
            const std::vector<int>& idxs = kv.second;

            resetOut(out);

            out.entry = (int)iev;
            out.eventID = eventID;
            out.trackID = tid;
            out.muonIndex = muonIndex++;
            out.nMuons = nMuons;

            int pdg_muon = 13;

            for (int idx : idxs) {
                if (std::abs((*pdg_br)[idx]) == 13) {
                    pdg_muon = (*pdg_br)[idx];
                    break;
                }
            }

            out.pdg = pdg_muon;

            // Per-muon RNG: deterministic seed from eventID+trackID
            // ensures smearing is identical whether running single-event or batch
            std::mt19937 muon_rng(42 + (uint32_t)(eventID * 10000 + tid));
            auto smearRadius = [&](double r_true) -> double {
                if (smear_mm <= 0.0) return r_true;
                std::normal_distribution<double> gaussR(0.0, smear_mm);
                double r = r_true + gaussR(muon_rng);
                if (r < 0.0) r = 0.0;
                if (r > tubeInnerRadius) r = tubeInnerRadius;
                return r;
            };

            // pTrue: total momentum |p| = sqrt(px²+py²+pz²) of this trackID's first hit.
            if (pz_br && !idxs.empty()) {
                int i0 = idxs[0];
                double ipx = (px_br && i0 < (int)px_br->size()) ? (*px_br)[i0] : 0.0;
                double ipy = (py_br && i0 < (int)py_br->size()) ? (*py_br)[i0] : 0.0;
                double ipz = (*pz_br)[i0];
                out.pTrue = std::sqrt(ipx*ipx + ipy*ipy + ipz*ipz) / 1000.0;
            }

            std::map<std::pair<int, int>, int> seen;
            std::set<int> stations;
            std::vector<Hit> hits;

            for (int idx : idxs) {
                auto key = std::make_pair((*stationID)[idx], (*planeID)[idx]);
                if (seen.count(key)) continue;
                seen[key] = 1;

                Hit h;
                h.wx = (tubeCenterX && idx < (int)tubeCenterX->size()) ? (*tubeCenterX)[idx] : 0.0;
                h.wy = (*tubeCenterY)[idx];
                h.wz = (*tubeCenterZ)[idx];
                h.rtrue = (*trueDriftRadius)[idx];
                h.r = smearRadius(h.rtrue);
                h.ty = (*trueY)[idx];
                h.tz = (*trueZ)[idx];
                h.hitX = (hitX_br && idx < (int)hitX_br->size()) ? (*hitX_br)[idx] : 0.0;
                h.sta = (*stationID)[idx];
                h.pln = (*planeID)[idx];

                hits.push_back(h);
                stations.insert(h.sta);
            }

            std::sort(hits.begin(), hits.end(),
                      [](const Hit& a, const Hit& b) {
                          if (a.sta != b.sta) return a.sta < b.sta;
                          return a.pln < b.pln;
                      });

            const int N = (int)hits.size();
            out.nHits = N;
            out.nStations = (int)stations.size();

            nMuonRows++;

            if (N < 5 || N > 20 || out.pTrue <= 0.0) {
                tout->Fill();
                continue;
            }

            // Warn: brute-force LR enumeration is O(2^N * N).
            // For N=20 this is ~20M iterations per track — can be slow.
            if (N > 14) {
                std::cerr << "WARNING: entry " << iev << " trackID " << tid
                          << " has N=" << N << " hits; 2^N=" << (1 << N)
                          << " LR combos — this will be slow.\n";
            }

            const double sigma = (smear_mm > 0.0) ? smear_mm : 0.001;

            double bestChi2 = 1e99;
            double bestP[5] = {0, 0, 0, 0, 0};
            int bestCombo = 0;
            int nComb = 1 << N;

            for (int combo = 0; combo < nComb; combo++) {
                double ATA[5][5] = {};
                double ATy[5] = {};

                for (int h = 0; h < N; h++) {
                    int sign = ((combo >> h) & 1) ? +1 : -1;
                    double y_assigned = hits[h].wy + sign * hits[h].r;

                    double A[5];
                    basisRow(hits[h].wz, A);

                    for (int i = 0; i < 5; i++) {
                        for (int j = 0; j < 5; j++) {
                            ATA[i][j] += A[i] * A[j];
                        }
                        ATy[i] += A[i] * y_assigned;
                    }
                }

                double p[5];
                if (!solve5(ATA, ATy, p)) continue;

                double chi2 = 0.0;
                for (int h = 0; h < N; h++) {
                    double yfit = trackY(hits[h].wz, p);
                    double dist = std::fabs(yfit - hits[h].wy);
                    double res = (dist - hits[h].r) / sigma;
                    chi2 += res * res;
                }

                if (chi2 < bestChi2) {
                    bestChi2 = chi2;
                    bestCombo = combo;
                    for (int i = 0; i < 5; i++) bestP[i] = p[i];
                }
            }

            if (bestChi2 >= 1e98) {
                tout->Fill();
                continue;
            }

            out.lrOK = 1;
            out.bestLRChi2 = bestChi2;
            out.bestLRNdf = N - 5;

            double totalBendAnalytic = bestP[2] + bestP[3] + bestP[4];
            out.bendAnalytic = totalBendAnalytic;
            //if (std::fabs(totalBendAnalytic) > 1e-8) {
            //    out.pAnalytic = 0.3 * B_T * totalMagLength_m / std::fabs(totalBendAnalytic);
            //}
            double signedBdl_Tm = 0.0;
            for (int m = 0; m < 3; ++m) {
                double zc_mm = zm_fit[m];
                // fitted y position at magnet center
                double y_mm = trackY(zc_mm, bestP);
                // for this simplified MDT fit, x is not reconstructed
                double x_mm = 0.0;
                double Bx_T = getBxTesla(x_mm, y_mm, zc_mm);
                double L_m = 2.0 * magnetHalfZ_mm * 1e-3; // 400 mm = 0.4 m
                signedBdl_Tm += Bx_T * L_m;
            }

            if (std::fabs(totalBendAnalytic) > 1e-8 &&
                std::fabs(signedBdl_Tm) > 1e-8) {
                out.pAnalytic =
                    0.3 * std::fabs(signedBdl_Tm) / std::fabs(totalBendAnalytic);
                // Sign convention may need one truth-based flip
                out.qAnalytic =
                    (totalBendAnalytic * signedBdl_Tm > 0.0) ? -1.0 : +1.0;
            }
            else {
                out.pAnalytic = -999.0;
                out.qAnalytic = 0.0;
            }

            std::vector<MDTMeas> meas;

            for (int h = 0; h < N; h++) {
                int sign = ((bestCombo >> h) & 1) ? +1 : -1;

                MDTMeas m;
                // Global X of hit = tube center X + local hitX along tube
                m.x_mm = hits[h].wx + hits[h].hitX;
                m.y_mm = hits[h].wy + sign * hits[h].r;
                m.z_mm = hits[h].wz;
                m.wireY_mm = hits[h].wy;
                m.wireZ_mm = hits[h].wz;
                m.r_meas_mm = hits[h].r;
                m.r_true_mm = hits[h].rtrue;
                m.trueY_mm = hits[h].ty;
                m.trueZ_mm = hits[h].tz;
                m.sta = hits[h].sta;
                m.pln = hits[h].pln;
                m.side = sign;

                meas.push_back(m);
            }

            genfit::AbsTrackRep* rep = nullptr;
            genfit::Track* fitTrack = nullptr;
            genfit::AbsKalmanFitter* fitter = nullptr;

            try {
                rep = new genfit::RKTrackRep(pdg_muon);

                TVector3 pos_seed(0.0,
                                  meas.front().y_mm / 10.0,
                                  meas.front().z_mm / 10.0);

                TVector3 direction(0.0, 0.0, 1.0);
                if (meas.size() > 1) {
                    TVector3 a(0.0, meas[0].y_mm, meas[0].z_mm);
                    TVector3 b(0.0, meas[1].y_mm, meas[1].z_mm);
                    TVector3 d = b - a;
                    if (d.Mag() > 0.0) direction = d.Unit();
                }

                //double p_seed = (out.pAnalytic > 0.0 && std::isfinite(out.pAnalytic))
                //              ? out.pAnalytic
                //              : 10.0;

                double p_seed = 10.0; //GeV/c

                TVector3 mom_seed = direction * p_seed;

                fitTrack = new genfit::Track(rep, pos_seed, mom_seed);

                int hitId = 0;
                int planeId = 0;

                for (const auto& m : meas) {
                    TVectorD hitCoords(2);
                    hitCoords[0] = m.x_mm / 10.0;   // global X in cm
                    hitCoords[1] = m.y_mm / 10.0;

                    TMatrixDSym hitCov(2);
                    hitCov.Zero();
                    hitCov(0,0) = sigmaX_cm * sigmaX_cm;
                    hitCov(1,1) = detRes_cm * detRes_cm;

                    genfit::PlanarMeasurement* pm =
                        new genfit::PlanarMeasurement(hitCoords, hitCov, 0, ++hitId, nullptr);

                    pm->setPlane(
                        genfit::SharedPlanePtr(new genfit::DetPlane(
                            TVector3(0.0, 0.0, m.z_mm / 10.0),
                            TVector3(1.0, 0.0, 0.0),
                            TVector3(0.0, 1.0, 0.0))),
                        ++planeId);

                    fitTrack->insertPoint(new genfit::TrackPoint(pm, fitTrack));
                }

                fitTrack->checkConsistency();

                fitter = new genfit::KalmanFitterRefTrack();
                fitter->setMaxIterations(10);
                fitter->setRelChi2Change(0.001);

                fitter->processTrackWithRep(fitTrack, rep);

                bool fitConverged = fitTrack->getFitStatus(rep)->isFitConverged();
                if (!fitConverged) {
                    delete fitter;
                    delete fitTrack;
                    tout->Fill();
                    continue;
                }

                genfit::FitStatus* status = fitTrack->getFitStatus(rep);

                out.chi2 = status->getChi2();
                out.ndf = status->getNdf();
                out.pval = status->getPVal();

                genfit::MeasuredStateOnPlane state = fitTrack->getFittedState(0, rep);
                TVector3 p_fit = state.getMom();

                out.pReco = p_fit.Mag();
                out.pErr = std::sqrt(state.getMomVar());
                out.charge = (int)std::round(state.getCharge());
                out.px = p_fit.X();
                out.py = p_fit.Y();
                out.pz = p_fit.Z();

                out.invPTrue = 1.0 / out.pTrue;
                out.invPReco = 1.0 / out.pReco;
                out.invPResidual = out.invPReco - out.invPTrue;

                double manualChi2 = 0.0;

                for (size_t i = 0; i < meas.size(); i++) {
                    genfit::MeasuredStateOnPlane st = fitTrack->getFittedState((int)i, rep);
                    TVector3 pos = st.getPos();

                    double yFit_mm = pos.Y() * 10.0;
                    double res_sig = (yFit_mm - meas[i].y_mm) / smear_mm;
                    manualChi2 += res_sig * res_sig;
                }

                out.manualChi2 = manualChi2;
                // Use GenFit's NDF (from the Kalman fit) not the analytic-fit NDF (N-5)
                out.manualChi2Ndf = (out.ndf > 0.0) ? manualChi2 / out.ndf : -999.;

                if (out.pTrue > 0.0 && out.pReco > 0.0) {
                    out.absRes = out.pReco - out.pTrue;
                    out.relRes = (out.pReco - out.pTrue) / out.pTrue;
                    out.relErr = out.pErr / out.pReco;
                }

                double q = state.getCharge();

                if (out.pReco > 0.0) {
                    out.qOverP = q / out.pReco;
                }

                out.qOverPErr = std::sqrt(state.getCov()(0,0));

                if (out.qOverPErr > 0.0 && std::isfinite(out.qOverPErr)) {
                    out.qOverPSignificance = std::fabs(out.qOverP) / out.qOverPErr;
                }

                out.fitOK = 1;
                nFitOK++;

                delete fitter;
                delete fitTrack;
            }
            catch (genfit::Exception& e) {
                std::cerr << "GenFit exception at entry " << iev
                          << " eventID " << eventID
                          << " trackID " << tid
                          << ": " << e.what() << "\n";

                if (fitter)   delete fitter;
                if (fitTrack) delete fitTrack;  // Track owns rep; only delete rep if Track wasn't created
                else          delete rep;
            }
            catch (std::exception& e) {
                std::cerr << "STD exception at entry " << iev
                          << " eventID " << eventID
                          << " trackID " << tid
                          << ": " << e.what() << "\n";

                if (fitter)   delete fitter;
                if (fitTrack) delete fitTrack;
                else          delete rep;
            }

            if (nDumped < nDebugDump) {
                dumpEventDebug(out, hits, meas, bestP, signedBdl_Tm);
                nDumped++;
            }

            tout->Fill();
        }

        if ((iev + 1) % 100 == 0) {
            std::cout << "  processed events " << (iev + 1)
                      << " / " << nEntries
                      << " | muon rows " << nMuonRows
                      << " | fits OK " << nFitOK << "\n";
        }
    }

    fout->cd();
    tout->Write();
    fout->Close();

    fin->Close();

    std::cout << "\nSaved output: " << outputFile << "\n";
    std::cout << "Muon rows: " << nMuonRows << "\n";
    std::cout << "Fits OK:   " << nFitOK << "\n";
    std::cout << "Done.\n";

    return 0;
}
