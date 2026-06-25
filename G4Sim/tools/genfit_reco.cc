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

    int bestCombo;

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

    out.bestCombo = -1;
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

void dumpEventDebug(const RecoOut& out, const std::vector<Hit>& hits, const std::vector<MDTMeas>& meas, const double* bestP, double signedBdl_Tm, double smearSigma_mm
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

    TCanvas* c = new TCanvas(Form("c_ev%d", out.eventID),
                            Form("ATLAS-MDT reco evt %d", out.eventID),
                            1200, 950);

    c->Divide(1, 2, 0.001, 0.001);

    c->cd(1);
    gPad->SetGrid();
    gPad->SetLeftMargin(0.07);
    gPad->SetRightMargin(0.22);

    double zlo = -950.;
    double zhi =  950.;

    double ylo = hits[0].ty;
    double yhi = hits[0].ty;

    for (auto& h : hits) {
        ylo = std::min(ylo, h.ty - 30.);
        yhi = std::max(yhi, h.ty + 30.);
    }

    for (auto& h : hits) {
        ylo = std::min(ylo, h.wy - 25.);
        yhi = std::max(yhi, h.wy + 25.);
    }

    double ypad = 0.15 * (yhi - ylo);
    ylo -= ypad;
    yhi += ypad;

    TH2F* hf1 = new TH2F("hf1",
    Form("Event %d | ATLAS-MDT drift-circle fit | p_{true}=%.0f GeV/c | p_{analytic}=%.0f GeV/c | p_{reco}=%.0f GeV/c | #sigma_{r}=%.1f #mum;Z [mm];Y [mm]",
         out.eventID, out.pTrue, out.pAnalytic, out.pReco, smearSigma_mm * 1000.),
        1, zlo, zhi,
        1, ylo, yhi);

    hf1->SetStats(0);
    hf1->GetYaxis()->SetTitleOffset(0.6);
    hf1->Draw();

    for (int m = 0; m < 3; m++) {
        TBox* mag = new TBox(zm_fit[m] - 200., ylo, zm_fit[m] + 200., yhi);
        mag->SetFillColorAlpha(kGray, 0.25);
        mag->SetLineColor(kGray + 1);
        mag->Draw("same");

        TLatex* mlbl = new TLatex(zm_fit[m], ylo + 0.05 * (yhi - ylo), Form("M%d", m + 1));
        mlbl->SetTextAlign(21);
        mlbl->SetTextSize(0.04);
        mlbl->SetTextColor(kGray + 2);
        mlbl->Draw("same");
    }
    int scols[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1};

    for (int h = 0; h < (int)hits.size(); h++) {
        TEllipse* tube = new TEllipse(hits[h].wz, hits[h].wy, 15., 15.);
        tube->SetLineColor(kGray + 1);
        tube->SetFillStyle(0);
        tube->Draw("same");

        TEllipse* dc = new TEllipse(hits[h].wz, hits[h].wy, hits[h].r, hits[h].r);
        dc->SetLineColor(scols[hits[h].sta - 1]);
        dc->SetLineWidth(2);
        dc->SetFillStyle(0);
        dc->Draw("same");

        TMarker* wire = new TMarker(hits[h].wz, hits[h].wy, 20);
        wire->SetMarkerSize(0.8);
        wire->Draw("same");

        TMarker* tr = new TMarker(hits[h].tz, hits[h].ty, 29);
        tr->SetMarkerColor(kRed);
        tr->SetMarkerSize(1.2);
        tr->Draw("same");
    }

    for (int sta = 1; sta <= 4; sta++) {
        double z0 = 1e9;
        double z1 = -1e9;

        for (auto& h : hits) {
        if (h.sta == sta) {
            z0 = std::min(z0, h.wz - 50.);
            z1 = std::max(z1, h.wz + 50.);
        }
        }

        if (z0 > z1) continue;

        TLine* sl = new TLine(z0, trackY(z0, bestP), z1, trackY(z1, bestP));
        sl->SetLineColor(kBlack);
        sl->SetLineWidth(2);
        sl->Draw("same");
    }
    for (int sta = 1; sta <= 4; sta++) {
        double zwire = 0.;
        int nc = 0;

        for (auto& h : hits) {
            if (h.sta == sta) {
                zwire += h.wz;
                nc++;
            }
        }

        if (nc) {
        zwire /= nc;

        TLatex* sl = new TLatex(zwire, yhi - 0.08 * (yhi - ylo), Form("Sta %d", sta));
        sl->SetTextAlign(21);
        sl->SetTextSize(0.04);
        sl->SetTextColor(scols[sta - 1]);
        sl->Draw("same");
        }
    }

    TLegend* leg1 = new TLegend(0.79, 0.55, 0.99, 0.92);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);

    //leg1->AddEntry((TObject*)0, Form("#chi^{2}/ndf=%.1f/%d", out.chi2, out.ndf), "");
    leg1->AddEntry((TObject*)0, Form("p_{true}=%.0f GeV/c", out.pTrue), "");
    leg1->AddEntry((TObject*)0, Form("p_{analytic}=%.0f GeV/c", out.pAnalytic), "");
    leg1->AddEntry((TObject*)0, Form("p_{reco}=%.0f GeV/c", out.pReco), "");
    leg1->AddEntry((TObject*)0, Form("#theta_{bend}=%.4f rad", out.bendAnalytic), "");

    TLine* ll1 = new TLine();
    ll1->SetLineColor(kBlack);
    ll1->SetLineWidth(2);
    leg1->AddEntry(ll1, "Analytic", "l");

    TMarker* mm1 = new TMarker();
    mm1->SetMarkerStyle(29);
    mm1->SetMarkerColor(kRed);
    leg1->AddEntry(mm1, "True hit", "p");

    leg1->AddEntry((TObject*)0, "Ellipses = MDT drift radius", "");
    leg1->Draw();

    c->cd(2);

    TPad* sub = new TPad("sub", "sub", 0, 0, 1, 1);
    sub->Divide(4, 1, 0.003, 0.003);
    sub->Draw();

    for (int sta = 1; sta <= 4; sta++) {
        sub->cd(sta);
        gPad->SetGrid();
        gPad->SetLeftMargin(0.14);
        gPad->SetBottomMargin(0.18);

        std::vector<Hit> shits;

        for (auto& h : hits) {
        if (h.sta == sta) shits.push_back(h);
        }

        if (shits.empty()) continue;

        double zlo2 = shits.front().wz - 55.;
        double zhi2 = shits.back().wz + 55.;

        double yl2 = shits[0].wy;
        double yh2 = shits[0].wy;

        for (auto& h : shits) {
        yl2 = std::min(yl2, h.wy - 22.);
        yh2 = std::max(yh2, h.wy + 22.);
        }

        for (auto& h : shits) {
        yl2 = std::min(yl2, h.ty - 5.);
        yh2 = std::max(yh2, h.ty + 5.);
        }

        yl2 -= 5.;
        yh2 += 5.;

        TH2F* hf = new TH2F(Form("hf_s%d_%d", sta, out.eventID),
                            Form("Sta %d;Z [mm];Y [mm]", sta),
                            1, zlo2, zhi2,
                            1, yl2, yh2);

        hf->SetStats(0);
        hf->GetXaxis()->SetTitleSize(0.07);
        hf->GetYaxis()->SetTitleSize(0.07);
        hf->GetXaxis()->SetLabelSize(0.06);
        hf->GetYaxis()->SetLabelSize(0.06);
        hf->GetYaxis()->SetTitleOffset(0.8);
        hf->Draw();

        for (int h = 0; h < int(hits.size()); h++) {
            if (hits[h].sta != sta) continue;

            TEllipse* tube = new TEllipse(hits[h].wz, hits[h].wy, 15., 15.);
            tube->SetLineColor(kGray + 1);
            tube->SetFillStyle(0);
            tube->Draw("same");

            TEllipse* dc = new TEllipse(hits[h].wz, hits[h].wy, hits[h].r, hits[h].r);
            dc->SetLineColor(scols[hits[h].sta - 1]);
            dc->SetLineWidth(2);
            dc->SetFillStyle(0);
            dc->Draw("same");

            TMarker* wire = new TMarker(hits[h].wz, hits[h].wy, 20);
            wire->SetMarkerSize(0.8);
            wire->Draw("same");

            TMarker* truth = new TMarker(hits[h].tz, hits[h].ty, 29);
            truth->SetMarkerColor(kRed);
            truth->SetMarkerSize(1.4);
            truth->Draw("same");
            }

            TLine* seg = new TLine(zlo2,
                                trackY(zlo2, bestP),
                                zhi2,
                                trackY(zhi2, bestP));

            seg->SetLineColor(kBlack);
            seg->SetLineWidth(2);
            seg->Draw("same");

            for (auto& h : shits) {
            TLatex* pl = new TLatex(h.wz,
                                    yh2 - 0.10 * (yh2 - yl2),
                                    Form("P%d", h.pln));

            pl->SetTextAlign(21);
            pl->SetTextSize(0.065);
            pl->SetTextColor(scols[h.sta - 1]);
            pl->Draw("same");
            }
        }

        c->Update();
        TString outname = Form("multistation_evt%d_ATLAS_MDT.png", out.eventID);
        c->SaveAs(outname);
        printf("  Plot saved: %s\n", outname.Data());

        // =============================================================================
        // Extra diagnostic plot: MDT hit-level reconstruction details
        // =============================================================================

        TCanvas* cdiag = new TCanvas(Form("cdiag_ev%d", out.eventID),
                                    Form("MDT hit diagnostics evt %d", out.eventID),
                                    1400, 900);

        cdiag->Divide(4, 3, 0.003, 0.003);

        for (int ih = 0; ih < int(hits.size()); ih++) {
        cdiag->cd(ih + 1);
        gPad->SetGrid();
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);

        const Hit& h = hits[ih];

        int sign = ((out.bestCombo >> ih) & 1) ? +1 : -1;

        double yRecoHit = h.wy + sign * h.r;
        double yFit     = trackY(h.wz, bestP);

        double dist     = fabs(yFit - h.wy);
        double mdtRes   = dist - h.r;
        double truthRes = yFit - h.ty;

        double zmin = h.wz - 35.;
        double zmax = h.wz + 35.;
        double ymin = std::min({h.wy - 20., h.ty - 10., yRecoHit - 10., yFit - 10.});
        double ymax = std::max({h.wy + 20., h.ty + 10., yRecoHit + 10., yFit + 10.});

        TH2F* frame = new TH2F(Form("diag_frame_%d_%d", out.eventID, ih),
                                Form("Sta %d Plane %d;Z [mm];Y [mm]", h.sta, h.pln),
                                1, zmin, zmax,
                                1, ymin, ymax);

        frame->SetStats(0);
        frame->GetXaxis()->SetTitleSize(0.06);
        frame->GetYaxis()->SetTitleSize(0.06);
        frame->GetXaxis()->SetLabelSize(0.05);
        frame->GetYaxis()->SetLabelSize(0.05);
        frame->Draw();

    // Tube outer radius
    TEllipse* tube = new TEllipse(h.wz, h.wy, 15., 15.);
    tube->SetLineColor(kGray + 1);
    tube->SetLineStyle(2);
    tube->SetFillStyle(0);
    tube->Draw("same");

    // Measured drift radius
    TEllipse* drift = new TEllipse(h.wz, h.wy, h.r, h.r);
    drift->SetLineColor(kBlue + 1);
    drift->SetLineWidth(2);
    drift->SetFillStyle(0);
    drift->Draw("same");

    // Wire centre
    TMarker* wire = new TMarker(h.wz, h.wy, 20);
    wire->SetMarkerSize(1.0);
    wire->SetMarkerColor(kBlack);
    wire->Draw("same");

    // True Geant4 hit
    TMarker* truth = new TMarker(h.tz, h.ty, 29);
    truth->SetMarkerSize(1.5);
    truth->SetMarkerColor(kRed + 1);
    truth->Draw("same");

    // Reconstructed left/right assigned point
    TMarker* reco = new TMarker(h.wz, yRecoHit, 24);
    reco->SetMarkerSize(1.3);
    reco->SetMarkerColor(kBlue + 2);
    reco->Draw("same");

    // Fitted track segment
    TLine* fitLine = new TLine(zmin, trackY(zmin, bestP),
                                zmax, trackY(zmax, bestP));
    fitLine->SetLineColor(kBlack);
    fitLine->SetLineWidth(2);
    fitLine->Draw("same");

    // Visual line from wire to assigned reco hit
    TLine* radLine = new TLine(h.wz, h.wy, h.wz, yRecoHit);
    radLine->SetLineColor(kBlue + 1);
    radLine->SetLineStyle(3);
    radLine->SetLineWidth(2);
    radLine->Draw("same");

    // Text summary
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.045);

    lat.DrawLatex(0.18, 0.88, Form("r_{true}=%.3f mm", h.rtrue));
    lat.DrawLatex(0.18, 0.82, Form("r_{meas}=%.3f mm", h.r));
    lat.DrawLatex(0.18, 0.76, Form("side=%s", sign > 0 ? "L/+" : "R/-"));
    lat.DrawLatex(0.18, 0.70, Form("y_{reco}=%.3f mm", yRecoHit));
    lat.DrawLatex(0.18, 0.64, Form("y_{true}=%.3f mm", h.ty));
    lat.DrawLatex(0.18, 0.58, Form("y_{fit}=%.3f mm", yFit));
    lat.DrawLatex(0.18, 0.52, Form("MDT res=%+.1f #mum", mdtRes * 1000.));
    lat.DrawLatex(0.18, 0.46, Form("fit-true=%+.1f #mum", truthRes * 1000.));
    }

    cdiag->Update();
    TString diagName = Form("multistation_evt%d_hitDiagnostics.png", out.eventID);
    cdiag->SaveAs(diagName);

    printf("  Hit diagnostic plot saved: %s\n", diagName.Data());

}


int main(int argc, char** argv)
{
    const char* inputFile  = (argc > 1) ? argv[1] : "mdt_hits.root";
    const char* gdmlFile   = (argc > 2) ? argv[2] : "detector_geometry.gdml";
    int    targetEvent = (argc > 3) ? std::atoi(argv[3]) : 4;
    double smear_mm        = (argc > 4) ? std::atof(argv[4]) : 0.080;
    int useMaterial        = (argc > 5) ? std::atoi(argv[5]) : 1;

    std::cout << "\n============================================================\n";
    std::cout << "  Batch MDT GenFit reconstruction, per trackID\n";
    std::cout << "  Input:       " << inputFile << "\n";
    std::cout << "  GDML:        " << gdmlFile << "\n";
    std::cout << "  Target event: " << targetEvent << "\n";
    std::cout << "  sigma_r:     " << smear_mm * 1000.0 << " um\n";
    std::cout << "  Material:    " << (useMaterial ? "ON" : "OFF") << "\n";
    std::cout << "============================================================\n\n";

    // -------------------------------------------------------------------------
    // Geometry initialization
    // -------------------------------------------------------------------------
    TGeoManager::Import(gdmlFile);
    if (!gGeoManager) {
        std::cerr << "ERROR: could not import GDML file: " << gdmlFile << "\n";
        return 1;
    }

    gGeoManager->CloseGeometry();
    gGeoManager->SetVerboseLevel(0);

    // -------------------------------------------------------------------------
    // GenFit initialization
    // -------------------------------------------------------------------------
    genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
    genfit::MaterialEffects::getInstance()->setDebugLvl(0);

    if (useMaterial) {
        genfit::MaterialEffects::getInstance()->setNoEffects(false);
        genfit::MaterialEffects::getInstance()->setMscModel("Highland");
    } else {
        genfit::MaterialEffects::getInstance()->setNoEffects(true);
    }

    genfit::FieldManager::getInstance()->init(new MDTMagneticField());

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
    // Prepare output ROOT file
    // -------------------------------------------------------------------------

    TFile* fout = new TFile("tmp.root", "RECREATE");
    TTree* tout = new TTree("Reco", "Per-muon MDT GenFit reconstruction results");

    RecoOut out;
    resetOut(out);

    tout->Branch("eventID", &out.eventID, "eventID/I");
    tout->Branch("entry", &out.entry, "entry/I");
    tout->Branch("trackID", &out.trackID, "trackID/I");
    tout->Branch("muonIndex", &out.muonIndex, "muonIndex/I");
    tout->Branch("nMuons", &out.nMuons, "nMuons/I");

    tout->Branch("bestCombo", &out.bestCombo, "bestCombo/I");

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

    // -------------------------------------------------------------------------
    // Collect MDT hits
    // -------------------------------------------------------------------------

    const double tubeInnerRadius = 14.6;
    const double detRes_cm = (smear_mm > 0.0) ? smear_mm / 10.0 : 0.001;
    const double sigmaX_cm = 1.0;
    // smearRadius is defined per-muon inside the loop using a deterministic seed

    Long64_t nEntries = tree->GetEntries();

    Long64_t nMuonRows = 0;
    Long64_t nFitOK = 0;

    //for (Long64_t iev = 0; iev < nEntries; iev++) {
        tree->GetEntry(targetEvent);
            std::cout << "\nEntry " << targetEvent
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

        if (!trueDriftRadius || !tubeCenterY || !tubeCenterZ || !tubeCenterX ||
            !trueY || !trueZ || !stationID || !planeID || !pdg_br || !trackID_br) {
            return 1;
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

            out.entry = (int)targetEvent;
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

            /// pTrue: total momentum |p| = sqrt(px²+py²+pz²) of this trackID's first hit.
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

            printf("\n============================================================\n");
            printf("  ATLAS-MDT-like GenFit reconstruction\n");
            printf("  Event: %d\n", eventID);
            printf("  p_true: %.2f GeV/c\n", out.pTrue);
            printf("  PDG used in RKTrackRep: %d\n", pdg_muon);
            printf("  Drift-radius smearing sigma: %.1f um\n", smear_mm * 1000.);
            printf("  Hits collected: %d\n", N);
            printf("============================================================\n\n");

            printf("  MDT hits:\n");
            printf("  Sta Pln  wireZ   wireY   trueR   measR   trueZ   trueY\n");
            printf("  --------------------------------------------------------\n");

            for (const auto& h : hits) {
                printf("   %d   %d  %7.1f %7.1f  %6.3f  %6.3f  %7.1f %7.3f\n",
                    h.sta, h.pln, h.wz, h.wy, h.rtrue, h.r, h.tz, h.ty);
            }

            nMuonRows++;

            if (N < 5 || N > 20 || out.pTrue <= 0.0) {
                tout->Fill();
                continue;
            }
            
            // -------------------------------------------------------------------------
            // Resolve L/R ambiguity using analytic broken-line fit
            // -------------------------------------------------------------------------
            // Warn: brute-force LR enumeration is O(2^N * N).
            // For N=20 this is ~20M iterations per track — can be slow.
            if (N > 14) {
                std::cerr << "WARNING: entry " << targetEvent << " trackID " << tid
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

            out.bestCombo = bestCombo;
            out.lrOK = 1;
            out.bestLRChi2 = bestChi2;
            out.bestLRNdf = N - 5;

            printf("\n  Best broken-line fit:\n");
            printf("    p0 intercept      = %+.4f mm\n", bestP[0]);
            printf("    p1 initial slope  = %+.6f rad\n", bestP[1]);
            printf("    k1 kick @ %+.0f mm = %+.6f rad\n", zm_fit[0], bestP[2]);
            printf("    k2 kick @ %+.0f mm = %+.6f rad\n", zm_fit[1], bestP[3]);
            printf("    k3 kick @ %+.0f mm = %+.6f rad\n", zm_fit[2], bestP[4]);

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

            // -------------------------------------------------------------------------
            // Convert MDT drift circles into selected pseudo-hits for GenFit
            // -------------------------------------------------------------------------
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

            printf("\n  MDT pseudo-hits passed to GenFit:\n");
            printf("  Sta Pln side   z_meas   y_meas   y_true   r_meas\n");
            printf("  --------------------------------------------------\n");

            for (const auto& m : meas) {
                printf("   %d   %d   %s  %8.2f %8.3f %8.3f %8.3f\n",
                    m.sta,
                    m.pln,
                    m.side > 0 ? "L/+" : "R/-",
                    m.z_mm,
                    m.y_mm,
                    m.trueY_mm,
                    m.r_meas_mm);
            }
            // -------------------------------------------------------------------------
            // Build GenFit track
            // -------------------------------------------------------------------------
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

                printf("\n  GenFit seed:\n");
                printf("    pos_seed = (%+.4f, %+.4f, %+.4f) cm\n",
                    pos_seed.X(), pos_seed.Y(), pos_seed.Z());
                printf("    mom_seed = (%+.4f, %+.4f, %+.4f) GeV/c, |p|=%.3f\n",
                    mom_seed.X(), mom_seed.Y(), mom_seed.Z(), mom_seed.Mag());
                printf("  --------------------------------------------------\n");
                
                fitTrack = new genfit::Track(rep, pos_seed, mom_seed);

                int hitId = 0;
                int planeId = 0;

                for (const auto& m : meas) {
                    TVectorD hitCoords(2);
                    //hitCoords[0] = 0.0;
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
                // -------------------------------------------------------------------------
                // Run GenFit Kalman fit
                // -------------------------------------------------------------------------
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
                // -------------------------------------------------------------------------
                // Extract GenFit results
                // -------------------------------------------------------------------------

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
                std::cerr << "GenFit exception at entry " << targetEvent
                          << " eventID " << eventID
                          << " trackID " << tid
                          << ": " << e.what() << "\n";

                if (fitter)   delete fitter;
                if (fitTrack) delete fitTrack;  // Track owns rep; only delete rep if Track wasn't created
                else          delete rep;
            }
            catch (std::exception& e) {
                std::cerr << "STD exception at entry " << targetEvent
                          << " eventID " << eventID
                          << " trackID " << tid
                          << ": " << e.what() << "\n";

                if (fitter)   delete fitter;
                if (fitTrack) delete fitTrack;
                else          delete rep;
            }

            if (nDumped) {
                dumpEventDebug(out, hits, meas, bestP, signedBdl_Tm,smear_mm);
                nDumped++;
            }

            tout->Fill();
        }
    //}

    fout->cd();
    tout->Write();
    fout->Close();

    fin->Close();

    std::cout << "Muon rows: " << nMuonRows << "\n";
    std::cout << "Fits OK:   " << nFitOK << "\n";
    std::cout << "Done.\n";
    return 0;
}