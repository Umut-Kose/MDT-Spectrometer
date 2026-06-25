// =============================================================================
// multistation_reco.C
// ATLAS-MDT-style multi-station reconstruction
//
// MDT treatment:
//   - tubeCenterY/tubeCenterZ are geometry, not event-smeared
//   - measured quantity is drift radius r
//   - smear only trueDriftRadius at reco level
//   - fit uses drift-circle residual: |y_fit - y_wire| - r
//
// Usage:
//   root -l 'multistation_reco.C(eventNumber, smearSigma_mm)'
//
// Examples:
//   root -l 'multistation_reco.C(4, 0.080)'
//   root -l 'multistation_reco.C(4, 0.0817)'
//   root -l 'multistation_reco.C(0, 0.0)'
// =============================================================================

#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <random>

#include "TFile.h"
#include "TTree.h"
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


// ── broken-line model ────────────────────────────────────────────────────────
// Magnet centre Z positions [mm]
const double zm[3] = {-526., 0., 526.};

double trackY(double z, const double* p)
{
  double y = p[0] + p[1] * z;

  for (int m = 0; m < 3; m++) {
    if (z > zm[m]) y += p[m + 2] * (z - zm[m]);
  }

  return y;
}

void basisRow(double z, double* A)
{
  A[0] = 1.;
  A[1] = z;

  for (int m = 0; m < 3; m++) {
    A[m + 2] = (z > zm[m]) ? (z - zm[m]) : 0.;
  }
}

// ── manual 5×5 Gauss–Jordan solve (ATA·p = ATy) ─────────────────────────────
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
      if (fabs(M[r][col]) > fabs(M[mx][col])) mx = r;
    }

    if (fabs(M[mx][col]) < 1e-15) return false;

    for (int j = 0; j <= 5; j++) std::swap(M[col][j], M[mx][j]);

    for (int r = 0; r < 5; r++) {
      if (r == col) continue;

      double f = M[r][col] / M[col][col];

      for (int j = col; j <= 5; j++) {
        M[r][j] -= f * M[col][j];
      }
    }
  }

  for (int i = 0; i < 5; i++) {
    p[i] = M[i][5] / M[i][i];
  }

  return true;
}


// ── main ─────────────────────────────────────────────────────────────────────
void multistation_reco(int targetEvent=0, double smearSigma_mm=0.080) 
{
  TFile* f = TFile::Open("mdt_hits.root");

  if (!f || f->IsZombie()) {
    printf("ERROR: cannot open mdt_hits.root\n");
    return;
  }

  TTree* tree = (TTree*)f->Get("Hits");

  if (!tree) {
    printf("ERROR: tree 'Hits' not found in mdt_hits.root\n");
    f->Close();
    return;
  }

  Long64_t nEntries = tree->GetEntries();

  if (targetEvent < 0 || targetEvent >= nEntries) {
    printf("ERROR: requested event %d, but tree has %lld entries\n",
           targetEvent, nEntries);
    f->Close();
    return;
  }

  int eventID = -1;

  std::vector<double>* trueDriftRadius = nullptr;
  std::vector<double>* tubeCenterY     = nullptr;
  std::vector<double>* tubeCenterZ     = nullptr;
  std::vector<double>* trueY           = nullptr;
  std::vector<double>* trueZ           = nullptr;
  std::vector<double>* pz_branch       = nullptr;

  std::vector<int>* stationID = nullptr;
  std::vector<int>* planeID   = nullptr;
  std::vector<int>* tubeID    = nullptr;
  std::vector<int>* pdg       = nullptr;

  tree->SetBranchAddress("eventID",         &eventID);
  tree->SetBranchAddress("trueDriftRadius", &trueDriftRadius);
  tree->SetBranchAddress("tubeCenterY",     &tubeCenterY);
  tree->SetBranchAddress("tubeCenterZ",     &tubeCenterZ);
  tree->SetBranchAddress("trueY",           &trueY);
  tree->SetBranchAddress("trueZ",           &trueZ);
  tree->SetBranchAddress("stationID",       &stationID);
  tree->SetBranchAddress("planeID",         &planeID);
  tree->SetBranchAddress("tubeID",          &tubeID);
  tree->SetBranchAddress("pdg",             &pdg);

  if (tree->GetBranch("pz")) {
    tree->SetBranchAddress("pz", &pz_branch);
  }

  tree->GetEntry(targetEvent);

  if (!trueDriftRadius || !tubeCenterY || !tubeCenterZ ||
      !trueY || !trueZ || !stationID || !planeID || !pdg) {
    printf("ERROR: one or more required branches are missing or null\n");
    f->Close();
    return;
  }

  std::mt19937 rng(42);

  const double tubeInnerRadius = 14.6; // mm, ATLAS MDT inner radius approximately 14.6 mm

  auto smearedR = [&](double r_true) -> double {
    if (smearSigma_mm <= 0.) return r_true;

    std::normal_distribution<double> gaussR(0.0, smearSigma_mm);

    double r = r_true + gaussR(rng);

    if (r < 0.) r = 0.;
    if (r > tubeInnerRadius) r = tubeInnerRadius;

    return r;
  };

  printf("\n========================================================\n");
  printf("  ATLAS-MDT-style multi-station reconstruction\n");
  printf("  Event: %d\n", eventID);
  printf("  Drift-radius smearing sigma: %.1f um\n", smearSigma_mm * 1000.);
  printf("  NOTE: only drift radius is smeared; wire geometry is not smeared\n");
  printf("========================================================\n");

  struct Hit {
    double wy;
    double wz;
    double r;
    double rtrue;
    double ty;
    double tz;
    int sta;
    int pln;
  };

  std::map<std::pair<int, int>, int> seen;
  std::vector<Hit> hits;

  for (int i = 0; i < (int)stationID->size(); i++) {
    if (abs((*pdg)[i]) != 13) continue;

    auto key = std::make_pair((*stationID)[i], (*planeID)[i]);

    if (seen.count(key)) continue;
    seen[key] = 1;

    Hit h;
    h.wy    = (*tubeCenterY)[i];
    h.wz    = (*tubeCenterZ)[i];
    h.rtrue = (*trueDriftRadius)[i];
    h.r     = smearedR(h.rtrue);
    h.ty    = (*trueY)[i];
    h.tz    = (*trueZ)[i];
    h.sta   = (*stationID)[i];
    h.pln   = (*planeID)[i];

    hits.push_back(h);
  }

  std::sort(hits.begin(), hits.end(),
            [](const Hit& a, const Hit& b) {
              if (a.sta != b.sta) return a.sta < b.sta;
              return a.pln < b.pln;
            });

  int N = (int)hits.size();

  printf("  Hits collected: %d  expected: 12 for 4 stations x 3 planes\n\n", N);

  if (N < 5) {
    printf("ERROR: not enough hits for 5-parameter fit\n");
    f->Close();
    return;
  }

  if (N > 20) {
    printf("ERROR: too many hits for brute-force left/right scan: N=%d\n", N);
    f->Close();
    return;
  }

  const double sigma = (smearSigma_mm > 0.) ? smearSigma_mm : 0.001;

  printf("  Sta Pln  wireZ   wireY   trueR   measR   trueZ   trueY\n");
  printf("  --------------------------------------------------------\n");

  for (int h = 0; h < N; h++) {
    printf("   %d   %d  %7.1f %7.1f  %6.3f  %6.3f  %7.1f %7.3f\n",
           hits[h].sta,
           hits[h].pln,
           hits[h].wz,
           hits[h].wy,
           hits[h].rtrue,
           hits[h].r,
           hits[h].tz,
           hits[h].ty);
  }

  double bestChi2 = 1e99;
  double bestP[5] = {0., 0., 0., 0., 0.};
  int bestCombo = 0;

  int nComb = 1 << N;

  for (int combo = 0; combo < nComb; combo++) {
    double ATA[5][5] = {};
    double ATy[5] = {};

    // Linearized fit for this left/right assignment:
    // y_assigned = wireY +/- measured drift radius
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

    // ATLAS-MDT-like residual:
    // residual = distance(line, wire) - measured drift radius.
    //
    // In this 2D Y-Z model, wire direction is perpendicular to Y-Z plane,
    // so distance is approximated by |y_fit - wireY|.
    double chi2 = 0.;

    for (int h = 0; h < N; h++) {
      double yfit = trackY(hits[h].wz, p);
      double dist = fabs(yfit - hits[h].wy);
      double res  = (dist - hits[h].r) / sigma;

      chi2 += res * res;
    }

    if (chi2 < bestChi2) {
      bestChi2 = chi2;
      for (int i = 0; i < 5; i++) bestP[i] = p[i];
      bestCombo = combo;
    }
  }

  int ndof = N - 5;

  printf("\n  Best broken-line fit:\n");
  printf("    p0 intercept      = %+.4f mm\n", bestP[0]);
  printf("    p1 initial slope  = %+.6f rad\n", bestP[1]);
  printf("    k1 kick @ %+.0f mm = %+.6f rad\n", zm[0], bestP[2]);
  printf("    k2 kick @ %+.0f mm = %+.6f rad\n", zm[1], bestP[3]);
  printf("    k3 kick @ %+.0f mm = %+.6f rad\n", zm[2], bestP[4]);


  if (ndof > 0) {
    printf("  chi2 / Ndof = %.3f / %d = %.3f\n\n",
           bestChi2, ndof, bestChi2 / ndof);
  } else {
    printf("  chi2 = %.3f, Ndof = %d\n\n", bestChi2, ndof);
  }

  printf("  Sta Pln  side  wireZ   y_fit   trueY   MDT_res [um]   fit-true [um]\n");
  printf("  --------------------------------------------------------------------\n");

  double sumRes2 = 0.;
  int nRes = 0;

  for (int h = 0; h < N; h++) {
    int sign = ((bestCombo >> h) & 1) ? +1 : -1;

    double yfit = trackY(hits[h].wz, bestP);
    double dist = fabs(yfit - hits[h].wy);
    double mdtRes_um = (dist - hits[h].r) * 1000.;
    double truthRes_um = (yfit - hits[h].ty) * 1000.;

    printf("   %d   %d   %s  %7.1f  %7.3f  %7.3f  %+9.1f   %+9.1f\n",
           hits[h].sta,
           hits[h].pln,
           sign > 0 ? "L" : "R",
           hits[h].wz,
           yfit,
           hits[h].ty,
           mdtRes_um,
           truthRes_um);

    sumRes2 += truthRes_um * truthRes_um;
    nRes++;
  }

  printf("\n  RMS fit-true residual: %.1f um\n", sqrt(sumRes2 / nRes));

  // Reco momentum from fitted kicks.
  // p[GeV/c] = 0.3 * B[T] * L[m] / theta[rad]
  // Here: three 0.4 m magnets, B = 1.5 T.
  double totalBendReco = bestP[2] + bestP[3] + bestP[4];

  printf("\n  Fitted total bending angle: %+.6f rad\n", totalBendReco);

  double pReco = 0.;

  if (fabs(totalBendReco) > 1e-8) {
    pReco = 0.3 * 1.5 * 1.2 / fabs(totalBendReco);
    printf("  Momentum estimate from fitted kicks: p_reco = %.2f GeV/c\n", pReco);
  } else {
    printf("  Momentum estimate unavailable: bending angle too small\n");
  }

  int chargeSign = (totalBendReco > 0.) ? +1 : -1;
  printf("  Bend sign proxy: %+d\n\n", chargeSign);

  double pTrue = 0.;

  if (pz_branch && !pz_branch->empty()) {
    pTrue = fabs((*pz_branch)[0]) / 1000.;
  }

  // Drawing
  TCanvas* c = new TCanvas(Form("c_ev%d", targetEvent),
                           Form("ATLAS-MDT reco evt %d", eventID),
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
    Form("Event %d | ATLAS-MDT drift-circle fit | p_{true}=%.0f GeV/c | p_{reco}=%.0f GeV/c | #sigma_{r}=%.1f #mum;Z [mm];Y [mm]",
         eventID, pTrue, pReco, smearSigma_mm * 1000.),
    1, zlo, zhi,
    1, ylo, yhi);

  hf1->SetStats(0);
  hf1->GetYaxis()->SetTitleOffset(0.6);
  hf1->Draw();

  for (int m = 0; m < 3; m++) {
    TBox* mag = new TBox(zm[m] - 200., ylo, zm[m] + 200., yhi);
    mag->SetFillColorAlpha(kGray, 0.25);
    mag->SetLineColor(kGray + 1);
    mag->Draw("same");

    TLatex* mlbl = new TLatex(zm[m], ylo + 0.05 * (yhi - ylo), Form("M%d", m + 1));
    mlbl->SetTextAlign(21);
    mlbl->SetTextSize(0.04);
    mlbl->SetTextColor(kGray + 2);
    mlbl->Draw("same");
  }

  int scols[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1};

  for (int h = 0; h < N; h++) {
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

  leg1->AddEntry((TObject*)0, Form("#chi^{2}/ndf=%.1f/%d", bestChi2, ndof), "");
  leg1->AddEntry((TObject*)0, Form("p_{true}=%.0f GeV/c", pTrue), "");
  leg1->AddEntry((TObject*)0, Form("p_{reco}=%.0f GeV/c", pReco), "");
  leg1->AddEntry((TObject*)0, Form("#theta_{bend}=%.4f rad", totalBendReco), "");

  TLine* ll1 = new TLine();
  ll1->SetLineColor(kBlack);
  ll1->SetLineWidth(2);
  leg1->AddEntry(ll1, "Reco broken-line", "l");

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

    TH2F* hf = new TH2F(Form("hf_s%d_%d", sta, targetEvent),
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

    for (int h = 0; h < N; h++) {
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

  TString outname = Form("multistation_evt%d_ATLAS_MDT.png", targetEvent);
  c->SaveAs(outname);

  printf("  Plot saved: %s\n", outname.Data());

// =============================================================================
// Extra diagnostic plot: MDT hit-level reconstruction details
// =============================================================================

TCanvas* cdiag = new TCanvas(Form("cdiag_ev%d", targetEvent),
                             Form("MDT hit diagnostics evt %d", eventID),
                             1400, 900);

cdiag->Divide(4, 3, 0.003, 0.003);

for (int ih = 0; ih < N; ih++) {
  cdiag->cd(ih + 1);
  gPad->SetGrid();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);

  const Hit& h = hits[ih];

  int sign = ((bestCombo >> ih) & 1) ? +1 : -1;

  double yRecoHit = h.wy + sign * h.r;
  double yFit     = trackY(h.wz, bestP);

  double dist     = fabs(yFit - h.wy);
  double mdtRes   = dist - h.r;
  double truthRes = yFit - h.ty;

  double zmin = h.wz - 35.;
  double zmax = h.wz + 35.;
  double ymin = std::min({h.wy - 20., h.ty - 10., yRecoHit - 10., yFit - 10.});
  double ymax = std::max({h.wy + 20., h.ty + 10., yRecoHit + 10., yFit + 10.});

  TH2F* frame = new TH2F(Form("diag_frame_%d_%d", targetEvent, ih),
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

TString diagName = Form("multistation_evt%d_hitDiagnostics.png", targetEvent);
cdiag->SaveAs(diagName);

printf("  Hit diagnostic plot saved: %s\n", diagName.Data());


  f->Close();
}


  