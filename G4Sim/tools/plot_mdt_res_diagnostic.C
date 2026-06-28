#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLine.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <TString.h>

void plot_mdt_res_diagnostic() {
        // Apply standard high-energy physics layout profiles
    gStyle->SetOptStat(0);
    gStyle->SetGridStyle(3);
    gStyle->SetGridColor(kGray+1);

    TFile* file = TFile::Open("results_mdt.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: Cannot open input file results_mdt.root" << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("Reco");
    if (!tree) {
        std::cerr << "ERROR: Cannot find TTree 'Reco'" << std::endl;
        file->Close();
        return;
    }

    // Allocate internal branches memory addresses
    Double_t pReco, pErr, chi2, ndf, pTrue, qOverP, qOverPErr;
    Int_t fitOK, nHits, charge;

    tree->SetBranchAddress("pReco",     &pReco);
    tree->SetBranchAddress("pErr",      &pErr);
    tree->SetBranchAddress("chi2",      &chi2);
    tree->SetBranchAddress("ndf",       &ndf);
    tree->SetBranchAddress("pTrue",     &pTrue);
    tree->SetBranchAddress("nHits",     &nHits);
    tree->SetBranchAddress("charge",    &charge);
    tree->SetBranchAddress("qOverP",    &qOverP);
    tree->SetBranchAddress("qOverPErr", &qOverPErr);
    tree->SetBranchAddress("fitOK",     &fitOK);

    // Define custom variable-width energy bin array requested
    Double_t eBinEdges[] = {10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 200., 400., 600., 800., 1000.};
    Int_t nEBins = sizeof(eBinEdges) / sizeof(Double_t) - 1;

    // Define intervals vector for the distinct momentum distributions
    std::vector<std::pair<double, double>> energyIntervals;
    for (int e = 10; e < 110; e += 10) 
    {    
        energyIntervals.push_back({double(e-0.5), double(e + 5.5)});
        std::cout << "Added energy interval: " << e-0.5 << " - " << (e + 5.5) << " GeV" << std::endl;
    }

    for (int e = 200; e < 1200; e += 200) {
        energyIntervals.push_back({double(e-10), double(e + 10)});
        std::cout << "Added energy interval: " << e-10 << " - " << (e + 10) << " GeV" << std::endl;
    }
    int nSeparatedBins = energyIntervals.size();

    // -------------------------------------------------------------------------
    // Book the 8 Diagnostic Framework Panels
    // -------------------------------------------------------------------------
    // Panel 1: 2D Correlation Matrix Plot
    TH2D* h2_reco_vs_true = new TH2D("h2_reco_vs_true", "Reconstructed vs True Momentum;p_{True} [GeV];p_{Reco} [GeV]", 50, 0, 1050, 50, 0, 1050);

    // Panel 2: Global Momentum Residuals
    TH1D* h_p_residual = new TH1D("h_p_residual", "Absolute Momentum Residuals;p_{Reco} - p_{True} [GeV];Tracks", 80, -50, 50);

    // Panel 3: Absolute 1/p Curvature Residuals
    TH1D* h_invp_residual = new TH1D("h_invp_residual", "Absolute Curvature Residuals (1/p);1/p_{Reco} - 1/p_{True} [1/GeV];Tracks", 80, -0.02, 0.02);

    // Panel 4: Mathematical Curvature Pull Distribution with standard Gaussian fit
    TH1D* h_pull_curvature = new TH1D("h_pull_curvature", "Track Parameter Pull (Curvature Space);Pull [#Delta(1/p) / #sigma_{1/p}];Probability Density", 60, -4, 4);

    // Panel 5: Track Fit Quality Distribution
    TH1D* h_reduced_chi2 = new TH1D("h_reduced_chi2", "Track Reduced #chi^{2} Distribution;#chi^{2} / NDF;Tracks", 60, 0, 5);

    // Panel 6: Structural Apparatus Multiplicities
    TH1D* h_n_hits = new TH1D("h_n_hits", "Track Hit Multiplicity Distribution;Number of Precision Hits (N_{Hits});Tracks", 8, 5.5, 13.5);

    // Panel 7: FIXED - Custom Variable-Width Resolution TProfile utilizing requested energy steps
    TProfile* prof_res_vs_energy = new TProfile("prof_res_vs_energy", "Relative Momentum Resolution vs Energy;p_{True} [GeV];#sigma(p)/p [%]", nEBins, eBinEdges, "s");

    // Panel 8: Relative Resolution Profile vs Hit Multiplicity
    TProfile* prof_res_vs_hits = new TProfile("prof_res_vs_hits", "Momentum Resolution vs Hit Multiplicity;Number of Hits (N_{Hits});#sigma(p)/p [%]", 8, 5.5, 13.5, "s");

    // -------------------------------------------------------------------------
    // Book the 14 Separate Momentum Distributions 
    // -------------------------------------------------------------------------
    std::vector<TH1D*> hSeparatedBins(nSeparatedBins, nullptr);
    for (int idx = 0; idx < nSeparatedBins; ++idx) {
        double eMin = energyIntervals[idx].first;
        double eMax = energyIntervals[idx].second;
        double TruthBinCenter;
        if(eMin<100.0)
        TruthBinCenter = eMin + 0.5;
        else
        TruthBinCenter = eMin + 10;

        TString histName = Form("h_sep_bin_%d", idx);
        TString histTitle = Form("For True Muon Energy %.0f GeV;p_{Reco} [GeV];Tracks", TruthBinCenter);
        
        // Outlier mitigation: limit upper bound view to max * 2.0
        double xMaxView = eMax * 2.0;
        hSeparatedBins[idx] = new TH1D(histName, histTitle, 50, 0.0, xMaxView);
        
        int baseColor = (idx >= 9) ? kTeal : kAzure;
        hSeparatedBins[idx]->SetLineColor(baseColor + (idx % 3) - 1);
        hSeparatedBins[idx]->SetFillColorAlpha(baseColor + (idx % 3) - 1, 0.35);
        hSeparatedBins[idx]->SetLineWidth(2);
    }


    // -------------------------------------------------------------------------
    // Processing Loop to Populate Diagnostic Plots
    // -------------------------------------------------------------------------
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (fitOK != 1) continue;

        double delta_p = pReco - pTrue;
        double rel_res_percent = (TMath::Abs(delta_p) / pTrue) * 100.0;

        double curvature_reco = qOverP; 
        double curvature_true = charge / pTrue;
        double delta_curvature = curvature_reco - curvature_true;

        // Populate histograms
        h2_reco_vs_true->Fill(pTrue, pReco);
        h_p_residual->Fill(delta_p);
        h_invp_residual->Fill(delta_curvature);
        h_reduced_chi2->Fill(chi2 / ndf);
        h_n_hits->Fill(nHits);

        if (qOverPErr > 0) {
            double pull = delta_curvature / qOverPErr;
            h_pull_curvature->Fill(pull);
        }

        // Fill profiles to automatically map resolution configurations downstream
        //prof_res_vs_energy->Fill(pTrue, rel_res_percent);
        //prof_res_vs_hits->Fill(nHits, rel_res_percent);
        // If a track's reconstructed error is larger than 100%, it's an unconstrained outlier
        if (rel_res_percent < 100.0) {
            prof_res_vs_energy->Fill(pTrue, rel_res_percent);
            prof_res_vs_hits->Fill(nHits, rel_res_percent);
        }

        // Fill Separated Bins
        for (int idx = 0; idx < nSeparatedBins; ++idx) {
            if (pTrue >= energyIntervals[idx].first && pTrue < energyIntervals[idx].second) {
                hSeparatedBins[idx]->Fill(pReco);
                break;
            }
        }
    }

    // -------------------------------------------------------------------------
    // Instantiate Master Canvas Layout Graphics
    // -------------------------------------------------------------------------
    TCanvas* cDiag = new TCanvas("cDiag", "ATLAS-MDT Spectrometer Global Tracker Performance Dashboard", 2200, 1100);
    cDiag->Divide(4, 2, 0.01, 0.01); 

    int colorBlue = kAzure + 2;
    int colorCrimson = kRed + 2;
    int colorGreen = kTeal + 2;

    // --- Pad 1: 2D Correlation Matrix Plot ---
    cDiag->cd(1); gPad->SetGrid(); gPad->SetLogz();
    h2_reco_vs_true->Draw("COLZ");
    TLine* diagonal = new TLine(0, 0, 1000, 1000);
    diagonal->SetLineColor(kRed); diagonal->SetLineStyle(2); diagonal->SetLineWidth(2);
    diagonal->Draw("SAME");

    // --- Pad 2: Absolute Momentum Residuals ---
    cDiag->cd(2); gPad->SetGrid();
    h_p_residual->SetFillColorAlpha(colorBlue, 0.4); h_p_residual->SetLineColor(colorBlue); h_p_residual->SetLineWidth(2);
    h_p_residual->Draw("HIST");

    // --- Pad 3: Curvature Residuals (1/p Profile) ---
    cDiag->cd(3); gPad->SetGrid();
    h_invp_residual->SetFillColorAlpha(colorGreen, 0.4); h_invp_residual->SetLineColor(colorGreen); h_invp_residual->SetLineWidth(2);
    h_invp_residual->Draw("HIST");

    // --- Pad 4: Curvature Pull Matrix with Real-Time Gaussian Fitting ---
    cDiag->cd(4); gPad->SetGrid();
    h_pull_curvature->SetMarkerStyle(20); h_pull_curvature->SetMarkerSize(0.8);
    h_pull_curvature->Sumw2();
    h_pull_curvature->Fit("gaus", "Q"); 
    TF1* gausFit = h_pull_curvature->GetFunction("gaus");
    if (gausFit) {
        gausFit->SetLineColor(kBlack);
        gausFit->SetLineWidth(2);
    }
    h_pull_curvature->Draw("E1");

    // --- Pad 5: Track Reduced Chi2 ---
    cDiag->cd(5); gPad->SetGrid();
    h_reduced_chi2->SetLineColor(kOrange + 7); h_reduced_chi2->SetLineWidth(2);
    h_reduced_chi2->Draw("HIST");

    // --- Pad 6: Hit Multiplicities Histogram ---
    cDiag->cd(6); gPad->SetGrid();
    h_n_hits->SetFillColorAlpha(kGray+2, 0.5); h_n_hits->SetLineColor(kBlack); h_n_hits->SetLineWidth(2);
    h_n_hits->Draw("HIST");

    // --- Pad 7: Resolution Scaling Profile vs Custom Energy Steps ---
    cDiag->cd(7); 
    gPad->SetGrid();
    gPad->SetLogx(); // <-- ADD THIS LINE to distribute the 10-100 GeV bins beautifully
    
    prof_res_vs_energy->SetMarkerStyle(21); 
    prof_res_vs_energy->SetMarkerColor(colorBlue);
    prof_res_vs_energy->SetLineColor(colorBlue); 
    prof_res_vs_energy->SetLineWidth(2);
    
    // Explicitly set the axis range to match your boundaries perfectly
    prof_res_vs_energy->GetXaxis()->SetMoreLogLabels(); // Adds clean intermediate numbers (20, 30, 50...)
    prof_res_vs_energy->Draw("E1");

    // --- Pad 8: Resolution Scaling Profile vs Hits ---
    cDiag->cd(8); gPad->SetGrid();
    prof_res_vs_hits->SetMarkerStyle(22); prof_res_vs_hits->SetMarkerColor(colorGreen);
    prof_res_vs_hits->SetLineColor(colorGreen); prof_res_vs_hits->SetLineWidth(2);
    prof_res_vs_hits->Draw("E1");

    cDiag->SaveAs("spectrometer_tracker_diagnostics.pdf");

// -------------------------------------------------------------------------
    // CANVAS 2: Dedicated Momentum Distributions Landscape (4 Columns x 4 Rows)
    // -------------------------------------------------------------------------
    TCanvas* cBins = new TCanvas("cBins", "Reconstructed Momentum Performance per Separated Energy Pad", 1800, 1400);
    cBins->Divide(4, 4, 0.01, 0.01);

    for (int idx = 0; idx < nSeparatedBins; ++idx) {
        cBins->cd(idx + 1); // Sequentially draw from pad 1 to 14
        gPad->SetGrid();
        gPad->SetBottomMargin(0.14);
        gPad->SetLeftMargin(0.14);

        if (hSeparatedBins[idx]->GetEntries() > 0) {
            hSeparatedBins[idx]->Draw("HIST");
            
            // Set red dashed truth midpoint indicator lines
            double TruthBinCenter;
            if(energyIntervals[idx].first<100.0)
                TruthBinCenter = energyIntervals[idx].first + 0.5;
            else
                TruthBinCenter = energyIntervals[idx].first + 10;

            TLine* line = new TLine(TruthBinCenter, 0.0, TruthBinCenter, hSeparatedBins[idx]->GetMaximum());
            line->SetLineColor(kRed + 1); line->SetLineStyle(2); line->SetLineWidth(2);
            line->Draw("SAME");
        }
    }
    // Clean up empty remaining subplots in the final row of the 4x4 matrix canvas
    for (int j = nSeparatedBins + 1; j <= 16; ++j) {
        cBins->cd(j);
        gPad->SetBorderMode(0);
    }

    cBins->SaveAs("separated_energy_bins_performance.pdf");

}