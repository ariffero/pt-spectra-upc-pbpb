// -----------------------------------------------------------------
// Headers
// -----------------------------------------------------------------

// C++ headers
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// ROOT headers
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

// RooFit headers
#include "RooAbsData.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGlobalFunc.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"

using namespace RooFit;
using namespace std;

// -----------------------------------------------------------------
// Read a TH2: X = pT, Y = DCAxy
// -----------------------------------------------------------------
TH2 *getDCAHistogram(TFile *file, string base, string particle, string category){
  string path = base + category + "_" + particle + "s";
  TH2 *h = dynamic_cast<TH2*>(file->Get(path.c_str()));
  if(!h){
    cout << "Missing TH2: " << path << endl;
  }
  return h;
}

// -----------------------------------------------------------------
// Check that data and MC have the same bin boundaries
// -----------------------------------------------------------------
bool sameAxis(TAxis *axis1, TAxis *axis2){
  if(axis1->GetNbins() != axis2->GetNbins()){
    return false;
  }
  for(int i = 1; i <= axis1->GetNbins() + 1; i++){
    if(abs(axis1->GetBinLowEdge(i) - axis2->GetBinLowEdge(i)) > 1.e-10){
      return false;
    }
  }
  return true;
}

// -----------------------------------------------------------------
// Project one pT bin, using only the selected DCA bins.
// Adjacent DCA bins are summed; their errors are added in quadrature.
// The caller owns the returned histogram.
// -----------------------------------------------------------------
TH1D *projectDCAHistogram(TH2 *h, string name, int ptBin, int firstDcaBin,
                         int lastDcaBin, int rebinFactor, const vector<double> &dcaEdges){
  TH1D *result = new TH1D(Form("%s_ptbin%d", name.c_str(), ptBin), h->GetTitle(),
                          static_cast<int>(dcaEdges.size() - 1), dcaEdges.data());
  result->SetDirectory(nullptr);
  result->GetXaxis()->SetTitle(h->GetYaxis()->GetTitle());
  result->GetYaxis()->SetTitle("Entries");

  int outputBin = 1;
  for(int first = firstDcaBin; first <= lastDcaBin; first += rebinFactor, outputBin++){
    int last = min(first + rebinFactor - 1, lastDcaBin);
    double content = 0.;
    double variance = 0.;
    for(int y = first; y <= last; y++){
      content += h->GetBinContent(ptBin, y);
      double error = h->GetBinError(ptBin, y);
      variance += error * error;
    }
    result->SetBinContent(outputBin, content);
    result->SetBinError(outputBin, sqrt(variance));
  }
  return result;
}

// -----------------------------------------------------------------
// Require nonnegative, finite contents and a positive integral
// -----------------------------------------------------------------
bool isUsableHistogram(TH1D *h){
  for(int i = 1; i <= h->GetNbinsX(); i++){
    if(!isfinite(h->GetBinContent(i)) || h->GetBinContent(i) < 0){
      return false;
    }
  }
  return h->Integral() > 0;
}

// -----------------------------------------------------------------
// Fit primary/secondary fractions in each existing pT bin.
// Proton secondary = weak decay; material is a separate category.
// -----------------------------------------------------------------
void fitDCASpectra(string particle = "pion", string detector = "TPC", int rebinFactor = 1){
  if((detector != "TPC" && detector != "TOF") || rebinFactor < 1){
    cerr << "Choose TPC or TOF and a positive rebin factor." << endl;
    return;
  }
  gROOT->SetBatch(kTRUE);
  const double requestedDcaMin = -0.2;
  const double requestedDcaMax = 0.2;

  if(particle != "pion" && particle != "kaon" && particle != "proton"){
    cerr << "Choose pion, kaon, or proton." << endl; return;
  }
  //-----Read the files for data and MC----------------------
  unique_ptr<TFile> fData(TFile::Open("AnalysisResults-spectra-data.root"));
  unique_ptr<TFile> fMC(TFile::Open("AnalysisResults-spectra-mc.root"));
  if(!fData || fData->IsZombie() || !fMC || fMC->IsZombie()){
    cerr << "Cannot open data or MC input." << endl; return;
  }
  //-----Read the histograms for the selected particle and PID------
  bool useMaterialSecondaries = particle == "proton";
  string base = "pt-spectra-inclusive-upc/DCAxy_" + detector + "_";
  TH2* data2 = getDCAHistogram(fData.get(), base, particle, "data");
  TH2* primary2 = getDCAHistogram(fMC.get(), base, particle, "primary");
  // All nonprimaries for pion/kaon; weak-decay secondaries ONLY for proton.
  TH2* secondary2 = getDCAHistogram(fMC.get(), base, particle, "secondary");
  TH2* material2 = useMaterialSecondaries ? getDCAHistogram(fMC.get(), base, particle, "material") : nullptr;
  if(!data2 || !primary2 || !secondary2 ||
      (useMaterialSecondaries && !material2)) return;
  for(TH2* h : vector<TH2*>{primary2, secondary2, material2}){
    if(h && (!sameAxis(data2->GetXaxis(), h->GetXaxis()) ||
              !sameAxis(data2->GetYaxis(), h->GetYaxis()))){
      cerr << "Input pT and DCA bin edges must match." << endl; return;
    }
  }
  //-----Select the DCA range before projection--------------------
  // Select complete source Y bins inside the requested DCA interval once.
  const TAxis* dcaAxis = data2->GetYaxis();
  int firstDcaBin = 1;
  while(firstDcaBin <= dcaAxis->GetNbins() &&
         dcaAxis->GetBinLowEdge(firstDcaBin) < requestedDcaMin - 1.e-10)
    ++firstDcaBin;
  int lastDcaBin = dcaAxis->GetNbins();
  while(lastDcaBin >= firstDcaBin &&
         dcaAxis->GetBinUpEdge(lastDcaBin) > requestedDcaMax + 1.e-10)
    --lastDcaBin;
  if(firstDcaBin > lastDcaBin){
    cerr << "No complete source DCA bins inside requested range." << endl;
    return;
  }
  // Group only selected bins. Keep a final smaller group if necessary.
  vector<double> dcaEdges;
  for(int i=firstDcaBin; i<=lastDcaBin; i+=rebinFactor)
    dcaEdges.push_back(dcaAxis->GetBinLowEdge(i));
  dcaEdges.push_back(dcaAxis->GetBinUpEdge(lastDcaBin));
  cout << "Selected DCA range [" << dcaEdges.front() << ", "
            << dcaEdges.back() << "] cm before projection; "
            << dcaEdges.size()-1 << " output DCA bins." << endl;
  //-----Prepare output file and fraction histograms---------------
  string outdir = "DCAxy_fit_" + particle + "_" + detector;
  gSystem->mkdir(outdir.c_str(), true);
  TFile output((outdir + "/results.root").c_str(), "RECREATE");
  if(output.IsZombie()) return;
  int nPtBins = data2->GetNbinsX();
  vector<double> edges(nPtBins+1);
  for(int i=0; i<=nPtBins; ++i) edges[i]=data2->GetXaxis()->GetBinLowEdge(i+1);
  TH1D hPrimaryFraction("primaryFraction", "Primary fraction;p_{T} (GeV/c);Fraction",
                       nPtBins, edges.data());
  TH1D hSecondaryFraction("secondaryFraction", "All secondary fraction;p_{T} (GeV/c);Fraction",
                         nPtBins, edges.data());
  TH1D hWeakFraction("weakDecayFraction", "Weak decay fraction;p_{T} (GeV/c);Fraction",
                    nPtBins, edges.data());
  TH1D hMaterialFraction("materialFraction", "Material fraction;p_{T} (GeV/c);Fraction",
                        nPtBins, edges.data());
  TH1D hStatus("fitStatus", "Fit status (-99 = empty/invalid, -98 = unsupported);p_{T} (GeV/c);Status",
               nPtBins, edges.data());
  TH1D hCovQual("covQual", "Covariance quality (-99 = skipped);p_{T} (GeV/c);Quality",
                nPtBins, edges.data());

  //-----Loop on the existing pT bins------------------------------
  for(int ptBin=1; ptBin<=nPtBins; ++ptBin){
    hStatus.SetBinContent(ptBin, -99);
    hCovQual.SetBinContent(ptBin, -99);
    // Restricted Y projection for one pT bin, with DCA rebinning.
    // Copy/sum only selected source bins, never the full DCA histogram.
    auto dataH = unique_ptr<TH1D>(projectDCAHistogram(data2, "data", ptBin, firstDcaBin, lastDcaBin, rebinFactor, dcaEdges));
    auto primaryH = unique_ptr<TH1D>(projectDCAHistogram(primary2, "primary", ptBin, firstDcaBin, lastDcaBin, rebinFactor, dcaEdges));
    auto secondaryH = unique_ptr<TH1D>(projectDCAHistogram(secondary2, "secondary", ptBin, firstDcaBin, lastDcaBin, rebinFactor, dcaEdges));
    unique_ptr<TH1D> materialH;
    if(useMaterialSecondaries) materialH = unique_ptr<TH1D>(projectDCAHistogram(material2, "material", ptBin, firstDcaBin, lastDcaBin, rebinFactor, dcaEdges));
    TH1D* hData = dataH.get();
    TH1D* hMCprimary = primaryH.get();
    TH1D* hMCsecondary = secondaryH.get();
    TH1D* hMCsecFromMaterial = materialH.get();
    // Save projections even when the subsequent checks reject the fit.
    output.cd();
    hData->Write(); hMCprimary->Write(); hMCsecondary->Write();
    if(useMaterialSecondaries) hMCsecFromMaterial->Write();

    if(!isUsableHistogram(hData) || !isUsableHistogram(hMCprimary) || !isUsableHistogram(hMCsecondary) ||
        (useMaterialSecondaries && !isUsableHistogram(hMCsecFromMaterial))){
      cerr << "Skipping empty/invalid pT bin " << ptBin << endl;
      continue;
    }

    //-----Check MC coverage in every populated data bin-------------
    int unsupportedBins = 0;
    double unsupportedEntries = 0.;
    for(int i=1; i<=hData->GetNbinsX(); ++i){
      const double templateContent = hMCprimary->GetBinContent(i)
          + hMCsecondary->GetBinContent(i)
          + (useMaterialSecondaries ? hMCsecFromMaterial->GetBinContent(i) : 0.);
      if(hData->GetBinContent(i)>0. && templateContent==0.){
        ++unsupportedBins;
        unsupportedEntries += hData->GetBinContent(i);
        cout << "UNSUPPORTED: pT bin " << ptBin << ", DCA ["
                  << hData->GetXaxis()->GetBinLowEdge(i) << ", "
                  << hData->GetXaxis()->GetBinUpEdge(i) << "] cm, data = "
                  << hData->GetBinContent(i) << endl;
      }
    }
    cout << "pT bin " << ptBin << ": unsupported DCA bins = "
              << unsupportedBins << ", data entries in them = "
              << unsupportedEntries << endl;
    if(unsupportedBins>0){
      // -98 distinguishes unsupported model coverage from empty/invalid input (-99).
      hStatus.SetBinContent(ptBin, -98);
      cerr << "Skipping pT bin " << ptBin
                << ": observed data outside MC template support. "
                << "Projections saved; consider coarser DCA binning or more MC."
                << endl;
      continue;
    }
    
    //-----Construct the binned RooFit data and MC templates----------

    // Fit the cropped DCA range; yields and fractions refer only to this interval.
    RooRealVar DCAxy("DCAxy", "DCA_{xy} (cm)",
                    hData->GetXaxis()->GetXmin(), hData->GetXaxis()->GetXmax());

    // Create RooDataHists
    RooDataHist *Hist_Data = NULL;
    RooDataHist *Hist_MCprimary = NULL;
    RooDataHist *Hist_MCsecondary = NULL;
    RooDataHist *Hist_MCsecFromMaterial = NULL;
    
    Hist_Data = new RooDataHist("Hist_Data", "Data DCAxy distribution", DCAxy, Import(*hData));
    Hist_MCprimary = new RooDataHist("Hist_MCprimary", "MC primary DCAxy distribution", DCAxy, Import(*hMCprimary));
    Hist_MCsecondary = new RooDataHist("Hist_MCsecondary", "MC secondary DCAxy distribution", DCAxy, Import(*hMCsecondary));
    if(useMaterialSecondaries){
      Hist_MCsecFromMaterial = new RooDataHist("Hist_MCsecFromMaterial", "MC secondary from material DCAxy distribution", DCAxy, Import(*hMCsecFromMaterial));
    }

    // Create MC templates
    RooHistPdf *MCprimary_PDF = new RooHistPdf("MCprimary_PDF", "PDF for primary particles", DCAxy, *Hist_MCprimary);
    RooHistPdf *MCsecondary_PDF = new RooHistPdf("MCsecondary_PDF", "PDF for secondary particles", DCAxy, *Hist_MCsecondary);
    RooHistPdf *MCsecFromMaterial_PDF = NULL;
    if(useMaterialSecondaries){
      MCsecFromMaterial_PDF = new RooHistPdf("MCsecFromMaterial_PDF", "PDF for secondary particles from material", DCAxy, *Hist_MCsecFromMaterial);
    }

    //-----Define component yields and build the model---------------
    double nData = Hist_Data->sumEntries();

    RooRealVar nPrimary(
        "nPrimary",
        "Number of primary particles",
        0.8 * nData,
        0.0,
        2.0 * nData
    );

    RooRealVar nSecondary(
        "nSecondary",
        "Number of secondary particles",
        (useMaterialSecondaries ? 0.1 : 0.2) * nData,
        0.0,
        2.0 * nData
    );

    RooRealVar nSecFromMaterial(
        "nSecFromMaterial",
        "Number of secondary particles from material",
        0.1 * nData,
        0.0,
        2.0 * nData
    );

    RooAddPdf *model = NULL;
    if(useMaterialSecondaries){
      model = new RooAddPdf("model3", "Model for DCAxy distribution",
                           RooArgList(*MCprimary_PDF, *MCsecondary_PDF, *MCsecFromMaterial_PDF),
                           RooArgList(nPrimary, nSecondary, nSecFromMaterial));
    }
    else{
      model = new RooAddPdf("model2", "Model for DCAxy distribution",
                           RooArgList(*MCprimary_PDF, *MCsecondary_PDF),
                           RooArgList(nPrimary, nSecondary));
    }

    // Fit the model to the data
    RooFitResult *fitResult = model->fitTo(*Hist_Data, Save(), Extended(kTRUE), SumW2Error(kTRUE));

    //-----Calculate fractions from fitted yields--------------------
    RooFormulaVar *primaryFraction = NULL;
    RooFormulaVar *secondaryFraction = NULL;
    RooFormulaVar *materialFraction = NULL;

    if(useMaterialSecondaries){
      primaryFraction = new RooFormulaVar("primaryFraction", "@0 / (@0 + @1 + @2)", RooArgList(nPrimary, nSecondary, nSecFromMaterial));
      secondaryFraction = new RooFormulaVar("secondaryFraction", "(@1 + @2) / (@0 + @1 + @2)", RooArgList(nPrimary, nSecondary, nSecFromMaterial));
      materialFraction = new RooFormulaVar("materialFraction", "@2 / (@0 + @1 + @2)", RooArgList(nPrimary, nSecondary, nSecFromMaterial));
    }
    else{
      primaryFraction = new RooFormulaVar("primaryFraction", "@0 / (@0 + @1)", RooArgList(nPrimary, nSecondary));
      secondaryFraction = new RooFormulaVar("secondaryFraction", "@1 / (@0 + @1)", RooArgList(nPrimary, nSecondary));
    }

    if(!fitResult){
      delete primaryFraction;
      delete secondaryFraction;
      delete materialFraction;
      delete model;
      delete MCsecFromMaterial_PDF;
      delete MCprimary_PDF;
      delete MCsecondary_PDF;
      delete Hist_MCsecFromMaterial;
      delete Hist_MCsecondary;
      delete Hist_MCprimary;
      delete Hist_Data;
      continue;
    }
    output.cd();
    fitResult->Write(Form("fitResult_ptbin%d", ptBin));
    hStatus.SetBinContent(ptBin, fitResult->status());
    hCovQual.SetBinContent(ptBin, fitResult->covQual());
    // Only converged fits with a fully accurate covariance enter the summaries.
    // Zero in other bins is a placeholder: always consult fitStatus/covQual.
    if(fitResult->status()==0 && fitResult->covQual()==3){
      hPrimaryFraction.SetBinContent(ptBin, primaryFraction->getVal());
      hPrimaryFraction.SetBinError(ptBin, primaryFraction->getPropagatedError(*fitResult));
      hSecondaryFraction.SetBinContent(ptBin, secondaryFraction->getVal());
      hSecondaryFraction.SetBinError(ptBin, secondaryFraction->getPropagatedError(*fitResult));
      if(useMaterialSecondaries){
        RooFormulaVar weakFraction("weakFraction", "@1/(@0+@1+@2)",
                                   RooArgList(nPrimary,nSecondary,nSecFromMaterial));
        hWeakFraction.SetBinContent(ptBin, weakFraction.getVal());
        hWeakFraction.SetBinError(ptBin, weakFraction.getPropagatedError(*fitResult));
        hMaterialFraction.SetBinContent(ptBin, materialFraction->getVal());
        hMaterialFraction.SetBinError(ptBin, materialFraction->getPropagatedError(*fitResult));
      }
    }
    cout << detector << " " << particle << ", pT bin " << ptBin
              << " [" << edges[ptBin-1] << ", " << edges[ptBin] << "] GeV/c" << endl;
    // Print the fit results
    fitResult->Print("v");

    cout << "\n=== Fitted yields and fractions ===" << endl;
    cout << "nData      = " << nData << endl;
    cout << "nPrimary   = " << nPrimary.getVal()
              << " +/- " << nPrimary.getError() << endl;
    cout << "nSecondary = " << nSecondary.getVal()
              << " +/- " << nSecondary.getError() << endl;
    if(useMaterialSecondaries){
      cout << "nSecFromMaterial = " << nSecFromMaterial.getVal()
                << " +/- " << nSecFromMaterial.getError() << endl;
    }
    cout << "Primary fraction   = " << primaryFraction->getVal()
              << " +/- " << primaryFraction->getPropagatedError(*fitResult)
              << endl;
    cout << "Secondary fraction = " << secondaryFraction->getVal()
              << " +/- " << secondaryFraction->getPropagatedError(*fitResult)
              << endl;
    if(useMaterialSecondaries){
      cout << "Material secondary fraction = " << materialFraction->getVal()
                << " +/- " << materialFraction->getPropagatedError(*fitResult)
                << endl;
    }

    // Plot the results
    RooPlot *frame = DCAxy.frame(Title(Form("%s %s: %.3g < p_{T} < %.3g GeV/c", particle.c_str(),
                 detector.c_str(), edges[ptBin-1], edges[ptBin])));
    Hist_Data->plotOn(frame, Name("Data"), DataError(RooAbsData::SumW2), LineColor(kBlack), MarkerStyle(20));
    model->plotOn(frame, Name("Total Fit"), LineColor(kRed));
    model->plotOn(frame, Components("MCprimary_PDF"), Name("Primary"), LineStyle(kDashed), LineColor(kGreen + 1));
    model->plotOn(frame, Components("MCsecondary_PDF"), Name("Secondary"), LineStyle(kDashed), LineColor(kBlue + 1));
    if(useMaterialSecondaries){
      model->plotOn(frame, Components("MCsecFromMaterial_PDF"), Name("Secondary from Material"), LineStyle(kDashed), LineColor(kMagenta + 1));
    }

    TCanvas *c = new TCanvas("c", "DCAxy Fit", 1920, 1080);
    frame->Draw();

    frame->SetMinimum(0.1);
    c->SetLogy();

    TLegend *legend = new TLegend(0.65, 0.65, 0.9, 0.9);
    legend->SetFillStyle(0);

    // Add entries for each component
    legend->AddEntry(frame->findObject("Data"), "Data", "P");
    legend->AddEntry(frame->findObject("Total Fit"), "Total Fit", "L");
    legend->AddEntry(frame->findObject("Primary"), "Primary component", "L");
    legend->AddEntry(frame->findObject("Secondary"), useMaterialSecondaries ? "Weak-decay component" : "All-secondary component", "L");
    if(particle == "proton"){
      legend->AddEntry(frame->findObject("Secondary from Material"), "Secondary from Material", "L");
    }
    legend->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.SetTextFont(42);
    latex.DrawLatex(0.12, 0.74, Form("Primaries = (%.3f #pm %.3f)%%", 100. * primaryFraction->getVal(), 100. * primaryFraction->getPropagatedError(*fitResult)));
    latex.DrawLatex(0.12, 0.69, Form("All secondaries = (%.3f #pm %.3f)%%", 100. * secondaryFraction->getVal(), 100. * secondaryFraction->getPropagatedError(*fitResult)));
    if(useMaterialSecondaries){
      latex.DrawLatex(0.12, 0.64, Form("Secondaries from material = (%.3f #pm %.3f)%%", 100. * materialFraction->getVal(), 100. * materialFraction->getPropagatedError(*fitResult)));
    }
    if(fitResult->status()!=0 || fitResult->covQual()!=3)
      latex.DrawLatex(0.12, 0.59, "Fit requires inspection: see fitStatus/covQual");
    c->SaveAs(Form("%s/fit_ptbin%03d.png", outdir.c_str(), ptBin));

    //-----Clean objects created for this pT bin---------------------
    delete c;
    delete legend;
    delete frame;
    delete primaryFraction;
    delete secondaryFraction;
    delete materialFraction;
    delete fitResult;
    delete model;
    delete MCsecFromMaterial_PDF;
    delete MCprimary_PDF;
    delete MCsecondary_PDF;
    delete Hist_MCsecFromMaterial;
    delete Hist_MCsecondary;
    delete Hist_MCprimary;
    delete Hist_Data;
  } // pT bin loop

  //-----Save fractions and fit diagnostics versus pT--------------
  output.cd();
  hPrimaryFraction.Write();
  hSecondaryFraction.Write();
  if(useMaterialSecondaries){
    hWeakFraction.Write();
    hMaterialFraction.Write();
  }
  hStatus.Write();
  hCovQual.Write();
  cout << "Saved results to " << outdir << "/results.root" << endl;
}
