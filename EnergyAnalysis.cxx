// filename: EnergyAnalysis.cxx

//----------------//
//   C++ StdLib   //
//----------------//
#include <iostream>
#include <string>

//----------//
//   ROOT   //
//----------//

#include <TTree.h>
#include <TFile.h>
#include <TObject.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TPolyLine3D.h>
#include <TView3D.h>
#include <TView.h>

//----------------//
// LDMX Framework //
//----------------//

#include "Framework/EventProcessor.h"
#include "Ecal/Event/EcalHit.h"
#include "Hcal/Event/HcalHit.h"

// Declare file and tree
TFile *file = new TFile("rec_data.root","RECREATE");
TTree *tree = new TTree("Muon Data","data from muon hits");
TFile *fileGraph = new TFile("graph.root", "RECREATE");


// Defines EnergyAnalysis as a class that is based on framework::Analyzer
// Analyzer processes a constant copy of the event which cannot be updated.
class EnergyAnalysis : public framework::Analyzer {
  public:
    EnergyAnalysis(const std::string& name, framework::Process& p) : framework::Analyzer(name, p) {}
    ~EnergyAnalysis() override = default;
    void onProcessStart() override; // Empty by default. Used to declare things like histograms when events start
    void analyze(const framework::Event& event) override; // Process the event and make histograms or summaries.
    void onProcessEnd() override; // Used to write the TTree

  private:
    // Declare branch variables
    double totalEnergy;
    double totalTime;

    std::vector<TH1D*> xHistograms;
    std::vector<TH1D*> yHistograms;
    std::vector<TH1D*> zHistograms;

    TGraph2D *graph;
};

void EnergyAnalysis::onProcessStart() {
  getHistoDirectory(); // forget this -> silently no output histogram

  // Creates a root 2D histogram of type TH2F. The F in TH2F classifies that it stores four bytes per cell (float).
  // Maximum precision 7 digits.
  histograms_.create(
    "total_hcal_rec_energy", "Total HCal Rec Energy [MeV]", 100, 0.0, 1000.0 // name, xlabel, bins, xmin, xmax
  );
  histograms_.create(
    "total_hcal_rec_time", "Total HCal Rec Time", 100, 0, 2200 
  );
  // Declare a branches for storing data
  tree->Branch("energy", &totalEnergy, "energy/D");
  tree->Branch("TOF", &totalTime, "TOF/D");
}

void EnergyAnalysis::analyze(const framework::Event& event) {

  // Declare a vector that contains the hits on the hcal
  // event.getCollection receives the data from the event bus, of type ldmx::HcalHit
  const auto& hcal_rec_hits{event.getCollection<ldmx::HcalHit>("HcalRecHits")};

  // hcal_rec_hits[0].Print(); // Prints data from a vector with HcalHit type

  //std::cout << "---------------------------------------------------------------------------------------------" << std::endl;

  // Calculates the total energy hits
  totalEnergy = 0.0;
  totalTime = 0.0;

  graph = new TGraph2D("graph.root","%lg %lg %lg",",");
  graph->SetMarkerColor(kRed);

  event.Print();

  int index = 0;
    for (const auto& hit : hcal_rec_hits) {
      totalEnergy += hit.getEnergy();
      totalTime += hit.getTime();

      graph->SetPoint(index, hit.getXPos(), hit.getZPos(), hit.getYPos());
      index++;

      //std::cout << "energy: " << hit.getEnergy() << " X-pos: " << hit.getXPos() << " Y-pos: " << hit.getYPos() << " Z-pos: " << hit.getZPos() << "ID: " << hit.getID() << std::endl;
    }

  // Fills the histogram with total val, i.e. the histogram will show the total energy deposited by the particle 
  // throughout the whole detector.
  if (totalEnergy != 0) { // If the total is zero the particle didn't pass through our detector or didn't interact, i.e. we don't care about it
    histograms_.fill("total_hcal_rec_energy", totalEnergy);
    histograms_.fill("total_hcal_rec_time", totalTime);

    tree->Fill();
    fileGraph->WriteTObject(graph, Form("graph_%d", 1));

  }


  int eventNumber = event.getEventNumber();
  std::cout << "event number: " << eventNumber << std::endl;

  delete graph;
}

void EnergyAnalysis::onProcessEnd(){
  // Fill a 3D histogram with all x,y,z values

  file->Write();
  fileGraph->Write();

  // Produce and save a canvas from data

  delete file;
  delete fileGraph;
}

DECLARE_ANALYZER(EnergyAnalysis);