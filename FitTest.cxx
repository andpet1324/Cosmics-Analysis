// filename: FitTest.cxx
//----------------//
//   C++ StdLib   //
//----------------//
#include <iostream> // Printing
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <chrono> // For time
#include <iomanip> // For setprecision
#include <stdexcept> // For throwing errors

//----------//
//   ROOT   //
//----------//

#include <TTree.h>
#include <TFile.h>
#include <TObject.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TLeaf.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TColor.h>
#include <TROOT.h>
#include <TBox.h>
#include <TH2F.h>
#include <TPolyLine3D.h> // Used to draw the fitted line
#include <TDirectory.h>
#include <TMinuit.h>
#include <TVector3.h>
#include <TList.h>
#include <TString.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <Fit/Fitter.h>
#include "TF1.h"
#include "TF1Convolution.h"
#include "TRandom.h"


//----------------//
// LDMX Framework //
//----------------//

#include "Framework/EventProcessor.h"
#include "Ecal/Event/EcalHit.h"
#include "Hcal/Event/HcalHit.h"
#include "SimCore/Event/SimParticle.h"
#include "SimCore/Event/SimCalorimeterHit.h"
#include "DetDescr/HcalGeometry.h"

class FitTest : public framework::Analyzer {
 public:
  FitTest(const std::string& name, framework::Process& p)
    : framework::Analyzer(name, p) {}
  ~FitTest() override = default;
  void onProcessStart() override;
  void analyze(const framework::Event& event) override;

  TFile *graphFile_;
  TCanvas *canvas_;
  TGraph2D *graph_;

  const int hcalZMax_ = 5544;
  const int hcalZMin_ = 840;
  const int hcalXMin_ = -1000;
  const int hcalXMax_ = 1000;
  const float hcalYMin_ = -1000 + 19.05; // added support box shift
  const float hcalYMax_ = 1000 + 19.05; // added support box shift
};

void FitTest::onProcessStart() {
    getHistoDirectory(); // forget this -> silently no output histogram
    
    graphFile_ = new TFile("fit_test_graphs.root", "RECREATE");
}

void FitTest::analyze(const framework::Event& event) {
  // this is where we will fill the histograms

  // Get the data collections from the event bus
  const auto& hcal_sim_hits{event.getCollection<ldmx::SimCalorimeterHit>("HcalSimHits")};
  const auto& sim_particles{event.getMap<int, ldmx::SimParticle>("SimParticles", "cosmics")};
  const auto& hcal_rec_hits{event.getCollection<ldmx::HcalHit>("HcalRecHits")};

  int eventNumber = event.getEventNumber();

    // Declare graph and canvas
    canvas_ = new TCanvas();
    graph_ = new TGraph2D();

    // graph_->SetDirectory(nullptr);

    graph_->SetMarkerColor(kCyan);
    graph_->SetMarkerStyle(20);
    graph_->SetTitle("Cosmics muon path; x-axis ; y-axis ; z-axis");

    std::vector<float> energies;

    // Iterate over the simulated hcal hits
    Int_t index = 0;
    int lastBar = 0;
    std::vector<std::vector<float>> currentHitPosition;
    std::vector<double> currentHitEnergy;
    for (const auto& hit : hcal_sim_hits) {
        std::vector<float> pos = hit.getPosition();
        auto hitEdeps = (hit.getEdeps())[0];
        float xHit = pos[0];
        float yHit = pos[1];
        float zHit = pos[2];

        double energyHit = 0;

        // // Filter away points that are outside the detector
        if (std::abs(pos[0]) > hcalXMax_ || pos[1] < hcalYMin_ || pos[1] > hcalYMax_ || pos[2] < hcalZMin_ || pos[2] > hcalZMax_) {
            continue;
        }

        // We don't care about products in the ecal or side-hcal for now
        // This is redundant as we actually use it in the sim-particle for-loop
        if (std::abs(pos[0]) < 1500 && std::abs(pos[1]) < 1519.05 && pos[2] < 840) {
            delete graph_;
            //   std::cout << "Event " << eventNumber << " has hits in the ECal or side-HCal and is being ignored" << std::endl;
            return;
        }

        // We want to sum the hits in a bar, then take an average position weighted by the energy
        // But let's just start with a normal average
        if (lastBar == 0) {
            currentHitPosition.push_back(pos);
            currentHitEnergy.push_back(hitEdeps);
            lastBar = hit.getID();
            continue;
        }
        if (hit.getID() == lastBar) {
            currentHitPosition.push_back(pos);
            currentHitEnergy.push_back(hitEdeps);
            continue;
        }
        else {
            // Get the average position
            float xSum = 0;
            float ySum = 0;
            float zSum = 0;
            double energySum = 0;
            for (int i = 0; i < currentHitPosition.size(); i++) {
                xSum += currentHitPosition[i][0] * currentHitEnergy[i];
                ySum += currentHitPosition[i][1] * currentHitEnergy[i];
                zSum += currentHitPosition[i][2] * currentHitEnergy[i];
            }

            for (const double energy : currentHitEnergy) {
                energySum += energy;
            }  
            energyHit = energySum/currentHitEnergy.size(); 

            xHit = xSum/energySum;
            yHit = ySum/energySum;
            zHit = zSum/energySum;

            // Clear for next bar
            currentHitPosition.clear();
            currentHitEnergy.clear();
            currentHitPosition.push_back(pos);
            currentHitEnergy.push_back(hitEdeps);

            lastBar = hit.getID();
        }
        // Last element
        if (index == hcal_sim_hits.size() && currentHitPosition.size() != 1){

            // Get the average position
            float xSum = 0;
            float ySum = 0; 
            float zSum = 0;
            double energySum = 0;

            for (int i = 0; i < currentHitPosition.size(); i++) {
                xSum += currentHitPosition[i][0] * currentHitEnergy[i];
                ySum += currentHitPosition[i][1] * currentHitEnergy[i];
                zSum += currentHitPosition[i][2] * currentHitEnergy[i];
            }

            for (const double energy : currentHitEnergy) {
                energySum += energy;
            }  
            energyHit = energySum/currentHitEnergy.size();  

            xHit = xSum/energySum;
            yHit = ySum/energySum;
            zHit = zSum/energySum;
        }

        // Use defined thresholds (see the constructor)
        graph_->SetPoint(index++, xHit, yHit, zHit);
        
    }

    if (graph_->GetN() != 0) {
        // Sets the corners of our detector. Just for visual representation
        graph_->SetPoint(index, -1000, -980.95, 0);
        graph_->SetPoint(index+1, 1000, 1019.15, 5544);

        // Draw and update to canvas
        graph_->Draw("p0");
        canvas_->Update();

        // Save the canvas
        graphFile_->cd();
        canvas_->Write(Form("canvas_%d", eventNumber));
    }

    // Delete so we don't get memory issues
    delete graph_;
    delete canvas_;
    
}

DECLARE_ANALYZER(FitTest);