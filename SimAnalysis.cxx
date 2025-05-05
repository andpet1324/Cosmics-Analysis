// filename: SimAnalysis.cxx

//----------------//
//   C++ StdLib   //
//----------------//
#include <iostream>
#include <string>
#include <map>
#include <vector>

//----------//
//   ROOT   //
//----------//

#include "TTree.h"
#include "TFile.h"
#include "TObject.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TView3D.h"
#include "TView.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TLinearFitter.h"
#include "TVector.h"

//----------------//
// LDMX Framework //
//----------------//

#include "Framework/EventProcessor.h"
#include "Ecal/Event/EcalHit.h"
#include "Hcal/Event/HcalHit.h"
#include "SimCore/Event/SimParticle.h"
#include "SimCore/Event/SimCalorimeterHit.h"

//----------------//
//   My Classes   //
//----------------//

// #include "MuonData.h"

// Defines SimAnalysis as a class that is based on framework::Analyzer
// Analyzer processes a constant copy of the event which cannot be updated.
class SimAnalysis : public framework::Analyzer {
  public:
    SimAnalysis(const std::string& name, framework::Process& p) : framework::Analyzer(name, p) {}
    ~SimAnalysis() override = default;
    void onProcessStart() override; // Empty by default. Used to declare things like histograms when events start
    void analyze(const framework::Event& event) override; // Process the event and make histograms or summaries.
    void onProcessEnd() override; // Used to write the TTree
    double getAngle(ldmx::SimParticle particle);
    std::pair<std::vector<float>,float> linearFit(const std::vector<ldmx::SimCalorimeterHit, std::allocator<ldmx::SimCalorimeterHit>>& sim_particles);
    std::vector<float> crossProduct(std::vector<float> vec1, std::vector<float> vec2);

  private:
    // Declare file and tree
    TFile *fileGraph;
    TTree *tree;
    TGraph2D *graph;
    TH3D *th3;
    TH2D *th2;

    TLinearFitter *lf;

    int pdgID;
    int decay_count = 0;

    std::map<const int, int> color_list;

    float edeps_upper_threshold_ = 10;
    float edeps_lower_threshold_ = 0.5;
    float angle_threshold_ = 2;

    // Muon variables
    double muon_energy_;
    double muon_angle_;
    double muon_kinetic_;
    double muon_init_mom_;
    double muon_final_mom_;
    double muon_mass_;
    double muon_mom_ratio_;

    // Daughter particle variables
    int daughter_pdgID_;
    double daughter_energy_;
    double daughter_mass_;
    double daughter_init_mom_;
    double daughter_final_mom_;
    double daughter_kinetic_;

    // Other
    int nr_events_w_hadrons_ = 0;
    std::vector<int> nr_events_miss_;
    int decayed_muons_ = 0;
    int total_events_ = 0;
    int straight_tracks_ = 0;

};

double SimAnalysis::getAngle(ldmx::SimParticle particle) {
  auto start_point = particle.getVertex();
  auto end_point = particle.getEndPoint();
  auto start_mom = particle.getMomentum();
  auto end_mom = particle.getEndPointMomentum();
  muon_init_mom_ = std::hypot(start_mom[0], start_mom[1], start_mom[2]);
  muon_final_mom_ = std::hypot(end_mom[0], end_mom[1], end_mom[2]);
  muon_mom_ratio_ = muon_final_mom_/muon_init_mom_;
  std::vector<double> vec_r = {(end_point[0] - start_point[0]), (end_point[1] - start_point[1]), (end_point[2] - start_point[2])};
  double mag_vec_r = std::hypot(vec_r[0], vec_r[1], vec_r[2]); 
  muon_angle_ = acos((vec_r[0]*start_mom[0] + vec_r[1]*start_mom[1] + vec_r[2]*start_mom[2])/(mag_vec_r * muon_init_mom_)) * 180/M_PI;

  return muon_angle_;
}

std::vector<float> SimAnalysis::crossProduct(std::vector<float> vec1, std::vector<float> vec2) {

  // Calculates the cross product of two vectors

  float v1 = vec1[0]; float v2 = vec1[1]; float v3 = vec1[2];
  float b1 = vec2[0]; float b2 = vec2[1]; float b3 = vec2[3];

  std::vector<float> product = {((v2*b3) - (v3*b2)), ((v3*b1) - (v1*b3)), ((v1*b2) - (v2*b1))};

  return product;
}

std::pair<std::vector<float>,float> SimAnalysis::linearFit(const std::vector<ldmx::SimCalorimeterHit, std::allocator<ldmx::SimCalorimeterHit>>& sim_particles) {
  // First and last position of the muon
  const auto origin = sim_particles[0].getPosition();
  std::vector<float> best_line;
  float best_distance;

  for (int i = 1; i <= 10; i++) {
    const auto last_pos = (sim_particles[sim_particles.size() - i]).getPosition();

    // Get the line and line magnitude
    std::vector<float> line = {(last_pos[0] - origin[0]), (last_pos[1] - origin[1]), (last_pos[2] - origin[2])};
    float line_magnitude = std::hypot(line[0], line[1], line[2]);

    float max_distance = 0;
    for (const auto& hit : sim_particles) {
      const std::vector<float> edeps = hit.getEdeps();
      const std::vector<float> pos = hit.getPosition();
      
      // Selection criteria for points
      if (edeps[0] > edeps_upper_threshold_ && edeps[0] < edeps_lower_threshold_ && muon_angle_ > angle_threshold_){
        continue;
      }

      if (edeps[0] <= 0.000000001) continue; // I don't think we actually have values that are zero
      //                                      but this is here just in case.
      //                                      0.000000001 is 0.001 eV

      // vector from origin to point
      std::vector<float> vector = {(pos[0] - origin[0]), (pos[1] - origin[1]), (pos[2] - origin[2])};

      // Check the distance between line and points
      // The distance is the magnitude of the cross product of vector and line divided by the line magnitude.

      std::vector<float> product = crossProduct(vector, line);
      float product_magnitude = std::hypot(product[0], product[1], product[2]);

      float distance = product_magnitude/line_magnitude;
      if (distance > max_distance){
        max_distance = distance;
      }
    }

    std::cout << "max distance: " << max_distance << std::endl;

    if (max_distance == 0) continue;

    if (i == 1){
      best_distance = max_distance;
      best_line = line;
    }
    else if (max_distance < best_distance) {
      best_distance = max_distance;
      best_line = line;
    }
  std::cout << "best distance: " << best_distance << std::endl;
  }

  return {best_line, best_distance};
}



void SimAnalysis::onProcessStart() {
  getHistoDirectory(); // forget this -> silently no output histogram

  fileGraph = new TFile("sim_graph.root", "RECREATE");
  tree = new TTree("Data","data from simulated muon hits");

  th2 = new TH2D("h2", "2D Hit Plot; X [mm], Y [mm]", 100, -1500, 1500, 100, -1500, 1500);
  th2->SetDirectory(fileGraph);

}

void SimAnalysis::analyze(const framework::Event& event) {

    // lf = new TLinearFitter(3, "x ++ y ++ x*x*y*y");

    // Define graph and set bins to max
    graph = new TGraph2D("sim_graph.root");
    graph->SetNpx(500);
    graph->SetNpy(500);
    graph->SetMarkerColor(kRed);
    graph->SetTitle("Cosmics muon path; x-axis ; y-axis ; z-axis");

    // th3 = new TH3D("h3", "3D Hit Plot;X [mm];Y [mm];Z [mm]", 100, -1000, 1000, 100, -1000, 1000, 100, 0, 5000);
    // th3->SetDirectory(fileGraph);

    // Print event number for analysis
    int eventNumber = event.getEventNumber();
    // std::cout << "event number: " << eventNumber << std::endl;

    // Get sim data
    const auto& hcal_sim_particles{event.getCollection<ldmx::SimCalorimeterHit>("HcalSimHits")};
    const auto& sim_particles{event.getMap<int, ldmx::SimParticle>("SimParticles", "cosmics")};

    // Get the nr of particles in the event
    auto nr_particles = sim_particles.size();

    // Get the muon
    auto muon = sim_particles.at(1);

    // Iterate over every particle
    for (const auto& [trackID, particle]: sim_particles){

        // Gets the angle of the muon
        if (trackID == 1) 
        {
          getAngle(particle);
        }

        // Check that particle doesn't end in the detector!
        auto end_point = particle.getEndPoint();
        if ((end_point[0] < 1000 && end_point[0] > -1000) && (end_point[1] < 1000 && end_point[1] > -1000) && (end_point[2] < 5000 && end_point[2] > 0) && trackID == 1) 
        {

            decay_count += 1;
        }

        pdgID = particle.getPdgID();

        // Checks the weird particles
        if (pdgID > 100 || pdgID < -100) 
        {
        }

        // Plot the particles with different color depending on their pdgID

        // Assign first particle to the map as 1
        if (color_list.empty())
        {
          color_list.insert({static_cast<const int>(pdgID), 1});
        }

        // Take the last value of the map
        auto var = std::prev(color_list.end())->second;
        
        // Check if the particle is already assigned to the map.
        // If yes, then plot. If not, assign and plot
        if (color_list.find(pdgID) != color_list.end()) 
        {
          // Plot this particle to a specific color
        }
        else
        {
          // Add to color_list and plot
          color_list.insert({static_cast<const int>(pdgID), (std::prev(color_list.end())->second)+1});
          
        }
    }

    const std::pair<std::vector<float>,float> line_fit = linearFit(hcal_sim_particles);

    std::vector<float> line = line_fit.first;
    float distance = line_fit.second;

    std::cout << "final distance:" << distance << std::endl;

    // Iterate over the hcal hits
    Int_t index = 0;
    for (const auto& hit : hcal_sim_particles) {
        std::vector<float> pos = hit.getPosition();
        auto edeps = hit.getEdeps();
        Double_t x_hit = pos[0];
        Double_t y_hit = pos[1];
        Double_t z_hit = pos[2];

        // We don't care about products in the ecal or side-hcal for now
        if (std::abs(pos[0]) < 1000 && std::abs(pos[1]) < 1000 && pos[2] < 840) {
          delete graph;
          return;
        }

        // Now we set a threshold for the energy in each hit
        
        float muon_mom_threshold = 0.5;
        // edeps[0] < edeps_upper_threshold_ && edeps[0] > edeps_lower_threshold_ && muon_angle_ < angle_threshold_
        if (edeps[0] < edeps_upper_threshold_ && edeps[0] > edeps_lower_threshold_ && muon_angle_ < angle_threshold_){
          graph->SetPoint(index++, x_hit, y_hit, z_hit);
        }
    }

    graph->SetPoint(index, -1000, -1000, 0);
    graph->SetPoint(index+1, 1000, 1000, 5000);


    // if (graph->GetN() != 0) {
    fileGraph->WriteTObject(graph, Form("graph_%d", 1));
    // fileGraph->WriteTObject(th3, Form("th3_%d", 1));
    // }

    delete graph;

    std::cout << "Event " << eventNumber << std::endl;

}

void SimAnalysis::onProcessEnd(){
  
  for (const auto& [key, value] : color_list) {
    // std::cout << "Key: " << key << ", Value: " << value << std::endl;
  }

  fileGraph->Write();
    
}

DECLARE_ANALYZER(SimAnalysis);