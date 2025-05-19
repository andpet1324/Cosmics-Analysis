// filename: CosmicsAnalysis.h

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
#include <TStyle.h> // used for gstyle
#include <TColor.h> // used for setting colors
#include <TROOT.h>
#include <TH2F.h>
#include <TPolyLine.h> // Used to draw the fitted line
#include <TPolyLine3D.h> // Used to draw the fitted line
#include <TDirectory.h>
#include <TMinuit.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <TLatex.h>
#include <Fit/Fitter.h>
#include <TF1.h>
#include <TF1Convolution.h>
#include <TMath.h>
#include <Math/Functor.h>
#include <Math/Vector3D.h>

//----------------//
// LDMX Framework //
//----------------//

#include "Framework/EventProcessor.h"
#include "Ecal/Event/EcalHit.h"
#include "Hcal/Event/HcalHit.h"
#include "SimCore/Event/SimParticle.h"
#include "SimCore/Event/SimCalorimeterHit.h"
// #include "DetDescr/DetectorID.h"
// #include "DetDescr/HcalDigiID.h"
#include "DetDescr/HcalGeometry.h"


// Notes on the code ------------------------------------------------------------------------------------------------

// Right now, the script throws events where there are hits in the side-ecal and hcal. This results in a bad analysis.
// For example, the azimuthal angle distribution is not homogenuous, since events at large angles aren't included.
// An implementation of the side-hcal and ecal is needed. For my analysis, we only cared about the back-hcal, and 
// results from the ecal could be misleading. For a side-hcal analysis, the ecal will have an effect on the results, and
// should be considered.

// ------------------------------------------------------------------------------------------------------------------

// Struct for information that we save in the detector map.
struct DetectorInfo {
  // Scintillator bar geometry
  int xMin;
  int xMax;
  float yMin;
  float yMax;
  int zMin;
  int zMax;
  // orientation of the bar (horizontal or vertical)
  std::string orientation;
  // energy deposited in the bar
  double energy;
  // number of hits in the bar
  int numberOfHits;
  // estimated energy deposited in the bar
  double estimatedEnergy;
};


// Defines CosmicsAnalysis as a class that is based on framework::Analyzer
// Analyzer processes a constant copy of the event which cannot be updated.
class CosmicsAnalysis : public framework::Analyzer {
  public:
    CosmicsAnalysis(const std::string& name, framework::Process& p) : framework::Analyzer(name, p) {}
    ~CosmicsAnalysis() override = default;
    void onProcessStart() override;
    void analyze(const framework::Event& event) override;
    void onProcessEnd() override; 
    // angle funcs
    void getAngle(const std::map<int, ldmx::SimParticle, std::less<int>, std::allocator<std::pair<const int, ldmx::SimParticle> > >& sim_particles, const std::vector<ldmx::SimCalorimeterHit, std::allocator<ldmx::SimCalorimeterHit>>& sim_hits);
    void getAngleEstimate(double& polarAngleEstimate, double& azimuthalAngleEstimate, const std::vector<float>& line);
    // plotting funcs
    void getPlot2D(std::pair<std::string, std::string> pair);
    void getPlot1D(std::pair<std::string, std::string> pair);
    // detector map funcs
    void buildMap(std::map<int,DetectorInfo>& map);
    void addToMap(std::map<int,DetectorInfo>& map, std::vector<float> position, double edep);
    void AddDaughterEnergy(std::map<int,DetectorInfo>& map, const std::vector<double>& vertex, const std::vector<double>& endPoint, const double& edep);
    bool IsInBar(std::vector<float> position);
    void AddEstimatedEnergy(std::map<int,DetectorInfo>& map, const std::vector<float>& position, const double& stepsize);
    // fitting funcs
    void lineFunc(double t, const double *p, double &x, double &y, double &z);
    void curveFunc(double t, const double *p, double &x, double &y, double &z);
    TVectorD fitSurfaceInitialGuess(TGraph2D* graph);

  private:
    // Declare files, trees and graphs
    TFile *file_;
    TFile *graphFile_;
    TTree *sim_muon_tree_;
    TTree *sim_daughter_tree_;
    TTree *hcal_hit_tree_;
    TGraph2D *graph_;
    TCanvas *canvas_;
    
    // Muon variables for storing quantities
    int muon_size_;
    std::vector<float> muon_origin_;
    double muon_energy_;
    double muon_deviation_angle_;
    double muon_polar_angle_;
    double muon_azimuthal_angle_;
    double muon_kinetic_;
    double muon_init_mom_;
    double muon_final_mom_;
    double muon_mom_diff_;
    double muon_mass_;
    double muon_mom_ratio_;
    std::vector<double> muon_edeps_;
    std::vector<double> muon_est_edeps_;
    std::vector<double>* edeps_vec_ptr_ = &muon_edeps_;
    std::vector<double>* est_edeps_vec_ptr_ = &muon_est_edeps_;
    double muon_edeps_total_event_;
    Int_t muon_bars_hit_;
    double muon_edeps_per_bar_;
    float muon_edeps_per_barmm_;
    double muon_est_polar_angle_;
    std::vector<double> edeps_est_edeps_ratio_;
    std::vector<double>* est_edeps_ratio_ptr_ = &edeps_est_edeps_ratio_;

    // Daughter particle quantities
    int daughter_pdgID_;
    double daughter_energy_;
    double daughter_mass_;
    double daughter_init_mom_;
    double daughter_final_mom_;
    double daughter_kinetic_;

    // Hcal sim quantities
    double hcal_hit_energy_;
    std::vector<float> recPE_;

    // Hcal map
    std::map<int, DetectorInfo> detectorMap_;


    // Other
    int nr_events_w_hadrons_ = 0;
    int nr_events_miss_ = 0;
    int nr_events_ecal_ = 0;
    int decayed_muons_ = 0;
    int total_events_ = 0;
    int straight_tracks_ = 0;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    int bar_orientation_;
    double largest_polar_angle_estimate_deviation_ = 0;
    int specific_event_;

    // Thresholds for defining straight tracks. The edeps thresholds weren't used in the end
    float edeps_upper_threshold_ = 10000;
    float edeps_lower_threshold_ = 0;
    float angle_threshold_ = 1;

    // If we want to run true, otherwise false
    // for 3D graphs
    bool graphcheck = false;
    // for line-fits
    bool estimatedcheck = false; // this requires graphcheck = true to run
  
    // bar orientation analysis
    int chose_orientation_ = 2; //0 for horizontal, 1 for vertical, 2 for both

    // back-hcal dimensions
    const int hcalZMax_ = 5544;
    const int hcalZMin_ = 840;
    const int hcalXMin_ = -1000;
    const int hcalXMax_ = 1000;
    const float hcalYMin_ = -1000 + 19.05; // added support box shift
    const float hcalYMax_ = 1000 + 19.05; // added support box shift
};