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
// #include "DetDescr/DetectorID.h"
// #include "DetDescr/HcalDigiID.h"
#include "DetDescr/HcalGeometry.h"


//----------------//
//   My Classes   //
//----------------//

// First define the struct
struct DetectorInfo {
  int xMin;
  int xMax;
  float yMin;
  float yMax;
  int zMin;
  int zMax;
  std::string orientation;
  double energy;
  double estimatedEnergy;
};


// Defines CosmicsAnalysis as a class that is based on framework::Analyzer
// Analyzer processes a constant copy of the event which cannot be updated.
class CosmicsAnalysis : public framework::Analyzer {
  public:
    CosmicsAnalysis(const std::string& name, framework::Process& p) : framework::Analyzer(name, p) {}
    ~CosmicsAnalysis() override = default;
    void onProcessStart() override; // Empty by default. Used to declare things like histograms when events start
    void analyze(const framework::Event& event) override; // Process the event and make histograms or summaries.
    void onProcessEnd() override; // Used to write the TTree
    void getAngle(const std::map<int, ldmx::SimParticle, std::less<int>, std::allocator<std::pair<const int, ldmx::SimParticle> > >& sim_particles, const std::vector<ldmx::SimCalorimeterHit, std::allocator<ldmx::SimCalorimeterHit>>& sim_hits);
    void getAngleEstimate(double& polarAngleEstimate, double& azimuthalAngleEstimate, std::vector<float> line);
    void getPlot2D(std::pair<std::string, std::string> pair);
    void getPlot2D3Entries(std::string str1, std::string str2 , std::string str3);
    void getPlot1D(std::pair<std::string, std::string> pair);
    std::pair<std::vector<float>,float>  linearFit(const std::vector<ldmx::SimCalorimeterHit, std::allocator<ldmx::SimCalorimeterHit>>& sim_hits, const std::map<int, ldmx::SimParticle, std::less<int>, std::allocator<std::pair<const int, ldmx::SimParticle>>>& sim_particles);
    void linearFitR2(const std::vector<ldmx::SimCalorimeterHit, std::allocator<ldmx::SimCalorimeterHit>>& sim_hits, const std::map<int, ldmx::SimParticle, std::less<int>, std::allocator<std::pair<const int, ldmx::SimParticle>>>& sim_particles, std::vector<float>& a, double& b);
    // void chi2_fcn(Int_t& npar, Double_t*, Double_t& fval, Double_t* par, Int_t);
    std::vector<float> crossProduct(std::vector<float> vec1, std::vector<float> vec2);
    std::pair<std::vector<std::vector<float>>, std::vector<double>> weightedAverage(const std::vector<ldmx::SimCalorimeterHit, std::allocator<ldmx::SimCalorimeterHit>>& simHits);
    std::pair<std::vector<std::vector<float>>, std::vector<double>> sort(std::pair<std::vector<std::vector<float>>, std::vector<double>> hits);
    void line(double t, const double *p, double &x, double &y, double &z);
    bool isInDiscreteRange(int value);
    void buildMap(std::map<int,DetectorInfo>& map);
    void addToMap(std::map<int,DetectorInfo>& map, std::vector<float> position, double edep);
    void AddDaughterEnergy(std::map<int,DetectorInfo>& map, const std::vector<double>& vertex, const std::vector<double>& endPoint, const double& edep);
    bool IsInBar(std::vector<float> position);
    void AddEstimatedEnergy(std::map<int,DetectorInfo>& map, const std::vector<float>& position, const double& stepsize);

  private:
    // Declare file, tree and graph
    TFile *file_;
    TFile *graphFile_;
    TTree *sim_muon_tree_;
    TTree *sim_daughter_tree_;
    TTree *hcal_hit_tree_;
    TTree *rec_tree_;
    TGraph2D *graph_;
    TCanvas *canvas_;
    TPolyLine3D *polyline_fit_;
    
    // Muon variables
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

    // Daughter particle variables
    int daughter_pdgID_;
    double daughter_energy_;
    double daughter_mass_;
    double daughter_init_mom_;
    double daughter_final_mom_;
    double daughter_kinetic_;

    // Hcal sim variables
    double hcal_hit_energy_;
    std::vector<float> recPE_;

    // Hcal map
    std::map<int, int> hcalBars_;

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

    // Thresholds for defining straight tracks
    float edeps_upper_threshold_ = 10000;
    float edeps_lower_threshold_ = 0;
    float angle_threshold_ = 1;

    // If we want to run true, otherwise false
    // for line-fits
    bool linecheck = true;
    // for 3D graphs
    bool graphcheck = true;
    // bar orientation analysis
    int chose_orientation_ = 1; //0 for horizontal, 1 for vertical, 2 for both

    const int hcalZMax_ = 5544;
    const int hcalZMin_ = 840;
    const int hcalXMin_ = -1000;
    const int hcalXMax_ = 1000;
    const float hcalYMin_ = -1000 + 19.05; // added support box shift
    const float hcalYMax_ = 1000 + 19.05; // added support box shift

    float allowedLineDistance_ = 10;

    int currentLayer_;
    int currentStrip_;

    std::map<int, DetectorInfo> detectorMap_;

};