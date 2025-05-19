// filename: CosmicsAnalysis.cxx

#include "CosmicsAnalysis.h"

// <--------------------------------------- SCINTILLATOR BAR MAP FUNCS ----------------------------------------------->

void CosmicsAnalysis::buildMap(std::map<int,DetectorInfo>& map) {

    // Function that builds the map of the scintillator bars. The map is a dictionary with the ID as key and the bar information as value.

    const float backHcaldz = 4704;
    const float sideHcalStartZ = 240;
    const float sideHcaldz = 600;
    const float back_hcal_startZ = sideHcaldz/2 - backHcaldz/2;

    const float airThick = 2;
    const float absoThick = 25;
    const float layerThick = 49;
    const float scintWidth = 50;
    const float scintThick = 20;

    const float supportBox = 19.05;

    const float hcal_mother_vol_center = ((backHcaldz + sideHcaldz) / 2) + sideHcalStartZ; // center of the hcal mother volume, including side and back-hcal. 240 is distance from target to start of side-hcal.
    // Positions with respect to hcal mother volume
    const float back_hcal_startZAbso = back_hcal_startZ + airThick + absoThick/2; // Center of the first absorber with respect to center of hcal mother volume
    const float back_hcal_startZScint = back_hcal_startZ + 2*airThick + absoThick + scintThick/2; // Center of the first absorber with respect to center of hcal mother volume
    // Positions with respect to target, i.e. z=0
    const float startZ = hcal_mother_vol_center + back_hcal_startZ; // Hcal start = 840 mm
    const float absoStartZ = hcal_mother_vol_center + back_hcal_startZAbso; // the center of the absorber
    const float horizontalScintStartZ = hcal_mother_vol_center + back_hcal_startZScint; // center of the first horizontal bars
    const float verticalScintStartZ = hcal_mother_vol_center + back_hcal_startZScint + layerThick; // center of the first vertical bars

    const int layers = 96;
    const int strips = 40;
    
    for (int i = 1; i <= layers; i += 2) { // Iterate through layers

        const float horizontalLayerPos = horizontalScintStartZ + (i-1)*layerThick;
        const float verticalLayerPos = verticalScintStartZ + (i-1)*layerThick;

        const int horizontalZMin = horizontalLayerPos - scintThick/2;
        const int horizontalZMax = horizontalLayerPos + scintThick/2;
        const int verticalZMin = verticalLayerPos - scintThick/2;
        const int verticalZMax = verticalLayerPos + scintThick/2;

        for (int l = 1; l <= strips; l++) { // Iterate through strips (scints) Strips go from -1000 to 1000 in x & y.

            const float horizontalStripPos = -975 + supportBox + (l-1)*scintWidth;
            const int verticalStripPos = -975 + (l-1)*scintWidth;

            const int horizontalXMin = -1000;
            const int horizontalXMax = 1000;
            const float horizontalYMin = horizontalStripPos - scintWidth/2;
            const float horizontalYMax = horizontalStripPos + scintWidth/2;
            const int verticalXMin = verticalStripPos - scintWidth/2;
            const int verticalXMax = verticalStripPos + scintWidth/2;
            const float verticalYMin = -1000 + supportBox;
            const float verticalYMax = 1000 + supportBox;

            std::string horizontalIDString = std::string("10") + (i < 10 ? "0" : "") + std::to_string(i) + (l < 10 ? "0" : "") + std::to_string(l); // 10 for horizontal + layer + strip
            std::string verticalIDString = std::string("20") + ((i+1) < 10 ? "0" : "") + std::to_string(i+1) + (l < 10 ? "0" : "") + std::to_string(l); // + 20 for vertical + (layer+1) + strip
            DetectorInfo horizontalData = {horizontalXMin, horizontalXMax, horizontalYMin, horizontalYMax, horizontalZMin, horizontalZMax, "horizontal" , 0.0, 0, 0.0};
            DetectorInfo verticalData = {verticalXMin, verticalXMax, verticalYMin, verticalYMax, verticalZMin, verticalZMax, "vertical" , 0.0, 0, 0.0};

            const int horizontalID = std::stoi(horizontalIDString);
            const int verticalID = std::stoi(verticalIDString);

            map.insert({horizontalID, horizontalData});
            map.insert({verticalID, verticalData});
        }
    }
}

void CosmicsAnalysis::addToMap(std::map<int,DetectorInfo>& map, std::vector<float> position, double edep) {
    // get detector id by position
    const float x = position[0];
    const float y = position[1];
    const float z = position[2];

    // iterates through every bar. This is probably a bad method.

    const int layers = 96;
    const int strips = 40;

    const float supportBox = 19.05;

    for (const auto& pair : map) {
        int ID = pair.first;
        DetectorInfo detectorData = pair.second;
        std::string orientation = detectorData.orientation;

        if ((x >= detectorData.xMin && x < detectorData.xMax)
         && (y >= detectorData.yMin && y < detectorData.yMax)
         && (z >= detectorData.zMin && z <= detectorData.zMax)) {

            std::string stringID = std::to_string(ID);

            int layer = std::stoi(stringID.substr(1,3));
            int strip = std::stoi(stringID.substr(4,5));

            map[ID].energy += edep; // add the energy to the bar
            map[ID].numberOfHits++;
            return;
         }
    }
}

void CosmicsAnalysis::AddDaughterEnergy(std::map<int,DetectorInfo>& map, const std::vector<double>& vertex, const std::vector<double>& endPoint, const double& edep) {

    // Function that adds the energy of a daughter particle to the map, in the bar where the particle was created.
    // The particle has to decay in an absorber.

    const float x = vertex[0];
    const float y = vertex[1];
    const float z = vertex[2];

    std::vector<float> v{x,y,z};

    const float endX = vertex[0];
    const float endY = vertex[1];
    const float endZ = vertex[2];

    const float backHcaldz = 4704;
    const float sideHcaldz = 600;
    const float sideHcalStartZ = 240;

    const float airThick = 2.;
    const float absoThick = 25.;
    const float layerThick = 49.;
    const float supportBox = 19.05;

    const int numLayers = 96;
    const int strips = 40;

    const float absoStartZ = ((backHcaldz + sideHcaldz) / 2) + sideHcalStartZ + sideHcaldz/2 - backHcaldz/2; + airThick + absoThick/2; // the center of the absorber

    for (int i = 1; i <= numLayers; i++) { // iterate through absorbers
        const float absoLayer = absoStartZ + (i-1)*layerThick;
        const float absoZMin = absoLayer - absoThick;
        const float absoZMax = absoLayer + absoThick;
        
        if (endZ > absoZMin && endZ < absoZMax) { // decay in absorber. Score!

            for (const auto& pair : map) {
                int ID = pair.first;
                DetectorInfo detectorData = pair.second;
                std::string orientation = detectorData.orientation;

                if ((x >= detectorData.xMin && x < detectorData.xMax)
                && (y >= detectorData.yMin && y < detectorData.yMax)
                && (z >= detectorData.zMin && z <= detectorData.zMax)) {

                    if (map[ID].energy == 0.0) return; // No thanks!
                    else map[ID].energy += edep; // add the energy to the bar

                    return;
                }
            }
        }
    }

}

bool CosmicsAnalysis::IsInBar(std::vector<float> position) {

    // This function checks if the position is in a bar.

    // get detector id by position
    const float x = position[0];
    const float y = position[1];
    const float z = position[2];

    auto map = detectorMap_;

    // iterate through ids, 010101, 020201, 010301, (orientation, layer, strip)

    int numLayers = 96;

    // Checks if the position is in the absorber
    for (int i = 1; i <= numLayers; i++ ) {
        std::string strID  = ((i % 2 == 0) ? std::string("20") : std::string("10")) + (i < 10 ? "0" : "") + std::to_string(i) + std::string("01");
        int ID = std::stoi(strID);
        DetectorInfo detectorData = map[ID];
        int absoZMin = detectorData.zMin - 2 - 25;
        int absoZMax = detectorData.zMin - 2;

        if (z > absoZMin && z < absoZMax) {
            return false; // the position is in an absorber, i.e. not in a bar
        }
    }
    return true;

    // This checks if the position is specifically in a bar
    // It would take less time to iterate through the absorbers.
    // for (const auto& pair : map) {
    //     DetectorInfo detectorData = pair.second;
    //     if ((x >= detectorData.xMin && x < detectorData.xMax)
    //      && (y >= detectorData.yMin && y < detectorData.yMax)
    //      && (z >= detectorData.zMin && z <= detectorData.zMax)) {
    //         return true;
    //      }
    // }
    // return false;
}

void CosmicsAnalysis::AddEstimatedEnergy(std::map<int,DetectorInfo>& map, const std::vector<float>& position, const double& stepsize) {
    // get detector id by position
    const float x = position[0];
    const float y = position[1];
    const float z = position[2];

    double BetheBloch = 0.233;  // MeV/mm

    for (const auto& pair : map) {
        const int ID = pair.first;
        const DetectorInfo detectorData = pair.second;

        if ((x >= detectorData.xMin && x < detectorData.xMax)
         && (y >= detectorData.yMin && y < detectorData.yMax)
         && (z >= detectorData.zMin && z <= detectorData.zMax)) {
            map[ID].estimatedEnergy += stepsize * BetheBloch; // add the estimated energy to the bar. Depends on the stepsize. For stepsize = 1, each addition is 0.233 MeV for 1 mm.
            return;
         }
    }
}

// <------------------------------------------- FITTING FUNCS ------------------------------------------------------->

void CosmicsAnalysis::lineFunc(double t, const double *p, double &x, double &y, double &z) {
    // a parametric line is define from 6 parameters but 4 are independent
    // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
    // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
    
    x = p[0] + p[1]*t;
    y = p[2] + p[3]*t;
    z = t;
}

void CosmicsAnalysis::curveFunc(double t, const double *p, double &x, double &y, double &z) {
    x = t;
    y = t;

    // Evaluate surface: z = a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
    z = p[0]*x*x + p[1]*y*y + p[2]*x*y + p[3]*x + p[4]*y + p[5];
}

struct SumDistance2 { // for line fitting
    // the TGraph is a data member of the object
    TGraph2D *fGraph;
  
    SumDistance2(TGraph2D *g) : fGraph(g) {}

    double huber(double r2, double delta = 1.0) {
        if (sqrt(r2) <= delta)
            return r2;
        else
            return 2*delta * std::sqrt(r2) - delta * delta;
    }
  
    // calculate distance line-point
    double distance2(double x,double y,double z, const double *p) {
       // distance line point is D= | (xp-x0) cross  ux |
       // where ux is direction of line and x0 is a point in the line (like t = 0)
       ROOT::Math::XYZVector xp(x,y,z);
       ROOT::Math::XYZVector x0(p[0], p[2], 0. );
       ROOT::Math::XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
       ROOT::Math::XYZVector u = (x1-x0).Unit();
       double d2 = ((xp-x0).Cross(u)).Mag2();
       return d2;
    }
  
    // implementation of the function to be minimized
    double operator() (const double *par) {
       assert(fGraph != nullptr);
       double * x = fGraph->GetX();
       double * y = fGraph->GetY();
       double * z = fGraph->GetZ();
       int npoints = fGraph->GetN();
       double sum = 0;
       for (int i  = 0; i < npoints; ++i) {
          double d = huber(distance2(x[i],y[i],z[i],par));
          sum += d;
       }
       return sum;
    }
  
};

struct SumDistance2Curve { // for curve fitting
    // the TGraph is a data member of the object
    TGraph2D *fGraph;
  
    SumDistance2Curve(TGraph2D *g) : fGraph(g) {}

    double huber(double r2, double delta = 1.0) {
        if (sqrt(r2) <= delta)
            return r2;
        else
            return 2*delta * std::sqrt(r2) - delta * delta;
    }

    // Calculate squared residual between point and model surface
    double distance2(double x, double y, double z, const double *p) const {

        double z_fit = p[0]*x*x + p[1]*y*y + p[2]*x*y + p[3]*x + p[4]*y + p[5];
        double res = z - z_fit;
        return res * res;
    }

    // Function to be minimized
    double operator() (const double *par) {
        assert(fGraph != nullptr);
        double *x = fGraph->GetX();
        double *y = fGraph->GetY();
        double *z = fGraph->GetZ();
        int npoints = fGraph->GetN();
        double sum = 0;
        for (int i = 0; i < npoints; ++i) {
            double d = distance2(x[i], y[i], z[i], par);
            sum += d;
        }
        return sum;
    }

};

TVectorD CosmicsAnalysis::fitSurfaceInitialGuess(TGraph2D* graph) {
    int n = graph->GetN();
    double* x = graph->GetX();
    double* y = graph->GetY();
    double* z = graph->GetZ();

    TMatrixD A(n, 6);
    TVectorD Z(n);

    for (int i = 0; i < n; ++i) {
        A(i, 0) = x[i] * x[i];    
        A(i, 1) = y[i] * y[i];    
        A(i, 2) = x[i] * y[i];    
        A(i, 3) = x[i];           
        A(i, 4) = y[i];        
        A(i, 5) = 1.0;             
        Z(i)    = z[i];            
    }

    // Solve A^T A p = A^T z
    TMatrixD At(TMatrixD::kTransposed, A);
    TMatrixD AtA = At * A;
    TVectorD AtZ = At * Z;

    TDecompSVD svd(AtA);
    Bool_t ok;
    TVectorD p = svd.Solve(AtZ, ok);
    if (!ok) {
        std::cerr << "Initial surface fit failed!" << std::endl;
    }

    return p;  
}  


// <-------------------------------------------- OTHER FUNCS -------------------------------------------------------->

void CosmicsAnalysis::getAngle(const std::map<int, ldmx::SimParticle, std::less<int>, std::allocator<std::pair<const int, ldmx::SimParticle> > >& sim_particles, const std::vector<ldmx::SimCalorimeterHit, std::allocator<ldmx::SimCalorimeterHit>>& sim_hits) {

    // This function finds the polar and azimuthal angle given some vector

    auto particle = sim_particles.at(1);
    auto start_point = particle.getVertex();
    std::vector<float> end_point;

    for (int i = 1; i < sim_hits.size(); i++){ // i has to be at least 1
        std::vector<float> pos = (sim_hits[sim_hits.size() - i]).getPosition();
        double x = pos[0], y = pos[1], z = pos[2];

        // Find the last point in the detector. This does not account for the side-hcal and ecal
        if (std::abs(x) < hcalXMax_ && y > hcalYMin_ && y < hcalYMax_ && z > hcalZMin_ && z < hcalZMax_) { // Inside the detector
            end_point = pos;
            break;
        }
        // Last iteration. No detector was found. This usually means that all hits were in the side-hcal or ecal
        if (i == sim_hits.size() - 1) {
            muon_deviation_angle_ = 365;
            return;
        }

    }

    // Vectors for reference
    auto start_mom = particle.getMomentum();
    auto end_mom = particle.getEndPointMomentum();
    muon_init_mom_ = std::hypot(start_mom[0], start_mom[1], start_mom[2]);
    muon_final_mom_ = std::hypot(end_mom[0], end_mom[1], end_mom[2]);
    muon_mom_ratio_ = muon_final_mom_/muon_init_mom_;

    // Calculates the 
    std::vector<double> vec_r = {(end_point[0] - start_point[0]), (end_point[1] - start_point[1]), (end_point[2] - start_point[2])};
    double mag_vec_r = std::hypot(vec_r[0], vec_r[1], vec_r[2]); 
    muon_deviation_angle_ = acos((vec_r[0]*start_mom[0] + vec_r[1]*start_mom[1] + vec_r[2]*start_mom[2])/(mag_vec_r * muon_init_mom_)) * 180/M_PI;

    // Azimuthal vector references
    std::vector<double> start_mom_azimuthal = {start_mom[0], 0, start_mom[2]};
    double mag_start_mom_azimuthal = std::hypot(start_mom_azimuthal[0], start_mom_azimuthal[1], start_mom_azimuthal[2]);

    // Gets the polar and azimuthal angles
    std::vector<double> vec_y = {0,-1,0};
    std::vector<double> vec_z = {0,0,1};

    muon_polar_angle_ = acos((vec_y[0]*start_mom[0] + vec_y[1]*start_mom[1] + vec_y[2]*start_mom[2])/(muon_init_mom_)) * 180/M_PI;  
    muon_azimuthal_angle_ = acos((vec_z[0]*start_mom_azimuthal[0] + vec_z[1]*start_mom_azimuthal[1] + vec_z[2]*start_mom_azimuthal[2])/(mag_start_mom_azimuthal)) * 180/M_PI;

    // For 360 degree angle
    if (start_mom_azimuthal[0] < 0){
        muon_azimuthal_angle_ = -muon_azimuthal_angle_;
    }
}

void CosmicsAnalysis::getAngleEstimate(double& polarAngleEstimate, double& azimuthalAngleEstimate, const std::vector<float>& line) {

    // This function finds the polar and azimuthal angle given a line-fit (i.e. estimated angle)

    // Vectors for reference
    std::vector<double> vec_y;
    std::vector<double> vec_z = {0,0,1};

    if (line[1] > 0) {
        vec_y = {0,1,0};
    }
    else {
        vec_y = {0,-1,0};
    }

    double lineMagnitude = std::hypot(line[0], line[1], line[2]);

    // polar angle
    polarAngleEstimate = acos((vec_y[0]*line[0] + vec_y[1]*line[1] + vec_y[2]*line[2])/(lineMagnitude)) * 180/M_PI;  

    // azimuthal angle
    std::vector<double> azimuthalLine{line[0], 0, line[2]};
    double azimuthalLineMagnitude = std::hypot(azimuthalLine[0], azimuthalLine[1], azimuthalLine[2]);
    azimuthalAngleEstimate = acos((azimuthalLine[0]*vec_z[0]) + (azimuthalLine[1]*vec_z[1]) + (azimuthalLine[2]*vec_z[2])/azimuthalLineMagnitude) * 180/M_PI;
}


// <------------------------------------------- PLOTTING FUNCS ------------------------------------------------------>

void CosmicsAnalysis::getPlot1D(std::pair<std::string, std::string> pair) {

    // this needs to be redefined (along with other things in this func) if we want to analyze another TTree
    TLeaf* xLeaf = sim_muon_tree_->GetLeaf((pair.first).c_str());
    TLeaf* yLeaf = sim_muon_tree_->GetLeaf((pair.second).c_str());
    Long64_t nEntries = sim_muon_tree_->GetEntries();

    if (!xLeaf || !yLeaf) {
        std::cerr << "Error: One or both leaves not found!" << std::endl;
        return;
    }

    // Get max value in the yLeaf
    double max_val_y = 0;
    double min_val_y = 1000000;
    double max_val_x = 0;
    double min_val_x = 1000000;

    // finds the max y and x vals
    for (Long64_t i = 0; i < nEntries; i++) {
        sim_muon_tree_->GetEntry(i);
        double y = yLeaf->GetValue(); 
        double x = xLeaf->GetValue();
        
        if (y > max_val_y) {
            max_val_y = y;
        }
        if (y < min_val_y) {
            min_val_y = y;
        }
        if (x > max_val_x) {
            max_val_x = x;
        }
        if (x < max_val_x) {
            min_val_x = x;
        }
    }

    // Define canvas, legend and histogram
    TCanvas* T1 = new TCanvas("t1");
    TLegend* leg = new TLegend(0.6,0.7,0.85,0.85);
    TH1F* H1 = new TH1F("Data", "", 90, 0, 2);

    // Resolution of the canvas
    T1->SetCanvasSize(1200, 800);

    // Adds upper and right ticks and sets a grid
    T1->SetTickx();
    T1->SetTicky();
    T1->SetGridx();
    // T1->SetGridy();

    // Adds margins
    T1->SetLeftMargin(0.15);
    T1->SetRightMargin(0.15);
    T1->SetBottomMargin(0.15);
    T1->SetTopMargin(0.1);

    // Sets title and label sizes
    H1->GetXaxis()->SetTitleSize(0.05);
    H1->GetXaxis()->SetLabelSize(0.04);
    H1->GetYaxis()->SetTitleSize(0.05);
    H1->GetYaxis()->SetLabelSize(0.04);

    // H1->GetYaxis()->SetRangeUser(0, 2000); // Sets the range for visualization
    
    // Set the histograms stat box location
    gStyle->SetStatX(0.85);
    gStyle->SetStatY(0.90);
    gStyle->SetOptStat(11);

    // Set the color of the colz
    // gStyle->SetStatBorderSize(0); // Removes border of the statistics from the hist

    H1->GetXaxis()->SetTitle("E_{muon} / E_{est}");
    H1->GetYaxis()->SetTitle("Nr of interactions");
    H1->GetYaxis()->SetNoExponent(kFALSE);
    H1->GetYaxis()->SetNdivisions(505, kTRUE);
    H1->GetXaxis()->SetNdivisions(505, kTRUE);

    // Set log scale on z
    // gPad->SetLogy(1);

    // Set log scale on x,y
    // gPad->SetLogx(1);
    // gPad->SetLogy(1);

    // If y is a vector
    /////////////////////////////COMMENT/UNCOMMENT/////////////////////////////

    double x;
    std::vector<float>* y = nullptr;

    sim_muon_tree_->SetBranchAddress((pair.first).c_str(), &x);
    sim_muon_tree_->SetBranchAddress((pair.second).c_str(), &y);
    // hcal_hit_tree_->SetBranchAddress((pair.second).c_str(), &y);


    for (Long64_t i = 0; i < nEntries; i++) {
        sim_muon_tree_->GetEntry(i);
        // hcal_hit_tree_->GetEntry(i);
    
        // Loop over elements in the vector
        for (float y_val : *y) {
            // if (y_val > 1000 || y_val < -1000) {continue;}
            if (y_val == 0) continue;
            H1->Fill(y_val, 1);
        }
    }
    /////////////////////////////COMMENT/UNCOMMENT/////////////////////////////
    
    // If y is not a vector
    /////////////////////////////COMMENT/UNCOMMENT/////////////////////////////

    // for (Long64_t i = 0; i < nEntries; i++) {
    //     sim_muon_tree_->GetEntry(i);
    //     double x = xLeaf->GetValue();
    //     double y = yLeaf->GetValue();
        
    //     // if (y <= 0.000000001) continue; // Not a hit. 0.000000001 is 0.001 eV
    //     if (y <= 0.0001) {
    //         continue;
    //     }
    //     H1->Fill(x,1);
    // }

    /////////////////////////////COMMENT/UNCOMMENT/////////////////////////////

    // Code below this line was used for fitting a landau-gaus convolution fit
    // -------------------------------------------------------------------------------------------------------
    // float minRange = 0.1;
    // float maxRange = 0.3;

    // TF1Convolution *conv = new TF1Convolution("landau", "gaus", minRange, maxRange, false);
    // conv->SetNofPointsFFT(10000);
    // conv->SetRange(minRange, maxRange);

    // TF1 *fitFunc = new TF1("fitFunc", *conv, minRange, maxRange, conv->GetNpar());
    // fitFunc->SetParNames("Norm", "MPV", "Gauss Mean", "Landau Sigma", "Gauss Sigma");

    // fitFunc->SetParameters(
    //     H1->GetMaximum(), 
    //     1.62791e-01, //H1->GetMean(),  
    //     1.67177e-01,
    //     3.82356e-03,  
    //     1.39512e-02
    // );
    
    // // fitFunc->SetParLimits(0, 10, 1e6); 
    // // fitFunc->SetParLimits(1, 0.15, 0.21);
    // // fitFunc->SetParLimits(2, 10, 1e6); 
    // // fitFunc->SetParLimits(3, 0.15, 0.21); 
    // // fitFunc->SetParLimits(4, 0.001, 0.3); 

    // H1->GetXaxis()->SetRangeUser(minRange, maxRange);
    // TF1 *landauFunc = new TF1("landauFunc", "landau", 0.0, 0.5);
    // TF1 *landauShifted = new TF1("landauShifted", "[0]*TMath::Landau(x, [1], [2]) + [3]", 0.0, 1);
    // landauShifted->SetParameters(H1->GetMaximum(), 0.233, 0.05, 100);
    // landauFunc->SetParameters(H1->GetMaximum(), 0.233, 0.05);  // Initial guesses
    // H1->Fit(landauShifted, "R");  // "R" restricts fit to function range
    // -------------------------------------------------------------------------------------------------------

    H1->Fit("gaus");
    H1->Draw();
    T1->SaveAs("output.pdf");
    H1->Write();
    delete H1;
    delete T1;
}

void CosmicsAnalysis::getPlot2D(std::pair<std::string, std::string> pair) {

    // this needs to be redefined (along with other things in this func) if we want to analyze another TTree
    TLeaf* xLeaf = sim_muon_tree_->GetLeaf((pair.first).c_str());
    TLeaf* yLeaf = sim_muon_tree_->GetLeaf((pair.second).c_str());
    Long64_t nEntries = sim_muon_tree_->GetEntries();

    if (!xLeaf || !yLeaf) {
        std::cerr << "Error: One or both leaves not found!" << std::endl;
        return;
    }

    // finds the max y and x vals
    double max_val_y = 0;
    double max_val_x = 0;
    for (Long64_t i = 0; i < nEntries; i++) {
        sim_muon_tree_->GetEntry(i);
        double y = yLeaf->GetValue(); 
        double x = xLeaf->GetValue();
        
        if (y > max_val_y) {
            max_val_y = y;
        }
        if (x > max_val_x) {
            max_val_x = x;
        }
    }

    // Define canvas, legend and histogram
    TCanvas* T1 = new TCanvas("t1");
    TLegend* leg = new TLegend(0.6,0.7,0.85,0.85);
    TH2F* H2 = new TH2F("Data", "", 100, 0, 75, 80, 0, 20); // use this to set the limits of the plot (xbin, xmin, xmax, ybin, ymin, ymax)

    // Resolution of the canvas
    T1->SetCanvasSize(1200, 800);

    // Adds upper and right ticks and sets a grid
    T1->SetTickx();
    T1->SetTicky();
    // T1->SetGridx();
    // T1->SetGridy();

    // Adds margins
    T1->SetLeftMargin(0.15);
    T1->SetRightMargin(0.15);
    T1->SetBottomMargin(0.15);
    T1->SetTopMargin(0.1);

    // Sets title and label sizes
    H2->GetXaxis()->SetTitleSize(0.05);
    H2->GetXaxis()->SetLabelSize(0.04);
    H2->GetYaxis()->SetTitleSize(0.05);
    H2->GetYaxis()->SetLabelSize(0.04);
    

    // Set the histogram stat box location
    gStyle->SetStatX(0.85);
    gStyle->SetStatY(0.90);
    gStyle->SetOptStat(11);

    // gStyle->SetStatBorderSize(0); // Removes border of the statistics from the hist

    // Set axis labels and number of divisions on y-axis
    H2->GetXaxis()->SetTitle("Polar angle #theta [deg]");
    H2->GetYaxis()->SetTitle("Estimated energies in bars [MeV]");
    H2->GetYaxis()->SetNdivisions(505, kTRUE);
    H2->GetXaxis()->SetNdivisions(505, kTRUE);

    // Set the axes to show values in exponential form
    // Only applies if val > 100000 I think
    H2->GetXaxis()->SetNoExponent(kFALSE);
    H2->GetYaxis()->SetNoExponent(kFALSE);

    // Set log scale on z
    gPad->SetLogz(1);

    // Choose one of the following methods below for analyzing y as a vector or as a number
    // If y is a vector
    /////////////////////////////COMMENT/UNCOMMENT/////////////////////////////

    double x;
    std::vector<double>* y = nullptr;

    sim_muon_tree_->SetBranchAddress((pair.first).c_str(), &x);
    sim_muon_tree_->SetBranchAddress((pair.second).c_str(), &y);

    for (Long64_t i = 0; i < nEntries; i++) {
        sim_muon_tree_->GetEntry(i);
    
        // Loop over elements in the vector
        for (float y_val : *y) {
            // std::cout << std::setprecision(30) << y_val << std::endl;
            if (y_val <= 0.000000001) continue; // Not a hit. 0.000000001 is 0.001 eV
            // if (y_val > 1000 || y_val < -1000) {continue;}
            // if (y_val < 2) continue;

            H2->Fill(x, y_val);
        }
    }
    /////////////////////////////COMMENT/UNCOMMENT/////////////////////////////
    
    // If y is not a vector
    /////////////////////////////COMMENT/UNCOMMENT/////////////////////////////

    // for (Long64_t i = 0; i < nEntries; i++) {
    //     sim_muon_tree_->GetEntry(i);
    //     double x = xLeaf->GetValue(); 
    //     double y = yLeaf->GetValue(); 

    //     if (y == NAN || y == -NAN) continue;
        
    //     if (y <= 0.000000001) {
    //         continue;
    //     }

    //     H2->Fill(x, y);        
    // }

    /////////////////////////////COMMENT/UNCOMMENT/////////////////////////////    

    // Draw the histogram
    H2->Draw("COLZ");
    T1->SaveAs("output.pdf");
    H2->Write();
    delete H2;
    delete T1;
}

// <---------------------------------------------- MAIN ANALYSIS ----------------------------------------------------->

void CosmicsAnalysis::onProcessStart() { // Runs once at start
    getHistoDirectory(); // forget this -> silently no output histogram

    file_ = new TFile("Analysis/cosmics_data.root", "RECREATE"); // File for storing TTrees

    sim_daughter_tree_ = new TTree("Sim_Daughter_Data","data from simulated daughter particles"); // For storing data about daughter particles from SimParticles
    sim_daughter_tree_->Branch("daughter_pdgID_", &daughter_pdgID_, "daughter_pdgID_/I");
    sim_daughter_tree_->Branch("daughter_energy_", &daughter_energy_, "daughter_energy_/D");
    sim_daughter_tree_->Branch("daughter_mass_", &daughter_mass_, "daughter_mass_/D");
    sim_daughter_tree_->Branch("daughter_kinetic_", &daughter_kinetic_, "daughter_kinetic_/D");

    sim_muon_tree_ = new TTree("Sim_Muon_Data", "data from simulated muons"); // For storing data from SimParticles
    sim_muon_tree_->Branch("muon_energy_", &muon_energy_, "muon_energy_/D");
    sim_muon_tree_->Branch("muon_deviation_angle_", &muon_deviation_angle_, "muon_deviation_angle_/D");
    sim_muon_tree_->Branch("muon_polar_angle_", &muon_polar_angle_, "muon_polar_angle_/D");
    sim_muon_tree_->Branch("muon_kinetic_", &muon_kinetic_, "muon_kinetic_/D");
    sim_muon_tree_->Branch("muon_init_mom_", &muon_init_mom_, "muon_init_mom_/D");
    sim_muon_tree_->Branch("muon_final_mom_", &muon_final_mom_, "muon_final_mom_/D");
    sim_muon_tree_->Branch("muon_mom_diff_", &muon_mom_diff_, "muon_mom_diff_/D");
    sim_muon_tree_->Branch("muon_mass_", &muon_mass_, "muon_mass_/D");
    sim_muon_tree_->Branch("muon_mom_ratio_", &muon_mom_ratio_, "muon_mom_ratio_/D");
    sim_muon_tree_->Branch("muon_edeps_", &edeps_vec_ptr_);
    sim_muon_tree_->Branch("muon_est_edeps_", &est_edeps_vec_ptr_);
    sim_muon_tree_->Branch("muon_size_", &muon_size_, "muon_size_/I");
    sim_muon_tree_->Branch("muon_edeps_total_event_", &muon_edeps_total_event_, "muon_edeps_total_event_/D");
    sim_muon_tree_->Branch("muon_bars_hit_", &muon_bars_hit_, "muon_bars_hit_/I");
    sim_muon_tree_->Branch("muon_edeps_per_bar_", &muon_edeps_per_bar_, "muon_edeps_per_bar_/D");
    sim_muon_tree_->Branch("muon_edeps_per_barmm_", &muon_edeps_per_barmm_, "muon_edeps_per_barmm_/F");
    sim_muon_tree_->Branch("muon_azimuthal_angle_", &muon_azimuthal_angle_, "muon_azimuthal_angle_/D");
    sim_muon_tree_->Branch("muon_est_polar_angle_", &muon_est_polar_angle_, "muon_est_polar_angle_/D");
    sim_muon_tree_->Branch("edeps_est_edeps_ratio_", &est_edeps_ratio_ptr_);

    hcal_hit_tree_ = new TTree("Hcal_Hit_Data", "data from simulated hits"); // For storing data from HcalSimHits
    hcal_hit_tree_->Branch("recPE_", &recPE_);

    graphFile_ = new TFile("Analysis/sim_graphs.root", "RECREATE"); // For saving plots

    // Get start-time
    start_time = std::chrono::high_resolution_clock::now();

    // Builds a map that contains information from every event (unless used otherwise)
    buildMap(detectorMap_);

}

void CosmicsAnalysis::analyze(const framework::Event& event) {

    // Get the data collections from the event bus
    const auto& hcal_sim_hits{event.getCollection<ldmx::SimCalorimeterHit>("HcalSimHits")};
    const auto& sim_particles{event.getMap<int, ldmx::SimParticle>("SimParticles", "cosmics")};
    // const auto& hcal_rec_hits{event.getCollection<ldmx::HcalHit>("HcalRecHits")}; // This is reconstructed information that was not used for the project, but it was analyzed in the beginning.

    // Get the event number
    int eventNumber = event.getEventNumber();

    if (eventNumber % 10000 == 0) { // print every 10,000 events
        std::cout << "event: " << eventNumber << " started" << std::endl;
    }

    // Some events do not interact with the detector. We want to find these events and store that event number
    // This shouldn't really happend with the cosmic generator. I assume that a particle is created close to the edge
    // and passes through a tiny part of the detector without a hit, or it only passes through absorber material.
    // Now for SOME reason, the code doesn't work with 1 hit either. It stores 1 hit but when I iterate through it, 
    // it's empty.
    if (hcal_sim_hits.size() <= 1) {
        nr_events_miss_++;
        return;
    }

    // line for fitting later
    std::vector<float> line; // linear fit
    std::vector<float> unitLine; //normalized line
    double R2; // coeff of determination

    std::map<int,DetectorInfo> eventMap; // Map of the current event
    buildMap(eventMap);

    // vectors for storing the energy deposits for the event
    std::vector<double> energies;
    std::vector<double> estimatedEnergies;

    // variables for storing angle estimates
    double polarAngleEstimate;
    double azimuthalAngleEstimate;

    // Finds the polar and azimuthal angles
    getAngle(sim_particles, hcal_sim_hits);

    if (graphcheck && muon_deviation_angle_ < angle_threshold_) { // selection on angle

        // Declare graph and canvas
        canvas_ = new TCanvas();
        graph_ = new TGraph2D();

        graph_->SetMarkerStyle(24);
        graph_->SetTitle("Cosmics muon path; x-axis ; y-axis ; z-axis");

        // determine the origin
        muon_origin_ = hcal_sim_hits[0].getPosition();

        // Iterate over the simulated hcal hits
        Int_t index = 0;
        int lastBar = 0;
        std::vector<std::vector<float>> currentHitPosition;
        std::vector<double> currentHitEnergy;

        // means for R²
        float xMean = 0;
        float yMean = 0;
        float zMean = 0;

        int its = 0;

        for (const auto& hit : hcal_sim_hits) {

            std::vector<float> pos = hit.getPosition();
            float xHit = pos[0];
            float yHit = pos[1];
            float zHit = pos[2];

            auto hitEdeps = (hit.getEdeps())[0];

            double energyHit = 0;

            // We don't care about products in the ecal or side-hcal for now
            if (std::abs(pos[0]) < 1500 && std::abs(pos[1]) < 1519.05 && pos[2] < 840) {
                delete graph_;
                nr_events_ecal_++;
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

                xMean += xHit;
                yMean += yHit;
                zMean += zHit;

                its++;

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

                xMean += xHit;
                yMean += yHit;
                zMean += zHit;

                its++;
            }

            graph_->SetPoint(index++, xHit, yHit, zHit);

        }

        // divide total values by number of hits
        xMean /= its;
        yMean /= its;
        zMean /= its;

        if (graph_->GetN() != 0) {
            
            // We have the graph, now we fit 
            ROOT::Fit::Fitter fitter;

            // make the functor objet
            SumDistance2 sdist(graph_);
            ROOT::Math::Functor fcn(sdist,4);

            auto muon = sim_particles.at(1);
            auto momentum = muon.getMomentum();
            
            // set the function and the initial parameter values
            double pStart[4] = {muon_origin_[0],(momentum[1]/momentum[2]),muon_origin_[1],(momentum[1]/momentum[2])}; //x_0, dx/dz, y_0, dy/dz
            fitter.SetFCN(fcn,pStart);

            // set step sizes different than default ones (0.3 times parameter values)
            for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);
            
            bool ok = fitter.FitFCN();
            if (!ok) { // I'm unsure why the fits fail. I assume bad initial param vals
                // Error("line3Dfit","Line3D Fit failed");
                return;
            }
            
            const ROOT::Fit::FitResult & result = fitter.Result();
            auto minFCN = result.MinFcnValue();

            if (minFCN > 400) return; // selection

            // result.Print(std::cout); // Here we can print the results

            // get fit parameters
            const double * parFit = result.GetParams();
            
            // draw the fitted line
            int n = 100000;
            double t0 = 0;
            double dt = 100000;
            std::vector<float> initPos;
            TPolyLine3D *l = new TPolyLine3D(n);
            for (int i = 0; i < n; ++i) {
                double t = t0+ dt*i/n;
                double x,y,z;
                lineFunc(t,parFit,x,y,z);
                float xf = float(x), yf = float(y), zf = float(z);
                if (i == 0) {
                    initPos = {xf,yf,zf};
                }
                if (i == 1) {
                    line = {(xf - initPos[0]),(yf - initPos[1]),(zf - initPos[2])};
                }
                l->SetPoint(i,x,y,z);
            }
            l->SetLineColor(kRed);

            float mag = std::hypot(line[0],line[1],line[2]);
            unitLine = {line[0]/mag, line[1]/mag, line[2]/mag};

            // Calculate R2
            double SSE = 0;
            double SST = 0;

            for (const auto& hit : hcal_sim_hits) {
                const std::vector<float>& pos = hit.getPosition();
                const float& x = pos[0];
                const float& y = pos[1];
                const float& z = pos[2];
                const float& x0 = initPos[0];
                const float& y0 = initPos[1];
                const float& z0 = initPos[2];

                const std::vector<float> pointLine{(x-x0), (y-y0), (z-z0)};

                auto t = ((pointLine[0]*unitLine[0]) + (pointLine[1]*unitLine[1]) + (pointLine[2]*unitLine[2]));

                std::vector<float> linePoint = {
                    initPos[0] + line[0]/mag * t,
                    initPos[1] + line[1]/mag * t,
                    initPos[2] + line[2]/mag * t
                };
    
                float xline = linePoint[0];
                float yline = linePoint[1];
                float zline = linePoint[2];

                SSE += std::pow((x-xline),2) + std::pow((y-yline),2) + std::pow((z-zline),2);
                SST += std::pow((x-xMean),2) + std::pow((y-yMean),2) + std::pow(z-zMean,2);
            }

            double R2 = 1 - (SSE/SST);

            if (R2 < 0.99) return;

            straight_tracks_ += 1; // event passed selections
            

            // Now do a curve fit and check if it is better than the line fit. If it is we throw the event
            // NOTE! This was never implemented fully. It finds a fit, but it is not good at all

            // ROOT::Fit::Fitter curveFitter;

            // SumDistance2Curve sCurve(graph_);
            // ROOT::Math::Functor fcnCurve(sCurve,6);

            // TVectorD pGuess = fitSurfaceInitialGuess(graph_);
            // double pStartCurve[6];
            // for (int i = 0; i < 6; ++i) pStartCurve[i] = pGuess[i];
            // // for (int i = 0; i < 6; ++i) std::cout << std::setprecision(5) << pStartCurve[i] << std::endl;
            // curveFitter.SetFCN(fcnCurve,pStartCurve);
            // for (int i = 0; i < 6; ++i) curveFitter.Config().ParSettings(i).SetStepSize(0.01);

            // bool okCurve = curveFitter.FitFCN();
            // if (!okCurve) {
            //     Error("line3Dfit","Line3D Fit failed");
            //     return;
            // }
            
            // const ROOT::Fit::FitResult & resultCurve = curveFitter.Result();
            // auto minFCNCurve = resultCurve.MinFcnValue();

            // // resultCurve.Print(std::cout);

            // const double * parFitC = resultCurve.GetParams();

            // TPolyLine3D *lc = new TPolyLine3D(n);
            // for (int i = 0; i < n; ++i) {
            //     double t = t0+ dt*i/n;
            //     double x,y,z;
            //     curveFunc(t,parFitC,x,y,z);
            //     float xf = float(x), yf = float(y), zf = float(z);
            //     lc->SetPoint(i,x,y,z);
            // }
            // lc->SetLineColor(kBlue);

            
            if (estimatedcheck) {
                // add estimated energies
                // This is VERY slow since the initial point from the fit can be far away from the detector.
                // Also the add method is VERY slow since it checks every bar

                float norm = std::hypot(line[0], line[1], line[2]); // norm of vec
                std::vector<float> lineNorm = {line[0]/norm, line[1]/norm, line[2]/norm};

                bool inside = false;

                const int stepsize = 1;
                int step = 0;

                while (true) {

                    float x = initPos[0] + lineNorm[0]*step;
                    float y = initPos[1] + lineNorm[1]*step;
                    float z = initPos[2] + lineNorm[2]*step;
                    step += stepsize;

                    std::vector<float> position{x,y,z};

                    if (abs(x) < hcalXMax_ && (y > hcalYMin_ && y < hcalYMax_) && (z > hcalZMin_ && z < hcalZMax_)) {
                        AddEstimatedEnergy(eventMap, position, stepsize);
                        inside = true;
                        continue;
                    }

                    if (inside) break; // breaks when we've left the detector
                }
            } 

            // Sets the corners of our detector. Just for visual representation
            // graph_->SetPoint(index, hcalXMin_, -hcalYMin_, 0);
            // graph_->SetPoint(index+1, hcalXMax_, hcalYMax_, hcalZMax_);

            graph_->GetXaxis()->SetTitleOffset(1.5);
            graph_->GetYaxis()->SetTitleOffset(2.0);
            graph_->GetZaxis()->SetTitleOffset(1.5);
            graph_->SetTitle("");

            graph_->Draw("P0"); // draw hit points
            l->Draw("same"); // draw linear fit
            // lc->Draw("same"); // draw curve fit
            canvas_->Update();

            // Save the canvas
            graphFile_->cd();
            canvas_->Write(Form("canvas_%d", eventNumber));

            if (eventNumber == 100000000) { // Here we can check an individual event and make 2D projections of it

                std::cout << "R²: " << R2 << std::endl;
                std::cout << "minFCN: " << minFCN << std::endl;

                canvas_->SaveAs("3D plot.pdf");

                int n = graph_->GetN();
                double* x = graph_->GetX();
                double* y = graph_->GetY();

                TCanvas* c2 = new TCanvas("c2", "", 800, 600);
                c2->SetTitle("");
                
                c2->SetCanvasSize(1200, 800);

                // Adds upper and right ticks and sets a grid
                c2->SetTickx();
                c2->SetTicky();
                c2->SetGridx();

                // Adds margins
                c2->SetLeftMargin(0.15);
                c2->SetRightMargin(0.15);
                c2->SetBottomMargin(0.15);
                c2->SetTopMargin(0.1);

                TGraph* g_yx = new TGraph(n, x, y);
                g_yx->SetTitle("");
                g_yx->SetMarkerStyle(24);
                g_yx->SetMarkerSize(2);

                g_yx->GetXaxis()->SetTitle("X [mm]");
                g_yx->GetYaxis()->SetTitle("Y [mm]");

                int nLinePts = l->GetN();
                TPolyLine* projLine = new TPolyLine(nLinePts);
                for (int i = 0; i < n; ++i) {
                    double t = t0+ dt*i/n;
                    double x,y,z;
                    lineFunc(t,parFit,x,y,z);
                    projLine->SetPoint(i,x,y);
                }
                projLine->SetLineColor(kRed);

                g_yx->Draw("AP");
                projLine->Draw("L same");

                TLegend* legend = new TLegend(0.38, 0.75, 0.62, 0.85); //x1, y1, x2, y2
                legend->AddEntry(g_yx, "Data", "p");
                legend->AddEntry(projLine, "Least squares fit", "l");
                legend->SetTextSize(0.035);
                legend->Draw();
                                
                c2->SaveAs("fit xy.pdf");

                delete c2;
                delete g_yx;
                delete legend;

                double* z = graph_->GetZ();

                TCanvas* c3 = new TCanvas("c3", "", 800, 600);
                
                c3->SetCanvasSize(1200, 800);

                // Adds upper and right ticks and sets a grid
                c3->SetTickx();
                c3->SetTicky();
                c3->SetGridx();

                // Adds margins
                c3->SetLeftMargin(0.15);
                c3->SetRightMargin(0.15);
                c3->SetBottomMargin(0.15);
                c3->SetTopMargin(0.1);

                TGraph* g_yz = new TGraph(n, z, y);
                g_yz->SetTitle("");
                g_yz->SetMarkerStyle(24);
                g_yz->SetMarkerSize(2);

                g_yz->GetXaxis()->SetTitle("Z [mm]");
                g_yz->GetYaxis()->SetTitle("Y [mm]");

                TPolyLine* projLine2 = new TPolyLine(nLinePts);
                for (int i = 0; i < n; ++i) {
                    double t = t0+ dt*i/n;
                    double x,y,z;
                    lineFunc(t,parFit,x,y,z);
                    projLine2->SetPoint(i,z,y);
                }
                projLine2->SetLineColor(kRed);

                g_yz->Draw("AP");
                projLine2->Draw("L same");

                TLegend* legend2 = new TLegend(0.38, 0.75, 0.62, 0.85); //x1, y1, x2, y2
                legend2->AddEntry(g_yx, "Data", "p");
                legend2->AddEntry(projLine, "Least squares fit", "l");
                legend2->SetTextSize(0.035);
                legend2->Draw();
                                
                c3->SaveAs("fit zy.pdf");

                delete c3;
                delete g_yz;
                delete legend2;
            }
        }

        // Delete so we don't get memory issues. For some reason, sometimes deleting doesn't clear the address?
        // Basically, to avoid memory issues, try to not create global ROOT members that store information from all events
        delete graph_;
        delete canvas_;
    }

    for (const auto& [trackID, particle]: sim_particles){

        if (trackID == 1) { // trackID == 1 is the muon

            muon_size_ = 1; // This is pretty pointless

            // Check that muon doesn't decay in the detector!
            // this can be commented out if we want to include that data
            auto end_point = particle.getEndPoint();
            if ((abs(end_point[0] < hcalXMax_)) && (end_point[1] < hcalYMax_ && end_point[1] > hcalYMin_) && (end_point[2] < hcalZMax_ && end_point[2] > hcalZMin_)) {
                decayed_muons_ += 1;
                return;
            }

            if (line.size() != 0) { // line.size() == 0 means that the line fit didn't work
                // get angle estimates
                getAngleEstimate(polarAngleEstimate, azimuthalAngleEstimate, line);

                // This saves the event with the largest deviation, for analysis.
                if (largest_polar_angle_estimate_deviation_ < abs(muon_polar_angle_ - polarAngleEstimate)) {
                    specific_event_ = eventNumber;
                    largest_polar_angle_estimate_deviation_ = abs(muon_polar_angle_ - polarAngleEstimate);
                }
            }
            
            // Assing these variables for filling TTrees.
            muon_energy_ = particle.getEnergy();
            muon_kinetic_ = particle.getKineticEnergy();
            muon_mass_ = particle.getMass();

            muon_edeps_total_event_ = 0;
            muon_bars_hit_ = 0;

            auto initPos = hcal_sim_hits[0].getPosition();
            auto finalPos = (hcal_sim_hits.back()).getPosition();
            std::vector<float> diff{(finalPos[0] - initPos[0]), (finalPos[1] - initPos[1]), (finalPos[2] - initPos[2])};
            float detectorLength = std::hypot(diff[0], diff[1], diff[2]);
            muon_mom_diff_ = (muon_init_mom_ - muon_final_mom_)/detectorLength;


            // Iterate over the hcal hits
            // One particle can leave energy in one bar several times. We want to store energy hits per bar
            int count = 0;
            for (const auto& hit : hcal_sim_hits) {
                count++;

                auto pdgIDs = hit.getPdgIds();
                // We only want to store hits from the muon!
                // This removes all the low energy hits (probably from electrons) so we can focus on the MIP peak
                if (pdgIDs[0] != -13 && pdgIDs[0] != 13) continue;

                auto position = hit.getPosition();
                auto length = hit.getPathLength();
                auto hitEdeps = (hit.getEdeps())[0]; // Get the hit

                auto detectorID = hit.getID();
                auto id = ldmx::HcalID(detectorID);

                // 1 for vertical, 0 for horizontal
                // even layers vertical, odd horizontal.
                bar_orientation_ = id.layer() % 2 == 0 ? 1 : 0;

                if (hitEdeps <= 0.000000001) continue; // I don't think we actually have values that are zero
                //                                      but this is here just in case.
                //                                      0.000000001 is 0.001 eV

                // add the energy to the corresponding bar
                if ((chose_orientation_ == 2 || bar_orientation_ == chose_orientation_)) {
                    addToMap(eventMap, position, hitEdeps); // map for current event
                    // addToMap(detectorMap_, position, hitEdeps); // map for all events
                }
                
                // Adds the total energy of the event
                muon_edeps_total_event_ += hitEdeps;
            }

            muon_edeps_per_bar_ = muon_edeps_total_event_/muon_bars_hit_;
        }
        
        daughter_pdgID_ = particle.getPdgID();
        daughter_energy_ = particle.getEnergy();
        daughter_mass_ = particle.getMass();

        sim_daughter_tree_->Fill();
    }

    std::vector<double> ratio;

    for (const auto& pair : eventMap) {
        int ID = pair.first;
        DetectorInfo data = pair.second;
        double energy = data.energy;
        double estEnergy = data.estimatedEnergy;
        if (energy == 0.0 && estEnergy == 0.0) continue; // unnessecary to add energies
    
        estimatedEnergies.push_back(estEnergy);
        energies.push_back(energy);

        ratio.push_back(energy/estEnergy);
    }
    
    muon_edeps_ = energies;
    muon_est_edeps_ = estimatedEnergies;
    edeps_est_edeps_ratio_ = ratio;
    sim_muon_tree_->Fill();

    total_events_ += 1; // This checks how many events passed
    
    // Print event and time
    auto current_time = std::chrono::high_resolution_clock::now();
    auto time_diff = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time) / 1000.0;
    // std::cout << "Event " << eventNumber << " done! Time elapsed: " << std::fixed << std::setprecision(1) << time_diff.count() << "s" << std::endl;
}

void CosmicsAnalysis::onProcessEnd(){

    std::cout << "--------------------------Finished---------------------------------" << std::endl;
    std::cout << "Largest estimated polar angle deviation: " << largest_polar_angle_estimate_deviation_ << " for event: " << specific_event_ << std::endl;
    std::cout << "total events analyzed: " << total_events_ << std::endl;
    std::cout << "events that miss: " << nr_events_miss_ << std::endl;

    float decay_ratio = static_cast<float>(decayed_muons_)/total_events_;
    std::cout << "Ratio of muons decaying: " << decay_ratio << std::endl;
    float straight_ratio = static_cast<float>(straight_tracks_)/total_events_;
    std::cout << "Ratio of straight tracks, given selections (R², d², cone): " << straight_ratio << std::endl;

    // Check how many hits we have in each scintillator bar
    int maxInt = 0;
    int minInt = 100000;
    for (const auto& [key, data] : detectorMap_) {
        int hits = data.numberOfHits;
        if (hits > maxInt) {
            maxInt = hits;
        }
        if (hits < minInt) {
            minInt = hits;
        }
    }

    std::cout << "max hits in a bar: " << maxInt << std::endl;
    std::cout << "min hits in a bar: " << minInt << std::endl;

    // Define leaf pairs to analyze
    std::vector<std::pair<std::string, std::string>> leafPairs = { // these were just used as different combinations for my plotting functions below
        {"muon_deviation_angle_", "muon_kinetic_"}, // 0
        {"muon_deviation_angle_", "muon_final_mom_"}, // 1
        {"muon_deviation_angle_", "muon_mass_"}, // 2
        {"muon_deviation_angle_", "muon_mom_ratio_"}, // 3
        {"muon_deviation_angle_", "muon_edeps_"}, // 4
        {"daughter_mass_", "daughter_energy_"}, // 5
        {"daughter_mass_", "daughter_kinetic_"}, // 6
        {"muon_polar_angle_", "muon_kinetic_"}, // 7
        {"muon_polar_angle_", "muon_edeps_"}, // 8
        {"muon_polar_angle_", "muon_est_edeps_"},// 9
        {"muon_polar_angle_", "muon_size_"},// 10
        {"muon_polar_angle_", "muon_edeps_total_event_"},// 11
        {"muon_energy_", "muon_edeps_total_event_"},// 12
        {"muon_energy_", "muon_size_"},// 13
        {"muon_bars_hit_", "muon_edeps_total_event_"},// 14
        {"muon_bars_hit_", "muon_polar_angle_"},// 15
        {"muon_polar_angle_", "muon_bars_hit_"},// 16
        {"muon_polar_angle_", "muon_edeps_per_bar_"},// 17
        {"muon_energy_", "muon_polar_angle_"},//18
        {"muon_polar_angle_", "muon_edeps_per_barmm_"},//19
        {"muon_azimuthal_angle_", "muon_edeps_"},//20
        {"muon_azimuthal_angle_", "muon_edeps_total_event_"},//21
        {"muon_azimuthal_angle_", "muon_bars_hit_"}, //22
        {"muon_azimuthal_angle_", "muon_edeps_per_bar_"}, //23
        {"muon_azimuthal_angle_", "muon_size_"},//24
        {"muon_polar_angle_", "recPE_"},//25
        {"muon_polar_angle_", "muon_mom_diff_"},//26
        {"muon_est_polar_angle_", "muon_est_edeps_"}, //27
        {"muon_energy_", "muon_edeps_"}, //28
        {"muon_polar_angle_", "edeps_est_edeps_ratio_"} // 29
    };

    TDirectory *Plots = file_->mkdir("Plots");
    Plots->cd();

    file_->Write();

    // Here we can make beutified plots 
    // These two funcs need to be configured manually to use. e.g. you might need to change which ttree it is taking data from
    // getPlot2D(leafPairs[8]);
    // getPlot1D(leafPairs[29]);

    if (file_) delete file_;
}

DECLARE_ANALYZER(CosmicsAnalysis);