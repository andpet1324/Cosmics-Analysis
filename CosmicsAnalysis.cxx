// filename: CosmicsAnalysis.cxx

#include "CosmicsAnalysis.h"

void CosmicsAnalysis::buildMap(std::map<int,DetectorInfo>& map) {

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
            DetectorInfo horizontalData = {horizontalXMin, horizontalXMax, horizontalYMin, horizontalYMax, horizontalZMin, horizontalZMax, "horizontal" , 0.0, 0.0};
            DetectorInfo verticalData = {verticalXMin, verticalXMax, verticalYMin, verticalYMax, verticalZMin, verticalZMax, "vertical" , 0.0, 0.0};

            const int horizontalID = std::stoi(horizontalIDString);
            const int verticalID = std::stoi(verticalIDString);

            map.insert({horizontalID, horizontalData});
            map.insert({verticalID, verticalData});
        }
    }
}

void CosmicsAnalysis::onProcessStart() {
    getHistoDirectory(); // forget this -> silently no output histogram

    file_ = new TFile("Analysis/cosmics_data.root", "RECREATE");

    sim_daughter_tree_ = new TTree("Sim_Daughter_Data","data from simulated daughter particles");
    sim_daughter_tree_->Branch("daughter_pdgID_", &daughter_pdgID_, "daughter_pdgID_/I");
    sim_daughter_tree_->Branch("daughter_energy_", &daughter_energy_, "daughter_energy_/D");
    sim_daughter_tree_->Branch("daughter_mass_", &daughter_mass_, "daughter_mass_/D");
    sim_daughter_tree_->Branch("daughter_kinetic_", &daughter_kinetic_, "daughter_kinetic_/D");

    sim_muon_tree_ = new TTree("Sim_Muon_Data", "data from simulated muons");
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

    hcal_hit_tree_ = new TTree("Hcal_Hit_Data", "data from simulated hits");
    hcal_hit_tree_->Branch("recPE_", &recPE_);

    graphFile_ = new TFile("Analysis/sim_graphs.root", "RECREATE");

    // Get start-time
    start_time = std::chrono::high_resolution_clock::now();

    buildMap(detectorMap_);

}

 
// define the parametric line equation
void line(double t, const double *p, double &x, double &y, double &z) {
    // a parametric line is define from 6 parameters but 4 are independent
    // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
    // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
    x = p[0] + p[1]*t;
    y = p[2] + p[3]*t;
    z = t;
 }
  

void CosmicsAnalysis::getAngle(const std::map<int, ldmx::SimParticle, std::less<int>, std::allocator<std::pair<const int, ldmx::SimParticle> > >& sim_particles, const std::vector<ldmx::SimCalorimeterHit, std::allocator<ldmx::SimCalorimeterHit>>& sim_hits) {
    auto particle = sim_particles.at(1);
    std::vector<double> vec_y = {0,-1,0};
    std::vector<double> vec_z = {0,0,1};
    auto start_point = particle.getVertex();
    std::vector<float> end_point;

    for (int i = 1; i < sim_hits.size(); i++){ // i has to be at least 1
        std::vector<float> pos = (sim_hits[sim_hits.size() - i]).getPosition();
        double x = pos[0], y = pos[1], z = pos[2];

        // Find the last point in the detector
        if (std::abs(x) < hcalXMax_ && y > hcalYMin_ && y < hcalYMax_ && (z > hcalZMin_ && z < hcalZMax_)) { // Inside the detector
            end_point = pos;
            break;
        }
        // Last iteration. No detector was found;
        if (i == sim_hits.size() - 1) {
            muon_deviation_angle_ = 365;
            return;
        }

    }

    auto start_mom = particle.getMomentum();
    auto end_mom = particle.getEndPointMomentum();
    muon_init_mom_ = std::hypot(start_mom[0], start_mom[1], start_mom[2]);
    muon_final_mom_ = std::hypot(end_mom[0], end_mom[1], end_mom[2]);
    muon_mom_ratio_ = muon_final_mom_/muon_init_mom_;
    std::vector<double> vec_r = {(end_point[0] - start_point[0]), (end_point[1] - start_point[1]), (end_point[2] - start_point[2])};
    double mag_vec_r = std::hypot(vec_r[0], vec_r[1], vec_r[2]); 
    muon_deviation_angle_ = acos((vec_r[0]*start_mom[0] + vec_r[1]*start_mom[1] + vec_r[2]*start_mom[2])/(mag_vec_r * muon_init_mom_)) * 180/M_PI;
    std::vector<double> start_mom_azimuthal = {start_mom[0], 0, start_mom[2]};
    double mag_start_mom_azimuthal = std::hypot(start_mom_azimuthal[0], start_mom_azimuthal[1], start_mom_azimuthal[2]);
    muon_polar_angle_ = acos((vec_y[0]*start_mom[0] + vec_y[1]*start_mom[1] + vec_y[2]*start_mom[2])/(muon_init_mom_)) * 180/M_PI;  
    muon_azimuthal_angle_ = acos((vec_z[0]*start_mom_azimuthal[0] + vec_z[1]*start_mom_azimuthal[1] + vec_z[2]*start_mom_azimuthal[2])/(mag_start_mom_azimuthal)) * 180/M_PI;

    if (start_mom_azimuthal[0] < 0){
        muon_azimuthal_angle_ = -muon_azimuthal_angle_;
    }
}
void CosmicsAnalysis::getAngleEstimate(double& polarAngleEstimate, double& azimuthalAngleEstimate, std::vector<float> line) {
    // Vectors for reference
    std::vector<double> vec_y = {0,-1,0};
    std::vector<double> vec_z = {0,0,1};

    double lineMagnitude = std::hypot(line[0], line[1], line[2]);
    polarAngleEstimate = acos((vec_y[0]*line[0] + vec_y[1]*line[1] + vec_y[2]*line[2])/(lineMagnitude)) * 180/M_PI;  
    std::vector<double> azimuthalLine{line[0], 0, line[2]};
    double azimuthalLineMagnitude = std::hypot(azimuthalLine[0], azimuthalLine[1], azimuthalLine[2]);
    azimuthalAngleEstimate = acos((azimuthalLine[0]*vec_z[0]) + (azimuthalLine[1]*vec_z[1]) + (azimuthalLine[2]*vec_z[2])/azimuthalLineMagnitude) * 180/M_PI;
}

std::vector<float> CosmicsAnalysis::crossProduct(std::vector<float> vec1, std::vector<float> vec2) {

    // Calculates the cross product of two vectors
  
    float v1 = vec1[0]; float v2 = vec1[1]; float v3 = vec1[2];
    float b1 = vec2[0]; float b2 = vec2[1]; float b3 = vec2[3];
  
    std::vector<float> product = {((v2*b3) - (v3*b2)), ((v3*b1) - (v1*b3)), ((v1*b2) - (v2*b1))};
  
    return product;
}

std::pair<std::vector<std::vector<float>>, std::vector<double>> CosmicsAnalysis::sort(std::pair<std::vector<std::vector<float>>, std::vector<double>> hits) {
    // Create a copy of the original vector
    std::vector<std::vector<float>> positions = hits.first;
    std::vector<double> edeps = hits.second;

    // Create index vector
    std::vector<size_t> indices(positions.size());
    std::iota(indices.begin(), indices.end(), 0);

    // Sort indices based on the y-value (index 1) of positions
    std::sort(indices.begin(), indices.end(),
            [&positions](size_t i1, size_t i2) {
                return positions[i1][1] > positions[i2][1]; // descending by y
            });

    // Create sorted vectors
    std::vector<std::vector<float>> sorted_positions;
    std::vector<double> sorted_edeps;
    sorted_positions.reserve(positions.size());
    sorted_edeps.reserve(edeps.size());

    for (size_t i : indices) {
        sorted_positions.push_back(positions[i]);
        sorted_edeps.push_back(edeps[i]);
    }

    return {sorted_positions, sorted_edeps};
}

std::pair<std::vector<std::vector<float>>, std::vector<double>> CosmicsAnalysis::weightedAverage(const std::vector<ldmx::SimCalorimeterHit, std::allocator<ldmx::SimCalorimeterHit>>& simHits) {
    
    std::vector<std::vector<float>> newPositions;
    std::vector<double> newEnergies;

    // Iterate over the simulated hcal hits
    Int_t index = 0;
    int lastBar = 0;
    std::vector<std::vector<float>> currentHitPosition;
    std::vector<double> currentHitEnergy;
    double energyHit = 0;
    for (const auto& hit : simHits) {
        std::vector<float> pos = hit.getPosition();
        double hitEnergy = (hit.getEdeps())[0];
        float xHit = pos[0];
        float yHit = pos[1];
        float zHit = pos[2];

        // Filter away points that are outside the detector
        if (std::abs(pos[0]) > hcalXMax_ || pos[1] < hcalYMin_ || pos[1] > hcalYMax_ || pos[2] < hcalZMin_ || pos[2] > hcalZMax_) {
            continue;
        }

        // We want to sum the hits in a bar, then take an average position weighted by the energy
        // But let's just start with a normal average
        if (lastBar == 0) {
            currentHitPosition.push_back(pos);
            currentHitEnergy.push_back(hitEnergy);
            lastBar = hit.getID();
            continue;
        }
        if (hit.getID() == lastBar) {
            currentHitPosition.push_back(pos);
            currentHitEnergy.push_back(hitEnergy);
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
        
            newPositions.push_back({xHit, yHit, zHit});
            newEnergies.push_back(energyHit);

            // Clear for next bar
            currentHitPosition.clear();
            currentHitEnergy.clear();
            currentHitPosition.push_back(pos);
            currentHitEnergy.push_back(hitEnergy);

            lastBar = hit.getID();
        }
        // Last element
        if (index == simHits.size() && currentHitPosition.size() != 1){

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

            newPositions.push_back({xHit, yHit, zHit});
            newEnergies.push_back(energyHit);
        }
    }

    return {newPositions, newEnergies};
}
  
std::pair<std::vector<float>,float> CosmicsAnalysis::linearFit(const std::vector<ldmx::SimCalorimeterHit, std::allocator<ldmx::SimCalorimeterHit>>& sim_hits, const std::map<int, ldmx::SimParticle, std::less<int>, std::allocator<std::pair<const int, ldmx::SimParticle>>>& sim_particles) {

    // Note: accepted hits are hits that have met selection criteria and are within the detector

    std::vector<float> bestLine;
    std::vector<float> bestPoint;
    float bestR2;
    
    std::pair<std::vector<std::vector<float>>, std::vector<double>> averageHits = weightedAverage(sim_hits);
    std::pair<std::vector<std::vector<float>>, std::vector<double>> sortedHits = sort(averageHits);
    auto pos = sortedHits.first;
    auto edeps = sortedHits.second;

    if (pos.size() == 0) {
        // std::cout << "Error: Missed the HCal for the currently defined limits, but had hits from simHits" << std::endl;
        return {{-1,-1,-1},1};
    }

    // Original position of the muon when it enters the detector
    for (int i = 0; i < pos.size(); i++) {
        if (edeps[i] > edeps_upper_threshold_ || edeps[i] < edeps_lower_threshold_ || muon_deviation_angle_ > angle_threshold_) {
            continue;
        }
        else {
            muon_origin_ = pos[i];
            break;
        }
    }

    // Checks that we have enough hits to perform our fit. Higher max iterator should give a better fit. 
    int iteratorMax = 3;
    if (pos.size() <= (iteratorMax + 1)){
        // std::cout << "Too few hits for a linear fit" << std::endl;
        return {{-1,-1,-1},1};
    }

    // Loop through the hits and find the last ones in the detector
    bool check = true;
    for (int i = 1; i < pos.size(); i++) {
        // Current position and energy deposit

        const auto ref_pos = pos[pos.size() - i];
        auto ref_edeps = edeps[edeps.size() - i];

        // Checks if we've reached the end of the sim hits
        if (i == pos.size() - 1) {
            // throw std::invalid_argument( "Bingbangbong" );
            // std::cout << "Error: iterated through all the positions" << std::endl;
            return{{-1,-1,-1}, 1};
        }

        // Ignore hit if selection criteria aren't met
        if (ref_edeps > edeps_upper_threshold_ || ref_edeps < edeps_lower_threshold_ || muon_deviation_angle_ > angle_threshold_){
            iteratorMax++; // We don't want this to count as an iteration of accepted hits
            continue;
        }
    
        // Get the line and line magnitude to current point
        std::vector<float> line = {(ref_pos[0] - muon_origin_[0]), (ref_pos[1] - muon_origin_[1]), (ref_pos[2] - muon_origin_[2])};
        float lineMagnitude = std::hypot(line[0], line[1], line[2]);
    
        // R2 to the furthest point from the line, of all accepted points
        float R2Max = 0;

        // Check the R2 to all accepted points
        // also check R².
        for (int k = 0; k < pos.size(); k++) {
            
            // Selection criteria for points
            if (edeps[k] > edeps_upper_threshold_ || edeps[k] < edeps_lower_threshold_ || muon_deviation_angle_ > angle_threshold_){
                continue;
            }
    
            if (edeps[k] <= 0.000000001) continue; // I don't think we actually have values that are zero
            //                                      but this is here just in case.
            //                                      0.000000001 is 0.001 eV
    
            // vector from origin to point
            std::vector<float> vector = {(pos[k][0] - muon_origin_[0]), (pos[k][1] - muon_origin_[1]), (pos[k][2] - muon_origin_[2])};
    
            // Check the R2 between line and points
            // The R2 is the magnitude of the cross product of vector and line divided by the line magnitude.
            std::vector<float> product = crossProduct(vector, line);
            float productMagnitude = std::hypot(product[0], product[1], product[2]);
            float R2 = productMagnitude/lineMagnitude;

            // Checks if current R2 is larger than the previously largest R2
            if (R2 > R2Max){
                R2Max = R2;
            }
        }
    
        // We don't except perfect lines
        if (R2Max == 0) {
            continue;
        }

        // For first iteration
        if (check){
            bestR2 = R2Max;
            bestLine = line;
            bestPoint = ref_pos;
            check = false;
        }
        // For second iteration and onwards
        else if (R2Max < bestR2) {
            bestR2 = R2Max;
            bestLine = line;
            bestPoint = ref_pos;
        }
        // For last iteration
        if (i >= iteratorMax) {
            break;
        }
    }

    if (bestR2 <= 0.0000001) {
        // std::cout << "The best R2 from the fitted line to a point is very small! It is: " << bestR2 << std::endl;
        return {{-1,-1,-1},1};
    }

    // Returns the line and R2, where 
    // line = (ax,by,cz) and R2 is the max R2 from the line to a point
    return {bestLine, bestR2};
}

void CosmicsAnalysis::linearFitR2(const std::vector<ldmx::SimCalorimeterHit, std::allocator<ldmx::SimCalorimeterHit>>& sim_hits, const std::map<int, ldmx::SimParticle, std::less<int>, std::allocator<std::pair<const int, ldmx::SimParticle>>>& sim_particles, std::vector<float>& a, double& b) {

    // Note: accepted hits are hits that have met selection criteria and are within the detector
    std::vector<float> bestPoint;
    
    std::pair<std::vector<std::vector<float>>, std::vector<double>> averageHits = weightedAverage(sim_hits);
    std::pair<std::vector<std::vector<float>>, std::vector<double>> sortedHits = sort(averageHits);
    auto pos = sortedHits.first;
    auto edeps = sortedHits.second;

    if (pos.size() == 0) {
        // std::cout << "Error: Missed the HCal for the currently defined limits, but had hits from simHits" << std::endl;
        return;
    }

    // Original position of the muon when it enters the detector
    for (int i = 0; i < pos.size(); i++) {
        if (edeps[i] > edeps_upper_threshold_ || edeps[i] < edeps_lower_threshold_ || muon_deviation_angle_ > angle_threshold_) {
            continue;
        }
        else {
            muon_origin_ = pos[i];
            break;
        }
    }

    // Less hits generally give very wide estimated polar angle deviations
    if (pos.size() <= 5) {
        // std::cout << "Too few hits for a linear fit" << std::endl;
        return;
    }

    // Checks that we have enough hits to perform our fit. Higher max iterator should give a better fit. 
    int iteratorMax = 10;
    if (pos.size() <= (iteratorMax + 1)){
        iteratorMax = pos.size()/2;
    }

    float bestSSE;
    float bestSST;

    // Loop through the hits and find the last ones in the detector
    bool check = true;
    for (int i = 1; i < pos.size(); i++) {
        // Current position and energy deposit

        const auto ref_pos = pos[pos.size() - i];
        auto ref_edeps = edeps[edeps.size() - i];

        // Checks if we've reached the end of the sim hits
        if (i == pos.size() - 1) {
            // std::cout << "Iterated through all the positions" << std::endl;
            // return{{-1,-1,-1}, 1};
            return;
        }

        // Ignore hit if selection criteria aren't met
        if (ref_edeps > edeps_upper_threshold_ || ref_edeps < edeps_lower_threshold_ || muon_deviation_angle_ > angle_threshold_){
            iteratorMax++; // We don't want this to count as an iteration of accepted hits
            continue;
        }
    
        // Get the line and line magnitude to current point
        std::vector<float> line = {(ref_pos[0] - muon_origin_[0]), (ref_pos[1] - muon_origin_[1]), (ref_pos[2] - muon_origin_[2])};
        float lineMagnitude = std::hypot(line[0], line[1], line[2]);

        double meanX = 0;
        double meanY = 0;
        double meanZ = 0;
        int n = 0;

        for (int l = 0; l < pos.size(); l++) {
            
            // Selection criteria for points
            if (edeps[l] > edeps_upper_threshold_ || edeps[l] < edeps_lower_threshold_ || muon_deviation_angle_ > angle_threshold_){
                continue;
            }
    
            if (edeps[l] <= 0.000000001) continue; // I don't think we actually have values that are zero
            //                                      but this is here just in case.
            //                                      0.000000001 is 0.001 eV
    
            // Data points
            float xPrime = pos[l][0];
            float yPrime = pos[l][1];
            float zPrime = pos[l][2];

            // Add to mean
            meanX += xPrime;
            meanY += yPrime;
            meanZ += zPrime;
            n += 1;
        }

        meanX = meanX/n;
        meanY = meanY/n;
        meanZ = meanZ/n;
    
        float SSE = 0;
        float SST = 0;

        std::vector<float> unitLine{line[0]/lineMagnitude, line[1]/lineMagnitude, line[2]/lineMagnitude};

        // Check R².
        for (int k = 0; k < pos.size(); k++) {
            
            // Selection criteria for points
            if (edeps[k] > edeps_upper_threshold_ || edeps[k] < edeps_lower_threshold_ || muon_deviation_angle_ > angle_threshold_){
                continue;
            }
    
            if (edeps[k] <= 0.000000001) continue; // I don't think we actually have values that are zero
            //                                      but this is here just in case.
            //                                      0.000000001 is 0.001 eV
    
            // Data points
            float xPrime = pos[k][0];
            float yPrime = pos[k][1];
            float zPrime = pos[k][2];

            // Line from origin to point
            std::vector<float> point = {(xPrime-muon_origin_[0]), (yPrime-muon_origin_[1]), (zPrime-muon_origin_[2])};

            // length unit t
            auto t = ((point[0]*unitLine[0]) + (point[1]*unitLine[1]) + (point[2]*unitLine[2]));

            if (t == 0) {
                continue;
            }

            // Projection of points on the line
            std::vector<float> linePoint = {
                muon_origin_[0] + unitLine[0] * t,
                muon_origin_[1] + unitLine[1] * t,
                muon_origin_[2] + unitLine[2] * t
            };

            float x = linePoint[0];
            float y = linePoint[1];
            float z = linePoint[2];

            float dist = (std::pow((xPrime-x),2) + std::pow((yPrime-y),2) + std::pow((zPrime-z),2));

            if (dist > 300) { // We don't want to include points that are to far away from the line
                continue;
            }

            SSE += dist;
            SST += (std::pow((xPrime-meanX),2) + std::pow((yPrime-meanY),2) + std::pow((zPrime-meanZ),2));
        }

        float R2 = 1 - (SSE/SST);

        // For first iteration
        if (check){
            b = R2;
            a = line;
            bestPoint = ref_pos;
            check = false;
            bestSSE = SSE;
            bestSST = SST;
        }
        // For second iteration and onwards
        else if (R2 > b) {
            b = R2;
            a = line;
            bestPoint = ref_pos;
            bestSSE = SSE;
            bestSST = SST;
        }
        // For last iteration
        if (i >= iteratorMax) {
            break;
        }
    
    }

    // std::cout << "R2: " << b << std::endl;
    // std::cout << "best SSE: " << bestSSE << std::endl;
    // std::cout << "best SST: " << bestSST << std::endl;

    // if (bestR2 <= 0.0000001) {
    //     std::cout << "The best R2 from the fitted line to a point is very small! It is: " << bestR2 << std::endl;
    //     return {{-1,-1,-1},1};
    // }
}

void CosmicsAnalysis::getPlot1D(std::pair<std::string, std::string> pair) {

    TLeaf* xLeaf = sim_muon_tree_->GetLeaf((pair.first).c_str());
    TLeaf* yLeaf = sim_muon_tree_->GetLeaf((pair.second).c_str());
    // TLeaf* yLeaf = hcal_hit_tree_->GetLeaf((pair.second).c_str());
    Long64_t nEntries = sim_muon_tree_->GetEntries();

    if (!xLeaf || !yLeaf) {
        std::cerr << "Error: One or both leaves not found!" << std::endl;
        return;
    }

    // If the yLeaf is a vector

    // Get max value in the yLeaf
    double max_val_y = 0;
    double min_val_y = 1000000;
    double max_val_x = 0;
    double min_val_x = 1000000;

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
    TH1F* H1 = new TH1F("Data", "", 100, 0, 75);
    // TH2F* H2 = new TH2F("Data", "", 50, 0.0, 75, 30, 0, 5000);

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
    

    // Set the histograms stat box location
    gStyle->SetStatX(0.85);
    gStyle->SetStatY(0.90);
    gStyle->SetOptStat(11);

    // Set the color of the colz
    // gStyle->SetStatBorderSize(0); // Removes border of the statistics from the hist

    H1->GetXaxis()->SetTitle("Polar angle #theta [deg]");
    H1->GetYaxis()->SetTitle("Counts");
    H1->GetYaxis()->SetNoExponent(kFALSE);
    H1->GetYaxis()->SetNdivisions(505, kTRUE);
    H1->GetXaxis()->SetNdivisions(505, kTRUE);

    // Set log scale on z
    // gPad->SetLogz(1);

    // Set log scale on x,y
    // gPad->SetLogx(1);
    // gPad->SetLogy(1);

    // If y is a vector

    // double x;
    // std::vector<float>* y = nullptr;

    // sim_muon_tree_->SetBranchAddress((pair.first).c_str(), &x);
    // sim_muon_tree_->SetBranchAddress((pair.second).c_str(), &y);
    // // hcal_hit_tree_->SetBranchAddress((pair.second).c_str(), &y);


    // for (Long64_t i = 0; i < nEntries; i++) {
    //     sim_muon_tree_->GetEntry(i);
    //     // hcal_hit_tree_->GetEntry(i);
    
    //     // Loop over elements in the vector
    //     for (float y_val : *y) {
    //         // if (y_val > 1000 || y_val < -1000) {continue;}
    //         // if (y_val > 1) continue;
    //         H1->Fill(y_val, 1);
    //     }
    // }
    
    // If y is not a vector

    for (Long64_t i = 0; i < nEntries; i++) {
        sim_muon_tree_->GetEntry(i);
        double x = xLeaf->GetValue();
        double y = yLeaf->GetValue();
        
        // if (y <= 0.000000001) continue; // Not a hit. 0.000000001 is 0.001 eV
        if (y <= 0.0001) {
            continue;
        }
        H1->Fill(x,1);
    }

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

    H1->Draw();
    T1->SaveAs("output.pdf");
    H1->Write();
    delete H1;
    delete T1;
}

void CosmicsAnalysis::getPlot2D(std::pair<std::string, std::string> pair) {

    TLeaf* xLeaf = sim_muon_tree_->GetLeaf((pair.first).c_str());
    TLeaf* yLeaf = sim_muon_tree_->GetLeaf((pair.second).c_str());
    Long64_t nEntries = sim_muon_tree_->GetEntries();

    if (!xLeaf || !yLeaf) {
        std::cerr << "Error: One or both leaves not found!" << std::endl;
        return;
    }


    // If the yLeaf is a vector

    // Get max value in the yLeaf
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
    TH2F* H2 = new TH2F("Data", "", 100, 0, 75, 75, 0, 20);

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
    H2->GetYaxis()->SetTitle("Estimated energy in bars [MeV]");
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
    /////////////////////////////UNCOMMENT/////////////////////////////

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
    ///////////////////////////UNCOMMENT/////////////////////////////
    
    // If y is not a vector
    /////////////////////////////UNCOMMENT/////////////////////////////

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

    /////////////////////////////UNCOMMENT/////////////////////////////    

    // Draw the histogram
    H2->Draw("COLZ");
    T1->SaveAs("output.pdf");
    H2->Write();
    delete H2;
    delete T1;
}

void CosmicsAnalysis::getPlot2D3Entries(std::string str1, std::string str2 , std::string str3) {

    TLeaf* xLeaf = sim_muon_tree_->GetLeaf((str1).c_str());
    TLeaf* yLeaf = sim_muon_tree_->GetLeaf((str3).c_str());
    TLeaf* thirdLeaf = sim_muon_tree_->GetLeaf((str2).c_str());

    Long64_t nEntries = sim_muon_tree_->GetEntries();

    if (!xLeaf || !yLeaf) {
        std::cerr << "Error: One or both leaves not found!" << std::endl;
        return;
    }

    // If the yLeaf is a vector

    // Get max value in the yLeaf
    double max_val_y = 0;
    double max_val_x = 0;

    for (Long64_t i = 0; i < nEntries; i++) {
        sim_muon_tree_->GetEntry(i);
        double y = yLeaf->GetValue(); 
        double x = xLeaf->GetValue();
        double third = thirdLeaf->GetValue();
        
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
    TH2F* H2 = new TH2F("Data", "", 100, 0, 75, 100, 0, 4);

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
    H2->GetYaxis()->SetTitle("Detected energy horizontal bars per mm [MeV/mm]");
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
    /////////////////////////////UNCOMMENT/////////////////////////////

    double x;
    double third;
    std::vector<double>* y = nullptr;


    sim_muon_tree_->SetBranchAddress((str1).c_str(), &x);
    sim_muon_tree_->SetBranchAddress((str3).c_str(), &y);
    sim_muon_tree_->SetBranchAddress((str2).c_str(), &third);

    for (Long64_t i = 0; i < nEntries; i++) {
        sim_muon_tree_->GetEntry(i);

        // Selection criteria on third
        if (abs(third) < 80 || abs(third) > 100) {
            continue;
        }
    
        // Loop over elements in the vector
        for (float y_val : *y) {
            // std::cout << std::setprecision(30) << y_val << std::endl;
            if (y_val <= 0.000000001) continue; // Not a hit. 0.000000001 is 0.001 eV
            // if (y_val > 1000 || y_val < -1000) {continue;}
            // if (y_val < 2) continue;

            H2->Fill(x, y_val);
        }
    }
    ///////////////////////////UNCOMMENT/////////////////////////////
    
    // If y is not a vector
    /////////////////////////////UNCOMMENT/////////////////////////////

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

    /////////////////////////////UNCOMMENT/////////////////////////////    

    // Draw the histogram
    H2->Draw("COLZ");
    T1->SaveAs("output.pdf");
    H2->Write();
    delete H2;
    delete T1;
}

bool CosmicsAnalysis::isInDiscreteRange(int value) {
    const int start = 1142;
    const int interval_gap = 49;
    const int interval_width = 25;
    const int max_value = 5544;

    for (int current_start = start; current_start <= max_value; current_start += (interval_gap + interval_width)) {
        int current_end = current_start + interval_width;
        if (value >= current_start && value <= current_end) {
            return true;
        }
    }
    return false;
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

            // if (currentLayer_ != layer || currentStrip_ != strip) {
            //     std::cout << "detectorID layer: " << currentLayer_ << " map layer: " << layer << std::endl;
            //     std::cout << "detectorID strip: " << currentStrip_ << " map strip: " << strip << std::endl;
            //     std::cout << "correlation error!" << std::endl;
            // }

            // std::cout << "x: " << x << " y: " << y - supportBox << " z: " << z << std::endl;
            // std::cout << ID << std::endl;
            map[ID].energy += edep; // add the energy to the bar
            // std::cout << "energy after adding: " << map[ID].energy << std::endl;
            // std::cout << "from ID: " << map[ID].orientation << std::endl;
            return;
         }
    }
}

void CosmicsAnalysis::AddDaughterEnergy(std::map<int,DetectorInfo>& map, const std::vector<double>& vertex, const std::vector<double>& endPoint, const double& edep) {

    const float x = vertex[0];
    const float y = vertex[1];
    const float z = vertex[2];

    std::vector<float> v{x,y,z};

    const float endX = vertex[0];
    const float endY = vertex[1];
    const float endZ = vertex[2];

    // Check if the particle gets absorbed in an absorber

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
    // get detector id by position
    const float x = position[0];
    const float y = position[1];
    const float z = position[2];

    auto map = detectorMap_;

    // iterate through ids, 010101, 020201, 010301, 

    int numLayers = 96;

    for (int i = 1; i <= numLayers; i++ ) {
        std::string strID  = ((i % 2 == 0) ? std::string("20") : std::string("10")) + (i < 10 ? "0" : "") + std::to_string(i) + std::string("01");
        int ID = std::stoi(strID);
        DetectorInfo detectorData = map[ID];
        int absoZMin = detectorData.zMin - 2 - 25;
        int absoZMax = detectorData.zMin - 2;

        if (z > absoZMin && z < absoZMax) {
            return false;
        }
    }
    return true;

    // method that checks if the pos is in an absorber. It sucks
    // for (auto it = map.begin(); it != map.end(); it++) {
    //     int ID = it->first;
    //     std::string strID = std::to_string(ID);
    //     int strip = std::stoi(strID.substr(4,5));

    //     if (strip != 1) {
    //         continue;
    //     }
    //     std::cout << count << std::endl;
    //     DetectorInfo detectorData = it->second;
    //     int absoZMin = detectorData.zMin - 2 - 25;
    //     int absoZMax = detectorData.zMin - 2;

    //     if (z > absoZMin && z < absoZMax) {
    //         return false;
    //     }
    // }
    // return true;

    // This checks if the position is specifically in a bar
    // It would take less time to iterate through the absorbers.
    for (const auto& pair : map) {
        DetectorInfo detectorData = pair.second;
        if ((x >= detectorData.xMin && x < detectorData.xMax)
         && (y >= detectorData.yMin && y < detectorData.yMax)
         && (z >= detectorData.zMin && z <= detectorData.zMax)) {
            return true;
         }
    }
    return false;
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


void CosmicsAnalysis::analyze(const framework::Event& event) {

    // Get the data collections from the event bus
    const auto& hcal_sim_hits{event.getCollection<ldmx::SimCalorimeterHit>("HcalSimHits")};
    const auto& sim_particles{event.getMap<int, ldmx::SimParticle>("SimParticles", "cosmics")};
    const auto& hcal_rec_hits{event.getCollection<ldmx::HcalHit>("HcalRecHits")};

    // Get the event number
    int eventNumber = event.getEventNumber();

    if (eventNumber % 10000 == 0) {
        std::cout << "event: " << eventNumber << " started" << std::endl;
    }

    std::cout << "event: " << eventNumber << std::endl;

    // Some events do not interact with the detector. We want to find these events and store that event number
    // todo: change the way this is handled. It is weird
    if (hcal_sim_hits.size() == 0) { // This event doesn't hit the HCal, so we ignore it for data collection
        nr_events_miss_++;
        // std::cout << "Event " << eventNumber << " missed the HCal" << std::endl;
        return;
    }

    // Add to total events
    total_events_ += 1;

    // Finds the amount of particles in the event (muon + daughters)
    auto nr_particles = sim_particles.size();

    // line for fitting later
    std::vector<float> line;
    double R2;

    // map with energies
    std::map<int, std::map<float, double> > hitMap;

    std::map<int,DetectorInfo> eventMap;

    // build detector map
    buildMap(eventMap);

    std::vector<double> energies;
    std::vector<double> estimatedEnergies;

    int missingParticles = 0;

    for (const auto& [trackID, particle]: sim_particles){
        
        // trackID == 1 is the initial muon
        if (trackID == 1) {
            muon_size_ = 1;
            // Check that muon doesn't END in the detector!
            auto end_point = particle.getEndPoint();
            if ((abs(end_point[0] < hcalXMax_)) && (end_point[1] < hcalYMax_ && end_point[1] > hcalYMin_) && (end_point[2] < hcalZMax_ && end_point[2] > 240)) {
                decayed_muons_ += 1;
                return;
            }

            // If size <= 1 then we dont have any hits in the detector or enough hits to conduct analysis
            if (hcal_sim_hits.size() <= 1) {
                // std::cout << "Event " << eventNumber << " missed the detector and is being ignored" << std::endl;
                return;
            }

            // Finds the angle between initial momentum and final position
            getAngle(sim_particles, hcal_sim_hits);
            
            // Now we choose which events to ignore
            // We define a "cone of straight tracks" to be within a certain angle threshold

            if (muon_deviation_angle_ < angle_threshold_) { // The track is straight per this definition.
                straight_tracks_ += 1;
            }
            
            muon_energy_ = particle.getEnergy();
            muon_kinetic_ = particle.getKineticEnergy();
            muon_mass_ = particle.getMass();

            // std::vector<double> energies;
            muon_edeps_total_event_ = 0;
            ldmx::HcalID lastDetectorID;
            muon_bars_hit_ = 0;

            double polarAngleEstimate;
            double azimuthalAngleEstimate;

            // Just a check if we want to run the linefit or not
            if (linecheck) {

                // const std::pair<std::vector<float>,float> line_fit = linearFit(hcal_sim_hits, sim_particles);
                linearFitR2(hcal_sim_hits, sim_particles, line, R2);

                if (R2 < 0.999) {
                    // std::cout << "R2 too small!" << std::endl;
                    return;
                }
    
                if (line[0] == -1 && line[1] == -1 && line[2] == -1){
                    // std::cout << "Event " << eventNumber << " has a line-fit that gives and error. No linefit will be done" << std::endl;
                    return;
                }
                
                if (R2 < 0) {
                    std::cout << R2 << std::endl;
                    throw std::invalid_argument( "received negative value" );
                }

                // Storing the largest deviation
                // if (largest_polar_angle_estimate_deviation_ < abs(muon_polar_angle_ - polarAngleEstimate)) {
                //     specific_event_ = eventNumber;
                //     largest_polar_angle_estimate_deviation_ = abs(muon_polar_angle_ - polarAngleEstimate);
                // }

                // Based on the line, we want to estimate how many bars were hit, and which ones.
                // Use the line to go through all the bars that it intersects. Assuming that the muon deposits energy in
                // every bar that it goes through, we can estimate how long it would take to hit each bar.

                // float stepSize = 0; // in mm
                // std::vector<float> initPos = muon_origin_;
                // while (true) {
                //     std::vector<float> pos = {
                //         initPos[0] + line[0] * stepSize,
                //         initPos[1] + line[1] * stepSize,
                //         initPos[2] + line[2] * stepSize
                //     };
                //     stepSize += 0.1;
                    
                //     if (abs(pos[0]) > 1000 || abs(pos[1]) > 1000 || pos[2] < 840 || pos[2] > hcalZMax_) break;

                // }
            }

            auto initPos = hcal_sim_hits[0].getPosition();
            auto finalPos = (hcal_sim_hits.back()).getPosition();
            std::vector<float> diff{(finalPos[0] - initPos[0]), (finalPos[1] - initPos[1]), (finalPos[2] - initPos[2])};
            float detectorLength = std::hypot(diff[0], diff[1], diff[2]);
            detectorLength = 2000;

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

                // We don't care about products in the ecal or side-hcal for now
                // todo: method is bad, fix!
                // REMEMBER! Side-hcal is 1500x1500, not 1000x1000
                if (std::abs(position[0]) < 1500 && std::abs(position[1]) < 1519.05 && position[2] < 840) {
                    nr_events_ecal_++;
                    // std::cout << "Event " << eventNumber << " has hits in the ECal or side-HCal and is being ignored" << std::endl;
                    return;
                }

                // Only saves point close to the line
                if (linecheck) {
                    std::vector<float> vector = {(position[0] - muon_origin_[0]), (position[1] - muon_origin_[1]), (position[2] - muon_origin_[2])};
                    std::vector<float> product = crossProduct(vector, line);
                    float productMagnitude = std::hypot(product[0], product[1], product[2]);
                    float lineMagnitude = std::hypot(line[0], line[1], line[2]);
                    float distance = productMagnitude/lineMagnitude;
                    if (distance > allowedLineDistance_) {
                        continue; // hit too far from the line.
                    }
                }

                auto detectorID = hit.getID();
                auto id = ldmx::HcalID(detectorID);

                // 1 for vertical, 0 for horizontal
                // even layers vertical, odd horizontal.
                bar_orientation_ = id.layer() % 2 == 0 ? 1 : 0;

                auto time = hit.getTime();

                auto hitEdeps = (hit.getEdeps())[0]; // Get the hit

                if (hitEdeps <= 0.000000001) continue; // I don't think we actually have values that are zero
                //                                      but this is here just in case.
                //                                      0.000000001 is 0.001 eV

                // add the energy to the right bar
                if (chose_orientation_ == 2 || bar_orientation_ == chose_orientation_) {
                    addToMap(eventMap, position, hitEdeps);
                }
                
                // Adds the total energy of the event
                muon_edeps_total_event_ += hitEdeps;

            }

            if (linecheck) {
                // get angle estimates
                getAngleEstimate(polarAngleEstimate, azimuthalAngleEstimate, line);

                if (largest_polar_angle_estimate_deviation_ < abs(muon_polar_angle_ - polarAngleEstimate)) {
                    specific_event_ = eventNumber;
                    largest_polar_angle_estimate_deviation_ = abs(muon_polar_angle_ - polarAngleEstimate);
                }

                float lineMagnitude = std::hypot(line[0], line[1], line[2]);
                std::vector<float> lineNorm{(line[0]/lineMagnitude), (line[1]/lineMagnitude), (line[2]/lineMagnitude)};

                // we don't actually use the angles for the energy estimation. We instead iterate over the line step-wise
                int step = 0;
                double stepsize = 1; // shorter stepsize will give more accurate estimation
                double length = 0;

                while (true) {
                    float x = muon_origin_[0] + lineNorm[0]*step;
                    float y = muon_origin_[1] + lineNorm[1]*step;
                    float z = muon_origin_[2] + lineNorm[2]*step;
                    std::vector<float> position{x,y,z};
                    step += stepsize;

                    if (y < hcalYMin_) {
                        break;
                    }

                    // check if the point is inside a scintillator
                    // this process is slow but works
                    // if (IsInBar(position)) {
                    //     length += stepsize;
                    // }

                    // method that adds the estimated energy to each bar so that we can plot it later.
                    // very slow
                    AddEstimatedEnergy(eventMap, position, stepsize);
                }
            }

            muon_mom_diff_ = (muon_init_mom_ - muon_final_mom_)/detectorLength;

            // if (abs(muon_mom_diff_) < 0.001) return;

            // We calculate the average length traveled in one bar from the incident angle and that the
            // bar width is 20mm.

            double angle_rad = muon_polar_angle_*(M_PI/180);
            double barLength = 20/(std::sin(angle_rad));

            muon_edeps_per_bar_ = muon_edeps_total_event_/muon_bars_hit_;
            muon_edeps_per_barmm_ = muon_edeps_per_bar_/barLength;

            // muon_edeps_ = energies;

            auto daughters = particle.getDaughters();

            // std::cout << daughters[0].getPdgID() << std::endl;

            for (const auto& daughter : daughters) {
                auto it = sim_particles.find(daughter);
                if (it != sim_particles.end()) {
                    const auto& particle = it->second;

                    auto vertex = particle.getVertex();
                    auto endPoint = particle.getEndPoint();
                    auto kinetic = particle.getKineticEnergy();

                    // AddDaughterEnergy(eventMap, vertex, endPoint, kinetic);
                } else {
                    missingParticles++;
                    continue;
                    // handle the case where the key is not found
                    // std::cerr << "SimParticle with ID " << daughters[0] << " not found.\n";
                }
            }



            // if (sim_particles.size() < daughters.size()) return;
            // try { // For some reason, the sim_particles map is smaller than the amount of daughter particles sometimes. This makes no sense at all and is presumably a bug in ldmx-sw
            //     auto particle = sim_particles.at(daughters.size());
            // }
            // catch (const std::out_of_range& e) {
            //     return;
            // }

            // for (const auto& daughter : daughters) { // want to add the daughter energy here
            //     // std::cout << daughter << std::endl;
            //     try { // some particles aren't defined I guess?? I don't understand why this is a problem
            //         auto particle = sim_particles.at(daughter);
            //     }
            //     catch (const std::out_of_range& e) {
            //         // std::cout << "error" << std::endl;
            //         continue;
            //     }
            //     auto vol = particle.getVertexVolume();
            //     if ((vol == "back_hcal_scintYVolume" || vol == "back_hcal_scintXVolume")) {
            //         auto time = particle.getTime(); // time the daughter particle was created. Should correlate to the time of a hit
            //         auto kinetic = particle.getEnergy();
            //         bool add = false;
            //         for (const auto& firstPair : hitMap) {
            //             if (add) break; 
            //             int it = 0;
            //             for (const auto& secondPair : firstPair.second){
            //                 if ( abs(time - secondPair.first) < 0.001) {
            //                     std::cout << it++ << std::endl;
            //                     std::cout << "before: " << hitMap[firstPair.first][secondPair.first] << std::endl;
            //                     std::cout << "kinetic: " << kinetic << std::endl;
            //                     // hitMap[firstPair.first][secondPair.first] += kinetic;
            //                     std::cout << "after: " << hitMap[firstPair.first][secondPair.first] << std::endl;
            //                     add = true;
            //                     break;
            //                 }
            //             }
            //         }
            //     }
            // }
            // std::vector<double> Edeps;
            // // sum up to energies
            // for (const auto& firstPair : hitMap) {
            //     float barEnergy = 0;
            //     for (const auto& secondPair : firstPair.second) {
            //         barEnergy += secondPair.second;
            //     }
            //     Edeps.push_back(barEnergy);
            // }
            
            // muon_edeps_ = energies;
            // sim_muon_tree_->Fill();
        }
        else { // Not a muon
            // auto vol = particle.getVertexVolume();
            // auto parents = particle.getParents();
            auto vertex = particle.getVertex();
            auto endPoint = particle.getEndPoint();
            auto kinetic = particle.getKineticEnergy();

            // AddDaughterEnergy(eventMap, vertex, endPoint, kinetic);

            // if (parents[0] == 1) {// parent is muon

            // }
            // if ((vol == "back_hcal_scintYVolume" || vol == "back_hcal_scintXVolume") && parents[0] == 1) { // this particle is produced in a scint bar, and its parent is a muon
            //     // Check that it decays in an absorber
            //     auto endZ = (particle.getEndPoint())[2];
            //     if (isInDiscreteRange(endZ)) { //checks that the particle is in an absorber layer
            //         auto time = particle.getTime();
            //         auto kinetic = particle.getEnergy();
            //         bool add = false;
            //         for (const auto& firstPair : hitMap) {
            //             auto timeEnergy = firstPair.second;
            //             if (add) break; 
            //             for (const auto& secondPair : timeEnergy){
            //                 if ( abs(time - secondPair.first) < 0.0001) {
            //                     hitMap[firstPair.first][secondPair.first] += kinetic;
            //                     add = true;
            //                     break;
            //                 }
            //             }
            //         }
            //         if (add != true) {
            //             std::cout << time << std::endl;
            //             throw "ERROR";
            //         }
            //     }
            // }
        }

        // std::vector<double> energies;
        // sum up to energies
        // for (const auto& firstPair : hitMap) {
        //     float barEnergy = 0;
        //     for (const auto& secondPair : firstPair.second) {
        //         barEnergy += secondPair.second;
        //     }
        //     energies.push_back(barEnergy);
        // }

        // for (const auto& pair : eventMap) {
        //     int ID = pair.first;
        //     DetectorInfo data = pair.second;
        //     double energy = data.energy;
        //     if (energy == 0.0) continue; // unnessecary to add energies

        //     energies.push_back(energy);
        // }

        // muon_edeps_ = energies;
        // sim_muon_tree_->Fill();
        
        daughter_pdgID_ = particle.getPdgID();
        daughter_energy_ = particle.getEnergy();
        daughter_mass_ = particle.getMass();
        
        // sim_daughter_tree_->Fill();
    }

    for (const auto& pair : eventMap) {
        int ID = pair.first;
        DetectorInfo data = pair.second;
        double energy = data.energy;
        double estEnergy = data.estimatedEnergy;
        if (energy == 0.0 && estEnergy == 0.0) continue; // unnessecary to add energies
    
        estimatedEnergies.push_back(estEnergy);
        energies.push_back(energy);
    }
    
    muon_est_edeps_ = estimatedEnergies;
    muon_edeps_ = energies;
    sim_muon_tree_->Fill();

    // std::cout << "the amount of missing particles in the sim_particles list are: " << missingParticles << std::endl;

    std::vector<float> PE;

    for (const auto& hit : hcal_rec_hits) {
        PE.push_back(hit.getPE());
    }

    recPE_ = PE;
    hcal_hit_tree_->Fill();

    if (graphcheck) {
        // Declare graph and canvas
        canvas_ = new TCanvas();
        graph_ = new TGraph2D();

        graph_->SetDirectory(nullptr);

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
                if (graphcheck) {
                    delete graph_;
                }
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
            if (energyHit < edeps_upper_threshold_ && energyHit > edeps_lower_threshold_ && muon_deviation_angle_ < angle_threshold_){
                graph_->SetPoint(index++, xHit, yHit, zHit);
            }
            // hcal_hit_tree_->Fill();
        }

        if (linecheck) {
            // Declare polyline
            polyline_fit_ = new TPolyLine3D(2);
            polyline_fit_->SetLineColor(kRed);

            // Variables for filling the line
            double x, y, z;
            double t_min = -10, t_max = 10;  // Range for t

            for (int i = 0; i < 2; i++) {
                double t = t_min + (t_max - t_min) * i; // Line parameter
                x = muon_origin_[0] + line[0]*t;
                y = muon_origin_[1] + line[1]*t;
                z = muon_origin_[2] + line[2]*t;
                polyline_fit_->SetPoint(i,x,y,z);
            }
            
        } 

        if (int nPoints = graph_->GetN() != 0) {
            // Sets the corners of our detector. Just for visual representation
            // graph_->SetPoint(index, hcalXMin_, -hcalYMin_, 0);
            // graph_->SetPoint(index+1, hcalXMax_, hcalYMax_, hcalZMax_);

            // Draw and update to canvas
            graph_->Draw("p0");
            if (linecheck) polyline_fit_->Draw("same");
            canvas_->Update();

            // Save the canvas
            graphFile_->cd();
            canvas_->Write(Form("canvas_%d", eventNumber));
        }

        // Delete so we don't get memory issues
        delete graph_;
        delete canvas_;
        // if (linecheck) delete polyline_fit_;
    }
    
    // Print event and time
    auto current_time = std::chrono::high_resolution_clock::now();
    auto time_diff = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time) / 1000.0;
    // std::cout << "Event " << eventNumber << " done! Time elapsed: " << std::fixed << std::setprecision(1) << time_diff.count() << "s" << std::endl;
}

void CosmicsAnalysis::onProcessEnd(){

    std::cout << "--------------------------Finished---------------------------------" << std::endl;
    std::cout << "Largest estimated polar angle deviation: " << largest_polar_angle_estimate_deviation_ << " for event: " << specific_event_ << std::endl;
    std::cout << "Total events (in the HCal): " << total_events_ << std::endl;
    std::cout << "Muons that don't hit the HCal: " << nr_events_miss_ << std::endl;
    std::cout << "Muons that hit the ECal or side-HCal: " << nr_events_ecal_ << std::endl;

    float decay_ratio = static_cast<float>(decayed_muons_)/total_events_;
    std::cout << "Ratio of muons decaying: " << decay_ratio << std::endl;
    float straight_ratio = static_cast<float>(straight_tracks_)/total_events_;
    std::cout << "Ratio of straight tracks, given angle threshold: " << straight_ratio << std::endl;

    // Check how many hits we have in each scintillator bar

    std::cout << "number of bars hit: " << hcalBars_.size() << std::endl;
    int maxInt = 0;
    int minInt = 100000;
    for (const auto& pair : hcalBars_) {
        if (pair.second > maxInt) {
            maxInt = pair.second;
        }
        if (pair.second < minInt) {
            minInt = pair.second;
        }
    }
    std::cout << "max hits in a bar: " << maxInt << std::endl;
    std::cout << "min hits in a bar: " << minInt << std::endl;

    // Define leaf pairs to analyze
    std::vector<std::pair<std::string, std::string>> leafPairs = {
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
        {"muon_polar_angle_", "muon_mom_diff_"}//26
    };

    TDirectory *Plots = file_->mkdir("Plots");
    Plots->cd();

    for (const auto& pair : leafPairs) {
        std::string xLeaf = pair.first;
        std::string yLeaf = pair.second;

        // Define and save to Canvas
        TCanvas *c1 = new TCanvas("Plot", "Scatter Plot", 1000, 800);
        c1->SetLeftMargin(0.15);
        if (sim_muon_tree_->GetLeaf(xLeaf.c_str()) && sim_muon_tree_->GetLeaf(yLeaf.c_str())){
            sim_muon_tree_->Draw((yLeaf + ":" + xLeaf).c_str(), "", "colz");
        }
        else if (sim_daughter_tree_->GetLeaf(xLeaf.c_str()) && sim_daughter_tree_->GetLeaf(yLeaf.c_str())){
            sim_daughter_tree_->Draw((yLeaf + ":" + xLeaf).c_str(), "", "colz");
        }

        c1->Draw();
        c1->Write();
        if (c1) delete c1;
    }

    file_->Write();

    // Here we can make beutified plots 
    getPlot2D(leafPairs[9]);
    // getPlot1D(leafPairs[0]);

    // 3 entries
    // getPlot2D3Entries("muon_polar_angle_", "muon_azimuthal_angle_", "muon_edeps_");

    if (file_) delete file_;
}

DECLARE_ANALYZER(CosmicsAnalysis);