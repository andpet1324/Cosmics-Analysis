# What the script does

This script analyzes cosmic muons in the back-HCal. An algorithm for the side-HCal (and ECal) is not implemented.

The whole analysis uses truth-level data, but should be implemented for reconstruction.

A least squares method fits a line to the track of cosmic muons. A quadratic fit method was started but never finished (it finds a fit but it is not aligned to the data).

Three selections exist for MIM-like (minimally ionizing muons) events. These are angle_threshold_, R^2 and MinFCN.

We can build a map that stores information for all scintillator bars in the back-HCal using the buildMap function.

With graphcheck enabled, the script plots the path of the muons in 3D plots.
With estimatedcheck enabled, the script calculates the estimated energy in each bar.

We can analyse horisontal and vertical bars individually by changing chose_orientation_.

# How to run

Fork this repository and cd into it.

## Install denv and change image
install denv as according to https://ldmx-software.github.io/

I used ldmx/dev:4.2.2. 

To change image:

denv config image ldmx/dev:4.2.2
denv config image pull

## Generate an LHE-file

The LHE-file contains initial information about the muon.

run:
denv python3 lheData/cosmic_muon_lhe_generator_v14.py --numEvents=xx --detector=backHcal

more options exist for the generator. See cosmic_muon_lhe_generator_v14.py.
An LHE-file is now created in lheData/

## Simulate in Geant4

cosmics_config.py contains an option for max allowed events. Change this to the number of events you want to run.

run:
denv fire cosmics_config.py lheData/name_of_file.lhe

This creates a ROOT-file in cosmicEvents/

## Analyse

config file uses CosmicsAnalysis.cxx & CosmicsAnalysis.h for analysis.

run:
denv fire Analysis/cosmics_ana_config.py cosmicEvents/name_of_file.root

Graphs and histograms are saved in Analysis/



