#!/usr/bin/env python

"""
file: cosmic_muon_lhe_generator.py
author: Reese Petersen, Tom Eichlersmith, Martin Meier
affiliations: University of Minnesota, LDMX Collaboration
email: pet00831@umn.edu, eichl008@umn.edu

date: 16 July 2018
updated: January 2021

This script writes lhe files which descibe cosmic muons 
with energy and momentum distributed according to a function provided by
A. Dar in "Atmospheric Neutrinos, Astrophysical Neutrons, and Proton-Decay Experiments",
18 July, 1983, Phys. Rev. Lett. Vol 51, No. 3, pg 227
Other papers necessary to understand this function in its entirety include reference 7 from A. Dar:
"The Muon Charge Ratio at Sea-Level" by A. Liland, 
16th International Cosmic Ray Conference, 
Conference Papers vol. 13 (Late Papers), Kyoto, Japan, August 1979

Examples
--------
    ./cosmic_muon_lhe_generator.py --numEvents 1000 
"""

import math
import random
import ROOT #2D functional spectrums
import ctypes #pass by reference doubles

class lhe_file :
    """Helper class for writing a file that is supposed to be in the LHE format.

    Parameters
    ----------
    name : str
        Name of File to write out

    Attributes
    ----------
    f : file
        File we are writing out
    """

    def __init__(self,name) :
        self.f = open(name,'w')

    def __del__(self) :
        self.f.close()

    def write_event(self,vertex,four_momentum,pdg_id,mass,muondelay) :
        """Write an Event to the LHE file
        
        Parameters
        ----------
        vertex : list of float
            space 3-point that the particle is originating from
        four_momentum : list of float
            four momentum of particle [GeV]
        pdg_id : int
            PDG ID of particle
        mass : float
            Mass of particle [GeV]
        muondelay : float
            Time shift of particle generation
        """

        self.f.write("<event>\n")
        self.f.write(" 1   0  0.3584989E-03  0.9118800E+02  0.7818608E-02  0.1180000E+00\n") #unused
        self.f.write("#vertex/time  %f  %f  %f  %f\n"%(vertex[0],vertex[1],vertex[2],muondelay))
        self.f.write("%9d    1    0    0    0    0  %.11E  %.11E  %.11E  %.11E  %.11E  0  1.0 \n"%(pdg_id,four_momentum[1],four_momentum[2],four_momentum[3],four_momentum[0],mass))
        self.f.write("</event>\n")

class cosmic_muon_spectrum :
    """Storage and use of the spectrum of cosmic muons

    ADar is the muon differential flux as a function of energy (in GeV) and zenith angle (in degrees) given by A. Dar (Eq. 6) in 
    "Atmospheric Neutrinos, Astrophysical Neutrons, and Proton-Decay Experiments", 18 July, 1983, Phys. Rev. Lett. Vol 51, No. 3, pg 227

    The muon spectrum has two parts, one for positive muons, and one for negative muons.
    The original function has units of muons/GeV/sr/s/cm^2, ADarplus and ADarminus are still expressed in these units, ADar has been normalized to 1 over the default range.

    This spectrum is fit to data taken at DESY, Hamburg, Germany (53.58 deg N, 9.88 deg E) at 0 degrees and 75 degrees in zenith angle at approximately sea-level altitude. (40 m)

    This spectrum matches data in the range: E = 0.3-1000 GeV, th = 0-75 GeV

    Effects from the environment of the detector, such as the building it will be housed in, have not been considered.

    The East-West charge ratio has not been accounted for. There is no azimuthal angular dependence. In other words, Earth's magnetic field has not been accounted for.
    The East-West charge ratio would have very little effect for muons with energy > 10 GeV and zenith angle < 75 degrees.

    ADar is used to simulate the energy and zenith angle of an incoming cosmic muon, ADar = ADarplus + ADarminus normalized to 1 over the default range (1-100 GeV, 0-75 degrees)
    ADarplus and ADarminus, which are the separate spectra for positive and negative muons, are used to determine the charge of each cosmic muon after the energy and angle are generated.

    Parameters
    ----------
    emin : float
        minimum energy for spectrum function [GeV]
    emax : float
        maximum energy for spectrum function [GeV]
    thmin : float
        minimum azimuthal angle for spectrum function [degrees]
    thmax : float
        maximum azimuthal angle for spectrum function [degrees]

    Attributes
    ----------
    cmmc_spec : ROOT.TF2
        2D Cosmic Muon Spectrum (energy vs theta)
    cmmc_spec_plus : ROOT.TF2
        2D Muon Plus Spectrum (energy vs theta)
    cmmc_spec_minus : ROOT.TF2
        2D Muon Minus Spectrum (energy vs theta)
    cmmc_spec_low : ROOT.TF2
        2D Low (<10GeV) Energy Muon Spectrum (energy vs theta)
    cmmc_spec_high : ROOT.TF2
        2D High (>10GeV) Energy Muon Spectrum (energy vs theta)
    low_frac : float
        Fraction of spectrum with E < 10 GeV
    """

    def __init__(self,emin,emax,thmin,thmax) :

        ADar = "(0.3403105255419549*TMath::Exp(2.172701227027019/(-1.0609+0.1236*TMath::Cos(TMath::Pi()/180*y)-x*TMath::Cos(TMath::Pi()/180*y)))*TMath::Power(TMath::Cos(TMath::Pi()/180*y),2.67-1.0106422229299041/(-1.0609+0.1236*TMath::Cos(TMath::Pi()/180*y)-x*TMath::Cos(TMath::Pi()/180*y)))*(-0.09553777560317707*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+6.68495653879148e-9*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3.5)+TMath::Sqrt(2.1218-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y))-2.3800027794354056e-8*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3)*TMath::Sqrt(2.1218-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y))+TMath::Power(x,2)*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3)*(3.2818816073967673e-7*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))-1.1684275435535303e-6*TMath::Sqrt(2.1218-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)))+x*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3)*(-8.112811333484808e-8*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+2.888352887664327e-7*TMath::Sqrt(2.1218-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)))+TMath::Power(x,3)*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3)*(-4.425406698215705e-7*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+1.5755495463235305e-6*TMath::Sqrt(2.1218-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)))+TMath::Power(TMath::Cos(TMath::Pi()/180*y),2)*(-6.315346819441329e-6*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+0.00005334447875055668*TMath::Sqrt(2.1218-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y))+x*(0.00005109503899224376*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))-0.0004315896339041803*TMath::Sqrt(2.1218-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)))+TMath::Power(x,2)*(-0.0001033475707771921*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+0.0008729563792560281*TMath::Sqrt(2.1218-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y))))+TMath::Cos(TMath::Pi()/180*y)*(0.0014548779564248476*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))-0.014359425845202903*TMath::Sqrt(2.1218-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y))+x*(-0.005885428626314109*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+0.05808829225405703*TMath::Sqrt(2.1218-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y))))))/((1-0.00040406601131876197*TMath::Cos(TMath::Pi()/180*y)+0.0016345712431988753*x*TMath::Cos(TMath::Pi()/180*y))*(1-0.0029397862016351797*TMath::Cos(TMath::Pi()/180*y)+0.011892339003378558*x*TMath::Cos(TMath::Pi()/180*y))*(1-0.005035070902595827*TMath::Cos(TMath::Pi()/180*y)+0.02036840980014493*x*TMath::Cos(TMath::Pi()/180*y))*(1-0.008226881320143641*TMath::Cos(TMath::Pi()/180*y)+0.03328026424006327*x*TMath::Cos(TMath::Pi()/180*y))*TMath::Power(2.1218-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y),2.67)*(-0.42008790004955787*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+TMath::Sqrt(2.1218-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y))))"
        ADarplus = "(0.04454993403165476*TMath::Exp(2.172701227027019/(0.1236*TMath::Cos(TMath::Pi()/180*y)-x*TMath::Cos(TMath::Pi()/180*y)-1.0609))*TMath::Power(TMath::Cos(TMath::Pi()/180*y),2.67-1.0106422229299041/(0.1236*TMath::Cos(TMath::Pi()/180*y)-x*TMath::Cos(TMath::Pi()/180*y)-1.0609))*(3.831797941255795e-9*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3.5)-2.5067089548238198e-8*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3)*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218)+TMath::Power(x,2)*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3)*(1.8811651375285701e-7*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))-1.230632086565563e-6*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218))+x*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3)*(-4.6502402199706266e-8*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+3.042122517990072e-7*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218))+TMath::Power(x,3)*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3)*(-2.5366304443481265e-7*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+1.6594283799427765e-6*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218))+(0.09269410920144934*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218))+TMath::Power(TMath::Cos(TMath::Pi()/180*y),2)*(3.2762610110392485e-6*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+0.00005369460580045414*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218)+x*(-0.000026506966108731775*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))-0.0004344223770263279*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218))+TMath::Power(x,2)*(0.00005361441365034748*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+0.0008786860376746114*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218)))+TMath::Cos(TMath::Pi()/180*y)*(-0.0012174109045183516*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))-0.014383269769682962*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218)+x*(0.004924801393682649*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+0.05818474825923528*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218)))))/((-0.00040406601131876197*TMath::Cos(TMath::Pi()/180*y)+0.0016345712431988753*x*TMath::Cos(TMath::Pi()/180*y)+1)*(-0.0029397862016351797*TMath::Cos(TMath::Pi()/180*y)+0.011892339003378558*x*TMath::Cos(TMath::Pi()/180*y)+1)*(-0.005035070902595827*TMath::Cos(TMath::Pi()/180*y)+0.02036840980014493*x*TMath::Cos(TMath::Pi()/180*y)+1)*(-0.008226881320143641*TMath::Cos(TMath::Pi()/180*y)+0.03328026424006327*x*TMath::Cos(TMath::Pi()/180*y)+1)*TMath::Power(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218,2.67)*(-0.42008790004955787*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218)))"
        ADarminus = "(0.031055887893745797*TMath::Exp(2.172701227027019/(0.1236*TMath::Cos(TMath::Pi()/180*y)-x*TMath::Cos(TMath::Pi()/180*y)-1.0609))*TMath::Power(TMath::Cos(TMath::Pi()/180*y),2.67-1.0106422229299041/(0.1236*TMath::Cos(TMath::Pi()/180*y)-x*TMath::Cos(TMath::Pi()/180*y)-1.0609))*(1.0777836695267741e-8*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3.5)-2.1982416984362076e-8*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3)*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218)+TMath::Power(x,2)*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3)*(5.29122123868283e-7*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))-1.0791946001214612e-6*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218))+x*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3)*(-1.3079898902023958e-7*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+2.667769051500252e-7*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218))+TMath::Power(x,3)*TMath::Power(TMath::Cos(TMath::Pi()/180*y),3)*(-7.134872220446103e-7*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+1.455224649570471e-6*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218))+(-0.3655580074957874*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218))+TMath::Power(TMath::Cos(TMath::Pi()/180*y),2)*(-0.0000200745894328675*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+0.00005284221853410446*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218)+x*(0.0001624157721105785*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))-0.000427526039919939*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218))+TMath::Power(x,2)*(-0.0003285108659194549*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+0.0008647371357603943*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218)))+TMath::Cos(TMath::Pi()/180*y)*(0.005288298944929909*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))-0.014325221533728156*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218)+x*(-0.021392795084667918*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+0.057949925298253054*TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218)))))/((-0.00040406601131876197*TMath::Cos(TMath::Pi()/180*y)+0.0016345712431988753*x*TMath::Cos(TMath::Pi()/180*y)+1)*(-0.0029397862016351797*TMath::Cos(TMath::Pi()/180*y)+0.011892339003378558*x*TMath::Cos(TMath::Pi()/180*y)+1)*(-0.005035070902595827*TMath::Cos(TMath::Pi()/180*y)+0.02036840980014493*x*TMath::Cos(TMath::Pi()/180*y)+1)*(-0.008226881320143641*TMath::Cos(TMath::Pi()/180*y)+0.03328026424006327*x*TMath::Cos(TMath::Pi()/180*y)+1)*TMath::Power(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218,2.67)*(-0.42008790004955787*TMath::Sqrt(TMath::Cos(TMath::Pi()/180*y))+TMath::Sqrt(-0.2472*TMath::Cos(TMath::Pi()/180*y)+x*TMath::Cos(TMath::Pi()/180*y)+2.1218)))"
        
        self.cmmc_spec = ROOT.TF2("cmmc_spec",ADar,emin,emax,thmin,thmax)
        self.cmmc_spec_plus = ROOT.TF2("cmmc_spec_plus",ADarplus,emin,emax,thmin,thmax)
        self.cmmc_spec_minus = ROOT.TF2("cmmc_spec_minus",ADarminus,emin,emax,thmin,thmax)

        # The spectrum is simulated in two parts, in the range E = 1(or emin) to 10 GeV and in the range E = 10 to 100(or emax) GeV. 
        # This allows better resolution of the steep peak at low energy and low zenith angle. Most muons have E < 10 GeV
        # full spectrum split along 10 GeV:
        # set the resolution over each energy range and over theta to 250 points.
        self.cmmc_spec_low = ROOT.TF2("cmmc_spec1",ADar,emin,10.0,thmin,thmax)
        self.cmmc_spec_low.SetNpy(250)
        self.cmmc_spec_low.SetNpx(250)

        self.cmmc_spec_high = ROOT.TF2("cmmc_spec2",ADar,10.0,emax,thmin,thmax)
        self.cmmc_spec_high.SetNpy(250)
        self.cmmc_spec_high.SetNpx(250)
        
        # Find the fraction of events in the low E spectrum
        full_int = self.cmmc_spec.Integral(emin,emax,thmin,thmax)
        self.low_frac = self.cmmc_spec_low.Integral(emin,10,thmin,thmax)/full_int

    def rand_energy_theta(self) :
        """Pull a random energy [GeV] and azimuthal angle [degrees] from the spectrum
        
        Returns
        -------
        two floats
            Energy [GeV] and Polar Angle [degrees] randomly sampled from spectrum
        """

        E = ctypes.c_double()
        thdeg = ctypes.c_double()

        rand = random.random()
        if rand <= self.low_frac: # E <= 10 GeV
            self.cmmc_spec_low.GetRandom2(E,thdeg)
        else: # E > 10 GeV
            self.cmmc_spec_high.GetRandom2(E,thdeg)

        return E.value,thdeg.value

    def rand_charge(self,energy_GeV,theta_deg) :
        """Pull (randomly) the charge of the muon using the already generated energy and angle

        Returns
        -------
        int
            PDG ID of Muon randomly pulled from spectrums
        """

        # determine the charge of the muon
        flux_plus = self.cmmc_spec_plus.Eval(energy_GeV,theta_deg)
        flux_minus = self.cmmc_spec_minus.Eval(energy_GeV,theta_deg)
        prob_plus = flux_plus/(flux_plus+flux_minus)
        rand_pm = random.random()
        if rand_pm < prob_plus:
          pdgID = -13
        else:
          pdgID = 13

        return pdgID

def main( numEvents , numFiles , outDir , energyMin , energyMax , thetaMin , thetaMax , hcalDepth, timeWindow, timeShift, detector ) :
    """Main generation of LHE files done here

    We use the helper classes defined above to generate LHE files containing single cosmic
    muons per event. Here we also generate a vertex for originating the cosmic and check
    that the generated cosmic will intersect some portion of the HCal volume.

    Parameters
    ----------
    numEvents : int
        Number of events to put into a single file (one muon per event)
    numFiles : int
        Number of files to generate
    outDir : str
        Directory to put generated files into
    energyMin : float
        Minimum energy of cosmic muons to generate [GeV]
    energyMax : float
        Maximum energy of cosmic muons to generate [GeV]
    thetaMin : float
        Minimum polar angle of cosmics (relative to vertical) [degrees]
    thetaMax : float
        Maximum polar angle of cosmics (relative to vertical) [degrees]
    hcalDepth : float
        Thickness of HCal along z direction [mm]
    """

    ROOT.gRandom.SetSeed(0) #seeds random number generator based off time
    # muon mass in GeV
    muon_mass = 0.1056583715


    ### detector parameters
    # detector regions are determined from the GDML files in the ldmx-det-v12 data directory

    #full detector region
    if detector == 'fullDetector': 
        xmin = -1550.
        xmax = 1550.
        ymin = -1550.
        ymax = 1550.
        zmin = -860.
        zmax = 220. + arg.hcalDepth

    #full Hcal
    #the front of the side HCal and ECal is located 220mm downstream from z=0
    if detector == 'fullHcal': 
        xmin = -1550.
        xmax = 1550.
        ymin = -1550.
        ymax = 1550.
        zmin = 220.
        zmax = zmin + arg.hcalDepth

    #side Hcal
    elif detector == 'sideHcal': 
        xmin = -1550.
        xmax = 1550.
        ymin = -1550.
        ymax = 1550.
        zmin = 220.
        zmax = 820.

    #back Hcal
    elif detector == 'backHcal': 
        xmin = -1550.
        xmax = 1550.
        ymin = -1550.
        ymax = 1550.
        zmin = 820.
        zmax = zmin + arg.hcalDepth - 600.   #subtract off side HCal thickness from zmin
        
    #Ecal
    elif detector == 'ecal': 
        xmin = -260.
        xmax = 260.
        ymin = -260.
        ymax = 260.
        zmin = 220.
        zmax = 666.

    #target area
    elif detector == 'target': 
        xmin = -20.
        xmax = 20.
        ymin = -50.
        ymax = 50.
        zmin = -4.677
        zmax = 4.677

    #tagging tracker
    elif detector == 'taggingTracker': 
        xmin = -212.5
        xmax = 212.5
        ymin = -172.5
        ymax = 172.5
        zmin = -860.
        zmax = -4.677

    #recoil tracker
    elif detector == 'recoilTracker': 
        xmin = -212.5
        xmax = 212.5
        ymin = -172.5
        ymax = 172.5
        zmin = 4.677
        zmax = 204.677

    #offset detector volume by 0.001mm clearance so that particle are not generate at the interface of two materials (can cause issues with Geant4)
    clearance=0.001
    Xmin = xmin-clearance
    Xmax = xmax+clearance
    Ymin = ymin-clearance
    Ymax = ymax+clearance
    Zmin = zmin-clearance
    Zmax = zmax+clearance

    the_cosmic_muon_spectrum = cosmic_muon_spectrum( energyMin, energyMax, thetaMin, thetaMax )

    for n in range(numFiles):
        filename = outDir+"/"+"cmmc_energy_"+str(int(energyMin))+"_"+str(int(energyMax))+"_GeV_theta_"+str(int(thetaMin))+"_"+str(int(thetaMax))+"_deg_" +str(numEvents)+"_events_%04d.lhe"%(n)
        new_lhe = lhe_file(filename)
        i=0
        while i < numEvents:
            # determine the energy and angle of the muon
            energy, theta_deg = the_cosmic_muon_spectrum.rand_energy_theta()
            # get charge
            pdgID = the_cosmic_muon_spectrum.rand_charge(energy,theta_deg)
            # get azimuthal angle
            phi = random.uniform(0,2*math.pi)
        
            # do math to get the kinematics into 3-momentum form
            p = math.sqrt(math.pow(energy,2)-math.pow(muon_mass,2))
            theta_rad = theta_deg*math.pi/180
            # get the muon momentum in detector cooridinates (y-axis up and z-axis down-beam)
            px = -p*math.sin(theta_rad)*math.sin(phi)
            py = -p*math.cos(theta_rad)
            pz = -p*math.sin(theta_rad)*math.cos(phi)
        
            # pick a random point inside the full hcal volume in detector coordinates, (0,0,0) is the center of the target
            z = random.uniform(zmin,zmax)
            if detector != 'fullDetector':              #anything besides the full detector
                x = random.uniform(xmin,xmax)
                y = random.uniform(ymin,ymax)
            else:
                if z < 220:                             #this z region encompasses the tagger/target region
                    x = random.uniform(-212.5,212.5)
                    y = random.uniform(-172.5,172.5)
                else:
                    x = random.uniform(xmin,xmax)
                    y = random.uniform(ymin,ymax)


            # find the intersection with the detector volume
            # xt, yt, and zt are where the incoming muon would intersect infinite planes containing the planes of the hcal volume
            # check for muons intersecting the top plane first
            xt = x + px/py*(Ymax-y)
            zt = z + pz/py*(Ymax-y)
            if Xmin < xt < Xmax and Zmin < zt < Zmax:
              xv = xt
              yv = Ymax
              zv = zt


            # now check for muons intersecting the positive x, negative x, positive z, and negative z planes
            else:
              if px < 0:
                yt = y + py/px*(Xmax-x)
                zt = z + pz/px*(Xmax-x)
                if Ymin < yt < Ymax and Zmin < zt < Zmax:
                  xv = Xmax
                  yv = yt
                  zv = zt

              elif px > 0:
                yt = y + py/px*(Xmin-x)
                zt = z + pz/px*(Xmin-x)
                if Ymin < yt < Ymax and Zmin < zt < Zmax:
                  xv = Xmin
                  yv = yt
                  zv = zt

              elif pz < 0:
                xt = x + px/pz*(Zmax-z)
                yt = y + py/pz*(Zmax-z)
                if Xmin < xt < Xmax and Ymin < yt < Ymax:
                  xv = xt
                  yv = yt
                  zv = Zmax

              elif pz > 0:
                xt = x + px/pz*(Zmin-z)
                yt = y + py/pz*(Zmin-z)
                if Xmin < xt < Xmax and Ymin < yt < Ymax:
                  xv = xt
                  yv = yt
                  zv = Zmin
                else:
                  continue
            

            # find muon TOF from [xv,xt,xz] to [x,y,z]
            dist=math.sqrt(math.pow(x - xv,2) + math.pow(y - yv,2) + math.pow(z - zv,2)) #distance from point of muon generation to random location in the detector volume (mm) 
            v=3*10**11 #(c in mm/s)
            muon_time=(dist/v)*10**9 #(nanoseconds)

            # find beam TOF from z=0 to z (randomly sampled Hcal depth from before)
            beamTOF=(z/v)*10**9 #(nanoseconds)

            muon_delay=(beamTOF-muon_time)+(random.uniform(-0.5,0.5)*timeWindow)+timeShift #time offset of muon generation

            new_lhe.write_event( [xv,yv,zv] , [energy,px,py,pz] , pdgID , muon_mass, muon_delay)

            i+=1
        #end loop over events
    #end loop over files
#end main def

if __name__ == '__main__' :
    """We are main, let's do this!"""

    import os
    import argparse

    # get number of events
    parser = argparse.ArgumentParser(
        description="Makes lhe files for atmospheric muons intersecting the HCal of the LDMX detector.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

    parser.add_argument("--numEvents" , dest = "numEvents" , help = "number of events per lhe file"                     , default = 10 , type = int)
    parser.add_argument("--numFiles"  , dest = "numFiles"  , help = "number of lhe files to make"                       , default = 1  , type = int)
    parser.add_argument("--outputDir" , dest = "outputDir" , help = "directory to put the output lhe files"             , default = os.getcwd())
    parser.add_argument("--energyMin" , dest = "energyMin" , help = "low energy limit for the muon spectrum in GeV"     , default = 1.    , type=float)
    parser.add_argument("--energyMax" , dest = "energyMax" , help = "high energy limit for the muon spectrum in GeV"    , default = 100.  , type=float)
    parser.add_argument("--thetaMin"  , dest = "thetaMin"  , help = "low theta limit for the muon specturm in degrees (theta is measured with respect to the vertical, or +y-axis in detector coordinates)"  , default = 0.    , type=float)
    parser.add_argument("--thetaMax"  , dest = "thetaMax"  , help = "high theta limit for the muon spectrum in degrees" , default = 75.   , type=float)
    parser.add_argument("--hcalDepth" , dest = "hcalDepth" , help = "thickness of the full hcal (side and back) in z"   , default = 5500. , type=float)
    parser.add_argument("--timeWindow", dest = "timeWindow", help = "window (half-width) of time around coinciding beam and muon events (ns)" , default = 15. , type=float)
    parser.add_argument("--timeShift" , dest = "timeShift" , help = "positive or negative time shift (ns) of the center of the time window to make the muon arrive earlier or later than the beam particle" , default = 0. , type=float)
    parser.add_argument("--detector"  , dest = "detector"  , help = "which detector (fullDetector, fullHcal, sideHcal, backHcal, ecal, target, taggingTracker, recoilTracker)"    , default = 'fullDetector', type=str)

    arg = parser.parse_args()
    
    # make output directory if it doesn't exist
    full_path_out_dir = os.path.realpath(arg.outputDir)
    if not os.path.exists(full_path_out_dir) :
        os.makedirs(full_path_out_dir)

    main( arg.numEvents , arg.numFiles , full_path_out_dir
            , arg.energyMin , arg.energyMax 
            , arg.thetaMin , arg.thetaMax 
            , arg.hcalDepth, arg.timeWindow
            , arg.timeShift, arg.detector )
