
from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('cosmics')

import LDMX.Ecal.ecal_hardcoded_conditions 
import LDMX.Hcal.HcalGeometry 
import LDMX.Ecal.EcalGeometry 
import LDMX.Hcal.hcal_hardcoded_conditions

# import processor templates
#import LDMX.Ecal.digi as ecal_digi
import LDMX.Hcal.digi as hcal_digi
import LDMX.Hcal.HcalGeometry

from datetime import datetime
current_date = datetime.now().strftime("%Y%m%d")  # Added current date to the filename
output_string = "cosmicEvents/cosmic_events_" + current_date + ".root"

import time

start = time.time()

p.outputFiles = [output_string]
p.run = 1
p.logFrequency = 1
p.termLogLevel = 1

from LDMX.SimCore import simulator

sim = simulator.simulator('lhe_cosmics_simulation')
sim.setDetector( 'ldmx-det-v14' )


import os,sys
full_lhe_file_path = os.path.realpath(sys.argv[1])
print(os.listdir())
from LDMX.SimCore import generators
from LDMX.Hcal.hcal import HcalOldDigiProducer

sim.generators = [ generators.lhe( 'cosmic_muons' , full_lhe_file_path ) ]

p.maxEvents = 10000

p.sequence = [ sim , hcal_digi.HcalDigiProducer(), hcal_digi.HcalRecProducer()]
# p.sequence = [ sim , hcal_digi.HcalDigiProducer(), hcal_digi.HcalSimpleDigiAndRecProducer()]
# p.sequence.append( HcalOldDigiProducer() )