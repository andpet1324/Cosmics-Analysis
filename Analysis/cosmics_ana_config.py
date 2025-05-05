
from LDMX.Framework import ldmxcfg

p = ldmxcfg.Process('ana')
p.sequence = [ ldmxcfg.Analyzer.from_file('CosmicsAnalysis.cxx') ]

import LDMX.Hcal.HcalGeometry
import LDMX.DetDescr.HcalGeometry
import LDMX.Hcal.hcal_hardcoded_conditions

import os,sys
cosmic_events_root = os.path.realpath(sys.argv[1]) # filepath to the root file

p.inputFiles = [ cosmic_events_root ]
p.histogramFile = 'Analysis/cosmics_hist.root'
