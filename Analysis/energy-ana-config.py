from LDMX.Framework import ldmxcfg
p = ldmxcfg.Process('ana')
p.sequence = [ ldmxcfg.Analyzer.from_file('EnergyAnalysis.cxx') ]

import os,sys
cosmic_events_root = os.path.realpath(sys.argv[1]) # filepath to the root file

#output_string = "anaData/energy.root"
#p.outputFiles = [output_string]

p.inputFiles = [ cosmic_events_root ]
p.histogramFile = 'hist.root'
