import argparse
import numpy.random as random
from icecube import icetray, dataio, simclasses, dataclasses

parser = argparse.ArgumentParser()

parser.add_argument("-g", "--gcdfile", default = "/mnt/home/dillonb5/cascades/gcdfile/PONE_800mGrid.i3.gz")
args = parser.parse_args()
gcdfile = args.gcdfile

tray = icetray.I3Tray()
tray.AddModule('I3Reader', 'reader', FilenameList = [gcdfile])
modules = []
def gcdprinter(frame):
	for key in frame['I3ModuleGeoMap'].keys():
		modules.append(key)

tray.AddModule('gcdprinter', Streams = [icetray.I3Frame.DAQ])
tray.Execute()
tray.Finish()

print(modules)	
