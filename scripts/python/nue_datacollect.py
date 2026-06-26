#!/usr/bin/env python

# We are ignoring the possibility of muon productions

import argparse
import os
import numpy as np
import numpy.random as random
from icecube import icetray, dataio, simclasses, dataclasses 
import time
start = time.time()
parser = argparse.ArgumentParser(description = "Takes I3photons and I3Electrons from simulation files to read out positional data")
parser.add_argument("-i", "--infile", default = "/mnt/home/dillonb5/cascades/nue_data/gen_001")
parser.add_argument("-r", "--runnumber", type=int, default = 1)
parser.add_argument("-g", "--gcdfile", default = "/mnt/home/dillonb5/cascades/gcdfile/PONE_800mGrid.i3.gz")

args = parser.parse_args()
runnumber = -999
if args.runnumber < 10:
	runnumber = '00'+str(args.runnumber)
elif args.runnumber < 100:
	runnumber = '0'+str(args.runnumber)
else:
	runnumber = str(args.runnumber)
infile = args.infile + runnumber+".i3.zst"
gcdfile = args.gcdfile

dt = []
dphilst = []
dr = []
Eanglelst = []
eNRGlst=[]
lambdalst = []
domlst = []

tray = icetray.I3Tray()
tray.AddModule('I3Reader', 'reader', FilenameList = [gcdfile, infile])
## will need to be able to find: dr, dt, angle from cherenkov peak to angle photon is observed from
c = 299792458
n = 1.33
def drdtfinder(frame):
	if (len(frame['I3MCTree']) != 0) and (len(frame['I3Photons']) != 0):	
		Epos = frame['I3MCTree'][1].pos # Assuming LeptonInjector Nu_e simulation
		doms = frame['I3ModuleGeoMap']
		appendr = dr.append
		appendt = dt.append
#		print(frame['I3MCTree'][1])
		omkeys = frame['I3Photons'].keys()
		photons = frame['I3Photons']
		for key in omkeys:
			modulekey = dataclasses.ModuleKey(key.string, key.om)
			dompos = doms[modulekey].pos
			for photon in photons[key]:
				flight = dompos + photon.pos - Epos
				appendr(flight.magnitude)
				offset = flight.magnitude * n / c
				appendt(photon.time - offset*10**9)

def dphifinder(frame):
	if len(frame['I3Photons']) != 0:
#		omkeys = frame['I3Photons'].keys()
		Tree = frame['I3Photons']
		append = dphilst.append
		for photons in Tree.values():
			for photon in photons:
				phi = photon.dir
				randx, randy, randz = random.uniform(-1, 1, 3)
				x = phi.x
				y = phi.y
				z = phi.z
				dx = x - randx
				dy = y - randy
				dz = z - randz
				dphi = dataclasses.I3Direction(dx, dy, dz)
				append((dphi.zenith, dphi.azimuth))

def Eanglefinder(frame):
	if len(frame['I3Photons']) != 0:
#		omkeys = frame['I3Photons'].keys()
		Tree1 = frame['I3Photons']
		Tree2 = frame['I3MCTree']
		append = Eanglelst.append
		electronangle = Tree2[1].dir
		for photons in Tree1.values():
			for photon in photons:
				append((electronangle.zenith, electronangle.azimuth))

def metadata(frame):
	if len(frame['I3Photons']) != 0:
#		omkeys = frame['I3Photons'].keys()
		Tree1 = frame['I3Photons']
		Tree2 = frame['I3MCTree']
		NRG = Tree2[1].energy
		for photons in Tree1.values():
			for photon in photons:
				eNRGlst.append(NRG)
				lambdalst.append(photon.wavelength)	

tray.AddModule(drdtfinder, Streams = [icetray.I3Frame.DAQ])
tray.AddModule(dphifinder, Streams = [icetray.I3Frame.DAQ])
tray.AddModule(Eanglefinder, Streams = [icetray.I3Frame.DAQ])
tray.AddModule(metadata, Streams=[icetray.I3Frame.DAQ])

print('Tray Populated')
tray.Execute()
tray.Finish()
print('Tray Finished')
ary = np.column_stack((dt, dphilst, dr, Eanglelst, eNRGlst, lambdalst))
cols=['dt', 'dZenith', 'dAzimuth', 'dr', 'EZenith', 'EAzimuth', 'eNRG', 'Lambda']
np.savez('/mnt/home/dillonb5/cascades/outfiles/ary_'+runnumber+'.npz', data=ary, columns= cols)
end = time.time()
print('Elapsed time: ' + str(end - start)) 
