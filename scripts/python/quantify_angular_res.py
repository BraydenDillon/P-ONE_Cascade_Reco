import sys
import os

sys.path.insert(0, os.path.abspath('..'))
import photospline
from icecube import dataio, dataclasses, icetray, simclasses
import numpy as np
import scipy.optimize as optimize
import argparse
import random


c = 299792458
n = 1.33

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

angular_resolution = []

tray = icetray.I3Tray()
tray.AddModule('I3Reader', 'reader', FilenameList = [gcdfile, infile])

string_subset = np.array([266, 199, 220, 275,  96, 113, 112, 286, 173, 116,   8, 240, 130,
       307, 281, 112, 324, 306, 183, 303, 314, 289,  72, 113, 147, 177,
       160,  57,  54, 142, 291, 204,  75, 215, 179, 143, 315, 201, 182,
        78,  60,   1, 326,  46, 272, 232, 134, 162, 268, 101, 139, 320,
       195,  11,  86, 300,  84,  67, 129,  63, 310,  59, 215, 312, 188,
       138,  71, 120, 139,  12,  26, 117, 129,  40,   6, 156,  79, 132,
       127, 161])

c = 299792458
n = 1.34
def datacollect(frame):
    dt = []
    t = []
    dphilst = []
    dr = []
    xyz = []
    

    Epos = frame['I3MCTree'][1].pos
    doms = frame['I3ModuleGeoMap']
    omkeys = frame['I3Photons'].keys()
    photons = frame['I3Photons']
    for key in omkeys:
        if key.string not in string_subset:
            continue
        modulekey = dataclasses.ModuleKey(key.string, key.om)
        dompos = doms[modulekey].pos
        for photon in photons[key]:
            photon_pos = dompos + photon.pos
            xyz.append([photon_pos.x, photon_pos.y, photon_pos.z])
            flight = dompos + photon.pos - Epos
            dr.append(flight.magnitude)
            offset = flight.magnitude * n / c
            dt.append(photon.time - offset*10**9)
            t.append(photon.time)
            phi = photon.dir
            # randx, randy, randz = random.uniform(-1, 1, 3)
            # x = phi.x
            # y = phi.y
            # z = phi.z
            # dx = x - randx
            # dy = y - randy
            # dz = z - randz
            # dphi = dataclasses.I3Direction(dx, dy, dz)
            # dphilst.append(dphi.zenith)
            # Etheta.append(flight.azimuth)
            # Ephi.append(flight.zenith)
            
    return np.column_stack([xyz, t])
    
    

def displacement_magnitude(pos1: np.array, pos2: np.array) -> float:
    vector = pos1 - pos2
    return np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)

splinefit_3d = photospline.SplineTable('/mnt/home/dillonb5/cascades/fits/splinelog_3D.fits')


def Likelihood_3d(coords: np.array, Event: np.array):
    L = 0
    # coords should have shape [x,y,z,theta,phi,t]
    # Event has shape [N, 6] cols:(x,y,z,t,dr,dt). We only use the first 4 here
    
    event_xyz = Event[:,0:3]
    event_t = Event[:,3]
    
    
    
    # Calculate Displacement Magnitude
    diff = coords[0:3] - event_xyz
    dr = np.linalg.norm(diff, axis=1)
    # Calculate Time Residual
    dt = abs(coords[5] - event_t) - (1.34*dr/c * 1e9)

    # Construct Electron direction unit vector from zenith and azimuth
    Ex = np.sin(coords[4])*np.cos(coords[3])
    Ey = np.sin(coords[4])*np.sin(coords[3])
    Ez = np.cos(coords[4])


    # Calculate angle between electron travel vector and displacement vector
    Eangle = np.array([Ex, Ey, Ez])
    Ephi = np.arccos(np.dot(diff, Eangle) / dr)
    #print(Ephi)
    #Calculate Likelihood from constructed coordinates
    params = np.array([dr, Ephi, dt])
    vals = splinefit_3d.evaluate_simple([params[0], params[1], params[2]])
    L = np.where(vals == 0, -30, vals)
    return -np.sum(L)


# Take absolute positions of PMTs and time 
# Electron is moving which changes photon paths based on fixed photon hits
# takes an electron position guess, then moves the electron around based on output from likelihood to optimize
def minimizer(guess, event, function=Likelihood_3d):
    # Guess should be xyzt
    minimized = optimize.minimize(function, 
                                  x0=guess,
                                  args=(event,),  
                                  method='Nelder-Mead', 
                                  tol=1e-3)
    return minimized

def angular_res(frame):
    if (len(frame['I3MCTree']) != 0) and (len(frame['I3Photons']) != 0):
        EventData = datacollect(frame) 
        if EventData.shape[0] == 0:
            return  # skip this frame, no photons on selected strings
        truth = np.array([frame['I3MCTree'][1].pos.x, frame['I3MCTree'][1].pos.y, frame['I3MCTree'][1].pos.z, frame['I3MCTree'][1].dir.azimuth, frame['I3MCTree'][1].dir.zenith, frame['I3MCTree'][1].time])
        Energy = frame['I3MCTree'][0].energy
        best_fit = minimizer(truth, EventData)
        fit_zenith = best_fit['x'][4]
        fit_azimuth = best_fit['x'][3]
        truth_zenith=truth[4]
        truth_azimuth=truth[3]

        space_angle = 2*np.arcsin(np.sqrt(np.sin((fit_zenith - truth_zenith)/2)**2 + np.cos(truth_zenith)*np.cos(fit_zenith)*np.sin((fit_azimuth - truth_azimuth)/2)**2))
        
        # zenith_res = np.min(abs(fit_zenith - truth[4]), 360 - abs(fit_zenith - truth[4]))
        # azimuth_res = np.min(abs(fit_azimuth - truth[3]), 360 - abs(fit_azimuth - truth[3]))
        
        angular_resolution.append((space_angle, Energy))


tray.AddModule(angular_res, Streams = [icetray.I3Frame.DAQ])

print('Tray Populated')
tray.Execute()
tray.Finish()
print('Tray Finished')
np.save('/mnt/scratch/dillonb5/resolutions_80str/resolution_'+runnumber+'.npy', angular_resolution)

