import sys
import os

sys.path.insert(0, os.path.abspath(".."))
import photospline
from icecube import dataio, dataclasses, icetray, simclasses, phys_services
import numpy as np
from scipy import optimize, integrate
import argparse
from scripts.python.SplineEval import evalPdf
import random

c = 299792458
n = 1.35557

parser = argparse.ArgumentParser(
    description="Takes I3photons and I3Electrons from simulation files to read out positional data"
)
parser.add_argument(
    "-i", "--infile", default="/mnt/home/dillonb5/cascades/nue_data/gen_001"
)
parser.add_argument("-r", "--runnumber", type=int, default=1)
parser.add_argument(
    "-g", "--gcdfile", default="/mnt/home/dillonb5/cascades/gcdfile/PONE_800mGrid.i3.gz"
)
parser.add_argument("-d", "--Dimensionality", type=str, default="2d")
args = parser.parse_args()
runnumber = -999
if args.runnumber < 10:
    runnumber = "00" + str(args.runnumber)
elif args.runnumber < 100:
    runnumber = "0" + str(args.runnumber)
else:
    runnumber = str(args.runnumber)
infile = args.infile + runnumber + ".i3.zst"
gcdfile = args.gcdfile
Ndim = args.Dimensionality

delta_logL = []

tray = icetray.I3Tray()
tray.AddModule("I3Reader", "reader", FilenameList=[gcdfile, infile])


def datacollect(frame):
    t = []
    xyz = []

    doms = frame["I3ModuleGeoMap"]
    omkeys = frame["new_photons"].keys()
    photons = frame["new_photons"]
    for key in omkeys:
        modulekey = dataclasses.ModuleKey(key.string, key.om)
        dompos = doms[modulekey].pos
        for photon in photons[key]:
            photon_pos = dompos + photon.pos
            xyz.append([photon_pos.x, photon_pos.y, photon_pos.z])
            t.append(photon.time)

    return np.column_stack([xyz, t])


def displacement_magnitude(pos1: np.array, pos2: np.array) -> float:
    vector = pos1 - pos2
    return np.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)


splinefit = photospline.SplineTable("/mnt/home/dillonb5/cascades/fits/splinelog.fits")
splinefit_3d = photospline.SplineTable(
    "/mnt/home/dillonb5/cascades/fits/splinelog_3D.fits"
)


def Likelihood(coords: np.array, Event):
    L = 0
    # coords should have shape [x,y,z,t]
    # Event has shape [x,y,z,t,dr,dt,dphi,Ephi]. We only use the first 4 here
    event_xyz = Event[:, 0:3]
    event_t = Event[:, 3]

    diff = coords[0:3] - event_xyz
    dr = np.linalg.norm(diff, axis=1)
    # Calculate Time Residual
    # dt = abs(coords[-1] - event_t) - (n*dr/c * 1e9)
    # For sampled data, t is stored as time residual
    dt = event_t

    params = np.array([dr, dt])

    vals = splinefit.evaluate_simple([params[0], params[1]])
    L = np.where(vals == 0, -30, vals)
    return -np.sum(L)


def Likelihood_3d(coords: np.array, Event: np.array):
    L = 0
    # coords should have shape [x,y,z,theta,phi,t]
    # Event has shape [N, 6] cols:(x,y,z,t,dr,dt). We only use the first 4 here

    event_xyz = Event[:, 0:3]
    event_t = Event[:, 3]

    # Calculate Displacement Magnitude
    diff = -coords[0:3] + event_xyz
    dr = np.linalg.norm(diff, axis=1)
    # Calculate Time Residual
    dt = event_t - coords[5] - (1.34*dr/c * 1e9)
    # dt = event_t

    # Construct Electron direction unit vector from zenith and azimuth
    Ex = np.sin(coords[4]) * np.cos(coords[3])
    Ey = np.sin(coords[4]) * np.sin(coords[3])
    Ez = np.cos(coords[4])

    # Calculate angle between electron travel vector and displacement vector
    Eangle = np.array([Ex, Ey, Ez])
    Ephi = np.arccos(np.dot(diff, Eangle) / dr)

    # print(Ephi)
    # Calculate Likelihood from constructed coordinates
    params = np.array([dr, Ephi, dt])
    vals = splinefit_3d.evaluate_simple([params[0], params[1], params[2]])
    L = np.where(vals == 0, -30, vals)
    return -np.sum(L)


# Take absolute positions of PMTs and time
# Electron is moving which changes photon paths based on fixed photon hits
# takes an electron position guess, then moves the electron around based on output from likelihood to optimize
def minimizer(guess, event, function=Likelihood):
    # Guess should be xyzt
    minimized = optimize.minimize(
        function, x0=guess, args=(event), method="Nelder-Mead", tol=1e-5
    )
    return minimized


def evaluate_frame(frame):
    if (len(frame["I3MCTree"]) != 0) and (len(frame["new_photons"]) != 0):
        EventData = datacollect(frame)
        if Ndim == "2d":
            truth = np.array(
                [
                    frame["I3MCTree"][1].pos.x,
                    frame["I3MCTree"][1].pos.y,
                    frame["I3MCTree"][1].pos.z,
                    frame["I3MCTree"][1].time,
                ]
            )
            func = Likelihood
        elif Ndim == "3d":
            truth = np.array(
                [
                    frame["I3MCTree"][1].pos.x,
                    frame["I3MCTree"][1].pos.y,
                    frame["I3MCTree"][1].pos.z,
                    frame["I3MCTree"][1].dir.azimuth,
                    frame["I3MCTree"][1].dir.zenith,
                    frame["I3MCTree"][1].time,
                ]
            )
            func = Likelihood_3d

        else:
            print("Error: Ndim must be either '3d' or '2d'")
            return True
        model = minimizer(truth, EventData, func)
        model_likelihood = frame['LLHFitFitParams']
        truth_likelihood = func(truth, EventData)
        delta_logL.append(
            [truth_likelihood - model_likelihood, len(frame["new_photons"])]
        )

    # return truth_likelihood - model_likelihood


tray.AddModule(evaluate_frame, Streams=[icetray.I3Frame.DAQ])


print("Tray Populated")
tray.Execute()
tray.Finish()
print("Tray Finished")
ary = np.array(delta_logL)
np.save("/mnt/scratch/dillonb5/cdf_per_photon_logL/delta_ary_" + runnumber + ".npy", ary)
