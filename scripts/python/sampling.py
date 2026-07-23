from icecube import dataclasses, phys_services, icetray, simclasses, dataio
import numpy as np
import sys
import os

sys.path.insert(0, os.path.abspath("/mnt/home/dillonb5/cascades/"))
from scripts.python.SplineEval import evalPdf
import photospline
import argparse
from scipy.stats import norm

# Re-written functions retrofitted for 4d functionality ########################################################################################################
LOG_FLOOR = -30.0
def evalLogPdf(spline, dr, phiE, dphi, dt):
    """PDF (density) at fixed scalar (dr, phiE) over a dt array-like.

    Clamps dr/phiE into extent, floors out-of-dt-support to LOG_FLOOR, then exp()."""
    # unpack all 3 ranges
    (drLo, drHi), (peLo, peHi), (dpLo, dpHi), (dtLo, dtHi) = spline.extents

    # clamp dr and phiE and dphi
    drC = float(np.clip(dr, drLo, drHi))
    peC = float(np.clip(phiE, peLo, peHi))
    dpC = float(np.clip(dphi, dpLo, dpHi))

    # normalize dt into a 1D float array.
    # Just as a safeguard because evaluate_simple takes an array (I think)
    dt = np.atleast_1d(np.asarray(dt, dtype=float))

    # clamp dt
    dtC = np.clip(dt, dtLo, dtHi)  # keep evaluate_simple supported

    # query the spline at the clamped coordinate
    # .squeeze drops extra dims
    logp = spline.evaluate_simple([[drC], [peC], [dpC], dtC]).squeeze()

    # for each dt point, ask if it was originally inside the fitted window
    # if yes, keep the spline's value
    # if no, floor it to -30
    logp = np.where((dt >= dtLo) & (dt <= dtHi), logp, LOG_FLOOR)

    # return linear density (convert from log space)
    return np.exp(logp)


def evalPdf(spline, dr, phiE, dphi, dt):
    """Unit-normalised PDF: evalLogPdf divided by its own dt-integral.

    dt must span the spline's dt support for the normalisation to be meaningful
    (e.g. np.linspace(dtLo, dtHi, N))."""
    dt = np.atleast_1d(np.asarray(dt, dtype=float))
    pdf = evalLogPdf(spline, dr, phiE, dphi, dt)
    area = np.trapz(pdf, dt)
    return pdf / area if area > 0 else pdf

#########################################################################################################

## Set up argument parsing 
parser = argparse.ArgumentParser(description="Takes a .i3 file and a gcd file")
parser.add_argument(
    "-i", "--infile", default="/mnt/home/dillonb5/cascades/sacrifice_data/gen_"
) # format should be the name of the file without the run number or .i3.zst extension, e.g. gen_ for gen_001.i3.zst
parser.add_argument(
    "-g", "--gcdfile", default="/mnt/home/dillonb5/cascades/gcdfile/PONE_800mGrid.i3.gz"
)
parser.add_argument("-r", "--runnumber", type=int, default=50) 
parser.add_argument(
    "-o", "--outfile", default="/mnt/scratch/dillonb5/sampled_4d_7-23/new_"
) 
args = parser.parse_args()
runnumber = -999
if args.runnumber < 10:
    runnumber = "000" + str(args.runnumber)
elif args.runnumber < 100:
    runnumber = "00" + str(args.runnumber)
else:
    runnumber = str(args.runnumber)
infile = args.infile + runnumber + ".i3.zst" # This is where input file is completed to readable format
gcdfile = args.gcdfile
outfile = args.outfile + runnumber + ".i3.zst" # Outfile is treated the same to copy formatting so that it can be tracked back to the infile


tray = icetray.I3Tray()
tray.AddModule("I3Reader", "reader", FilenameList=[gcdfile, infile]) # Reads through all frames in FilenameList,
                                                                     # All Tray Modules will be executed on each frame read by I3Reader


# spline = photospline.SplineTable("/mnt/home/dillonb5/cascades/fits/splinelog_3D.fits") # Load in spline fit
spline = photospline.SplineTable("/mnt/research/IceCube/jalabadz/iter_6.0_4D_I3Photons/spline_result/splinelog.fits") # Load in spline fit
tgrid = np.linspace(spline.extents[2][0], spline.extents[2][1], 10000) # Define tgrid for interpolation of spline fit for sampling
# tgrid = np.linspace(-5, 5, 1000)
N_GROUP = 1.34
N_PHASE = 1.35557

rng = np.random.default_rng() # Initialize random number generator for sampling from spline fit


def resample(frame): # Sampling function. For each frame samples random value from cdf of spline fit and adds new object with new time
    stats = dict(events=0, doms=0, pulses=0, doms_clamped=0, doms_skipped=0, bad_doms=0)
    if "I3Photons" not in frame or "I3MCTree" not in frame:
        return True # Returns True and does not execute if there are not photons or no particles
    electron = dataclasses.I3Particle(frame["I3MCTree"][0]) # should change to neutrino if sampling in 4d
    electron.shape = dataclasses.I3Particle.Cascade # Not sure if needed, but explicitely sets the shape of the particle
    geo = frame["I3ModuleGeoMap"] # Grabs geo module map from gcd file
    old = frame["I3Photons"] # Old is the original I3Photons object, which is iterated over and resampled from cdf
    new = simclasses.I3CompressedPhotonSeriesMap() # Initializes new compressed photon series map
    stats["events"] += 1

    for omkey, series in old.items(): # Iterates over old photons, using both the omkey and the series of photons associated
        stats["doms"] += 1
        if omkey not in geo:
            new[omkey] = series # If omkey is not found in geo module map, keep original times
            stats["bad_doms"] += 1
            continue
        dompos = geo[omkey].pos # Grabs position of omkey, I3Position vector
        Epos = electron.pos # Grabs position of electron, I3Position vector
        diff = dompos - Epos # Displacement vector from optical module to electron
        
        dist = np.linalg.norm(diff) # Magnitude of displacement vector
        # Construction of electron direction unit vector and calculation of angle between electron travel vector and displacement vector
        # zenith = electron.dir.zenith
        # azimuth = electron.dir.azimuth
        # Ex = np.sin(zenith) * np.cos(azimuth)
        # Ey = np.sin(zenith) * np.sin(azimuth)
        # Ez = np.cos(zenith)
        # Eangle = np.array([Ex, Ey, Ez])
        # phiE = np.arccos(np.dot(diff, Eangle) / dist)
        phiE = np.arccos((diff * electron.dir) / dist) # I3Direction object has built in operator for dot product

        
        # pdf = evalPdf(spline, dist, phiE, tgrid) # Evaluates the pdf (true pdf, not log) of the spline in 4d
        
        # # pdf = evalPdf(spline, dist, phiE, tgrid) # Evaluates the pdf (true pdf, not log) of the spline in 3d
        # pdf = np.clip(pdf, 0, None) # Clips the pdf to be non-negative, as negative values are not valid for a pdf
        # tot = pdf.sum()
        # if not np.isfinite(tot) or tot <= 0:
        #     stats["doms_skipped"] += 1
        #     new[omkey] = series  # Keeps original times, nothing to sample from
        #     continue
        # cdf = np.cumsum(pdf) # Calculates the cumulative distribution function
        # cdf /= cdf[-1] # Normalizes the cdf to be between 0 and 1

        ns = simclasses.I3CompressedPhotonSeries() # Initializes a new photonseries object to hold the resampled photons
        rows = []
        for p in series: # Iterates over the original series of photons and samples new tres values from the cdf
            photonpos = p.pos
            photondir = p.dir

            dphi = max(-1.0, min(1.0, np.arccos((-photonpos*photondir) / photonpos.magnitude)))

            pdf = evalPdf(spline, dist, phiE, dphi, tgrid) # Evaluates the pdf (true pdf, not log) of the spline in 4d
                   
            # pdf = evalPdf(spline, dist, phiE, tgrid) # Evaluates the pdf (true pdf, not log) of the spline in 3d
            pdf = np.clip(pdf, 0, None) # Clips the pdf to be non-negative, as negative values are not valid for a pdf
            tot = pdf.sum()
            if not np.isfinite(tot) or tot <= 0:
                stats["doms_skipped"] += 1
                new[omkey] = series  # Keeps original times, nothing to sample from
                continue
            cdf = np.cumsum(pdf) # Calculates the cumulative distribution function
            cdf /= cdf[-1] # Normalizes the cdf to be between 0 and 1

            old_tres = p.time - electron.time - 1.34*dist / dataclasses.I3Constants.c # Uses the old photons to calculate the actual time residual
            new_tres = float(np.interp(rng.random(), cdf, tgrid)) # Ransomly samples a time residual from the cdf of the spline
            rows.append( # Adds a new photon with the desired information to the rows list
                (
                    p.time - old_tres + new_tres,
                    p.weight,
                    p.wavelength,
                    p.dir.zenith,
                    p.dir.azimuth,
                    p.pos,
                )
            )
            stats["pulses"] += 1
        rows.sort(key=lambda r: r[0])
        for t, w, lamb, zen, az, ph_pos in rows: # Creates new compressed photon with the information in rows, including new tres
            np_ = simclasses.I3CompressedPhoton()
            np_.time, np_.weight, np_.wavelength, np_.dir, np_.pos = (
                t,
                w,
                lamb,
                dataclasses.I3Direction(zen, az),
                ph_pos,
            )
            ns.append(np_) # adds the new photons to the new series
        new[omkey] = ns # assigns new series object to the same omkey in the new map
        # new[omkey] = rows

    # del frame['new_photons']
    frame["new_photons"] = new # saves the new map to the frame
    return True


tray.AddModule(resample, streams=[icetray.I3Frame.DAQ]) # Adds the resample module to the tray
tray.AddModule( # Adds writer module, which actually writes new frames to the output file
    "I3Writer",
    filename=outfile,
    Streams=[icetray.I3Frame.DAQ],
    DropOrphanStreams=[icetray.I3Frame.Calibration],
)
tray.Execute()
tray.Finish()
