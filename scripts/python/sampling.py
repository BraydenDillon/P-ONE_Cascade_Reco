from icecube import dataclasses, phys_services, icetray, simclasses, dataio
import numpy as np
import sys
import os

sys.path.insert(0, os.path.abspath(".."))
from scripts.python.SplineEval import evalPdf
import photospline
import argparse
from scipy.stats import norm

parser = argparse.ArgumentParser(description="Takes a .i3 file and a gcd file")
parser.add_argument(
    "-i", "--infile", default="/mnt/home/dillonb5/cascades/sacrifice_data/gen_"
)
parser.add_argument(
    "-g", "--gcdfile", default="/mnt/home/dillonb5/cascades/gcdfile/PONE_800mGrid.i3.gz"
)
parser.add_argument("-r", "--runnumber", type=int, default=50)
parser.add_argument(
    "-o", "--outfile", default="/mnt/scratch/dillonb5/cdf_per_photon/new_"
)
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
outfile = args.outfile + runnumber + ".i3.zst"


tray = icetray.I3Tray()
tray.AddModule("I3Reader", "reader", FilenameList=[gcdfile, infile])

c = 299792458

spline = photospline.SplineTable("/mnt/home/dillonb5/cascades/fits/splinelog_3D.fits")
tgrid = np.linspace(min(spline.knots[-1]), max(spline.knots[-1]), 1000)
# tgrid = np.linspace(-5, 5, 1000)
N_GROUP = 1.34
N_PHASE = 1.35557

rng = np.random.default_rng()


def resample(frame):
    stats = dict(events=0, doms=0, pulses=0, doms_clamped=0, doms_skipped=0, bad_doms=0)
    if "I3Photons" not in frame or "I3MCTree" not in frame:
        return True
    electron = dataclasses.I3Particle(frame["I3MCTree"][1])
    electron.shape = dataclasses.I3Particle.Cascade
    # mu.shape = dataclasses.I3Particle.InfiniteTrack
    geo = frame["I3ModuleGeoMap"]
    old = frame["I3Photons"]
    omkeys = frame["I3Photons"].keys()
    new = simclasses.I3CompressedPhotonSeriesMap()
    stats["events"] += 1

    for omkey, series in old.items():
        stats["doms"] += 1
        if omkey not in geo:
            new[omkey] = series
            stats["bad_doms"] += 1
            continue
        pos = geo[omkey].pos
        Epos = electron.pos
        # diff = Epos - pos
        # # dist = phys_services.I3Calculator.cherenkov_distance(
        # #     electron, pos, N_GROUP, N_PHASE
        # # )
        # dist = np.linalg.norm(diff)
        # zenith = electron.dir.zenith
        # azimuth = electron.dir.azimuth
        # Ex = np.sin(zenith) * np.cos(azimuth)
        # Ey = np.sin(zenith) * np.sin(azimuth)
        # Ez = np.cos(zenith)
        # Eangle = np.array([Ex, Ey, Ez])
        # phiE = np.arccos(np.dot(diff, Eangle) / dist)

        # pdf = evalPdf(spline, dist, phiE, tgrid)
        # # dist_for_pdf = 20.0
        # # phiE_for_pdf = 0.3
        # # pdf = evalPdf(spline, dist_for_pdf, phiE_for_pdf, tgrid)
        # # pdf = norm.pdf(tgrid)
        # pdf = np.clip(pdf, 0, None)
        # tot = pdf.sum()
        # if not np.isfinite(tot) or tot <= 0:
        #     stats["doms_skipped"] += 1
        #     new[omkey] = series  # keep original times, nothing to sample from
        #     continue
        # cdf = np.cumsum(pdf)
        # cdf /= cdf[-1]

        ns = simclasses.I3CompressedPhotonSeries()
        rows = []
        for p in series:
            # p_pos = dataclasses.I3Position(Epos.x + 20, Epos.y, Epos.z)
            # p_angle = dataclasses.I3Direction(zenith + 0.3, azimuth)
            # old_tres = phys_services.I3Calculator.time_residual(
            #     electron, pos, p.time, N_GROUP, N_PHASE
            # )
            diff = Epos - pos + p.pos
        # dist = phys_services.I3Calculator.cherenkov_distance(
        #     electron, pos, N_GROUP, N_PHASE
        # )
            dist = np.linalg.norm(diff)
            zenith = electron.dir.zenith
            azimuth = electron.dir.azimuth
            Ex = np.sin(zenith) * np.cos(azimuth)
            Ey = np.sin(zenith) * np.sin(azimuth)
            Ez = np.cos(zenith)
            Eangle = np.array([Ex, Ey, Ez])
            phiE = np.arccos(np.dot(diff, Eangle) / dist)

            pdf = evalPdf(spline, dist, phiE, tgrid)
        # dist_for_pdf = 20.0
        # phiE_for_pdf = 0.3
        # pdf = evalPdf(spline, dist_for_pdf, phiE_for_pdf, tgrid)
        # pdf = norm.pdf(tgrid)
            pdf = np.clip(pdf, 0, None)
            tot = pdf.sum()
            if not np.isfinite(tot) or tot <= 0:
                stats["doms_skipped"] += 1
                new[omkey] = series  # keep original times, nothing to sample from
                continue
            cdf = np.cumsum(pdf)
            cdf /= cdf[-1]

            old_tres = p.time - electron.time - 1.34*dist / c * 1e9
            new_tres = float(np.interp(rng.random(), cdf, tgrid))
            rows.append(
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
        for t, w, lamb, zen, az, ph_pos in rows:
            np_ = simclasses.I3CompressedPhoton()
            np_.time, np_.weight, np_.wavelength, np_.dir, np_.pos = (
                t,
                w,
                lamb,
                dataclasses.I3Direction(zen, az),
                ph_pos,
            )
            ns.append(np_)
        new[omkey] = ns
        # new[omkey] = rows

    # del frame['new_photons']
    frame["new_photons"] = new
    return True


tray.AddModule(resample, streams=[icetray.I3Frame.DAQ])
tray.AddModule(
    "I3Writer",
    filename=outfile,
    Streams=[icetray.I3Frame.DAQ],
    DropOrphanStreams=[icetray.I3Frame.Calibration],
)
tray.Execute()
tray.Finish()
