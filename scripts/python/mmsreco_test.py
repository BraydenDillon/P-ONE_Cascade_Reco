from icecube import icetray, dataio, dataclasses, phys_services
from icecube import gulliver, lilliput, gulliver_modules
from icecube.icetray import I3Tray
from I3Tray import *
import argparse

# import numpy as np
icetray.load("mmsreco")

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--infile")
parser.add_argument(
    "-g", "--gcdfile", default="/mnt/home/dillonb5/cascades/gcdfile/PONE_800mGrid.i3.gz"
)

args = parser.parse_args()
infile = args.infile
gcdfile = args.gcdfile

tray = I3Tray()
tray.AddModule("I3Reader", "reader", Filenamelist=[gcdfile, infile])
tray.AddModule("I3NullSplitter", "Splitter")


def mctruth(fr):
    fr["MCTruth"] = fr["I3MCTree"][1]
    fr["MCTruth"].fit_status = fr["MCTruth"].OK


tray.Add(mctruth, Streams=[icetray.I3Frame.DAQ])


pulses = "I3Photons"
seed = "MCTruth"
tray.AddService("I3BasicSeedServiceFactory", "seed", FirstGuess=seed)
tray.AddService(
    "I3GSLSimplexFactory",
    "minimizeit",
    Tolerance=0.01,
    SimplexTolerance=0.01,
    MaxIterations=50000,
)
tray.AddService(
    "I3SimpleParametrizationFactory",
    "simpleparam",
    StepX=10 * I3Units.m,
    StepY=10 * I3Units.m,
    StepZ=10 * I3Units.m,
    StepT=10 * I3Units.ns,
    StepZenith=0.3 * I3Units.deg,
    StepAzimuth=0.3 * I3Units.deg,
)  # , BoundsZenith=[0, 180*I3Units.degree], BoundsAzimuth=[0, 360*I3Units.degree])
tray.AddService(
    "MMSLikelihoodFactory",
    "mms",
    InputPhotons=pulses,
    SplineTablePath="/mnt/home/dillonb5/cascades/fits/splinelog_3D.fits",
    ExpectNoise=False,
)
tray.AddModule(
    "I3SimpleFitter",
    "LLHFit",
    SeedService="seed",
    Parametrization="simpleparam",
    LogLikelihood="mms",
    Minimizer="minimizeit",
    If=lambda frame: pulses in frame and seed in frame,
    OutputName="LLHFit",
)

tray.AddModule("I3Writer", Filename="/mnt/scratch/dillonb5/mmsreco_test/test6.i3.gz")
tray.Execute(50)
tray.Finish()
