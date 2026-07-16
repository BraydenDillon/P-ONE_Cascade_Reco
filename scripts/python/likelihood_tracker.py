from icecube import icetray, dataio, dataclasses, phys_services
from icecube import gulliver, lilliput, gulliver_modules
from icecube.icetray import I3Tray
from I3Tray import *
import argparse

icetray.load("mmsreco")

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--infile", default="/mnt/scratch/dillonb5/const_tres/new_080.i3.zst")
# parser.add_argument("-r", "--runnumber", type=int, default=1)
parser.add_argument(
    "-g", "--gcdfile", default="/mnt/home/dillonb5/cascades/gcdfile/PONE_800mGrid.i3.gz"
)

parser.add_argument('-o', "--outfile", default = "/mnt/scratch/dillonb5/throwaway.i3.zst")

args = parser.parse_args()
# runnumber = -999
# if args.runnumber < 10:
#     runnumber = "00" + str(args.runnumber)    
# elif args.runnumber < 100:
#     runnumber = "0" + str(args.runnumber)
# else:
#     runnumber = str(args.runnumber)
infile = args.infile #+ runnumber + ".i3.zst"
gcdfile = args.gcdfile
outfile = args.outfile #+ runnumber + ".i3.zst"

tray = I3Tray()
tray.AddModule("I3Reader", "reader", Filenamelist=[gcdfile, infile])
tray.AddModule("I3NullSplitter", "Splitter")


splinepath = "/mnt/home/dillonb5/cascades/fits/splinelog_3D.fits"
#splinepath = "/mnt/research/IceCube/jalabadz/iter_6.0_4D_I3Photons/spline_result/splinelog.fits"

def lazyfilter(frame):
    if frame['I3EventHeader'].event_id != 22:
        return False
    return frame['I3EventHeader'].event_id == 22
def mctruth(fr):
    fr["MCTruth"] = fr["I3MCTree"][1]
    fr["MCTruth"].fit_status = fr["MCTruth"].OK

tray.AddModule(mctruth, streams = [icetray.I3Frame.DAQ])

tray.AddModule(lazyfilter, streams = [icetray.I3Frame.Physics])

pulses = "new_photons"
seed = "MCTruth"
# tray.AddService("I3BasicSeedServiceFactory", "seed1", FirstGuess=seed)
# tray.AddService(
#     "I3GSLSimplexFactory",
#     "minimizeit",
#     Tolerance=0.001,
#     SimplexTolerance=0.001,
#     MaxIterations=50000,)

# tray.AddService("I3SimpleParametrizationFactory", "simpleparam",
#     StepX=10 * I3Units.m,
#     StepY=10 * I3Units.m,
#     StepZ=10 * I3Units.m,
#     StepT=10 * I3Units.ns,
#     StepZenith=0.7 * I3Units.deg,
#     StepAzimuth=0.7 * I3Units.deg,)

# tray.AddService(
#     "MMSLikelihoodFactory",
#     "mms_step1",
#     InputPhotons=pulses,
#     SplineTablePath=splinepath,
#     ExpectNoise=False,
#     ConvolutionWidth=0
# )

# tray.AddModule(
#     "I3SimpleFitter",
#     "LLHFit_step1",
#     SeedService="seed1",
#     Parametrization="simpleparam",
#     LogLikelihood="mms_step1",
#     Minimizer="minimizeit",
#     If=lambda frame: pulses in frame and seed in frame,
#     OutputName="LLHFit",
# )

tray.AddService("MMSLikelihoodFactory", "mmsreco_truth",
                InputPhotons=pulses,  ExpectNoise=False, ConvolutionWidth=0.0,
                SplineTablePath=splinepath)

tray.AddModule("I3LogLikelihoodCalculator", "LLHFit_mctruth",
        LogLikelihoodService = "mmsreco_truth",
        FitName              = "MCTruth",
        NFreeParameters      = 6,
        If = lambda frame: pulses in frame and "MCTruth" in frame)




tray.AddModule("I3Writer", Filename = outfile)
tray.Execute()
tray.Finish()