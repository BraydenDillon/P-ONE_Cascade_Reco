from icecube import icetray, dataio, dataclasses, phys_services
from icecube import gulliver, lilliput, gulliver_modules
from icecube.icetray import I3Tray
from I3Tray import *
import argparse

icetray.load("mmsreco")

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--infile", default="/mnt/scratch/dillonb5/mmsreco_positive/llhfit_fixedDiff_080.i3.zst")
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

def lazyfilter(frame):
    if frame['I3EventHeader'].event_id != 22:
        return False
    return frame['I3EventHeader'].event_id == 22

tray.AddModule(lazyfilter, streams = [icetray.I3Frame.Physics])

pulses = "new_photons"
seed = "MCTruth"

tray.AddService("MMSLikelihoodFactory", "mmsreco_truth",
                InputPhotons=pulses,  ExpectNoise=False, ConvolutionWidth=0.0,
                SplineTablePath="/mnt/home/dillonb5/cascades/fits/splinelog_3D.fits")

tray.AddModule("I3LogLikelihoodCalculator", "LLHFit_mctruth_2",
        LogLikelihoodService = "mmsreco_truth",
        FitName              = "MCTruth",
        NFreeParameters      = 6,
        If = lambda frame: pulses in frame and "MCTruth" in frame)


tray.AddModule("I3Writer", Filename=outfile)
tray.Execute()
tray.Finish()
# check what channels are involved and how they contribute to the likelihood