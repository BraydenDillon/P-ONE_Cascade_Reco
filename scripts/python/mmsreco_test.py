from icecube import icetray, dataio, dataclasses, phys_services
from icecube import gulliver, lilliput, gulliver_modules
from icecube.icetray import I3Tray
from I3Tray import *
import argparse

# import numpy as np
icetray.load("mmsreco")

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--infile", default="/mnt/scratch/dillonb5/sampled_7-21_fine_bin/new_") # I change the input file here, so that I don't need to change two scripts
parser.add_argument("-r", "--runnumber", type=int, default=1)
parser.add_argument(
    "-g", "--gcdfile", default="/mnt/home/dillonb5/cascades/gcdfile/PONE_800mGrid.i3.gz"
)
parser.add_argument('-o', "--outfile", default = "/mnt/scratch/dillonb5/mmsreco_7-21_fine_bin/llhfit_") # Output file name, will be appended with run number

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

tray = I3Tray()
tray.AddModule("I3Reader", "reader", Filenamelist=[gcdfile, infile])
tray.AddModule("I3NullSplitter", "Splitter")
# splinepath = '/mnt/research/IceCube/jalabadz/iter_6.0_4D_I3Photons/spline_result/splinelog.fits' # Self-explanatory but this is the path to 4d spline fit. Needs to switch for 3d fits
splinepath = '/mnt/home/dillonb5/cascades/fits/splinelog_3D.fits' # Self-explanatory but this is the path to 3d spline fit. Needs to switch for 4d fits

def mctruth(fr): # Function that will be added to the tray. Extracts the MCTruth from the MCtree and adds it to the frame independently
    fr["MCTruth"] = fr["I3MCTree"][1] # For 3d we treat the electron (second entry in mctree) as truth particle. For 4d we treat the neutrino itself as the truth which is the first particle.
    fr["MCTruth"].fit_status = fr["MCTruth"].OK

# def mchadrons(fr):
#     fr["MCHadrons"] = fr["I3MCTree"][2]
#     fr["MCHadrons"].fit_status = fr["MCHadrons"].OK


tray.Add(mctruth, Streams=[icetray.I3Frame.DAQ]) # adds mctruth function to tray


pulses = "new_photons" # 
seed = "MCTruth"
tray.AddService("I3BasicSeedServiceFactory", "seed1", FirstGuess=seed) # First step of convolution. Seed is MCTruth
tray.AddService( # Initializes simplex minimizer to be used in fitter module
    "I3GSLSimplexFactory",
    "minimizeit",
    Tolerance=0.001,
    SimplexTolerance=0.001,
    MaxIterations=50000,
)
tray.AddService( # Parameterizes the space for fitter module. Essentially translates coordinates from hypothesis to likelihood space
    "I3SimpleParametrizationFactory",
    "simpleparam",
    StepX=10 * I3Units.m,
    StepY=10 * I3Units.m,
    StepZ=10 * I3Units.m,
    StepT=10 * I3Units.ns,
    StepZenith=0.7 * I3Units.deg,
    StepAzimuth=0.7 * I3Units.deg,
)  # , BoundsZenith=[0, 180*I3Units.degree], BoundsAzimuth=[0, 360*I3Units.degree])

# tray.AddService("I3HalfSphereParametrizationFactory", "simpleparam",
#         DirectionStepsize = 0.3 * I3Units.deg,
#         VertexStepsize = 5 * I3Units.m,
#         TimeStepsize = 1 * I3Units.ns)

tray.AddService( # calls our edited mmsreco likelihood service
    "MMSLikelihoodFactory",
    "mms_step1",
    InputPhotons=pulses,
    SplineTablePath=splinepath,
    ExpectNoise=False,
    ConvolutionWidth=20.0 # Start with wide convolution width to find the minimum efficiently, then later narrow down
)
tray.AddModule( # Compiles all services into fitter module to find best fit for likelihood 
    "I3SimpleFitter",
    "LLHFit_step1",
    SeedService="seed1",
    Parametrization="simpleparam",
    LogLikelihood="mms_step1",
    Minimizer="minimizeit",
    If=lambda frame: pulses in frame and seed in frame,
    OutputName="LLHFit_step1",
)

#########################################################################################
# Step 2
tray.AddService(
    "I3BasicSeedServiceFactory", "seed2", firstguess = "LLHFit_step1"
)
tray.AddService(
    "MMSLikelihoodFactory",
    "mms_step2",
    InputPhotons=pulses,
    SplineTablePath=splinepath,
    ExpectNoise=False,
    ConvolutionWidth=15.0
)

tray.AddModule(
    "I3SimpleFitter",
    "LLHFit_step2",
    SeedService="seed2",
    Parametrization="simpleparam",
    LogLikelihood="mms_step2",
    Minimizer="minimizeit",
    If=lambda frame: pulses in frame and 'LLHFit_step1' in frame,
    OutputName="LLHFit_step2"
)

########################################################################################
# Step 3
tray.AddService(
    "I3BasicSeedServiceFactory", "seed3", firstguess = "LLHFit_step2"
)
tray.AddService(
    "MMSLikelihoodFactory",
    "mms_step3",
    InputPhotons=pulses,
    SplineTablePath=splinepath,
    ExpectNoise=False,
    ConvolutionWidth=10.0
)

tray.AddModule(
    "I3SimpleFitter",
    "LLHFit_step3",
    SeedService="seed3",
    Parametrization="simpleparam",
    LogLikelihood="mms_step3",
    Minimizer="minimizeit",
    If=lambda frame: pulses in frame and 'LLHFit_step2' in frame,
    OutputName="LLHFit_step3"
)

#####################################################################################
# Step 4
tray.AddService(
    "I3BasicSeedServiceFactory", "seed4", firstguess = "LLHFit_step3"
)
tray.AddService(
    "MMSLikelihoodFactory",
    "mms_step4",
    InputPhotons=pulses,
    SplineTablePath=splinepath,
    ExpectNoise=False,
    ConvolutionWidth=5.0
)

tray.AddModule(
    "I3SimpleFitter",
    "LLHFit_step4",
    SeedService="seed4",
    Parametrization="simpleparam",
    LogLikelihood="mms_step4",
    Minimizer="minimizeit",
    If=lambda frame: pulses in frame and 'LLHFit_step3' in frame,
    OutputName="LLHFit_step4"
)

#####################################################################################
# Step 5, conv width zero to find exact minimum
tray.AddService(
    "I3BasicSeedServiceFactory", "seed5", firstguess = "LLHFit_step4"
)
tray.AddService(
    "MMSLikelihoodFactory",
    "mms_step5",
    InputPhotons=pulses,
    SplineTablePath=splinepath,
    ExpectNoise=False,
    ConvolutionWidth=0
)

tray.AddModule(
    "I3SimpleFitter",
    "LLHFit_step5",
    SeedService="seed5",
    Parametrization="simpleparam",
    LogLikelihood="mms_step5",
    Minimizer="minimizeit",
    If=lambda frame: pulses in frame and 'LLHFit_step4' in frame,
    OutputName="LLHFit_step5"
)

#################################################################################################
# No fit, finds the likelihood of the MCTruth Hypothesis
tray.AddService("MMSLikelihoodFactory", "mmsreco_truth",
                InputPhotons=pulses,  ExpectNoise=False, ConvolutionWidth=0.0,
                SplineTablePath=splinepath)

tray.AddModule("I3LogLikelihoodCalculator", "LLHFit_mctruth",
        LogLikelihoodService = "mmsreco_truth",
        FitName              = "MCTruth",
        NFreeParameters      = 6,
        If = lambda frame: pulses in frame and "MCTruth" in frame)


tray.AddModule("I3Writer", Filename=outfile)
tray.Execute()
tray.Finish()
