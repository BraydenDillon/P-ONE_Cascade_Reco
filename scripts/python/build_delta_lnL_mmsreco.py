import sys
import os

sys.path.insert(0, os.path.abspath(".."))
import photospline
from icecube import dataio, dataclasses, icetray, gulliver, simclasses, phys_services
import numpy as np
from scipy import optimize, integrate
import argparse
from scripts.python.SplineEval import evalPdf
import random

# This file takes a file with log likelihood information and calculates the delta log likelihood between the best fit and the MCTruth
# It then compiles the data into .npy files for easy analysis in an interactive notebook

parser = argparse.ArgumentParser(
    description="Takes I3photons and I3Electrons from simulation files to read out positional data"
)
parser.add_argument(
    "-i", "--infile", default="/mnt/scratch/dillonb5/mmsreco_7-20/llhfit_" # This is the line I change when I run on different files
        # Change here instead of in submission script so that submit script never needs to change as it doesn't have input file arguments 
)
parser.add_argument("-r", "--runnumber", type=int, default=1)

args = parser.parse_args()
runnumber = -999
if args.runnumber < 10:
    runnumber = "00" + str(args.runnumber)
elif args.runnumber < 100:
    runnumber = "0" + str(args.runnumber)
else:
    runnumber = str(args.runnumber)
infile = args.infile + runnumber + ".i3.zst"
# No outfile argument assembly this time, output is npy file not i3 file
delta_logL = [] # list to be appended to 

tray = icetray.I3Tray()
tray.AddModule("I3Reader", "reader", FilenameList=[infile]) # no need for a gcd file here. 


def calculate_dlnL(frame):
    if (len(frame['new_photons']) != 0) and (len(frame['I3MCTree']) != 0):
        if not np.isnan(frame['LLHFit_step5FitParams'].logl):
            bestfit = frame['LLHFit_step5FitParams'].logl
        elif not np.isnan(frame['LLHFit_step4FitParams'].logl):
            bestfit = frame['LLHFit_step4FitParams'].logl
        elif not np.isnan(frame['LLHFit_step3FitParams'].logl):
            bestfit = frame['LLHFit_step3FitParams'].logl
        elif not np.isnan(frame['LLHFit_step2FitParams'].logl):
            bestfit = frame['LLHFit_step2FitParams'].logl
        if not np.isnan(frame['LLHFit_step1FitParams'].logl):
            bestfit = frame['LLHFit_step1FitParams'].logl
        else:
            return

        delta_logL.append(frame['LLHFit_mctruth'].logl - bestfit)

tray.AddModule(calculate_dlnL, Streams = [icetray.I3Frame.Physics])

tray.Execute()
tray.Finish()
ary = np.array(delta_logL)
np.save("/mnt/scratch/dillonb5/7-20_logL/delta_ary_" + runnumber +".npy", ary) # saves list to npy file
        