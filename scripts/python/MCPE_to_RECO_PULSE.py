#!/usr/bin/env python3
from icecube import icetray, dataclasses, dataio, simclasses
from I3Tray import I3Tray
import argparse

"""
From running a separate script I understand that
I3MCPE
+ Has 3 attributes
++ time (float, ns)
++ npe (int) - describes number of photoelectrons in this hit
++ ID - particle ID for monte carlo?

I3RecoPulse
+ Has 5 attributes
++ time (double, ns)
++ charge (float, PE) - how much charge in photoelectron units
++ width (float, ns) - pulse width = duration of the pulse
++ flags (int) - this is a bit mask I don't need it
++ PulseFlags (enum) - for flag values
"""

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--input-prefix",
    required=True,
    help="Path prefix; full input path is <prefix><run>.i3.zst",
)
parser.add_argument("-r", "--run", type=int, required=True, help="Run / array index")
args = parser.parse_args()

# GCD file
GCDFILE = "/mnt/home/jalabadz/myWork/myScripts/DATA/gcdfile/PONE_800mGrid.i3.gz"

# Build the single input file this task is responsible for
INFILE = f"{args.input_prefix}{args.run:03d}.i3.zst"

# Where to save the output
OUTFILE = f"/mnt/home/jalabadz/myWork/myScripts/RESULTS/mcpe2pulses_RESULT/output_recopulses_{args.run:03d}.i3.zst"

# Just for cleanness I put this here
MCPE_KEY = "Accepted_MCPEMap"


def mcpe_to_pulses(frame):
    # Guard against mistakes
    if MCPE_KEY not in frame:
        return

    # empty container for output. Map of OMKey -> list of pulses
    pulses = dataclasses.I3RecoPulseSeriesMap()

    """
    We have a dict with key = OMKey and val = lists of MCPE hits 
    for a particular DOM.
    
    Loop over each DOM. 
    
    omkey identifies which DOM we're on.
    
    mcpes is the list of photoelectron hits for this DOM
    """
    for omkey, mcpes in frame[MCPE_KEY].items():
        # for THIS DOM, make a list of I3RecoPulseSeries. Fill it one pulse per MCPE
        series = dataclasses.I3RecoPulseSeries()

        # for every MCPE hit on this DOM
        # make an I3RecoPulse
        # copy the time across (both in ns)
        # copy the charge across
        # set the pulse width to a default of 1
        # append the series for the DOM
        for hit in mcpes:
            pulse = dataclasses.I3RecoPulse()
            pulse.time = hit.time
            pulse.charge = float(
                hit.npe
            )  # npe is int, charge is float so have to convert
            pulse.width = 1.0  # default
            series.append(pulse)

        # Add this DOM to the dict if it has pulses
        if len(series) > 0:
            pulses[omkey] = series

    # The finished dict gets to be in the frame with key "I3RecoPulses"
    frame["MCPulses"] = (
        pulses  # SHOULD HAVE A DIFFERENT NAME, THE NAME IS HOW WE IDENTIFY WHICH PULSE SERIES IT'S GOING TO BE. IN THIS CASE CALL THEM MCPulses or something
    )


tray = I3Tray()
tray.Add("I3Reader", FilenameList=[GCDFILE, INFILE])
tray.Add(mcpe_to_pulses, Streams=[icetray.I3Frame.DAQ])
tray.Add("I3Writer", Filename=OUTFILE)
tray.Execute()
