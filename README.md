This repo is used to work on cascade event reconstruction for the Pacific Ocean Neutrino Experiment

## Basic Framework ##
There are multiple scripts in this repo that are used for various purposes. Some build on each other and require different scripts to be run in series to get a result. I'll get into that but what you need to know at a basic level is that all of the scripts in the scripts/python directory are source code, and they are executed via slurm workload manager using submission scipts in scripts/submit_scripts. These were written to be run on Michigan State University's High Performance Computing Cluster (HPCC) which uses slurm as its workload manager. I'll be talking about each python file individually, and in each section there will be a reference to the submission script used to run the program. 

## MMSReco ## 
First and formost, this is the most important script written in this repo, as it is the (hopefully) working reconstruction for cascade events. The python script is **mmsreco_test.py** and it is run with **run_mmsreco.sb**. This script takes an input of a .i3 file containing monte-carlo (or real, if people are using this when the detector is taking data) event information and takes it through five convolution steps of calculating the log-likelihood and optimizing to the best fit electron. Inputs also include a gcd file, a run number, and optionally an outfile if you want to change from the default. 

The script is set up to read a set of files within a directory, so the infile should take the form */path-to-directory/file_name_* with the number of the file intentionally left out (The standard for MC event files has followed the syntax of *file_name_{file_number}*). The runnumber argument will then be appended to the input file name for fetching within the script, as well as the file suffixes (e.g. .i3.zst). This may make running on single files inconvenient but makes it much easier to run on many files with a single array job submission. 

The GCD and output file arguments have defaults, so are not required but can be changed when the function is called or the defaults can be edited. So if you wanted to run this script as-is on say, file gen_100.i3.gz, you'd run the command **python3 mmsreco_test.py -i /path-to-directory/gen_ -r 100**. 

mmsreco_test sets up the required icetray services to be called by the I3SimpleFitter module. BasicSeedServiceSeedservice creates the seed, SimpleParameterizationFactory converts the parameters from hypothesis coordinates to likelihood coordinates, LikelihoodFactory is what actually computes the likelihood with mmsreco, and finally SimpleFitter iterates until the likelihood is optimized. Everything but mmsreco is untouched from the public version of icetray.

These fits are stored in new .i3 files as llhfit_step{1-5} for the best fit particle hypotheses at each step of convolution, as well as in llhfit_step{1-5}FitParams storing the log likelihood and other metadata like the number of iterations the fitter takes to converge. The truth particle is also stored in MCTruth and the MCTruth log likelihood is stored in LLHFit_mctruth. 

I have "written" two versions of the mmsreco script that I keep in the **Edited_Icetray** directory. The two versions are 3d_mmsreco and 4d_mmsreco, appended by the date they were last updated. These scripts are not called or referenced anywhere, and only serve as alternate versions of the mmsreco script that I can look at and copy into the source directory when I want to edit the functionality of the script.

All versions of mmsreco also include a commented line with print statements tracking the likelihood contributions. These are being used to diagnose the problems with the reconstruction and can be disregarded if needed. 

## Sampling ##

This will not be needed in the actual use case for this reconstruction, but as a diagnostic I wrote a script to sample a new time residual from the spline to attempt to reconstruct it. By construction, the fit using the spline should be a good fit to data that was sampled from that same spline. 

The sampling used is a cdf sampling of the spline distribution of time residual values. First the pdf for a set of coordinates is calculated, then the cdf is calculated on a tgrid array using that pdf. A random number generator is used to sample the inverse of the cdf to obtain a time residual value. Then the actual time residual value is calculated and subtracted off of the photon time, then the sampled time residual is added back on. This new value is put into a new photon object (aptly called new_photons) along with the same information from the original photon save the original time. This creates a new absolute time that, when calculated with a reco, will return the sampled time resolution value.

The sampling source script is called **sampling.py** and it is run using **run_sampling.sb**. The input data, output file, run number, and gcd file are customizeable as inputs but as of now the spline needs to be changed within the actual function. 

Sampling is also currently configured to sample the 4d spline, which involves calculating an individual cdf for every photon due to the nature of the 4th dimension of the spline, impact angle (dphi). 

## Check Sampling ##

This has been used a few times for a few different checks that the sampling function atcually accurately samples from the spline. The most important test that has stayed the longest is the cdf test. The script steps through a t array and compiles pdfs to calculate a cumulative distribution function and produce a plot of the cdf calculated at every t value in the array. This histogram should, by construction, be a completely uniform distribution since it was sampled from the cdf originally. There is some commented code in there as well for tests where nonsense values were sampled instead to make sure sampling was doing anything at all. Of note there is code to re-assemble a gaussian sample, and code to graph a pdf when sampled data was all taken from the same set of coordinates (except fot tres). Unless there is uncertainty about a future sampling method, these can be ignored/discarded. 

Following the standard format, the python file is called **check_sampling.py** and the batch file that runs it is called **run_check.sb**. Another file called run_check_norm.sb exists but was only for the aformentioned gaussian sampling check. 

