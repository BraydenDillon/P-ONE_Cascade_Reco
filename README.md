This repo is used to work on cascade event reconstruction for the Pacific Ocean Neutrino Experiment

## Basic Framework ##
There are multiple scripts in this repo that are used for various purposes. Some build on each other and require different scripts to be run in series to get a result. I'll get into that but what you need to know at a basic level is that all of the scripts in the scripts/python directory are source code, and they are executed via slurm workload manager using submission scipts in scripts/submit_scripts. These were written to be run on Michigan State University's High Performance Computing Cluster (HPCC) which uses slurm as its workload manager. I'll be talking about each python file individually, and in each section there will be a reference to the submission script used to run the program. 

## MMSReco ## 
First and formost, this is the most important script written in this repo, as it is the (hopefully) working reconstruction for cascade events. The python script is **mmsreco_test.py** and it is run with **run_mmsreco.sb**. This script takes an input of a .i3 file containing monte-carlo (or real, if people are using this when the detector is taking data) event information and takes it through five convolution steps of calculating the log-likelihood and optimizing to the best fit electron. Inputs also include a gcd file, a run number, and optionally an outfile if you want to change from the default. 

The script is set up to read a set of files within a directory, so the infile should take the form */path-to-directory/file_name_* with the number of the file intentionally left out (The standard for MC event files has followed the syntax of *file_name_{file_number}*). The runnumber argument will then be appended to the input file name for fetching within the script, as well as the file suffixes (e.g. .i3.zst). This may make running on single files inconvenient but makes it much easier to run on many files with a single array job submission. 

The GCD and output file arguments have defaults, so are not required but can be changed when the function is called or the defaults can be edited. So if you wanted to run this script as-is on say, file gen_100.i3.gz, you'd run the command **python3 mmsreco_test.py -i /path-to-directory/gen_ -r 100**. 

mmsreco_test sets up the required icetray services to be called by the I3SimpleFitter module. BasicSeedServiceSeedservice creates the seed, SimpleParameterizationFactory converts the parameters from hypothesis coordinates to likelihood coordinates, LikelihoodFactory is what actually computes the likelihood with mmsreco, and finally SimpleFitter iterates until the likelihood is optimized. Everything but mmsreco is untouched from the public version of icetray.

I have "written" two versions of the mmsreco script that I keep in the **Edited_Icetray** directory. The two versions are 3d_mmsreco and 4d_mmsreco, appended by the date they were last updated. These scripts are not called or referenced anywhere, and only serve as alternate versions of the mmsreco script that I can look at and copy into the source directory when I want to edit the functionality of the script.

All versions of mmsreco also include a commented line with print statements tracking the likelihood contributions. These are being used to diagnose the problems with the reconstruction and can be disregarded if needed. 

