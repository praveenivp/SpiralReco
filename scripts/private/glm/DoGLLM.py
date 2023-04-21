#!/usr/bin/env python

import os,sys,itertools,subprocess,glob,re


# script to generate fsl design files with custom paradigm
# input parameters are nifti volume file, template design file, custom paradigm file and smoothfac
# if smoothfac is not passed it is set to 1

#usage: runFSL.sh <moco_nifiti_volume>  <paradigm folder> <smoothfac>

# check whether user at least 3 arguments
if len(sys.argv) < 2:
    print "Usage: %s <fullpath to nifti volume>  <fullpath to custom paradigm file>" % sys.argv[0]
    sys.exit(1)

# check whether the input files exist
if not os.path.isdir(sys.argv[1]):
    print "Error: %s does not exist" % sys.argv[1]
    sys.exit(1)
if not os.path.isdir(sys.argv[2]):
    print "Error: %s does not exist" % sys.argv[2]
    sys.exit(1)






def RunGLM(fullfile,templatefile,stimfile,smoothfac=1.0):


	# Parse file name components
	pn = os.path.dirname(fullfile)
	fn = os.path.basename(fullfile)
	ext = os.path.splitext(fn)[1][1:]
	basefn = os.path.splitext(fn)[0]

	# Get TR, resolution, smooth size, number of volumes, and highpass filter value
	FSLContainer="singularity exec -B /ptmp /home/pvalsala/ptmp/MyContainers/fsl_ubuntu4_devel-2022-06-09-e06711f0139f.sif"
	#print(FSLContainer+ " fslval "+ fullfile+ " pixdim4")
	os.system("module load singularity")
	TR = float(subprocess.check_output(FSLContainer+ " fslval "+ fullfile+ " pixdim4", shell=True))
	res = float(subprocess.check_output(FSLContainer+ " fslval "+ fullfile+ " pixdim1",shell=True))
	smoothsize = smoothfac * res
	nVol = int(subprocess.check_output(FSLContainer+ " fslval "+ fullfile+ " dim4",shell=True))
	Highpass = 30 # 30seconds is stimulus specific



	smoothlabel = "_{:.1f}_".format(smoothsize).replace('.', 'p')
	outputdir = os.path.join(os.getcwd(), "smooth{}{}".format(smoothlabel, basefn))
	designFile = os.path.join(os.getcwd(), "smooth{}{}.fsf".format(smoothlabel, basefn))
	#os.system(" ".join(["cp",templatefile,designFile]))



	#List all find and replace pairs
	pattern=["set\sfmri\(tr\).*","set\sfmri\(smooth\).*", "set\sfmri\(npts\).*","set\sfmri\(custom1\).*","set\sfeat_files\(1\).*","set\sfmri\(outputdir\).*" ]
	subst=["set fmri(tr) {:f}".format(TR),"set fmri(smooth) {:f}".format(smoothsize) , "set fmri(npts) {:d}".format(nVol), "set fmri(custom1) \"" + stimfile + "\"","set feat_files(1) \"" + fullfile + "\"",   "set fmri(outputdir) \"" + outputdir + "\""]

	#print("\n".join(subst)) #debug

	old_file=open(templatefile)
	new_file=open(designFile,'w+')
	for line in old_file:
	  writeLine=True
	  for cSubst,cPattern in itertools.izip_longest(subst,pattern):
	    if re.search(cPattern, line):
	      #print("writing " + cSubst + "\n\n")
	      new_file.write(cSubst)
	      writeLine=False
	      break
	  if(writeLine):
	    new_file.write(line)
	old_file.close()
	new_file.close()




	# Print modified lines
	print("Modified lines in deignFile : " + designFile)

	os.system(" ".join(["grep", "-E", "\"set fmri\(tr\)*\"", designFile]))
	os.system(" ".join(["grep", "-E", "\"set fmri\(npts*\"", designFile]))
	os.system(" ".join(["grep", "-E", "\"set feat_files(1)*\"", designFile]))
	os.system(" ".join(["grep", "-E", "\"set fmri\(outputdir*\"", designFile]))
	os.system(" ".join(["grep", "-E", "\"set fmri\(custom1*\"", designFile]))


	#run feat in background
	print(" ".join(["Running:",FSLContainer, "feat", designFile, "&"]))
	os.system(" ".join([FSLContainer, "feat", designFile, "&"]))

############################################################################################################
############################################################################################################



data_folder = sys.argv[1]
stim_folder = sys.argv[2]

#iterate through the files in the data folder with *.nii extension
for filename in os.listdir(data_folder):
    if filename.endswith(".nii"):

        data_filepath = os.path.join(data_folder, filename)
        # Extract the number between "M" and "_" in the filename
        MeasUID = filename.split("M")[1].split("_")[0]
        #print(MeasUID)
        # Build the filename of the corresponding file in the target folder with "*M<MEASUID>_*.nii" pattern
        stim_filename = "*M" + MeasUID + "_*.txt"
        stim_filepath = os.path.join(stim_folder, stim_filename)
        stim_filepath = glob.glob(stim_filepath)[0]
        # Check if the file exists in the target folder
        if os.path.exists(stim_filepath):

	    RunGLM(data_filepath,"template.fsf",stim_filepath,0.0)
	    RunGLM(data_filepath,"template.fsf",stim_filepath,1.0)
        else:
            # Move on to the next file in the source folder
            print("expected stim file :" + stim_filepath)
            continue








