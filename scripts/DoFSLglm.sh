#!/bin/sh
# script to generate fsl design files with custom paradigm
# input parameters are nifti volume file, template design file, custom paradigm file and smoothfac
# if smoothfac is not passed it is set to 1

# check whether user atleast 3 arguments
if [ $# -lt 3 ]
then
    echo "Usage: $0 <fullpath to nifti volume> <fullpath to template design file> <fullpath to custom paradigm file>"
    exit 1
fi
#check whether the input files exist
if [ ! -f "$1" ]
then
    echo "Error: $1 does not exist"
    exit 1
fi
if [ ! -f "$2" ]
then
    echo "Error: $2 does not exist"
    exit 1
fi
if [ ! -f "$3" ]
then
    echo "Error: $3 does not exist"
    exit 1
fi

#if fourth parameter is passed set smothfac to that otherwise set it to 0
if [ $# -eq 4 ]
then
    smoothfac=$4
else
    smoothfac=1
fi


fullfile=$1
templatefile=$2
stimfile=$3

#helper function to convert paths ro sed compactible string
sedConvert(){
    <<< "$1" sed -e 's`[][\\/.*^$]`\\&`g'
}

module load singularity
FSLContianer="singularity exec  /home/pvalsala/ptmp/MyContainers/fsl_ubuntu4_devel-2022-06-09-e06711f0139f.sif"

 #split the fullfile into path,filename and extension
 # extension can have two dot like .nii.gz 
pn=${fullfile%/*}
fn=${fullfile##*/}
ext=${fn#*.}
basefn=$(basename "$fn" ".$ext")


TR=$($FSLContianer fslval  "$fullfile" pixdim4)
res=$($FSLContianer fslval  "$fullfile" pixdim1)
smoothsize=$(echo "$smoothfac*$res" | bc)
nVol=$($FSLContianer fslval  "$fullfile" dim4)
Highpass=30 #seconds

smoothlabel=$(printf "_%0.1f_" $smoothsize)
#replace all . with p  in the smoothlabel
smoothlabel=${smoothlabel//./p}


#make output directory name with $(pwd)/smooth$smmothfac_$basefn
outputdir=$(pwd)/smooth$smoothlabel$basefn

#make design file name with $(pwd)/smooth$smmothfac_$basefn.fsf
designFile=$(pwd)/smooth$smoothlabel$basefn.fsf
cp $templatefile $designFile


# set fmri(tr) 1.960000 
sed -i "s/set fmri(tr) [0-9]*\.[0-9]*/set fmri(tr) $TR/g" $designFile
#set fmri(smooth) 0
sed -i "s/set fmri(smooth) [0-9]*/set fmri(smooth) $smoothsize/g" $designFile
# set fmri(npts) 168
sed -i "s/set fmri(npts) [0-9]*/set fmri(npts) $nVol/g" $designFile
# set fmri(paradigm_hp)
sed -i "s/set fmri(paradigm_hp) [0-9]*/set fmri(paradigm_hp) $Highpass/g" $designFile

# set fmri(custom1) "/ptmp/pvalsala/HRXA/EXPDATA/stim_p8.txt"
stimfile=$(sedConvert "$stimfile")
sed -i -E "s@^set\sfmri\(custom1\)\s*".*"@set fmri(custom1) \"$stimfile\"@g" $designFile
# set feat_files(1) "/ptmp/pvalsala/HRXA/moco/corr/moco_m68_B0MTI_DCFJackson"
featfile=$(sedConvert "$fullfile")
sed -i -E "s@^set\sfeat_files\(1\)\s*".*"@set feat_files(1) \"$featfile\"@g" $designFile
#set fmri(outputdir) "/ptmp/pvalsala/HRXA/blob/p8_smooth0"
outputdir=$(sedConvert "$outputdir")
sed -i -E "s@^set\sfmri\(outputdir\)\s*".*"@set fmri(outputdir) \"$outputdir\"@g" $designFile

#print all modified lines 
echo "Modified lines in $designFile"
grep -E "set feat_files(1)*" $designFile
grep -E "set fmri\(outputdir*" $designFile
grep -E "set fmri\(custom1*" $designFile
grep -E "set fmri\(tr\)*" $designFile
grep -E "set fmri\(smooth*" $designFile
grep -E "set fmri\(npts*" $designFile
grep -E "set fmri\(paradigm_hp*" $designFile

#run feat in background
echo "Running $FSLContianer feat $designFile &"
$FSLContianer feat $designFile &

