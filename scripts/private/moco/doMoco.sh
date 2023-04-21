#!/bin/sh
# usage : ./doMocoWithBaseVol Vol1.nii

module load singularity

for filename in *.nii; do
    echo "working on "  ${filename}
    fn=$(echo "$filename" | cut -f 1 -d '.')
    fullFilename=$(pwd)/$filename
    FileName1D=$(pwd)/1D_$fn
    OUTfile=moco_$fn

    echo "OUTfile=" $OUTfile
    # motion correction time series

    if [ $# -eq 0 ]
    then 
 echo "using the first volume of the time series as base\n"

    singularity run -B $pwd:/home/afni_user/work  \
    /home/pvalsala/ptmp/MyContainers/afni_afni_make_build_latest-2022-09-02-223e89c19492.sif \
    3dvolreg -verbose -zpad 1 -prefix $OUTfile -1Dfile $FileName1D -cubic $fullFilename
else
    singularity run -B $pwd:/home/afni_user/work  \
    /home/pvalsala/ptmp/MyContainers/afni_afni_make_build_latest-2022-09-02-223e89c19492.sif \
    3dvolreg -verbose -zpad 1 -prefix $OUTfile -1Dfile $FileName1D -base $1 -cubic $fullFilename
fi



    # convert AFNI files to NIFTI
    BIRKFile=$OUTfile+orig;
    singularity run  -B $pwd:/home/afni_user/work  \
    /home/pvalsala/ptmp/MyContainers/afni_afni_make_build_latest-2022-09-02-223e89c19492.sif \
    3dAFNItoNIFTI $BIRKFile
    # plot the motion paramters for 
    singularity run  -B $pwd:/home/afni_user/work  \
    /home/pvalsala/ptmp/MyContainers/afni_afni_make_build_latest-2022-09-02-223e89c19492.sif \
     1dplot -volreg   -jpgs 1024  $FileName1D $FileName1D

done
echo "Cleaning AFNI files"
rm *.BRIK
rm *.HEAD
