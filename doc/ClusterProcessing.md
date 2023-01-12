# List of steps

## Preparation
copy all data from irods to `/ptmp/pvalsala/<EXP_PSEUDONYM>/raw/`

## Get all csm for all the spiral data
As `bart` is not inside matlab container, we have to write `cfl` files of the calibration data and call `bart` container outside using the generated `runcsm2.sh` script.If we are using the reference scan data, there is no need for additional registration. 

Make sure the B0map volume has same center,Normal and inplane rotation for proper interpolation results. Ideally it should have same FOV in-plane and more FOV long partition direction. 
<details><summary>MATLAB scripts for coilmaps and fieldmap</summary>

```matlab
%% Create calib data, script for bart and registered fieldmap
addpath(genpath('/ptmp/pvalsala/MATLAB'))
Measpath='/ptmp/pvalsala/HRXA';
dir_st=dir(fullfile(Measpath,'raw','*peSpiral_R4*.dat'));
clear commands;
for i=1:length(dir_st)
    [~,measID]=regexp(dir_st(i).name,'\S*#M(\d{2,}+)\S*','match','tokens');
    fm_fn=findClosestFm(fullfile(Measpath,'raw'),str2double(measID{1}));
    ra=SpiralReco(fullfile(dir_st(i).folder,dir_st(i).name),'doCoilCombine','sos','RepSel',1,'CompMode','CPU3D');
    fmobj=B0map(fm_fn,'UnwrapMode','SpatialUnwrap');
    %only works if the orientation matches
    fmobj.PerformResampling(ra.twix);
    
    fmobj2=B0map(fm_fn,'UnwrapMode','UMPIRE');
    fmobj2.PerformResampling(ra.twix);
    outfile=fullfile(Measpath,'dep',sprintf('UMPIRE_FM_MeasUID%d.mat',str2double(measID{1})));
    fmobj.saveFmap(outfile); 
    
%   as(cat(5,squeeze(ra.img),fmobj.regIm,fmobj.Fmap_registered/(2*pi)));
    outfile=fullfile(Measpath,'dep',sprintf('fm_csm_MeasUID%d.mat',str2double(measID{1})));
    fmobj.saveFmap(outfile);  
    commands{i}=GetCoilMaps(ra,[],'Cluster');
    fprintf(' %s uses fm:  %s \n',measID{1}{:},fm_fn);
end
 save(fullfile(Measpath,'dep','allCommands.mat'),'commands')
fileID = fopen(fullfile(Measpath,'dep','runcsm2.sh'),'w');
fprintf(fileID,'module load singularity\n');
for i=1:length(commands)
fprintf(fileID,commands{i}{1}(4:end-1));
fprintf(fileID,'\n');
end
fclose(fileID);

%% load csm and put it back into the fm_csm_MeasUIDxxx.mat
load(fullfile(Measpath,'/dep/allCommands.mat'))
for i=1:length(dir_st)
    [~,measID]=regexp(dir_st(i).name,'\S*#M(\d{2,}+)\S*','match','tokens');  
    outfile=fullfile(Measpath,'dep',sprintf('fm_csm_MeasUID%d.mat',str2double(measID{1})));
    %load Coilmaps
    evalc(commands{i}{2});
    %select 1st map and save it to our fm_csm_MeasUIDxxx.mat
    csm=permute(Coilmaps(:,:,:,:,1),[4 2 1 3 5]);
    save(outfile,'csm','-append')
    
end

%% If the data is small you can do reconstruction here
for i=1:length(dir_st)
        [~,measID]=regexp(dir_st(i).name,'\S*#M(\d{2,}+)\S*','match','tokens');
        outfile=fullfile(Measpath,'dep',sprintf('fm_csm_MeasUID%d.mat',str2double(measID{1})));
load(outfile,'fm_interp','csm')

    ra=SpiralReco(fullfile(dir_st(i).folder,dir_st(i).name),'doCoilCombine','sos','RepSel',1,'CompMode','CPU3D');
rB0=SpiralReco(fullfile(dir_st(i).folder,dir_st(i).name),'RepSel',1,...
    'doPAT','CGSENSE','csm',csm,'maxit',10,'reg','Tikhonov','reg_lambda',1e-3,'compMode','CPU3D',...
    'fm',-1*fm_interp,'doB0Corr','MTI','doDCF','Jackson');
   r=SpiralReco(fullfile(dir_st(i).folder,dir_st(i).name),'RepSel',1,...
    'doPAT','CGSENSE','csm',csm,'maxit',10,'reg','Tikhonov','reg_lambda',1e-3,'compMode','CPU3D',...
    'fm',-1*fm_interp,'doB0Corr','none','doDCF','Jackson'); 
%   as(cat(5,squeeze(rB0.img),squeeze(r.img),squeeze(ra.img)));
    
end
%% supporting fucntions
function fm_fn=findClosestFm(folderpath,mid)
% simple function to get the closest fieldmap file from folder
% 
dir_st=dir(fullfile(folderpath,'*B0*.dat'));
[~,measID]=regexp({dir_st.name}','\S*#M(\d{2,}+)\S*','match','tokens');
measID=cellfun(@(x)str2double(x{:}),measID);
[~,idx]=min(abs(measID-mid));
fm_fn=fullfile(dir_st(idx).folder,dir_st(idx).name);
end
```
</details>


## Start reconstruction
Go the root folder of the measurement(`/ptmp/pvalsala/<EXP_PSEUDONYM>`) and Make sure all four files are in there. The following lines will give to niftis from raw data if the dependicies are already there. Cancel the whole job array and restart if some runtime exceptions happen in the begining. 

```bash
recojobid=$(sbatch runArrayjob.sh raw/*M<MEASUID>* <B0MODE{none,MFI,MTI}>)
recojobid2=$(sbatch runArrayjob.sh raw/*M<MEASUID>* <B0MODE{none,MFI,MTI}>)
sbatch --dependency=afterok:$recojobid:$recojobid2 runMergeJob.sh raw/*M68*
```

<details><summary>ArrayScript.m</summary>

```matlab
%%
addpath(genpath('/ptmp/pvalsala/MATLAB'))
filepattern=getenv('SPIRAL_FN');
b0mode=getenv('B0MODE');
fprintf("inoput file=%s\n",filepattern);
fprintf("B0mode=%s\n",b0mode);

if(~isfolder(fullfile(pwd,'proc')))
    mkdir(fullfile(pwd,'proc'))
end

%twix without extension
dirst=dir(filepattern);
 twix=mapVBVD(fullfile(dirst(1).folder,dirst(1).name));
 NVol=twix.image.NRep;
[pn,filename,~]=fileparts(fullfile(dirst(1).folder,dirst(1).name));
mid=cellfun(@(x)str2double(x{1}),(regexp(filename,'\S*#M(\d{2,}+)#\S*','tokens')));

protName=strsplit(filename,'#');
protName=sprintf('M%d_%s',mid,protName{end});
if(~isfolder(fullfile(pwd,'proc',protName)))
    mkdir(fullfile(pwd,'proc',protName))
end
cd(fullfile(pwd,'proc',protName))
% load csm,fm and start recon
load(sprintf('../../dep/fm_csm_MeasUID%d.mat',mid))
%% Spliting data for Array job
minTask=str2double(getenv('SLURM_ARRAY_TASK_MIN'));
maxTask=str2double(getenv('SLURM_ARRAY_TASK_MAX'));
nTask=str2double(getenv('SLURM_ARRAY_TASK_COUNT'));
cTask=str2double(getenv('SLURM_ARRAY_TASK_ID')); 
fprintf('\n minTASK : %d, ,maxTask: %d , nTask: %d, cTask : %d : reconstruction following volumes  \n',minTask,maxTask,nTask,cTask);
idx=floor(linspace(1,NVol+1,nTask+1));

%% do Reco
 fprintf('\n pwd is %s \n',pwd)
 fprintf(strcat(datestr(datetime),': Starting job \n'))
 fprintf('\n %s',filename)
 fprintf('\n nTask: %d, cTask : %d : reconstruction following volumes  \n',nTask,cTask);
 fprintf('%d, ',idx(cTask):(idx(cTask+1)-1));fprintf('\n ');
% %change filename, leave Repsel,1 as it is
r=SpiralReco(twix,'RepSel',idx(cTask):(idx(cTask+1)-1),...
    'doPAT','CGSENSE','csm',csm,'maxit',10,'reg','Tikhonov','reg_lambda',1e-3,'compMode','CPU2DHybrid',...
    'fm',-1*fm_interp,'doB0Corr',b0mode,'doDCF','Jackson','precision','double');

fprintf(strcat(datestr(datetime),': Finishing job \n'))
% save results
im=squeeze(r.img);
im=im(:,:,:,idx(cTask):(idx(cTask+1)-1));
    
flags=r.flags;
OutFile=sprintf('%s_B0%s_DCF%s_part%d.mat',protName,flags.doB0Corr,flags.doDCF,cTask);
if(cTask==1)
    sp=r.SpiralPara;
    fn=r.filename;
    ro=(2*sp.ADCLength*sp.DwellTime)/1e6; % ms
    vTR=(sp.TR*sp.Ninterleaves*sp.NPartitions)/(sp.R_PE*sp.R_3D*1e6); %s
    descrip=(sprintf('R%dx%dC%d TR=%.1fms RO=%.2fms vTR=%.1fs',sp.R_PE,sp.R_3D,sp.CAIPIShift,sp.TR/1e3,ro,vTR));
    descrip_reco=sprintf('%s PAT=%s coilcomb=%s B0=%s DCF=%s CompMode=%s',flags.CompMode,flags.doPAT, flags.doCoilCombine, flags.doB0Corr,flags.doDCF,flags.CompMode);   
    save(OutFile,'im','sp','flags','descrip','descrip_reco','fn','-v7.3')
else
    save(OutFile,'im','-v7.3')
end


```
</details>

<details><summary>runArrayjob.sh</summary>

```bash
#!/bin/bash -l
# A comment

#SBATCH -o ./job.out.%A_%a
#SBATCH -e ./job.err.%A_%a
#SBATCH -D ./
#SBATCH -J SpiralReco
# --- resource specification (which resources for how long) ---
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120000        # memory in MB required by the job
#SBATCH --time=24:00:00   # run time in h:m:s, up to 24h possible
#SBATCH --array=1-5
#SBATCH --requeue
# --- start from a clean state and load necessary environment modules ---
module purge
module load singularity

export OMP_NUM_THREADS=16
export SINGULARITYENV_OMP_NUM_THREADS=16
export SINGULARITYENV_SLURM_ARRAY_TASK_MIN=$SLURM_ARRAY_TASK_MIN
export SINGULARITYENV_SLURM_ARRAY_TASK_MAX=$SLURM_ARRAY_TASK_MAX
export SINGULARITYENV_SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID
export SINGULARITYENV_SLURM_ARRAY_TASK_COUNT=$SLURM_ARRAY_TASK_COUNT
export SINGULARITYENV_SPIRAL_FN=$1 
export SINGULARITYENV_B0MODE=$2

# put the command into a variable and run the matlab container with it
scriptloc=$(pwd)
command="cd $scriptloc; ArrayScript;"
srun singularity run -B /ptmp /ptmp/containers/matlab-r2020b.sif -nodisplay -batch "$command"

# move log files
mkdir -p logs/$SLURM_ARRAY_JOB_ID
mv job.*$SLURM_ARRAY_JOB_ID*$SLURM_ARRAY_TASK_ID logs/$SLURM_ARRAY_JOB_ID/

```
</details>

<details><summary>MergeScript.m</summary>

```matlab
%% Merge and make nifti
addpath(genpath('/ptmp/pvalsala/MATLAB'))
filepattern=getenv('SPIRAL_FN');
b0mode=getenv('B0MODE');

fprintf("inoput file= %s\n",filepattern);
fprintf("B0mode= %s\n",b0mode);


%twix without extension
dirst=dir(filepattern);
twix=mapVBVD(fullfile(dirst(1).folder,dirst(1).name));
NVol=twix.image.NRep;
[pn,filename,~]=fileparts(fullfile(dirst(1).folder,dirst(1).name));
MeasUID=cellfun(@(x)str2double(x{1}),(regexp(filename,'\S*#M(\d{2,}+)#\S*','tokens')));

protName=strsplit(filename,'#');
protName=sprintf('M%d_%s',MeasUID,protName{end});
if(~isfolder(fullfile(pwd,'proc',protName)))   
    error('Processed Data directoty not found: %s',fullfile(pwd,'proc',protName))
end
cd(fullfile(pwd,'proc',protName))

%%
dirst_part1=dir(fullfile(pwd,sprintf('%s*part1.mat',protName)));
for i=1:length(dirst_part1)
    FilePat=regexprep(dirst_part1(i).name,'part[0-9]+','*');
    dirst=dir(fullfile(pwd,FilePat));
    im_all=[];
    for i=1:length(dirst)
        fprintf('Loading file : %s \n',dirst(i).name)
        load(dirst(i).name)
        im_all=cat(4,im_all,im);
    end
    
    %some error checking
    if(NVol~=size(im_all,4))
        error('%d volumes expected only %d volumes are reconstructed',Nvol,size(im_all,4))
    end
    
    %do NORDIC
    im_all_Nordic=NORDIC_imComplex(im_all);
    outFile=sprintf('m%d_B0%s_merged_%s.mat',MeasUID,flags.doB0Corr,datetime('now','Format','d_MMM_y_HH_mm_ss'));
    save(outFile,'im_all','sp','flags','descrip_reco','descrip','fn','im_all_Nordic','-v7.3')
    delete(dirst.name)
    
    % write nifti file
    niiFile=sprintf('%s_B0%s_DCF%s',protName,flags.doB0Corr,flags.doDCF);
    MyNIFTIWriteSpiral(single(abs(im_all)),twix,niiFile);
    MyNIFTIWriteSpiral(single(abs(im_all_Nordic)),twix,strcat('NORDIC_',niiFile));
end

```
</details>

<details><summary>runMergeJob.sh</summary>

```matlab
#!/bin/bash -l
# A comment

#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
#SBATCH -D ./
#SBATCH -J MergeScript
# --- resource specification (which resources for how long) ---
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400000        # memory in MB required by the job
#SBATCH --time=24:00:00   # run time in h:m:s, up to 24h possible
# --- start from a clean state and load necessary environment modules ---
module purge
module load singularity

export OMP_NUM_THREADS=64
export SINGULARITYENV_OMP_NUM_THREADS=64
export SINGULARITYENV_SLURM_ARRAY_TASK_MIN=$SLURM_ARRAY_TASK_MIN
export SINGULARITYENV_SLURM_ARRAY_TASK_MAX=$SLURM_ARRAY_TASK_MAX
export SINGULARITYENV_SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID
export SINGULARITYENV_SLURM_ARRAY_TASK_COUNT=$SLURM_ARRAY_TASK_COUNT
export SINGULARITYENV_SPIRAL_FN=$1 
export SINGULARITYENV_B0MODE=$2

# put the command into a variable and run the matlab container with it
scriptloc=$(pwd)
command="cd $scriptloc; MergeScript;"
srun singularity run -B /ptmp /ptmp/containers/matlab-r2020b.sif -nodisplay -batch "$command"

# move log files
mkdir -p logs/$SLURM_JOB_ID
mv job.*$SLURM_JOB_ID* logs/$SLURM_JOB_ID/

```
</details>

