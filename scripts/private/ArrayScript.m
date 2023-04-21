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
OutFile=sprintf('%s_B0%s_DCF%s_part%02d.mat',protName,flags.doB0Corr,flags.doDCF,cTask);
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
