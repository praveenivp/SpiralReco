%% Create calib data, script for bart and registered fieldmap
addpath(genpath('/ptmp/pvalsala/MATLAB'))
Measpath='/ptmp/pvalsala/YU3S-VKP3';
dir_st=dir(fullfile(Measpath,'TWIX','*peSpiral_R4*.dat'));
clear commands;
for i=1:length(dir_st)
    [~,measID]=regexp(dir_st(i).name,'\S*#M(\d{2,}+)\S*','match','tokens');
    fm_fn=findClosestFm(fullfile(Measpath,'TWIX'),str2double(measID{1}));
    ra=SpiralReco(fullfile(dir_st(i).folder,dir_st(i).name),'doCoilCombine','sos','RepSel',1,'CompMode','CPU3D');
    fmobj=B0map(fm_fn,'UnwrapMode','SpatialUnwrap');
    %only works if the orientation matches
    fmobj.PerformResampling(ra.twix);
    
    fmobj2=B0map(fm_fn,'UnwrapMode','UMPIRE');
    fmobj2.PerformResampling(ra.twix);
    outfile=fullfile(Measpath,'dep',sprintf('UMPIRE_FM_MeasUID%d.mat',str2double(measID{1})));
    fmobj2.saveFmap(outfile); 
    
%   as(cat(5,squeeze(ra.img),fmobj.regIm,fmobj.Fmap_registered/(2*pi)));
    outfile=fullfile(Measpath,'dep',sprintf('fm_csm_MeasUID%d.mat',str2double(measID{1})));
    fmobj.saveFmap(outfile);  
    commands{i}=GetCoilMaps(ra,[],'Cluster');
    fprintf(' %s uses fm:  %s \n',measID{1}{:},fm_fn);
end
 save(fullfile(Measpath,'dep','allCommands.mat'),'commands')
fileID = fopen(fullfile(Measpath,'dep','runcsm2.sh'),'w');
fileattrib(fullfile(Measpath,'dep','runcsm2.sh'),'+x'); 
fprintf(fileID,'module load singularity\n');
for i=1:length(commands)
fprintf(fileID,commands{i}{1}(4:end-1));
fprintf(fileID,' & \n');
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

%% after running the BART
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