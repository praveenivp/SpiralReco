%% generate paradigm file with the logfile from psychopy
% Only works if the folders in /proc and logfiles in /EXP_DATA/stim matches
data_dirst=dir('../proc/M*R4_*');
stimulus_dirst=dir('../EXPDATA/stim/fmri*.txt');
mkdir('../EXPDATA/paradigm')
if(length(data_dirst)~=length(stimulus_dirst)); error('ip files and out put doesn''t match');  end
for i=1:length(stimulus_dirst)
dat=textscan(fopen(fullfile(stimulus_dirst(i).folder,stimulus_dirst(i).name)),'%d %f %d','HeaderLines',11);
 vTR=mean(diff(dat{2}-dat{2}(1))); %s
 stim= [dat{2}(diff([0; dat{3}])>0)-dat{2}(1) vTR+dat{2}(diff([dat{3};0])<0)-dat{2}(diff([0; dat{3}])>0)];
 stim=padarray(stim,[0 1],1,'post');
 fprintf(fopen(['../EXPDATA/paradigm/',data_dirst(i).name,'.txt'],'w+'),'%.4f %.4f %d\n',stim');
end