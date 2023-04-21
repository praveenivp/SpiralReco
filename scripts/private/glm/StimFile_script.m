%%
op_dirst=dir('/ptmp/pvalsala/JIJP-YO7X/proc/M*');
ip_dirst=dir('/ptmp/pvalsala/JIJP-YO7X/EXPDATA/stim/fmri*.txt');
if(length(op_dirst)~=length(ip_dirst)); error('ip files and out put doesn''t match');  end
for i=1:length(ip_dirst)
dat=textscan(fopen(fullfile(ip_dirst(i).folder,ip_dirst(i).name)),'%d %f %d','HeaderLines',11);
 vTR=mean(diff(dat{2}-dat{2}(1))); %s
 stim= [dat{2}(diff([0; dat{3}])>0)-dat{2}(1) vTR+dat{2}(diff([dat{3};0])<0)-dat{2}(diff([0; dat{3}])>0)];
 stim=padarray(stim,[0 1],1,'post');
 fprintf(fopen([op_dirst(i).name,'.txt'],'w+'),'%.4f %.4f %d\n',stim');
end
