 Subjects={'JEL7-26IB','DKSR-UFYK','H734-SPCU','JIJP-YO7X','YU3S-VKP3','FX6C-SZ37'};
% Subjects={'H734-SPCU'};

%% check the order and filter
clearvars -except Subjects
for cSub=1:length(Subjects)
        pn=sprintf('/ptmp/pvalsala/%s/moco/allmoco',Subjects{cSub});
        pn2=sprintf('/ptmp/pvalsala/%s/glm',Subjects{cSub});
        pn3=sprintf('/ptmp/pvalsala/%s/reg',Subjects{cSub});
%  %TE     
%  im_dir_TE= dir(fullfile(pn,'moco*TR10*.nii*'));
%  feat_dir_TE=dir(fullfile(pn2,'smooth_1p0_moco_M*TR10*.feat'));
% %TR
% im_dir_TR= cat(1,dir(fullfile(pn,'moco*TR6*.nii*')),dir(fullfile(pn,'moco*TR8*.nii*')),dir(fullfile(pn,'moco*TE2_TR10*.nii*')));
% feat_dir_TR=cat(1,dir(fullfile(pn2,'smooth_1p0_moco_M*TR6*.feat')),dir(fullfile(pn2,'smooth_1p0_moco_M*TR8*.feat')),dir(fullfile(pn2,'smooth_1p0_moco_M*TE2_TR10*.feat')));
% % high res
im_dir_HR= cat(1,dir(fullfile(pn,'moco*p8*.nii*')),dir(fullfile(pn,'moco*p6*.nii*')));
feat_dir_HR=cat(1,dir(fullfile(pn2,'smooth_0p8_moco_M*p8*.feat')),dir(fullfile(pn2,'smooth_0p6_moco_M*p6*.feat')));
reg_dir_HR=cat(1,dir(fullfile(pn3,'moco_M*p8*')),dir(fullfile(pn3,'moco_M*p6*')));
% 

fprintf('TE: %s\n',Subjects{cSub})
disp({im_dir_HR.name})
disp({feat_dir_HR.name})
fprintf('\n\n')

im_dir_all{cSub}=[im_dir_HR];
feat_dir_all{cSub}=[feat_dir_HR];
reg_dir_all{cSub}=[reg_dir_HR];
end

% clearvars -except im_dir_all feat_dir_all Subjects
%%
clc
% SlcSel=3:16; % 2+2 slice removed
mypermute=@(imx) permute(imx,[2 1 3 4]);
clear s1 s2 im1  im volTR im_info sig_chang1 sig_chang2 T1 ribbon brain
for cSub=1:length(Subjects)
    fprintf('TE: %s\n',Subjects{cSub})

for i=1:length(im_dir_HR)
    

    

 im_all=niftiread(fullfile(im_dir_all{cSub}(i).folder,im_dir_all{cSub}(i).name));
 im_info{i,cSub}=niftiinfo(fullfile(im_dir_all{cSub}(i).folder,im_dir_all{cSub}(i).name));

volTR{i,cSub}=im_info{i,cSub}.PixelDimensions(4);
 im1{i,cSub}=mypermute(mean(im_all(:,:,:,5,1),4));
im{i,cSub}=mypermute(mean(im_all(:,:,:,5:60,1),4));
im_ts{i,cSub}=mypermute(im_all);

reg_folder=fullfile(reg_dir_all{cSub}(i).folder,reg_dir_all{cSub}(i).name);
T1_file=dir(fullfile(reg_folder,'coregT1*.nii'));
T1{i,cSub}=mypermute( niftiread(fullfile(reg_folder,T1_file(1).name)));
brain_file=dir(fullfile(reg_folder,'coregBM*.nii'));
brain{i,cSub}=mypermute( niftiread(fullfile(reg_folder,brain_file(1).name)));
ribbon_file=dir(fullfile(reg_folder,'coregR*.nii'));
ribbon{i,cSub}=mypermute(niftiread(fullfile(reg_folder,ribbon_file(1).name)));

s1{i,cSub}=mypermute(niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'rendered_thresh_zstat1.nii.gz')));
s2{i,cSub}=mypermute(niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'rendered_thresh_zstat2.nii.gz')));

% spatial signal change
imf=niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'mean_func.nii.gz'));
beta1=niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'stats/pe1.nii.gz'));
beta2=niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'stats/pe2.nii.gz'));
sig_chang1{i,cSub}=mypermute(100*(beta1./imf(:,:,:)));
 sig_chang2{i}=mypermute(100*(beta2./imf(:,:,:))); %intercept
end
end

save 'dataHR_smooth.mat' s1 s2 im im1 volTR im_info sig_chang1 im_dir_all feat_dir_all Subjects 

%% get feieldmaps
clear fm_all
clc
for cds=1:size(im_ts,1) % dataset
    for cSub=1:size(im_ts,2)
        
        measID= regexp(im_dir_all{cSub}(cds).name,'moco_M([0-9]+).*','tokens');
        measID=str2double(measID{1});
        mat_file=sprintf('fm_csm_MeasUID%d.mat',measID);
        fprintf ('%d | loading %s \n',measID,mat_file)
        pn=sprintf('/ptmp/pvalsala/%s/dep',Subjects{cSub});
        load(fullfile(pn,mat_file),'fm_interp');
        fm_all{cds,cSub}=fm_interp;
    end
    fprintf('\n');
end



%%  TSNR B0
clc
clear TSNR_mean TSNR_std B0_hz_mean B0_hz_std
for cds=1:size(im_ts,1)
    for cSub=1:size(im_ts,2)
        cTSNR=mean(im_ts{cds,cSub}(:,:,:,10:end),4)./std(im_ts{cds,cSub}(:,:,:,10:end),[],4);
        cMask=brain{cds,cSub}>0;
        cMask=imerode(cMask,strel('sphere',5));
        % set 20 % of slice on each size to zezro
        cMask(:,:,1:round(0.1*size(cMask,3)))=0;
        cMask(:,:,end-round(0.1*size(cMask,3)):end)=0;
        
        TSNR_mean(cds,cSub)=mean(cTSNR(cMask),'omitnan');
        TSNR_std(cds,cSub)=std(cTSNR(cMask),'omitnan');
        B0_hz_mean(cds,cSub)=mean(fm_all{cds,cSub}(cMask),'omitnan')./(2*pi);
        B0_hz_std(cds,cSub)=std(fm_all{cds,cSub}(cMask),'omitnan')./(2*pi);
        
        fprintf ('%s | mean TSNR %.2f | std TSNR : %.2f |',Subjects{cSub},TSNR_mean(cds,cSub),TSNR_std(cds,cSub))
        fprintf ('mean B0 %.2f | std B0 %.2f Hz | \n',B0_hz_mean(cds,cSub),B0_hz_std(cds,cSub))
    end
    fprintf('\n');
end

t=table(TSNR_mean,TSNR_std,B0_hz_mean,B0_hz_std);
writetable(t,'TNSR_B0.xlsx')

%% fmri stats
clc
clear nvoxels_all zstat_median sigch_median zstat_prc95 sigch_prc95

for cds=1:size(im_ts,1)
    for cSub=1:size(im_ts,2)
        act_mask=s1{cds,cSub}>4;
% %         
%                 act_mask(:,:,1:round(0.1*size(act_mask,3)))=0;
%         act_mask(:,:,end-round(0.1*size(act_mask,3)):end)=0;
        
        nvoxels_all(cds,cSub)=sum(act_mask,'all');
        zstat_median(cds,cSub)=median(s1{cds,cSub}(act_mask));
        sigch_median(cds,cSub)=median(sig_chang1{cds,cSub}(act_mask));
        
        zstat_prc95(cds,cSub)=prctile(s1{cds,cSub}(act_mask),95);
        sigch_prc95(cds,cSub)=prctile(sig_chang1{cds,cSub}(act_mask),95);
        
        
        fprintf ('%s |  %d  voxels| zstat median %.2f | sig_change %.2f |%.2f |%.2f | \n',Subjects{cSub},nvoxels_all(cds,cSub),zstat_median(cds,cSub),sigch_median(cds,cSub), zstat_prc95(cds,cSub),sigch_prc95(cds,cSub))

    end
    fprintf('\n');
end
t=table(TSNR_mean,TSNR_std,B0_hz_mean,B0_hz_std,nvoxels_all, zstat_median, sigch_median, zstat_prc95, sigch_prc95);
writetable(t,'TNSR_B0_fmri_stats.xlsx')

% save allstats.mat nvoxels_all zstat_median sigch_median zstat_prc95 sigch_prc95 TSNR_mean TSNR_std B0_hz_mean B0_hz_std


%% better table
% t=table(TSNR_mean,TSNR_std,B0_hz_mean,B0_hz_std,nvoxels_all, , , , );
header={'measures','S1','S2','S3','S4','S5','S6'};
t3=table();
nm=10;
measures={'mean tSNR','SD tSNR','mean B0 [Hz]','SD B0 [Hz]','activated voxels (Z>4)','median Z score','P_95% Zscore','median \Delta S' ,'P_95% \Delta S' };

for cds=1:2
t3(1 +nm*(cds-1),:)=[measures{1},num2cell(TSNR_mean(cds,:))];
t3(2 +nm*(cds-1),:)=[measures{2},num2cell(TSNR_std(cds,:))];
t3(3 +nm*(cds-1),:)=[measures{3},num2cell(B0_hz_mean(cds,:))];
t3(4 +nm*(cds-1),:)=[measures{4},num2cell(B0_hz_std(cds,:))];
t3(5 +nm*(cds-1),:)=[measures{5},num2cell(nvoxels_all(cds,:))];
t3(6 +nm*(cds-1),:)=[measures{6},num2cell(zstat_median(cds,:))];
t3(7 +nm*(cds-1),:)=[measures{7},num2cell(zstat_prc95(cds,:))];
t3(8 +nm*(cds-1),:)=[measures{8},num2cell(sigch_median(cds,:))];
t3(9 +nm*(cds-1),:)=[measures{9},num2cell(sigch_prc95(cds,:))];
    
end

t3.Properties.VariableNames=string(header);
writetable(t3,'TNSR_B0_fmri_stats_better.xlsx')
%% signal change if 
Thres_P=4;Thres_N=4;
clear nVoxels Zscores sig_change_vec Contrasts_idx;
% TEs=[];
% Zscores=[];
% Contrasts_idx=[];
% sig_change_vec=[];
SlcSel1={5:24,3:16};
clc
for cSub=1:length(Subjects)
    fprintf('TR: %s\n',Subjects{cSub})

for i=1:length(im_dir_all{1})
    SlcSel=SlcSel1{i};
Pmask=(s1{i,cSub}>Thres_P);
Pmask(:,:,1:(SlcSel(1)-1))=0;Pmask(:,:,(SlcSel(end)+1):end)=0;Pmask(isnan(Pmask))=0;
nVoxels{i,cSub,1}=sum(Pmask,'all');
% TEs{i,cSub,1} =[TE_list{i}*ones(nVoxels{i,cSub,1},1)];
Zscores{i,cSub,1}=[s1{i,cSub}(Pmask)];
sig_change_vec{i,cSub,1}=[sig_chang1{i,cSub}(Pmask) ];
Contrasts_idx{i,cSub,1}=[ones(nVoxels{i,cSub,1},1)];

Nmask=(s2{i}>Thres_N);
Nmask(:,:,1:(SlcSel(1)-1))=0;Nmask(:,:,(SlcSel(end)+1):end)=0;Nmask(isnan(Nmask))=0;
nVoxels{i,cSub,2}=sum(Nmask,'all');
Zscores{i,cSub,2}=[ s2{i}(Nmask)];
sig_change_vec{i,cSub,2}=[sig_chang1{i,cSub}(Nmask) ];
Contrasts_idx{i,cSub,2}=[ 2*ones(nVoxels{i,cSub,1},1)];

end

end

%%
function m=getmask(x)
cutoff=min([5000 length(x)]);
[xs,idx]=sort(x(:),1,'descend');
disp(xs(1))
m=zeros(size(x));
m(idx(1:cutoff))=1;
m=(m>0);
end


function [cim,Ic1,crop_fun,composed]=cropim(im,Other,mytransform,SlcFrac,padSize)

if(~exist('SlcFrac','var'))
SlcFrac=[0.5, 0.8 0.7 0.3];
end
if(~exist('padSize','var'))
padSize=[0 5];
end
im=mytransform(im);
thres=0.1;
   [idx1]=find(squeeze(mean(im,[2,3])) > thres*max(squeeze(mean(im,[2,3]))),1,'first');
   [idx2]=find(squeeze(mean(im,[2,3])) > thres*max(squeeze(mean(im,[2,3]))),1,'last');
   [idy1]=find(squeeze(mean(im,[1,3])) > thres*max(squeeze(mean(im,[1,3]))),1,'first');
   [idy2]=find(squeeze(mean(im,[1,3])) > thres*max(squeeze(mean(im,[1,3]))),1,'last');
   [idz1]=find(squeeze(mean(im,[1,2])) > thres*max(squeeze(mean(im,[1,2]))),1,'first');
   [idz2]=find(squeeze(mean(im,[1,2])) > thres*max(squeeze(mean(im,[1,2]))),1,'last');
   disp([idx1,idx2,idy1,idy2,idz1,idz2])
   crop_fun=@(im3D)im3D(idx1:idx2,idy1:idy2,idz1:idz2);
   
   
   cim=crop_fun(im);
   ACS= floor([size(cim,3) size(cim,1) size(cim,2) size(cim,2)].*SlcFrac);
   
   Ic1={cim(:,:,ACS(1)),squeeze(cim(ACS(2),:,:))',[],squeeze(cim(:,ACS(3),:)),squeeze(cim(:,ACS(4),:))};
   width=floor(size(Ic1{2},2)/2);
   cfac=10; % cut more in the middle
   Ic1{3}=cat(1,Ic1{4}((end-width-cfac):(end-cfac),:),flip(Ic1{5}(end-width-cfac:end-cfac,:),1));
   Ic1{3}=Ic1{3}(1:size(Ic1{2},2),:)';
   
   if(padSize(1)<0)
   pz1=padSize(1);
   padSize(1)=0;
   Ic1{1}=Ic1{1}(1:end+pz1,:);
   else 
       pz1=0;
   end
   
   composed{1}=cat(1,Ic1{1},zeros(padSize(1),size(Ic1{1},2)),Ic1{2},zeros(padSize(2),size(Ic1{1},2)),Ic1{3});
   for i=1:length(Other)
       
       tmp=crop_fun(mytransform(Other{i}));
       tmp2={tmp(:,:,ACS(1)),squeeze(tmp(ACS(2),:,:))',[],squeeze(tmp(:,ACS(3),:)),squeeze(tmp(:,ACS(4),:))};
          tmp2{3}=cat(1,tmp2{4}((end-width-cfac):(end-cfac),:),flip(tmp2{5}(end-width-cfac:end-cfac,:),1));
   tmp2{3}=tmp2{3}(1:size(tmp2{2},2),:)';
   tmp2{1}=tmp2{1}(1:end+pz1,:);
   composed{i+1}=cat(1,tmp2{1},zeros(padSize(1),size(Ic1{1},2)),tmp2{2},zeros(padSize(2),size(tmp2{1},2)),tmp2{3});

   
   end
   
end

function plotString(cim,cax)
if(~exist('cax','var'))
cax=gca;
end
text(cax,size(cim,2)/2,size(cim,1)+12,'Coronal',...
'Interpreter','latex','FontSize',18,'FontWeight','bold','Color',[1 1 1],...
    'EdgeColor','none','HorizontalAlignment','center');
text(cax,size(cim,2)/2,size(cim,1)+60,'L~~~~Sagital~~~~R',...
'Interpreter','latex','FontSize',18,'FontWeight','bold','Color',[1 1 1],...
    'EdgeColor','none','HorizontalAlignment','center');
set(gcf,'Position',[429 149 1228 860])
end

function plotStringp6(cim,cax)
if(~exist('cax','var'))
cax=gca;
end
text(cax,size(cim,2)/2,size(cim,1)+12,'Coronal',...
'Interpreter','latex','FontSize',18,'FontWeight','bold','Color',[1 1 1],...
    'EdgeColor','none','HorizontalAlignment','center');
text(cax,size(cim,2)/2,size(cim,1)+55,'L~~~~Sagital~~~~R',...
'Interpreter','latex','FontSize',18,'FontWeight','bold','Color',[1 1 1],...
    'EdgeColor','none','HorizontalAlignment','center');
set(gcf,'Position',[429 149 1228 860])
end
