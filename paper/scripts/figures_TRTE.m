%% set data paths
Subjects={'JEL7-26IB','DKSR-UFYK','H734-SPCU','JIJP-YO7X','YU3S-VKP3','FX6C-SZ37'};
for cSub=1:length(Subjects)

        pn=sprintf('/ptmp/pvalsala/%s/moco/allmoco',Subjects{cSub});
        pn2=sprintf('/ptmp/pvalsala/%s/glm',Subjects{cSub});
 %TE     
 im_dir_TE= dir(fullfile(pn,'moco*TR10*.nii*'));
 feat_dir_TE=dir(fullfile(pn2,'smooth_1p0_moco_M*TR10*.feat'));
%TR
im_dir_TR= cat(1,dir(fullfile(pn,'moco*TR6*.nii*')),dir(fullfile(pn,'moco*TR8*.nii*')),dir(fullfile(pn,'moco*TE2_TR10*.nii*')));
feat_dir_TR=cat(1,dir(fullfile(pn2,'smooth_1p0_moco_M*TR6*.feat')),dir(fullfile(pn2,'smooth_1p0_moco_M*TR8*.feat')),dir(fullfile(pn2,'smooth_1p0_moco_M*TE2_TR10*.feat')));
% high res
im_dir_HR= cat(1,dir(fullfile(pn,'moco*p8*.nii*')),dir(fullfile(pn,'moco*p6*.nii*')));
feat_dir_HR=cat(1,dir(fullfile(pn2,'smooth_0p0_moco_M*p8*.feat')),dir(fullfile(pn2,'smooth_0p0_moco_M*p6*.feat')));
fprintf('TE: %s\n',Subjects{cSub})
disp({im_dir_TR.name})
disp({feat_dir_TR.name})
fprintf('\n\n')

im_dir_all{cSub}=[im_dir_TE;im_dir_TR;];
feat_dir_all{cSub}=[feat_dir_TE;feat_dir_TR;];
end

clearvars -except im_dir_all feat_dir_all Subjects
%% load data
clc
load('/ptmp/pvalsala/paper/tSNR_new/AllMask.mat')
SlcSel=3:18; % 2+2 slice removed
mypermute=@(imx) permute(imx,[2 1 3 4]);
% clear s1 s2 im volTR im_info sig_chang1 sig_chang2 TSNR_mean mask_all im_ts
for cSub=1:length(Subjects)
    fprintf('TE/TR: %s\n',Subjects{cSub})

for i=1:length(im_dir_all{1})
    

 im_all=niftiread(fullfile(im_dir_all{cSub}(i).folder,im_dir_all{cSub}(i).name));

im_info{i,cSub}=niftiinfo(fullfile(im_dir_all{cSub}(i).folder,im_dir_all{cSub}(i).name));
 im_ts{i,cSub}=im_all;
 volTR{i,cSub}=im_info{i}.PixelDimensions(4);
 im{i,cSub}=mypermute(mean(im_all(:,:,:,5:end,1),4));
% cMask= im{i,cSub}> prctile(col(im{i,cSub}),50);
% cMask(:,:,[1 2 19 20])=0;
% cMask=imerode(cMask,strel('sphere',5));
% mask_all{i,cSub}= cMask;

cMask=AllMask_cut{1,cSub};
 cMask=imerode(cMask,strel('sphere',3));
cTSNR=mean(im_ts{i,cSub}(:,:,:,5:end),4)./std(im_ts{i,cSub}(:,:,:,5:end),[],4);
TSNR_mean(i,cSub)=mean(cTSNR(cMask),'omitnan');
TSNR_std(i,cSub)=std(cTSNR(cMask),'omitnan');

s1{i,cSub}=mypermute(niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'rendered_thresh_zstat1.nii.gz')));
s2{i,cSub}=mypermute(niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'rendered_thresh_zstat2.nii.gz')));

% spatial signal change
imf=niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'filtered_func_data.nii.gz'));
beta1=niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'stats/pe1.nii.gz'));
beta2=niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'stats/pe2.nii.gz'));
 imf1=mean(imf(:,:,:,10:end),4);
 sig_chang1{i,cSub}=mypermute(100*(beta1./imf1(:,:,:)));
  sig_chang2{i}=mypermute(100*(beta2./imf1(:,:,:))); %intercept
end
end

save 'data_smooth2.mat' s1 s2 im volTR im_info sig_chang1 im_dir_all feat_dir_all Subjects 

%% signal change calc
Thres_P=4;Thres_N=4;
clear nVoxels Zscores sig_change_vec Contrasts_idx;
TR_list={6,8,10};
clc
for cSub=1:length(Subjects)
    fprintf('TR: %s\n',Subjects{cSub})

for i=1:length(im_dir_all{1})
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
% TEs{i,cSub,2} =[TE_list{i}*ones(nVoxels{i,cSub,1},1)];
Zscores{i,cSub,2}=[ s2{i}(Nmask)];
sig_change_vec{i,cSub,2}=[sig_chang1{i,cSub}(Nmask) ];
Contrasts_idx{i,cSub,2}=[ 2*ones(nVoxels{i,cSub,1},1)];

end

end


%% TR plot : figure 4
clc
plotSLC=[6,8,12,  10,6,8];
Thres_P=[4 20];Thres_N=[4 8];
figure(12),clf,
t = tiledlayout(4,4,'TileSpacing','compact');
% TR_list={2,4,6};
TR_list={6,8,10};
crop_RL=20:180;
for i=1:3
% subplot(4,1,i*1)
nexttile([1 4])

im_plot=[];
s1_plot=[];
s2_plot=[];
mask_plot=[];
 title(sprintf('Activation maps : TE %d ms',TR_list{i})...
     ,'Interpreter','latex','FontSize',16,'FontWeight','bold','Color',[1 1 1])
for cSub=1:6
im_plot=cat(3,im_plot,im{i+3,cSub}(:,crop_RL,plotSLC(cSub)));
mask_plot=cat(2,mask_plot,AllMask_cut{1,cSub}(:,crop_RL,plotSLC(cSub)));
im_plot(:,:,cSub)=im_plot(:,:,cSub)./prctile(col(im_plot(:,:,cSub)),90);
s1_plot=cat(3,s1_plot,s1{i+3,cSub}(:,crop_RL,plotSLC(cSub)));
s2_plot=cat(3,s2_plot,s2{i+3,cSub}(:,crop_RL,plotSLC(cSub)));
end

clim_im=[0,2];
cmap_P=hot(4096);
[allAx]=makeBlobPlot(flip(abs(im_plot),1),flip(s1_plot,1),'Thres',Thres_P,'SlcSel',1:6,'im_horz',6, ...
    'title_im',sprintf('Activation maps : TR %d ms',TR_list{i}),'colorbar',i==3,'caxis_im',clim_im, ...
    'alpha_blobs',0.8);

 alpha_blobs=0.9;
 cmap_N=flip(cmap_P,2);

makeBlobPlot([],flip(s2_plot,1),'Thres',Thres_N,'SlcSel',1:6,'im_horz',6,'title_im',sprintf('Activation maps : TR %d ms',TR_list{i}),...
    'negBold',true,'cmap',cmap_N,'colorbar',i==3,'allAx',allAx,'alpha_blobs',0.8)

% hold on
% contour(flip(mask_plot,1),'color','blue')
end


% summary stats
TR_sel=4:6;
MSC=cellfun(@(x,y)median(x(getmask(y))),sig_change_vec(TR_sel,:,1),Zscores(TR_sel,:,1));
MZ=cellfun(@(x,y)median(x(getmask(y))),Zscores(TR_sel,:,1),Zscores(TR_sel,:,1));
AV=cell2mat(nVoxels(TR_sel,:,1));
TSNR_TR=TSNR_mean(TR_sel,:)./sqrt([volTR{TR_sel,1}])';

% figure,
nexttile(t,[1 1]),scatter([TR_list{:}],MSC,100,'filled','LineWidth',1,'MarkerFaceAlpha',0.9),xlabel('TR [ms]'),title('median sig. change [%]','Color','w'),
hold on,plot([TR_list{:}],mean(MSC,2),'x-','LineWidth',3,'Markersize',15,'Color',[0.6350 0.0780 0.1840])
set(gca,'Box','on','YColor',[1 1 1],'Xcolor',[1 1 1],'Color',[0 0 0],'fontweight','bold','FontSize',10);
xticks([6 8 10]),xlim([5.5 10.5]),ylim([2 5]),grid on

nexttile,scatter([TR_list{:}],MZ,100,'filled','LineWidth',1,'MarkerFaceAlpha',1),xlabel('TR [ms]'),title('median z-score','Color','w')
hold on,plot([TR_list{:}],mean(MZ,2),'x-','LineWidth',3,'Markersize',15,'Color',[0.6350 0.0780 0.1840])
set(gca,'Box','on','YColor',[1 1 1],'Xcolor',[1 1 1],'Color',[0 0 0],'fontweight','bold','FontSize',10);
xticks([6 8 10]),xlim([5.5 10.5]),ylim([5 13]),grid on

nexttile,scatter([TR_list{:}],AV./1000,100,'filled','LineWidth',1,'MarkerFaceAlpha',1),xlabel('TR [ms]'),
title('activated voxels x1000','Color','w')

hold on,plot([TR_list{:}],mean(AV./1000,2),'x-','LineWidth',3,'Markersize',15,'Color',[0.6350 0.0780 0.1840])
set(gca,'Box','on','YColor',[1 1 1],'Xcolor',[1 1 1],'Color',[0 0 0],'fontweight','bold','FontSize',10);
xticks([6 8 10]),xlim([5.5 10.5]),ylim([0 18]),grid on



nexttile,scatter([TR_list{:}],TSNR_TR,100,'filled','LineWidth',1,'MarkerFaceAlpha',1),xlabel('TR [ms]'),title('mean tSNR eff. [1/√s]','Color','w')
hold on,plot([TR_list{:}],mean(TSNR_TR,2),'x-','LineWidth',3,'Markersize',15,'Color',[0.6350 0.0780 0.1840])
set(gca,'Box','on','YColor',[1 1 1],'Xcolor',[1 1 1],'Color',[0 0 0],'fontweight','bold','FontSize',10);
xticks([6 8 10]),xlim([5.5 10.5]),ylim([11 21]),grid on 

legend({'sub1','sub2','sub3','sub4','sub5','sub6','mean'},'textcolor','w','Location','eastoutside') 

% title('$Z_{scores}$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color',[1 1 1])
set(gcf,'Position',[[293 112 1100 1000]])

%%
savefig(gcf,'TR_dependece_allsub_Smooth_masked')
FigH=gcf;
FigH.InvertHardcopy = 'off';
FigH.Color='black';
print(FigH,'TR_dependece_allsub_Smooth_masked.tiff','-dtiff','-r600')

%% TE plot : figure 3
set(0,'defaulttextinterpreter','tex')
clc
plotSLC=[6,8,12,10,6,8];
Thres_P=[4 16];Thres_N=[4 8];
figure(13),clf,
t = tiledlayout(4,4,'TileSpacing','compact');
 TE_list={2,4,6};
for i=1:3
% subplot(4,1,i*1)
nexttile([1 4])

im_plot=[];
s1_plot=[];
s2_plot=[];
mask_plot=[];
 title(sprintf('Activation maps : TE %d ms',TE_list{i})...
     ,'Interpreter','none','FontSize',14,'FontWeight','bold','Color',[1 1 1])
for cSub=1:6
im_plot=cat(3,im_plot,im{i,cSub}(:,crop_RL,plotSLC(cSub)));
im_plot(:,:,cSub)=im_plot(:,:,cSub)./prctile(col(im_plot(:,:,cSub)),90);
mask_plot=cat(2,mask_plot,AllMask_cut{1,cSub}(:,crop_RL,plotSLC(cSub)));
s1_plot=cat(3,s1_plot,s1{i,cSub}(:,crop_RL,plotSLC(cSub)));
s2_plot=cat(3,s2_plot,s2{i,cSub}(:,crop_RL,plotSLC(cSub)));
end
% as(im_plot)
% plotSLC=1:6;

% colormap_start_P=round(Thres_P(1)/max(s1{i}(:))*4096);
colormap_start_P=round(Thres_P(1)/Thres_P(2)*4096);
clim_im=[0,2];
alpha_blobs=0.9;
cmap_P=hot(4096);
[allAx]=makeBlobPlot(flip(abs(im_plot),1),flip(s1_plot,1),'Thres',Thres_P,'SlcSel',1:6,'im_horz',6, ...
    'title_im',sprintf('Activation maps : TE %d ms',TE_list{i}),'colorbar',i==3, ...
    'caxis_im',clim_im,'alpha_blobs',0.8);


%  colormap_start_N=round(Thres_N(1)/Thres_N(2)*4096);
 alpha_blobs=0.9;
 cmap_N=flip(cmap_P,2);
% % cmap=flip(cool(4096),1);
% makeBlobPlot([],flip(s2_plot,1),Thres_N,plotSLC,clim_im,alpha_blobs,cmap_N)

makeBlobPlot([],flip(s2_plot,1),'Thres',Thres_N,'SlcSel',1:6,'im_horz',6,'title_im',sprintf('Activation maps : TE %d ms',TE_list{i}),...
    'negBold',true,'cmap',cmap_N,'colorbar',i==3,'allAx',allAx,'alpha_blobs',0.8)
%  hold on
%  contour(flip(mask_plot,1),'color','blue')
end


% summary stats
TE_sel=1:3;

MSC=cellfun(@(x,y)median(x(getmask(y))),sig_change_vec(TE_sel,:,1),Zscores(TE_sel,:,1));
MZ=cellfun(@(x,y)median(x(getmask(y))),Zscores(TE_sel,:,1),Zscores(TE_sel,:,1));
AV=cell2mat(nVoxels(TE_sel,:,1));
TSNR_TE=TSNR_mean(TE_sel,:)./sqrt([volTR{TE_sel,1}])';



% figure,
nexttile(t,[1 1]),scatter([TE_list{:}],MSC,100,'filled','LineWidth',1,'MarkerFaceAlpha',0.9),xlabel('TE [ms]'),title('median sig. change [%]','Color','w'),
hold on,plot([TE_list{:}],mean(MSC,2),'x-','LineWidth',3,'Markersize',15,'Color',[0.6350 0.0780 0.1840])
set(gca,'Box','on','YColor',[1 1 1],'Xcolor',[1 1 1],'Color',[0 0 0],'fontweight','bold','FontSize',10);
xticks([2 4 6]),xlim([1.5 6.5]),ylim([2 7]),grid on

nexttile,scatter([TE_list{:}],MZ,100,'filled','LineWidth',1,'MarkerFaceAlpha',0.9),xlabel('TE [ms]'),title('median z-score','Color','w')
hold on,plot([TE_list{:}],mean(MZ,2),'x-','LineWidth',3,'Markersize',15,'Color',[0.6350 0.0780 0.1840])
set(gca,'Box','on','YColor',[1 1 1],'Xcolor',[1 1 1],'Color',[0 0 0],'fontweight','bold','FontSize',10);
xticks([2 4 6]),xlim([1.5 6.5]),ylim([5.5 7.5]),grid on

nexttile,scatter([TE_list{:}],AV./1000,100,'filled','LineWidth',1,'MarkerFaceAlpha',0.9),xlabel('TE [ms]'),title('activated voxels x1000','Color','w')
hold on,plot([TE_list{:}],mean(AV./1000,2),'x-','LineWidth',3,'Markersize',15,'Color',[0.6350 0.0780 0.1840])
set(gca,'Box','on','YColor',[1 1 1],'Xcolor',[1 1 1],'Color',[0 0 0],'fontweight','bold','FontSize',10);
xticks([2 4 6]),xlim([1.5 6.5]),ylim([2 11]),grid on

nexttile,scatter([TE_list{:}],TSNR_TE,100,'filled','LineWidth',1,'MarkerFaceAlpha',0.9),xlabel('TE [ms]'),title('mean tSNR eff. [1/√s]','Color','w')
hold on,plot([TE_list{:}],mean(TSNR_TE,2),'x-','LineWidth',3,'Markersize',15,'Color',[0.6350 0.0780 0.1840])
set(gca,'Box','on','YColor',[1 1 1],'Xcolor',[1 1 1],'Color',[0 0 0],'fontweight','bold','FontSize',10);
xticks([2 4 6]),xlim([1.5 6.5]),ylim([9 17]),grid on


legend({'sub1','sub2','sub3','sub4','sub5','sub6','mean'},'textcolor','w','Location','eastoutside') 
set(gcf,'Position',[293 69 1100 1000])

%%
savefig(gcf,'TE_dependece_allsub_Smooth_masked')
FigH=gcf;
FigH.InvertHardcopy = 'off';
FigH.Color='black';
print(FigH,'TE_dependece_allsub_Smooth_masked.tiff','-dtiff','-r600')



%% register blobs to GRE image
gre_dir=dir('/ptmp/pvalsala/Spiral_bSSFP/sub-0*');

SlcSel=3:18; % 2+2 slice removed
mypermute=@(imx) permute(imx,[2 1 3 4]);
for cSub=1:length(Subjects)
    fprintf('TE: %s\n',Subjects{cSub})
    gunzip(fullfile(gre_dir(cSub).folder,gre_dir(cSub).name,'anat/GRE_TE16.nii.gz'));
gre_fn=dir(fullfile(gre_dir(cSub).folder,gre_dir(cSub).name,'anat/GRE_TE16.nii'));

for i=1:length(im_dir_all{1})   
z1_fn=dir(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'rendered_thresh_zstat*.nii.gz'));
gunzip(fullfile(z1_fn(1).folder,'rendered_thresh_zstat*.nii.gz'));
z1_fn=dir(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'rendered_thresh_zstat*.nii'));

resliced_vol=myspm_reslice(gre_fn, z1_fn, 'nearest','rG_');

end
resliced_vol2=myspm_reslice(z1_fn(1),gre_fn, 'nearest','rF_');
end


%% load and plot
clear gre_vol Z1vol_G 
plotSLC2=round([6,8,12,  10,6,8]*(88*0.4/(20*1)))+22;

for cSub=1:length(Subjects)
    fprintf('TE: %s\n',Subjects{cSub})
Seldata=@(imx) flip(permute(imx(50:430,20:ceil(end*0.45),plotSLC2(cSub)),[2 1 3 4 5]),2);
gre_vol{cSub}=Seldata(niftiread(fullfile(gre_dir(cSub).folder,gre_dir(cSub).name,'anat/GRE_TE16.nii')));

for i=1:length(im_dir_all{1})   
Z1vol_G{i,cSub}=Seldata(niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'rG_rendered_thresh_zstat1.nii')));

end
end
%%
figure(26),clf
tt=tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile()
 Z1vol_G_sel=Z1vol_G(1:3,:);
%   Z1vol_G_sel=Z1vol_G(4:6,:); %TR
im_plot=reshape(repmat(squeeze(cell2mat(gre_vol)),[1 1]),[],size(gre_vol{1},2),6);
s1_plot=squeeze(reshape(cat(3,Z1vol_G_sel{:}),[],size(gre_vol{1},2),3,6));
s1_mask=s1_plot>(4);
m_all=sum(s1_mask,3)>1; %all active

s1_mask=cat(3,s1_mask(:,:,1,:), s1_mask(:,:,2,:)&(~m_all), ...
    s1_mask(:,:,3,:)&(~m_all),m_all);
s1_im=squeeze(sum(s1_mask.*permute((1:4),[1 3 2]),3));
Thres_P=[0 5];

clim_im=[0,500];
alpha_blobs=0.9;
% cmap_P=jet(4096);

cmap_P=cat(1,[0 0 0],[1 0 0],[0 1 0],[0 0 1],[1 1 1]);

[allAx]=makeBlobPlot(flip(abs(im_plot),1),1.*flip(s1_im,1),'Thres',Thres_P,'SlcSel',1:6,'im_horz',3, ...
    'title_im',sprintf('TE scans (Z>%d)',(ifig+2)),'colorbar',false, ...
    'caxis_im',clim_im,'alpha_blobs',0.5,'cmap',cmap_P);

for i=1:3
text(1+400*(i-1),180,sprintf('S%d',i),'FontSize',16,'color','w')
text(1+400*(i-1),380,sprintf('S%d',i+3),'FontSize',16,'color','w')
end


%%%% TR

nexttile()

 Z1vol_G_sel=Z1vol_G(4:6,:); %TR
im_plot=reshape(repmat(squeeze(cell2mat(gre_vol)),[1 1]),[],size(gre_vol{1},2),6);
s1_plot=squeeze(reshape(cat(3,Z1vol_G_sel{:}),[],size(gre_vol{1},2),3,6));
s1_mask=s1_plot>(4);
m_all=sum(s1_mask,3)>1; %all active

s1_mask=cat(3,s1_mask(:,:,1,:), s1_mask(:,:,2,:)&(~m_all), ...
    s1_mask(:,:,3,:)&(~m_all),m_all);
s1_im=squeeze(sum(s1_mask.*permute((1:4),[1 3 2]),3));
Thres_P=[0 5];

clim_im=[0,500];
alpha_blobs=0.9;
% cmap_P=jet(4096);

cmap_P=cat(1,[0 0 0],[1 0 0],[0 1 0],[0 0 1],[1 1 1]);

[allAx]=makeBlobPlot(flip(abs(im_plot),1),1.*flip(s1_im,1),'Thres',Thres_P,'SlcSel',1:6,'im_horz',3, ...
    'title_im',sprintf('TR scans(Z>%d)',(ifig+2)),'colorbar',false, ...
    'caxis_im',clim_im,'alpha_blobs',0.5,'cmap',cmap_P);

for i=1:3
text(1+400*(i-1),180,sprintf('S%d',i),'FontSize',16,'color','w')
text(1+400*(i-1),380,sprintf('S%d',i+3),'FontSize',16,'color','w')
end

set(gcf,'Position', [401 101 1268 980],'InvertHardcopy','off')

%% calculate venn diagram table
figure(4),clf
tt=tiledlayout(1,6,'TileSpacing','none','Padding','tight');
colors = [
    255, 138, 101;...
    129, 199, 132;...
    79, 195, 247;....

]./255;
for cSub=1:6
   nexttile()
   data_cSub=cat(2,s1{1,cSub}(cMask),s1{2,cSub}(cMask),s1{3,cSub}(cMask))>=4;
% data_cSub=cat(4,s1{1:3,cSub})>4;
% data_cSub(end/2:end,:,:,:)=0;
% data_cSub(:,:,1:2,:)=0;
% data_cSub(:,:,19:20,:)=0;
data_cSub=reshape(data_cSub,[],3);
% all

% A, B, C, A&B, A&C, B&C and A&B&C.
ABC=sum((data_cSub(:,1)+data_cSub(:,2)+data_cSub(:,3)) ==3 ,'all');
AB=sum((data_cSub(:,1)+data_cSub(:,2)) ==2 ,'all')-ABC;
AC=sum((data_cSub(:,1)+data_cSub(:,3)) ==2 ,'all')-ABC;
BC=sum((data_cSub(:,2)+data_cSub(:,3)) ==2 ,'all')-ABC;

A=sum(data_cSub(:,1) ==1,'all')-ABC-AC-AB;
B=sum(data_cSub(:,2) ==1,'all')-ABC-BC-AB;
C=sum(data_cSub(:,3) ==1,'all')-ABC-AC-BC;

assert(A+AB+AC+ABC==sum(data_cSub(:,1) ==1,'all'),'something wrong With A');
assert(B+AB+BC+ABC==sum(data_cSub(:,2) ==1,'all'),'something wrong With B');
assert(C+AC+BC+ABC==sum(data_cSub(:,3) ==1,'all'),'something wrong With C');

venn_fields=[A, B, C, AB, AC, BC, ABC];
% venn_fields={'A', 'B', 'C', 'AB', 'AC', 'BC', 'ABC'};
venn(3,'labels',venn_fields,'sets',{'TE2','TE4','TE6'},'colors',linspecer(3),'labelC','w','alpha',0.6);
end
set(gcf,'color','black','Position',[445 662 809 164])

sgtitle('Postive activation distribution','Interpreter','tex','FontSize',14,'FontWeight','bold','Color',[1 1 1])

%%
figure(4),clf
tt=tiledlayout(1,6,'TileSpacing','none','Padding','tight');

for cSub=1:6
   nexttile()
   cMask=AllMask_cut{1,cSub};
data_cSub=cat(2,s1{4,cSub}(cMask),s1{5,cSub}(cMask),s1{6,cSub}(cMask))>=4;
% data_cSub(end/2:end,:,:,:)=0;
% data_cSub(:,:,1:2,:)=0;
% data_cSub(:,:,19:20,:)=0;
data_cSub=reshape(data_cSub,[],3);

% all

% A, B, C, A&B, A&C, B&C and A&B&C.
ABC=sum((data_cSub(:,1)+data_cSub(:,2)+data_cSub(:,3)) ==3 ,'all');
AB=sum((data_cSub(:,1)+data_cSub(:,2)) ==2 ,'all')-ABC;
AC=sum((data_cSub(:,1)+data_cSub(:,3)) ==2 ,'all')-ABC;
BC=sum((data_cSub(:,2)+data_cSub(:,3)) ==2 ,'all')-ABC;

A=sum(data_cSub(:,1) ==1,'all')-ABC-AC-AB;
B=sum(data_cSub(:,2) ==1,'all')-ABC-BC-AB;
C=sum(data_cSub(:,3) ==1,'all')-ABC-AC-BC;
venn_fields=[A, B, C, AB, AC, BC, ABC];
 venn(3,'labels',venn_fields,'sets',{'TR6','TR8','TR10'},'colors',linspecer(3),'labelC','w','alpha',0.6);
% venn(3,'labels',venn_fields,'sets',{'TE2','TE4','TE6'},'colors',linspecer(3),'labelC','w','alpha',0.6);
end
sgtitle('Postive activation distribution','Interpreter','tex','FontSize',14,'FontWeight','bold','Color',[1 1 1])

set(gcf,'color','black','Position',[445 662 809 164],'InvertHardcopy',off)

%%
function m=getmask(x)
cutoff=min([5000 length(x)]);
[xs,idx]=sort(x(:),1,'descend');
disp(xs(1))
m=zeros(size(x));
m(idx(1:cutoff))=1;
m=(m>0);
end



