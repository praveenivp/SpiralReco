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
imf=niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'filtered_func_data.nii.gz'));
beta1=niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'stats/pe1.nii.gz'));
beta2=niftiread(fullfile(feat_dir_all{cSub}(i).folder,feat_dir_all{cSub}(i).name,'stats/pe2.nii.gz'));
imf1=mean(imf(:,:,:,10:end),4);
sig_chang1{i,cSub}=mypermute(100*(beta1./imf1(:,:,:)));
 sig_chang2{i}=mypermute(100*(beta2./imf1(:,:,:))); %intercept
end
end

save 'dataHR_smooth.mat' s1 s2 im im1 volTR im_info sig_chang1 im_dir_all feat_dir_all Subjects 

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

%% p8 figure

figure(26),
clf
title_str={'Mean functional','Structural','Structural overlay','Subject 1','Subject 2','Subject 3','Subject 6'};
cSub=4;
ds=1; %0.8 mm
slcFrac={[0.6, 0.85 0.4 0.54],[0.6, 0.82 0.4 0.54 ],[0.5, 0.8 0.4 0.6],[0.7, 0.82 0.35 0.62 ],[0.5, 0.8 0.3 0.7],[0.65, 0.88 0.38 0.67]};
cax_imall={[0 1.2e-4],[0 1e-4] , [0 1.2e-4] , [5e-6 1.e-4] , [0 1.2e-4] ,[0 200]};
mypermute=@(imx) ndflip(permute(imx,[1 2 3 4]),[1,2,3]);

[cim,acs,cf,composed]=cropim(im{ds,cSub},{s1{ds,cSub},s2{ds,cSub},brain{ds,cSub},T1{ds,cSub}},mypermute,slcFrac{cSub},[20 20]);
t=tiledlayout(4,12,'TileSpacing','tight');
 
nexttile(t,1,[3 4]), 
 cmap_N=flip(hot(4096),2);
[allax_cb1]=makeBlobPlot(composed{1},composed{2}.*double(composed{4}>5),'Thres',[4 16],'SlcSel',1,'im_horz',1,'title_im','','colorbar',true,'caxis_im',cax_imall{cSub});
plotString(cim,allax_cb1{1})
allax_cb2=makeBlobPlot([],composed{3}.*double(composed{4}>0),'Thres',[4 8],'SlcSel',1,'im_horz',1,'title_im',title_str{1},'colorbar',true,...
    'AllAx',allax_cb1,'negBold',true,'cmap', cmap_N);

%overlay
nexttile(t,5,[3 4])
allax=makeBlobPlot(composed{4},composed{2}.*double(composed{4}>5),'Thres',[4 16],'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'caxis_im',[1 2e2]);
makeBlobPlot([],composed{3}.*double(composed{4}>5),'Thres',[4 8],'SlcSel',1,'im_horz',1,'title_im',title_str{2},'colorbar',false,...
    'AllAx',allax,'negBold',true,'cmap', cmap_N);
plotString(cim)

nexttile(t,9,[3 4])
 imagesc(composed{1}),colormap gray,axis image,box off,caxis(cax_imall{cSub})
 ax1=gca;
 ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';   % remove x-axis

hold on
contour(composed{4},'-r','LevelList',[55 100],'LineWidth',1.3)
title(title_str{3},'Interpreter','latex','FontSize',18,'FontWeight','bold','Color',[1 1 1])
plotString(cim)

sgtitle(' Subject 4 : $$0.8~mm^{3}$$ protocol','Interpreter','latex','FontSize',20,'FontWeight','bold','Color',[1 1 1])
%
% makeBlobPlot(composed{4},composed{2},'Thres',[4 20],'SlcSel',1,'im_horz',1,'title_im','','colorbar',i==3,'caxis_im',[1 2e2]);
set(allax_cb1{2},'Position',[0.94 0.54 0.02 0.2])
set(allax_cb2{2},'Position',[0.94 0.29 0.02 0.2])
cropFrac=0.55;
picked={1,2,3,6};
for i=1:4
cSub=picked{i};
nexttile(t,37+(i-1)*3,[1 3])

[cim,acs,cf,composed]=cropim(im{ds,cSub},{s1{ds,cSub},s2{ds,cSub},brain{ds,cSub},T1{ds,cSub}},mypermute,slcFrac{cSub},[-10 5]);
cropSel=floor(cropFrac*size(composed{1},1)):size(composed{1},1);
allax=makeBlobPlot(composed{4}(cropSel,:),composed{2}(cropSel,:).*double(composed{4}(cropSel,:)>5),'Thres',[4 16],'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'caxis_im',[1 2e2]);
makeBlobPlot([],composed{3}(cropSel,:).*double(composed{4}(cropSel,:)>5),'Thres',[4 8],'SlcSel',1,'im_horz',1,'title_im',title_str{i+3},'colorbar',false,...
    'AllAx',allax,'negBold',true,'cmap', cmap_N);
end

set(gcf,'Position',[429 149 1228 860])

%%
savefig(gcf,'Highres_p8_Sub4_2')
FigH=gcf;
FigH.InvertHardcopy = 'off';
FigH.Color='black';
print(FigH,'highres_p8_sub4_2.tiff','-dtiff','-r600')

%% p6 figure
Pthres=[4 16];
Nthres=[4 8];

  figure(25),
clf
title_str={'Mean functional','Structural','Structural overlay','Subject 1','Subject 3','Subject 4','Subject 6'};
cSub=2;
ds=2; %0.6 mm
slcFrac={[0.6, 0.8 0.4 0.54],[0.6, 0.82 0.46 0.56],[0.5, 0.8 0.46 0.56 ],[0.4, 0.83 0.4 0.6],[0.5, 0.82 0.3 0.7 ],[0.65, 0.88 0.38 0.67]};
cax_imall={[0 .5e-3],[3e-5 4.4e-4] , [0 .5e-3],[1e-5 .4e-3],[0 .5e-3],[0 .5e-3],[0 .5e-3]};
mypermute=@(imx) ndflip(permute(imx,[1 2 3 4]),[1,2,3]);

[cim,acs,cf,composed]=cropim(im{ds,cSub},{s1{ds,cSub},s2{ds,cSub},brain{ds,cSub},T1{ds,cSub}},mypermute,slcFrac{cSub},[22 25]);
t=tiledlayout(4,12,'TileSpacing','tight');
  
nexttile(t,1,[3 4]), 
% imagesc(cat(2,composed{i})),colormap gray,axis image
 cmap_N=flip(hot(4096),2);
[allax_cb1]=makeBlobPlot(composed{1},composed{2},'Thres',Pthres,'SlcSel',1,'im_horz',1,'title_im','','colorbar',true,'caxis_im',cax_imall{cSub});
allax_cb2=makeBlobPlot([],composed{3},'Thres',Nthres,'SlcSel',1,'im_horz',1,'title_im',title_str{1},'colorbar',true,...
    'AllAx',allax_cb1,'negBold',true,'cmap', cmap_N);
 plotStringp6(cim,allax_cb1{1})

nexttile(t,5,[3 4])
% imagesc(cat(2,composed{i})),colormap gray,axis image
allax=makeBlobPlot(composed{4},composed{2}.*double(composed{4}>0),'Thres',Pthres,'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'caxis_im',[1 2e2]);
makeBlobPlot([],composed{3}.*double(composed{4}>0),'Thres',Nthres,'SlcSel',1,'im_horz',1,'title_im',title_str{2},'colorbar',false,...
    'AllAx',allax,'negBold',true,'cmap', cmap_N);
plotStringp6(cim)
nexttile(t,9,[3 4])
 imagesc(composed{1}),colormap gray,axis image,box off,caxis(cax_imall{cSub})
 ax1=gca;
 ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';   % remove x-axis

hold on
contour(composed{4},'-r','LevelList',[55 100],'LineWidth',1.2)
title(title_str{3},'Interpreter','latex','FontSize',18,'FontWeight','bold','Color',[1 1 1])
plotStringp6(cim)

sgtitle(' Subject 2 : $$0.6~mm^{3}$$ protocol','Interpreter','latex','FontSize',20,'FontWeight','bold','Color',[1 1 1])
% makeBlobPlot(composed{4},composed{2},'Thres',[4 20],'SlcSel',1,'im_horz',1,'title_im','','colorbar',i==3,'caxis_im',[1 2e2]);
set(allax_cb1{2},'Position',[0.94 0.54 0.02 0.2])
set(allax_cb2{2},'Position',[0.94 0.29 0.02 0.2])
cropFrac=0.58;
picked={1,3,4,6};
%
for i=1:4
cSub=picked{i};
nexttile(t,37+(i-1)*3,[1 3]),cla

[cim,acs,cf,composed]=cropim(im{ds,cSub},{s1{ds,cSub},s2{ds,cSub},brain{ds,cSub},T1{ds,cSub}},mypermute,slcFrac{cSub},[-10 5]);
cropSel=floor(cropFrac*size(composed{1},1)):size(composed{1},1);
allax=makeBlobPlot(composed{4}(cropSel,:),composed{2}(cropSel,:).*double(composed{4}(cropSel,:)>0),'Thres',Pthres,'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'caxis_im',[1 2e2]);
makeBlobPlot([],composed{3}(cropSel,:).*double(composed{4}(cropSel,:)>0),'Thres',Nthres,'SlcSel',1,'im_horz',1,'title_im',title_str{i+3},'colorbar',false,...
    'AllAx',allax,'negBold',true,'cmap', cmap_N);
end

set(gcf,'Position',[429 149 1228 860])

%%
savefig(gcf,'Highres_p6_sub2_2')
FigH=gcf;
FigH.InvertHardcopy = 'off';
FigH.Color='black';
print(FigH,'highres_p6_sub2_2.tiff','-dtiff','-r600')



%% plot a slice of all subjects (p8)
figure(72),clf
tiledlayout(6,2,'TileSpacing','none','Padding','compact')
for i=1:6
    nexttile(i*2-1)
    
    imslc=createImMontage(mypermute(im{1,i}(:,:,4:5:25)),5);
    imslc=imslc./prctile(imslc,95,'all');
    imagesc(imslc),colormap('gray'),caxis([0 2]), axis image,
    xticks([]),yticks([])
            ax=gca;
    ax.XAxis.Color='black';
    ax.YAxis.Color='black';
    ylabel(sprintf('subject %d',i),'FontSize',12,'FontWeight','bold','Color','w')
    if(i==1), title('mean spiral images at 0.8 mm^3','FontSize',20,'Color','w'); end
    box off
end
ax.XAxis.Color='white';
ax.XAxis.FontSize=12;
sz=size(im{1,i},1);
xticks(sz/2:sz:sz*6),xticklabels(4:5:25)
xlabel('slice number','FontSize',16,'FontWeight','bold','Color','w')
set(gcf,'Color','black')
for i=1:6
    nexttile(i*2)
    
    imslc=createImMontage(mypermute(im{2,i}(:,:,round(2:3.25:16))),5);
    imslc=imslc./prctile(imslc,95,'all');
    imagesc(imslc),colormap('gray'),caxis([0 2]), axis image,
    xticks([]),yticks([])
        ax=gca;
    ax.XAxis.Color='black';
    ax.YAxis.Color='black';
    ylabel(sprintf('subject %d',i),'FontSize',12,'FontWeight','bold','Color','w')
    if(i==1), title('mean spiral images at 0.6 mm^3','FontSize',20,'Color','w'); end
    box off
end
ax.XAxis.Color='white';
ax.XAxis.FontSize=12;
sz=size(im{2,i},1);
xticks(sz/2:sz:sz*6),xticklabels(round(3:3.25:16))



xlabel('slice number','FontSize',16,'FontWeight','bold','Color','w')


%%
savefig(gcf,'data_overview')
FigH=gcf;
FigH.InvertHardcopy = 'off';
FigH.Color='black';
print(FigH,'data_overview.tiff','-dtiff','-r300')
%% plot a slice of all subjects (p6)
figure(72)
% tiledlayout(6,2,'TileSpacing','none','Padding','none')
for i=1:6
    nexttile(i*2)
    imslc=createImMontage(mypermute(im{2,i}(:,:,3:4:16)),4);
    imslc=imslc./prctile(imslc,95,'all');
    imagesc(imslc),colormap('gray'),caxis([0 2]), axis image,xticks([]),yticks([])
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
