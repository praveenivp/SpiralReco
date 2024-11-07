%% read all
dirst=dir('/ptmp/pvalsala/DKSR-UFYK/SNR/*.mat');
clearvars snr_all g_all
dirst(2)=[];
for i=1:7
    load(fullfile(dirst(i).folder,dirst(i).name),'snr','g')
    g_all{i}=g;
    snr_all{i}=snr;
end

%%
% [tsnr,panel_obj,mask]=plotTSNR(Vols);


%% get all TSNR

Subjects={'JEL7-26IB','DKSR-UFYK','H734-SPCU','JIJP-YO7X','YU3S-VKP3','FX6C-SZ37'};

% check the order and filter
clc
for cSub=1:length(Subjects)
% no smooth
        pn=sprintf('/ptmp/pvalsala/%s/moco/allmoco',Subjects{cSub});
        pn2=sprintf('/ptmp/pvalsala/%s/glm',Subjects{cSub});
         pn3=sprintf('/ptmp/pvalsala/%s/reg',Subjects{cSub});
 %TE     
 im_dir_TE= dir(fullfile(pn,'moco*TR10*.nii*'));
  feat_dir_TE=dir(fullfile(pn2,'smooth_0p0_moco_M*TR10*.feat'));
  reg_dir_TE=cat(1,dir(fullfile(pn3,'moco_M*TR10*')));
%TR
im_dir_TR= cat(1,dir(fullfile(pn,'moco*TR6*.nii*')),dir(fullfile(pn,'moco*TR8*.nii*')),dir(fullfile(pn,'moco*TE2_TR10*.nii*')));
feat_dir_TR=cat(1,dir(fullfile(pn2,'smooth_0p0_moco_M*TR6*.feat')),dir(fullfile(pn2,'smooth_0p0_moco_M*TR8*.feat')),dir(fullfile(pn2,'smooth_0p0_moco_M*TE2_TR10*.feat')));
reg_dir_TR=cat(1,dir(fullfile(pn3,'moco_M*TR6*')),dir(fullfile(pn3,'moco_M*TR8*')),dir(fullfile(pn3,'moco_M*TE2_TR10*')));

% high res
im_dir_HR= cat(1,dir(fullfile(pn,'moco*p8*.nii*')),dir(fullfile(pn,'moco*p6*.nii*')));
feat_dir_HR=cat(1,dir(fullfile(pn2,'smooth_0p0_moco_M*p8*.feat')),dir(fullfile(pn2,'smooth_0p0_moco_M*p6*.feat')));
reg_dir_HR=cat(1,dir(fullfile(pn3,'moco_M*p8*')),dir(fullfile(pn3,'moco_M*p6*')));

fprintf('TE: %s\n',Subjects{cSub})
disp({im_dir_TR.name})
disp({feat_dir_TR.name})
fprintf('\n\n')

im_dir_all{cSub}=[im_dir_TE;im_dir_TR;im_dir_HR];
feat_dir_all{cSub}=[feat_dir_TE;feat_dir_TR;feat_dir_HR];
reg_dir_all{cSub}=[reg_dir_TE;reg_dir_TR; reg_dir_HR];
end

volSel= {getVolSel(7,4,132),getVolSel(7,4,132),getVolSel(7,4,132),...TE
    getVolSel(12,6,216),getVolSel(9,4,156),getVolSel(7,4,132),...TR
    getVolSel(10,5,180),getVolSel(10,5,180)};


clc
% SlcSel=3:16; % 2+2 slice removed
mypermute=@(imx) permute(imx,[2 1 3 4]);
clear im1  im volTR im_info   T1 ribbon brain tsnr_all mask_all
for cSub=2
    fprintf('TE: %s\n',Subjects{cSub})

for i=1:length(im_dir_all{cSub})
    
 im_all=niftiread(fullfile(im_dir_all{cSub}(i).folder,im_dir_all{cSub}(i).name));
 
 % calculate TSNR
 im_all_picked=im_all(:,:,:,volSel{i});
 tsnr_all{i,cSub}=mean(im_all_picked,4)./std(im_all_picked,[],4);
 im_info{i,cSub}=niftiinfo(fullfile(im_dir_all{cSub}(i).folder,im_dir_all{cSub}(i).name));

volTR{i,cSub}=im_info{i,cSub}.PixelDimensions(4);
 im1{i,cSub}=mypermute(mean(im_all(:,:,:,5,1),4));
im{i,cSub}=mypermute(mean(im_all(:,:,:,5:end,1),4));

reg_folder=fullfile(reg_dir_all{cSub}(i).folder,reg_dir_all{cSub}(i).name);
brain_file=dir(fullfile(reg_folder,'coregBM*.nii'));
brain{i,cSub}=mypermute( niftiread(fullfile(reg_folder,brain_file(1).name)));
 mask_all{i,cSub}=brain{i,cSub}>0.1;
ribbon_file=dir(fullfile(reg_folder,'coregR*.nii'));
ribbon{i,cSub}=mypermute(niftiread(fullfile(reg_folder,ribbon_file(1).name)));


end
end

 im_all=niftiread('/ptmp/pvalsala/PMNU-4ZUF/moco/allmoco/moco_M377_peSpiral_R4x2_1p2iso_ss_sWE_run_B0MTI_DCFJackson.nii');

  im_all_picked=im_all(:,:,:,getVolSel(7,4,132));
 tsnr_all{i+1,cSub}=mean(im_all_picked,4)./std(im_all_picked,[],4);


% save 'TNSR_allsub.mat' im im1 volTR im_info  im_dir_all feat_dir_all reg_dir_all  Subjects tsnr_all mask_all

%% 1mm protocols
csub=2;
slc=12;
figure(77),clf, 
t=tiledlayout(3,8,'TileSpacing','compact','padding','compact');

mycolorbar= @(ax)colorbar('FontSize',12,'FontWeight','bold');
nexttile(1,[1,4])
snr_1mm=ndflip(abs(cat(4,snr_all{1:5})),[1 2 ])./sqrt(2.8);
imagesc(createImMontage(squeeze(snr_1mm(:,:,slc,:)),5)),colormap hot,caxis([0 120]),mycolorbar(gca),axis image
mask_1mm=ndflip(repmat(AllMask_cut{1,2},[1 1 1 5]),[1 2]);
hold on
contour(createImMontage(squeeze(mask_1mm(:,:,slc,:)),5),'Color','blue')

cMask=imerode(AllMask_cut{1,2},strel('sphere',3));
snr1mm_av=cellfun(@(x)mean(abs(x(cMask)),'omitnan'),snr_all(1:5));
snr1mm_std=cellfun(@(x)std(abs(x(cMask)),'omitnan'),snr_all(1:5)); % this is okay

for i=1:5
text(140+(i-1)*200,25,{sprintf('%.1f ±',snr1mm_av(i));sprintf('%.1f',snr1mm_std(i))},'Color','blue','FontSize',12,'FontWeight','bold')
end

xticks(100:200:1000),
% xticklabels(strrep({dirst(1).name(18:25),dirst(2).name(18:25),dirst(3).name(18:25),dirst(4).name(18:24),dirst(5).name(18:24)},'_',','))
xticklabels({'TE2/TR10','TE4/TR10','TE6/TR10','TE2/TR8','TE2/TR6'})
ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;
title('SNR_0 efficiency','Interpreter','tex')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;


nexttile(9,[1,4])
g_1mm=ndflip(abs(cat(4,g_all{1:5})),[1 2 ]);
imagesc(createImMontage(squeeze(g_1mm(:,:,slc,:)),5)),colormap(gca,'jet'),caxis([0.8 1.2 ]),mycolorbar(gca),axis image
xticks(100:200:1000),
xticklabels({'TE2/TR10','TE4/TR10','TE6/TR10','TE2/TR8','TE2/TR6'})
% xticklabels(strrep({dirst(1).name(18:25),dirst(2).name(18:25),dirst(3).name(18:25),dirst(4).name(18:24),dirst(5).name(18:24)},'_',','))
title('g-factor')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;

nexttile(17,[1,4])
tsnr_1mm=permute(ndflip(abs(cat(4,tsnr_all{1:3,csub},tsnr_all{5,csub},tsnr_all{4,csub})),[1 2 3]),[2 1 3 4])./sqrt(2.8);
imagesc(createImMontage(squeeze(tsnr_1mm(:,:,slc,:)),5)),colormap(gca,'hot'),caxis([0 40]),mycolorbar(gca),axis image
hold on
contour(createImMontage(squeeze(mask_1mm(:,:,slc,:)),5),'Color','blue')
xticks(100:200:1000),
% xticklabels(strrep({dirst(1).name(18:25),dirst(2).name(18:25),dirst(3).name(18:25),dirst(4).name(18:24),dirst(5).name(18:24)},'_',','))
xticklabels({'TE2/TR10','TE4/TR10','TE6/TR10','TE2/TR8','TE2/TR6'})
title('tSNR efficiency')


tsnr_trans=@(xxx) permute(ndflip(abs(xxx),[20 ]),[2 1 3 4]);
cMask=tsnr_trans(imerode(AllMask_cut{1,2},strel('sphere',3)))==1;
tsnr1mm_av=cellfun(@(xx)mean(abs(xx(cMask)),'omitnan'),tsnr_all(1:5,csub))./sqrt(2.8);
tsnr1mm_std=cellfun(@(xx)std(abs(xx(cMask)),'omitnan'),tsnr_all(1:5,csub))./sqrt(2.8);

for i=1:5
text(140+(i-1)*200,25,{sprintf('%.1f ±',tsnr1mm_av(i));sprintf('%.1f',tsnr1mm_std(i))},'Color','blue','FontSize',12,'FontWeight','bold')
end



yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;
set(gcf,'color','w')
 set(gcf,'Position', [-11 154 1921 868])


%0.8 mm
nexttile(5)
snr_hs=ndflip(abs(cat(4,snr_all{7})),[1 2])./(sqrt(1.96));
slc=floor(size(snr_hs,3)/2);
imagesc((squeeze(snr_hs(:,:,slc,:)))),colormap(gca,'hot'),caxis([0 100]),mycolorbar(gca),axis image
hold on

cMask=ndflip(imerode(AllMask_cut{2,csub},strel('sphere',3)),[1 2]);
contour(cMask(:,:,slc),'Color','blue')
text(175,35,{sprintf('%.1f ±',mean(snr_hs(cMask)));sprintf('%.1f',std(snr_hs(cMask)))},'Color','blue','FontSize',12,'FontWeight','bold')

xticks(size(snr_hs,1)*0.5:size(snr_hs,1):1000),
xticklabels(('0.8 mm'))
title('SNR_0 eff.','Interpreter','tex')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;


nexttile(13)
snr_hs=ndflip(abs(cat(4,g_all{7})),[1 2]);
slc=floor(size(snr_hs,3)/2);
imagesc((squeeze(snr_hs(:,:,slc,:)))),colormap(gca,'jet'),caxis([0.8 1.2 ]),mycolorbar(gca),axis image
xticks(size(snr_hs,1)*0.5:size(snr_hs,1):1000),
xticklabels(('0.8 mm'))
title('g-factor')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;

% TNSR
nexttile(21)
tsnr_p8=permute(ndflip(tsnr_all{7,cSub},[1 2]),[2 1 3 4])./(sqrt(1.96));
slc=floor(size(tsnr_p8,3)/2);
imagesc(squeeze(tsnr_p8(:,:,slc,:))),colormap(gca,'hot'),caxis([0 30]),mycolorbar(gca),axis image
xticks(100:200:1000),
xticklabels('0.8 mm')
title('tSNR eff.')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;
hold on

cMask=permute(ndflip(imerode(AllMask_cut{2,csub},strel('sphere',3)),[1 2]),[1 2 3 4]);
contour(cMask(:,:,slc),'Color','blue')
text(175,35,{sprintf('%.1f±',mean(tsnr_p8(cMask)));sprintf('%.1f',std(tsnr_p8(cMask)))},'Color','blue','FontSize',12,'FontWeight','bold')



%%%%%%%%%%%%%%%%%%%%0.6 mm

nexttile(6)
snr_hs=ndflip(abs(cat(4,snr_all{6})),[1 2])./(sqrt(2.16));
slc=floor(size(snr_hs,3)/2);
imagesc((squeeze(snr_hs(:,:,slc,:)))),colormap(gca,'hot'),caxis([0 50]),mycolorbar(gca),axis image
xticks(size(snr_hs,1)*0.5:size(snr_hs,1):1000),
xticklabels(('0.6 mm'))
title('SNR_0 eff.','Interpreter','tex')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;

hold on

cMask=ndflip(imerode(AllMask_cut{3,csub},strel('sphere',3)),[1 2]);
contour(cMask(:,:,slc),'Color','blue')
text(230,45,{sprintf('%.1f ±',mean(snr_hs(cMask)));sprintf('%.1f',std(snr_hs(cMask)))},'Color','blue','FontSize',12,'FontWeight','bold')



nexttile(14)
snr_hs=ndflip(abs(cat(4,g_all{6})),[1 2]);
slc=floor(size(snr_hs,3)/2);
imagesc((squeeze(snr_hs(:,:,slc,:)))),colormap(gca,'jet'),caxis([0.8 1.2 ]),mycolorbar(gca),axis image
xticks(size(snr_hs,1)*0.5:size(snr_hs,1):1000),
xticklabels(('0.6 mm '))
title('g-factor')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;

% TNSR
nexttile(22)
tsnr_p6=permute(ndflip(tsnr_all{8,cSub},[1 2]),[2 1 3 4])./(sqrt(2.16));
slc=floor(size(tsnr_p6,3)/2);
imagesc(squeeze(tsnr_p6(:,:,slc,:))),colormap(gca,'hot'),caxis([0 15]),mycolorbar(gca),axis image
xticks(size(tsnr_p6,2)* 0.5:size(tsnr_p6,2):1000),
xticklabels('0.6 mm')
title('tSNR eff.')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;

hold on

cMask=permute(ndflip(imerode(AllMask_cut{3,csub},strel('sphere',3)),[1 2]),[1 2 3 4]);
contour(cMask(:,:,slc),'Color','blue')
text(230,45,{sprintf('%.1f±',mean(tsnr_p6(cMask)));sprintf('%.1f',std(tsnr_p6(cMask)))},'Color','blue','FontSize',12,'FontWeight','bold')



%%%%%%%%%%%%%%%%%%%% whole head

 nexttile(7,[1 2]),cla
load('/ptmp/pvalsala/PMNU-4ZUF/SNR/M377_peSpiral_R4x2_1p2iso_ss_sWE_run1_highmem.mat','g','snr')

 snr_wh=abs(padarray(permute(ndflip(snr,[1 2]),[1 2 3 4])./(sqrt(2.88)),[0 0 20]));
slc=floor(size(snr_wh,3)/2)*[1 1];
imagesc([(squeeze(snr_wh(:,slc(2),:))') snr_wh(:,:,slc(1)) ]),colormap(gca,'hot'),caxis([0 200]),mycolorbar(gca),axis image
 hold on
cMask=AllMask_cut{1,7};
cMask=padarray(permute(ndflip(cMask,[1 2]),[1 2 3 4]),[0 0 20]);
contour([flipud(squeeze(cMask(:,slc(2),:))') cMask(:,:,slc(1)) ],'color','blue')

xticks(size(snr_wh,2)* 1:size(snr_wh,2):1000),
xticklabels(('1.2 mm'))
text(160,20,{sprintf('%.1f ±',mean(snr_wh(cMask)));sprintf('%.1f',std(snr_wh(cMask)))},'Color','blue','FontSize',12,'FontWeight','bold')
title('SNR_0 efficiency','Interpreter','tex')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;
 
 nexttile(15,[1 2])
 g_wh=abs(padarray(permute(ndflip(g,[1 2]),[1 2 3 4]),[0 0 20]));
slc=floor(size(g_wh,3)/2)*[1 1];
imagesc([(squeeze(g_wh(:,slc(2),:))') g_wh(:,:,slc(1)) ]),colormap(gca,'jet'),caxis([0.8 1.2]),mycolorbar(gca),axis image
xticks(size(snr_wh,1)*1:size(snr_wh,1):1000),
 xticklabels(('1.2 mm'))
 
title('tSNR efficiency')
title('g-factor')
 yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;
% 
 nexttile(23,[1 2]),cla
tsnr_wh=padarray(permute(ndflip(tsnr_all{9,cSub},[1 2]),[2 1 3 4])./(sqrt(2.16)),[0 0 20]);
slc=floor(size(tsnr_wh,3)/2)*[1 1];
imagesc([flipud(squeeze(tsnr_wh(:,slc(2),:))') tsnr_wh(:,:,slc(1)) ]),colormap(gca,'hot'),caxis([0 40]),mycolorbar(gca),axis image
hold on
cMask=AllMask_cut{1,7};
cMask=padarray(permute(ndflip(cMask,[1 2]),[1 2 3 4]),[0 0 20]);
contour([flipud(squeeze(cMask(:,slc(2),:))') cMask(:,:,slc(1)) ],'color','blue')

xticks(size(snr_wh,2)* 1:size(snr_wh,2):1000),
xticklabels('1.2 mm')
text(160,20,{sprintf('%.1f ±',mean(tsnr_wh(cMask)));sprintf('%.1f',std(tsnr_wh(cMask)))},'Color','blue','FontSize',12,'FontWeight','bold')
title('tSNR efficiency')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;


ax=gcf;
for iii=[1:6:24, 5:6:24]
tb=annotation(gcf,'textbox',...
[ax.Children.Children(iii).Position(1:2) 0 0]+ [-0.009 -0.04 0.03 0.03],...
'String',{'1/√s'},'FontSize',12,'Fontweight','bold','EdgeColor','w');
end
%%
savefig(gcf,'TNSR_all_masked')
FigH=gcf;
FigH.InvertHardcopy = 'off';
% FigH.Color='black';
print(FigH,'TNSR_all_masked.tiff','-dtiff','-r600')


%% other stuff
figure(24),clf

t2=tiledlayout(3,3,'TileSpacing','tight')

% %
nexttile
snr_hs=ndflip(abs(cat(4,snr_all{6})),[1 2])./(sqrt(1.96));
slc=floor(size(snr_hs,3)/2);
imagesc((squeeze(snr_hs(:,:,slc,:)))),colormap hot,caxis([0 100]),colorbar,axis image
xticks(size(snr_hs,1)*0.5:size(snr_hs,1):1000),
xticklabels(('0.8mm TR10'))
title('SNR_0/sqrt(1.96)')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;

nexttile
snr_hs=ndflip(abs(cat(4,snr_all{7})),[1 2])./(sqrt(2.16));
slc=floor(size(snr_hs,3)/2);
imagesc((squeeze(snr_hs(:,:,slc,:)))),colormap hot,caxis([0 100]),colorbar,axis image
xticks(size(snr_hs,1)*0.5:size(snr_hs,1):1000),
xticklabels(('0.6mm TR10'))
title('SNR_0/sqrt(2.16)')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;

nexttile
snr_hs=ndflip(abs(cat(4,snr_all{8})),[1 2])./(sqrt(2.88));
slc=floor(size(snr_hs,3)/2);
imagesc((squeeze(snr_hs(:,:,70,:)))),colormap hot,caxis([0 200]),colorbar,axis image
xticks(size(snr_hs,1)*0.5:size(snr_hs,1):1000),
xticklabels(('1.2 mm TR5'))
title('SNR_0/sqrt(2.88)')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;

nexttile
snr_hs=ndflip(abs(cat(4,g_all{6})),[1 2]);
slc=floor(size(snr_hs,3)/2);
imagesc((squeeze(snr_hs(:,:,slc,:)))),colormap(gca,'jet'),caxis([0.8 1.5]),colorbar,axis image
xticks(size(snr_hs,1)*0.5:size(snr_hs,1):1000),
xticklabels(('0.8mm TR10'))
title('g-factor')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;

nexttile
snr_hs=ndflip(abs(cat(4,g_all{7})),[1 2]);
slc=floor(size(snr_hs,3)/2);
imagesc((squeeze(snr_hs(:,:,slc,:)))),colormap(gca,'jet'),caxis([0.8 1.5]),colorbar,axis image
xticks(size(snr_hs,1)*0.5:size(snr_hs,1):1000),
xticklabels(('0.6mm TR10'))
title('g-factor')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;

nexttile
snr_hs=ndflip(abs(cat(4,g_all{8})),[1 2]);
slc=floor(size(snr_hs,3)/2);
imagesc((squeeze(snr_hs(:,:,70,:)))),colormap(gca,'jet'),caxis([0.8 1.5]),colorbar,axis image
xticks(size(snr_hs,1)*0.5:size(snr_hs,1):1000),
xticklabels(('1.2 mm TR5'))
title('g-factor')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;

% TNSR
nexttile
tsnr_p8=permute(ndflip(tsnr_all{6,5},[1 2]),[2 1 3 4])./(sqrt(1.96));
slc=floor(size(tsnr_p8,3)/2);
imagesc(squeeze(tsnr_p8(:,:,slc,:))),colormap(gca,'hot'),caxis([0 40]),colorbar,axis image
xticks(100:200:1000),
xticklabels('0.8mm')
title('TSNR/qrt(1.96)')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;


% TNSR
nexttile
tsnr_p6=permute(ndflip(tsnr_all{7,5},[1 2]),[2 1 3 4])./(sqrt(2.16));
slc=floor(size(tsnr_p6,3)/2);
imagesc(squeeze(tsnr_p6(:,:,slc,:))),colormap(gca,'hot'),caxis([0 20]),colorbar,axis image
xticks(100:200:1000),
xticklabels('0.6mm')
title('TSNR/sqrt(2.16)')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;

nexttile
% im_wh=niftiread('/ptmp/pvalsala/HWEZ/wh_blobology/moco_m1412_B0none_DCFJackson.nii');
tsnr_wh=permute(ndflip(mean(im_wh(:,:,:,5:end),4)./std(im_wh(:,:,:,5:end),[],4),[1 2]),[2 1 3 4]);
tsnr_wh=tsnr_wh./sqrt(2.88);
slc=floor(size(tsnr_wh,3)/2);
% as(tsnr_wh(:,:,slc))

imagesc(squeeze(tsnr_wh(:,:,slc,:))),colormap(gca,'hot'),caxis([0 40]),colorbar,axis image
xticks(100:200:1000),
xticklabels('WB 1.2mm')
title('TSNR (Subject: HWEZ) /sqrt(2.88)')
yticks([]), ax=gca;ax.XAxis.FontSize=12; ax.Title.FontSize=14;



set(gcf,'color','w')
set(gcf,'Position', [552 282 1232 833])



%%
savefig(gcf,'TNSR_Other')
FigH=gcf;
FigH.InvertHardcopy = 'off';
% FigH.Color='black';
print(FigH,'TNSR_other.tiff','-dtiff','-r300')


%%
VolSel=getVolSel(7,4,132)

%%
% r='allData#S94Tuebingen#F47946#M140#D190423#T171532#peSpiral_R1_p6.dat';

Measpath='/ptmp/pvalsala/YU3S-VKP3';
dir_st=dir(fullfile(Measpath,'TWIX','*peSpiral_R1*.dat'));
for i=6;%3:(length(dir_st)-1)
     SpObj=SpiralReco(fullfile(dir_st(i).folder,dir_st(i).name),'RepSel',1,...
    'doCoilCombine','sos','compMode','CPU3D'...
    ,'doDCF','Jackson','precision','double','NormNoiseData',false);
end

data=ones(size(SpObj))


function VolSel=getVolSel(rest,stim,Nvol)
% get all really resting
VTR=30/(rest+stim); %rough vTR
skip=round(6/VTR); %6s
bl=rest+stim; %block length

st=  skip+((bl+1):bl:Nvol); % start
Nblocks= rest-skip;
% skip first block and pick last (rest-skip) volumes until Nvol
VolSel=repmat(st,[Nblocks 1])+repmat((0:(Nblocks-1))',[1 length(st)]);

VolSel=VolSel(:);
end