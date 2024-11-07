

%% HWEZ
pn='/ptmp/pvalsala/PMNU-4ZUF/reg/moco_M377_peSpiral_R4x2_1p2iso_ss_sWE_run_B0MTI_DCFJackson/';

anat=niftiread(fullfile(pn,'coreg_moco_M377_peSpiral_R4x2_1p2iso_ss_sWE_run_B0MTI_DCFJackson_vol1_anat.nii'));
anat2=niftiread(fullfile(pn,'moco_M377_peSpiral_R4x2_1p2iso_ss_sWE_run_B0MTI_DCFJackson_vol1.nii'));
anat3=niftiread('/ptmp/pvalsala/PMNU-4ZUF/moco/allmoco/moco_M377_peSpiral_R4x2_1p2iso_ss_sWE_run_B0MTI_DCFJackson.nii');
%  func=niftiread('/ptmp/pvalsala/HWEZ/wh_blobology/moco_s1p2_negBOLD.feat/rendered_thresh_zstat1.nii.gz');
%   func=niftiread('/ptmp/pvalsala/PMNU-4ZUF/glm/smooth_1p2_moco_M377_peSpiral_R4x2_1p2iso_ss_sWE_run_B0MTI_DCFJackson.feat/rendered_thresh_zstat1.nii.gz');
   func=niftiread('/ptmp/pvalsala/PMNU-4ZUF/glm/smooth_0p0_moco_M377_peSpiral_R4x2_1p2iso_ss_sWE_run_B0MTI_DCFJackson.feat/rendered_thresh_zstat1.nii.gz');
  func2=niftiread('/ptmp/pvalsala/PMNU-4ZUF/glm/smooth_0p0_moco_M377_peSpiral_R4x2_1p2iso_ss_sWE_run_B0MTI_DCFJackson.feat/rendered_thresh_zstat2.nii.gz');


%% wh figure

title_str={'axial'};
 slcFrac=[0.62, 0.85 0.59 0.44];

mypermute=@(imx) ndflip(permute(imx,[2 1 3 4]),[1,2,3]);

[cim,acs,cf,composed]=cropim(anat2,{func,func2},mypermute,slcFrac);
fac=1.5;
anat_crop=imresize3(cf(mypermute(anat)),fac);
   anat2_crop=imresize3(cf(mypermute(anat2)),fac);
   func_crop=imresize3(cf(mypermute(func)),fac);
   func2_crop=imresize3(cf(mypermute(func2)),fac);
%
  figure(24),
clf
t=tiledlayout(2,3,'TileSpacing','tight');
 cmap_N=flip(hot(4096),2);

slc_range=@(slc) slc-1:slc+1;
cax_im=[40 7e2];
% imagesc(cat(2,composed{i})),colormap gray,axis image
 nexttile(),
 sl=round(slcFrac(1)*size(anat_crop,3));
 %axial
[allax_cb1]=makeBlobPlot(anat2_crop(:,:,sl),max(func_crop(:,:,slc_range(sl)),[],3),'Thres',[4,15],'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'caxis_im',cax_im);
[allax_cb1]=makeBlobPlot([],max(func2_crop(:,:,slc_range(sl)),[],3),'Thres',[4,8],'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'caxis_im',cax_im,'negBold',true,'cmap', cmap_N);
%coronal
nexttile(),
sl=round(slcFrac(2)*size(anat_crop,1));
[allax_cb1]=makeBlobPlot(squeeze(anat2_crop(sl,:,:))',squeeze(max(func_crop(slc_range(sl),:,:),[],1))','Thres',[4,15],'SlcSel',1,'im_horz',1,'title_im','bSSFP functional volume','colorbar',false,'caxis_im',cax_im);
[allax_cb2]=makeBlobPlot([],squeeze(max(func2_crop(slc_range(sl),:,:),[],1))','Thres',[4,8],'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'caxis_im',cax_im,'negBold',true,'cmap', cmap_N);
title(allax_cb2{1},'Functional volume','FontSize',20)

nexttile(),
sl=round(slcFrac(3)*size(anat_crop,2));
 [allax_cb1]=makeBlobPlot(squeeze(anat2_crop(:,sl,:))',squeeze(max(func_crop(:,slc_range(sl),:),[],2))','Thres',[4,15],'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'caxis_im',cax_im);
 [allax_cb1]=makeBlobPlot([],squeeze(max(func2_crop(:,slc_range(sl),:),[],2))','Thres',[4,8],'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'negBold',true,'cmap', cmap_N);
cax_im=[0 200];%[0 1.2e-4];
 nexttile(),
sl=round(slcFrac(1)*size(anat_crop,3));
 %axial
[allax_cb1]=makeBlobPlot(anat_crop(:,:,sl),max(func_crop(:,:,slc_range(sl)),[],3),'Thres',[4,15],'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'caxis_im',cax_im);
[allax_cb1]=makeBlobPlot([],max(func2_crop(:,:,slc_range(sl)),[],3),'Thres',[4,8],'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'caxis_im',cax_im,'negBold',true,'cmap', cmap_N);
setxlabel(allax_cb1{1},'Transversal')

%coronal
nexttile(),
sl=round(slcFrac(2)*size(anat_crop,1));
[allax_cb1]=makeBlobPlot(squeeze(anat_crop(sl,:,:))',squeeze(max(func_crop(slc_range(sl),:,:),[],1))','Thres',[4,15],'SlcSel',1,'im_horz',1,'title_im','anatomical reference','colorbar',false,'caxis_im',cax_im);
[allax_cb2]=makeBlobPlot([],squeeze(max(func2_crop(slc_range(sl),:,:),[],1))','Thres',[4,8],'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'caxis_im',cax_im,'negBold',true,'cmap', cmap_N);

title(allax_cb2{1},'Structural volume','FontSize',20)
setxlabel(allax_cb1{1},'Coronal')

nexttile(),

sl=round(slcFrac(3)*size(anat_crop,2));
 [allax_cb1]=makeBlobPlot(squeeze(anat_crop(:,sl,:))',squeeze(max(func_crop(:,slc_range(sl),:),[],2))','Thres',[4,15],'SlcSel',1,'im_horz',1,'title_im','','colorbar',true,'caxis_im',cax_im);
[allax_cb2]=makeBlobPlot([],squeeze(max(func2_crop(:,slc_range(sl),:),[],2))','Thres',[4,8],'SlcSel',1,'im_horz',1,'title_im','','colorbar',true,'negBold',true,'Allax',allax_cb1,'cmap', cmap_N);
allax_cb1{2}.Position=allax_cb1{2}.Position+[0.02,-0.1,0,0];
allax_cb2{2}.Position=allax_cb2{2}.Position+[0.02,-0.1,0,0];

setxlabel(allax_cb1{1},'Sagital')
 sgtitle('Whole brain 1.2 mm^3 protocol','Color','w','FontSize',20)
% set(gcf,'Postion',[207 146 1345 896])
%%
nexttile(t,5,[3 4])
% imagesc(cat(2,composed{i})),colormap gray,axis image
allax=makeBlobPlot(composed{4},composed{2},'Thres',[4 20],'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'caxis_im',[1 2e2]);
makeBlobPlot([],composed{3},'Thres',[4 6],'SlcSel',1,'im_horz',1,'title_im',title_str{2},'colorbar',false,...
    'AllAx',allax,'negBold',true,'cmap', cmap_N);

nexttile(t,9,[3 4])
 imagesc(composed{1}),colormap gray,axis image,box off,caxis(cax_im{cSub})
 ax1=gca;
 ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';   % remove x-axis

hold on
contour(composed{4},'-r','LevelList',[55 100],'LineWidth',2)
title(title_str{3},'Interpreter','latex','FontSize',18,'FontWeight','bold','Color',[1 1 1])

% makeBlobPlot(composed{4},composed{2},'Thres',[4 20],'SlcSel',1,'im_horz',1,'title_im','','colorbar',i==3,'caxis_im',[1 2e2]);
set(allax_cb1{2},'Position',[0.94 0.52 0.02 0.2])
set(allax_cb2{2},'Position',[0.94 0.29 0.02 0.2])
cropFrac=0.6;
picked={1,3,4,6};
for i=1:4
cSub=picked{i};
nexttile(t,37+(i-1)*3,[1 3])

[cim,acs,cf,composed]=cropim(im{ds,cSub},{s1{ds,cSub},s2{ds,cSub},brain{ds,cSub},T1{ds,cSub}},mypermute,slcFrac{cSub});
cropSel=floor(cropFrac*size(composed{1},1)):size(composed{1},1);
allax=makeBlobPlot(composed{4}(cropSel,:),composed{2}(cropSel,:),'Thres',[4 20],'SlcSel',1,'im_horz',1,'title_im','','colorbar',false,'caxis_im',[1 2e2]);
makeBlobPlot([],composed{3}(cropSel,:),'Thres',[4 6],'SlcSel',1,'im_horz',1,'title_im',title_str{i+3},'colorbar',false,...
    'AllAx',allax,'negBold',true,'cmap', cmap_N);
end

set(gcf,'Position',[427 269 750*1.414 750])

%%
savefig(gcf,'wh_threeaxis_subPMNU_377')
FigH=gcf;
FigH.InvertHardcopy = 'off';
FigH.Color='black';
print(FigH,'wh_threeaxis_subPMNU_377.tiff','-dtiff','-r300')

%%

function setxlabel(ax,textstr)
axes(ax);
ax.XAxis.Visible='on';
ax.XAxis.Color='black';
xlabel(textstr,'Color','w','FontSize',14)
end

function m=getmask(x)
cutoff=min([5000 length(x)]);
[xs,idx]=sort(x(:),1,'descend');
disp(xs(1))
m=zeros(size(x));
m(idx(1:cutoff))=1;
m=(m>0);
end


function [cim,Ic1,crop_fun,composed]=cropim(im,Other,mytransform,SlcFrac)

if(~exist('SlcFrac','var'))
SlcFrac=[0.5, 0.8 0.7 0.3];
end
im=mytransform(im);
thres=0.1;
   [idx1]=find(squeeze(mean(im,[2,3])) > thres*max(squeeze(mean(im,[2,3]))),1,'first')+6;
   [idx2]=find(squeeze(mean(im,[2,3])) > thres*max(squeeze(mean(im,[2,3]))),1,'last');
   [idy1]=find(squeeze(mean(im,[1,3])) > thres*max(squeeze(mean(im,[1,3]))),1,'first');
   [idy2]=find(squeeze(mean(im,[1,3])) > thres*max(squeeze(mean(im,[1,3]))),1,'last');
   [idz1]=find(squeeze(mean(im,[1,2])) > thres*max(squeeze(mean(im,[1,2]))),1,'first');
   [idz2]=find(squeeze(mean(im,[1,2])) > thres*max(squeeze(mean(im,[1,2]))),1,'last');
   disp([idx1,idx2,idy1,idy2,idz1,idz2])
   crop_fun=@(im3D)padarray(im3D(idx1:idx2,idy1:idy2,idz1:idz2),round(0.5*([164 164 164]-[160 164 124])),0,'both');
   
   
   cim=crop_fun(im);
   ACS= floor([size(cim,3) size(cim,1) size(cim,2) size(cim,2)].*SlcFrac);
   
   Ic1={cim(:,:,ACS(1)),squeeze(cim(ACS(2),:,:))',[],squeeze(cim(:,ACS(3),:)),squeeze(cim(:,ACS(4),:))};
   width=floor(size(Ic1{2},2)/2);
   cfac=10; % cut more in the middle
   padsz=5;
   Ic1{3}=cat(1,Ic1{4}((end-width-cfac):(end-cfac),:),flip(Ic1{5}(end-width-cfac:end-cfac,:),1));
   Ic1{3}=Ic1{3}(1:size(Ic1{2},2),:)';
   
   composed{1}=cat(1,Ic1{1},Ic1{2},zeros(padsz,size(Ic1{1},2)),Ic1{3});
   for i=1:length(Other)
%        Other{i}=mytransform(Other{i});
%        tmp=crop_fun(mytransform(Other{i}));
%        tmp2={tmp(:,:,ACS(1)),squeeze(tmp(ACS(2),:,:))',[],squeeze(tmp(:,ACS(3),:)),squeeze(tmp(:,ACS(4),:))};
%           tmp2{3}=cat(1,tmp2{4}((end-width-cfac):(end-cfac),:),flip(tmp2{5}(end-width-cfac:end-cfac,:),1));
%    tmp2{3}=tmp2{3}(1:size(tmp2{2},2),:)';
   composed{i+1}=0;%cat(1,tmp2{1},tmp2{2},zeros(padsz,size(tmp2{1},2)),tmp2{3});

   
   end
   
end
