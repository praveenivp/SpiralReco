%% spm segmentation
% addpath('/ptmp/pvalsala/Packages/spm12/')
% addpath(genpath('/ptmp/pvalsala/MATLAB'))
gre_vol=fullfile(pwd,'DICOM_pe_gre_CP_0p4iso_TE16_20230414102640_20.nii');


matlabbatch{1}.spm.spatial.preproc.channel.vols = {[gre_vol,',1']};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/ptmp/pvalsala/Packages/spm12/tpm/TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/ptmp/pvalsala/Packages/spm12/tpm/TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/ptmp/pvalsala/Packages/spm12/tpm/TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/ptmp/pvalsala/Packages/spm12/tpm/TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/ptmp/pvalsala/Packages/spm12/tpm/TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/ptmp/pvalsala/Packages/spm12/tpm/TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];

spm_jobman('run',matlabbatch)

%%
% skull strip after SPM segmentation
seg_st=dir('c*.nii');
mask1=niftiread(seg_st(1).name)|niftiread(seg_st(2).name)|niftiread(seg_st(3).name);
mask1=imdilate(imfill(mask1,'holes'),strel('sphere',9));

%% erode the mask to get bit more tight fit
mask2=imerode(imfill(mask1,'holes'),strel('sphere',9));
brain=double(niftiread(gre_vol));
as(cat(4,xor(mask2,mask1),brain.*mask2,brain))
%%
all_slc=[1 3 5];
for ijk=1:3
%% minimum intensity projection across sliding slices
brain=double(niftiread(gre_vol));
mIP_slc=all_slc(ijk); % only odd 1,3,5,
slc_res=0.4;%mm
brain_mip=im2col(padarray(reshape(brain,[],size(brain,3)),[0 floor(mIP_slc/2)],'replicate','both'),[1 mIP_slc]);
disp([size(brain_mip) numel(brain)])
brain_mip=reshape(min(brain_mip,[],1),size(brain));
brain_mip_masked=brain_mip.*mask2;
% as(brain_mip_masked)

%% frangi filt
clear J;
for jj=1:size(brain_mip_masked,3)
[J(:,:,jj)] = FrangiFilter2D(brain_mip_masked(:,:,jj), struct('verbose',false,'FrangiScaleRatio',1,'FrangiBetaOne',10,'FrangiBetaTwo',30));
end

%% visual
% as(cat(4,(mask2.*J>0.6),brain_mip_masked,brain),'select',':,:,3')

%% vessel distance with imdilate

ves_mask2=(mask2.*J>0.6);
upscale_fact=1;
ves_mask2_up=imresize3(double(ves_mask2),size(ves_mask2)*upscale_fact)>0.5;
nIter=40*upscale_fact;
dil_iter=zeros([size(ves_mask2) ,nIter]);
dil_iter(:,:,:,1)=ves_mask2;
dist=zeros(size(ves_mask2_up))*Inf;
dist(ves_mask2_up)=0;
old=ves_mask2_up;
for i=2:nIter
dil_iter(:,:,:,i)=(imdilate(dil_iter(:,:,:,i-1),strel('sphere',1)));
% dilated mask to dist
dist(xor(dil_iter(:,:,:,i-1),dil_iter(:,:,:,i)))=slc_res.*(i-1)/upscale_fact;

% new=(imdilate(old,strel('sphere',1)));
% dist(xor(old,new))=slc_res.*(i-1)/upscale_fact;
% old=new;
% fprintf('%.0f %%',100*(i/nIter))
end
dist_dila=imresize3(double(dist),size(ves_mask2));
dist_dila=dist_dila.*mask2;
% as(cat(4,dist_dila,ves_mask2,brain_mip_masked))


%% vessel distance map(slow Eucledian)
[Vx,Vy,Vz]=ind2sub(size(ves_mask2),find(ves_mask2>0));
ves_idx=cat(2,Vx(:),Vy(:),Vz(:));
[Xq,Yq,Zq]=ndgrid(1:size(ves_mask2,1),1:size(ves_mask2,2),1:size(ves_mask2,3));
voxel_idx=cat(2,Xq(mask2),Yq(mask2),Zq(mask2));
tic
[k,dist]=dsearchn(ves_idx,voxel_idx(:,:));
toc
dist_eucl=zeros(size(ves_mask2));
dist_eucl(mask2)=dist;

save(sprintf('out_%d.mat',ijk),'dist_eucl','dist_dila','mIP_slc','mask2','J','ves_mask2','all_slc')
end

%% write nifti files
nii_hdr=niftiinfo(gre_vol);
% nii_hdr.Datatype='double';
% nii_hdr.BitsPerPixel=64;
% nii_hdr.raw.datatype=64; %https://brainder.org/2012/09/23/the-nifti-file-format/
nii_hdr.Datatype='single';
nii_hdr.BitsPerPixel=32;
% nii_hdr.raw.datatype=16; %https://brainder.org/2012/09/23/the-nifti-file-format/
nii_hdr.Info="Vessel Distance in mm";
for i=1:length(all_slc)
    load(sprintf('out_%d.mat',i),'dist_eucl')
    niftiwrite(single(dist_eucl),sprintf('Vessel_dist_mm_mipSlc%d',all_slc(i)),nii_hdr);
end
% as(cat(4,dist,dist2,ves_mask2,ves))
%% try region grwoing
slc=22;
[poly,region_mask] = regionGrowing(squeeze(J(:,:,slc)),[],.4,400,false,true,false);
% [xq,yq]=meshgrid(1:size(J,1));
% [in,on]=inpolygon(xq,yq,poly(:,1),poly(:,2));
as((mask2(:,:,5:2:84).*J));

%%
as(cat(5,niftiread('c1p4_GRE.nii'),niftiread('c2p4_GRE.nii'),niftiread('c3p4_GRE.nii'),niftiread('c4p4_GRE.nii'),niftiread('c1p5_GRE.nii'),niftiread('c7p4_GRE.nii')))

%% some plots
Myori=@(x) ndflip(permute(x,[2 1 3 4]),1);
figure(2),clf,imagesc(createImMontage(Myori(brain(:,:,[10:12:80])),3))

title('GRE 0.4mm','Interpreter','latex','FontSize',24,'FontWeight','bold','Color',[1 1 1])
ax=gca;
set(ax,'Color',[0 0 0])
ax.YAxis.Color=[0 0 0];
ax.XAxis.Color=[0 0 0];

set(gcf,'Color',[0,0,0])
set(gcf,'InvertHardcopy','off')
axis image,colormap gray

%%
figure(3),clf,imagesc(createImMontage(Myori(brain_mip_masked(:,:,[10:12:80])),3))

title('mIP 3 slices + masked','Interpreter','latex','FontSize',24,'FontWeight','bold','Color',[1 1 1])
ax=gca;
set(ax,'Color',[0 0 0])
ax.YAxis.Color=[0 0 0];
ax.XAxis.Color=[0 0 0];
set(gcf,'Color',[0,0,0])
set(gcf,'InvertHardcopy','off')
axis image,colormap gray

%%
figure(4),clf,imagesc(createImMontage(Myori(J(:,:,[10:12:80])),3))

title('Frangi(2D) filter','Interpreter','latex','FontSize',24,'FontWeight','bold','Color',[1 1 1])
ax=gca;
set(ax,'Color',[0 0 0])
ax.YAxis.Color=[0 0 0];
ax.XAxis.Color=[0 0 0];
set(gcf,'Color',[0,0,0])
set(gcf,'InvertHardcopy','off')
axis image,colormap gray
%%
figure(5),clf,imagesc(createImMontage(Myori(ves_mask2(:,:,[10:12:80])),3))

title('Frangi(2D) (thres+masked)','Interpreter','latex','FontSize',24,'FontWeight','bold','Color',[1 1 1])
ax=gca;
set(ax,'Color',[0 0 0])
ax.YAxis.Color=[0 0 0];
ax.XAxis.Color=[0 0 0];
set(gcf,'Color',[0,0,0])
set(gcf,'InvertHardcopy','off')
axis image,colormap gray

%%
dist_eucl2=dist_eucl;
dist_eucl2(~mask2)=nan;
cim12=createImMontage(Myori(dist_eucl2(:,:,[10:12:80])),3);
% cim12(cim12==0)=nan;
 imAlpha=ones(size(cim12));
 imAlpha(isnan(cim12))=0;
figure(6),clf,imagesc(cim12,'AlphaData',imAlpha)

title('Vessel Distance(mm)','Interpreter','latex','FontSize',24,'FontWeight','bold','Color',[1 1 1])
ax=gca;
set(ax,'Color',[0 0 0])
ax.YAxis.Color=[0 0 0];
ax.XAxis.Color=[0 0 0];
set(gcf,'Color',[0,0,0])
set(gcf,'InvertHardcopy','off')
axis image,colormap(flip(summer)),%
gf=gcf;
 colorbar
 gf.Children(1).Color=[1,1,1];
 gf.Children(1).FontSize=20;
 
caxis([0,8])