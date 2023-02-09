%%
[fn,pn]=uigetfile(fullfile(pwd,'*.nii*'));
as(niftiread(fullfile(pn,fn)))

%% SSFP
cd('/ptmp/pvalsala/HRXA/reg/moco_m73_B0MTI_DCFJackson')
anatVol=niftiread('coreg_moco_m73_B0MTI_DCFJackson_vol1_anat.nii');
funcVol=niftiread('moco_m73_B0MTI_DCFJackson_vol1.nii');
 figure(6),clf
%  FuncVol=abs(mean(V(:,:,:,1),4));
  transform_tra=@(x) ndflip(permute(x,[2 1 3 4]),[1 2 3]);
 SlcSel=6:15;
 im_horz_tra=5;
Thres_P=[3 -2];%[4 18];
Thres_N=[4 -3];%[4 18];
  [allax_func,st]=makeBlobPlot(funcVol,[],'Thres',Thres_P,'SlcSel',SlcSel,...
     'im_horz',im_horz_tra,'transform',transform_tra,'caxis_im',[0 5e-4],...
     'title_im','Mean functional volume with Activations');
 
  figure(7),clf
  [allax,st]=makeBlobPlot(anatVol,[],'Thres',Thres_P,'SlcSel',SlcSel,...
     'im_horz',im_horz_tra,'transform',transform_tra,'caxis_im',[0 2e2],...
     'title_im','Mean functional volume with Activations');
 
 %
 ax=allax{1}.Children;
 ax_func=allax_func{1};
 
nblines=1;
linestyle='r-';
linewidth=1;
% ax=gca;
     CData = sqrt(sum(get(ax,'CData').^2, 3));
    CData(isinf(CData)) = NaN;
    CData(isnan(CData)) = 0;
        
        [C,lh] = ...
            contour(ax_func,CData,...
            nblines,linestyle,'LineWidth',linewidth,'LabelSpacing',1000);
                [C2,lh2] = ...
            contour(allax{1},CData,...
            nblines,linestyle,'LineWidth',linewidth,'LabelSpacing',200);
%% GRE

cd('/ptmp/pvalsala/HRXA/reg/moco_DICOM_dzne_ep3d_Pat4_0p8iso_run1_20221221115340_23')
anatVol=niftiread('coreg_moco_DICOM_dzne_ep3d_Pat4_0p8iso_run1_20221221115340_23_vol1_anat.nii');
funcVol=niftiread('moco_DICOM_dzne_ep3d_Pat4_0p8iso_run1_20221221115340_23_vol1.nii');

%
 figure(61),clf
%  FuncVol=abs(mean(V(:,:,:,1),4));
  transform_tra=@(x) ndflip(permute(x,[2 1 3 4]),[1 2 3]);
 SlcSel=6:15;
 im_horz_tra=5;
Thres_P=[3 -2];%[4 18];
Thres_N=[4 -3];%[4 18];
  [allax_func,st]=makeBlobPlot(funcVol,[],'Thres',Thres_P,'SlcSel',SlcSel,...
     'im_horz',im_horz_tra,'transform',transform_tra,'caxis_im',[0 2.5e2],...
     'title_im','Mean functional volume with Activations');
 
  figure(71),clf
  [allax,st]=makeBlobPlot(anatVol,[],'Thres',Thres_P,'SlcSel',SlcSel,...
     'im_horz',im_horz_tra,'transform',transform_tra,'caxis_im',[0 2e2],...
     'title_im','Mean functional volume with Activations');
 
 %
 ax=allax{1}.Children;
 ax_func=allax_func{1};
 
nblines=1;
linestyle='r-';
linewidth=1;
% ax=gca;
     CData = sqrt(sum(get(ax,'CData').^2, 3));
    CData(isinf(CData)) = NaN;
    CData(isnan(CData)) = 0;
        
        [C,lh] = ...
            contour(ax_func,CData,...
            nblines,linestyle,'LineWidth',linewidth,'LabelSpacing',1000);
                [C2,lh2] = ...
            contour(allax{1},CData,...
            nblines,linestyle,'LineWidth',linewidth,'LabelSpacing',200);
%%
lh.LevelList