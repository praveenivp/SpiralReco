%% plot higher order terms

mhd=simobj.ParseMeasHeader(ImagingEventBlocks(1));

q=quaternion(mhd{1}.sSD.aflQuaternion);
center_SCT=[mhd{1}.sSD.sSlicePosVec.flSag; mhd{1}.sSD.sSlicePosVec.flCor; mhd{1}.sSD.sSlicePosVec.flTra; ];
center_PRS=quatrotate(q,center_SCT');



res_PRS= r.SpiralPara.res_PRS;

im=ndflip(abs(squeeze(rg.img)),[1 2 3]);
load('E:\test\sub-04\sub4-p8_R4.mat','mean_im')
mean_im=ndflip(abs(squeeze(mean_im)),[1 2 3]);
mean_im=mean_im./max(mean_im,[],'all');
%%
voxel=[200,194,8];
figure(24),clf,
tt=tiledlayout(3,4,"TileSpacing","compact","Padding","compact");
im_h=nexttile(tt,6,[2 2 ]);
imagesc(mean_im(:,:,voxel(3))),colormap(gca,'gray'),axis image,clim([0 0.3])
title('reference 0.8mm image')

% not sure about the polarity of the 


tile_Nr=[9 5 1 2 3 4  8  12];
color_all=[lines(7);0 0.5 0.5];
% clear voxel_all;
for i=1:length(tile_Nr)


% pt_handle=drawpoint(im_h,"color",color_all(i,:));
% voxel=[pt_handle.Position(2) pt_handle.Position(1) 14];
% voxel_all{i}=voxel;

 voxel=voxel_all{i};
voxel_PRS_mm=(voxel - size(im)/2).*(-1*res_PRS)+center_PRS;
voxel_xyz_mm=quatrotate(q',voxel_PRS_mm);

rectangle(im_h,'Position',[[voxel(2) voxel(1)]-[1 1]  2 2],'EdgeColor',color_all(i,:),'LineWidth',5)

  [H0field]=getHigherOrderField(G_corr_SPH_hist(:,end-1000:end,end-190:end),voxel_xyz_mm*1e-3);
H0field=cumsum(H0field,1)*42.567e3*2*pi*10e-6; %[rad]
taxis=(0:1:(size(H0field,1)-1))*10e-3; %[ms]
ax_handle=nexttile(tt,tile_Nr(i));
plot(taxis(1:end),rad2deg(H0field(1:end,:)),'LineWidth',1.1);
grid on

%  [H0field]=getHigherOrderField(G_corr_SPH_nohist(:,:,:),voxel_xyz_mm*1e-3);
% H0field=cumsum(H0field,1)*42.567e3*2*pi*10e-6; %[rad]
% taxis=(0:1:(size(H0field,1)-1))*10e-3; %[ms]
% ax_handle=nexttile(tt,tile_Nr(i));
% plot(taxis(1:end),rad2deg(H0field(:,:)),'LineWidth',1.1);


%   [confield,header]= getConcomitantfield(G_corr_SPH_hist(2:4,19e3:end,1:28),voxel_xyz_mm*1e-3,9.4);

xlabel('time [ms]'),ylabel('phase error [deg]')
title(sprintf('(x,y,z)=(%.0f,%.0f,%.0f) mm',voxel_xyz_mm)) 
 set(gca,'ColorOrder',linspecer(size(H0field,2)));

set(ax_handle,'XColor',color_all(i,:),'YColor',color_all(i,:),'LineWidth',1.5)


end


cb=colorbar;colormap(gca,linspecer(size(H0field,2)))
fontsize(gcf,"increase")
fontsize(gcf,"increase")
fontsize(gcf,"increase")
set(gcf,'color','w')
cb.Ticks=[];
cb.Label.String='Acquistion order';

sgtitle('Phase error due to higher order eddy currents','fontsize',24,'fontweight','bold');

