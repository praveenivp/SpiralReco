%%
% addpath(genpath('C:\Users\pvalsala\Documents\Packages2\IDEAsim'))
%  simobj=IDEASim('D:\VM Shared\simulations\subject\sub-04\sub04_p8_r1');
 simobj=IDEASim('D:\VM Shared\simulations\subject\sub-04\sub04_p8_20rep');
 simobj.performGIRF();
 %

% simobj_nocorr=IDEASim('D:\VM Shared\simulations\subject\sub-04\sub04_p8_r1');
 simobj_nocorr=IDEASim('D:\VM Shared\simulations\subject\sub-04\sub04_p8_20rep');
EvalMask=simobj_nocorr.getEvalInfoMask();
 ImagingEventBlocks=find(EvalMask.MDH_IMASCAN|EvalMask.MDH_ACQEND);

 %%
 figure,simobj_nocorr.plot(ImagingEventBlocks(100),{'GRX','GRY','GRZ'})
simobj.plot(ImagingEventBlocks(100),{'GRX','GRY','GRZ'})
legend('Uncorrected','GIRF','location','north')

%%
grad_diff=simobj.EventBlocks(ImagingEventBlocks(100)).GRX(:)   - ...
   1.* simobj_nocorr.EventBlocks(ImagingEventBlocks(100)).GRX(:);

k=cumsum(grad_diff).*10e-6*42.567e3*2*pi;
figure,plot(k)

%%
figure,simobj.plot(ImagingEventBlocks(100:120),{'B0Siemens','B0GIRF'}),ylabel('B0 modulation (mT)')


%% 
load('E:\test\sub-04\fm_csm_MeasUID386.mat','fm_interp','csm');
rg=SpiralReco('X:\mrdata\echtdata\studies\44\experiments\JIJP-YO7X\TWIX\allData#S94Tuebingen#F47323#M385#D140423#T112053#peSpiral_R1_p8.dat', ...
    'CompMode','CPU2DHybrid','csm',csm,'fm',-1*fm_interp,'doB0corr','MTI','doCoilCombine','adapt2');


%% load R4 data
%  twix=mapVBVD('X:\mrdata\echtdata\studies\44\experiments\JIJP-YO7X\TWIX\allData#S94Tuebingen#F47324#M386#D140423#T112305#peSpiral_R4_p8.dat');
r4=SpiralReco(twix, ...
    'CompMode','CPU3D','csm',csm,'fm',-1*fm_interp,'doB0corr','none','doCoilCombine','adapt2','doPAT','CGSENSE','RepSel',3);


%% 
mhd=simobj.ParseMeasHeader(ImagingEventBlocks(1));

q=quaternion(mhd{1}.sSD.aflQuaternion);
pos_SCT=[mhd{1}.sSD.sSlicePosVec.flSag; mhd{1}.sSD.sSlicePosVec.flCor; mhd{1}.sSD.sSlicePosVec.flTra; ];
quatrotate(q,pos_SCT')

%% get all eventblocks for the rep;
% 7*28=196
[evntblock,st]=getEventblocks(simobj_nocorr,3);
%%
history_time=500e-3; %[s]
hist_point=round(history_time/10e-6)+3;
adc_time=[2660-930 (2660+6550)]; %[us]
adc_samp=adc_time(1)/10 :adc_time(2)/10 ;
Grad_hist=zeros(3,10+hist_point+length(adc_samp),length(evntblock));

for iii=1:length(evntblock)
    cblk=evntblock(iii);

    adc_samp2=simobj_nocorr.EventBlocks(cblk).StartTime/10+adc_samp;
if(adc_samp2(1)-hist_point>1)
    adc_samp2=adc_samp2(1)-hist_point:adc_samp2(end);
else
    adc_samp2=1:adc_samp2(end);
end
Grad_hist(1,1+end-length(adc_samp2):end,iii)=simobj_nocorr.Traces.GRX(adc_samp2);
Grad_hist(2,1+end-length(adc_samp2):end,iii)=simobj_nocorr.Traces.GRY(adc_samp2);
Grad_hist(3,1+end-length(adc_samp2):end,iii)=simobj_nocorr.Traces.GRZ(adc_samp2);
end

[Grad_hist]= shift5us(Grad_hist);
figure(10),clf,plot(squeeze(Grad_hist(1:3,:,1)'))

%%
history_time=0.01e-3; %[S]
hist_point=history_time/10e-6;
adc_time=[2660 (2660+6550)]; %[us]
adc_samp=adc_time(1)/10 :adc_time(2)/10 ;
Grad_nohist=zeros(3,hist_point+length(adc_samp),length(evntblock));

for iii=1:length(evntblock)
    cblk=evntblock(iii);

    adc_samp2=simobj_nocorr.EventBlocks(cblk).StartTime/10+adc_samp;
if(adc_samp2(1)-hist_point>1)
    adc_samp2=adc_samp2(1)-hist_point:adc_samp2(end);
else
    adc_samp2=1:adc_samp2(end);
end
Grad_nohist(1,1+end-length(adc_samp2):end,iii)=simobj_nocorr.Traces.GRX(adc_samp2);
Grad_nohist(2,1+end-length(adc_samp2):end,iii)=simobj_nocorr.Traces.GRY(adc_samp2);
Grad_nohist(3,1+end-length(adc_samp2):end,iii)=simobj_nocorr.Traces.GRZ(adc_samp2);
end
 Grad_nohist(:,1:2,:)=0;
figure(10),clf,plot(squeeze(Grad_nohist(1:3,:,8)).')

[Grad_nohist]= shift5us(Grad_nohist);



%% rotate to PRS both and apply GIRF 
grad_PRS=quatrotate(q,Grad_hist(:,:).');
grad_PRS=reshape(grad_PRS.',size(Grad_hist));
% apply GIRF
load(fullfile('C:\Users\pvalsala\Documents\Packages2\SpiralReco','kspace','PSF_MAR2022'),'PSF')
G_corr_SPH_hist=GIRF_correction_Freq(permute(Grad_hist,[2 1 3]),PSF,'isCrossTermPSFCorr',true,'isHigherOrder',false);
G_corr_SPH_hist=permute(G_corr_SPH_hist,[2 1 3]);
grad_PRS_corr=quatrotate(q,G_corr_SPH_hist(2:4,:).');
grad_PRS_corr=reshape(grad_PRS_corr.',size(Grad_hist));

%% no history
G_corr_SPH_nohist=GIRF_correction_Freq(permute(Grad_nohist,[2 1 3]),PSF,'isCrossTermPSFCorr',true,'isHigherOrder',false);
G_corr_SPH_nohist=permute(G_corr_SPH_nohist,[2 1 3]);
grad_PRS_corr_nohist=quatrotate(q,G_corr_SPH_nohist(2:4,:).');
grad_PRS_corr_nohist=reshape(grad_PRS_corr_nohist.',size(Grad_nohist));



st_grad=655;
figure,plot(grad_PRS_corr(1:3,end-st_grad:end,2).'-grad_PRS(1:3,end-st_grad:end,2).')
hold on
plot(grad_PRS_corr(1:3,end-st_grad:end,2).'-grad_PRS_corr_nohist(1:3,end-st_grad:end,2).')


%% Set interval to integrate and calc k
st_enc=655+93;
gammaH=2*pi*42.575575e3; %[Hz/mT]
k_PRS_corr=cumsum(grad_PRS_corr(:,end-st_enc:end,:),2)*gammaH*10e-6;
k_PRS_corr_nohist=cumsum(grad_PRS_corr_nohist(:,end-st_grad:end,:),2)*gammaH*10e-6;
k_PRS=cumsum(grad_PRS(:,end-st_enc:end,:),2)*gammaH*10e-6;


% cut only ADC part

k_PRS_corr(:,1:(end-st_grad-1),:)=[];
% k_PRS_corr_nohist(:,1:(end-st_grad-1),:)=[];
k_PRS(:,1:(end-st_grad-1),:)=[];

% Kz encoding needs to be added

k_PRS_corr_nohist(3,:,:)=k_PRS_corr_nohist(3,:,:)+k_PRS(3,:,:);

%% Kspace figure

taxis=(0:1:st_grad)*10e-3; %[ms

figure(9),clf
tt=tiledlayout(4,2,'TileSpacing','tight','Padding','compact');

nexttile()
plot(taxis,cumsum( squeeze(0-G_corr_SPH_nohist(1,end-st_grad:end,:)),1)*gammaH*10e-6)
set(gca,'ColorOrder',linspecer(size(k_PRS,3)));
xlabel('time [ms]'),ylabel('\Deltak_0 [rad]')
title('Trajectory deviations with only spiral waveform','Interpreter','none')
grid on
ylim([-1 1]*1)

nexttile()
plot(taxis,cumsum(squeeze(G_corr_SPH_hist(1,end-st_grad:end,:)-G_corr_SPH_nohist(1,end-st_grad:end,:)),1)*gammaH*10e-6)
set(gca,'ColorOrder',linspecer(size(k_PRS,3)));
title('Additional trajectory deviations with 0.5s gradient history','Interpreter','none')
% title('k_PRS-k_PRS_corr_nohist','Interpreter','none')
xlabel('time [ms]'),ylabel('\Deltak_0 [rad]')
ylim([-1 1]*1)
grid on

ylabel_str={'\Deltak_{x} [rad/m]','\Deltak_{y} [rad/m]','\Deltak_{z} [rad/m]'};

for i=1:3

nexttile()
plot(taxis,squeeze(k_PRS(i,:,:)-k_PRS_corr_nohist(i,:,:)))
set(gca,'ColorOrder',linspecer(size(k_PRS,3)));
% title('k_PRS-k_PRS_corr_nohist','Interpreter','none')
 ylim([-35 35]),yticks(-30:10:30)
xlabel('time [ms]'),ylabel(ylabel_str{i})
grid on

ax2=nexttile();

plot(taxis,squeeze(k_PRS_corr(i,:,:)-k_PRS_corr_nohist(i,:,:)));
set(gca,'ColorOrder',linspecer(size(k_PRS,3)));

xlabel('time [ms]'),ylabel(ylabel_str{i})
 ylim([-35 35]),yticks(-30:10:30)
grid on
end


fontsize(gcf,"increase")
fontsize(gcf,"increase")
fontsize(gcf,"increase")
set(gcf,'color','w')
cb=colorbar;colormap(gca,linspecer(size(k_PRS,3)))
cb.Ticks=[];
cb.Label.String='Acquistion order';
% need K0


%% gradient figure
taxis=(0:1:st_grad)*10e-3; %[ms
figure(10),clf
tt=tiledlayout(4,2);

nexttile()
plot(taxis,1e3*squeeze(0-G_corr_SPH_nohist(1,end-st_grad:end,:)))
set(gca,'ColorOrder',jet(size(k_PRS,3)));
title('Field deviations with only spiral waveform','Interpreter','none')
grid on
xlabel('time [ms]'),ylabel('B_0 [\muT]')
ylim([-15 15])

nexttile()
plot(taxis,1e3*squeeze(G_corr_SPH_hist(1,end-st_grad:end,:)-G_corr_SPH_nohist(1,end-st_grad:end,:)))
set(gca,'ColorOrder',jet(size(k_PRS,3)));
title('Additional Field deviations with 0.5s gradient history','Interpreter','none')
grid on
xlabel('time [ms]'),ylabel('B_0 [\muT]')
ylim([-15 15])

ylabel_str={'\DeltaG_{x} [mT/m]','\DeltaG_{y} [mT/m]','\DeltaG_{z} [mT/m]'};
for i=1:3

nexttile()
plot(taxis,squeeze(grad_PRS(i,end-st_grad:end,:)-grad_PRS_corr_nohist(i,end-st_grad:end,:)))
set(gca,'ColorOrder',jet(size(k_PRS,3)));

yticks([-0.5,-0.3,-0.1,0,0.1,0.3,0.5])
ylim([-0.5 0.5])
xlabel('time [ms]'),ylabel(ylabel_str{i})
grid on

nexttile()
plot(taxis,squeeze(grad_PRS_corr(i,end-st_grad:end,:)-grad_PRS_corr_nohist(i,end-st_grad:end,:)))
set(gca,'ColorOrder',jet(size(k_PRS,3)));

ylim([-0.5 0.5]),yticks([-0.5,-0.3,-0.1,0,0.1,0.3,0.5])
xlabel('time [ms]'),ylabel(ylabel_str{i})
grid on
end
cb=colorbar;colormap('jet')
fontsize(gcf,"increase")
fontsize(gcf,"increase")
fontsize(gcf,"increase")
set(gcf,'color','w')
cb.Ticks=[];
cb.Label.String='Acquistion order';
%% get siemens terms
Grad_nohist1=Grad_nohist;
Grad_nohist1(:,1,:)=0;

Grad_hist1=Grad_hist;
Grad_hist1(:,1,:)=0;

B0_siemens_nohist= getEddyB0driftIIR(permute(Grad_nohist1,[2 1 3]),r.twix,'IIR');
B0_siemens_hist= getEddyB0driftIIR(permute(Grad_hist,[2 1 3]),r.twix,'IIR');



figure,
tt=tiledlayout("flow");
nexttile()
plot(B0_siemens_nohist)
nexttile()
plot(B0_siemens_hist)
%%

figure,plotk(k_PRS_corr)

%% resample to dwell time grid

%% use it for reco
figure(56),clf
% plot(r.time-r.time(1),r.KTraj(1:3,:,1).')
intlv=1;
hold on
inGrid=0:10:(size(k_PRS_corr_nohist,2)-1)*10; %us
outGrid=rg.time-rg.time(1)+(2.9023+rg.SpiralPara.DwellTime*1.4486*1e-3)+5.1;
k_interp=ipermute(interp1(inGrid,permute(k_PRS ,[2 1 3]),outGrid,'linear','extrap'),[2 1 3]);
k_interp=[-1;1;-1].*k_interp;
plot(rg.time-rg.time(1),rg.KTraj(1:2,:,intlv).'-(k_interp(1:2,:,intlv))')       
%  ylim([-40 40])


k0_interp=cumsum( squeeze(0-G_corr_SPH_nohist(1,end-st_grad:end,:)),1)*gammaH*10e-6;
k0_interp=interp1(inGrid,k0_interp,outGrid,'linear','extrap');
kmax=pi./(rg.SpiralPara.res_PRS(:)*1e-3);
adc_time=rg.time*1e-6;
DCF=rg.DCF;
%   sig=rg.sig;
%  save('sub_p8_R4_nom.mat','k_interp','kmax','adc_time','DCF','k0_interp','sig');

%%
k_interp=ipermute(interp1(inGrid,permute(k_PRS_corr_nohist,[2 1 3]),outGrid,'linear',0),[2 1 3]);
k_interp=[-1;1;-1].*k_interp;
kmax=pi./(r4.SpiralPara.res_PRS(:)*1e-3);
% SOSobj=StackofSpirals(reshape(k_interp./(2*kmax),3,[],7,28),r4.DCF,[264 264 28],permute(csm,[2 3 4 1]),'CompMode','CPU2DHybrid');
tic
SOSB0obj=StackofSpiralsB0(reshape(k_interp./(2*kmax),3,[],7,28),r4.DCF,[264 264 28],permute(csm,[2 3 4 1]), -1*fm_interp,r4.time*1e-6,'Method','MTI',...
    'CompMode','GPU3D');
% im=SOSB0obj'*permute(r4.sig,[2 3 4 1]);
im_pat_knom=spiralCGSENSE(SOSB0obj,permute(r4.sig,[2,3,4,1]),...
                                    'maxit',10,'tol',1e-6,'reg','Tikhonov',...
                                   'lambda',1e-1);
toc
 as(sos(im_pat_knom,4),'title','hist','select',':,:,24')

%%
k_interp=ipermute(interp1(inGrid,permute(k_PRS_corr  ,[2 1 3]),outGrid,'linear',0),[2 1 3]);
k_interp=[-1;1;1].*k_interp;
kmax=pi./(r4.SpiralPara.res_PRS(:)*1e-3);
% SOSobj=StackofSpirals(reshape(k_interp./(2*kmax),3,[],7,28),r4.DCF,[264 264 28],permute(csm,[2 3 4 1]),'CompMode','CPU2DHybrid');
SOSB0obj=StackofSpiralsB0(reshape(k_interp./(2*kmax),3,[],7,28),r4.DCF,[264 264 28],permute(csm,[2 3 4 1]), -1*fm_interp,r4.time*1e-6,'Method','MTI',...
    'CompMode','GPU3D');
% im=SOSB0obj'*permute(r4.sig,[2 3 4 1]);
im_pat_kcorr=spiralCGSENSE(SOSB0obj,permute(r4.sig,[2,3,4,1]),...
                                    'maxit',10,'tol',1e-6,'reg','Tikhonov',...
                                    'lambda',1e-3);


k_interp=ipermute(interp1(inGrid,permute(k_PRS_corr_nohist  ,[2 1 3]),outGrid,'linear',0),[2 1 3]);
k_interp=[-1;1;1].*k_interp;
kmax=pi./(r4.SpiralPara.res_PRS(:)*1e-3);
% SOSobj=StackofSpirals(reshape(k_interp./(2*kmax),3,[],7,28),r4.DCF,[264 264 28],permute(csm,[2 3 4 1]),'CompMode','CPU2DHybrid');
SOSB0obj=StackofSpiralsB0(reshape(k_interp./(2*kmax),3,[],7,28),r4.DCF,[264 264 28],permute(csm,[2 3 4 1]), -1*fm_interp,r4.time*1e-6,'Method','MTI',...
    'CompMode','CPU3D');
% im=SOSB0obj'*permute(r4.sig,[2 3 4 1]);
im_pat_kcorr_nohist=spiralCGSENSE(SOSB0obj,permute(r4.sig,[2,3,4,1]),...
                                    'maxit',10,'tol',1e-6,'reg','Tikhonov',...
                                    'lambda',1e-3);

%%
as(sos(im,4)-squeeze(rg.img))


function [shifted_waveform]= shift5us(Grad)
sz=size(Grad);
taxis1=(0:(sz(2)-1))*10e-6;
taxis2=(0:(sz(2)-1))*10e-6 - 5e-6;

shifted_waveform=ipermute(interp1(taxis1,permute(Grad,[2 1 3]),taxis2,'linear',0),[2 1 3]);


end
