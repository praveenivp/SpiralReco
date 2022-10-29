%%
SimObj=IDEASim('D:\VM Shared\simulations\nonsel_3d\peSpiral_R4x2_1p2iso_ns');
 SimObj.performGIRF();
EvalMask=SimObj.getEvalInfoMask();
ImagingEventBlocks=find(EvalMask.MDH_IMASCAN|EvalMask.MDH_ACQEND);
[~,Gxyz]=SimObj.getkTraj(ImagingEventBlocks);


%% get Siemens and GIRF B0
ECC=zeros([2 size(Gxyz,2) size(Gxyz,3)]);
for idx=1:length(ImagingEventBlocks)
    ECC(1,:,idx)=SimObj.EventBlocks(ImagingEventBlocks(idx)).B0Siemens;
    ECC(2,:,idx)=SimObj.EventBlocks(ImagingEventBlocks(idx)).B0GIRF;
end


%% get the middle of RF pulse
% [~,st]=max(SimObj.EventBlocks(ImagingEventBlocks(1)).RFD,'last');
st=find(SimObj.EventBlocks(ImagingEventBlocks(1)).RFD==max(SimObj.EventBlocks(ImagingEventBlocks(1)).RFD),1,'last');
st=round(st*(SimObj.Trace_headers.RFD.HORIDELTA/SimObj.Trace_headers.GRX.HORIDELTA));
G_0xyz=cat(1,-1*ECC(1,:,:)+0*ECC(2,:,:),Gxyz);
G_0xyz(:,1:st,:)=0;
k_0xyz=cumsum(G_0xyz,2); %(mT/m)*(10 us)
% GammaH1=42.577478518e3; %Hz/mT
GammaH1=42.575575e3;
k_0xyz=k_0xyz*GammaH1*10e-6; % 1 and (1/m)
%% get dwell time and ADC samples
taxis=+5e-6+linspace(0,size(G_0xyz,2)-1,size(G_0xyz,2))*10e-6; %s
cEventInfo=SimObj.ParseEventInfo(ImagingEventBlocks(1));
mhd=SimObj.ParseMeasHeader(ImagingEventBlocks(1));
cellIdx=find(cellfun(@(x) isfield(x,'ADC'),cEventInfo),1,'last');
ADCsamples=cEventInfo{cellIdx}.ADC.points;
ADCDurarion=cEventInfo{cellIdx}.ADC.Dur*1e-6; %s
ADCDwell=ADCDurarion/ADCsamples;
ADCstart=cEventInfo{cellIdx}.RelTime*1e-6;
taxis_adc=linspace(ADCstart,ADCstart+ADCDurarion-ADCDwell,ADCsamples);


% rotate gradient to Encoding space
%Encoding space to patient co-ordnidate to  (PRS -> SCT)
Q=quaternion(mhd{2}.sSD.aflQuaternion);
% Gradient to patient co-ordinate (XYZ*[1 -1 -1] -> SCT )
k_PRS=quat2rotm(Q)'*(k_0xyz(2:end,:).*[1; -1; -1;]);
k_PRS=reshape(k_PRS,3,size(k_0xyz,2),size(k_0xyz,3));
k_PRS=cat(1,k_0xyz(1,:,:),k_PRS(1,:,:),k_PRS(2,:,:),k_PRS(3,:,:));
%slice position
pos_SCT=[mhd{2}.sSD.sSlicePosVec.flSag ;mhd{2}.sSD.sSlicePosVec.flCor;mhd{2}.sSD.sSlicePosVec.flTra];
pos_PRS=quat2rotm(Q)'*(pos_SCT*1e-3); %m

% comapre Ktraj
%-5us is not required as Gradient are already 5us shifted in simulation
GradDelay=(2.9023e-6+ADCDwell*1.4486)-5e-6; %s
k_0PRS_interp=interp1(taxis,permute(k_PRS,[2 3 1]),taxis_adc+GradDelay);
k_0PRS_interp=permute(k_0PRS_interp,[3 1 2]);
k_PRS_interp=reshape(k_0PRS_interp(2:end,:,:),[3 sz(2:end)]);
%% FOV shift, undo Siemens ECC, add B0 from GIRF
B0_mod=exp(-1i*sum(bsxfun(@times,2*pi*k_0PRS_interp(2:4,:,:),[pos_PRS(1:2);0]),1));
B0_mod=B0_mod.*exp(-1i*1*pi*reshape(k_0PRS_interp(1,:,:),size(B0_mod)));


% %%
% intlv=1;
% par=77;
% % figure(67),
% clf,plotk(kFinal2(:,:,intlv,par))
% hold on,plotk(ra.KTraj(:,:,intlv,par)/(2*pi))

% do simple reco

% load('m1301_csm_k6.mat')
% ra=SpiralReco('meas_MID01301_FID30329_peSpiral_R4x2_1p2iso.dat','csm',csm);
SpiralPara=ra.SpiralPara;
flags=ra.flags;
 DCF=ra.DCF;

%  DCF=ones(size(ra.DCF));
%get sig
sig=ra.twix.image(:,:,:,:,1);
acq_sel=(ra.twix.image.Rep==1&ra.twix.image.Sli==1);
Lin_ordering=reshape(ra.twix.image.Lin(acq_sel),round(SpiralPara.Ninterleaves/SpiralPara.R_PE),[]);
Par_ordering=reshape(ra.twix.image.Par(acq_sel),round(SpiralPara.Ninterleaves/SpiralPara.R_PE),[]);
sig=sig(:,:,sub2ind([size(sig,3) size(sig,4)],Lin_ordering,Par_ordering));
sig=permute(sig,[2 1 3 4 5]);
sz=[size(sig,1),size(sig,2),size(Lin_ordering)];
sig= reshape(ra.D*sig(:,:),sz);


B0_mod=reshape(B0_mod,[1 sz(2:end)]);
B0_mod=B0_mod.*reshape(sqrt(DCF),[1 sz(2:end)]);
sig=bsxfun(@times,sig,B0_mod);

N=round(SpiralPara.FOV(1)/SpiralPara.Resolution);
% DCF=jacksonDCF2(squeeze(complex(k_PRS_interp(1,:,:,1),k_PRS_interp(2,:,:,1))),SpiralPara);

kmax=2*pi*(0.5/(SpiralPara.Resolution*1e-3));
kmax=[kmax;kmax;(2*pi*0.5)/(1e-3*SpiralPara.slice{1}.FOV_PRS(3)/SpiralPara.NPartitions)];

save('ForNUFFT.mat','k_0PRS_interp','kmax','DCF','SpiralPara')

csm2=permute(csm,[2 3 4 1]);
NUFFT_obj= StackofSpirals((2*pi)*k_PRS_interp./(2*kmax),(DCF),[N N SpiralPara.NPartitions],csm2,...
    'CompMode',flags.CompMode,'precision',flags.precision);

img = permute(NUFFT_obj'*(permute(sig,[2,3,4,1])),[4,1,2,3]);
as(cat(5,squeeze(img),squeeze(ra.img)),'')

%% CGSNESE
tic
im_pat=spiralCGSENSE(ra.NUFFT_obj,permute(ra.sig,[2,3,4,1]),...
    'maxit',flags.maxit,'tol',flags.tol,'reg',flags.reg,...
    'lambda',flags.reg_lambda);
toc
as(cat(5,squeeze(im_pat),squeeze(im_pat_sim)),'sel',':,:,37,1,1')
        
        %%
        figure(8),clf,plotk((2*pi)*0.999*k_PRS_interp(:,1:1:end,:,1:4:end)./kmax),view([0 90]),title('Simulated')
         figure(7), clf,plotk(ra.KTraj(:,1:1:end,:,1:4:end)./kmax),view([0 90])

% xlim([156.5344  215.0621]),
% zlim(   [-420 -410])



