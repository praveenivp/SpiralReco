function [B0_drift,phase_drift]=getEddyB0driftIIR(G_xyz,twix_obj,mode)
%[B0,phase]=getEddyB0drift(G_xyz,twix_obj)
% 
% Function calculates temporal B0 drifts due the eddy currents induced by
% the gradients. It uses decaying exponential model to fit the B0
% fluctuation. It is implemented at hardware level in siemens scanner.
%
%G_xyz: ND matrix [time x grad axis x interleaves x rep .....]
%       Sampling rate should be 10us
%twix_obj: raw data object from mapvbvd
%
%OUTPUTS:
%B0 [mT] : Same time scale as G_xyz, tail is cropped
%          ND-matrix [time x interleaves x rep ...] 
%phase [radians]: same time scale as G_xyz, tail is cropped
%                 ND-matrix [time x interleaves x rep ...] 

if(nargin<3)
    mode='IIR';
end
if(nargin<2||~isfield(twix_obj,'hdr'))
    %load default B0 terms (9T at 30.03.2022)
    %time constants are in seconds
    %amplitude/100 is in uT (Siemens for some reason give them in percentage of 1uT) 
CoeffsX=struct("Amp",[-0.0335807017982 0.0215377993882 0.044569298625],"Tau",[0.100905001163 0.0276404991746 0.00049977801973]);
CoeffsY=struct("Amp",[0.0207749009132 -0.206313997507 0.141762003303],"Tau",[1.99937999249 0.236003994942 0.169891998172]);
CoeffsZ=struct("Amp",[0.120907001197 0.216192007065 0.0163282006979],"Tau",[1.23278999329 0.445013999939 0.0538129992783]);

else
    CoeffsX.Amp = cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.asGPAData{1}.sB0CompensationX.aflAmplitude(1:3));
CoeffsX.Tau = cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.asGPAData{1}.sB0CompensationX.aflTimeConstant(1:3));
CoeffsY.Amp = cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.asGPAData{1}.sB0CompensationY.aflAmplitude(1:3));
CoeffsY.Tau = cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.asGPAData{1}.sB0CompensationY.aflTimeConstant(1:3));
CoeffsZ.Amp = cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.asGPAData{1}.sB0CompensationZ.aflAmplitude(1:3));
CoeffsZ.Tau = cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.asGPAData{1}.sB0CompensationZ.aflTimeConstant(1:3));
    
end

 gammaH=42.575575e3; %Hz/mT
 dt=10e-6; %s gradient raster time
sz=[size(G_xyz) 1];
G_xyz=reshape(G_xyz,sz(1),sz(2),prod(sz(3:end)));

%plot impulse response function
% max_tc=max([CoeffsX.Tau CoeffsY.Tau CoeffsZ.Tau]); % seconds
% %After 5 time constants there is usually nothing left(~0.67%) 
% time_axis=(0:dt:ceil(max_tc*5)); %1s
% [h]=plot_impulseFunction(twix_obj,time_axis)

 SR_xyz=cat(1,zeros([1 sz(2) ,prod(sz(3:end))]) ,(diff(G_xyz,1,1)./(10e-6*1e3))); % slew rate mT/m/ms

 B0x = RCfilter_timedomain(squeeze(SR_xyz(:,1,:)), CoeffsX.Amp/100, CoeffsX.Tau, dt,mode); %uT
 B0y = RCfilter_timedomain(squeeze(SR_xyz(:,2,:)), CoeffsY.Amp/100, CoeffsY.Tau, dt,mode); %uT
 B0z = RCfilter_timedomain(squeeze(SR_xyz(:,3,:)), CoeffsZ.Amp/100, CoeffsZ.Tau, dt,mode); %uT
 B0_drift=1e-3*(B0x+B0y+B0z); %mT

  %rehape back to original size
  B0_drift=reshape(B0_drift,[sz(1) ,sz(3:end)]);
  if(nargout>1)
   phase_drift=cumsum(B0_drift,1)*dt*(2*pi*gammaH); %radians
  phase_drift=reshape(phase_drift,[sz(1) ,sz(3:end)]);
  %interpolate to adc time grid
  taxis_grad=0:dt:dt*(size(phase_drift,1)-1);
  dw=twix_obj.hdr.Phoenix.sRXSPEC.alDwellTime{1}*1e-9; %s
  ADC_points=twix_obj.image.NCol;
  taxis_adc=0:dw:(ADC_points-1)*dw;
  phase_drift=interp1(taxis_grad,phase_drift,taxis_adc,'linear',0);
  end
  
end

%%

function fw = RCfilter_timedomain(dgdt, amp, tau, dt,mode)
% function fw = RCfilter_timedomain(dgdt, tau, dt)
switch mode
    case 'FIR' % bit slow only for testing
        t=0:dt:(10*max(tau(:))); %5*tau -10*tau
        % impulseFunctionsum=@(t,amp,tc) sum(repmat(amp(:),[ 1 length(t)]).*exp(-1*repmat(t(:),[ 1 length(tc)])'./ repmat(tc(:),[1 length(t)])),1);
        % psf=impulseFunctionsum(t,amp,tau);
        psf=exp(-1*(t+dt*0.5).*(1./tau(:)));
        psf=amp(:).*psf;
        psf=sum(psf,1);
        fw=zeros([size(dgdt,1)+length(t)-1 size(dgdt,2)]);
        for cIntlv=1:size(dgdt,2)
            fw(:,cIntlv)=fw(:,cIntlv)+conv(dgdt(:,cIntlv),psf(:));
        end
        fw((size(dgdt,1)+1):end,:)=[];
      case 'IIR'
        %  adaoted from :https://github.com/filip-szczepankiewicz/safe_pns_prediction.git
         fw = zeros([length(tau) size(dgdt,2) size(dgdt,1)]);
            alpha =1+((dt)./(tau(:)));
            fw(:,:,1) = repmat(dgdt(1,:),[length(tau) 1]);
            for s = 2:length(dgdt)
                fw(:,:,s) = (alpha(:) .* dgdt(s,:) + (2-alpha(:)).* fw(:,:,s-1));
            end
            fw=permute(sum(amp(:).*fw,1),[3 2 1]);
end

end

function [h]=plot_impulseFunction(twix_obj,t)
impulseFunction=@(t,amp,tc) (repmat(amp(:),[ 1 length(t)]).*exp(-1*repmat(t(:),[ 1 length(tc)])'./ repmat(tc(:),[1 length(t)]))).';
impulseFunctionsum=@(t,amp,tc) sum(repmat(amp(:),[ 1 length(t)]).*exp(-1*repmat(t(:),[ 1 length(tc)])'./ repmat(tc(:),[1 length(t)])),1);

    
fprintf('\nAmplitude GRAD X: [%1.6f %1.6f %1.6f ] \n',cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.sB0CompensationX.aflAmplitude(1:3)))
fprintf('Decay constant(s) : [%1.6f %1.6f %1.6f ] \n',cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.sB0CompensationX.aflTimeConstant))
fprintf('Amplitude GRAD Y: [%1.6f %1.6f %1.6f ] \n',cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.sB0CompensationY.aflAmplitude(1:3)))
fprintf('Decay constant(s) : [%1.6f %1.6f %1.6f ] \n',cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.sB0CompensationY.aflTimeConstant))
fprintf('Amplitude GRAD Z: [%1.6f %1.6f %1.6f ] \n',cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.sB0CompensationZ.aflAmplitude(1:3)))
fprintf('Decay constant(s) : [%1.6f %1.6f %1.6f ] \n',cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.sB0CompensationZ.aflTimeConstant))
 % x gradient  
 if(nargin<2)
 t=linspace(0,5,2048); %s
 end
amp=cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.sB0CompensationX.aflAmplitude);
amp(amp==0)=[];
tc=cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.sB0CompensationX.aflTimeConstant);
tc(tc==0)=[];
hx=impulseFunctionsum(t,amp,tc);
hx_ind=impulseFunction(t,amp,tc);


legend_String=cell(1,1+length(tc));
legend_String{1}='All Components';
for i=1:length(tc)
    legend_String{i+1}=sprintf('%d*exp(-t/%d)',amp(i),tc(i));
end 
figure,subplot(311),plot(t,hx,'LineWidth',2),hold on ,plot (t,hx_ind,'--'),title('B_0 Compensation : response function X'),xlabel('Seconds'),ylabel('% per 1 mT/m/ms'),legend(legend_String)

%grad Y
amp=cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.sB0CompensationY.aflAmplitude);
amp(amp==0)=[];
tc=cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.sB0CompensationY.aflTimeConstant);
tc(tc==0)=[];
hy=impulseFunctionsum(t,amp,tc);
hy_ind=impulseFunction(t,amp,tc);


legend_String=cell(1,1+length(tc));
legend_String{1}='All Components';
for i=1:length(tc)
    legend_String{i+1}=sprintf('%d*exp(-t/%d)',amp(i),tc(i));
end 
subplot(312),plot(t,hy,'LineWidth',2),hold on ,plot (t,hy_ind,'--'),title('B_0 Compensation : response function Y'),xlabel('Seconds'),ylabel('% per 1 mT/m/ms'),legend(legend_String)

%grad Y
amp=cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.sB0CompensationZ.aflAmplitude);
amp(amp==0)=[];
tc=cell2mat(twix_obj.hdr.Phoenix.sGRADSPEC.sB0CompensationZ.aflTimeConstant);
tc(tc==0)=[];
hz=impulseFunctionsum(t,amp,tc);
hz_ind=impulseFunction(t,amp,tc);


legend_String=cell(1,1+length(tc));
legend_String{1}='All Components';
for i=1:length(tc)
    legend_String{i+1}=sprintf('%d*exp(-t/%d)',amp(i),tc(i));
end 
subplot(313),plot(t,hz,'LineWidth',2),hold on ,plot (t,hz_ind,'--'),title('B_0 Compensation : response function Z'),xlabel('Seconds'),ylabel('% per 1 mT/m/ms '),legend(legend_String)

h=[hx;hy;hz];

end