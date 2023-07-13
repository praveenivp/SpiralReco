function [fwhm]=plotPSF(SpRecoObj,T2_star,B0)
% [fwhm]=plotPSF(SpRecoObj,T2_star,B0)
% fucntion to plot point spread function of SpiralReco object
% matches the prediction of theoretical FWHM values
%10.1016/j.mri.2012.04.017
% plotPSF(SpiralObj,Inf,0)./resolution_PRS should give ~1.4 for circular
% and 1.2 for cartesian
%%
if(nargin<2)
T2_star=18e-3; %s
B0=0; %Hz
end

sz=size(SpRecoObj.sig);
sig=ones(sz(2:end));

% apply sqrt(DCF)
sig=sig.*sqrt(SpRecoObj.DCF);


% add T2 start
taxis=SpRecoObj.time*1e-6; %s
sig_decay=exp(-taxis./T2_star);
% add additional B0 modulation
phase_evol=exp(1i*2*pi*B0*taxis);

sig=sig.*(sig_decay.*phase_evol).';

PSF=SpRecoObj.NUFFT_obj'*sig;

figure,clf
st_title={'phase','read','slice'};
for i=1:3
  
    [~,Idx]=max(abs(PSF),[],'all','linear');
    [psf_mid(1),psf_mid(2),psf_mid(3)]=ind2sub(size(PSF),Idx);
    switch(i)
        case 1
            psf_1d=squeeze(abs(PSF(:,psf_mid(2),psf_mid(3))));
        case 2
            psf_1d=squeeze(abs(PSF(psf_mid(1),:,psf_mid(3))));
        case 3
            psf_1d=squeeze(abs(PSF(psf_mid(1),psf_mid(2),:)));
            
    end
    psf_1d_interp=ifft(fftshift(padarray(ifftshift(fft(psf_1d(:))) ,round(0.5*(2^18-length(psf_1d))),'both')));

    
    xaxis=linspace(0,SpRecoObj.SpiralPara.FOV_PRS(i),length(psf_1d_interp)); %mm
    xaxis=xaxis-mean(xaxis);
subplot(3,1,i),plot(xaxis,real(psf_1d_interp),'LineWidth',2),xlabel('distance [mm]')
[fwhm(i),points]=getFWHM(xaxis,psf_1d_interp);
hold on
plot(points.x,points.y,'LineWidth',2)
title(sprintf('%s PSF: FWHM=%2.2f mm | Nominal=%1.2f mm',st_title{i},fwhm(i),SpRecoObj.SpiralPara.res_PRS(i)))
ylabel('Amplitude [a.u]')
 xlim([-1 1]*5)
end
sgtitle(sprintf('PSF T2*= %2.2f ms & B0= %d Hz',T2_star*1e3,B0))
set(gcf,'Color','w')

end


function [fwhm,points]=getFWHM(x,y)
x=x(:);
y=y(:);
halfMax = (max(y)-0*min(y)) *(0.5);
index1 = find(y >= halfMax, 1, 'first');
index2 = find(y >= halfMax, 1, 'last');
fwhm = x(index2) - x(index1);

points.x=[x(index1); x(index2)];
points.y=[y(index1);y(index2)];
end