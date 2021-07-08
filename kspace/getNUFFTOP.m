function [FT,sig,adcTime]=getNUFFTOP(twix,RXY,RZ,CAIPI)
% [FT,sig,adcTime]=getNUFFTOP(twix,RXY,RZ,CAIPI)
% get cell array of @D NUFFFT operators for 3D stack of spiral with caipiu
% shift
cSet=1;
sig=squeeze(twix.image(:,:,:,:,1,1,1,1,1,cSet,:,:,:,:) );%first rep


SpiralPara=getSpiralPara(twix);
SpiralPara.GradDelay=15.4;%[1; 1; 1]*(SpiralPara.GRAD_RASTER_TIME_DEFAULT-4.5); 
soda_obj=getSoda(twix,SpiralPara);


 [~,G_xyz,grad_MOM]=GetGradients(twix,SpiralPara,soda_obj,1);
 load('GIRF_20200210_reg500.mat','PSF_time')
 G_corr=(GIRF_Correction(G_xyz,PSF_time,'isCrossTermPSFCorr',true));
 Grad=GradientXYZ2PRS(G_corr(:,2:4,:),soda_obj);
SpiralPara.grad_MOM=grad_MOM;
[KTraj,adcTime]=Grad2Traj(Grad,SpiralPara,'my');


            kmax=2*pi*(0.5/(SpiralPara.Resolution*1e-3));
            k_scaled=KTraj./(2*kmax);
            

N=SpiralPara.FOV(1)/SpiralPara.Resolution;
DCF=jacksonDCF2(KTraj,SpiralPara);

if(CAIPI==0)
    cintlv=1:RXY:SpiralPara.Ninterleaves;
FT= NUFFT((k_scaled(:,cintlv)),((DCF(:,cintlv))),1,0,[N,N], 2);
sig=sig(:,:,cintlv,:,:);
sig=performFOVShift(sig,KTraj(:,1:RXY:end),SpiralPara,soda_obj,sqrt(DCF(:,cintlv)));
else
    FT=cell(1,CAIPI+1);
    for i=0:(CAIPI)
        cintlv=twix.image.Lin((RXY*i)+(1:floor(SpiralPara.Ninterleaves/RXY)));
        FT{i+1}= NUFFT((k_scaled(:,cintlv)),(DCF(:,cintlv)),1,0,[N,N], 2);
    end
    sig_fov=zeros(ceil(size(sig)./[1 1 RXY 1]));
    for cpar=1:size(sig,4)
        cintlv=twix.image.Lin((RXY*(cpar-1))+(1:floor(SpiralPara.Ninterleaves/RXY)));
        sig_fov(:,:,:,cpar)=performFOVShift(sig(:,:,cintlv,cpar),KTraj(:,cintlv),SpiralPara,soda_obj,sqrt(DCF(:,cintlv)));
    end
     sig=sig_fov;
end            

end
        function sig_FOV=performFOVShift(sig,KTraj,SpiralPara,soda_obj,DCF)
            if(any(SpiralPara.slice{1}.Position ~=0))
                %             [~,idx]=sort(obj.SpiralPara.slice{1}.Normal);
                %             posi=1e-3*obj.SpiralPara.slice{1}.Position(idx); %m
                pos_PRS=GradientXYZ2PRS(1e-3*[1 -1 -1].*SpiralPara.slice{1}.Position,soda_obj,1); %only work for head first-supine
                
%                 if(obj.flags.doB0Driftcorr)
%                     [kHO]=Grad2TrajHigherorder(obj.Grad,obj.SpiralPara);
%                     B0_mod=exp(-1i*(real(obj.KTraj).*pos_PRS(1)+imag(obj.KTraj).*pos_PRS(2)+squeeze(kHO(:,3,:)).*pos_PRS(3) ));
%                 else
                B0_mod=exp(-1i*(real(KTraj).*pos_PRS(1)+imag(KTraj).*pos_PRS(2)));
%                 end
 
                B0_mod=B0_mod.*reshape(DCF,size(B0_mod));
                
                
                

                % negative displacement added to real part of kTRaj moves up down in array show
                sig_FOV=bsxfun(@times, sig,reshape(B0_mod,size(sig,1),1,size(sig,3)));
            end
        end


        function soda_obj=getSoda(twix,spiralPara)
            soda_obj=   SODA_OBJ( 'mrprot',twix.hdr);
            %             ReadoutFOV
            
            soda_obj.NPixelReadout=spiralPara.FOV(1)/spiralPara.Resolution;
            soda_obj.NPixelPhase=spiralPara.FOV(1)/spiralPara.Resolution;
            soda_obj.PixelSizePhase= soda_obj.ReadoutFOV/soda_obj.NPixelPhase;
            soda_obj.PixelSizeReadout= soda_obj.PhaseFOV/soda_obj.NPixelReadout;
            soda_obj=soda_obj.calcPixLoc();
            
            %             PhaseFOV
            
        end