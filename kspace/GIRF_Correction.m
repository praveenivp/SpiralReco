function [G_corr,parameters]=GIRF_Correction(G_xyz,PSF,varargin)
% [G_corr,parameters]=GIRF_Correction(G_xyz,PSF,varargin)
%This function performs GIRF Correction in Time domain
%G_xyz - uncorrected gradient waveform
%        (N x Axis x NIntlv) matrix correspond to each axis
%        G and PSF should have same sampling time.
%        Gradients should be in physical XYZ axis or provide proper rotation matrix
%PSF is in time domain (N x SpHa x Axis)
%    
%Parameters: Struct or name value pairs see below.

 %% Input parameter check
  switch(nargin)
     case 1
         error('need two input parameters Gradient[Samples x Axis x ...] and PSF [Samples x SPH x Axis] ');
     case 2
          warning('Using default Parmeters')
          parameters=struct("RotMat",[],"isHigherOrderPSFCorr",0,"isCrossTermPSFCorr",0,"isB0Term",0,"verbose",0,"ConvType","same");
          disp(parameters)
     case 3
         if(isstruct(varargin{1}))
            parameters=varargin{1};
         else
            error('Third paramter is optional struct or Name-value pair') 
         end
      otherwise
          p=inputParser;
         addParameter(p,'RotMat',[],@(x) ismatrix(x) && all(size(x,1)==[3,3]));
         addParameter(p,'isHigherOrderPSFCorr',false,@(x) islogical(x));
         addParameter(p,'isCrossTermPSFCorr',true,@(x) islogical(x));
         addParameter(p,'isB0Term',true,@(x) islogical(x));
         addParameter(p,'verbose',false,@(x) islogical(x));
          %use 'full' for complete conv output or 'same' for truncated output
         addParameter(p,'ConvType','same',@(x) any(validatestring(x,{'same','full'})));
          parse(p,varargin{:});
          parameters=p.Results;
           
  end
          

 %make input gradient matrix into a 3D matrix
 G_xyz=reshape(G_xyz,size(G_xyz,1),size(G_xyz,2),[]);

%Calculate physical gradients from XYZ gradient
if(~isempty(parameters.RotMat))
    sz=size(G_xyz);
    G_xyz=parameters.RotMat*reshape(permute(G_xyz,[2 1 3]),sz(2),[]);
    G_xyz=ipermute(reshape(G_xyz,sz(2),sz(1),sz(3)),[2 1 3]);
end

%Spherical harmonics selector:
SPH_Sel=1:4; % B0,X,Y,Z by default

if(parameters.isHigherOrderPSFCorr) 
    SPH_Sel =1:size(PSF,2); 
end
if(strcmpi(parameters.ConvType,'same'))
    gc=zeros([size(G_xyz,1) max(SPH_Sel) size(G_xyz,3) 3]);
else
    gc=zeros([size(G_xyz,1)+size(PSF,1)-1 max(SPH_Sel) size(G_xyz,3) 3]);
end
 
for cintlv=1:size(G_xyz,3)
    for axis=1:3
     for sph=SPH_Sel
         gc(:,sph,cintlv,axis)=conv(G_xyz(:,axis,cintlv),PSF(:,sph,axis),parameters.ConvType);
     end
    end
end

% %Different Convolution mode:
% if(strcmpi(parameters.ConvType,'same'))
%    idx=101; %pad_size in GIRF computation
% gc=gc(idx:idx+size(G_xyz,1)-1,:,:,:);
% end


%remote all cross term compenents when cross term effects are not required
%Only removing Cross effects for X,Y,Z axis
if(~parameters.isCrossTermPSFCorr)
    gc(:,[3 4],:,1)=0;
    gc(:,[2, 4],:,2)=0;
    gc(:,[2, 3],:,3)=0;   
end

gc=sum(gc,4);

if(parameters.verbose)
figure,plot(gc(:,2:4,1)),hold on,plot(G_xyz(:,1:3))
end


if(parameters.isB0Term)
    G_corr=gc;
else
    G_corr=gc(:,2:4,:);
end

end

%old stuff 
% %Frequency domain
% function [G_corr,parameters]=GIRF_Correction_freq(G,PSF,parameters,verbose)
% 
% %G - uncorrected gradient waveform
%     %G: (N x Axis x NIntlv) matrix correspond to each axis
%  %PSF is in time domain
%     
%  %Parameters: Structure with following parameter
%   parameters.rotmat=eye(3);
%   parameters.isHigherOrderPSFCorr= false;
%   parameters.isCrossTermPSFCorr=true;
%   parameters.isB0Term=true;
%   
%   %warning
%   if(size(G,1)>1500)
%    warn('PSF FT might be unstable')
%    %Need higher PSF of large time points or reduce input gradient length
%   end
%  
%   %Calculate physical gradients from XYZ gradient
% %   G_physical=parameters.rotmat*G;
%   
% ZeroPadLength=800;%floor((4095-size(G,1))/2);  %800; % 8ms length of PSF measurement
% FFTLength= size(G,1)+ZeroPadLength;
% G_pad=padarray(G,[ZeroPadLength 0 0],0,'both');
% 
% %   PSF_FT=fft(PSF,FFTLength,1);
% %   PSF_self=cat(2,PSF_FT(:,2,1),PSF_FT(:,3,2),PSF_FT(:,4,3));
% %   freq_axis=linspace(-0.5/10e-6,0.5/10e-6,size(PSF_self,1));
% %   pl(freq_axis,fftshift(PSF_self),[-50e3 50e3])
% 
% 
% % %calculate PSF
% load('D:\Data\GIRF_fieldcamera\20190718\GIRF_mat17082020_2','cal_sig','rec_sig')
% 
% win_size=  FFTLength; 
% time_sel=1:800;
% in_sig=fft(cal_sig(time_sel,:,:,:,:,:),win_size,1);
% out_sig=fft(rec_sig(time_sel,:,:,:,:,:),win_size,1);
%  num=squeeze(sum(conj(in_sig).*out_sig,[2,4,5]));
%  denom=squeeze(sum(conj(in_sig).*in_sig,[2,4,5]));
%  PSF=num./denom;
% PSF(isnan(PSF))=0;
% PSF(isinf(PSF))=0;
% freq_axis=linspace(-1/1e-6,1/1e-6,size(PSF,1));
%  PSF_self=cat(2,PSF(:,2,1),PSF(:,3,2),PSF(:,4,3));
%   pl(freq_axis,fftshift(PSF_self),[-1.e6 1.e6])
% 
% PSF_FT=PSF;
% 
% 
% G_FT=fft(G_pad,FFTLength,1);
% 
% %Spherical harmonics selector:
% SPH_Sel=2:4; % X,Y,Z by default
% if(parameters.isB0Term) 
%     SPH_Sel = 1:4; 
% end
% if(parameters.isHigherOrderPSFCorr) 
%     SPH_Sel =1:size(PSF,2); 
% end
% 
% G_corr=zeros([size(G,1) length(SPH_Sel)   size(G,3)]);
% for cintlv=1:size(G,3)
% 
% GX_effects=ifft(bsxfun(@times,G_FT(:,1,cintlv),PSF_FT(:,SPH_Sel,1)),[],1);
% GY_effects=ifft(bsxfun(@times,G_FT(:,2,cintlv),PSF_FT(:,SPH_Sel,2)),[],1); 
% GZ_effects=ifft(bsxfun(@times,G_FT(:,3,cintlv),PSF_FT(:,SPH_Sel,3)),[],1);
% 
% %remote all cross term compenents when cross term effects are not required
% if(~parameters.isCrossTermPSFCorr)
%     GX_effects(:,[2,3]+parameters.isB0Term)=0;
%     GY_effects(:,[1 3]+parameters.isB0Term)=0;
%     GZ_effects(:,[1, 2]+parameters.isB0Term)=0;   
% end
% %remove padding
% % G_corr(:,:,cintlv)=real(GX_effects(1:size(G,1),:)+GY_effects(1:size(G,1),:)+GZ_effects(1:size(G,1),:));
% G_corr(:,:,cintlv)=real(GX_effects(ZeroPadLength+(1:size(G,1)),:)+GY_effects(ZeroPadLength+(1:size(G,1)),:)+GZ_effects(ZeroPadLength+(1:size(G,1)),:));
% end
% end
