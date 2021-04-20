
function [G_corr,parameters]=GIRF_Correction(G_xyz,PSF,varargin)

%G - uncorrected gradient waveform
    %G: (N x Axis x NIntlv) matrix correspond to each axis
    %G and PSF should have same sampling time.
 %PSF is in time domain
%    
%Parameters: Structure with following parameter
 
 
 % TASKS
% 0. Rotate gradients interleaves(will do in GradtoK)
% 1. calculate PSF with appropriate zeropadding(not required)
% 2. apply PSF(done)
% 3. Crossterm correction(done)
% 4. give higher order kspace(done)
% 5. verbose mode
% 6. Rotation matrix/quatinions support
% 7. write parameter parser and help text


%time domain
 
 %% Input parameter check
  switch(nargin)
     case 1
         error('need two input parameters Gradient[Samples x Axis x ...] and PSF [Samples x SPH x Axis] ');
     case 2
          warning('Using default Parmeters')
          parameters.RotMat=eye(3);
          parameters.isHigherOrderPSFCorr= false;
          parameters.isCrossTermPSFCorr=true;
          parameters.isB0Term=true;
          parameters.verbose=false;
          parameters.ConvType='same';
          disp(parameters)
     case 3
         if(isstruct(varargin{1}))
            parameters=varargin{1};
         else
            error('Third paramter is optional struct or Name-value pair') 
         end
      otherwise
          p=inputParser;
         addParameter(p,'RotMat',eye(3),@(x) ismatrix(x) && all(size(x,1)==[3,3]));
         addParameter(p,'isHigherOrderPSFCorr',false,@(x) islogical(x));
         addParameter(p,'isCrossTermPSFCorr',true,@(x) islogical(x));
         addParameter(p,'isB0Term',true,@(x) islogical(x));
         addParameter(p,'verbose',false,@(x) islogical(x));
          %use 'full' for complete conv output or 'same' for truncated
          %output
         addParameter(p,'ConvType','same',@(x) any(validatestring(x,{'same','full'})));
          parse(p,varargin{:});
          parameters=p.Results;
           
  end
          
%%

%make input gradient matrix into a 3D matrix
 G_xyz=reshape(G_xyz,size(G_xyz,1),size(G_xyz,2),[]);
 
  %Calculate physical gradients from XYZ gradient
%   G_physical=parameters.rotmat*G;


%Spherical harmonics selector:
SPH_Sel=1:4; % B0,X,Y,Z by default

if(parameters.isHigherOrderPSFCorr) 
    SPH_Sel =1:size(PSF,2); 
end
% PSF=flip(PSF,1);
 gc=zeros([size(G_xyz,1)+size(PSF,1)-1  max(SPH_Sel) size(G_xyz,3) 3]);
for cintlv=1:size(G_xyz,3)
    for axis=1:3
     for sph=SPH_Sel
         gc(:,sph,cintlv,axis)=conv(G_xyz(:,axis,cintlv),PSF(:,sph,axis));
     end
    end
end


% %sanity check: Normally center/maximum position of all PSF function(self term) is the same 
% [~,idx]=max(PSF(:,[2 3+size(PSF,2) 4+size(PSF,2)*2]));
% if(sum(diff(idx))~=0)
%      warning('PSF has different Max idx: [%d %d %d]',idx(1),idx(2),idx(3)) 
% end
% 
% idx(1)=idx(1)+1;
% %Different Convolution mode:
% if(strcmpi(parameters.ConvType,'same'))
% gc=gc(idx(1):idx(1)+size(G_xyz,1)-1,:,:,:);
% end
% 




%Different Convolution mode:
if(strcmpi(parameters.ConvType,'same'))
   idx=101; %pad_size in GIRF computation
gc=gc(idx:idx+size(G_xyz,1)-1,:,:,:);
end


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


%Frequency domain
function [G_corr,parameters]=GIRF_Correction_freq(G,PSF,parameters,verbose)

%G - uncorrected gradient waveform
    %G: (N x Axis x NIntlv) matrix correspond to each axis
 %PSF is in time domain
    
 %Parameters: Structure with following parameter
  parameters.rotmat=eye(3);
  parameters.isHigherOrderPSFCorr= false;
  parameters.isCrossTermPSFCorr=true;
  parameters.isB0Term=true;
  
  %warning
  if(size(G,1)>1500)
   warn('PSF FT might be unstable')
   %Need higher PSF of large time points or reduce input gradient length
  end
 
  %Calculate physical gradients from XYZ gradient
%   G_physical=parameters.rotmat*G;
  
ZeroPadLength=800;%floor((4095-size(G,1))/2);  %800; % 8ms length of PSF measurement
FFTLength= size(G,1)+ZeroPadLength;
G_pad=padarray(G,[ZeroPadLength 0 0],0,'both');

%   PSF_FT=fft(PSF,FFTLength,1);
%   PSF_self=cat(2,PSF_FT(:,2,1),PSF_FT(:,3,2),PSF_FT(:,4,3));
%   freq_axis=linspace(-0.5/10e-6,0.5/10e-6,size(PSF_self,1));
%   pl(freq_axis,fftshift(PSF_self),[-50e3 50e3])


% %calculate PSF
load('D:\Data\GIRF_fieldcamera\20190718\GIRF_mat17082020_2','cal_sig','rec_sig')

win_size=  FFTLength; 
time_sel=1:800;
in_sig=fft(cal_sig(time_sel,:,:,:,:,:),win_size,1);
out_sig=fft(rec_sig(time_sel,:,:,:,:,:),win_size,1);
 num=squeeze(sum(conj(in_sig).*out_sig,[2,4,5]));
 denom=squeeze(sum(conj(in_sig).*in_sig,[2,4,5]));
 PSF=num./denom;
PSF(isnan(PSF))=0;
PSF(isinf(PSF))=0;
freq_axis=linspace(-1/1e-6,1/1e-6,size(PSF,1));
 PSF_self=cat(2,PSF(:,2,1),PSF(:,3,2),PSF(:,4,3));
  pl(freq_axis,fftshift(PSF_self),[-1.e6 1.e6])

PSF_FT=PSF;


G_FT=fft(G_pad,FFTLength,1);

%Spherical harmonics selector:
SPH_Sel=2:4; % X,Y,Z by default
if(parameters.isB0Term) 
    SPH_Sel = 1:4; 
end
if(parameters.isHigherOrderPSFCorr) 
    SPH_Sel =1:size(PSF,2); 
end

G_corr=zeros([size(G,1) length(SPH_Sel)   size(G,3)]);
for cintlv=1:size(G,3)

GX_effects=ifft(bsxfun(@times,G_FT(:,1,cintlv),PSF_FT(:,SPH_Sel,1)),[],1);
GY_effects=ifft(bsxfun(@times,G_FT(:,2,cintlv),PSF_FT(:,SPH_Sel,2)),[],1); 
GZ_effects=ifft(bsxfun(@times,G_FT(:,3,cintlv),PSF_FT(:,SPH_Sel,3)),[],1);

%remote all cross term compenents when cross term effects are not required
if(~parameters.isCrossTermPSFCorr)
    GX_effects(:,[2,3]+parameters.isB0Term)=0;
    GY_effects(:,[1 3]+parameters.isB0Term)=0;
    GZ_effects(:,[1, 2]+parameters.isB0Term)=0;   
end
%remove padding
% G_corr(:,:,cintlv)=real(GX_effects(1:size(G,1),:)+GY_effects(1:size(G,1),:)+GZ_effects(1:size(G,1),:));
G_corr(:,:,cintlv)=real(GX_effects(ZeroPadLength+(1:size(G,1)),:)+GY_effects(ZeroPadLength+(1:size(G,1)),:)+GZ_effects(ZeroPadLength+(1:size(G,1)),:));
end
end
