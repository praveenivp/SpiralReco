function [G_corr] = GIRF_correction_Freq(G_xyz,PSF,varargin)
%G - uncorrected gradient waveform
    %G_xyz: (N x Axis x NIntlv) matrix correspond to each axis
    %G and PSF should have same sampling time.
 %PSF is in Frequency domain 
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
         
%make input gradient matrix into a 3D matrix
G_xyz=reshape(G_xyz,size(G_xyz,1),size(G_xyz,2),[]);
 
%Spherical harmonics selector:
SPH_Sel=1:4; % B0,X,Y,Z by default
if(parameters.isHigherOrderPSFCorr) 
    SPH_Sel =1:size(PSF,2); 
end



% adapted from % https://github.com/MRI-gradient/girf : 
%  Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
dt=10e-6; %always 10us gradient raster
fmax=1/dt;
dFreq=fmax/size(PSF,1);
faxis_PSF = dFreq*((0:(size(PSF,1)-1))-floor(size(PSF,1)/2)) ;

% Set length of prediction
TH = 1/dFreq; 
TO = size(G_xyz,1)*dt;
TZF = min(TH,TO);

% Prepare input
nZF = ceil(TZF/dt);
Gxyz_pad=padarray(G_xyz,[nZF 0 0],0,'both');
tIn=dt*(0:(size(G_xyz,1)-1));
tIn = [tIn(1)+dt*(-nZF:-1)'; tIn'; tIn(end)+dt*(1:nZF)'];
fIn = time2freq(tIn);
IN = fftshift(fft(Gxyz_pad,[],1),1)*dt;
nsIn = length(tIn);
% T_In = tIn(end) - tIn(1);

% Interpolate GIRF onto input grid
HIp = interp1(faxis_PSF,PSF,fIn,'pchip');
HIp(isnan(HIp)) = 0;


%     hTime = real(ifft(ifftshift(HIp,1)));
%     tStart = 1e-3; % impulse response startup time to use
%     tSettle = TZF; % impulse response settling time to use
%     tEndInd = ceil(tSettle/dt);
%     tStartInd = ceil(tStart/dt);
%     hTime([tEndInd+1:nsIn-tStartInd-1],:,:) = 0;
%     HNew = fftshift(fft(hTime),1);
%     HNew(isnan(HNew)) = 0;
%     HIp = HNew;

 gc=zeros([length(fIn)  max(SPH_Sel) size(G_xyz,3) 3]);

    for axis=1:3
     for sph=SPH_Sel
         for cintlv=1:size(G_xyz,3)
            gc(:,sph,cintlv,axis)=HIp(:,sph,axis).*IN(:,axis,cintlv);
         end
     end
    end
     
gc = ifft(ifftshift(gc,1))/dt;
gc = real(gc);

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

if(strcmpi(parameters.ConvType,'same'))
gc(1:nZF,:,:,:)=[];
gc(end-nZF+1:end,:,:,:)=[];
end

if(parameters.isB0Term)
    G_corr=gc;
else
    G_corr=gc(:,2:end,:);
end
end


function [f, df, f_max] = time2freq(t)
% Function to compute frequency vector corresponding to a time-vector 
% (assuming fftshift)
%
% USE
% [f, df, f_max] = time2freq(t)
%
% IN
% t         [n_samples x 1] time vector [s]
% 
% OUT
% f         [n_samples x 1] frequency vector [Hz]
% df        frequency resolution [Hz]
% f_max     bandwidth [Hz]
%
% Author:   Johanna Vannesjo (johanna.vannesjo@gmail.com)
% Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
%               2016 FMRIB centre, University of Oxford
%
% This file is part of a code package for GIRF computation and application. 
% The package is available under a BSD 3-clause license. Further info see:
% https://github.com/MRI-gradient/girf
%

nrs = length(t);
dt = t(2)-t(1);
f_max = 1/dt;
df = f_max/nrs;
f = ([0:nrs-1]'-floor(nrs/2))*df;
end
