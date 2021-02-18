function [CoilMaps]=GetCoilMaps(ref,mode,imSize)
% Calculate coils maps


CoilMaps=CalESpiritCoilMaps(ref,imSize);



end


function CoilMaps=CalESpiritCoilMaps(ref,imSize)

if(nargin<2)
    sz=size(ref);
    imSize=[sz(1:3)];
end

[sx,sy,sz,Nc] = size(ref);

fprintf("Cpar(%d):",sz)
if(sz>1) % make it seperable along PArtition direction
    ref=fftshift(fft(ref,imSize(3),3),3);
end

CoilMaps=zeros([imSize Nc],'double');
for Cpar=1:sz
    fprintf(1,'\b%d',Cpar);
ncalib = 24; % use 24 calibration lines to compute compression
ksize = [6,6]; % kernel size

DATA=squeeze(ref(:,:,Cpar,:));
% Threshold for picking singular vercors of the calibration matrix
% (relative to largest singlular value.
eigThresh_1 = 0.02;
% threshold of eigen vector decomposition in image space. 
eigThresh_2 = 0.95;

% crop a calibration area
calib = crop(DATA,[ncalib,ncalib,Nc]);

% compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels
[k,S] = dat2Kernel(calib,ksize);
idx = max(find(S >= S(1)*eigThresh_1));

[M,~] = kernelEig(k(:,:,:,1:idx),imSize(1:2));

CoilMaps(:,:,Cpar,:)=permute(M(:,:,:,end),[1 2 4 3]);
end
 fprintf('\n')
end