function [allkernel,allrawkernel] = calibSPIRiT3D(kCalib, kSize, CalibTyk,sampling)


Ncha=size(kCalib,4);
if nargin < 4
	sampling = ones([kSize,Ncha]);
end

[AtA] = dat2AtA3D(kCalib,kSize);
%AtA is correlation matrix of all kernel points across all coils calculted
%by all all slided kernels
%size(AtA) = [prod(kSize)*nCha   x    prod(kSize)*nCha ]



allkernel=zeros([kSize Ncha Ncha]);
if(nargout==2)
allrawkernel=zeros([kSize Ncha Ncha]);
end

for coil=1:Ncha
sampling = ones([kSize,Ncha]);
dummyK = zeros(kSize(1),kSize(2),kSize(3),Ncha); dummyK((end+1)/2,(end+1)/2,(end+1)/2,coil) = 1;
idxY = find(dummyK);
sampling(idxY) = 0;
idxA = find(sampling);

Aty = AtA(:,idxY); Aty = Aty(idxA);
AtAtmp = AtA(idxA,:); AtAtmp =  AtAtmp(:,idxA);

kernel = zeros(size(sampling));

CalibTyk_curr = norm(AtAtmp,'fro')/size(AtAtmp,1)*CalibTyk;

rawkernel = mldivide((AtAtmp + eye(size(AtAtmp))*CalibTyk_curr),Aty);
kernel(idxA) = rawkernel; 

allkernel(:,:,:,:,coil) = kernel;
if(nargout==2)
allrawkernel(:,:,:,:,coil)=rawkernel;
end

end
end

