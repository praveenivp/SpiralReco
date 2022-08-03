%%
cd('S:\Phantom\20220215_2dtest')
r1=SpiralReco('meas_MID00559_FID08369_peSpiral_R1.dat','doDCF','voronoi');
R=[1 1];
csm_chaLast=permute(csm,[2 3 4 1]);
[sosop,sig]=undersampleSpiral(r1,csm_chaLast,R);
im=sosop'*sig;
as(im)
%%
R=[3 1];
csm_chaLast=permute(csm,[2 3 4 1]);
[sosop,sig]=undersampleSpiral(r1,csm_chaLast,R);
im_formation_fun=@(sig)  spiralCGSENSE(sosop,sig,...
    'maxit',10,'tol',1e-6,'reg','none',...
    'lambda',1e-2);
im=im_formation_fun(sig);
as(im./sqrt(prod(R)))



%% pseudo replica
[snr,g,noise_psf,img_noise_rep] = ismrm_pseudo_replica(sig*1e3, im_formation_fun, 50);
as(snr)

%%
function [sosop,sig,img_noise_rep]=undersampleSpiral(r1,csm,R)

sig=r1.sig(:,:,1:R(1):end,1:R(2):end,1);
sig=permute(sig,[2 3 4 1]);
ktraj=r1.KTraj(:,:,1:R(1):end,1:R(2):end,1);
DCF=r1.DCF(:,1:R(1):end,1:R(2):end,1);

            kmax=2*pi*(0.5/(r1.SpiralPara.Resolution*1e-3));
%             kmax=[kmax;kmax;(2*pi*0.5)/(1e-3*obj.SpiralPara.slice{1}.FOV_PRS(3)/obj.SpiralPara.NPartitions)];

sosop=StackofSpirals(ktraj./(2*kmax),DCF,[r1.NUFFT_obj.imSize(1:2) 0],csm,...
'CompMode','CPU2DHybrid','precision','single');
end

function [snr,g,noise_psf] = ismrm_pseudo_replica(in, image_formation_func, reps)
baseline = image_formation_func(in);

parfor r=1:reps,
    fprintf('Running pseudo replica %d/%d\n',r,reps);
    n = complex(randn(size(in)),randn(size(in)));
    s = in + n;
    tmp = image_formation_func(s);
    img_noise_rep(:,r) = tmp(:);
end

img_noise_rep = reshape(img_noise_rep,[size(baseline),reps]);
rep_dim = length(size(img_noise_rep));

g = std(abs(img_noise_rep + max(abs(img_noise_rep(:)))),[],rep_dim); %Measure variation, but add offset to create "high snr condition"
g(g < eps) = 1;
snr = mean(img_noise_rep,3)./g;

img_noise_rep = img_noise_rep - repmat(baseline,[ones(1,length(size(baseline))) reps]);


noise_psf=0;
% ftdims = 1:(rep_dim-1);
% noise_psf = ismrm_transform_kspace_to_image(mean(abs(ismrm_transform_image_to_kspace(img_noise_rep,ftdims)).^2,3));

end