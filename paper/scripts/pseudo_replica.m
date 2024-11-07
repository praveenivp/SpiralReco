
%% Create calib data, script for bart and registered fieldmap
addpath(genpath('/ptmp/pvalsala/MATLAB'))
Measpath='/ptmp/pvalsala/DKSR-UFYK';
dir_st=dir(fullfile(Measpath,'TWIX','*peSpiral_R4*.dat'));
clear commands;
for i=7%1:(length(dir_st)-2)
    disp(dir_st(i).name)
    [~,measID]=regexp(dir_st(i).name,'\S*#M(\d{2,}+)\S*','match','tokens');
     csmfile=fullfile('/ptmp/pvalsala/DKSR-UFYK/SNR/dep_hansen',sprintf('fm_csm_MeasUID%d.mat',str2double(measID{1})));
      load(csmfile,'csm','fm_interp')
     SpObj=SpiralReco(fullfile(dir_st(i).folder,dir_st(i).name),'RepSel',2,...
    'doPAT','CGSENSE','csm',csm,'maxit',10,'reg','Tikhonov','reg_lambda',1e-3,'compMode','CPU3D',...
    'fm',-1*fm_interp,'doB0Corr','MTI','doDCF','Jackson','precision','double','NormNoiseData',false);


R=[1 1];
csm_chaLast=permute(csm,[2 3 4 1]);
[sosop,sig]=undersampleSpiral(SpObj,csm_chaLast,R);
im_formation_fun=@(sig)  spiralCGSENSE(sosop,sig,...
    'maxit',10,'tol',1e-6,'reg','Tikhonov','lambda',1e-3);
% im=im_formation_fun(permute(sig,[1 2 3 4]));
% as(im./sqrt(prod(R)))
   noise                = permute(SpObj.twix.noise(''),[1,3,4,5,6,2]);
noise                = reshape(noise,[],size(noise,ndims(noise)));
R                    = (noise.'*conj(noise))./(2*size(noise,1));
R=R./mean(abs(diag(R)));
noise_cov_matrix=0.5*(R+R.');

[snr,g,noise_psf,~,basline] = ismrm_pseudo_replica(sig, im_formation_fun, 128,noise_cov_matrix);
Nrep=128;
flags=SpObj.flags;
    save(sprintf('M%d_%s.mat',str2double(measID{1}),SpObj.twix.hdr.Config.ProtocolName),'snr','g','noise_psf','basline','Nrep','flags')
end

%%

[snr,g,noise_psf,~,basline] = ismrm_pseudo_replica(sig, im_formation_fun, 32,noise_cov_matrix);
    save(sprintf('pseudo_replica_M%d.mat',measID),'snr','g','noise_psf','basline')
%% reco
csm=permute(Coilmaps(:,:,:,:,1),[4 2 1 3 5]);
r=SpiralReco(twix,'RepSel',2,...
    'doPAT','CGSENSE','csm',csm,'maxit',10,'reg','Tikhonov','reg_lambda',1e-3,'compMode','CPU2DHybrid',...
    'fm',[],'doB0Corr','none','doDCF','Jackson','precision','double');
%% reco
%  load('coilmaps_hansen.mat')
csm=permute(Coilmaps(:,:,:,:,1),[4 2 1 3 5]);

%%
R=[1 1];
csm_chaLast=permute(csm,[2 3 4 1]);
[sosop,sig]=undersampleSpiral(rh,csm_chaLast,R);
im=sosop'*sig;
as(im)
%%
R=[1 1];
csm_chaLast=permute(csm,[2 3 4 1]);
[sosop,sig]=undersampleSpiral(rh,csm_chaLast,R);
im_formation_fun=@(sig)  spiralCGSENSE(sosop,sig,...
    'maxit',10,'tol',1e-6,'reg','Tikhonov','lambda',1e-3);
im=im_formation_fun(permute(sig,[1 2 3 4]));
as(im./sqrt(prod(R)))

%%

noise                = permute(twix.noise(''),[1,3,4,5,6,2]);
noise                = reshape(noise,[],size(noise,ndims(noise)));
R                    = (noise.'*conj(noise))./(2*size(noise,1));
R=R./mean(abs(diag(R)));
noise_cov_matrix=0.5*(R+R.');
n_white = complex(randn(size(sig)),randn(size(sig)));

n_color=chol(noise_cov_matrix)*reshape(n_white,[],size(n_white,ndims(n_white))).';

n_color=reshape(n_color.',size(n_white));

%% pseudo replica
[snr,g,noise_psf,img_noise_rep] = ismrm_pseudo_replica(sig, im_formation_fun, 10,noise_cov_matrix);
as(snr)

%%
for cIter=1:2:40
    fprintf('Curr iter :%d\n',cIter)
    im_formation_fun=@(sig)  spiralCGSENSE(sosop,sig,...
    'maxit',cIter,'tol',1e-6,'reg','Tikhonov','lambda',1e-3);
    [snr,g,noise_psf,img_noise_rep,basline] = ismrm_pseudo_replica(sig, im_formation_fun, 128,noise_cov_matrix);
    save(sprintf('M113_cIter_%d.mat',cIter),'snr','g','noise_psf','basline')
end


%% plotting
mask=imerode(squeeze(sos(csm,1)>0),strel('sphere',10));
snr_all=[];g_all=[]; baseline_all=[];
for cIter=1:2:40
    fprintf('Curr iter :%d\n',cIter)
  
    load(sprintf('M113_cIter_%d.mat',cIter),'snr','g','basline')
    
    snr_all=cat(4,snr_all,snr);
    g_all=cat(4,g_all,g);
    baseline_all=cat(4,baseline_all,basline);
    
end

as(snr,'title','iter 39','select',':,:,10')

%%
slcSel=4:12;

nanMask=single(mask);
nanMask(nanMask==0)=nan;

g_plot=g_all.*nanMask;
g_plot=g_plot(:,:,slcSel,:);
g_plot=reshape(g_plot(~isnan(g_plot)),[],size(g_plot,4));
snr_plot=snr_all.*nanMask;
snr_plot=snr_plot(:,:,slcSel,:);
snr_plot=reshape(snr_plot(~isnan(snr_plot)),[],size(snr_plot,4));


imdiff_plot=diff(baseline_all.*nanMask,4);
imdiff_plot=imdiff_plot(:,:,slcSel,:);
imdiff_plot=reshape(imdiff_plot(~isnan(imdiff_plot)),[],size(imdiff_plot,4));

figure(3),clf,subplot(411),violin(g_plot);
xticks(1:20);
xticklabels(strsplit(num2str(1:2:40),' '))

subplot(412),
violin(abs(snr_plot)./sqrt(2.8));
xticks(1:20);
xticklabels(strsplit(num2str(1:2:40),' '))
ylim([0 300])
hold on, plot(sum(abs(snr_plot)./sqrt(2.8),1),'LineWidth',2)
xlabel('CG Iteration')



subplot(413),
plot(1:2:40,sum(abs(imdiff_plot),1))
% violin(abs(imdiff_plot));
% xticks(1:20);
% xticklabels(strsplit(num2str(1:2:40),' '))
% % ylim([0 300])
xlabel('CG Iteration')

subplot(414),plot(1:2:40,mean(abs(snr_plot)./sqrt(2.8),1),'LineWidth',2)

%%
function [sosop,sig]=undersampleSpiral(r1,csm,R)

sig=r1.sig(:,:,1:R(1):end,1:R(2):end,1);
sig=permute(sig,[2 3 4 1]);
ktraj=r1.KTraj(:,:,1:R(1):end,1:R(2):end,1);
DCF=r1.DCF(:,1:R(1):end,1:R(2):end,1);

            kmax=2*pi*(0.5/(r1.SpiralPara.Resolution*1e-3));
%             kmax=[kmax;kmax;(2*pi*0.5)/(1e-3*obj.SpiralPara.slice{1}.FOV_PRS(3)/obj.SpiralPara.NPartitions)];
sosop=r1.NUFFT_obj;
% sosop=StackofSpirals(ktraj./(2*kmax),DCF,[r1.NUFFT_obj.imSize],csm,...
% 'CompMode','CPU3D','precision','single');
end

function [snr,g,noise_psf,img_noise_rep,baseline] = ismrm_pseudo_replica(in, image_formation_func, reps,Noise_cov)
baseline = image_formation_func(in);

parfor r=1:reps,
    fprintf('Running pseudo replica %d/%d\n',r,reps);
    n_white = complex(randn(size(in)),randn(size(in)));    
    n_color=chol(Noise_cov)*reshape(n_white,[],size(n_white,ndims(n_white))).';
    n_color=reshape(n_color.',size(n_white));
    
    s = in + n_color;
    tmp = image_formation_func(s);
    img_noise_rep(:,r) = tmp(:);
end

img_noise_rep = reshape(img_noise_rep,[size(baseline),reps]);
rep_dim = length(size(img_noise_rep));

g = std(abs(img_noise_rep + max(abs(img_noise_rep(:)))),[],rep_dim); %Measure variation, but add offset to create "high snr condition"
g(g < eps) = 1;
snr = mean(img_noise_rep,ndims(img_noise_rep))./g;

img_noise_rep = img_noise_rep - repmat(baseline,[ones(1,length(size(baseline))) reps]);


noise_psf=0;
% ftdims = 1:(rep_dim-1);
% noise_psf = ismrm_transform_kspace_to_image(mean(abs(ismrm_transform_image_to_kspace(img_noise_rep,ftdims)).^2,3));

end

%
function [snr,g,noise_psf] = ismrm_pseudo_replica2(in, image_formation_func, reps,Noise_cov)
baseline = image_formation_func(in);

for r=1:reps,
    fprintf('Running pseudo replica %d/%d\n',r,reps);
%     n = complex(randn(size(in)),randn(size(in)));
        n_white = 1e-5*complex(randn(size(in)),randn(size(in)));    
    n_color=chol(Noise_cov)*reshape(n_white,[],size(n_white,ndims(n_white))).';
    n_color=reshape(n_color.',size(n_white));
    s = in + n_color;
    tmp = image_formation_func(s);
    img_noise_rep(:,r) = tmp(:);
end

img_noise_rep = reshape(img_noise_rep,[size(tmp),reps]);
rep_dim = length(size(img_noise_rep));

g = std(abs(img_noise_rep + max(abs(img_noise_rep(:)))),[],rep_dim); %Measure variation, but add offset to create "high snr condition"
g(g < eps) = 1;
snr = mean(img_noise_rep,3)./g;

img_noise_rep = img_noise_rep - repmat(baseline,[ones(1,length(size(baseline))) reps]);

ftdims = 1:(rep_dim-1);
noise_psf = ismrm_transform_kspace_to_image(mean(abs(ismrm_transform_image_to_kspace(img_noise_rep,ftdims)).^2,3));

end