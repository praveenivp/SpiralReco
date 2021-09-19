%% 
clc
% cd('S:\Subject_scan\20210521_subject3477_bssfp_sfmri')
FTOP=NUFFT(complex(1,1),1,0,[1 1]); %NUFFT class works
load('S:\Subject_scan\20210521_subject3477_bssfp_sfmri\processeddata\testdata_R3_set2_mid447_20210521.mat')

% load twix and undersampling 
% set 2 is in the test data
% twix=mapVBVD('meas_MID447_so_BSSFP_i21_R1_fmri_FID38133.dat');
% [FTOP,sig,adcTime]=getNUFFTOP(twix,3,1,0);

%% zeropadded recon
sig2=permute(sig,[1,3,4,2]);
im=FTOP'*double(sig2);
im=fftshift(ifft(ifftshift(im,3),[],3),3);
as(sum(im.*conj(permute(csm,[2 3 4 1])),4))
%% MFI
MFIOP=MFI3D(-1*fm,adcTime(:)*1e-6,FTOP,double(csm));
% MFIOP=MFI3D(0,adcTime(:)*1e-6,FTOP,double(csm)); %no B0
im_MFI=MFIOP'*double(permute(sig,[2,1,3,4]));
%   as(cat(4,im_MFI,sum(im.*conj(permute(flip(csm,40),[2 3 4 1])),4)))
%% MFI cgsense
clc
nIterCG=4;
limit=1e-6;
csm_sq = squeeze(sum(csm .* conj(csm),1)); csm_sq(csm_sq < eps) = 1;
M = spdiag(sqrt(csm_sq)); %Preconditioner

E_MFI=@(x,transp) myCGSENSE3D_MFI(x,MFIOP,double(csm),transp);
tic
[img_cgsense_MFI,flag,relres,iter,resvec] = lsqr(E_MFI, double(col(permute(sig,[2 1 3 4]))), limit,nIterCG,M);
rtime=toc
as(cat(4,reshape(img_cgsense_MFI,FTOP.imSize(1),FTOP.imSize(2),[]),im_MFI,sum(im.*conj(permute(csm,[2 3 4 1])),4)))

%% Try the same with Stackofspiral class
kxy=reshape(FTOP.st.om.',2,[],size(sig,3))/(2*pi);
kxy=repmat(kxy,[1 1 1 size(sig,4)]);
kz=-0.5:1/size(sig,4):0.5-1/size(sig,4);%linspace(-0.5,0.5,obj.SpiralPara.NPartitions);
kz=repmat(permute(kz(:),[2 3 4 1]),[1 size(kxy,2) size(kxy,3) 1]);
kxyz=cat(1,kxy,kz);
DCF3d=repmat(FTOP.w.^2,[1  1 size(sig,4)]);
N=FTOP.imSize(1);
SOSOP= StackofSpirals(kxyz,DCF3d,[N N size(sig,4)],permute(csm,[2 3 4 1]),...
                                'CompMode','CPU2DHybrid','precision','single');
SOSOPB0=StackofSpiralsB0(kxyz,DCF3d,[N N size(sig,4)],permute(csm,[2 3 4 1]),-1*fm,adcTime(:)*1e-6,...
                     'CompMode','CPU2DHybrid','precision','single',...
                     'Method','MFI');                            
im_zp=SOSOP'*permute(sig,[1 3 4 2]);
imcg=spiralCGSENSE(SOSOP,permute(sig,[1 3 4 2]),'maxit',10,'reg','none'); %cgsense
imcgB0=spiralCGSENSE(SOSOPB0,permute(sig,[1 3 4 2]),'maxit',10,'reg','none'); %cgsense
as(cat(4,im_zp,imcg,imcgB0,reshape(img_cgsense_MFI,FTOP.imSize(1),FTOP.imSize(2),[])));
%% MTI
% MTIOP=MTI_3D(-1*fm,adcTime(:)*1e-6,FTOP,double(csm));
% MTIOP=MTI_3D(0,adcTime(:)*1e-6,FTOP,flip(double(csm),4)); %no B0
im_MTI=MTIOP'*double(permute(sig,[2 1 3 4]));
  as(cat(4,im_MTI,sum(im.*conj(permute(csm,[2 3 4 1])),4)))
 %% MTI-cgsense
clc
nIterCG=4;
E_MTI=@(x,transp) myCGSENSE3D_MFI(x,MTIOP,double(csm),transp);
tic
[img_cgsense_MTI,flag,relres,iter,resvec] = lsqr(E_MTI, double(col(permute(sig,[2 1 3 4]))), limit,nIterCG,M);

rtime=toc
as(cat(4,reshape(img_cgsense_MTI,FTOP.imSize(1),FTOP.imSize(2),[]),im_MTI,sum(im.*conj(permute(csm,[2 3 4 1])),4)))


%%  cgsense no regularization 
E_FT=@(x,transp) myCGSENSE3D(x,FTOP,double(csm),transp);
reg_out=[];
limit=1e-9;

csm_sq = squeeze(sum(csm .* conj(csm),1)); csm_sq(csm_sq < eps) = 1;
M = spdiag(sqrt(csm_sq)); %Preconditioner

sig2=permute(sig,[1,3,2,4]);
nIterCG=20;
tic
[img_cgsense_NUFFT,flag,relres,iter,resvec] = lsqr(E_FT, double([sig2(:); reg_out]), limit,nIterCG,M);
img_cgsense_NUFFT=reshape(img_cgsense_NUFFT,FTOP.imSize(1),FTOP.imSize(2),[]);
rtime=toc
as(img_cgsense_NUFFT)

%% tikhonov regularization
E_reg=@(x,transp) myCGSENSE3D_reg(x,FTOP,double(csm),1000,transp);
reg_out=col(zeros([192 192 20]));
tic
[img_cgsense_NUFFT_reg,flag,relres,iter,resvec_] = lsqr(E_reg, double([sig2(:); reg_out]), limit,nIterCG,M);
img_cgsense_NUFFT_reg = reshape(img_cgsense_NUFFT_reg,size(csm,2),size(csm,3),size(csm,4));
rtime=toc
as(cat(4,img_cgsense_NUFFT,img_cgsense_NUFFT_reg,img_cgsense_NUFFT-img_cgsense_NUFFT_reg))

figure,plot(resvec),hold on,plot(resvec_tk),legend('No regularization','Tikonov reg'),title('residual vector CGSENSE @ lambda=1000'),xlabel('Iterations')
%%

function outp =  myCGSENSE3D_reg(inp,NUFFT_obj,Coilsens,lambda,transpose_indicator)

% Coilsens [cha x colx Lin x Par]
% inp[nFE x nIntlv x nCha x npar ]
% scale = sqrt(prod(prod(nufft_st.Kd))/numel(weights(:)));
scale=1; %acceleration factor
Npar=size(Coilsens,4);
if (strcmp(transpose_indicator,'transp'))
       imsize=[NUFFT_obj.dataSize(1),NUFFT_obj.dataSize(2),size(Coilsens,1),size(Coilsens,4)];
      inp=reshape(inp(1:prod(imsize)),imsize);
      outp = NUFFT_obj'*double(inp);
      outp=(fftshift(ifft(ifftshift(outp,4),[],4),4))*(sqrt(Npar).*sqrt(scale));
      outp= double(col(squeeze(sum(outp.*conj(permute(Coilsens,[2,3,1,4])),3))));
      if(numel(inp)>prod(imsize)) %reg
         outp=outp+lambda*double(inp((prod(imsize)+1):end));
      end
elseif (strcmp(transpose_indicator, 'notransp'))
    inp=reshape(inp,[NUFFT_obj.imSize(1),NUFFT_obj.imSize(2),size(Coilsens,4)]);
    reg=double(inp(:).*conj(inp(:)))*lambda;
    inp=permute(bsxfun(@times,inp,permute(Coilsens,[2,3,4,1])),[1 2,4,3]); % colxlinxCoilxpar
    inp=(fftshift(fft(ifftshift(inp,4),[],4),4))/(sqrt(Npar).*sqrt(scale));
    outp = NUFFT_obj*double(inp);
    outp=[double(outp(:)) ;reg];
else
    error('Transpose flag not appropriately defined');
end

end

function outp =  myCGSENSE3D(inp,NUFFT_obj,Coilsens,transpose_indicator)

% Coilsens [cha x colx Lin x Par]
% inp[nFE x nIntlv x nCha x npar ]
% scale = sqrt(prod(prod(nufft_st.Kd))/numel(weights(:)));
scale=1; %acceleration factor
Npar=size(Coilsens,4);
if (strcmp(transpose_indicator,'transp'))
      inp=reshape(inp,NUFFT_obj.dataSize(1),NUFFT_obj.dataSize(2),size(Coilsens,1),size(Coilsens,4));
      outp = NUFFT_obj'*double(inp);
      outp=(fftshift(ifft(ifftshift(outp,4),[],4),4))*(sqrt(Npar).*sqrt(scale));
      outp= double(col(squeeze(sum(outp.*conj(permute(Coilsens,[2,3,1,4])),3))));

elseif (strcmp(transpose_indicator, 'notransp'))
    inp=reshape(inp,[NUFFT_obj.imSize(1),NUFFT_obj.imSize(2),size(Coilsens,4)]);
    inp=permute(bsxfun(@times,inp,permute(Coilsens,[2,3,4,1])),[1 2,4,3]); % colxlinxCoilxpar
    inp=(fftshift(fft(ifftshift(inp,4),[],4),4))/(sqrt(Npar).*sqrt(scale));
    outp = NUFFT_obj*double(inp);
    outp=double(outp(:));
else
    error('Transpose flag not appropriately defined');
end

end

function outp =  myCGSENSE3D_MFI(inp,MFIOP,Coilsens,transpose_indicator)

% Coilsens [cha x colx Lin x Par]
% inp[nFE x nIntlv x nCha x npar ]
if (strcmp(transpose_indicator,'transp'))
      inp=reshape(inp,size(Coilsens,1),MFIOP.NUFFTOP.dataSize(1),MFIOP.NUFFTOP.dataSize(2),size(Coilsens,4));
      outp = MFIOP'*double(inp);
        outp=col(double(outp));
elseif (strcmp(transpose_indicator, 'notransp'))
    inp=reshape(inp,[MFIOP.NUFFTOP.imSize(1),MFIOP.NUFFTOP.imSize(2),size(Coilsens,4)]);
    outp = MFIOP*double(inp);
    outp=double(outp(:));
else
    error('Transpose flag not appropriately defined');
end

end
