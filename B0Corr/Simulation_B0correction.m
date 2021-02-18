%% Simulate and test gridding without B0 correction
  parameter=struct('Ninterleaves', 1,'Resolution',1,'FOV',[1 1 1 1].*32,'MaxGradAmp',4,...
      'MaxSlewrate',20,'DwellTime',2e-6,'B0effects',0,'T2star',1);%fov,resolution
a=MriSim( parameter);
tic
FT = NUFFT(col(a.k.k),col(a.k.wi),1,0,[a.MatSize,a.MatSize], 2);
im=FT'*a.data(:);
FT_time=toc;
as(im)

%% Use MFI operator
tic
MFIOP=MFI(flip(flip(-2*pi*a.b0_map,1),2),a.t,FT,1);
imc_MFI=MFIOP'*a.data(:);
imc_MFI=imc_MFI*0.0125;% ~1/80
as(imc_MFI)
MFI_time=toc;
%% MTI test
tic
 MTIOP=MTI(flip(flip(-2*pi*a.b0_map,1),2),a.t,FT,1);
imc_MTI=squeeze(MTIOP'*(a.data(:).*FT.w).');
imc_MTI=imc_MTI*10;
MTI_time=toc;
as(imc_MTI)
%%
figure,
subplot(2,3,[1,4]),imagesc([],[],real(a.im_true),[0 2]),title('True image'),axis square,colorbar
subplot(232),imagesc([],[],real(im).*0.01,[0 2]),title(sprintf('Simple Gridding(%3.4f sec)',FT_time)),axis square,colorbar
subplot(233),imagesc(a.b0_map),title('B0 map(Hz)'),axis square,colorbar
subplot(235),imagesc([],[],real(imc_MFI),[0 2]),title(sprintf('MFI correction(%3.4f sec)',MFI_time)),axis square,colorbar
subplot(236),imagesc([],[],real(imc_MTI),[0 2]),title(sprintf('MTI correction(%3.4f sec)',MTI_time)),axis square,colorbar

%% use newly implemented forward operator
% W=ones(size(FT.w.'));
W=FT.w.';
nch=1;
beta=1e-6;
Fope=FOP(FT,MTIOP);
M=ones([prod(FT.imSize) 1]);
tic
%data
y=a.data.';
norm_old=-1;
clear normvec
% %initial estimate
  f0=MTIOP'*(y.*W);
%     f0=zeros(size(f0));
 f=zeros(numel(f0),20);
f(:,1)=f0(:);
for i=1:20
r=(y-Fope*reshape(f(:,i),nch,FT.imSize(1),FT.imSize(2))); % residual
% figure(11),plot(real(r)),hold on
g_new=col(Fope.*(r.*W))-col(beta*getRoughness(f(:,i),FT.imSize));

% normvec(i)=norm(g_new);
%   figure(12),plot(real(g_new)),hold on
  if (norm_old*1.1<norm(g_new)&& i>3)
%        break;
  else
      norm_old=norm(g_new);
  end

if(i==1)
    gamma=0;
    d=M.*g_new;
    q=Fope*reshape(d,1,FT.imSize(1),[]);
    alpha=(d'*g_new/((q.*W)*q')) *prod(FT.imSize);
else
    gamma=(g_new'*(M.*g_new))/(g_old'*(M.*g_old));
    d=M.*g_new+gamma*d;
%     d=d.*1e3;
    q=Fope*reshape(d,1,FT.imSize(1),[]);
    alpha=(d'*g_new)/((q.*W)*q'+beta*d'*col(getRoughness(d,FT.imSize)))*prod(FT.imSize);
end

f(:,i+1)=f(:,i)+alpha(:)*d;
r=r-alpha*q;
g_old=g_new;
end
f(:,i+1:end)=[];
iterative_time= toc;
as(reshape(f,FT.imSize(1),FT.imSize(2),[]),'complexSelect','m','colormap','parula')

% figure,plot(normvec)
%%
figure,
subplot(231),imagesc([],[],real(im)*0.01,[0 2]),title('Uncorrected image'),axis square,colorbar
subplot(232),imagesc(a.b0_map),title('B0 map(Hz)'),axis square,colorbar
subplot(233),imagesc([],[],real(imc_MFI),[0 2]),title(sprintf('MFI correction(%3.4f sec)',MFI_time)),axis square,colorbar
subplot(234),imagesc([],[],real(imc_MTI),[0 2]),title(sprintf('MTI correction(%3.4f sec)',MTI_time)),axis square,colorbar
subplot(235),imagesc([],[],(real(reshape(im_exact1,FT.imSize))),[0 2]),title('PCG exact(~310 sec)'),axis square,colorbar
subplot(236),imagesc([],[],(real(reshape(f(:,end-1),FT.imSize))),[0 2]),title(sprintf('Fast Iterative(~%2.2f+%2.2f sec)',MTI_time,iterative_time)),axis square,colorbar



%% Exact iterative and PCG 
E=a.getEncodingmatrix();
% don't do this more than imsize=100x100 
% tic,im_pcg=pcg(E'*E,E'*a.data(:));,toc
% as(reshape(im_pcg,[100 100]))


tic
clc
beta=1e-1;
%  M=ones([prod(FT.imSize) 1]);
%  M=FT.w(:);
  M=double(col(a.im_true>0));

%data
y=a.data(:);

%initial estimate
% f0=FT'*(y.*FT.w);
%   f0=E'*(y) -col(beta*getRoughness(E'*(y),FT.imSize));
    f0=zeros(FT.imSize);
f=zeros(numel(f0),50);
f(:,1)=f0(:);
for i=1:50
r=y-E*f(:,i);
g_new=E'*r(:)-col(beta*getRoughness(f(:,i),FT.imSize));
if(i==1)
    gamma=0;
    d=M.*g_new;
    q=E*d;
    alpha=d'*g_new/(q'*q);
else
    gamma=(g_new'*(M.*g_new))/(g_old'*(M.*g_old));
    d=M.*g_new+gamma*d;
    q=E*d;
    alpha=d'*g_new/(q'*(q)+beta*d'*col(getRoughness(d,FT.imSize)));
end
 
f(:,i+1)=f(:,i)+alpha(:)*d;
r=r-alpha*q;
g_old=g_new;
end
toc
as(reshape(f,FT.imSize(1),FT.imSize(2),[]))


%% Calculate MFI_weights
sig=permute(a.data(:),[2 1]);
[im_b0corr,MFI_Weights]=B0correctNUFFT(sig,-2*pi*a.b0_map,a.t*1e6,FT,1,[]);
as(im_b0corr)

%% Use the calculated MFI weights
[im_b0corr]=B0correctNUFFT(sig,-2*pi*a.b0_map,a.t*1e6,FT,1,MFI_Weights);
as(im_b0corr)

%% suporting fucntions
function R=getRoughness(omega,im_size) %
%get roughness by calculating the finite difference of neighboring pixel in
%all dimension
% for 2D the kernel h is
%  0 -1  0
% -1  4 -1
%  0 -1  0
% The kerenel is seperable with multiple 1D  kernel [-1 2 -1]

h=[-1;2;-1];
omega=reshape(omega,im_size);
R1=imfilter(omega,h,'replicate','same','corr');
R2=imfilter(omega,permute(h,[2 1 3]),'replicate','same','corr');
% if(im_size(3)>1)
%     R3=imfilter(omega,permute(h,[2 3 1]),'replicate','same','corr');
%     R=0.5*(R1.^2+R2.^2+R3.^2);
% else
       R=0.5*(conj(R1).*R1+conj(R2).*R2);
%             R=0.5*(conj(R1+R2).*(R1+R2));
%       R=0.5*(R1+R2);
% end
% R=0.5*(R1(:)'*R2(:));
end



%% old trash
% %% MTI test: all steps
% b0map=-2*pi*a.b0_map;
% tk=a.t;
% obj.nLevels=ceil((2*max(abs(b0map(:)))*max(tk)/pi));
% obj.tau=max(tk)/(obj.nLevels+1);
% tic
% img_tau0 =FT'*(a.data(:));
% 
% 
% wn_tau_l=exp(-1i*obj.tau.* (b0map(:)*(0:obj.nLevels)));
% % MTI_weights=zeros([numel(b0map) obj.nLevels]);
% 
% actual=exp(1i*b0map(:)*tk');
% 
%  MTI_weights= mldivide(wn_tau_l'*wn_tau_l,wn_tau_l'*actual);
% % MTI_weights= mldivide(wn_tau_l,actual);
% 
% ksp_MTI=zeros([length(tk) obj.nLevels+1]);
% img_MTI=zeros([FT.imSize obj.nLevels+1]);
% for idx_freq=0:obj.nLevels
%     ksp_MTI(:,idx_freq+1) =FT*(img_tau0.*exp(1i*b0map*obj.tau*idx_freq));
%     img_MTI(:,:,idx_freq+1)= FT'*ksp_MTI(:,idx_freq+1);
% end
% 
% % sig_MTI=sum(MTI_weights'.*ksp_MTI,2);
% sig_MTI=sum(MTI_weights'.*ksp_MTI,2);
% imc_MTI=FT'*sig_MTI(:);
% toc
% as(imc_MTI)
% 
% 
% 
% %% memory hog way
% E=a.getEncodingmatrix();
% sig=E*a.im_true(:);
% % FT = NUFFT(col(a.k.k),col(a.k.wi),[a.MatSize,a.MatSize]*0.,[a.MatSize,a.MatSize]);
% im=FT'*sig;
% 
% %%
% %gpuNUFFT
% % FT = NUFFT(col(a.k.k),col(a.k.wi),1,0,[a.MatSize,a.MatSize], 2);
% 
% %Lustig NUFFT
% FT = NUFFT(col(a.k.k),col(a.k.wi),[a.MatSize,a.MatSize]*0.,[a.MatSize,a.MatSize]);
% 
% %% Iterative with MFI operator
% % here y is kspace data 
% % therefore r
% % 
% % and 
% % g_new,g_old should be im _size
% 
% clc
% beta=1e-1;
% 
% M=ones([prod(FT.imSize) 1]);
% 
% %data
% y=a.data(:);
% 
% %initial estimate
% f0=FT'*y;
% f0=zeros(size(f0));
% f=zeros(numel(f0),20);
% f(:,1)=f0(:);
% for i=1:20
% r=y-FT*f(:,i);
% g_new=col(FT'*r-beta*getRoughness(f(:,i),FT.imSize));
% if(i==1)
%     gamma=0;
%     d=M.*g_new;
%     q=FT*d;
%     alpha=d'*g_new/(q'*q);
% else
%     gamma=(g_new'*(M.*g_new))/(g_old'*(M.*g_old));
%     d=M.*g_new+gamma*d;
%     q=FT*d;
%     alpha=d'*g_new/(q'*q+beta*d'*col(getRoughness(d,FT.imSize)));
% end
% 
% f(:,i+1)=f(:,i)+alpha(:)*d;
% r=r-alpha*q;
% g_old=g_new;
% end
% 
% as(reshape(f,FT.imSize(1),FT.imSize(2),[]))
% %% Check the transpose  and forward operator
% MTIOP=MTI(2*pi*a.b0_map,a.t,FT,1);
% y=a.data(:);
% 
% Exact_tranpose=flip(flip(reshape(E'*(y.*FT.w),FT.imSize),4),5);
% MTI_tranpose= (MTIOP'*y);
% % as(cat(3,Exact_tranpose./100,MTI_tranpose))
% figure,subplot(211),plot(real([Exact_tranpose(:)./100 MTI_tranpose(:)]) ),legend('Exact','MTI'),xlim([1000 2500])
% subplot(212),plot(imag([Exact_tranpose(:)./100 MTI_tranpose(:)]) ),legend('Exact','MTI'),title('imag'),xlim([1000 2500])
% 
% disp(sum(Exact_tranpose./100-MTI_tranpose,'all'))
% 
% %% forward
% 
% sig_exact=E*a.im_true(:);
% sig_MTI= MTIOP*a.im_true(:);
% figure, plot(real([sig_exact./99 sig_MTI])),legend('excat','MTI')
% 
% figure, plot(imag([sig_exact./99  sig_MTI])),legend('excat','MTI'),title('imag')
% % disp(sum((sig_exact./99-sig_MTI_forww).^2,'all'))
% disp(sum((sig_exact./99-sig_MTI).^2,'all'))
% 
% %% Iterative with MTI operator
% % here y is kspace data 
% % therefore r
% % 
% % and 
% % g_new,g_old should be im _size
% 
% clc
% beta=1e-1;
% FT1=FT;
% FT1.w=ones(size(FT.w));
% MTIOP=MTI(flip(flip(-2*pi*a.b0_map,1),2),a.t,FT1,1);
% 
% M=ones([prod(FT.imSize) 1]);
% 
% %data
% y=a.data(:);
% 
% %initial estimate
% % f0=MTIOP'*y;
%  f0=zeros(size(f0));
% f=zeros(numel(f0),25);
% f(:,1)=f0(:);
% for i=1:25
% r=y-MTIOP*f(:,i);
% g_new=col(MTIOP'*r-beta*getRoughness(f(:,i),FT.imSize));
% if(i==1)
%     gamma=0;
%     d=M.*g_new;
%     q=MTIOP*d;
%     alpha=d'*g_new/(q'*q);
% else
%     gamma=(g_new'*(M.*g_new))/(g_old'*(M.*g_old));
%     d=M.*g_new+gamma*d;
%     q=MTIOP*d;
%     alpha=d'*g_new/(q'*q+beta*d'*col(getRoughness(d,FT.imSize)));
% end
% 
% f(:,i+1)=f(:,i)+alpha(:)*d;
% r=r-alpha*q;
% g_old=g_new;
% end
% 
% as(reshape(f,FT.imSize(1),FT.imSize(2),[]))
% %%
% % here y is kspace data 
% % therefore r
% % 
% % and 
% % g_new,g_old should be im _size
% 
% clc
% beta=1e0;
% 
% M=ones([prod(MFIOP.NUFFTOP.imSize) 1]);
% 
% %data
% y=a.data(:);
% 
% %initial estimate
% f0=MFIOP'*y;
% % f0=zeros(size(f0));
% f=zeros(numel(f0),20);
% f(:,1)=f0(:);
% for i=1:5
% r=y-MFIOP*f(:,i);
% g_new=col(MFIOP'*r-beta*getRoughness(f(:,i),MFIOP.NUFFTOP.imSize));
% if(i==1)
%     gamma=0;
%     d=M.*g_new;
%     q=MFIOP*d;
%     alpha=d'*g_new/(q'*q);
% else
%     gamma=(g_new'*(M.*g_new))/(g_old'*(M.*g_old));
%     d=M.*g_new+gamma*d;
%     q=MFIOP*d;
%     alpha=d'*g_new/(q'*q+beta*d'*col(getRoughness(d,MFIOP.NUFFTOP.imSize)));
% end
% 
% f(:,i+1)=f(:,i)+alpha(:)*d;
% r=r-alpha*q;
% g_old=g_new;
% end
% 
% as(reshape(f,MFIOP.NUFFTOP.imSize(1),MFIOP.NUFFTOP.imSize(2),[]))