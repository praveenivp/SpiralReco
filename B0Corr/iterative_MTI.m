function [im,im_coil]=iterative_MTI(MTIOP,coilData,beta,niter)
%
%[im,im_coil]=iterative_MTI(MTIOP,coilData,beta,niter);
%
%Funciton to perform iterative B0 correction with Multi-time interpolation.
%only works for 2D at the moment.
%
%
%INPUTS:
% MTIOP- Multi-time interpolation object from MTI.m class
%CoilData- 
%
%
%Dependecies:
% ESPIRIT package for NUFFT
% MTI operator from SpiralReco package
%
%Reference:
%Sutton, B. P., Noll, D. C., & Fessler, J. A. (2003). 
% Fast, iterative image reconstruction for MRI in the presence of field inhomogeneities.
% IEEE Transactions on Medical Imaging, 22(2), 178–188. 
% https://doi.org/10.1109/TMI.2002.808360
%
%praveenivp

switch nargin
    case {0,1}
        error('need MTIOP and CoilData')
    case 2
        beta=1e-1;
        niter=5;
    case 3
        niter=5;
    case 4
    otherwise
        error('Too many input arguments: [im,im_coil]=iterative_MTI(MTIOP,coilData,beta,niter)')
end

nCh=size(coilData,1);


FT=MTIOP.NUFFTOP;

M=ones([nCh prod(FT.imSize) ]);
tic

norm_old=-1;
% %initial estimate
MTIOP1=MTIOP;
MTIOP1.CoilSens=1;
f0=MTIOP1'*(bsxfun(@times,coilData,(MTIOP.NUFFTOP.w).'));
im_coil=zeros(nCh,prod(MTIOP.NUFFTOP.imSize),niter+1,'single');
im_coil(:,:,1)=f0;

for i=1:niter
    resi=coilData-MTIOP1*im_coil(:,:,i); % residual
    % figure(11),plot(real(r)),hold on
    temp2=getRoughnessMC(im_coil(:,:,i),[nCh FT.imSize]);
    g_new=MTIOP1'*resi-beta*temp2;
    % normvec(i)=norm(g_new);
    %   figure(12),plot(real(g_new)),hold on
    if (norm_old*1.1<norm(g_new(:))&& i>3)
        %        break;
    else
        norm_old=norm(g_new(:));
    end
    if(i==1)
        gamma=0;
        d=M.*g_new;
        q=MTIOP1*reshape(d,[nCh FT.imSize]);
        alpha=(sum(conj(d).*g_new,2)./sum(conj(q).*q,2)) *prod(FT.imSize);
    else
        gamma=(col(g_new)'*col(M.*g_new))/(col(g_old)'*col(M.*g_old));
        gamma=sum(conj(g_new).*(M.*g_new),2)./sum(conj(g_old).*(M.*g_old),2);
        d=M.*g_new+bsxfun(@times,d,gamma);
        %     d=d.*1e3;
        q=MTIOP1*d;
        
        alpha=(sum(conj(d).*g_new,2)./(sum(conj(q).*q,2)+  ...
            beta*sum(conj(d).*getRoughnessMC(im_coil(:,:,i),[nCh FT.imSize]),2)))*prod(FT.imSize);
        %     alpha=(d'*g_new)/(q*q'+beta*d'*col(getRoughness(d,FT.imSize)))*prod(FT.imSize);
    end
    
    im_coil(:,:,i+1)=im_coil(:,:,i)+bsxfun(@times,d,alpha);
    resi=resi-bsxfun(@times,q,alpha);
    g_old=g_new;
end

im_coil=reshape(im_coil,[],MTIOP.NUFFTOP.imSize(1),MTIOP.NUFFTOP.imSize(1),niter+1);
if(~isempty(MTIOP.CoilSens))
   
    im=squeeze(sum(bsxfun(@times,MTIOP.CoilSens,im_coil),1));
else
    im=im_coil;
end

fprintf('Finished %d iterations in %.2f seconds\n',niter,toc)
end

function R=getRoughnessMC(omega,im_size) %
%get roughness by calculating the finite difference of neighboring pixel in
%all dimension
% for 2D the kernel h is
%  0 -1  0
% -1  4 -1
%  0 -1  0
% The kerenel is seperable with multiple 1D  kernel [-1 2 -1]

omega=reshape(omega,im_size);
omega=permute(omega,[2,3,1]);
h=[-1;2;-1];

R1=imfilter(omega,h,'replicate','same','corr');
R2=imfilter(omega,permute(h,[2 1 3]),'replicate','same','corr');

R=0.5*(conj(R1).*R1+conj(R2).*R2);
R=reshape(R,[],im_size(1)).';
end
