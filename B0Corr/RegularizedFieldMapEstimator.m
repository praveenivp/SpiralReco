function [FieldMap,Fm_steps,g_steps]=RegularizedFieldMapEstimator(im,TE,Fm_initial,beta,maxIter)
%function [FieldMap,Fm_steps,g_steps]=RegularizedFieldMapEstimator(im,TE,Fm_initial,beta,maxIter)
%
%INPUTS:
%im: complex Echo images
%    4D matrix: first three are physical dimension and 4th is echo dim
%TE:  Echo times corresponds to 4th dimension of im
%      1D array : in Seconds
%Fm_initial: Noisy field map in rad/s (optional)
%     calculted from Echo images (im) when not provided.
%beta: Regularization paramter (optional)
%      scalar value
%maxIter: maximum iterations  (optional)
%         scalar
%
%OUTPUTS:
%FieldMap: 3D matrix with B0 fieldmap (in rad/s)
% Following outputs are only for debugging and need large memory
% Fm_steps : 4D matrix with field maps at each iterations(rad/s)
% g_steps  : 4D matrix with calculated gradients for each iteration.
%
% author:praveen.ivp@gmail.com
%
%Reference
% Funai, A. K., Fessler, J. A., Yeo, D. T. B., Noll, D. C., & Olafsson, V. T. (2008).
% Regularized field map estimation in MRI. IEEE Transactions on Medical Imaging,
% 27(10), 1484–1494.
% https://doi.org/10.1109/TMI.2008.923956

% reshaping inputs
sz=size(im);
TE=TE(:);
im=reshape(im,[],length(TE));
Fm_initial=Fm_initial(:);


%Input parsing
switch(nargin)
    case {1,2}
        error('first 2 arguments necessary');
    case 3
        Fm_initial=[];
        beta=1e-1;
        maxIter=10;
    case 4
        beta=1e-1;
        maxIter=10;
    case 5
    otherwise
        error('Too many input arguments')
end

if(sz(4)~=length(TE))
    error('#TE and echo image dim didn''t match')
end

if(isempty(Fm_initial))
    if(exist('UMPIRE_unwrapp_3D.m','file') && length(TE)>2)
        Fm_initial=UMPIRE_unwrapp_3D(reshape(im(:,1:3),[sz(1:3) 3]),TE(1:3));
        Fm_initial=Fm_initial(:);
    else
        Fm_initial=angle(im(:,2).*conj(im(:,1)))./(diff(TE(1:2)));
    end
end

%calculate diag{1/(d_j+beta .c)} here we use array multiplication.
d_j=get_dj(im,TE);
factor=1./(d_j+16*beta);


%Calculte the fieldmaps iteratively
switch(nargout)
    case {0,1}
        FieldMap=Fm_initial;
        for i=1:maxIter-1
            grad=calcgradients(FieldMap(:),im,TE,sz(1:3),beta);
            FieldMap=FieldMap-factor.*grad;
        end
        FieldMap=reshape(FieldMap,sz(1:3));
    case 2 % store intermediate field maps: need more memory
        Fm_steps=zeros([length(Fm_initial) maxIter]);
        Fm_steps(:,1)=Fm_initial;
        for i=1:maxIter-1
            grad=calcgradients(Fm_steps(:,i),im,TE,sz(1:3),beta);
            Fm_steps(:,i+1)=Fm_steps(:,i)-factor.*grad;
        end
        Fm_steps=reshape(Fm_steps,[sz(1:3) maxIter]);
        FieldMap=Fm_steps(:,:,end);
    case 3 % store intermediate field maps and gradients: need even more memory
        Fm_steps=zeros([length(Fm_initial) maxIter]);
        g_steps=zeros([length(Fm_initial) maxIter-1]);
        Fm_steps(:,1)=Fm_initial;
        for i=1:maxIter-1
            g_steps(:,i)=calcgradients(Fm_steps(:,i),im,TE,sz(1:3),beta);
            Fm_steps(:,i+1)=Fm_steps(:,i)-factor.*g_steps(:,i);
        end
        Fm_steps=reshape(Fm_steps,[sz(1:3) maxIter]);
        g_steps=reshape(g_steps,[sz(1:3) maxIter-1]);
        FieldMap=Fm_steps(:,:,end);
    otherwise
        error ('only three output arguments supported')
end
%for debugging: plot cost function and gradients of first pixel
% plotCostfunction(,im(1,N),TE)
end

%% supporting fucntions
function w_j=get_w_j(y,m,n)
%refer paper
denom=sum(abs(y).^2,2);
w_j= abs(y(:,m)).*abs(y(:,n))./(denom+1e-3*(mean(denom(:))));

if(sum(isnan(w_j(:)))>0)
    warning('Setting %d values of w_j to zero',sum(isnan(w_j(:))))
    w_j(isnan(w_j))=0;
end
end
function d_j=get_dj(im,TE)
d_j=zeros([ size(im,1) 1]);
for m=1:(length(TE)-1)
    for n=(m+1):length(TE)
        w_j=get_w_j(im,m,n);
        d_j= d_j+abs(im(:,n).*im(:,m)).*w_j.*(TE(n)-TE(m)).^2;
        
    end
end

end

function grad=calcgradients(omega,y,TE,im_size,beta)
%gradient for next iteration
grad1=zeros(numel(omega),1);
for m=1:(length(TE)-1)
    for n=(m+1):length(TE)
        w_j=get_w_j(y,m,n);
        grad1=grad1+ abs(y(:,m).*y(:,n)).*w_j.*(-1*(TE(n)-TE(m))*sin(angle(y(:,n).*conj(y(:,m)))-omega*(TE(n)-TE(m))));
    end
end
grad1(isnan(grad1))=1e-3;
%   figure,imagesc(reshape(grad1,im_size))
grad=grad1+beta*(getRoughness_grad(omega,im_size));

if(sum(isnan(grad(:))) >0)
    warning('NaN in gradient calculation')
end
%  figure,imagesc(reshape(grad,im_size))
end

function R=getRoughness_grad(omega,im_size) %

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
if(im_size(3)>1)
    R3=imfilter(omega,permute(h,[2 3 1]),'replicate','same','corr');
    R=0.5*(R1+R2+R3);
else
    R=0.5*(R1+R2);
end
R=R(:);
end


function plotCostfunction(Fm_initial,im,TE)
%bvec is the aray of B0 values in Hz to be tested
%im_pxl is the array of all echoes of an invidual pixel
%TE: echo times in seconds

bvec=[min(Fm_initial(:)):max(Fm_initial(:))]*2;



c=zeros([ length(bvec) nchoosek(length(TE),2) ]);
g=zeros([ length(bvec) nchoosek(length(TE),2) ]);
for i=1:length(bvec)
    [~,~,c(i,:),g(i,:),all_legend]=costfunc1D(bvec(i),im_pxl,TE);
end
figure,subplot(211),plot(bvec./(2*pi),reshape( c,length(bvec),[])),title('Cost fucntion')
legend(all_legend)
subplot(212),plot(bvec./(2*pi),reshape( g,length(bvec),[])),title('Gradient')
legend(all_legend)
figure,subplot(211),plot(bvec./(2*pi),sum(reshape( c,length(bvec),[]),2)),title('cumulative Cost fucntion')
subplot(212),plot(bvec./(2*pi),sum(reshape( g,length(bvec),[]),2)),title('Cumulative Gradient')

end
function [cost,grad,cost1,grad1,alllegend]=costfunc1D(omega,im_pxl,TE)
%just to plot cost function for a single voxel without roughness.
%to see the convergence for a given TE values
%
%omega: is in rad/s

L=length(TE);
cost1=zeros([1 nchoosek(L,2)]);
idx=1;
alllegend=cell(1,nchoosek(L,2));
for m=1:(length(TE)-1)
    for n=(m+1):length(TE)
        denom=sum(abs(im_pxl).^2,'all');
        w_j= abs(im_pxl(m)).*abs(im_pxl(n))./denom;
        cost1(idx)= abs(im_pxl(m).*im_pxl(n)).*w_j.*(1-cos(angle(im_pxl(n).*conj(im_pxl(m)))-omega*(TE(n)-TE(m))));
        %         figure,imagesc(reshape( abs(y(:,m).*y(:,n)).*w_j.*(1-cos(angle(y(:,n).*conj(y(:,m)))-omega*(TE(n)-TE(m)))),im_size))
        alllegend{idx}=sprintf('Echo(%d ms,%d ms)',TE(m)*1e3,TE(n)*1e3);
        idx=idx+1;
        
    end
end
cost=sum(cost1,'all');

grad1=zeros([1 nchoosek(L,2)]);
idx=1;
for m=1:(length(TE)-1)
    for n=(m+1):length(TE)
        denom=sum(abs(im_pxl).^2,'all');
        w_j= abs(im_pxl(m)).*abs(im_pxl(n))./denom;
        grad1(idx)= -1*(TE(n)-TE(m).*abs(im_pxl(m).*im_pxl(n))).*w_j.*(sin(angle(im_pxl(n).*conj(im_pxl(m)))-omega*(TE(n)-TE(m))));
        %         figure,imagesc(reshape( abs(y(:,m).*y(:,n)).*w_j.*(1-cos(angle(y(:,n).*conj(y(:,m)))-omega*(TE(n)-TE(m)))),im_size))
        idx=idx+1;
    end
end
% grad=grad1;
grad=sum(grad1,'all');
end
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
if(im_size(3)>1)
    R3=imfilter(omega,permute(h,[2 3 1]),'replicate','same','corr');
    R=0.5*(R1.^2+R2.^2+R3.^2);
else
    R=0.5*(R1.^2+R2.^2);
end
R=R(:);
end
