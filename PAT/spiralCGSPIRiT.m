function [img_cg,flag,relres,iter,resvec]=spiralCGSPIRiT(SOSOP,GOP3D,Datau,varargin)
%  [img_cgsense,flag,relres,iter,resvec]=spiralCGSENSE(SOSOP,Datau,varargin)
% 
% INPUTS:
% StackofSpiral or StackofSpiralB0 operator
% undersampled data: [ColxLinxParxCha]
% varargin : Name-Value pairs
% 'tol'    : scalar lsqr tolerance
% 'maxit'  : Integer lsqr max interations
% 'reg'    :{'none','Tikhonov'}
% 'lambda' : scalar as regularization parameter

p=inputParser;
p.addParameter('tol',1e-6,@(x)isscalar(x));
p.addParameter('maxit',10,@(x) isscalar(x));
p.addParameter('reg','none',@(x) any(validatestring(x,{'none','Tikhonov'})));
p.addParameter('lambda',1,@(x) isscalar(x));
parse(p,varargin{:});

if(~isa(Datau,SOSOP.precision))
    if(strcmp(SOSOP.precision,'single'))
        Datau=single(Datau);
    else
        Datau=double(Datau);
    end
end

switch(p.Results.reg)
    case 'none'
        E = @(x,tflag) mySPIRiT(x,SOSOP,GOP3D,p.Results.lambda,tflag);
    case 'Tikhonov'
        error('not implemented')
%         E_FT=@(x,transp) myCGSENSE3D_reg(x,SOSOP,p.Results.lambda,transp);
%         reg_out=zeros(SOSOP.imSize,'like',Datau);
%         reg_out=reg_out(:);
end
imSize=size(GOP3D.KERNEL); imSize(end)=[];        
b=zeros(imSize,'like',Datau);

if(strcmp(SOSOP.precision,'double'))
%     csm_sq = sum(SOSOP.op.sens.^2,[1 3]); csm_sq(csm_sq < eps) = 1;
%     M=spdiags(double(csm_sq(:)), 0, numel(csm_sq), numel(csm_sq));
     [img_cg,flag,relres,iter,resvec] = lsqr(E, [Datau(:);b(:)], p.Results.tol,p.Results.maxit);
 else
    [img_cg,flag,relres,iter,resvec] = lsqr(E, [Datau(:);b(:)], p.Results.tol,p.Results.maxit);
end
img_cg=reshape(img_cg,imSize);

end


function outp =  mySPIRiT(inp,SOSOP,GOP,lambda,transpose_indicator)
% inp[Ncha nFE nIntlv nPar]
% w_scale =sum(sqrt(SOSOP.w(:)))/numel(SOSOP.w(:));
ImSize=size(GOP.KERNEL); ImSize(end)=[];  
Ncha = size(GOP.KERNEL,4);
dataSize=[SOSOP.dataSize, Ncha ];
if (strcmp(transpose_indicator,'transp'))
    inp1 = reshape(inp(1:prod(dataSize)),dataSize);
    inp2 = reshape(inp((prod(dataSize)+1):end),ImSize);
    outp1=SOSOP'*inp1;
    outp2=GOP'*inp2;
    outp = (outp1(:) + lambda*outp2(:));
    
elseif (strcmp(transpose_indicator, 'notransp'))
    inp = reshape(inp,ImSize);
    outp1=SOSOP*inp;
    outp2=GOP*inp;
    outp =  [outp1(:);lambda*outp2(:)];
else
    error('Transpose flag not appropriately defined');
end

end