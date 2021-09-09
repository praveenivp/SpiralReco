function [img_cgsense,flag,relres,iter,resvec]=spiralCGSENSE(SOSOP,Datau,varargin)
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
if(SOSOP.op.sensChn==0)
    error('Input StackofSpiral Operator doesnt have coil senstivies');
end

p=inputParser;
p.addParameter('tol',1e-6,@(x)isscalar(x));
p.addParameter('maxit',10,@(x) isscalar(x));
p.addParameter('reg','none',@(x) any(validatestring(x,{'none','Tikhonov'})));
p.addParameter('lambda',1e-3,@(x) isscalar(x));
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
        E_FT=@(x,transp) myCGSENSE3D(x,SOSOP,transp);
        reg_out=[];
    case 'Tikhonov'
        E_FT=@(x,transp) myCGSENSE3D_reg(x,SOSOP,p.Results.lambda,transp);
        reg_out=zeros(SOSOP.imSize(SOSOP.imSize>0),'like',Datau);
        reg_out=reg_out(:);
end
        

if(strcmp(SOSOP.precision,'double'))
    if(ndims(SOSOP.op.sens)==3)
        csm_sq = sum(SOSOP.op.sens.^2,[1,3]);
    else
        csm_sq = sum(SOSOP.op.sens.^2,4);
    end
    csm_sq(csm_sq < eps) = 1;
    M=spdiags(double(csm_sq(:)), 0, numel(csm_sq), numel(csm_sq));
    [img_cgsense,flag,relres,iter,resvec] = lsqr(E_FT, ([Datau(:); reg_out]),p.Results.tol,p.Results.maxit,M);
    else
    [img_cgsense,flag,relres,iter,resvec] = lsqr(E_FT, ([Datau(:); reg_out]),p.Results.tol,p.Results.maxit);
end
img_cgsense=reshape(img_cgsense,SOSOP.imSize(1),SOSOP.imSize(2),[]);

end

function outp =  myCGSENSE3D(inp,NUFFT_obj,transpose_indicator)
% inp[nFE x nIntlv x nCha x npar ]
% scale = sqrt(prod(prod(nufft_st.Kd))/numel(weights(:)));
Ncha=NUFFT_obj.op.sensChn;
if (strcmp(transpose_indicator,'transp'))
      inp=reshape(inp,[NUFFT_obj.dataSize,Ncha ]);
      outp = NUFFT_obj'*inp;
      outp= outp(:);

elseif (strcmp(transpose_indicator, 'notransp'))
    inp=reshape(inp,NUFFT_obj.imSize(1),NUFFT_obj.imSize(2),[]);
    outp = NUFFT_obj*inp;
    outp=outp(:);
else
    error('Transpose flag not appropriately defined');
end

end

function outp =  myCGSENSE3D_reg(inp,NUFFT_obj,lambda,transpose_indicator)

% Coilsens [cha x colx Lin x Par]
% inp[nFE x nIntlv x nCha x npar ]
Ncha=NUFFT_obj.op.sensChn;
if (strcmp(transpose_indicator,'transp'))
       inpSize=[NUFFT_obj.dataSize,Ncha];
      outp = NUFFT_obj'*reshape(inp(1:prod(inpSize)),inpSize);
      if(numel(inp)>prod(inpSize)) %reg
         outp=outp(:)+lambda*double(inp((prod(inpSize)+1):end));
      end
elseif (strcmp(transpose_indicator, 'notransp'))
    inp=reshape(inp,NUFFT_obj.imSize(1),NUFFT_obj.imSize(2),[]);
    reg=double(inp(:).*conj(inp(:)));
    outp = NUFFT_obj*double(inp);
    outp=[outp(:) ;reg];
else
    error('Transpose flag not appropriately defined');
end

end