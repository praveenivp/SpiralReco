function [img_cgsense,flag,relres,iter,resvec]=spiralCGSENSE(SOSOP,Datau)

E_FT=@(x,transp) myCGSENSE3D(x,SOSOP,transp);
reg_out=[];
limit=1e-9;

% csm_sq = squeeze(sum(csm .* conj(csm),1)); csm_sq(csm_sq < eps) = 1;
% M = spdiag(sqrt(csm_sq)); %Preconditioner
nIterCG=5; 
[img_cgsense,flag,relres,iter,resvec] = lsqr(E_FT, ([Datau(:); reg_out]), limit,nIterCG);
img_cgsense=reshape(img_cgsense,SOSOP.imSize(1),SOSOP.imSize(2),[]);

end

function outp =  myCGSENSE3D(inp,NUFFT_obj,transpose_indicator)
% inp[nFE x nIntlv x nCha x npar ]
% scale = sqrt(prod(prod(nufft_st.Kd))/numel(weights(:)));
% scale=1; %acceleration factor
% Npar=NUFFT_obj.dataSize(3);
Ncha=NUFFT_obj.op.sensChn;
if (strcmp(transpose_indicator,'transp'))
      inp=reshape(inp,[NUFFT_obj.dataSize,Ncha ]);
      outp = NUFFT_obj'*inp;
      outp= outp(:);

elseif (strcmp(transpose_indicator, 'notransp'))
    inp=reshape(inp,NUFFT_obj.imSize);
    outp = NUFFT_obj*inp;
    outp=outp(:);
else
    error('Transpose flag not appropriately defined');
end

end