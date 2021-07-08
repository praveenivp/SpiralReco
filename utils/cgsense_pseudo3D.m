function [im_cg,im_zp,rtime]=cgsense_pseudo3D(inpdata,FTOP,Coilmaps,B0map_rad_s,adcTime,CGpara,SetSel)
% This function does fft along partition and split the 3D SENSE problem
% into multiple 2D problem which seem working. only works for stack of
% spiral trajectories with no CAIPI shift
% for B0 correction it uses MFI
%
%INPUTS:
%inp: Coildata or twix (FOV shifted and sqrt(DCF compensated) [CHAxCOLxLINxPAR]
%      twix from mapVBVD
%FTOP: NUFFT operator
% Coilmaps: from ESPIRIT [CHAxCOLxLINxPAR]
%B0map_rad_s: B0map in rad/s [COLxLINxPAR]
%CGPARA=struct('limit',1e-6,'nIterCG',10)
%do fft along partition to split into 2D problem
%
%praveen.ivp@gmail.com
if(isa(inpdata,'struct'))
    SpiralPara=getSpiralPara(inpdata);
    soda_obj=getSoda(inpdata,SpiralPara);

    nPar=inpdata.image.NPar;
else
    sig=permute(inpdata,[2 3 1 4 5]); %[CHAxCOLxLINxPAR]-> [COLxLINxCHAxPAR]
    nPar=size(sig,4);
end


im_cg=zeros([FTOP.imSize nPar numel(SetSel)] ,'single');
im_zp=zeros([FTOP.imSize nPar numel(SetSel)] ,'single');
reg_out=[];
rtime=zeros(1,3);
MFIOP2D=cell(1,nPar);
for cPar=1:nPar
    csm2d=Coilmaps(:,:,:,cPar);
    
    tic
    if(B0map_rad_s(1)==0 || isempty(B0map_rad_s))
        MFIOP2D{cPar}=MFI(0,adcTime,FTOP,double(csm2d));
    else
        MFIOP2D{cPar}=MFI(B0map_rad_s(:,:,cPar),adcTime,FTOP,double(csm2d));
    end
    rtime(3)= rtime(3)+toc;
end
fprintf('MFI weights calculated\n')
print_str ='';
for cSet=1:numel(SetSel)
    [sig]=getSig(inpdata,SetSel(cSet),SpiralPara,soda_obj,FTOP);
    sig=permute(sig,[1 3 2 4 5]); %[COLxLINxCHAxPAR]
    sig=(fftshift(ifft(ifftshift(sig,4),[],4),4))/sqrt(nPar);
    for cPar=1:nPar
        
        sig2D=sig(:,:,:,cPar);
        inp=reshape(permute(sig2D,[3,1,2]),size(Coilmaps,1),[]);
        %zeropadding recon
        tic
        im_zp(:,:,cPar,cSet)=single(MFIOP2D{cPar}'*inp);
        rtime(2)= rtime(2)+toc;
        
        %CGSENSE recon
        tic
        E_CGSENSE=@(x,transp) myMFI_SENSEop(x,MFIOP2D{cPar},transp);
        [img_cgsense,flag,relres,iter,resvec] = lsqr(E_CGSENSE, [double(inp(:)); reg_out], CGpara.limit,CGpara.nIterCG);
        im_cg(:,:,cPar,cSet) = single(reshape(img_cgsense,size(Coilmaps,2),size(Coilmaps,3)));
        rtime(1)= rtime(1)+toc;
        
        %     disp(resvec)
        fprintf(repmat('\b',1,numel(print_str)));
        print_str = sprintf(['Par = %2.0f/%2.0f'...
            ' | Set = %2.0f/%2.0f'...
            ' | time = %6.1f min'], cPar,nPar,cSet,numel(SetSel),sum(rtime(1:2))/60);
        fprintf(print_str);
    end
    
    
    
end
rtime(1:2)= rtime(1:2)./numel(SetSel);
fprintf('\nMFI calc time for %d Par is %3.4f seonds\n',nPar,rtime(3))
fprintf('zeropad NUFFT  took  %3.4f seonds\n',rtime(2))
fprintf('CGSENSE took %3.4f seonds for %d Iteration\n',rtime(1),CGpara.nIterCG)

end

function outp =  myMFI_SENSEop(inp,MFIOP,transpose_indicator)

% scale = sqrt(prod(prod(nufft_st.Kd))/numel(weights(:)));
if (strcmp(transpose_indicator,'transp'))
    inp=reshape(inp,31,[]);
    outp=col(MFIOP'*inp);
elseif (strcmp(transpose_indicator, 'notransp'))
    inp=reshape(inp,MFIOP.NUFFTOP.imSize(1),MFIOP.NUFFTOP.imSize(2));
    outp=MFIOP*double(inp);
    outp=outp(:);
else
    error('Transpose flag not appropriately defined');
end

end


function [sig]=getSig(twix,cSet,SpiralPara,soda_obj,FTOP)
cintlv=1:SpiralPara.R_PE:SpiralPara.Ninterleaves;
sig=squeeze(twix.image(:,:,cintlv,:,1,1,1,1,:,cSet,:,:,:,:) );%first rep
% sig=sig(:,:,cintlv,:,:);
KTraj=reshape(complex(FTOP.st.om(:,1),FTOP.st.om(:,2)),size(FTOP.w)) .*(1000/SpiralPara.Resolution);
sig=performFOVShift(sig,KTraj,SpiralPara,soda_obj,FTOP.w);
end

function sig_FOV=performFOVShift(sig,KTraj,SpiralPara,soda_obj,DCF)
if(any(SpiralPara.slice{1}.Position ~=0))
    %             [~,idx]=sort(obj.SpiralPara.slice{1}.Normal);
    %             posi=1e-3*obj.SpiralPara.slice{1}.Position(idx); %m
    pos_PRS=GradientXYZ2PRS(1e-3*[1 -1 -1].*SpiralPara.slice{1}.Position,soda_obj,1); %only work for head first-supine
    
    %                 if(obj.flags.doB0Driftcorr)
    %                     [kHO]=Grad2TrajHigherorder(obj.Grad,obj.SpiralPara);
    %                     B0_mod=exp(-1i*(real(obj.KTraj).*pos_PRS(1)+imag(obj.KTraj).*pos_PRS(2)+squeeze(kHO(:,3,:)).*pos_PRS(3) ));
    %                 else
    B0_mod=exp(-1i*(real(KTraj).*pos_PRS(1)+imag(KTraj).*pos_PRS(2)));
    %                 end
    
    B0_mod=B0_mod.*reshape(DCF,size(B0_mod));
    
    
    
    
    % negative displacement added to real part of kTRaj moves up down in array show
    sig_FOV=bsxfun(@times, sig,reshape(B0_mod,size(sig,1),1,size(sig,3)));
end
end


function soda_obj=getSoda(twix,spiralPara)
soda_obj=   SODA_OBJ( 'mrprot',twix.hdr);
%             ReadoutFOV

soda_obj.NPixelReadout=spiralPara.FOV(1)/spiralPara.Resolution;
soda_obj.NPixelPhase=spiralPara.FOV(1)/spiralPara.Resolution;
soda_obj.PixelSizePhase= soda_obj.ReadoutFOV/soda_obj.NPixelPhase;
soda_obj.PixelSizeReadout= soda_obj.PhaseFOV/soda_obj.NPixelReadout;
soda_obj=soda_obj.calcPixLoc();

%             PhaseFOV

end