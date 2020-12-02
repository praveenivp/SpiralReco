function [im_b0corr]=B0correctNUFFT(CoilData,b0Map,tk,FT,CoilSens)
%[im_b0corr]=B0coreect(CoilData,b0Map,k,senseMap,parameter)
%function to correct the B0 inhomgenity
% Inputs
% CoilData- kspace data (Readout*Interleaves) X coils
% k - kspace Trajectory scaled [-0.5 0.5 ]
% b0Map- Field map in rad/s
% senseMap- Coils sensitivity map 4D matrix [image_dim1 image_dim2 image_dim3 coils ]
% paramter: struct
%        parameter.time= time in s for 1  interleave
%        parameter.psotion= off center postion of slice [x y].*(2*kmax) ||||| factor is beacuse of scaled k

%steps for MFI Least square approach.
% get max B0 offresonance frequency
% calculate number of quantisations required(L).
tk=tk(:)*1e-6;
  Levels=ceil((0.8*max(abs(b0Map(:)))*max(tk)/pi));
%  Levels=100;
% Levels=ceil(2.2*(2*pi*600*max(tk))/pi);
wi=linspace(-0.8*max(abs(b0Map(:))),0.8*max(abs(b0Map(:))),Levels);
% wi=linspace(-600*2*pi,600*2*pi,Levels);

% [MFI_weights]=CalcWeights('LeastSquares');

%hard coded para
N=FT.imSize(1);
[nCh,nFE,nSlc,nRep]=size(CoilData);
nIntlv=nFE/length(tk);
img_MFI=zeros(nCh,N,N,nSlc,nRep,length(wi));
for idx_freq=1:length(wi)
    %tk is in s
     b0term=exp(-1i.*(repmat(tk.*wi(idx_freq),[1 nIntlv]))); % add the frequency levels
for rep=1:nRep
    for slc=1:nSlc
        for ii=1:nCh
            img_MFI(ii,:,:,slc,rep,idx_freq) = FT'*(col(CoilData(ii,:,slc,rep)).*b0term(:));
        end
    end
end    
end
% img_MFI= (flip(flip(img_MFI,2),3));
% if(~exist('CoilSens','var'))
%     img_sos=sum(conj(img_MFI(:,:,:,:,1,ceil(length(wi)/2))).*img_MFI(:,:,:,:,1,ceil(length(wi)/2)),1);
% CoilSens=img_MFI(:,:,:,:,1,ceil(length(wi)/2))./
% end

img_MFI_combined=zeros(N,N,nSlc,nRep,Levels);
for idx_freq=1:length(wi)
for rep=1:nRep
    for slc=1:nSlc
    img_MFI_combined(:,:,slc,rep,idx_freq)=squeeze(sum((CoilSens).* permute(img_MFI(:,:,:,slc,rep,idx_freq),[ 1 2 3 ]),1));
%      img_MFI_combined(:,:,slc,rep,idx_freq)=squeeze(adaptiveCombine(permute(img_MFI(:,:,:,slc,rep,idx_freq),[3 1 2])...
%          ,[9 9 1],true,true,[1 1 1],1,true));
%     img_MFI_combined(:,:,slc,rep,idx_freq)=squeeze(adaptiveCombine2(permute(img_MFI(:,:,:,slc,rep,idx_freq),[3 1 2])));
    end
end    
end
 img_MFI_combined=squeeze(img_MFI_combined);
%  img_MFI_combined= (flip(flip(img_MFI_combined,1),2));
%  img_MFI_combined=squeeze(sqrt(sum(abs(img_MFI).^2,3)));

% for every value of omega in the field map get the weights to combine the
% multi frequency images.(pixel*pixel*weights 3d matrix)

[MFI_weights]=CalcWeights('LeastSquares');
% MFI_weights_shift1=flip(MFI_weights,3);
im_b0corr=sum(flip(MFI_weights,3).*img_MFI_combined ,3);
%multiply and add everything
toc

    %nested funciton: can share varibles of the parent fucntion
    function [MFI_weights]=CalcWeights(Mode)
        if(nargin<1)
            Mode='LeastSquares';
        end
        Levels=length(wi);
        %         wi=linspace(-1*max(abs(b0Map(:))),max(abs(b0Map(:))),Levels);
        %         tk=time(1:33068/4);
        witk=exp(1i.*tk*wi);
        MFI_weights=zeros([size(b0Map) Levels]);
        % calculate weights for all values in b0maps
        if(strcmp(Mode,'LeastSquares'))
             warning('turn off:MATLAB:nearlySingularMatrix');
             warning('off','MATLAB:nearlySingularMatrix');
            for i=1:size(MFI_weights,1)
                for j=1:size(MFI_weights,2)
                    B=exp(-1i*b0Map(i,j)*tk);
%                     if(abs(det(witk'*witk))>1e-16)
                    MFI_weights(i,j,:)= mldivide((witk'*witk),(witk'*B));
%                     else
%                          MFI_weights(i,j,:)=0;% mldivide((witk'*witk),witk'*B);
%                     end
                    %         MFI_weights(i,j,:)= mldivide(witk,B);
                end
            end
            
             warning('on','MATLAB:nearlySingularMatrix');
        elseif(strcmp(Mode,'NearestNeighbour') )
            for i=1:size(MFI_weights,1)
                for j=1:size(MFI_weights,2)
                    [~,I]=min(abs(wi-b0Map(i,j)));
                    MFI_weights(i,j,I)=1;
                end
            end
        end
        
    end
end