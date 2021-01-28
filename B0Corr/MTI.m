classdef MTI
    % An operator for correcting B0 inhomogenity using Multi frequency
    % interpolation for spiral imaging
    %
    %
    %
    %Usage:
    % MFIOP=MFI(b0Map,tk,FT,CoilSens)
    % im_B0corr=MFI*CoilData;
    % CoilData should be of the form [nCh,nFE,nSlc,nRep]=size(CoilData);
    %
    %
    %Dependencies:
    % NUFFT operator from MIRT toolbox(by fressler) and
    % NUFFT operator from ESPIRIT toolbox (by LUSTIG)
    %
    %
    %Reference:
    %Man, L. C., Pauly, J. M., & Macovski, A. (19m 97).
    %Multifrequency interpolation for fast off-resonance correction.
    %Magnetic Resonance in Medicine, 37(5), 785–792.
    %DOI: https://doi.org/10.1002/mrm.1910370523
    
    properties
        MTI_weights %Coeffients for interpolation
        MTI_weights_Forward %Coeffients for interpolation
        tflag %transpose flag
        B0map %in HZ
        NUFFTOP %FRESSLER/LUSTIG toolboxes
        tk % %time points in seconds
        CoilSens
        mode % Interpolation mode {'LeastSquares','NearestNeighbour'}
        nLevels % Number of frequency levels for MFI
        tau
    end
    methods
        function obj=MTI(varargin)
            if(nargin==4)
                obj.B0map=varargin{1};
                obj.tk=varargin{2};
                obj.NUFFTOP=varargin{3};
                obj.CoilSens=varargin{4};
                obj.tflag=0;
                
                %calculate some parameters
                obj.nLevels=ceil((2*max(abs(obj.B0map(:)))*max(obj.tk)/pi));
                obj.tau=max(obj.tk)/(obj.nLevels+1);
                obj.mode='LeastSquares'; % {'LeastSquares','NearestNeighbour'}
                
                % calculate weights
                obj=obj.CalcWeights();
                
            else
                error('Need several input parameters')
            end
        end
        
        function out=mtimes(obj,InData)
            if(obj.tflag==1) % Reverse operator: kspace data to image
                %out: is corrected image
                %InData: is kspace data
                
                %hard coded para
                N=obj.NUFFTOP.imSize(1);
                [nFE,nCh,nSlc,nRep]=size(InData);
                nIntlv=nFE/length(obj.tk);
                ksp_MTI=zeros(nCh,nFE,nSlc,nRep,obj.nLevels+1);
                img_MFI_tau0 =obj.NUFFTOP'*(col(InData(:,1,1,1)));
                for idx_freq=0:obj.nLevels
                    %tk is in s
                    
                    for rep=1:nRep
                        for slc=1:nSlc
                            for ii=1:nCh
                                ksp_MTI(ii,:,slc,rep,idx_freq+1) =obj.NUFFTOP*(((img_MFI_tau0)).*exp(1i*obj.B0map*obj.tau*idx_freq));
                            end
                        end
                    end
                end
                ksp_MTI=squeeze(ksp_MTI);
                % img_MFI= (flip(flip(img_MFI,2),3));
                % if(~exist('CoilSens','var'))
                %     img_sos=sum(conj(img_MFI(:,:,:,:,1,ceil(length(wi)/2))).*img_MFI(:,:,:,:,1,ceil(length(wi)/2)),1);
                % CoilSens=img_MFI(:,:,:,:,1,ceil(length(wi)/2))./
                % end
                
                %                 img_MFI_combined=zeros(N,N,nSlc,nRep,obj.nLevels);
                %                 for idx_freq=1:length(obj.wi)
                %                     for rep=1:nRep
                %                         for slc=1:nSlc
                %                             img_MFI_combined(:,:,slc,rep,idx_freq)=squeeze(sum((obj.CoilSens).* permute(img_MFI(:,:,:,slc,rep,idx_freq),[ 1 2 3 ]),1));
                %                             %      img_MFI_combined(:,:,slc,rep,idx_freq)=squeeze(adaptiveCombine(permute(img_MFI(:,:,:,slc,rep,idx_freq),[3 1 2])...
                %                             %          ,[9 9 1],true,true,[1 1 1],1,true));
                %                             %     img_MFI_combined(:,:,slc,rep,idx_freq)=squeeze(adaptiveCombine2(permute(img_MFI(:,:,:,slc,rep,idx_freq),[3 1 2])));
                %                         end
                %                     end
                %                 end
                %                 img_MFI_combined=squeeze(img_MFI_combined);
                
                k_interp=sum(obj.MTI_weights'.* ksp_MTI ,2);
                 out=conj(obj.NUFFTOP'*k_interp(:));
                
%                 out=img_MFI_tau0;
                
                %just trying
                
                
                
            else % forward operation image to kspace
                
                
                
                
                %                  nIntlv=round(length(obj.NUFFTOP.w)/length(obj.tk));
                nCh=1;
                ksp_MTI=zeros(nCh,length(obj.tk),obj.nLevels+1);
                
                for idx_freq=0:obj.nLevels
                    %tk is in s
                    
                    for ii=1:nCh
                        ksp_MTI(ii,:,idx_freq+1) =obj.NUFFTOP*(reshape(InData,obj.NUFFTOP.imSize).*exp(-1i*obj.B0map*obj.tau*idx_freq));
                    end
                end
                
                ksp_MTI=squeeze(ksp_MTI);
                
                %using tranpose weigths doesn't seem to work, may be i'm stupid
%                 out=conj(sum(flip(obj.MTI_weights.',2).* ksp_MTI ,2)); %
%                 
                out=(sum(obj.MTI_weights_Forward'.* ksp_MTI,2)); %Using forward weights has 3 order of magnitude less error
                
               
            end
        end
        function obj=ctranspose(obj)
            if(obj.tflag==0)
                obj.tflag=1;
            end
        end
        
        function obj=CalcWeights(obj)
%             [W,ai,bi]=unique(obj.B0map(:));
%             wn_tau_l=exp(-1i*obj.tau.*(W(:)*(0:obj.nLevels)));
%             actual=exp(1i*W(:)*obj.tk');
            wn_tau_l=exp(-1i*(obj.B0map(:)*(obj.tau.* (0:obj.nLevels))));
%             obj.MTI_weights=zeros([numel(obj.B0map) obj.nLevels]);
            
            % calculate weights for all values in b0maps
            if(strcmp(obj.mode,'LeastSquares'))
                
                %                 for i=1:size(obj.MTI_weights,1)
                %                     for j=1:size(obj.MTI_weights,2)
                %                         B=exp(-1i*obj.B0map(i,j)*obj.tk);
                actual=exp(1i*obj.B0map(:)*obj.tk');
                %                         obj.MTI_weights= mldivide(wn_tau_l,actual);
                %much mfater
                obj.MTI_weights= mldivide(wn_tau_l'*wn_tau_l,wn_tau_l'*actual);
                 actual=exp(-1i*col(obj.B0map)*obj.tk');
                 obj.MTI_weights_Forward=mldivide(wn_tau_l'*wn_tau_l,wn_tau_l'*actual);
                %                     end
                %                 end
                
                warning('on','MATLAB:nearlySingularMatrix');
            elseif(strcmp(obj.mode,'NearestNeighbour') )
                error('Not implemented')
                %                 for i=1:size(obj.MFI_weights,1)
                %                     for j=1:size(obj.MFI_weights,2)
                %                         [~,I]=min(abs(obj.wi-obj.B0map(i,j)));
                %                         obj.MFI_weights(i,j,I)=1;
                %                     end
                %                 end
            end
            
        end
    end
    
    
end