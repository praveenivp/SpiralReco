classdef MTI_3D
    % An operator for correcting B0 inhomogenity using Multi frequency
    % interpolation for spiral imaging
    %
    %
    %
    %Usage:
    % MTIOP=MTI_3D(b0Map,tk,FT,CoilSens)
    % im_B0corr=MTIOP*CoilData;
    % CoilData should be of the form [nCh,nFE,nIntlv,nPar]=size(CoilData);
    % B0Map is in rad/s
    % tk is in seconds of one interleave
    % CoilSens: [CHAxCOLxLINxPAR]
    %
    %Dependencies:
    % NUFFT functions from MIRT toolbox(by fressler) and
    % NUFFT operator from ESPIRIT toolbox (by LUSTIG)
    %
    %
    %Reference:
    %Man, L. C., Pauly, J. M., & Macovski, A. (19m 97).
    %Multifrequency interpolation for fast off-resonance correction.
    %Magnetic Resonance in Medicine, 37(5), 785–792.
    %DOI: https://doi.org/10.1002/mrm.1910370523
    %
    %Author: praveenivp
    properties
        MTI_weights %Coeffients for interpolation
        tflag %transpose flag
        B0map % Nd matrix in HZ
        NUFFTOP %Gridder from FRESSLER/LUSTIG toolbox
        tk % %time points in seconds
        CoilSens
        mode % Interpolation mode {'LeastSquares','NearestNeighbour'}
        nLevels % Number of frequency levels for MFI
        tau
    end
    methods
        function obj=MTI_3D(varargin)
            if(nargin==4)
                obj.B0map=varargin{1};
                obj.tk=varargin{2};
                obj.NUFFTOP=varargin{3};
                obj.CoilSens=varargin{4};
                obj.tflag=0;
                
                %calculate some parameters
                obj.nLevels=round((max(abs(obj.B0map(:)))*max(obj.tk)/pi));
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
                %InData: is kspace data  size:  [size(obj.CoilSens,1) obj.NUFFTOP.dataSize size(obj.CoilSens,4)]
                
                %
                [nCh,nFE,nIntlv,nPar]=size(InData);
                InData =permute(InData,[2 3 4 1]);
                %                 InData=reshape(InData,[nCh,nFE,nIntlv,nPar]);
                sig_MTI=bsxfun(@times,permute(conj(obj.MTI_weights),[2 3 4 5 1]),InData);
                
                all_images=zeros(obj.NUFFTOP.imSize(1),obj.NUFFTOP.imSize(2),nPar,obj.nLevels+1);
                
                for idx_freq=0:obj.nLevels
                    
                    ch_images=  obj.NUFFTOP'*double(sig_MTI(:,:,:,:,idx_freq+1));
                    ch_images=fftshift(ifft(ifftshift(ch_images,3),[],3),3)*sqrt(nPar);
                    
                    ch_images=bsxfun(@times,ch_images,exp(1i*obj.B0map*obj.tau*idx_freq));
                    %coil combination
                    all_images(:,:,:,idx_freq+1)=sum(ch_images.*permute(conj(obj.CoilSens),[2, 3, 4,1]),4);
                end
                out=sum(all_images,4);
            else % forward operation image to kspace
                
                [nCh,~,~,nPar]=size(obj.CoilSens);
                %                  nIntlv=round(length(obj.NUFFTOP.w)/length(obj.tk));
                InData=reshape(InData,[],obj.NUFFTOP.imSize(1),obj.NUFFTOP.imSize(2),nPar);
                InData=bsxfun(@times,InData,(obj.CoilSens));
                InData=permute(InData,[2,3,4,1]); % COLxLINxPARxCHA
                
                %                  nIntlv=prod(obj.NUFFTOP.dataSize)/length(obj.tk);
                ksp_MTI=zeros(obj.NUFFTOP.dataSize(1),obj.NUFFTOP.dataSize(2),size(obj.CoilSens,4),nCh,obj.nLevels+1);
                
                for idx_freq=0:obj.nLevels
                    %tk is in s
                    
                    temp=double(bsxfun(@times,InData,exp(-1i*obj.B0map*obj.tau*idx_freq)));
                    temp=(fftshift(fft(ifftshift(temp,3),[],3),3))/sqrt(nPar);
                    ksp_MTI(:,:,:,:,idx_freq+1) =obj.NUFFTOP*temp;
                    
                end
                
                %                 ksp_MTI=squeeze(ksp_MTI);
                out=sum(bsxfun(@times,permute((obj.MTI_weights),[2 3 4 5 1]), ksp_MTI),5); %Using forward weights has 3 order of magnitude less error
                out=permute(out,[4,1,2,3]);
            end
        end
        function obj=ctranspose(obj)
            if(obj.tflag==0)
                obj.tflag=1;
            end
        end
        
        function obj=CalcWeights(obj)
            if(sum(obj.B0map,'all')==0 ||isempty(obj.B0map))
                obj.MTI_weights=1;
            else
                %                 h=histogram(obj.B0map,"NumBins",prod(0.5*obj.NUFFTOP.imSize));
                %                 bins=h.BinEdges(h.BinCounts>0)+0.5*h.BinWidth;
                [Ncount,edges] = histcounts(obj.B0map(:),prod(0.5*obj.NUFFTOP.imSize));
                bins=edges(Ncount>0)+0.5*diff(edges(1:2));
                wn_tau_l=exp(1i*(bins(:)*(obj.tau.* (0:obj.nLevels))));
                if(strcmp(obj.mode,'LeastSquares'))
                    %                 actual=exp(1i*col(obj.B0map)*obj.tk');
                    actual=exp(1i*col(bins)*obj.tk');
                    obj.MTI_weights=mldivide(wn_tau_l'*wn_tau_l,wn_tau_l'*actual);
                elseif(strcmp(obj.mode,'NearestNeighbour') )
                    error('Not implemented')
                end
            end
        end
    end
    
    
end