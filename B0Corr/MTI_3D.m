classdef MTI_3D
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
                obj.nLevels=ceil(2*(max(abs(obj.B0map(:)))*max(obj.tk)/pi));
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
                %                 N=obj.NUFFTOP.imSize(1);
                [nCh,nFE,nPar,nRep]=size(InData);
                nIntlv=nFE/length(obj.tk);
                
                all_images=zeros(obj.NUFFTOP.imSize(1),obj.NUFFTOP.imSize(2),nPar,nRep,obj.nLevels+1);
                ch_images=zeros(nCh,obj.NUFFTOP.imSize(1),obj.NUFFTOP.imSize(2),nPar);
                sig_MTI=bsxfun(@times,obj.MTI_weights',reshape(permute(InData,[2 1 3 4 5]),nFE/nIntlv,1,[]));
                sig_MTI=reshape(sig_MTI,nFE/nIntlv,obj.nLevels+1,nIntlv,nCh,nPar,nRep);
                sig_MTI=reshape(permute(sig_MTI,[1 3 2 4 5 6 7]),nFE,obj.nLevels+1,nCh,nPar,nRep);
                
                sig_MTI=bsxfun(@times,sig_MTI,(obj.NUFFTOP.w));
                
                for idx_freq=0:obj.nLevels
                    for rep=1:nRep
                        for ii=1:nCh
                        for par=1:nPar
                                temp=  obj.NUFFTOP.st.p'*double(sig_MTI(:,idx_freq+1,ii,par,rep));
                                temp=reshape(temp,obj.NUFFTOP.imSize(1)*2,obj.NUFFTOP.imSize(2)*2);
                                
                                temp=ifft2(temp);
                                temp=temp(1:obj.NUFFTOP.imSize(1),1:obj.NUFFTOP.imSize(2));
                                temp=temp.*(obj.NUFFTOP.st.sn');
                                %                         temp=temp.*(1./obj.FT.w);
                                
                                ch_images(ii,:,:,par) =temp;
                                %                             all_images(ii,:,:,slc,rep,idx_freq+1) =(obj.NUFFTOP'*sig_MTI(nFE,idx_freq+1,ii,slc,rep)).*exp(-1i*obj.B0map*obj.tau*idx_freq);
                        end 
                            
                        end
                        ch_images=fftshift(fft(ch_images,[],4),4);
                        ch_images=bsxfun(@times,ch_images,permute(exp(-1i*obj.B0map*obj.tau*idx_freq),[4 1 2 3]));
                        all_images(:,:,:,rep,idx_freq+1)=squeeze(sum(ch_images.*obj.CoilSens,1));
                    end
                end
                out=sum(all_images,5);
            else % forward operation image to kspace
                % just for testing for single channel sim data not properly implemented
                
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
                out=(sum(obj.MTI_weights'.* ksp_MTI,2)); %Using forward weights has 3 order of magnitude less error
                
                
            end
        end
        function obj=ctranspose(obj)
            if(obj.tflag==0)
                obj.tflag=1;
            end
        end
        
        function obj=CalcWeights(obj)
               h=histogram(obj.B0map,"NumBins",prod(10*obj.NUFFTOP.imSize));
                bins=h.BinEdges(h.BinCounts>0)+0.5*h.BinWidth;
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