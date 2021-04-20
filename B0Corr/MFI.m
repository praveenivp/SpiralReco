classdef MFI
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
    %Man, L. C., Pauly, J. M., & Macovski, A. (1997).
    %Multifrequency interpolation for fast off-resonance correction. 
    %Magnetic Resonance in Medicine, 37(5), 785–792. 
    %DOI: https://doi.org/10.1002/mrm.1910370523
   
    properties
        MFI_weights %Coeffients for interpolation
        tflag %transpose flag
        B0map %in HZ
        NUFFTOP %FRESSLER/LUSTIG toolboxes
        tk % %time points in seconds
        CoilSens
        mode % Interpolation mode {'LeastSquares','NearestNeighbour'}
        nLevels % Number of frequency levels for MFI
        wi % frequencies in rad/s for MFI
    end
    methods
        function obj=MFI(varargin)
            if(nargin==4)
                obj.B0map=varargin{1};
                obj.tk=varargin{2};
                obj.NUFFTOP=varargin{3};
                obj.CoilSens=varargin{4};
                obj.tflag=0;
                
                %calculate some parameters
                obj.nLevels=ceil((1*max(abs(obj.B0map(:)))*max(obj.tk)/pi));
                obj.wi=linspace(-1*max(abs(obj.B0map(:))),1*max(abs(obj.B0map(:))),obj.nLevels);
                
                %if field map is zero, it becomes normal NUFFT operator
                if(obj.nLevels==0)
                        obj.nLevels=1;
                        obj.wi=0;
                end
                obj.mode='LeastSquares'; % {'LeastSquares','NearestNeighbour'}
                
                % calculate weights
                obj=obj.CalcWeights();
                
            else
                error('Need several input parameters')
            end
        end
        
        function out=mtimes(obj,InData)
            if(obj.tflag==1) % Reverse operator: kspace data to image
                %hard coded para
                N=obj.NUFFTOP.imSize(1);
                [nCh,nFE]=size(InData);
                nIntlv=nFE/length(obj.tk);
                img_MFI=zeros(nCh,N,N,length(obj.wi));
                b0term=exp(-1i.*col(repmat(obj.tk(:),[1 nIntlv]))*obj.wi(:).');
                InData=double(bsxfun(@times,InData,permute(b0term,[3, 1,2])));
                for idx_freq=1:length(obj.wi)
                    %tk is in s
%                     b0term=exp(-1i.*(repmat(obj.tk.*obj.wi(idx_freq),[1 nIntlv]))); % add the frequency levels
                    
                    for ii=1:nCh
                        	b = col(InData(ii,:,idx_freq)).*(obj.NUFFTOP.w(:));
                            res = nufft_adj(b, obj.NUFFTOP.st)/sqrt(prod(obj.NUFFTOP.imSize));
                            img_MFI(ii,:,:,idx_freq)= reshape(res, obj.NUFFTOP.imSize);
%                         img_MFI(ii,:,:,idx_freq) =obj.NUFFTOP'*col(InData(ii,:,idx_freq));
                    end
                end

                out=sum(bsxfun(@times,permute(flip(obj.MFI_weights,3),[4, 1, 2, 3]),img_MFI),4);
                
               if(~(isscalar(obj.CoilSens)||isempty(obj.CoilSens)))
                     %try to combine coils
                    out=squeeze(sum(out.*obj.CoilSens,1));
                end
                
%                 img_MFI_combined=zeros(N,N,obj.nLevels);
%                 for idx_freq=1:length(obj.wi)
%                      img_MFI_combined(:,:,idx_freq)=squeeze(sum((obj.CoilSens).* permute(img_MFI(:,:,:,slc,rep,idx_freq),[ 1 2 3 ]),1));
%                 end
%                 img_MFI_combined=squeeze(img_MFI_combined);
%                 
%                 out=sum(flip(obj.MFI_weights,3).*img_MFI_combined ,3);
                
            else % forward operation image to kspace
                
               if(~(isscalar(obj.CoilSens)||isempty(obj.CoilSens)))
                     %try to combine coils
                    out=bsxfun(@times, permute(InData,[3, 1, 2]),conj(obj.CoilSens));
               end
                
                out=bsxfun(@times,permute(flip(conj(obj.MFI_weights),3),[4, 1, 2, 3]),out);
                 nIntlv=round(length(obj.NUFFTOP.w)/length(obj.tk));
                kData=zeros(length(obj.tk),length(obj.wi));
                
                
               out1=zeros([prod(obj.NUFFTOP.dataSize) length(obj.wi) size(out,1)]);
                for idx_freq=1:length(obj.wi)
                for ch=1:size(out,1)
                    
                    	b = reshape(out(ch,:,:,idx_freq),obj.NUFFTOP.imSize(1),obj.NUFFTOP.imSize(2));
                        res = nufft(b, obj.NUFFTOP.st)/sqrt(prod(obj.NUFFTOP.imSize)).*(obj.NUFFTOP.w(:));
                    	out1(:,idx_freq,ch) = res(:);
                   
                end
                end
                b0term=exp(1i.*col(repmat(obj.tk(:),[1 nIntlv]))*obj.wi(:).');
                
                out=squeeze(sum(bsxfun(@times,b0term,out1),2)).';
%                 out=  mean(kData*reshape(flip(obj.MFI_weights,4),prod(obj.NUFFTOP.imSize) ,[])',2); %kspace
            end
        end
        function obj=ctranspose(obj)
            if(obj.tflag==0)
                obj.tflag=1;
            end
        end
        
        function obj=CalcWeights(obj)

            witk=exp(1i.*obj.tk*obj.wi);
            obj.MFI_weights=zeros([size(obj.B0map) obj.nLevels]);
            
            % calculate weights for all values in b0maps
            if(strcmp(obj.mode,'LeastSquares'))
%                 warning('turn off:MATLAB:nearlySingularMatrix');
%                 warning('off','MATLAB:nearlySingularMatrix');
%                 for i=1:size(obj.MFI_weights,1)
%                     for j=1:size(obj.MFI_weights,2)
%                         B=exp(-1i*obj.B0map(i,j)*obj.tk);
%                         obj.MFI_weights(i,j,:)= mldivide((witk'*witk),(witk'*B));
%                     end
%                 end
%                  warning('on','MATLAB:nearlySingularMatrix');
                %memory hog way but faster 
                 B=exp(-1i*obj.B0map(:)*obj.tk.').';
                obj.MFI_weights=mldivide((witk'*witk),(witk'*B));
                obj.MFI_weights=permute(reshape(obj.MFI_weights,[],size(obj.B0map,1),size(obj.B0map,2)),[2 3 1]);
               
            elseif(strcmp(obj.mode,'NearestNeighbour') )
                for i=1:size(obj.MFI_weights,1)
                    for j=1:size(obj.MFI_weights,2)
                        [~,I]=min(abs(obj.wi-obj.B0map(i,j)));
                        obj.MFI_weights(i,j,I)=1;
                    end
                end
            end
            
        end
    end
    
    
end