classdef MFI3D
    % An operator for correcting B0 inhomogenity using Multi frequency
    % interpolation for spiral imaging
    %
    %
    %
    %Usage:
    % MFIOP=MFI3D(b0Map,tk,FT,CoilSens)
    % im_B0corr=MFI*CoilData;
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
    %Man, L. C., Pauly, J. M., & Macovski, A. (1997).
    %Multifrequency interpolation for fast off-resonance correction.
    %Magnetic Resonance in Medicine, 37(5), 785–792.
    %DOI: https://doi.org/10.1002/mrm.1910370523
    %
    %Author: praveenivp
    properties
        MFI_weights %Coeffients for interpolation
        tflag %transpose flag
        B0map %in rad/s
        B0map_thres %in rad/s
        NUFFTOP %FRESSLER/LUSTIG toolboxes
        tk % %time points in seconds
        CoilSens
        mode % Interpolation mode {'LeastSquares','NearestNeighbour'}
        nLevels % Number of frequency levels for MFI
        wi % frequencies in rad/s for MFI
        max_size
        precision
    end
    methods
        function obj=MFI3D(varargin)
            if(nargin==4)
                obj.B0map=varargin{1};
                obj.tk=varargin{2};
                obj.NUFFTOP=varargin{3};
                obj.CoilSens=varargin{4};
                obj.tflag=0;
                
                % hard coded settings
                obj.max_size=4; %GB
                obj.precision='single';
                obj.mode='LeastSquares'; % {'LeastSquares','NearestNeighbour'}
                
                if(sum(obj.B0map(:))==0||isempty(obj.B0map))
                    %if field map is zero, it becomes normal NUFFT operator
                
                    obj.nLevels=1;
                    obj.wi=0;
                    obj.MFI_weights=1;
                else
                %calculate some parameters
                if(strcmpi(obj.precision,'single'))
                    obj.B0map_thres=single(obj.B0map(squeeze(abs(obj.CoilSens(1,:,:,:)))>0));
                    obj.nLevels=single(ceil((1*max(abs(obj.B0map_thres(:)))*max(obj.tk)/pi)));
                    obj.wi=single(linspace(min((obj.B0map_thres(:))),1*max((obj.B0map_thres(:))),obj.nLevels));
                else
                    obj.B0map_thres=double(obj.B0map(squeeze(abs(obj.CoilSens(1,:,:,:)))>0));
                    obj.nLevels=double(ceil((1*max(abs(obj.B0map_thres(:)))*max(obj.tk)/pi)));
                    obj.wi=double(linspace(min(obj.B0map_thres(:)),1*max(abs(obj.B0map_thres(:))),obj.nLevels));
                    
                end
                                % calculate weights
                obj=obj.CalcWeights();
 
                end
                
            else
                error('Need several input parameters')
            end
        end
        
        function [out]=mtimes(obj,InData)
            if(obj.tflag==1) % Reverse operator: kspace data to image
                %hard coded para
                N=obj.NUFFTOP.imSize(1);
                [nCh,nFE,nIntlv,nPar]=size(InData);
                InData =permute(InData,[2 3 4 1]);
                img_MFI=zeros(N,N,nPar,nCh,length(obj.wi),obj.precision);
                b0term=exp(1i.*obj.tk(:)*obj.wi(:).');

                for idx_freq=1:length(obj.wi)
                        img_MFI(:,:,:,:,idx_freq)= obj.NUFFTOP'*double(bsxfun(@times,InData,b0term(:,idx_freq)));
                end
                img_MFI=fftshift(ifft(ifftshift(img_MFI,3),[],3),3).*(sqrt(nPar));
                out=sum(bsxfun(@times,permute((obj.MFI_weights),[2 3 4 5 1]),permute(img_MFI,[1 2 3 4 5])),5);
                if(~(isscalar(obj.CoilSens)||isempty(obj.CoilSens)))
                    %try to combine coils
                     out=sum(bsxfun(@times,out,permute(conj(obj.CoilSens),[2 3 4 1])),4);
                end
                
            else % forward operation image to kspace
                
                if(~(isscalar(obj.CoilSens)||isempty(obj.CoilSens)))
                    %try to combine coils
                    CoilIm=bsxfun(@times, InData,permute(obj.CoilSens,[2,3,4,1]));
                end

%                 kData=zeros(length(obj.tk),length(obj.wi));

                out=zeros([obj.NUFFTOP.dataSize size(obj.MFI_weights,4) size(obj.CoilSens,1)]);
                b0term=exp(-1i.*obj.tk(:)*obj.wi(:).');
                for idx_freq=1:length(obj.wi)
                        temp1=bsxfun(@times,squeeze(conj(obj.MFI_weights(idx_freq,:,:,:))),CoilIm);
                        temp1=fftshift(fft(ifftshift(temp1,3),[],3),3)./(sqrt(size(temp1,3)));
                        temp2= obj.NUFFTOP*double(temp1);
                        temp2=bsxfun(@times,b0term(:,idx_freq),temp2);
                        out=out+temp2;
                end
                out=permute(out,[4 1 2 3]); %CHAxCOLxLINxPAR
            end
        end
        function obj=ctranspose(obj)
            if(obj.tflag==0)
                obj.tflag=1;
            end
        end
        
        function obj=CalcWeights(obj)
            
            
            obj.MFI_weights=zeros([size(obj.B0map) obj.nLevels]);
            
            
            % calculate weights for all values in b0maps
            if(strcmp(obj.mode,'LeastSquares'))
                idx=find(squeeze(sum(abs(obj.CoilSens),1))>0);
                n_split=ceil(8*length(obj.tk)*numel(obj.B0map_thres)/1e9/obj.max_size);
                warning('Too large matrix to invert: spliting into %d  parts',n_split);
                obj.MFI_weights=zeros(numel(obj.wi),numel(obj.B0map),obj.precision);
                witk=single(exp(1i.*obj.tk*obj.wi));
                for i=1:n_split
                    B=exp(1i*obj.B0map_thres(i:n_split:end)*obj.tk.').'; %this is a HUGE matrix limited by max_size
                    obj.MFI_weights(:,idx(i:n_split:end))=mldivide((witk'*witk),(witk'*B));
                end
                obj.MFI_weights=reshape(obj.MFI_weights,[numel(obj.wi),size(obj.B0map)]);
 
                
                
            elseif(strcmp(obj.mode,'LeastSquares_HighMemory'))
                witk=exp(1i.*obj.tk*obj.wi);
                %memory hog way but faster
                B=exp(1i*obj.B0map(:)*obj.tk.').';
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