classdef StackofSpiralsB0<StackofSpirals
    
    properties
        B0map
        B0Para
        tk
        Interp_weights
        mask
        
    end
    methods
        function obj=StackofSpiralsB0(k,w,imSize,sens,B0map,adcTime,varargin)
            obj@StackofSpirals(k,w,imSize,sens,varargin);
            obj.B0Para=struct('Method','MFI','Levels',[],'nLevels',0,'max_mem_GB',4,'mode','LeastSquares','precision','single');
            obj.tk=adcTime;
            
            obj.B0map=B0map;
            if(isempty(obj.op.sens)||isscalar(isempty(obj.op.sens)))
                obj.mask=ones(size(B0map),'logical');
            else
                switch(obj.Mode)
                    case 'GPU3D'
                       obj.mask = reshape(obj.op.sens(1,:,1),[imSize(1) imSize(2) max(1,imSize(3))])>0;         
                    otherwise
                        obj.mask=abs(obj.op.sens(:,:,:,1))>0;
                end
            end
            
            if(isscalar(obj.B0map)||isempty(obj.B0map))
                %if field map is scalar or empty, it becomes normal NUFFT operator
                
                obj.B0Para.nLevels=1;
                obj.B0Para.Levels=0;
                obj.Interp_weights=1;
                obj.B0map=0;
            else
                if(strcmpi(obj.B0Para.precision,'single'))
                    %                     obj.B0map=single(obj.B0map(obj.mask)); %(squeeze(abs(obj.CoilSens(1,:,:,:)))>0)
                    obj.B0Para.nLevels=single(ceil((1*max(abs(obj.B0map(:)))*max(obj.tk)/pi)));
                    % calculate weights
                    switch(obj.B0Para.Method)
                        case 'MFI'
                            obj.B0Para.Levels=single(linspace(min(obj.B0map(obj.mask)),1*max(obj.B0map(obj.mask)),obj.B0Para.nLevels));
                            obj=obj.CalcWeightsMFI();
                        case 'MTI'
                            obj.B0Para.Levels=linspace(0, max(obj.tk),obj.B0Para.nLevels);
                            obj=obj.CalcWeightsMTI();
                    end
                    
                else
                    obj.B0map=double(obj.B0map); %(squeeze(abs(obj.CoilSens(1,:,:,:)))>0)
                    obj.B0Para.nLevels=double(ceil((1*max(abs(obj.B0map(:)))*max(obj.tk)/pi)));
                    obj.B0Para.Levels=double(linspace(min(obj.B0map(:)),1*max(abs(obj.B0map(:))),obj.B0Para.nLevels));
                    
                end
                
                
            end
            
            
        end
        function [out]=mtimes(obj,bb)
            %             bb=[nFE,nIntlv,nPar,nCha]% adj case
            %             bb=[ImX,Imy,Imz,nCha] % forward case
            
            if(obj.adjoint==1) % Reverse operator: kspace data to image
                 all_images=zeros([obj.imSize(1:2),max(1,obj.imSize(3)),max(1,size(bb,4)-obj.op.sensChn),length(obj.B0Para.Levels)],obj.B0Para.precision);
                switch(obj.B0Para.Method)
                    case 'MFI'
                        b0term=exp(1i.*obj.tk(:)*obj.B0Para.Levels(:).');
                        for idx_freq=1:length(obj.B0Para.Levels)
                            all_images(:,:,:,:,idx_freq)= mtimes@StackofSpirals(obj,bsxfun(@times,bb,b0term(:,idx_freq)));
                        end
                        out=sum(bsxfun(@times,permute((obj.Interp_weights),[2 3 4 5 1]),all_images),5);
                    case 'MTI'
                        for idx_tau=1:obj.B0Para.nLevels
                            sig_MTI=bsxfun(@times,obj.Interp_weights(idx_tau,:).',bb);
                            all_images(:,:,:,:,idx_tau)= mtimes@StackofSpirals(obj,sig_MTI);
                            all_images(:,:,:,:,idx_tau)=bsxfun(@times,all_images(:,:,:,:,idx_tau),exp(1i*obj.B0map*obj.B0Para.Levels(idx_tau)));
                        end
                        out=sum(all_images,5);
                end
            else % forward operation image to kspace
                out=zeros([obj.dataSize max(size(bb,4),obj.op.sensChn)]);
                switch(obj.B0Para.Method)
                    case 'MFI'
                        b0term=exp(-1i.*obj.tk(:)*obj.B0Para.Levels(:).');
                        for idx_freq=1:obj.B0Para.nLevels
                            temp1=bsxfun(@times,squeeze(conj(obj.Interp_weights(idx_freq,:,:,:))),bb);
                            temp2= mtimes@StackofSpirals(obj,temp1);
                            temp2=bsxfun(@times,b0term(:,idx_freq),temp2);
                            out=out+temp2;
                        end
                    case 'MTI'
                        for idx_tau=1:obj.B0Para.nLevels
                            temp1=bsxfun(@times,bb,exp(-1i*obj.B0map*obj.B0Para.Levels(idx_tau)));
                            temp2= mtimes@StackofSpirals(obj,temp1);
                            temp2=bsxfun(@times,temp2,conj(obj.Interp_weights(idx_tau,:).'));
                            out=out+temp2;
                        end   
                end
            end
        end
        function obj=ctranspose(obj)
            if(obj.adjoint==0)
                obj.adjoint=1;
            end
        end
        
        function obj=CalcWeightsMFI(obj)
            obj.Interp_weights=zeros([size(obj.B0map) obj.B0Para.nLevels]);
            % calculate weights for all values in b0maps
            if(strcmp(obj.B0Para.mode,'LeastSquares'))
                idx=find(obj.mask);
                n_split=ceil(8*length(obj.tk)*numel(obj.B0map)/1e9/obj.B0Para.max_mem_GB);
                warning('Too large matrix to invert: spliting into %d  parts',n_split);
                obj.Interp_weights=zeros(numel(obj.B0Para.Levels),numel(obj.B0map),obj.B0Para.precision);
                witk=single(exp(1i.*obj.tk*obj.B0Para.Levels));
                for i=1:n_split
                    B=exp(1i*obj.tk*obj.B0map(i:n_split:end)); %this is a HUGE matrix limited by max_size
                    obj.Interp_weights(:,idx(i:n_split:end))=mldivide((witk'*witk),(witk'*B));
                end
                obj.Interp_weights=reshape(obj.Interp_weights,[numel(obj.B0Para.Levels),size(obj.B0map)]);
                
            elseif(strcmp(obj.B0Para.mode,'LeastSquares_HighMemory'))
                witk=exp(1i.*obj.tk*obj.B0Para.Levels);
                %memory hog way but faster
                B=exp(1i*obj.B0map(:)*obj.tk.').';
                obj.Interp_weights=mldivide((witk'*witk),(witk'*B));
                obj.Interp_weights=permute(reshape(obj.Interp_weights,[],size(obj.B0map,1),size(obj.B0map,2)),[2 3 1]);
                
            elseif(strcmp(obj.B0Para.mode,'NearestNeighbour') )
                for i=1:size(obj.Interp_weights,1)
                    for j=1:size(obj.Interp_weights,2)
                        [~,I]=min(abs(obj.B0Para.Levels-obj.B0map(i,j)));
                        obj.Interp_weights(i,j,I)=1;
                    end
                end
            end
            
        end
        function obj=CalcWeightsMTI(obj)
            if(sum(obj.B0map,'all')==0 ||isempty(obj.B0map))
                obj.Interp_weights=1;
            else
                %                 h=histogram(obj.B0map,"NumBins",prod(0.5*obj.NUFFTOP.imSize));
                %                 bins=h.BinEdges(h.BinCounts>0)+0.5*h.BinWidth;
                [Ncount,edges] = histcounts(obj.B0map(obj.mask),prod(0.5*obj.imSize(1:2))*max(1,obj.imSize(3)));
                bins=edges(Ncount>0)+0.5*diff(edges(1:2));
                wn_tau_l=exp(1i*(bins(:)*(obj.B0Para.Levels(:).')));
                if(strcmp(obj.B0Para.mode,'LeastSquares'))
                    actual=exp(1i*col(bins)*obj.tk(:)');
                    obj.Interp_weights=mldivide(wn_tau_l'*wn_tau_l,wn_tau_l'*actual);
                elseif(strcmp(obj.B0Para.mode,'NearestNeighbour') )
                    error('Not implemented')
                end
            end
        end
    end
    
end


