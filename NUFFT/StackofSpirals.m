classdef StackofSpirals

% Exmaple:
%     sig1=permute(r.sig,[2,3,1,4]);
%     t=MygpuNUFFT(r.NUFFT_obj.st.om'./(2*pi),col(r.NUFFT_obj.w).^2,2,5,12,r.NUFFT_obj.imSize,[],true);
%     imgg=t'*sig1;
%     imcoil=permute(fftshift(fft(imgg,[],4),4),[3,1,2,4]);
%     as(adaptiveCombine2(imcoil))
    
    properties
        op
        adjoint
        imSize
        dataSize
        w
        GridPara
        is2D
        Mode
    end
    methods
        function [obj] = StackofSpirals(k,w,imSize,sens,varargin)
            % function m = gpuNUFFT(k,w,osf,wg,sw,imageDim,sens,varargin)
            %
            %     k -- k-trajectory, scaled -0.5 to 0.5
            %          dims: 2 (3) ... x y (z)
            %                N ... # sample points
            %     w -- k-space weighting, density compensation
            %     osf -- oversampling factor (usually between 1 and 2)
            %     wg -- interpolation kernel width (usually 3 to 7)
            %     sw -- sector width to use
            %     imageDim -- image dimensions [n n n]
            %     sens -- coil sensitivity data
            %     varargin
            %        opt  -- true/false for atomic operation (default true)
            %             -- true/false for using textures on gpu (default true)
            %             -- true/false for balanced operation (default true)
            %
            %  res -- gpuNUFFT operator
            
            obj.GridPara.atomic = true;
            obj.GridPara.use_textures = true;
            obj.GridPara.balance_workload = true;
            obj.GridPara.SectorWidth=12;
            obj.GridPara.osf=1.5;
            obj.GridPara.KernelWidth=5;
            obj.Mode='CPU3D';
            
            
            
            obj.is2D = false;
            if (length(imSize) < 3)
                imSize(3) = 0;
                obj.is2D = true;
            end
            obj.imSize=imSize;
            obj.dataSize=ones(1,3);
            obj.dataSize(1:length(size(w)))=size(w);
            obj.adjoint = false;
            

            
            obj.w= sqrt(w(:));
            k=k(:,:);
            switch(obj.Mode)
                case 'GPU3D'
                    w = w(:);
                    obj.op.params.img_dims = uint32(imSize);
                    obj.op.params.osr = single(obj.GridPara.osf);
                    obj.op.params.kernel_width = uint32(obj.GridPara.KernelWidth);
                    obj.op.params.sector_width = uint32(obj.GridPara.SectorWidth);
                    obj.op.params.trajectory_length = uint32(size(k,2));
                    obj.op.params.use_textures = obj.GridPara.use_textures;
                    obj.op.params.balance_workload = obj.GridPara.balance_workload;
                    obj.op.params.is2d_processing = obj.is2D;
                    
                    [obj.op.dataIndices,obj.op.sectorDataCount,...
                        obj.op.densSorted,obj.op.coords,...
                        obj.op.sectorCenters,obj.op.sectorProcessingOrder,...
                        obj.op.deapoFunction] = mex_gpuNUFFT_precomp_f(single(k)',single(w(:))',obj.op.params);
                    obj.op.atomic = obj.GridPara.atomic;
                    obj.op.verbose = false;
                case 'CPU2D'
                    %use Fressler's CPU NUFFT
                    obj.GridPara.osf=2;
                    n_shift=obj.imSize(1:2)/2;
                    Kd=floor(obj.imSize(1:2)*obj.GridPara.osf);
                    Jd=ones(1,2)*obj.GridPara.KernelWidth;
                    obj.op = nufft_init((2*pi)*k(1:2,:)', obj.imSize(1:2), Jd, Kd, n_shift, 'minmax:kb');
                case 'CPU3D'
                    %use Fressler's CPU NUFFT
                    obj.GridPara.osf=2;
                    n_shift=obj.imSize/2;
                    Kd=floor(obj.imSize*obj.GridPara.osf);
                    Jd=ones(1,3)*obj.GridPara.KernelWidth;
                    obj.op = nufft_init((2*pi)*k', obj.imSize, Jd, Kd, n_shift, 'minmax:kb');
                  
            end
            
            
            if ~isempty(sens)
                obj.op.sensChn = size(sens,4);
                switch(obj.Mode)
                case 'GPU3D'
                         obj.op.sens = [real(sens(:))'; imag(sens(:))'];
                        obj.op.sens = reshape(obj.op.sens,[2 imSize(1)*imSize(2)*max(1,imSize(3)) obj.op.sensChn]);
                otherwise
                        obj.op.sens=sens;
                end
            else
                obj.op.sens = sens;
                obj.op.sensChn = 0;
            end
           
        end
        function res = mtimes(obj,bb)
            %             bb=[nFE,nIntlv,nPar,nCha]% adj case
            %             bb=[ImX,Imy,Imz,nCha] % forward case
            
            if (obj.adjoint)
                res=zeros([obj.imSize(1),size(bb,4)-obj.op.sensChn]);
                switch (obj.Mode)
                    case 'GPU3D'
                        res = gpuNUFFT_adj(obj.op,reshape(bb,[],size(bb,4)));
                    case 'CPU3D'
                        bb=bsxfun(@times,reshape(bb,[],size(bb,4)),obj.w(:));
                        res = nufft_adj(double(reshape(bb,[],size(bb,2))),obj.op);
                        if(~isempty(obj.op.sens))
                            res=sum(conj(obj.op.sens).*res,4);
                        end
                  
                end
            else
                
                switch (obj.Mode)
                    case 'GPU3D'
                        res = gpuNUFFT_forw(obj.op,squeeze(bb));
                    case 'CPU3D'
                        if(~isempty(obj.op.sens))
                            bb=bsxfun(@times,obj.op.sens,bb);
                        end
                        res = nufft(bb,obj.op);
                        res=bsxfun(@times,reshape(res,[],size(bb,4)),obj.w(:));
                end
                res=reshape(res,obj.dataSize(1),obj.dataSize(2),obj.dataSize(3),[]);
            end
            
        end
        
        function res = ctranspose(a)
            a.adjoint = xor(a.adjoint,1);
            res = a;
        end
        function res = transpose(a)
            a.adjoint = xor(a.adjoint,1);
            res = a;
        end
    end
end
