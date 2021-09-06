classdef StackofSpirals
%[obj] = StackofSpirals(k,w,imSize,sens,varargin)
%Interfaceclass for CPU and GPU-NUFFT
%INPUTS:
% k - kspace trajectory,scaled -0.5 to 0.5 [Axis x #points x NInterleaves x NPartitions]
% w - density compensation function [ #points x NInterleaves x NPartitions]
% imSize - [ X Y] or [X Y Z]
% sens - Coil sensitvities(can be empty) [ColxLinxParxCha]  or []
%
% varargin            : Name-Value pair
% 'CompMode'          : Computation mode {'GPU3D','CPU3D','CPU2DHybrid'} (default GPU3D)
% 'precision'         : {'single','double'} (default single)
% 'KernelWidth'       : scalar for GPU gridding kernel size/#neighbours in CPU case (default 5)
% 'osf'               : oversampling factor (default 1.5(GPU),2(CPU))
%ONLY relavent for GPU
% 'SectorWidth'       : default(12) (recommended: 8 for 3D,16 for 2D)
% 'atomic'            : true/false for atomic operation (default true)
% 'use_textures'      : true/false for using textures on gpu (default true)
% 'balance_workload' : true/false for balanced operation (default true)
%
%
%Example USAGE:
%SOSOP=StackofSpirals(kxyz,DCF3d,[FTOP.imSize 20],(permute(csm,[2 3 4 1])),...
%'CompMode','CPU2DHybrid','precision','single');
%im3d=SOSOP'*sig3d;%[ColxLinxParxCha] %Adjoint
%sigr=SOSOP*im3d;  %[ColxLinxPar] %forward
%im=spiralCGSENSE(SOSOP,sig3d,'maxit',10,'reg','none'); %cgsense
%
%
%
%Dependencies: 
%CPU: NUFFT functions from https://github.com/JeffFessler/mirt
%GPU: gpuNUFFT https://github.com/andyschwarzl/gpuNUFFT
%
%Author: praveen.ivp@gmail.com

    properties
        op
        adjoint
        imSize
        dataSize
        w
        GridPara
        is2D
        CompMode
        precision
    end
    methods
        function [obj] = StackofSpirals(k,w,imSize,sens,varargin)
            %constructor
            obj=obj.ParseInputParam(varargin{:});
            
            
            obj.is2D = false;
            if (length(imSize) < 3)
                imSize(3) = 0;
                obj.is2D = true;
            end
            obj.imSize=imSize;
            obj.dataSize=ones(1,3);
            obj.dataSize(1:length(size(w)))=size(w);
            obj.adjoint = false;
               
            obj.w= sqrt(w);
            switch(obj.CompMode)
                case 'GPU3D'
                    w = w(:);
                    k=k(:,:);
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
                case 'CPU2DHybrid'
                    obj.w= sqrt(w(:,:,1));
                    k=k(1:2,:,:,1);
%                     obj.precision='double';
                    %use Fressler's CPU NUFFT
                    obj.GridPara.osf=2;
                    n_shift=obj.imSize(1:2)/2;
                    Kd=floor(obj.imSize(1:2)*obj.GridPara.osf);
                    Jd=ones(1,2)*obj.GridPara.KernelWidth;
                    obj.op = nufft_init((2*pi)*k(:,:)', obj.imSize(1:2), Jd, Kd, n_shift, 'minmax:kb');
                case 'CPU3D'
                    k=k(:,:);
%                     obj.precision='double';
                    %use Fressler's CPU NUFFT
                    obj.GridPara.osf=2;
                    n_shift=obj.imSize/2;
                    Kd=floor(obj.imSize*obj.GridPara.osf);
                    Jd=ones(1,3)*obj.GridPara.KernelWidth;
                    obj.op = nufft_init((2*pi)*k', obj.imSize, Jd, Kd, n_shift, 'minmax:kb');
                  
            end

            
           if ~isempty(sens)
                obj.op.sensChn = size(sens,4);
                switch(obj.CompMode)
                case 'GPU3D'
                         obj.op.sens = [real(sens(:))'; imag(sens(:))'];
                        obj.op.sens = reshape(obj.op.sens,[2 imSize(1)*imSize(2)*max(1,imSize(3)) obj.op.sensChn]);
                otherwise
                        obj.op.sens=double(sens);
                end
            else
                obj.op.sens = sens;
                obj.op.sensChn = 0;
            end
           
        end
     function obj=ParseInputParam(obj,varargin)
            p=inputParser;
            p.KeepUnmatched=1;
             p.PartialMatching = true;
            addParameter(p,'CompMode','GPU3D',@(x) any(strcmp(x,{'GPU3D','CPU3D','CPU2DHybrid'})));
            addParameter(p,'precision','single',@(x) any(validatestring(x,{'single','double'})));
            parse(p,varargin{:});
            obj.precision=p.Results.precision; 
            obj.CompMode=p.Results.CompMode;
            p=inputParser;
            p.KeepUnmatched=1;
            addParameter(p,'atomic',true,@(x)islogical(x));
            addParameter(p,'use_textures',true,@(x) islogical(x));
            addParameter(p,'balance_workload',true,@(x) islogical(x));
            addParameter(p,'SectorWidth',12,@(x) isscalar(x));
            addParameter(p,'KernelWidth',5,@(x) isscalar(x));
            addParameter(p,'osf',1.5,@(x) isscalar(x));
            parse(p,varargin{:});
            obj.GridPara=p.Results;            
        end
        
        function res = mtimes(obj,bb)
            %             bb=[nFE,nIntlv,nPar,nCha]% adj case
            %             bb=[ImX,Imy,Imz,nCha] % forward case
              
            
            if (obj.adjoint)
                switch (obj.CompMode)
                    case 'GPU3D'
                        sz=size(bb);sz(4)=max(1,size(bb,4)-obj.op.sensChn);
                        res = gpuNUFFT_adj(obj.op,reshape(bb,[],size(bb,4)));
                        res=reshape(res,size(res,1),size(res,2),[],sz(4));
                    case 'CPU3D'
                        bb=bsxfun(@times,bb,obj.w);
                        res = nufft_adj(double(reshape(bb,[],size(bb,4))),obj.op);
                        res=res/(sqrt(numel(obj.w)));
                        if(~isempty(obj.op.sens))
                            res=sum(conj(obj.op.sens).*res,4);
                        end
                    case 'CPU2DHybrid'
                        sz=size(bb);
                        bb=bsxfun(@times,bb,obj.w);
                        res = nufft_adj(double(reshape(bb,prod(sz(1:2)),[])),obj.op);
                        res=res/(numel(obj.w));
                        res=reshape(res,[size(res,1),size(res,2),sz(3:end)]);
                        res=fftshift(ifft(ifftshift(res,3),[],3),3).*(sqrt(sz(3)));
                        if(~isempty(obj.op.sens))
                            res=sum(conj(obj.op.sens).*res,4);
                        end     
                end
            else
                
                switch (obj.CompMode)
                    case 'GPU3D'
                        res = gpuNUFFT_forw(obj.op,squeeze(bb));
                    case 'CPU3D'
                        if(~isempty(obj.op.sens))
                            bb=bsxfun(@times,double(obj.op.sens),double(bb));
                        end
                        res = nufft(bb,obj.op);
                        res=res/(sqrt(numel(obj.w)));
                        res=bsxfun(@times,reshape(res,[],size(bb,4)),obj.w(:));
                    case 'CPU2DHybrid'
                      if(~isempty(obj.op.sens))
                            bb=bsxfun(@times,double(obj.op.sens),double(bb));  
                      end
                      sz=size(bb);
                      bb=fftshift(fft(ifftshift(bb,3),[],3),3)./(sqrt(size(bb,3)));
                      res = nufft(reshape(bb,sz(1),sz(2),prod(sz(3:4))),obj.op);
                      res=res/sqrt(numel(obj.w));
                      res=reshape(res,obj.dataSize(1),obj.dataSize(2),obj.dataSize(3),[]); 
                      res=bsxfun(@times,res,obj.w);
                      
                        
                end
                res=reshape(res,obj.dataSize(1),obj.dataSize(2),obj.dataSize(3),[]);
            end
            if(~isa(res,obj.precision))
                if(strcmp(obj.precision,'single'))
                    res=single(res);
                else
                    res=double(res);
                end
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


%% construct SOSOP from LUSTIG's NUFFT

