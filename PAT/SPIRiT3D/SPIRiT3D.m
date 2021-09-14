classdef SPIRiT3D
    %Kernel consitency operator for SPIRIT reconstruction. It is at the
    %moment tested only for 3D and works only on image domain. Therefore,
    %it could be a memory hog for large matrix size and coils.
    %Constructor: SPIRiT3D(InpArg1,imSize,method)
    %
    %INPUTS:
    %InpArg1: SpiralReco obj or calculated kernel[kSize1 x KSize2 x Ksize3 x coil x coil]. 
    %imSize : final image size [COLxLINxPAR].optinal if InpArg1 is SpiralReco obj
    %method: default('image'), 'conv' coming soon
    % 
    %
    %Author:praveen.ivp@gmail.com
    properties
        kernel
        KERNEL
        method
        imSize
        tflag
    end
    methods
        function obj= SPIRiT3D(InpArg1,imSize,method)
            if(isnumeric(InpArg1))
            obj.kernel=InpArg1;
            obj.imSize=imSize;
            else
                if(~isfield(InpArg1.twix,'refscan'))
                    error('No Calibration data for the SPIRIT Kernel:(.')
                end
                ref=removeOS(squeeze(InpArg1.twix.refscan(:,InpArg1.flags.CoilSel,:,:,:,:)),1,2);
                ref=permute(ref,[1 3 4 2]);
                ref_sz=size(ref);
                Ncha=size(ref,4);
                kSize=[7 7 7];
                CalibSize=[32, 32 16];
                CalibSize=min(CalibSize,ref_sz(1:3));
                kSize=min(CalibSize,kSize);
                ref=crop(ref,[CalibSize Ncha]);
                ref=flip(flip(ref,1),2);
                obj.kernel = calibSPIRiT3D(ref, kSize, 0.02);
                obj.imSize=max(1,InpArg1.NUFFT_obj.imSize);
            end
            
            if(nargin<3)
                obj.method='image';
            else
            obj.method=method;
            end
           
            obj.tflag=1;
            obj=obj.precompute();
        end
        function obj=precompute(obj)
            if(strcmpi(obj.method,'image'))
                ksz=size(obj.kernel);           
                obj.KERNEL=zeros([obj.imSize ksz(4:5)],'single');
                for n=1:size(obj.kernel,5)
                    obj.KERNEL(:,:,:,:,n)=zpad(obj.kernel(end:-1:1,end:-1:1,end:-1:1,:,n)*sqrt(prod(obj.imSize)), [obj.imSize ,ksz(4)]);
                end

                    obj.KERNEL=myifft(obj.KERNEL,[1 2 3],true,obj.imSize);

            else
                error('not implemented: use only image')
            end
            
        end
        
        
        function out=mtimes(obj,InData)
            Ncha=size(obj.kernel,4);
%             kSize=size(obj.kernel); kSize=kSize(1:3);          
            switch obj.method
                case {'image'}
                    out = zeros(size(InData));
                    if ~(obj.tflag)
                        for n=1:Ncha
                            tmpk = permute(conj(obj.KERNEL(:,:,:,n,:)),[1 2 3 5 4]);
                            out(:,:,:,n) = sum(tmpk.*InData,4);
                        end
                        out = out - InData;
                    else
                        for n=1:Ncha
                            tmpk = obj.KERNEL(:,:,:,:,n);
                            out(:,:,:,n) = sum(tmpk.*InData,4);
                        end
                        out = out-InData;
                    end
            end
        end
    function obj=ctranspose(obj)
    if(obj.tflag==0)
        obj.tflag=1;
    end
    end
    
    
end

end


%% supportiing fucntions
function [im]=myfft(kdata1,dim,shift,sz)
switch (nargin)
    case 1
        dim=1;
        shift=true;
        sz=size(kdata1);
    case 2
        shift=true;
        sz=size(kdata1);
    case 3      
        sz=size(kdata1);
end

im=kdata1;
if(shift)
    for i=dim
        im=fftshift(fft(ifftshift(im,i),sz(i),i),i)./sqrt(sz(i));
    end
else
    for i=dim
        im=fft(im,sz(i),i);
    end
end

end

function [kdata]=myifft(im,dim,shift,sz)
switch (nargin)
    case 1
        dim=1;
        shift=true;
        sz=size(im);
    case 2
        shift=true;
        sz=size(im);
    case 3    
        sz=size(im);  
end

kdata=im;
if(shift)
    for i=dim
        kdata=fftshift(ifft(ifftshift(kdata,i),sz(i),i),i).*sqrt(sz(i));
    end
else
    for i=dim
        kdata=ifft(kdata,sz(i),i);
    end
end

end
