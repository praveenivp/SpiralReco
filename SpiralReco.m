classdef SpiralReco<handle
    properties
        KTraj
        Grad
        DCF
        NUFFT_obj
        soda_obj
        time
        
        B0Drift %temporal
        B0Map %spatial
        
        twix
        SpiralPara
        flags
        
        LoopCounter
        
        filename
        
        
        sig % raw data
        img % reconstructed image
        coilSens
        
        D % noise decoraltion matrix
        V % Coil compression matrix
        
    end
    methods
        
        %constructor: get full path of dat file or no arguments
        function obj=SpiralReco(varargin)
            if(nargin==0)
                %path='D:\Data\Spiral\20200918_B0test_sameres';
                path='M:\Subject_data\20201020_subject4847';
                [fn, pathname, ~] = uigetfile(strcat(path,'\*.dat'), 'Pick a DATA file');
                obj.filename=fullfile(pathname,fn);
                obj.twix = mapVBVD(obj.filename,'IgnoreSeg');
            else
                if(isfile(varargin{1}))
                    obj.filename=varargin{1};
                elseif(isfolder(varargin{1}))
                    path=varargin{1};
                    [fn, pathname, ~] = uigetfile(strcat(path,'\*.dat'), 'Pick a DATA file');
                    obj.filename=fullfile(pathname,fn);
                end
                obj.twix = mapVBVD(obj.filename, 'IgnoreSeg');
                if(iscell(obj.twix)) %VE
                    obj.twix=obj.twix{2};
                end
            end
            obj.SpiralPara=getSpiralPara(obj.twix);
            obj.flags=obj.getflags(varargin{2:end});
            obj.getSoda();
            obj.getNUFFTobj();
            obj.performRecon();
            
        end
        function flags=getflags(obj,varargin)
            switch nargin    
                case 2
                    flags=varargin{1};
                otherwise %parse name valur pairs
                    p=inputParser;
                    addParameter(p,'doGIRF',true,@(x) islogical(x));
                    addParameter(p,'doB0Driftcorr',false,@(x) islogical(x));
                    addParameter(p,'doCoilCombine','adapt2',@(x) any(validatestring(x,{'none','sos','adapt2'})));
                    addParameter(p,'doDCF','Jackson',@(x) any(validatestring(x,{'none','Jackson','voronoi'})));
                    addParameter(p,'doCoilCompression',false,@(x) (islogical(x)||islogical(x)));
                    addParameter(p,'doNoiseDecorr',true,@(x) islogical(x))
                    addParameter(p,'CoilSel',1:obj.twix.image.NCha, @(x) isvector(x));
                    addParameter(p,'RepSel',1:obj.twix.image.NRep, @(x) isvector(x));
                    addParameter(p,'SlcSel',1:obj.twix.image.NSli, @(x) (isvector(x)&& all(x<=obj.twix.image.NSli)));
                    addParameter(p,'is3D',(obj.twix.image.NPar>1),@(x)islogical(x));
                    addParameter(p,'isGIRFCorrected',false,@(x) islogical(x));
                    addParameter(p,'doB0Corr','none',@(x) any(validatestring(x,{'none','MTI','MFI','Iterative'})));
                    parse(p,varargin{:});
                    flags=p.Results;             
            end
            obj.SpiralPara.GradDelay=[1; 1; 1]*(obj.SpiralPara.GRAD_RASTER_TIME_DEFAULT-4.5); % 4.4us is filter delay of the ADC(2.4 us dwell)
        end
        
        function getNUFFTobj(obj)
            [obj.Grad,G_xyz,grad_MOM]=GetGradients(obj.twix,obj.SpiralPara,obj.soda_obj);
            
            if(obj.flags.doGIRF)
%                 load('PSF_time.mat','PSF_time')
                 load('S:\KYB\AGKS\pvalsala\Fieldcamera\20200210_GIRF_lowtrigdelay\processeddata\sid_5_PSF_full_reg500.mat','PSF_time')
                G_corr=(GIRF_Correction(G_xyz,PSF_time,'isCrossTermPSFCorr',true));
                obj.Grad=GradientXYZ2PRS(G_corr(:,2:4,:),obj.twix);
                obj.B0Drift=squeeze(G_corr(:,1,:));
                obj.flags.isGIRFCorrected=true;
            end
            obj.SpiralPara.grad_MOM=grad_MOM;
            [obj.KTraj,obj.time]=Grad2Traj(obj.Grad,obj.SpiralPara,'my');
            
%             obj.KTraj=obj.KTraj(2:(2*obj.SpiralPara.ADCLength+1),:);
            switch(obj.flags.doDCF)
                case 'none'
                    obj.DCF= ones(size(obj.KTraj));
                case 'Jackson'
                    obj.DCF=jacksonDCF2(obj.KTraj,obj.SpiralPara);
                case 'voronoi'
                    warning('not implemented')
                    obj.DCF=ones(size(obj.KTraj));
            end
            
            %scale kspace to -0.5 to 0.5 range
            kmax=2*pi*(0.5/(obj.SpiralPara.Resolution*1e-3));
            k_scaled=obj.KTraj./(2*kmax);
            
            N=obj.SpiralPara.FOV(1)/obj.SpiralPara.Resolution;
   % Lustif NUFFT operator splits the DCF in Forward and reverse
            % operator
%             obj.sig=bsxfun(@times, obj.sig,reshape(sqrt(obj.DCF),1,[]));
            obj.NUFFT_obj= NUFFT(col(k_scaled),(col(obj.DCF)),1,0,[N,N], 2);
            
        end
        function performFOVShift(obj)
            if(any(obj.SpiralPara.slice{obj.LoopCounter.cSlc}.Position ~=0))
                %             [~,idx]=sort(obj.SpiralPara.slice{1}.Normal);
                %             posi=1e-3*obj.SpiralPara.slice{1}.Position(idx); %m
                pos_PRS=GradientXYZ2PRS(1e-3*[1 -1 -1].*obj.SpiralPara.slice{obj.LoopCounter.cSlc}.Position,obj.twix); %only work for head first-supine
                
                if(obj.flags.doB0Driftcorr)
                    [kHO]=Grad2TrajHigherorder(obj.Grad,obj.SpiralPara);
                    B0_mod=exp(-1i*(real(obj.KTraj).*pos_PRS(1)+imag(obj.KTraj).*pos_PRS(2)+squeeze(kHO(:,3,:)).*pos_PRS(3) ));
                else
                B0_mod=exp(-1i*(real(obj.KTraj).*pos_PRS(1)+imag(obj.KTraj).*pos_PRS(2)));
                end
             %when using NUFFT operator from ESPIRIT toolbox check DCF in
             %all steps
                %B0_mod=B0_mod.*reshape(sqrt(obj.NUFFT_obj.w),size(B0_mod));
                
                
                

                % negative displacement added to real part of kTRaj moves up down in array show
                obj.sig=bsxfun(@times, obj.sig,reshape(B0_mod,1,[]));
            end
        end
        
        function performNoiseDecorr(obj)
            
            if(obj.flags.doNoiseDecorr)
                if(isempty(obj.D))
                    obj.calcNoiseDecorrMatrix();
                end
                if(ismatrix(obj.D) )
                    sz     = size(obj.sig);
                    obj.sig   = obj.D*obj.sig(:,:);
                    obj.sig    = reshape(obj.sig,sz);
                end
            end
        end
        
        function calcNoiseDecorrMatrix(obj)
            if isfield(obj.twix,'noise')
                noise                = permute(obj.twix.noise(:,obj.flags.CoilSel,:),[2,1,3]);
                noise                = noise(:,:).';
                R                    = cov(noise);
                R                    = R./mean(abs(diag(R)));
                R(eye(size(R,1))==1) = abs(diag(R));
                obj.D               = sqrtm(inv(R)).';
            else
                obj.D = 1;
            end
        end
        
        function performCoilCompression(obj)
            if(obj.flags.doCoilCompression>0)
                
                data=obj.sig(:,:).';
                [~,sval,vec]=svd(data,'econ');
                if(obj.flags.doCoilCompression<1 && obj.flags.doCoilCompression>0)
                    tol= obj.flags.doCoilCompression;
                else
                    tol = 0.95;
                end
                NcCha = find((cumsum(diag(sval))/sum(diag(sval))) >tol,1,'first');
                obj.V=vec(1:NcCha,:);
                fprintf('Compressing %d coils into %d coils with %d%% signal energy\n',size(sval,1),NcCha,tol*100);
                sz    = size(obj.sig);
                obj.sig   = obj.V*obj.sig(:,:);
                obj.sig   = reshape(obj.sig,[NcCha sz(2:end)]);
                
            end
        end
        
        function performRecon(obj)
            tic;
            print_str='';
            fprintf('\n');
            N=obj.SpiralPara.FOV(1)/obj.SpiralPara.Resolution;
            obj.img=zeros(length(obj.flags.CoilSel),N,N,obj.twix.image.NPar,max(obj.flags.SlcSel),max(obj.flags.RepSel));
            obj.coilSens=zeros(length(obj.flags.CoilSel),N,N,obj.twix.image.NPar,max(obj.flags.SlcSel));
            for cRep=obj.flags.RepSel
                for cSlc=obj.flags.SlcSel
                    
                    obj.LoopCounter.cRep=cRep;
                    obj.LoopCounter.cSlc=cSlc;
                    
                    %             {'Col','Cha','Lin','Par','Sli','Ave','Phs','Eco','Rep',
                    %     'Set','Seg','Ida','Idb','Idc','Idd','Ide'}
                    
                    obj.sig = obj.twix.image(:, obj.flags.CoilSel,:,:, obj.LoopCounter.cSlc, 1,1,1,obj.LoopCounter.cRep);
                    obj.sig=permute(obj.sig,[2 1 3 4]);
                    [nCh,nFE,nIntlv,nPar]=size(obj.sig);
%                     obj.sig(:,end-128:end,:)=0;  %stupid line
                    obj.sig = reshape(obj.sig,[nCh,nFE*nIntlv,nPar]);
                    
                    obj.performNoiseDecorr();
                    obj.performCoilCompression();
                    obj.performFOVShift();
                    obj.performNUFFT();
                    obj.performCoilCombination();
                    
                    fprintf(repmat('\b',1,numel(print_str)));
                    print_str = sprintf(['slice = %2.0f /%2.0f)'...
                                         ' | Rep = %2.0f/%2.0f'...
                                         ' | time = %6.1f s'], cSlc,max(obj.flags.SlcSel),cRep,max(obj.flags.RepSel),toc);
                    fprintf(print_str);
                    
                end
            end
            
            if(~strcmpi(obj.flags.doCoilCombine,'none'))
                obj.img(2:end,:,:,:,:,:,:)=[];
            end
            fprintf('\n');
        end
        
        
        function performNUFFT(obj)
            
            if(obj.flags.is3D)
                for par=1:obj.twix.image.NPar
                    for cha=1:size(obj.sig,1)
                        obj.img(cha,:,:,par,obj.LoopCounter.cSlc,obj.LoopCounter.cRep) = obj.NUFFT_obj'*double(obj.sig(cha,:,par));
                    end
                end
                obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep)=(fftshift(fft(fftshift(obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep),4),[],4),4));
                
            else
                
                for cha=1:size(obj.sig,1)
                    obj.img(cha,:,:,1,obj.LoopCounter.cSlc,obj.LoopCounter.cRep) = obj.NUFFT_obj'*double(obj.sig(cha,:));
                end
            end
        end
        
        function performCoilCombination(obj)
            if(~(size(obj.img,1)>1))
                obj.flags.doCoilCombine='none';
            end
            %coil combination
            switch(obj.flags.doCoilCombine) %{'none','sos','adapt','adapt2'}
                case 'sos'
                    obj.img(1,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep) = sqrt(sum(abs(obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep)).^2,1));
                case 'adapt2'
                    if(obj.LoopCounter.cRep==1 && ~any(obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc),'all') )
                    [obj.img(1,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep), ...
                        obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc)]...
                        = adaptiveCombine2(obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep));
                    else
                        obj.img(1,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep)=sum(obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep).*...
                                                                                    obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc),1);                                                                 
                    end
                    
                otherwise
                    warning('Exiting without coil combination')
            end
        end
        
        function obj=getSoda(obj)
            obj.soda_obj=   SODA_OBJ( 'mrprot',obj.twix.hdr);
            %             ReadoutFOV
            
            obj.soda_obj.NPixelReadout=obj.SpiralPara.FOV(1)/obj.SpiralPara.Resolution;
            obj.soda_obj.NPixelPhase=obj.SpiralPara.FOV(1)/obj.SpiralPara.Resolution;
            obj.soda_obj.PixelSizePhase= obj.soda_obj.ReadoutFOV/obj.soda_obj.NPixelPhase;
            obj.soda_obj.PixelSizeReadout= obj.soda_obj.PhaseFOV/obj.soda_obj.NPixelReadout;
            obj.soda_obj=obj.soda_obj.calcPixLoc();
            
            %             PhaseFOV
            
        end
        
        function obj=setDataDelay(obj,datadelay_us)
            %usually in the range 3-10 us  
            obj.SpiralPara.GradDelay=[1; 1; 1].*datadelay_us;
            obj.getNUFFTobj();
            obj.performRecon();
        end
        function WriteImages(obj)
            [fPath,fn,~]=fileparts(obj.filename);
            if(~isfolder(fullfile(fPath,'processeddata')))
                 mkdir(fullfile(fPath,'processeddata'))
            end
            description=strcat(obj.SpiralPara.SpiralTypeName,'_i',num2str(obj.SpiralPara.Ninterleaves));
            MyWriteNIFTI(squeeze(abs(obj.img)),fullfile(fPath,'\processeddata\',fn),description);
        
            
        end
        
    end
end


