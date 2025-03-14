classdef SpiralReco<matlab.mixin.Copyable
    properties
        KTraj % [rad/m]
        Grad  % [mT/m]
        DCF
        NUFFT_obj
        SPIRIT3D_obj %kernel consitency operator
        soda_obj
        time  %[s]
        
        B0Drift %temporal k0 [rad]
        B0Map %spatial [rad/s]
        
        twix
        SpiralPara
        flags
        
        LoopCounter
        
        filename
        
        
        sig % raw data
        img % reconstructed image [CHAxCOLxLINxPARxSLCxREP]
        coilSens %%[CHAxCOLxLINxPARxSLC]
        coilNormMat %%[COLxLINxPARxSLC]
        
        D % noise decoraltion matrix
        V % Coil compression matrix
        
    end
    methods
        
        %constructor: get full path of dat file or no arguments
        function obj=SpiralReco(varargin)
            if(nargin==0)
                %path='D:\Data\Spiral\20200918_B0test_sameres';
                path=userpath;
                [fn, pathname, ~] = uigetfile(strcat(path,'\*.dat'), 'Pick a DATA file');
                obj.filename=fullfile(pathname,fn);
                obj.twix = mapVBVD(obj.filename,'IgnoreSeg');
            else
                if(isstruct(varargin{1})||iscell(varargin{1}))
                    obj.twix =varargin{1};
                elseif(isfile(varargin{1}))
                    obj.filename=varargin{1};
                    obj.twix = mapVBVD(obj.filename, 'IgnoreSeg');
                elseif(isfolder(varargin{1}))
                    path=varargin{1};
                    [fn, pathname, ~] = uigetfile(strcat(path,'\*.dat'), 'Pick a DATA file');
                    obj.filename=fullfile(pathname,fn);
                    obj.twix = mapVBVD(obj.filename, 'IgnoreSeg');
                else
                    error('File not  found or invalid input');
                end
            end
            
            obj.coilSens=[];
            obj.SPIRIT3D_obj=[];
            if(iscell(obj.twix)) %VE  
                obj.SpiralPara=getSpiralPara(obj.twix{end});
                TW1=obj.twix{1};
                obj.twix=obj.twix{end};
                obj.getflags(varargin{2:end});
                obj.calcNoiseDecorrMatrix(TW1);
                
            else
                obj.SpiralPara=getSpiralPara(obj.twix);
                obj.getflags(varargin{2:end});
            end
            
            obj.getSoda();
            obj.performRecon();
            
        end
        function getflags(obj,varargin)
            switch nargin    
                case 2
                    obj.flags=varargin{1};
                otherwise %parse name valur pairs
                    p=inputParser;
                    p.KeepUnmatched=1;
                    addParameter(p,'doGIRF',true,@(x) islogical(x));
                    addParameter(p,'doB0Driftcorr',false,@(x) islogical(x));
                    addParameter(p,'doCoilCombine','sos',@(x) any(strcmp(x,{'none','sos','adapt2'})));
                    addParameter(p,'doDCF','Jackson',@(x) any(strcmp(x,{'none','Jackson','voronoi'})));
                    addParameter(p,'doCoilCompression',false,@(x) (islogical(x)||isscalar(x)));
                    addParameter(p,'doNoiseDecorr',true,@(x) islogical(x))
                    addParameter(p,'NormNoiseData',false,@(x) islogical(x))
                    addParameter(p,'CoilSel',1:obj.twix.image.NCha, @(x) isvector(x));
                    addParameter(p,'RepSel',1: (obj.twix.image.NRep*obj.twix.image.NSet*obj.twix.image.NAve), @(x) isvector(x));
                    addParameter(p,'SlcSel',1:obj.twix.image.NSli, @(x) (isvector(x)&& all(x<=obj.twix.image.NSli)));
                    addParameter(p,'is3D',(obj.twix.image.NPar>1),@(x)islogical(x));
                    addParameter(p,'isGIRFCorrected',false,@(x) islogical(x));
                    addParameter(p,'doB0Corr','none',@(x) any(strcmp(x,{'none','MTI','MFI'})));
                    addParameter(p,'doPAT','none',@(x) any(strcmp(x,{'none','CGSENSE','SPIRIT'})))
                    addParameter(p,'precision','single',@(x) any(strcmp(x,{'single','double'})));
                    addParameter(p,'CompMode','CPU2DHybrid',@(x) any(strcmp(x,{'GPU3D','CPU3D','CPU2DHybrid'})));
                    addParameter(p,'maxit',10,@(x)isscalar(x));
                    addParameter(p,'tol',1e-6,@(x)isscalar(x));
                    addParameter(p,'reg','none',@(x) any(strcmp(x,{'none','Tikhonov'})));
                    addParameter(p,'reg_lambda',0,@(x)isscalar(x));

                    
                    parse(p,varargin{:});
                    
                    obj.flags=p.Results;   
                    if(isfield(p.Unmatched,'csm')) %[CHAxCOLxLINxPARxSLC]
                        obj.coilSens=p.Unmatched.csm;
                        obj.flags.doCoilCombine='none';
%                         obj.SpiralPara.R_PE=2;
%                         obj.flags.doPAT='CGSENSE';
                    end
                    if(isfield(p.Unmatched,'fm'))
                        obj.B0Map=p.Unmatched.fm;
%                         obj.flags.doB0Corr='MTI';
                    end
                    if(isfield(p.Unmatched,'GradDelay'))
                        obj.SpiralPara.GradDelay=ones(1,3)*p.Unmatched.GradDelay;
%                         obj.flags.doB0Corr='MTI';
                    end
            end
            
        end
        
        function getNUFFTobj(obj)
            [obj.Grad,G_xyz,grad_MOM]=GetGradients(obj.twix,obj.SpiralPara,obj.soda_obj,obj.LoopCounter.cSlc);
            obj.SpiralPara.grad_MOM=grad_MOM;
            if(obj.flags.doGIRF)
                [SpiralRecopath,~,~]=fileparts(mfilename('fullpath'));
                if(obj.twix.hdr.Dicom.lFrequency<300e6) %should calc PSF in freq for low field!!
                    load(fullfile(SpiralRecopath,'kspace','PSF_time_oct2021_3T.mat'),'PSF_time')
                    G_corr_SPH=GIRF_Correction(G_xyz,PSF_time,'isCrossTermPSFCorr',true);
                else %9.4T
                    load(fullfile(SpiralRecopath,'kspace','PSF_MAR2022'),'PSF')
                    G_corr_SPH=GIRF_correction_Freq(G_xyz,PSF,'isCrossTermPSFCorr',true);
                end
                obj.Grad=GradientXYZ2PRS(G_corr_SPH(:,2:4,:),obj.soda_obj);
                [k_SPH,obj.time]=Grad2TrajHigherorder(G_corr_SPH,obj.SpiralPara);
                k_PRS=GradientXYZ2PRS(k_SPH(:,2:4,:),obj.soda_obj);
                [~,siemensB0_ECC]=getEddyB0driftIIR(G_xyz,obj.twix,'IIR');
                obj.B0Drift=squeeze(k_SPH(:,1,:))-siemensB0_ECC;
                obj.flags.isGIRFCorrected=true;
            else
                [k_PRS,obj.time]=Grad2TrajHigherorder(obj.Grad,obj.SpiralPara);
            end
            switch(obj.flags.doDCF)
                case 'none'
                    obj.DCF= ones([size(k_PRS,1) size(k_PRS,3)]);
                case 'Jackson'
                    obj.DCF=jacksonDCF2(squeeze(complex(k_PRS(:,1,:),k_PRS(:,2,:))),obj.SpiralPara);
                case 'voronoi'
                    obj.DCF=VoronoiDCF2D(squeeze(complex(k_PRS(:,1,:),k_PRS(:,2,:))));
            end

            %scale kspace to -0.5 to 0.5 range
            kmax=2*pi*(0.5/(obj.SpiralPara.Resolution*1e-3));
            kmax=[kmax;kmax;(2*pi*0.5)/(1e-3*obj.SpiralPara.FOV_PRS(3)/obj.SpiralPara.NPartitions)];
            k_PRS=permute(k_PRS,[2 1 3]);
            k_PRS=repmat(k_PRS,[1 1 1 obj.SpiralPara.NPartitions]);
            kz=(-0.5:1/obj.SpiralPara.NPartitions:(0.5-1/obj.SpiralPara.NPartitions))*(2*kmax(3));
            kz=repmat(permute(kz(:),[2 3 4 1]),[1 size(k_PRS,2) size(k_PRS,3) 1]);
            k_PRS(3,:,:,:)=k_PRS(3,:,:,:)+kz;

            N=round(obj.SpiralPara.FOV(1)/obj.SpiralPara.Resolution);

            if (isempty(obj.coilSens)|| ~any(col(abs(obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc)))>0))
                csm=[];
            else
                csm=permute(obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc),[2 3 4 1]);
            end
            acq_sel=(obj.twix.image.Rep==obj.twix.image.NRep&obj.twix.image.Sli==1);
            Lin_ordering=reshape(obj.twix.image.Lin(acq_sel),round(obj.SpiralPara.Ninterleaves/obj.SpiralPara.R_PE),[]);
            Par_ordering=reshape(obj.twix.image.Par(acq_sel),round(obj.SpiralPara.Ninterleaves/obj.SpiralPara.R_PE),[]);
            %              Lin_ordering=Lin_ordering(1:3:end,:);Par_ordering=Par_ordering(1:3:end,:);

            obj.KTraj=k_PRS(:,:,sub2ind([size(k_PRS,3) size(k_PRS,4)],Lin_ordering,Par_ordering));
            obj.KTraj=reshape(obj.KTraj,[size(obj.KTraj,1),size(obj.KTraj,2),size(Lin_ordering)]);
            obj.DCF=repmat(obj.DCF,[1  1 obj.SpiralPara.NPartitions]);
            obj.DCF=obj.DCF(:,sub2ind([size(obj.DCF,2) size(obj.DCF,3)],Lin_ordering,Par_ordering));
            obj.DCF=reshape(obj.DCF,[size(obj.DCF,1),size(Lin_ordering)]);

            %                 if(obj.flags.is3D)
            switch(obj.flags.doB0Corr)
                case 'none'
                    obj.NUFFT_obj= StackofSpirals(obj.KTraj./(2*kmax),obj.DCF,[N N obj.SpiralPara.NPartitions],csm,...
                        'CompMode',obj.flags.CompMode,'precision',obj.flags.precision);
                case 'MFI'
                    obj.NUFFT_obj= StackofSpiralsB0(obj.KTraj./(2*kmax),obj.DCF,[N N obj.SpiralPara.NPartitions],...
                        csm,obj.B0Map,obj.time*1e-6,'Method','MFI','CompMode',obj.flags.CompMode,...
                        'precision',obj.flags.precision);
                case 'MTI'
                    obj.NUFFT_obj= StackofSpiralsB0(obj.KTraj./(2*kmax),obj.DCF,[N N obj.SpiralPara.NPartitions],...
                        csm,obj.B0Map,obj.time*1e-6,'Method','MTI','CompMode',obj.flags.CompMode,...
                        'precision',obj.flags.precision);
            end

            %                 else
            %                      obj.NUFFT_obj= StackofSpirals(obj.KTraj./(2*kmax),obj.DCF,[N N],csm,...
            %                          'CompMode',obj.flags.CompMode,'precision',obj.flags.precision);
            %                 end
        end
        function performFOVShift(obj)
            if(any(obj.SpiralPara.slice{obj.LoopCounter.cSlc}.Position ~=0))
                pos_PRS=GradientXYZ2PRS(1e-3*[-1 1 1].*obj.SpiralPara.slice{obj.LoopCounter.cSlc}.Position,obj.soda_obj,obj.LoopCounter.cSlc); %only work for head first-supine
                pos_PRS=[pos_PRS(1); pos_PRS(2);0*pos_PRS(3)];
                %                 if(strcmpi(obj.flags.CompMode,'GPU3D')&& obj.flags.is3D)
                %                     pos_PRS(3)=1e-3*obj.twix.hdr.Phoenix.sKSpace.dSliceResolution;%m
                %                 end
                B0_mod=exp(1i*sum(bsxfun(@times,obj.KTraj,pos_PRS(:)),1));
            else
                B0_mod=ones([1 size(obj.DCF)]);
            end
            acq_sel=(obj.twix.image.Rep==obj.twix.image.NRep&obj.twix.image.Sli==1);
            Lin_ordering=reshape(obj.twix.image.Lin(acq_sel),round(obj.SpiralPara.Ninterleaves/obj.SpiralPara.R_PE),[]);
            Par_ordering=reshape(obj.twix.image.Par(acq_sel),round(obj.SpiralPara.Ninterleaves/obj.SpiralPara.R_PE),[]);
            %             Lin_ordering=Lin_ordering(1:3:end,:);Par_ordering=Par_ordering(1:3:end,:);

            if(obj.flags.doB0Driftcorr)
                obj.B0Drift=repmat(obj.B0Drift,[1  1 obj.SpiralPara.NPartitions]);
                obj.B0Drift=obj.B0Drift(:,sub2ind([size(obj.B0Drift,2) size(obj.B0Drift,3)],Lin_ordering,Par_ordering));
                obj.B0Drift=reshape(obj.B0Drift,[1 size(obj.B0Drift,1),size(Lin_ordering)]);
                B0_mod=B0_mod.*exp(-1i*obj.B0Drift);
            end
            %do FOVshift and DCF compensation together
            B0_mod=B0_mod.*reshape(sqrt(obj.DCF),size(B0_mod));


            obj.sig=obj.sig(:,:,sub2ind([size(obj.sig,3) size(obj.sig,4)],Lin_ordering,Par_ordering));
            obj.sig=reshape(obj.sig,[size(obj.sig,1) ,size(obj.sig,2),size(Lin_ordering)]);
            obj.sig=bsxfun(@times, obj.sig,B0_mod);

        end

        function performNoiseDecorr(obj)

            if(obj.flags.doNoiseDecorr)
                if(isempty(obj.D))
                    obj.calcNoiseDecorrMatrix(obj.twix);
                end
                if(ismatrix(obj.D) )
                    sz     = size(obj.sig);
                    obj.sig   = obj.D*obj.sig(:,:);
                    obj.sig    = reshape(obj.sig,sz);
                end
            end
        end

        function calcNoiseDecorrMatrix(obj,twix)
            if isfield(twix,'noise')
                noise                = permute(twix.noise(:,obj.flags.CoilSel,:),[2,1,3]);
                noise                = noise(:,:).';
                R                    = cov(noise);
                R(eye(size(R,1))==1) = abs(diag(R));
                if(obj.flags.NormNoiseData)
                    R= R./mean(abs(diag(R)));
                    obj.D               = sqrtm(inv(R)).';
                else
                    scale_factor=1; %dwell time are the same
                    Rinv = inv(chol(R,'lower')).';
                    obj.D = Rinv*sqrt(2)*sqrt(scale_factor);

                end
            else
                obj.D = 1;
            end
        end
        
        function performCoilCompression(obj)
            if(obj.flags.doCoilCompression>0)
                if(isempty(obj.V))
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
                else
                   NcCha =size(obj.V,1); 
                end
                sz    = size(obj.sig);
                obj.sig   = obj.V*obj.sig(:,:);
                obj.sig   = reshape(obj.sig,[NcCha sz(2:end)]);
                
                %fix coil dimension
                obj.img(NcCha+1:end,:,:,:,:,:,:)=[];
                obj.coilSens(NcCha+1:end,:,:,:,:,:,:)=[];
            end
        end
        
        function performRecon(obj)
            tic;
            print_str='';
            fprintf('\n');
            N=round(obj.SpiralPara.FOV(1)/obj.SpiralPara.Resolution);
            
            if(isempty(obj.coilSens))
            obj.img=zeros(length(obj.flags.CoilSel),N,N,obj.SpiralPara.NPartitions,max(obj.flags.SlcSel),max(obj.flags.RepSel),obj.flags.precision);
            obj.coilSens=zeros(length(obj.flags.CoilSel),N,N,obj.SpiralPara.NPartitions,max(obj.flags.SlcSel),obj.flags.precision);
            else
                obj.img=zeros(1,N,N,obj.SpiralPara.NPartitions,max(obj.flags.SlcSel),max(obj.flags.RepSel),obj.flags.precision);
            end
            [aveidx,repidx,setidx]=ndgrid(1:obj.twix.image.NAve,1:obj.twix.image.NRep,1:obj.twix.image.NSet);
            repidx=int16(repidx(:));setidx=int16(setidx(:)); aveidx=int16(aveidx(:));
            
            
            for cSlc=obj.flags.SlcSel
                obj.LoopCounter.cSlc=cSlc;
                obj.getNUFFTobj();
                for cRep=obj.flags.RepSel
                    obj.LoopCounter.cRep=cRep;
                    
                    %             {'Col','Cha','Lin','Par','Sli','Ave','Phs','Eco','Rep',
                    %     'Set','Seg','Ida','Idb','Idc','Idd','Ide'}
                    obj.sig = obj.twix.image(:, obj.flags.CoilSel,:,:, obj.LoopCounter.cSlc, aveidx(cRep),1,1,repidx(cRep),setidx(cRep));
                    
                    
                    obj.sig=permute(obj.sig,[2 1 3 4]);         
                    obj.performNoiseDecorr();
                    obj.performCoilCompression();
                    obj.performFOVShift();
                    obj.performNUFFT();
                    obj.performCoilCombination();
                    
                    fprintf(repmat('\b',1,numel(print_str)));
                    print_str = sprintf(['slice = %2.0f/%2.0f'...
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
%                     if(obj.SpiralPara.CAIPIShift==0)
                        switch(obj.flags.doPAT)
                            case 'none'
                                 obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep) = permute(obj.NUFFT_obj'*(permute(obj.sig,[2,3,4,1])),[4,1,2,3]);
                            case 'CGSENSE'
                                im_pat=spiralCGSENSE(obj.NUFFT_obj,permute(obj.sig,[2,3,4,1]),...
                                    'maxit',obj.flags.maxit,'tol',obj.flags.tol,'reg',obj.flags.reg,...
                                    'lambda',obj.flags.reg_lambda);
                                obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep) = permute(im_pat,[4,1,2,3]);
                            case 'SPIRIT'
                                if(isempty(obj.SPIRIT3D_obj))
                                  obj.SPIRIT3D_obj=SPIRiT3D(obj);
                                end
                                img_spirit=spiralCGSPIRiT(obj.NUFFT_obj,obj.SPIRIT3D_obj,permute(obj.sig,[2,3,4,1])...
                                    ,'maxit',obj.flags.maxit,'lambda',1);
                                obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep) = permute(img_spirit,[4,1,2,3]);
                                
                        end
%                     else
%                         curr_CAIPI=mod(cpar-1,obj.SpiralPara.R_PE)+1;
%                         obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep) =permute(obj.NUFFT_obj{curr_CAIPI}'*double(permute(obj.sig,[2,3,4,1])),[4,1,2,3]); %sig is sqrt(DCF) pre-compensated
%                     end
                
            else
                         switch(obj.flags.doPAT)
                            case 'none'
                                 obj.img(:,:,:,1,obj.LoopCounter.cSlc,obj.LoopCounter.cRep) = permute(obj.NUFFT_obj'*permute(double(obj.sig),[2, 3, 4,1]),[4,1,2,3]);
                            case 'CGSENSE'
                                im_pat=spiralCGSENSE(obj.NUFFT_obj,permute(obj.sig,[2,3,4,1]),...
                                    'maxit',obj.flags.maxit,'tol',obj.flags.tol,'reg',obj.flags.reg,...
                                    'lambda',obj.flags.reg_lambda);
                                obj.img(:,:,:,1,obj.LoopCounter.cSlc,obj.LoopCounter.cRep) = permute(im_pat,[4,1,2,3]);
                            case 'SPIRIT'
                                if(isempty(obj.SPIRIT3D_obj))
                                  obj.SPIRIT3D_obj=SPIRiT3D(obj);
                                end
                                img_spirit=spiralCGSPIRiT(obj.NUFFT_obj,obj.SPIRIT3D_obj,permute(obj.sig,[2,3,4,1])...
                                    ,'maxit',obj.flags.maxit,'lambda',1);
                                obj.img(:,:,:,1,obj.LoopCounter.cSlc,obj.LoopCounter.cRep) = permute(img_spirit,[4,1,2,3]);
                                
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
                        obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc),obj.coilNormMat]...
                        = adaptiveCombine2(obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep),[],false);
                    obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc)=conj(obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc));
                    else
                        obj.img(1,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep)=sum(obj.img(:,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep).*...
                                                                                    conj(obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc)),1);                                                                 
                    end
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
            obj.performRecon();
        end
      
        function WriteImages(obj)
            [fPath,fn,~]=fileparts(obj.filename);
            try
                if(~isfolder(fullfile(fPath,'processeddata')))
                    mkdir(fullfile(fPath,'processeddata'))
                end
            catch
                fPath=pwd;
                if(~isfolder(fullfile(fPath,'processeddata')))
                    mkdir(fullfile(fPath,'processeddata'))
                end
            end
            description=strcat(obj.SpiralPara.SpiralTypeName,'_i',num2str(obj.SpiralPara.Ninterleaves));
            MyNIFTIWriteSpiral(squeeze(single(abs(obj.img))),obj.twix,fullfile(fPath,'processeddata',fn),description);
        
            
        end
        function SaveResults(obj)
            % save results toa MAT file
            im=squeeze(obj.img);
            
            flags=obj.flags;
            
            OutFile=sprintf('m%d_B0%s_DCF%s.mat',obj.twix.hdr.Config.MeasUID,obj.flags.doB0Corr,obj.flags.doDCF);
            sp=obj.SpiralPara;
            fn=obj.filename;
            ro=(2*sp.ADCLength*sp.DwellTime)/1e6; % ms
            vTR=(sp.TR*sp.Ninterleaves*sp.NPartitions)/(sp.R_PE*sp.R_3D*1e6); %s
            descrip=(sprintf('R%dx%dC%d TR=%.1fms RO=%.2fms vTR=%.1fs',sp.R_PE,sp.R_3D,sp.CAIPIShift,sp.TR/1e3,ro,vTR));
            descrip_reco=sprintf('%s PAT=%s coilcomb=%s B0=%s DCF=%s CompMode=%s',flags.CompMode,flags.doPAT, flags.doCoilCombine, flags.doB0Corr,flags.doDCF,flags.CompMode);
            save(OutFile,'im','sp','flags','descrip','descrip_reco','fn','-v7.3')
            
            
        end
        
    end
end

%% old crap
%         function obj=performB0corr(obj)
%             if(~strcmpi(obj.flags.doB0Corr,'none'))
%             if( isempty(obj.B0Map) || ~isa(obj.B0Map,'B0map'))
%                 [path,~]=fileparts(obj.filename);
%                 obj.B0Map=B0map(path);
%             end
%             if(obj.LoopCounter.cRep==1)
% %                 obj.B0Map=obj.B0Map.PerformSliceSelection(obj.soda_obj.Coords{obj.LoopCounter.cSlc},squeeze(abs(obj.img(1,:,:,:,obj.LoopCounter.cSlc,1,1))));
%                 if(obj.flags.is3D)
%                     obj.B0OPerator=MTI_3D(((2*pi).^2)*obj.B0Map.Fmap_registered, obj.time(:)*1e-6,obj.NUFFT_obj,obj.coilSens(:,:,:,:,obj.LoopCounter.cSlc));
%                 else
%                     obj.B0OPerator=MTI(((2*pi).^2)*obj.B0Map.Fmap_registered, obj.time(:)*1e-6,obj.NUFFT_obj,obj.coilSens(:,:,:,1,obj.LoopCounter.cSlc));
%                 end
%                 
%             end
%              if(obj.flags.is3D)
%                   obj.img(2,:,:,:,obj.LoopCounter.cSlc,obj.LoopCounter.cRep)=obj.B0OPerator'*double(obj.sig);
%              else
%                 obj.img(2,:,:,1,obj.LoopCounter.cSlc,obj.LoopCounter.cRep) = obj.B0OPerator'*double(obj.sig);
%              end
%             end
%         end


