classdef MriSim<matlab.mixin.Copyable
    properties
        Type='Analytical';
        k%kspace trajectory object
        t % time in s
        b0_map % in rad/s
        parameter
        
        data % kspace samples
        im_true
        im_recon
        CoilSens % [PxRx1xCoils]
        DCF % density compensation
        RecoSettings=struct('CompMode','CPU2DHybrid','doB0corr','MTI','Oversampling',1.5,'KernelSize',5) %gridding reconstruction
        NUFFT_obj
    end
    methods
        function obj=MriSim(varargin) %  constructor
              
                
                %units FOV and Resolution are both in mm, maxGrad: mT/m,
                obj.parameter= struct('Ninterleaves', 4,'Resolution',1,'FOV',[1 1 1 1]*192,...
            'MaxGradAmp',42,'MinRiseTime',5,'B0effects',1,'T2star',20e-3);

                obj.parameter.nCoils=4;
                 obj.parameter.AccelFactor=1;
                obj.parameter.MatSize=obj.parameter.FOV(1)/obj.parameter.Resolution;
                
                obj.parameter.NoiseFactor=1e-3;
                obj.parameter.TE=2000; %[us]
                obj.parameter.DwellTime=2500; %[ns]

                obj.parameter.GradDelay=[0 0 0];
                obj.parameter.GRAD_RASTER_TIME_DEFAULT=10; %[us]
                obj.parameter.gammaH=42.567e6; %[Hz/T]
                obj.parameter.SpiralType=1;
                
                SpiralType= {'SpiralOut','SpiralIn', 'DoubleSpiral', 'SpiralInAndOut' }; %selSpiraldir
                obj.parameter.SpiralTypeName=SpiralType{obj.parameter.SpiralType};

            if(mod(nargin,2)==0)
                % allow users to chnage the parameters
                fn=fieldnames(obj.parameter);
                fn2=fieldnames(obj.RecoSettings);
                for i=1:2:nargin
                    idx1=strcmpi(varargin(i),fn);
                    if(sum(idx1)==1); obj.parameter.(fn{idx1})=varargin{i+1};end

                    idx2=strcmpi(varargin(i),fn2);
                    if(sum(idx2)==1); obj.RecoSettings.(fn2{idx2})=varargin{i+1};end

                    if(~any([idx1; idx2])); warning('Unknown parameter : %s ',varargin{i}); end
                end
                obj.parameter.Resolution=obj.parameter.FOV(1)/obj.parameter.MatSize;
               
            end
            %functions
                obj.InitializePhantom();
                obj.getKtraj();
                obj.getCoilSens();
                obj.b0_map=obj.setB0map(obj.parameter.B0effects);
                obj.data=obj.getsig();
                obj.performUndersampling();
                obj.PerformRecon();
        end

        function obj=recalculate(obj)
                obj.im_true=obj.InitializePhantom();
                obj=obj.getKtraj();
                obj=obj.getCoilSens();
                obj.b0_map=obj.setB0map(obj.parameter.B0effects);
       
                 obj.b0_map=[];
       
                obj.data=obj.getsig();
                obj.performUndersampling();
                obj=obj.PerformRecon();
        end
        function showImages(obj)
            mask=obj.im_true>0;
            figure,subplot(121),imagesc(obj.im_true,[0 1.8]),colorbar,title('True image')
            subplot(122),imagesc(real(obj.im_recon)-min(real(obj.im_recon))),colorbar,title('Simulated')
        end

        function phantom=InitializePhantom(obj)
            if(strcmpi(obj.Type,'Analytical')) % we need analytical dft as well
                [x ,y]=ndgrid(linspace(-1,1,obj.parameter.MatSize));
                maskCircle=@(x,y,center,radius) ((x-center(1)).^2+(y-center(2)).^2 < radius.^2);
                maskRect=@(x,y,center,length,breadth) (abs(x-center(1))<=length/2 & abs(y-center(2)) <=breadth/2);
                % centerx,centery, length, breadth,amplitude
                rect_dim=[0 ,0.05,0.85*2 ,0.05,0.5;...
                    0 ,-0.05,0.85*2 ,0.05,0.5;...
                    0.1000, -0.4, 0.1000, 0.4,0.8;...
                    0.2075, -0.4, 0.0750, 0.4,.80;...
                    0.2900, -0.4, 0.0500, 0.4,.80;...
                    0.3475, -0.4, 0.0250, 0.4,.80;...
                    0.3862, -0.4, 0.0125, 0.4,.80;...
                    0.4156, -0.4, 0.0063, 0.4,.80;...
                    -0.4, -0.5500, 0.4, 0.1000,.80;...
                    -0.4, -0.4425, 0.4, 0.0750,.80;...
                    -0.4, -0.3600, 0.4, 0.0500,.80;...
                    -0.4, -0.3025, 0.4, 0.0250,.80;...
                    -0.4, -0.2638, 0.4, 0.0125,.80;...
                    -0.4, -0.2344, 0.4, 0.0063,.80;...
                    ];
                %centerx,centery, radius,amplitude
                circ_dim=[0,0,0.9,1;...
                    -0.5 ,0.4,0.15,0.7;...
                    0 ,0.4,0.15,0.7;...
                    0.5 ,0.4,0.15,0.7;...
                    ];
                % circ_amp=[1 0.7 0.7 0.7];
                % rect_amp=[0.5000    0.5000    0.8000    0.8000    0.8000    0.8000    0.8000    0.8000    0.8000    0.8000    0.8000    0.8000    0.8000    0.8];
                phantom=zeros(obj.parameter.MatSize);
                for i=1:size(circ_dim,1)
                    p=zeros(obj.parameter.MatSize);
                    p(maskCircle(x,y,circ_dim(i,1:2),circ_dim(i,3)))= circ_dim(i,4);
                    phantom=phantom+p;
                end
                for i=1:size(rect_dim,1)
                    p=zeros(obj.parameter.MatSize);
                    p(maskRect(x,y,rect_dim(i,1:2),rect_dim(i,3),rect_dim(i,4)))= rect_dim(i,5);
                    phantom=phantom+p;
                end
            end
             obj.im_true=phantom;
        end % Initialze phantom
        function obj=getKtraj(obj)

            p=obj.parameter;
            [G] = getSpiralTraj(p.Ninterleaves, p.Resolution,p.FOV,[0,0.2,0.3,1],p.MaxGradAmp,p.MinRiseTime,p.SpiralType);
            G_int=zeros(size(G,1), size(G,2) ,p.Ninterleaves);
            phase=-2*(0:(p.Ninterleaves-1))*pi/p.Ninterleaves;

            for m=1:p.Ninterleaves
                RotMat=[cos(phase(m)) -sin(phase(m)) 0; sin(phase(m)) cos(phase(m)) 0; 0 0 1 ];
                G_int(:,:,m)=RotMat*G;
            end

            %add some more parameters for kspace calculation
            obj.parameter.ADCLength=round((size(G_int,2)*10)/(2*p.DwellTime*1e-3));

            %calculate kspace
            [obj.k,obj.t]=Grad2TrajHigherorder(permute(G_int,[2 1 3]),obj.parameter);
            obj.k=permute(obj.k,[2 1 3]);
            kmax=2*pi*(0.5/(obj.parameter.Resolution*1e-3));
            obj.k=obj.k/(2*kmax);
            obj.t=obj.t.*1e-6; %[us -> s]
            % calculate the sampling density compenation function
            obj.DCF=jacksonDCF2(squeeze(complex(obj.k(1,:,:),obj.k(2,:,:))),obj.parameter);
%             obj.DCF=VoronoiDCF2D(squeeze(complex(obj.k(1,:,:),obj.k(2,:,:))));
        end

        function mrsig=getsig(obj)
%             warning('Inverse Crime!!!!!!!!!!!! Need analytical formula')
            ms=obj.parameter.MatSize;
            [rx,ry]=ndgrid(((0:(ms-1))-floor(ms/2)).*2*pi);
            k1=obj.k;
            t1=obj.t(:);
            mrsig=zeros([size(k1,2) size(k1,3) 1 size(obj.CoilSens,4)]);
            dxdy=(obj.parameter.FOV(1)/obj.parameter.MatSize).^2; %pixel area in m^2       
                for tp=1:size(k1,2)
                     if(obj.parameter.B0effects>0)
                         b0mod=exp(-1i*obj.b0_map*t1(tp));
                     else
                         b0mod=1;
                     end
                    T2starTerm= exp(-t1(tp)/obj.parameter.T2star);
                    for intlv=1:size(k1,3)
                        basis=exp(-1i*(rx*k1(1,tp,intlv)+ry*k1(2,tp,intlv)));
                        for coil=1:size(obj.CoilSens,4)
                            mrsig(tp,intlv,1,coil)=sum((obj.CoilSens(:,:,1,coil)).*obj.im_true.*basis.*T2starTerm.*b0mod*dxdy,'all');
                        end
                    end
                end
                mrsig=mrsig+obj.parameter.NoiseFactor*complex(randn(size(mrsig)),randn(size(mrsig)));
        end%getsig()

        function b0map_true=setB0map(obj,Mapselect)
            %generate some field map in Hz
            [X,Y]=meshgrid(linspace(-1,1,obj.parameter.MatSize));
            switch Mapselect
                case 0
                    b0map_true=[];
                case 1 % smooth maps
                    b0map_true=(50/(2*pi*0.1))*exp(-0.5* (X.^2 + Y.^2)./(2*0.05))- ...
                        (200/(2*pi*0.2))*exp(-0.5* ((X+0.2).^2 + (Y+0.1).^2)./(2*0.1))-...
                        (600/(2*pi*0.2))*exp(-0.5* ((X-0.2).^2 + (Y-0.1).^2)./(2*0.1));
                case 2 % Discontinous maps
                    b0map_true= -60+((0.01*(X-0.25).^2 + 0.04*(Y-0.25).^2)<(0.01*0.04))*120+...
                        ((0.01*(X+0.25).^2 + 0.04*(Y+0.25).^2)<(0.01*0.04))*120;
            end


        end

        function E=getEncodingmatrix(obj)
            
            [rx,ry]=ndgrid( ((0:(obj.parameter.MatSize-1)) -floor(obj.parameter.MatSize/2)).*2*pi );
            k1=obj.k;
            t1=obj.t(:);
            warning('Need %.2f GB memory',prod([size(rx) size(k1,2)*size(k1,3)*size(obj.CoilSens,4)])*16/(2^30));
            if((prod([size(rx) size(k1,2)*size(k1,3)*size(obj.CoilSens,4)])*16/(2^30)) >16 ); error('Need more than 16 GB'); end
              idx=kron(eye(2)>0,ones([numel(rx) size(k1,2)*size(k1,3)],'logical'))>0;

            E=zeros([numel(rx) size(k1,2) size(k1,3) size(obj.CoilSens,4)]);
            dxdy=(obj.parameter.FOV(1)/obj.parameter.MatSize).^2; %pixel area in m^2
            %             dxdy=(obj.parameter.Resolution*1e-2).^2; %pixel area in m^2

                for tp=1:size(k1,2)
                     if(obj.parameter.B0effects>0)
                         b0mod=exp(-1i*obj.b0_map*t1(tp));
                     else
                         b0mod=1;
                     end
                    for intlv=1:size(k1,3)
                        basis=exp(-1i*(rx*k1(1,tp,intlv)+ry*k1(2,tp,intlv)));
                        for coil=1:size(obj.CoilSens,4)
                            tmp=(obj.CoilSens(:,:,1,coil)).*basis.*b0mod*dxdy;
                            E(:,tp,intlv,coil)=tmp(:);
                        end
                    end
                end
                E=reshape(E,[numel(rx) size(k1,2)*size(k1,3)*size(obj.CoilSens,4)]).';

%                 % mutli coil bad memory hog
%                 idx=find(kron(eye(size(obj.CoilSens,4),'single'),ones([numel(rx) size(k1,2)*size(k1,3)],'single'))>0);
%                 [row,col] = ind2sub([numel(rx)*size(obj.CoilSens,4) size(k1,2)*size(k1,3)*size(obj.CoilSens,4) ],(idx));
%                 E2=sparse(col,row,E(:));
%                 %dense
% %                 E1=zeros([numel(rx)*size(obj.CoilSens,4) size(k1,2)*size(k1,3)*size(obj.CoilSens,4) ]);
% %                  E1(idx)=E(:);
                
              
        end

        function obj=getCoilSens(obj)
            obj.CoilSens=ones(obj.parameter.MatSize,obj.parameter.MatSize,1,obj.parameter.nCoils);
            if(obj.parameter.nCoils>1)
                sigma=eye(2)*0.5;
                nCoils= obj.parameter.nCoils;
                for nc=1:nCoils
                    mu=[real(exp(1i*(2*pi*nc)/nCoils)) imag(exp(-1i*(2*pi*nc)/nCoils)) ];
                    % Define the grid of points
                    x = linspace(-0.5, 0.5, obj.parameter.MatSize);
                    y = linspace(-0.5, 0.5, obj.parameter.MatSize);
                    [X, Y] = meshgrid(x, y);


                    % Calculate the Gaussian values at each point
                    Z = 1 / (2 * pi * det(sigma)^0.5) * ...
                        exp(-0.5 * ((X - mu(1)).^2 / sigma(1,1) + ...
                        (Y - mu(2)).^2 / sigma(2,2)));

                    obj.CoilSens(:,:,1,nc)=Z.*exp(1i*2*pi*(X.^2+Y.^2));

                end

                obj.CoilSens=obj.CoilSens./sum(abs(obj.CoilSens),4);

            end

        end
        function obj=performUndersampling(obj)
            
            R=obj.parameter.AccelFactor;
            if(R>1 && size(obj.data,2)==obj.parameter.Ninterleaves)
            obj.DCF=obj.DCF(:,1:R:end,:);
            obj.k=obj.k(:,:,1:R:end,:);
            obj.data=obj.data(:,1:R:end,:,:);
            end

        end
        function obj=PerformRecon(obj)
    
            if(size(obj.CoilSens,4)==1)
                csm=[];
            else
                csm=obj.CoilSens;
            end




            switch(obj.RecoSettings.doB0corr)
                case {'MTI','MFI'}
            obj.NUFFT_obj=StackofSpiralsB0(obj.k,obj.DCF,[obj.parameter.MatSize, obj.parameter.MatSize,1],csm,obj.b0_map,obj.t,'CompMode',obj.RecoSettings.CompMode,'Method',obj.RecoSettings.doB0corr);
                otherwise
            obj.NUFFT_obj=StackofSpirals(obj.k,obj.DCF,[obj.parameter.MatSize, obj.parameter.MatSize,1],csm,'CompMode',obj.RecoSettings.CompMode);
            
            end
            sig=obj.data.*sqrt(obj.DCF);
            obj.im_recon=obj.NUFFT_obj'*sig;
            %             imr_oversamp=gridkb(obj.data,obj.k.k,obj.k.wi,obj.parameter.FOV(1),obj.RecoSettings.Oversampling,obj.RecoSettings.KernelSize,'image');
            %             imr=image_removeOS(imr_oversamp,obj.RecoSettings.Oversampling);
            %             imr=0;
        end

        function ShowReco(obj)
            figure
            subplot(1,2,1)
            imagesc(abs(obj.im_recon)),
            title('Magnitude image'),colorbar,axis square
            
            subplot(1,2,2)
            imagesc((angle(obj.im_recon))),
            title('Phase image'),colorbar,axis square
        end


    end
end
