classdef MriSim
    properties
        Type='Analytical';
        k%kspace trajectory object
        t % time in s
        b0_map % in rad/s
        %k_dist % distorted traj
        parameter=struct('Ninterleaves', 4,'Resolution',1,'FOV',[128 128 128 128],...
            'MaxGradAmp',4,'MaxSlewrate',20,'DwellTime',5e-6,'B0effects',1,...
            'T2star',10e-3);%fov,resolution
        %units FOV and Resolution are both in cm, maxGrad: G/cm,
        %maxslew is G/cm/s
        
        MatSize% phantom Resolution
        data % kspace samples
        im_true
        im_recon
        %wi % density compensation
        GridSetting=struct('Oversampling',1.5,'KernelSize',5) %gridding reconstruction
    end
    methods
        function obj= MriSim(varargin) %  constructor
            if(nargin==0)
                obj.MatSize=128;
                obj.im_true=obj.InitializePhantom();
                obj.k= kspace(obj.parameter);
                obj.t=repmat((0:(size(obj.k.k,1)-1))*obj.parameter.DwellTime,[ size(obj.k.k,2) 1]).';
                if(obj.parameter.B0effects>0)
                     obj.b0_map=obj.setB0map(obj.parameter.B0effects);
                else
                    obj.b0_map=[];
                end
%                 w=jacksonDCF2(obj.k,obj.parameter);
                %obj.wi=w./max(w(:));
                obj.data=obj.getsig();
                obj.im_recon=obj.Griddingrecon();
            elseif(nargin==1)
                parameter=varargin{1};
                obj.MatSize=ceil(parameter.FOV(1)./parameter.Resolution);
                obj.parameter=parameter;
                obj.im_true=obj.InitializePhantom();
                obj.k=kspace(obj.parameter);
                obj.t=repmat((0:(size(obj.k.k,1)-1))*obj.parameter.DwellTime,[ size(obj.k.k,2) 1]).';
                if(obj.parameter.B0effects>0)
                     obj.b0_map=obj.setB0map(obj.parameter.B0effects);
                else
                    obj.b0_map=zeros(size(obj.im_true));
                end
                obj.data=obj.getsig();
                obj.im_recon=obj.Griddingrecon();
                obj.t=repmat((0:(size(obj.k.k,1)-1))*obj.parameter.DwellTime,[ size(obj.k.k,2) 1]).';
           
            elseif(nargin==2)
                parameter=varargin{1};
                GridSetting=varargin{2}; % Default constructor
                obj.MatSize=parameter.FOV(1);
                %             obj.parameter=parameter;
                obj.parameter=parameter;
                obj.GridSetting=GridSetting;
                obj.im_true=obj.InitializePhantom();
                obj.k=kspace(obj.parameter);
                obj.t=repmat((0:(size(obj.k.k,1)-1))*obj.parameter.DwellTime,[ size(obj.k.k,2) 1]).';
                if(obj.parameter.B0effects>0)
                     obj.b0_map=obj.setB0map(obj.parameter.B0effects);
                else
                    obj.b0_map=zeros(size(obj.im_true));
                end
                obj.data=obj.getsig();
                obj.im_recon=obj.Griddingrecon();
                obj.t=repmat((0:(size(obj.k.k,1)-1))*obj.parameter.DwellTime,[ size(obj.k.k,2) 1]).';
           
                
            end
        end
        
        function obj=recalculate(obj)
                obj.im_true=obj.InitializePhantom();
                obj.k=kspace(obj.parameter);
                obj.data=obj.getsig();
                obj.im_recon=obj.Griddingrecon();
            
        end
        function showImages(obj)
            mask=obj.im_true>0;
            figure,subplot(121),imagesc(obj.im_true,[0 1.8]),colorbar,title('True image')
            subplot(122),imagesc(real(obj.im_recon)-min(real(obj.im_recon))),colorbar,title('Simulated')
        end
        
        function phantom=InitializePhantom(this)
            if(strcmpi(this.Type,'Analytical'))
                [x ,y]=ndgrid(linspace(-1,1,this.MatSize));
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
                phantom=zeros(this.MatSize);
                for i=1:size(circ_dim,1)
                    p=zeros(this.MatSize);
                    p(maskCircle(x,y,circ_dim(i,1:2),circ_dim(i,3)))= circ_dim(i,4);
                    phantom=phantom+p;
                end
                for i=1:size(rect_dim,1)
                    p=zeros(this.MatSize);
                    p(maskRect(x,y,rect_dim(i,1:2),rect_dim(i,3),rect_dim(i,4)))= rect_dim(i,5);
                    phantom=phantom+p;
                end
            end
        end % Initialze phantom

        function mrsig=getsig(this)
            [rx,ry]=ndgrid( ((0:(this.MatSize-1)) -floor(this.MatSize/2)).*2*pi );
            k1=this.k.k(:);
            t1=this.t(:);
            mrsig=zeros(size(k1(:)));
            dxdy=(this.parameter.FOV(1)/this.MatSize).^2; %pixel area in m^2
            %             dxdy=(this.parameter.Resolution*1e-2).^2; %pixel area in m^2
            
            if(this.parameter.B0effects>0)
                for i=1:length(k1)
                    basis=exp(-1i*(rx*real(k1(i))+ry*imag(k1(i))));
                    b0mod=exp(flip(flip(this.b0_map,1),2)*-2i*pi*t1(i));
                    T2starTerm= exp(-t1(i)/this.parameter.T2star);
%                     mask=(this.im_true>0);
                    mrsig(i)=sum(this.im_true.*basis.*T2starTerm.*b0mod*dxdy,'all');
                end
            else
                for i=1:length(k1)
                    T2starTerm= exp(-t1(i)/this.parameter.T2star);
                    basis=exp(-1i*(rx*real(k1(i))+ry*imag(k1(i))));
%                     mask=(this.im_true>0);
                    mrsig(i)=sum(this.im_true.*basis*dxdy.*T2starTerm,'all');
                end
            end
            %filter
%             mrsig_r=(myfilterPhaseData((reshape(real(mrsig),[],this.parameter.Ninterleaves)),this.parameter.DwellTime,2*this.parameter.DwellTime));
%             mrsig_i=(myfilterPhaseData((reshape(imag(mrsig),[],this.parameter.Ninterleaves)),this.parameter.DwellTime,2*this.parameter.DwellTime));
%             mrsig=complex(mrsig_r(:),mrsig_i(:));
        end%getsig()
        
        function b0map_true=setB0map(this,Mapselect)
            %generate some field map in Hz
            [X,Y]=meshgrid(linspace(-1,1,this.MatSize));
            switch Mapselect
                case 1 % smooth maps
            b0map_true=(200/(2*pi*0.1))*exp(-0.5* (X.^2 + Y.^2)./(2*0.05))- ...
                (200/(2*pi*0.2))*exp(-0.5* ((X+0.2).^2 + (Y+0.1).^2)./(2*0.1))-...
                (300/(2*pi*0.2))*exp(-0.5* ((X-0.2).^2 + (Y-0.1).^2)./(2*0.1));
                case 2 % Discontinous maps
                    b0map_true= -60+((0.01*(X-0.25).^2 + 0.04*(Y-0.25).^2)<(0.01*0.04))*120+...
                                    ((0.01*(X+0.25).^2 + 0.04*(Y+0.25).^2)<(0.01*0.04))*120;
            end
            
            
        end
        
        function E=getEncodingmatrix(this,os)
            if(nargin<2)
            os=1; %oversampling factor
            end
            [rx,ry]=ndgrid( ((0:(os*this.MatSize-1)) -floor(os*this.MatSize/2)).*2*pi );
            k1=this.k.k(:);
            t1=this.t(:);
%             mrsig=zeros(size(k1(:)));
            dxdy=(this.parameter.FOV(1)/(os*this.MatSize)).^2; %pixel area in m^2
            %             dxdy=(this.parameter.Resolution*1e-2).^2; %pixel area in m^2
            E=zeros(length(k1),numel(rx));
            if(this.parameter.B0effects>0)
                for i=1:length(k1)
                    basis=exp(-1i*(rx*real(k1(i))+ry*imag(k1(i))));
                    b0mod=exp(flip(flip(this.b0_map,1),2)*-2i*pi*t1(i));
                    T2starTerm= exp(-t1(i)/this.parameter.T2star);
                    E(i,:)=basis(:).*T2starTerm.*b0mod(:)*dxdy;
%                     mask=(this.im_true>0);
%                     mrsig(i)=sum(this.im_true.*basis.*T2starTerm.*b0mod*dxdy,'all');
                end
            else
                for i=1:length(k1)
                    T2starTerm= exp(-t1(i)/this.parameter.T2star);
                    basis=exp(-1i*(rx*real(k1(i))+ry*imag(k1(i))));
%                     mask=(this.im_true>0);
                     E(i,:)=basis(:).*T2starTerm.*dxdy;
%                     mrsig(i)=sum(this.im_true.*basis*dxdy.*T2starTerm,'all');
                end
            end
        end
        function imr=Griddingrecon(this)
            imr_oversamp=gridkb(this.data,this.k.k,this.k.wi,this.parameter.FOV(1),this.GridSetting.Oversampling,this.GridSetting.KernelSize,'image');
            imr=image_removeOS(imr_oversamp,this.GridSetting.Oversampling);
        end

    end
end


%% old function moved to kspace objects
%         function obj=getk(obj)
%             p=obj.parameter;
%             [~,G] = vdSpiralDesign(p.Ninterleaves, 4, p.Resolution*10,p.FOV,[0,0.2,0.3,1],p.MaxGradAmp,p.MaxSlewrate,p.DwellTime*1e3,[],'cubic');
%             %when resolution is 1cm, then fov becomes L cm
%             gammaH=4.2576e3; %Hz/G
%             
%             G=complex(G(:,1),G(:,2)); %G/cm
%             Gr=zeros(length(G)+2,p.Ninterleaves);
%             for intlv=1:p.Ninterleaves
%                 Gr(:,intlv)=[0; G*exp(-1i*2*pi*(intlv/p.Ninterleaves)); 0;];
%             end
%             
%             %calculate kspace
%             k1= cumsum(2*pi*gammaH*Gr*p.DwellTime,1); % rad/cm
%             obj.k.k=k1./(2*pi)*(p.Resolution); %scale to  kmax=0.5    
%         end %getk
%         function k1=distortk(this)
%                
%             p=this.parameter;
%             [~,G] = vdSpiralDesign(p.Ninterleaves, 4, p.Resolution*10,p.FOV,[0,0.2,0.3,1],p.MaxGradAmp,p.MaxSlewrate,p.DwellTime*1e3,[],'cubic');
%             %when resolution is 1cm, then fov becomes L cm
%             gammaH=4.2576e3; %Hz/G
%             
%       
%             
%             %PSF
%             win_size=length(G)+2;
%             load('GIRF_0612_gammaH.mat','cal_sig','rec_sig');
%             in_sig=fft(cal_sig(:,:,:,:,:),win_size,1);
%             out_sig=fft(rec_sig(:,:,:,:,:),win_size,1);
%             num=squeeze(sum(sum(sum(conj(in_sig).*out_sig,2),4),5));
%             denom=squeeze(sum(sum(sum((abs(in_sig).^2),2),4),5));
%             PSF=num./denom;
%             PSF(isnan(PSF))=0;
%             PSF(isinf(PSF))=0;
%             
%             G=complex(G(:,1),G(:,2)); %G/cm
%             Gr=zeros(length(G)+2,p.Ninterleaves);
%             G_distorted=zeros(length(G)+2,p.Ninterleaves);
%             for intlv=1:p.Ninterleaves
%                 Gr(:,intlv)=[0; G*exp(-1i*2*pi*(intlv/p.Ninterleaves)); 0;];
%                 
%                 
%             end
% 
%             %calculate kspace
%             k1= cumsum(2*pi*gammaH*Gr*p.DwellTime,1); % rad/cm
%             k1=k1./(2*pi)*(p.Resolution); %scale to  kmax=0.5
%             
%         end