classdef kspace
    properties
        Type;
        Gradient % G/cm
        G_dist % G/cm
        k%kspace trajectory
        k_dist % distorted traj
        parameter=struct('Ninterleaves', 32,'Resolution',0.3,'FOV',[128 128 128 128],'MaxGradAmp',4,'MaxSlewrate',20,'DwellTime',5e-6)%fov,resolution
        ZeroPadLength
        %units FOV and Resolution are both in cm, maxGrad: G/cm,
        %maxslew is G/cm/s
        wi % density compensation
        PSF
        
    end
    methods
        function obj= kspace(varargin) %  constructor
            if(nargin==0)
                obj.Type='SpiralOut';
                obj.ZeroPadLength=512;
                obj=obj.getk();
                w=jacksonDCF2(obj.k,obj.parameter);
                obj.wi=w./max(w(:));
                obj=obj.calcPSF();
                obj=distortk(obj);
            elseif(nargin==1)
                obj.parameter=varargin{1};
                obj.Type='SpiralOut';
                obj.ZeroPadLength=512;
                obj=obj.getk();
                w=jacksonDCF2(obj.k,obj.parameter);
                obj.wi=w./max(w(:));
                obj=obj.calcPSF();
                obj=distortk(obj);
            end
        end
        
        function obj=getk(obj)
            p=obj.parameter;
            [~,G] = vdSpiralDesign(p.Ninterleaves, 4, p.Resolution*10,p.FOV,[0,0.2,0.3,1],p.MaxGradAmp,p.MaxSlewrate,p.DwellTime*1e3,[],'cubic');
            %when resolution is 1cm, then fov becomes L cm
            gammaH=4.2576e3; %Hz/G
            
            G=complex(G(:,1),G(:,2)); %G/cm
            obj.Gradient=zeros(length(G)+2,p.Ninterleaves);
            for intlv=1:p.Ninterleaves
                obj.Gradient(:,intlv)=[0; G*exp(-1i*2*pi*(intlv/p.Ninterleaves)); 0;];
            end
            
            %calculate kspace
            obj.k= cumsum(2*pi*gammaH*obj.Gradient*p.DwellTime,1); % rad/cm
            obj.k=obj.k./(2*pi)*(p.Resolution); %scale to  kmax=0.5
            
        end %getk
        function obj=calcPSF(obj)
            %Loads meaured triangular pulses to  
            win_size=length(obj.Gradient)+obj.ZeroPadLength*2;
            load('GIRF_0612_gammaH.mat','cal_sig','rec_sig');
            in_sig=fft(cal_sig(:,:,:,:,:),win_size,1);
            out_sig=fft(rec_sig(:,:,:,:,:),win_size,1);
            num=squeeze(sum(sum(sum(conj(in_sig).*out_sig,2),4),5));
            denom=squeeze(sum(sum(sum((abs(in_sig).^2),2),4),5));
            obj.PSF=num./denom;
            obj.PSF(isnan(obj.PSF) | isinf(obj.PSF))=0;
        end
        
        function obj=distortk(obj)
            win_size=length(obj.Gradient)+obj.ZeroPadLength*2;
            gammaH=4.2576e3; %Hz/G
            
            
            if (obj.ZeroPadLength>0)
                gx=[zeros(obj.ZeroPadLength,obj.parameter.Ninterleaves); real(obj.Gradient); zeros(obj.ZeroPadLength,obj.parameter.Ninterleaves)];
                gy=[zeros(obj.ZeroPadLength,obj.parameter.Ninterleaves); imag(obj.Gradient); zeros(obj.ZeroPadLength,obj.parameter.Ninterleaves)];
            end
            gxf=fft(gx,win_size,1);
            gyf=fft(gy,win_size,1);
            
            gdist=complex(real(ifft(gxf.*obj.PSF(:,1))),real(ifft(gyf.*obj.PSF(:,2))));
            
            obj.G_dist=gdist(obj.ZeroPadLength+1:end-obj.ZeroPadLength,:);
            %calculate kspace
            obj.k_dist= cumsum(2*pi*gammaH*obj.G_dist*obj.parameter.DwellTime,1); % rad/cm
            obj.k_dist=obj.k_dist./(2*pi)*(obj.parameter.Resolution); %scale to  kmax=0.5
            
        end
        function obj=fixk(obj) 
            % Dirty function: don't use
            %replace the first and last few points of distorted Gradient
            %from original gradient waveform
            
            gammaH=4.2576e3; %Hz/G
            %             timeshift=0; %apparently after applying the PSF, gradients are shifted by 3 points
%             obj.G_dist(1:end-timeshift+1,:)=obj.G_dist(timeshift:end,:);
%             obj.G_dist(end-timeshift:end,:)=0;
            
            %replace first 10 and last ten values with orginal G
            obj.G_dist(1:2,:)=obj.Gradient(1:2,:);
            obj.G_dist(end-2:end,:)=obj.Gradient(end-2:end,:);
            
            %calculate kspace
            obj.k_dist= cumsum(2*pi*gammaH*obj.G_dist*obj.parameter.DwellTime,1); % rad/cm
            obj.k_dist=obj.k_dist./(2*pi)*(obj.parameter.Resolution); %scale to  kmax=0.5
        end
        
        function obj = voronoidens(obj)
            %
            % function area = voronoidens(k);
            %
            % input:  k = kx + i ky is the  k-space trajectory
            % output: area of cells for each point
            %           (if point doesn't have neighbors the area is NaN
           
            kxy = [real(obj.k(:)),imag(obj.k(:))];
            % returns vertices and cells of voronoi diagram
            [V,C] = voronoin(kxy);
            area = zeros(size(C));
            for j = 1:length(kxy)
                
                try
                    x = V(C{j},1);
                    y = V(C{j},2);
                    lxy = length(x);
                    A = abs(sum( 0.5*(x([2:lxy 1]) - x(:)).*(y([2:lxy 1]) + y(:))));
                catch
                    warning('substituting area=0 for index %d',j);
                    A=NaN;
                end
                area(j) =A;
            end
            area=reshape(area,size(obj.k));
            
            %remove last n% samples with samples:large area near the
            %boundary
            clip_point=round(0.9*size(obj.k,1));
            area(clip_point:end,:)=repmat(area(clip_point,:),[length(clip_point:size(obj.k,1)) 1]);
            
            obj.wi =area./max(area(:)) ;
        end
        
        function plotGradient(obj,intlv)
            if(nargin<2)
                intlv=1;
            end
            t=linspace(0,(size(obj.Gradient,1)-1)*10,size(obj.Gradient,1)); %us
            %nominal gradients
            gx=real(obj.Gradient(:,intlv))*10; % G/cm to mT/m
            gy=imag(obj.Gradient(:,intlv))*10; % G/cm to mT/m
            %distorted with PSF
            gxr=real(obj.G_dist(:,intlv))*10; % G/cm to mT/m
            gyr=imag(obj.G_dist(:,intlv))*10; % G/cm to mT/m
            
            figure,
            plot(t,gxr,'color',	'#EDB120'),hold on,plot(t,gx,'color', 	'#D95319')
            plot(t,gyr, 'LineStyle','--' ,   'color',	'#EDB120'),hold on,plot(t,gy,'LineStyle','--' ,  'color', 	'#D95319')
            
            legend('Reconstructed','Nominal'),xlabel('Time(us)'),ylabel('Grad Amplitude(mT/m)'),title('Gradients: GIRF results')
            figure,
            plot(t,(gx-gxr),'color',	'#EDB120'),hold on,plot(t,ones(size(gx)).*0.1),plot(t,ones(size(gx)).*-0.1),plot(t,zeros(size(t)))
            plot(t,(gy-gyr), 'LineStyle','--' ,   'color',	'#EDB120')
            legend('nominal-distorted','0.1mT/m','-0.1mT/m'),xlabel('time(us)'),ylabel('Deviations(mT/m)')
        end
        function plotKTraj(obj,intlv)
            %plot the kspace trajectories of the interleave
            if(nargin<2)
                intlv=1;
            end
            t=linspace(0,(size(obj.Gradient,1)-1)*10,size(obj.Gradient,1)); %us
            %nominal gradients
            kx=real(obj.k(:,intlv));
            ky=imag(obj.k(:,intlv));
            %distorted with PSF
            kxr=real(obj.k_dist(:,intlv));
            kyr=imag(obj.k_dist(:,intlv));
            
            figure,
            plot(t,kxr,'color',	'#EDB120'),hold on,plot(t,kx,'color', 	'#D95319')
            plot(t,kyr, 'LineStyle','--' ,   'color',	'#EDB120'),hold on,plot(t,ky,'LineStyle','--' ,  'color', 	'#D95319')
            
            legend('Reconstructed','Nominal'),xlabel('Time(us)'),ylabel('Amplitude(Normalised)'),title('Kspace Trajectories: GIRF results')
%             figure,
%             plot(t,(kx-kxr),'color',	'#EDB120'),hold on,
%             plot(t,(ky-kyr), 'LineStyle','--' ,   'color',	'#EDB120')
%             legend('nominal-distorted'),xlabel('time(us)'),ylabel('Deviations(mT/m)')
        end
    end
end
