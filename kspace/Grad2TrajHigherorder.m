function [k2,t]=Grad2TrajHigherorder(G_SPH,parameter)
% [k2,t]=Grad2TrajHigherorder(G_SPH,parameter)
%Function to get kspace trajectories from the gradient waveforms
%gradient units is mT/m sampled at 10 us(typ)
%3D matrix: Timepoints x Gradient Axis x interleaves
% ONLY FOR SPIRALOUT at the moment

%Get all parameters
Gyro=42.575575e3;%Hz/mT %parameter.GyroH;

if(contains(parameter.SpiralTypeName,'SpiralInAndOut'))
    parameter.ADCLength=parameter.ADCLength*4; % read over sampling and 2 ADC
else
    parameter.ADCLength=parameter.ADCLength*2; % read over sampling
end

    G_SPH=padarray(G_SPH,[1 0 0 0],0,'Pre');
%     G_SPH=0.5*(G_SPH(1:end-1,:,:,:)+G_SPH(2:end,:,:,:));

    k_cal=(parameter.gammaH*2*pi)*(parameter.GRAD_RASTER_TIME_DEFAULT*1e-6)*(cumsum(G_SPH,1)*1e-3);
    
   t_grad=(0:(size(k_cal,1)-1))*parameter.GRAD_RASTER_TIME_DEFAULT-5;
   
    if(contains(parameter.SpiralTypeName,{'DoubleSpiral','SpiralInAndOut'},'IgnoreCase',true))
        premom=(k_cal(floor(size(k_cal,1)/2),2:4,:));       
        k_cal(:,2:4,:)=-premom+k_cal(:,2:4,:);   
        parameter.GradDelay=parameter.GradDelay+1*parameter.DwellTime*1e-3; 
    elseif(contains(parameter.SpiralTypeName,'SpiralIn','IgnoreCase',true))
        premom=(k_cal(end,2:4,:));
%          kmax=2*pi*(0.5/(parameter.Resolution*1e-3));
%         phase=2*(0:(parameter.Ninterleaves-1))*pi/parameter.Ninterleaves;
%         for m=1:parameter.Ninterleaves
%             RotMat=[cos(phase(m)) -sin(phase(m)) 0; sin(phase(m)) cos(phase(m)) 0; 0 0 1 ];
%             premom2(1,:,m)=RotMat*[-1*kmax;0;0];
%         end
%         premom2(1,2,:)=premom2(1,2,:)*-1;
        k_cal(:,2:4,:)=-premom+k_cal(:,2:4,:);
        parameter.GradDelay=parameter.GradDelay;
    elseif(contains(parameter.SpiralTypeName,'RIO','IgnoreCase',true))
        premom=(k_cal(floor(size(k_cal,1)/2)+1,2:4,:));
        k_cal(:,2:4,:)=-premom+k_cal(:,2:4,:);
        parameter.GradDelay=parameter.GradDelay+ 2*parameter.DwellTime*1e-3; 
        
     elseif(contains(parameter.SpiralTypeName,'ROI','IgnoreCase',true))
%          premom=k_cal(1,2:4,:);
%          k_cal(:,2:4,:)=-premom+k_cal(:,2:4,:);
        parameter.GradDelay=parameter.GradDelay+2*parameter.DwellTime*1e-3; 
    end
    dwelltime=parameter.DwellTime*1e-3; %us
    for i=1:size(k_cal,2)
        if(size(k_cal,2)==4)
            if(i~=1)
            t_k= [0:dwelltime:dwelltime*(parameter.ADCLength)-1 ]-parameter.GradDelay(i-1);
            else
                t_k= [0:dwelltime:dwelltime*(parameter.ADCLength)-1 ]-parameter.GradDelay(1);
            end
        else
            t_k= [0:dwelltime:dwelltime*(parameter.ADCLength)-1 ]-parameter.GradDelay(1); 
        end
        
    %extrapolation error handling
    t_k(t_k<min(t_grad))=min(t_grad);
    t_k(t_k>max(t_grad))= max(t_grad);
    k2(:,i,:)=interp1(t_grad,k_cal(:,i,:),t_k);
    end
     
    t=t_k;%us
    
    if(contains(parameter.SpiralTypeName,'SpiralOut'))
        t=t+parameter.TE(1);%us
    elseif(contains(parameter.SpiralTypeName,'DoubleSpiral','IgnoreCase',true))
        t=t+parameter.TE(1)-max(t)/2;%us
    else
        t=t+parameter.TE(1); %us
    end
end





