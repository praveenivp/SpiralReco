function [k2,t]=Grad2Traj(G_PRS,parameter,mode)
%Function to get kspace trajectories from the gradient waveforms
% only works for 2D traj (i.e slice gradient is ommited)
%gradient units is mT/m sampled at 10 us(typ)
%3D matrix: Timepoints x Gradient Axis x interleaves

% mode: mode of gradient integration
% PE or My Both are same at the moment

%k2 is calculated kspace
% real is phase and imaginary part is read gradients
% unit is (rad/m)
%t is the time axis in us

g=squeeze(complex(G_PRS(:,1,:),G_PRS(:,2,:)));
if(nargin<3)
    mode='PE';
end

%Get all parameters
Gyro=42.575575e3;%Hz/mT %parameter.GyroH;

if(contains(parameter.SpiralTypeName,'SpiralInAndOut'))
    parameter.ADCLength=parameter.ADCLength*4; % read over sampling and 2 ADC
    dwelltime=(size(g,1)*parameter.GRAD_RASTER_TIME_DEFAULT+parameter.ADCShift)/parameter.ADCLength;
else
    parameter.ADCLength=parameter.ADCLength*2; % read over sampling
    dwelltime=(size(g,1)*parameter.GRAD_RASTER_TIME_DEFAULT+parameter.ADCShift)/parameter.ADCLength;
end
if(strcmpi(mode,'PE'))
    %%
    if(contains(parameter.SpiralTypeName,'SpiralOut','IgnoreCase',true))
        nFillPre=2+floor(parameter.ADCShift/parameter.GRAD_RASTER_TIME_DEFAULT);
        nFillPost=2;
    elseif(contains(parameter.SpiralTypeName,'DoubleSpiral','IgnoreCase',true))
        nFillPre=2+floor(parameter.ADCShift/parameter.GRAD_RASTER_TIME_DEFAULT);
        nFillPost=2;
    else
        nFillPost=2+floor(parameter.ADCShift/parameter.GRAD_RASTER_TIME_DEFAULT);
        nFillPre=2;
    end
    lFilledSamples=size(g,1)+nFillPre+nFillPost;
    
    
    % Calculate time axis for gradient and adc for interpolation
    t_grad=((0:(lFilledSamples-1))-nFillPre+1)*parameter.GRAD_RASTER_TIME_DEFAULT;
    %      t_relative= ((0:parameter.ADCLength-1)+0.5)*dwelltime+parameter.GradDelay;
    if(contains(parameter.SpiralTypeName,'SpiralOut'))
        t_relative= ((0:parameter.ADCLength-1)+0.5)*dwelltime-parameter.ADCShift; %+parameter.GradDelay
    elseif(contains(parameter.SpiralTypeName,'DoubleSpiral','IgnoreCase',true))
        t_relative= ((0:parameter.ADCLength-1)+0.5)*dwelltime-parameter.ADCShift;%+parameter.GradDelay
    else
        t_relative= ((0:parameter.ADCLength-1)+0.5)*dwelltime+parameter.ADCShift;%+parameter.GradDelay
    end
    % calculte k
    k=zeros(lFilledSamples,1);
    k2=zeros(length(t_relative),size(g,2)); %dwell time matched
    for i=1:size(g,2)
        k((2:length(g))+nFillPre)=cumsum((g(1:end-1,i))+(g(2:end,i)))*0.5;
        %postfill Looks wrong + nprefill to get the last sample
        k(length(g)+nFillPre+(1:nFillPost))=[k(length(g)+nFillPre)+0.5*(g(end,i))  k(length(g)+nFillPre)+1*(g(end,i))];
        
        k=k.*parameter.GRAD_RASTER_TIME_DEFAULT*1e-6;
        
        if(contains(parameter.SpiralTypeName,'DoubleSpiral','IgnoreCase',true))
            phase=(parameter.Ninterleaves-i+1)*pi/parameter.Ninterleaves;
            premom=(k(ceil(length(k)/2)));
            %         premom1=(parameter.grad_MOM.phase.PreMoment+1i*parameter.grad_MOM.read.PreMoment).*(1e-6*exp(1i*phase)); %(mT*us)/m -> (mT*s)/m
            %         if(all(abs(premom-premom1)>1e-10))
            %         disp([premom;premom1])
            %         disp(abs(premom-premom1))
            %         warning('all(abs(premom-premom1)>1e-10)')
            %         end
            k=-1*(k-premom);
        end
        k2(:,i)=interp1(t_grad,k,t_relative).';
        
        if (isfield(parameter,'grad_MOM'))
            premom1=(parameter.grad_MOM.phase.PreMoment+1i*parameter.grad_MOM.read.PreMoment).*1e-3; %mT-> T %units mT/m*us
            %     else
            premom=sum(g(1:floor(size(g,1)/2),i).*parameter.GRAD_RASTER_TIME_DEFAULT*1e-6); %mT/m*us
        end
        
        
        
        if(contains(parameter.SpiralTypeName,'SpiralInAndOut'))
            k2(:,i)=k2(:,i)-k2(ceil(length(t_relative)/2),i);
        end
    end
    k2=k2.*(Gyro*2*pi);
    t=0:dwelltime:(size(k2,1)-1)*dwelltime;
    
else
    k_cal=(parameter.gammaH*2*pi)*(parameter.GRAD_RASTER_TIME_DEFAULT*1e-6)*(cumsum(G_PRS,1)*1e-3);
    k1=squeeze(complex(k_cal(:,1,:),k_cal(:,2,:)));
    
    if(contains(parameter.SpiralTypeName,'DoubleSpiral','IgnoreCase',true))
        %         phase=(parameter.Ninterleaves-i+1)*pi/parameter.Ninterleaves;
        premom=(k1(floor(length(k1)/2),:));
        %         premom1=(parameter.grad_MOM.phase.PreMoment+1i*parameter.grad_MOM.read.PreMoment).*(1e-6*exp(1i*phase)); %(mT*us)/m -> (mT*s)/m
        %         if(all(abs(premom-premom1)>1e-10))
        %         disp([premom;premom1])
        %         disp(abs(premom-premom1))
        %         warning('all(abs(premom-premom1)>1e-10)')
        %         end
        k1=(premom-k1);
        
        
    elseif(contains(parameter.SpiralTypeName,'SpiralInAndOut'))
        k1=k1-k1(ceil(length(k1)/2),:);
        
    end
    
   t_grad=(0:(size(k1,1)-1))*parameter.GRAD_RASTER_TIME_DEFAULT;
   
%     dwelltime=((size(k1,1))*parameter.GRAD_RASTER_TIME_DEFAULT+parameter.ADCShift)/parameter.ADCLength;    
%     t_kx= [0 0 0:dwelltime:max(t_grad) max(t_grad)  max(t_grad)]-parameter.GradDelay(1); % delay is moved to GetGradients();
%     t_ky= [0 0 0:dwelltime:max(t_grad) max(t_grad)  max(t_grad)]-parameter.GradDelay(1);
%    
    dwelltime=parameter.DwellTime*1e-3; %us
    t_kx= [0:dwelltime:dwelltime*(parameter.ADCLength)-1 ]-parameter.GradDelay(1); % delay is moved to GetGradients();
    t_ky= [0:dwelltime:dwelltime*(parameter.ADCLength)-1 ]-parameter.GradDelay(1);
    %extrapolation error handling
    t_kx(t_kx<0)=0;
    t_kx(t_kx>max(t_grad))= max(t_grad);
    t_ky(t_ky<0)=0;
    t_ky(t_ky>max(t_grad))= max(t_grad);
    %((0:parameter.ADCLength*2-1)+0.5)*dwelltime+parameter.GradDelay-parameter.ADCShift;
    k2=complex(interp1(t_grad,real(k1),t_kx),interp1(t_grad,imag(k1),t_ky));
    
    
    
    
    
    t=t_ky;%us
    
    if(contains(parameter.SpiralTypeName,'SpiralOut'))
        t=t+parameter.TE(1);%us
    elseif(contains(parameter.SpiralTypeName,'DoubleSpiral','IgnoreCase',true))
        t=t+parameter.TE(1)-max(t)/2;%us
    else
        t=t+parameter.TE(1); %us
    end
    
    %  k2=k2.*(parameter.gammaH*2*pi);
    
end

end




% %% my old calc
%    %% my calculation
%    % Calculate time axis for gradient and adc for interpolation
%     t_grad=((0:(size(g,1)-1+2))-1)*parameter.GRAD_RASTER_TIME_DEFAULT;
% %      t_relative= ((0:parameter.ADCLength-1)+0.5)*dwelltime+parameter.GradDelay;
%     if(contains(parameter.SpiralTypeName,'SpiralOut'))
%         t_relative= ((0:parameter.ADCLength-1)+0.5)*dwelltime+parameter.GradDelay-parameter.ADCShift;
%     else
%         t_relative= ((0:parameter.ADCLength-1))*dwelltime+parameter.GradDelay+parameter.ADCShift;
%     end
%    % calculte k
%    k=zeros(size(g,1)+2,1);
% k2=zeros(length(t_relative),size(g,2)); %dwell time matched
% for i=1:size(g,2)
%     k((2:length(g))+2)=cumsum((g(1:end-1,i)+g(2:end,i)))*0.5;
%     k(end)=k(length(g))+g(end,i);
%     k(end-1)=k(length(g))+0.5*g(end,i);
% %      k(end-2)=k(length(g))+0.25*g(end,i);
% %      k(end-3)=k(length(g))+0.125*g(end,i);
%     k=k.*parameter.GRAD_RASTER_TIME_DEFAULT*1e-6;
%     k2(:,i)=interp1(t_grad,k,t_relative,'linear').';
%
%     if (isfield(parameter,'grad_MOM'))
%          premom1=(parameter.grad_MOM.phase.PreMoment+1i*parameter.grad_MOM.read.PreMoment).*1e-3; %mT-> T
% %     else
%       premom=sum(g(1:floor(size(g,1)/2),i).*parameter.GRAD_RASTER_TIME_DEFAULT*1e-6);
%    end
%
%     if(contains(parameter.SpiralTypeName,'DoubleSpiral'))
%         k2(:,i)=k2(:,i)-premom;
%     end
% end
%  k2=k2.*(Gyro*2*pi);
