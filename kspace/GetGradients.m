function [Grad_PRS,grad_XYZ,grad_MOM]=GetGradients(twix_obj,SpiralPara,soda_obj,cSlc)
% [Grad_PRS,grad_XYZ,grad_MOM]=GetGradients(twix_obj)
%INPUT:
%twix_obj: from mapVBVD of peNC_spiral Sequence
%SpiralPara : struct with seq parameters from getSpiralPara()
%soda_obj : slice orientation data from SODA_OBJ()
%cSlc: current slice when slices have different orientation
%
% OUTPUT:
% GRAD_PRS:   Gradient  in Encoding space
%             3D matrix (time X gradient Axis X interleaves)
% GRAD_XYZ:   Physical Gradients
%             3D matrix (time X gradient Axis X interleaves)


if(nargin<3)
soda_obj=SODA_OBJ('mrprot',twix_obj.hdr);
SpiralPara=getSpiralPara(twix_obj);
cSlc=1;
end


% [G_PRS1,grad_MOM]=CalculateGradient(SpiralPara);//old
[G_PRS,grad_MOM]=getSpiralTraj(SpiralPara); %just gradients


%when the two angles are used to define the plane orientation,
% an offset needs to be added to the interleave dimension 
euler_angle=rotm2eul(soda_obj.RotMatrix{cSlc},"XYZ");
signmat=[1 ; -1 ;-1]; %looks like patient postion matrix
[~,Idx]=max(soda_obj.Normal{cSlc}); % main oriantaion 1->S 2->COR 3->TRA
offset=pi+signmat(Idx)*euler_angle(Idx);%deg2rad(-5.1237);
G_PRS_int=CalculateInterleaves(G_PRS,SpiralPara,offset);
Grad_PRS=permute(G_PRS_int,[2 1 3]); %second dim is axis

% % %%test
% if(contains(parameter.SpiralTypeName,'SpiralOut'))
%     nFillPre=2+floor(parameter.ADCShift/parameter.GRAD_RASTER_TIME_DEFAULT);
%     nFillPost=2;
% else
%     nFillPost=2+floor(parameter.ADCShift/parameter.GRAD_RASTER_TIME_DEFAULT);
%     nFillPre=2;
% end
% lFilledSamples=size(Grad_PRS,2)+nFillPre+nFillPost;
% k=zeros(size(Grad_PRS,1),lFilledSamples,size(Grad_PRS,3));
% for i=1:3
%     for j=1:size(Grad_PRS,3)
%    
%        k(i,2:size(Grad_PRS,2),j)=cumsum((Grad_PRS(i,1:end-1,j))+(Grad_PRS(i,2:end,j)))*0.5;
%         k(size(Grad_PRS,2)+nFillPre+(1:nFillPost))=...
%             [k(size(Grad_PRS,2)+nFillPre)+0.5*(Grad_PRS(i,end,j))...
%             k(size(Grad_PRS,2)+nFillPre)+1*(Grad_PRS(i,end,j))];
%     end
% end
% 
% 
% Grad_PRS=diff(k,1,2);
% 
% test_mat=[    0.8829    0.1984    0.4255;...
%    -0.0000    0.9063   -0.4226;...
%    -0.4695    0.3731    0.8002];

%rotate the gradients
grad_XYZ=zeros(size(G_PRS_int));
for i=1:SpiralPara.Ninterleaves
G1=(soda_obj.InPlaneRotMatrix{1}*G_PRS_int(:,:,i));
%[1 0 0; 0 -1 0 ; 0 0 -1] head-first supine [SCT]->XYZ
rc=[1 0 0; 0 -1 0 ; 0 0 -1]*soda_obj.RotMatrix{1};
% rc=[1 0 0; 0 -1 0 ; 0 0 -1]*test_mat;
 

grad_XYZ(:,:,i)=(rc*G1);
end
grad_XYZ=permute(grad_XYZ,[2 1 3]);

% grad_XYZ=DelayGradients(grad_XYZ,SpiralPara.GradDelay,SpiralPara.GRAD_RASTER_TIME_DEFAULT);
end


function [G_xyz_delayed]=DelayGradients(G_xyz,Delay_us,RasterTime_us)
%G_xyz should be 3D matrix with time x 3 axis x interleaves
%Gradient raster time is assumed to be 10 us
% Delay_us: delay value for three physical axis

if(nargin<3)
    RasterTime_us=10;
end

if(numel(Delay_us)<3)
    Delay_us=ones(3,1)*Delay_us(1);
end

t_nom=0:RasterTime_us:(size(G_xyz,1)-1)*RasterTime_us;
G_xyz_delayed=zeros(size(G_xyz));
for intlv=1:size(G_xyz,3)
    for Gaxis=1:size(G_xyz,2)
        t_delay=t_nom+Delay_us(Gaxis);
        %avoid extrapolation NAN
        t_delay(t_delay<0)=0;
        t_delay(t_delay>max(t_nom)) =max(t_nom);
        G_xyz_delayed(:,Gaxis,intlv)=interp1(t_nom,G_xyz(:,Gaxis,intlv),t_delay,'linear');
        
    end
end
end

function [G_PRS,grad_MOM]=CalculateGradient(parameter)
% GetPESpiralTraj(Nitlv, dResolution, fov, nfov, radius, dGradMaxAmpl, Smax, typ, GRAD_RASTER_TIME);
a1=GetPESpiralTraj(parameter.Ninterleaves, parameter.Resolution,length(parameter.FOV),parameter.FOV, parameter.radius , parameter.MaxGradAmp,parameter.MaxSlewrate, parameter.SpiralType, parameter.GRAD_RASTER_TIME_DEFAULT);

if( parameter.SpiralType==4)
    a1=reshape(a1,[numel(a1)/2 2]);
end

grad.read.Amplitude=a1(1,1);
grad.phase.Amplitude=a1(2,1);
grad.slice.Amplitude=a1(3,1);

grad_MOM.read.Moment=a1(4,1);
grad_MOM.phase.Moment=a1(5,1);
grad_MOM.slice.Moment=a1(6,1);

grad_MOM.read.PreMoment=a1(7,1);
grad_MOM.phase.PreMoment=a1(8,1);
grad_MOM.slice.PreMoment=a1(9,1);

grad_MOM.read.PostMoment=a1(10,1);
grad_MOM.phase.PostMoment=a1(11,1);
grad_MOM.slice.PostMoment=a1(12,1);
a1(1:12,:)=[];

grad.read.shape= reshape(fliplr(a1(1:end/2,:)), [numel(a1)/2 1]); %gx
grad.phase.shape=reshape(fliplr(a1(end/2+1:end,:)), [numel(a1)/2 1]); %gy

%  READ and phase gradients
G_PRS=[ grad.phase.Amplitude *grad.phase.shape ...
    grad.read.Amplitude *grad.read.shape ...
    zeros(size(grad.phase.shape))];
G_PRS=G_PRS.';
%     figure,plot(G(:,1:3,1)),legend('G_x','G_y','G_z');



end

function G_int=CalculateInterleaves(G,parameter,offset)

if(nargin<3)
    offset=0;
end
if(size(G,3)~=parameter.Ninterleaves)
    G_int=zeros(size(G,1), size(G,2) ,parameter.Ninterleaves);
    
    %calculate inplane rotation angle
    if(parameter.SpiralType==3) %double spiral
        phase=offset+1*(0:(parameter.Ninterleaves-1))*pi/parameter.Ninterleaves; %interleaves are rotated only in 180 deg range
    else %else rotate in the 2*pi range
        phase=offset+2*(0:(parameter.Ninterleaves-1))*pi/parameter.Ninterleaves;
    end
%     phase= ;%[phase(end) phase(1:end-1)];
    
    for m=1:parameter.Ninterleaves
        RotMat=[cos(phase(m)) -sin(phase(m)) 0; sin(phase(m)) cos(phase(m)) 0; 0 0 1 ];
        %   k1(:,m)=k(:,1).*exp(-1i*pi*((m)/(M)));
        G_int(:,:,m)=RotMat*G;
        %   if(SpiralType==3)
        %       k1(:,m)=k1(:,m)-k1(ceil(length(k1)/2),m);
        %   end
    end
else
    G_int=G;
end

end


function [q_final,rotMat]=Normal2quaternions(normal_vec,verbose)


if(nargin<2)
    verbose=false;
end

%the main orientation should have postive coefficient
%all angles has to be acute

%
[~,sort_idx]=sort((normal_vec),'descend');
%main orientation vector
  v1=zeros(3,1); 
  v1(sort_idx(1))=1;

% project normal vector to the plane perpendicular to the rAxis 
rAxis1=zeros(3,1);  %not the final one
rAxis1(sort_idx(3))=1;
proj=(normal_vec.*(rAxis1==0))./norm(normal_vec.*(rAxis1==0));


%angle between normal_vec and vector projected
dp=dot( v1,proj );
euAngle= -0.5*acos(dp); 
rAxis1=cross(proj,v1); %find perpendicular
rAxis1=rAxis1/norm(rAxis1); % make it a unit vector
q1=quaternion(cos(euAngle),sin(euAngle)*rAxis1(1),sin(euAngle)*rAxis1(2),sin(euAngle)*rAxis1(3));
v1_rot=rotatepoint(q1,v1'); %this should be same as proj



% The axis close to v_90 is the second rotation axis
v_90=cross(normal_vec,v1_rot); % find perpendicular
% rAxis2=(abs(v_90)==max(abs(v_90))).*(1-2*(v_90<0));
rAxis2=v_90./norm(v_90);
euAngle2=real(-0.5*(acos(dot(v1_rot,normal_vec))));



q2=quaternion(cos(euAngle2),sin(euAngle2)*rAxis2(1),sin(euAngle2)*rAxis2(2),sin(euAngle2)*rAxis2(3));
%   v1_rot=quatrotate(v1,rAxis,euAngle);
q_final=q1*q2;
rotMat=quat2rotm(q_final);


%plot everything
if(verbose)
    svec=[1;0;0];
cvec=[0;1;0];
tvec=[0;0;1];
v1_final=rotatepoint(q_final,v1');
   figure,axis square,hold on,
  title(num2str(rad2deg(euAngle))),
  
  plotVector(normal_vec,'V_{orig}',1)
    plotVector(svec,'SAG',1)
      plotVector(cvec,'COR',1)
      plotVector(tvec,'TRA',1)
      plotVector(v1,'V_{initial}')
       plotVector(proj,'V_{Proj}',0.75)
      plotVector(v1_rot','V_{rot}',0.25)
      plotVector(v1_final','V_{final}')
      
  
   xlabel('SAG (X)'),ylabel('COR(-Y)'),zlabel('TRA(-Z)')
    lim=[-1 1];
  xlim(lim),ylim(lim),zlim(lim)
  grid on,view([175 15]) 
end

if(all(normal_vec~=rotatepoint(q_final,v1')))
    disp(normal_vec)
    disp(rotatepoint(q_final,v1'))
   warn('Result is wrong')
end

end

function []=plotVector(vec,label,offset)
if(nargin<3)
    offset=0.5; %position of label 0-> bottom 1-> top
end
%plot the postion vector i.e origin(0,0,0) and point
 quiver3(0,0,0,vec(1),vec(2),vec(3),'LineWidth',2)
  text(offset*vec(1),offset*vec(2),offset*vec(3),label)
end