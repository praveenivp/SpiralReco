function [G_PRS]=GradientXYZ2PRS(G_xyz,soda,cSlc)
%[G_PRS]=GradientXYZ2PRS(G_xyz,twix_obj)
%
% Function to convert Physical gradients to reconstruction space PRS
% Needs SODA_OBJ.m for Rotation matrix calculation
%G_XYZ   Time X 3 axis X interleaves

if(size(G_xyz,2)~=3)
    error('G_XYZ  matrix dimension: Time X 3 aixs X interleaves')
end

%make axis as first dim
G_xyz=permute(G_xyz,[2 1 3]); 
if(~isa(soda,'SODA_OBJ'))
soda=SODA_OBJ('mrprot',soda.hdr);
end
if(nargin<3)
    cSlc=1;
end
% parameter=getSpiralPara(twix_obj);



%rotate the gradients
G_PRS=zeros(size(G_xyz));
 for i=1:size(G_xyz,3)%parameter.Ninterleaves

%[1 0 0; 0 -1 0 ; 0 0 -1] head-first supine [SCT]->XYZ
rc=([1 0 0; 0 -1 0 ; 0 0 -1]*soda.RotMatrix{cSlc});

G1=(rc\G_xyz(:,:,i));
G_PRS(:,:,i)=((soda.InPlaneRotMatrix{cSlc})\G1);
end

G_PRS=permute(G_PRS,[2,1,3]);
% %UNDO OFFET
% %when the two angles are used to define the plane orientation,
% % an offset needs to be added to the interleave dimension 
% euler_angle=rotm2eul(s.RotMatrix{1},"XYZ");
% signmat=[1 ; -1 ;-1]; %looks like patient postion matrix
% [~,Idx]=max(s.Normal{1}); % main oriantaion 1->S 2->COR 3->TRA
% offset= -1*signmat(Idx)*euler_angle(Idx);%deg2rad(-5.1237);
%  RotMat=[cos(offset) -sin(offset) 0; sin(offset) cos(offset) 0; 0 0 1 ];
%  for m=1:parameter.Ninterleaves
%         G_PRS(:,:,m)=RotMat*G_PRS(:,:,m);
% end
% 

end