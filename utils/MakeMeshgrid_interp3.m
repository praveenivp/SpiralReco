function Vq = MakeMeshgrid_interp3(X,Y,Z,V,Xq,Yq,Zq,varargin)
% This function is an overload of interp3 to get sub volume from the
% SAG,COR and TRA voxel location matrix from SODA objects
X=permute(X,[2 1 3]);
Y=permute(Y,[2 1 3]);
Z=permute(Z,[2 1 3]);
V=permute(V,[2 1 3 ]);

Xq=permute(Xq,[2 1 3]);
Yq=permute(Yq,[2 1 3]);
Zq=permute(Zq,[2 1 3]);


XYZ=cat(4,X,Y,Z);
XYZq=cat(4,Xq,Yq,Zq);

% %find the order of XYZ to make it a valid mesh grid
 order_idx=[1 2 3];
temp=abs([ sum(diff(X,1,2),'all') ,sum(diff(X,1,1),'all') sum(diff(X,1,3),'all')]);
order_idx(1)=find(temp==max(temp),1);
temp=abs([ sum(diff(Y,1,2),'all') ,sum(diff(Y,1,1),'all') sum(diff(Y,1,3),'all')]);
order_idx(2)=find(temp==max(temp),1);
temp=abs([ sum(diff(Z,1,2),'all') ,sum(diff(Z,1,1),'all') sum(diff(Z,1,3),'all')]);
order_idx(3)=find(temp==max(temp),1);

% V=permute(V,order_idx);

% % find the order in 
% orderq_idx=[1 2 3];
% idx=1;
% temp=[ sum(diff(Xq,1,2),'all') ,sum(diff(Xq,1,1),'all') sum(diff(Xq,1,3),'all')];
% 
% if(sum(temp)~=0)
% orderq_idx(idx)=find(temp==max(temp),1);
% idx=idx+1;
% end
% temp=[ sum(diff(Yq,1,2),'all') ,sum(diff(Yq,1,1),'all') sum(diff(Yq,1,3),'all')];
% 
% if(sum(temp)~=0)
% orderq_idx(idx)=find(temp==max(temp),1);
% idx=idx+1;
% end
% 
% temp=[ sum(diff(Zq,1,2),'all') ,sum(diff(Zq,1,1),'all') sum(diff(Zq,1,3),'all')];
% if(sum(temp)~=0)
% orderq_idx(idx)=find(temp==max(temp),1);
% end
% orderq_idx(3)=setdiff([1,2,3],orderq_idx(1:2));

Vq=interp3(XYZ(:,:,:,order_idx(1)),XYZ(:,:,:,order_idx(2)),XYZ(:,:,:,order_idx(3)) ,V,...
            XYZq(:,:,:,order_idx(1)),XYZq(:,:,:,order_idx(2)),XYZq(:,:,:,order_idx(3)),varargin{:});
        
%undo the stupid permute
%   [~,reorder_idx]=sort(order_idx,'ascend');
    Vq=squeeze(permute(Vq,[2 1 3]));
   Vq=fliplr(rot90(Vq));
end

