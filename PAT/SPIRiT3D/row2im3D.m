function [res,W] = row2im3D(ImRow,imSize, KernelSize)
%[res,W] = row2im(mtx,imSize, winSize);
%
%W counts number of times the pixel overlaps

NCha = size(ImRow,3);
sx = imSize(1); sy = imSize(2); sz=imSize(3);
res = zeros(imSize(1),imSize(2),imSize(3),NCha);
W = res;



count=0;
for z=1:KernelSize(3)
for y=1:KernelSize(2)
    for x=1:KernelSize(1)
        count = count+1;
        res(x:sx-KernelSize(1)+x,y:sy-KernelSize(2)+y,z:sz-KernelSize(3)+z,:) = res(x:sx-KernelSize(1)+x,y:sy-KernelSize(2)+y,z:sz-KernelSize(3)+z,:) + reshape(ImRow(:,count,:),(sx-KernelSize(1)+1),(sy-KernelSize(2)+1),(sz-KernelSize(3)+1),NCha);
        W(x:sx-KernelSize(1)+x,y:sy-KernelSize(2)+y,z:sz-KernelSize(3)+z,:) = W(x:sx-KernelSize(1)+x,y:sy-KernelSize(2)+y,z:sz-KernelSize(3)+z,:)+1;
    end
end
end

res = res./W;
