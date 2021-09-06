function ImRow=im2row3D(im,winSize)

%INPUTS:
% im: 4D matrix [CHAxLINxPARxNCHA]
%KernelSize: 3 element vector
%OUTPUTS:
%imRow:[Nkernels x KernelPoints x NCHA]
%praveenivp


[sx,sy,sz,NCha] = size(im);
Nkernel=(sx-winSize(1)+1)*(sy-winSize(2)+1)*(sz-winSize(3)+1);
ImRow = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1)*(sz-winSize(3)+1),prod(winSize),NCha);
count=0;
for z=1:winSize(3)
for y=1:winSize(2)
    for x=1:winSize(1)
        count = count+1;
        ImRow(:,count,:) = reshape(im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,z:sz-winSize(3)+z,:),...
            Nkernel,1,NCha);
    end
end
end

end