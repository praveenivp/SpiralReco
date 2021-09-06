function [AtA,A,kernel] = dat2AtA3D(data, KernelSize)

% [AtA,A,kernel] = dat2AtA3D(data, kSize)
%
% Function computes the calibration matrix from calibration data. 
% INPUTS
% data: 4D calib data [COLxLINxPARxCHA]
% KernelSize : 3 element vector [Kx Ky Kz]
% praveenivp



tmp = im2row3D(data,KernelSize); [NKer,NKerPt,Ncha] = size(tmp);
A = reshape(tmp,NKer,NKerPt*Ncha);

AtA = A'*A;
if(nargout>2)
kernel = AtA;
kernel = reshape(kernel,KernelSize(1),KernelSize(2),size(data,4),size(kernel,2));
end
end