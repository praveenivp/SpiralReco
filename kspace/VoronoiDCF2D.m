function DCF = VoronoiDCF2D(kTRaj_complex)
%
% DCF = VoronoiDCF(k_PRS)
% 
% Function calulated density compensation function for 2D trajectory. 
%The area of the peripheral kspace points exponentially decays to zero.    
%
% input: kTRaj_complex kx+ iky
% output: DCF : scaled area 
%
% adapted from :http://web.stanford.edu/class/ee369c/mfiles/voronoidens.m

kx=real(kTRaj_complex)';
ky=imag(kTRaj_complex)';

[row,column] = size(kx);

% uncomment these to plot voronoi diagram
%[vx, vy] = voronoi(kx,ky);
%plot(kx,ky,'r.',vx,vy,'b-'); axis equal

kxy = [kx(:),ky(:)];
kmax=max(abs(kxy(:)));
% returns vertices and cells of voronoi diagram
[V,C] = voronoin(kxy); 
DCF = zeros(size(kxy,1),1);
for j = 1:length(C)
  x = V(C{j},1);
  y = V(C{j},2); 
  lxy = length(x);
    
  %if vextex is outside the Trajectory support
  if(any(sqrt(x.^2+y.^2)>=kmax))
  DCF(j)=0;
  continue;
  end
  
  DCF(j) = abs(sum( 0.5*(x([2:lxy 1]) - x(:)).*(y([2:lxy 1]) + y(:))));
end

DCF = reshape(DCF,row,column);
DCF=DCF';
%if outer kspace- area smoothly decays to zero
for intlv=1:size(DCF,2)
    start_idx=find(DCF(:,intlv)==0,1,'first')-5;
    L=size(DCF,1)-start_idx+1;
    DCF(start_idx:end,intlv)=DCF(start_idx,intlv)*exp(-1*(1:L)/(L/5));
    
end
%scaling to 1?
DCF=DCF./max(DCF(:));

end