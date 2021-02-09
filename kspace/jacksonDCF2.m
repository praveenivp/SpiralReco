 function [wi]=jacksonDCF2(k,parameter)

 % compute weighting function for non-cartesian traj
 % k - complex space
 %parameter is struct with
%  parameter.FOV 
%  parameter.Resolution 
%  parameter.Ninterleaves
 
gridsize=4*floor(parameter.FOV(1)/(4.*parameter.Resolution) + 0.5);
nsamples = size(k,1);%/parameter.Ninterleaves;
zeta=1;


kmax=max(abs(k(:,1)));
zeta= zeta*2*kmax*(2/gridsize(1));

[~,I]=min(abs(abs(k(:,1))-0.85*kmax));
tmp=abs(k(I,1));

I2=I+(nsamples-I)/4+1;

zeta_sq=zeta*zeta;
wi=ones(size(k));
for ns=1:nsamples
    goal=0;
    kxk=real(k(ns,1));
    kyk=imag(k(ns,1));
    for intlv=1:parameter.Ninterleaves
        for m=1:nsamples
            dx=(real(k(m,intlv))-kxk)^2;
            if(dx<zeta_sq)
                dr=dx+(imag(k(m,intlv))-kyk)^2;
                if(dr<zeta_sq)
                    dr=sqrt(dr);
                    kern=0.5-0.5*cos(2*pi*(1+dr/zeta)/2);
                    goal=goal+kern;
                end
                
                
            end
            
        end
    end
    wi(ns,1)=1/goal;
end
    

cutoff_val=sum(wi(round(I):round(I2),1));
cutoff_val=cutoff_val/max([1 I2-I]);

wi(:,1)=wi(:,1)./cutoff_val;

wi=repmat(wi(:,1),[1 parameter.Ninterleaves]);
%cut DCF at 0.85 of kmax

 end