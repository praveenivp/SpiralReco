function data_noOS=removeOS(data,dim,OSfactor)
%function to remove Oversampling along particular dim. typ readout
%data_noOS=removeOS(data,dim=1,OSfactor=2)
%praveenivp
if(nargin<2)
    dim=1;
    OSfactor=2;
end
sz=size(data);
data_noOS=data;
for cdim=dim
    till_this=round(0.5*(sz(cdim)/OSfactor));
    from_this=sz(cdim)-till_this+1 ;
    keepOS       = [1:till_this, from_this:sz(cdim)]; % data points  to keep
    data_noOS=ifft( data_noOS,[],cdim);
    data_noOS=permute(data_noOS,unique([cdim 1:ndims(data)],'stable'));
    data_noOS=  fft( data_noOS(keepOS,:,:,:,:),[],1);
    data_noOS=ipermute(data_noOS,unique([cdim 1:ndims(data)],'stable'));
end

end