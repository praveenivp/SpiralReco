function B0map=UMPIRE_unwrapp_3D(im,TE)
%function B0map=UMPIRE_unwrapp_3D(im,TE)
%
%INPUTS:
%im: 4D matrix 
%    Complex multi-Echo images
%    First three are physical dimension and 4th should be Echo
%TE : 1D array
%     all echo times should be in seconds.
%OUTPUT:
%B0map :3D matrix
%       The B0map is in Hz
%
%REFERENCE:
% Robinson, S., Schödl, H., & Trattnig, S. (2014).
% A method for unwrapping highly wrapped multi-echo phase images at 
% very high field: UMPIRE. Magnetic Resonance in Medicine, 72(1), 80–92. 
% https://doi.org/10.1002/mrm.24897

%error checking
TE=TE(:);
if(length(TE)~=size(im,4))
error('Number of Echo images and echo times didn''t mactch')
end
deltaTE=( TE(3)-2*TE(2)+TE(1));
if(abs(deltaTE)<0.5e-3)
    warning('DeltaTE is too low %d : The double Phase difference image will be noisy')
else
    fprintf('The maximum range for unwarpping is [%4.2f %4.2f] Hz\n',-0.5/abs(deltaTE),0.5/abs(deltaTE))
end




%The step number corresponds to Steps mentioned in tha paper.
% Step 2: Calculate Phase difference(PD) image
PD=angle(im(:,:,:,2:end).*conj(im(:,:,:,1:2)));
% Step 3: get Double phase difference (DPD) images
DPD= angle(exp(1i*(angle(im(:,:,:,3))-2*angle(im(:,:,:,2))+angle(im(:,:,:,1)) )));
% step 4: Calculate the field map from the DPD images
omega_map=DPD./deltaTE;
%step 5 : calulate n_i: number of wraps in PD images
n_i=round((PD-bsxfun(@times,permute(diff(TE),[4, 3, 2, 1]),omega_map))./(2*pi));
%Step 6 : unwrapped PD images with n_i
PD_unwrap=PD-2*pi*n_i;
%step 7: Calculate the field map from the PD images and get mean
omega_star=bsxfun(@times,PD_unwrap,permute(1./diff(TE),[4, 3, 2, 1]));
omega_star=mean(omega_star,4);
%Step 8: calculate number of unwraps in the Echo images(Phase).
P=angle(im);
n_i=round((P-bsxfun(@times,permute((TE),[4, 3, 2, 1]),omega_star))./(2*pi));
%Step 9: unwrap the Phase images with n_i
P_unwrap=P-2*pi*n_i;
% Step 10: calculate Receiver phase offset.
R=(TE(1)*P_unwrap(:,:,:,2)-TE(2)*P_unwrap(:,:,:,1))./(TE(1)-TE(2));
%Step 11: remove the phase offset from the phase images
P_i0=angle(exp(1i*(P-R)));
%Step 12: calculate number of phase wraps in offset removed phase images
n_i=round((P_i0-bsxfun(@times,permute((TE),[4, 3, 2, 1]),omega_star))./(2*pi));
%Step 13: calculate number of phase wraps in offset removed phase images
P_final= P_i0-2*pi*n_i;

% calculate field map from offset removed unwrapped 
omega_final=mean(bsxfun(@times,diff(P_final,1,4),permute(1./diff(TE),[4, 3, 2, 1])),4);
% convert to Hz
B0map=omega_final./(2*pi);
end