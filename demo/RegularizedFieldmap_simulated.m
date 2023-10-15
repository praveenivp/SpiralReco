%% Regularized field map demo
% generate true field map in Hz
MS=64; %matrix size
[X,Y]=meshgrid(linspace(-1,1,MS));

b0map_true=(200/(2*pi*0.1))*exp(-0.5* (X.^2 + Y.^2)./(2*0.05))- ...
    (200/(2*pi*0.2))*exp(-0.5* ((X+0.2).^2 + (Y+0.1).^2)./(2*0.1))-...
    (300/(2*pi*0.2))*exp(-0.5* ((X-0.2).^2 + (Y-0.1).^2)./(2*0.1));
b0map_true=b0map_true*4; %lot of warps
figure,imagesc(b0map_true),colorbar,title(' field map [Hz]'),axis image

%% generate multi-echo images with complex noise
m = phantom('Modified Shepp-Logan',MS);
T2star=m.*20e-3;
TE=[2,4,6.5]*1e-3; %second
L=length(TE);
im_size=[MS MS];
im=100*m.*... Proton density
    exp(2i*pi*bsxfun(@times,b0map_true,permute(TE,[3 1 2]))).*... Phase
    exp(-1*(bsxfun(@rdivide,permute(TE,[3 1 2]) ,T2star)))... %relaxation effects
    +1e-3*complex(randn([MS MS length(TE)]),randn([MS MS length(TE)])); %noise
im(isinf(im))=0;

figure,
for i=1:L
subplot(2,L,i)
imagesc(abs(im(:,:,i)),[min(abs(im(:))) max(abs(im(:)))]),
colorbar,title(sprintf('Magnitude image(TE=%d ms)',TE(i)*1e3)),axis image
subplot(2,L,i+L)
imagesc(angle(im(:,:,i))),
colorbar,title(sprintf('Phase image(TE=%d ms)',TE(i)*1e3)),axis image
end

im=permute(im,[1 2 4 3]);
%% temporal unwarping and denoising with RegularizedFieldMapEstimator
fm_initial=UMPIRE_unwrapp_3D(im,TE);
[FieldMap,Fm_steps,g_steps]=RegularizedFieldMapEstimator(im,TE,fm_initial,1e-4,10e3);
 range=[min(b0map_true(:)) max(b0map_true(:))];
 % plot
 figure,subplot(231),imagesc([],[],b0map_true,range),
 title('True value [Hz]'),colormap('parula'),colorbar,axis image
 subplot(232),imagesc([],[],fm_initial./(2*pi),range),
 title('UMPIRE B0map [Hz]'),colormap('parula'),colorbar,axis image
 subplot(233),imagesc([],[],Fm_steps(:,:,end)./(2*pi),range),
 title('denoised B0map [Hz]'),colormap('parula'),colorbar,axis image
  subplot(235),imagesc([],[],fm_initial./(2*pi)-b0map_true,[-10 10]),
 title('UMPIRE B0map-true [Hz]'),colormap('parula'),colorbar,axis image
 subplot(236),imagesc([],[],Fm_steps(:,:,end)./(2*pi)-b0map_true,[-10 10]),
 title('denoised B0map-true [Hz]'),colormap('parula'),colorbar,axis image