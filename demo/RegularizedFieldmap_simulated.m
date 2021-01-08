%% Regularized field map demo
% generate true field map in Hz
MS=64; %matrix size
[X,Y]=meshgrid(linspace(-1,1,MS));

b0map_true=(200/(2*pi*0.1))*exp(-0.5* (X.^2 + Y.^2)./(2*0.05))- ...
    (200/(2*pi*0.2))*exp(-0.5* ((X+0.2).^2 + (Y+0.1).^2)./(2*0.1))-...
    (300/(2*pi*0.2))*exp(-0.5* ((X-0.2).^2 + (Y-0.1).^2)./(2*0.1));

figure,imagesc(b0map_true),colorbar,title(' field map(in Hz)');

%% generate echoe images with complex noise
m = phantom('Modified Shepp-Logan',MS);
T2star=m.*100e-3;
TE=[2,3,7]*1e-3; %second
L=length(TE);
im_size=[MS MS];
im=100*m.*... Proton density
    exp(2i*pi*mod(bsxfun(@times,b0map_true,permute(TE,[3 1 2])),0.5)).*... Phase
    exp(-1*(bsxfun(@rdivide,permute(TE,[3 1 2]) ,T2star)))... %relaxation effects
    +1e-3*complex(randn([MS MS length(TE)]),randn([MS MS length(TE)])); %noise
im(isinf(im))=0;

figure,
for i=1:L
subplot(2,L,i)
imagesc(abs(im(:,:,i)),[min(abs(im(:))) max(abs(im(:)))]),
colorbar,title(sprintf('Magnitude image(TE=%d ms)',TE(i)*1e3))
subplot(2,L,i+L)
imagesc(angle(im(:,:,i))),
colorbar,title(sprintf('Phase image(TE=%d ms)',TE(i)*1e3))
end

im=permute(im,[1 2 4 3]);
%%
[FieldMap,Fm_steps,g_steps]=RegularizedFieldMapEstimator(im,TE,[],1e-3,200);
 range=[min(b0map_true(:)) max(b0map_true(:))];
 figure,subplot(221),imagesc([],[],b0map_true,range),title('True value(Hz)'),colormap('parula'),colorbar
 subplot(222),imagesc([],[],Fm_steps(:,:,1)./(2*pi),range),title('Initial B0map(Hz)'),colormap('parula'),colorbar
 subplot(223),imagesc([],[],Fm_steps(:,:,end)./(2*pi),range),title('final B0map(Hz)'),colormap('parula'),colorbar
 subplot(224),imagesc([],[],Fm_steps(:,:,end)./(2*pi)-b0map_true,[-10 10]),title('final B0map-true(Hz)'),colormap('parula'),colorbar