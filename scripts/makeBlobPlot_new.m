function [allax,st]= makeBlobPlot(im,blob_nifti,varargin)
st=Parseinput(im,varargin{:});
allax={};
if(ischar(blob_nifti))
    s1=niftiread(blob_nifti);
else
    s1=blob_nifti;
end

%
s1=st.transform(s1);
im=st.transform(im);



if(~isempty(im))
    if (size(im,3)>1)
    im2=createImMontage(abs(im(:,:,st.SlcSel)),st.im_horz);
    else
        im2=abs(im);
    end
    imagesc(im2,st.caxis_im),colormap(gca,'gray'),axis image
allax{1}=gca;
end

hold on

if(~isempty(s1)) % do blob plot only if not empty
    if (size(s1,3)>1)
        blobs=createImMontage(s1(:,:,st.SlcSel),st.im_horz);
    else
        blobs=s1;
    end
    
    blob_mask=double((blobs)>st.Thres(1));
    blobs(~blob_mask)=0;
    
    if(st.Thres(2)<=0)
        st.Thres(2)=max(blobs(:))+st.Thres(2);
    end

%%convert gray blobs to colorfull blobs
bs0=int16((blobs./st.Thres(2))*(size(st.cmap,1)));
bs_rgb=ind2rgb(bs0,st.cmap);
    
    
    
if(~isempty(st.AllAx))
    h2=image(st.AllAx{1},bs_rgb);
    ax1 = st.AllAx{1};                   % gca = get current axis
else
    h2=image(bs_rgb);
    ax1=gca;
    allax{1}=gca;
end
set(h2,'AlphaData',st.alpha_blobs*double(blob_mask)),axis(ax1,'image')
title(ax1,st.title_im,st.title_im_format{:})
ax1.YAxis.Visible = 'off';   % remove y-axis
ax1.XAxis.Visible = 'off';   % remove x-axis
ax1.Box='off';
set(gcf,'Color',[0 0 0]);



%% make colorbar
if(st.colorbar)
if(~st.negBold)
        cb=colorbar();
        set(cb,'Ticks',[],'visible','off')

%     ax=cb.get('axes');
%     cla(ax)
 ax=axes('Position',cb.Position);
cmap_start=(size(st.cmap,1)./st.Thres(2))*st.Thres(1);
image(ax,[],linspace(1,0,size(st.cmap,1)),flip(repmat(permute(st.cmap(cmap_start:end,:),[1 3 2]),[1 100 1]),10))
xticks([]),
yticks([0 0.25 0.5 0.75 1]),
% yticklabels({-1*Thres_N,-1*round(max(s2(:)),1),-1*round(max(s2(:)),1)}),
yticklabels(round(linspace(st.Thres(2),st.Thres(1),5),1)),
ax.YAxisLocation='right';
set(ax,'YAxisLocation','right','Box','off','YColor',0.*[1 1 1]);
ax.YAxis.FontSize=12;
% title('$Z_{scores}$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color',[1 1 1])
else
    cb=colorbar(gca);

    ax=cb.get('axes');
    cla(ax)
% ax=axes('Position',[0.91 0.40 0.02 0.2]);
cmap_start=(size(st.cmap,1)./st.Thres(2))*st.Thres(1);
image(ax,[],linspace(1,0,size(st.cmap,1)),flip(repmat(permute(st.cmap(cmap_start:end,:),[1 3 2]),[1 100 1]),1))
xticks([]),
yticks([0 0.25 0.5 0.75 1]),
% yticklabels({-1*Thres_N,-1*round(max(s2(:)),1),-1*round(max(s2(:)),1)}),
yticklabels(round(linspace(-1*st.Thres(1),-1*st.Thres(2),5),1)),
% ax.YAxis.Visible='off';
% set(ax,'YAxisLocation','right','Box','off','YColor',[1 1 1]);
% title('z-score','Interpreter','latex','FontSize',18,'FontWeight','bold','Color',[1 1 1])
ax.YAxis.FontSize=12;
% image(ax,[],linspace(1,0,size(st.cmap,1)),flip(repmat(permute(st.cmap,[1 3 2]),[1 100 1]),10))
% xticks([]),
% yticks([0 0.25 0.5 0.75 1]),
% % yticklabels({-1*Thres_N,-1*round(max(s2(:)),1),-1*round(max(s2(:)),1)}),
% yticklabels(round(linspace(st.Thres(2),st.Thres(1),5),1)),
% ax.YAxisLocation='right';
% set(ax,'YAxisLocation','right','Box','off','YColor',[1 1 1]);
% title('$Z_{scores}$','Interpreter','latex','FontSize',16,'FontWeight','bold','Color',[1 1 1])
end
end
allax{2}=gca;

end











end

function st=Parseinput(im,varargin)
    
    Nslc=size(im,3);
    im_max=prctile(im(:),98);
    

    p=inputParser;
    p.KeepUnmatched=1;
    addParameter(p,'Thres',[3.1 15],@(x) isvector(x));
    addParameter(p,'SlcSel',floor(0.1*Nslc):ceil(0.9*Nslc),@(x) isvector(x));
    addParameter(p,'caxis_im',[0,im_max],@(x) isvector(x));
    addParameter(p,'alpha_blobs',0.95,@(x) isscalar(x));
    addParameter(p,'cmap',hot(4096),@(x) ismatrix(x));
    addParameter(p,'im_horz',floor(sqrt(Nslc)*(16/9)),@(x) isscalar(x));
    addParameter(p,'transform',@(x)x);
    addParameter(p,'title_im','',@(x) ischar(x))
    addParameter(p,'title_im_format',{'Interpreter','latex','FontSize',18,'FontWeight','bold','Color',[1 1 1]},@(x) iscell(x))
    addParameter(p,'negBold',false,@(x)islogical(x))
    addParameter(p,'colorbar',false,@(x)islogical(x))
    addParameter(p,'AllAx',{})
    
   
    
    
%     addParameter(p,'doDCF','Jackson',@(x) any(strcmp(x,{'none','Jackson','voronoi'})));
    
    parse(p,varargin{:});
    
    st=p.Results;   
end