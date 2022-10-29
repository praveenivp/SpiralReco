function [Coilmaps,im,ref]=GetCoilMaps(recoObj,imSize,mode)
if(nargin<3)
    mode='ESPIRIT';
end
if(isnumeric(recoObj))
    ref=recoObj;%[COLxLINxPARxCHA]
    recoObj=struct;
    recoObj.hdr.Meas.MeasUID=6666;
else
    if(isfield(recoObj.twix,'refscan'))
        recoObj.twix.refscan.flagRemoveOS=1;
        ref=recoObj.twix.refscan(:,:,:,:,1);
        ref=performNoiseDecorrAndCoilCompression(permute(ref,[2 1 3 4 5]),recoObj.D,recoObj.V);
        ref=permute(ref,[2,3,4,1]);
        
        if(recoObj.twix.hdr.Phoenix.sKSpace.ucDimension==2)
            imSize=[recoObj.twix.hdr.Phoenix.sKSpace.lBaseResolution recoObj.twix.hdr.Phoenix.sKSpace.lBaseResolution 1];
        else
            imSize=[recoObj.twix.hdr.Phoenix.sKSpace.lBaseResolution recoObj.twix.hdr.Phoenix.sKSpace.lBaseResolution recoObj.twix.hdr.Phoenix.sKSpace.lPartitions];
        end
    else
        error('no reference data');
    end
    
end
ncoil=size(ref,4);
ref1=zpad(flip(flip(flip(ref,1),2),3),[imSize(1),imSize(2),imSize(3),ncoil]);
% im=fftshift(fftshift(fftshift(fft(fft(fft(ref1,[],1),[],3),[],3),1),2),3);
% im=myfft3d(ref1,[imSize,ncoil]);
im=myfft(ref1,[1 2 3]);
switch (mode)
    case 'sos'
        
        im_sos=sum(abs(im),4);
        Coilmaps=bsxfun(@rdivide,im,im_sos);
    case 'adapt2'
        
        %         im=fftshift(fftshift(fftshift(fft(fft(fft(ref1,imSize(1),1),imSize(2),3),imSize(3),3),1),3),3);
        [~,Coilmaps]=adaptiveCombine2(permute(im,[4,1,2,3]),9,true);
        Coilmaps=conj(permute(Coilmaps,[2,3,4,1]));
        
    case 'ESPIRIT'
        cs=24;
        calibSize=[min(cs,size(ref,1)) min(cs,size(ref,2)) min(cs,size(ref,3))];
        shift=[0,0,0,0];
        calib=crop(ref,[calibSize ncoil],[],[],[],shift);
        kSize=[1 1 1]*7; % kernelsize
        if(exist('imSize','var'))
            imSize=[imSize ncoil];
            calib=zpad(calib,imSize+[0 0 0 0]);
        end
        
        %         [im,Coilmaps]=twix2calibcfl(twix,calibSize,kSize,imSize,shift);
        %% Calculate using bart
        % calib_zp=zpad(calib,imSize);
        
        src=sprintf('/mnt/d/Programs/test/calib_measID%d_%s',recoObj.twix.hdr.Meas.MeasUID,datetime('now','Format','HHmmss'));
        targ=regexprep(src,'calib_','sensWSL_');
        writecfl(wsl2winpath(src),single(calib));
        fprintf('\n **** Now executing it WSL*****\n')
        tic
        if(size(ref,3)==1)
            c2=sprintf(' "~/packages/bart/bart ecalib  -k %d -r %d -m 2  %s %s"',kSize(1),calibSize(1),src,targ); %-t 0.00001
        else
            
            c2=sprintf(' "~/packages/bart/bart ecalib  -k %d:%d:%d -r %d:%d:%d -m 2  %s %s"',(kSize),calibSize,src,targ); %t 0.00001
        end
        system(strcat('bash -c ',c2));
        Coilmaps=readcfl(wsl2winpath(targ));
        
        fprintf('It took %.4f s\n',toc);
        
        
            case 'Cluster'
        cs=24;
        calibSize=[min(cs,size(ref,1)) min(cs,size(ref,2)) min(cs,size(ref,3))];
        shift=[0,0,0,0];
        calib=crop(ref,[calibSize ncoil],[],[],[],shift);
        kSize=[1 1 1]*7; % kernelsize
        if(exist('imSize','var'))
            imSize=[imSize ncoil];
            calib=zpad(calib,imSize+[0 0 0 0]);
        end
        
        %         [im,Coilmaps]=twix2calibcfl(twix,calibSize,kSize,imSize,shift);
        %% Calculate using bart
        % calib_zp=zpad(calib,imSize);
        
        src=sprintf('~/ptmp/bart/calib_measID%d_%s',recoObj.twix.hdr.Meas.MeasUID,datetime('now','Format','HHmmss'));
        targ=regexprep(src,'calib_','sensWSL_');
        writecfl(src,single(calib));
        tic
        if(size(ref,3)==1)
            c2=sprintf(' "~/packages/bart/bart ecalib  -k %d -r %d -m 2  %s %s"',kSize(1),calibSize(1),src,targ); %-t 0.00001
        else
            
            c2=sprintf(' " singularity exec -B ~/ptmp/:/ptmp ~/ptmp/gadgetron.sif bart ecalib -k %d:%d:%d -r %d:%d:%d -m 2  %s %s"',(kSize),calibSize,src(2:end),targ(2:end)); %t 0.00001
        end
%         system(strcat('bash -c ',c2));
        Coilmaps{1}=c2;
        Coilmaps{2}=sprintf('Coilmaps=readcfl(''%s'');',targ);
        
        fprintf('It took %.4f s\n',toc);
        % fprintf(' etime=$SECONDS && ~/packages/bart/bart ecalib  -k %d -r %d %s %s && expr $SECONDS - $etime\n',kSize(1),calibSize(1),src,targ);
        
    case 'Gadgetron'
        cs=24;
        calibSize=[min(cs,size(ref,1)) min(cs,size(ref,2)) min(cs,size(ref,3))];
        shift=[0,0,0,0];
        calib=crop(ref,[calibSize ncoil],[],[],[],shift);
        kSize=[1 1 1]*7; % kernelsize
        if(exist('imSize','var'))
            imSize=[imSize ncoil];
            calib=zpad(calib,imSize+[0 0 0 0]);
        end
        calibdata=permute(calib,[4 1 2 3]);
        calibfilename=fullfile(getenv('TEMP'),'Calibdata.h5');
        outfile=fullfile(getenv('TEMP'),'test_csm.h5');
        BARTcmd='ecalib -m 2 -k 6:6:6 -r 24:24:24';
        [pathFolder,calibfilename]=data2ismrmrd(calibdata,BARTcmd, calibfilename);
        
        
        cmdStr{1}=fullfile(pathFolder,'..\IsmrmrdClient-win10-x64-Release\gadgetron_ismrmrd_client ');
        cmdStr{2}=sprintf(' -f %s ',calibfilename);
        cmdStr{3}=sprintf(' -C %s\\..\\gadgetron\\ecalib.xml ',pathFolder);
        cmdStr{4}=' -a 10.41.60.157 -p 9002 ';
        cmdStr{5}=sprintf(' -o %s ',outfile);
        tic
        [status,cmdout] = system(strcat(cmdStr{:}));
        fprintf('It took %.4f s\n',toc);
        %% read the coil maps back
        [Coilmaps,header,file_info]=readH5File(outfile);
end
end
function winpath= wsl2winpath(wslpath)

temp=strsplit(wslpath,'/');
winpath =strcat(temp{3},':\',strjoin(temp(4:end),'\'));
end
function data_out=performNoiseDecorrAndCoilCompression(data,D,V)
%V is coil compression matrix
sz     = size(data);
data_out   = D*data(:,:);
if(~isempty(V))
    data_out=V*data_out;
    sz(1)=size(V,1);
end

data_out    = reshape(data_out,sz);
end

function writecfl(filenameBase,data)
%WRITECFL  Write complex data to file.
%   WRITECFL(filenameBase, data) writes reconstruction data to 
%   filenameBase.cfl (complex float) and its dimensions to filenameBase.hdr.
%
%   Written to edit data with the Berkeley Advanced Reconstruction Toolbox (BART).
%
%   Parameters:
%       filenameBase:   path and filename of cfl file (without extension)
%       data:           array/matrix to be written
%
% Copyright 2013. Joseph Y Cheng.
% Copyright 2016. CBClab, Maastricht University.
% 2012 Joseph Y Cheng (jycheng@mrsrl.stanford.edu).
% 2016 Tim Loderhose (t.loderhose@student.maastrichtuniversity.nl).

    dims = size(data);
    writeReconHeader(filenameBase,dims);

    filename = strcat(filenameBase,'.cfl');
    fid = fopen(filename,'w');
    
    data = data(:);
    
    fwrite(fid,[real(data)'; imag(data)'],'float32');
    fclose(fid);
end

function writeReconHeader(filenameBase,dims)
    filename = strcat(filenameBase,'.hdr');
    fid = fopen(filename,'w');
    fprintf(fid,'# Dimensions\n');
    for N=1:length(dims)
        fprintf(fid,'%d ',dims(N));
    end
    if length(dims) < 5
        for N=1:(5-length(dims))
            fprintf(fid,'1 ');
        end
    end
    fprintf(fid,'\n');
    
    fclose(fid);
end


function CoilMaps=CalESpiritCoilMaps(ref,imSize)

if(nargin<2)
    sz=size(ref);
    imSize=[sz(1:3)];
end

[sx,sy,sz,Nc] = size(ref);

fprintf("Cpar(%d):",sz)
if(sz>1) % make it seperable along PArtition direction
    ref=fftshift(fft(ref,imSize(3),3),3);
end

CoilMaps=zeros([imSize Nc],'double');
for Cpar=1:sz
    fprintf(1,'\b%d',Cpar);
ncalib = 24; % use 24 calibration lines to compute compression
ksize = [6,6]; % kernel size

DATA=squeeze(ref(:,:,Cpar,:));
% Threshold for picking singular vercors of the calibration matrix
% (relative to largest singlular value.
eigThresh_1 = 0.02;
% threshold of eigen vector decomposition in image space. 
eigThresh_2 = 0.95;

% crop a calibration area
calib = crop(DATA,[ncalib,ncalib,Nc]);

% compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels
[k,S] = dat2Kernel(calib,ksize);
idx = max(find(S >= S(1)*eigThresh_1));

[M,~] = kernelEig(k(:,:,:,1:idx),imSize(1:2));

CoilMaps(:,:,Cpar,:)=permute(M(:,:,:,end),[1 2 4 3]);
end
 fprintf('\n')
end