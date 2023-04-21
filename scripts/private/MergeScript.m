%% Goes though \proc folder and merge part files for nifti creation and Nordic
addpath(genpath('/ptmp/pvalsala/MATLAB'))

filepattern=getenv('SPIRAL_FN');
if(isempty(filepattern)); filepattern='*M*'; end
fprintf('running mergeScript with file pattern : %s \n' ,filepattern)

AllFolders=dir(fullfile('proc',filepattern));

procPath=fullfile(pwd,'proc');
rawPath=fullfile(pwd,'raw');
if(~isfolder(procPath)&& ~isfolder(rawPath)); error('Bad starting path: it should have \raw \proc'); end

for cFolder=1:length(AllFolders)

cFolderName= AllFolders(cFolder).name;
% check whether files needs to be merged
dirst_part1=dir(fullfile(procPath,cFolderName,sprintf('%s*part01.mat',cFolderName)));


if(isempty(dirst_part1))   
    warning('Processed  directory : %s has no part1 file',cFolderName)
    continue;
end

%twix without extension
MeasUID=cellfun(@(x)str2double(x{1}),(regexp(cFolderName,'M(\d{2,}+)_peSpiral\S*','tokens')));
dir_raw=dir(fullfile(rawPath,sprintf('*#M%d#*.dat',MeasUID)));
twix=mapVBVD(fullfile(rawPath,dir_raw(1).name));
NVol=twix.image.NRep;

cd(fullfile(procPath,cFolderName))
% merge and read

for i=1:length(dirst_part1)
    FilePat=regexprep(dirst_part1(i).name,'part[0-9]+','*');
    dirst=dir(fullfile(pwd,FilePat));
    im_all=[];
    for i=1:length(dirst)
        fprintf('Loading file : %s \n',dirst(i).name)
        load(dirst(i).name)
        im_all=cat(4,im_all,im);
    end
    
    %some error checking
    if(NVol~=size(im_all,4))
        warning('%d volumes expected only %d volumes are reconstructed',NVol,size(im_all,4))
        continue;
    end
    
    %do NORDIC
    im_all_Nordic=NORDIC_imComplex(im_all);
    outFile=sprintf('m%d_B0%s_merged_%s.mat',MeasUID,flags.doB0Corr,datetime('now','Format','d_MMM_y_HH_mm_ss'));
    save(outFile,'im_all','sp','flags','descrip_reco','descrip','fn','im_all_Nordic','-v7.3')
    delete(dirst.name)
    
    % write nifti file
    niiFile=sprintf('%s_B0%s_DCF%s',cFolderName,flags.doB0Corr,flags.doDCF);
    MyNIFTIWriteSpiral(single(abs(im_all)),twix,niiFile);
    MyNIFTIWriteSpiral(single(abs(im_all_Nordic)),twix,strcat('NORDIC_',niiFile));
end

end
