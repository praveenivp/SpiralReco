
function resliced_vol=myspm_reslice(def_space, reslice_vols, interp_method,prefix)
% myspm_reslice(def_space, reslice_vols, interp_method,prefix)
%example
% myspm_reslice(def_space, reslice_vols, 'nearest','r')


if(~exist('interp_method','var'))
    interp_method=4; % 4th deg spline
elseif(ischar(interp_method))
    validatestring(interp_method,{'nearest','linear','2nd Spline','3rd Spline','4th Spline','5th Spline'});
    interp_method=find(contains({'nearest','linear','2nd Spline','3rd Spline','4th Spline','5th Spline'},interp_method,'IgnoreCase',true),1)-1;
end
if(~exist('prefix','var'))
    prefix='r'; % 4th deg spline
end

if(~exist('spm_jobman','file'))
    warning('Trying to add spm12');
    addpath('/ptmp/pvalsala/Packages/spm12/')
end


% def_space=dir(def_space);
 def_space={sprintf('%s,1',fullfile(def_space(1).folder,def_space(1).name))};

% collect all volumes to reslice

% reslice_vols=dir(reslice_vols);
fprintf('reslicing %d volumes: %s\n%s\n%s\n%s\n', length(reslice_vols),reslice_vols.name);
source={};idx=1;
for i=1:length(reslice_vols)
    NIIinfo=niftiinfo(fullfile(reslice_vols(i).folder,reslice_vols(i).name));
    for cv=1:NIIinfo.raw.dim(5)
        source{idx}=sprintf('%s,%d',fullfile(reslice_vols(i).folder,reslice_vols(i).name),cv);
        idx=idx+1;
    end
end
source=source';

spm('defaults','fmri');
spm_jobman('initcfg');

matlabbatch{1}.spm.spatial.coreg.write.ref = def_space;
matlabbatch{1}.spm.spatial.coreg.write.source = source;
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = interp_method;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = prefix;
spm_jobman('run',matlabbatch)


% load all volumes back

resliced_vol=[];

for i=1:length(reslice_vols)
    resliced_vol_fn=sprintf('%s',fullfile(reslice_vols(i).folder,strcat(prefix,reslice_vols(i).name)));
    resliced_vol{i}=niftiread(resliced_vol_fn);
end
sz_all=cellfun(@(x1)size(x1), resliced_vol,'UniformOutput',false);
% if(isequal(sz_all{:}))
%    resliced_vol=cat(ndims(resliced_vol{1})+1, resliced_vol{:});
% end

end