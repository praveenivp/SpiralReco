%% extract and copy 1 volume
% fn='/ptmp/pvalsala/HRXA/moco/moco_DICOM_dzne_ep3d_Pat4_0p8iso_run1_20221221115340_23.nii';
% addpath('/ptmp/pvalsala/Packages/spm12/')

[pathmfile,~]=fileparts(mfilename('fullpath'));
cd(pathmfile)

%%

% [fn,pn]=uigetfile('/ptmp/pvalsala/JIJP-YO7X/moco/allmoco/*.nii');
dir_st=dir('/ptmp/pvalsala/JIJP-YO7X/moco/allmoco/*.nii');
for ii=3:length(dir_st)
fn=fullfile(dir_st(ii).folder,dir_st(ii).name);
% Select the volume that you want to extract (e.g. the 5th volume)
volume_number = 1;

[~,fname,~]=fileparts(fn);
mkdir(fname);
cd(fname);


% Load the 4D time series
V = spm_vol(fn);

% Load the selected volume
data = spm_read_vols(V(volume_number));

% Save the selected volume as a separate 3D image
V_out = V(volume_number);
V_out.n=1;

V_out.fname = sprintf('%s_vol%d.nii',fname,volume_number);
spm_write_vol(V_out, data);


% try image co-registration
!cp ../anat.nii .
!cp ../ves3.nii .
!cp ../GRE.nii .
 coreg_resl3('anat.nii',V_out.fname)
  coreg_resl2('GRE.nii',V_out.fname,'ves3.nii')
  cd('..')
end
  
%% if that fails, try reslicing

  
%%

function coreg_resl3(anat, func)

    VF = spm_vol(func);
    VG = spm_vol(anat);
    flags = struct('sep', [ 4 2 1 0.6], 'cost_fun', 'ecc');


    x = spm_coreg(VG, VF, flags);
    M = spm_matrix(x);

%%%%%%%%%%%%%%

    [~,func_img,~] = fileparts(func);

%     spm_get_space(anat, M*VG.mat);
    files = {func, anat};
    res_flags = struct('mean', false, 'which', 1, 'wrap', [1 1 0], 'interp', 1, 'prefix', ['coreg_' char(func_img) '_']);
    spm_reslice(files, res_flags)
    

end

function coreg_resl2(GRE, func,OV)

    VF = spm_vol(func);
    VG = spm_vol(GRE);
    VG1=spm_vol(OV); % bad contrast volume
    flags = struct('sep', [4 2 1 0.4], 'cost_fun', 'ecc');


    x = spm_coreg(VG, VF, flags);
    M = spm_matrix(x);
    save('res.mat','M')
%%%%%%%%%%%%%%

    spl = split(func, '.');
    func_img = spl(1);

    spm_get_space(GRE, M*VG.mat);
    files = {func, GRE};
    res_flags = struct('mean', false, 'which', 1, 'wrap', [1 1 0], 'interp', 1, 'prefix', ['coregG_' char(func_img) '_']);
    spm_reslice(files, res_flags)
    
     spm_get_space(OV, M*VG1.mat);
    files = {func, OV};
    res_flags = struct('mean', false, 'which', 1, 'wrap', [1 1 0], 'interp', 1, 'prefix', ['coregV_' char(func_img) '_']);
    spm_reslice(files, res_flags)
    

end

