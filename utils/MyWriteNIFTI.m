%%
function [nifti_Info]=MyWriteNIFTI(V,filename,Description)
%Write Nifiti file with description. 
%
%
%strip file extension
[sPath, sFilename, ~] = fileparts( filename );
filename = fullfile( sPath, sFilename);

%setup default struct
nifti_Info=struct("Combined",true,"Compressed",false,"Endian","little","Version","NIfTI1");
nifti_Info.Info = images.internal.nifti.niftiImage.niftiDefaultHeader(...
    V, nifti_Info.Combined, nifti_Info.Version);
nifti_Info.Info.descrip=Description;

if strcmp(nifti_Info.Endian, 'little')
    machineFmt = 'ieee-le';
else
    machineFmt = 'ieee-be';
end
NV = images.internal.nifti.niftiImage(nifti_Info.Info);
fid = fopen([filename '.nii'], 'w', machineFmt);
% write header.
[fid, headerBytes] = NV.writeHeader(fid, machineFmt);
assert(headerBytes == 348||headerBytes == 540);
% Write empty data until vox_offset
skipBytes = double(nifti_Info.Info.vox_offset) - headerBytes;
fwrite(fid, zeros(1,skipBytes), 'uint8');
% write image data.
fid = NV.writeImage(V, fid, machineFmt);
fclose(fid);

if nifti_Info.Compressed
    gzip([filename '.nii'], path);
    delete([filename '.nii']);
end