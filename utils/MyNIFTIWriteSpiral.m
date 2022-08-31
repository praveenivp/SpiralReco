function [nifti_header]=MyNIFTIWriteSpiral(Vol,twix,filename,Description)
%[nifti_header]=MyNIFTIWriteSpiral(Vol,twix,filename,Description)
%
%[INPUTS]:
% Vol : upto 7D Input volume
%       Please make sure the first three physical dim are in this order: PHASExREADxPARTITION
% twix : twix object from MAPVBVD
% filename : string (optinal)
% Description : string (optional) 
%
%Example:
%   MyNIFTIWrite(Image_volume,twix_obj);
%   MyNIFTIWrite(Image_volume,twix_obj,'test.nii','Some Documentaion here');
%
%References:
% https://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
% https://brainder.org/2012/09/23/the-nifti-file-format/
%
%praveen.ivp
%

%input paramter handling
if(nargin<3)
    filename=[twix.hdr.Config.ProtocolName '.nii'];
end
if(nargin<4)
    Description='';
end


sa=twix.hdr.Phoenix.sSliceArray.asSlice{1};
kspst=twix.hdr.Phoenix.sKSpace;
R_PE=twix.hdr.MeasYaps.sPat.lAccelFactPE;
R_3D=twix.hdr.MeasYaps.sPat.lAccelFact3D;


if(kspst.ucTrajectory~=1) % then Spiral
    volTR=kspst.lRadialViews*kspst.lPartitions*twix.hdr.Phoenix.alTR{1}*1e-6/(R_3D*R_PE); %s
    sp=getSpiralPara(twix);
    ro=(2*sp.ADCLength*sp.DwellTime)/1e6; % ms
    if(isempty(Description))
        Description=sprintf('R%dx%dC%d TR=%.1fms RO=%.2fms vTR=%.1fs',sp.R_PE,sp.R_3D,sp.CAIPIShift,sp.TR/1e3,ro,volTR);
    end
else
    volTR=kspst.lPhaseEncodingLines*kspst.lPartitions*twix.hdr.Phoenix.alTR{1}*1e-6/(R_3D*R_PE);
end

FOV_PRS=[sa.dPhaseFOV sa.dReadoutFOV sa.dThickness];
% include slice oversampling in the FOV
FOV_PRS(3)=sa.dThickness.*(kspst.lPartitions/kspst.lImagesPerSlab); 
res=sa.dReadoutFOV/kspst.lBaseResolution;
Res_PRS=[kspst.dPhaseResolution*res   res FOV_PRS(3)/kspst.lPartitions];

pos_SCT=mygetfield(twix.hdr.Phoenix.sSliceArray.asSlice{1}.sPosition);
normal_SCT=mygetfield(twix.hdr.Phoenix.sSliceArray.asSlice{1}.sNormal);

if(isfield(sa,'dInPlaneRot'))
    inplanerot=sa.dInPlaneRot; %rad
else
    inplanerot= 0;
end
%from the above five it should be possible to find the affine matrix
[~,MainOrientation] = max(abs(Normal));
if(MainOrientation~=3)
    warning('Nifti Affine matrix might not work properly')
end
Vol=flip(permute(Vol,[2 1 3 4 5]), 3);
% calculate rotation  matrix
[RM,iRM]=getRotmat(normal_SCT,inplanerot);

% Roatation matrix to Affine transform
%patient coordinates(deined SCT here ) is LPS:  +x is Left, +y is posterior and +z is Superior
%NIFTI  needs RAS : +x is Right, +y is Anterior and +z is Superior
LPS2RAS=[-1 -1 1];

AffineMat=zeros(4,4);
cRM=(RM);
AffineMat(1:3,1:3)=(RM*(diag(Res_PRS(:))) )'*diag(LPS2RAS);
inplane_flip=1-2*[any(iRM(:,1)<0) any(iRM(:,2)<0)   any(iRM(:,3)<0)];
AffineMat(1:3,1:3)=AffineMat(1:3,1:3).*(inplane_flip');

%   AffineMat= [ -0.5882    0.0000    0.0000         0  ;...
%     0.0000    0.5415    0.2298         0;...
%    -0.0000   -0.2344    0.5523         0;...
%   100.5882  -82.5184  -49.1938    1.0000];

 % predict translation vector
 pos_RAS=LPS2RAS.*pos_SCT;
  CenterVoxel=floor(0.5*[size(Vol,1) size(Vol,2) size(Vol,3) ]-[1 1 1]); 
%  RM'*CenterVoxel(:)+trans_vec(:) %gives position of center voxel in Patient coordinates (but in RAS convention)
 trans_vec= pos_RAS(:)-AffineMat(1:3,1:3)'*CenterVoxel(:);

  AffineMat(4,:)=[trans_vec' 1];
 



%% try to write the nifti
deafult_header=images.internal.nifti.niftiImage.niftiDefaultHeader(Vol, true, 'NIfTI1');
%modify
deafult_header.pixdim=[1 Res_PRS volTR];
% deafult_header.pixdim=[1 0.8333 0.8000 0.8333 0.0070 0 0 0];
deafult_header.sform_code=1;
deafult_header.srow_x=AffineMat(:,1)';
deafult_header.srow_y=AffineMat(:,2)';
deafult_header.srow_z=AffineMat(:,3)';
deafult_header.descrip=Description;
deafult_header.dim_info=bin2dec(strcat('11','10','01'));% (slice,phase,read) -> 3 2 1
deafult_header.xyzt_units=setSpaceTimeUnits('Millimeter', 'Second');
deafult_header.cal_max=max(Vol(:));
deafult_header.cal_max=min(Vol(:));

% qform not tested but should work
deafult_header.qform_code=0;
quat=rotm2quat(cRM);
deafult_header.quatern_b=-1*quat(2);
deafult_header.quatern_c=quat(3);
deafult_header.quatern_d=quat(4);
deafult_header.qoffset_x=trans_vec(1);
deafult_header.qoffset_y=trans_vec(2);
deafult_header.qoffset_z=trans_vec(3);

deafult_header.scl_slope=1;

%simplify header
NV = images.internal.nifti.niftiImage(deafult_header);
nifti_header=NV.simplifyStruct();
nifti_header.raw=deafult_header;
% 
niftiwrite(Vol,filename,nifti_header,'Compressed',false)
end

%% supporting fucntion
function [RotMatrix,InPlaneRotMatrix]=getRotmat(Normal,InplaneRot)

% Determine dominant MainOrientation
[~,MainOrientation] = max(abs(Normal));

%% Transformation from logical to patient coordinate system
if(MainOrientation==1) % SAG
    InPlaneRotMatrix = [    0                    0                     1; ...
        cos(InplaneRot)  sin(InplaneRot)   0;...
        -sin(InplaneRot)  cos(InplaneRot)   0];
    initNormal = [1; 0; 0];
elseif(MainOrientation==2) %COR
    InPlaneRotMatrix = [    cos(InplaneRot)  sin(InplaneRot)   0;...
        0                    0                     1; ...
        sin(InplaneRot)  -cos(InplaneRot)  0];
    initNormal = [0; 1; 0];
elseif(MainOrientation==3) %TRA
    InPlaneRotMatrix = [    sin(InplaneRot)  -cos(InplaneRot)  0;...
        cos(InplaneRot)  sin(InplaneRot)   0; ...
        0                    0                     1];
    initNormal = [0; 0; 1];
end


%% Compute rotation matrix to align the logical coordinate normal vector with normal vector in the patient coordiate system [SAG, COR, TRA]
if(norm(Normal)>1.001 || norm(Normal)<0.999 )
    error('SODA.Normal must be unit vector')
end

v = cross(initNormal,Normal);
s = norm(v);     % sine of angle
c = dot(initNormal,Normal);   % cosine of angle

if(s == 0)
    RotMatrix = eye(3)*c;
else
    V = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    RotMatrix = eye(3)  + V + V*V*(1-c)/s^2;
end

end

function field=mygetfield(Struct,Fieldname)

if(nargin<2)
    field=zeros(1,3);
    if(isfield(Struct,'dCor'))
        field(2)=Struct.('dCor');
    end
    if(isfield(Struct,'dSag'))
        field(1)=Struct.('dSag');
    end
    if(isfield(Struct,'dTra'))
        field(3)=Struct.('dTra');
    end
else
    field=0;
    if(isfield(Struct,Fieldname))
        field=struct.(Fieldname);
    end
end
end

        function xyztCode = setSpaceTimeUnits(spaceUnitText, timeUnitText)
        %setSpaceTimeUnits: convert simplified space time units to standard
        %form.
        %    This is a helper function to convert simplified space and time
        %    units to the raw format as specified in the NIfTI header.
        
            spaceKey   = {0, 1, 2, 3};
            spaceValue = {'Unknown', 'Meter', 'Millimeter', 'Micron'};

            spaceMap = containers.Map(spaceValue, spaceKey);
            
            if isempty(find(strcmp(spaceValue, spaceUnitText), 1))
               error(message('images:nifti:spaceUnitNotSupported')); 
            end
            
            spaceUnits = spaceMap(spaceUnitText);

            timeValue = {'None', 'Second', 'Millisecond', 'Microsecond', 'Hertz', 'PartsPerMillion', 'Radian'};
            timeKey = {0, 8, 16, 24, 32, 40, 48};

            timeMap = containers.Map(timeValue, timeKey);
            
            if isempty(find(strcmp(timeValue, timeUnitText), 1))
               error(message('images:nifti:timeUnitNotSupported')); 
            end
            
            timeUnits = timeMap(timeUnitText);

            spaceUnitCode = bitand(uint8(spaceUnits),uint8(7));
            timeUnitCode  = bitand(uint8(timeUnits),uint8(56)); % 0x38

            xyztCode = bitor(spaceUnitCode, timeUnitCode);

        end



%%
% function [nifti_Info]=MyWriteNIFTI_old(V,filename,Description)
% %Write Nifiti file with description. 
% %
% %
% %strip file extension
% [sPath, sFilename, ~] = fileparts( filename );
% filename = fullfile( sPath, sFilename);
% 
% %setup default struct
% nifti_Info=struct("Combined",true,"Compressed",false,"Endian","little","Version","NIfTI1");
% nifti_Info.Info = images.internal.nifti.niftiImage.niftiDefaultHeader(...
%     V, nifti_Info.Combined, nifti_Info.Version);
% nifti_Info.Info.descrip=Description;
% 
% if strcmp(nifti_Info.Endian, 'little')
%     machineFmt = 'ieee-le';
% else
%     machineFmt = 'ieee-be';
% end
% NV = images.internal.nifti.niftiImage(nifti_Info.Info);
% fid = fopen([filename '.nii'], 'w', machineFmt);
% % write header.
% [fid, headerBytes] = NV.writeHeader(fid, machineFmt);
% assert(headerBytes == 348||headerBytes == 540);
% % Write empty data until vox_offset
% skipBytes = double(nifti_Info.Info.vox_offset) - headerBytes;
% fwrite(fid, zeros(1,skipBytes), 'uint8');
% % write image data.
% fid = NV.writeImage(V, fid, machineFmt);
% fclose(fid);
% 
% if nifti_Info.Compressed
%     gzip([filename '.nii'], path);
%     delete([filename '.nii']);
% end