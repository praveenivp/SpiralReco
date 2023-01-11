% extract essential parameter from special card and phoenix of the twix_obj
%[parameter struct]=getSpiralPara(twix_obj)
function [para]= getSpiralPara(twix_obj)


if(nargin<1)
    parameter.Ninterleaves=8;
    parameter.SpiralTypeName='DoubleSpiral';
parameter.ADCShift=0;
parameter.GRAD_RASTER_TIME_DEFAULT=10;%us
parameter.ADCLength=1024; 
parameter.GradDelay=0;
parameter.GammaH=42.575575e6; %Hz/T
para=parameter;
return;
end


if(isfield(twix_obj.hdr.Phoenix,'sWipMemBlock')) 
    %This is VE twix;)
    [para]= getSpiralParaVE(twix_obj);
    return;
end

%2^16=65636 difference 50397444 50462980
LoopOrder={'LininPar','ParinLin'};
%2^8=256 difference   16843009 ; 16843 265; 16843521 ; 16843777             
SpiralType= {'SpiralOut','SpiralIn', 'DoubleSpiral', 'SpiralInAndOut' }; %selSpiraldir

%2^0= 1 shift; 16843777;1684378;  16843779; 16843780;   16843781
SeqType={'GRE','Echo-shifted GRE','PSIF','bSSFP',    'modified ES-GRE'};

% compact_selection=(twix_obj.hdr.Phoenix.sWiPMemBlock.alFree{1});
% para.Spiraltype=SpiralType{bi2de(compact_selection(9:10))+1};
para.LoopOrder=LoopOrder{1}; 

try
sel1=mod((compact_selection-16843009),2^3); %0:GRE +1:shifted +2:PASIF +3 bssfp 
para.SeqType=SeqType{1+sel1};
sel3=mod((compact_selection-16843009)-sel1-sel2*256,2^17);
para.LoopOrder=LoopOrder{1+sel3};


end

%twix_obj.hdr.Phoenix.sNavigatorPara has all para required for recon

para.FOV=cell2mat (twix_obj.hdr.Phoenix.sNavigatorPara.adFree(6:9));
para.OSCenter=twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{2};
para.FOV(1)=para.FOV(1)*para.OSCenter;
para.radius=cell2mat (twix_obj.hdr.Phoenix.sNavigatorPara.adFree(10:13));
if(length(para.radius)==(length(para.FOV)-1))
    para.radius=[0 para.radius];
else
    warn('Dimension of radius and FOV didn''t match')
end
para.PATmode=twix_obj.hdr.Meas.ucPATMode;
para.MinRiseTime=twix_obj.hdr.Phoenix.sNavigatorPara.adFree{2};
para.MaxSlewrate=1000/para.MinRiseTime; %mT/m/ms
para.MaxGradAmp=twix_obj.hdr.Phoenix.sNavigatorPara.adFree{1}; %mT/m
para.Resolution=twix_obj.hdr.Phoenix.sNavigatorPara.adFree{3}; %mm
para.SpiralType= twix_obj.hdr.Phoenix.sNavigatorPara.alFree{4};
para.SpiralTypeName=SpiralType{para.SpiralType};
para.GRAD_RASTER_TIME_DEFAULT=10; %10us
para.gammaH=42.575575e6; %Hz/T
para.DwellTime=twix_obj.hdr.Meas.alDwellTime(1);

para.Ninterleaves= twix_obj.hdr.Phoenix.sKSpace.lRadialViews;
para.NRepetitions=twix_obj.hdr.Meas.lRepetitions;
para.NAverages=twix_obj.hdr.Meas.lAverages;
para.NPartitions=twix_obj.hdr.Meas.NPar;
para.NDummyScans=twix_obj.hdr.Phoenix.sWiPMemBlock.alFree{5}; %number of shots
% para.ADCLength=twix_obj.hdr.Meas.lColSlopeLength;
para.ADCLength=twix_obj.image.NCol/2;
if(isempty(twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{6})||(twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{6})>5)
    para.ADCShift=0;
else
    para.ADCShift=twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{6};
end
%para.ADCLength=twix_obj.hdr.Phoenix.sKSpace.lBaseResolution/2;

% para.Resolution= twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{1};  %mm
para.Koversampling= twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{2};  %frac
para.MaxGradFrac=twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{3};  %fraction
para.SlewrateFrac=twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{4}; %fraction

if(isempty(twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{5}))
    para.GradDelay=0; %us
else
    para.GradDelay= ones(3,1)*twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{5};
end

%Acceleration parameter
para.R_PE=twix_obj.hdr.MeasYaps.sPat.lAccelFactPE;
para.R_3D=twix_obj.hdr.MeasYaps.sPat.lAccelFact3D;
para.CAIPIShift=getCAIPIShift(twix_obj,para);


para.RFPulDur=twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{7};%us
para.FlipAngle=twix_obj.hdr.Meas.FlipAngle; %deg
para.RFTBWproduct=twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{9};
try
para.RFSpoilmoment=twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{10};
catch
    warning('para.RFSpoilmoment=twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{10}');
end
para.TR=twix_obj.hdr.Phoenix.alTR{1};
para.TE(1)=twix_obj.hdr.Phoenix.alTE{1};
try 
    para.TE(2)=twix_obj.hdr.Phoenix.alTE{2};
catch 
   para.TE(2)=0;
end

para.slice=getSlicePosition(twix_obj);

end


function [para]= getSpiralParaVE(twix_obj)

%
%function to get parameters for VE
%


%2^16=65636 difference 50397444 50462980
LoopOrder={'LIN_IN_PAR','PAR_IN_LIN','SQUARE_SPIRAL'};
%2^8=256 difference   16843009 ; 16843 265; 16843521 ; 16843777             
SpiralType= {'SpiralOut','SpiralIn', 'DoubleSpiral', 'ROI','RIO','SpiralDouble' }; 

%2^0= 1 shift; 16843777;1684378;  16843779; 16843780;   16843781
SeqType={'GRE_STRONG','GRE_WEAK','bSSFP'};

para.SeqType=SeqType{1+getCellValue(twix_obj.hdr.Phoenix.sNavigatorPara.alFree,2)};
para.PerformFOVShift=getCellValue(twix_obj.hdr.Phoenix.sNavigatorPara.alFree,4);
para.LoopOrder=LoopOrder{1+getCellValue(twix_obj.hdr.Phoenix.sNavigatorPara.alFree,3)};


%twix_obj.hdr.Phoenix.sNavigatorPara has all para required for recon

para.FOV=ones(1,4)*cell2mat (twix_obj.hdr.Phoenix.sNavigatorPara.adFree(4));
para.OSCenter=cell2mat (twix_obj.hdr.Phoenix.sNavigatorPara.adFree(5));
para.FOV=para.FOV./para.OSCenter;
% para.radius=zeros(size(para.FOV));
para.radius=[0,0.18,.25,1];%cell2mat (twix_obj.hdr.Phoenix.sNavigatorPara.adFree(8:10));
para.GReadDeph.Amplitude=getCellValue( twix_obj.hdr.Phoenix.sNavigatorPara.adFree,1);
para.PhaseCycle=mygetfield(twix_obj.hdr.Phoenix.sFastImaging,'dTrufiPhaseCycle');
% para.GReadReph.Amplitude=getCellValue( twix_obj.hdr.Phoenix.sNavigatorPara.adFree,2);
% para.GReadReph.FlatTopTime=getCellValue( twix_obj.hdr.Phoenix.sNavigatorPara.alFree,2);
% para.GReadReph.RampupTime=getCellValue( twix_obj.hdr.Phoenix.sNavigatorPara.alFree,1);
para.PATmode=twix_obj.hdr.Phoenix.sPat.ucPATMode;
para.MinRiseTime=twix_obj.hdr.Dicom.flGradAmplGp;
para.MaxSlewrate=1000/para.MinRiseTime; %mT/m/ms
para.MaxGradAmp=twix_obj.hdr.Dicom.flGradAmplGr; %mT/m
para.Resolution=para.FOV(1)/twix_obj.hdr.Phoenix.sKSpace.lBaseResolution; %mm
para.SpiralType= twix_obj.hdr.Phoenix.sNavigatorPara.alFree{1}; %twix_obj{2}.hdr.Phoenix.aulServicePara{1}
para.SpiralTypeName=SpiralType{para.SpiralType};
para.GRAD_RASTER_TIME_DEFAULT=10; %10us
para.gammaH=42.575575e6; %Hz/T
para.DwellTime=twix_obj.hdr.Phoenix.sRXSPEC.alDwellTime{1};

para.Ninterleaves= twix_obj.hdr.Phoenix.sKSpace.lRadialViews;
para.NRepetitions=twix_obj.hdr.Config.NRepMeas;
para.NAverages=twix_obj.hdr.Config.Averages;
if(twix_obj.hdr.Phoenix.sKSpace.ucSlicePartialFourier~=16)
    para.NPartitions=twix_obj.hdr.Phoenix.sKSpace.lPartitions;
else
    para.NPartitions=twix_obj.hdr.Config.NPar;
end
para.NDummyScans=twix_obj.hdr.Phoenix.sWipMemBlock.alFree{3}; %number of shots
% para.ADCLength=twix_obj.hdr.Meas.lColSlopeLength;
para.ADCLength=twix_obj.image.NCol/2;
para.ADCShift=0; %this is parameter is removed from seq


% para.Resolution= twix_obj.hdr.Phoenix.sWipMemBlock.adFree{1};  %mm
para.Koversampling= twix_obj.hdr.Phoenix.sWipMemBlock.adFree{2};  %frac
para.MaxGradFrac=twix_obj.hdr.Phoenix.sWipMemBlock.adFree{2};  %fraction
para.SlewrateFrac=twix_obj.hdr.Phoenix.sWipMemBlock.adFree{3}; %fraction

try
para.skope.Gradienfree_us=twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{6};
para.skope.skopeAcq_us=twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{7};
para.skope.TrigADC_us=twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{8};
end

%Acceleration parameter
para.R_PE=twix_obj.hdr.MeasYaps.sPat.lAccelFactPE;
para.R_3D=twix_obj.hdr.MeasYaps.sPat.lAccelFact3D;
para.CAIPIShift=getCAIPIShift(twix_obj,para);

para.TR=twix_obj.hdr.Phoenix.alTR{1};
para.TE(1)=twix_obj.hdr.Phoenix.alTE{1};
try 
    para.TE(2)=twix_obj.hdr.Phoenix.alTE{2};
catch 
   para.TE(2)=0;
end

para.ReadoutLength_ms=(2*para.ADCLength*para.DwellTime)/1e6; % ms
para.vTR_s=(para.TR*para.Ninterleaves*para.NPartitions)/(para.R_PE*para.R_3D*1e6); %s
para.description=(sprintf('R%dx%dC%d TR=%.1fms RO=%.2fms vTR=%.1fs',para.R_PE,para.R_3D,para.CAIPIShift,para.TR/1e3,para.ReadoutLength_ms,para.vTR_s));

%Mostly data delay
if(twix_obj.hdr.Dicom.lFrequency/para.gammaH <8) %3T
    para.GradDelay=ones(3,1)*(3.85+para.DwellTime/2e3);
else
    para.GradDelay=-1*ones(3,1)*5+(2.9023+para.DwellTime*1.4486*1e-3);%-1*ones(3,1)*(2+para.DwellTime/2e3)-[2;0;2];
end
 
 
 
para.RFPulDur=twix_obj.hdr.Phoenix.sWipMemBlock.alFree{4};%ms
para.FlipAngle=twix_obj.hdr.Phoenix.adFlipAngleDegree; %deg
para.RFTBWproduct=twix_obj.hdr.Phoenix.sWipMemBlock.adFree{4};
try
para.RFSpoilmoment=twix_obj.hdr.Phoenix.sWipMemBlock.adFree{10};
catch
    para.RFSpoilmoment=0;
    %warning('para.RFSpoilmoment=twix_obj.hdr.Phoenix.sWipMemBlock.adFree{10} Failed');
end
sa=twix_obj.hdr.Phoenix.sSliceArray.asSlice{1};
kp=twix_obj.hdr.Phoenix.sKSpace;

try
para.FOV_PRS=[sa.dReadoutFOV sa.dPhaseFOV sa.dThickness + sa.dThickness*kp.dSliceOversamplingForDialog]; %mm
catch
    para.FOV_PRS=[sa.dReadoutFOV sa.dPhaseFOV sa.dThickness]; %mm
end
para.MatrixSize=[kp.lBaseResolution  kp.lBaseResolution kp.lPartitions];
para.res_PRS=para.FOV_PRS./para.MatrixSize; %mm


para.slice=getSlicePosition(twix_obj);

end


function slice=getSlicePosition(twix_obj)

p=twix_obj.hdr.Phoenix.sSliceArray;
slice=cell(1,size(p.asSlice,2));
for i=1:size(p.asSlice,2)
if(isfield(p.asSlice{i},'sPosition'))
    slice{i}.PosFieldName=['.dSag';'.dCor';'.dTra';];
    slice{i}.Position=mygetfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{i}.sPosition);
else
    slice{i}.Position=[0;0;0];
         fprintf('Position: Isocenter\n');
end
if(isfield(p.asSlice{i},'sNormal'))
    slice{i}.Normal=mygetfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{i}.sNormal);
end
if(isfield(twix_obj.hdr.Meas,'SliceThickness'))
  slice{i}.thickness=twix_obj.hdr.Meas.SliceThickness; %mm
end
if(isfield(p.asSlice{i},'dThickness') &&isfield(p.asSlice{i},'dReadoutFOV'))
  slice{i}.FOV_PRS=[p.asSlice{i}.dReadoutFOV p.asSlice{i}.dPhaseFOV  p.asSlice{i}.dThickness]; %mm
end

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
        field=Struct.(Fieldname);
    end
end
end


function val=getCellValue(ip_cell,idx)
try
if(isempty(ip_cell{idx}))
    val=0;
else
    val=ip_cell{idx};
end
catch
    val=0;
end

end


function CAIPIShift =getCAIPIShift(twix,para)
is3D=twix.image.NPar>1;
if(is3D)
 try
 CAIPIShift=abs(mean(twix.image.Lin(twix.image.Rep==1 & twix.image.Set==1 & twix.image.Par==1)-...
 twix.image.Lin(twix.image.Rep==1 & twix.image.Set==1 & twix.image.Par==(1+para.R_3D))));
 catch
     % LIN_PAR ordering or acceleration along PAR will break calc
     warning('getSpiralPara():Check the CAIPIShift');
     CAIPIShift=0;
 end

else
CAIPIShift=0;
end

end
%Tried to decode compact selection value


% compact_selection=(twix_obj.hdr.Phoenix.sWiPMemBlock.alFree{1});
% %para.Spiraltype=SpiralType{bi2de(compact_selection(9:10))+1};
% %para.LoopOrder=LoopOrder{1}; 
% 
% 
% sel1=mod((compact_selection-16843009),2^3); %0:GRE +1:shifted +2:PASIF +3 bssfp 
% para.SeqType=SeqType{1+sel1};
% sel2=mod((compact_selection-16843009)-sel1,2^11); %+256*0 out +1 in 2*256:double +3: in &out
% try
%     para.Spiraltype=SpiralType{sel2+1};
% catch
%      para.SpiraltypeNAME=SpiralType{1};
% end
% sel3=mod((compact_selection-16843009)-sel1-sel2*256,2^17);
% para.LoopOrder=LoopOrder{1+sel3};