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

function slice=getSlicePosition(twix_obj)
p=twix_obj.hdr.Phoenix.sSliceArray;
slice=cell(1,size(p.asSlice,2));
for i=1:size(p.asSlice,2)
if(isfield(p.asSlice{i},'sPosition'))
    slice{i}.PosFieldName=['.dSag';'.dCor';'.dTra';];
    slice{i}.Position=getfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{i}.sPosition);
else
         fprintf('Position: Isocenter\n');
end
if(isfield(p.asSlice{i},'sNormal'))
    slice{i}.Normal=getfield(twix_obj.hdr.Phoenix.sSliceArray.asSlice{i}.sNormal);
end
if(isfield(twix_obj.hdr.Meas,'SliceThickness'))
  slice{i}.thickness=twix_obj.hdr.Meas.SliceThickness; %mm
end

end
end



function field=getfield(structName,Fieldname)

if(nargin<2)
    field=zeros(1,3);
    if(isfield(structName,'dCor'))
        field(2)=eval(strcat('structName','.dCor'));
    end
    if(isfield(structName,'dSag'))
        field(1)=eval(strcat('structName','.dSag'));
    end
    if(isfield(structName,'dTra'))
        field(3)=eval(strcat('structName','.dTra'));
    end
else
    field=0;
    if(isfield(structName,Fieldname))
        field=eval(strcat('structName','.',Fieldname));
    end
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