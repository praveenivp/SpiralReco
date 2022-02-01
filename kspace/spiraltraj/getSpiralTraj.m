function [G_PRS,st]=getSpiralTraj(Nintlv,Res_mm,fov,radius,MaxGrad_mT_m,minRiseTime_ms,SpiralType)
%[G_PRS,st]=getSpiralTraj(Nintlv,Res_mm,fov,radius,MaxGrad_mT_m,minRiseTime_ms,SpiralType)
%[G_PRS,st]=mex_SpiralTraj(SpiralPara);
%
%Wrapper function for a MEX binary(mex_SpiralTraj). MEX binary calulates
%the 2D spiral gradients (only the first interleaf) in the Encoding domain.
%It support variable density spiral waveforms. Check the readme.md
% INPUTS:
% Nintlv- Number of Spiral Interleaves
% Res_mm- Resolution in mm
% radius- fraction to divide encoding space for variable density spiral(4 elements 1D vector )
% fov- Field of view in mm in the corresponding region(4 elements 1D vector )
% MaxGrad_mT_m-   Maximum gradient Amplitude in mT/m
% minRiseTime_ms- Minimum Rise time in ms
% SpiralType - Integer scalar
%         SpiralOut = 1,
%         SpiralIn = 2,
%         DoubleSpiral = 3,
%         InAndOut = 4,
%         RIO = 5
%
%OUTPUTS:
% G_PRS- Calulated gradients in encoding space [Grad Axis x Timepoints]
%st - struct with momentum of prephaser and rephaser gradients
%
%eg
%[g,st]=getSpiralTraj(12,1,col([192 192 192 192]),col([0,0.25,.5,1]),42,5,1);
%
%Mex compilation:
% mex mex_SpiralTraj.cpp vdspiral.cpp nonCartesianTraj.cpp
%
%praveen.ivp@gmail.com

if(nargin==1)
    sp=Nintlv;
    FOV=sp.FOV;
    FOV(1:2)=FOV(1:2)*sp.OSCenter;
    [g,st]=mex_SpiralTraj(sp.Ninterleaves,sp.Resolution,FOV,sp.radius,sp.MaxGradAmp,sp.MinRiseTime,sp.SpiralType);
    G_PRS=zeros(3,size(g,1));
    G_PRS(1,:)=g(:,1);
    G_PRS(2,:)=g(:,2);
    return
end

[g,st]=mex_SpiralTraj(Nintlv,Res_mm,fov,radius,MaxGrad_mT_m,minRiseTime_ms,SpiralType);

G_PRS=zeros(3,size(g,1));
G_PRS(1,:)=g(:,1);
G_PRS(2,:)=g(:,2);
end