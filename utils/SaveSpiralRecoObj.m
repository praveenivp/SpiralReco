function outFilename=SaveSpiralRecoObj(r)
sp=r.SpiralPara;
flags=r.flags;
fn=r.filename;
ro=(2*sp.ADCLength*sp.DwellTime)/1e6; % ms
 vTR=(sp.TR*sp.Ninterleaves*sp.NPartitions)/(sp.R_3D*sp.R_3D*1e6); %s
 descrip=(sprintf('R%dx%dC%d TR=%.1fms RO=%.2fms vTR=%.1fs',sp.R_PE,sp.R_3D,sp.CAIPIShift,sp.TR/1e3,ro,vTR));
 descrip_reco=sprintf('%s PAT=%s coilcomb=%s B0=%s',flags.CompMode,flags.doPAT, flags.doCoilCombine, flags.doB0Corr);

 if(strcmpi(r.flags.doPAT,'CGSENSE'))
 coilSens=r.coilSens;
 else
     coilSens=[];
 end
 
 if(~strcmpi(r.flags.doB0Corr,'none'))
 B0map=r.B0Map;
 else
     B0map=[];
 end
 
 im=squeeze(r.img);
 dtStr = datetime('now','TimeZone','local','Format','d_MMM_y_HH_mm_ss');
outFilename=sprintf('m%d_reco_%s.mat',r.twix.hdr.Config.MeasUID,dtStr);
 save(outFilename,'sp','flags','im','descrip','descrip_reco','coilSens','B0map','fn');
 
end