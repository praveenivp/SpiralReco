classdef B0map<matlab.mixin.Copyable
    
    %Use this class only for 3D B0 maps.
    %
    %Dependencies:
    %recoVBVD
    %UnWrap_mex
    %SODA_OBJ
    
    properties
        reco_obj %reco
        filename
        soda_obj % slice/slab orientation data
        TE_s
        
        flags
        
        Fmap %B0 map in rad/s
        Fmap_registered
        regIm%fieldmap,cartim ,,Spiralim
        mask
    end
    methods
        function obj=B0map(varargin)
            if(nargin==0)
                path='D:\Data\Spiral\20200918_B0test_sameres';
                [fn, pathname, ~] = uigetfile(strcat(path,'\*.dat'), 'Pick a B0map dat file');
                obj.filename=fullfile(pathname,fn);
                obj.reco_obj=recoVBVD(obj.filename,'coilcombine','adapt3D');
                
            else
                
                if(isfile(varargin{1}))
                    obj.filename=varargin{1};
                elseif(isfolder(varargin{1}))
                    path=varargin{1};
                    [fn, pathname, ~] = uigetfile(strcat(path,'\*B0*.dat'), 'Pick a B0map .dat file');
                    obj.filename=fullfile(pathname,fn);
                end
                obj.reco_obj=recoVBVD(obj.filename,'coilcombine','adapt3D');
            end
            obj.flags=obj.getFlags(varargin{2:end});
            obj.soda_obj=SODA_OBJ( 'mrprot',obj.reco_obj.twix.hdr);
            
            obj.TE_s=[obj.reco_obj.twix.hdr.Phoenix.alTE{1:obj.reco_obj.twix.hdr.Phoenix.lContrasts}]*1e-6; %us
            obj.performB0mapping();
            obj.CalcMask();
        end
        
        function flags =getFlags(obj,varargin)
            switch nargin    
                case 2
                    flags=varargin{1};
                otherwise %parse name valur pairs
                    p=inputParser;
                    addParameter(p,'UnwrapMode','UMPIRE',@(x) any(validatestring(x,{'none','3EchoFit','SpatialUnwrap','UMPIRE'})));
                    addParameter(p,'doRegularization',true,@(x) islogical(x));
                    addParameter(p,'RegIteration',100,@(x) isscalar(x));
                    addParameter(p,'RegBeta',1e0,@(x) issclar(x));
                    addParameter(p,'doMasking',true,@(x) islogical(x)); %just during resampling
                    addParameter(p,'Interpmode','linear',@(x) any(validatestring(x,{'linear','pchip','spline'})));
                    addParameter(p,'doRegistration',false,@(x)islogical(x));
                    
                    addParameter(p,'is3D',(obj.reco_obj.twix.image.NPar>1),@(x)islogical(x));

                    parse(p,varargin{:});
                    flags=p.Results;             
            end
            
            
            
%             flags.UnwrapMode='UMPIRE'; %{'3EchoFit','SpatialUnwrap','UMPIRE'};
%             flags.doRegularization=true; %Boolean
%             flags.RegIteration=20; %scaler
%             flags.RegBeta=1e1;
%             flags.doMasking=true; %boolean: greythresh mag image and fill the image
%             flags.Interpmode='linear';
%             flags.doImageReg=true; %this uses elastics and 1st echo image to do image registration.
%             
%             %auto flags don't touch
%             flags.is3D=(obj.reco_obj.twix.image.NPar>1);       
            
        
        end
        function performB0mapping(obj)
            im=permute(obj.reco_obj.img,[2 3 4 7 1 5 6]);
%             if(size(im,4)<3 && ~strcmpi(obj.flags.UnwrapMode,'SpatialUnwrap'))
%                 warning('Only 2 echoes: Performing spatial Unwrap');
%                 obj.flags.UnwrapMode='SpatialUnwrap';
%             end
            switch obj.flags.UnwrapMode
                case 'none'
                    deltaPhi=single(diff(angle(im),1,4));
                    obj.Fmap=bsxfun(@rdivide,deltaPhi,permute(diff(obj.TE_s(:)),[2 3 4 1]));
                case '3EchoFit'
                    %             if(strcmpi(mode,all_modes(1)))          
                    obj.Fmap=fit_fieldmap_3echo(im(:,:,:,1:3),obj.TE_s(1:3)); % Plilips code                   
                case 'SpatialUnwrap'
%                     Reference Applied Optics, Vol. 46, No. 26, pp. 6623-6635, 2007.
                    deltaPhi=single(diff(angle(im),1,4));
                    for cEcho=1:size(deltaPhi,4)
                        deltaPhi(:,:,:,cEcho)=UnWrap_mex(deltaPhi(:,:,:,cEcho));
                    end
                    obj.Fmap=bsxfun(@rdivide,deltaPhi,permute(diff(obj.TE_s(:)),[2 3 4 1]));
                case 'UMPIRE'
                    obj.Fmap=UMPIRE_unwrapp_3D(im(:,:,:,1:3),obj.TE_s(1:3));
                otherwise
                    error('Unkown fieldmap computation method\n Supported modes: none,3EchoFit,SpatialUnwrap,UMPIRE')         
            end
            
            if(obj.flags.doRegularization)
                obj.Fmap=obj.Fmap(:,:,:,1);
                obj.Fmap=RegularizedFieldMapEstimator(im,obj.TE_s,obj.Fmap,obj.flags.RegBeta,obj.flags.RegIteration);  
            end
        end
        
                    
        function CalcMask(obj)
            im=abs(permute(obj.reco_obj.img(:,:,:,:,1),[2 3 4 7 1 5 6]));
            level=graythresh(im(:,:,:,1));
            obj.mask=imbinarize(im(:,:,:,1),level);
            obj.mask=imclose(obj.mask, strel('disk', 20, 0));
        end
        
        function PerformResampling(obj,twix)
            %twix of spiral data
            obj.flags.sp_st=obj.getResolution(twix,0);
            obj.flags.fm_st=obj.getResolution(obj.reco_obj.twix,0);
            obj.flags.sp_st.pos_PRS=obj.flags.sp_st.pos_PRS-obj.flags.sp_st.res_PRS;
            fm_im=abs(squeeze(obj.reco_obj.img(1,:,:,:,end)));
            
            % center of the volume, normal and inplane rotation should match!!!!!!!!!
            bitf=obj.checkSoda(twix,obj.reco_obj.twix,true);

            [P,R,S]= obj.getMeshgrid(obj.flags.sp_st);
            [p,r,s]= obj.getMeshgrid(obj.flags.fm_st);
            
            
            obj.regIm= interp3(p,r,s,fm_im,P,R,S);
            obj.regIm=ndflip(permute(obj.regIm,[2 1 3 4]),[1 2 3]);
            % fm_interp is in rad/s ready to use with SpiralReco
            if(obj.flags.doMasking)
                obj.Fmap_registered= interp3(p,r,s,obj.Fmap.*obj.mask,P,R,S,'linear',0);
            else
                obj.Fmap_registered= interp3(p,r,s,obj.Fmap,P,R,S,'linear',0);
            end
            obj.Fmap_registered=ndflip(permute(obj.Fmap_registered,[2 1 3 4]),[1 2 3]);
            

        end
        
        function saveFmap(obj,outfile)
            if(~exist('outfile','var'))
            [pathFolder,fn,~]=fileparts(obj.filename);
            outfile=fullfile(pathFolder,strcat(fn,'.mat'));
            end
            fm_fn=obj.filename;
            fm_interp=obj.Fmap_registered;
            pflags=obj.flags;
            fm_im_interp=obj.regIm;
            fm=obj.Fmap;
            Info='All fieldmap are in rad/s';
            save(outfile,'fm_interp','fm','pflags','fm_im_interp','Info','fm_fn');
                      
            
        end
        
    end
    
    
    
    methods(Access=private)
        function st=getResolution(~,twix,rmos)
            sa=twix.hdr.Phoenix.sSliceArray.asSlice{1};
            kp=twix.hdr.Phoenix.sKSpace;
            soda_obj2=SODA_OBJ('mrprot',twix.hdr);
            if(rmos)
                st.FOV_PRS=[sa.dReadoutFOV sa.dPhaseFOV sa.dThickness ];
                st.Npixels=[soda_obj2.NPixelPhase soda_obj2.NPixelReadout kp.lImagesPerSlab];
                st.res_PRS=st.FOV_PRS./st.Npixels;
                
            else
                if(~isfield(kp,'dSliceOversamplingForDialog')) kp.dSliceOversamplingForDialog=0; end
                st.FOV_PRS=[sa.dReadoutFOV sa.dPhaseFOV sa.dThickness + sa.dThickness*kp.dSliceOversamplingForDialog];
                st.Npixels=[soda_obj2.NPixelPhase soda_obj2.NPixelReadout soda_obj2.NPixelSlice];
                st.res_PRS=st.FOV_PRS./st.Npixels;
                
            end
            st.sSlice=sa;
            st.sKspace=kp;
            st.pos_PRS= GradientXYZ2PRS(soda_obj2.Position{1}',twix);
        end
        
        function [P,R,S]= getMeshgrid(~,st)
            [P,R,S]=meshgrid( linspace(-0.5*st.FOV_PRS(1),0.5*st.FOV_PRS(1)-st.res_PRS(1),st.Npixels(1)) -st.pos_PRS(1),...
                linspace(-0.5*st.FOV_PRS(2),0.5*st.FOV_PRS(2)-st.res_PRS(2),st.Npixels(2))-st.pos_PRS(2),...
                linspace(-0.5*st.FOV_PRS(3),0.5*st.FOV_PRS(3)-st.res_PRS(3),st.Npixels(3))-st.pos_PRS(3));
            
            
        end

        function [output]= checkSoda(~,twix1,twix2,verbose)
            
            if(~isvar('verbose'))
                verbose=true;
            end
            s1=SODA_OBJ('mrprot',twix1.hdr);
            s2=SODA_OBJ('mrprot',twix2.hdr);
            % bit          16    15    14    13    12    11    10     9     8     7      6  5  4           3  2  1
            %output is bit map [               |                               Position (S  C  T)| Normal (S  C  T)    ]
            
            bitfield=zeros(1,16);
            %
            bitfield(1:3)=(s1.Normal{1}==s2.Normal{1});
            bitfield(4:6)=(s1.Position{1}==s2.Position{1});
            bitfield(7)=(s1.InplaneRot{1}==s2.InplaneRot{1});
            bitfield(8:10)=[s1.ReadoutFOV==s2.ReadoutFOV s1.PhaseFOV==s2.PhaseFOV s1.Thickness==s2.Thickness];
            bitfield(11:13)=[s1.NPixelReadout==s2.NPixelReadout s1.NPixelPhase==s2.NPixelPhase s1.NPixelSlice==s2.NPixelSlice];
            
            bitfield(14)=(twix1.hdr.Dicom.lFrequency==twix2.hdr.Dicom.lFrequency);
            
            if(verbose)
                if(any(bitfield(1:3)==0))
                    fprintf('normal vector are not equal : \n  (%.2f %.2f %.2f) ~= (%.2f %.2f %.2f)',s1.Normal{1},s2.Normal{1})
                %else
                  %  disp('Normal of both orientation are the same')
                end
                
                if(any(bitfield(4:6)==0))
                    warning(sprintf('position vector are not equal : \n  (%.2f %.2f %.2f) ~= (%.2f %.2f %.2f)',s1.Position{1},s2.Position{1}));
                %else
                    %disp('position of both orientation are the same')
                end
                if(bitfield(7)==0)
                    warning('Inplane rotation of the volumes do not match %.2f rad!=%.2f',s1.InplaneRot{1},s2.InplaneRot{1});
                end
            end
            
            output=(bitfield);
        end
       
        
    end
    
end



% old code
%         function obj=PerformSliceSelection(obj,slice_coord,spiralIm)
%             if(isempty(obj.soda_obj.Coords))
%                 obj.soda_obj=obj.soda_obj.calcPixLoc();
%             end
%             
%             if(size(slice_coord.SAG,1)~=size(slice_coord.SAG,2))
%                 error('you forgot to set proper Matrix size');
%             end
%             
%             
%              
%             
%             %much faster Gridded interpolation: but is your data meshgrid?
% %             fieldmap=MakeMeshgrid_interp3(obj.soda_obj.Coords{1}.SAG,obj.soda_obj.Coords{1}.COR,obj.soda_obj.Coords{1}.TRA,...
% %                 obj.Fmap(:,:,:,1),...
% %                 slice_coord.SAG,slice_coord.COR,slice_coord.TRA,'cubic',0);
% %             Magim=MakeMeshgrid_interp3(obj.soda_obj.Coords{1}.SAG,obj.soda_obj.Coords{1}.COR,obj.soda_obj.Coords{1}.TRA,...
% %                 squeeze(abs(obj.reco_obj.img(1,:,:,:,1,1,1))),...
% %                 slice_coord.SAG,slice_coord.COR,slice_coord.TRA);
%             fieldmap= obj.Fmap(:,:,:,1);
%             Magim=squeeze(abs(obj.reco_obj.img(1,:,:,:,1,1,1)));
%             if(obj.flags.doMasking)
%                 obj.mask=(imfill(Magim>graythresh(Magim),'holes'));
%             else
%                 obj.mask=1;
%             end
%             
%             
%             
%             % Try to convert the slice orientation to the Fielmap co-oridanated to use gridded interpolant.
%             
%             
%             
%             %super slow scattered interpolant
%             %              F_fm=scatteredInterpolant(obj.soda_obj.Coords{1}.SAG(:),obj.soda_obj.Coords{1}.COR(:),obj.soda_obj.Coords{1}.TRA(:),double(col(obj.Fmap(:,:,:,1))));
%             %              F_magim=scatteredInterpolant(obj.soda_obj.Coords{1}.SAG(:),obj.soda_obj.Coords{1}.COR(:),obj.soda_obj.Coords{1}.TRA(:),double(col(abs(squeeze(obj.reco_obj.img(1,:,:,:,1,1,1))))));
%             %             fieldmap=F_fm(slice_coord.SAG,slice_coord.COR,slice_coord.TRA);
%             %             Magim=F_magim(slice_coord.SAG,slice_coord.COR,slice_coord.TRA);
%             %             Magim=   fliplr(rot90(Magim));
%             %             fieldmap=fliplr(rot90(fieldmap));
%             %
%             if(obj.flags.doRegistration)
%                 [obj.Fmap_registered,obj.regIm]=RegisterFieldMap2D(spiralIm,Magim,fieldmap);
%             else
%                 obj.Fmap_registered=fieldmap;
%                 obj.regIm=cat(4,fieldmap,Magim,spiralIm);
%             end
%             obj.mask=(imfill(obj.regIm(:,:,2)>graythresh(obj.regIm(:,:,2)),'holes'));
%             %             obj.Fmap_registered=obj.Fmap_registered.*obj.mask;
%             
%             %   obj.Fmap_registered=interp3(obj.soda_obj.Coords{1}.COR(:),obj.soda_obj.Coords{1}.SAG(:),obj.soda_obj.Coords{1}.TRA(:),...
%             %                           col(obj.Fmap(:,:,:,1)),...
%             %                                         slice_coord.COR(:),slice_coord.SAG(:),slice_coord.TRA(:));
%             
%             %             [SAG_im,COR_im,TRA_im]=ndgrid(linspace(-0.5*FOV(1),0.5*FOV(1),MatSize(1))-imCenter(1),...
%             %                 linspace(-0.5*FOV(1),0.5*FOV(3),MatSize(2))-imCenter(2),...
%             %                 -1*imCenter(3));
%             
