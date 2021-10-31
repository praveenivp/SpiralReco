classdef B0map
    
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
        
        Fmap %B0 map in Hz
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
                obj.reco_obj=recoVBVD(obj.filename,'coilcombine','adapt3D','removeos',false);
                
            else
                
                if(isfile(varargin{1}))
                    obj.filename=varargin{1};
                elseif(isfolder(varargin{1}))
                    path=varargin{1};
                    [fn, pathname, ~] = uigetfile(strcat(path,'\*B0*.dat'), 'Pick a B0map .dat file');
                    obj.filename=fullfile(pathname,fn);
                end
                obj.reco_obj=recoVBVD(obj.filename,'coilcombine','adapt3D','removeos',false);
            end
            obj.flags=obj.getFlags();
            obj.soda_obj=SODA_OBJ( 'mrprot',obj.reco_obj.twix.hdr);
            
            obj.TE_s=[obj.reco_obj.twix.hdr.Phoenix.alTE{1:3}]*1e-6; %us
            obj=obj.performB0mapping();
        end
        
        function flags =getFlags(obj,varargin)
            switch nargin    
                case 2
                    flags=varargin{1};
                otherwise %parse name valur pairs
                    p=inputParser;
                    addParameter(p,'UnwrapMode','UMPIRE',@(x) any(validatestring(x,{'none','3EchoFit','SpatialUnwrap','UMPIRE'})));
                    addParameter(p,'doRegularization',true,@(x) islogical(x));
                    addParameter(p,'RegIteration',40,@(x) issclar(x));
                    addParameter(p,'RegBeta',1e1,@(x) issclar(x));
                    addParameter(p,'doMasking',false,@(x) islogical(x));
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
        function obj=performB0mapping(obj)
            im=squeeze(obj.reco_obj.img);
            
            switch obj.flags.UnwrapMode
                case '3EchoFit'
                    %             if(strcmpi(mode,all_modes(1)))          
                    fmap_original=fit_fieldmap_3echo(im(:,:,:,1:3),obj.TE_s(1:3)); % Plilips code
                    
                    obj.Fmap=fmap_original./(2*pi); %Hz
                    
                case 'SpatialUnwrap'
                    
                    t=zeros(size(im));
                    for i=1:size(im,4)
                        t(:,:,:,i)=UnWrap_mex(single(angle(im(:,:,:,i))));
                    end
                    fmap_original=diff((t),1,4);
                    fmap_original(:,:,:,1)=fmap_original(:,:,:,1)./diff(obj.TE_s(1:2));
                    fmap_original(:,:,:,2)=fmap_original(:,:,:,2)./diff(obj.TE_s(2:3));
                    obj.Fmap=fmap_original./(2*pi); %Hz
                case 'UMPIRE'
                    obj.Fmap=UMPIRE_unwrapp_3D(im(:,:,:,1:3),obj.TE_s(1:3));
                otherwise
                    error('Unkown fieldmap computation method\n Supported modes: 3EchoFit,SpatialUnwrap,UMPIRE')
                    
                    
            end
            
            if(obj.flags.doRegularization)
                obj.Fmap=RegularizedFieldMapEstimator(squeeze(obj.reco_obj.img),obj.TE_s,obj.Fmap,obj.flags.RegBeta,obj.flags.RegIteration);  
            end
        end
        
        function obj=PerformSliceSelection(obj,slice_coord,spiralIm)
            if(isempty(obj.soda_obj.Coords))
                obj.soda_obj=obj.soda_obj.calcPixLoc();
            end
            
            if(size(slice_coord.SAG,1)~=size(slice_coord.SAG,2))
                error('you forgot to set proper Matrix size');
            end
            
            
            
            
            
            
            
            %much faster Gridded interpolation: but is your data meshgrid?
%             fieldmap=MakeMeshgrid_interp3(obj.soda_obj.Coords{1}.SAG,obj.soda_obj.Coords{1}.COR,obj.soda_obj.Coords{1}.TRA,...
%                 obj.Fmap(:,:,:,1),...
%                 slice_coord.SAG,slice_coord.COR,slice_coord.TRA,'cubic',0);
%             Magim=MakeMeshgrid_interp3(obj.soda_obj.Coords{1}.SAG,obj.soda_obj.Coords{1}.COR,obj.soda_obj.Coords{1}.TRA,...
%                 squeeze(abs(obj.reco_obj.img(1,:,:,:,1,1,1))),...
%                 slice_coord.SAG,slice_coord.COR,slice_coord.TRA);
            fieldmap= obj.Fmap(:,:,:,1);
            Magim=squeeze(abs(obj.reco_obj.img(1,:,:,:,1,1,1)));
            if(obj.flags.doMasking)
                obj.mask=(imfill(Magim>graythresh(Magim),'holes'));
            else
                obj.mask=1;
            end
            
            
            
            % Try to convert the slice orientation to the Fielmap co-oridanated to use gridded interpolant.
            
            
            
            %super slow scattered interpolant
            %              F_fm=scatteredInterpolant(obj.soda_obj.Coords{1}.SAG(:),obj.soda_obj.Coords{1}.COR(:),obj.soda_obj.Coords{1}.TRA(:),double(col(obj.Fmap(:,:,:,1))));
            %              F_magim=scatteredInterpolant(obj.soda_obj.Coords{1}.SAG(:),obj.soda_obj.Coords{1}.COR(:),obj.soda_obj.Coords{1}.TRA(:),double(col(abs(squeeze(obj.reco_obj.img(1,:,:,:,1,1,1))))));
            %             fieldmap=F_fm(slice_coord.SAG,slice_coord.COR,slice_coord.TRA);
            %             Magim=F_magim(slice_coord.SAG,slice_coord.COR,slice_coord.TRA);
            %             Magim=   fliplr(rot90(Magim));
            %             fieldmap=fliplr(rot90(fieldmap));
            %
            if(obj.flags.doRegistration)
                [obj.Fmap_registered,obj.regIm]=RegisterFieldMap2D(spiralIm,Magim,fieldmap);
            else
                obj.Fmap_registered=fieldmap;
                obj.regIm=cat(4,fieldmap,Magim,spiralIm);
            end
            obj.mask=(imfill(obj.regIm(:,:,2)>graythresh(obj.regIm(:,:,2)),'holes'));
            %             obj.Fmap_registered=obj.Fmap_registered.*obj.mask;
            
            %   obj.Fmap_registered=interp3(obj.soda_obj.Coords{1}.COR(:),obj.soda_obj.Coords{1}.SAG(:),obj.soda_obj.Coords{1}.TRA(:),...
            %                           col(obj.Fmap(:,:,:,1)),...
            %                                         slice_coord.COR(:),slice_coord.SAG(:),slice_coord.TRA(:));
            
            %             [SAG_im,COR_im,TRA_im]=ndgrid(linspace(-0.5*FOV(1),0.5*FOV(1),MatSize(1))-imCenter(1),...
            %                 linspace(-0.5*FOV(1),0.5*FOV(3),MatSize(2))-imCenter(2),...
            %                 -1*imCenter(3));
            
        end
        
    end
    
    
    
end