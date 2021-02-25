classdef B0map<handle
    
    %Use this class only for 3D B0 maps. It takes into account only the
    %
    %Dependencies:
    %recoVBVD
    %UnWrap_mex
    %SODA_OBJ
    
    properties
        reco_obj %reco
        filename
        soda_obj % slice/slab orientation data
        TE
        mode % field map compuataion mode
        
        Fmap %B0 map in Hz
        Fmap_registered
        regIm%fieldmap,cartim ,,Spiralim
        mask
    end
    methods
        function obj=B0map(varargin)
            if(nargin==0)
                path='D:\Data\Spiral\20200918_B0test_sameres';
                [fn, pathname, ~] = uigetfile(strcat(path,'\*.dat'), 'Pick a DATA file');
                obj.filename=fullfile(pathname,fn);
                obj.reco_obj=recoVBVD(obj.filename,'coilcombine','adapt3D');
                
            else
                obj.filename=varargin{1};
                obj.reco_obj=recoVBVD(obj.filename,'coilcombine','adapt3D');
            end
            obj.mode='UMPIRE';
            obj.soda_obj=SODA_OBJ( 'mrprot',obj.reco_obj.twix.hdr);
            
            obj.TE=[obj.reco_obj.twix.hdr.Phoenix.alTE{1:3}]; %us
            obj.performB0mapping();
        end
        
        function obj=performB0mapping(obj)
            all_modes = {'Philip','Chris','UMPIRE'};
            im=squeeze(obj.reco_obj.img);
            
            
            switch obj.mode
                case 'Philip'
%             if(strcmpi(mode,all_modes(1)))
                
                fmap_original=fit_fieldmap_3echo(im(:,:,:,1:3),obj.TE(1:3).*1e-6); % Plilips code
                
                obj.Fmap=fmap_original./(2*pi); %Hz
                
                case 'Chris'
                
                t=zeros(size(im));
                for i=1:size(im,4)
                    t(:,:,:,i)=UnWrap_mex(single(angle(im(:,:,:,i))));
                end
                fmap_original=diff((t),1,4);
                fmap_original(:,:,:,1)=fmap_original(:,:,:,1)./diff(obj.TE(1:2)*1e-6);
                fmap_original(:,:,:,2)=fmap_original(:,:,:,2)./diff(obj.TE(2:3)*1e-6);
                obj.Fmap=fmap_original./(2*pi); %Hz
                case 'UMPIRE'
                    obj.Fmap=UMPIRE_unwrapp_3D(im(:,:,:,1:3),obj.TE(1:3).*1e-6);
                otherwise
                    error('Unkown fieldmap computation method\n Supported modes: Philip,Chris,UMPIRE')
                    
                
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
             fieldmap=MakeMeshgrid_interp3(obj.soda_obj.Coords{1}.SAG,obj.soda_obj.Coords{1}.COR,obj.soda_obj.Coords{1}.TRA,...
     obj.Fmap(:,:,:,1),...
    slice_coord.SAG,slice_coord.COR,slice_coord.TRA,'cubic',0);
  Magim=MakeMeshgrid_interp3(obj.soda_obj.Coords{1}.SAG,obj.soda_obj.Coords{1}.COR,obj.soda_obj.Coords{1}.TRA,...
     squeeze(abs(obj.reco_obj.img(1,:,:,:,1,1,1))),...
     slice_coord.SAG,slice_coord.COR,slice_coord.TRA);
%  obj.mask=(imfill(Magim>graythresh(Magim),'holes'));
            
              
            
%             %super slow scattered interpolant
%              F_fm=scatteredInterpolant(obj.soda_obj.Coords{1}.SAG(:),obj.soda_obj.Coords{1}.COR(:),obj.soda_obj.Coords{1}.TRA(:),double(col(obj.Fmap(:,:,:,1))));
%              F_magim=scatteredInterpolant(obj.soda_obj.Coords{1}.SAG(:),obj.soda_obj.Coords{1}.COR(:),obj.soda_obj.Coords{1}.TRA(:),double(col(abs(squeeze(obj.reco_obj.img(1,:,:,:,1,1,1))))));
%             fieldmap=F_fm(slice_coord.SAG,slice_coord.COR,slice_coord.TRA);
%             Magim=F_magim(slice_coord.SAG,slice_coord.COR,slice_coord.TRA);
%             Magim=   fliplr(rot90(Magim));
%             fieldmap=fliplr(rot90(fieldmap));
            
            if(false)
            [obj.Fmap_registered,obj.regIm]=RegisterFieldMap2D(spiralIm,Magim,fieldmap);
            else
                obj.Fmap_registered=fieldmap;
                obj.regIm=cat(3,fieldmap,Magim,spiralIm);
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