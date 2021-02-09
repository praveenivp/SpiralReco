classdef SODA_OBJ
    %SODA Slice orientation date
    %   All values are either in units of mm or radians if not specified
    %   otherwise
    
    %EXAMPLES
%     s=SODA_OBJ('mrprot',twix_obj.hdr); passing twix header
%     s=SODA_OBJ('File',fullfilepath); %passing filename
    
    properties
            Dim
            MainOrientation
            
            ReadoutFOV
            PhaseFOV
            Thickness
            
            Position           
            Normal             
            
            NPixelReadout
            NPixelPhase
            NPixelSlice
            
            NPixelReadoutInterp = 0
            NPixelPhaseInterp   = 0
            NPixelSliceInterp   = 0
            
            PixelSizeReadout
            PixelSizePhase
            PixelSizeSlice   
            
            PixelSizeReadoutInterp = 0
            PixelSizePhaseInterp   = 0
            PixelSizeSliceInterp   = 0 
            
            NSlices
            Coords
            RotMatrix            
            Sliceshift
            InplaneRot
            InPlaneRotMatrix

            LarmorConst = 42.5756 % [MHz/T]
            DefaultFileName = 'Y:\9T_Chris\SODA\SODA.ini';
    end
    
    methods
        function obj = SODA_OBJ(varargin)
            
            if(mod(nargin,2))
                error('Number of arguments must be even: "File",File or "mrprot", mrprot')
            end
            
            filename = [];
            mrprot = [];
            ext = [];
            Name = [];
            
            
            for i=1:2:nargin
                switch varargin{i}
                    case 'File'
                        [pathname,fname,ext] = fileparts(varargin{i+1});
                        filename = [pathname,'\',fname,ext];
                    case 'mrprot'
                        mrprot = varargin{i+1};
                end
            end
            
            
            if(~isempty(mrprot))     
                    obj = extract(obj,mrprot); 
                    obj = calcPixLoc(obj);
                    for i=1:obj.NSlices
                        obj.Sliceshift{i} = dot( obj.Normal{i},obj.Position{i}); % [mm]
                    end   
            else
                 if(strcmp(ext,'.dat'))
                     mrprot = readVB17Header(filename);
                     obj = extract(obj,mrprot); 
                     obj = calcPixLoc(obj);
                     for i=1:obj.NSlices
                        obj.Sliceshift{i} = dot( obj.Normal{1},obj.Position{i}); % [mm]
                     end
                 elseif(strcmp(ext,'.hdr')||strcmp(ext,'.ini'))
                    obj = readSODAfromFile(obj,filename);
                    obj.NSlices = size(obj.Position,2);
                    obj = calcPixLoc(obj);

                     for i=1:obj.NSlices
                        obj.Sliceshift{i} = dot( obj.Normal{i},obj.Position{i}); % [mm]
                     end
                 else
                     error('Wrong file extension. Expecting "*.hdr" or "*.ini"')
                 end
                
            end
            
        end
        
        function obj = setCSIVoxelsize(obj,NRead, NPhase, NSlice)
            obj.NPixelReadoutInterp = NRead;
            obj.NPixelPhaseInterp   = NPhase;
            obj.NPixelSliceInterp   = NSlice;
        end
       
        function obj = readSODAfromFile(obj,filename)
            
            if(isempty(filename))
                warning(['Reading from ',obj.DefaultFileName])
                filename = obj.DefaultFileName;
            end

            %% Get data
            % Slice
            obj.NSlices = cell2mat(inifile(filename,'read',{ 'Parameters','','NSlices','d','none'}));
            
            % Dimension
            obj.Dim = cell2mat(inifile(filename,'read',{ 'Parameters','','Dim','d','none'}));
             
            % FoV
            obj.ReadoutFOV	= cell2mat(inifile(filename,'read',{ 'FOV','','Read','d','none'}));
            obj.PhaseFOV    = cell2mat(inifile(filename,'read',{ 'FOV','','Phase','d','none'}));
            obj.Thickness   = cell2mat(inifile(filename,'read',{ 'FOV','','Slice','d','none'}));
            
            % Matrix Dimension
            obj.NPixelReadout	= cell2mat(inifile(filename,'read',{ 'MatrixSize','','Read','d','none'}));
            obj.NPixelPhase     = cell2mat(inifile(filename,'read',{ 'MatrixSize','','Phase','d','none'}));
            obj.NPixelSlice     = cell2mat(inifile(filename,'read',{ 'MatrixSize','','Slice','d','none'})); 
            
            % Pixelsize            
            obj.PixelSizeReadout   = cell2mat(inifile(filename,'read',{ 'PixelSize','','Read','d','none'}));
            obj.PixelSizePhase     = cell2mat(inifile(filename,'read',{ 'PixelSize','','Phase','d','none'}));
            obj.PixelSizeSlice     = cell2mat(inifile(filename,'read',{ 'PixelSize','','Slice','d','none'}));
             
            for cSli = 1:obj.NSlices
               
                % Rotation
                obj.InplaneRot{cSli} = cell2mat(inifile(filename,'read',{ 'Parameters','','InPlaneRot','d','none'}));

                % Normal
                obj.Normal{cSli}(1) = cell2mat(inifile(filename,'read',{ 'Normal','','SAG','d','none'}));
                obj.Normal{cSli}(2) = cell2mat(inifile(filename,'read',{ 'Normal','','COR','d','none'}));
                obj.Normal{cSli}(3) = cell2mat(inifile(filename,'read',{ 'Normal','','TRA','d','none'}));

                % Position
                obj.Position{cSli} = cell2mat(inifile(filename,'read',{ 'Position','',['POS[',num2str(cSli-1),']'],'d','none'}));
            end
            
            
        end
                
        function obj = extract(obj,mrprot)
            obj.NSlices = length(mrprot.MeasYaps.sSliceArray.asSlice);
            
            for cSlice =1:obj.NSlices

                %% Determine normal
                try
                    NSag = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sNormal.dSag;
                catch
                    NSag = 0;
                end
                try
                    NTra = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sNormal.dTra;
                catch
                    NTra = 0;
                end
                try
                    NCor = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sNormal.dCor;
                catch
                    NCor = 0;
                end

                obj.Normal{cSlice} = [NSag;NCor;NTra];
 
                %% Determine centre of slice 
                try
                    if isstr(mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dSag)
                        dSag = str2num(mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dSag);
                    else
                        dSag = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dSag;
                    end
                catch
                    dSag = 0;
                end
                try
                    if isstr(mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dCor)
                        dCor = str2num(mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dCor);
                    else
                        dCor = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dCor;
                    end
                catch
                dCor = 0;
                end
                try
                    if isstr(mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dTra)
                        dTra = str2num(mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dTra);
                    else
                        dTra = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.sPosition.dTra;
                    end
                    catch
                    dTra = 0;
                end
                obj.Position{cSlice} = [dSag;dCor;dTra];

                %% Inplane rotation
                try
                    if isstr(mrprot.MeasYaps.sSliceArray.asSlice{1, 1}.dInPlaneRot)
                        dRot = str2num(mrprot.MeasYaps.sSliceArray.asSlice{1, 1}.dInPlaneRot);
                    else
                        dRot = mrprot.MeasYaps.sSliceArray.asSlice{1, 1}.dInPlaneRot;
                    end
                catch
                    dRot = 0;
                end
                obj.InplaneRot{cSlice} = dRot;

                %% Slice dimensions
                obj.Thickness  = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dThickness;
                obj.PhaseFOV   = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dPhaseFOV;
                obj.ReadoutFOV = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dReadoutFOV;

                %% Determine dimension
                if(strcmp(mrprot.MeasYaps.sKSpace.ucDimension,'0x2') ||mrprot.MeasYaps.sKSpace.ucDimension==2)  % 2D
                    obj.Dim = 2;
                elseif(strcmp(mrprot.MeasYaps.sKSpace.ucDimension,'0x4') ||mrprot.MeasYaps.sKSpace.ucDimension==4) % 3D
                    obj.Dim = 3;
                end                
                
                %% Pixels
                obj.NPixelPhase   = mrprot.MeasYaps.sKSpace.lPhaseEncodingLines;
                obj.NPixelReadout = mrprot.MeasYaps.sKSpace.lBaseResolution;
                if(obj.Dim==2)
                    obj.NPixelSlice   = 1;
                else
                    obj.NPixelSlice   = mrprot.MeasYaps.sKSpace.lPartitions;
                end

                %% Pixels size
                obj.PixelSizePhase   = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dPhaseFOV/mrprot.MeasYaps.sKSpace.lPhaseEncodingLines;
                obj.PixelSizeReadout = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dReadoutFOV/mrprot.MeasYaps.sKSpace.lBaseResolution;
                
                if(strcmp(mrprot.MeasYaps.sKSpace.ucDimension,'0x4')||mrprot.MeasYaps.sKSpace.ucDimension==4)  % 3D
                    obj.PixelSizeSlice   = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dThickness/mrprot.MeasYaps.sKSpace.lPartitions;
                elseif(strcmp(mrprot.MeasYaps.sKSpace.ucDimension,'0x2') ||mrprot.MeasYaps.sKSpace.ucDimension==2) % 2D
                    obj.PixelSizeSlice   = mrprot.MeasYaps.sSliceArray.asSlice{cSlice}.dThickness;
                else
                    warning('Unknown dimension');
                end

            end
            
        end
        
        
        function obj = calcPixLoc(obj)
            
          
            
            %% Loop over slices
            for cSlice =1:obj.NSlices
                
                % Determine dominant obj.MainOrientation
                [~,obj.MainOrientation{cSlice}] = max(abs(obj.Normal{cSlice}));
                
                %% Calculate pixel positions for slice in logical coordinate system [PHASE, READ, SLICE]
                if(obj.Dim==2) % 2D
                   [PHASE, READ, SLICE] = ndgrid(   linspace(   -obj.PhaseFOV/2, obj.PhaseFOV/2   -obj.PhaseFOV/  obj.NPixelPhase,        obj.NPixelPhase       ),... 
                                                    linspace(    obj.ReadoutFOV/2,-obj.ReadoutFOV/2 +obj.ReadoutFOV/obj.NPixelReadout,    obj.NPixelReadout     ),...
                                                    linspace(   -obj.Thickness/2, obj.Thickness/2  ,       obj.NPixelSlice    )      );
%                 [PHASE, READ, SLICE] = ndgrid(   linspace(   -1, 1   ,        100      ),... 
%                                                 linspace(    -1,1 ,    100    ),...
%                                                  linspace(   -0.05, 0.05  ,       10    )      );
                 
               elseif(obj.Dim==3) % 3D
                   [PHASE, READ, SLICE] = ndgrid(   linspace(   -obj.PhaseFOV/2, obj.PhaseFOV/2   -obj.PhaseFOV/  obj.NPixelPhase,        obj.NPixelPhase       ),... 
                                                    linspace(    obj.ReadoutFOV/2,-obj.ReadoutFOV/2 +obj.ReadoutFOV/obj.NPixelReadout,    obj.NPixelReadout     ),...
                                                    linspace(   -obj.Thickness/2, obj.Thickness/2  -obj.Thickness/ obj.NPixelSlice,       obj.NPixelSlice       )      );
                   
                end
                
                %% Transformation from logical to patient coordinate system           
                if(obj.MainOrientation{cSlice}==1) % SAG
                    obj.InPlaneRotMatrix{cSlice} = [    0                    0                     1; ... 
                                                        cos(obj.InplaneRot{cSlice})  sin(obj.InplaneRot{cSlice})   0;...
                                                        -sin(obj.InplaneRot{cSlice})  cos(obj.InplaneRot{cSlice})   0];
                    initNormal = [1; 0; 0];  
                elseif(obj.MainOrientation{cSlice}==2) %COR
                    obj.InPlaneRotMatrix{cSlice} = [    cos(obj.InplaneRot{cSlice})  sin(obj.InplaneRot{cSlice})   0;...
                                                        0                    0                     1; ... 
                                                        sin(obj.InplaneRot{cSlice})  -cos(obj.InplaneRot{cSlice})  0];
                    initNormal = [0; 1; 0];  
                elseif(obj.MainOrientation{cSlice}==3) %TRA
                    obj.InPlaneRotMatrix{cSlice} = [    sin(obj.InplaneRot{cSlice})  -cos(obj.InplaneRot{cSlice})  0;...
                                                        cos(obj.InplaneRot{cSlice})  sin(obj.InplaneRot{cSlice})   0; ... 
                                                        0                    0                     1];
                    initNormal = [0; 0; 1];                       
                end
                
                SAG = obj.InPlaneRotMatrix{cSlice}(1,1)*PHASE + obj.InPlaneRotMatrix{cSlice}(1,2)*READ + obj.InPlaneRotMatrix{cSlice}(1,3)*SLICE;
                COR = obj.InPlaneRotMatrix{cSlice}(2,1)*PHASE + obj.InPlaneRotMatrix{cSlice}(2,2)*READ + obj.InPlaneRotMatrix{cSlice}(2,3)*SLICE;
                TRA = obj.InPlaneRotMatrix{cSlice}(3,1)*PHASE + obj.InPlaneRotMatrix{cSlice}(3,2)*READ + obj.InPlaneRotMatrix{cSlice}(3,3)*SLICE;

                %% Compute rotation matrix to align the logical coordinate normal vector with normal vector in the patient coordiate system [SAG, COR, TRA]   
                if(norm(obj.Normal{cSlice})>1.001 || norm(obj.Normal{cSlice})<0.999 )
                   error('SODA.Normal must be unit vector') 
                end

                v = cross(initNormal,obj.Normal{cSlice});
                s = norm(v);     % sine of angle
                c = dot(initNormal,obj.Normal{cSlice});   % cosine of angle

                if(s == 0)
                   obj.RotMatrix{cSlice} = eye(3)*c;
                else
                   V = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
                   obj.RotMatrix{cSlice} = eye(3)  + V + V*V*(1-c)/s^2;
                end

                %% Apply rotation matrix, add shift ans set READ as first dimension
                obj.Coords{cSlice}.SAG = permute(obj.RotMatrix{cSlice}(1,1)*SAG + obj.RotMatrix{cSlice}(1,2)*COR + obj.RotMatrix{cSlice}(1,3)*TRA + obj.Position{cSlice}(1),[2 1 3]);
                obj.Coords{cSlice}.COR = permute(obj.RotMatrix{cSlice}(2,1)*SAG + obj.RotMatrix{cSlice}(2,2)*COR + obj.RotMatrix{cSlice}(2,3)*TRA + obj.Position{cSlice}(2),[2 1 3]);
                obj.Coords{cSlice}.TRA = permute(obj.RotMatrix{cSlice}(3,1)*SAG + obj.RotMatrix{cSlice}(3,2)*COR + obj.RotMatrix{cSlice}(3,3)*TRA + obj.Position{cSlice}(3),[2 1 3]);

                obj.Coords{cSlice}.SAG = obj.Coords{cSlice}.SAG(end:-1:1,:,:);
                obj.Coords{cSlice}.COR = obj.Coords{cSlice}.COR(end:-1:1,:,:);
                obj.Coords{cSlice}.TRA = obj.Coords{cSlice}.TRA(end:-1:1,:,:);
                
            end
            
            
            
        end
        
        function obj = calcPixLocInterp(obj)
            
            % Check
            if(obj.NPixelReadoutInterp==0 || obj.NPixelPhaseInterp==0 || obj.NPixelSliceInterp==0)
                error('NPixelInterp have to be set with setCSIVoxelsize() before calling calcPixLocInterp()')
            end
            
            % Calculate new pixel sizes
            obj.PixelSizeReadoutInterp = obj.ReadoutFOV/obj.NPixelReadoutInterp;
            obj.PixelSizePhaseInterp   = obj.PhaseFOV/obj.NPixelPhaseInterp;
            obj.PixelSizeSliceInterp   = obj.Thickness/obj.NPixelSliceInterp;
            
            % Determine dominant obj.MainOrientation
            [~,obj.MainOrientation] = max(abs(obj.Normal{1}));
            
            %% Loop over slices
            for cSlice =1:obj.NSlices
                %% Calculate pixel positions for slice in logical coordinate system [PHASE, READ, SLICE]
                if(obj.Dim==2) % 2D
                         [PHASE, READ, SLICE] = ndgrid(   linspace(   -obj.PhaseFOV/2, obj.PhaseFOV/2   -obj.PhaseFOV/  obj.NPixelPhase,        obj.NPixelPhaseInterp       ) ,... + obj.PixelSizePhase/2 
                                                    linspace(    obj.ReadoutFOV/2,-obj.ReadoutFOV/2 +obj.ReadoutFOV/obj.NPixelReadout,    obj.NPixelReadoutInterp     ) ,...+ obj.PixelSizeReadout/2
                                                    linspace(   0, 0,        1      )    ); 
                    
                    
               elseif(obj.Dim==3) % 3D / For CSI the the center of the voxel defines the position
                   [PHASE, READ, SLICE] = ndgrid(   linspace(   -obj.PhaseFOV/2, obj.PhaseFOV/2   -obj.PhaseFOV/  obj.NPixelPhase,        obj.NPixelPhaseInterp       ) ,... + obj.PixelSizePhase/2 
                                                    linspace(    obj.ReadoutFOV/2,-obj.ReadoutFOV/2 +obj.ReadoutFOV/obj.NPixelReadout,    obj.NPixelReadoutInterp     ) ,...+ obj.PixelSizeReadout/2
                                                    linspace(   -obj.Thickness/2, obj.Thickness/2-obj.Thickness/obj.NPixelSlice,        obj.NPixelSliceInterp       )    ); %+ obj.PixelSizeSlice/2   
                   
                end
                
                %% Transformation from logical to patient coordinate system           
                if(obj.MainOrientation==1) % SAG
                    obj.InPlaneRotMatrix = [  0                    0                     1; ... 
                                              cos(obj.InplaneRot{1})  sin(obj.InplaneRot{1})   0;...
                                             -sin(obj.InplaneRot{1})  cos(obj.InplaneRot{1})   0];
                    initNormal = [1; 0; 0];  
                elseif(obj.MainOrientation==2) %COR
                    obj.InPlaneRotMatrix = [  cos(obj.InplaneRot{1})  sin(obj.InplaneRot{1})   0;...
                                              0                    0                     1; ... 
                                              sin(obj.InplaneRot{1})  -cos(obj.InplaneRot{1})  0];
                    initNormal = [0; 1; 0];  
                elseif(obj.MainOrientation==3) %TRA
                    obj.InPlaneRotMatrix = [  sin(obj.InplaneRot{1})  -cos(obj.InplaneRot{1})  0;...
                                              cos(obj.InplaneRot{1})  sin(obj.InplaneRot{1})   0; ... 
                                              0                    0                     1];
                    initNormal = [0; 0; 1];                       
                end
                
                SAG = obj.InPlaneRotMatrix(1,1)*PHASE + obj.InPlaneRotMatrix(1,2)*READ + obj.InPlaneRotMatrix(1,3)*SLICE;
                COR = obj.InPlaneRotMatrix(2,1)*PHASE + obj.InPlaneRotMatrix(2,2)*READ + obj.InPlaneRotMatrix(2,3)*SLICE;
                TRA = obj.InPlaneRotMatrix(3,1)*PHASE + obj.InPlaneRotMatrix(3,2)*READ + obj.InPlaneRotMatrix(3,3)*SLICE;

                %% Compute rotation matrix to align the logical coordinate normal vector with normal vector in the patient coordiate system [SAG, COR, TRA]   
                if(norm(obj.Normal{1})>1.001 || norm(obj.Normal{1})<0.999 )
                   error('SODA.Normal must be unit vector') 
                end

                v = cross(initNormal,obj.Normal{1}); %axis of rotation
                s = norm(v);     % sine of angle
                c = dot(initNormal,obj.Normal{1});   % cosine of angle

                if(s == 0)
                   obj.RotMatrix = eye(3)*c;
                else
                   V = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
                   obj.RotMatrix = eye(3)  + V + V*V*(1-c)/s^2;
                   %https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
                end

                %% Apply rotation matrix, add shift and set READ as first dimension
                obj.Coords{cSlice}.SAG = permute(obj.RotMatrix(1,1)*SAG + obj.RotMatrix(1,2)*COR + obj.RotMatrix(1,3)*TRA + obj.Position{cSlice}(1),[2 1 3]);
                obj.Coords{cSlice}.COR = permute(obj.RotMatrix(2,1)*SAG + obj.RotMatrix(2,2)*COR + obj.RotMatrix(2,3)*TRA + obj.Position{cSlice}(2),[2 1 3]);
                obj.Coords{cSlice}.TRA = permute(obj.RotMatrix(3,1)*SAG + obj.RotMatrix(3,2)*COR + obj.RotMatrix(3,3)*TRA + obj.Position{cSlice}(3),[2 1 3]);

                obj.Coords{cSlice}.SAG = obj.Coords{cSlice}.SAG(end:-1:1,:,:);
                obj.Coords{cSlice}.COR = obj.Coords{cSlice}.COR(end:-1:1,:,:);
                obj.Coords{cSlice}.TRA = obj.Coords{cSlice}.TRA(end:-1:1,:,:);
                
            end
            
            
            
        end
        
        function [S,C,T] = applyRotations(obj,PHASE,READ,SLICE)
            
            %% Transformation from logical to patient coordinate system           
            SAG = obj.InPlaneRotMatrix{1}(1,1)*PHASE + obj.InPlaneRotMatrix{1}(1,2)*READ + obj.InPlaneRotMatrix{1}(1,3)*SLICE;
            COR = obj.InPlaneRotMatrix{1}(2,1)*PHASE + obj.InPlaneRotMatrix{1}(2,2)*READ + obj.InPlaneRotMatrix{1}(2,3)*SLICE;
            TRA = obj.InPlaneRotMatrix{1}(3,1)*PHASE + obj.InPlaneRotMatrix{1}(3,2)*READ + obj.InPlaneRotMatrix{1}(3,3)*SLICE;
            
            %% Apply rotation matrix
            S = obj.RotMatrix{1}(1,1)*SAG + obj.RotMatrix{1}(1,2)*COR + obj.RotMatrix{1}(1,3)*TRA;
            C = obj.RotMatrix{1}(2,1)*SAG + obj.RotMatrix{1}(2,2)*COR + obj.RotMatrix{1}(2,3)*TRA;
            T = obj.RotMatrix{1}(3,1)*SAG + obj.RotMatrix{1}(3,2)*COR + obj.RotMatrix{1}(3,3)*TRA;

        end
        
        function [PHASE,READ,SLICE] = unapplyRotations(obj,S,C,T)
            %warning('Using first rot matrix for first slice')
            %% Unapply rotation matrix
            InvMat = inv(obj.RotMatrix{1});
            
            SAG = InvMat(1,1)*S + InvMat(1,2)*C + InvMat(1,3)*T;
            COR = InvMat(2,1)*S + InvMat(2,2)*C + InvMat(2,3)*T;
            TRA = InvMat(3,1)*S + InvMat(3,2)*C + InvMat(3,3)*T;
                        
            %% Transformation from logical to patient coordinate system  
            InvInplaneMat = inv(obj.InPlaneRotMatrix{1});
            
            PHASE = InvInplaneMat(1,1)*SAG + InvInplaneMat(1,2)*COR + InvInplaneMat(1,3)*TRA;
            READ  = InvInplaneMat(2,1)*SAG + InvInplaneMat(2,2)*COR + InvInplaneMat(2,3)*TRA;
            SLICE = InvInplaneMat(3,1)*SAG + InvInplaneMat(3,2)*COR + InvInplaneMat(3,3)*TRA;
            
        end
    end
end

