classdef MTI
    % An operator for correcting B0 inhomogenity using Multi time
    % interpolation for spiral imaging.
    %
    %
    %
    %Usage:
    % MTIOP=MTI(b0Map,tk,FT,CoilSens)
    % im_B0corr=MFI*CoilData;
    % CoilData should be of the form [nCh,nFE,nSlc,nRep]=size(CoilData);
    %
    %
    %Dependencies:
    % NUFFT operator from MIRT toolbox(by fressler) and
    % NUFFT operator from ESPIRIT toolbox (by LUSTIG)
    %
    %
    %Reference:
    % Sutton, B. P., Noll, D. C., & Fessler, J. A. (2003).
    % Fast, iterative image reconstruction for MRI in the presence of field inhomogeneities.
    % IEEE Transactions on Medical Imaging, 22(2), 178–188.
    % https://doi.org/10.1109/TMI.2002.808360
    
    properties
        MTI_weights %Coeffients for interpolation
        tflag %transpose flag
        B0map %in rad/s
        NUFFTOP %Gridder from FRESSLER/LUSTIG toolbox
        tk % %time points in seconds
        CoilSens
        mode % Interpolation mode {'LeastSquares','NearestNeighbour'}
        nLevels % Number of frequency levels for MFI
        tau
        
        flags
    end
    methods
        function obj=MTI(varargin)
            if(nargin==4)
                obj.B0map=varargin{1};
                obj.tk=varargin{2};
                obj.NUFFTOP=varargin{3};
                obj.CoilSens=varargin{4};
                obj.tflag=0;
                
                %calculate some parameters
                obj.nLevels=ceil(2*(max(abs(obj.B0map(:)))*max(obj.tk)/pi));
                obj.tau=max(obj.tk)/(obj.nLevels+1);
                obj.mode='LeastSquares'; % {'LeastSquares','NearestNeighbour'}
                
                % calculate weights
                obj=obj.CalcWeights();
                
                %set initial flags
                flags.DCFApplied=false;
                flags.CoilCombine='adapt2';
                
            else
                error('Need several input parameters')
            end
        end
        
        function out=mtimes(obj,InData)
            if(obj.tflag==1) % Inverse operator: kspace data to image
                %out: is corrected image
                %InData: is kspace data

                [nCh,nFE]=size(InData);
                nIntlv=nFE/length(obj.tk);
                
                %modulate the signal with MTI weights
                sig_MTI=bsxfun(@times,conj(obj.MTI_weights),reshape(InData.',1,nFE/nIntlv,[]));
                sig_MTI=reshape(sig_MTI,obj.nLevels+1,nFE,nCh);
                sig_MTI=reshape(permute(sig_MTI,[2 3 1]),nFE,nCh*(obj.nLevels+1));
                
                %Inverse NUFFT on the modulated signal
                temp=  obj.NUFFTOP.st.p'*double(sig_MTI);
                temp=reshape(temp,obj.NUFFTOP.st.Kd(1),obj.NUFFTOP.st.Kd(2),[]);
                temp=ifft2(temp);
                all_images=temp(1:obj.NUFFTOP.imSize(1),1:obj.NUFFTOP.imSize(2),:);
                
                %apply NUFFT scaling factor and add the phase shifts due BO
                fac=bsxfun(@times,reshape(exp(-1i*obj.B0map(:)*(obj.tau*(0:obj.nLevels))),obj.NUFFTOP.imSize(1),obj.NUFFTOP.imSize(2),[]),conj(obj.NUFFTOP.st.sn));
                all_images=bsxfun(@times,reshape(all_images,[obj.NUFFTOP.imSize, nCh,(obj.nLevels+1)]) , permute(fac,[1,2, 4, 3]));
                all_images=sum(all_images,4);
                
                % Combine coils and combine MTI components
                if((isscalar(obj.CoilSens)||isempty(obj.CoilSens)))
                    out=permute(all_images,[3,1,2]);
                    out=reshape(out,nCh,[]);
                else %try to combine coils
                    all_images=all_images.*permute(obj.CoilSens,[2,3,1]);
                    out=sum(all_images,3);
                end
                
            else % forward operation coil images to kspace
            [nCh,~,~]=size(InData);
            InData=reshape(InData,[nCh,obj.NUFFTOP.imSize]);
            nIntlv=size(obj.NUFFTOP.st.p,1)/length(obj.tk);
%             ksp_MTI=zeros(size(obj.NUFFTOP.st.p,1),obj.nLevels+1);
%             sig=zeros(nCh,size(obj.NUFFTOP.st.p,1),nSlc,nRep);
                
             fac=bsxfun(@times,reshape(exp(1i*obj.B0map(:)*(obj.tau*(0:obj.nLevels))),...
                 obj.NUFFTOP.imSize(1),obj.NUFFTOP.imSize(2),[]),obj.NUFFTOP.st.sn);
            temp=bsxfun(@times,permute(InData,[2, 3, 1]),permute(fac,[1, 2, 4, 3]));
            temp=fft2(temp,obj.NUFFTOP.st.Kd(1),obj.NUFFTOP.st.Kd(2));
            sig=obj.NUFFTOP.st.p*double(reshape(temp,prod(obj.NUFFTOP.st.Kd),[]));
            sig=reshape(sig,length(obj.tk),nIntlv,nCh,[]);
            sig=bsxfun(@times,sig,permute(obj.MTI_weights,[2 3 4 1]));
            out=sum(sig,4);
            out=reshape(out,[],nCh).';
            end
        end
        function obj=ctranspose(obj)
            if(obj.tflag==0)
                obj.tflag=1;
            end
        end
        
        function obj=CalcWeights(obj)
            wn_tau_l=exp(1i*(obj.B0map(:)*(obj.tau.* (0:obj.nLevels))));
            if(strcmp(obj.mode,'LeastSquares'))
                actual=exp(1i*col(obj.B0map)*obj.tk');
                obj.MTI_weights=mldivide(wn_tau_l'*wn_tau_l,wn_tau_l'*actual);
            elseif(strcmp(obj.mode,'NearestNeighbour') )
                error('Not implemented')
            end
            
        end
    end
    
    
end

%% old stuff
%         function obj=CalcWeights(obj)
% %             [W,ai,bi]=unique(obj.B0map(:));
% %             wn_tau_l=exp(-1i*obj.tau.*(W(:)*(0:obj.nLevels)));
% %             actual=exp(1i*W(:)*obj.tk');
%             wn_tau_l=exp(-1i*(obj.B0map(:)*(obj.tau.* (0:obj.nLevels))));
% %             obj.MTI_weights=zeros([numel(obj.B0map) obj.nLevels]);
%
%             % calculate weights for all values in b0maps
%             if(strcmp(obj.mode,'LeastSquares'))
%
%                 %                 for i=1:size(obj.MTI_weights,1)
%                 %                     for j=1:size(obj.MTI_weights,2)
%                 %                         B=exp(-1i*obj.B0map(i,j)*obj.tk);
% %                 actual=exp(1i*obj.B0map(:)*obj.tk');
%                 %                         obj.MTI_weights= mldivide(wn_tau_l,actual);
%                 %much mfater
% %                 obj.MTI_weights= mldivide(wn_tau_l'*wn_tau_l,wn_tau_l'*actual);
%                  actual=exp(-1i*col(obj.B0map)*obj.tk');
%                  obj.MTI_weights=mldivide(wn_tau_l'*wn_tau_l,wn_tau_l'*actual);
%                 %                     end
%                 %                 end
%
%                 warning('on','MATLAB:nearlySingularMatrix');
%             elseif(strcmp(obj.mode,'NearestNeighbour') )
%                 error('Not implemented')
%                 %                 for i=1:size(obj.MFI_weights,1)
%                 %                     for j=1:size(obj.MFI_weights,2)
%                 %                         [~,I]=min(abs(obj.wi-obj.B0map(i,j)));
%                 %                         obj.MFI_weights(i,j,I)=1;
%                 %                     end
%                 %                 end
%             end
%
%         end
%     end