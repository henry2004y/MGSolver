classdef BilinearInterpolator
   %BILINEARINTERPOLATOR Bi-linear interpolation of corrections.
   % Executes a second-order interpolation from level L to L+1.
   
   %======================== MEMBERS =================================
   properties (GetAccess = private, SetAccess = private)
      method  % 1 for 1D line, 2 for 3D matrix
   end
   %======================== CONSTRUCTORS ============================
   methods
      function obj = BilinearInterpolator(method)
         obj.method = method;
      end
   end
   
   %======================== METHODS =================================
   methods
      function u = interpolate(obj, coarseLevel, fineLevel, uc)
         % Restrict the fine-level function F at level FINELEVEL
         % to the function FC at level COARSELEVEL.
         
         if obj.method == 1
            % Interpolate along one dimension at a time
            u1 = obj.interpInX(coarseLevel, fineLevel, uc);
            u2 = obj.interpInY(coarseLevel, fineLevel, u1);
            u  = obj.interpInZ(coarseLevel, fineLevel, u2);
         elseif obj.method == 2
            % It is also possible to do interpolation all at once using
            % convn or interp3.
            u = obj.interp3D(coarseLevel, fineLevel, uc);
         else
            error('unknown interpolation method!')
         end
      end
   end
   
   %======================== PRIVATE METHODS =========================
   methods (Access = private)
      function u = interpInX(obj, coarseLevel, fineLevel, uc)
         % Linear interpolation in x
         % Aliases, allocate output array
         nf = fineLevel.n;
         nc = coarseLevel.n;
         u = zeros(nf(1),nc(2),nc(3));
         % Inject coarse points into the respective fine points
         u(1:2:nf(1),:,:) = uc;
         % Linearly interpolate into in-between fine-level points
         i_ = 1:nc(1)-1;
         u(2*i_,:,:) = 0.5*(uc(i_,:,:) + uc(i_+1,:,:));
      end
      
      function u = interpInY(obj, coarseLevel, fineLevel, uc)
         % Linear interpolation in y
         % Aliases, allocate output array
         nf = fineLevel.n;
         nc = coarseLevel.n;
         u = zeros(nf(1),nf(2),nc(3));
         % Inject coarse points into the respective fine points
         u(:,1:2:nf(2),:) = uc;
         % Linearly interpolate into in-between fine-level points
         j_ = 2:nc(2)-1;
         u(:,2*j_,:) = 0.5*(uc(:,j_,:) + uc(:,j_+1,:));
      end
      
      function u = interpInZ(obj, coarseLevel, fineLevel, uc)
         % Linear interpolation in z
         % Aliases, allocate output array
         nf = fineLevel.n;
         nc = coarseLevel.n;
         u = zeros(nf);
         % Inject coarse points into the respective fine points
         u(:,:,1:2:nf(3)) = uc;
         % Linearly interpolate into in-between fine-level points
         k_ = 1:nc(3)-1;
         u(:,:,2*k_) = 0.5*(uc(:,:,k_) + uc(:,:,k_+1));
      end
      
      function u = interp3D(obj, coarseLevel, fineLevel, uc)
         % Linear interpolation in 3D as a whole
         % Alias, allocate output array
         nf = fineLevel.n;
         u = zeros(nf);
         u(1:2:end ,1:2:end ,1:2:end ) = uc;
         % Define interpolation operator
         KER = zeros(3,3);
         KER(:,:,1) = 1/8*[1 2 1;2 4 2;1 2 1];
         KER(:,:,2) = 1/8*[2 4 2;4 8 4;2 4 2];
         KER(:,:,3) = 1/8*[1 2 1;2 4 2;1 2 1];
         u = convn(u,KER,'same');
      end
   end
   
end

