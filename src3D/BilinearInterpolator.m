classdef BilinearInterpolator   
   %BILINEARINTERPOLATOR Bi-linear interpolation of corrections.
   % Executes a second-order interpolation from level L to L+1.
   
   %======================== METHODS =================================
   methods
      function u = interpolate(obj, coarseLevel, fineLevel, uc)
         % Restrict the fine-level function F at level FINELEVEL
         % to the function FC at level COARSELEVEL.
         
         % Interpolate along one dimension at a time
         u1 = obj.interpInX(coarseLevel, fineLevel, uc);
         u2 = obj.interpInY(coarseLevel, fineLevel, u1);
         u  = obj.interpInZ(coarseLevel, fineLevel, u2);
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
         for i=1:nc(1)-1
            u(2*i,:,:) = 0.5*(uc(i,:,:) + uc(i+1,:,:));
         end
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
         for j=1:nc(2)-1
            u(:,2*j,:) = 0.5*(uc(:,j,:) + uc(:,j+1,:));
         end
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
         for k=1:nc(3)-1
            u(:,:,2*k) = 0.5*(uc(:,:,k) + uc(:,:,k+1));
         end
      end
   end
   
end

