classdef BilinearInterpolator   
   %BILINEARINTERPOLATOR Bi-linear interpolation of corrections.
   % Executes a second-order interpolation from level L to L+1.
   
   %======================== METHODS =================================
   methods
      function u = interpolate(obj, coarseLevel, fineLevel, uc)
         % Restrict the fine-level function F at level FINELEVEL
         % to the function FC at level COARSELEVEL.
         
         % Interpolate along one dimension at a time
         u = obj.interpInX(coarseLevel, fineLevel, uc);
      end
   end

   %======================== PRIVATE METHODS =========================
   methods (Access = private)
      function u = interpInX(obj, coarseLevel, fineLevel, uc)
         % Linear interpolation in x
         % Aliases, allocate output array
         nf = fineLevel.n;
         nc = coarseLevel.n;
         u = zeros(nf,1);
         % Inject coarse points into the respective fine points
         u(1:2:nf(1)) = uc;
         % Linearly interpolate into in-between fine-level points
         for i1=1:nc(1)-1
            u(2*i1) = 0.5*(uc(i1) + uc(i1+1));
         end
      end
      
   end
   
end

