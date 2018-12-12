classdef FwLinearRestrictor
   %FULLWEIGHTINGRESTRICTOR Second-order full weighting of residuals.
   % Executes a second-order full weighting from level L+1 to L.
   
   %======================== METHODS =================================
   methods
      function fc = restrict(obj, coarseLevel, fineLevel, f)
         % Restrict the fine-level function F at level FINELEVEL
         % to level COARSELEVEL.
         
         % Interpolate along one dimension at a time
         fc = obj.restrictInX(coarseLevel, fineLevel, f);
      end
   end
   
   %======================== PRIVATE METHODS =========================
   methods (Access = private)
      function fc = restrictInX(obj, coarseLevel, fineLevel, f)
         % Full-weighting in x
         
         % Aliases, allocate output array
         nf = fineLevel.n;
         nc = coarseLevel.n;
         fc = zeros(nc,1);
         % Full-weighting of boundary residuals
         fc([1 nc]) = f([1 nf]);
         % Full-weighting of interior residuals
         for i=2:nc-1
            fc(i) = 0.25*(f(2*i-2) + 2*f(2*i-1) + f(2*i));
         end
      end
      
   end
   
end