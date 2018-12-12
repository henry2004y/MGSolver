classdef FwLinearRestrictor
   %FULLWEIGHTINGRESTRICTOR Second-order full weighting of residuals.
   % Executes a second-order full weighting from level L+1 to L.
   
   %======================== METHODS =================================
   methods
      function fc = restrict(obj, coarseLevel, fineLevel, f)
         % Restrict the fine-level function F at level FINELEVEL
         % to level COARSELEVEL.
         
         % Interpolate along one dimension at a time
         f1 = obj.restrictInX(coarseLevel, fineLevel, f);
         fc = obj.restrictInY(coarseLevel, fineLevel, f1);
      end
   end
   
   %======================== PRIVATE METHODS =========================
   methods (Access = private)
      function fc = restrictInX(obj, coarseLevel, fineLevel, f)
         % Full-weighting in x
         % Aliases, allocate output array
         nf = fineLevel.n;
         nc = coarseLevel.n;
         fc = zeros(nc(1),nf(2));
         % Full-weighting of boundary residuals
         fc([1 nc(1)],:) = f([1 nf(1)],:);
         % Full-weighting of interior residuals
         for i1=2:nc(1)-1
            fc(i1,:) = 0.25*(f(2*i1-2,:) + 2*f(2*i1-1,:) + f(2*i1,:));
         end
      end
      
      function fc = restrictInY(obj, coarseLevel, fineLevel, f)
         % Full-weighting in y
         % Aliases, allocate output array
         nf = fineLevel.n;
         nc = coarseLevel.n;
         fc = zeros(nc);
         % Full-weighting of boundary residuals
         fc(:,[1 nc(2)]) = f(:,[1 nf(2)]);
         % Full-weighting of interior residuals
         for i2=2:nc(2)-1
            fc(:,i2) = 0.25*(f(:,2*i2-2) + 2*f(:,2*i2-1) + f(:,2*i2));
         end
      end
   end
   
end