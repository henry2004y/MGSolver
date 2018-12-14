classdef FwLinearRestrictor
   %FULLWEIGHTINGRESTRICTOR Second-order full weighting of residuals.
   % Executes a second-order full weighting from level L+1 to L.
   
   %======================== MEMBERS =================================
   properties (GetAccess = private, SetAccess = private)
      method  % 1 for 1D line, 2 for 3D matrix
   end
   %======================== CONSTRUCTORS ============================
   methods
      function obj = FwLinearRestrictor(method)
         obj.method = method;
      end
   end
   
   %======================== METHODS =================================
   methods
      function fc = restrict(obj, coarseLevel, fineLevel, f)
         % Restrict the fine-level function F at level FINELEVEL
         % to level COARSELEVEL.
         
         if obj.method == 1
            % Interpolate along one dimension at a time
            f1 = obj.restrictInX(coarseLevel, fineLevel, f);
            f2 = obj.restrictInY(coarseLevel, fineLevel, f1);
            fc = obj.restrictInZ(coarseLevel, fineLevel, f2);
         elseif obj.method == 2
            fc = obj.restrict3D(coarseLevel, fineLevel, f);
         else
            error('unknown restriction method!')
         end
      end
   end
   
   %======================== PRIVATE METHODS =========================
   methods (Access = private)
      function fc = restrictInX(obj, coarseLevel, fineLevel, f)
         % Full-weighting in x
         % Aliases, allocate output array
         nf = fineLevel.n;
         nc = coarseLevel.n;
         fc = zeros(nc(1),nf(2),nf(3));
         % Full-weighting of boundary residuals
         fc([1 nc(1)],:,:) = f([1 nf(1)],:,:);
         % Full-weighting of interior residuals
         for i=2:nc(1)-1
            fc(i,:,:) = 0.25*(f(2*i-2,:,:) + 2*f(2*i-1,:,:) + f(2*i,:,:));
         end
      end
      
      function fc = restrictInY(obj, coarseLevel, fineLevel, f)
         % Full-weighting in y
         % Aliases, allocate output array
         nf = fineLevel.n;
         nc = coarseLevel.n;
         fc = zeros(nc(1),nc(2),nf(3));
         % Full-weighting of boundary residuals
         fc(:,[1 nc(2)],:) = f(:,[1 nf(2)],:);
         % Full-weighting of interior residuals
         for j=2:nc(2)-1
            fc(:,j,:) = 0.25*(f(:,2*j-2,:) + 2*f(:,2*j-1,:) + f(:,2*j,:));
         end
      end
      
      function fc = restrictInZ(obj, coarseLevel, fineLevel, f)
         % Full-weighting in z
         % Aliases, allocate output array
         nf = fineLevel.n;
         nc = coarseLevel.n;
         fc = zeros(nc);
         % Full-weighting of boundary residuals
         fc(:,:,[1 nc(3)]) = f(:,:,[1 nf(3)]);
         % Full-weighting of interior residuals
         for k=2:nc(3)-1
            fc(:,:,k) = 0.25*(f(:,:,2*k-2) + 2*f(:,:,2*k-1) + f(:,:,2*k));
         end
      end
      
      function fc = restrict3D(obj, coarseLevel, fineLevel, f)
         % Full-weighting in 3D as a whole
         % Define restriction operator
         KER = zeros(3,3,3);
         KER(:,:,1) = 1/28*[0 1 0;1 2 1;0 1 0];
         KER(:,:,2) = 1/28*[1 2 1;2 4 2;1 2 1];
         KER(:,:,3) = 1/28*[0 1 0;1 2 1;0 1 0];
         fc = convn(f,KER,'same');
         % Aliases
         nf = fineLevel.n;
         nc = coarseLevel.n;
         % Full-weighting of boundary residuals
         fc([1 nf(1)],:,:) = f([1 nf(1)],:,:);
         fc(:,[1 nf(2)],:) = f(:,[1 nf(2)],:);
         fc(:,:,[1 nf(3)]) = f(:,:,[1 nf(3)]);
         % Actual value for the coarse grid
         fc = fc(1:2:end,1:2:end,1:2:end);
      end
   end
   
end