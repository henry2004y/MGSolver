classdef (Sealed) JacobiSmoother < Smoother
   %JACOBISMOOTHER Gauss-Seidel relaxation scheme.
   % This class executes weighted Jacobi relaxation.
   % It can be applied at any level.
   
   %======================== MEMBERS =================================
   properties (GetAccess = private, SetAccess = private)
      w % weight
   end
   %======================== CONSTRUCTORS ============================
   methods
      function obj = JacobiSmoother(w)
         obj.w = w;
      end
   end
   %======================== METHODS =================================
   methods
      function uNew = relax(obj, level, u)
         % Weighted Jacobi successive displacement.
         
         % Initialization
         uNew = Inf(size(u));
         % Useful aliases
         h2 = level.h^2;
         f = level.f;
         % Impose B.C.
         i = [1 level.n(1)]; uNew(i,:) = f(i,:);
         j = [1 level.n(2)]; uNew(:,j) = f(:,j);
         % Relax in the internal domain
         for j=2:level.n(2)-1
            for i=2:level.n(1)-1
               uNew(i,j) = obj.w*0.25*(h2*f(i,j) ...
                  + u(i  ,j-1) + u(i  ,j+1) ...
                  + u(i-1,j )  + u(i+1,j  ) ) + (1 - obj.w)*u(i,j);
            end
         end
         % A relaxation sweep is counted as one flop per internal
         % gridpoint
         addflops(prod(level.n-1));
      end
   end
   
end