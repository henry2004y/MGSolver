classdef (Sealed) JacobiSmoother < Smoother
   %JACOBISMOOTHER Relaxation scheme.
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
         % Because MATLAB passes array parameters by value, this does not
         % override the original U array.
         
         % Initialization
         uNew = Inf(size(u));
         
         % Useful aliases
         h2 = level.h^2;
         f = level.f;
         
         % Impose B.C.
         i = [1 level.n]; uNew(i) = f(i);
         % Relax in the internal domain
         for i=2:level.n-1
            uNew(i) = obj.w*0.5*(h2*f(i) + u(i-1) + u(i+1)) + ...
               (1-obj.w)*u(i);
         end
         % A relaxation sweep is counted as one flop per internal
         % gridpoint
         addflops(prod(level.n-1));
      end
   end
   
end