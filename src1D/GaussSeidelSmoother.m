classdef (Sealed) GaussSeidelSmoother < Smoother
   %GAUSSSEIDALSMOOTHER Relaxation scheme.
   % This class executes Gauss-Seidel relaxation sweeps in lexicographic
   % order. It can be applied at any level.
   
   %======================== MEMBERS =================================
   properties (GetAccess = private, SetAccess = private)
      numColors % Number of colors (1=LEX, 2=RB)
   end
   %======================== CONSTRUCTORS ============================
   methods
      function obj = GaussSeidelSmoother(numColors)
         obj.numColors = numColors;
      end
   end
   %======================== METHODS =================================
   methods
      function u = relax(obj, level, u)
         % Gauss-Seidel successive displacement.
         % Because MATLAB passes array parameters by value, this does not
         % override the original U array.
         
         % Useful aliases
         h2 = level.h^2;
         f = level.f;
         % Impose B.C.
         i = [1 level.n(1)]; u(i) = f(i);
         % Relax in the internal domain
         for i=2:level.n(1)-1
            u(i) = 0.5*(h2*f(i) + u(i-1) + u(i+1));
         end
         % A relaxation sweep is counted as one flop per internal
         % gridpoint
         addflops(prod(level.n-1));
      end
      
   end
   
end