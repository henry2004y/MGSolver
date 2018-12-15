classdef (Sealed) GaussSeidelSmoother < Smoother
   %GAUSSSEIDELSMOOTHER Gauss-Seidel relaxation scheme.
   % This class executes Gauss-Seidel relaxation sweeps in lexicographic
   % order. It can be applied at any level.
   
   %======================== MEMBERS =================================
   properties (GetAccess = private, SetAccess = private)
      numColors  double {mustBeInteger} % Number of colors (1=LEX, 2=RB)
   end
   %======================== CONSTRUCTORS ============================
   methods
      function obj = GaussSeidelSmoother(numColors)
         obj@Smoother('GS');
         obj.numColors = numColors;
      end
   end
   %======================== METHODS =================================
   methods
      function u = relax(obj, level, u)
         % Gauss-Seidel successive displacement in lexicographic ordering.
         % Because MATLAB passes array parameters by value, this does not
         % override the original U array.
         
         % Useful aliases
         h2 = level.h^2;
         f = level.f;
         % Impose B.C.
         i = [1 level.n(1)]; u(i,:,:) = f(i,:,:);
         j = [1 level.n(2)]; u(:,j,:) = f(:,j,:);
         k = [1 level.n(3)]; u(:,:,k) = f(:,:,k);
         % Relax in the internal domain
         for colorSweep=0:obj.numColors-1
            for k=2:level.n(3)-1
               for j=2:level.n(2)-1
                  for i=2:level.n(1)-1
                     if mod(i+j+k, obj.numColors) == colorSweep
                        u(i,j,k) = 0.25*(h2*f(i,j,k) ...
                           + u(i-1,j,k)/3 + u(i+1,j,k)/3 ...
                           + u(i,j-1,k)/3 + u(i,j+1,k)/3 ...
                           + u(i,j,k-1)/3 + u(i,j,k+1)/3 ...
                           + u(i-1,j-1,k)/6 + u(i+1,j-1,k)/6 ...
                           + u(i-1,j+1,k)/6 + u(i+1,j+1,k)/6 ...
                           + u(i,j-1,k-1)/6 + u(i,j-1,k+1)/6 ...
                           + u(i,j+1,k-1)/6 + u(i,j+1,k+1)/6 ...
                           + u(i-1,j,k-1)/6 + u(i+1,j,k-1)/6 ...
                           + u(i-1,j,k+1)/6 + u(i+1,j,k+1)/6 );
                     end
                  end
               end
            end
            % A relaxation sweep is counted as one flop per internal
            % gridpoint
            addflops(prod(level.n-1));
         end
      end
   end
   
end
