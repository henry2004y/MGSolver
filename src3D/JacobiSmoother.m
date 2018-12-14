classdef (Sealed) JacobiSmoother < Smoother
   %JACOBISMOOTHER Gauss-Seidel relaxation scheme.
   % This class executes weighted Jacobi relaxation.
   % It can be applied at any level.
   
   %======================== MEMBERS =================================
   properties (GetAccess = private, SetAccess = private)
      w  {double} % weight
   end
   %======================== CONSTRUCTORS ============================
   methods
      function obj = JacobiSmoother(w)
         obj@Smoother('Jacobi');
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
         i = [1 level.n(1)]; uNew(i,:,:) = f(i,:,:);
         j = [1 level.n(2)]; uNew(:,j,:) = f(:,j,:);
         k = [1 level.n(3)]; uNew(:,:,k) = f(:,:,k);
         % Relax in the internal domain
         for k=2:level.n(3)-1
            for j=2:level.n(2)-1
               for i=2:level.n(1)-1
                  uNew(i,j,k) = obj.w*0.25*(h2*f(i,j,k) ...
                     + u(i-1,j,k)/3 + u(i+1,j,k)/3 ...
                     + u(i,j-1,k)/3 + u(i,j+1,k)/3 ...
                     + u(i,j,k-1)/3 + u(i,j,k+1)/3 ...
                     + u(i-1,j-1,k)/6 + u(i+1,j-1,k)/6 ...
                     + u(i-1,j+1,k)/6 + u(i+1,j+1,k)/6 ...
                     + u(i,j-1,k-1)/6 + u(i,j-1,k+1)/6 ...
                     + u(i,j+1,k-1)/6 + u(i,j+1,k+1)/6 ...
                     + u(i-1,j,k-1)/6 + u(i+1,j,k-1)/6 ...
                     + u(i-1,j,k+1)/6 + u(i+1,j,k+1)/6 )...
                     + (1 - obj.w)*u(i,j,k);
               end
            end
         end
         % A relaxation sweep is counted as one flop per internal
         % gridpoint
         addflops(prod(level.n-1));
      end
      
      function u = relax_conv(obj, level, u)
         % Weighted Jacobi successive displacement using convolution.
         % This is faster than the explicit index-based method on a large
         % grid, but slower on a small grid. An optimal choice may be a
         % combination of these two. My current test shows that for cell
         % numbers less than 64, I should use the brute-force method.
         
         % Stencil
         L = zeros(3,3,3);
         L(:,:,1) = [ 0 -1 0; -1 -2 -1; 0 -1 0];
         L(:,:,2) = [-1 -2 -1; -2 24 -2; -1 -2 -1];
         L(:,:,3) = [ 0 -1 0; -1 -2 -1; 0 -1 0];
         L = (1/6)*L;
         
         % Useful aliases
         h2 = level.h^2;
         f = level.f;
         % Relax in the internal domain
         L_sm = L;
         L_sm(2,2,2) = 0;
         u = obj.w/L(2,2,2)*(h2*f - convn(u,L_sm,'same')) + (1 - obj.w)*u;
         
         % Impose B.C.
         i = [1 level.n(1)]; u(i,:,:) = f(i,:,:);
         j = [1 level.n(2)]; u(:,j,:) = f(:,j,:);
         k = [1 level.n(3)]; u(:,:,k) = f(:,:,k);
         
         % A relaxation sweep is counted as one flop per internal
         % gridpoint
         addflops(prod(level.n-1));
      end
   end
   
end