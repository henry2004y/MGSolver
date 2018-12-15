classdef (Sealed) Operator < handle
   %OPERATOR Discrete operator computer.
   % This class computes the discrete operator L(U) of a function U at
   % a certain level in the multi-level algorithm.
   
   %======================== MEMBERS =================================
   properties (GetAccess = private, SetAccess = private)
      level % A data structure holding convenient level variables
   end
   %======================== CONSTRUCTORS ============================
   methods
      function obj = Operator(level)
         % Initializes an operator computer at level LEVEL.
         obj.level = level;
      end
   end
   %======================== METHODS =================================
   methods
      function result = L(obj, u)
         % Apply the discrete operator L to a function U. This is the
         % 7-point Laplacian with Dirichlet boundary conditions.
         
         % Brute-force way, equivalent to the convolution method
         % You can choose to do this or use convolution instead.
%          % Allocate output array
%          result = zeros(obj.level.n);
%          % Set Dirichlet boundary conditions
%          i = [1 obj.level.n(1)]; result(i,:,:) = u(i,:,:);
%          j = [1 obj.level.n(2)]; result(:,j,:) = u(:,j,:);
%          k = [1 obj.level.n(3)]; result(:,:,k) = u(:,:,k);
%          % 19-point discrete Laplacian in the interior domain
%          rh2 = 1/obj.level.h^2;
%          for k=2:obj.level.n(3)-1
%             for j=2:obj.level.n(2)-1
%                for i=2:obj.level.n(1)-1
%                   result(i,j,k) = rh2/6*(...
%                      24*u(i,j,k) ...
%                      - 2*u(i-1,j,k) - 2*u(i+1,j,k) ...
%                      - 2*u(i,j-1,k) - 2*u(i,j+1,k) ...
%                      - 2*u(i,j,k-1) - 2*u(i,j,k+1) ...
%                      - u(i-1,j-1,k) - u(i+1,j-1,k) ...
%                      - u(i-1,j+1,k) - u(i+1,j+1,k) ...
%                      - u(i,j-1,k-1) - u(i,j-1,k+1) ...
%                      - u(i,j+1,k-1) - u(i,j+1,k+1) ...
%                      - u(i-1,j,k-1) - u(i+1,j,k-1) ...
%                      - u(i-1,j,k+1) - u(i+1,j,k+1) );
%                end
%             end
%          end

         % Stencil
         L = zeros(3,3,3);
         L(:,:,1) = [ 0 -1 0; -1 -2 -1; 0 -1 0];
         L(:,:,2) = [-1 -2 -1; -2 24 -2; -1 -2 -1];
         L(:,:,3) = [ 0 -1 0; -1 -2 -1; 0 -1 0];
         L = (1/6)*L;
         rh2 = 1/obj.level.h^2;
         result = rh2*convn(u,L,'same');
         % Set Dirichlet boundary conditions
         i = [1 obj.level.n(1)]; result(i,:,:) = u(i,:,:);
         j = [1 obj.level.n(2)]; result(:,j,:) = u(:,j,:);
         k = [1 obj.level.n(3)]; result(:,:,k) = u(:,:,k);
      end
   end
            
end