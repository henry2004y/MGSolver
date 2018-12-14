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
         
         % Allocate output array
         result = zeros(obj.level.n);
         % Set Dirichlet boundary conditions
         i = [1 obj.level.n(1)]; result(i,:,:) = u(i,:,:);
         j = [1 obj.level.n(2)]; result(:,j,:) = u(:,j,:);
         k = [1 obj.level.n(3)]; result(:,:,k) = u(:,:,k);
         % 7-point discrete Laplacian in the interior domain
         rh2 = 1/obj.level.h^2;
         for k=2:obj.level.n(3)-1
            for j=2:obj.level.n(2)-1
               for i=2:obj.level.n(1)-1
                  result(i,j,k) = rh2/6*(...
                     24*u(i,j,k) ...
                     - 2*u(i-1,j,k) - 2*u(i+1,j,k) ...
                     - 2*u(i,j-1,k) - 2*u(i,j+1,k) ...
                     - 2*u(i,j,k-1) - 2*u(i,j,k+1) ...
                     - u(i-1,j-1,k) - u(i+1,j-1,k) ...
                     - u(i-1,j+1,k) - u(i+1,j+1,k) ...
                     - u(i,j-1,k-1) - u(i,j-1,k+1) ...
                     - u(i,j+1,k-1) - u(i,j+1,k+1) ...
                     - u(i-1,j,k-1) - u(i+1,j,k-1) ...
                     - u(i-1,j,k+1) - u(i+1,j,k+1) );
               end
            end
         end
      end
   end
            
end