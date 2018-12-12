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
         % 5-point Laplacian with Dirichlet boundary conditions.
         
         % Allocate output array
         result = zeros(obj.level.n);
         % Set Dirichlet boundary conditions
         i1 = [1 obj.level.n(1)]; result(i1,:) = u(i1,:);
         i2 = [1 obj.level.n(2)]; result(:,i2) = u(:,i2);
         % 5-point discrete Laplacian in the interior domain
         rh2 = 1/obj.level.h^2;
         for i1 = 2:obj.level.n(1)-1
            for i2 = 2:obj.level.n(2)-1
               result(i1,i2) = rh2*(...
                  4*u(i1,i2) ...
                  -u(i1 ,i2-1) - u(i1 ,i2+1) ...
                  -u(i1-1,i2 ) - u(i1+1,i2 ) );
            end
         end
      end
   end
            
end