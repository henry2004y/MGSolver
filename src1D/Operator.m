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
         % 3-point Laplacian with Dirichlet boundary conditions.
         
         % Allocate output array
         result = zeros(obj.level.n,1);
         % Set Dirichlet boundary conditions
         i = [1 obj.level.n]; result(i) = u(i);
         % 3-point discrete Laplacian in the interior domain
         rh2 = 1/obj.level.h^2;
         for i=2:obj.level.n-1
            result(i) = rh2*( 2*u(i) - u(i-1) - u(i+1) );
         end
      end
   end
   
end