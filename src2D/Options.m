classdef (Sealed) Options < handle
   %OPTIONS Multi-level algorithm options.
   % Includes both model parameters and cycle parameters. Sets default
   % values for parameters that can be overriden by the user.
   
   %======================== MEMBERS =================================
   properties (Constant)
      % Model parameters
      domainSize = [2.0 3.0] % Domain size in all directions
      f = @(x,y)(sin(x.^2+y)+0.5) % Right-hand-side
      g = @(x,y)(cos(2*x+y)+0.5) % Dirichlet boundary condition
      % Known solution u = (2*pi^(-2))*sin(pi*x).*sin(pi*y)
      %f = @(x,y)(sin(pi*x).*sin(pi*y)) % Right-hand-side
      %g = @(x,y)(sin(pi*x).*sin(pi*y)) % Dirichlet boundary condition
      % To debug the cycle error, set f=g=0 so that u=error
      %f = @(x,y)(zeros(size(x))) % Right-hand-side
      %g = @(x,y)(zeros(size(x))) % Dirichlet boundary condition
      
      % Discretization
      nCoarsest = [2 3] % #coarsest grid intervals
      numLevels = 6 % #levels. numLevels = finest level
      % Discrete operator at level LEVEL
      operator = @(level)(Operator(level)) 
      
      % Relaxation parameters
      % Gauss-Seidel relaxation (1=LEX, 2=RB)
      smoother = GaussSeidelSmoother(1) 
      % Weighted Jacobi relaxation
      %smoother = JacobiSmoother(2/3)
      % Inter-grid transfers
      interpolator = BilinearInterpolator % Interpolation of corrections
      restrictor = FwLinearRestrictor % Residual transfer
      % Cycle parameters
      maxCycleLevels = 100 % # levels to employ in the cycle
      cycleIndex = 1 % V-cycle/W-cycle/etc.
      numCoarsestSweeps = 5 % # relaxation sweeps at coarsest level
      numPreSweeps = 2 % # pre-CGC relaxation sweeps
      numPostSweeps = 1 % # post-CGC relaxation sweeps
      % Multi-grid run
      numCycles = 1 % #cycles to run
      
      % Plotting
      DoPlot = true
   end
   
   properties
      % Miscellaneous
      logLevel = 1 % Cycle logging level
   end
      
end