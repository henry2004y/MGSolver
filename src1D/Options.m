classdef (Sealed) Options < handle
   %OPTIONS Multi-level algorithm options.
   % Includes both model parameters and cycle parameters. Sets default
   % values for parameters that can be overriden by the user.
   
   %======================== MEMBERS =================================
   properties (Constant) % Model parameters
      domainSize = 1.5 % Domain size in all directions (start from 0 point)
      % Known solution u = sin(k*pi*x/L)
      k = 2
      u = @(x)(sin(Options.k*pi*x/Options.domainSize))
      f = @(x)(Options.k*pi/Options.domainSize)^2*...
         sin(Options.k*pi*x/Options.domainSize) % Right-hand-side
      g = @(x)(0) % Dirichlet boundary condition
      % Known solution u = (2*pi^(-2))*sin(pi*x)
      %f = @(x)(sin(pi*x)) % Right-hand-side
      %g = @(x)(sin(pi*x)) % Dirichlet boundary condition
      % To debug the cycle error, set f=g=0 so that u=error
      %f = @(x)(zeros(size(x))) % Right-hand-side
      %g = @(x)(zeros(size(x))) % Dirichlet boundary condition
      
      % Discretization
      nCoarsest = 4 % #coarsest grid intervals
      numLevels = 3 % #levels. numLevels = finest level
      % Discrete operator at level LEVEL
      operator = @(level)(Operator(level)) 
      
      % Relaxation parameters
      % Gauss-Seidel relaxation (1=LEX, 2=RB)
      %smoother = GaussSeidelSmoother(1)
      % Weighted Jacobi relaxation
      smoother = JacobiSmoother(2/3)
      % Inter-grid transfers
      interpolator = BilinearInterpolator % Interpolation of corrections
      restrictor = FwLinearRestrictor % Residual transfer
      % Cycle parameters
      maxCycleLevels = 100 % # levels to employ in the cycle
      cycleIndex = 1 % V-cycle/W-cycle/etc.
      numCoarsestSweeps = 5 % # relaxation sweeps at coarsest level
      numPreSweeps = 1 % # pre-CGC relaxation sweeps
      numPostSweeps = 1 % # post-CGC relaxation sweeps
      % Multi-grid run
      numCycles = 10 % #cycles to run
      
      % Plotting
      DoPlot = true
   end
   
   properties
      % Miscellaneous
      logLevel = 1 % Cycle logging level
   end
      
end