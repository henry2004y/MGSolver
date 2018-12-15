classdef (Sealed) Options < handle
   %OPTIONS Multi-level algorithm options.
   % Includes both model parameters and cycle parameters. Sets default
   % values for parameters that can be overriden by the user.
   
   %======================== MEMBERS =================================
   properties (Constant)
      % Model parameters
      domainSize = [1.0 1.0 1.0] % Domain size in all directions
      % Define the initial guess and the RHS functions:
      kx = 5; ky = 5; kz = 5;
      % Right-hand-side
      f = @(x,y,z)...
         (sin(Options.kx*x*pi).*sin(Options.ky*y*pi).*sin(Options.kz*z*pi))
      g = @(x,y,z)(0)
      %g = @(x,y)(cos(2*x+y)+0.5) % Dirichlet boundary condition
      %u = @(x,y,z) (f(x,y,z)/(pi^2*kx^2 + pi^2*ky^2 + pi^2*kz^2));
      % Known solution u = (2*pi^(-2))*sin(pi*x).*sin(pi*y)
      %f = @(x,y)(sin(pi*x).*sin(pi*y)) % Right-hand-side
      %g = @(x,y)(sin(pi*x).*sin(pi*y)) % Dirichlet boundary condition
      % To debug the cycle error, set f=g=0 so that u=error
      %f = @(x,y,z)(zeros(size(x))) % Right-hand-side
      %g = @(x,y,z)(zeros(size(x))) % Dirichlet boundary condition
      
      % Discretization
      nCoarsest = [2 2 2] % #coarsest grid intervals
      numLevels = 6 % #levels. numLevels = finest level
      % Discrete operator at level LEVEL
      operator = @(level)(Operator(level)) 
      
      % Relaxation parameters
      % Gauss-Seidel relaxation (1=LEX, 2=RB)
      smoother = GaussSeidelSmoother(1) 
      % Weighted Jacobi relaxation
      %smoother = JacobiSmoother(6/7)
      % Inter-grid transfers
      % Interpolation of corrections, 1 for line, 2 for 3D matrix
      interpolator = BilinearInterpolator(1)
      % Residual transfer, 1 for line, 2 for 3D matrix
      restrictor = FwLinearRestrictor(1)
      % Cycle parameters
      maxCycleLevels = 100 % # levels to employ in the cycle
      cycleIndex = 1 % V-cycle/W-cycle/etc.
      numCoarsestSweeps = 10 % # relaxation sweeps at coarsest level
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