classdef TestCycle
   %TESTCYCLE Test the multigrid cycle for the 1D Poisson equation.
   % This class iteratively runs multigrid V-cycles and measures their
   % convergence factor.
   %
   % See also: ERROR_NORM, CYCLE.
   
   %======================== METHODS =================================
   methods(Static)
      function [u, finestLevel] = run
         % Initialize objects
         flops(0); % Reset flop count
         options = Options;
         levels = MultilevelBuilder().build(options);
         cycle = Cycle(options, levels);
         finest = length(levels);
         finestLevel = levels{finest};
         
         % Initial guess
         u = rand(finestLevel.n,1);
         eNew = errornorm(finestLevel, u);
         % Run cycles
         for numCycle=1:options.numCycles
            % Print debugging lines only for the first few cycles
            if numCycle <= 3
               options.logLevel = 1;
               fprintf('############# CYCLE #%d #############\n',numCycle);
            else
               options.logLevel = 0;
            end
            eOld = eNew;
            u = cycle.cycle(finest, u);
            eNew = errornorm(finestLevel, u);
            fprintf('CYCLE %#1d CONVERGENCE FACTOR = %.3f\n', ...
               numCycle, eNew/eOld);
         end
      end
      
      function plot(u)
         % Visualize the final output and compare with analytical solution.
         figure
         options = Options;
         fplot(options.u,[0 options.domainSize])
         hold on
         plot(linspace(0,1.5,numel(u)),u,'--+')
      end
   end
   
end