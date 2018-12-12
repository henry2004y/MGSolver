classdef (Sealed) Cycle < handle
   %CYCLE Multigrid cycle.
   % This class holds the entire multi-level data structure and executes
   % multigrid cycles. A cycling strategy with an integer index is
   % implemented (gamma=1: V-cycle; gamma=2: W-cycle).
   
   %======================== MEMBERS =================================
   properties (GetAccess = private, SetAccess = private)
      levels          % List of levels (1=finest, end=coarsest)
      options         % Contains cycle parameters
      finestRelaxWork % Estimated cost of finest-level relaxation sweep
   end
   %======================== CONSTRUCTORS ============================
   methods
      function obj = Cycle(options, levels)
         % Create a cycle executor with options OPTIONS, to act on the
         % level list LEVELS.
         obj.options = options;
         obj.levels = levels;
         obj.finestRelaxWork = prod(levels{end}.n-1);
      end
   end
   %======================== METHODS =================================
   methods
      function u = cycle(obj, finest, u)
         % The main call that executes a cycle at level FINEST.
         obj.printErrorNormHeader();
         u = obj.cycleAtLevel(finest, finest, u);
      end
   end
   %======================== PRIVATE METHODS =========================
   methods (Access = private)
      function u = cycleAtLevel(obj, l, finest, u)
         % Execute a cycle at level L. FINEST is the index of the finest
         % level in the cycle. Recursively calls itself with the
         % next-coarser level until NUM_LEVELS is reached.
         obj.printErrorNorm(l, 'Initial', u);
         if l == max(1, finest-obj.options.maxCycleLevels+1)
            % Coarsest level
            u = obj.relax(l, obj.options.numCoarsestSweeps, u, false);
         else
            %--- Pre-relaxation ---
            u = obj.relax(l, obj.options.numPreSweeps, u, true);
            %--- Coarse-grid correction ---
            c = l - 1;
            fineLevel = obj.levels{l};
            coarseLevel = obj.levels{c};
            % Transfer fine-level residuals
            r = fineLevel.residual(u);
            coarseLevel.f = fineLevel.restrict(r);
            % Solve residual equation at coarse level Correction
            % scheme: start from vc=0
            vc = zeros(coarseLevel.n);
            for i=1:obj.options.cycleIndex
               vc = obj.cycleAtLevel(c, finest, vc);
            end
            % Interpolate coarse-level correction and add it
            v = fineLevel.interpolate(vc);
            u = u + v;
            obj.printErrorNorm(l, 'Coarse-grid correction', u);
            %--- Post-relaxation ---
            u = obj.relax(l, obj.options.numPostSweeps, u, true);
         end
      end
      
      function u = relax(obj, l, nu, u, printEverySweep)
         % Perform NU relaxation sweeps on U at level LEVEL. If
         % PRINTEVERYSWEEP is true, generates debugging printouts after
         % every sweep; otherwise, only after the last sweep.
         for i=1:nu
            u = obj.levels{l}.relax(u);
            if printEverySweep
               obj.printErrorNorm(l, sprintf('Relaxation sweep %d', i), u);
            end
         end
         if ~printEverySweep
            obj.printErrorNorm(l, sprintf('Relaxation sweep %d', i), u);
         end
      end
      
      function printErrorNormHeader(obj)
         % Print a header line for cycle debugging printouts.
         if obj.options.logLevel >= 1
            fprintf('%-5s %-25s %-13s %-9s\n', 'LEVEL', 'ACTION', ...
               'ERROR NORM', 'WORK');
         end
      end
      
      function u = printErrorNorm(obj, l, action, u)
         % A debugging printout of the error norm at level L after a
         % certain action has been applied. The work per finest-level
         % relaxation sweep is also printed.
         if obj.options.logLevel >= 1
            fprintf('%-5d %-25s %.3e %6.2f\n', l, action, ...
               errornorm(obj.levels{l}, u), flops/obj.finestRelaxWork);
         end
      end
   end
   
end

